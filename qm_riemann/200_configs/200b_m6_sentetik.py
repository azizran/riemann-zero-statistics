"""
200b_m6_sentetik.py — KALEM 200-B: Kapı M6b (sentetik güç + SE kalibrasyonu + seçim gücü)
====================================================================================================

Hiçbir gerçek Z/M/δ̃/log(M/δ̃²) verisi kullanılmaz. Üç aşama (KALEM Kapı M6b):

  (i)   Sonlu-ε düzeltme tablosu: CUE Monte Carlo (200t_mc.py MANTIĞININ YENİDEN
        KULLANIMI — Haar QR, 41 noktalı ızgara + parabolik rafinasyon, TAM |Λ_N|
        üzerinden M_n; 200t_mc.py DEĞİŞTİRİLMEDİ, yalnız aynı yöntem burada N∈{9,10,
        11,12} için ölçeklendi), ε∈{0.05,0.1,0.15,0.2,0.3} bantlarında Δκ_r(ε,N) =
        örneklem k_r − kap(r,N,2) (kap: 200a_ortak.kap_vec, analitik kapalı form);
        ayrıca CUE Pearson(δ̃,x) ve ε/2 altındaki oran. Her N'de ε̃<0.2'de ≥4×10⁵ olay.
  (ii)  f-kalibrasyonu: GERÇEK olay sayıları (yalnız sayı — 200b_sayimlar.json,
        pozisyon-türevli, körlük ihlali değil) ile iki gerçekleştirilebilir örnekleyici
        (saf eğik-CUE(b=2) ürün yasası; sağa-çarpık = + bağımsız Euler p=2,3), pencere
        içinde düzgün rastgele t'ye yerleştirilmiş olaylar, 64-bloklu jackknife, 200
        tekrar ⇒ f = SD(tekrarlar)/ortalama(SE_jk) (pencere,ε,kümülant,örnekleyici);
        elementwise max ile birleştirilir.
  (iii) Seçim gücü: her hipotez (SABİT/GEOMETRİK/lin/RMT/Gauss) altında havuzlanmış
        D̄₃ (ve düzeltilmiş havuz κ₃) için Gauss yaklaşıklığıyla karar kurallarının
        (H-200B-1, H-200B-1b, H-200B-2) çıktı olasılıkları ve H-200B-3 sınıf
        karışıklık matrisi; kalibre SE'ler (iki örnekleyici senaryosu) kullanılır.

Çıktı: M6b_200B.json. --stage ile (i/ii/iii/all) ayrı ayrı çalıştırılabilir (her biri
tek başına <10 dk; 'all' toplamda daha uzun sürebilir, ayrı ayrı çalıştırılması
ÖNERİLİR).
"""
import argparse
import importlib.util
import json
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def _load(modname):
    spec = importlib.util.spec_from_file_location(modname, HERE / f"{modname}.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


ortak_b = _load("200b_ortak")
ortak_a = _load("200a_ortak")
orn = _load("200a_orneklem")

N_LIST = (9, 10, 11, 12)
EPS_BANDS = (0.05, 0.1, 0.15, 0.2, 0.3)
EPS_CUSHION = 0.3  # aday üst eşiği (yalnız ihtiyacımız olan en büyük eps'e kadar)
TARGET_EVENTS_AT_02 = 400_000
BATCH = 50_000
GRID_G = 41  # 200t_mc.py ile aynı (41 noktalı, iç 39 nokta + parabol)


# =====================================================================================
# STAGE (i): sonlu-ε CUE Monte Carlo tablosu (200t_mc.py mantığının yeniden kullanımı)
# =====================================================================================
def _mc_batch_events(rng, N, batch, eps_cushion):
    """200t_mc.py'nin AYNI mantığı: Haar QR, N-N kompleks matris; TÜM (n,n+1) ardışık
    çiftler için TAM |Lambda_N| ızgara-maks (41 nokta, parabolik rafinasyon, log
    domaininde — 200t_mc.py'nin kendi yöntemi). su<eps_cushion olan adaylar döner:
    (su, x=log(M/s^2)) burada s HAM (unfold edilmemiş) açısal aralık."""
    Z = (rng.standard_normal((batch, N, N)) + 1j * rng.standard_normal((batch, N, N))) / np.sqrt(2)
    Q, R = np.linalg.qr(Z)
    d = np.diagonal(R, axis1=1, axis2=2)
    U = Q * (d / np.abs(d))[:, None, :]
    th = np.sort(np.angle(np.linalg.eigvals(U)), axis=1)
    nxt = np.concatenate([th[:, 1:], th[:, :1] + 2 * np.pi], axis=1)
    s = nxt - th
    su = s * N / (2 * np.pi)
    b_idx, i_idx = np.nonzero(su < eps_cushion)
    if len(b_idx) == 0:
        return np.empty(0), np.empty(0)
    a = th[b_idx, i_idx]
    ss = s[b_idx, i_idx]
    G = GRID_G
    u = np.linspace(0, 1, G)[1:-1]
    grid = a[:, None] + ss[:, None] * u[None, :]
    ll = np.log(np.abs(2 * np.sin((grid[:, :, None] - th[b_idx][:, None, :]) / 2))).sum(-1)
    k = np.clip(ll.argmax(1), 1, G - 4)
    r = np.arange(len(k))
    y0, y1, y2 = ll[r, k - 1], ll[r, k], ll[r, k + 1]
    den = y0 - 2 * y1 + y2
    logM = y1 - 0.125 * (y2 - y0) ** 2 / np.where(den == 0, -1e-300, den)
    x = logM - 2 * np.log(ss)
    return su[b_idx, i_idx], x


def _stage_i_one_N(args):
    N, seed, target02, batch, eps_cushion, max_batches = args
    rng = np.random.default_rng(seed)
    # eps bandı başına blok-güç-toplamları (jackknife) + korelasyon toplamları
    # satır: [n, Sx1,Sx2,Sx3,Sx4, Ssu1,Ssu2,Ssux]
    block_rows = {e: [] for e in EPS_BANDS}
    n_at_02 = 0
    nmat_total = 0
    nbatches = 0
    t0 = time.time()
    while n_at_02 < target02 and nbatches < max_batches:
        su, x = _mc_batch_events(rng, N, batch, eps_cushion)
        nmat_total += batch
        nbatches += 1
        for e in EPS_BANDS:
            m = su < e
            nn = int(m.sum())
            if nn == 0:
                block_rows[e].append(np.zeros(8))
                continue
            xe, sue = x[m], su[m]
            row = np.array([nn, xe.sum(), (xe ** 2).sum(), (xe ** 3).sum(), (xe ** 4).sum(),
                             sue.sum(), (sue ** 2).sum(), (sue * xe).sum()])
            block_rows[e].append(row)
        n_at_02 += int((su < 0.2).sum())
    dt = time.time() - t0

    out = {"N": N, "nmat": nmat_total, "nbatches": nbatches, "runtime_s": dt, "bands": {}}
    for e in EPS_BANDS:
        blocks = np.array(block_rows[e])  # (B,8)
        xblocks = blocks[:, :5]  # (B,5) n,Sx1..4 -> jackknife_2 uyumlu
        point, cov, reps = ortak_b.jackknife_2(xblocks)
        se_jk = np.sqrt(np.diag(cov))
        k2, k3 = float(point[0]), float(point[1])
        se_k2, se_k3 = float(se_jk[0]), float(se_jk[1])

        total = blocks.sum(axis=0)
        n_tot, Sx1, Sx2, _, _, Ssu1, Ssu2, Ssux = total
        if n_tot > 1:
            mean_x = Sx1 / n_tot
            mean_su = Ssu1 / n_tot
            cov_xsu = Ssux / n_tot - mean_x * mean_su
            var_x = Sx2 / n_tot - mean_x ** 2
            var_su = Ssu2 / n_tot - mean_su ** 2
            pearson = float(cov_xsu / np.sqrt(max(var_x * var_su, 1e-300)))
        else:
            pearson = float("nan")

        n_e = int(n_tot)
        kap2_exact = float(ortak_a.kap_vec(2, N, 2))
        kap3_exact = float(ortak_a.kap_vec(3, N, 2))
        out["bands"][f"{e}"] = {
            "n": n_e, "k2": k2, "k3": k3, "k2_se_jk": se_k2, "k3_se_jk": se_k3,
            "kap2_exact": kap2_exact, "kap3_exact": kap3_exact,
            "delta_k2": k2 - kap2_exact, "delta_k3": k3 - kap3_exact,
            "pearson_su_x": pearson,
        }
    for e in EPS_BANDS:
        n_half = out["bands"].get(f"{e/2}", {}).get("n")
        out["bands"][f"{e}"]["frac_below_half"] = (
            float(n_half / out["bands"][f"{e}"]["n"]) if (n_half is not None and out["bands"][f"{e}"]["n"] > 0)
            else None)
    return out


def stage_i(N_list=N_LIST, target02=TARGET_EVENTS_AT_02, batch=BATCH, n_workers=4,
            max_batches=200, seed0=900000):
    jobs = [(N, seed0 + N, target02, batch, EPS_CUSHION, max_batches) for N in N_list]
    t0 = time.time()
    with Pool(n_workers) as pool:
        results = pool.map(_stage_i_one_N, jobs, chunksize=1)
    dt = time.time() - t0
    print(f"  stage (i): {len(N_list)} N degeri, toplam {dt:.0f} s")
    out = {}
    for r in results:
        out[str(r["N"])] = r
        print(f"    N={r['N']}: nmat={r['nmat']} nbatches={r['nbatches']} rt={r['runtime_s']:.0f}s "
              f"n(eps=0.2)={r['bands']['0.2']['n']} dk2={r['bands']['0.2']['delta_k2']:+.5f} "
              f"dk3={r['bands']['0.2']['delta_k3']:+.5f}")
    return out, dt


# =====================================================================================
# STAGE (ii): f-kalibrasyonu (gerçek olay SAYILARIYLA sentetik x; körlük ihlali değil)
# =====================================================================================
def _job_ii(args):
    w, eps, sampler, Lbar, n_events, seed, nblocks, t_a, t_b = args
    rng = np.random.default_rng(seed)
    if sampler == "cue":
        x = orn.sample_log_tilted_product(Lbar, ortak_b.B_RUNG, n_events, rng)
    else:  # "skew": + bağımsız Euler p=2,3
        x = orn.sample_log_tilted_product(Lbar, ortak_b.B_RUNG, n_events, rng)
        x = x + orn.sample_euler_factor([2, 3], n_events, rng)
    t_pos = rng.uniform(t_a, t_b, n_events)
    blk = ortak_b.block_index(t_pos, t_a, t_b, nblocks)
    blocks = ortak_b.power_sums_blocked(x, blk, nblocks)
    point, cov, _ = ortak_b.jackknife_2(blocks)
    se_jk = np.sqrt(np.diag(cov))
    return w, eps, sampler, point, se_jk


def stage_ii(n_reps=200, nblocks=64, n_workers=8, seed0=910000):
    sayimlar = ortak_b.load_sayimlar()
    real_windows = ortak_b.load_real_windows()
    jobs = []
    for w in ortak_b.WINDOWS:
        t_a, t_b = real_windows[w]["t_a"], real_windows[w]["t_b"]
        for eps in ortak_b.EPSILONS:
            rec = sayimlar[w][f"{eps:.1f}"]
            n_events = int(rec["n"])
            Lbar = float(rec["L_mean"])
            for sampler in ("cue", "skew"):
                for rep in range(n_reps):
                    seed = (seed0 + hash((w, eps, sampler)) % 100000 * 1000 + rep) % (2 ** 31 - 1)
                    jobs.append((w, eps, sampler, Lbar, n_events, seed, nblocks, t_a, t_b))
    t0 = time.time()
    with Pool(n_workers) as pool:
        results = pool.map(_job_ii, jobs, chunksize=8)
    dt = time.time() - t0
    print(f"  stage (ii): {len(jobs)} tekrar (2 ornekleyici x 3 pencere x 3 eps x {n_reps} tekrar), {dt:.0f} s")

    by_key = {}
    for w, eps, sampler, point, se_jk in results:
        key = (w, eps, sampler)
        by_key.setdefault(key, {"points": [], "se_jk": []})
        by_key[key]["points"].append(point)
        by_key[key]["se_jk"].append(se_jk)

    f_detail = {}
    f_combined = {}
    for w in ortak_b.WINDOWS:
        for eps in ortak_b.EPSILONS:
            f_per_sampler = {}
            for sampler in ("cue", "skew"):
                pts = np.array(by_key[(w, eps, sampler)]["points"])
                ses = np.array(by_key[(w, eps, sampler)]["se_jk"])
                sd_over_reps = pts.std(axis=0, ddof=1)
                mean_se = ses.mean(axis=0)
                f2 = sd_over_reps / np.maximum(mean_se, 1e-12)
                f_per_sampler[sampler] = f2
                f_detail[f"{w}_{eps}_{sampler}"] = {"SD_over_reps": sd_over_reps.tolist(),
                                                      "mean_SE_jk": mean_se.tolist(), "f": f2.tolist()}
            fmax = np.maximum(f_per_sampler["cue"], f_per_sampler["skew"])
            f_combined[f"{w}_{eps}"] = {"k2": float(fmax[0]), "k3": float(fmax[1])}
    return f_combined, f_detail, dt


# =====================================================================================
# STAGE (iii): seçim gücü (Gauss yaklaşıklığı) — H-200B-1/1b/2/3
# =====================================================================================
# Hipotez tahminleri (KALEM tablosu; 200b_A_kaymalari.json'dan tam-hassas)
def load_hypotheses():
    with open(HERE / "200b_A_kaymalari.json") as f:
        A = json.load(f)
    A2, A3 = ortak_b.A2_ARITH, ortak_b.A3_ARITH
    hyps = {
        "SABIT": {"D3": A3, "D3_sd": 0.0, "D2": A2, "D2_sd": 0.0},
        "GEOMETRIK": {"D3": A["k3"]["geo"], "D3_sd": A["k3"]["geo_sd"],
                      "D2": A["k2"]["geo"], "D2_sd": A["k2"]["geo_sd"]},
        "LIN": {"D3": A["k3"]["lin"], "D3_sd": A["k3"]["lin_sd"],
                "D2": A["k2"]["lin"], "D2_sd": A["k2"]["lin_sd"]},
        "RMT": {"D3": 0.0, "D3_sd": 0.0, "D2": 0.0, "D2_sd": 0.0},
        "GAUSS": {"D3": None, "D3_sd": 0.0, "D2": None, "D2_sd": 0.0},  # pencereye bağlı, ayrı ele alınır
    }
    return hyps


def pooled_se(se_per_window):
    """3 pencerenin ters-varyans agirlikli havuzlanmis SE'si."""
    se_per_window = np.asarray(se_per_window, dtype=np.float64)
    w = 1.0 / se_per_window ** 2
    return float(1.0 / np.sqrt(w.sum()))


def stage_iii(f_combined, sayimlar, n_draws=20000, seed=920000):
    """Gauss yaklaşıklığı: pencere başına SE_D3(eps=0.2) = f * SE_jk_iid_approx(n) —
    kesin jackknife SE_jk yerine, gerçek n ve tilted-CUE kappa2(L) ile NAIF SE formülü
    kullanılır (M6b(ii)'nin f'i zaten jackknife/naif oranını taşıyor; burada yalnız
    ÖLÇEK için n'ye bağlı bir taban SE gerekir — bu stage yalnız GÜÇ TAHMİNİ, karar
    hesaplarının kendisi 200b_analiz.py'de GERÇEK jackknife ile yapılacak)."""
    import importlib.util
    rng = np.random.default_rng(seed)
    hyps = load_hypotheses()

    se_D3_per_window = {}
    se_D2_per_window = {}
    for w in ortak_b.WINDOWS:
        rec = sayimlar[w]["0.2"]
        n_events = rec["n"]
        Lbar = rec["L_mean"]
        kap2 = float(ortak_a.kap_vec(2, Lbar, ortak_b.B_RUNG))
        kap3 = float(ortak_a.kap_vec(3, Lbar, ortak_b.B_RUNG))
        se_iid_k2 = np.sqrt(2 * kap2 ** 2 / n_events)
        se_iid_k3 = np.sqrt(6 * max(kap2, 1e-9) ** 3 / n_events)
        f = f_combined[f"{w}_0.2"]
        se_D3_per_window[w] = f["k3"] * se_iid_k3
        se_D2_per_window[w] = f["k2"] * se_iid_k2
        hyps.setdefault("_GAUSS_D2_per_window", {})[w] = 0.233653 - kap2  # Gauss D2 hyp (pencereye bağlı)

    SE_D3_pool = pooled_se(list(se_D3_per_window.values()))
    SE_D2_pool = pooled_se(list(se_D2_per_window.values()))

    classes = {
        "SABIT": {"pred": 0.233653, "sigma_h": 0.0},
        "GEOMETRIK": {"pred": hyps["GEOMETRIK"]["D3"], "sigma_h": hyps["GEOMETRIK"]["D3_sd"]},
        "LIN": {"pred": hyps["LIN"]["D3"], "sigma_h": hyps["LIN"]["D3_sd"]},
        "RMT": {"pred": 0.0, "sigma_h": 0.0},
    }

    def classify(d3_draw):
        chi2 = {}
        for cname, c in classes.items():
            chi2[cname] = (d3_draw - c["pred"]) ** 2 / (SE_D3_pool ** 2 + c["sigma_h"] ** 2 + 0.01 ** 2)
        chi2_sonmus = min(chi2["LIN"], chi2["RMT"])
        chi2_group = {"SABIT": chi2["SABIT"], "GEOMETRIK": chi2["GEOMETRIK"], "SONMUS": chi2_sonmus}
        best = min(chi2_group, key=chi2_group.get)
        order = sorted(chi2_group.values())
        delta = order[1] - order[0]
        if chi2_group[best] > 9:
            return "HICBIRI"
        if delta >= 9:
            return best
        return "BELIRSIZ(" + best + ")"

    true_d3 = {"SABIT": 0.233653, "GEOMETRIK": hyps["GEOMETRIK"]["D3"], "LIN": hyps["LIN"]["D3"], "RMT": 0.0}
    confusion = {}
    h1_outcomes = {}
    h1b_outcomes = {}
    for true_name, true_val in true_d3.items():
        draws = rng.normal(true_val, SE_D3_pool, size=n_draws)
        labels = [classify(d) for d in draws]
        uniq, counts = np.unique(labels, return_counts=True)
        confusion[true_name] = {u: float(c / n_draws) for u, c in zip(uniq, counts)}
        p_correct = confusion[true_name].get(true_name if true_name != "LIN" and true_name != "RMT"
                                              else true_name, 0.0)
        h1_outcomes[true_name] = {
            "P_VAR": float(np.mean((draws >= 5 * SE_D3_pool) & (draws >= 0.05))),
            "P_YOK": float(np.mean(draws <= 2 * SE_D3_pool)),
        }
        h1_outcomes[true_name]["P_ZAYIF"] = 1.0 - h1_outcomes[true_name]["P_VAR"] - h1_outcomes[true_name]["P_YOK"]

    low_power = []
    for true_name in true_d3:
        pcorr_class = confusion[true_name].get(
            "SABIT" if true_name == "SABIT" else ("GEOMETRIK" if true_name == "GEOMETRIK" else "SONMUS"), 0.0)
        if pcorr_class < 0.80:
            low_power.append(f"H-200B-3 true={true_name}: P(doğru sınıf)={pcorr_class:.2f}")

    return {
        "SE_D3_pool": SE_D3_pool, "SE_D2_pool": SE_D2_pool,
        "SE_D3_per_window": se_D3_per_window, "SE_D2_per_window": se_D2_per_window,
        "classes": classes, "confusion_H200B3": confusion, "H200B1_outcomes": h1_outcomes,
        "low_power_flagged_lt_0.80": low_power,
        "gauss_D2_per_window": hyps.get("_GAUSS_D2_per_window", {}),
        "note": ("Gauss-yaklasik guc tahmini; SE_D3/D2 tabani NAIF iid formulunden (n gercek "
                 "olay sayisi, kappa2(Lbar) tilted-CUE) f-kalibrasyonuyla olceklenir. GERCEK "
                 "karar SE'leri 200b_analiz.py'de blok-jackknife'tan gelir (M8b), bu yalnizca "
                 "ON GUC tahminidir."),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["i", "ii", "iii", "all"], default="all")
    ap.add_argument("--target02", type=int, default=TARGET_EVENTS_AT_02)
    ap.add_argument("--max-batches", type=int, default=200)
    ap.add_argument("--workers-i", type=int, default=4)
    ap.add_argument("--n-reps-ii", type=int, default=200)
    ap.add_argument("--workers-ii", type=int, default=8)
    ap.add_argument("--n-draws-iii", type=int, default=20000)
    ap.add_argument("--out", default=str(HERE / "M6b_200B.json"))
    args = ap.parse_args()

    out = {}
    if Path(args.out).exists():
        with open(args.out) as f:
            out = json.load(f)

    if args.stage in ("i", "all"):
        print("=== M6b (i): sonlu-ε CUE Monte Carlo tablosu ===")
        res_i, dt_i = stage_i(target02=args.target02, n_workers=args.workers_i,
                               max_batches=args.max_batches)
        out["stage_i_finite_eps_table"] = res_i
        out["stage_i_runtime_s"] = dt_i
        with open(args.out, "w") as f:
            json.dump(out, f, indent=1)
        print(f"  Kaydedildi ({args.out}), stage (i) tamam.")

    if args.stage in ("ii", "all"):
        print("=== M6b (ii): f-kalibrasyonu ===")
        sayimlar = ortak_b.load_sayimlar()
        f_combined, f_detail, dt_ii = stage_ii(n_reps=args.n_reps_ii, n_workers=args.workers_ii)
        out["stage_ii_f_calibration"] = f_combined
        out["stage_ii_f_detail"] = f_detail
        out["stage_ii_runtime_s"] = dt_ii
        with open(args.out, "w") as f:
            json.dump(out, f, indent=1)
        print(f"  f (birlesik, max) ornekler: {json.dumps({k: v for k, v in list(f_combined.items())[:3]})}")
        print(f"  Kaydedildi ({args.out}), stage (ii) tamam.")

    if args.stage in ("iii", "all"):
        print("=== M6b (iii): seçim gücü ===")
        if "stage_ii_f_calibration" not in out:
            raise RuntimeError("stage (iii) icin once stage (ii) calistirilmali (f_calibration gerekli)")
        sayimlar = ortak_b.load_sayimlar()
        res_iii = stage_iii(out["stage_ii_f_calibration"], sayimlar, n_draws=args.n_draws_iii)
        out["stage_iii_power"] = res_iii
        with open(args.out, "w") as f:
            json.dump(out, f, indent=1)
        print(f"  SE_D3_pool={res_iii['SE_D3_pool']:.5f}  SE_D2_pool={res_iii['SE_D2_pool']:.5f}")
        print(f"  P(dogru<0.80) sayisi: {len(res_iii['low_power_flagged_lt_0.80'])}")
        for lp in res_iii['low_power_flagged_lt_0.80']:
            print("   ", lp)
        print(f"  Kaydedildi ({args.out}), stage (iii) tamam.")


if __name__ == "__main__":
    main()
