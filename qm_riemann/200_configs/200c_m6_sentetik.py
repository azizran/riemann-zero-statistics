"""
200c_m6_sentetik.py — KALEM 200-C: Kapı M6c (sentetik güç + SE kalibrasyonu + eşik donma)
================================================================================================

Hiçbir gerçek Z/Z'/M verisi kullanılmaz (b0/b1/b2 için) — yalnız GERÇEK olay
SAYILARI (pozisyon-türevli: b0=1e6, b1=gerçek/alt-örnek sıfır sayısı, b2=gerçek
δ̃<eps olay sayısı) ve GERÇEK t-aralıkları (pencere sınırları, konum) kullanılır.

Üç aşama:
  (i)   Sonlu-ε CUE Monte Carlo tablosu, N ∈ {14,17,19,22}: 200b_m6_sentetik.py'nin
        AYNI mantığı/kodu YENİDEN KULLANILIR (İÇE AKTARILIR, DEĞİŞTİRİLMEZ) —
        yalnız N_list ve target farklı (KALEM 200-C: N ∈ {14,17,19,22}, ≥4e5
        olay ε̃<0.2'de).
  (ii)  f-kalibrasyonu: GERÇEK per-pencere sayılarla (b0=1e6, b1=200c_olcum'un
        karar verdiği sayı, b2≈olay sayısı @ eps=0.2), İKİ gerçekleştirilebilir
        örnekleyici (saf eğik-CUE(N≈L̄,b) ürün yasası; sağa-çarpık = +bağımsız
        Euler p=2,3), 64-blok jackknife, 200 tekrar ⇒ f (pencere,basamak,kümülant).
  (iii) H-200C-1 gücü: kalibre SE + σ_sys ile 12 kaymayı (4 pencere x 3 basamak,
        yalnız κ3) H_C / H_Sγ / H_S1 altında simüle et (200c_tahmin.json
        tablolarından), GLS ile γ̂ uydur, eşik ∈{3.0,3.5,4.0,5.0} için H_C altında
        yanlış-pozitif oranı ve H_Sγ altında güç raporla; KALEM kuralıyla eşiği
        DONDUR (en küçük eşik: FP<%1 (H_C) VE güç>=0.80 (H_Sγ); yoksa 3.5 tut).
        Ayrıca H-200C-2 (sınıf seçimi) karışıklık matrisi.

Çıktı: M6c_200C.json.
"""
import argparse
import importlib.util
import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))


def _load(modname):
    if modname in sys.modules:
        return sys.modules[modname]
    spec = importlib.util.spec_from_file_location(modname, HERE / f"{modname}.py")
    m = importlib.util.module_from_spec(spec)
    sys.modules[modname] = m
    spec.loader.exec_module(m)
    return m


ortak = _load("200c_ortak")
ortak_a = _load("200a_ortak")
ortak_b = _load("200b_ortak")
m6b_mod = _load("200b_m6_sentetik")   # REUSABLE — stage_i DEGISTIRILMEDI, yalniz cagirilir
orn = _load("200a_orneklem")
kap = ortak_a.kap_vec

N_LIST_200C = (14, 17, 19, 22)
TARGET02_200C = 400_000
A2, A3 = ortak_b.A2_ARITH, ortak_b.A3_ARITH
THRESH_CANDIDATES = (3.0, 3.5, 4.0, 5.0)
SIGMA_SYS = 0.005


# =====================================================================================
# STAGE (i): 200b_m6_sentetik.stage_i'nin AYNEN yeniden kullanımı (N_list farkli)
# =====================================================================================
def stage_i(n_workers=8, target02=TARGET02_200C, max_batches=400, batch=None):
    batch_use = batch or m6b_mod.BATCH
    return m6b_mod.stage_i(N_list=N_LIST_200C, target02=target02, batch=batch_use,
                            n_workers=n_workers, max_batches=max_batches, seed0=950000)


def m6c_table(m6_stage_i, N, eps, field):
    vals = [m6_stage_i[str(int(Nv))]["bands"][f"{eps}"][field] for Nv in N_LIST_200C]
    return float(np.interp(N, N_LIST_200C, vals))


# =====================================================================================
# STAGE (ii): f-kalibrasyonu — 3 basamak (b=0,1,2) x 4 pencere, GERCEK n ve t-araligi
# =====================================================================================
def _job_ii(args):
    w, rung, sampler, Lbar, n_events, seed, nblocks, t_a, t_b = args
    rng = np.random.default_rng(seed)
    if sampler == "cue":
        x = orn.sample_log_tilted_product(Lbar, rung, n_events, rng)
    else:  # "skew": +bagimsiz Euler p=2,3
        x = orn.sample_log_tilted_product(Lbar, rung, n_events, rng)
        x = x + orn.sample_euler_factor([2, 3], n_events, rng)
    t_pos = rng.uniform(t_a, t_b, n_events)
    blk = ortak_b.block_index(t_pos, t_a, t_b, nblocks)
    blocks = ortak_b.power_sums_blocked(x, blk, nblocks)
    point, cov, _ = ortak_b.jackknife_2(blocks)
    se_jk = np.sqrt(np.diag(cov))
    return w, rung, sampler, point, se_jk


def real_counts_and_L(b1_decision):
    """Pozisyon-turevli GERCEK sayimlar + Lbar (rung0: t-agirlikli, rung1/2: olay-
    agirlikli) — hicbir Z/Z'/M degeri hesaplanmaz, yalniz sifir KONUMLARI (M0c/M7c
    ile ayni korluk kurali). b1_decision: "all"/"subsample" (tum pencerelere) YA DA
    200c_olcum.B1_DECISION gibi bir dict (pencere->"all"/"subsample", TUTARLILIK
    icin tercih edilir)."""
    olcum = _load("200c_olcum")
    out = {}
    for w in ortak.WINDOWS:
        win = ortak.load_window(w)
        t_a, t_b = win.t_a_off, win.t_b_off
        L_grid = np.log(np.linspace(win.T0 + t_a, win.T0 + t_b, 200001) / ortak.TWO_PI)
        Lbar0 = float(np.mean(L_grid))
        L1_all = ortak.L_of_offset(win, win.off)
        Lbar1 = float(np.mean(L1_all))
        this_decision = b1_decision[w] if isinstance(b1_decision, dict) else b1_decision
        if this_decision == "all":
            assert this_decision == olcum.B1_DECISION[w], (
                f"{w}: real_counts_and_L 'all' diyor ama 200c_olcum.B1_DECISION "
                f"'{olcum.B1_DECISION[w]}' diyor -- f-kalibrasyonu 200c_olcum.py'nin "
                f"GERCEK kararıyla TUTARSIZ olur")
        n1 = win.n_read if this_decision == "all" else min(1_000_000, win.n_read)

        ev02 = ortak.events_from_window(win, eps_max=0.3)
        rec_eps = {}
        for eps in (0.1, 0.2, 0.3):
            mask = ev02["delta_tilde"] < eps
            n_ev = int(mask.sum())
            Lbar2 = float(np.mean(ev02["L"][mask])) if n_ev else Lbar1
            rec_eps[f"{eps}"] = {"n": n_ev, "Lbar": Lbar2}

        out[w] = {"t_a": win.T0 + t_a, "t_b": win.T0 + t_b, "n0": 1_000_000, "Lbar0": Lbar0,
                  "n1": int(n1), "Lbar1": Lbar1, "b2": rec_eps}
    return out


def stage_ii(counts, n_reps=200, nblocks=64, n_workers=8, seed0=960000):
    jobs = []
    for w in ortak.WINDOWS:
        c = counts[w]
        rung_specs = [(0, c["Lbar0"], c["n0"]), (1, c["Lbar1"], c["n1"]),
                      (2, c["b2"]["0.2"]["Lbar"], c["b2"]["0.2"]["n"])]
        for rung, Lbar, n_events in rung_specs:
            for sampler in ("cue", "skew"):
                for rep in range(n_reps):
                    seed = (seed0 + hash((w, rung, sampler)) % 100000 * 1000 + rep) % (2 ** 31 - 1)
                    jobs.append((w, rung, sampler, Lbar, n_events, seed, nblocks, c["t_a"], c["t_b"]))
    t0 = time.time()
    with Pool(n_workers) as pool:
        results = pool.map(_job_ii, jobs, chunksize=8)
    dt = time.time() - t0
    print(f"  stage (ii): {len(jobs)} tekrar (3 basamak x 2 ornekleyici x 4 pencere x {n_reps} tekrar), {dt:.0f} s")

    by_key = {}
    for w, rung, sampler, point, se_jk in results:
        by_key.setdefault((w, rung, sampler), {"points": [], "se_jk": []})
        by_key[(w, rung, sampler)]["points"].append(point)
        by_key[(w, rung, sampler)]["se_jk"].append(se_jk)

    f_combined = {}
    f_detail = {}
    for w in ortak.WINDOWS:
        for rung in (0, 1, 2):
            f_per_sampler = {}
            for sampler in ("cue", "skew"):
                pts = np.array(by_key[(w, rung, sampler)]["points"])
                ses = np.array(by_key[(w, rung, sampler)]["se_jk"])
                sd = pts.std(axis=0, ddof=1)
                mean_se = ses.mean(axis=0)
                f2 = sd / np.maximum(mean_se, 1e-12)
                f_per_sampler[sampler] = f2
                f_detail[f"{w}_b{rung}_{sampler}"] = {"SD_over_reps": sd.tolist(),
                                                        "mean_SE_jk": mean_se.tolist(), "f": f2.tolist()}
            fmax = np.maximum(f_per_sampler["cue"], f_per_sampler["skew"])
            f_combined[f"{w}_b{rung}"] = {"k2": float(fmax[0]), "k3": float(fmax[1])}
    return f_combined, f_detail, dt


# =====================================================================================
# STAGE (iii): H-200C-1 gucu (GLS gamma uydurma) + esik donma + H-200C-2 karisiklik
# =====================================================================================
def gls_fit_gamma(L_by_b, s_by_b, cov_full, a_k=A3):
    """12-nokta (4 pencere x 3 basamak, STAT sirasi: pencere-major, sonra b=0,1,2)
    s_b(L) = a_k + c_b*L^-gamma ortak gamma, c_b serbest (3 parametre + gamma).
    cov_full: 12x12 (Cholesky ile beyazlatilir). Doner: gamma_hat, sigma_gamma."""
    L_flat = L_by_b   # (12,) pencere-major, b=0,1,2 sirasiyla tekrarli
    s_flat = s_by_b
    Linv = np.linalg.cholesky(np.linalg.inv(cov_full))  # whitening: w = Linv @ resid

    def resid(p):
        gamma = p[0]
        c = p[1:4]
        pred = np.empty(12)
        for i in range(12):
            bi = i % 3
            pred[i] = a_k + c[bi] * L_flat[i] ** (-gamma)
        r = s_flat - pred
        return Linv.T @ r

    fit = least_squares(resid, x0=[0.5, 0.15, -0.2, -0.3], max_nfev=2000)
    J = fit.jac
    try:
        cov_p = np.linalg.inv(J.T @ J)
    except np.linalg.LinAlgError:
        cov_p = np.full((4, 4), np.nan)
    gamma_hat = float(fit.x[0])
    sigma_gamma = float(np.sqrt(max(cov_p[0, 0], 0.0))) if np.isfinite(cov_p[0, 0]) else float("nan")
    return gamma_hat, sigma_gamma, fit.x[1:].tolist()


def build_cov_12(se_by_wb, corr_cross_rung=0.3, sigma_sys=SIGMA_SYS):
    """12x12 kovaryans: pencere basina 3x3 blok (SE_jk*f kosegeni + basit ortak-
    korelasyon off-diagonal, cross-rung block jackknife'in kabaca modellenmesi
    icin varsayilan bir korelasyon katsayisi -- M6c GUC SIMULASYONU icin yeterli,
    GERCEK analiz 200c_analiz.py'de GERCEK ortak-blok jackknife kullanacak) +
    sigma_sys^2 kosegen."""
    C = np.zeros((12, 12))
    for wi in range(4):
        se3 = se_by_wb[wi]  # (3,) k3 SE'leri (b=0,1,2)
        block = np.outer(se3, se3) * corr_cross_rung
        np.fill_diagonal(block, se3 ** 2)
        C[3 * wi:3 * wi + 3, 3 * wi:3 * wi + 3] = block
    C += np.eye(12) * sigma_sys ** 2
    return C


def stage_iii(counts, tahmin, f_combined, n_sims=4000, seed=970000):
    rng = np.random.default_rng(seed)
    L_by_w = {}
    se_by_w = {}
    for w in ortak.WINDOWS:
        c = counts[w]
        Ls = [c["Lbar0"], c["Lbar1"], c["b2"]["0.2"]["Lbar"]]
        L_by_w[w] = Ls
        # naif iid SE tabani (f-kalibreli), b0/b1 icin b=0,1 kap; b2 icin b=2
        ses = []
        for bi, Lb in enumerate(Ls):
            kap2 = float(kap(2, Lb, bi))
            se_iid = np.sqrt(6 * max(kap2, 1e-9) ** 3 / (c["n0"] if bi == 0 else (c["n1"] if bi == 1 else c["b2"]["0.2"]["n"])))
            f = f_combined[f"{w}_b{bi}"]["k3"]
            ses.append(f * se_iid)
        se_by_w[w] = np.array(ses)

    L_flat = np.array([L_by_w[w][bi] for w in ortak.WINDOWS for bi in range(3)])
    se_stack = np.array([se_by_w[w] for w in ortak.WINDOWS])  # (4,3)
    cov12 = build_cov_12(se_stack)

    tab_k3 = tahmin["k3"]
    HC_pred = np.array([tab_k3["H_C"][str(bi)] for w in ortak.WINDOWS for bi in range(3)])
    # H_Sg / H_S1 tablo tahminleri: en yakin L_nominal'a gore (tahmin.json 4 nominal L icin tablolu)
    L_nom_order = ["14.194", "16.577", "18.884", "22.306"]

    def table_pred(hyp_key):
        out = np.empty(12)
        for wi, w in enumerate(ortak.WINDOWS):
            Lnom = L_nom_order[wi]
            for bi in range(3):
                out[3 * wi + bi] = tab_k3["tablo"][str(bi)][Lnom][hyp_key][0]
        return out

    HSg_pred = table_pred("H_Sg")
    HS1_pred = table_pred("H_S1")

    # Sinif karari (200c_analiz.py'nin GERCEK H-200C-1 mantigiyla AYNI): z=gamma_hat/
    # sigma_gamma >= esik -> SONLU_YUKSEKLIK; |gamma_hat|<=2*sigma_gamma -> SABIT;
    # aksi -> BELIRSIZ. M6c(iii) yalniz SONLU_YUKSEKLIK oranini (FP/guc) izler.
    # NOT (verimlilik): gamma_hat/sigma_gamma her cekilis icin YALNIZ BIR KEZ
    # uydurulur (esikten BAGIMSIZ); sinif kararı (classify) her esik icin bu AYNI
    # (gamma_hat,sigma_gamma) uzerinden ucuzca tekrarlanir -- ilk surum 4 esigi
    # AYRI AYRI simule ediyordu (4x gereksiz least_squares maliyeti); bu duzeltme
    # M6c(iii)'nin toplam suresini ~4x azaltir.
    rates_by_scenario = {}
    for true_pred, counter_name in ((HC_pred, "fp"), (HSg_pred, "power"), (HS1_pred, "power_s1")):
        z_vals = np.empty(n_sims)
        gamma_vals = np.empty(n_sims)
        for i in range(n_sims):
            draw = rng.multivariate_normal(true_pred, cov12)
            g_hat, sig_g, _ = gls_fit_gamma(L_flat, draw, cov12, a_k=A3)
            gamma_vals[i] = g_hat
            z_vals[i] = g_hat / sig_g if (np.isfinite(sig_g) and sig_g > 0) else np.nan
        rates_by_scenario[counter_name] = (gamma_vals, z_vals)

    results_by_thresh = {}
    for thresh in THRESH_CANDIDATES:
        rates = {}
        for counter_name in ("fp", "power", "power_s1"):
            z_vals = rates_by_scenario[counter_name][1]
            finite = z_vals[np.isfinite(z_vals)]
            rate = float(np.mean(finite >= thresh)) if len(finite) else float("nan")
            rates[counter_name] = rate
        results_by_thresh[thresh] = {"FP_rate_HC": rates["fp"], "power_HSg": rates["power"],
                                      "power_HS1": rates["power_s1"]}

    # KALEM esik dondurma kurali: en kucuk esik [3,5] icinde FP<0.01 (HC) VE power>=0.80 (HSg); yoksa 3.5
    frozen = None
    for thresh in sorted(THRESH_CANDIDATES):
        r = results_by_thresh[thresh]
        if r["FP_rate_HC"] < 0.01 and r["power_HSg"] >= 0.80:
            frozen = thresh
            break
    if frozen is None:
        frozen = 3.5

    # H-200C-2 karisiklik matrisi (sinif secimi: HC/HSg/HS1 arasinda chi2 en iyisi)
    def chi2_of(pred, draw, cov):
        r = draw - pred
        return float(r @ np.linalg.inv(cov) @ r)

    classes = {"H_C": HC_pred, "H_Sg": HSg_pred, "H_S1": HS1_pred}
    confusion = {}
    for true_name, true_pred in classes.items():
        counts_cls = {k: 0 for k in classes}
        for _ in range(n_sims):
            draw = rng.multivariate_normal(true_pred, cov12)
            chis = {k: chi2_of(v, draw, cov12) for k, v in classes.items()}
            best = min(chis, key=chis.get)
            counts_cls[best] += 1
        confusion[true_name] = {k: v / n_sims for k, v in counts_cls.items()}

    return {
        "n_sims": n_sims, "threshold_candidates": list(THRESH_CANDIDATES),
        "results_by_threshold": {str(k): v for k, v in results_by_thresh.items()},
        "frozen_threshold": frozen,
        "H200C2_confusion": confusion,
        "cov12_diag": np.diag(cov12).tolist(),
        "L_flat": L_flat.tolist(),
        "note": ("GUC SIMULASYONU icin cov12, pencere basina 3x3 blok (basit sabit "
                 "cross-rung korelasyon varsayimi, corr=0.3) + sigma_sys^2 kosegen "
                 "kullanir -- GERCEK H-200C-1 karari 200c_analiz.py'de GERCEK ortak-"
                 "blok jackknife'tan gelen kovaryansla yapilacak (bu yalnizca ON "
                 "esik-donma/guc tahmini)."),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["i", "ii", "iii", "all"], default="all")
    ap.add_argument("--target02", type=int, default=TARGET02_200C)
    ap.add_argument("--workers-i", type=int, default=8)
    ap.add_argument("--max-batches-i", type=int, default=400)
    ap.add_argument("--n-reps-ii", type=int, default=200)
    ap.add_argument("--workers-ii", type=int, default=8)
    ap.add_argument("--n-sims-iii", type=int, default=4000)
    ap.add_argument("--b1-decision", choices=["all", "subsample", "olcum"], default="olcum",
                     help="'olcum' (varsayilan): 200c_olcum.B1_DECISION dict'ini pencere-basina "
                          "kullanir (TUTARLILIK icin dogru secim); 'all'/'subsample' tum "
                          "pencerelere zorlar (yalniz test/duyarlilik).")
    ap.add_argument("--out", default=str(HERE / "M6c_200C.json"))
    args = ap.parse_args()

    out = {}
    if Path(args.out).exists():
        with open(args.out) as f:
            out = json.load(f)

    if args.stage in ("i", "all"):
        print("=== M6c (i): sonlu-eps CUE MC tablosu, N in {14,17,19,22} ===")
        res_i, dt_i = stage_i(n_workers=args.workers_i, target02=args.target02,
                               max_batches=args.max_batches_i)
        out["stage_i_finite_eps_table"] = res_i
        out["stage_i_runtime_s"] = dt_i
        with open(args.out, "w") as f:
            json.dump(out, f, indent=1)
        print(f"  Kaydedildi ({args.out}), stage (i) tamam ({dt_i:.0f}s).")

    if args.stage in ("ii", "all"):
        print("=== M6c (ii): f-kalibrasyonu (3 basamak x 4 pencere) ===")
        b1_arg = _load("200c_olcum").B1_DECISION if args.b1_decision == "olcum" else args.b1_decision
        counts = real_counts_and_L(b1_arg)
        out["real_counts"] = counts
        f_combined, f_detail, dt_ii = stage_ii(counts, n_reps=args.n_reps_ii, n_workers=args.workers_ii)
        out["stage_ii_f_calibration"] = f_combined
        out["stage_ii_f_detail"] = f_detail
        out["stage_ii_runtime_s"] = dt_ii
        with open(args.out, "w") as f:
            json.dump(out, f, indent=1)
        print(f"  Kaydedildi ({args.out}), stage (ii) tamam ({dt_ii:.0f}s).")

    if args.stage in ("iii", "all"):
        print("=== M6c (iii): H-200C-1 gucu + esik donma + H-200C-2 karisiklik ===")
        if "real_counts" not in out or "stage_ii_f_calibration" not in out:
            raise RuntimeError("stage (iii) icin once stage (ii) calistirilmali")
        tahmin = ortak.load_tahmin_200c()
        res_iii = stage_iii(out["real_counts"], tahmin, out["stage_ii_f_calibration"],
                             n_sims=args.n_sims_iii)
        out["stage_iii_power"] = res_iii
        with open(args.out, "w") as f:
            json.dump(out, f, indent=1)
        print(f"  Dondurulan esik: {res_iii['frozen_threshold']}")
        for th, r in res_iii["results_by_threshold"].items():
            print(f"    esik={th}: FP(H_C)={r['FP_rate_HC']:.3f} guc(H_Sg)={r['power_HSg']:.3f} "
                  f"guc(H_S1)={r['power_HS1']:.3f}")
        print(f"  Kaydedildi ({args.out}), stage (iii) tamam.")


if __name__ == "__main__":
    main()
