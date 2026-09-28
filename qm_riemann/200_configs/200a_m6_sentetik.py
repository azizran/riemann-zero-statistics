"""
200a_m6_sentetik.py — KALEM 200-A: Kapı M6 (sentetik güç + SE kalibrasyonu + seçim gücü)
====================================================================================================

Hiçbir gerçek Z/Z' verisi kullanılmaz — hepsi model tanımlarından (F ailesi) çekilen
sentetik örneklerdir (bkz 200a_orneklem.py). Gerçek pencere BOYUTLARI (n0=1e6, gerçek
sıfır sayıları n_zeros[W]) yalnızca örnekleyicinin ne kadar örnek çekeceğini belirler —
bunlar zaten 200a_tahmin.json'da açık pencere metadata'sı, gerçek Z/Z' DEĞİL.

Üç aşama:
  (a) Gerçekleştirilebilir modellerin (CUE, hyb2..11) "tam" (yüksek-N Monte Carlo)
      kümülantları — a_k analitik (200t_tam_yasa.py Bulgu-1: sonlu-L olasılık yasası
      yok), Gauss kapalı-form.
  (b) f-kalibrasyon çarpanları: hyb3 ve CUE örnekleyicileri, pencere/basamak başına
      GERÇEK n (b0: 1e6, b1: gerçek sıfır sayısı) ile 30 tekrar, her biri tam
      blok/jackknife hattından (64 blok, iid sıra) geçirilir; f = SD(tekrarlar) /
      ortalama(SE_jk).
  (c) Seçim gücü: F'nin her modeli için Gauss yaklaşıklığıyla (N(kappa^m, C_synth))
      20 000 çekiliş, tam khi-kare (sigma_teori dahil) ile seçim hattı çalıştırılır;
      iki senaryo (yalnız gürültü / gürültü + model-hatası). Karışıklık matrisi,
      P(M*=doğru), P(doğru in S*) raporlanır.

Çıktı: M6_200A.json.
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


ortak = _load("200a_ortak")
orn = _load("200a_orneklem")

WINDOWS = ortak.WINDOWS
STAT_ORDER = ortak.STAT_ORDER
MODELS = ortak.MODELS
REALIZABLE = ("CUE", "hyb2", "hyb3", "hyb5", "hyb7", "hyb11")  # a_k, Gauss ayrı ele alınır
GAUSS_CLOSED = {0: (np.pi ** 2 / 8, -7 * 1.2020569031595942 / 4),
                1: (np.pi ** 2 / 24, -1.2020569031595942 / 4)}

N0_REAL = 1_000_000


def n_zeros_real(tahmin):
    return {w: tahmin[w]["n_zeros"] for w in WINDOWS}


# ---------------------------------------------------------------------------
# (a) "Tam" (yüksek-N MC) kümülantlar, gerçekleştirilebilir modeller
# ---------------------------------------------------------------------------
def _job_a(args):
    w, b, model, L, n_exact, seed = args
    rng = np.random.default_rng(seed)
    if model == "CUE":
        x = orn.sample_log_tilted_product(L, b, n_exact, rng)
    else:
        X = int(model[3:])
        x, _ = orn.sample_hybrid(X, L, b, n_exact, rng)
    k2, k3 = ortak.kstat_from_sums(*ortak.power_sums_of(x)[:4])
    return w, b, model, float(k2), float(k3)


def stage_a(tahmin, n_exact, n_workers=8):
    jobs = []
    seed0 = 606000
    for wi, w in enumerate(WINDOWS):
        for b in (0, 1):
            L = tahmin[w][f"L_rung{b}"]
            for mi, model in enumerate(REALIZABLE):
                seed = seed0 + wi * 1000 + b * 100 + mi
                jobs.append((w, b, model, L, n_exact, seed))
    t0 = time.time()
    with Pool(n_workers) as pool:
        results = pool.map(_job_a, jobs, chunksize=1)
    dt = time.time() - t0
    print(f"  stage (a): {len(jobs)} MC kümülant işi, {dt:.0f} s (n_exact={n_exact})")

    out = {}
    for w, b, model, k2, k3 in results:
        out.setdefault(w, {}).setdefault(f"b{b}", {})[model] = [k2, k3]
    # Gauss (kapalı form) + a_k (analitik, tahmin.json'dan) ekle
    for w in WINDOWS:
        for b in (0, 1):
            out[w][f"b{b}"]["Gauss"] = list(GAUSS_CLOSED[b])
            out[w][f"b{b}"]["a_k"] = tahmin[w]["pred"][f"b{b}"]["a_k"]
    return out, dt


def kappa_true_vector(model, exact_cum, tahmin):
    """STAT_ORDER sırasıyla 12-vektör: gerçekleştirilebilir/Gauss -> MC/kapalı form;
    a_k -> analitik (tahmin.json, aynı zamanda exact_cum içinde de kopyalanmıştır)."""
    out = np.empty(12)
    for i, (w, b, r) in enumerate(STAT_ORDER):
        pair = exact_cum[w][f"b{b}"][model]
        out[i] = pair[0] if r == 2 else pair[1]
    return out


# ---------------------------------------------------------------------------
# (b) f-kalibrasyon çarpanları
# ---------------------------------------------------------------------------
def _job_b(args):
    w, model, L0, L1, n0, n1, seed, nblocks = args
    rng = np.random.default_rng(seed)
    if model == "CUE":
        x0 = orn.sample_log_tilted_product(L0, 0, n0, rng)
        x1 = orn.sample_log_tilted_product(L1, 1, n1, rng)
    else:
        X = int(model[3:])
        x0, _ = orn.sample_hybrid(X, L0, 0, n0, rng)
        x1, _ = orn.sample_hybrid(X, L1, 1, n1, rng)
    blocks = np.zeros((nblocks, 2, 5))
    for x, rung, n in ((x0, 0, n0), (x1, 1, n1)):
        edges = np.linspace(0, n, nblocks + 1).astype(int)
        for bi in range(nblocks):
            blocks[bi, rung] = ortak.power_sums_of(x[edges[bi]:edges[bi + 1]])
    point, cov, _ = ortak.jackknife_4x4(blocks)
    se_jk = np.sqrt(np.diag(cov))
    return w, model, point, se_jk, cov


def stage_b(tahmin, n_reps=30, nblocks=64, n_workers=8):
    nz = n_zeros_real(tahmin)
    jobs = []
    seed0 = 707000
    for wi, w in enumerate(WINDOWS):
        L0 = tahmin[w]["L_rung0"]
        L1 = tahmin[w]["L_rung1"]
        for model in ("CUE", "hyb3"):
            for rep in range(n_reps):
                seed = seed0 + wi * 100000 + hash(model) % 1000 * 1000 + rep
                seed = seed % (2 ** 31 - 1)
                jobs.append((w, model, L0, L1, N0_REAL, nz[w], seed, nblocks))
    t0 = time.time()
    with Pool(n_workers) as pool:
        results = pool.map(_job_b, jobs, chunksize=1)
    dt = time.time() - t0
    print(f"  stage (b): {len(jobs)} tekrar (30 rep x 2 model x 3 pencere), {dt:.0f} s")

    by_wm = {}
    for w, model, point, se_jk, cov in results:
        by_wm.setdefault((w, model), {"points": [], "se_jk": [], "cov": []})
        by_wm[(w, model)]["points"].append(point)
        by_wm[(w, model)]["se_jk"].append(se_jk)
        by_wm[(w, model)]["cov"].append(cov)

    f_per_model = {"CUE": np.empty(12), "hyb3": np.empty(12)}
    cov_avg_per_model = {"CUE": [None] * 3, "hyb3": [None] * 3}
    detail = {}
    for wi, w in enumerate(WINDOWS):
        for model in ("CUE", "hyb3"):
            pts = np.array(by_wm[(w, model)]["points"])       # (n_reps,4)
            ses = np.array(by_wm[(w, model)]["se_jk"])         # (n_reps,4)
            covs = np.array(by_wm[(w, model)]["cov"])          # (n_reps,4,4)
            sd_over_reps = pts.std(axis=0, ddof=1)
            mean_se_jk = ses.mean(axis=0)
            f4 = sd_over_reps / np.maximum(mean_se_jk, 1e-12)
            f_per_model[model][4 * wi:4 * wi + 4] = f4
            cov_avg_per_model[model][wi] = covs.mean(axis=0)
            detail[f"{w}_{model}"] = {"SD_over_reps": sd_over_reps.tolist(),
                                       "mean_SE_jk": mean_se_jk.tolist(), "f": f4.tolist()}

    f_combined = np.maximum(f_per_model["CUE"], f_per_model["hyb3"])
    cov_avg_combined = [0.5 * (cov_avg_per_model["CUE"][wi] + cov_avg_per_model["hyb3"][wi])
                         for wi in range(3)]
    return f_per_model, f_combined, cov_avg_combined, detail, dt


# ---------------------------------------------------------------------------
# (c) Seçim gücü (Gauss yaklaşıklığı)
# ---------------------------------------------------------------------------
def stage_c(tahmin, exact_cum, f_combined, cov_avg_combined, n_draws=20000, seed=808000):
    C_synth = ortak.build_full_cov(cov_avg_combined, f=f_combined)
    sigma_teori = ortak.sigma_teori_vector()
    C_full = C_synth + np.diag(sigma_teori ** 2)
    Cinv = np.linalg.inv(C_full)

    pred_mat = np.array([ortak.model_prediction_vector(tahmin, m) for m in MODELS])  # (8,12)
    kappa_true = {}
    for m in MODELS:
        if m in REALIZABLE or m == "Gauss":
            kappa_true[m] = kappa_true_vector(m, exact_cum, tahmin) if m in REALIZABLE else pred_mat[MODELS.index("Gauss")]
        else:  # a_k
            kappa_true[m] = pred_mat[MODELS.index("a_k")]

    rng = np.random.default_rng(seed)
    results = {}
    low_power_pairs = []
    for scenario in ("noise_only", "noise_plus_model_error"):
        conf = np.zeros((len(MODELS), len(MODELS)))
        p_correct = {}
        p_in_sstar = {}
        for mi, m in enumerate(MODELS):
            mean = kappa_true[m]
            draws = rng.multivariate_normal(mean, C_synth, size=n_draws)
            if scenario == "noise_plus_model_error":
                draws = draws + rng.multivariate_normal(np.zeros(12), np.diag(sigma_teori ** 2), size=n_draws)
            resid = draws[:, None, :] - pred_mat[None, :, :]
            chi2 = np.einsum("nmi,ij,nmj->nm", resid, Cinv, resid)
            mstar_idx = np.argmin(chi2, axis=1)
            chi2min = chi2.min(axis=1, keepdims=True)
            in_sstar = (chi2 - chi2min) < 9
            for j in range(len(MODELS)):
                conf[mi, j] = np.mean(mstar_idx == j)
            p_correct[m] = float(conf[mi, mi])
            p_in_sstar[m] = float(np.mean(in_sstar[:, mi]))
        results[scenario] = {
            "confusion_matrix": conf.tolist(),
            "models_order": list(MODELS),
            "P_Mstar_eq_true": p_correct,
            "P_true_in_Sstar": p_in_sstar,
        }
        for m in MODELS:
            if p_correct[m] < 0.80:
                low_power_pairs.append(f"{scenario}: {m} (P(doğru)={p_correct[m]:.2f})")
    return results, low_power_pairs, C_synth, C_full


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-exact", type=int, default=5_000_000)
    ap.add_argument("--n-reps", type=int, default=30)
    ap.add_argument("--n-draws", type=int, default=20_000)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--out", default=str(HERE / "M6_200A.json"))
    args = ap.parse_args()

    tahmin = ortak.load_tahmin()

    print("=== M6 (a): gerçekleştirilebilir modellerin 'tam' kümülantları (MC) ===")
    exact_cum, dt_a = stage_a(tahmin, args.n_exact, n_workers=args.workers)

    print("\n=== M6 (b): f-kalibrasyon çarpanları (CUE, hyb3; gerçek n; 30 tekrar; 64 blok) ===")
    f_per_model, f_combined, cov_avg_combined, detail_b, dt_b = stage_b(
        tahmin, n_reps=args.n_reps, nblocks=64, n_workers=args.workers)
    print("  f (CUE):   ", np.round(f_per_model["CUE"], 3).tolist())
    print("  f (hyb3):  ", np.round(f_per_model["hyb3"], 3).tolist())
    print("  f (birlesik, max):", np.round(f_combined, 3).tolist())

    print("\n=== M6 (c): seçim gücü (Gauss yaklaşıklığı, 20000 çekiliş) ===")
    results_c, low_power, C_synth, C_full = stage_c(
        tahmin, exact_cum, f_combined, cov_avg_combined, n_draws=args.n_draws)
    for scenario, r in results_c.items():
        print(f"  [{scenario}] P(M*=doğru):", {m: round(v, 3) for m, v in r["P_Mstar_eq_true"].items()})
        print(f"  [{scenario}] P(doğru in S*):", {m: round(v, 3) for m, v in r["P_true_in_Sstar"].items()})
    print("\n  P(doğru)<0.80 olan (model,senaryo) çiftleri (KALEM kuralı: bunlar S*'ta ÖNCEDEN birleştirilir):")
    for lp in low_power:
        print("   ", lp)
    if not low_power:
        print("    (yok — F'nin her üyesi her iki senaryoda da >=0.80 doğru seçim olasılığına sahip)")

    out = {
        "n_exact_MC": args.n_exact,
        "n_reps_f": args.n_reps,
        "n_draws_selection": args.n_draws,
        "runtime_s": {"stage_a": dt_a, "stage_b": dt_b},
        "exact_cumulants_realizable_and_Gauss_and_ak": exact_cum,
        "f_calibration_CUE": f_per_model["CUE"].tolist(),
        "f_calibration_hyb3": f_per_model["hyb3"].tolist(),
        "f_calibration_combined": f_combined.tolist(),
        "f_detail": detail_b,
        "C_synth_diag": np.diag(C_synth).tolist(),
        "selection_power": results_c,
        "low_power_flagged_lt_0.80": low_power,
        "note": ("C_synth = f-kalibreli, CUE ve hyb3 replikalarının ORTALAMA jackknife "
                 "kovaryansı (pencere-blok-köşegen 12x12); a_k icin kappa^true = analitik "
                 "tahmin (sonlu-L olasilik yasasi yok, KALEM Bulgu-1); Gauss kapali-formdan."),
    }
    with open(args.out, "w") as f:
        json.dump(out, f, indent=1)
    print(f"\nKaydedildi: {args.out}")


if __name__ == "__main__":
    main()
