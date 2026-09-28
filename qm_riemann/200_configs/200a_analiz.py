"""
200a_analiz.py — KALEM 200-A: ANALİZ (blok güç toplamlarından karar) → HUKUM_200A.json
=================================================================================================

*** BU BETİK GERÇEK VERİ ÜZERİNDE ÇALIŞTIRILMADI *** — yalnız 200a_bloklar.npz'nin
ŞEMASIYLA UYUMLU SENTETİK dosyalar üzerinde uçtan uca test edilmiştir (bkz. --bloklar /
--tahmin argümanları ve 200a_test_sentetik_analiz.py). Girdi yalnızca BLOK GÜÇ
TOPLAMLARI (n, Sx, Sx2, Sx3, Sx4) — hiçbir örnek-başına dizi okunmaz/işlenmez (körlük
zaten 200a_olcum.py aşamasında sağlanmıştı; bu betik onu bozamaz çünkü girdi formatı
buna izin vermiyor).

Akış: k2,k3 (yansız) <- güç toplamları -> bir-blok-dışı jackknife (64 blok, KALEM
varsayılanı; M8 32/128 duyarlılığı kontrol edilir) -> pencere başı 4x4, blok-köşegen
12x12 kovaryans -> f-kalibrasyonu (M6'dan) -> + sigma_teori^2 köşegen -> khi-kare (F
modelleri) -> H-200A-1 (Delta1), H-200A-2 (M*/S*), G (uyum kalitesi) -> ikincil
sınavlar (basamak-başı seçim, sürekli-X uyumu, tam-korelasyonlu sigma_teori varyantı,
Var_jk/Var_iid) -> 200-B toleransı.
"""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def _load(modname):
    spec = importlib.util.spec_from_file_location(modname, HERE / f"{modname}.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


ortak = _load("200a_ortak")
WINDOWS = ortak.WINDOWS
STAT_ORDER = ortak.STAT_ORDER
MODELS = ortak.MODELS
ARITH = ortak.ARITHMETIC_MODELS


# ---------------------------------------------------------------------------
# Sürekli-X hibrit tahmini (ikincil sınav): X in [2,13] gerçel, N_X = L/(e^gamma log X)
# ---------------------------------------------------------------------------
def hybrid_pred_continuous(L, b, X):
    NX = L / (np.e ** np.euler_gamma * np.log(X))
    k2 = float(ortak.kap_vec(2, NX, b))
    k3 = float(ortak.kap_vec(3, NX, b))
    for p in ortak.primes_leq(X):
        pk2, pk3 = ortak.prime_cumulant(p)
        k2 += pk2
        k3 += pk3
    return k2, k3, NX


# ---------------------------------------------------------------------------
# f-kalibrasyon vektörü: M6_200A.json'dan (varsa); yoksa 1.0 (kalibrasyonsuz, UYARI)
# ---------------------------------------------------------------------------
def load_f_vector(m6_path):
    m6_path = Path(m6_path)
    if not m6_path.exists():
        return np.ones(12), False
    with open(m6_path) as f:
        m6 = json.load(f)
    f = m6.get("f_calibration_combined")
    if f is None or len(f) != 12:
        return np.ones(12), False
    return np.array(f, dtype=np.float64), True


# ---------------------------------------------------------------------------
# Pencere başına jackknife (bir blok sayısı için)
# ---------------------------------------------------------------------------
def per_window_jk(blok_guc, nblocks):
    """blok_guc: (3,2,B0,5); nblocks<=B0 ve B0 % nblocks==0 olmalı. Döner:
    point (3,4), cov_list (3 adet 4x4), reps (3,nblocks,4)."""
    B0 = blok_guc.shape[2]
    assert B0 % nblocks == 0, f"{B0} bloktan {nblocks}'e birleştirilemez"
    factor = B0 // nblocks
    merged = ortak.merge_adjacent_blocks(blok_guc, factor) if factor > 1 else blok_guc
    points = np.empty((3, 4))
    covs = []
    reps_all = []
    for wi in range(3):
        # merged[wi] şekli (2 basamak, nblocks, 5); jackknife_4x4 (nblocks, 2, 5) bekliyor
        pt, cov, reps = ortak.jackknife_4x4(merged[wi].transpose(1, 0, 2))
        points[wi] = pt
        covs.append(cov)
        reps_all.append(reps)
    return points, covs, reps_all


def flatten_point(points):
    """(3,4) -> 12-vektör STAT_ORDER sırasıyla (zaten aynı sıra: pencere,sonra b0k2,b0k3,b1k2,b1k3)."""
    return points.reshape(-1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bloklar", default=str(HERE / "200a_bloklar.npz"))
    ap.add_argument("--tahmin", default=str(HERE / "200a_tahmin.json"))
    ap.add_argument("--m6", default=str(HERE / "M6_200A.json"))
    ap.add_argument("--m7", default=str(HERE / "M7_200A.json"))
    ap.add_argument("--out", default=str(HERE / "HUKUM_200A.json"))
    ap.add_argument("--primary-nblocks", type=int, default=64)
    ap.add_argument("--label", default="")
    args = ap.parse_args()

    d = np.load(args.bloklar, allow_pickle=True)
    blok_guc = d["blok_guc"]
    win_names = [str(x) for x in d["pencere_adlari"]]
    assert blok_guc.shape[0] == 3 and blok_guc.shape[1] == 2 and blok_guc.shape[3] == 5, (
        f"beklenmeyen blok_guc şekli {blok_guc.shape}")
    B0 = blok_guc.shape[2]

    tahmin = ortak.load_tahmin(args.tahmin)
    f_vec, f_from_m6 = load_f_vector(args.m6)
    if not f_from_m6:
        print("UYARI: M6_200A.json bulunamadı / f_calibration_combined eksik -> f=1 (kalibrasyonsuz) kullanılıyor.")

    # ---- M8: blok sayısı duyarlılığı (32/64/128) ----
    levels = [lv for lv in (32, 64, 128) if B0 % lv == 0 and lv <= B0]
    se_by_level = {}
    pts_by_level = {}
    covs_by_level = {}
    for lv in levels:
        pts, covs, _ = per_window_jk(blok_guc, lv)
        se = np.concatenate([np.sqrt(np.diag(c)) for c in covs])
        se_by_level[lv] = se
        pts_by_level[lv] = pts
        covs_by_level[lv] = covs

    primary = args.primary_nblocks if args.primary_nblocks in levels else max(levels)
    m8_switch = False
    m8_report = {}
    if 64 in se_by_level and 128 in se_by_level:
        rel_change = np.abs(se_by_level[128] - se_by_level[64]) / np.maximum(se_by_level[64], 1e-12)
        m8_switch = bool(np.any(rel_change > 0.20))
        m8_report = {"rel_change_64_to_128": rel_change.tolist(), "max_rel_change": float(rel_change.max()),
                      "switch_to_128": m8_switch}
        if m8_switch:
            primary = 128
    m8_report["levels_available"] = levels
    m8_report["primary_used"] = primary

    points, covs, reps = per_window_jk(blok_guc, primary)
    k_obs = flatten_point(points)          # 12-vektör (gözlenen), STAT_ORDER sırasıyla
    C_jk = ortak.build_full_cov(covs, f=None)      # ham jackknife (kalibrasyonsuz), tanı için
    C_f = ortak.build_full_cov(covs, f=f_vec)      # f-kalibreli
    sigma_teori = ortak.sigma_teori_vector()
    C_full = C_f + np.diag(sigma_teori ** 2)

    # Var_jk / Var_iid tanısı (primary blok sayısında)
    _merged_primary = ortak.merge_adjacent_blocks(blok_guc, B0 // primary) if B0 // primary > 1 else blok_guc
    var_iid_per_window = [ortak.jackknife_var_iid(_merged_primary[wi].transpose(1, 0, 2))
                           for wi in range(3)]
    var_iid = np.concatenate(var_iid_per_window)
    var_jk = np.diag(C_jk)
    var_ratio = var_jk / np.maximum(var_iid, 1e-300)

    # ---- khi-kare, F modelleri ----
    chi2_all = {}
    resid_all = {}
    for m in MODELS:
        pred = ortak.model_prediction_vector(tahmin, m)
        r = k_obs - pred
        chi2_all[m] = ortak.chi2(r, C_full)
        resid_all[m] = r

    chi2_arith = {m: chi2_all[m] for m in ARITH}
    m_arith_best = min(chi2_arith, key=chi2_arith.get)
    delta1 = chi2_all["CUE"] - chi2_arith[m_arith_best]
    if delta1 >= 25:
        h1_decision = "TUTAR"
    elif delta1 < 4:
        h1_decision = "OLU"
    else:
        h1_decision = "BELIRSIZ"
    delta1_gauss = chi2_all["Gauss"] - chi2_arith[m_arith_best]   # KALEM: aynı ölçüt Gauss için (Gauss − en iyi aritmetik)

    order = sorted(MODELS, key=lambda m: chi2_all[m])
    m_star = order[0]
    chi2_min = chi2_all[m_star]
    runner_up = order[1]
    delta_runner = chi2_all[runner_up] - chi2_min
    kesin = delta_runner >= 9
    s_star = [m for m in MODELS if chi2_all[m] - chi2_min < 9]

    G = chi2_min / 12.0
    if G <= 3:
        g_band = "IYI"
    elif G <= 10:
        g_band = "KABA"
    else:
        g_band = "YETERSIZ"

    # ---- ikincil: basamak-başı ayrı seçim ----
    idx_b0 = [i for i, (w, b, r) in enumerate(STAT_ORDER) if b == 0]
    idx_b1 = [i for i, (w, b, r) in enumerate(STAT_ORDER) if b == 1]
    sub_results = {}
    for label, idxset in (("b0", idx_b0), ("b1", idx_b1)):
        Csub = C_full[np.ix_(idxset, idxset)]
        chi2_sub = {}
        for m in MODELS:
            pred = ortak.model_prediction_vector(tahmin, m)
            r = (k_obs - pred)[idxset]
            chi2_sub[m] = ortak.chi2(r, Csub)
        best = min(chi2_sub, key=chi2_sub.get)
        sub_results[label] = {"chi2": chi2_sub, "M_star": best}
    ladder_consistent = sub_results["b0"]["M_star"] == sub_results["b1"]["M_star"]

    # ---- ikincil: sürekli-X uyumu ----
    Xgrid = np.linspace(2.0, 13.0, 1101)
    chi2_of_X = np.empty_like(Xgrid)
    for i, X in enumerate(Xgrid):
        pred = np.empty(12)
        for j, (w, b, r) in enumerate(STAT_ORDER):
            L = tahmin[w][f"L_rung{b}"]
            k2c, k3c, _ = hybrid_pred_continuous(L, b, X)
            pred[j] = k2c if r == 2 else k3c
        resid = k_obs - pred
        chi2_of_X[i] = ortak.chi2(resid, C_full)
    ibest = int(np.argmin(chi2_of_X))
    X_star = float(Xgrid[ibest])
    chi2_X_star = float(chi2_of_X[ibest])

    # ---- ikincil: tam-korelasyonlu sigma_teori varyantı ----
    # yorum: (basamak,kümülant) çifti başına PENCERELER ARASI ortak sistematik (korelasyon=1);
    # farklı (basamak,kümülant) çiftleri bağımsız. bkz KALEM "pencereler arası ortak sistematik".
    Sigma_corr = np.zeros((12, 12))
    for i, (wi, bi, ri) in enumerate(STAT_ORDER):
        for j, (wj, bj, rj) in enumerate(STAT_ORDER):
            if bi == bj and ri == rj:
                Sigma_corr[i, j] = ortak.SIGMA_TEORI[ri] * ortak.SIGMA_TEORI[rj]
    C_full_corr = C_f + Sigma_corr
    chi2_all_corr = {}
    for m in MODELS:
        pred = ortak.model_prediction_vector(tahmin, m)
        r = k_obs - pred
        chi2_all_corr[m] = ortak.chi2(r, C_full_corr)
    order_corr = sorted(MODELS, key=lambda m: chi2_all_corr[m])
    m_star_corr = order_corr[0]

    # ---- 200-B toleransı: kümülant başına ----
    tol = {}
    for r in (2, 3):
        idxset = [i for i, (w, b, rr) in enumerate(STAT_ORDER) if rr == r]
        pred = ortak.model_prediction_vector(tahmin, m_star)
        res_r = (k_obs - pred)[idxset]
        rms = float(np.sqrt(np.mean(res_r ** 2)))
        tol[r] = max(ortak.SIGMA_TEORI[r], 2 * rms)

    out = {
        "label": args.label,
        "input": {"bloklar": str(args.bloklar), "n_blocks_file": int(B0)},
        "f_calibration": {"vector": f_vec.tolist(), "from_M6": f_from_m6},
        "M8": m8_report,
        "diagnostics": {"Var_jk": var_jk.tolist(), "Var_iid_naive": var_iid.tolist(),
                         "Var_jk_over_Var_iid": var_ratio.tolist()},
        "k_obs": k_obs.tolist(),
        "STAT_ORDER": [f"{w}_b{b}_k{r}" for (w, b, r) in STAT_ORDER],
        "chi2": chi2_all,
        "H_200A_1": {"delta1_arith_best": m_arith_best, "delta1": delta1, "decision": h1_decision,
                     "delta1_gauss": delta1_gauss},
        "H_200A_2": {"M_star": m_star, "chi2_min": chi2_min, "runner_up": runner_up,
                     "delta_runner_up": delta_runner, "KESIN": kesin, "S_star": s_star},
        "G": {"value": G, "band": g_band},
        "secondary": {
            "per_rung_selection": sub_results,
            "ladder_consistent_M0_eq_M1": ladder_consistent,
            "continuous_X_fit": {"X_star": X_star, "chi2_X_star": chi2_X_star},
            "fully_correlated_sigma_teori": {"chi2": chi2_all_corr, "M_star": m_star_corr},
        },
        "tolerance_200B": {"kappa2": tol[2], "kappa3": tol[3]},
    }
    with open(args.out, "w") as f:
        json.dump(out, f, indent=1)

    print(f"=== 200a_analiz.py sonucu ({'GERÇEK' if 'W1' in win_names and args.label=='' else args.label}) ===")
    print(f"M8: kullanılan blok sayısı = {primary} (mevcut düzeyler {levels}); "
          f"64->128 geçiş tetiklendi mi: {m8_switch}")
    print(f"f-kalibrasyonu M6'dan mı: {f_from_m6}")
    print("khi-kare (F):", {m: round(v, 2) for m, v in chi2_all.items()})
    print(f"H-200A-1: Delta1={delta1:.2f} ({m_arith_best} en iyi aritmetik) -> {h1_decision}  "
          f"(Gauss icin Delta1_Gauss={delta1_gauss:.2f})")
    print(f"H-200A-2: M*={m_star} chi2_min={chi2_min:.2f}  runner-up={runner_up} "
          f"Delta_runnerup={delta_runner:.2f} -> {'KESIN' if kesin else 'S* bandı: ' + str(s_star)}")
    print(f"G={G:.3f} -> {g_band}")
    print(f"Basamak-başı: M*_0={sub_results['b0']['M_star']} M*_1={sub_results['b1']['M_star']} "
          f"tutarlı={ladder_consistent}")
    print(f"Sürekli-X uyumu: X*={X_star:.3f} chi2={chi2_X_star:.2f}")
    print(f"Tam-korelasyonlu sigma_teori varyantı: M*={m_star_corr}")
    print(f"200-B toleransı: kappa2 tol={tol[2]:.4f}  kappa3 tol={tol[3]:.4f}")
    print(f"Kaydedildi: {args.out}")


if __name__ == "__main__":
    main()
