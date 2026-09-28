"""
200c_analiz.py — KALEM 200-C: ANALIZ (200c_olcum.npz -> HUKUM_200C.json)
================================================================================

*** BU BETIK GERCEK VERIDE CALISTIRILMADI *** — yalniz 200c_olcum.npz'nin
SEMASIYLA UYUMLU SENTETIK dosyalarla uctan uca test edilmistir (korluk zaten
200c_olcum.py asamasinda saglandi; girdi yalniz blok GUC TOPLAMLARI (b0/b1) ve
olay-basina KONUM+M dizileri (b2) -- hicbir ornek-basina dagilim ozeti
onceden hesaplanmadi).

Akis: b0/b1 blok-guc-toplamlarindan + b2 olaylarindan (ayni 128 t-blogu ile
binlenerek) ORTAK 3-basamakli (b0,b1,b2) blok-jackknife -> pencere basina 6x6
kovaryans (k2,k3 x b0,b1,b2) -> f-kalibrasyonu (M6c) -> kayma s_b = k3-
kappa3^CUE,b(Lbar) (b2 icin ayrica -Delta-kappa3^CUE(eps=0.2,Lbar), M6c) ->
H-200C-1 (GLS gamma, 3x3 kappa3-alt-blok + sigma_sys^2), H-200C-2 (khi-kare
sinif), H-200C-3 (oranlar), H-200C-4 (7-pencere ortak A/gamma) -> ikincil
(kappa2 kaymalari, b2 bagimsizlik/aralik yasasi).
"""
import argparse
import importlib.util
import json
import sys
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


ortak_a = _load("200a_ortak")
ortak_b = _load("200b_ortak")
kap = ortak_a.kap_vec
kstat = ortak_a.kstat_from_sums
power_sums = ortak_a.power_sums_of

WINDOWS = ("C1", "C2", "C3", "C4")
A2, A3 = ortak_b.A2_ARITH, ortak_b.A3_ARITH
SIGMA_SYS = 0.005
L_NOMINAL = {"C1": "14.194", "C2": "16.577", "C3": "18.884", "C4": "22.306"}

# 200-A/B'den (MUHURLU) ESKI pencere verisi -- 200c_tahmin.py'nin gomulu D
# sozlugunden BIREBIR kopyalandi (kaynak: 200a_merdiven_alt_RAPOR.md /
# 200b_yakin_cift_RAPOR.md tablolari). H-200C-4'un 7-pencereli capraz denetimi
# icin GEREKLI (kalem: "A/B degerleri + yeni 4 pencere").
D_OLD_K3 = {
    0: [(9.313, 0.2983, 0.0181), (10.582, 0.2729, 0.0195), (11.650, 0.2791, 0.0130)],
    1: [(9.343, 0.1439, 0.0027), (10.589, 0.1457, 0.0021), (11.658, 0.1518, 0.0010)],
    2: [(9.355, 0.0764, 0.0059), (10.594, 0.0844, 0.0052), (11.662, 0.0928, 0.0046)],
}
D_OLD_K2 = {
    0: [(9.313, -0.0857, 0.0045), (10.582, -0.0787, 0.0044), (11.650, -0.0818, 0.0039)],
    1: [(9.343, -0.0733, 0.0022), (10.589, -0.0716, 0.0014), (11.658, -0.0731, 0.0011)],
    2: [(9.355, -0.1001, 0.0078), (10.594, -0.1047, 0.0061), (11.662, -0.0938, 0.0052)],
}


# ---------------------------------------------------------------------------
# 3-basamakli (b0,b1,b2) ortak-blok jackknife
# ---------------------------------------------------------------------------
def jackknife_3rung(blocks_win):
    """blocks_win: (B,3,5) [blok, basamak(0,1,2), (n,S1,S2,S3,S4)] -> point (6,)
    [k2_0,k3_0,k2_1,k3_1,k2_2,k3_2], cov (6,6), reps (B,6)."""
    blocks_win = np.asarray(blocks_win, dtype=np.float64)
    B = blocks_win.shape[0]
    total = blocks_win.sum(axis=0)
    point = np.empty(6)
    for ri in range(3):
        n, s1, s2, s3, s4 = total[ri]
        k2, k3 = kstat(n, s1, s2, s3)
        point[2 * ri], point[2 * ri + 1] = k2, k3
    reps = np.empty((B, 6))
    for i in range(B):
        loo = total - blocks_win[i]
        for ri in range(3):
            n, s1, s2, s3, s4 = loo[ri]
            k2, k3 = kstat(n, s1, s2, s3)
            reps[i, 2 * ri], reps[i, 2 * ri + 1] = k2, k3
    rep_bar = reps.mean(axis=0)
    d = reps - rep_bar[None, :]
    cov = (B - 1) / B * (d.T @ d)
    return point, cov, reps


def block_index(t, t_a, t_b, nblocks):
    idx = ((t - t_a) / (t_b - t_a) * nblocks).astype(np.int64)
    return np.clip(idx, 0, nblocks - 1)


def merge_blocks_by_factor(blocks_B_3_5, factor):
    """(B,3,5) -> (B/factor,3,5), bitisik `factor` blogu toplayarak."""
    B = blocks_B_3_5.shape[0]
    assert B % factor == 0, f"blok sayisi {B} {factor}'e bolunmeli"
    return blocks_B_3_5.reshape(B // factor, factor, 3, 5).sum(axis=1)


def m8c_sensitivity(blocks_3_B0):
    """KALEM M8c (200-B kurali): 32/64/128 blok duyarliligi -- yalniz TANI/rapor
    (KALEM'in kendisi 200c_olcum.py icin nblocks=128'i ZATEN SABITLEMIS; burada
    128'in 64'e/32'ye gore yakinsadigini dogrulariz, birincil blok sayisi
    DEGISTIRILMEZ). Doner: {level: SE_k3(3,)} ve 64->128 en buyuk goreli fark."""
    B0 = blocks_3_B0.shape[0]
    se_by_level = {}
    for level in (32, 64, 128):
        if B0 % level != 0:
            continue
        factor = B0 // level
        merged = merge_blocks_by_factor(blocks_3_B0, factor) if factor > 1 else blocks_3_B0
        _, cov6, _ = jackknife_3rung(merged)
        se_by_level[level] = np.sqrt(np.diag(cov6))[[1, 3, 5]]  # k3, basamak 0,1,2
    max_rel_change_64_128 = None
    if 64 in se_by_level and 128 in se_by_level:
        rel = np.abs(se_by_level[128] - se_by_level[64]) / np.maximum(se_by_level[64], 1e-300)
        max_rel_change_64_128 = float(np.max(rel))
    return {lv: se.tolist() for lv, se in se_by_level.items()}, max_rel_change_64_128


def m6c_interp(m6, N, eps, field):
    N_LIST = (14, 17, 19, 22)
    vals = [m6["stage_i_finite_eps_table"][str(Nv)]["bands"][f"{eps}"][field] for Nv in N_LIST]
    return float(np.interp(N, N_LIST, vals))


def h200c2_classify(shifts12, SE12, sigma_h_by_class, preds_by_class):
    chi2 = {}
    for cname, pred in preds_by_class.items():
        sigma_h = sigma_h_by_class[cname]
        var = SE12 ** 2 + sigma_h ** 2 + SIGMA_SYS ** 2
        chi2[cname] = float(np.sum((shifts12 - pred) ** 2 / var))
    order = sorted(chi2.items(), key=lambda kv: kv[1])
    best_name, best_chi2 = order[0]
    delta_next = order[1][1] - best_chi2
    if best_chi2 > 36:
        decision = "HICBIRI"
    elif delta_next >= 9:
        decision = "KESIN:" + best_name
    else:
        decision = "BELIRSIZ(en_iyi=" + best_name + ")"
    return {"chi2": {k: v for k, v in chi2.items()}, "best": best_name, "best_chi2": best_chi2,
            "delta_next": delta_next, "decision": decision}


def gls_fit(L_flat, s_flat, rung_idx, cov_full, free_A=False, a_k=A3):
    Linv = np.linalg.cholesky(np.linalg.inv(cov_full))
    rung_idx = np.asarray(rung_idx)

    if free_A:
        def resid(p):
            gamma, A = p[0], p[1]
            c = p[2:5]
            pred = A + c[rung_idx] * np.asarray(L_flat) ** (-gamma)
            return Linv.T @ (np.asarray(s_flat) - pred)
        x0 = [0.5, a_k, 0.15, -0.2, -0.3]
    else:
        def resid(p):
            gamma = p[0]
            c = p[1:4]
            pred = a_k + c[rung_idx] * np.asarray(L_flat) ** (-gamma)
            return Linv.T @ (np.asarray(s_flat) - pred)
        x0 = [0.5, 0.15, -0.2, -0.3]

    fit = least_squares(resid, x0, max_nfev=5000)
    J = fit.jac
    try:
        cov_p = np.linalg.inv(J.T @ J)
    except np.linalg.LinAlgError:
        cov_p = np.full((len(x0), len(x0)), np.nan)
    chi2 = float(np.sum(fit.fun ** 2))
    return fit.x, cov_p, chi2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--olcum", default=str(HERE / "200c_olcum.npz"))
    ap.add_argument("--m6c", default=str(HERE / "M6c_200C.json"))
    ap.add_argument("--tahmin", default=str(HERE / "200c_tahmin.json"))
    ap.add_argument("--out", default=str(HERE / "HUKUM_200C.json"))
    ap.add_argument("--label", default="")
    ap.add_argument("--windows-meta", default=None,
                     help="test modunda ozel pencere T0/t_a/t_b JSON'u; gercekte "
                          "200c_ortak.load_window() ile otomatik.")
    args = ap.parse_args()

    d = np.load(args.olcum, allow_pickle=True)
    blok_guc = d["blok_guc"]  # (nwin,2,B0,5)
    win_names = [str(x) for x in d["pencere_adlari"]]
    nwin = blok_guc.shape[0]
    B0 = blok_guc.shape[2]
    events = {k: d[k] for k in ("window", "n", "m", "L", "delta_tilde", "M")}

    with open(args.m6c) as f:
        m6 = json.load(f)
    f_calib = m6.get("stage_ii_f_calibration", {})
    with open(args.tahmin) as f:
        tahmin = json.load(f)

    if args.windows_meta:
        with open(args.windows_meta) as f:
            win_meta = json.load(f)
    else:
        # GERCEK MOD: Lbar0 (t-agirlikli, 200001-nokta izgara) ve Lbar1 (sifir-
        # agirlikli, GERCEK sifir konumlari) M7c/M6c ile AYNI yontemle hesaplanir
        # -- yalniz KONUM, korluk ihlali degil.
        ortak = _load("200c_ortak")
        win_meta = {}
        for w in win_names:
            if w in ortak.WINDOWS:
                win = ortak.load_window(w)
                L_grid0 = np.log((win.T0 + np.linspace(win.t_a_off, win.t_b_off, 200001)) / (2 * np.pi))
                L_grid1 = ortak.L_of_offset(win, win.off)
                win_meta[w] = {"t_a": win.t_a_off, "t_b": win.t_b_off,
                                "Lbar0": float(np.mean(L_grid0)), "Lbar1": float(np.mean(L_grid1))}

    shifts_k3 = np.full(3 * nwin, np.nan)
    shifts_k2 = np.full(3 * nwin, np.nan)
    se_k3 = np.full(3 * nwin, np.nan)
    se_k2 = np.full(3 * nwin, np.nan)
    per_window_cov3 = []  # kappa3 3x3 alt-blok, pencere basina
    per_window_record = {}
    L_by_w = {}
    pearson_by_w = {}
    frac_by_w = {}

    for wi, w in enumerate(win_names):
        meta = win_meta.get(w, {"t_a": float(d["t_kenar"][wi][0]), "t_b": float(d["t_kenar"][wi][-1])})
        t_a, t_b = meta["t_a"], meta["t_b"]
        b2mask = events["window"] == w

        m_ev = events["m"][b2mask]
        L_ev = events["L"][b2mask]
        dt_ev = events["delta_tilde"][b2mask]
        M_ev = events["M"][b2mask]
        mask02 = dt_ev < 0.2
        m2, L2, dt2, M2 = m_ev[mask02], L_ev[mask02], dt_ev[mask02], M_ev[mask02]
        n2 = len(m2)
        if n2 < 8:
            continue
        x2 = np.log(M2 / dt2 ** 2)
        Lbar2 = float(np.mean(L2))

        blk2 = block_index(m2, t_a, t_b, B0)
        b2_blocks = ortak_b.power_sums_blocked(x2, blk2, B0)
        blocks_3 = np.stack([blok_guc[wi, 0], blok_guc[wi, 1], b2_blocks], axis=1)  # (B0,3,5)
        point, cov6, reps = jackknife_3rung(blocks_3)

        # ---- M8c (200-B kurali): 32/64/128 blok duyarliligi, YALNIZ TANI ----
        # (KALEM 200c_olcum.py icin nblocks=128'i zaten SABITLIYOR -- burada
        # 128'in kendisinin 64'e gore yakinsadigini DOGRULARIZ; primary blok
        # sayisi degismez, yalniz buyuk bir sapma flag'lenir.)
        m8c_se_by_level, m8c_max_rel = m8c_sensitivity(blocks_3)
        m8c_flag = bool(m8c_max_rel is not None and m8c_max_rel > 0.20)

        n0_tot = blok_guc[wi, 0, :, 0].sum()
        n1_tot = blok_guc[wi, 1, :, 0].sum()
        # Lbar0/Lbar1: sentetik testte dosyadan gelmeyebilir -> t-agirlikli/sifir-agirlikli
        # yaklasik olarak L(t_a..t_b) araligindan (gercek modda 200c_ortak ile tam
        # hesaplanir; burada YALNIZCA yapı dogrulugu icin kaba bir tahmin yeterli).
        Lbar0 = meta.get("Lbar0", float(np.log(0.5 * (t_a + t_b) / (2 * np.pi) + 1)))
        Lbar1 = meta.get("Lbar1", Lbar0)
        L_by_w[w] = (Lbar0, Lbar1, Lbar2)

        f = np.ones(6)
        for ri in range(3):
            key = f"{w}_b{ri}"
            if key in f_calib:
                f[2 * ri] = f_calib[key]["k2"]
                f[2 * ri + 1] = f_calib[key]["k3"]
        se_jk = np.sqrt(np.diag(cov6))
        se_calib = f * se_jk

        have_m6_table = "stage_i_finite_eps_table" in m6
        dkap2 = m6c_interp(m6, Lbar2, 0.2, "delta_k2") if have_m6_table else 0.0
        dkap3 = m6c_interp(m6, Lbar2, 0.2, "delta_k3") if have_m6_table else 0.0

        base_k2 = [float(kap(2, Lbar0, 0)), float(kap(2, Lbar1, 1)), float(kap(2, Lbar2, 2))]
        base_k3 = [float(kap(3, Lbar0, 0)), float(kap(3, Lbar1, 1)), float(kap(3, Lbar2, 2))]

        s_k2 = [point[0] - base_k2[0], point[2] - base_k2[1], point[4] - base_k2[2] - dkap2]
        s_k3 = [point[1] - base_k3[0], point[3] - base_k3[1], point[5] - base_k3[2] - dkap3]

        shifts_k3[3 * wi:3 * wi + 3] = s_k3
        shifts_k2[3 * wi:3 * wi + 3] = s_k2
        se_k3[3 * wi:3 * wi + 3] = [se_calib[1], se_calib[3], np.sqrt(se_calib[5] ** 2 + dkap3 ** 2)]
        se_k2[3 * wi:3 * wi + 3] = [se_calib[0], se_calib[2], np.sqrt(se_calib[4] ** 2 + dkap2 ** 2)]

        idx3 = [1, 3, 5]
        cov3 = cov6[np.ix_(idx3, idx3)] * np.outer([f[1], f[3], f[5]], [f[1], f[3], f[5]])
        per_window_cov3.append(cov3)

        pearson_real = float(np.corrcoef(dt2, x2)[0, 1]) if n2 > 2 else float("nan")
        frac_real = float(np.mean(dt2 < 0.1))
        pearson_cue = m6c_interp(m6, Lbar2, 0.2, "pearson_su_x") if have_m6_table else float("nan")
        frac_cue = m6c_interp(m6, Lbar2, 0.2, "frac_below_half") if have_m6_table else float("nan")
        pearson_by_w[w] = {"real": pearson_real, "CUE_MC": pearson_cue}
        frac_by_w[w] = {"real": frac_real, "CUE_MC": frac_cue, "U13_theory": 1.0 / 8.0}

        per_window_record[w] = {
            "n_events_eps02": n2, "n0": float(n0_tot), "n1": float(n1_tot),
            "Lbar0": Lbar0, "Lbar1": Lbar1, "Lbar2": Lbar2,
            "point_k2k3_b0b1b2": point.tolist(),
            "shift_k3": s_k3, "shift_k2": s_k2,
            "SE_k3_calib": se_k3[3 * wi:3 * wi + 3].tolist(),
            "SE_k2_calib": se_k2[3 * wi:3 * wi + 3].tolist(),
            "M8c": {"SE_k3_by_level": m8c_se_by_level, "max_rel_change_64_to_128": m8c_max_rel,
                     "flag_gt_20pct": m8c_flag, "nblocks_used": B0},
        }

    valid = ~np.isnan(shifts_k3)
    n_valid_windows = int(valid.sum() // 3)

    # ---- H-200C-1: GLS gamma (yalniz kappa3, 3x3-per-pencere + sigma_sys^2) ----
    cov12 = np.zeros((3 * nwin, 3 * nwin))
    for wi in range(nwin):
        if wi < len(per_window_cov3):
            cov12[3 * wi:3 * wi + 3, 3 * wi:3 * wi + 3] = per_window_cov3[wi]
    cov12 += np.eye(3 * nwin) * SIGMA_SYS ** 2

    L_flat = np.array([L_by_w[w][ri] for w in win_names for ri in range(3)])
    rung_idx = np.array([ri for _ in win_names for ri in range(3)])

    frozen_threshold = m6.get("stage_iii_power", {}).get("frozen_threshold", 3.5)

    h1 = None
    if np.all(valid):
        params, cov_p, chi2_h1 = gls_fit(L_flat, shifts_k3, rung_idx, cov12, free_A=False, a_k=A3)
        gamma_hat = float(params[0])
        sigma_gamma = float(np.sqrt(max(cov_p[0, 0], 0.0))) if np.isfinite(cov_p[0, 0]) else float("nan")
        z = gamma_hat / sigma_gamma if sigma_gamma > 0 else float("nan")
        if np.isfinite(z) and z >= frozen_threshold:
            h1_decision = "SONLU_YUKSEKLIK"
        elif np.isfinite(sigma_gamma) and abs(gamma_hat) <= 2 * sigma_gamma:
            h1_decision = "SABIT_YAPI"
        else:
            h1_decision = "BELIRSIZ"
        h1 = {"gamma_hat": gamma_hat, "sigma_gamma": sigma_gamma, "z": z,
              "c_b": params[1:4].tolist(), "chi2": chi2_h1, "threshold_used": frozen_threshold,
              "decision": h1_decision}

    # ---- H-200C-2: sinif secimi (12 kayma vs H_C/H_Sg/H_S1 tablolari) ----
    h2 = None
    if np.all(valid):
        tab_k3 = tahmin["k3"]
        preds = {}
        sigma_h = {}
        for hyp in ("H_C", "H_Sg", "H_S1"):
            arr = np.empty(3 * nwin)
            sig = np.empty(3 * nwin)
            for wi, w in enumerate(win_names):
                Lnom = L_NOMINAL.get(w)
                for ri in range(3):
                    val, sd = tab_k3["tablo"][str(ri)][Lnom][hyp]
                    arr[3 * wi + ri] = val
                    sig[3 * wi + ri] = sd
            preds[hyp] = arr
            sigma_h[hyp] = sig
        # sigma_h TABLODAKI GERCEK NOKTA-BASINA yayilim (KALEM: "sigma_h: tablodaki
        # tahmin yayilimi") -- skaler ortalama DEGIL, 12-elemanli vektor olarak kullanilir.
        h2 = h200c2_classify(shifts_k3, se_k3, sigma_h, preds)

    # ---- H-200C-3: oranlar r1=s1/s0, r2=s2/s1 pencere basina ----
    h3 = {}
    for wi, w in enumerate(win_names):
        if w not in per_window_record:
            continue
        s0, s1, s2 = per_window_record[w]["shift_k3"]
        se0, se1, se2 = per_window_record[w]["SE_k3_calib"]
        r1 = s1 / s0 if s0 else float("nan")
        r2 = s2 / s1 if s1 else float("nan")
        r1_se = abs(r1) * np.sqrt((se1 / s1) ** 2 + (se0 / s0) ** 2) if s0 and s1 else float("nan")
        r2_se = abs(r2) * np.sqrt((se2 / s2) ** 2 + (se1 / s1) ** 2) if s1 and s2 else float("nan")
        h3[w] = {"r1_s1_s0": r1, "r1_SE": r1_se, "r2_s2_s1": r2, "r2_SE": r2_se}

    # ---- H-200C-4: 7-pencere ortak fit (A serbest, gamma ortak) ----
    h4 = None
    if np.all(valid):
        L7, s7, rung7 = [], [], []
        se7 = []
        for ri in range(3):
            for (Lv, sv, ev) in D_OLD_K3[ri]:
                L7.append(Lv)
                s7.append(sv)
                rung7.append(ri)
                se7.append(ev)
        for wi, w in enumerate(win_names):
            for ri in range(3):
                L7.append(L_by_w[w][ri])
                s7.append(shifts_k3[3 * wi + ri])
                rung7.append(ri)
                se7.append(se_k3[3 * wi + ri])
        L7 = np.array(L7)
        s7 = np.array(s7)
        rung7 = np.array(rung7)
        cov21 = np.diag(np.array(se7) ** 2) + np.eye(len(se7)) * SIGMA_SYS ** 2
        params4, covp4, chi2_h4 = gls_fit(L7, s7, rung7, cov21, free_A=True, a_k=A3)
        gamma4, A4 = float(params4[0]), float(params4[1])
        sigma_gamma4 = float(np.sqrt(max(covp4[0, 0], 0))) if np.isfinite(covp4[0, 0]) else float("nan")
        sigma_A4 = float(np.sqrt(max(covp4[1, 1], 0))) if np.isfinite(covp4[1, 1]) else float("nan")
        tutar = (abs(A4 - A3) <= 2 * sigma_A4) and (gamma4 / sigma_gamma4 >= 3 if sigma_gamma4 > 0 else False)
        conflict = None
        if h1 and h1["decision"] == "SONLU_YUKSEKLIK" and gamma4 <= 2 * sigma_gamma4:
            conflict = "H-200C-1 SONLU-YUKSEKLIK diyor ama H-200C-4 gamma<=2sigma -> CELISKI, manset BELIRSIZ'e iner"
        h4 = {"gamma_hat": gamma4, "sigma_gamma": sigma_gamma4, "A_hat": A4, "sigma_A": sigma_A4,
              "c_b": params4[2:5].tolist(), "chi2": chi2_h4, "n_points": len(s7),
              "TUTAR": bool(tutar), "conflict_note": conflict}

    out = {
        "label": args.label,
        "input": {"olcum": str(args.olcum), "n_blocks_file": int(B0), "n_windows": nwin,
                   "windows": win_names},
        "frozen_threshold_H200C1": frozen_threshold,
        "per_window": per_window_record,
        "H_200C_1": h1,
        "H_200C_2": h2,
        "H_200C_3": h3,
        "H_200C_4": h4,
        "secondary": {
            "shift_k2_per_point": shifts_k2.tolist(),
            "SE_k2_per_point": se_k2.tolist(),
            "pearson_delta_x": pearson_by_w,
            "frac_below_eps_half": frac_by_w,
        },
    }
    with open(args.out, "w") as f:
        json.dump(out, f, indent=1, default=lambda o: None)

    print(f"=== 200c_analiz.py sonucu ({'GERCEK' if args.label=='' and 'C1' in win_names else args.label}) ===")
    if h1:
        print(f"H-200C-1: gamma_hat={h1['gamma_hat']:.3f} +/- {h1['sigma_gamma']:.3f} "
              f"(z={h1['z']:.2f}, esik={h1['threshold_used']}) -> {h1['decision']}")
    if h2:
        print(f"H-200C-2: en iyi={h2['best']} chi2={h2['best_chi2']:.2f} delta_next={h2['delta_next']:.2f} "
              f"-> {h2['decision']}")
    if h4:
        print(f"H-200C-4: A_hat={h4['A_hat']:+.4f}+/-{h4['sigma_A']:.4f} (a_k={A3}) "
              f"gamma_hat={h4['gamma_hat']:.3f}+/-{h4['sigma_gamma']:.3f} -> TUTAR={h4['TUTAR']}")
        if h4["conflict_note"]:
            print("  CELISKI:", h4["conflict_note"])
    for w, r in h3.items():
        print(f"H-200C-3 {w}: r1(s1/s0)={r['r1_s1_s0']:.3f}+/-{r['r1_SE']:.3f}  "
              f"r2(s2/s1)={r['r2_s2_s1']:.3f}+/-{r['r2_SE']:.3f}")
    print(f"Kaydedildi: {args.out}")


if __name__ == "__main__":
    main()
