"""
200b_analiz.py — KALEM 200-B: ANALİZ (olay dosyasından karar) → HUKUM_200B.json
====================================================================================

*** BU BETİK GERÇEK VERİ ÜZERİNDE ÇALIŞTIRILMADI *** — yalnız 200b_olaylar.npz'nin
ŞEMASIYLA UYUMLU SENTETİK dosyalar üzerinde uçtan uca test edilmiştir (körlük zaten
200b_olcum.py aşamasında sağlanmıştı; bu betik onu bozamaz çünkü girdi zaten olay
başına (pencere,n,m,L,δ̃,M) — bu betiğin KENDİSİ çalıştırıldığında k-istatistiği/
korelasyon/χ² hesaplar, ki bu KALEM'in ta kendisidir; makine-inşa ajanı bu betiği
gerçek veriyle ÇALIŞTIRMAZ, yalnız inşa edip sentetik veriyle test eder).

Akış (KALEM "Karar kuralları" + Görev 2): ε∈{0.1,0.2,0.3}, pencere başına x=log(M/δ̃²)
-> yansız k2,k3 -> M8b blok-jackknife SE (32/64/128, %20 kuralı) -> f-kalibrasyonu
(M6b) -> baseline κ_r^{CUE,2}(L̄_W) -> sonlu-ε düzeltmesi Δκ_r(ε,L̄_W) (M6b tablosu,
N'de doğrusal ara değer) -> D_r=k_r−baseline−Δκ_r, SE_D²=(f·SE_jk)²+Δκ_r² -> havuz
D̄_r (ters-varyans) -> H-200B-1, H-200B-1b (ε=0.2, birincil), H-200B-2, H-200B-3 ->
ikincil (D̄₂ uzaklıkları, Pearson, δ̃<ε/2 oranı, tüm ε'ler).
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


ortak_b = _load("200b_ortak")
ortak_a = _load("200a_ortak")
kap = ortak_a.kap_vec

N_GRID_M6 = [9.0, 10.0, 11.0, 12.0]
A2, A3 = ortak_b.A2_ARITH, ortak_b.A3_ARITH  # -0.088124, +0.233653
GAUSS_KAPPA2, GAUSS_KAPPA3 = 0.233653, -0.1036  # R_Gauss (L-bagimsiz, TEOREM §8)


def load_hyps():
    with open(HERE / "200b_A_kaymalari.json") as f:
        return json.load(f)


def load_tol_200b():
    with open(HERE / "HUKUM_200A.json") as f:
        h = json.load(f)
    return h["tolerance_200B"]


def m6b_table(m6, N, eps, field):
    vals = [m6["stage_i_finite_eps_table"][str(int(Nv))]["bands"][f"{eps}"][field] for Nv in N_GRID_M6]
    return float(np.interp(N, N_GRID_M6, vals))


def per_window_analysis(events, w, eps, t_a, t_b, m6, f_calib):
    mask = (events["window"] == w) & (events["delta_tilde"] < eps)
    n_events = int(mask.sum())
    if n_events < 8:
        return None  # cok az olay (sentetik test icin guvenlik; gercek veri hicbir zaman bu kadar az degil)

    m_ev = events["m"][mask]
    L_ev = events["L"][mask]
    dt_ev = events["delta_tilde"][mask]
    M_ev = events["M"][mask]
    x = np.log(M_ev / dt_ev ** 2)

    Lbar_w = float(np.mean(L_ev))
    n_tot, s1, s2, s3, _ = ortak_a.power_sums_of(x)
    k2, k3 = ortak_a.kstat_from_sums(n_tot, s1, s2, s3)
    k2, k3 = float(k2), float(k3)

    baseline2 = float(kap(2, Lbar_w, ortak_b.B_RUNG))
    baseline3 = float(kap(3, Lbar_w, ortak_b.B_RUNG))
    dkap2 = m6b_table(m6, Lbar_w, eps, "delta_k2")
    dkap3 = m6b_table(m6, Lbar_w, eps, "delta_k3")
    pearson_cue = m6b_table(m6, Lbar_w, eps, "pearson_su_x")
    frac_cue = m6b_table(m6, Lbar_w, eps, "frac_below_half")

    D2 = k2 - baseline2 - dkap2
    D3 = k3 - baseline3 - dkap3

    # ---- M8b: 32/64/128 blok jackknife, %20 kurali ----
    se_by_level = {}
    for nb in (32, 64, 128):
        blk = ortak_b.block_index(m_ev, t_a, t_b, nb)
        blocks = ortak_b.power_sums_blocked(x, blk, nb)
        _, cov, _ = ortak_b.jackknife_2(blocks)
        se_by_level[nb] = np.sqrt(np.diag(cov))
    se32, se64, se128 = se_by_level[32], se_by_level[64], se_by_level[128]
    rel_change = np.abs(se32 - se64) / np.maximum(se64, 1e-300)
    m8b_trigger = rel_change > 0.20
    se_jk_final = np.where(m8b_trigger, np.maximum(se32, se64), se64)

    f_k2 = f_calib[f"{w}_{eps}"]["k2"]
    f_k3 = f_calib[f"{w}_{eps}"]["k3"]
    SE_D2 = float(np.sqrt((f_k2 * se_jk_final[0]) ** 2 + dkap2 ** 2))
    SE_D3 = float(np.sqrt((f_k3 * se_jk_final[1]) ** 2 + dkap3 ** 2))

    pearson_real = float(np.corrcoef(dt_ev, x)[0, 1])
    frac_real = float(np.mean(dt_ev < eps / 2))

    return {
        "n_events": n_events, "Lbar_w": Lbar_w, "k2": k2, "k3": k3,
        "baseline2_CUE": baseline2, "baseline3_CUE": baseline3,
        "delta_kappa2_finite_eps": dkap2, "delta_kappa3_finite_eps": dkap3,
        "D2": D2, "D3": D3, "SE_D2": SE_D2, "SE_D3": SE_D3,
        "M8b": {"SE32": se32.tolist(), "SE64": se64.tolist(), "SE128": se128.tolist(),
                "rel_change_32_64": rel_change.tolist(), "trigger": m8b_trigger.tolist(),
                "SE_jk_final": se_jk_final.tolist()},
        "f_calibration": {"k2": f_k2, "k3": f_k3},
        "pearson_delta_x_real": pearson_real, "pearson_delta_x_CUE_MC": pearson_cue,
        "frac_below_half_real": frac_real, "frac_below_half_CUE_MC": frac_cue,
        "frac_below_half_U13_theory": 1.0 / 8.0,
        "k3_minus_finite_eps_corr": k3 - dkap3,  # H-200B-1b icin
    }


def pooled(values, ses):
    values = np.asarray(values, dtype=np.float64)
    ses = np.asarray(ses, dtype=np.float64)
    w = 1.0 / ses ** 2
    val = float(np.sum(w * values) / np.sum(w))
    se = float(1.0 / np.sqrt(np.sum(w)))
    return val, se


def h200b3_classify(D3bar, SE, hyps):
    classes = {
        "SABIT": {"pred": A3, "sigma_h": 0.0},
        "GEOMETRIK": {"pred": hyps["k3"]["geo"], "sigma_h": hyps["k3"]["geo_sd"]},
    }
    lin = {"pred": hyps["k3"]["lin"], "sigma_h": hyps["k3"]["lin_sd"]}
    rmt = {"pred": 0.0, "sigma_h": 0.0}

    def chi2(pred, sigma_h):
        return (D3bar - pred) ** 2 / (SE ** 2 + sigma_h ** 2 + 0.01 ** 2)

    chi2_sabit = chi2(**classes["SABIT"])
    chi2_geo = chi2(**classes["GEOMETRIK"])
    chi2_lin = chi2(**lin)
    chi2_rmt = chi2(**rmt)
    chi2_sonmus = min(chi2_lin, chi2_rmt)

    group = {"SABIT": chi2_sabit, "GEOMETRIK": chi2_geo, "SONMUS": chi2_sonmus}
    order = sorted(group.items(), key=lambda kv: kv[1])
    best_name, best_chi2 = order[0]
    delta_next = order[1][1] - best_chi2

    if best_chi2 > 9:
        decision = "HICBIRI"
    elif delta_next >= 9:
        decision = "KESIN:" + best_name
    else:
        decision = "BELIRSIZ(en_iyi=" + best_name + ")"

    return {"chi2": {"SABIT": chi2_sabit, "GEOMETRIK": chi2_geo, "SONMUS": chi2_sonmus,
                      "SONMUS_lin": chi2_lin, "SONMUS_rmt": chi2_rmt},
            "best_class": best_name, "best_chi2": best_chi2, "delta_next": delta_next,
            "decision": decision}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--olaylar", default=str(HERE / "200b_olaylar.npz"))
    ap.add_argument("--m6b", default=str(HERE / "M6b_200B.json"))
    ap.add_argument("--out", default=str(HERE / "HUKUM_200B.json"))
    ap.add_argument("--label", default="")
    ap.add_argument("--windows-json", default=None,
                     help="test modunda ozel pencere t_a/t_b/isim JSON'u (gercekte kullanilmaz, "
                          "gercek pencereler otomatik 200b_ortak.load_real_windows()'tan gelir)")
    args = ap.parse_args()

    d = np.load(args.olaylar, allow_pickle=True)
    events = {k: d[k] for k in d.files}

    with open(args.m6b) as f:
        m6 = json.load(f)
    f_calib = m6["stage_ii_f_calibration"]
    hyps = load_hyps()
    tol = load_tol_200b()

    if args.windows_json:
        with open(args.windows_json) as f:
            win_bounds = json.load(f)
    else:
        rw = ortak_b.load_real_windows()
        win_bounds = {w: {"t_a": rw[w]["t_a"], "t_b": rw[w]["t_b"]} for w in ortak_b.WINDOWS}

    win_names = sorted(set(events["window"].tolist()))

    all_eps_results = {}
    for eps in ortak_b.EPSILONS:
        per_window = {}
        for w in win_names:
            if w not in win_bounds:
                continue
            r = per_window_analysis(events, w, eps, win_bounds[w]["t_a"], win_bounds[w]["t_b"], m6, f_calib)
            if r is not None:
                per_window[w] = r
        if not per_window:
            continue

        D2bar, SE_D2bar = pooled([r["D2"] for r in per_window.values()],
                                  [r["SE_D2"] for r in per_window.values()])
        D3bar, SE_D3bar = pooled([r["D3"] for r in per_window.values()],
                                  [r["SE_D3"] for r in per_window.values()])
        k3corr_pooled, se_k3corr = pooled([r["k3_minus_finite_eps_corr"] for r in per_window.values()],
                                           [r["SE_D3"] for r in per_window.values()])

        all_eps_results[eps] = {"per_window": per_window, "D2bar": D2bar, "SE_D2bar": SE_D2bar,
                                 "D3bar": D3bar, "SE_D3bar": SE_D3bar,
                                 "k3_corrected_pooled": k3corr_pooled, "SE_k3_corrected_pooled": se_k3corr}

    if ortak_b.EPS_PRIMARY not in all_eps_results:
        raise RuntimeError(f"birincil eps={ortak_b.EPS_PRIMARY} icin yeterli olay yok")
    prim = all_eps_results[ortak_b.EPS_PRIMARY]
    D3bar, SE_D3bar = prim["D3bar"], prim["SE_D3bar"]
    D2bar, SE_D2bar = prim["D2bar"], prim["SE_D2bar"]

    # ---- H-200B-1 ----
    if D3bar >= 5 * SE_D3bar and D3bar >= 0.05:
        h1 = "VAR"
    elif D3bar <= 2 * SE_D3bar:
        h1 = "YOK"
    else:
        h1 = "ZAYIF"

    # ---- H-200B-1b ----
    k3corr, se_k3corr = prim["k3_corrected_pooled"], prim["SE_k3_corrected_pooled"]
    if k3corr > 0 and abs(k3corr) >= 3 * se_k3corr:
        h1b = "DONDU"
    elif k3corr < 0 and abs(k3corr) >= 3 * se_k3corr:
        h1b = "DONMEDI"
    else:
        h1b = "BELIRSIZ"

    # ---- H-200B-2 ----
    k3_tutar = abs(D3bar - A3) <= tol["kappa3"]
    k2_tutar = abs(D2bar - A2) <= tol["kappa2"]
    if k3_tutar and k2_tutar:
        h2 = "TUTAR"
    elif k3_tutar or k2_tutar:
        h2 = "KISMI"
    else:
        h2 = "OLU"

    # ---- H-200B-3 ----
    h3 = h200b3_classify(D3bar, SE_D3bar, hyps)

    # ---- Onceklik/celiski notu ----
    contradiction = None
    if h2 in ("TUTAR", "KISMI") and not h3["decision"].startswith("KESIN:SABIT"):
        contradiction = (f"H-200B-2 a_k={h2} (D3bar={D3bar:.4f} tol={tol['kappa3']:.4f} bandinda) "
                          f"ama H-200B-3 manset sinifi={h3['decision']} (SABIT degil) -> manset H-200B-3'e gore.")

    # ---- ikincil: D2bar uzaklıklari ----
    d2_distances = {
        "a_k": D2bar - A2, "GEOMETRIK": D2bar - hyps["k2"]["geo"], "LIN": D2bar - hyps["k2"]["lin"],
        "RMT": D2bar - 0.0,
    }
    gauss_d2_per_window = {}
    for w, r in prim["per_window"].items():
        gauss_pred_w = GAUSS_KAPPA2 - r["baseline2_CUE"]
        gauss_d2_per_window[w] = {"pred": gauss_pred_w, "D2_w": r["D2"], "dist": r["D2"] - gauss_pred_w}

    out = {
        "label": args.label,
        "primary_eps": ortak_b.EPS_PRIMARY,
        "tolerance_200A_source": tol,
        "by_eps": {str(eps): {
            "D2bar": res["D2bar"], "SE_D2bar": res["SE_D2bar"],
            "D3bar": res["D3bar"], "SE_D3bar": res["SE_D3bar"],
            "k3_corrected_pooled": res["k3_corrected_pooled"],
            "SE_k3_corrected_pooled": res["SE_k3_corrected_pooled"],
            "per_window": res["per_window"],
        } for eps, res in all_eps_results.items()},
        "H_200B_1": {"decision": h1, "D3bar": D3bar, "SE_D3bar": SE_D3bar,
                     "threshold_var": [5 * SE_D3bar, 0.05], "threshold_yok": 2 * SE_D3bar},
        "H_200B_1b": {"decision": h1b, "k3_corrected_pooled": k3corr, "SE": se_k3corr,
                      "n_sigma": abs(k3corr) / se_k3corr if se_k3corr > 0 else None},
        "H_200B_2": {"decision": h2, "k3_tutar": k3_tutar, "k2_tutar": k2_tutar,
                     "D3bar": D3bar, "D2bar": D2bar, "tol_kappa3": tol["kappa3"], "tol_kappa2": tol["kappa2"]},
        "H_200B_3": h3,
        "contradiction_note": contradiction,
        "secondary": {
            "D2bar_distances_to_hypotheses": d2_distances,
            "D2bar_distance_to_Gauss_per_window": gauss_d2_per_window,
            "per_window_D3_trend": {w: r["D3"] for w, r in prim["per_window"].items()},
            "per_window_D2_trend": {w: r["D2"] for w, r in prim["per_window"].items()},
        },
    }
    with open(args.out, "w") as f:
        json.dump(out, f, indent=1)

    print(f"=== 200b_analiz.py sonucu (eps_primary={ortak_b.EPS_PRIMARY}) ===")
    print(f"D3bar={D3bar:+.5f} SE={SE_D3bar:.5f} ({D3bar/SE_D3bar:.1f} sigma) -> H-200B-1: {h1}")
    print(f"H-200B-1b: kappa3_corrected={k3corr:+.5f} SE={se_k3corr:.5f} -> {h1b}")
    print(f"H-200B-2: D3bar dist to a_k={D3bar-A3:+.5f} (tol {tol['kappa3']:.4f}) "
          f"D2bar dist to a_k={D2bar-A2:+.5f} (tol {tol['kappa2']:.4f}) -> {h2}")
    print(f"H-200B-3: best={h3['best_class']} chi2={h3['best_chi2']:.2f} "
          f"delta_next={h3['delta_next']:.2f} -> {h3['decision']}")
    if contradiction:
        print("  CELISKI:", contradiction)
    print(f"Kaydedildi: {args.out}")


if __name__ == "__main__":
    main()
