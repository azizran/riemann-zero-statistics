"""
200c_m7_karisim.py — KALEM 200-C: Kapı M7c (karışım düzeltmesi)
====================================================================

YALNIZCA sıfırların KONUMLARI (pozisyon) ve model kümülant fonksiyonları
(kap(r,N,b) analitik, 200a_ortak'tan İÇE AKTARILIR; sonlu-ε tablosu M6c_200C.json'dan
ARA DEĞER ile) kullanılır — hiçbir Z/Z'/M değeri hesaplanmaz.

Karışım düzeltmesi (baseline):
  b0 (rastgele t): E_{t~Uniform[t_a,t_b]}[kappa_r(L(t),0)] - kappa_r(L0bar,0)
  b1 (sıfırlar):   E_{gamma_n in pencere}[kappa_r(L(gamma_n),1)] - kappa_r(L1bar,1)
  b2 (yakın çift, eps sabit): E_{olay, L_n}[kappa_r(L_n,2)] - kappa_r(L2bar,2)

Karışım düzeltmesi (sonlu-eps Δκ, yalnız b2): E_{olay}[Δκ_r^CUE(eps,L_n)] -
Δκ_r^CUE(eps, L2bar)  (M6c'nin ara-deger tablosundan, N=L kullanılarak).

KALEM: C1'de ΔL ≈ 0.072 (diğerlerinde daha küçük) — bu, düzeltmelerin küçük
olması BEKLENDİĞİ anlamına gelir (200a/200b M7 örüntüsüyle tutarlı); > 0.002
olan her (pencere,basamak,kümülant[,eps]) flag'lenir.

Çıktı: M7c_200C.json.
"""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

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
kap = ortak.kap

FLAG_THRESH = 0.002


def m6c_interp(m6, N, eps, field):
    N_LIST = (14, 17, 19, 22)
    vals = [m6["stage_i_finite_eps_table"][str(Nv)]["bands"][f"{eps}"][field] for Nv in N_LIST]
    return float(np.interp(N, N_LIST, vals))


def m6c_interp_vec(m6, N_arr, eps, field):
    """Vektorize surum (m6c_interp'in dizi girdili hali) -- buyuk olay
    dizilerinde (C4'te ~40k) Python-basina-cagri maliyetinden kacinmak icin."""
    N_LIST = (14, 17, 19, 22)
    vals = [m6["stage_i_finite_eps_table"][str(Nv)]["bands"][f"{eps}"][field] for Nv in N_LIST]
    return np.interp(N_arr, N_LIST, vals)


def main():
    m6_path = HERE / "M6c_200C.json"
    have_m6 = m6_path.exists()
    m6 = None
    if have_m6:
        with open(m6_path) as f:
            m6 = json.load(f)
        have_m6 = "stage_i_finite_eps_table" in m6

    out = {"threshold": FLAG_THRESH, "windows": {}, "flags": [], "m6c_available": have_m6}

    for w in ortak.WINDOWS:
        win = ortak.load_window(w)
        t_a, t_b = win.t_a_off, win.t_b_off

        # rung 0: t duzgun izgara (200001 nokta)
        L_grid0 = np.log((win.T0 + np.linspace(t_a, t_b, 200001)) / ortak.TWO_PI)
        Lbar0 = float(np.mean(L_grid0))

        # rung 1: gercek sifirlarin L'si
        L_grid1 = ortak.L_of_offset(win, win.off)
        Lbar1 = float(np.mean(L_grid1))

        rec = {"n_t_grid": int(L_grid0.size), "n_zeros": int(L_grid1.size),
               "Lbar_rung0": Lbar0, "Lbar_rung1": Lbar1, "baseline": {}, "delta_eps": {}}

        for rung, L_grid, Lbar in ((0, L_grid0, Lbar0), (1, L_grid1, Lbar1)):
            for r in (2, 3):
                mean_k = float(np.mean(kap(r, L_grid, rung)))
                k_at_mean = float(kap(r, Lbar, rung))
                corr = mean_k - k_at_mean
                rec["baseline"][f"b{rung}_kappa{r}"] = {"correction": corr, "kappa_at_Lbar": k_at_mean,
                                                          "flag": abs(corr) > FLAG_THRESH}
                if abs(corr) > FLAG_THRESH:
                    out["flags"].append(f"{w} b{rung} kappa{r} baseline: correction={corr:+.5f}")

        # rung 2: eps basina (olay L dagilimi)
        ev = ortak.events_from_window(win, eps_max=0.3)
        for eps in (0.1, 0.2, 0.3):
            mask = ev["delta_tilde"] < eps
            L_ev = ev["L"][mask]
            n_ev = int(len(L_ev))
            Lbar2 = float(np.mean(L_ev)) if n_ev else float("nan")
            eps_rec = {"n_events": n_ev, "Lbar": Lbar2}
            for r in (2, 3):
                mean_k = float(np.mean(kap(r, L_ev, 2))) if n_ev else float("nan")
                k_at_mean = float(kap(r, Lbar2, 2)) if n_ev else float("nan")
                corr = mean_k - k_at_mean if n_ev else float("nan")
                eps_rec[f"baseline_kappa{r}"] = {"correction": corr, "flag": bool(abs(corr) > FLAG_THRESH) if n_ev else False}
                if n_ev and abs(corr) > FLAG_THRESH:
                    out["flags"].append(f"{w} b2 eps={eps} kappa{r} baseline: correction={corr:+.5f}")

                if have_m6:
                    mean_dk = float(np.mean(m6c_interp_vec(m6, L_ev, eps, f"delta_k{r}"))) if n_ev else float("nan")
                    dk_at_mean = m6c_interp(m6, Lbar2, eps, f"delta_k{r}") if n_ev else float("nan")
                    corr_dk = mean_dk - dk_at_mean if n_ev else float("nan")
                    eps_rec[f"delta_kappa{r}"] = {"correction": corr_dk,
                                                   "flag": bool(abs(corr_dk) > FLAG_THRESH) if n_ev else False}
                    if n_ev and abs(corr_dk) > FLAG_THRESH:
                        out["flags"].append(f"{w} b2 eps={eps} DELTA-kappa{r}: correction={corr_dk:+.5f}")
            rec["delta_eps"][f"{eps}"] = eps_rec

        out["windows"][w] = rec

    with open(HERE / "M7c_200C.json", "w") as f:
        json.dump(out, f, indent=1)

    print("=== M7c (karisim duzeltmesi) ===")
    print(f"M6c mevcut (Delta-kappa duzeltmesi hesaplanabilir): {have_m6}")
    for w in ortak.WINDOWS:
        b = out["windows"][w]["baseline"]
        print(f"{w}: b0 k2={b['b0_kappa2']['correction']:+.5f} k3={b['b0_kappa3']['correction']:+.5f} | "
              f"b1 k2={b['b1_kappa2']['correction']:+.5f} k3={b['b1_kappa3']['correction']:+.5f}")
    print(f"\nEsik ({FLAG_THRESH}) asan duzeltme sayisi: {len(out['flags'])}")
    for fl in out["flags"]:
        print("  FLAG:", fl)
    if not out["flags"]:
        print("  (hicbiri esigi asmadi -> M7c: hicbir duzeltme uygulanmaz)")


if __name__ == "__main__":
    main()
