"""
200b_m7_karisim.py — KALEM 200-B: Kapı M7b (karışım düzeltmesi)
====================================================================

YALNIZCA GERÇEK SIFIRLARIN KONUMLARI (KALEM Görülmüşlük beyanı: "Sayımlar yalnız
sıfır konumlarından") ve model kümülant fonksiyonu kap(r,N,b=2) (analitik,
200a_ortak.kap_vec, 200-A'dan İÇE AKTARILIR, DEĞİŞTİRİLMEZ) kullanılır — hiçbir
Z/M/log(M/δ̃²) değeri HESAPLANMAZ. Bu, 200a_m7_karisim.py'nin b=2 / olay-tabanlı
sürümüdür (aynı körlük-güvenli gerekçe: 200a_m7_karisim.py §1 docstring'ine bkz.).

Karışım düzeltmesi = E_{L ~ pencere-içi GERÇEK olayların L_n dağılımı}[kappa_r(L,2)]
                      − kappa_r(L̄_W, 2)
  (pencere w, eşik ε için; olay = δ̃_n < ε, L_n = log(m_n/2π), m_n=(γ_n+γ_{n+1})/2 —
  hepsi yalnız KONUM.)

Çıktı: M7b_200B.json. |düzeltme| > 0.002 olan her (pencere,ε,kümülant) flag'lenir
(KALEM Kapı M7b eşiği).
"""
import importlib.util
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
FLAG_THRESH = 0.002


def _load(modname):
    spec = importlib.util.spec_from_file_location(modname, HERE / f"{modname}.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


ortak_b = _load("200b_ortak")
ortak_a = _load("200a_ortak")
kap = ortak_a.kap_vec  # kap(r, N, b) — vektörize, b=2 burada kullanılacak


def main():
    real_windows = ortak_b.load_real_windows()
    zeros = np.load(HERE.parent / "128_odl_zeros6_2e6_zeros.npz")["zeros"]
    sayimlar = ortak_b.load_sayimlar()

    out = {"threshold": FLAG_THRESH, "windows": {}, "flags": []}

    for w in ortak_b.WINDOWS:
        i0, i1 = real_windows[w]["idx"]
        zero_t = zeros[i0:i1]
        out["windows"][w] = {}
        for eps in ortak_b.EPSILONS:
            ev = ortak_b.events_from_zeros(zero_t, i0, eps_max=eps)
            L_events = ev["L"]
            n_events = int(len(L_events))
            Lbar_direct = float(np.mean(L_events)) if n_events else float("nan")
            Lbar_sayimlar = sayimlar[w][f"{eps:.1f}"]["L_mean"]

            rec = {"n_events": n_events, "Lbar_direct": Lbar_direct,
                   "Lbar_sayimlar_200b": Lbar_sayimlar, "kappa": {}}
            for r in (2, 3):
                mean_kap = float(np.mean(kap(r, L_events, ortak_b.B_RUNG))) if n_events else float("nan")
                kap_at_mean = float(kap(r, Lbar_direct, ortak_b.B_RUNG)) if n_events else float("nan")
                correction = mean_kap - kap_at_mean
                flag = bool(abs(correction) > FLAG_THRESH)
                rec["kappa"][f"kappa{r}"] = {"correction": correction, "kappa_at_Lbar": kap_at_mean,
                                              "mean_kappa_over_events": mean_kap, "flag": flag}
                if flag:
                    out["flags"].append(f"{w} eps={eps} kappa{r}: correction={correction:+.5f}")
            out["windows"][w][f"eps_{eps}"] = rec

    with open(HERE / "M7b_200B.json", "w") as f:
        json.dump(out, f, indent=1)

    print("=== M7b (karışım düzeltmesi, b=2, yalnız konum) ===")
    for w in ortak_b.WINDOWS:
        for eps in ortak_b.EPSILONS:
            rec = out["windows"][w][f"eps_{eps}"]
            k2c = rec["kappa"]["kappa2"]["correction"]
            k3c = rec["kappa"]["kappa3"]["correction"]
            print(f"{w} eps={eps}: n={rec['n_events']} Lbar={rec['Lbar_direct']:.4f} "
                  f"k2 duzeltme={k2c:+.5f} k3 duzeltme={k3c:+.5f}")
    print(f"\nEşik ({FLAG_THRESH}) aşan düzeltme sayısı: {len(out['flags'])}")
    for fl in out["flags"]:
        print("  FLAG:", fl)
    if not out["flags"]:
        print("  (hiçbiri eşiği aşmadı -> M7b: hiçbir düzeltme uygulanmaz, nokta-L baseline "
              "kap(r,L̄_W,2) direkt kullanılır)")


if __name__ == "__main__":
    main()
