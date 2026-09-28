"""
200a_m7_karisim.py — KALEM 200-A: Kapı M7 (karışım düzeltmesi)
====================================================================

YALNIZCA MODEL KÜMÜLANT FONKSİYONLARI kullanılır (kap(r,N,b), analitik) — hiçbir
Z/Z' değeri hesaplanmaz. Kullanılan tek "veri", pencere tanımları (t aralığı,
sıfır İNDEKSLERİ ve pozisyonları — bunlar zaten 200a_tahmin.json / zeros6 dosyasında
açık) — bu körlük ihlali DEĞİLDİR (KALEM: "hiçbir istatistik gerçek Z/Z' dağılımından
hesaplanmaz"; burada hiç Z/Z' yok, yalnız sıfırların KONUMLARI, pencere tanımının bir
parçası zaten).

Karışım düzeltmesi = E_{L ~ pencere-içi dağılım}[kappa_r(L)] − kappa_r(L̄)
  rung 0: L ~ t düzgün dağılımı üzerinden L(t)=log(t/2pi) (ince ızgara ile integral)
  rung 1: L ~ pencere içindeki GERÇEK sıfırların L(gamma_n) dağılımı (eşit ağırlık)

M_CUE ve M_ak için hesaplanır (M_ak = kap + sabit A_r; A_r sabiti L'den bağımsız
olduğundan karışım düzeltmesi M_CUE ile TAM AYNI çıkmalı — bu betik ikisini de açıkça
hesaplayıp eşitliği bir iç-tutarlılık kontrolü olarak raporlar).

Çıktı: M7_200A.json. |düzeltme| > 0.002 olan her (pencere,basamak,kümülant,model) flag'lenir.
"""
import importlib.util
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
A2, A3 = -0.088124, 0.233653  # kappa_r^{arit} (KALEM / 200t_aritmetik_kesin.py)
FLAG_THRESH = 0.002


def _load(modname):
    spec = importlib.util.spec_from_file_location(modname, HERE / f"{modname}.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


ortak = _load("200a_ortak")
kap = ortak.kap_vec  # vektörize kap(r,N,b) — bkz 200a_ortak.py (scipy.special, mpmath ile ~1e-14 örtüşür)


def main():
    with open(HERE / "200a_tahmin.json") as f:
        tahmin = json.load(f)
    zeros = np.load(HERE.parent / "128_odl_zeros6_2e6_zeros.npz")["zeros"]
    Lz_all = np.log(zeros / TWO_PI)

    out = {"threshold": FLAG_THRESH, "windows": {}, "flags": []}

    for w in ("W1", "W2", "W3"):
        idx0, idx1 = tahmin[w]["idx"]
        t_a, t_b = tahmin[w]["t"]
        Lbar0 = tahmin[w]["L_rung0"]
        Lbar1 = tahmin[w]["L_rung1"]

        # rung0: t düzgün ızgara (200001 nokta, pencereler_tahmin.py ile aynı çözünürlük)
        t_grid = np.linspace(t_a, t_b, 200001)
        L_grid0 = np.log(t_grid / TWO_PI)

        # rung1: pencere içindeki GERÇEK sıfırların L'si (eşit ağırlık, körlük ihlali değil:
        # yalnız pozisyon, Z/Z' yok)
        L_grid1 = Lz_all[idx0:idx1]

        rec = {"n_t_grid": int(t_grid.size), "n_zeros": int(L_grid1.size),
               "Lbar_rung0": Lbar0, "Lbar_rung1": Lbar1, "models": {}}

        for model, const in (("CUE", {2: 0.0, 3: 0.0}), ("a_k", {2: A2, 3: A3})):
            mrec = {}
            for r in (2, 3):
                mean_k0 = float(np.mean(kap(r, L_grid0, 0))) + const[r]
                k_at_mean0 = float(kap(r, Lbar0, 0)) + const[r]
                corr0 = mean_k0 - k_at_mean0

                mean_k1 = float(np.mean(kap(r, L_grid1, 1))) + const[r]
                k_at_mean1 = float(kap(r, Lbar1, 1)) + const[r]
                corr1 = mean_k1 - k_at_mean1

                mrec[f"kappa{r}"] = {
                    "rung0": {"correction": corr0, "kappa_at_Lbar": k_at_mean0,
                              "flag": abs(corr0) > FLAG_THRESH},
                    "rung1": {"correction": corr1, "kappa_at_Lbar": k_at_mean1,
                              "flag": abs(corr1) > FLAG_THRESH},
                }
                for rung, cc in (("rung0", corr0), ("rung1", corr1)):
                    if abs(cc) > FLAG_THRESH:
                        out["flags"].append(f"{w} {model} kappa{r} {rung}: correction={cc:+.5f}")
            mrec_check = abs(mrec["kappa2"]["rung0"]["correction"]
                              - (float(np.mean(kap(2, L_grid0, 0))) - float(kap(2, Lbar0, 0))))
            mrec["_note"] = "M_ak duzeltmesi M_CUE ile ayni olmali (A_r sabiti L-bagimsiz)"
            out["windows"].setdefault(w, {})
            out["windows"][w][model] = mrec

    with open(HERE / "M7_200A.json", "w") as f:
        json.dump(out, f, indent=1)

    print("=== M7 (karışım düzeltmesi) ===")
    for w in ("W1", "W2", "W3"):
        for model in ("CUE", "a_k"):
            m = out["windows"][w][model]
            print(f"{w} {model}: k2 rung0={m['kappa2']['rung0']['correction']:+.5f} "
                  f"rung1={m['kappa2']['rung1']['correction']:+.5f} | "
                  f"k3 rung0={m['kappa3']['rung0']['correction']:+.5f} "
                  f"rung1={m['kappa3']['rung1']['correction']:+.5f}")
    print(f"\nEşik ({FLAG_THRESH}) aşan düzeltme sayısı: {len(out['flags'])}")
    for fl in out["flags"]:
        print("  FLAG:", fl)
    if not out["flags"]:
        print("  (hiçbiri eşiği aşmadı -> KALEM kuralı: ihmal edilebilir, tahmin.json'daki nokta-L tahminleri kullanılır)")


if __name__ == "__main__":
    main()
