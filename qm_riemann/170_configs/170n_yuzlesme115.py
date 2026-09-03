"""
170 — EK KAPI (K2′) YÜZLEŞME: ön kayıt ↔ L115'in ölçülen c'si
==============================================================
Ölçüm parçası kopyalanmaz: `170g_yuzlesme.olc` (o da `169_k1.band_jk/
saglikli` → `166_T1.bant_agg/_jk`) aynen import edilir.
Ön kayıt DOSYADAN okunur (`ONKAYIT_L115.json`, 21:16:11); ölçüm
`170k_olcum115.py` 21:16:21'de başlatıldı.

Çıktı: scratchpad/170/K2p_yuzlesme.json
"""
import importlib
import json
import sys
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
for _p in ("170_configs", "169_configs", "167_configs", "166_configs",
           "165_configs", "163_configs", "160_configs", "159_configs",
           "155_configs", "154_configs"):
    sys.path.insert(0, str(QM / _p))
YZ = importlib.import_module("170g_yuzlesme")

SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
C_HIP = 4 / np.pi ** 2
ET = {"Q1": "H-D1 (T-yasası, oranlı)          [İLK-İLKE]",
      "Q2": "H-C1  c = 4/π²                   [İLK-İLKE]",
      "Q3": "derin limite ALTTAN yaklaşma      [İLK-İLKE]",
      "Q4": "λ-kuadratik (4 nokta)             [EMPİRİK]",
      "Q5": "W_X^ν üyesi (ν = 0.830)           [EMPİRİK]",
      "Q6": "KALİB λ-değişmez (168 §A3)        [ölmüştü]"}


def main():
    OK = json.load(open(SCR / "170/ONKAYIT_L115.json"))
    M = YZ.olc("L115")
    H = YZ.olc("Hkeskin")
    print("=" * 96)
    print("K2′ — EK KAPI, ÖN KAYITLI SINAV: λ = 1.15 (L115)")
    print("=" * 96)
    print(f"  ÖN KAYIT: ONKAYIT_L115.json (21:16:11)   "
          f"ÖLÇÜM: C_L115.json (21:16:21 başladı, 4.2 dk)")
    print(f"  marjinaller: σ_ds={OK['sigds']:.5f} σ_X̃={OK['sigX']:.5f} "
          f"σ_Ĉ={OK['sigC']:.5f}")
    print("\n  bant  τ_eff    KALİB_u2      W_X      c_WX ± jk")
    for b in M["bant"]:
        print(f"  {b['lo']:.2f}  {b['tau_eff']:.4f}  {b['KALIB']:.4f}     "
              f"{b['WX']:.4f}   {b['c']:.4f} ± {b['s']:.4f}")
    print(f"\n  **c(L115) = {M['c']:.4f} ± {M['s_jk']:.4f}(jk) "
          f"± {M['s_bant']:.4f}(bant) = ±{M['s_tot']:.4f}**  [{M['nb']} bant]")
    print(f"  çapa c(Hkeskin) = {H['c']:.4f} ± {H['s_tot']:.4f}")

    print("\n  " + "-" * 92)
    print(f"  {'ön-kayıtlı öngörü':46s} {'değer':>8s} {'ölçüm−öngörü':>13s} "
          f"{'σ_tot':>10s}")
    print("  " + "-" * 92)
    res = {}
    for k in ("Q1", "Q2", "Q3", "Q4", "Q5", "Q6"):
        p = OK["ongoru"][k]
        z = (M["c"] - p) / M["s_tot"]
        res[k] = dict(ongoru=p, yuzde=100 * (M["c"] / p - 1), z=z)
        ek = ""
        if k == "Q3":
            ek = ("  [aralık (0.4035, 0.4053]: ÖLÇÜM ARALIĞIN "
                  + ("İÇİNDE]" if 0.4035 < M["c"] <= C_HIP else "ÜSTÜNDE]"))
        print(f"  {ET[k]:46s} {p:8.4f} {100*(M['c']/p-1):+12.2f}% "
              f"{z:+9.2f}σ{ek}")
    print("  " + "-" * 92)

    # KALİB λ-değişmezliği (Q6'nın içeriği) doğrudan
    print("\n  KALİB_u2'nin λ=1.00 → 1.15 değişimi (Q6'nın içeriği):")
    hb = {b["lo"]: b for b in H["bant"]}
    for b in M["bant"]:
        r = b["KALIB"] / hb[b["lo"]]["KALIB"]
        print(f"    lo={b['lo']:.2f}:  {hb[b['lo']]['KALIB']:.4f} → "
              f"{b['KALIB']:.4f}   ({100*(r-1):+.2f}%)     "
              f"W_X: {hb[b['lo']]['WX']:.4f} → {b['WX']:.4f}   "
              f"({100*(b['WX']/hb[b['lo']]['WX']-1):+.2f}%)")
    rk = float(np.mean([b["KALIB"] / hb[b["lo"]]["KALIB"] for b in M["bant"]]))
    rw = float(np.mean([b["WX"] / hb[b["lo"]]["WX"] for b in M["bant"]]))
    print(f"    ortalama: KALİB {100*(rk-1):+.2f}%   W_X {100*(rw-1):+.2f}%"
          f"   ⇒ ν(1.00→1.15) = {np.log(rk)/np.log(rw):+.4f}")

    json.dump(dict(olcum=M, onkayit=OK["ongoru"], sonuc=res,
                   KALIB_oran=rk, WX_oran=rw,
                   nu=float(np.log(rk) / np.log(rw))),
              open(SCR / "170/K2p_yuzlesme.json", "w"), indent=1, default=float)
    print(f"\n-> {SCR}/170/K2p_yuzlesme.json")


if __name__ == "__main__":
    main()
