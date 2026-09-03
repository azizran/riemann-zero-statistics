"""
170 — K2 YÜZLEŞME: ÖN KAYITLI ÖNGÖRÜLER ↔ L060'IN ÖLÇÜLEN c'Sİ
===============================================================
Ölçüm parçası KOPYALANMAZ: `169_k1.band_jk/saglikli` (o da
`166_T1.bant_agg/_jk`) aynen import edilir; ölçülen dosya
`scratchpad/167/C_L060.json` (`170d_olcum60.py` = `167_olcum.kos`).

ÖN KAYIT: `scratchpad/170/ONKAYIT_L060.json`, yazılma saati 20:48:26;
ölçüm 20:48:52'de başlatıldı (3 Eylül 2026). Bu betik ön kaydı DOSYADAN
okur, yeniden hesaplamaz.

Ayrıca λ = 0.85 ve 0.70 için T-yasasının İLK-İLKE değerleri (170b'de,
L060 ölçülmeden ÖNCE hesaplandı) ölçülenlerle yüzleştirilir.

Çıktı: scratchpad/170/K2_yuzlesme.json
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
K1 = importlib.import_module("169_k1")

SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
C_HIP = 4.0 / np.pi ** 2
ETIKET = {
    "P1a": "H-D1 MUTLAK   Π_b T(σ_b)                [İLK-İLKE]",
    "P1b": "H-D1 ORANLI   c(Hk)·ΠT(λ)/ΠT(1.00)      [İLK-İLKE]",
    "P2":  "H-C1          c = 4/π² evrensel         [İLK-İLKE]",
    "P3":  "λ-kuadratik   3 ölçülen λ noktası       [EMPİRİK]",
    "P4":  "W_X^ν üyesi   ν = 0.823                 [EMPİRİK]",
    "P5":  "KALİB λ-değişmez (168 §A3)              [ölmüştü]",
}


def olc(ad):
    d = json.load(open(SCR / f"167/C_{ad}.json"))
    B = K1.saglikli(d)
    rows = [K1.band_jk(b["cizgi"]) for b in B]
    c = np.array([a["c_WX"] for a in rows])
    s = np.array([a["sc_WX"] for a in rows])
    lg = np.log(c)
    cg = float(np.exp(lg.mean()))
    s_jk = float(np.sqrt(np.sum((s / c) ** 2)) / len(c)) * cg
    s_sc = float(np.std(lg, ddof=1) / np.sqrt(len(lg))) * cg
    return dict(c=cg, s_jk=s_jk, s_bant=s_sc,
                s_tot=float(np.hypot(s_jk, s_sc)), nb=len(c),
                bant=[dict(lo=b["lo"], tau_eff=a["tau_eff"], c=a["c_WX"],
                           s=a["sc_WX"], KALIB=a["KALIB_u2"], WX=a["W_X"])
                      for b, a in zip(B, rows)])


def main():
    OK = json.load(open(SCR / "170/ONKAYIT_L060.json"))
    M = olc("L060")
    print("=" * 100)
    print("K2 — ÖN KAYITLI ÖRNEKLEM-DIŞI SINAV:  λ = 0.60  (L060)")
    print("=" * 100)
    print(f"  ÖN KAYIT dosyası : ONKAYIT_L060.json  (20:48:26)")
    print(f"  ÖLÇÜM            : C_L060.json        (20:48:52 başladı, "
          f"3.7 dk)")
    print(f"  gaz marjinalleri : σ_ds={OK['sigds']:.5f}  "
          f"σ_X̃={OK['sigX']:.5f}  σ_Ĉ={OK['sigC']:.5f}")
    print(f"\n  bant  τ_eff    KALİB_u2      W_X      c_WX ± jk")
    for b in M["bant"]:
        print(f"  {b['lo']:.2f}  {b['tau_eff']:.4f}  {b['KALIB']:.4f}     "
              f"{b['WX']:.4f}   {b['c']:.4f} ± {b['s']:.4f}")
    print(f"\n  **c(L060) = {M['c']:.4f} ± {M['s_jk']:.4f}(jk) "
          f"± {M['s_bant']:.4f}(bant) = ±{M['s_tot']:.4f}**  "
          f"[{M['nb']} bant]")

    print("\n  " + "-" * 96)
    print(f"  {'öngörü':52s} {'değer':>8s} {'ölçüm−öngörü':>13s} "
          f"{'σ_tot cinsinden':>16s}")
    print("  " + "-" * 96)
    res = {}
    for k in ("P1a", "P1b", "P2", "P3", "P4", "P5"):
        p = OK["ongoru"][k]
        d = M["c"] - p
        z = d / M["s_tot"]
        res[k] = dict(ongoru=p, fark=d, yuzde=100 * (M["c"] / p - 1), z=z)
        print(f"  {ETIKET[k]:52s} {p:8.4f} {100*(M['c']/p-1):+12.2f}% "
              f"{z:+15.2f}σ")
    print("  " + "-" * 96)
    print(f"  {'4/π²':52s} {C_HIP:8.4f} {100*(M['c']/C_HIP-1):+12.2f}% "
          f"{(M['c']-C_HIP)/M['s_tot']:+15.2f}σ")

    # --- λ dizisi: T-yasasının ilk-ilke değerleri ----------------------
    KT = json.load(open(SCR / "170/K1_T.json"))["gaz"]
    print("\n" + "=" * 100)
    print("K2(b) — T-YASASININ İLK-İLKE DEĞERLERİ ↔ ÖLÇÜLEN c  (λ dizisi)")
    print("=" * 100)
    OLC = {a: olc(a) for a in ("Hkeskin", "L085", "L070")}
    OLC["L060"] = M
    piT0 = KT["Hkeskin"]["PiT"]
    c0 = OLC["Hkeskin"]["c"]
    print(f"  {'gaz':9s} {'λ':>5s} {'σ_Ĉ':>8s} {'Π_b T(σ_b)':>11s} "
          f"{'H-D1 oranlı':>12s} {'ÖLÇÜLEN c':>10s} {'±σ_tot':>7s} "
          f"{'fark':>8s} {'σ':>8s}")
    lam_tab = []
    for ad, lam in (("Hkeskin", 1.00), ("L085", 0.85), ("L070", 0.70),
                    ("L060", 0.60)):
        pT = (KT[ad]["PiT"] if ad in KT
              else float(json.load(open(SCR / "170/ONKAYIT_L060.json"))
                         ["ongoru"]["P1a"]))
        sC = (KT[ad]["sigC"] if ad in KT else OK["sigC"])
        pred = c0 * pT / piT0
        o = OLC[ad]
        z = (o["c"] - pred) / o["s_tot"]
        lam_tab.append(dict(gaz=ad, lam=lam, sigC=sC, PiT=pT, hd1=pred,
                            c=o["c"], s=o["s_tot"], z=z))
        print(f"  {ad:9s} {lam:5.2f} {sC:8.5f} {pT:11.5f} {pred:12.4f} "
              f"{o['c']:10.4f} {o['s_tot']:7.4f} "
              f"{100*(o['c']/pred-1):+7.2f}% {z:+8.2f}σ")

    json.dump(dict(olcum=M, onkayit=OK["ongoru"], sonuc=res,
                   lam=lam_tab, C_HIP=C_HIP),
              open(SCR / "170/K2_yuzlesme.json", "w"), indent=1, default=float)
    print(f"\n-> {SCR}/170/K2_yuzlesme.json")


if __name__ == "__main__":
    main()
