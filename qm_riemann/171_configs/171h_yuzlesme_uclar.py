"""
171h — T2c YÜZLEŞME: λ = 0.50 ve λ = 1.30 (ön kayıt ↔ ölçüm)
=============================================================
Yeni ölçüm YOK; `171f`in zaman damgalı ön kayıtları ile `171g`in
ölçümlerinin yüzleşmesi. `169_k1` AYNEN import edilir (üye tanımı
c_WX = KALİB_u2/W_X — 170/171 boyunca aynı).

ÖN KAYITLAR:  ONKAYIT_L050.json (2026-09-04 00:00:05)
              ONKAYIT_L130.json (2026-09-04 00:03:56)

DÜRÜSTLÜK NOTU: `171g`in ekran çıktısında L050'nin gaz-düzeyi **c_ampX**'i
(0.4369, 6 bant — 167'nin kendi penceresi ve üyesi) görüldü; c_WX ≈
c_ampX·W_amp olduğundan bu betiğin ön-mührü artık KÖR DEĞİLDİR ve
öyle işaretlenir. Hükmü belirleyen ön kayıt (00:00:05 / 00:03:56)
bundan ÖNCE mühürlenmiştir ve dokunulmamıştır.

ÖN-MÜHÜR (kör DEĞİL; 4 Eylül 2026 00:15):
 V1 c_WX(L050) ≈ 0.371 (c_ampX·W_amp kestirimi) ⇒ M(0.50) ≈ 1.20:
    M2 (1.143) ve M3 (1.135) **ALTTAN ISKALAR** (≈ −%5), M9 (1.224)
    yakın. Yani eşik ailelerinin düşük-λ kolu FAZLA SIĞ.
 V2 λ=1.30 için kestirim yok (kör): M2/M3 tam 1.0000 der; M1 0.9396,
    M7 0.9674, M9 1.0175.
 V3 Yeni gazların ŞEKLİ A(τ)'ya uyar: 5-bant S/A2Hk sapması rms ≤ %2
    (L050) ve ≤ %2.5 (L130); gaz-başına σ*/2 λ-ailesi bandına
    (0.2269–0.2493) yakın ama L130'da alt ucun ALTINA (≈0.220) düşebilir.

SONUÇ yalnız gerçek koşudan.   Çıktı: scratchpad/171/T2C_YUZLESME.json

SONUÇ (gerçek koşudan):
 c(L050) = 0.3800 ± 0.0072   M(L050) = 1.2266 ± 0.0069
 c(L130) = 0.4540 ± 0.0136   M(L130) = 0.9832 ± 0.0157
 V1 ✓ (kör değil) M(0.50) ölçülen 1.2266: M2 +12.0σ, M3 +13.2σ, M8 +13.8σ,
    M1 +15.3σ, M7 +9.4σ, M0 +32.6σ ile ÖLDÜ; yalnız EMPİRİK M9 tuttu
    (+0.4σ). Türetimli M kimliği YOK.
 V2 ✓ M(1.30) = 0.9832, 1'den −1.1σ ⇒ TEK-TARAFLILIK AYAKTA; M1 (0.9396)
    +2.8σ ile bir kez daha öldü.
 V3 ~ L050 şekilde A(τ)'nun üyesi (rms %1.50, σ*/2 = 0.22779 — bandın
    içinde ✓); L130 DEĞİL (rms %3.62 > öngörü %2.5; σ*/2 = 0.18661,
    bandın %18 ALTINDA — öngörü ≈0.220 idi, ISKA).
 EK: c(0.50) = 0.3800 > c(0.60) = 0.3689 ⇒ c'nin TABANI YOK (bkz. 171k).
"""
import importlib
import json
import os
import sys
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
for _p in ("169_configs", "167_configs", "166_configs", "165_configs",
           "163_configs", "160_configs", "159_configs", "155_configs",
           "154_configs"):
    sys.path.insert(0, str(QM / _p))
K1 = importlib.import_module("169_k1")

SCR = ("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
       "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
SCR171 = SCR + "/171"
LAM5 = ["L115", "Hkeskin", "L085", "L070", "L060"]
C_HIP = 4.0 / np.pi ** 2


def gaz_c(D, g):
    B = K1.saglikli(D[g], lo_max=0.68)
    rows = [K1.band_jk(b["cizgi"]) for b in B]
    c = np.array([r["c_WX"] for r in rows])
    s = np.array([r["sc_WX"] for r in rows])
    lg = np.log(c)
    cg = float(np.exp(lg.mean()))
    s_jk = float(np.sqrt(np.sum((s / c) ** 2)) / len(c)) * cg
    s_sc = float(np.std(lg, ddof=1) / np.sqrt(len(lg))) * cg
    return dict(rows=rows, c=cg, s_jk=s_jk, s_bant=s_sc,
                s_tot=float(np.hypot(s_jk, s_sc)),
                c_bant=[float(x) for x in c],
                K=np.array([r["KALIB_u2"] for r in rows]),
                sK=np.array([r["sKALIB_u2"] for r in rows]),
                W=np.array([r["W_X"] for r in rows]),
                t=np.array([r["tau_eff"] for r in rows]))


def main():
    os.makedirs(SCR171, exist_ok=True)
    D = K1.yukle()
    MO = json.load(open(SCR171 + "/M_ONKAYIT.json"))
    sX_hk = D["Hkeskin"]["sigX"]
    A2 = lambda t: np.exp(-2 * np.pi ** 2 * t ** 2 * sX_hk ** 2)  # noqa: E731
    H = gaz_c(D, "Hkeskin")
    out = {}

    for ad in ("L050", "L130"):
        ON = json.load(open(SCR171 + "/ONKAYIT_%s.json" % ad))
        G = gaz_c(D, ad)
        print("=" * 104)
        print("171h — T2c YÜZLEŞME: %s (λ = %.2f)   ön kayıt %s"
              % (ad, ON["lam"], ON["zaman"]))
        print("=" * 104)
        print("  MARJİNALLER: σ_ds = %.5f (ön kayıt %.5f), σ_X̃ = %.5f "
              "(%.5f), σ_Ĉ = %.5f (%.5f)"
              % (D[ad]["sigds"], ON["sigds"], D[ad]["sigX"], ON["sigX"],
                 D[ad]["sigC"], ON["sigC"]))
        print("\n  bant  τ_eff   KALİB_u2±jk      W_X(ölç)  W_X(ön)   "
              "**c_WX±jk**")
        for i in range(len(G["t"])):
            print("   %d   %.4f  %.4f±%.4f   %.4f    %.4f   %.4f"
                  % (i + 1, G["t"][i], G["K"][i], G["sK"][i], G["W"][i],
                     ON["W_X"][i], G["c_bant"][i]))
        print("\n  **c(%s) = %.4f ± %.4f(jk) ± %.4f(bant) = ±%.4f**  "
              "[5 bant]   4/π² farkı %+.2f%%"
              % (ad, G["c"], G["s_jk"], G["s_bant"], G["s_tot"],
                 100 * (G["c"] / C_HIP - 1)))

        # ölçülen M
        Mo = float(np.exp(np.mean(np.log(G["K"] / H["K"]))))
        sj = np.sqrt(np.sum((G["sK"] / G["K"]) ** 2
                            + (H["sK"] / H["K"]) ** 2)) / len(G["K"]) * Mo
        ss = float(np.std(np.log(G["K"] / H["K"]), ddof=1)
                   / np.sqrt(len(G["K"]))) * Mo
        sMo = float(np.hypot(sj, ss))
        print("  **M(%s) = %.4f ± %.4f**   (bant bant: %s)"
              % (ad, Mo, sMo, "  ".join("%.4f" % x
                                        for x in G["K"] / H["K"])))

        print("\n  ÖN-KAYITLI ÖNGÖRÜLERLE YÜZLEŞME")
        print("  %-34s %8s %8s | %9s %8s | %8s %8s"
              % ("aday", "M(ön)", "c(ön)", "c ölç−ön", "σ_tot", "M ölç−ön",
                 "σ_M"))
        hk = {}
        for nm, r in ON["ongoru"].items():
            dc = 100 * (G["c"] / r["c"] - 1)
            zc = (G["c"] - r["c"]) / G["s_tot"]
            dm = 100 * (Mo / r["M"] - 1)
            zm = (Mo - r["M"]) / sMo
            hk[nm] = dict(M_on=r["M"], c_on=r["c"], dc=dc, zc=float(zc),
                          dm=dm, zm=float(zm))
            print("  %-34s %8.4f %8.4f | %+8.2f%% %+7.1fσ | %+7.2f%% %+7.1fσ"
                  % (nm, r["M"], r["c"], dc, zc, dm, zm))
        print("  %-34s %8s %8.4f | %+8.2f%% %+7.1fσ |"
              % ("(çapa) 4/π²", "—", C_HIP, 100 * (G["c"] / C_HIP - 1),
                 (G["c"] - C_HIP) / G["s_tot"]))

        # şekil: A(τ)'ya uyuyor mu
        S = G["K"] / G["K"][2]
        Ap = A2(G["t"]) / A2(G["t"][2])
        sap = 100 * (S / Ap - 1)
        a = -np.polyfit(G["t"] ** 2, np.log(G["K"]), 1)[0]
        se = float(np.sqrt(a / (2 * np.pi ** 2)))
        ge5 = [MO["genislik"][g]["sX_eff"] for g in LAM5]
        print("\n  ŞEKİL: S/A2Hk − 1 (%%) = %s   rms %.2f%%"
              % ("  ".join("%+6.2f" % x for x in sap),
                 float(np.sqrt((sap ** 2).mean()))))
        print("  gaz-başına α = %.4f ⇒ σ*/2 = %.5f   "
              "(λ-ailesi bandı [%.5f, %.5f]; kendi σ_X̃ = %.5f)"
              % (a, se, min(ge5), max(ge5), D[ad]["sigX"]))

        out[ad] = dict(zaman_onkayit=ON["zaman"], lam=ON["lam"],
                       c=G["c"], s_tot=G["s_tot"], s_jk=G["s_jk"],
                       s_bant=G["s_bant"], c_bant=G["c_bant"],
                       M=Mo, sM=sMo, hukum=hk, sekil=[float(x) for x in sap],
                       alfa=float(a), sX_eff=se,
                       tau=[float(x) for x in G["t"]],
                       KALIB=[float(x) for x in G["K"]],
                       W_X=[float(x) for x in G["W"]])
        print()

    # ---- yedi noktalı M ve c defteri -----------------------------------
    print("=" * 104)
    print("T2c.3 — YEDİ NOKTALI λ DEFTERİ (c ve M)")
    print("=" * 104)
    print("  %-8s %6s | %8s %8s | %8s %8s | %8s"
          % ("gaz", "λ", "c_WX", "±", "M", "±", "σ_X̃"))
    defter = []
    for g, l in (("L130", 1.30), ("L115", 1.15), ("Hkeskin", 1.00),
                 ("L085", 0.85), ("L070", 0.70), ("L060", 0.60),
                 ("L050", 0.50), ("son", None)):
        G = gaz_c(D, g)
        Mv = float(np.exp(np.mean(np.log(G["K"] / H["K"]))))
        print("  %-8s %6s | %8.4f %8.4f | %8.4f %8s | %8.5f"
              % (g, ("%.2f" % l) if l else "—", G["c"], G["s_tot"], Mv,
                 "", D[g]["sigX"]))
        defter.append(dict(gaz=g, lam=l, c=G["c"], s=G["s_tot"], M=Mv,
                           sigX=D[g]["sigX"]))
    out["defter"] = defter
    json.dump(out, open(SCR171 + "/T2C_YUZLESME.json", "w"), indent=1,
              default=float)
    print("\n-> %s/T2C_YUZLESME.json" % SCR171)


if __name__ == "__main__":
    main()
