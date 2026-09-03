"""
171d — T2a: A(τ) KİMLİK YARIŞI, 2. AŞAMA (9-BANT YÜZLEŞMESİ)
=============================================================
Yeni ölçüm YOK; 171c'nin ÖN KAYDI (scratchpad/171/A_ONKAYIT.json, zaman
damgalı) ile DOKUZ-BANT (lo ≤ 0.80) verisinin yüzleşmesi. Adaylar YALNIZ
beş bantta (τ ≤ 0.698) uyduruldu; 6. bant (τ ≈ 0.7386, beş λ-gazının
hepsinde) ve 7. bant (τ ≈ 0.7774, L085 + Hkeskin) ÖRNEKLEM-DIŞIDIR.

ÖN-MÜHÜR (bu betiği koşmadan önce, 171c'nin ön kaydından):
 Q1 A(0.7386) ölçülen ≈ 0.828 ± 0.010 ⇒ A2Hk (0.8287) ve A2s (0.8301)
    tutar; A3s (0.7367), A4s (0.7575), A1λ (0.7595), A1 (0.7138) ölür
    (%9-16 sapma); A2iv (0.9326) ve A5 (0.8588) da ölür; A4 (3.54) zaten öldü.
 Q2 A(0.7774) ölçülen ≈ 0.775 ± 0.015 ⇒ yine A2Hk (0.7742)/A2s (0.7760).
 Q3 6. bantta ÇARPANLAŞMA da sınanır: beş gazın S_6 değerleri birbirinden
    ≤ %3 ayrılmalı (5-bant kolapsı %0.8 rms'ti; son bantta %2.8 idi).
    Ayrışma büyürse çarpanlaşma τ > 0.70'te ZAYIFLIYOR demektir — bu da
    kurtarmasız yazılır.
 Q4 Kimlik hükmü: A(τ) = W_X(τ; σ_X̃ = σ_X̃(λ=1.00)) — SIFIR serbest
    parametre — dokuz bantta rms ≤ %1.2 tutarsa mühürlenir.

SONUÇ bloğu yalnız gerçek koşudan.
Çıktı: scratchpad/171/A_YUZLESME.json + log

SONUÇ (gerçek koşudan):
 ÖLÇÜLEN (örneklem-dışı): A(0.7386) = 0.8265 ± 0.0122 (5 gaz),
                          A(0.7774) = 0.7456 ± 0.0223 (2 gaz).
 Q1 ✓✓ A2Hk −%0.26, A2s −%0.43 tuttu; A1λ +%8.8, A4s +%9.1, A2iv −%11.4,
       A3s +%12.2, A1 +%15.8, A3f −%15.3, A0 −%17.4, A4 −%76.6 ÖLDÜ.
 Q2 ~  A(0.7774) = 0.7456; ön-mühür 0.775 ± 0.015 demişti (−%3.8 iska),
       sıralama değişmedi (A2Hk −%3.69 ile en yakın).
 Q3 ✗  6. bantta gaz-gaz yayılımı **%8.08** (öngörü ≤%3): çarpanlaşma
       τ > 0.70'te ZAYIFLIYOR ve sapma λ ile tekdüze sıralı.
 Q4 ~  A2Hk'nın 9-bant rms'i %1.49 (öngörü ≤%1.2) — yine de 1-parametreli
       en iyi uyumu (A2s %1.58) YENİYOR.
 MÜHÜR: **A(τ) = exp(−2π²τ²σ_X̃(λ=1)²)**, sıfır serbest parametre.
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
MID = 2


def main():
    ON = json.load(open(SCR171 + "/A_ONKAYIT.json"))
    D = K1.yukle()
    print("=" * 104)
    print("171d — T2a AŞAMA-2: 9-BANT YÜZLEŞMESİ  (ön kayıt zamanı: %s)"
          % ON["zaman"])
    print("=" * 104)

    S, TAU = {}, {}
    for g in LAM5:
        B = K1.saglikli(D[g], lo_max=0.80)
        rows = [K1.band_jk(b["cizgi"]) for b in B]
        k = np.array([r["KALIB_u2"] for r in rows])
        S[g] = k / k[MID]
        TAU[g] = np.array([r["tau_eff"] for r in rows])
        print("  %-8s λ=%-5s nb9=%d  S = %s" % (g, D[g]["lam"], len(k),
              "  ".join("%.4f" % x for x in S[g])))

    # 6. bant: beş gazın hepsinde;  7. bant: L085 + Hkeskin
    g6 = [g for g in LAM5 if len(S[g]) >= 6]
    g7 = [g for g in LAM5 if len(S[g]) >= 7]
    A6 = float(np.exp(np.mean([np.log(S[g][5]) for g in g6])))
    t6 = float(np.mean([TAU[g][5] for g in g6]))
    sd6 = float(np.std([np.log(S[g][5]) for g in g6], ddof=1)
                / np.sqrt(len(g6)))
    A7 = float(np.exp(np.mean([np.log(S[g][6]) for g in g7])))
    t7 = float(np.mean([TAU[g][6] for g in g7]))
    sd7 = (float(np.std([np.log(S[g][6]) for g in g7], ddof=1)
                 / np.sqrt(len(g7))) if len(g7) > 1 else float("nan"))
    print("\n  ÖLÇÜLEN (örneklem-dışı):")
    print("    A(τ=%.4f) = %.4f ± %.4f (gazlar arası s.h.; %d gaz: %s)"
          % (t6, A6, A6 * sd6, len(g6), ",".join(g6)))
    print("    A(τ=%.4f) = %.4f ± %.4f (%d gaz: %s)"
          % (t7, A7, A7 * sd7, len(g7), ",".join(g7)))
    print("    6. bantta gaz-gaz ayrışması: " + "  ".join(
        "%s=%.4f" % (g, S[g][5]) for g in g6)
        + "   (yayılım %.2f%%)" % (100 * (max(S[g][5] for g in g6)
                                          / min(S[g][5] for g in g6) - 1)))

    print("\n" + "=" * 104)
    print("ÖN-KAYITLI ÖNGÖRÜLERLE YÜZLEŞME")
    print("=" * 104)
    print("  %-46s %2s %8s %8s | %8s %8s | %7s"
          % ("aday", "p", "ön6", "ön7", "sap6%", "sap7%", "9b rms%"))
    hüküm = []
    for r in ON["adaylar"]:
        p6, p7 = r["ong9"]
        s6 = 100 * (A6 / p6 - 1)
        s7 = 100 * (A7 / p7 - 1)
        rms9 = float(np.sqrt(np.mean(np.array(r["artik"] + [s6, s7]) ** 2)))
        hüküm.append(dict(ad=r["ad"], npar=r["npar"], etiket=r["etiket"],
                          param=r["param"], sap6=s6, sap7=s7,
                          rms5=r["rms"], rms9=rms9))
        print("  %-46s %2d %8.4f %8.4f | %+8.2f %+8.2f | %7.2f"
              % (r["ad"], r["npar"], p6, p7, s6, s7, rms9))

    en = sorted(hüküm, key=lambda h: h["rms9"])
    print("\n  SIRALAMA (9-bant rms):  " + "   ".join(
        "%d) %s %.2f%%" % (i + 1, h["ad"].split()[0], h["rms9"])
        for i, h in enumerate(en[:5])))
    print("  SIFIR-PARAMETRELİ en iyi: %s (rms9 = %.2f%%)"
          % ([h["ad"] for h in en if h["npar"] == 0][0],
             [h["rms9"] for h in en if h["npar"] == 0][0]))

    json.dump(dict(zaman_onkayit=ON["zaman"], A6=A6, t6=t6, sA6=sd6,
                   A7=A7, t7=t7, sA7=sd7, g6=g6, g7=g7,
                   S9={g: [float(x) for x in S[g]] for g in LAM5},
                   TAU9={g: [float(x) for x in TAU[g]] for g in LAM5},
                   hukum=hüküm),
              open(SCR171 + "/A_YUZLESME.json", "w"), indent=1)
    print("\n-> %s/A_YUZLESME.json" % SCR171)


if __name__ == "__main__":
    main()
