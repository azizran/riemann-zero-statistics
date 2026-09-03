"""
171i — T3: GERÇEK GAZIN DİLİ — ΔM = +%5.6'nın adresi
=====================================================
Yeni ölçüm YOK. Girdi: scratchpad/167/C_<gaz>.json + 171e'nin M defteri +
170 §K3'ün ÖLÇÜLEN bacak kırpma aktarımları (rapordaki sayılar, kaynak
`170e_bacak_aktarim.py` → `169/K2b_*.json`).

SORULAR (kalem):
  (a) Gerçek gaz A(τ)'nun üyesi mi — 9 bantta ÖRNEKLEM-DIŞI dahil?
  (b) ΔM'yi hangi ÖLÇÜLÜR marjinal öngörür?
  (c) 170 §K3'ün ölçülen bacak-kırpma kapanışı (c dilinde %86) A·M
      dilinde neye karşılık gelir?

ÖN-MÜHÜR (koşudan ÖNCE, 4 Eylül 2026 00:10):
 S1 M(son) = 1.0561 ± 0.0060; ayrışma g_E +%2.49, g_X² +%0.15, θ +%2.90
    (≈ 45/55 paylaşım).
 S2 A-ÜYELİĞİ TUTAR: son'un gaz-başına Gauss genişliği σ*/2 = 0.24493,
    σ_X̃(λ=1)'den +%1.19 — λ-ailesinin ±%5 bandının İÇİNDE ve Hkeskin'in
    kendisinden (−%2.57) DAHA YAKIN. 9-bant örneklem-dışı sapması ≤ %4.
 S3 KIRPMA KURALININ M DİLİ: M_kırpma = (ρ₃/ρ₃⁰)^{4/3}·W_X(son)/W_X(Hk)
    = 1.0182·1.0333 = 1.0521 ⇒ ΔM'nin **%93'ünü** kapatır (c dilindeki
    %86'dan iyi; çünkü W_X payı ayrıldı).
 S4 HİÇBİR MARJİNAL ΔM'yi öngörmez: σ_X̃ +%5.02 (+8.4σ), σ_Ĉ ≈ +%5.3,
    σ_ds ≈ +%4.6 açık bırakır.

SONUÇ yalnız gerçek koşudan.   Çıktı: scratchpad/171/T3.json

SONUÇ (gerçek koşudan):
 S1 ✓✓ M(son) = 1.0561 ± 0.0060 (ΔM = +%5.61 = +9.4σ);
       g_E +%2.49, g_X² +%0.15, θ +%2.90; log-paylar %45 / %3 / %52.
 S2 ~  A-üyeliği TUTTU: 5-bant rms %0.77 (λ-ailesi 0.53–1.24%),
       σ*/2 = 0.24493 (+%1.19). Ama örneklem-dışı 7. bant −%7.12
       (öngörü ≤%4 — ISKA); L085 aynı bantta −%6.53, yani τ > 0.77
       bütün gazlarda A'nın dışı.
 S3 ✓✓ M_kırpma = 1.0182 · 1.0333 = 1.0521 ⇒ ΔM'nin **%93'ü** kapandı,
       kalan +%0.37 = +0.7σ.
 S4 ✓  Hiçbir marjinal öngörmüyor: σ_ds +7.9σ, σ_X̃ +8.4σ, σ_Ĉ +8.9σ.
       (g_X satırı KÖTÜ KOŞULLU — λ-ailesindeki yayılımı ±%0.3 ve tekdüze
       bile değil; kimlik sayılmaz.)
 EK: aynı kırpma kuralı λ ekseninde M dilinde −10.6σ / −24.4σ / −43.3σ
     çöküyor (c dilindeki −6.7σ/−10.1σ/−13.4σ'dan 3 kat sert).
"""
import importlib
import json
import os
import sys

import numpy as np
from pathlib import Path

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
# 170 §K3.2/K3.3'ün ÖLÇÜLEN kırpma aktarımları (169_k2b ile)
RHO = {"Hkeskin": (0.8293, 0.6586, 0.5284), "son": (0.8321, 0.6585, 0.5356),
       "L085": (0.8306, 0.6502, 0.5285), "L070": (0.8448, 0.6641, 0.5542),
       "L060": (0.8611, 0.6894, 0.5904)}


def main():
    os.makedirs(SCR171, exist_ok=True)
    D = K1.yukle()
    MO = json.load(open(SCR171 + "/M_ONKAYIT.json"))
    M, sM = MO["M"], MO["sM"]
    print("=" * 100)
    print("171i — T3: GERÇEK GAZIN DİLİ (ΔM = M(son) − 1)")
    print("=" * 100)
    ay = MO["ayrisim"]["son"]
    print("  M(son) = %.4f ± %.4f   ⇒  ΔM = %+.2f%%  (%.1fσ)"
          % (M["son"], sM["son"], 100 * (M["son"] - 1),
             (M["son"] - 1) / sM["son"]))
    print("  ayrışma:  g_E/g₀ = %.4f (%+.2f%%)   (g_X/g₀)² = %.4f (%+.2f%%)"
          "   θ/θ₀ = %.4f (%+.2f%%)"
          % (ay["gE"], 100 * (ay["gE"] - 1), ay["gX2"], 100 * (ay["gX2"] - 1),
             ay["th"], 100 * (ay["th"] - 1)))
    pay_gE = np.log(ay["gE"]) / np.log(M["son"])
    print("  ΔM'nin log-payları:  g_E %.0f%%   g_X² %.0f%%   θ %.0f%%"
          % (100 * pay_gE, 100 * np.log(ay["gX2"]) / np.log(M["son"]),
             100 * np.log(ay["th"]) / np.log(M["son"])))

    # ---------- (a) A-üyeliği, 9 bant örneklem-dışı dahil --------------
    print("\n" + "=" * 100)
    print("T3.1 — GERÇEK GAZ A(τ)'NUN ÜYESİ Mİ?  (şekil sınavı)")
    print("=" * 100)
    sX_hk = D["Hkeskin"]["sigX"]
    A2 = lambda t: np.exp(-2 * np.pi ** 2 * t ** 2 * sX_hk ** 2)  # noqa: E731
    tab = {}
    for g in LAM5 + ["son"]:
        B = K1.saglikli(D[g], lo_max=0.80)
        rows = [K1.band_jk(b["cizgi"]) for b in B]
        k = np.array([r["KALIB_u2"] for r in rows])
        t = np.array([r["tau_eff"] for r in rows])
        S = k / k[2]
        Ap = A2(t) / A2(t[2])
        tab[g] = dict(S=[float(x) for x in S], tau=[float(x) for x in t],
                      sap=[float(100 * (S[b] / Ap[b] - 1))
                           for b in range(len(S))])
        print("  %-8s nb9=%d  S/A2Hk − 1 (%%) = %s"
              % (g, len(S), "  ".join("%+6.2f" % x for x in tab[g]["sap"])))
    s5 = np.array(tab["son"]["sap"][:5])
    print("\n  son, 5 bant: rms %.2f%%, en kötü %.2f%%   "
          "(λ-ailesi 5-bant rms'leri: %s)"
          % (np.sqrt((s5 ** 2).mean()), np.abs(s5).max(),
             "  ".join("%.2f%%" % np.sqrt((np.array(tab[g]["sap"][:5]) ** 2)
                                          .mean()) for g in LAM5)))
    print("  son, ÖRNEKLEM-DIŞI 6. ve 7. bant sapması: %+.2f%%  %+.2f%%"
          % (tab["son"]["sap"][5], tab["son"]["sap"][6]))
    gen = MO["genislik"]
    print("  gaz-başına genişlik σ*/2:  son %.5f (%+.2f%% σ_X̃(λ=1)'den);  "
          "λ-ailesi %s"
          % (gen["son"]["sX_eff"], 100 * (gen["son"]["oran_hk"] - 1),
             "  ".join("%.5f" % gen[g]["sX_eff"] for g in LAM5)))

    # ---------- (b) hangi marjinal ΔM'yi öngörür ----------------------
    print("\n" + "=" * 100)
    print("T3.2 — ΔM'yi HANGİ ÖLÇÜLÜR BÜYÜKLÜK ÖNGÖRÜR?")
    print("=" * 100)
    lam = np.array([D[g]["lam"] for g in LAM5], float)
    Mv = np.array([M[g] for g in LAM5])
    print("  %-28s %10s %10s %10s %10s"
          % ("öngörücü (λ-ailesinden)", "değer(son)", "M öngörü", "ölç−öng",
             "σ"))
    ong = {}
    for et, key in (("σ_ds", "sigds"), ("σ_X̃", "sigX"), ("σ_Ĉ", "sigC"),
                    ("g_E", None), ("g_X", None)):
        if key:
            arr = np.array([D[g][key] for g in LAM5])
            v = D["son"][key]
        else:
            arr = np.array([D[g]["artik"][et.replace("_", "")]
                            for g in LAM5])
            v = D["son"]["artik"][et.replace("_", "")]
        p = np.polyfit(np.log(arr), np.log(Mv), 2)
        mo = float(np.exp(np.polyval(p, np.log(v))))
        z = (M["son"] - mo) / sM["son"]
        ong[et] = dict(deger=float(v), M_ong=mo,
                       sapma=float(100 * (M["son"] / mo - 1)), z=float(z))
        print("  %-28s %10.5f %10.4f %+9.2f%% %+9.1f"
              % (et, v, mo, 100 * (M["son"] / mo - 1), z))
    # ölçülen kırpma aktarımı (169_k2b) — marjinal DEĞİL, korelatör düzeyi
    WXs = np.mean([K1.band_jk(b["cizgi"])["W_X"]
                   for b in K1.saglikli(D["son"], 0.68)])
    WXh = np.mean([K1.band_jk(b["cizgi"])["W_X"]
                   for b in K1.saglikli(D["Hkeskin"], 0.68)])
    r3 = (RHO["son"][2] / RHO["Hkeskin"][2]) ** (4.0 / 3.0)
    M_kirp = float(r3 * WXs / WXh)
    print("\n  [korelatör düzeyi, marjinal DEĞİL] 3-bacak kırpma aktarımı:")
    print("     ρ₃(son)/ρ₃(Hk) = %.4f  ⇒ (…)^{4/3} = %.4f  (170 §K3.2 kuralı (e))"
          % (RHO["son"][2] / RHO["Hkeskin"][2], r3))
    print("     W_X(son)/W_X(Hk) = %.4f  ⇒  M_kırpma = %.4f"
          % (WXs / WXh, M_kirp))
    kap = (M_kirp - 1) / (M["son"] - 1)
    print("     ⇒ ΔM'nin kapanan payı = %.0f%%   (kalan açık %+.2f%% = %+.1fσ)"
          % (100 * kap, 100 * (M["son"] / M_kirp - 1),
             (M["son"] - M_kirp) / sM["son"]))
    # aynı kural λ ekseninde ne yapıyor?
    print("\n  AYNI KURAL λ EKSENİNDE (170 §K3'ün en sert kısıtı):")
    for g in ("L085", "L070", "L060"):
        WXg = np.mean([K1.band_jk(b["cizgi"])["W_X"]
                       for b in K1.saglikli(D[g], 0.68)])
        mk = (RHO[g][2] / RHO["Hkeskin"][2]) ** (4.0 / 3.0) * WXg / WXh
        print("     %-6s M_kırpma = %.4f   M ölçülen = %.4f   sapma %+.2f%% "
              "(%+.1fσ)" % (g, mk, M[g], 100 * (M[g] / mk - 1),
                            (M[g] - mk) / sM[g]))

    json.dump(dict(M_son=M["son"], sM_son=sM["son"], ayrisim=ay,
                   sekil=tab, ongorucu=ong, M_kirpma=M_kirp,
                   kapanan=float(kap)),
              open(SCR171 + "/T3.json", "w"), indent=1, default=float)
    print("\n-> %s/T3.json" % SCR171)


if __name__ == "__main__":
    main()
