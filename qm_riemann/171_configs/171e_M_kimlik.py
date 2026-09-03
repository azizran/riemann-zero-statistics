"""
171e — T2b: M(λ)'NIN KİMLİĞİ (tek-değişkenlik + aday aileler + A-genişliği)
===========================================================================
Yeni ölçüm YOK. Girdi: scratchpad/167/C_<gaz>.json. `169_k1` (o da
`166_T1`) AYNEN import edilir.

M TANIMI (çarpanlaşmadan, A'yı ÖZDEŞ olarak eleyen biçim):
    M(g) = geo.ort_b [ KALİB_b(g) / KALİB_b(Hkeskin) ]      (5 bant)
Hata: bantlar arası saçılım/√5 (bant-içi jackknife ile birleştirilir).

SINAV SIRASI (kalem): (1) TEK-DEĞİŞKENLİK, (2) TEK-TARAFLILIK, (3) uyum.

ÖN-MÜHÜR (koşudan ÖNCE, 3 Eylül 2026 23:55):
 N1 M = 0.999 / 1.000 / 1.018 / 1.067 / 1.129 (λ = 1.15→0.60), ±0.005;
    M(son) = 1.056.
 N2 ÖZDEŞLİK M = (g_E/g₀)(g_X/g₀)²(θ/θ₀) ‰1 içinde kapanır; g_X λ-değişmez
    (≤‰3); M'nin yükselişini g_E ile θ NEREDEYSE EŞİT paylaşır
    (λ=0.60'ta +%6.3 ve +%4.9).
 N3 TEK-DEĞİŞKENLİK: λ-ailesi İÇİNDE λ, σ_ds, σ_X̃, σ_Ĉ birebir dejeneredir
    (hepsi birbirinin tekdüze fonksiyonu) ⇒ "hangi değişken" sorusu λ
    ailesinden CEVAPLANAMAZ. Dejenerasyonu KIRAN tek nokta GERÇEK gazdır:
    σ_X̃(son) = 0.23384 ⇒ λ-eşdeğeri ≈ 0.955 ⇒ σ_X̃-tek-değişkenli M
    öngörüsü ≈ 1.005; ölçülen 1.056 ⇒ **+%5.1 (≈ +5σ) ile M σ_X̃'in
    tek-değişkenli fonksiyonu DEĞİL.** Aynı sonuç σ_ds ve σ_Ĉ için de.
 N4 TEK-TARAFLILIK: M(1.15) − 1 = −0.001 ± 0.005 ⇒ 1'den ayırt edilemez
    (ama TEK nokta; λ=1.30 hakem). Öngörü: türetimli aday M1 (Gram-doyum,
    ∝λ) λ=1.15'te ≈ −%2, M7 (merdiven-sürüşlü kesir) ≈ −%4 verip
    ölürler; eşikli aileler (M2, M3) tanımı gereği geçer.
 N5 A-GENİŞLİĞİ λ ile kayıyor: gaz-başına σ*_g/2 ≈ 0.226 / 0.237 / 0.247 /
    0.252 / 0.249 (λ = 1.15→0.60); ortalama 0.242 = σ_X̃(λ=1) ama
    yayılım ±%5 ve SİSTEMATİK ⇒ §T2a.4'ün 6. bant yayılımının kaynağı.

SONUÇ bloğu yalnız gerçek koşudan.
Çıktı: scratchpad/171/M_ONKAYIT.json + log

SONUÇ (gerçek koşudan, 23:52:46):
 N1 ~  M = 1.0022/1.0000/1.0144/1.0595/1.1205 (λ = 1.15→0.60),
       M(son) = 1.0561 ± 0.0060. Orta-bant ön-mührünün %0.3–0.8 altında
       (tanım farkı: bant-ortalaması ↔ orta bant). Tam isabet DEĞİL.
 N2 ✓  Özdeşlik ‰0'da kapandı. λ ekseninde g_E ve θ neredeyse eşit paylı
       (L060: +%6.31 ve +%4.93), g_X² ≤ ‰5. Kesim → yalnız θ (g_E ARTAR),
       pencere → yalnız g_X² (HkT4b 0.7833).
 N3 ✓✓ TEK-DEĞİŞKENLİK ÖLDÜ: aynı σ_X̃'te dokuz gaz −%17.8 … +%5.0
       (−17.7σ … +8.4σ). σ_ds ve σ_Ĉ için de aynı.
 N4 ✓  M(1.15) = 1.0022 ± 0.0054 (1'den +0.4σ). Tek-taraflılık M1'i
       (−%3.3), M1b'yi (−%4.0) ve M7'yi (−%2.1) öldürdü. Ayakta: M2
       (rms %1.19) ve M3 (%1.24) — eşik λ_c = 1.0109 SABİT, uydurulmadı.
 N5 ✓  σ*_g/2 = 0.22685/0.23582/0.24438/0.24926/0.24596 (λ = 1.15→0.60);
       ort 0.24045 = σ_X̃(λ=1) − %0.66, yayılım %9.9. Kesim gazlarında
       +%2.7…+%41.0, pencere gazlarında −%7.5…−%21.7.
"""
import importlib
import json
import os
import sys
import time
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
HEP = LAM5 + ["son", "K090", "K070", "HA4", "E060",
              "HkT2a", "HkT2b", "HkT4a", "HkT4b"]
LAM_C = 1.9147 / 1.894        # rms S'(λ) = N̄' doyum noktası = 1.01094


def bantlar(D, g, lo_max=0.68):
    B = K1.saglikli(D[g], lo_max=lo_max)
    return [K1.band_jk(b["cizgi"]) for b in B]


def main():
    os.makedirs(SCR171, exist_ok=True)
    D = K1.yukle()
    R = {g: bantlar(D, g) for g in HEP}
    KH = np.array([r["KALIB_u2"] for r in R["Hkeskin"]])

    # ---------------- M(g) ---------------------------------------------
    print("=" * 104)
    print("171e — T2b: M(λ) = geo.ort_b [KALİB_b(g)/KALİB_b(Hkeskin)]")
    print("=" * 104)
    M, sM = {}, {}
    print("  %-8s %6s | %-42s | %8s %8s" % ("gaz", "λ", "bant bant oran",
                                            "M", "±"))
    for g in HEP:
        k = np.array([r["KALIB_u2"] for r in R[g]])
        if len(k) != len(KH):
            continue
        o = k / KH
        lg = np.log(o)
        M[g] = float(np.exp(lg.mean()))
        sj = np.array([r["sKALIB_u2"] / r["KALIB_u2"] for r in R[g]])
        sh = np.array([r["sKALIB_u2"] / r["KALIB_u2"] for r in R["Hkeskin"]])
        s_jk = float(np.sqrt(np.sum(sj ** 2 + sh ** 2)) / len(o)) * M[g]
        s_sc = float(np.std(lg, ddof=1) / np.sqrt(len(lg))) * M[g]
        sM[g] = float(np.hypot(s_jk, s_sc))
        print("  %-8s %6s | %s | %8.4f %8.4f"
              % (g, D[g]["lam"], "  ".join("%.4f" % x for x in o),
                 M[g], sM[g]))

    # ---------------- ÖZDEŞLİK: M = (g_E/g₀)(g_X/g₀)²(θ/θ₀) ------------
    print("\n" + "=" * 104)
    print("T2b.1 — M'NİN ÇARPAN AYRIŞMASI  (KALİB_u2 = g_E·g_X²·θ)")
    print("=" * 104)
    def art(g, k):
        return D[g]["artik"][k]
    gE0, gX0 = art("Hkeskin", "gE"), art("Hkeskin", "gX")
    kal0 = float(np.exp(np.mean(np.log(KH))))
    th0 = kal0 / (gE0 * gX0 ** 2)
    print("  %-8s %6s | %8s %8s %8s | %9s %9s | %7s"
          % ("gaz", "λ", "g_E/g₀", "(g_X/g₀)²", "θ/θ₀", "çarpım", "M(ölçülen)",
             "fark%"))
    ayr = {}
    for g in HEP:
        if g not in M:
            continue
        kal = float(np.exp(np.mean(np.log(
            [r["KALIB_u2"] for r in R[g]]))))
        gE, gX = art(g, "gE"), art(g, "gX")
        th = kal / (gE * gX ** 2)
        c = (gE / gE0) * (gX / gX0) ** 2 * (th / th0)
        ayr[g] = dict(gE=gE / gE0, gX2=(gX / gX0) ** 2, th=th / th0,
                      carp=float(c))
        print("  %-8s %6s | %8.4f %8.4f %8.4f | %9.4f %9.4f | %+7.3f"
              % (g, D[g]["lam"], gE / gE0, (gX / gX0) ** 2, th / th0,
                 c, M[g], 100 * (c / M[g] - 1)))

    # ---------------- TEK-DEĞİŞKENLİK ----------------------------------
    print("\n" + "=" * 104)
    print("T2b.2 — TEK-DEĞİŞKENLİK: M yalnız hangi değişkenin fonksiyonu?")
    print("=" * 104)
    lam = np.array([D[g]["lam"] for g in LAM5], float)
    sX = np.array([D[g]["sigX"] for g in LAM5])
    sd = np.array([D[g]["sigds"] for g in LAM5])
    sC = np.array([D[g]["sigC"] for g in LAM5])
    Mv = np.array([M[g] for g in LAM5])
    print("  λ-ailesi içinde dejenerasyon (hepsi tekdüze artan):")
    print("     λ      : " + "  ".join("%.4f" % x for x in lam))
    print("     σ_ds   : " + "  ".join("%.4f" % x for x in sd))
    print("     σ_X̃    : " + "  ".join("%.4f" % x for x in sX))
    print("     σ_Ĉ    : " + "  ".join("%.4f" % x for x in sC))
    print("     M      : " + "  ".join("%.4f" % x for x in Mv))
    print("  ⇒ λ-ailesi TEK BOYUTLUDUR: 'hangi değişken' sorusu bu aileden"
          " cevaplanamaz.")
    print("\n  DEJENERASYONU KIRAN GAZLAR (λ = 1.00 ama σ farklı):")
    print("  %-8s %8s %8s %8s | %9s %9s %9s | %9s | %7s"
          % ("gaz", "σ_ds", "σ_X̃", "σ_Ĉ", "λ_eş(ds)", "λ_eş(X̃)", "λ_eş(Ĉ)",
             "M_öngörü", "ölç−öng"))
    kir = {}
    for g in ("son", "HkT2a", "HkT2b", "HkT4a", "HkT4b",
              "K090", "K070", "HA4", "E060"):
        if g not in M:
            continue
        le = {}
        for et, v, arr in (("ds", D[g]["sigds"], sd), ("X", D[g]["sigX"], sX),
                           ("C", D[g]["sigC"], sC)):
            p = np.polyfit(np.log(arr), np.log(lam), 2)
            le[et] = float(np.exp(np.polyval(p, np.log(v))))
        pM = np.polyfit(np.log(sX), np.log(Mv), 2)
        Mong = float(np.exp(np.polyval(pM, np.log(D[g]["sigX"]))))
        z = (M[g] - Mong) / sM[g]
        kir[g] = dict(lam_ds=le["ds"], lam_X=le["X"], lam_C=le["C"],
                      M_ong=Mong, M=M[g], sM=sM[g], z=float(z))
        print("  %-8s %8.4f %8.4f %8.4f | %9.4f %9.4f %9.4f | %9.4f | "
              "%+6.2f%% (%+.1fσ)"
              % (g, D[g]["sigds"], D[g]["sigX"], D[g]["sigC"],
                 le["ds"], le["X"], le["C"], Mong,
                 100 * (M[g] / Mong - 1), z))

    # ---------------- ADAY M AİLELERİ ----------------------------------
    print("\n" + "=" * 104)
    print("T2b.3 — ADAY M AİLELERİ  (1. sınav: TEK-TARAFLILIK, λ≥1'de düz)")
    print("=" * 104)
    AX, BX = np.polyfit(lam ** 2, sX ** 2, 1)
    u = sX ** 2
    uc = float(AX * LAM_C ** 2 + BX)
    fL = AX * lam ** 2 / (AX * lam ** 2 + BX)
    gEr = np.array([ayr[g]["gE"] for g in LAM5])
    gX2r = np.array([ayr[g]["gX2"] for g in LAM5])

    aday = []

    def kayit(ad, npar, et, f, sinir=None):
        """f(lam, sX) -> M; sinir verilirse 1 parametre uydurulur."""
        if sinir is None:
            v = f(lam, sX)
            p = None
        else:
            gr = np.linspace(sinir[0], sinir[1], 8001)
            best = (np.inf, None)
            for q in gr:
                r = np.log(Mv) - np.log(np.abs(f(lam, sX, q)))
                s = float((r ** 2).sum())
                if s < best[0]:
                    best = (s, q)
            p = float(best[1])
            v = f(lam, sX, p)
        r = 100 * (Mv / v - 1)
        rec = dict(ad=ad, npar=npar, etiket=et, param=p,
                   rms=float(np.sqrt((r ** 2).mean())),
                   artik=[float(x) for x in r],
                   ong=[float(x) for x in v],
                   tek_taraf=float(v[0]))     # λ=1.15 öngörüsü
        aday.append(rec)
        print("  %-44s %2d %-9s %6.2f | λ=1.15 öngörü %.4f (ölç %.4f) | %s"
              % (ad, npar, ("%.4g" % p) if p is not None else "—",
                 rec["rms"], v[0], Mv[0],
                 "  ".join("%+5.2f" % x for x in r)))
        return rec

    kayit("M0 düz (168 §A3: KALİB λ-değişmez)", 0, "türetimli",
          lambda l, s: np.ones_like(l))
    kayit("M8 ÖLÇÜLEN (g_E/g₀)(g_X/g₀)² [θ λ-değişmez]", 0, "türetimli-ölçülen",
          lambda l, s: gEr * gX2r)
    kayit("M1 Gram-doyum 1/g_E = 1+Cλ", 1, "türetimli",
          lambda l, s, C: (1 + C) / (1 + C * l), (0.001, 3.0))
    kayit("M1b Gram-doyum 1/g_E = 1+Cλ²", 1, "türetimli",
          lambda l, s, C: (1 + C) / (1 + C * l ** 2), (0.001, 3.0))
    kayit("M2 kırık kuvvet, λ_c = %.4f SABİT" % LAM_C, 1, "türetimli",
          lambda l, s, p: np.where(l >= LAM_C, 1.0, (LAM_C / l) ** p),
          (0.01, 2.0))
    kayit("M3 varyans-açığı 1+κ(u_c−u)_+/u_c", 1, "türetimli",
          lambda l, s, k: 1 + k * np.maximum(0.0, (uc - s ** 2) / uc),
          (0.01, 3.0))
    kayit("M7 merdiven-sürüşlü kesir f_L^{−q}", 1, "türetimli",
          lambda l, s, q: (fL / fL[1]) ** (-q), (0.0, 3.0))
    kayit("M9 log-λ kuadratik", 2, "EMPİRİK-ÇIPA",
          lambda l, s: np.exp(np.polyval(np.polyfit(np.log(lam),
                                                    np.log(Mv), 2),
                                         np.log(l))))
    # 170 §K3.3'ün ÖLÜ referansı (aday DEĞİL; ölçek için)
    PiT = np.array([0.73806, 0.77072, 0.80875, 0.85336, 0.88654])
    v = PiT / PiT[1]
    r = 100 * (Mv / v - 1)
    print("  %-44s %2s %-9s %6.2f | λ=1.15 öngörü %.4f (ölç %.4f) | ÖLÜ "
          "(170), aday DEĞİL | %s"
          % ("[ref] ΠT/ΠT₀ (170 §K3.3)", "-", "—",
             float(np.sqrt((r ** 2).mean())), v[0], Mv[0],
             "  ".join("%+5.2f" % x for x in r)))

    # ---------------- A'nın gaz-başına genişliği -----------------------
    print("\n" + "=" * 104)
    print("T2b.4 — A(τ)'NUN GAZ-BAŞINA GENİŞLİĞİ (λ-değişmezlik ne kadar tam?)")
    print("=" * 104)
    print("  %-8s %6s | %8s %8s | %10s %10s" % ("gaz", "λ", "α_g", "σ*_g/2",
                                                "/σ_X̃(Hk)", "/σ_X̃(kendi)"))
    gen = {}
    for g in HEP:
        if g not in M:
            continue
        t = np.array([r["tau_eff"] for r in R[g]])
        y = np.log([r["KALIB_u2"] for r in R[g]])
        a = -np.polyfit(t ** 2, y, 1)[0]
        se = float(np.sqrt(a / (2 * np.pi ** 2)))
        gen[g] = dict(alfa=float(a), sX_eff=se,
                      oran_hk=float(se / D["Hkeskin"]["sigX"]),
                      oran_ken=float(se / D[g]["sigX"]))
        print("  %-8s %6s | %8.4f %8.5f | %+9.2f%% %+9.2f%%"
              % (g, D[g]["lam"], a, se,
                 100 * (se / D["Hkeskin"]["sigX"] - 1),
                 100 * (se / D[g]["sigX"] - 1)))
    ge5 = np.array([gen[g]["sX_eff"] for g in LAM5])
    print("  λ-ailesi: ort = %.5f  (σ_X̃(Hk) = %.5f, fark %+.2f%%)  "
          "yayılım %.1f%%"
          % (ge5.mean(), D["Hkeskin"]["sigX"],
             100 * (ge5.mean() / D["Hkeskin"]["sigX"] - 1),
             100 * (ge5.max() / ge5.min() - 1)))

    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    json.dump(dict(zaman=ts, M=M, sM=sM, ayrisim=ayr, kiran=kir,
                   adaylar=aday, genislik=gen, LAM_C=LAM_C,
                   AX=float(AX), BX=float(BX), uc=uc),
              open(SCR171 + "/M_ONKAYIT.json", "w"), indent=1, default=float)
    print("\n-> %s/M_ONKAYIT.json  (%s)" % (SCR171, ts))


if __name__ == "__main__":
    main()
