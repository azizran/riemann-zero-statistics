# -*- coding: utf-8 -*-
"""
172d — G2: θ'NIN BİRLEŞİK YASASI ARANIYOR (4 Eylül 2026)
Veri: 167/C_*.json + 170/K0.json (φ) + 171/M_ONKAYIT.json (α_g).
YENİ KOŞU YOK.

══════════════════════════════════════════════════════════════════════
TÜRETİM (uyumdan ÖNCE)
══════════════════════════════════════════════════════════════════════
θ := KALİB_u2/g_cal ,  KALİB_u2 = Mu2/Pu2 (166_T1.bant_agg),
Mu2 ∝ Re[h̄_Q ⟨e1·x1²·e^{−iWs}⟩] , Pu2 ∝ Re[h̄_Q ⟨E·X²·e^{−iWs}⟩].

EKK artıkları ε := e1 − g_E E (⟨εE⟩ = 0) ve ξ := x1 − g_X X (⟨ξX⟩ = 0)
ile ölçülen üçlü korelatör TAM olarak açılır:

  ⟨e1 x1² e^{−iWs}⟩ = g_E g_X²⟨E X² e^{−iWs}⟩                 [= model]
      + 2g_E g_X⟨E X ξ e^{−iWs}⟩ + g_E⟨E ξ² e^{−iWs}⟩
      + g_X²⟨ε X² e^{−iWs}⟩ + 2g_X⟨ε X ξ e^{−iWs}⟩
      + ⟨ε ξ² e^{−iWs}⟩

  ⇒ **θ = 1 + (beş artık teriminin toplamı)/(g_E g_X²·model)**   (ÖZDEŞ)

Yani θ, g'lerin çözemediği ARTIK ALAN EŞLEŞMESİDİR: her terim ε ya da ξ
eklentisi taşır, en düşük mertebeliler ξ² (∝ varX_res) ve εX² (∝ √varE_res).
Bu, θ'nın taşıyıcısının hangi büyüklükler OLABİLECEĞİNİ belirler ve iki
adayı ÖLÇÜMDEN ÖNCE eler:

 (E1) σ_ds, σ_X̃, σ_Ĉ, W_X — pencere ekseninde bunlar ‰4 içinde DONUK
      (170 §K0 tablosu) ama θ orada %11'e kadar oynuyor ⇒ taşıyamazlar.
 (E2) d/N — λ VE kesim eksenlerinde TAM SABİT (172c Ö4) ⇒ λ eğrisi
      kurulamaz ⇒ aday biçimsel olarak tanımsız.

BİRLEŞİK ADAY (kaptanın G2 formu, β 168'den DONDURULMUŞ):

    **θ/θ₀ = (1 − βφ)·(x/x₀)^p ,  β = 0.2175 SABİT, p YALNIZ λ
    ailesinden (7 gaz, φ = 0) uydurulur; pencere (4) + kesim (4) +
    gerçek gaz (1) ÖRNEKLEM DIŞIDIR.**

══════════════════════════════════════════════════════════════════════
ÖN-MÜHÜR (koşudan önce; 4 Eyl, 172d)
══════════════════════════════════════════════════════════════════════
 Ö1  (E1) adayları pencere ekseninde ölür: örneklem-dışı pencere rms > %7.
 Ö2  (E2) d/N adayı λ ekseninde tanımsız (p sonsuz/NaN) ⇒ ölü.
 Ö3  rE ve rX: λ'da p < 0 (θ↓, r↑) ama pencerede r ARTIYOR ve θ da
     ARTIYOR ⇒ **İŞARET ÇELİŞKİSİ**; pencere öngörüsü θ/θ₀ < 1 verecek,
     ölçüm > 1 ⇒ ÖLÜM. Aynı çelişki 1−kor² için de geçerli.
 Ö4  Q_X üç eksende de İŞARETİ doğru olan tek aday (λ: ikisi de ↓;
     pencere: ikisi de ↑; kesim: ikisi de ↓). Ama 172c'nin tablosundan
     λ ekseni içindeki yerel üslerin 0.09–1.43 arasında oynadığını
     görüyorum ⇒ **ön-kaydım: Q_X da tek güç yasası olarak ölür**
     (λ uyum rms > %1.5, örneklem-dışı rms > %5). Yani beklentim:
     **HİÇBİR TEK DEĞİŞKEN θ'yı üç eksende taşımaz** — 171 §T2b.2'nin
     M için verdiği hükmün θ'ya taşınması.
 Ö5  g_cal adayı SİREN'dir (θ ≡ KALİB/g_cal, tanım gereği bağlı);
     ölçülür ama kimlik sayılmaz.
 Ö6  **τ-EĞİMİ TÜRETİMİ (sıfır parametre).** 171 §T4.2'nin β-artığındaki
     tekdüze τ-eğimi (−15.2 → −141.7 %/τ) θ'nın SEVİYE yasasının borcu
     DEĞİL, gazın kendi A-genişliğinin işidir: A_g(τ) = e^{−α_g τ²}
     ⇒ artık(τ) = 100[(θ̄_g/(θ̄₀(1−βφ)))·e^{−Δα(τ²−τ̄²)} − 1],
     Δα := α_g − α₀, τ̄² := ⟨τ_b²⟩. α_g 171'in gaz-başına uyumundan,
     β 168'den, θ̄ ölçümden — HİÇBİR YENİ PARAMETRE YOK. Ön-kayıt:
     dört kesim gazında öngörülen ve ölçülen τ-eğimi **%20 içinde**.
 Ö7  Aynı türetim λ ekseninde de tutmalı (φ = 0): L060'ın β-artığı
     171'de "DÜZ" görünüyordu; öngörü −12 ± 3 %/τ (Δα = 0.0964),
     ölçüm ≈ −13.5 %/τ olmalı. Tutmazsa Ö6 de şüpheli sayılır.
 Ö8  α_g DENETİMİ: 172d'nin kendi uyumu (log KALİB_b ~ τ_b²) 171'in
     `M_ONKAYIT.genislik.alfa` değerlerini %1 içinde yeniden üretmeli
     (14 gaz); L050/L130 için yeni değerler eklenir.

══════════════════════════════════════════════════════════════════════
SONUÇ (yalnız gerçek koşudan, `172/log_172d.txt`)
══════════════════════════════════════════════════════════════════════
 Ö8 ✓ TAM. α_g 14 gazda 171'in değerlerini **+0.00%** ile yeniden üretti;
      θ, 170 §K0'ın değerlerini **2.2e−16** ile. Yeni: α(L050) = 1.0242,
      α(L130) = 0.6874 (L130 şekil bakımından da aykırı, λ-ailesi
      1.016-1.226 bandının çok altında).
 Ö1 ✗(kısmi) σ-adayları pencerede rms **%6.0-6.2** — öldüler (eşik %5)
      ama ön-kaydım ">%7" demişti, ıskaladı. W_X de aynı yerde (%6.21).
 Ö2 ✓ d/N λ ekseninde donuk ⇒ aday tanımsız.
 Ö3 ✓ rE (%8.74), rX (%12.68), Q_E (%9.17), 1−g_E (%8.30) — hepsi
      pencerede öngörülen İŞARET ÇELİŞKİSİYLE öldü.
 Ö4 ✗ **BENİM ÖN-KAYDIM ÖLDÜ.** "Hiçbir tek değişken taşımaz" dedim;
      **ρ_X = varX_mod/varX_olc AYAKTA KALDI**: p = +0.651 (yalnız λ'dan),
      λ rms %1.37, PENCERE %4.59, KESİM %2.56, son −%1.91.
      Q_X (%17.01 pencere) ise öngördüğüm gibi öldü.
      DÜRÜSTLÜK: %5'lik eşik kaba bir ağdır; θ'nın kendi jk hatası ‰5
      olduğu için ρ_X yasası T/4 gazlarında hâlâ 5-10σ dışındadır
      (gaz-başına artıklar: T2a +0.68, T2b +4.77, T4a +6.25, T4b +4.69;
      kesim: K090 +0.78, K070 −2.90, HA4 −1.23, E060 −3.95).
      **Kimlik değil, en iyi tek-değişkenli taşıyıcı.**
      (μ̂²_X ≡ 0 olduğu için ρ_X ile ρ_var_X aynı sayıdır — özdeşlik
      denetiminin ikinci kez tutması.)
 Ö5 ✓ g_cal SİREN olarak işaretlendi; zaten pencerede %17.87 ile öldü.
 Ö6 ✓✓✓ **TAM İSABET, SIFIR PARAMETRE.** β-artığının τ-eğimi, gazın
      kendi A-genişliğinden (α_g) türetilen −2Δα τ̄ ile ON İKİ gazda
      tutuyor. Kesim: K090 −15.0/−15.2 (−0.8%), K070 −40.1/−33.3
      (+20.2%), HA4 −74.4/−71.6 (+3.9%), E060 −148.9/−141.7 (+5.1%).
      Yalnız K070 ön-kayıtlı %20 sınırının tam üstünde.
 Ö7 ✓ λ ekseninde de tutuyor: L060 öngörü −12.5, ölçüm −11.8 (+6.3%;
      ön-kayıt −12 ± 3). Dahası eğim İŞARET DEĞİŞTİRDİĞİ yerlerde de
      tutuyor: L050 +9.9/+10.6, L130 +49.6/+49.2, son −11.0/−10.7,
      HkT4b +53.3/+54.4. **Ölçülen −149'dan +54'e uzanan bütün τ-eğimi
      yelpazesi TEK bir sayıyla (α_g) açıklanıyor.**
"""
import json
import os

import numpy as np

SCR = ("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
       "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/")
BETA = 0.2175                       # 168 §A2.7 — DONDURULMUŞ
LAM7 = ["L050", "L060", "L070", "L085", "Hkeskin", "L115", "L130"]
KESIM = ["K090", "K070", "HA4", "E060"]
PENC = ["HkT2a", "HkT2b", "HkT4a", "HkT4b"]
HEPSI = LAM7 + ["son"] + KESIM + PENC
LAM = {"L050": .50, "L060": .60, "L070": .70, "L085": .85, "Hkeskin": 1.00,
       "L115": 1.15, "L130": 1.30}

G1 = json.load(open(SCR + "172/G1.json"))
K0 = json.load(open(SCR + "170/K0.json"))["gaz"]
MO = json.load(open(SCR + "171/M_ONKAYIT.json"))["genislik"]


def bantlar(g):
    d = json.load(open(SCR + "167/C_%s.json" % g))
    B = [b for b in d["bant"] if b.get("olculdu")
         and 0.52 - 1e-9 <= b["lo"] <= 0.68 + 1e-9]
    return (np.array([b["KALIB_u2"] for b in B]),
            np.array([b["tau_eff"] for b in B]),
            np.array([b.get("sKALIB_u2", np.nan) for b in B]))


D = {}
for g in HEPSI:
    k, t, sk = bantlar(g)
    alfa = -np.polyfit(t ** 2, np.log(k), 1)[0]
    D[g] = dict(k=k, t=t, sk=sk, KAL=float(np.exp(np.log(k).mean())),
                gcal=G1[g]["E"] and None or None, alfa=float(alfa))
    D[g]["gcal"] = json.load(open(SCR + "167/C_%s.json" % g))["artik"]["gcal"]
    D[g]["th"] = D[g]["KAL"] / D[g]["gcal"]
    D[g]["phi"] = K0[g]["phi"] if g in K0 else 0.0

print("== Ö8: θ SEVİYESİ ve α_g DENETİMİ ==")
print("gaz       λ      KALİB(5b)  g_cal    **θ**     θ/θ₀    φ       "
      "α_g(172d)  α_g(171)  fark%")
th0 = D["Hkeskin"]["th"]
for g in HEPSI:
    a171 = MO[g]["alfa"] if g in MO else float("nan")
    print("%-8s %-6s %.5f    %.5f  %.5f  %.4f  %.4f  %.4f     %s"
          % (g, ("%.2f" % LAM[g]) if g in LAM else "-", D[g]["KAL"],
             D[g]["gcal"], D[g]["th"], D[g]["th"] / th0, D[g]["phi"],
             D[g]["alfa"],
             ("%.4f    %+.2f" % (a171, 100 * (D[g]["alfa"] / a171 - 1)))
             if a171 == a171 else "   —         (yeni)"))
if "theta" in K0.get("Hkeskin", {}):
    sap = max(abs(D[g]["th"] / K0[g]["theta"] - 1) for g in HEPSI if g in K0)
    print("  170/K0 ile θ denetimi: maks fark %.2e" % sap)

# ---------------------------------------------------------------- adaylar
ADAY = {
    "rE": lambda g: G1[g]["E"]["r"], "rX": lambda g: G1[g]["X"]["r"],
    "Q_E": lambda g: G1[g]["E"]["Q"], "Q_X": lambda g: G1[g]["X"]["Q"],
    "rho_E": lambda g: G1[g]["E"]["rho"], "rho_X": lambda g: G1[g]["X"]["rho"],
    "1-gE": lambda g: 1 - G1[g]["gE"], "1-gX": lambda g: 1 - G1[g]["gX"],
    "sig_ds": lambda g: G1[g]["sigds"], "sig_X": lambda g: G1[g]["sigX"],
    "sig_C": lambda g: G1[g]["sigC"], "W_X": lambda g: G1[g]["WX"],
    "d/N": lambda g: 2.0 * G1[g]["nline"] / G1[g]["N"],
    "g_cal(SİREN)": lambda g: D[g]["gcal"],
}
# --- İKİNCİ TUR (ÖN-MÜHÜRSÜZ, etiketli): ρ_X'in ayakta kalmasından SONRA
#     eklendi. Bunlar KİMLİK sayılmaz; yalnız "ρ_X mi kor_X mi" sorusunu
#     açık tutmak için ölçülür (aynı örneklem-dışı eksenler harcanmıştır).
IKINCI = {
    "kor_X^-2*": lambda g: 1.0 / (G1[g]["X"]["pi"] ** 2 / G1[g]["X"]["rho_var"]),
    "kor_E^-2*": lambda g: 1.0 / (G1[g]["E"]["pi"] ** 2 / G1[g]["E"]["rho_var"]),
    "rho_var_X*": lambda g: G1[g]["X"]["rho_var"],
    "rho_E*rho_X²*": lambda g: G1[g]["E"]["rho"] * G1[g]["X"]["rho"] ** 2,
}
print()
print("== Ö1-Ö5: ÇAPRAZ SINAV — θ/θ₀ = (1−βφ)(x/x₀)^p, p YALNIZ λ'dan ==")
print("%-14s %8s | %7s | %7s %7s %7s | %s"
      % ("aday x", "p", "λ rms%", "PENCERE", "KESİM", "son%", "hüküm"))
sonuc = {}
for ad, f in list(ADAY.items()) + [("--- ikinci tur (ön-mühürsüz) ---",
                                    None)] + list(IKINCI.items()):
    if f is None:
        print(ad)
        continue
    x0 = f("Hkeskin")
    xl = np.array([f(g) for g in LAM7])
    yl = np.log(np.array([D[g]["th"] for g in LAM7]) / th0)
    if np.ptp(xl) < 1e-12:
        print("%-14s %8s | %7s | %7s %7s %7s | (E2) λ'da DONUK ⇒ TANIMSIZ"
              % (ad, "—", "—", "—", "—", "—"))
        sonuc[ad] = dict(p=None, olu="tanimsiz")
        continue
    lx = np.log(xl / x0)
    p = float(np.sum(lx * yl) / np.sum(lx * lx))       # kesişimsiz EKK
    rl = 100 * (np.exp(p * lx) - np.exp(yl))
    rmsl = float(np.sqrt((rl ** 2).mean()))

    def ogr(g):
        return (1 - BETA * D[g]["phi"]) * (f(g) / x0) ** p

    rp = np.array([100 * (ogr(g) / (D[g]["th"] / th0) - 1) for g in PENC])
    rk = np.array([100 * (ogr(g) / (D[g]["th"] / th0) - 1) for g in KESIM])
    rs = 100 * (ogr("son") / (D["son"]["th"] / th0) - 1)
    rmsp, rmsk = float(np.sqrt((rp ** 2).mean())), float(np.sqrt((rk ** 2).mean()))
    huk = "AYAKTA" if max(rmsl, rmsp, rmsk) < 5 else "ÖLDÜ"
    print("%-14s %+8.3f | %7.2f | %7.2f %7.2f %+7.2f | %s"
          % (ad, p, rmsl, rmsp, rmsk, rs, huk))
    sonuc[ad] = dict(p=p, lam_rms=rmsl, pen_rms=rmsp, kes_rms=rmsk, son=rs,
                     pen=[float(v) for v in rp], kes=[float(v) for v in rk])

print()
print("== AYAKTA KALANLARIN GAZ-BAŞINA ARTIKLARI (öngörü/ölçüm − 1, %) ==")
for ad in sorted(sonuc, key=lambda a: max(sonuc[a].get("pen_rms", 9e9),
                                          sonuc[a].get("kes_rms", 9e9)))[:3]:
    s = sonuc[ad]
    print("  %-14s p=%+.3f | pencere %s | kesim %s | son %+.2f"
          % (ad, s["p"], "  ".join("%+6.2f" % x for x in s["pen"]),
             "  ".join("%+6.2f" % x for x in s["kes"]), s["son"]))
    print("     (%s)  (%s)" % ("  ".join("%6s" % g for g in PENC),
                               "  ".join("%6s" % g for g in KESIM)))

print()
print("== Ö6/Ö7: β-ARTIĞININ τ-EĞİMİ = A-GENİŞLİĞİ (sıfır parametre) ==")
KH, TH = D["Hkeskin"]["k"], D["Hkeskin"]["t"]
a0 = D["Hkeskin"]["alfa"]
print("%-8s %7s %8s | %-36s | %8s %8s %7s"
      % ("gaz", "φ", "Δα", "ölçülen β-artığı (%)", "eğim ölç", "eğim ÖNG",
         "fark%"))
egim = {}
for g in KESIM + ["L060", "L050", "L130", "son"] + PENC:
    fac = (D[g]["gcal"] / D["Hkeskin"]["gcal"]) * (1 - BETA * D[g]["phi"])
    r = 100 * (D[g]["k"] / (KH * fac) - 1)
    eg = float(np.polyfit(D[g]["t"], r, 1)[0])
    da = D[g]["alfa"] - a0
    tb2 = float((D[g]["t"] ** 2).mean())
    seviye = (D[g]["th"] / th0) / (1 - BETA * D[g]["phi"])
    rp = 100 * (seviye * np.exp(-da * (D[g]["t"] ** 2 - tb2)) - 1)
    egp = float(np.polyfit(D[g]["t"], rp, 1)[0])
    egim[g] = dict(olc=eg, ong=egp, da=float(da))
    print("%-8s %7.4f %+8.4f | %s | %+8.1f %+8.1f %+7.1f"
          % (g, D[g]["phi"], da, "  ".join("%+6.2f" % x for x in r), eg, egp,
             100 * (egp / eg - 1) if abs(eg) > 1e-9 else float("nan")))

os.makedirs(SCR + "172", exist_ok=True)
json.dump(dict(theta={g: D[g]["th"] for g in HEPSI},
               KAL={g: D[g]["KAL"] for g in HEPSI},
               alfa={g: D[g]["alfa"] for g in HEPSI},
               phi={g: D[g]["phi"] for g in HEPSI},
               tau={g: [float(x) for x in D[g]["t"]] for g in HEPSI},
               kband={g: [float(x) for x in D[g]["k"]] for g in HEPSI},
               adaylar=sonuc, egim=egim),
          open(SCR + "172/G2.json", "w"), indent=1)
print("\n-> %s172/G2.json" % SCR)
