# -*- coding: utf-8 -*-
"""
172b — G1: g_E'NİN YASASI, ÖLÇÜM KODUNDAN TÜRETİLDİ (4 Eylül 2026)
Veri: 167/C_*.json önbellekleri (17 gaz). YENİ KOŞU YOK, YENİ GAZ YOK.

══════════════════════════════════════════════════════════════════════
TÜRETİM (koda bakarak, kalemle — uyumdan ÖNCE)
══════════════════════════════════════════════════════════════════════
`165_cekirdek.Model165.alanlar` (satır 272-300) tam olarak şunu yapar:

    E   = sentez(s, ω, hp[msk])      hp_q = 2⟨η_{n+1} e^{−iω_q s}⟩
    e1  = Y.e1 − ort(Y.e1)                       (ORTALAMASI TAM SIFIR)
    gE  = np.dot(E,e1) / np.dot(E,E)             (KESİŞİMSİZ EKK eğimi)
    varE_mod = Var(E)   varE_olc = Var(e1)   varE_res = Var(e1−E)
    korE = kor(E, e1)

Ham ikinci momentleri ⟨fg⟩ := (1/N)Σ f g ile yazalım. ort(e1)=0 olduğu
için ⟨E e1⟩ = Kov(E,e1) =: K ve ⟨E²⟩ = Var(E) + μ², μ := ort(E). Varyans
cebiri (hiçbir varsayım yok, ÖZDEŞLİK):

    (Ö-A)  K = ½·[ varE_olc + varE_mod − varE_res ]      = korE·√(V_M V_O)
    (Ö-B)  g_E = K / (varE_mod + μ²)   ⇒   μ̂² := (K/g_E − V_M)/V_O ≥ 0

Yani ÜÇ varyans + korE + g_E arasında İKİ özdeşlik var; ikincisi μ² için
ÇÖZÜLEBİLİR (μ = model alanının DC kaçağı: κ(ω_q) asal rezonansları).

**PROJEKSİYON KUSURU.** S_M := ⟨E²⟩ = K/g_E. Boyutsuzlaştır:
    π := K/V_O        (yakalanan güç payı)
    ρ := S_M/V_O      (model gücü payı, DC dahil)
    **Q := ρ − π = ⟨E·(E−e1)⟩ / V_O**
E gerçek bir DİK İZDÜŞÜM olsaydı (E−e1) ⊥ E olurdu ⇒ **Q ≡ 0** ve
g_E ≡ 1. Demek ki Q, merdiven çizgilerinin DİKGEN OLMAMASININ (Gram
sızıntısı) tek sayılık ölçüsüdür ve

    (Ö-C)  **g_E = 1 − Q/ρ**        (ÖZDEŞLİK, sıfır parametre)
    (Ö-D)  rE := varE_res/varE_olc  ⇒  ρ = 1 − rE + 2Q − μ̂²

Gram diliyle: u_k ∈ {cos ω_q s, sin ω_q s}, Γ = U^T U, v = U^T e1,
E = 2Uv ⇒ K = 2|v|², S_M = 4v^TΓv, Γ = ½(I+Δ) ⇒ g_E = 1/(1+δ),
δ = v^TΔv/|v|², **Q = πδ = 2v^TΔv/V_O**.

ÖLÇEK-DEĞİŞMEZLİĞİ: e1 → αe1 altında v → αv, V_O → α²V_O ⇒ **Q, π, ρ,
g_E'nin HEPSİ DEĞİŞMEZ**. λ ekseni alanı (yaklaşık) tek bir çarpanla
küçültüyorsa Q λ-DEĞİŞMEZ olmalıdır. Bu, türetimin kendi öngörüsüdür.

══════════════════════════════════════════════════════════════════════
ÖN-MÜHÜR (koşudan önce; 4 Eyl, 172b)
══════════════════════════════════════════════════════════════════════
 Ö1  ÖZDEŞLİK DENETİMİ: |K_var − K_kor|/K < 1e-6, 17 gazın hepsinde.
     (Bu bir öngörü değil, defter denetimi; tutmazsa JSON bozuk demektir.)
 Ö2  μ̂² > 0 her gazda (μ² ≥ 0 zorunlu) VE μ̂²/ρ < 0.10 (DC kaçağı küçük
     bir paydır). μ̂² < 0 çıkarsa BÜTÜN resim ölür.
 Ö3  **YASA:** Q, λ merdiveninde (0.50→1.30) DEĞİŞMEZ — yayılım
     (maks−min)/ort ≤ %5. [türetim: ölçek-değişmezliği]
 Ö4  Ö3 doğruysa g_E'nin rE'ye eğimi TÜRETİLİR: dg_E/drE = −Q/ρ²;
     Hkeskin'in sayılarıyla (172a'da zaten basılı) bu **−0.189**'dur.
     172a'nın koklaması **−0.373** demişti. **ÇELİŞKİ ÖNGÖRÜLÜ:** ikisi
     birden doğru olamaz. Eğer ölçülen eğim ≈ −0.373 ise Ö3 ÖLÜR ve o
     zaman ZORUNLU olarak dQ/drE ≈ +2.41 çıkmalıdır (aynı cebirden:
     dg/dr = −q'/ρ + Q(−1+2q')/ρ²). Bu sayı ön-kayıtlıdır; ölçüm
     tutmazsa hem Ö3 hem bu çıkarım ölür, kurtarma yok.
 Ö5  PENCERE EKSENİ: T küçüldükçe gürültü izdüşümü (d/N ∝ 1/T) büyür;
     bu Δ'nın izsel katkısını (σ_n²‖Δ‖_F²/N) artırır ⇒ **Q(T/4) > Q(T/2)
     > Q(T)** (yalnız İŞARET ön-kayıtlı; büyüklük değil).
 Ö6  KESİM EKSENİ: kesim yüksek-τ çizgilerini söndürür; yoğunluk orada
     en büyüktür ⇒ eş-güçlü çizgi çiftleri azalır ⇒ **Q(kesim) < Q(Hk)**
     ve kesim şiddetiyle (φ) tekdüze düşer.
 Ö7  GERÇEK GAZ: `son` λ-ailesinin Q bandının İÇİNDE (172a: g_E doğrusu
     üstünde %0.1-0.3) ⇒ |Q(son) − ⟨Q⟩_λ| ≤ λ-ailesinin yayılımı.
 Ö8  g_X için AYNI cebir aynen geçerli (aynı satırlar, X ve x1 ile).
     Q_X'in λ-yayılımı Q_E'ninkinden KÜÇÜK olmalı (171: g_X λ-değişmez).

══════════════════════════════════════════════════════════════════════
SONUÇ (yalnız gerçek koşudan, 4 Eyl — `172/log_172b.txt`)
══════════════════════════════════════════════════════════════════════
 Ö1 ✓ TAM. |ΔK|/K ≤ 5.2e−15 (16 gazın hepsinde, çoğu 1e−16). Özdeşlik
      MAKİNE HASSASİYETİNDE. ⇒ g_E serbest bir sayı DEĞİL; üç varyans +
      μ² tarafından belirlenir.
 Ö2 ✓ μ̂²_E > 0 her gazda (0.0055…0.0627), μ̂²/ρ ≤ 0.029 (< 0.10 öngörüsü).
      BONUS (öngörülmemişti, koddan çıkıyor): X alanı `alanlar` içinde
      AÇIKÇA ortalanıyor (`X = X − X.mean()`), E ORTALANMIYOR ⇒ μ̂²_X ≡ 0
      olmalı: ölçüm |μ̂²_X| ≤ 1e−17 veriyor. Cebrin kod-düzeyi mührü.
      (Betikteki "μ̂²<0 olan gaz" listesi bu ±0'ları sayıyor — okuma
      hatası, gerçek negatif YOK.)
 Ö3 ✗ **ÖLDÜ.** Q_E λ-yayılımı **%26.8** (öngörü ≤ %5). Q λ-DEĞİŞMEZ
      DEĞİL. Üstelik TEK-YÖNLÜ de değil: Q_E(λ) bir **TÜMSEK** —
      0.6903 (0.50) ↗ 0.9156 (1.00) ↘ 0.8390 (1.30); tepe λ ≈ 1.0-1.06.
 Ö4 ✗ **ÖLDÜ.** Ölçülen eğim −0.3550 (172a'nın koklaması −0.373 ile
      uyumlu), sabit-Q türetimi −0.1893 diyordu: **1.9 kat.** Ön-kayıtlı
      gereklilik dQ/drE = +2.41 idi; ÖLÇÜLEN +1.12. O da tutmadı — çünkü
      ön-kayıt μ̂²'nin de sürüklendiğini ihmal ediyordu (dμ̂²/drE ≈ +0.25).
      Kurtarma yok: iki ön-kayıtlı sayı da ıskaladı; ayakta kalan yalnız
      ÖZDEŞLİK g_E = 1 − Q/ρ'dır.
 Ö5 ✓ **TAM İSABET.** Q_E: T 0.9156 < T/2 {0.9513, 0.9761} < T/4
      {1.0824, 1.0977}. Q_X: 0.4621 < {0.5242, 0.5370} < {0.7152, 0.7379}.
      İşaret ön-kayıtlıydı; büyüklük 172c'de d/N yasasıyla sınanır.
 Ö6 ✗/✓ K090 TERS gitti (0.9411 > 0.9156, +%2.8); K070/HA4/E060 tekdüze
      düştü (0.8305/0.7959/0.6548). φ tek değişken değil (E060 φ=0.960 <
      K070 φ=0.927 sıralamasına rağmen Q'da en düşük).
 Ö7 ✓ Q_E(son) = 0.8758, λ-ailesi [0.6903, 0.9156] — içinde; g_E–rE
      doğrusunda +%0.45 (ailenin kendi artık rms'i %0.45).
 Ö8 ✓ Q_X λ-yayılımı %12.5 < Q_E'nin %26.8'i. Ve Q_X λ'da TEKDÜZE
      AZALAN (0.4920→0.4363), Q_E gibi tümsek DEĞİL.
"""
import json
import math
import os
from math import erfc

import numpy as np

SCR = ("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
       "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/")
K167 = SCR + "167/"
OUT = SCR + "172/"

LAM7 = ["L050", "L060", "L070", "L085", "Hkeskin", "L115", "L130"]
KESIM = ["K090", "K070", "HA4", "E060"]
PENC = ["HkT2a", "HkT2b", "HkT4a", "HkT4b"]
HEPSI = LAM7 + ["son"] + KESIM + PENC
LAM = {"L050": .50, "L060": .60, "L070": .70, "L085": .85, "Hkeskin": 1.00,
       "L115": 1.15, "L130": 1.30}


def oku(g):
    d = json.load(open(K167 + "C_%s.json" % g))
    a = d["artik"]
    orta = [b for b in d["bant"] if abs(b["lo"] - 0.60) < 1e-9][0]
    r = dict(a)
    r.update(N=d["N"], T=d["T"], L=d["L"], sigX=d["sigX"], dres=d["dres"],
             KAL=orta["KALIB_u2"], WX=orta["W_X"], Wamp=orta["W_amp"],
             th=orta["KALIB_u2"] / a["gcal"], pen=d["pen"],
             tau_ust=d["tau_ust"] or 1.0, lam=d["lam"])
    return r


def kusur(V_M, V_O, V_R, kor, g):
    """Türetimin (Ö-A..Ö-D) tek fonksiyonu."""
    Kv = 0.5 * (V_O + V_M - V_R)          # Kov(E,e1) — varyans yolundan
    Kk = kor * math.sqrt(V_M * V_O)       # Kov — korelasyon yolundan
    S = Kv / g                            # ⟨E²⟩
    return dict(K=Kv, dK=abs(Kv - Kk) / abs(Kv), S=S,
                mu2=(S - V_M) / V_O, pi=Kv / V_O, rho=S / V_O,
                Q=(S - Kv) / V_O, r=V_R / V_O, rho_var=V_M / V_O)


def main():
    V = {g: oku(g) for g in HEPSI}
    # φ (168_olcek ile aynı tanım) — merdiven τ'ları C_*.json'da yok, ama
    # w_q yalnız kesim biçimine bağlı; τ ızgarasını Hkeskin'in çizgilerinden
    # değil merdivenden almak gerekir ⇒ φ'yi rapordan DEĞİL, aşağıda
    # 172c'de merdivenden hesaplayacağız. Burada yalnız etiket tutuyoruz.
    R = {}
    print("== Ö1/Ö2: ÖZDEŞLİK DENETİMİ ve DC KAÇAĞI ==")
    print("gaz       nline   N      |ΔK|/K     g_E     g_E*(=π/ρ_var)  "
          "μ̂²      μ̂²/ρ    | g_X     |ΔK_X|/K   μ̂²_X")
    kotu = 0
    for g in HEPSI:
        v = V[g]
        E = kusur(v["varE_mod"], v["varE_olc"], v["varE_res"], v["korE"],
                  v["gE"])
        X = kusur(v["varX_mod"], v["varX_olc"], v["varX_res"], v["korX"],
                  v["gX"])
        R[g] = dict(E=E, X=X, lam=LAM.get(g), N=v["N"], T=v["T"],
                    nline=v["nline"], th=v["th"], gE=v["gE"], gX=v["gX"],
                    KAL=v["KAL"], WX=v["WX"], sigX=v["sigX"],
                    sigds=v["sigds"], sigC=v["sigC"], gcal=v["gcal"])
        gs = E["pi"] / E["rho_var"]
        if E["dK"] > 1e-6 or X["dK"] > 1e-6:
            kotu += 1
        print("%-8s %5d %7d %.2e   %.4f  %.4f        %+.5f %+.4f | "
              "%.4f  %.2e  %+.5f"
              % (g, v["nline"], v["N"], E["dK"], v["gE"], gs, E["mu2"],
                 E["mu2"] / E["rho"], v["gX"], X["dK"], X["mu2"]))
    print("  Ö1: %d/%d gazda özdeşlik 1e-6'yı aştı" % (kotu, len(HEPSI)))
    neg = [g for g in HEPSI if R[g]["E"]["mu2"] < 0 or R[g]["X"]["mu2"] < 0]
    print("  Ö2: μ̂² < 0 olan gaz: %s" % (neg or "YOK"))

    print()
    print("== Ö3/Ö8: PROJEKSİYON KUSURU Q — BÜTÜN EKSENLER ==")
    print("gaz       λ      rE      ρ_E     **Q_E**   g_E(=1−Q/ρ) | "
          "rX      ρ_X     **Q_X**   g_X")
    for g in HEPSI:
        E, X = R[g]["E"], R[g]["X"]
        print("%-8s %-6s %.4f  %.4f  %.5f  %.4f       | %.4f  %.4f  "
              "%.5f  %.4f"
              % (g, ("%.2f" % LAM[g]) if g in LAM else "-", E["r"], E["rho"],
                 E["Q"], 1 - E["Q"] / E["rho"], X["r"], X["rho"], X["Q"],
                 1 - X["Q"] / X["rho"]))

    for et, ky in (("Q_E", "E"), ("Q_X", "X")):
        q = np.array([R[g][ky]["Q"] for g in LAM7])
        print("  %s λ-merdiveni: ort %.5f, yayılım (maks−min)/ort = %.2f%%"
              % (et, q.mean(), 100 * (q.max() - q.min()) / q.mean()))
    qk = np.array([R[g]["E"]["Q"] for g in KESIM])
    qp = np.array([R[g]["E"]["Q"] for g in PENC])
    print("  Ö6 kesim Q_E: %s (Hk = %.5f)"
          % (np.array2string(qk, precision=5), R["Hkeskin"]["E"]["Q"]))
    print("  Ö5 pencere Q_E: T/2 %.5f %.5f  T/4 %.5f %.5f  (T %.5f)"
          % (R["HkT2a"]["E"]["Q"], R["HkT2b"]["E"]["Q"], R["HkT4a"]["E"]["Q"],
             R["HkT4b"]["E"]["Q"], R["Hkeskin"]["E"]["Q"]))
    print("  Ö7 son: Q_E = %.5f, λ-ailesi [%.5f, %.5f]"
          % (R["son"]["E"]["Q"], min(R[g]["E"]["Q"] for g in LAM7),
             max(R[g]["E"]["Q"] for g in LAM7)))

    print()
    print("== Ö4: g_E–rE EĞİMİ — TÜRETİM vs ÖLÇÜM ==")
    rr = np.array([R[g]["E"]["r"] for g in LAM7])
    gg = np.array([R[g]["gE"] for g in LAM7])
    egim, kes = np.polyfit(rr, gg, 1)
    art = gg - (egim * rr + kes)
    H = R["Hkeskin"]["E"]
    tur = -H["Q"] / H["rho"] ** 2
    print("  ölçülen (7 λ gazı):   g_E = %+.5f %+.5f·rE   (artık rms %.2e)"
          % (kes, egim, art.std()))
    print("  TÜRETİM (sabit Q):    dg_E/drE = −Q/ρ² = %+.5f  (Hkeskin'de)"
          % tur)
    # Ö4'ün ikinci yarısı: gerekli dQ/drE
    qq = np.array([R[g]["E"]["Q"] for g in LAM7])
    qeg = np.polyfit(rr, qq, 1)[0]
    gerek = (egim + H["Q"] / H["rho"] ** 2) / (
        -1.0 / H["rho"] + 2 * H["Q"] / H["rho"] ** 2)
    print("  ön-kayıtlı gereklilik: dQ/drE = %+.3f   ÖLÇÜLEN dQ/drE = %+.3f"
          % (2.41, qeg))
    print("  (aynı cebirin ölçülen eğimden çözdüğü dQ/drE = %+.3f)" % gerek)
    for g in ["son"] + KESIM + PENC:
        E = R[g]["E"]
        print("    %-8s rE=%.4f  g_E ölç %.4f  doğru %.4f  fark %+.2f%%"
              % (g, E["r"], R[g]["gE"], egim * E["r"] + kes,
                 100 * ((egim * E["r"] + kes) / R[g]["gE"] - 1)))

    os.makedirs(OUT, exist_ok=True)
    json.dump({g: {"lam": R[g]["lam"], "N": R[g]["N"], "T": R[g]["T"],
                   "nline": R[g]["nline"], "gE": R[g]["gE"],
                   "gX": R[g]["gX"], "th": R[g]["th"], "KAL": R[g]["KAL"],
                   "WX": R[g]["WX"], "gcal": R[g]["gcal"],
                   "sigX": R[g]["sigX"], "sigds": R[g]["sigds"],
                   "sigC": R[g]["sigC"],
                   "E": R[g]["E"], "X": R[g]["X"]} for g in HEPSI},
              open(OUT + "G1.json", "w"), indent=1)
    print("\n-> %sG1.json" % OUT)


if __name__ == "__main__":
    main()
