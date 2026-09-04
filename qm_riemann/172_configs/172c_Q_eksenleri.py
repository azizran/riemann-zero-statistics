# -*- coding: utf-8 -*-
"""
172c — G1 KAPANIŞI + G5: PROJEKSİYON KUSURU Q'NUN EKSEN YASALARI
Veri: 172/G1.json (172b'nin ürettiği). YENİ KOŞU YOK.

══════════════════════════════════════════════════════════════════════
TÜRETİM (172b'nin Gram cebrinin devamı — uyumdan ÖNCE)
══════════════════════════════════════════════════════════════════════
Q = 2 v^TΔv / V_O ,  v = U^T e1 ,  Δ = 2Γ − I  (Γ = Gram, U = merdiven
cos/sin tabanı, d = 2·nline boyut). e1 = f + n (f = koherent çizgi
içeriği, n = geri kalan) ayrıştırmasıyla İKİ katkı:

  (i)  KOHERENT SIZINTI    Q_koh = 2 v_f^TΔv_f/V_O
  (ii) GÜRÜLTÜ İZDÜŞÜMÜ    E[2n^TUΔU^Tn]/V_O = (σ_n²/V_O)·(‖Δ‖_F²/d)·(d/N)

(ii) 4Γ² = (I+Δ)² açılımından: E[⟨E,E⟩_n] = (σ_n²/N)(d + 2trΔ + trΔ²),
E[⟨E,e1⟩_n] = σ_n² d/N, trΔ ≈ 0 ⇒ fark = σ_n²‖Δ‖_F²/N.

⇒ **Q = Q_koh + B·(d/N)** ,  B = (σ_n²/V_O)·(‖Δ‖_F²/d).

İKİ MEKANİZMA, İKİ AYIRT EDİCİ ORAN. Δ YALNIZ tabana (merdiven + s
ızgarası) bağlıdır; E ve X AYNI tabanı kullanır ⇒ ‖Δ‖_F²/d ORTAKTIR:

  (M1) SAF GÜRÜLTÜ-AŞIRI UYUMU:  B_E/B_X = (σ_n,E²/V_O,E)/(σ_n,X²/V_O,X)
       ≈ (1−kor_E²)/(1−kor_X²)  = **1.025**  (Hkeskin'in kendi sayıları)
  (M2) SAF KOHERENT κ-GENİŞLEMESİ (pencere kısalınca κ 1/T genişler ⇒
       çözünmemiş çizgi çifti sayısı ∝ 1/T ∝ d/N): bu durumda B ∝ Q_koh
       ⇒ B_E/B_X = Q_E/Q_X (Hkeskin) = **1.981**, ve dahası Q'nun TAMAMI
       1/T ile ölçeklenmeli: Q(T/4)/Q(T) = **4.00**.

══════════════════════════════════════════════════════════════════════
ÖN-MÜHÜR (koşudan önce; 4 Eyl, 172c)
══════════════════════════════════════════════════════════════════════
 Ö1  BİÇİM: Q, pencere ekseninin BEŞ noktasında (Hk, T2a, T2b, T4a, T4b)
     d/N'de DOĞRUSAL; artık rms ≤ Q'nun %1'i. (5 nokta, 2 parametre ⇒
     3 serbestlik; gerçek bir sınav.) Hem E hem X için.
 Ö2  KESİŞİM: Q_koh = Q(d/N → 0) > 0 ve Q(Hk)'nin ALTINDA (aşırı uyum
     terimi pozitif). E: Q_koh ∈ (0.50, 0.9156); X: Q_koh ∈ (0.20, 0.4621).
 Ö3  MEKANİZMA HAKEMİ: B_E/B_X = 1.025 (M1) ya da 1.981 (M2).
     DÜRÜSTLÜK NOTU: 172b'nin basılı tablosundan iki-nokta kaba eğimle
     bu oranın ≈ 0.66 çıktığını görüyorum ⇒ **İKİ MEKANİZMA DA ÖLECEK**
     diye bekliyorum; ön-kayıt yine de aynen sınanır ve kurtarılmaz.
     Ayrıca M2'nin ikinci ayağı Q(T/4)/Q(T) = 4.00 ölçülür.
 Ö4  λ EKSENİ TEMİZ: λ ailesinde d/N SABİT (nline = 8981, N = 299998)
     ⇒ Q'nun λ hareketi TAMAMEN Q_koh'tur, aşırı uyumla ilgisi yoktur.
     (Kesim ekseni de aynı: d/N sabit.) Bu bir öngörü değil, ayrıştırma.
 Ö5  TÜMSEĞİN YERİ: log Q_E(λ)'ye beş orta noktadan (0.60…1.15) parabol;
     tepe **λ* = λ_c = 1.0109** (170 §K0.2'nin inşa doyum eşiği, SIFIR
     yeni parametre) ± 0.05. Yedi noktalı parabol da aynı yeri vermeli.
 Ö6  G5 — g_X TEK ÇERÇEVE: g_X = 1 − Q_X/ρ_X ÖZDEŞ. λ = 0.50 uyanışı
     (Δlog g_X² = +0.0157) ve T/4 sızıntısı (g_X² −%22) AYNI Q_X
     değişkeninde okunmalı; ikisinin de Q_X/ρ_X'teki payı hesaplanır.
     Ön-kayıt: λ = 0.50 adımında g_X'in artışı Q_X'ten DEĞİL ρ_X'ten
     gelir (Q_X 0.60→0.50'de ‰7 oynuyor, ρ_X %1.3), yani "uyanış" bir
     MODEL GÜCÜ olayıdır, sızıntı olayı değildir.
 Ö7  GERÇEK GAZ: Q_E(son) ve Q_X(son), λ-ailesinin Q(λ) eğrisinden
     okunan λ_eş değerleri TUTARLI olmalı (|λ_eş^E − λ_eş^X| ≤ 0.10);
     tutmazsa gerçek gaz "tek bir eşdeğer λ" ile tarif edilemez.

══════════════════════════════════════════════════════════════════════
SONUÇ (yalnız gerçek koşudan, `172/log_172c.txt`)
══════════════════════════════════════════════════════════════════════
 Ö1 ✓/✗ E için artık rms %0.88 (≤%1 ✓); X için %2.40 (✗). Doğrusal biçim
      E'de tutuyor, X'te tutmuyor (X'in eğrisi d/N'de hafif konveks).
 Ö2 ✓ Q_koh(E) = 0.8439 ∈ (0.50, 0.9156); Q_koh(X) = 0.3504 ∈ (0.20,
      0.4621). Aşırı uyum terimi pozitif ve Q(Hk)'nin sırasıyla %8 ve
      %24'ü. B_E = 1.052, B_X = 1.603.
 Ö3 ✗✗ **İKİ MEKANİZMA DA ÖLDÜ** (beklendiği gibi, kurtarmasız):
      B_E/B_X = 0.656 (M1: 1.025, M2: 1.981); M2'nin ikinci ayağı
      Q(T/4)/Q(T) = 1.191 (E) / 1.572 (X), öngörü 4.00 — çok uzak.
      Yani T ekseni ne saf gürültü-aşırı uyumu ne de saf κ-genişlemesi.
      KOKLAMA (ön-mühürsüz, etiketli): B·r ≈ 0.59 iki alanda da
      (1.052·0.575 = 0.605, 1.603·0.362 = 0.580) — B ∝ 1/artık payı.
      Bu bir iddia DEĞİL, 173'e not.
 Ö4 ✓ d/N = 0.059874 λ ailesinde ve kesim ailesinde TAM SABİT (nline
      8981, N 299998 hepsinde). ⇒ **λ ve kesim eksenlerinde Q'nun bütün
      hareketi KOHERENT SIZINTIDIR**; aşırı uyum yalnız T ekseninde.
 Ö5 ✓✓ **TAM İSABET (sıfır yeni parametre).** Q_E tümseğinin tepesi
      λ* = **1.0191** (5 orta nokta) ve **1.0009** (7 nokta); ön-kayıt
      λ_c = 1.0109 ± 0.05 — ikisi de |Δ| ≤ 0.011 içinde. **Projeksiyon
      kusuru, inşa denkleminin doyum eşiğinde ZİRVE yapıyor.**
      Q_X'in parabol tepesi 0.24-0.30, yani ölçüm aralığının dışında:
      Q_X ölçülen bütün λ'larda TEKDÜZE azalıyor (Ö8/172b ile uyumlu).
 Ö6 ✓/kısmi (G5): λ=0.60→0.50 adımında Δlog g_X² = +0.01566'nın **%66'sı
      ρ_X'ten** (+0.01253), %34'ü Q_X'ten (−0.00651). T/4'te ise TERS:
      Δlog g_X² = −0.24419'un baskın kaynağı **Q_X** (+0.46795), ρ_X
      (+0.22606) onu kısmen telafi ediyor. ⇒ **HÜKÜM: λ=0.50 uyanışı ile
      T/4 sızıntısı AYNI ODA DEĞİL.** Biri model gücü (ρ_X), öteki
      projeksiyon kusuru (Q_X).
 Ö7 ✗ **ÖLDÜ.** λ_eş(Q_E) = 0.7845, λ_eş(Q_X) = 0.9463 — fark 0.162
      (eşik 0.10). **Gerçek gaz tek bir eşdeğer λ ile tarif edilemez**;
      E ve X kanalları farklı λ'lara işaret ediyor.
"""
import json
import numpy as np

SCR = ("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
       "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/172/")
G = json.load(open(SCR + "G1.json"))

LAM7 = ["L050", "L060", "L070", "L085", "Hkeskin", "L115", "L130"]
PENC = ["Hkeskin", "HkT2a", "HkT2b", "HkT4a", "HkT4b"]
KESIM = ["K090", "K070", "HA4", "E060"]
LAMV = np.array([.50, .60, .70, .85, 1.00, 1.15, 1.30])
LAM_C = 1.0109


def dN(g):
    return 2.0 * G[g]["nline"] / G[g]["N"]


print("== Ö1/Ö2/Ö3: PENCERE EKSENİ — Q = Q_koh + B·(d/N) ==")
x = np.array([dN(g) for g in PENC])
sonuc = {}
for ky, ad in (("E", "Q_E"), ("X", "Q_X")):
    y = np.array([G[g][ky]["Q"] for g in PENC])
    B, Q0 = np.polyfit(x, y, 1)
    art = y - (B * x + Q0)
    print("  %s: d/N = %s" % (ad, np.array2string(x, precision=5)))
    print("     Q  = %s" % np.array2string(y, precision=5))
    print("     ** Q_koh = %+.5f,  B = %+.5f **  artık rms = %.5f "
          "(= Q(Hk)'nin %%%.2f'si)" % (Q0, B, art.std(), 100 * art.std() / y[0]))
    print("     artıklar: %s" % np.array2string(art, precision=5))
    sonuc[ad] = dict(Q0=Q0, B=B, rms=float(art.std()))
oran = sonuc["Q_E"]["B"] / sonuc["Q_X"]["B"]
print("  ÖLÇÜLEN B_E/B_X = %.3f   |  M1 öngörü 1.025   M2 öngörü 1.981"
      % oran)
for ky in ("E", "X"):
    q4 = 0.5 * (G["HkT4a"][ky]["Q"] + G["HkT4b"][ky]["Q"])
    print("  M2 ikinci ayak %s: Q(T/4)/Q(T) = %.3f  (öngörü 4.00)"
          % (ky, q4 / G["Hkeskin"][ky]["Q"]))

print()
print("== Ö4/Ö5: λ EKSENİ — d/N SABİT mi, TÜMSEK NEREDE ==")
print("  d/N: %s" % np.array2string(np.array([dN(g) for g in LAM7]),
                                    precision=6))
print("  kesim d/N: %s" % np.array2string(np.array([dN(g) for g in KESIM]),
                                          precision=6))
qE = np.array([G[g]["E"]["Q"] for g in LAM7])
qX = np.array([G[g]["X"]["Q"] for g in LAM7])
for ad, q in (("Q_E", qE), ("Q_X", qX)):
    for et, sl in (("5 orta (0.60-1.15)", slice(1, 6)), ("7 nokta", slice(None))):
        c = np.polyfit(LAMV[sl], np.log(q[sl]), 2)
        tepe = -c[1] / (2 * c[0]) if c[0] != 0 else float("nan")
        print("  %s %-20s parabol tepesi λ* = %+.4f  (a2 = %+.4f)"
              % (ad, et, tepe, c[0]))
print("  ÖN-KAYIT λ_c = %.4f (170 §K0.2)" % LAM_C)

print()
print("== Ö6 (G5): g_X TEK ÇERÇEVEDE ==")
print("gaz       λ      Q_X      ρ_X      Q_X/ρ_X   g_X      Δlog g_X   "
       "  pay Q_X    pay ρ_X")
ref = G["Hkeskin"]["X"]
for g in LAM7 + ["son"] + KESIM + ["HkT2a", "HkT4b"]:
    X = G[g]["X"]
    dg = np.log(G[g]["gX"] / G["Hkeskin"]["gX"])
    # d log g = -(Q/ρ)/(1-Q/ρ) * (dlogQ - dlogρ)  → payları ayır
    f = (X["Q"] / X["rho"])
    f0 = ref["Q"] / ref["rho"]
    payQ = -(f0 / (1 - f0)) * np.log(X["Q"] / ref["Q"])
    payR = +(f0 / (1 - f0)) * np.log(X["rho"] / ref["rho"])
    print("%-8s %-6s %.5f  %.5f  %.5f  %.5f  %+.5f    %+.5f   %+.5f"
          % (g, ("%.2f" % G[g]["lam"]) if G[g]["lam"] else "-", X["Q"],
             X["rho"], f, G[g]["gX"], dg, payQ, payR))
print("  λ=0.60→0.50 adımı: Δlog Q_X = %+.5f, Δlog ρ_X = %+.5f, "
      "Δlog g_X² = %+.5f"
      % (np.log(G["L050"]["X"]["Q"] / G["L060"]["X"]["Q"]),
         np.log(G["L050"]["X"]["rho"] / G["L060"]["X"]["rho"]),
         2 * np.log(G["L050"]["gX"] / G["L060"]["gX"])))
print("  T/4 (HkT4b vs Hk):  Δlog Q_X = %+.5f, Δlog ρ_X = %+.5f, "
      "Δlog g_X² = %+.5f"
      % (np.log(G["HkT4b"]["X"]["Q"] / ref["Q"]),
         np.log(G["HkT4b"]["X"]["rho"] / ref["rho"]),
         2 * np.log(G["HkT4b"]["gX"] / G["Hkeskin"]["gX"])))

print()
print("== Ö7: GERÇEK GAZIN EŞDEĞER λ'sı (Q'dan) ==")
for ad, q in (("Q_E", qE), ("Q_X", qX)):
    v = G["son"]["E" if ad == "Q_E" else "X"]["Q"]
    qq, ll = q[:5], LAMV[:5]              # λ ≤ 1 kolu (tekdüze kol)
    if qq[0] > qq[-1]:                    # azalan ⇒ np.interp için çevir
        qq, ll = qq[::-1], ll[::-1]
    le = float(np.interp(v, qq, ll))
    print("  %s(son) = %.5f  ⇒  λ_eş = %.4f  (λ ≤ 1 kolunda, %s)"
          % (ad, v, le, "artan" if q[0] < q[-1] else "azalan"))

print()
print("== EK (ön-mühürsüz, Ö4'ün ölümünden SONRA): g_E–rE doğrusunun "
      "gerçek üyeleri ve L130 ==")
rE = np.array([G[g]["E"]["r"] for g in LAM7])
gE = np.array([G[g]["gE"] for g in LAM7])
for et, sl in (("7 λ gazı", slice(None)), ("L130 HARİÇ (6 gaz)", slice(0, 6))):
    a, b = np.polyfit(rE[sl], gE[sl], 1)
    art = gE[sl] - (a * rE[sl] + b)
    print("  %-20s g_E = %+.5f %+.5f·rE   artık rms %.2e" % (et, b, a, art.std()))
    if sl.stop == 6:
        for g in LAM7 + ["son"]:
            E = G[g]["E"]
            print("      %-8s fark %+.2f%%"
                  % (g, 100 * ((a * E["r"] + b) / G[g]["gE"] - 1)))
print("  L115→L130 yerel eğimler: ΔrE = %+.5f, ΔQ = %+.5f ⇒ dQ/drE = %+.2f"
      "  (Δg_E = %+.5f ⇒ dg/dr = %+.2f, İŞARET DÖNÜYOR)"
      % (rE[6] - rE[5], G["L130"]["E"]["Q"] - G["L115"]["E"]["Q"],
         (G["L130"]["E"]["Q"] - G["L115"]["E"]["Q"]) / (rE[6] - rE[5]),
         gE[6] - gE[5], (gE[6] - gE[5]) / (rE[6] - rE[5])))
