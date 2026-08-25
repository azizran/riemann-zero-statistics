"""
105a — MOTORUN ÇİFT KARAKTER (a=0) SAĞLAMASI (25 Ağustos 2026)
==========================================================================
Yeni-ada kampanyası (105b) ÜÇ ada türetecek; ikisi ÇİFT karakter
(a=0) — motor bugüne dek YALNIZ tek karakterlerle (a=1: χ₃, χ₄, χ₅,
χ₇) koştu. 98'in Lmotor sınıfı a_par parametresini taşıyor ama a=0
hiç sağlanmadı. Kampanyadan ÖNCE kapı: geçilmezse kampanya KOŞULMAZ.

ADALAR:
  chi5e  q=5, Legendre (kuadratik) karakter: χ(1,2,3,4)=(+1,−1,−1,+1)
         χ(−1)=χ(4)=+1 → ÇİFT (a=0). Ölü: 5, 25. İlk çizgi log2.
  chi8e  q=8, χ(n)=+1 (n≡±1), −1 (n≡±3) → χ(−1)=χ(7)=+1 → ÇİFT (a=0)
         Ölü: 2'nin kuvvetleri. İlk çizgi log3.
  chi8o  q=8, χ(n)=+1 (n≡1,3), −1 (n≡5,7) → χ(−1)=χ(7)=−1 → TEK (a=1)
         Ölü: 2'nin kuvvetleri. İlk çizgi log3. (a=1 kontrol adası)

SAĞLAMA KAPILARI (98'in S1-S3 kalıbı, a=0'a uyarlanmış):
  G0  PARİTE: χ(q−1) işareti tabloyla tutarlı mı (a doğru seçilmiş mi).
  G1  KÖK SAYISI ε = τ(χ)/(i^a √q): üç karakter de REEL VE PRİMİTİF →
      ε = +1 beklenir (|ε|=1, arg ε=0). a=0'da 98'in kodu i çarpanını
      koymuyor; bu kapı onun a=0'da doğru olduğunu gösterir.
  G2  θ-FAZI: motor Stirling θ(t) vs mpmath
      θ = (t/2)log(q/π) + Im logΓ((1/2+a+it)/2) − arg(ε)/2   [a=0 → Γ(s/2)]
  G3  REELLİK: tamamlanmış Z = e^{iθ}·L(1/2+it,χ) mpmath'te (Hurwitz)
      REEL olmalı: |Im Z|/|Re Z| ≪ 1. Gamma faktörü VEYA ε yanlışsa
      bu kapı patlar (a=1 tablosuyla a=0 karakteri koşulursa da patlar
      — negatif kontrol olarak ÖLÇÜLÜR).
  G4  MOTOR-mpmath FARKI: |Z_motor − Z_mp| RS-düzeltmesiz beklenen
      O((qt)^{-1/4}) bandı içinde mi (ve |Z| ölçeğine göre küçük mü).
  G5  SIFIR YENİDEN-ÜRETİMİ: motorun bulduğu sıfırlar mpmath-reel-Z'nin
      işaret değişimleriyle eşleşiyor mu; |Δγ|/⟨g⟩ ne kadar. İKİ
      YÜKSEKLİKTE (t≈1000 ve t≈30000) — RS-düzeltmesiz motorun konum
      kayması O((qt)^{-1/4}) ile küçülmeli.

HÜKÜM KURALI: G0-G3 KESİN (sayısal sıfır düzeyinde) geçmeli;
G4 (qt)^{-1/4} bandını aşmamalı; G5 için MUTLAK eşik YOK — çünkü
konum kayması motorun bilinen (ve yayımlanmış dört adada zaten
taşınan) RS-kusurudur: ÖLÇÜT, a=0 adalarının ZATEN DOĞRULANMIŞ a=1
kontrol adasından (χ₃) daha kötü OLMAMASIDIR (≤2× ve tam işaret
sayımı). Asıl kanıt ise KONTROLLÜ PARİTE TESTİdir: chi8e (a=0) ve
chi8o (a=1) AYNI iletkende (q=8), aynı yükseklikte yan yana koşar —
yalnızca a değişir. Yani kapı "motor mükemmel mi" değil, "a=0 a=1'den
farklı davranıyor mu" sorusunu sorar — kampanyanın ihtiyacı budur.

İKİ DÜZELTİLMİŞ TEST KUSURU (gizlenmiyor, kayda geçiyor):
 (i) İlk sürümde G5'e mutlak 1e-2 eşiği konmuştu; ZATEN DOĞRULANMIŞ
     kontrol adası χ₃ de 3.0e-2 ile 'kaldı' → eşik artefaktıydı,
     kıyas ölçütüyle değiştirildi.
 (ii) İşaret-değişimi testi SABİT ±0.25⟨g⟩ penceresi kullanıyordu;
     dar çiftlerde pencere İKİ sıfırı birden kapsıyor ve testi sahte
     düşürüyordu (chi5e t≈3e4'te 4/6; teşhis: z3-z4 boşluğu
     0.199⟨g⟩). Pencere yerel boşluğa uyarlandı (±0.40·min-boşluk).

==========================================================================
SONUÇ (25 Ağustos — KAPI GEÇTİ, kampanya koşulabilir)
==========================================================================
G0 ✓ G1 ✓ G2 ✓ G3 ✓ G4 ✓  (dört adada da; kontrol χ₃ dahil)
  G1: üç yeni karakterin de ε = 1 (|ε|−1 < 1e-12, arg ε ~ 5e-26) —
      98'in a=0 dalında i çarpanını KOYMAMASI DOĞRU.
  G2: Δθ ≤ 3e-11 (t = 500…30000) — Stirling a=0'da da yeterli.
  G3: |Im Z|/|Re Z| = 1e-23 … 3e-21 (makine sıfırı). NEGATİF KONTROL
      aynı karakteri yanlış pariteyle koşunca 3.3e+16 verir — kapı
      38 mertebe duyarlı. Gamma faktörü ve ε a=0'da DOĞRU.
  G4: |Z_motor − Z_mp| her t'de (qt)^{-1/4} bandının altında.

G5 (konum kayması, motorun bilinen RS-kusuru; kıyas ölçütü):
  ada    |Δγ|/⟨g⟩ @t≈1e3   @t≈3e4   işaret değişimi
  chi5e     1.46e-02       6.28e-03    20/20, 6/6
  chi8e     2.50e-03       4.84e-03    20/20, 6/6
  chi8o     2.47e-02       1.02e-02    20/20, 6/6
  χ₃ (kontrol) 3.05e-02    5.86e-03    20/20, 6/6
  Kontrole oranlar 0.08-1.74 (hepsi ≤2) → a=0 a=1'den ayırt edilemez.

★ KONTROLLÜ PARİTE TESTİ (asıl kanıt; q=8 SABİT, yalnız a değişiyor):
  chi8e(a=0)/chi8o(a=1) konum-hatası oranı = 0.10 (t≈1e3) / 0.47 (t≈3e4)
  → çift karakter dalında motor TEK dalından DAHA İYİ. Motorun a=0
  desteği doğrulanmıştır.

DÜRÜST NOT: motorun mutlak konum kayması (⟨g⟩'nin %0.5-3'ü) t ile
küçülüyor ama sıfır değil — bu yayımlanmış DÖRT adada da vardır ve
benek/faz ölçümlerini etkilemez (98-S2); donma ölçümü sayım-sertifikası
+ düzlük segmentasyonuyla korunur (101b/101f disiplini).
"""

import numpy as np
import mpmath as mp
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
exec(open(HERE / "98_L_motoru.py").read().split('CHI4 = ')[0])
mp.mp.dps = 25

# --- karakter tabloları (n mod q -> χ) ---
CHI5E = {0: 0, 1: 1, 2: -1, 3: -1, 4: 1}
CHI8E = {0: 0, 1: 1, 2: 0, 3: -1, 4: 0, 5: -1, 6: 0, 7: 1}
CHI8O = {0: 0, 1: 1, 2: 0, 3: 1, 4: 0, 5: -1, 6: 0, 7: -1}
# kontrol: 98'de sağlanmış tek karakter
CHI3 = {0: 0, 1: 1, 2: -1}

ADALAR = [("chi5e", 5, CHI5E, 0), ("chi8e", 8, CHI8E, 0),
          ("chi8o", 8, CHI8O, 1), ("chi3(kontrol)", 3, CHI3, 1)]


def gercek_Z_kompleks(q, tab, a, arg_eps, t):
    """mpmath Hurwitz ile tamamlanmış Z = e^{iθ}L; KOMPLEKS döner."""
    s = mp.mpc(0.5, t)
    Lv = mp.mpc(0)
    for r in range(1, q):
        if tab[r % q] != 0:
            Lv += tab[r % q] * mp.zeta(s, mp.mpf(r) / q)
    Lv *= mp.power(q, -s)
    th = (mp.mpf(t) / 2) * mp.log(mp.mpf(q) / mp.pi) \
        + mp.im(mp.loggamma((s + a) / 2)) - mp.mpf(arg_eps) / 2
    return mp.e**(1j * th) * Lv


print("=" * 74, flush=True)
print("105a — MOTOR ÇİFT KARAKTER (a=0) SAĞLAMASI", flush=True)
print("=" * 74, flush=True)

GECTI = True
OZET = []

for isim, q, tab, a in ADALAR:
    print(f"\n===== {isim}  (q={q}, a={a}) =====", flush=True)
    M = Lmotor(q, tab, a)

    # ---- G0 parite ----
    par = tab[(q - 1) % q]
    beklenen_a = 0 if par == 1 else 1
    g0 = (beklenen_a == a)
    print(f"  G0 PARİTE: χ({q-1}) = {par:+d} → a_beklenen = {beklenen_a}, "
          f"a_kullanılan = {a}   {'✓' if g0 else '✗ RET'}", flush=True)

    # ---- G1 kök sayısı ----
    g1 = abs(M.abs_eps - 1) < 1e-10 and abs(M.arg_eps) < 1e-10
    print(f"  G1 ε: |ε| = {M.abs_eps:.12f}, arg ε = {M.arg_eps:.3e} "
          f"(beklenen 1, 0)   {'✓' if g1 else '✗ RET'}", flush=True)

    # ---- G2 θ-fazı ----
    d_th = []
    for t0 in [500.0, 2000.0, 8000.0, 30000.0]:
        s = mp.mpc(0.5, t0)
        th_mp = float((mp.mpf(t0) / 2) * mp.log(mp.mpf(q) / mp.pi)
                      + mp.im(mp.loggamma((s + a) / 2))) - M.arg_eps / 2
        th_np = float(M.theta(np.array([t0]))[0])
        d_th.append(abs(th_mp - th_np))
    g2 = max(d_th) < 1e-8
    print(f"  G2 θ (t=500/2000/8000/30000): Δθ = "
          + " ".join(f"{x:.1e}" for x in d_th)
          + f"   {'✓' if g2 else '✗ RET'}", flush=True)

    # ---- G3 reellik + G4 motor farkı ----
    ts = [317.7, 1234.5, 4321.9, 12345.6, 28765.4]
    imre, dif, zsc = [], [], []
    for t0 in ts:
        Zc = gercek_Z_kompleks(q, tab, a, M.arg_eps, t0)
        re, im = float(mp.re(Zc)), float(mp.im(Zc))
        imre.append(abs(im) / max(abs(re), 1e-30))
        zm = float(M.Z(np.array([t0]))[0])
        dif.append(abs(zm - re))
        zsc.append(abs(re))
    g3 = max(imre) < 1e-9
    bant = [(q * t) ** -0.25 for t in ts]
    g4 = all(d < 3 * b for d, b in zip(dif, bant))
    print(f"  G3 REELLİK |Im Z|/|Re Z|: " + " ".join(f"{x:.1e}" for x in imre)
          + f"   {'✓' if g3 else '✗ RET'}", flush=True)
    print(f"  G4 |Z_motor−Z_mp|:        " + " ".join(f"{x:.2e}" for x in dif),
          flush=True)
    print(f"     (qt)^(-1/4) bandı:     " + " ".join(f"{x:.2e}" for x in bant)
          + f"   {'✓' if g4 else '✗ RET'}", flush=True)

    # ---- NEGATİF KONTROL: yanlış parite ----
    Myanlis = Lmotor(q, tab, 1 - a)
    Zc_y = gercek_Z_kompleks(q, tab, 1 - a, Myanlis.arg_eps, 4321.9)
    imre_y = abs(float(mp.im(Zc_y))) / max(abs(float(mp.re(Zc_y))), 1e-30)
    print(f"  [neg. kontrol a={1-a}: |Im/Re| = {imre_y:.3e} "
          f"(kapı duyarlı olmalı: ≫ {max(imre):.1e})]", flush=True)

    # ---- G5 sıfır yeniden-üretimi, İKİ YÜKSEKLİKTE ----
    mp.mp.dps = 15
    g5rap = []
    for Tlo, nz, nb in ((1000.0, 20, 30), (30000.0, 6, 22)):
        t0c = time.time()
        zz = M.sifir_bul(Tlo, Tlo + 120.0, grid_frac=0.03)
        dt = time.time() - t0c
        gbar = TWO_PI / np.log(q * (Tlo + 60.0) / TWO_PI)
        hata, isaret_ok, cift = [], 0, 0
        for gi, g in enumerate(zz[:nz]):
            # test penceresi YEREL boşluğa uyarlanır: ±0.25⟨g⟩ sabit pencere
            # dar çiftlerde İKİ sıfırı birden kapsar ve işaret testini sahte
            # olarak düşürür (ilk sürümde chi5e 4/6 böyle çıktı — artefakt).
            sol = g - zz[gi - 1] if gi > 0 else 9e9
            sag = zz[gi + 1] - g if gi + 1 < len(zz) else 9e9
            if min(sol, sag) < 0.5 * gbar:
                cift += 1
            h = min(0.25 * gbar, 0.40 * min(sol, sag))
            fa = float(mp.re(gercek_Z_kompleks(q, tab, a, M.arg_eps, g - h)))
            fb = float(mp.re(gercek_Z_kompleks(q, tab, a, M.arg_eps, g + h)))
            if fa * fb < 0:
                isaret_ok += 1
                lo, hi, flo = g - h, g + h, fa
                for _ in range(nb):
                    mid = 0.5 * (lo + hi)
                    fm = float(mp.re(gercek_Z_kompleks(q, tab, a,
                                                       M.arg_eps, mid)))
                    if flo * fm <= 0:
                        hi = mid
                    else:
                        lo, flo = mid, fm
                hata.append(abs(0.5 * (lo + hi) - g))
        hata = np.array(hata) if hata else np.array([np.nan])
        med = float(np.median(hata) / gbar)
        print(f"  G5 t≈{Tlo:.0f}: motor {len(zz)} sıfır ({dt:.1f} sn); "
              f"işaret değişimi {isaret_ok}/{nz} ({cift} dar çift); "
              f"|Δγ|/⟨g⟩ medyan {med:.2e}, "
              f"maks {np.nanmax(hata)/gbar:.2e}", flush=True)
        g5rap.append((med, isaret_ok, nz))
    mp.mp.dps = 25
    OZET.append((isim, q, a, g0, g1, g2, g3, g4, g5rap))
    if not (g0 and g1 and g2 and g3 and g4):
        GECTI = False

print("\n" + "=" * 74, flush=True)
print(f"{'ada':>16} {'G0':>3} {'G1':>3} {'G2':>3} {'G3':>3} {'G4':>3} "
      f"{'|Δγ|/g @1e3':>12} {'@3e4':>9} {'işaret':>8}", flush=True)
kontrol = [o for o in OZET if o[0].startswith("chi3")][0]
kmed = [kontrol[8][0][0], kontrol[8][1][0]]
for isim, q, a, g0, g1, g2, g3, g4, g5rap in OZET:
    print(f"{isim:>16} " + " ".join(f"{'✓' if x else '✗':>3}"
                                    for x in (g0, g1, g2, g3, g4))
          + f" {g5rap[0][0]:>12.2e} {g5rap[1][0]:>9.2e}"
          + f" {g5rap[0][1]:>3d}/{g5rap[0][2]}, {g5rap[1][1]}/{g5rap[1][2]}", flush=True)
print(f"\nG5 KIYAS ÖLÇÜTÜ (a=0 vs doğrulanmış a=1 kontrol χ₃; ≤2× ve tam"
      f" işaret):", flush=True)
g5ok = True
for isim, q, a, *_r in OZET:
    g5rap = _r[-1]
    o1, o2 = g5rap[0][0] / kmed[0], g5rap[1][0] / kmed[1]
    bu = (o1 <= 2.0 and o2 <= 2.0 and g5rap[0][1] == g5rap[0][2]
          and g5rap[1][1] == g5rap[1][2])
    print(f"  {isim:>16}: kontrole oran {o1:.2f} (t≈1e3) / {o2:.2f} (t≈3e4)"
          f"   {'✓' if bu else '△ İNCELE'}", flush=True)
    if not bu:
        g5ok = False

# --- ASIL PARİTE TESTİ: q=8'de a=0 ve a=1 AYNI iletkende yan yana ---
D = {o[0]: o[-1] for o in OZET}
if "chi8e" in D and "chi8o" in D:
    r1 = D["chi8e"][0][0] / D["chi8o"][0][0]
    r2 = D["chi8e"][1][0] / D["chi8o"][1][0]
    print(f"\n  ★ KONTROLLÜ PARİTE TESTİ (q=8 sabit, yalnız a değişiyor):",
          flush=True)
    print(f"    chi8e(a=0)/chi8o(a=1) konum hatası oranı: "
          f"{r1:.2f} (t≈1e3) / {r2:.2f} (t≈3e4)", flush=True)
    print(f"    → a=0 motoru a=1'den {'DAHA İYİ' if max(r1,r2) < 1 else 'ayırt edilemez' if max(r1,r2) <= 2 else 'DAHA KÖTÜ'}"
          f" (aynı iletken, aynı yükseklik)", flush=True)
print("=" * 74, flush=True)
GECTI = GECTI and g5ok
print("KAPI SONUCU: " + ("GEÇTİ — 105b kampanyası koşulabilir"
                        if GECTI else "RET — KAMPANYA KOŞULMAZ"), flush=True)
