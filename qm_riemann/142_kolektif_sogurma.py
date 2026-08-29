"""
142 — KOLEKTİF SOĞURMA BİRLEŞİMİ: EKRAN ve ⅓-KUYRUK TEK EĞRİ Mİ? (29 Ağu)
==========================================================================
141: yerinde ekran r(τ)≈1.033−0.642τ (τ≤0.37) + dip kalıntısı 0.93.
Lab (eski): örgü koherent kuyruk gücünün ~⅔'ünü soğurur (⅓ hayatta).
GÖZLEM: ekran doğrusunun karesi derin τ'da tam ⅓'e iniyor —
soğurma = ekranın derin ucu OLABİLİR (tek perde hipotezi).
A) Ekran haritası τ≤0.70'e uzatılır (basis q≤e^{0.7L}, ~660 çizgi,
   1300+ sütun chunked regresyon; yüksek-τ çizgiler τ-kutulu,
   Rice-yansızlaştırılmış: A1²−2σ_c²).
B) Varyans defteri: V_zero = Σ_{τ≤0.7} r²w + ρ_tail·Σ_{τ>0.7}w + σ_η²
   → ρ_tail veriden çözülür (σ_η²=0.022 [109]; duyarlılık ±0.004).
ÖN-MÜHÜR:
  S1  ρ_tail ≈ ⟨r_uzatma²⟩_kuyruk (≈0.29) ≈ ⅓ (lab) → TEK PERDE:
      ekran ve soğurma aynı eğri; "⅔-soğurma" ekranın derin ucudur.
  S2  ρ_tail belirgin altında → kuyrukta ekstra kolektif soğurma
      (iki ayrı olgu; dip kalıntısı 0.93'ün büyük kardeşi).
  S3  Harita 0.37-0.70'te doğrusal mı sürer, kırılır mı — kayıt.

SONUÇ (29 Ağustos; ana koşu + üç teşhis koşusu) — S1 DÜŞTÜ, YERİNE
MEKANİZMA MÜHÜRLENDİ: "TEK PERDE" DEĞİL, ÇARPIMSAL ORKESTRA-BAĞLAŞIMI.
  (i)  TABAN-AKIŞI: tek-çizgi genlikler taban-bağımlı — r(2) =
       1.001 (τ≤0.52) → 0.950 (τ≤0.60) → 0.899 (τ≤0.70); BÜTÜN
       çizgiler tekdüze ~0.897 çarpanıyla iniyor (fazlar hep ≤0.002).
       "Bir çizginin genliği" örgüde tek başına iyi-tanımlı değil;
       141'in 0.93 kalıntısı KONVANSİYONA bağlıymış — asıl nesne
       kolektif. (0.756 çıplak-oran konvansiyonsuz ve sağlam.)
  (ii) GRAM ÇAPRAZ-GÜCÜ: V_fit'in %23'ü (0.52 taban) → %37'si (0.70)
       çizgi ÇİFTLERİNE yazılı; ΣA1²/2 ≠ V_fit (0.101 vs 0.160).
  (iii) TEMİZ KUYRUK SPEKTROSKOPİSİ (çözünür bant 0.70-0.78, 822 çizgi,
       ara-nokta referanslı; ilk tasarımın çözünürlük hatası düzeltildi):
       ρ_tail = 0.127 — doğrusal ekran-uzatması (0.31) ve lab ⅓ (0.33)
       ÇOK ALTINDA. Derin kuyruk neredeyse tamamen örgüye yutuluyor.
  (iv) MEKANİZMA (Gram'ın kendisinde): çizgi-çifti bağlaşımı ÇARPIMSAL
       üçlülerde yoğunlaşıyor — q₂=q₁·q₃ (üçü de örgüde) çiftlerinde
       ⟨|G|⟩ kontrol çiftlerinin 9.6 KATI; en güçlü bağlar kule
       basamakları (3↔9, 2↔4, 5↔25, 49↔343 — aracı hep aynı asal).
       141'in kule-güçlenmesi (r(4),r(8),r(16)>1) da aynı makine.
  ⇒ Ekran ve soğurma iki olgu değil; TEK mekanizmanın iki yüzü:
  örgünün kendi modlarına çarpım-tablosu üstünden kurduğu üç-dalga
  bağlaşımı (Hecke dejenerasyonu — MHBH'nin dışladığı, BK95'in
  Hardy–Littlewood'la beslediği yapı). K(τ) çekirdeğinin aday kaynağı
  da bu — bütün kapılar aynı odaya açılıyor. Fazlar süreçte kusursuz
  sıfır: bağlaşım salt genlik/soğurucu, faz-bükümsüz (teori ipucu).
  AÇIK: bağlaşımın nicel teorisi (G'yi aracı-dalga genliğinden öngör),
  ρ_tail(τ) profili, C(n)-çukur genliğinin kolektif hesabı.
"""

import numpy as np
from pathlib import Path
from sympy import primerange, factorint

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi

d = np.load(HERE / "128_odl_zeros6_2e6_zeros.npz")
Z = np.sort(np.asarray(d["zeros"], dtype=float))
zz = Z[len(Z) - 300000:]
g = np.diff(zz)
mid = 0.5 * (zz[:-1] + zz[1:])
Lw = np.log(mid / TWO_PI)
L = float(Lw.mean())
ds = g * Lw / TWO_PI - 1
Nn = len(ds)

TAU_B = 0.70
lim_b = int(np.exp(TAU_B * L))
qs = []
for p in primerange(2, lim_b + 1):
    q = p
    while q <= lim_b:
        qs.append(q); q *= p
qs = sorted(set(qs))
freqs = [np.log(q) for q in qs]
C = 3 + 2 * len(freqs)
print(f"L={L:.3f}  basis: {len(qs)} çizgi (τ≤{TAU_B}), {C} sütun", flush=True)

tt = (mid - mid.mean()) / (mid[-1] - mid[0])
XtX = np.zeros((C, C)); Xty = np.zeros(C)
fr = np.array(freqs)
chunk = 40000
for s0 in range(0, Nn, chunk):
    sl = slice(s0, min(s0 + chunk, Nn))
    arg = np.outer(mid[sl], fr)
    Xc = np.empty((sl.stop - sl.start, C))
    Xc[:, 0] = 1.0; Xc[:, 1] = tt[sl]; Xc[:, 2] = tt[sl]**2
    Xc[:, 3::2] = np.cos(arg); Xc[:, 4::2] = np.sin(arg)
    XtX += Xc.T @ Xc; Xty += Xc.T @ ds[sl]
    del Xc, arg
b = np.linalg.solve(XtX, Xty)

# artık varyansı (σ_c² için)
res_ss = 0.0
for s0 in range(0, Nn, chunk):
    sl = slice(s0, min(s0 + chunk, Nn))
    arg = np.outer(mid[sl], fr)
    fit = (b[0] + b[1]*tt[sl] + b[2]*tt[sl]**2 +
           np.cos(arg) @ b[3::2] + np.sin(arg) @ b[4::2])
    res_ss += float(np.sum((ds[sl] - fit)**2))
    del arg
sig_res2 = res_ss / Nn
sig_c2 = 2 * sig_res2 / Nn
print(f"artık varyans (τ≤0.7 tabanla) = {sig_res2:.4f}  σ_c²={sig_c2:.2e}",
      flush=True)

# A) çizgi başına r; kutulu harita (Rice-yansız)
taus = fr / L
A1sq = b[3::2]**2 + b[4::2]**2 - 2 * sig_c2
pred = np.array([(2 * np.sin(np.pi * t) /
                  (np.pi * factorint(q).popitem()[1] * np.sqrt(q)))**2
                 for q, t in zip(qs, taus)])
print("\nA) ekran haritası (kutulu, yansız):")
print("   τ-kutu    ⟨r²⟩      r      doğru(1.033−0.642τ)")
edges = [0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.62, 0.70]
rmap = []
for lo, hi in zip(edges[:-1], edges[1:]):
    m = (taus >= lo) & (taus < hi)
    if m.sum() == 0:
        continue
    r2 = float(np.sum(A1sq[m]) / np.sum(pred[m]))
    tmid_b = float(np.average(taus[m], weights=pred[m]))
    rmap.append((tmid_b, r2))
    print(f"  {lo:.2f}-{hi:.2f}  {r2:.3f}   {np.sqrt(max(r2,0)):.3f}   "
          f"{1.033-0.642*tmid_b:.3f}", flush=True)

# B) varyans defteri → ρ_tail çözümü
V_zero = float(np.var(ds))
S_in = float(np.sum(A1sq) / 2)          # ölçülen çizgilerin katkısı Σ A1²/2
lim_t = int(np.exp(L))
S_out = 0.0
for p in primerange(2, lim_t + 1):
    q, mm = p, 1
    while q <= lim_t:
        if q > lim_b:
            t = np.log(q) / L
            S_out += 2 * np.sin(np.pi * t)**2 / (np.pi**2 * mm**2 * q)
        q *= p; mm += 1
print(f"\nB) defter: V_zero={V_zero:.4f}  Σ_iç(ölçülen)={S_in:.4f}  "
      f"artık(τ≤0.7 taban)={sig_res2:.4f}  Σ_dış(çıplak)={S_out:.4f}")
tt_ = np.linspace(TAU_B, 1.0, 500)
wt = 2 / np.pi**2 * np.sin(np.pi * tt_)**2 / tt_
rext = np.clip(1.033 - 0.642 * tt_, 0.0, None)
r2_ext_kuyruk = float(np.sum(wt * rext**2) / np.sum(wt))
# artık = ρ_tail·S_out + σ_inc²  (σ_inc: gerçek-inkoherent; 132c'ye göre ≤0.01)
for si in (0.000, 0.005, 0.010):
    rho = (sig_res2 - si) / S_out
    print(f"  σ_inc²={si:.3f} → ρ_tail = {rho:.3f}")
print(f"  karşılaştırma: ekran-uzatma ⟨r²⟩_kuyruk = {r2_ext_kuyruk:.3f}  "
      f"| lab ⅓ = 0.333")
