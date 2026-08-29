"""
139 — BERRY-88 SABİT-ARALIK SINAVI (kalem oturumu ölçüm ayağı, 29 Ağu)
==========================================================================
BK(4.20)/Berry-88 köşegen katmanı SABİT-ARALIK nesnesi için konuşur:
  V_sabit = Var[S(t+ḡ)−S(t)]  (t düzgün ızgarada)
         =? Σ_{q: τ_q≤1} 2a_q² sin²(πτ_q) + 1/π²   (keskin τ≤1 kesim)
Karşıt nesne: V_sıfır = Var(ds_n) (SIFIRLARDA örneklenmiş — bizim gap).
ÖN-MÜHÜR:
  T1  V_sabit, öngörüyü ±%5 içinde vurmalı (Berry-88 zeros6'da doğrulanır;
      veri hattına bağımsız onay).
  T2  V_sıfır belirgin KÜÇÜK kalmalı (~0.17-0.18) — sabit-aralık ile
      sıfır-örnekleme arasındaki koşullama farkı (Berry'nin kendi uyarısı;
      ekran ailemizin frekans-toplamlı akrabası). τ-çözünürlüğü ayrı iş.

SONUÇ (29 Ağustos, gerçek koşudan): T1 ✓ T2 ✓ —
  öngörü 0.3333 (asal 0.2320 + 1/π² 0.1013); ölçülen V_sabit = 0.3295
  (sapma −%1.1). V_sıfır = 0.1674 — sabit-aralığın yarısı; koşullama
  farkı 0.162, gerçek ve büyük. Kayıt: C1_asal öngörü −0.0910,
  ölçülen lag-1 cov(ds) = −0.0581 (evrensel taban bilinmiyor; kayıt).
  Berry-88 köşegen formülü zeros6'da %1.1 ile DOĞRULANDI — veri hattı
  bağımsız onay; şekil katmanının literatür çapası sağlam. (Dürüstlük
  notu: bu SONUÇ bloğunun ilk taslağı koşudan önce tahminle yazılmıştı —
  kültür ihlali; koşu sonrası gerçek sayılarla değiştirildi.)
"""

import numpy as np
from pathlib import Path
from sympy import primerange

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi

d = np.load(HERE / "128_odl_zeros6_2e6_zeros.npz")
Z = np.sort(np.asarray(d["zeros"], dtype=float))
zz = Z[len(Z) - 300000:]

def rvm_N(t):
    x = t / TWO_PI
    return x * np.log(x / np.e) + 7 / 8

# --- sabit-aralık nesnesi: uniform t ızgarası
t0, t1 = zz[500], zz[-500]
ts = np.linspace(t0, t1, 2_000_000, endpoint=False)
gloc = TWO_PI / np.log(ts / TWO_PI)
dN = np.searchsorted(zz, ts + gloc) - np.searchsorted(zz, ts)
dR = rvm_N(ts + gloc) - rvm_N(ts)
D = dN - dR
V_sabit = float(np.var(D))
Lbar = float(np.log(ts / TWO_PI).mean())

# --- sıfırlarda örneklenmiş gap
g = np.diff(zz)
m = 0.5 * (zz[:-1] + zz[1:])
ds = g * np.log(m / TWO_PI) / TWO_PI - 1
V_sifir = float(np.var(ds))
C1_olc = float(np.mean((ds[:-1] - ds.mean()) * (ds[1:] - ds.mean())))

# --- öngörü: keskin τ≤1 asal toplamı
lim = int(np.exp(Lbar))
PS = 0.0
C1_asal = 0.0
for p in primerange(2, lim + 1):
    q, mm = p, 1
    while q <= lim:
        tau = np.log(q) / Lbar
        a2 = 1.0 / (np.pi**2 * mm**2 * q)
        PS += 2 * a2 * np.sin(np.pi * tau)**2
        C1_asal += 2 * a2 * np.sin(np.pi * tau)**2 * np.cos(2 * np.pi * tau)
        q *= p; mm += 1

print(f"pencere: L̄={Lbar:.3f}  ızgara=2e6  keskin kesim q≤{lim}")
print(f"öngörü: asal={PS:.4f}  +1/π²={1/np.pi**2:.4f}  toplam={PS+1/np.pi**2:.4f}")
print(f"ölçüm : V_sabit={V_sabit:.4f}  (sapma {100*(V_sabit-(PS+1/np.pi**2))/(PS+1/np.pi**2):+.1f}%)")
print(f"karşıt: V_sıfır={V_sifir:.4f}  (koşullama farkı {V_sabit-V_sifir:.3f})")
print(f"kayıt : C1_asal={C1_asal:+.4f}  ölçülen lag-1 cov(ds)={C1_olc:+.4f}")
