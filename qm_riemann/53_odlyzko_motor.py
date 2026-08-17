"""
53 — ODLYZKO YÜKSEKLİĞİ: ÇAPALI RIEMANN-SIEGEL MOTORU (17 Ağustos 2026)
=========================================================================

Hedef: t ≈ 2.6765×10¹¹ (10¹²'inci sıfır civarı, L=24.48) — zeros3 tablosunun
9,999 aralığında max|Z| hesabı.

Sorun: çıplak float64'te faz t·log n ~ 3×10¹² → hata ~10⁻³ rad (sınırda),
θ(t) benzer. Çözüm: ÇAPA AÇILIMI —
  faz(n, t0+dt) = [θ(t0) − t0·ln n mod 2π]  (mpmath, dps=30, bir kez)
               + (θ'(t0) − ln n)·dt          (float64, |·| ≤ 3×10⁴ → tam)
               + dt²/(4t0)                   (ortak skaler, ≤ 6×10⁻⁶)

t0 = 267653395647 (tablo tabanı, TAM SAYI — mod 2π hesabı temiz).

Doğrulama: tablodaki sıfır konumlarında |Z| ≈ 0 çıkmalı (tablo = hakikat).
Çıktı: 53_odlyzko_amps.npz (gaps, max_amps, t_mid ofsetleri).
"""

import numpy as np
import mpmath as mp
from pathlib import Path
import time

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
T0 = 267653395647  # tablo tabanı (tam sayı)

# ---------------------------------------------------------------
# 1) Tabloyu oku
# ---------------------------------------------------------------
offs = []
with open(HERE / "odlyzko_zeros3.txt") as f:
    for line in f:
        s = line.strip()
        try:
            offs.append(float(s))
        except ValueError:
            continue
offs = np.array(offs)
print(f"{len(offs)} sıfır ofseti okundu: [{offs[0]:.4f}, {offs[-1]:.4f}]")
L_here = float(np.log((T0 + offs.mean()) / TWO_PI))
print(f"L = {L_here:.4f},  ort. aralık = {np.diff(offs).mean():.4f} (2π/L = {TWO_PI/L_here:.4f})")

# ---------------------------------------------------------------
# 2) Çapa hesapları (mpmath, bir kez)
# ---------------------------------------------------------------
mp.mp.dps = 40
t0m = mp.mpf(T0)
x0 = t0m / (2 * mp.pi)          # t0/2π
a0m = mp.sqrt(x0)               # RS kesim uzunluğu
N_RS = int(a0m)
print(f"RS terim sayısı N = {N_RS}")

theta0 = t0m / 2 * mp.log(x0) - t0m / 2 - mp.pi / 8 + 1 / (48 * t0m) + 7 / (5760 * t0m**3)
theta0_mod = float(mp.fmod(theta0, 2 * mp.pi))
thp = float(mp.log(x0) / 2)     # θ'(t0) = ½ log(t0/2π)
a0_frac_base = a0m - N_RS       # p tabanı
a_slope = float(1 / (4 * mp.pi * a0m))  # da/dt

print("Taban fazlar hesaplanıyor (206k × mpmath)...")
t1 = time.time()
base = np.empty(N_RS)
mp.mp.dps = 30
for n in range(1, N_RS + 1):
    # (θ(t0) − t0·ln n) mod 2π ;  θ kısmını float'ta ekleyeceğiz
    base[n - 1] = float(mp.fmod(t0m * mp.log(n), 2 * mp.pi))
print(f"  {time.time()-t1:.0f} s")
base = (theta0_mod - base) % TWO_PI     # taban faz ∈ [0, 2π)
lnn = np.log(np.arange(1, N_RS + 1))
wts = np.arange(1, N_RS + 1) ** -0.5
coef = thp - lnn                        # dt katsayısı (float64 tam)
sign_c0 = (-1.0) ** (N_RS - 1)
p_base = float(a0_frac_base)

def Z_anchor(dt, chunk=40):
    """Z(t0+dt); dt: 1B array (0..~2600)."""
    dt = np.asarray(dt, dtype=np.float64)
    out = np.empty_like(dt)
    for s in range(0, len(dt), chunk):
        d = dt[s:s + chunk]
        q = d * d / (4.0 * T0)          # ortak kuadratik terim
        ph = base[None, :] + d[:, None] * coef[None, :] + q[:, None]
        z = 2.0 * (np.cos(ph) @ wts)
        # C0 düzeltmesi
        p = (p_base + a_slope * d) % 1.0
        cp = np.cos(TWO_PI * p)
        cp = np.where(np.abs(cp) < 1e-8, 1e-8, cp)
        psi = np.cos(TWO_PI * (p * p - p - 1.0 / 16.0)) / cp
        z += sign_c0 * float((x0) ** -0.25) * psi
        out[s:s + chunk] = z
    return out

# ---------------------------------------------------------------
# 3) DOĞRULAMA: tablo sıfırlarında |Z| ≈ 0 mı?
# ---------------------------------------------------------------
print("\nDoğrulama: 300 tablo sıfırında |Z|...")
idx = np.linspace(0, len(offs) - 1, 300).astype(int)
Zz = np.abs(Z_anchor(offs[idx]))
print(f"  medyan |Z(sıfır)| = {np.median(Zz):.2e},  %95'lik = {np.quantile(Zz, 0.95):.2e}")
mid = 0.5 * (offs[:-1] + offs[1:])
Zm = np.abs(Z_anchor(mid[idx[:150]]))
print(f"  kıyas: aralık ortası medyan |Z| = {np.median(Zm):.3f}")
print(f"  → sıfır/tepe oranı ~ {np.median(Zz)/np.median(Zm):.1e}  (küçükse motor sağlam)")

# ---------------------------------------------------------------
# 4) TÜM ARALIKLARDA MAX |Z|
# ---------------------------------------------------------------
print("\n9,999 aralıkta max|Z| (24 nokta + parabolik)...")
t1 = time.time()
g_lo, g_hi = offs[:-1], offs[1:]
gaps = g_hi - g_lo
u = np.arange(1, 25) / 25
tt = g_lo[:, None] + gaps[:, None] * u[None, :]
Zg = np.abs(Z_anchor(tt.ravel(), chunk=60)).reshape(tt.shape)
j = np.argmax(Zg, axis=1)
rows = np.arange(len(j))
jc = np.clip(j, 1, 22)
y0, y1, y2 = Zg[rows, jc - 1], Zg[rows, jc], Zg[rows, jc + 1]
den = y0 - 2 * y1 + y2
dpos = np.where(np.abs(den) > 1e-12, 0.5 * (y0 - y2) / den, 0.0)
t_pk = tt[rows, jc] + dpos * gaps / 25
A_pk = np.abs(Z_anchor(t_pk, chunk=60))
max_amps = np.maximum(Zg[rows, j], A_pk)
print(f"  {(time.time()-t1)/60:.1f} dk")

np.savez(HERE / "53_odlyzko_amps.npz",
         gaps=gaps, max_amps=max_amps, t_mid=mid, L=L_here)
print(f"Kaydedildi: 53_odlyzko_amps.npz  (max|Z| ort = {max_amps.mean():.3f})")
