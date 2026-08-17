"""
57 — 10¹¹ PENCERESİ, ÇAPALI MOTORLA YENİDEN (17 Ağustos 2026, gece)
=====================================================================

56'nın kalite bayrağı: f64 taramasında 10¹¹ penceresi 158 "kurtarılan
yakın çift" üretti (diğer pencereler 2'şer) ve p=2 noktası Odlyzko-tabanlı
komşusundan 2.8σ yukarıda. Şüphe: 1.5e-3'lük faz hata tabanı sahte sıfır
bölünmeleri yarattı.

Çözüm: 53'ün ÇAPA AÇILIMI aynen — T0 = 10¹¹ (tam sayı), taban fazlar
mpmath'te bir kez, artanlar float64. Beklenen hata ~1e-8 (C0 kalıntısı).

Öngörü: temiz fazlarla 158 kurtarmanın büyük kısmı BUHARLAŞMALI
(sahtelerse) ve p=2 w noktası Odlyzko'ya yaklaşmalı.

Çıktı: 55_win_1e+11.npz YENİLENİR (eski f64 sürümü .bak olarak saklanır).
"""

import numpy as np
import mpmath as mp
from pathlib import Path
import time

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
T0 = 100_000_000_000  # tam sayı çapa

# ---------------------------------------------------------------
# Çapa kurulumu
# ---------------------------------------------------------------
mp.mp.dps = 40
t0m = mp.mpf(T0)
x0 = t0m / (2 * mp.pi)
a0m = mp.sqrt(x0)
N_RS = int(a0m)
frac_a0 = float(a0m - N_RS)
L0 = float(mp.log(x0))
print(f"T0=1e11: N_RS={N_RS}, frac(a0)={frac_a0:.6f}, L={L0:.4f}", flush=True)
assert 0.005 < frac_a0 < 0.995, "frac(a0) sınıra yakın — T0 kaydır!"

theta0 = t0m / 2 * mp.log(x0) - t0m / 2 - mp.pi / 8 + 1 / (48 * t0m) + 7 / (5760 * t0m**3)
theta0_mod = float(mp.fmod(theta0, 2 * mp.pi))
thp = float(mp.log(x0) / 2)
a_slope = float(1 / (4 * mp.pi * a0m))

t1 = time.time()
base = np.empty(N_RS)
mp.mp.dps = 30
for n in range(1, N_RS + 1):
    base[n - 1] = float(mp.fmod(t0m * mp.log(n), 2 * mp.pi))
base = (theta0_mod - base) % TWO_PI
print(f"taban fazlar: {time.time()-t1:.0f} s", flush=True)
lnn = np.log(np.arange(1, N_RS + 1))
wts = np.arange(1, N_RS + 1) ** -0.5
coef = thp - lnn
sign_c0 = (-1.0) ** (N_RS - 1)
c0_amp = float(x0 ** -0.25)

def Z_anchor(dt, chunk=400):
    dt = np.asarray(dt, dtype=np.float64)
    out = np.empty_like(dt)
    for s in range(0, len(dt), chunk):
        d = dt[s:s + chunk]
        q = d * d / (4.0 * T0)
        ph = base[None, :] + d[:, None] * coef[None, :] + q[:, None]
        z = 2.0 * (np.cos(ph) @ wts)
        p = (frac_a0 + a_slope * d) % 1.0
        cp = np.cos(TWO_PI * p)
        cp = np.where(np.abs(cp) < 1e-8, 1e-8, cp)
        psi = np.cos(TWO_PI * (p * p - p - 1.0 / 16.0)) / cp
        z += sign_c0 * c0_amp * psi
        out[s:s + chunk] = z
    return out

# doğrulama (mpmath)
mp.mp.dps = 15
rngv = np.random.default_rng(57)
dts = rngv.uniform(-50, 50, 5)
za = Z_anchor(dts)
zm = np.array([float(mp.siegelz(T0 + float(d))) for d in dts])
err = np.abs(za - zm).max()
print(f"doğrulama (5 nokta): max hata = {err:.2e}", flush=True)

# ---------------------------------------------------------------
# Tarama (dt uzayında)
# ---------------------------------------------------------------
NG = 20000
L = L0
half = NG * (TWO_PI / L) / 2
step = (TWO_PI / L) / 18
print(f"tarama: dt ∈ [{-half:.0f}, {half:.0f}], adım {step:.4f}", flush=True)
t1 = time.time()
grid = np.arange(-half, half, step)
Zg = Z_anchor(grid)
sc = np.where(np.sign(Zg[:-1]) != np.sign(Zg[1:]))[0]
lo, hi = grid[sc].copy(), grid[sc + 1].copy()
flo = Zg[sc].copy()
for _ in range(18):
    mid = 0.5 * (lo + hi)
    fm = Z_anchor(mid)
    left = np.sign(fm) == np.sign(flo)
    lo = np.where(left, mid, lo)
    flo = np.where(left, fm, flo)
    hi = np.where(left, hi, mid)
zeros = 0.5 * (lo + hi)
print(f"kaba sıfır: {len(zeros)} ({(time.time()-t1)/60:.1f} dk)", flush=True)

extra = []
for i in range(len(zeros) - 1):
    a_, b_ = zeros[i], zeros[i + 1]
    mask = (grid > a_) & (grid < b_)
    if mask.sum() and np.abs(Zg[mask]).min() < 0.1:
        tf = np.linspace(a_ + step / 50, b_ - step / 50, 200)
        Zf = Z_anchor(tf)
        sc2 = np.where(np.sign(Zf[:-1]) != np.sign(Zf[1:]))[0]
        for jj in sc2:
            l2, h2, fl2 = tf[jj], tf[jj + 1], Zf[jj]
            for _ in range(30):
                m2 = 0.5 * (l2 + h2)
                f2 = float(Z_anchor(np.array([m2]))[0])
                if np.sign(f2) == np.sign(fl2):
                    l2, fl2 = m2, f2
                else:
                    h2 = m2
            extra.append(0.5 * (l2 + h2))
if extra:
    zeros = np.sort(np.concatenate([zeros, np.array(extra)]))
print(f"kurtarılan yakın-çift üyesi: {len(extra)}  "
      f"(f64 taramada 158 idi — buharlaşma testi!)", flush=True)

g_lo, g_hi = zeros[:-1], zeros[1:]
gaps = g_hi - g_lo
u = np.arange(1, 25) / 25
tt = g_lo[:, None] + gaps[:, None] * u[None, :]
Am = np.abs(Z_anchor(tt.ravel())).reshape(tt.shape)
j = np.argmax(Am, axis=1)
rows = np.arange(len(j))
jc = np.clip(j, 1, 22)
y0, y1, y2 = Am[rows, jc - 1], Am[rows, jc], Am[rows, jc + 1]
den = y0 - 2 * y1 + y2
dpos = np.where(np.abs(den) > 1e-12, 0.5 * (y0 - y2) / den, 0.0)
t_pk = tt[rows, jc] + dpos * gaps / 25
A_pk = np.abs(Z_anchor(t_pk))
max_amps = np.maximum(Am[rows, j], A_pk)
print(f"genlikler tamam ({(time.time()-t1)/60:.1f} dk toplam)", flush=True)

# ---------------------------------------------------------------
# Eski f64 penceresiyle kıyas + kayıt
# ---------------------------------------------------------------
old = np.load(HERE / "55_win_1e+11.npz")
print(f"\nKIYAS: eski {len(old['gaps'])} aralık, yeni {len(gaps)} aralık")
tmid_abs = T0 + 0.5 * (g_lo + g_hi)
import shutil
shutil.move(HERE / "55_win_1e+11.npz", HERE / "55_win_1e11_eski_f64.bak")
np.savez(HERE / "55_win_1e+11.npz", gaps=gaps, amps=max_amps, tmid=tmid_abs)
print("Kaydedildi: 55_win_1e+11.npz (eski → 55_win_1e11_eski_f64.bak)", flush=True)
