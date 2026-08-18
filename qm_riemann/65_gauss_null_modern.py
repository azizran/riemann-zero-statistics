"""
65 — MODERN GAUSSIAN NULL: ANA PENCEREDE SURROGATE TESTİ (S3) (18 Ağustos)
============================================================================

Eski test (29): t ∈ [50, 500], ~260 aralık, 30 deneme → r_null = 0.505,
+6.7σ. Denetim S3: künye yanlıştı ve deney zayıftı.

Modern versiyon: L=9.86 penceresi (t ∈ [107k, 133k], ~40k aralık).
  1. Z(t) ince uniform gridde (adım 0.02, 1.27M nokta; RS 138 terim — hızlı)
  2. Faz-randomize surrogate'ler (|FFT| korunur, fazlar rastgele → aynı
     güç spektrumlu Gaussian süreç)
  3. Gerçek ve surrogate İÇİN AYNI grid-tabanlı çıkarım: işaret değişimi
     sıfırları, aralık başına grid-max, unfold, Pearson r
  4. 40 deneme → null dağılımı → gerçek ζ kaç σ dışarıda?
"""

import numpy as np
from pathlib import Path
from scipy.stats import pearsonr
import time

rng = np.random.default_rng(65)
HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi

def theta(t):
    return t / 2 * np.log(t / TWO_PI) - t / 2 - np.pi / 8 + 1 / (48 * t)

def Z_rs(t, chunk=200000):
    t = np.asarray(t, dtype=np.float64)
    out = np.empty_like(t)
    for s in range(0, len(t), chunk):
        tt = t[s:s + chunk]
        a = np.sqrt(tt / TWO_PI)
        N = a.astype(np.int64)
        th = theta(tt)
        z = np.zeros_like(tt)
        for Nv in np.unique(N):
            m = N == Nv
            n = np.arange(1, Nv + 1)
            ph = th[m, None] - tt[m, None] * np.log(n)[None, :]
            z[m] = 2 * (np.cos(ph) @ (n**-0.5))
        p = a - N
        cp = np.cos(TWO_PI * p)
        cp = np.where(np.abs(cp) < 1e-8, 1e-8, cp)
        psi = np.cos(TWO_PI * (p**2 - p - 1 / 16)) / cp
        z += (-1) ** (N - 1) * (tt / TWO_PI) ** -0.25 * psi
        out[s:s + chunk] = z
    return out

def gap_max_r(series, grid_t):
    """Grid serisinden: sıfırlar (işaret değişimi + lineer interp),
    aralık başına max|·|, unfold, Pearson r."""
    sgn = np.sign(series)
    sc = np.where(sgn[:-1] * sgn[1:] < 0)[0]
    if len(sc) < 100:
        return np.nan, 0
    # lineer interpolasyonla sıfır konumu
    t0 = grid_t[sc] - series[sc] * (grid_t[sc + 1] - grid_t[sc]) / (series[sc + 1] - series[sc])
    gaps = np.diff(t0)
    # aralık başına grid-max: reduceat (kesim indeksleri sc+1)
    absv = np.abs(series)
    mx = np.maximum.reduceat(absv, sc + 1)[:-1]  # son parça açık uçlu → at
    ok = gaps > 0
    gaps, mx = gaps[ok], mx[ok]
    Lw = np.log(0.5 * (t0[:-1] + t0[1:])[ok] / TWO_PI)
    g_u = gaps * Lw / TWO_PI
    a_u = mx / np.sqrt(np.maximum(Lw, 1.0))  # kaba ölçek; pencere içinde ~sabit
    r, _ = pearsonr(g_u, a_u)
    return r, len(gaps)

# 1) gerçek Z gridde
T_LO, T_HI, STEP = 107252.0, 132748.0, 0.02
print("Z(t) grid değerlendirmesi...")
t1 = time.time()
grid_t = np.arange(T_LO, T_HI, STEP)
Z = Z_rs(grid_t)
print(f"  {len(grid_t)} nokta, {time.time()-t1:.0f} s")

r_real, n_real = gap_max_r(Z, grid_t)
print(f"GERÇEK ζ (grid-tabanlı): r = {r_real:.4f} ({n_real} aralık)")

# 2) surrogate döngüsü
print("Surrogate'ler...")
F = np.fft.rfft(Z - Z.mean())
mag = np.abs(F)
rs = []
t1 = time.time()
N_TRIALS = 40
for k in range(N_TRIALS):
    ph = rng.uniform(0, TWO_PI, len(F))
    ph[0] = 0.0
    Fs = mag * np.exp(1j * ph)
    Zs = np.fft.irfft(Fs, n=len(Z))
    r_s, n_s = gap_max_r(Zs, grid_t)
    rs.append(r_s)
rs = np.array(rs)
print(f"  {N_TRIALS} deneme, {time.time()-t1:.0f} s")
print(f"\nNULL: r = {rs.mean():.4f} ± {rs.std(ddof=1):.4f}  "
      f"(aralık [{rs.min():.4f}, {rs.max():.4f}])")
exc = (r_real - rs.mean()) / rs.std(ddof=1)
print(f"GERÇEK ζ FAZLASI: ({r_real:.4f} − {rs.mean():.4f}) / {rs.std(ddof=1):.4f} "
      f"= +{exc:.1f}σ")
print(f"\n(Eski test: t∈[50,500], 260 aralık, +6.7σ idi. "
      f"Bu: L=9.86, {n_real} aralık.)")
