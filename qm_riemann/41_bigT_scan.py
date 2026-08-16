"""
41 — BÜYÜK T TARAMASI: VEKTÖRİZE RIEMANN-SIEGEL (16 Ağustos 2026)
===================================================================

Amaç (Plan A): N_eff − L ≈ +0.84 sabitliğini L=12'ye kadar test etmek.
36'nın mpmath yolu T=100k'da 2-4 saatti; burada RS ana toplamı numpy ile
vektörize → T=1.6M'a kadar 6 pencere, dakikalar mertebesinde.

Z(t) = 2 Σ_{n≤N} n^{-1/2} cos(θ(t) − t·log n) + R(t),  N = ⌊√(t/2π)⌋
R(t) ≈ (−1)^{N−1} (t/2π)^{-1/4} Ψ(p),  p = √(t/2π) − N
Ψ(p) = cos(2π(p² − p − 1/16)) / cos(2πp)      (C0 düzeltmesi)

Hata ~ t^{-3/4}: t ≥ 10⁵'te ~1e-4 → genlik (O(1-30)) ve sıfır konumu
(aralık ~0.5) için fazlasıyla yeterli. Doğrulama: (a) mpmath örneklem,
(b) 36'nın çifte-doğrulanmış verisiyle örtüşme bölgesi.

Kaçırılan yakın çiftler: grid adımı ort_aralık/20 + şüpheli aralıklarda
(iç min|Z| < 0.1) 10× ince yeniden tarama.
"""

import numpy as np
from pathlib import Path
import mpmath as mp
import time

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi

# ---------------------------------------------------------------
# Vektörize Riemann-Siegel
# ---------------------------------------------------------------
def theta(t):
    return t / 2 * np.log(t / TWO_PI) - t / 2 - np.pi / 8 + 1 / (48 * t) + 7 / (5760 * t**3)

def Z_rs(t, chunk=20000):
    """Vektörize Z(t) (C0 düzeltmeli). t: 1B array (pozitif, büyük)."""
    t = np.asarray(t, dtype=np.float64)
    out = np.empty_like(t)
    for s in range(0, len(t), chunk):
        tt = t[s:s + chunk]
        a = np.sqrt(tt / TWO_PI)
        N = a.astype(np.int64)
        th = theta(tt)
        z = np.zeros_like(tt)
        # N gruplarına böl (pencere içinde N az değişir)
        for Nv in np.unique(N):
            m = N == Nv
            n = np.arange(1, Nv + 1)
            ph = th[m, None] - tt[m, None] * np.log(n)[None, :]
            z[m] = 2 * (np.cos(ph) @ (n**-0.5))
        # C0 düzeltmesi
        p = a - N
        cp = np.cos(TWO_PI * p)
        cp = np.where(np.abs(cp) < 1e-8, 1e-8 * np.sign(cp + 1e-300), cp)  # tekillik koruması
        psi = np.cos(TWO_PI * (p**2 - p - 1 / 16)) / cp
        z += (-1) ** (N - 1) * (tt / TWO_PI) ** -0.25 * psi
        out[s:s + chunk] = z
    return out

# ---------------------------------------------------------------
# Pencere işleyici: sıfırlar + aralık başına max|Z|
# ---------------------------------------------------------------
def scan_window(t_lo, t_hi, verbose=True):
    L = np.log(0.5 * (t_lo + t_hi) / TWO_PI)
    step = (TWO_PI / L) / 20  # ort. aralığın 1/20'si
    grid = np.arange(t_lo, t_hi, step)
    Zg = Z_rs(grid)

    # işaret değişimleri → bisection (18 iter)
    sc = np.where(np.sign(Zg[:-1]) != np.sign(Zg[1:]))[0]
    lo, hi = grid[sc].copy(), grid[sc + 1].copy()
    flo = Zg[sc].copy()
    for _ in range(18):
        mid = 0.5 * (lo + hi)
        fm = Z_rs(mid)
        left = np.sign(fm) == np.sign(flo)
        lo = np.where(left, mid, lo)
        flo = np.where(left, fm, flo)
        hi = np.where(left, hi, mid)
    zeros = 0.5 * (lo + hi)

    # şüpheli aralıklar: iç |Z| min'i küçükse ince yeniden tara (kaçmış çift?)
    extra = []
    gaps0 = np.diff(zeros)
    for i in np.where(gaps0 > 0)[0]:
        a_, b_ = zeros[i], zeros[i + 1]
        mask = (grid > a_) & (grid < b_)
        if mask.sum() and np.abs(Zg[mask]).min() < 0.1:
            tf = np.linspace(a_ + step / 50, b_ - step / 50, 200)
            Zf = Z_rs(tf)
            sc2 = np.where(np.sign(Zf[:-1]) != np.sign(Zf[1:]))[0]
            for j in sc2:
                l2, h2, fl2 = tf[j], tf[j + 1], Zf[j]
                for _ in range(30):
                    m2 = 0.5 * (l2 + h2)
                    f2 = float(Z_rs(np.array([m2]))[0])
                    if np.sign(f2) == np.sign(fl2):
                        l2, fl2 = m2, f2
                    else:
                        h2 = m2
                extra.append(0.5 * (l2 + h2))
    if extra:
        zeros = np.sort(np.concatenate([zeros, np.array(extra)]))

    # aralık başına max|Z|: 24 iç nokta + parabolik incelik
    g_lo, g_hi = zeros[:-1], zeros[1:]
    u = np.arange(1, 25) / 25
    tt = g_lo[:, None] + (g_hi - g_lo)[:, None] * u[None, :]
    A = np.abs(Z_rs(tt.ravel()).reshape(tt.shape))
    j = np.argmax(A, axis=1)
    rows = np.arange(len(j))
    jc = np.clip(j, 1, 23)
    y0, y1, y2 = A[rows, jc - 1], A[rows, jc], A[rows, jc + 1]
    denom = y0 - 2 * y1 + y2
    dpos = np.where(np.abs(denom) > 1e-12, 0.5 * (y0 - y2) / denom, 0.0)
    t_pk = tt[rows, jc] + dpos * (g_hi - g_lo) / 25
    A_pk = np.abs(Z_rs(t_pk))
    max_amps = np.maximum(A[rows, j], A_pk)

    if verbose:
        print(f"  [{t_lo:.0f}, {t_hi:.0f}] L={L:.2f}: {len(zeros)} sıfır, "
              f"{len(extra)} kurtarılan yakın-çift üyesi")
    return zeros, np.diff(zeros), max_amps, 0.5 * (g_lo + g_hi)

# ---------------------------------------------------------------
# 0) DOĞRULAMA
# ---------------------------------------------------------------
print("=== DOĞRULAMA 1: mpmath örneklem (200 nokta) ===")
mp.mp.dps = 20
rng = np.random.default_rng(41)
ts = rng.uniform(1e5, 1.6e6, 200)
z_rs = Z_rs(ts)
z_mp = np.array([float(mp.siegelz(t)) for t in ts])
err = np.abs(z_rs - z_mp)
print(f"  max hata = {err.max():.2e}, medyan = {np.median(err):.2e}  "
      f"({'OK' if err.max() < 1e-3 else 'SORUN!'})")

print("\n=== DOĞRULAMA 2: 36 verisiyle örtüşme (t ∈ [70000, 74000]) ===")
d36 = np.load(HERE / "36_T100k.npz")
zeros_v, gaps_v, amps_v, tm_v = scan_window(70000, 74000)
m36 = (d36["t_mid"] > zeros_v[0]) & (d36["t_mid"] < zeros_v[-1])
g36, a36, t36 = d36["intervals"][m36], d36["max_amps"][m36], d36["t_mid"][m36]
if len(g36) == len(gaps_v):
    print(f"  aralık sayısı eşleşti: {len(g36)}")
    print(f"  gap farkı   : max {np.abs(g36 - gaps_v).max():.2e}")
    rel = amps_v / a36 - 1
    print(f"  amp oranı−1 : medyan {np.median(rel):+.4f}, "
          f"%95'lik |·| {np.quantile(np.abs(rel), 0.95):.4f}")
    print("  (36 kaba optimizer kullanmıştı; bizim ince grid hafif YÜKSEK olmalı)")
else:
    print(f"  UYARI: sayı farkı! 36: {len(g36)}, RS: {len(gaps_v)} → incele")

# ---------------------------------------------------------------
# 1) YENİ PENCERELER
# ---------------------------------------------------------------
print("\n=== YENİ PENCERELER (her biri ~40k aralık) ===")
centers = [1.2e5, 2.0e5, 3.5e5, 6.0e5, 1.0e6, 1.6e6]
N_GAPS = 40000
out = {}
t0 = time.time()
for c in centers:
    L = np.log(c / TWO_PI)
    half = N_GAPS * (TWO_PI / L) / 2
    zeros, gaps, amps, tmid = scan_window(c - half, c + half)
    key = f"{int(c/1000)}k"
    out[f"gaps_{key}"] = gaps
    out[f"amps_{key}"] = amps
    out[f"tmid_{key}"] = tmid
print(f"Toplam süre: {(time.time()-t0)/60:.1f} dk")

np.savez(HERE / "41_bigT_windows.npz", **out)
print("Kaydedildi: 41_bigT_windows.npz")
