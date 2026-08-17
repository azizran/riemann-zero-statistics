"""
55 — KÜÇÜK-τ SEFERİ: ARA YÜKSEKLİK PENCERELERİ (17 Ağustos 2026)
==================================================================

w(τ→0) sorusu: perde tam saydamlaşıyor mu (w→1) yoksa ~0.94'te mi duruyor?
Ekstrapolasyonu beslemek için p=2'nin küçük-τ noktalarını yoğunlaştırıyoruz:

  t = 10⁸  (L=16.58, τ₂=0.0418)   30k aralık
  t = 10⁹  (L=18.88, τ₂=0.0367)   30k aralık
  t = 10¹⁰ (L=21.18, τ₂=0.0327)   25k aralık
  t = 10¹¹ (L=23.48, τ₂=0.0295)   20k aralık

Motor: 41'in düz float64 RS'i (bu aralıkta faz hatası ≤3×10⁻⁴ rad — güvenli;
mpmath örneklemiyle her yükseklikte doğrulanır). Çıktı: 55_kucuk_tau.npz
"""

import numpy as np
import mpmath as mp
from pathlib import Path
import time

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi

def theta(t):
    return t / 2 * np.log(t / TWO_PI) - t / 2 - np.pi / 8 + 1 / (48 * t) + 7 / (5760 * t**3)

def Z_rs(t, chunk=None):
    t = np.asarray(t, dtype=np.float64)
    if chunk is None:
        Nmax = int(np.sqrt(t.max() / TWO_PI)) + 1
        chunk = max(200, int(1.0e8 / Nmax))
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
        cp = np.where(np.abs(cp) < 1e-8, 1e-8 * np.sign(cp + 1e-300), cp)
        psi = np.cos(TWO_PI * (p**2 - p - 1 / 16)) / cp
        z += (-1) ** (N - 1) * (tt / TWO_PI) ** -0.25 * psi
        out[s:s + chunk] = z
    return out

def scan_window(t_lo, t_hi):
    L = np.log(0.5 * (t_lo + t_hi) / TWO_PI)
    step = (TWO_PI / L) / 18
    grid = np.arange(t_lo, t_hi, step)
    Zg = Z_rs(grid)
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
    # şüpheli aralık (kaçmış yakın çift) taraması
    extra = []
    for i in range(len(zeros) - 1):
        a_, b_ = zeros[i], zeros[i + 1]
        mask = (grid > a_) & (grid < b_)
        if mask.sum() and np.abs(Zg[mask]).min() < 0.1:
            tf = np.linspace(a_ + step / 50, b_ - step / 50, 200)
            Zf = Z_rs(tf)
            sc2 = np.where(np.sign(Zf[:-1]) != np.sign(Zf[1:]))[0]
            for jj in sc2:
                l2, h2, fl2 = tf[jj], tf[jj + 1], Zf[jj]
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
    g_lo, g_hi = zeros[:-1], zeros[1:]
    gaps = g_hi - g_lo
    u = np.arange(1, 25) / 25
    tt = g_lo[:, None] + gaps[:, None] * u[None, :]
    Am = np.abs(Z_rs(tt.ravel())).reshape(tt.shape)
    j = np.argmax(Am, axis=1)
    rows = np.arange(len(j))
    jc = np.clip(j, 1, 23)
    y0, y1, y2 = Am[rows, jc - 1], Am[rows, jc], Am[rows, jc + 1]
    den = y0 - 2 * y1 + y2
    dpos = np.where(np.abs(den) > 1e-12, 0.5 * (y0 - y2) / den, 0.0)
    t_pk = tt[rows, jc] + dpos * gaps / 25
    A_pk = np.abs(Z_rs(t_pk))
    return gaps, np.maximum(Am[rows, j], A_pk), 0.5 * (g_lo + g_hi), len(extra)

CENTERS = [(1e8, 30000), (1e9, 30000), (1e10, 25000), (1e11, 20000)]
NVAL = {1e8: 10, 1e9: 8, 1e10: 5, 1e11: 3}

out = {}
for c, ng in CENTERS:
    L = np.log(c / TWO_PI)
    # mpmath doğrulama
    mp.mp.dps = 15
    rngv = np.random.default_rng(int(c % 97) + 7)
    ts = c + rngv.uniform(0, 50, NVAL[c])
    zr = Z_rs(ts)
    t1 = time.time()
    zm = np.array([float(mp.siegelz(t)) for t in ts])
    err = np.abs(zr - zm).max()
    print(f"t={c:.0e}: doğrulama max hata = {err:.2e} ({time.time()-t1:.0f}s)", flush=True)
    if err > 5e-3:
        print("  UYARI: hata büyük — bu pencere şüpheli!", flush=True)
    half = ng * (TWO_PI / L) / 2
    t1 = time.time()
    gaps, amps, tmid, nex = scan_window(c - half, c + half)
    key = f"{c:.0e}"
    out[f"gaps_{key}"] = gaps
    out[f"amps_{key}"] = amps
    out[f"tmid_{key}"] = tmid
    print(f"  L={L:.3f}: {len(gaps)} aralık, {nex} kurtarılan, "
          f"{(time.time()-t1)/60:.1f} dk", flush=True)

np.savez(HERE / "55_kucuk_tau.npz", **out)
print("Kaydedildi: 55_kucuk_tau.npz", flush=True)
