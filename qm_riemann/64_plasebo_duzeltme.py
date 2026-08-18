"""
64 — PLASEBO TASARIM DÜZELTMESİ + TABAN YENİDEN ÖLÇÜMÜ (K2) (18 Ağustos)
==========================================================================

Denetim K2: 52/54'te plasebo çekim aralığı (log 320'ye dek) kaçınma
listesini (p^k ≤ Q) aşıyordu → (Q, ~360] bandındaki GERÇEK asal çizgileri
plaseboya karışabiliyordu.

Düzeltme: kaçınma listesi = çekim aralığındaki TÜM asal kuvvetleri
(p^k ≤ 360), min mesafe 0.015. Yeniden ölçüm: 6 pencere (41) × Q ∈
{13, 50, 150, 300}, her biri 5 tekrar → temiz plasebo tabanları
(r-kayması ve soyulan varyans).
"""

import numpy as np
from pathlib import Path
from scipy.stats import pearsonr
from sympy import primerange

rng = np.random.default_rng(64)
HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
A = (np.e**2 - 5) / 2
B0, B1 = 2.7580, -0.0543

def prime_powers_upto(Q):
    out = []
    for p in primerange(2, Q + 1):
        pk = p
        while pk <= Q:
            out.append(float(pk)); pk *= p
    return sorted(out)

AVOID = np.log(np.array(prime_powers_upto(360)))  # TÜM asal kuvvetleri ≤360

def placebo_freqs(n_freq, lo, hi):
    out = []
    while len(out) < n_freq:
        c = rng.uniform(lo, hi)
        if np.abs(AVOID - c).min() > 0.015 and all(abs(o - c) > 0.015 for o in out):
            out.append(c)
    return np.array(out)

def strip(y, tvals, freqs):
    cols = [np.ones_like(y)]
    for om in freqs:
        cols += [np.cos(tvals * om), np.sin(tvals * om)]
    C = np.vstack(cols).T
    b, *_ = np.linalg.lstsq(C, y, rcond=None)
    part = C[:, 1:] @ b[1:]
    return y - part, float(part.var())

d41 = np.load(HERE / "41_bigT_windows.npz")
QS = [13, 50, 150, 300]
print(f"{'L':>6} | " + " ".join(f"Q={q}: Δr(maks)/Var".rjust(20) for q in QS))
summary = {q: {"dr": [], "va": []} for q in QS}
for wkey in sorted({x.split('_')[1] for x in d41.files}, key=lambda s: int(s[:-1])):
    gaps = d41[f"gaps_{wkey}"]; tmid = d41[f"tmid_{wkey}"]; amps = d41[f"amps_{wkey}"]
    Lw = np.log(tmid / TWO_PI); L = float(Lw.mean())
    g_u = gaps * Lw / TWO_PI
    a_u = amps / np.sqrt(A * Lw + B0 + B1 / Lw)
    a_u /= np.sqrt((a_u**2).mean())
    y_a, y_g = np.log(a_u), np.log(g_u)
    r0, _ = pearsonr(g_u, a_u)
    row = []
    for Q in QS:
        nf = len(prime_powers_upto(Q))
        drs, vas = [], []
        for rep in range(5):
            pf = placebo_freqs(nf, np.log(2) * 0.9, np.log(320))
            ya, va = strip(y_a, tmid, pf)
            yg, _ = strip(y_g, tmid, pf)
            ar = np.exp(ya); ar /= np.sqrt((ar**2).mean())
            gr = np.exp(yg); gr *= g_u.mean() / gr.mean()
            r1, _ = pearsonr(gr, ar)
            drs.append(r1 - r0); vas.append(va)
        dmax = max(abs(x) for x in drs)
        vmax = max(vas)
        summary[Q]["dr"].append(dmax); summary[Q]["va"].append(vmax)
        row.append(f"{dmax:.4f}/{vmax:.4f}")
    print(f"{L:>6.2f} | " + " ".join(x.rjust(20) for x in row))

print("\n=== DÜZELTİLMİŞ PLASEBO TABANLARI (6 pencere üstünden en kötü) ===")
for Q in QS:
    nf = len(prime_powers_upto(Q))
    print(f"  Q={Q:>3} ({nf:>2} frekans): |Δr| ≤ {max(summary[Q]['dr']):.4f}, "
          f"soyulan varyans ≤ {max(summary[Q]['va']):.4f}")
print("\n(Not 2'deki taban cümleleri bu sayılarla güncellenmeli.)")
