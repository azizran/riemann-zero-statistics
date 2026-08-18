"""
71 — V_res İNŞAATI: KALINTININ KİMLİĞİ (18 Ağustos 2026)
==========================================================

V_res(L) = 0.163→0.244 (Q=300 soyması sonrası). Üç aday bileşen:
  (A) soyulmamış kuyruk asalları (p>300) — ölçülmüş w(τ)'dan HESAPLANIR
  (B) boşluk-modülasyonlu asal kanalı — 43'ün w(g̃) bulgusu; sabit-katsayılı
      soyma bunu ıskalar → g̃-etkileşimli soymayla ÖLÇÜLÜR
  (C) gerçek öz-gürültü — artakalan

A: V_kuyruk(L) = ½ Σ_{p>300} w(τ_p)²/p ; w(τ): pozitif dal (yasa) +
   negatif plato (−0.167, τ>0.40; τ_end duyarlılığıyla)
B: taban = [p^k ≤ 300 cos/sin] + [aynıları × (g̃−1)] (316 sütun);
   eş-yapılı plasebo (tam kaçınma listesi) ile taban düşülür.
"""

import numpy as np
from pathlib import Path
from scipy.stats import pearsonr
from sympy import primerange

rng = np.random.default_rng(71)
HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
A = (np.e**2 - 5) / 2
B0, B1 = 2.7580, -0.0543
B_LAW = 1.0683
W_PLATO = -0.167
TAU0 = 0.40

def w_curve(tau):
    """Ölçülmüş tam w(τ): pozitif dal (yasa+v-eğrisi yaklaşıklığı) + plato."""
    # pozitif dal: v ≈ 2.014τ küçük-τ; v-doyumu ~0.575; √w = 1 − 1.0683 v
    v = np.minimum(2.014 * tau, 0.575)
    pos = np.clip(1 - B_LAW * v, 0, None) ** 2
    return np.where(tau < TAU0, pos, W_PLATO**2 * 0 + W_PLATO * np.ones_like(tau) * 0 + W_PLATO)  # düz plato

def prime_powers_upto(Q):
    out = []
    for p in primerange(2, Q + 1):
        pk = p
        while pk <= Q:
            out.append(float(pk)); pk *= p
    return sorted(out)

AVOID = np.log(np.array(prime_powers_upto(360)))

# ---------------- A: kuyruk muhasebesi ----------------
print("A — KUYRUK ASALLARI (hesap, ölçülmüş w(τ) ile):")
for L, vres_meas in [(9.86, 0.163), (12.45, 0.244)]:
    for tau_end, tag in [(0.7, "τ_end=0.7"), (1.0, "τ_end=1.0"), (1.5, "τ_end=1.5")]:
        pmax = int(np.exp(tau_end * L))
        ps = np.array(list(primerange(301, pmax)))
        taus = np.log(ps) / L
        wt = w_curve(taus)
        V_tail = 0.5 * np.sum(wt**2 / ps)
        if tag == "τ_end=1.0":
            print(f"  L={L}: V_kuyruk = {V_tail:.4f}  (ölçülen V_res {vres_meas})"
                  f"  → pay %{100*V_tail/vres_meas:.1f}   [{tag}, p≤{pmax}]")

# ---------------- B: g̃-etkileşimli soyma ----------------
def strip_vres(y_a, y_g, g_u, tmid, freqs, interact):
    cols = [np.ones_like(y_a)]
    gc = g_u - 1.0
    for om in freqs:
        c, s = np.cos(tmid * om), np.sin(tmid * om)
        cols += [c, s]
        if interact:
            cols += [c * gc, s * gc]
    X = np.vstack(cols).T
    ba, *_ = np.linalg.lstsq(X, y_a, rcond=None)
    bg, *_ = np.linalg.lstsq(X, y_g, rcond=None)
    ya = y_a - X[:, 1:] @ ba[1:]
    yg = y_g - X[:, 1:] @ bg[1:]
    ar = np.exp(ya); ar /= np.sqrt((ar**2).mean())
    gr = np.exp(yg); gr *= g_u.mean() / gr.mean()
    return float(ya.var()), pearsonr(gr, ar)[0]

def placebo_freqs(n):
    out = []
    while len(out) < n:
        c = rng.uniform(np.log(2) * 0.9, np.log(320))
        if np.abs(AVOID - c).min() > 0.015 and all(abs(o - c) > 0.015 for o in out):
            out.append(c)
    return out

FREQS = [np.log(q) for q in prime_powers_upto(300)]
d41 = np.load(HERE / "41_bigT_windows.npz")
print(f"\nB — G̃-ETKİLEŞİMLİ SOYMA (taban {len(FREQS)} frekans; "
      f"etkileşimli {4*len(FREQS)+1} sütun):")
print(f"{'L':>6} {'V_res(std)':>10} {'V_res(mod)':>10} {'yenen':>7} "
      f"{'plasebo-yeme':>12} {'NET-mod':>8} {'r*(std)':>8} {'r*(mod)':>8}")
for k in sorted({x.split("_")[1] for x in d41.files}, key=lambda s: int(s[:-1])):
    gaps, tmid, amps = d41[f"gaps_{k}"], d41[f"tmid_{k}"], d41[f"amps_{k}"]
    Lw = np.log(tmid / TWO_PI); L = float(Lw.mean())
    g_u = gaps * Lw / TWO_PI
    a_u = amps / np.sqrt(A * Lw + B0 + B1 / Lw)
    a_u /= np.sqrt((a_u**2).mean())
    y_a, y_g = np.log(a_u), np.log(g_u)
    v_std, r_std = strip_vres(y_a, y_g, g_u, tmid, FREQS, False)
    v_mod, r_mod = strip_vres(y_a, y_g, g_u, tmid, FREQS, True)
    pf = placebo_freqs(len(FREQS))
    v_p0, _ = strip_vres(y_a, y_g, g_u, tmid, pf, False)
    v_p1, _ = strip_vres(y_a, y_g, g_u, tmid, pf, True)
    eaten = v_std - v_mod
    plc_eat = v_p0 - v_p1
    print(f"{L:>6.2f} {v_std:>10.4f} {v_mod:>10.4f} {eaten:>7.4f} "
          f"{plc_eat:>12.4f} {eaten-plc_eat:>8.4f} {r_std:>8.4f} {r_mod:>8.4f}")
