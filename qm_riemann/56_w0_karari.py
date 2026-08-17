"""
56 — w(τ→0) KARARI: TAM SAYDAMLIK MI, PLATO MU? (17 Ağustos 2026)
===================================================================

Tüm pencereler (41: L=9.9-12.5, 55: L=16.6-23.5, 53: L=24.5) tek
konvansiyonla — w kanalı, tmid fazları, g̃+g̃² kontrollü — yeniden ölçülür.
p=2 ve p=3 noktaları (τ ≤ 0.12) kuadratik fit ile τ=0'a taşınır:

    w(τ) = w0 + c1·τ + c2·τ²

Karar: w0 ≈ 1 (perde sonsuz yükseklikte tam saydam — explicit formulanın
naif beklentisi) vs w0 ≈ 0.94 (kalıcı plato — yeni bir sabit demek olur).
Model kıyası: serbest w0 vs w0=1 sabitli fit, Δχ².
"""

import numpy as np
import mpmath as mp
import matplotlib.pyplot as plt
from pathlib import Path

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
A = (np.e**2 - 5) / 2
B0, B1 = 2.7580, -0.0543
T0_ODL = 267653395647
PRIMES_FIT = [2, 3]          # küçük-τ fit noktaları
PRIMES_ALL = [2, 3, 5, 7]    # ölçülen (şekil için)

mp.mp.dps = 30

def measure_w(gaps, amps, tmid, phase_const=None):
    """tmid-fazlı, kontrollü w_p ölçümü. phase_const: {p: sabit faz} (çapa)."""
    Lw = np.log(tmid / TWO_PI) if phase_const is None else \
         np.log((T0_ODL + tmid) / TWO_PI)
    g_u = gaps * Lw / TWO_PI
    a_u = amps / np.sqrt(A * Lw + B0 + B1 / Lw)
    a_u /= np.sqrt((a_u**2).mean())
    y = np.log(a_u)
    cols = [np.ones_like(y), g_u, g_u**2]
    for p in PRIMES_ALL:
        ph0 = phase_const[p] if phase_const else 0.0
        arg = ph0 + tmid * np.log(p)
        cols += [np.cos(arg), np.sin(arg)]
    X = np.vstack(cols).T
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ b
    se = np.sqrt(res.var() * np.diag(np.linalg.inv(X.T @ X)))
    L = float(Lw.mean())
    out = []
    for i, p in enumerate(PRIMES_ALL):
        w = b[3 + 2 * i] / p**-0.5
        s = se[3 + 2 * i] / p**-0.5
        out.append((p, np.log(p) / L, w, s))
    return L, out

points = []  # (p, tau, w, se, L)

# 41 pencereleri
d41 = np.load(HERE / "41_bigT_windows.npz")
for k in sorted({x.split("_")[1] for x in d41.files}, key=lambda s: int(s[:-1])):
    L, rows = measure_w(d41[f"gaps_{k}"], d41[f"amps_{k}"], d41[f"tmid_{k}"])
    points += [(p, tau, w, s, L) for p, tau, w, s in rows]

# 55 pencereleri
d55 = np.load(HERE / "55_kucuk_tau.npz")
for k in sorted({x.split("_")[1] for x in d55.files}, key=lambda s: float(s)):
    L, rows = measure_w(d55[f"gaps_{k}"], d55[f"amps_{k}"], d55[f"tmid_{k}"])
    points += [(p, tau, w, s, L) for p, tau, w, s in rows]

# 53 Odlyzko (çapalı sabit fazlar)
d53 = np.load(HERE / "53_odlyzko_amps.npz")
pc = {p: float(mp.fmod(T0_ODL * mp.log(p), 2 * mp.pi)) for p in PRIMES_ALL}
L, rows = measure_w(d53["gaps"], d53["max_amps"], d53["t_mid"], phase_const=pc)
points += [(p, tau, w, s, L) for p, tau, w, s in rows]

pts = np.array([(p, t, w, s, L) for p, t, w, s, L in points])
print(f"{'p':>3} {'L':>7} {'τ':>7} {'w':>8} {'±':>7}")
for p, t, w, s, L in sorted(points, key=lambda r: r[1]):
    print(f"{int(p):>3} {L:>7.2f} {t:>7.4f} {w:>8.4f} {s:>7.4f}")

# ---------------------------------------------------------------
# FIT: p=2,3 noktaları, τ ≤ 0.12
# ---------------------------------------------------------------
m = (np.isin(pts[:, 0], PRIMES_FIT)) & (pts[:, 1] <= 0.12)
tau, w, se = pts[m, 1], pts[m, 2], pts[m, 3]
W = 1 / se**2

def chi2(params, fix_w0=None):
    if fix_w0 is None:
        w0, c1, c2 = params
    else:
        w0 = fix_w0; c1, c2 = params
    return np.sum(W * (w - (w0 + c1 * tau + c2 * tau**2))**2)

# serbest fit (ağırlıklı polinom)
X = np.vstack([np.ones_like(tau), tau, tau**2]).T
sw = np.sqrt(W)
beta, *_ = np.linalg.lstsq(X * sw[:, None], w * sw, rcond=None)
cov = np.linalg.inv((X * W[:, None]).T @ X)
w0_free, se_w0 = beta[0], np.sqrt(cov[0, 0])
chi_free = chi2(beta)

# w0 = 1 sabitli
X1 = np.vstack([tau, tau**2]).T
b1, *_ = np.linalg.lstsq(X1 * sw[:, None], (w - 1.0) * sw, rcond=None)
chi_1 = chi2(b1, fix_w0=1.0)

# w0 = 0.94 sabitli
b94, *_ = np.linalg.lstsq(X1 * sw[:, None], (w - 0.94) * sw, rcond=None)
chi_94 = chi2(b94, fix_w0=0.94)

dof = m.sum() - 3
print(f"\nFIT ({m.sum()} nokta, p∈{{2,3}}, τ≤0.12):")
print(f"  SERBEST : w0 = {w0_free:.4f} ± {se_w0:.4f}   χ²/dof = {chi_free/dof:.2f}")
print(f"  w0 = 1.00 sabit: χ² = {chi_1:.1f}  (serbest: {chi_free:.1f}, Δχ² = {chi_1-chi_free:.1f})")
print(f"  w0 = 0.94 sabit: χ² = {chi_94:.1f}  (Δχ² = {chi_94-chi_free:.1f})")
print(f"\n  1'den uzaklık: {(1-w0_free)/se_w0:.1f}σ   0.94'ten: {abs(0.94-w0_free)/se_w0:.1f}σ")

# ---------------------------------------------------------------
# GRAFİK
# ---------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8.5, 5.6))
for p, c in zip(PRIMES_ALL, ["firebrick", "steelblue", "darkorange", "teal"]):
    mm = pts[:, 0] == p
    ax.errorbar(pts[mm, 1], pts[mm, 2], yerr=pts[mm, 3], fmt="o", ms=5, c=c,
                label=f"p={int(p)}")
xx = np.linspace(0, 0.14, 100)
ax.plot(xx, beta[0] + beta[1] * xx + beta[2] * xx**2, "k-", lw=1.4,
        label=f"serbest fit: w0={w0_free:.3f}±{se_w0:.3f}")
ax.plot(xx, 1.0 + b1[0] * xx + b1[1] * xx**2, "g--", lw=1.2, label="w0=1 sabitli")
ax.axhline(1, color="gray", lw=0.6)
ax.set_xlabel("τ = log p / L"); ax.set_ylabel("w — genlik kanalı iletimi")
ax.set_title("Perde τ→0'da tam saydamlaşıyor mu?")
ax.legend(fontsize=9); ax.grid(alpha=0.3)
plt.tight_layout()
out = HERE / "56_w0_karari.png"
plt.savefig(out, dpi=110)
print(f"\nGrafik: {out.name}")
