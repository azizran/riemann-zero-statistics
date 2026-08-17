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

# 55b pencereleri (pencere-başına dosyalar)
for f in sorted(HERE.glob("55_win_*.npz")):
    d55 = np.load(f)
    L, rows = measure_w(d55["gaps"], d55["amps"], d55["tmid"])
    points += [(p, tau, w, s, L) for p, tau, w, s in rows]
    print(f"  yüklendi: {f.name} (L={L:.2f})")

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
# FIT VARYANTLARI (ön koşudaki ders: saf τ-çökmesi yetmiyor —
# aynı τ'da büyük p daha az geçiriyor. Sistematik pay için 4 model.)
# ---------------------------------------------------------------
def wfit(mask, cols_extra, label):
    tau = pts[mask, 1]; w = pts[mask, 2]; se = pts[mask, 3]
    pv = pts[mask, 0]; Lv = pts[mask, 4]
    W = 1 / se**2; sw = np.sqrt(W)
    cols = [np.ones_like(tau), tau, tau**2]
    for ce in cols_extra:
        cols.append({"invp": 1 / pv, "invL": 1 / Lv}[ce])
    X = np.vstack(cols).T
    beta, *_ = np.linalg.lstsq(X * sw[:, None], w * sw, rcond=None)
    cov = np.linalg.inv((X * W[:, None]).T @ X)
    resid = w - X @ beta
    chi = float(np.sum(W * resid**2)); dof = mask.sum() - X.shape[1]
    # w0: τ→0 limiti. 1/p ve 1/L terimleri τ→0'da p→∞/L→∞ limitinde düşer
    w0, s0 = beta[0], np.sqrt(cov[0, 0])
    print(f"  {label:<28} w0 = {w0:.4f} ± {s0:.4f}   χ²/dof = {chi/dof:.2f}"
          f"   (1'e uzaklık {abs(1-w0)/s0:.1f}σ)")
    return w0, s0, chi / dof, beta

m23 = (np.isin(pts[:, 0], [2, 3])) & (pts[:, 1] <= 0.12)
m2 = (pts[:, 0] == 2)
m235 = (np.isin(pts[:, 0], [2, 3, 5])) & (pts[:, 1] <= 0.14)

print(f"\nFIT VARYANTLARI:")
r1 = wfit(m23, [], "A: p∈{2,3} saf kuadratik")
r2 = wfit(m23, ["invp"], "B: p∈{2,3} + d/p terimi")
r3 = wfit(m23, ["invL"], "C: p∈{2,3} + e/L terimi")
r4 = wfit(m2, [], "D: yalnız p=2, kuadratik")
r5 = wfit(m235, ["invp"], "E: p∈{2,3,5} + d/p (geniş)")

w0s = [r[0] for r in (r1, r2, r3, r4, r5)]
print(f"\n  w0 aralığı (5 model): [{min(w0s):.3f}, {max(w0s):.3f}]")
print(f"  → sistematik pay istatistikten büyükse aralığı raporla, tek sayıyı değil")
beta = r1[3]; w0_free, se_w0 = r1[0], r1[1]

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
        label=f"model A: w0={w0_free:.3f}±{se_w0:.3f}")
ax.axhline(1, color="gray", lw=0.6)
ax.set_xlabel("τ = log p / L"); ax.set_ylabel("w — genlik kanalı iletimi")
ax.set_title("Perde τ→0'da tam saydamlaşıyor mu?")
ax.legend(fontsize=9); ax.grid(alpha=0.3)
plt.tight_layout()
out = HERE / "56_w0_karari.png"
plt.savefig(out, dpi=110)
print(f"\nGrafik: {out.name}")
