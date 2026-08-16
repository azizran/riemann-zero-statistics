"""
42b — N_eff SAĞLAMLIK: PÜRÜZSÜZ CUE + SPEARMAN (16 Ağustos 2026)
==================================================================

42'deki zayıf trend (Δχ²=7.4) iki artefakttan olabilir:
(a) CUE r(N) tablosundaki Monte-Carlo zikzağı (±0.0013 → N_eff'e ±0.1)
    → çare: r_CUE(N) = r_∞ + c₁/N + c₂/N² pürüzsüz fiti
(b) Pearson'ın marjinal şekle duyarlılığı (ζ ve CUE genlik dağılımları
    farklı kuyruklara sahip olabilir)
    → çare: Spearman (rank) ile aynı analiz — marjinallerden bağımsız
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import pearsonr, spearmanr
from scipy.optimize import brentq

rng = np.random.default_rng(423)
HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
A = (np.e**2 - 5) / 2
B0, B1 = 2.7580, -0.0543

# ---------------------------------------------------------------
# 1) CUE: her N için HEM Pearson HEM Spearman (200k aralık/N)
# ---------------------------------------------------------------
def haar_unitary_batch(M, N):
    G = (rng.standard_normal((M, N, N)) + 1j * rng.standard_normal((M, N, N))) / np.sqrt(2)
    Q, R = np.linalg.qr(G)
    diag = np.einsum("mii->mi", R)
    return Q * (diag / np.abs(diag))[:, None, :]

def cue_both(N, n_gaps=200000, grid=48, chunk=1200):
    M = int(np.ceil(n_gaps / N))
    gs, ms = [], []
    for s in range(0, M, chunk):
        m = min(chunk, M - s)
        U = haar_unitary_batch(m, N)
        ph = np.sort(np.angle(np.linalg.eigvals(U)), axis=1)
        gaps = np.diff(np.concatenate([ph, ph[:, :1] + TWO_PI], axis=1), axis=1)
        u = np.arange(1, grid + 1) / (grid + 1)
        th = ph[:, :, None] + gaps[:, :, None] * u[None, None, :]
        df = th[:, :, :, None] - ph[:, None, None, :]
        amp = np.prod(2 * np.abs(np.sin(df / 2)), axis=-1)
        gs.append(gaps.ravel()); ms.append(amp.max(axis=-1).ravel())
    g = np.concatenate(gs); a = np.concatenate(ms)
    return pearsonr(g, a)[0], spearmanr(g, a)[0], len(g)

print("CUE (200k aralık/N):")
Ns = np.arange(5, 20)
cueP, cueS = [], []
for N in Ns:
    rp, rs, n = cue_both(N)
    cueP.append(rp); cueS.append(rs)
    print(f"  N={N:>2}: Pearson {rp:.4f}  Spearman {rs:.4f}")
cueP, cueS = np.array(cueP), np.array(cueS)

# pürüzsüz fit: r(N) = r∞ + c1/N + c2/N²
def smooth_fit(rvals):
    X = np.vstack([np.ones_like(Ns, dtype=float), 1 / Ns, 1 / Ns**2]).T
    beta, res, *_ = np.linalg.lstsq(X, rvals, rcond=None)
    fit = X @ beta
    rms = np.sqrt(np.mean((rvals - fit)**2))
    return beta, rms

bP, rmsP = smooth_fit(cueP)
bS, rmsS = smooth_fit(cueS)
print(f"\nPürüzsüz fit r(N) = r∞ + c₁/N + c₂/N²:")
print(f"  Pearson : r∞ = {bP[0]:.4f}, c₁ = {bP[1]:+.3f}, c₂ = {bP[2]:+.3f}  (fit RMS {rmsP:.4f})")
print(f"  Spearman: r∞ = {bS[0]:.4f}, c₁ = {bS[1]:+.3f}, c₂ = {bS[2]:+.3f}  (fit RMS {rmsS:.4f})")

def inv_smooth(beta, r_target):
    # fit N* = −2c₂/c₁ ≈ 4.6'da tepe yapar; monoton azalan dala (N>5) bak
    f = lambda N: beta[0] + beta[1] / N + beta[2] / N**2 - r_target
    return brentq(f, 5.0, 80)

# ---------------------------------------------------------------
# 2) ZETA pencereleri: Pearson + Spearman
# ---------------------------------------------------------------
def win_stats(gaps, amps, tmid):
    L = np.log(tmid / TWO_PI)
    g = gaps * L / TWO_PI
    a = amps / np.sqrt(A * L + B0 + B1 / L)
    rp, _ = pearsonr(g, a)
    rs, _ = spearmanr(g, a)
    return L.mean(), rp, rs, len(g)

rows = []
d36 = np.load(HERE / "36_T100k.npz")
edges = np.geomspace(d36["t_mid"][0], d36["t_mid"][-1] * 1.0001, 13)
for i in range(12):
    m = (d36["t_mid"] >= edges[i]) & (d36["t_mid"] < edges[i + 1])
    if m.sum() < 500:
        continue
    rows.append(win_stats(d36["intervals"][m], d36["max_amps"][m], d36["t_mid"][m]))
d41 = np.load(HERE / "41_bigT_windows.npz")
keys = sorted({k.split("_")[1] for k in d41.files}, key=lambda s: int(s[:-1]))
for k in keys:
    rows.append(win_stats(d41[f"gaps_{k}"], d41[f"amps_{k}"], d41[f"tmid_{k}"]))

zw_L = np.array([r[0] for r in rows])
zw_P = np.array([r[1] for r in rows])
zw_S = np.array([r[2] for r in rows])
zw_n = np.array([r[3] for r in rows])
seP = (1 - zw_P**2) / np.sqrt(zw_n)
seS = (1 - zw_S**2) / np.sqrt(zw_n)

# ---------------------------------------------------------------
# 3) N_eff (iki metrik) + trend testi
# ---------------------------------------------------------------
def neff_analysis(zr, se, beta, label):
    N_eff = np.array([inv_smooth(beta, r) for r in zr])
    dN = 1e-3
    drdN = np.abs(np.array([(beta[1] * (-1 / N**2) + beta[2] * (-2 / N**3)) for N in N_eff]))
    se_N = se / drdN
    shift = N_eff - zw_L
    w = 1 / se_N**2
    c_w = np.sum(w * shift) / np.sum(w)
    chi2_c = np.sum(w * (shift - c_w)**2)
    p = np.polyfit(zw_L, shift, 1, w=np.sqrt(w))
    chi2_l = np.sum(w * (shift - np.polyval(p, zw_L))**2)
    print(f"\n{label}:")
    print(f"  {'L':>5} {'N_eff−L':>8} {'±':>5}")
    for i in range(len(zw_L)):
        print(f"  {zw_L[i]:>5.2f} {shift[i]:>+8.2f} {se_N[i]:>5.2f}")
    print(f"  sabit: {c_w:+.3f} (χ²/dof {chi2_c/(len(shift)-1):.2f}) | "
          f"lineer eğim: {p[0]:+.4f}/L (χ²/dof {chi2_l/(len(shift)-2):.2f}) | "
          f"Δχ² = {chi2_c-chi2_l:.1f}")
    return shift, se_N, c_w, p

shP, seNP, cP, pP = neff_analysis(zw_P, seP, bP, "PEARSON N_eff − L")
shS, seNS, cS, pS = neff_analysis(zw_S, seS, bS, "SPEARMAN N_eff − L")

# ---------------------------------------------------------------
# 4) GRAFİK
# ---------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
for ax, sh, seN, c_w, p, ttl in [
    (axes[0], shP, seNP, cP, pP, "Pearson"),
    (axes[1], shS, seNS, cS, pS, "Spearman"),
]:
    ax.errorbar(zw_L, sh, yerr=seN, fmt="o", ms=5,
                c="firebrick" if ttl == "Pearson" else "teal")
    ax.axhline(c_w, color="gray", ls="--", label=f"sabit {c_w:+.2f}")
    xs = np.linspace(zw_L.min(), zw_L.max(), 50)
    ax.plot(xs, np.polyval(p, xs), "g-", lw=1, label=f"lineer ({p[0]:+.3f}/L)")
    ax.set_xlabel("L"); ax.set_ylabel("N_eff − L")
    ax.set_title(f"{ttl} tabanlı N_eff"); ax.legend(fontsize=9); ax.grid(alpha=0.3)

plt.tight_layout()
out = HERE / "42b_robustluk.png"
plt.savefig(out, dpi=110)
print(f"\nGrafik: {out.name}")
