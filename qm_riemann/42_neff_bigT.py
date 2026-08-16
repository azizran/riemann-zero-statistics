"""
42 — N_eff BÜYÜK T ANALİZİ (16 Ağustos 2026)
==============================================

Veri: 36 (100k aralık, L=5.6–9.1) + 41 (240k aralık, L=9.9–12.5).
Soru: N_eff − L ≈ +0.84 kayması L=12.5'e kadar SABİT mi?

  Sabit  → "ζ, N=L+c boyutlu CUE gibi korele" (temiz, tek-parametreli gözlem)
  Trendli→ kayma yapısal değil, başka bir şey oluyor (dürüstçe raporlanır)
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import pearsonr

rng = np.random.default_rng(42)
HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
A = (np.e**2 - 5) / 2
B0, B1 = 2.7580, -0.0543  # HLP-C yerel ölçek (38'de doğrulandı)

# ---------------------------------------------------------------
# 1) ZETA PENCERELERİ (unfold + Pearson r)
# ---------------------------------------------------------------
def unfolded_r(gaps, amps, tmid):
    L = np.log(tmid / TWO_PI)
    g = gaps * L / TWO_PI
    a = amps / np.sqrt(A * L + B0 + B1 / L)
    r, _ = pearsonr(g, a)
    return L.mean(), r, (1 - r**2) / np.sqrt(len(g)), len(g)

rows = []
# eski veri: 12 log-pencere
d36 = np.load(HERE / "36_T100k.npz")
edges = np.geomspace(d36["t_mid"][0], d36["t_mid"][-1] * 1.0001, 13)
for i in range(12):
    m = (d36["t_mid"] >= edges[i]) & (d36["t_mid"] < edges[i + 1])
    if m.sum() < 500:
        continue
    rows.append(unfolded_r(d36["intervals"][m], d36["max_amps"][m], d36["t_mid"][m])
                + ("36",))
# yeni pencereler
d41 = np.load(HERE / "41_bigT_windows.npz")
keys = sorted({k.split("_")[1] for k in d41.files}, key=lambda s: int(s[:-1]))
for k in keys:
    rows.append(unfolded_r(d41[f"gaps_{k}"], d41[f"amps_{k}"], d41[f"tmid_{k}"])
                + ("41",))

zw_L = np.array([r[0] for r in rows])
zw_r = np.array([r[1] for r in rows])
zw_se = np.array([r[2] for r in rows])
src = np.array([r[4] for r in rows])

print("ZETA pencereleri (unfold):")
for Lv, rv, sev, nv, sv in rows:
    print(f"  L={Lv:5.2f}  r={rv:.4f} ±{sev:.4f}  (n={nv}, {sv})")

# ---------------------------------------------------------------
# 2) CUE r(N), N=5..19
# ---------------------------------------------------------------
def haar_unitary_batch(M, N):
    G = (rng.standard_normal((M, N, N)) + 1j * rng.standard_normal((M, N, N))) / np.sqrt(2)
    Q, R = np.linalg.qr(G)
    diag = np.einsum("mii->mi", R)
    return Q * (diag / np.abs(diag))[:, None, :]

def cue_r(N, n_gaps=150000, grid=48, chunk=1200):
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
    r, _ = pearsonr(g, a)
    return r, (1 - r**2) / np.sqrt(len(g))

print("\nCUE r(N):")
Ns = np.arange(5, 20)
cue = np.array([cue_r(N) for N in Ns])
for N, (r, se) in zip(Ns, cue):
    print(f"  N={N:>2}: r = {r:.4f} ± {se:.4f}")

# ---------------------------------------------------------------
# 3) N_eff ve TREND
# ---------------------------------------------------------------
N_eff = np.interp(-zw_r, -cue[:, 0], Ns.astype(float))
# se(N_eff) = se(r) / |dr/dN| (yerel eğim)
drdN = np.abs(np.gradient(cue[:, 0], Ns))
se_Neff = zw_se / np.interp(N_eff, Ns, drdN)
shift = N_eff - zw_L

print("\nN_eff analizi:")
print(f"  {'L':>5} {'r_ζ':>8} {'N_eff':>7} {'N_eff−L':>8} {'±':>6}")
for i in range(len(zw_L)):
    print(f"  {zw_L[i]:>5.2f} {zw_r[i]:>8.4f} {N_eff[i]:>7.2f} "
          f"{shift[i]:>+8.2f} {se_Neff[i]:>6.2f}")

# ağırlıklı sabit + trend fiti
w = 1 / se_Neff**2
c_w = np.sum(w * shift) / np.sum(w)
chi2_const = np.sum(w * (shift - c_w)**2)
p = np.polyfit(zw_L, shift, 1, w=np.sqrt(w))
chi2_lin = np.sum(w * (shift - np.polyval(p, zw_L))**2)
dof_c, dof_l = len(shift) - 1, len(shift) - 2
print(f"\n  SABİT model : N_eff−L = {c_w:+.3f},  χ²/dof = {chi2_const/dof_c:.2f}")
print(f"  LİNEER model: eğim = {p[0]:+.4f}/L,  χ²/dof = {chi2_lin/dof_l:.2f}")
print(f"  Eğim anlamlılığı: Δχ² = {chi2_const - chi2_lin:.1f} (1 dof; >4 ise trend anlamlı)")

# ---------------------------------------------------------------
# 4) GRAFİK
# ---------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
ax = axes[0]
for sv, c, lbl in [("36", "firebrick", "ζ (36 verisi)"), ("41", "darkorange", "ζ (41, yeni)")]:
    m = src == sv
    ax.errorbar(zw_L[m], zw_r[m], yerr=2 * zw_se[m], fmt="o", ms=5, c=c, label=lbl)
ax.errorbar(Ns, cue[:, 0], yerr=2 * cue[:, 1], fmt="s-", ms=4, c="steelblue", label="CUE")
ax.set_xlabel("L  ↔  N"); ax.set_ylabel("Pearson r")
ax.set_title("ζ vs CUE — L=12.5'e kadar"); ax.legend(fontsize=9); ax.grid(alpha=0.3)

ax = axes[1]
for sv, c in [("36", "firebrick"), ("41", "darkorange")]:
    m = src == sv
    ax.errorbar(zw_L[m], shift[m], yerr=se_Neff[m], fmt="o", ms=5, c=c)
ax.axhline(c_w, color="gray", ls="--", label=f"sabit fit {c_w:+.3f}")
xs = np.linspace(zw_L.min(), zw_L.max(), 50)
ax.plot(xs, np.polyval(p, xs), "g-", lw=1,
        label=f"lineer fit (eğim {p[0]:+.3f}/L)")
ax.set_xlabel("L = log(t/2π)"); ax.set_ylabel("N_eff − L")
ax.set_title("Kayma sabit mi?"); ax.legend(fontsize=9); ax.grid(alpha=0.3)

plt.tight_layout()
out = HERE / "42_neff_bigT.png"
plt.savefig(out, dpi=110)
print(f"\nGrafik: {out.name}")
