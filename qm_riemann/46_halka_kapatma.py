"""
46 — HALKAYI KAPATMA: w(p)'DEN AÇIĞI TÜRETMEK (17 Ağustos 2026)
=================================================================

Ölçülen w spektrumundan −0.019 r açığını ve 43'ün varyans profilini
TAHMİN etmeye çalışıyoruz. İki aday model (ikisi de 0 serbest parametre):

  EKLEME : log ã_ζ = log m̃_CUE + S       (asal dalga CUE üstüne biner)
  İKAME  : log ã_ζ = α·log m̃_CUE + S     (matris kısmı sessizleşir;
           α² = [Var(log ã_ζ) − σ_S²] / Var(log m̃_CUE) — ölçümden)

S = Σ_p w_p p^{-1/2} cos φ_p, φ bağımsız uniform (tek-aralık istatistiği
için yeterli), w'ler 45'te ölçüldü.

KİLİT TEŞHİS (varyans bütçesi): Var(log ã_ζ) ölçülür;
  Var(log m̃_CUE) + σ_S² ile karşılaştırılır.
  ≈ eşitse → EKLEME. Küçükse → İKAME (sıfır gazı rastgelelik devrediyor).

Sonra iki model simüle edilir → r ve koşullu varyans profili → gerçek ζ
(L=12.45 penceresi) ile kıyas.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import pearsonr

rng = np.random.default_rng(46)
HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
A = (np.e**2 - 5) / 2
B0, B1 = 2.7580, -0.0543
PRIMES = [2, 3, 5, 7, 11, 13]

def theta(t):
    return t / 2 * np.log(t / TWO_PI) - t / 2 - np.pi / 8 + 1 / (48 * t) + 7 / (5760 * t**3)

def Z_rs(t, chunk=20000):
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
        cp = np.where(np.abs(cp) < 1e-8, 1e-8 * np.sign(cp + 1e-300), cp)
        psi = np.cos(TWO_PI * (p**2 - p - 1 / 16)) / cp
        z += (-1) ** (N - 1) * (tt / TWO_PI) ** -0.25 * psi
        out[s:s + chunk] = z
    return out

BINS = np.concatenate([[0.05], np.linspace(0.25, 2.2, 14), [3.5]])
BC = 0.5 * (BINS[:-1] + BINS[1:])

def bin_var(g, a, min_n=200):
    var = np.full(len(BC), np.nan)
    idx = np.digitize(g, BINS) - 1
    for i in range(len(BC)):
        m = idx == i
        if m.sum() >= min_n:
            var[i] = a[m].var()
    return var

# ---------------------------------------------------------------
# 1) ZETA (L=12.45): g̃, ã, t_pk, w regresyonu
# ---------------------------------------------------------------
d41 = np.load(HERE / "41_bigT_windows.npz")
gaps = d41["gaps_1600k"]; tmid = d41["tmid_1600k"]
Lw = np.log(tmid / TWO_PI); L = float(Lw.mean())
g_lo = tmid - gaps / 2
u = np.arange(1, 25) / 25
tt = g_lo[:, None] + gaps[:, None] * u[None, :]
print("Z değerlendiriliyor...")
Zg = np.abs(Z_rs(tt.ravel())).reshape(tt.shape)
j = np.argmax(Zg, axis=1)
t_pk = tt[np.arange(len(j)), j]
amps = Zg[np.arange(len(j)), j]

g_u = gaps * Lw / TWO_PI
a_u = amps / np.sqrt(A * Lw + B0 + B1 / Lw)
a_u /= np.sqrt((a_u**2).mean())
y = np.log(a_u)

cols = [np.ones_like(y), g_u, g_u**2]
for p in PRIMES:
    cols += [np.cos(t_pk * np.log(p)), np.sin(t_pk * np.log(p))]
X = np.vstack(cols).T
beta, *_ = np.linalg.lstsq(X, y, rcond=None)
w_meas = [beta[3 + 2 * i] / p**-0.5 for i, p in enumerate(PRIMES)]
amp_S = np.array([w * p**-0.5 for w, p in zip(w_meas, PRIMES)])
sigma_S2 = float(np.sum(amp_S**2) / 2)

r_zeta, _ = pearsonr(g_u, a_u)
V_zeta_log = float(y.var())
print(f"\nζ (L={L:.2f}): r = {r_zeta:.4f}, Var(log ã) = {V_zeta_log:.4f}")
print(f"w'ler: " + " ".join(f"{w:.3f}" for w in w_meas))
print(f"σ_S² (p≤13) = {sigma_S2:.4f}")

# ---------------------------------------------------------------
# 2) CUE (N=12,13 karışımı → N≈L etkisi)
# ---------------------------------------------------------------
def haar_unitary_batch(M, N):
    G = (rng.standard_normal((M, N, N)) + 1j * rng.standard_normal((M, N, N))) / np.sqrt(2)
    Q, R = np.linalg.qr(G)
    diag = np.einsum("mii->mi", R)
    return Q * (diag / np.abs(diag))[:, None, :]

def cue_sample(N, n_gaps=300000, grid=24, chunk=1200):
    M = int(np.ceil(n_gaps / N))
    gs, ms = [], []
    for s in range(0, M, chunk):
        m = min(chunk, M - s)
        U = haar_unitary_batch(m, N)
        ph = np.sort(np.angle(np.linalg.eigvals(U)), axis=1)
        gp = np.diff(np.concatenate([ph, ph[:, :1] + TWO_PI], axis=1), axis=1)
        uu = np.arange(1, grid + 1) / (grid + 1)
        th = ph[:, :, None] + gp[:, :, None] * uu[None, None, :]
        df = th[:, :, :, None] - ph[:, None, None, :]
        amp = np.prod(2 * np.abs(np.sin(df / 2)), axis=-1)
        gs.append(gp.ravel()); ms.append(amp.max(axis=-1).ravel())
    g = np.concatenate(gs); a = np.concatenate(ms)
    return g, a

print("\nCUE örnekleniyor (N=12 ve 13)...")
frac = L - 12  # 0.45
g12, m12 = cue_sample(12); g13, m13 = cue_sample(13)
n13 = int(frac * len(g13))
gC = np.concatenate([g12[:int((1 - frac) * len(g12))] * 12 / TWO_PI,
                     g13[:n13] * 13 / TWO_PI])
mC = np.concatenate([m12[:int((1 - frac) * len(g12))], m13[:n13]])
mC /= np.sqrt((mC**2).mean())
r_cue, _ = pearsonr(gC, mC)
lC = np.log(mC)
V_cue_log = float(lC.var())
print(f"CUE(N≈{L:.2f}): r = {r_cue:.4f}, Var(log m̃) = {V_cue_log:.4f}")

# ---------------------------------------------------------------
# 3) VARYANS BÜTÇESİ — kilit teşhis
# ---------------------------------------------------------------
print("\n=== VARYANS BÜTÇESİ (log domeni) ===")
print(f"  Var(log ã_ζ)              = {V_zeta_log:.4f}")
print(f"  Var(log m̃_CUE) + σ_S²     = {V_cue_log:.4f} + {sigma_S2:.4f} = {V_cue_log + sigma_S2:.4f}")
print(f"  → EKLEME modeli öngörüsü {V_cue_log + sigma_S2:.3f}; ölçülen {V_zeta_log:.3f}")
alpha2 = (V_zeta_log - sigma_S2) / V_cue_log
print(f"  → İKAME modeli: α² = {alpha2:.4f}, α = {np.sqrt(max(alpha2, 0)):.4f}")
print(f"    (α<1: sıfır gazı matris-gürültüsünün %{100*(1-alpha2):.0f}'ini asal düzenine devretmiş)")

# ---------------------------------------------------------------
# 4) İKİ MODELİN SİMÜLASYONU → r ve varyans profili
# ---------------------------------------------------------------
S = np.zeros(len(mC))
for a_j in amp_S:
    S += a_j * np.cos(rng.uniform(0, TWO_PI, len(mC)))

sims = {}
# EKLEME
Xa = np.exp(lC + S); Xa /= np.sqrt((Xa**2).mean())
sims["EKLEME"] = Xa
# İKAME
alpha = np.sqrt(max(alpha2, 0))
Xi = np.exp(alpha * (lC - lC.mean()) + lC.mean() + S)
Xi /= np.sqrt((Xi**2).mean())
sims["İKAME"] = Xi

vC = bin_var(gC, mC)
vZ = bin_var(g_u, a_u)
ok = ~np.isnan(vC) & ~np.isnan(vZ) & (BC < 2.2)
big = ok & (BC > 1.3)

print("\n=== MODEL → GÖZLEM KIYASI ===")
print(f"  {'':>8} {'r':>8} {'r−r_CUE':>9} {'ΔVar(g̃>1.3)':>13}")
print(f"  {'ζ gerçek':>8} {r_zeta:>8.4f} {r_zeta - r_cue:>+9.4f} "
      f"{np.nanmean((vZ - vC)[big]):>+13.4f}")
for name, Xs in sims.items():
    r_s, _ = pearsonr(gC, Xs)
    vS = bin_var(gC, Xs)
    print(f"  {name:>8} {r_s:>8.4f} {r_s - r_cue:>+9.4f} "
          f"{np.nanmean((vS - vC)[big]):>+13.4f}")

# profil grafiği
fig, ax = plt.subplots(figsize=(8, 5.2))
ax.plot(BC[ok], (vZ - vC)[ok], "o-", ms=5, c="firebrick", label="ζ gözlem")
for name, Xs, c in [("EKLEME", sims["EKLEME"], "steelblue"),
                    ("İKAME", sims["İKAME"], "darkgreen")]:
    vS = bin_var(gC, Xs)
    ax.plot(BC[ok], (vS - vC)[ok], "s--", ms=4, c=c, label=f"{name} modeli")
ax.axhline(0, color="k", lw=0.7)
ax.set_xlabel("g̃"); ax.set_ylabel("Var − Var_CUE")
ax.set_title(f"Koşullu varyans fazlası: gözlem vs iki model (L={L:.2f})")
ax.legend(); ax.grid(alpha=0.3)
plt.tight_layout()
out = HERE / "46_halka_kapatma.png"
plt.savefig(out, dpi=110)
print(f"\nGrafik: {out.name}")
