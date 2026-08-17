"""
48 — ORTAK-FAZLI TAM MODEL: HALKA KAPANIYOR MU? (17 Ağustos 2026)
===================================================================

46 negatifti (bağımsız faz). 47 eksik kanalı buldu (v: boşluk sürüşü).
Şimdi tam model — her aralıkta TEK dalga fazı φ_p, iki kanala birden:

  log g̃' = log g̃_C + S_v            S_v = Σ_p v_p p^{-1/2} cos φ_p
  log ã' = log m̃_C + β·S_v + S_w    S_w = Σ_p w_p p^{-1/2} cos φ_p

β = CUE'da log-genliğin log-boşluğa regresyon eğimi (esneyen boşluk
genliği matris bağıntısı boyunca taşır; 45'in g̃-kontrollü w'si tam bu
tanımla uyumlu — w SAF genlik kanalıdır).

Varyantlar:
  V1 esneme          : yukarıdaki gibi
  V2 esneme+sessizlik: matris kısmı α ile ölçeklenir, α toplam
                       Var(log ã_ζ)'yi tam tutturacak şekilde (ölçümden)

Karşılaştırma: r, ΔVar(g̃>1.3) profili, log-varyans bütçeleri (genlik VE
boşluk), + öz-tutarlılık (simüle veriden w₂ ve v₂ geri ölçümü).
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import pearsonr

rng = np.random.default_rng(48)
HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
A = (np.e**2 - 5) / 2
B0, B1 = 2.7580, -0.0543
PRIMES = np.array([2, 3, 5, 7, 11, 13], dtype=float)
# L=12.45 penceresi kanal ölçümleri (45 ve 47)
W_CH = np.array([0.782, 0.674, 0.544, 0.454, 0.345, 0.311])
V_CH = np.array([0.111, 0.179, 0.263, 0.312, 0.381, 0.408])
AMP_W = W_CH * PRIMES**-0.5
AMP_V = V_CH * PRIMES**-0.5

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
# 1) ZETA gözlemleri (L=12.45, saklı veriden)
# ---------------------------------------------------------------
d41 = np.load(HERE / "41_bigT_windows.npz")
gaps = d41["gaps_1600k"]; tmid = d41["tmid_1600k"]; amps = d41["amps_1600k"]
Lw = np.log(tmid / TWO_PI); L = float(Lw.mean())
g_z = gaps * Lw / TWO_PI
a_z = amps / np.sqrt(A * Lw + B0 + B1 / Lw)
a_z /= np.sqrt((a_z**2).mean())
r_z, _ = pearsonr(g_z, a_z)
Vla_z = float(np.log(a_z).var())
Vlg_z = float(np.log(g_z).var())

# ---------------------------------------------------------------
# 2) CUE tabanı (N=12/13 karışımı)
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
    return np.concatenate(gs), np.concatenate(ms)

print("CUE örnekleniyor...")
frac = L - 12
g12, m12 = cue_sample(12); g13, m13 = cue_sample(13)
k12 = int((1 - frac) * len(g12)); k13 = int(frac * len(g13))
gC = np.concatenate([g12[:k12] * 12 / TWO_PI, g13[:k13] * 13 / TWO_PI])
mC = np.concatenate([m12[:k12], m13[:k13]])
mC /= np.sqrt((mC**2).mean())
r_C, _ = pearsonr(gC, mC)
lg, lm = np.log(gC), np.log(mC)
beta = float(np.cov(lg, lm)[0, 1] / lg.var())
print(f"CUE: r={r_C:.4f}, Var(log m̃)={lm.var():.4f}, Var(log g̃)={lg.var():.4f}, β={beta:.4f}")

# ---------------------------------------------------------------
# 3) SİMÜLASYON — ortak faz
# ---------------------------------------------------------------
phi = rng.uniform(0, TWO_PI, (len(gC), len(PRIMES)))
S_v = (np.cos(phi) * AMP_V).sum(axis=1)
S_w = (np.cos(phi) * AMP_W).sum(axis=1)

def make(alpha):
    gS = np.exp(lg + S_v)
    aS = np.exp(alpha * (lm - lm.mean()) + lm.mean() + beta * S_v + S_w)
    aS /= np.sqrt((aS**2).mean())
    return gS, aS

# V2 için α: Var(log ã) hedefe otursun (β kanalı + S_w + α²·matris)
var_common = np.var(beta * S_v + S_w)
alpha2 = (Vla_z - var_common) / lm.var()
alpha = np.sqrt(max(alpha2, 0.0))

variants = {"V1 esneme": make(1.0), "V2 esneme+sessizlik": make(alpha)}

# ---------------------------------------------------------------
# 4) KARŞILAŞTIRMA
# ---------------------------------------------------------------
vZ = bin_var(g_z, a_z); vC = bin_var(gC, mC)
ok = ~np.isnan(vZ) & ~np.isnan(vC) & (BC < 2.2)
big = ok & (BC > 1.3)

print(f"\nα² = {alpha2:.4f} → α = {alpha:.4f} (V2)")
print(f"\n{'':>22} {'r':>8} {'r−r_CUE':>9} {'ΔVar_büyük':>11} {'Var(log ã)':>11} {'Var(log g̃)':>11}")
print(f"{'ζ gerçek':>22} {r_z:>8.4f} {r_z - r_C:>+9.4f} "
      f"{np.nanmean((vZ - vC)[big]):>+11.4f} {Vla_z:>11.4f} {Vlg_z:>11.4f}")
print(f"{'CUE (taban)':>22} {r_C:>8.4f} {'—':>9} {'—':>11} {lm.var():>11.4f} {lg.var():>11.4f}")
for name, (gS, aS) in variants.items():
    r_s, _ = pearsonr(gS, aS)
    vS = bin_var(gS, aS)
    print(f"{name:>22} {r_s:>8.4f} {r_s - r_C:>+9.4f} "
          f"{np.nanmean((vS - vC)[big]):>+11.4f} {np.log(aS).var():>11.4f} {np.log(gS).var():>11.4f}")

# öz-tutarlılık: simüle veriden w2, v2 geri ölçümü (V1)
gS, aS = variants["V1 esneme"]
X = np.column_stack([np.ones(len(gS)), gS, gS**2,
                     np.cos(phi[:, 0]), np.sin(phi[:, 0])])
b_ = np.linalg.lstsq(X, np.log(aS), rcond=None)[0]
w2_back = b_[3] / 2**-0.5
Xg = np.column_stack([np.ones(len(gS)), np.cos(phi[:, 0]), np.sin(phi[:, 0])])
bg_ = np.linalg.lstsq(Xg, np.log(gS), rcond=None)[0]
v2_back = np.hypot(bg_[1], bg_[2]) / 2**-0.5
print(f"\nÖz-tutarlılık (V1): girilen w₂={W_CH[0]:.3f} → geri ölçülen {w2_back:.3f}; "
      f"girilen v₂={V_CH[0]:.3f} → geri ölçülen {v2_back:.3f}")

# profil grafiği
fig, ax = plt.subplots(figsize=(8.5, 5.2))
ax.plot(BC[ok], (vZ - vC)[ok], "o-", ms=5, c="firebrick", label="ζ gözlem")
for (name, (gS, aS)), c in zip(variants.items(), ["steelblue", "darkgreen"]):
    vS = bin_var(gS, aS)
    okS = ok & ~np.isnan(vS)
    ax.plot(BC[okS], (vS - vC)[okS], "s--", ms=4, c=c, label=name)
ax.axhline(0, color="k", lw=0.7)
ax.set_xlabel("g̃"); ax.set_ylabel("Var − Var_CUE")
ax.set_title(f"Ortak-fazlı model vs gözlem (L={L:.2f})")
ax.legend(); ax.grid(alpha=0.3)
plt.tight_layout()
out = HERE / "48_ortak_faz.png"
plt.savefig(out, dpi=110)
print(f"\nGrafik: {out.name}")
