"""
75 — AYNA KANALININ SESİ: ÇADIR ŞARKISI (18 Ağustos 2026, gece)
=================================================================

74'ün çözümü: noktasal |Z|² modülasyonu kaydırılmış-ikinci-moment ana
terimini izlemeli (BCHB-tipi):
  R_tam(τ) = 2·(L − log p + 2γ − 1) / (L + 2γ − 1)   [≈ 2(1−τ)]
— katta kesilmez, p ≈ t/2π'de (TAM ufuk) ölür.

AYRIŞIM:
  direkt çiftler = 2·S₁(N/p)/S₁(N)  [≈ 2(1−2τ), katta keskin ölür — 74 türetimi]
  AYNA kanalı  = R_tam − direkt     [≈ ÇADIR: 2τ ↑ katta 1.0 ↓ 2(1−τ)]

Test: iki pencere (L=9.86 ve 12.45), parametresiz; ek çapraz-kontrol:
c₀ serbest fit → 2γ−1 = 0.154'e yakın mı?
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
GAMMA = 0.5772156649
C0_PRED = 2 * GAMMA - 1  # 0.1544

def theta(t):
    return t / 2 * np.log(t / TWO_PI) - t / 2 - np.pi / 8 + 1 / (48 * t)

def Z_rs(t, chunk=120000):
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
        cp = np.where(np.abs(cp) < 1e-8, 1e-8, cp)
        psi = np.cos(TWO_PI * (p**2 - p - 1 / 16)) / cp
        z += (-1) ** (N - 1) * (tt / TWO_PI) ** -0.25 * psi
        out[s:s + chunk] = z
    return out

S1 = lambda M: np.sum(1.0 / np.arange(1, max(int(M), 1) + 1))

def window_R(T_LO, T_HI, PS, step=0.02):
    grid = np.arange(T_LO, T_HI, step)
    Y = Z_rs(grid) ** 2
    Y = Y / Y.mean()
    L = float(np.log(0.5 * (T_LO + T_HI) / TWO_PI))
    N_RS = int(np.sqrt(0.5 * (T_LO + T_HI) / TWO_PI))
    cols = [np.ones_like(Y)]
    for p in PS:
        arg = grid * np.log(p)
        cols += [np.cos(arg), np.sin(arg)]
    X = np.vstack(cols).T
    b, *_ = np.linalg.lstsq(X, Y, rcond=None)
    se = np.sqrt((Y - X @ b).var() * np.diag(np.linalg.inv(X.T @ X)))
    out = []
    for i, p in enumerate(PS):
        out.append((np.log(p) / L, b[1 + 2*i] / p**-0.5, se[1 + 2*i] / p**-0.5, p))
    return L, N_RS, out

PS1 = [2, 3, 5, 7, 11, 13, 17, 23, 31, 43, 59, 79, 101, 127, 137, 139, 149,
       163, 181, 199, 251, 307, 397, 499, 641, 809]
print("Pencere 1 (L≈9.86)...")
L1, N1, R1 = window_R(107252.0, 132748.0, PS1)
PS2 = [2, 3, 5, 7, 13, 23, 43, 79, 149, 251, 439, 761, 1327, 2311, 4021,
       6997, 12163, 21169]
print("Pencere 2 (L≈12.45)...")
L2, N2, R2 = window_R(1589905.0, 1610095.0, PS2, step=0.015)

def r_full(tau, L):
    return 2 * (L * (1 - tau) + C0_PRED) / (L + C0_PRED)

def r_direct(p, L, N):
    return 2 * S1(N / p) / S1(N) if p <= N else 0.0

print(f"\n{'pencere':>8} {'p':>6} {'τ':>6} {'ölçülen':>8} {'±':>6} "
      f"{'TAM(BCHB)':>9} {'direkt':>7} {'AYNA(ölç−dir)':>13} {'çadır-tahmin':>12}")
tent_res = []
for (L, N, RR), tag in [((L1, N1, R1), "9.86"), ((L2, N2, R2), "12.45")]:
    for tau, meas, err, p in RR:
        full = r_full(tau, L)
        dirc = r_direct(p, L, N)
        mirror_meas = meas - dirc
        mirror_pred = full - dirc
        tent_res.append((tau, mirror_meas, err, mirror_pred, L))
        if p in (2, 13, 43, 137, 139, 199, 499, 809, 2311, 21169):
            print(f"{tag:>8} {p:>6} {tau:>6.3f} {meas:>8.3f} {err:>6.3f} "
                  f"{full:>9.3f} {dirc:>7.3f} {mirror_meas:>13.3f} {mirror_pred:>12.3f}")

tent = np.array([(t, m, e, pr) for t, m, e, pr, L in tent_res])
resid = tent[:, 1] - tent[:, 3]
print(f"\nÇADIR uyumu ({len(tent)} nokta, İKİ pencere, PARAMETRESİZ): "
      f"artık RMS = {np.sqrt(np.mean(resid**2)):.3f} (sinyal ölçeği ~1)")

# c0 serbest çapraz-kontrol (tam-formül üzerinden)
from scipy.optimize import least_squares
allpts = [(t, m, e, L) for t, m, e, pr, L in tent_res]
meas_full = np.array([(t, m + r_direct(np.exp(t*L), L, int(np.exp(L/2))), e, L)
                      for t, m, e, L in allpts])
def res_c0(c):
    return (2*(meas_full[:,3]*(1-meas_full[:,0])+c[0])/(meas_full[:,3]+c[0])
            - meas_full[:,1]) / meas_full[:,2]
fc = least_squares(res_c0, [0.15])
print(f"c₀ serbest fit = {fc.x[0]:.3f}   (tahmin 2γ−1 = {C0_PRED:.3f})")

# grafik: çadır
fig, ax = plt.subplots(figsize=(9, 5.6))
for L, c, lab in [(L1, "firebrick", "L=9.86"), (L2, "teal", "L=12.45")]:
    m = np.array([(t, mm, e) for t, mm, e, pr, LL in tent_res if abs(LL-L) < 0.1])
    ax.errorbar(m[:, 0], m[:, 1], yerr=m[:, 2], fmt="o", ms=4, c=c,
                label=f"ayna kanalı (ölçüm − direkt), {lab}")
tt = np.linspace(0.02, 0.95, 300)
tent_c = np.where(tt < 0.5, 2*tt, 2*(1-tt))
ax.plot(tt, tent_c, "k-", lw=1.4, label="ÇADIR: 2τ ↑ kat ↓ 2(1−τ)")
ax.axvline(0.5, color="gray", ls="--", lw=1, label="kat çizgisi")
ax.axhline(0, color="gray", lw=0.6)
ax.set_xlabel("τ = log p / L"); ax.set_ylabel("ayna-kanalı genliği")
ax.set_title("Aynanın sesi: çadır şarkısı — parametresiz")
ax.legend(fontsize=9); ax.grid(alpha=0.3)
plt.tight_layout()
out = HERE / "75_ayna_sesi.png"
plt.savefig(out, dpi=110)
print(f"\nGrafik: {out.name}")
