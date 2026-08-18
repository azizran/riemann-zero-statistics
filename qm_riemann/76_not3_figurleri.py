"""
76 — NOT 3 FİGÜRLERİ (İngilizce) (19 Ağustos 2026)
====================================================
fig_law_en: √w vs v, yasa doğrusu, out-of-sample p=11,13
fig_crossing_en: w(τ) işaret geçişi + termal-Bragg modeli
fig_tent_en: ayna kanalının çadırı (75 verisi, 2 pencere)
"""

import numpy as np
import mpmath as mp
import matplotlib.pyplot as plt
from pathlib import Path

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
A_CG = (np.e**2 - 5) / 2
B0, B1 = 2.7580, -0.0543
T0_ODL = 267653395647
mp.mp.dps = 30

def channels(gaps, amps, tmid, primes, anchored=False, controls=True):
    Lw = np.log(((T0_ODL + tmid) if anchored else tmid) / TWO_PI)
    g_u = gaps * Lw / TWO_PI
    a_u = amps / np.sqrt(A_CG * Lw + B0 + B1 / Lw)
    a_u /= np.sqrt((a_u**2).mean())
    ya, yg = np.log(a_u), np.log(g_u)
    cw = ([np.ones_like(ya), g_u, g_u**2] if controls else [np.ones_like(ya)])
    cv = [np.ones_like(yg)]
    for p in primes:
        ph0 = float(mp.fmod(T0_ODL * mp.log(p), 2 * mp.pi)) if anchored else 0.0
        arg = ph0 + tmid * np.log(p)
        cw += [np.cos(arg), np.sin(arg)]
        cv += [np.cos(arg), np.sin(arg)]
    o = 3 if controls else 1
    Xw = np.vstack(cw).T
    bw, *_ = np.linalg.lstsq(Xw, ya, rcond=None)
    sew = np.sqrt((ya - Xw @ bw).var() * np.diag(np.linalg.inv(Xw.T @ Xw)))
    Xv = np.vstack(cv).T
    bv, *_ = np.linalg.lstsq(Xv, yg, rcond=None)
    L = float(Lw.mean())
    out = []
    for i, p in enumerate(primes):
        out.append((np.log(p) / L, bw[o + 2*i] / p**-0.5, sew[o + 2*i] / p**-0.5,
                    float(np.hypot(bv[1 + 2*i], bv[2 + 2*i]) / p**-0.5), p, L))
    return out

W41, W36 = [], []
d41 = np.load(HERE / "41_bigT_windows.npz")
for k in sorted({x.split("_")[1] for x in d41.files}, key=lambda s: int(s[:-1])):
    W41.append((d41[f"gaps_{k}"], d41[f"amps_{k}"], d41[f"tmid_{k}"], False))
for f in ["55_win_1e+08.npz", "55_win_1e+09.npz", "55_win_1e+10.npz", "55_win_1e+11.npz"]:
    d = np.load(HERE / f)
    W41.append((d["gaps"], d["amps"], d["tmid"], False))
d53 = np.load(HERE / "53_odlyzko_amps.npz")
W41.append((d53["gaps"], d53["max_amps"], d53["t_mid"], True))
d36 = np.load(HERE / "36_T100k.npz")
edges = np.geomspace(d36["t_mid"][0], d36["t_mid"][-1] * 1.0001, 13)
for i in range(12):
    m = (d36["t_mid"] >= edges[i]) & (d36["t_mid"] < edges[i + 1])
    if m.sum() >= 500:
        W36.append((d36["intervals"][m], d36["max_amps"][m], d36["t_mid"][m], False))

# ---- Fig 1: the law ----
P4, P2 = [2, 3, 5, 7], [11, 13]
tr, oos = [], []
for gaps, amps, tmid, anch in W41:
    for tau, w_, s_, v_, p, L in channels(gaps, amps, tmid, P4 + P2, anch):
        if w_ > 0:
            (tr if p in P4 else oos).append((v_, np.sqrt(w_)))
tr = np.array(tr); oos = np.array(oos)
fig, ax = plt.subplots(figsize=(8.2, 5.4))
ax.plot(tr[:, 0], tr[:, 1], "o", ms=4, c="steelblue", alpha=0.6,
        label=r"$p \in \{2,3,5,7\}$ (44 pairs, law derived here)")
ax.plot(oos[:, 0], oos[:, 1], "s", ms=6, c="firebrick", zorder=5,
        label=r"$p \in \{11,13\}$ (out of sample)")
xx = np.linspace(0, 0.62, 50)
ax.plot(xx, 1 - 1.0683 * xx, "k--", lw=1.2, label=r"$\sqrt{w} = 1 - 1.068\,v$")
ax.set_xlabel(r"$v$ — zero-displacement channel")
ax.set_ylabel(r"$\sqrt{w}$ — amplitude transmission")
ax.legend(fontsize=9); ax.grid(alpha=0.3)
plt.tight_layout(); plt.savefig(HERE / "fig_law_en.png", dpi=110); plt.close()
print("fig_law_en ok")

# ---- Fig 2: crossing + thermal model ----
P11 = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
pts = []
cs = []
def rvm_N(t):
    x = t / TWO_PI
    return x * np.log(x / np.e) + 7 / 8
for gaps, amps, tmid, anch in (W36 + W41[:6]):
    row = channels(gaps, amps, tmid, P11, anch)
    pts += [(t, w_, s_) for t, w_, s_, v_, p, L in row]
    t0 = tmid[0] - gaps[0] / 2
    kk = np.arange(len(tmid)) + 0.5
    t = t0 + kk * TWO_PI / np.log(t0 / TWO_PI)
    for _ in range(6):
        fdel = rvm_N(t) - rvm_N(t0) - kk
        t = t - fdel / (np.log(t / TWO_PI) / TWO_PI)
    L = float(np.log(tmid / TWO_PI).mean())
    cs.append((L * (tmid - t).std()) ** 2 / 2)
pts = np.array(sorted(pts))
tau, w, sw = pts.T
c_jit = float(np.mean(cs))
B_PL = float(np.mean(w[tau > 0.5]))
from scipy.optimize import least_squares
fT = least_squares(lambda p: ((p[0]*(1-2*tau)*np.exp(-c_jit*tau**2)*(tau<0.5) + B_PL) - w) / sw, [1.05])
fig, ax = plt.subplots(figsize=(8.6, 5.4))
ax.errorbar(tau, w, yerr=sw, fmt="o", ms=3, c="firebrick", alpha=0.55,
            label=r"$w(\tau)$ (132 measurements)")
tt = np.linspace(0.001, 0.62, 500)
ax.plot(tt, fT.x[0]*(1-2*tt)*np.exp(-c_jit*tt**2)*(tt < 0.5) + B_PL, "k-", lw=1.4,
        label=rf"thermal Bragg: $A(1-2\tau)e^{{-c\tau^2}} + B$ ($A$={fT.x[0]:.2f})")
ax.axvline(0.5, color="gray", ls="--", lw=1, label=r"ideal mirror $\tau = 1/2$")
ax.axhline(0, color="gray", lw=0.6)
ax.set_xlabel(r"$\tau = \log p / L$"); ax.set_ylabel(r"$w$")
ax.legend(fontsize=9); ax.grid(alpha=0.3)
plt.tight_layout(); plt.savefig(HERE / "fig_crossing_en.png", dpi=110); plt.close()
print(f"fig_crossing_en ok (c_jit={c_jit:.2f}, B={B_PL:+.3f}, A={fT.x[0]:.3f})")

# ---- Fig 3: the tent (75'in makinesi, hafif grid) ----
def theta(t):
    return t / 2 * np.log(t / TWO_PI) - t / 2 - np.pi / 8 + 1 / (48 * t)
def Z_rs(t, chunk=120000):
    t = np.asarray(t, dtype=np.float64)
    out = np.empty_like(t)
    for s in range(0, len(t), chunk):
        tt2 = t[s:s + chunk]
        a = np.sqrt(tt2 / TWO_PI)
        N = a.astype(np.int64)
        th = theta(tt2)
        z = np.zeros_like(tt2)
        for Nv in np.unique(N):
            m = N == Nv
            n = np.arange(1, Nv + 1)
            ph = th[m, None] - tt2[m, None] * np.log(n)[None, :]
            z[m] = 2 * (np.cos(ph) @ (n**-0.5))
        pfr = a - N
        cp = np.cos(TWO_PI * pfr)
        cp = np.where(np.abs(cp) < 1e-8, 1e-8, cp)
        z += (-1.0)**(N - 1) * (tt2 / TWO_PI)**-0.25 * np.cos(TWO_PI*(pfr*pfr - pfr - 1/16)) / cp
        out[s:s + chunk] = z
    return out
S1f = lambda M: np.sum(1.0 / np.arange(1, max(int(M), 1) + 1))
def windowR(T_LO, T_HI, PS, step):
    grid = np.arange(T_LO, T_HI, step)
    Y = Z_rs(grid)**2; Y = Y / Y.mean()
    L = float(np.log(0.5*(T_LO+T_HI)/TWO_PI))
    N_RS = int(np.sqrt(0.5*(T_LO+T_HI)/TWO_PI))
    cols = [np.ones_like(Y)]
    for p in PS:
        arg = grid*np.log(p); cols += [np.cos(arg), np.sin(arg)]
    X = np.vstack(cols).T
    b, *_ = np.linalg.lstsq(X, Y, rcond=None)
    se = np.sqrt((Y - X@b).var()*np.diag(np.linalg.inv(X.T@X)))
    return [(np.log(p)/L, b[1+2*i]/p**-0.5, se[1+2*i]/p**-0.5,
             (2*S1f(N_RS/p)/S1f(N_RS) if p <= N_RS else 0.0)) for i, p in enumerate(PS)]
PSa = [2,3,5,7,11,13,17,23,31,43,59,79,101,127,137,139,149,163,181,199,251,307,397,499]
Ra = windowR(107252.0, 132748.0, PSa, 0.02)
PSb = [2,3,5,7,13,23,43,79,149,251,439,761,1327,2311,4021]
Rb = windowR(1589905.0, 1610095.0, PSb, 0.015)
fig, ax = plt.subplots(figsize=(8.6, 5.4))
for R, c, lab in [(Ra, "firebrick", r"$L=9.86$"), (Rb, "teal", r"$L=12.45$")]:
    m = np.array([(t, meas - dirc, e) for t, meas, e, dirc in R if t < 0.7])
    ax.errorbar(m[:, 0], m[:, 1], yerr=m[:, 2], fmt="o", ms=4, c=c,
                label=f"mirror channel (measured − direct), {lab}")
tt = np.linspace(0.02, 0.75, 300)
ax.plot(tt, np.where(tt < 0.5, 2*tt, 2*(1-tt)), "k-", lw=1.3,
        label=r"tent: $2\tau$ up to the fold, $2(1-\tau)$ beyond")
ax.axvline(0.5, color="gray", ls="--", lw=1, label=r"the fold $\tau = 1/2$")
ax.axhline(0, color="gray", lw=0.6)
ax.set_xlabel(r"$\tau = \log p / L$"); ax.set_ylabel("mirror-channel amplitude")
ax.legend(fontsize=9); ax.grid(alpha=0.3)
plt.tight_layout(); plt.savefig(HERE / "fig_tent_en.png", dpi=110); plt.close()
print("fig_tent_en ok")
