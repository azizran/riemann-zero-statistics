"""
94 — NOT 4 FİGÜRLERİ (İngilizce) (21 Ağustos 2026)
==========================================================================
fig_diffraction_en : benekler (mutlak EF öngörüsü vs ölçüm) + karanlık
                     alan (üç ızgara) — 86/89 makinesi
fig_geometry_en    : benek geometrisi πτ·cot(πτ) çökmesi — 88-M1
fig_bridge_en      : v-toplam kuralı + mutlak benek yasası — 89
fig_twobranch_en   : iki-dal yanıtı (asal perdesi / CUE / spontane) —
                     92-93 (CUE M=1200 hızlı tekrar; D_spont 93 ızgarası)
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
P11 = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]

def Ghat(t, omegas, chunk=30000):
    out = np.zeros(len(omegas), dtype=complex)
    for s0 in range(0, len(t), chunk):
        tt = t[s0:s0 + chunk]
        out += np.exp(1j * np.outer(omegas, tt)).sum(axis=1)
    return out / len(t)

def rvm_N(t):
    x = t / TWO_PI
    return x * np.log(x / np.e) + 7 / 8

def sigma_t(gaps, tmid):
    t0 = tmid[0] - gaps[0] / 2
    kk = np.arange(len(tmid)) + 0.5
    ts = t0 + kk * TWO_PI / np.log(t0 / TWO_PI)
    for _ in range(6):
        ts = ts - (rvm_N(ts) - rvm_N(t0) - kk) / (np.log(ts / TWO_PI) / TWO_PI)
    return float((tmid - ts).std()), ts

d41 = np.load(HERE / "41_bigT_windows.npz")
K41 = sorted({x.split("_")[1] for x in d41.files}, key=lambda s: int(s[:-1]))
WNDS = [(d41[f"gaps_{k}"], d41[f"amps_{k}"], d41[f"tmid_{k}"]) for k in K41[:6]]
d83 = np.load(HERE / "83_tam_taban_egri.npz")

# ============ A: fig_diffraction_en ============
meas, pred, taus = [], [], []
for p in P11:
    ms, prs, ts_ = [], [], []
    for gaps, amps, tmid in WNDS:
        L = float(np.log(tmid / TWO_PI).mean())
        tau = np.log(p) / L
        sig, _ = sigma_t(gaps, tmid)
        dw = np.exp(-np.log(p)**2 * sig**2 / 2)
        ms.append(abs(Ghat(tmid, np.array([np.log(p)]))[0]))
        prs.append(tau * p**-0.5 * np.cos(np.pi * tau) * dw)
        ts_.append(tau)
    meas.append(np.mean(ms)); pred.append(np.mean(prs)); taus.append(np.mean(ts_))
gaps, amps, tmid = WNDS[1]
n = len(tmid)
rng = np.random.default_rng(94)
g_sh = rng.permutation(gaps)
t_sh = tmid[0] + np.cumsum(g_sh) - g_sh / 2
_, ts_smooth = sigma_t(gaps, tmid)
lines = [np.log(q) for q in
         [2,3,4,5,7,8,9,11,13,16,17,19,23,25,27,29,31,32,37,41,43,47,49,53,59,61,64]]
om_scan = np.linspace(0.15, 4.2, 1200)
om_scan = np.array([o for o in om_scan if min(abs(o - l) for l in lines) > 0.015])
G_real = np.abs(Ghat(tmid, om_scan))
G_shuf = np.abs(Ghat(t_sh, om_scan))
G_rvm = np.abs(Ghat(ts_smooth, om_scan))

fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.0))
ax = axes[0]
ax.plot(taus, pred, "s-", c="steelblue", ms=6,
        label=r"explicit formula: $\tau p^{-1/2}\cos(\pi\tau)\cdot\mathrm{DW}$")
ax.plot(taus, meas, "o", c="firebrick", ms=6, label="measured spot $|\\hat G(\\log p)|$")
ax.set_xlabel(r"$\tau = \log p / L$"); ax.set_ylabel(r"$|\hat G|$")
ax.set_title("Arithmetic spots: parameter-free prediction")
ax.legend(fontsize=9); ax.grid(alpha=0.3)
ax = axes[1]
ax.semilogy(om_scan, G_real, ".", ms=2, c="firebrick", alpha=0.6, label="real grid")
ax.semilogy(om_scan, G_shuf, ".", ms=2, c="steelblue", alpha=0.5, label="shuffled-gap surrogate")
ax.semilogy(om_scan, G_rvm, ".", ms=2, c="gray", alpha=0.5, label="smooth (RvM) grid")
ax.axhline(1/np.sqrt(n), color="k", lw=0.8, ls="--", label=r"$1/\sqrt{n}$")
ax.set_xlabel(r"$\omega$ (off-line)"); ax.set_ylabel(r"$|\hat G(\omega)|$")
ax.set_title("Dark field: hyperuniformity in diffraction")
ax.legend(fontsize=9); ax.grid(alpha=0.3)
plt.tight_layout(); plt.savefig(HERE / "fig_diffraction_en.png", dpi=110); plt.close()
print("fig_diffraction_en ok")

# ============ B: fig_geometry_en ============
rats, fss, tts = [], [], []
for p in P11:
    rs, fs2, t2 = [], [], []
    for gaps_, amps_, tmid_ in WNDS:
        L = float(np.log(tmid_ / TWO_PI).mean())
        m = (np.abs(d83["L"] - L) < 0.01) & (np.abs(d83["tau"] - np.log(p)/L) < 1e-9)
        if m.sum() != 1:
            continue
        v = float(d83["v"][m][0])
        tau = np.log(p) / L
        rs.append(abs(Ghat(tmid_, np.array([np.log(p)]))[0]) / (v * p**-0.5 / 2))
        fs2.append(np.pi * tau / np.tan(np.pi * tau)); t2.append(tau)
    rats.append(np.mean(rs)); fss.append(np.mean(fs2)); tts.append(np.mean(t2))
fig, ax = plt.subplots(figsize=(8.2, 5.2))
tt = np.linspace(0.02, 0.49, 200)
ax.plot(tt, np.pi * tt / np.tan(np.pi * tt), "k-", lw=1.3,
        label=r"discrete geometry: $f(\tau) = \pi\tau\cot(\pi\tau)$")
ax.plot(tts, rats, "o", c="firebrick", ms=6, label="measured spot / naive $v$-prediction")
ax.axvline(0.5, color="gray", ls="--", lw=1, label=r"Nyquist: $f(1/2)=0$")
ax.set_xlabel(r"$\tau$"); ax.set_ylabel("ratio")
ax.set_title(r"Spot attenuation is pure sampling geometry ($\times\,k_0 = 1.083$)")
ax.legend(fontsize=9); ax.grid(alpha=0.3)
plt.tight_layout(); plt.savefig(HERE / "fig_geometry_en.png", dpi=110); plt.close()
print("fig_geometry_en ok")

# ============ C: fig_bridge_en ============
tauv, vv = d83["tau"], d83["v"]
A_v = np.sum(vv * np.sin(np.pi * tauv)) / np.sum(np.sin(np.pi * tauv)**2)
fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.0))
ax = axes[0]
ax.plot(tauv, vv, "o", ms=4, c="steelblue", alpha=0.55, label="measured $v$ (132 points)")
tt = np.linspace(0, 0.62, 200)
ax.plot(tt, (2/np.pi) * np.sin(np.pi * tt), "k-", lw=1.4,
        label=r"derived: $v = \frac{2}{\pi}\sin(\pi\tau)$")
ax.plot(tt, A_v * np.sin(np.pi * tt), "k--", lw=1,
        label=rf"free amplitude: ${A_v:.3f}\,\sin(\pi\tau)$")
ax.plot(tt, 2 * tt, ":", c="gray", lw=1, label=r"derived onset $2\tau$")
ax.axhline(2/np.pi, color="firebrick", lw=0.8, ls=":", label=r"saturation $2/\pi$")
ax.set_xlabel(r"$\tau$"); ax.set_ylabel("$v$")
ax.set_title("The gap-channel sum rule")
ax.legend(fontsize=8.5); ax.grid(alpha=0.3)
ax = axes[1]
ax.plot(taus, pred, "s-", c="steelblue", ms=6, label="explicit-formula absolute prediction")
ax.plot(taus, meas, "o", c="firebrick", ms=6, label="measured spot")
ax.set_xlabel(r"$\tau$"); ax.set_ylabel(r"$|\hat G(\log p)|$")
ax.set_title("Absolute spot law (no channel input)")
ax.legend(fontsize=9); ax.grid(alpha=0.3)
plt.tight_layout(); plt.savefig(HERE / "fig_bridge_en.png", dpi=110); plt.close()
print("fig_bridge_en ok")

# ============ D: fig_twobranch_en ============
gaps, amps, tmid = WNDS[1]
L0 = float(np.log(tmid / TWO_PI).mean())
z = np.empty(len(gaps) + 1)
z[0] = tmid[0] - gaps[0] / 2
z[1:] = z[0] + np.cumsum(gaps)
tau_cal, av_cal = [], []
for taustar in [0.06, 0.12, 0.20, 0.30, 0.40, 0.50]:
    om = taustar * L0
    U = 0.02
    zp = z + U * np.cos(om * z)
    gp = np.diff(zp); mp_ = 0.5 * (zp[:-1] + zp[1:])
    yg = np.log(gp * np.log(mp_ / TWO_PI) / TWO_PI)
    X = np.vstack([np.ones_like(mp_), np.cos(om * mp_), np.sin(om * mp_)]).T
    b, *_ = np.linalg.lstsq(X, yg, rcond=None)
    yg0 = np.log(gaps * np.log(tmid / TWO_PI) / TWO_PI)
    X0 = np.vstack([np.ones_like(tmid), np.cos(om * tmid), np.sin(om * tmid)]).T
    b0, *_ = np.linalg.lstsq(X0, yg0, rcond=None)
    tau_cal.append(taustar)
    av_cal.append(np.hypot(b[1]-b0[1], b[2]-b0[2]) / (U * om))
D_emp = vv / (2 * tauv * np.interp(tauv, tau_cal, av_cal))

N, M = 256, 1200
ms = [8, 13, 26, 51, 77, 102, 115, 128]
rho = np.zeros((M, len(ms)), dtype=complex)
Gm = np.zeros((M, len(ms)), dtype=complex)
for s in range(M):
    A = (rng.normal(size=(N, N)) + 1j * rng.normal(size=(N, N))) / np.sqrt(2)
    Q, R = np.linalg.qr(A)
    Q = Q * (np.diagonal(R) / np.abs(np.diagonal(R)))
    th = np.sort(np.angle(np.linalg.eigvals(Q)))
    dth = np.diff(np.concatenate([th, [th[0] + TWO_PI]]))
    mid = th + dth / 2
    ds = dth * N / TWO_PI - 1
    for j, m_ in enumerate(ms):
        rho[s, j] = np.exp(1j * m_ * th).sum()
        Gm[s, j] = (ds * np.exp(1j * m_ * mid)).sum()
D_cue = []
for j, m_ in enumerate(ms):
    kap = TWO_PI * m_ / N
    D_cue.append((m_ / N, abs(np.mean(Gm[:, j] * np.conj(rho[:, j])))
                  / ((2 * np.sin(kap/2) / kap) * np.mean(np.abs(rho[:, j])**2))))
D_cue = np.array(D_cue)

lines_all = lines + [np.log(q) for q in [67,71,73,79,81,83,89,97,101,103,107,109,113,121,125,127,128]]
D_sp = []
for tau0 in [0.04, 0.06, 0.08, 0.10, 0.125, 0.15, 0.20, 0.30, 0.40, 0.50]:
    num = 0.0 + 0j; den = 0.0
    for k in K41[:6]:
        gz, tm = d41[f"gaps_{k}"], d41[f"tmid_{k}"]
        Lz = float(np.log(tm / TWO_PI).mean())
        zz = np.empty(len(gz) + 1)
        zz[0] = tm[0] - gz[0] / 2
        zz[1:] = zz[0] + np.cumsum(gz)
        dsz = gz * np.log(tm / TWO_PI) / TWO_PI - 1
        oms = tau0 * Lz + np.linspace(-0.02 * Lz, 0.02 * Lz, 160)
        oms = np.array([o for o in oms
                        if min(abs(o - l) for l in lines_all) > 0.01])
        for s0 in range(0, len(oms), 40):
            ob = oms[s0:s0+40]
            rr = np.exp(1j * np.outer(ob, zz)).sum(axis=1)
            GG = (np.exp(1j * np.outer(ob, tm)) * dsz[None, :]).sum(axis=1)
            num += (GG * np.conj(rr)).sum()
            den += (np.abs(rr)**2).sum()
    kap = 2 * np.pi * tau0
    D_sp.append((tau0, abs(num) / ((2 * np.sin(kap/2) / kap) * den)))
D_sp = np.array(D_sp)

fig, ax = plt.subplots(figsize=(8.8, 5.8))
ax.plot(tauv, D_emp, "o", ms=4, c="steelblue", alpha=0.5,
        label="prime branch: $D = v/(2\\tau a_v)$ (132 points)")
ax.plot(D_cue[:, 0], D_cue[:, 1], "s-", c="firebrick", ms=6, lw=1.3,
        label="CUE conditional response")
ax.plot(D_sp[:, 0], D_sp[:, 1], "^-", c="#2E8B57", ms=7, lw=1.3,
        label="spontaneous branch (off-line bands)")
ax.axhline(1, color="gray", lw=0.8, ls=":", label="additive-core null ($D=1$)")
ax.axvline(np.log(2)/10.4, color="gray", lw=0.8, ls="--",
           label=r"first prime line $\tau_2 = \log 2/L$")
ax.set_xlabel(r"$\tau = \kappa/2\pi$")
ax.set_ylabel("$D$ — conditional gap response")
ax.set_title("The two-branch response of the zero gas")
ax.set_ylim(0, 1.28)
ax.legend(fontsize=8.5); ax.grid(alpha=0.3)
plt.tight_layout(); plt.savefig(HERE / "fig_twobranch_en.png", dpi=110); plt.close()
print("fig_twobranch_en ok")
