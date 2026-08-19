"""
95 — NOT 4 DENETİM EKLERİ (21 Ağustos 2026, gece)
==========================================================================
Denetçinin işaret ettiği üç eksik/kaba ölçümün numaralı kaydı:
  T1  KUVVET-FAZI KILL-TESTİ (S4): gerçek kuvvet-beneği fazları +
      istatistiksel faz belirsizlikleri (91'in öldürme testi kod olarak).
  T2  DAL KESİŞMESİ (K1): asal dalı D_emp (83 + 6-nokta a_v) × spontane
      dal (93 makinesi, ince bantlar) → τ* kesişmesi. (92'nin kaba
      τ≈0.17'si eski-ızgara artefaktıydı.)
  T3  ENJEKSİYON DOĞRUSALLIĞI (S9): U = 0.02 vs 0.10, τ* = 0.25.
"""

import numpy as np
from pathlib import Path

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi

def Ghat(t, omegas, chunk=30000):
    out = np.zeros(len(omegas), dtype=complex)
    for s0 in range(0, len(t), chunk):
        tt = t[s0:s0 + chunk]
        out += np.exp(1j * np.outer(omegas, tt)).sum(axis=1)
    return out / len(t)

d41 = np.load(HERE / "41_bigT_windows.npz")
K41 = sorted({x.split("_")[1] for x in d41.files}, key=lambda s: int(s[:-1]))
WNDS = [(d41[f"gaps_{k}"], d41[f"amps_{k}"], d41[f"tmid_{k}"]) for k in K41[:6]]

# ============ T1: kuvvet-fazı kill-testi ============
print("T1 — KUVVET-BENEĞİ FAZLARI (öldürme testi; lab öngörüsü 10-30° dönme):")
print(f"{'q':>3} {'faz ort':>9} {'sapma':>7} {'ist.belirsizlik':>15}")
for q in [5, 13, 4, 8, 9, 16, 25, 27, 32, 49]:
    phs, mags = [], []
    for gaps, amps, tmid in WNDS:
        G = Ghat(tmid, np.array([np.log(q)]))[0]
        ph = np.degrees(np.angle(G))
        phs.append(ph if ph > 0 else ph + 360)
        mags.append(abs(G))
        n = len(tmid)
    phs = np.array(phs); mags = np.array(mags)
    sig_ph = np.degrees((1/np.sqrt(2*n)) / mags.mean()) / np.sqrt(6)
    print(f"{q:>3} {phs.mean():>8.2f}° {abs(phs.mean()-180):>6.2f}° {sig_ph:>13.1f}°")

# ============ T2: dal kesişmesi ============
print("\nT2 — DAL KESİŞMESİ:")
d83 = np.load(HERE / "83_tam_taban_egri.npz")
tauv, vv = d83["tau"], d83["v"]
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
srt = np.argsort(tauv)
tt_s, DD_s = tauv[srt], D_emp[srt]
grid = np.linspace(0.08, 0.25, 60)
Dp_g = np.array([np.median(DD_s[np.abs(tt_s - g) < 0.035]) for g in grid])

lines_all = [np.log(q) for q in
    [2,3,4,5,7,8,9,11,13,16,17,19,23,25,27,29,31,32,37,41,43,47,49,53,59,
     61,64,67,71,73,79,81,83,89,97,101,103,107,109,113,121,125,127,128]]
sp = []
for tau0 in [0.10, 0.125, 0.15, 0.175, 0.20]:
    num = 0.0 + 0j; den = 0.0
    for k in K41[:6]:
        gz, tm = d41[f"gaps_{k}"], d41[f"tmid_{k}"]
        Lz = float(np.log(tm / TWO_PI).mean())
        zz = np.empty(len(gz) + 1)
        zz[0] = tm[0] - gz[0] / 2
        zz[1:] = zz[0] + np.cumsum(gz)
        ds_loc = gz * np.log(tm / TWO_PI) / TWO_PI - 1
        tt2 = (tm - tm.mean()) / (tm[-1] - tm[0])
        P = np.vstack([np.ones_like(tt2), tt2, tt2**2, tt2**3]).T
        dsz = ds_loc - P @ np.linalg.lstsq(P, ds_loc, rcond=None)[0]
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
    sp.append((tau0, abs(num) / ((2 * np.sin(kap/2) / kap) * den)))
sp = np.array(sp)
Dsp_g = np.interp(grid, sp[:, 0], sp[:, 1])
fark = Dp_g - Dsp_g
ix = np.argmax(fark < 0)
tau_star = grid[ix-1] + (grid[ix]-grid[ix-1]) * fark[ix-1]/(fark[ix-1]-fark[ix])
print(f"  spontane dal: " + "  ".join(f"D({t:.3f})={d:.3f}" for t, d in sp))
print(f"  KESİŞME: τ* = {tau_star:.3f}  (aralık taraması 0.13-0.15; "
      f"çapraz-estimatör sistematiği dahil değil)")

# ============ T3: enjeksiyon doğrusallığı ============
print("\nT3 — ENJEKSİYON DOĞRUSALLIĞI (τ*=0.25):")
om = 0.25 * L0
for U in [0.02, 0.10]:
    zp = z + U * np.cos(om * z)
    gp = np.diff(zp); mp_ = 0.5 * (zp[:-1] + zp[1:])
    yg = np.log(gp * np.log(mp_ / TWO_PI) / TWO_PI)
    X = np.vstack([np.ones_like(mp_), np.cos(om * mp_), np.sin(om * mp_)]).T
    b, *_ = np.linalg.lstsq(X, yg, rcond=None)
    yg0 = np.log(gaps * np.log(tmid / TWO_PI) / TWO_PI)
    X0 = np.vstack([np.ones_like(tmid), np.cos(om * tmid), np.sin(om * tmid)]).T
    b0, *_ = np.linalg.lstsq(X0, yg0, rcond=None)
    a_v = np.hypot(b[1]-b0[1], b[2]-b0[2]) / (U * om)
    dG = abs(Ghat(mp_, np.array([om]))[0] - Ghat(tmid, np.array([om]))[0])
    a_G = dG / (om * U / 2)
    print(f"  U = {U:.2f}: a_v = {a_v:.3f}  a_G = {a_G:.3f}")
