"""
77 — YASA MANŞET FİTİ (Not 3 §2'nin sayıları, tek script) (19 Ağustos)
========================================================================
Denetim S5: 44-çift fitinin numaralı scripti yoktu. Bu script Not 3'teki
manşet sayıları üretir: serbest fit, a≡1 fiti, pencere-blok bootstrap
eğim hatası, serbest-üs (blok bootstrap) ve form yarışı.
"""

import numpy as np
import mpmath as mp
from pathlib import Path
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
A_CG = (np.e**2 - 5) / 2
B0, B1 = 2.7580, -0.0543
T0_ODL = 267653395647
mp.mp.dps = 30

def channels(gaps, amps, tmid, primes, anchored=False):
    Lw = np.log(((T0_ODL + tmid) if anchored else tmid) / TWO_PI)
    g_u = gaps * Lw / TWO_PI
    a_u = amps / np.sqrt(A_CG * Lw + B0 + B1 / Lw)
    a_u /= np.sqrt((a_u**2).mean())
    ya, yg = np.log(a_u), np.log(g_u)
    cw = [np.ones_like(ya), g_u, g_u**2]
    cv = [np.ones_like(yg)]
    for p in primes:
        ph0 = float(mp.fmod(T0_ODL * mp.log(p), 2 * mp.pi)) if anchored else 0.0
        arg = ph0 + tmid * np.log(p)
        cw += [np.cos(arg), np.sin(arg)]
        cv += [np.cos(arg), np.sin(arg)]
    Xw = np.vstack(cw).T
    bw, *_ = np.linalg.lstsq(Xw, ya, rcond=None)
    Xv = np.vstack(cv).T
    bv, *_ = np.linalg.lstsq(Xv, yg, rcond=None)
    return [(bw[3 + 2*i] / p**-0.5,
             float(np.hypot(bv[1 + 2*i], bv[2 + 2*i]) / p**-0.5)) for i, p in enumerate(primes)]

P4 = [2, 3, 5, 7]
WNDS = []
d41 = np.load(HERE / "41_bigT_windows.npz")
for k in sorted({x.split("_")[1] for x in d41.files}, key=lambda s: int(s[:-1])):
    WNDS.append((d41[f"gaps_{k}"], d41[f"amps_{k}"], d41[f"tmid_{k}"], False))
for f in ["55_win_1e+08.npz", "55_win_1e+09.npz", "55_win_1e+10.npz", "55_win_1e+11.npz"]:
    d = np.load(HERE / f)
    WNDS.append((d["gaps"], d["amps"], d["tmid"], False))
d53 = np.load(HERE / "53_odlyzko_amps.npz")
WNDS.append((d53["gaps"], d53["max_amps"], d53["t_mid"], True))

pairs, widx = [], []
for wi, (gaps, amps, tmid, anch) in enumerate(WNDS):
    for w_, v_ in channels(gaps, amps, tmid, P4, anch):
        if w_ > 0:
            pairs.append((w_, v_)); widx.append(wi)
pairs = np.array(pairs); widx = np.array(widx)
wv, vv = pairs.T
sq = np.sqrt(wv)
print(f"{len(pairs)} eşleşmiş çift (4 asal × 11 pencere)")

X = np.vstack([np.ones_like(vv), vv]).T
c, *_ = np.linalg.lstsq(X, sq, rcond=None)
print(f"SERBEST: √w = {c[0]:.4f} − {-c[1]:.4f}·v   RMS = "
      f"{np.sqrt(np.mean((sq - X@c)**2)):.4f}")
b1 = np.sum((1 - sq) * vv) / np.sum(vv * vv)
print(f"a≡1   : b = {b1:.4f}   RMS = {np.sqrt(np.mean((sq - 1 + b1*vv)**2)):.4f}")

rng = np.random.default_rng(77)
slopes, ks = [], []
uw = np.unique(widx)
for _ in range(400):
    sel = rng.choice(uw, len(uw), replace=True)
    ii = np.concatenate([np.where(widx == s)[0] for s in sel])
    cc, *_ = np.linalg.lstsq(X[ii], sq[ii], rcond=None)
    slopes.append(-cc[1])
    try:
        r = least_squares(lambda p: np.clip(1 - p[0]*vv[ii], 1e-6, None)**p[1] - wv[ii],
                          [1.07, 2.0], bounds=([0.3, 0.5], [3.0, 6.0]))
        ks.append(r.x[1])
    except Exception:
        pass
ks = np.array(ks)
print(f"pencere-blok bootstrap: eğim = {np.mean(slopes):.4f} ± {np.std(slopes):.4f}")

# serbest üs: ortak (b,k) fiti dejenere bir vadide gezer — dürüst rapor:
rj = least_squares(lambda p: np.clip(1 - p[0]*vv, 1e-6, None)**p[1] - wv,
                   [1.07, 2.0], bounds=([0.3, 0.5], [3.0, 6.0]))
rk = least_squares(lambda p: np.clip(1 - b1*vv, 1e-6, None)**p[0] - wv, [2.0])
print(f"serbest üs — ortak fit: b={rj.x[0]:.3f}, k={rj.x[1]:.3f} (dejenere vadi; "
      f"blok-bootstrap medyan {np.median(ks):.2f}, %16-84 "
      f"[{np.percentile(ks,16):.2f},{np.percentile(ks,84):.2f}])")
print(f"serbest üs — koşullu (b={b1:.3f} doğrusal yasadan): k = {rk.x[0]:.4f}")

cq = np.sum((1 - wv) * vv**2) / np.sum(vv**4)
c2, *_ = np.linalg.lstsq(X, wv, rcond=None)
r_sqrt = np.sqrt(np.mean(((X @ c)**2 - wv)**2))
r_lin = np.sqrt(np.mean((X @ c2 - wv)**2))
r_uni = np.sqrt(np.mean(((1 - cq*vv**2) - wv)**2))
print(f"form yarışı (w-ölçeği): √w-lineer {r_sqrt:.4f} | w-lineer {r_lin:.4f} | "
      f"şiddet-ünitaritesi {r_uni:.4f} (×{r_uni/r_sqrt:.1f})")
