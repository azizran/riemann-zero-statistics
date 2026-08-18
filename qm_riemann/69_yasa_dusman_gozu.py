"""
69 — YASANIN DÜŞMAN-GÖZ TURU (18 Ağustos 2026, gece nefesi)
=============================================================

Aday: √w + b·v = 1 (b ≈ 1.07). Üç saldırı:

S1 — ORTAK-MOD TESTİ (yarı-bölme): her pencerede boşluklar çift/tek
     indekslere ayrılır; w yarı-A'dan, v yarı-B'den (ve tersi) ölçülür →
     gürültüler bağımsız. Çapraz-yarı RMS ≈ aynı-yarı RMS ise ortak-mod
     hata açıklaması ÖLÜR.

S2 — ÜS ÖLÇÜMÜ: w = (1 − b·v)^k, k SERBEST. k = 2'ye ne kadar yakın?
     (bootstrap hatasıyla)

S3 — FORM YARIŞI (eşit parametre): √w-lineer vs w-lineer vs
     şiddet-ünitaritesi (w + c·v² = 1) vs üstel (w = e^{−cv}).
"""

import numpy as np
import mpmath as mp
from pathlib import Path
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
A = (np.e**2 - 5) / 2
B0, B1 = 2.7580, -0.0543
T0_ODL = 267653395647
PRIMES6 = [2, 3, 5, 7, 11, 13]
mp.mp.dps = 30

def channels(gaps, amps, tmid, anchored, idx=None):
    if idx is not None:
        gaps, amps, tmid = gaps[idx], amps[idx], tmid[idx]
    Lw = np.log(((T0_ODL + tmid) if anchored else tmid) / TWO_PI)
    g_u = gaps * Lw / TWO_PI
    a_u = amps / np.sqrt(A * Lw + B0 + B1 / Lw)
    a_u /= np.sqrt((a_u**2).mean())
    ya, yg = np.log(a_u), np.log(g_u)
    cw = [np.ones_like(ya), g_u, g_u**2]
    cv = [np.ones_like(yg)]
    for p in PRIMES6:
        ph0 = float(mp.fmod(T0_ODL * mp.log(p), 2 * mp.pi)) if anchored else 0.0
        arg = ph0 + tmid * np.log(p)
        cw += [np.cos(arg), np.sin(arg)]
        cv += [np.cos(arg), np.sin(arg)]
    Xw = np.vstack(cw).T
    bw, *_ = np.linalg.lstsq(Xw, ya, rcond=None)
    Xv = np.vstack(cv).T
    bv, *_ = np.linalg.lstsq(Xv, yg, rcond=None)
    ws = [bw[3 + 2*i] / p**-0.5 for i, p in enumerate(PRIMES6)]
    vs = [float(np.hypot(bv[1 + 2*i], bv[2 + 2*i]) / p**-0.5) for i, p in enumerate(PRIMES6)]
    return ws, vs

WINDOWS = []
d41 = np.load(HERE / "41_bigT_windows.npz")
for k in sorted({x.split("_")[1] for x in d41.files}, key=lambda s: int(s[:-1])):
    WINDOWS.append((d41[f"gaps_{k}"], d41[f"amps_{k}"], d41[f"tmid_{k}"], False))
for f in ["55_win_1e+08.npz", "55_win_1e+09.npz", "55_win_1e+10.npz", "55_win_1e+11.npz"]:
    d = np.load(HERE / f)
    WINDOWS.append((d["gaps"], d["amps"], d["tmid"], False))
d53 = np.load(HERE / "53_odlyzko_amps.npz")
WINDOWS.append((d53["gaps"], d53["max_amps"], d53["t_mid"], True))

# ---- S1: yarı-bölme ----
same, cross = [], []
for gaps, amps, tmid, anch in WINDOWS:
    n = len(gaps)
    iA, iB = np.arange(0, n, 2), np.arange(1, n, 2)
    wA, vA = channels(gaps, amps, tmid, anch, iA)
    wB, vB = channels(gaps, amps, tmid, anch, iB)
    for i in range(len(PRIMES6)):
        if wA[i] > 0 and wB[i] > 0:
            same.append((np.sqrt(wA[i]), vA[i]))
            same.append((np.sqrt(wB[i]), vB[i]))
            cross.append((np.sqrt(wA[i]), vB[i]))
            cross.append((np.sqrt(wB[i]), vA[i]))

def law_rms(pairs):
    sq = np.array([p[0] for p in pairs]); vv = np.array([p[1] for p in pairs])
    Amat = np.vstack([np.ones_like(vv), vv]).T
    c, *_ = np.linalg.lstsq(Amat, sq, rcond=None)
    return np.sqrt(np.mean((sq - Amat @ c) ** 2)), c

rms_same, c_same = law_rms(same)
rms_cross, c_cross = law_rms(cross)
print("S1 — ORTAK-MOD (yarı-bölme):")
print(f"  aynı-yarı  çiftler ({len(same)}): RMS = {rms_same:.4f}  "
      f"(a={c_same[0]:.3f}, b={-c_same[1]:.3f})")
print(f"  ÇAPRAZ-yarı çiftler ({len(cross)}): RMS = {rms_cross:.4f}  "
      f"(a={c_cross[0]:.3f}, b={-c_cross[1]:.3f})")
print("  → çapraz ≈ aynı ise ortak-mod hata açıklaması ölür\n")

# ---- tam-veri 66 çifti (68'deki ölçüm, tam pencerelerle) ----
full = []
for gaps, amps, tmid, anch in WINDOWS:
    ws, vs = channels(gaps, amps, tmid, anch)
    for i in range(len(PRIMES6)):
        if ws[i] > 0:
            full.append((ws[i], vs[i]))
wf = np.array([p[0] for p in full]); vf = np.array([p[1] for p in full])
sqf = np.sqrt(wf)

# ---- S2: üs ölçümü ----
def model_k(par):
    b, k = par
    base = np.clip(1 - b * vf, 1e-6, None)
    return base**k - wf

fit = least_squares(model_k, [1.07, 2.0], bounds=([0.3, 0.5], [3.0, 6.0]))
b_fit, k_fit = fit.x
rng = np.random.default_rng(69)
ks = []
for _ in range(300):
    ii = rng.integers(0, len(wf), len(wf))
    try:
        r = least_squares(lambda p: (np.clip(1 - p[0]*vf[ii],1e-6,None)**p[1] - wf[ii]),
                          [b_fit, k_fit], bounds=([0.3, 0.5], [3.0, 6.0]))
        ks.append(r.x[1])
    except Exception:
        pass
print(f"S2 — ÜS ÖLÇÜMÜ: w = (1 − b·v)^k")
print(f"  b = {b_fit:.4f},  k = {k_fit:.3f} ± {np.std(ks):.3f}  (bootstrap)")
print(f"  → k = 2 (genlik-korunumu) {abs(k_fit-2)/np.std(ks):.1f}σ uzakta\n")

# ---- S3: form yarışı ----
def rms_of(resid):
    return np.sqrt(np.mean(resid**2))

Am = np.vstack([np.ones_like(vf), vf]).T
c1, *_ = np.linalg.lstsq(Am, sqf, rcond=None)
c2, *_ = np.linalg.lstsq(Am, wf, rcond=None)
cq = np.sum((1 - wf) * vf**2) / np.sum(vf**4)
ce = np.sum(-np.log(wf) * vf) / np.sum(vf**2)
print("S3 — FORM YARIŞI (w-ölçeğinde artık RMS):")
print(f"  √w = a − b·v (2p)      : {rms_of((Am @ c1)**2 - wf):.4f}")
print(f"  w  = a − b·v (2p)      : {rms_of(Am @ c2 - wf):.4f}")
print(f"  w + c·v² = 1 (1p)      : {rms_of((1 - cq*vf**2) - wf):.4f}")
print(f"  w = exp(−c·v) (1p)     : {rms_of(np.exp(-ce*vf) - wf):.4f}")
print(f"  w = (1−b·v)^k (2p)     : {rms_of(np.clip(1-b_fit*vf,1e-6,None)**k_fit - wf):.4f}")
