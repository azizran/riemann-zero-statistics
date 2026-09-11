# -*- coding: utf-8 -*-
"""
GAP 2: 88_S_modeli.py M2 bloğu BİREBİR (yazarın tanımı) — tarak parlaklığı.
  x = rvm_N(tmid); P1 = exp(2πi x); delta = angle(P1·conj(mean yönü))/2π
  sig = delta.std();  meas = |<exp(2πi k x)>|;  gauss = exp(−(2πk)² sig²/2)
"""
import numpy as np
TWO_PI = 2*np.pi
def rvm_N(t):
    x = t/TWO_PI
    return x*np.log(x/np.e) + 7/8

d = np.load("../qm_riemann/41_bigT_windows.npz")
K = sorted({f.split("_")[1] for f in d.files}, key=lambda s: int(s[:-1]))
WNDS = [(d[f"gaps_{k}"], d[f"amps_{k}"], d[f"tmid_{k}"]) for k in K[:6]]
print(f"{'L':>7} {'k':>2} {'ölçüm':>8} {'Gauss':>8} {'oran':>7} {'σ_δ':>7} {'kurt':>7}")
oranlar=[]
for gaps, amps, tmid in WNDS[:6]:
    L = float(np.log(tmid/TWO_PI).mean())
    x = rvm_N(tmid)
    P1 = np.exp(2j*np.pi*x)
    delta = np.angle(P1*np.conj(P1.mean()/abs(P1.mean())))/TWO_PI
    sig = delta.std()
    k4 = float(((delta-delta.mean())**4).mean()/sig**4 - 3)
    for k in [1,2]:
        meas = abs(np.exp(2j*np.pi*k*x).mean())
        kk = TWO_PI*k
        gauss = np.exp(-kk**2*sig**2/2)
        edge = gauss*(1 + k4*sig**4*kk**4/24)
        o = meas/gauss
        if k==1: oranlar.append(o)
        print(f"{L:7.3f} {k:2d} {meas:8.4f} {gauss:8.4f} {o:7.4f} {sig:7.4f} {k4:+7.2f}")
print(f"\nk=1 oran ortalaması = {np.mean(oranlar):.4f} ± {np.std(oranlar):.4f}   (makale: 0.88–0.89 sabit)")

# ---------- İKİNCİ TABAN: unfold edilmiş jitter (sigma_t → σ_u) ----------
def sigma_t(gaps, tmid):
    t0 = tmid[0] - gaps[0]/2; kk = np.arange(len(tmid)) + 0.5
    ts = t0 + kk*TWO_PI/np.log(t0/TWO_PI)
    for _ in range(6):
        ts = ts - (rvm_N(ts) - rvm_N(t0) - kk)/(np.log(ts/TWO_PI)/TWO_PI)
    return float((tmid - ts).std()), ts

print("\n=== İKİNCİ TABAN: sigma_t (mutlak) → σ_u = σ_t·L/2π  →  DW = exp(−(2πσ_u)²/2) ===")
print(f"{'L':>7} {'σ_t':>8} {'σ_u':>8} {'ölçüm':>8} {'DW':>8} {'oran':>7}")
oran2=[]
for gaps, amps, tmid in WNDS[:6]:
    L = float(np.log(tmid/TWO_PI).mean())
    x = rvm_N(tmid)
    meas = abs(np.exp(2j*np.pi*x).mean())
    sig_t,_ = sigma_t(gaps, tmid)
    sig_u = sig_t*L/TWO_PI
    dw = np.exp(-0.5*(TWO_PI*sig_u)**2)
    o = meas/dw; oran2.append(o)
    print(f"{L:7.3f} {sig_t:8.5f} {sig_u:8.5f} {meas:8.4f} {dw:8.4f} {o:7.4f}")
print(f"\nunfold-jitter tabanıyla oran = {np.mean(oran2):.4f} ± {np.std(oran2):.4f}   (makale: 0.88–0.89 sabit)")
