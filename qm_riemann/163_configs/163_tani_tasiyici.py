"""163 probe5 — TASIYICI TANISI: s2'yi X~'nin CIZGILERI mi, artigi mi tasiyor?"""
import importlib, sys, time
from pathlib import Path
import numpy as np
sys.path.insert(0, "/Users/ugursezen/Desktop/arin/deney/qm_riemann/163_configs")
K = importlib.import_module("163_cekirdek")
veri = sys.argv[1]; TU = 0.95
Y = K.gaz(veri, 0.40); L=Y.L
d = np.load(f"sp_{veri}_{TU}.npz"); W=d["W"]; h1=d["h"]; x1=d["x"]; TAU=d["tau"]
m = Y.m0; N=len(m)
def sentez(amp, W, blok=20000, fblok=400):
    out = np.zeros(N)
    for f0 in range(0,len(W),fblok):
        fs=slice(f0,min(f0+fblok,len(W))); Wc=W[fs]; A=amp[fs]
        for s in range(0,N,blok):
            sl=slice(s,min(s+blok,N))
            arg=np.outer(m[sl],Wc)
            out[sl]+= np.cos(arg)@A.real - np.sin(arg)@A.imag
            del arg
    return out
t0=time.time()
Xl = sentez(x1, W); Xr = Y.Xtil0 - Xl
El = sentez(h1, W); Er = Y.e0 - El
print(f"{veri}: sentez {time.time()-t0:.0f}s  Var(X~0)={np.var(Y.Xtil0):.5f} "
      f"Var(Xlin)={np.var(Xl):.5f} Var(Xart)={np.var(Xr):.5f} kov={2*np.mean(Xl*Xr):.5f}", flush=True)
print(f"   Var(eta)={np.var(Y.e0):.5f} Var(eta_lin)={np.var(El):.5f} Var(eta_art)={np.var(Er):.5f}", flush=True)
print(f"{'q':>7} {'tau':>7} | {'s2':>10} {'ll':>10} {'lr':>10} {'rr':>10} | {'s1':>10} {'l':>10} {'r':>10} | {'u1':>10} {'l':>10} {'r':>10}")
for tq in (0.46,0.50,0.54,0.58,0.62,0.66,0.70,0.78):
    i=int(np.argmin(np.abs(TAU-tq))); Wf=W[i]
    cw,sw=np.cos(Wf*m),np.sin(Wf*m)
    zc=complex(2*np.mean(Y.e0*cw),-2*np.mean(Y.e0*sw))
    cw1,sw1=np.cos(Wf*Y.mid[1:]),np.sin(Wf*Y.mid[1:])
    c1r,c1i=Y.e1*cw1,-Y.e1*sw1; zb=zc/2
    rho=c1r*zb.real+c1i*zb.imag; sig=c1i*zb.real-c1r*zb.imag; rm=rho.mean()
    s2=np.mean(sig*Y.Xtil0**2)/rm
    ll=np.mean(sig*Xl*Xl)/rm; lr=2*np.mean(sig*Xl*Xr)/rm; rr=np.mean(sig*Xr*Xr)/rm
    s1=np.mean(sig*Y.Xtil0)/rm; s1l=np.mean(sig*Xl)/rm; s1r=np.mean(sig*Xr)/rm
    u1=np.mean(rho*Y.Xtil0)/rm; u1l=np.mean(rho*Xl)/rm; u1r=np.mean(rho*Xr)/rm
    print(f"{int(d['q'][i]):>7d} {TAU[i]:.4f} | {s2:+10.5f} {ll:+10.5f} {lr:+10.5f} {rr:+10.5f} | "
          f"{s1:+10.5f} {s1l:+10.5f} {s1r:+10.5f} | {u1:+10.5f} {u1l:+10.5f} {u1r:+10.5f}", flush=True)
