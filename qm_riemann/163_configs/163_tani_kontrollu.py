"""163 probe6 — KONTROLLU SINAV: alanlari KENDI cizgilerinden yeniden kur,
   uclu-toplam formulunu ORADA sina (butun genlik/fazlar bilinir)."""
import importlib, sys, time
from pathlib import Path
import numpy as np
sys.path.insert(0, "/Users/ugursezen/Desktop/arin/deney/qm_riemann/163_configs")
K = importlib.import_module("163_cekirdek")
veri = sys.argv[1]; TU=0.95; UNIF = len(sys.argv)>2 and sys.argv[2]=="unif"
Y = K.gaz(veri, 0.40); L=Y.L
d=np.load(f"sp_{veri}_{TU}.npz"); W=d["W"]; h1=d["h"]; x1=d["x"]; TAU=d["tau"]; QQ=d["q"]
m = Y.m0.copy(); N=len(m)
if UNIF:                      # DUZGUN tarak: G-katmani olmez
    m = m[0] + (2*np.pi/L)*np.arange(N)
def sentez(amp, Wv, mm, blok=20000, fblok=500):
    out=np.zeros(len(mm))
    for f0 in range(0,len(Wv),fblok):
        fs=slice(f0,min(f0+fblok,len(Wv))); Wc=Wv[fs]; A=amp[fs]
        for s in range(0,len(mm),blok):
            sl=slice(s,min(s+blok,len(mm)))
            arg=np.outer(mm[sl],Wc)
            out[sl]+= np.cos(arg)@A.real - np.sin(arg)@A.imag
            del arg
    return out
t0=time.time()
Et = sentez(h1, W, m); Xt = sentez(x1, W, m)
Xt = Xt - Xt.mean()
print(f"{veri}{' UNIF' if UNIF else ''}: sentez {time.time()-t0:.0f}s Var(eta*)={np.var(Et):.5f} Var(X*)={np.var(Xt):.5f}", flush=True)
# yapay gazin kendi Phi tablosu (X* marjinali)
class Yb: pass
Yb.Xtil = Xt; Yb.Xtil0 = Xt - Xt.mean()
PT = K.PhiTablo(Yb)
Et1 = np.concatenate((Et[1:], Et[:1]))     # eta*_{n+1}
mid1 = np.concatenate((m[1:], m[:1]+ (m[-1]-m[-2])))
print(f"{'q':>7} {'tau':>7} | {'s2*':>10} {'s2*_pred':>10} {'oran':>7} | {'u1*':>10} {'u1*p':>10} | {'u0*p':>7}")
for tq in (0.46,0.50,0.54,0.58,0.62,0.66,0.70,0.78):
    i=int(np.argmin(np.abs(TAU-tq))); Wf=W[i]
    cw,sw=np.cos(Wf*m),np.sin(Wf*m)
    zc=complex(2*np.mean(Et*cw),-2*np.mean(Et*sw))
    cw1,sw1=np.cos(Wf*mid1),np.sin(Wf*mid1)
    c1r,c1i=Et1*cw1,-Et1*sw1; zb=zc/2
    rho=c1r*zb.real+c1i*zb.imag; sig=c1i*zb.real-c1r*zb.imag; rm=rho.mean()
    X0=Yb.Xtil0
    s2=np.mean(sig*X0**2)/rm; u1=np.mean(rho*X0)/rm; 
    xq=complex(2*np.mean(X0*cw),-2*np.mean(X0*sw))
    msk=np.abs(TAU-TAU[i])>1e-12
    pr=K.ongor_cizgi(zc, xq, TAU[i], TAU[msk], h1[msk], x1[msk], PT, 2)
    o = s2/pr['s2_pred'] if pr['s2_pred'] else float('nan')
    print(f"{int(QQ[i]):>7d} {TAU[i]:.4f} | {s2:+10.5f} {pr['s2_pred']:+10.5f} {o:7.2f} | "
          f"{u1:+10.5f} {pr['u1_pred']:+10.5f} | {pr['u0_pred']:7.4f}", flush=True)
