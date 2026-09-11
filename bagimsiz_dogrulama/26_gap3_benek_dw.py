# -*- coding: utf-8 -*-
"""
GAP 3: benek yasasının TAM formu (DW dahil) — makale: pred = τ·q^{-1/2}·cos(πτ)·DW
DW = exp(−(log q)² σ²/2),  σ = std(orta nokta − pürüzsüz kafes)   [89/103 ile aynı]
Benim önceki testim DW'siz idi → o yüzden oran 0.99→0.79 diye düşüyordu.
"""
import numpy as np, json
TWO_PI=2*np.pi
zeros=np.sort(np.asarray(np.load("../qm_riemann/128_odl_zeros6_2e6_zeros.npz")["zeros"],float))
i0=int(np.argmin(np.abs(zeros-1.2e5))); g=zeros[i0-20000:i0+20001]
mid=0.5*(g[:-1]+g[1:]); L=float(np.log(np.mean(g)/TWO_PI))

def rvm_N(t):
    x=t/TWO_PI
    return x*np.log(x/np.e)+7/8

def sigma_t(gaps, tmid):
    t0=tmid[0]-gaps[0]/2; kk=np.arange(len(tmid))+0.5
    ts=t0+kk*TWO_PI/np.log(t0/TWO_PI)
    for _ in range(6):
        ts = ts - (rvm_N(ts)-rvm_N(t0)-kk)/(np.log(ts/TWO_PI)/TWO_PI)
    return float((tmid-ts).std()), ts

sig,_=sigma_t(np.diff(g), mid)
print(f"pencere: L={L:.5f}  σ_t (mutlak birim) = {sig:.5f}   → σ/ortalama_aralık = {sig/np.mean(np.diff(g)):.3f}")

def asal_kuv(lim):
    ps=[p for p in range(2,lim+1) if all(p%d for d in range(2,int(p**.5)+1))]
    o=[]
    for p in ps:
        q=p;k=1
        while q<=lim: o.append((q,p,k)); q*=p;k+=1
    return sorted(o)

print(f"\n{'q':>5} {'τ':>6} {'ölçülen':>9} {'pred(DW)':>9} {'pred(DWsiz)':>11} {'oran(DW)':>9} {'oran(DWsiz)':>11}")
sat=[]
for (q,p,k) in asal_kuv(84):
    tau=np.log(q)/L
    Gm=abs(np.exp(1j*np.log(q)*mid).mean())
    dw=np.exp(-np.log(q)**2*sig**2/2)
    pred_dw  = tau*q**-0.5*np.cos(np.pi*tau)*dw
    pred_yok = tau*q**-0.5*np.cos(np.pi*tau)
    sat.append((q,k,Gm,pred_dw,pred_yok))
    print(f"{q:5d} {tau:6.3f} {Gm:9.5f} {pred_dw:9.5f} {pred_yok:11.5f} {Gm/pred_dw:9.3f} {Gm/pred_yok:11.3f}")
a=np.array([s[2]/s[3] for s in sat]); b=np.array([s[2]/s[4] for s in sat])
print(f"\nDW DAHİL   : oran ort = {a.mean():.4f} ± {a.std():.4f}   (min {a.min():.3f}, maks {a.max():.3f})")
print(f"DW'siz     : oran ort = {b.mean():.4f} ± {b.std():.4f}   (min {b.min():.3f}, maks {b.max():.3f})")
asal=[s for s in sat if s[0]==s[0] and s[1]==1]
if asal:
    x=np.array([s[2]/s[3] for s in asal]); print(f"yalnız asallar (DW dahil): {x.mean():.4f} ± {x.std():.4f}, saçılma %{100*x.std()/x.mean():.1f}")
json.dump(dict(sigma_t=sig, L=L, satir=[(int(q),int(k),float(m),float(pd),float(py)) for q,k,m,pd,py in sat]),
          open("26_gap3.json","w"), indent=1)
