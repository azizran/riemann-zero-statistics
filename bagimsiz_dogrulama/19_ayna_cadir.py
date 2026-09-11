# -*- coding: utf-8 -*-
"""
NOT 3 (iv): |Z|^2 modulasyonu ve AYNA (cadir) kanali — parametresiz test.
  R_pred(q) = 2 q^{-1/2} S1(N/q)/S1(N),  N=sqrt(t/2pi),  S1(x)=sum_{n<=x} 1/n
  tau>1/2  (q>N)  =>  TAM SIFIR (keskin ufuk)
  R_tam(tau) = 2(L - log q + 2g - 1)/(L + 2g - 1)   (~2(1-tau))
  AYNA = R_tam - R_pred                             (~cadir, foldda tepe)
"""
import numpy as np, importlib.util, json
def yukle(ad,yol):
    s=importlib.util.spec_from_file_location(ad,yol); m=importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
s3=yukle("s3","03_surrogate.py"); TWO_PI=2*np.pi; GAMMA=0.5772156649015329

t0,t1,h=107252.0,132748.0,0.01992
ts=t0+h*np.arange(int((t1-t0)/h))
Zg=s3.Z(ts); I=np.abs(Zg)**2
n=I.size; tc=0.5*(t0+t1); L=float(np.log(tc/TWO_PI)); N=np.sqrt(tc/TWO_PI)
S1N=sum(1.0/k for k in range(1,int(N)+1))
DC=I.mean()
print(f"izgara {n} nokta, L={L:.4f}, N=√(t/2π)={N:.2f}, DC=mean|Z|²={DC:.4f}")
print(f"keskin ufuk: τ=1/2  <=>  q=N={N:.1f}\n")
def asal_mi(x): return x>1 and all(x%d for d in range(2,int(x**0.5)+1))
print(f"{'q':>5} {'τ':>6} {'ölçülen':>9} {'R_pred':>8} {'R_tam':>8} {'ayna(ölç)':>10} {'ayna(pred)':>11}")
kayit=[]
for q in [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,67,71,79,83,97,101,103,107,113,127,131,137,139,149,151,157]:
    tau=np.log(q)/L
    c=(I*np.exp(-1j*np.log(q)*ts)).mean()
    olc=2*abs(c)/DC                                     # bağıl modülasyon genliği
    if q<=N:
        S1q=sum(1.0/k for k in range(1,int(N//q)+1)) if N//q>=1 else 0.0
        R_pred=2*q**-0.5*S1q/S1N
    else:
        R_pred=0.0
    R_tam=2*(L-np.log(q)+2*GAMMA-1)/(L+2*GAMMA-1)
    ayna_olc=olc-R_pred; ayna_pred=R_tam-R_pred
    kayit.append((q,tau,olc,R_pred,R_tam,ayna_olc,ayna_pred))
    print(f"{q:5d} {tau:6.3f} {olc:9.4f} {R_pred:8.4f} {R_tam:8.4f} {ayna_olc:10.4f} {ayna_pred:11.4f}")

alt=[k for k in kayit if k[1]<0.45]
ust=[k for k in kayit if k[1]>0.55]
if alt:
    a=np.array([k[2] for k in alt]); p=np.array([k[3] for k in alt])
    print(f"\nfold ALTINDA (τ<0.45, n={len(alt)}): ölçülen/pred oranı = {np.mean(a/p):.4f} ± {np.std(a/p):.4f}, r={np.corrcoef(a,p)[0,1]:.4f}")
if ust:
    a=np.array([k[2] for k in ust])
    print(f"fold ÜSTÜNDE (τ>0.55, n={len(ust)}): ölçülen ort = {a.mean():.4f} (pred = 0.0000, keskin ufuk)")
    print(f"  (örneklem gürültüsü mertebesi ~ {2/np.sqrt(n)/DC*0+0.002:.4f} civarı beklenir)")
json.dump(dict(L=L,n=n,N=N,kayit=kayit), open("19_ayna.json","w"), indent=1)
