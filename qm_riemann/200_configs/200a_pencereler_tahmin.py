# KALEM 200-A: windows (zero positions only; no Z, no M) + parameter-free model predictions for rungs b=0,1 (and b=2 for reference)
import numpy as np, mpmath as mp, importlib, json
from sympy import isprime
_k=importlib.import_module('200t_merdiven').kap
kap=lambda r,N,b: float(_k(r,N,b))
z=np.load('../128_odl_zeros6_2e6_zeros.npz')['zeros']
Lz=np.log(z/(2*np.pi))
W={'W1':(8.0,10.0),'W2':(10.0,11.0),'W3':(11.0,float(Lz.max()))}
def S12(x):
    s=mp.mpf(0); H=mp.mpf(0); m=1; t=x
    while True:
        m+=1; H+=mp.mpf(1)/(m-1); t*=x; term=H*t/m**2; s+=term
        if term<mp.mpf(10)**-25: return s
def euler(X):
    ps=[p for p in range(2,X+1) if isprime(p)]
    return float(mp.fsum(mp.polylog(2,mp.mpf(1)/p)/2 for p in ps)), float(mp.fsum(1.5*S12(mp.mpf(1)/p) for p in ps))
A2,A3=-0.088124,0.233653
eg=float(mp.e**mp.euler)
out={}
for w,(a,b) in W.items():
    i0,i1=np.searchsorted(Lz,a),np.searchsorted(Lz,b,side='right')
    ta,tb=float(z[i0]),float(z[i1-1])
    # rung 0: t uniform on [ta,tb] -> mean L = mean of log(t/2pi) over uniform t
    L0=float(np.mean(np.log(np.linspace(ta,tb,200001)/(2*np.pi))))
    L1=float(Lz[i0:i1].mean())                    # rung 1: zero-weighted
    rec={'idx':[int(i0),int(i1)],'t':[ta,tb],'n_zeros':int(i1-i0),'L_rung0':L0,'L_rung1':L1,'pred':{}}
    for bb,L in ((0,L0),(1,L1),(2,L1)):
        P={'CUE':(kap(2,L,bb),kap(3,L,bb)),'a_k':(kap(2,L,bb)+A2,kap(3,L,bb)+A3)}
        for X in (2,3,5,7,11):
            NX=L/(eg*np.log(X))
            if NX-bb>0.3:
                e2,e3=euler(X); P[f'hyb{X}']=(kap(2,NX,bb)+e2,kap(3,NX,bb)+e3)
        gs={0:(np.pi**2/8,-14*1.2020569031595942/8),1:(np.pi**2/24,-2*1.2020569031595942/8),2:(float(mp.psi(1,1.5))/4,float(mp.psi(2,1.5))/8)}
        P['Gauss']=gs[bb]
        rec['pred'][f'b{bb}']={k:[round(v[0],4),round(v[1],4)] for k,v in P.items()}
    out[w]=rec
json.dump(out,open('200a_tahmin.json','w'),indent=1)
for w,r in out.items():
    print(w,'zeros',r['idx'],'n',r['n_zeros'],'t [%.0f, %.0f]'%tuple(r['t']),'L0 %.3f L1 %.3f'%(r['L_rung0'],r['L_rung1']))
    for bb in ('b0','b1','b2'):
        print('  ',bb,' '.join(f'{k}:{v[0]:.3f}/{v[1]:+.3f}' for k,v in r['pred'][bb].items()))
