# -*- coding: utf-8 -*-
"""
NOT 3 (ii): w isaret degistiriyor mu? tau0 ~ 0.447, faz 0->pi donmeden atliyor mu?
Dusuk-t pencereleri (zeros_100k, t~1.7e3..5.6e4) + yuksek-t pencereler.
ISARETLI kosinus katsayisi (a) ve kuadratar (b) ayri ayri olculur.
"""
import numpy as np, importlib.util, json
def yukle(ad,yol):
    s=importlib.util.spec_from_file_location(ad,yol); m=importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
o=yukle("o","02_olcum.py"); s7=yukle("s7","07_wv.py"); TWO_PI=2*np.pi

z100=np.sort(np.asarray([float(l) for l in open("../zeros_100k.txt")],float))
z6=np.sort(np.asarray(np.load("../qm_riemann/128_odl_zeros6_2e6_zeros.npz")["zeros"],float))

def isaretli(d,M,tmid,Q,gap_kontrol=True):
    L=float(np.log(np.mean(tmid)/TWO_PI)); dt=d*L/TWO_PI
    Mt=M/np.sqrt(o.A_CG*L+o.ALFA_M1+2*o.A_CG+(o.ALFA_0+o.ALFA_M1)/L)
    cols=[np.ones_like(dt)]
    for (q,p) in Q: cols += [np.cos(tmid*np.log(q)), np.sin(tmid*np.log(q))]
    Xw=np.column_stack(cols)
    X=Xw if not gap_kontrol else np.column_stack([Xw,np.log(dt)])
    c,*_=np.linalg.lstsq(X,np.log(Mt),rcond=None)
    out={}
    for i,(q,p) in enumerate(Q):
        k=int(round(np.log(q)/np.log(p)))
        out[q]=(float(c[1+2*i]*np.sqrt(q)*k), float(c[2+2*i]*np.sqrt(q)*k))   # k-duzeltmeli isaretli
    return L,out

def pencere(zeros, mc, n):
    i0=int(np.argmin(np.abs(zeros-mc))); lo=max(0,i0-n//2)
    g=zeros[lo:lo+n+1]
    d=np.diff(g); tmid=0.5*(g[:-1]+g[1:])
    M=np.array([o.M_aralik(g[i],g[i+1],24) for i in range(len(d))])
    return d,M,tmid

kayit=[]
# dusuk-t pencereleri (zeros_100k)
for mc,n in [(1700,1000),(3400,2300),(6900,5200),(14000,11000),(28000,24000),(56000,40000)]:
    d,M,tmid=pencere(z100,mc,n)
    L0=float(np.log(np.mean(tmid)/TWO_PI)); Q=s7.taban_olustur(L0)
    L,c=isaretli(d,M,tmid,Q)
    kayit.append((L,c)); print(f"L={L:.3f} n={len(d):6d} taban={len(Q):3d} hesaplandı", flush=True)
# yuksek-t pencereleri (zeros6)
for mc in [1.2e5,3.5e5,1.0e6]:
    d,M,tmid=pencere(z6,mc,40000)
    L0=float(np.log(np.mean(tmid)/TWO_PI)); Q=s7.taban_olustur(L0)
    L,c=isaretli(d,M,tmid,Q)
    kayit.append((L,c)); print(f"L={L:.3f} n={len(d):6d} taban={len(Q):3d} hesaplandı", flush=True)

print("\n=== ISARETLI w (a) vs tau : butun pencereler ===")
print(f"{'tau':>6} {'a(cos)':>9} {'b(sin)':>9} {'q':>4} {'L':>6}")
satir=[]
for L,c in kayit:
    for q,(a,b) in c.items():
        tau=np.log(q)/L
        if 0.35<=tau<=0.62: satir.append((tau,a,b,q,L))
for tau,a,b,q,L in sorted(satir):
    print(f"{tau:6.3f} {a:9.4f} {b:9.4f} {q:4d} {L:6.2f}")
if satir:
    t=np.array([s[0] for s in satir]); a=np.array([s[1] for s in satir])
    m=(t>0.38)&(t<0.52)
    if m.sum()>3:
        c=np.polyfit(t[m],a[m],1); kok=-c[1]/c[0]
        print(f"\ntau[0.38,0.52] doğrusal fit: a = {c[1]:.4f} + {c[0]:.4f}·tau  ->  sıfır geçişi tau0 = {kok:.3f}   (makale: 0.447 ± 0.005)")
    print(f"kuadratür oranı |b|/|a| ortalaması (güçlü çizgiler): {np.mean([abs(s[2]/s[1]) for s in satir if abs(s[1])>0.05]):.4f}")
json.dump([(L,{str(k):v for k,v in c.items()}) for L,c in kayit], open("16_isaret.json","w"), indent=1)
