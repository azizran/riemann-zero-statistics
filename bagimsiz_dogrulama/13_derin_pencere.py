# -*- coding: utf-8 -*-
"""
10^12 PENCERESI: v kanali (yalniz araliklardan) — tau-yasasi 2.5e6x yuksek t'de tutuyor mu?
Veri: Odlyzko zeros3 (taban 267653395647 + degerler), 10.000 sifir.
"""
import numpy as np, importlib.util, json
def yukle(ad,yol):
    s=importlib.util.spec_from_file_location(ad,yol); m=importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
o=yukle("o","02_olcum.py"); s7=yukle("s7","07_wv.py")
TWO_PI=2*np.pi

satirlar=[l.strip() for l in open("../zeros3.txt") if l.strip()]
vals=np.array([float(x) for x in satirlar if x.replace('.','',1).replace('-','',1).isdigit()])
zeros=267653395647.0+vals
print(f"10^12 penceresi: {zeros.size} sifir, t={zeros[0]:.3f} .. {zeros[-1]:.3f}")
d=np.diff(zeros); tmid=0.5*(zeros[:-1]+zeros[1:])
L=float(np.log(np.mean(tmid)/TWO_PI)); dt=d*L/TWO_PI
print(f"L={L:.4f}  ortalama aralik={d.mean():.6f}  ({tmid.mean():.4g})")
Q=s7.taban_olustur(L, tau_max=0.45)
print(f"taban: {len(Q)} asal kuvveti (tau<=0.45, q<=720)")

# v kanali: log dt ~ taban
cols=[np.ones_like(dt)]
for (q,p) in Q: cols += [np.cos(tmid*np.log(q)), np.sin(tmid*np.log(q))]
X=np.column_stack(cols)
cd,*_=np.linalg.lstsq(X,np.log(dt),rcond=None)
v={q: float(np.hypot(cd[1+2*i],cd[2+2*i])*np.sqrt(q)) for i,(q,p) in enumerate(Q)}
print(f"\n{'q':>4} {'tau':>6} {'v(10^12)':>9} | {'v(L=9.86)':>9} {'v(L=11.98)':>10}   makale |v|")
makale={2:0.11,3:0.18,5:0.25,7:0.30,11:0.35,13:0.37}
# karsilastirma icin diger pencereler
dz=np.load("../qm_riemann/128_odl_zeros6_2e6_zeros.npz"); z6=np.sort(np.asarray(dz["zeros"],float))
def v_pencere(mc):
    i0=int(np.argmin(np.abs(z6-mc))); g=z6[i0-20000:i0+20001]
    dd=np.diff(g); tm=0.5*(g[:-1]+g[1:]); LL=float(np.log(np.mean(tm)/TWO_PI))
    QQ=s7.taban_olustur(LL); cc=[np.ones_like(tm)]
    for (q,p) in QQ: cc += [np.cos(tm*np.log(q)), np.sin(tm*np.log(q))]
    c,*_=np.linalg.lstsq(np.column_stack(cc),np.log(dd*LL/TWO_PI),rcond=None)
    return LL,{q: float(np.hypot(c[1+2*i],c[2+2*i])*np.sqrt(q)) for i,(q,p) in enumerate(QQ)}
L1,v1=v_pencere(1.2e5); L2,v2=v_pencere(1.0e6)
for q in [2,3,5,7,11,13]:
    print(f"{q:4d} {np.log(q)/L:6.3f} {v[q]:9.4f} | {v1[q]:9.4f} {v2[q]:10.4f}   {makale[q]:.2f}")
print(f"\nTAU KOLAPSI: 10^12'de tau(2..13) = " + ", ".join(f"{np.log(q)/L:.3f}" for q in [2,3,5,7,11,13]))
json.dump(dict(L=L, v={str(k):v_ for k,v_ in v.items()}, L1=L1, v1={str(k):v_ for k,v_ in v1.items()},
               L2=L2, v2={str(k):v_ for k,v_ in v2.items()}), open("13_derin.json","w"), indent=1)
