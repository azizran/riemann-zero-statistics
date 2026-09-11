# -*- coding: utf-8 -*-
"""
PLASEBO: toplam kurali (u_q ~ 1) aritmetik mi, regresyon artefakti mi?
Ayni guc tayfli Gaussian surecler (faz rastgele) uzerinde AYNI regresyonu kos.
Gercek denizde u~1 cikiyorsa ve plaseboda u~0 ise: icerik aritmetik.
"""
import numpy as np, json, importlib.util
def yukle(ad,yol):
    s=importlib.util.spec_from_file_location(ad,yol); m=importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
s3=yukle("s3","03_surrogate.py"); s7=yukle("s7","07_wv.py")
TWO_PI=s3.TWO_PI

t0,t1,h=107252.0,132748.0,0.01992
ts=t0+h*np.arange(int((t1-t0)/h))
yz=s3.Z(ts)
L0=float(np.log(np.mean(ts)/TWO_PI)); Q=s7.taban_olustur(L0)
print(f"izgara {ts.size} nokta, L={L0:.3f}, taban {len(Q)} satır")

def regres_u(ys):
    d,M,tmid=s3.sifir_ve_maks(ts,ys)
    L=float(np.log(np.mean(tmid)/TWO_PI)); dt=d*L/TWO_PI
    Mt=M/np.sqrt(s7.A*L+s7.AM1+2*s7.A+(s7.A0+s7.AM1)/L)
    cols=[np.ones_like(dt)]
    for (q,p) in Q: cols += [np.cos(tmid*np.log(q)), np.sin(tmid*np.log(q))]
    X=np.column_stack(cols); c,*_=np.linalg.lstsq(X,np.log(Mt),rcond=None)
    return {q: float(np.hypot(c[1+2*i],c[2+2*i])*np.sqrt(q)) for i,(q,p) in enumerate(Q)}, L

u_gercek,Lg=regres_u(yz)
rng=np.random.default_rng(4242); F=np.fft.rfft(yz); amp=np.abs(F); n=ts.size
n_sur=20; u_null={q:[] for q,_ in Q}
for k in range(n_sur):
    faz=rng.uniform(0,2*np.pi,amp.size); faz[0]=0.0
    if n%2==0: faz[-1]=0.0
    ys=np.fft.irfft(amp*np.exp(1j*faz),n=n)
    uk,_=regres_u(ys)
    for q in u_null: u_null[q].append(uk[q])
print(f"\n{'q':>4} {'tau':>6} | {'u_gercek':>9} {'u_null ort':>11} {'u_null std':>11} {'ayrim':>8}")
for (q,p) in Q:
    if p!=q or q>13: continue
    arr=np.array(u_null[q])
    print(f"{q:4d} {np.log(q)/Lg:6.3f} | {u_gercek[q]:9.4f} {arr.mean():11.4f} {arr.std(ddof=1):11.4f} {(u_gercek[q]-arr.mean())/arr.std(ddof=1):7.0f}s")
allq=[q for (q,p) in Q if q>13][:0]
A=np.array([np.array(u_null[q]).mean() for (q,p) in Q if q<=13 and p==q])
print(f"\nplasebo tabanı (ilk 6 asal, ortalama): {A.mean():.4f}")
print(f"gerçek içerik        (ilk 6 asal, ortalama): {np.mean([u_gercek[q] for q in [2,3,5,7,11,13]]):.4f}")
json.dump(dict(u_gercek={str(k):v for k,v in u_gercek.items()},
               u_null={str(k):list(map(float,v)) for k,v in u_null.items()}, L=Lg),
          open("11_placebo.json","w"), indent=1)
