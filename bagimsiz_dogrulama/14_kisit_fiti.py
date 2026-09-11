# -*- coding: utf-8 -*-
"""
√w-v KISITI: tum satirlar (asal + asal-kuvveti) uzerinden yeniden fit.
Makale: sqrt(w) = 1.017 - 0.884 v (RMS 0.0075); alternatifler:
  (a) sqrt-w dogrusal   (b) w dogrusal   (c) yogunluk-birimlik w = 1 - v^2
"""
import numpy as np, importlib.util, json
def yukle(ad,yol):
    s=importlib.util.spec_from_file_location(ad,yol); m=importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
o=yukle("o","02_olcum.py"); s7=yukle("s7","07_wv.py")
zeros=np.sort(np.asarray(np.load("../qm_riemann/128_odl_zeros6_2e6_zeros.npz")["zeros"],float))
noktalar=[]
for mc in [1.2e5,2.0e5,3.5e5,6.0e5,1.0e6]:
    d,M,tmid=s7.pencere_verisi(zeros,mc,40000)
    L0=float(np.log(np.mean(tmid)/o.TWO_PI)); Q=s7.taban_olustur(L0)
    L,w,v,cM,cd=s7.regresyon(d,M,tmid,Q,gap_kontrol=True)
    for (q,p) in Q:
        noktalar.append((float(np.log(q)/L), float(w[q]), float(v[q]), int(q), int(p), float(L)))
print(f"toplam nokta: {len(noktalar)}  (tau aralığı {min(n[0] for n in noktalar):.3f}..{max(n[0] for n in noktalar):.3f})")
tau=np.array([n[0] for n in noktalar]); W=np.array([n[1] for n in noktalar]); V=np.array([n[2] for n in noktalar])
asal=np.array([n[3]==n[4] for n in noktalar])
def rapor(ad, W_, V_, maske=None):
    if maske is None: maske=np.ones(len(W_),bool)
    w_,v_=W_[maske],V_[maske]
    y=np.sqrt(w_); c=np.polyfit(v_,y,1); rms=np.sqrt(np.mean((y-np.polyval(c,v_))**2))
    print(f"  {ad:<22} sqrt(w) = {c[1]:.3f} + {c[0]:.3f} v   RMS={rms:.4f}  (n={maske.sum()})")
    return c,rms
print("\n=== TÜM SATIRLAR ===")
c_all,r_all=rapor("tüm satırlar", W, V)
c_pr,r_pr=rapor("yalnız asallar", W, V, asal)
print("\n=== Alternatif formlar (tüm satırlar) ===")
# (b) w dogrusal
cb=np.polyfit(V,W,1); rms_b=np.sqrt(np.mean((W-np.polyval(cb,V))**2)); print(f"  w-lineer               w = {cb[1]:.3f} + {cb[0]:.3f} v   RMS={rms_b:.4f}")
# (c) yogunluk-birimlik: w = 1 - v^2 (parametresiz)
rms_c=np.sqrt(np.mean((W-(1-V**2))**2)); print(f"  yoğunluk-birimlik      w = 1 - v² (parametresiz)  RMS={rms_c:.4f}")
rms_c2=np.sqrt(np.mean((W-(1-V)**2))**2); print(f"  genlik-birimlik        w = (1-v)² (parametresiz)  RMS={rms_c2:.4f}")
print("\n=== MAKALE: sqrt(w) = 1.017 - 0.884 v, RMS 0.0075 ===")
# makale katsayilariyla benim veride RMS
y=np.sqrt(W); pred=1.017-0.884*V
print(f"  makale katsayılarıyla benim veride RMS = {np.sqrt(np.mean((y-pred)**2)):.4f}")
# tau bantlarina gore egim
print("\n=== tau bandı başına eğim (tüm satırlar) ===")
for lo,hi in [(0.05,0.15),(0.15,0.25),(0.25,0.35),(0.35,0.46)]:
    m=(tau>=lo)&(tau<hi)
    if m.sum()>5:
        c=np.polyfit(V[m],np.sqrt(W[m]),1); rr=np.corrcoef(V[m],np.sqrt(W[m]))[0,1]
        print(f"  tau[{lo},{hi}): n={m.sum():4d}  eğim={c[0]:7.3f}  kesişim={c[1]:6.3f}  r={rr:6.3f}")
json.dump(dict(noktalar=noktalar, c_all=list(map(float,c_all)), rms_all=float(r_all),
               c_pr=list(map(float,c_pr)), rms_pr=float(r_pr), rms_makale_kats=float(np.sqrt(np.mean((y-pred)**2)))),
          open("14_kisit.json","w"), indent=1)
