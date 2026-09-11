# -*- coding: utf-8 -*-
"""
v(τ) KOLAPSININ DERİN ZİNCİRİ: 10^5 → 10^12 → 10^21 → 10^22
Yalnız aralık verisi (M gerekmez). Eşleşen τ'lerde asal↔asal karşılaştırma.
"""
import numpy as np, importlib.util, json
def yukle(ad,yol):
    s=importlib.util.spec_from_file_location(ad,yol); m=importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
o=yukle("o","02_olcum.py"); s7=yukle("s7","07_wv.py"); TWO_PI=2*np.pi

def asal(n): return n>1 and all(n%d for d in range(2,int(n**0.5)+1))

def v_egri(zeros):
    """zeros: MUTLAK gamma (float64'te guvenli oldugu yukseklikler icin).
    Faz sutunlari icin ortalamadan arindirilmis t kullanilir (kaydirma-invaryant)."""
    d=np.diff(zeros); tmid=0.5*(zeros[:-1]+zeros[1:])
    L=float(np.log(np.mean(tmid)/TWO_PI)); dt=d*L/TWO_PI
    return _fit(dt, tmid, L)

def v_egri_ofset(ofset, taban):
    """Derin tablolar: gamma = taban + ofset. float64 tabani yutar; faz icin
    ofset kullanilir (sabit faz kaymasi cos/sin ciftince yutulur => esdeger)."""
    d=np.diff(ofset); t_true=taban+0.5*(ofset[:-1]+ofset[1:])
    L=float(np.log(np.mean(t_true)/TWO_PI)); dt=d*L/TWO_PI
    return _fit(dt, 0.5*(ofset[:-1]+ofset[1:]), L)

def _fit(dt, faz_zamani, L):
    Q=s7.taban_olustur(L)
    cols=[np.ones_like(dt)]
    z=faz_zamani-faz_zamani.mean()          # sayisal kararlilik
    for (q,p) in Q: cols += [np.cos(z*np.log(q)), np.sin(z*np.log(q))]
    c,*_=np.linalg.lstsq(np.column_stack(cols),np.log(dt),rcond=None)
    v={q: float(np.hypot(c[1+2*i],c[2+2*i])*np.sqrt(q)) for i,(q,p) in enumerate(Q)}
    return L,{q:x for q,x in v.items() if asal(q)}

def tablo_oku(dosya, taban):
    sat=[l.strip() for l in open(dosya) if l.strip()]
    vals=np.array([float(x) for x in sat if x.replace('.','',1).replace('-','',1).isdigit()])
    return taban+vals

# referans: L=9.86 penceresi (zeros6)
z6=np.sort(np.asarray(np.load("../qm_riemann/128_odl_zeros6_2e6_zeros.npz")["zeros"],float))
i0=int(np.argmin(np.abs(z6-1.2e5))); L_ref,v_ref=v_egri(z6[i0-20000:i0+20001])
def ofset_oku(dosya):
    sat=[l.strip() for l in open(dosya) if l.strip()]
    return np.array([float(x) for x in sat if x.replace('.','',1).replace('-','',1).isdigit()])
B12,B21,B22 = 267653395647.0, 144176897509546973000.0, 1370919909931995300000.0
a12,a21,a22 = ofset_oku("../zeros3.txt"), ofset_oku("../zeros4.txt"), ofset_oku("../zeros5.txt")
L12,v12=v_egri_ofset(a12,B12); L21,v21=v_egri_ofset(a21,B21); L22,v22=v_egri_ofset(a22,B22)
z12,z21,z22 = B12+a12, B21+a21, B22+a22
print(f"pencereler: L_ref={L_ref:.3f} (t≈1.2e5) | L12={L12:.3f} (t≈{z12[0]:.3e}) | "
      f"L21={L21:.3f} (t≈{z21[0]:.3e}) | L22={L22:.3f} (t≈{z22[0]:.3e})")

def karsilastir(ad, La, va, Lb, vb):
    """va'yı vb'nin τ'larında interpole edip oran verir."""
    ta=np.array([np.log(q)/La for q in sorted(va)]); xa=np.array([va[q] for q in sorted(va)])
    o_=np.argsort(ta); ta,xa=ta[o_],xa[o_]
    oranlar=[]
    for q in sorted(vb):
        t=np.log(q)/Lb
        if ta.min()<=t<=ta.max():
            oranlar.append(vb[q]/np.interp(t,ta,xa))
    if oranlar:
        a=np.array(oranlar)
        print(f"  {ad:<26} n={len(a):3d}  oran ort={a.mean():.3f} ± {a.std():.3f}   "
              f"[τ aralığı {min(np.log(q)/Lb for q in vb if ta.min()<=np.log(q)/Lb<=ta.max()):.3f}–"
              f"{max(np.log(q)/Lb for q in vb if ta.min()<=np.log(q)/Lb<=ta.max()):.3f}]")
    else:
        print(f"  {ad:<26} örtüşme yok")
    return oranlar

print("\nv(derin)/v(referans veya 10^12), eşleşen τ'lerde (yalnız asallar):")
karsilastir("10^12 / 10^5",  L_ref,v_ref, L12,v12)
karsilastir("10^21 / 10^12", L12,v12,    L21,v21)
karsilastir("10^22 / 10^12", L12,v12,    L22,v22)
karsilastir("10^22 / 10^21", L21,v21,    L22,v22)
karsilastir("10^21 / 10^5",  L_ref,v_ref, L21,v21)
karsilastir("10^22 / 10^5",  L_ref,v_ref, L22,v22)
print(f"\nt oranı 10^22/10^5 = {z22[0]/1.2e5:.2e}  (~{np.log10(z22[0]/1.2e5):.1f} mertebe)")
json.dump(dict(L_ref=L_ref,L12=L12,L21=L21,L22=L22,
               v_ref={str(k):v for k,v in v_ref.items()}, v12={str(k):v for k,v in v12.items()},
               v21={str(k):v for k,v in v21.items()}, v22={str(k):v for k,v in v22.items()}),
          open("20_zincir.json","w"), indent=1)
