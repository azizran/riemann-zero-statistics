# -*- coding: utf-8 -*-
"""w icin gap-kosullama konvansiyonlari — normal denklemlerle (hizli).
 V1: duz log-gap kontrolu   V2: ayni-frekans gap kontrolu   V3(ref): kontrolsuz (=u)"""
import numpy as np, importlib.util, json, sys
def yukle(ad,yol):
    s=importlib.util.spec_from_file_location(ad,yol); m=importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
o=yukle("o","02_olcum.py"); s7=yukle("s7","07_wv.py")
zeros=np.sort(np.asarray(np.load("../qm_riemann/128_odl_zeros6_2e6_zeros.npz")["zeros"],float))
makale={2:0.83,3:0.72,5:0.60,7:0.52,11:0.42,13:0.38}

def cozum(G,b):
    return np.linalg.solve(G,b)

def isle(d,M,tmid,Q):
    L=float(np.log(np.mean(tmid)/o.TWO_PI)); dt=d*L/o.TWO_PI
    Mt=M/np.sqrt(o.A_CG*L+o.ALFA_M1+2*o.A_CG+(o.ALFA_0+o.ALFA_M1)/L)
    yM=np.log(Mt); yd=np.log(dt)
    cols=[np.ones_like(dt)]
    for (q,p) in Q: cols += [np.cos(tmid*np.log(q)), np.sin(tmid*np.log(q))]
    Xb=np.column_stack(cols)                      # (n, 1+2Q)
    G=Xb.T@Xb; b=Xb.T@yM
    # V1: duz gap kontrolu
    z=yd.reshape(-1,1); G1=np.block([[G, Xb.T@z],[z.T@Xb, z.T@z]]); b1=np.r_[b, z.T@yM]
    c1=cozum(G1,b1)
    # V2: her q icin ayni-frekans gap sutunlari
    out={}
    for i,(q,p) in enumerate(Q):
        Z=np.column_stack([np.cos(tmid*np.log(q)), np.sin(tmid*np.log(q))])
        G2=np.block([[G, Xb.T@Z],[Z.T@Xb, Z.T@Z]]); b2=np.r_[b, Z.T@yM]
        c2=cozum(G2,b2)
        out[q]=float(np.hypot(c2[1+2*i],c2[2+2*i])*np.sqrt(q))     # dalga katsayisi (kontrol sutunlari sona eklendi)
    w1={q: float(np.hypot(c1[1+2*i],c1[2+2*i])*np.sqrt(q)) for i,(q,p) in enumerate(Q)}
    u ={q: float(np.hypot(cozum(G,b)[1+2*i],cozum(G,b)[2+2*i])*np.sqrt(q)) for i,(q,p) in enumerate(Q)}
    return L,w1,out,u

merkezler=[1.2e5,2.0e5,3.5e5,6.0e5,1.0e6]
top={k:{q:[] for q in makale} for k in ["V1","V2","u"]}
for mc in merkezler:
    d,M,tmid=s7.pencere_verisi(zeros,mc,40000)
    L0=float(np.log(np.mean(tmid)/o.TWO_PI)); Q=s7.taban_olustur(L0)
    L,w1,w2,u=isle(d,M,tmid,Q)
    for q in makale: top["V1"][q].append(w1[q]); top["V2"][q].append(w2[q]); top["u"][q].append(u[q])
    print(f"L={L:.3f} tamam", flush=True)
print(f"\n{'varyant':>8} | " + " ".join(f"w({p})".rjust(7) for p in [2,3,5,7,11,13]))
for k in ["V1","V2","u"]:
    print(f"{k:>8} | " + " ".join(f"{np.mean(top[k][q]):7.3f}" for q in [2,3,5,7,11,13]))
print(f"{'MAKALE':>8} | " + " ".join(f"{makale[q]:7.2f}" for q in [2,3,5,7,11,13]))
json.dump({k:{str(q):v for q,v in d2.items()} for k,d2 in top.items()}, open("12_w_konv.json","w"), indent=1)
