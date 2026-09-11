"""
BAGIMSIZ w-v REGRESYONU (Not 2/3):
  unfold:  d_tilde = d*L/2pi ,  M_tilde = M/sqrt(AL+a_m1+2A+(a0+a_m1)/L)
  taban Q: log q/L <= 0.45, q<=720 olan TUM asal kuvvetleri
  log M_tilde ~ 1 + sum[a cos(t_n log q) + b sin(...)] + c*log d_tilde   -> w_q = |A_q|*sqrt(q)
  log d_tilde ~ ayni taban (kontrolsuz)                                    -> v_q = A_q*sqrt(q)
"""
import numpy as np, json, importlib.util, sys
spec=importlib.util.spec_from_file_location("o","02_olcum.py"); o=importlib.util.module_from_spec(spec); spec.loader.exec_module(o)
TWO_PI=o.TWO_PI; A=o.A_CG; AM1=o.ALFA_M1; A0=o.ALFA_0

def asal_kuvvetleri(limit):
    ps=[p for p in range(2,limit+1) if all(p%d for d in range(2,int(p**0.5)+1))]
    out=[]
    for p in ps:
        q=p
        while q<=limit: out.append((q,p)); q*=p
    return sorted(out)

def pencere_verisi(zeros, merkez, n_gap, k=24):
    i0=int(np.argmin(np.abs(zeros-merkez))); lo=max(0,i0-n_gap//2)
    g=zeros[lo:lo+n_gap+1]
    d=np.diff(g); tmid=0.5*(g[:-1]+g[1:])
    M=np.array([o.M_aralik(g[i],g[i+1],k) for i in range(len(d))])
    return d, M, tmid

def taban_olustur(L, qmax=720, tau_max=0.45):
    Q=[(q,p) for (q,p) in asal_kuvvetleri(qmax) if np.log(q)/L<=tau_max]
    return Q

def regresyon(d, M, tmid, Q, gap_kontrol=False):
    L=float(np.log(np.mean(tmid)/TWO_PI))
    dt=d*L/TWO_PI
    Mt=M/np.sqrt(A*L+AM1+2*A+(A0+AM1)/L)
    cols=[np.ones_like(dt)]
    for (q,p) in Q:
        cols.append(np.cos(tmid*np.log(q))); cols.append(np.sin(tmid*np.log(q)))
    yM=np.log(Mt); yd=np.log(dt)
    Xw=np.column_stack(cols)                      # dalga tabani (her iki kanal icin ortak)
    cd,*_=np.linalg.lstsq(Xw,yd,rcond=None)       # v kanali: gap ~ taban
    X=Xw if not gap_kontrol else np.column_stack([Xw, yd])
    cM,*_=np.linalg.lstsq(X,yM,rcond=None)        # w kanali: max ~ taban + kendi gap kontrolu
    w={}; u={}; v={}
    for i,(q,p) in enumerate(Q):
        aM,bM=cM[1+2*i],cM[2+2*i]; ad,bd=cd[1+2*i],cd[2+2*i]
        w[q]=float(np.hypot(aM,bM)*np.sqrt(q)); u[q]=float(np.hypot(aM,bM)*np.sqrt(q))
        v[q]=float(np.hypot(ad,bd)*np.sqrt(q))
    return L, w, v, cM, cd

if __name__=="__main__":
    zeros=np.sort(np.asarray(np.load("../qm_riemann/128_odl_zeros6_2e6_zeros.npz")["zeros"],float))
    merkezler=[1.2e5,2.0e5,3.5e5,6.0e5,1.0e6]
    kayit=[]
    for mc in merkezler:
        d,M,tmid=pencere_verisi(zeros,mc,40000)
        L0=float(np.log(np.mean(tmid)/TWO_PI)); Q=taban_olustur(L0)
        L,w,v,cM,cd=regresyon(d,M,tmid,Q,gap_kontrol=True)
        # kontrolsuz (kosulsuz icerik u_q) icin gap kontrolu olmadan M regresyonu
        _,u,_,_,_=regresyon(d,M,tmid,Q,gap_kontrol=False)
        kayit.append(dict(L=L,n=len(d),Q=len(Q),w=w,v=v,u=u))
        print(f"L={L:.3f} n={len(d)} taban={len(Q)} satir | w(2)={w[2]:.3f} w(3)={w[3]:.3f} w(5)={w[5]:.3f} w(7)={w[7]:.3f} w(11)={w[11]:.3f} w(13)={w[13]:.3f} | u(2)={u[2]:.3f} u(3)={u[3]:.3f}")
    print("\n--- w (gap-kosullu) ortalama, makale tablosu: 0.83 0.72 0.60 0.52 0.42 0.38 ---")
    for p in [2,3,5,7,11,13]:
        print(f"  p={p:2d}: w={np.mean([k['w'][p] for k in kayit]):.3f}   v={np.mean([k['v'][p] for k in kayit]):.3f}  (makale |v|: {dict(zip([2,3,5,7,11,13],[0.11,0.18,0.25,0.30,0.35,0.37]))[p]})")
    print("\n--- u (kosulsuz icerik, toplam kural: asallarda ~1) ---")
    for p in [2,3,5,7,11,13]: print(f"  p={p:2d}: u={np.mean([k['u'][p] for k in kayit]):.3f}")
    # w-v kisiti: tum (tau, w, v) noktalari
    pts=[]
    for k in kayit:
        for q in k["w"]:
            tau=np.log(q)/k["L"]
            pts.append((tau,k["w"][q],k["v"][q],q))
    import numpy as np
    wv=np.array([[p[1],p[2]] for p in pts if p[3] in (2,3,5,7,11,13)])
    sw=np.sqrt(wv[:,0]); vv=wv[:,1]
    A_fit=np.polyfit(vv,sw,1)
    print(f"\n--- KISIT: sqrt(w) = {A_fit[1]:.3f} + {A_fit[0]:.3f}*v   (makale: 1.017 - 0.884 v, yani v isaretli negatifse egim +0.884)")
    r=np.corrcoef(vv,sw)[0,1]; print(f"    korelasyon r={r:.4f}, RMS artik={np.std(sw-np.polyval(A_fit,vv)):.5f}")
    json.dump(kayit, open("07_wv_sonuc.json","w"), indent=1)
