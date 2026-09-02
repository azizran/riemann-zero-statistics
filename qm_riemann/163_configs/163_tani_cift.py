"""163 probe7 — T4 cekirdek: hangi (q_a,q_b) CIFTI isareti tasiyor?
   143'un verteks yasasi: G(q_a,q_b) yalnizca q_a q_b^{+-1} asal-kuvvetse buyuk."""
import importlib, sys
from pathlib import Path
import numpy as np
from sympy import factorint
sys.path.insert(0, "/Users/ugursezen/Desktop/arin/deney/qm_riemann/163_configs")
K = importlib.import_module("163_cekirdek")
TWO_PI=2*np.pi
veri=sys.argv[1]; HED=float(sys.argv[2]) if len(sys.argv)>2 else 0.54
Y=K.gaz(veri,0.40); L=Y.L; m=Y.m0
d=np.load(f"sp_{veri}_0.95.npz"); Wl,x1,TAUl,Ql = d["W"],d["x"],d["tau"],d["q"]
QS=[int(q) for q in Ql[:14]]                    # en kucuk 14 cizgi
XI=[]
for q in QS:
    i=int(np.where(Ql==q)[0][0]); w=Wl[i]; xq=x1[i]
    v=xq.real*np.cos(w*m)-xq.imag*np.sin(w*m); XI.append(v-v.mean())
Xr = Y.Xtil0 - sum(XI); Xr-=Xr.mean(); XI.append(Xr); QS.append(0)
n=len(XI)
BA=K.bant_adaylari(Y,[b for b in K.IZGARA_T1 if b[0]>=0.40-1e-9])
bd=[b for b in BA if abs(b["tau"]-HED)<1e-6][0]
acc=np.zeros((n,n)); pN=[];pO=[];Aon=[];Aof=[]; tots=[]
MON=np.zeros((n,n)); MOF=np.zeros((n,n)); son=[];sof=[]
for cz in bd["cizgi"]:
    for Wf,et in ((cz["w"],"on"),(cz["w"]+cz["gap"]/2,"off")):
        cw,sw=np.cos(Wf*m),np.sin(Wf*m)
        zc=complex(2*np.mean(Y.e0*cw),-2*np.mean(Y.e0*sw))
        cw1,sw1=np.cos(Wf*Y.mid[1:]),np.sin(Wf*Y.mid[1:])
        c1r,c1i=Y.e1*cw1,-Y.e1*sw1; zb=zc/2
        rho=c1r*zb.real+c1i*zb.imag; sig=c1i*zb.real-c1r*zb.imag; rm=rho.mean()
        A=TWO_PI*Wf/L
        M=np.zeros((n,n))
        for a in range(n):
            sa=sig*XI[a]
            for b in range(a,n):
                M[a,b]=(1.0 if a==b else 2.0)*np.mean(sa*XI[b])/rm
        if et=="on": pN.append(abs(zc)**2); Aon.append(A); MON+=abs(zc)**2*A*A*M; son.append(np.mean(sig*Y.Xtil0**2)/rm*abs(zc)**2*A*A)
        else: pO.append(abs(zc)**2); Aof.append(A); MOF+=abs(zc)**2*A*A*M; sof.append(np.mean(sig*Y.Xtil0**2)/rm*abs(zc)**2*A*A)
payda=float(np.sum(np.array(pN)-np.array(pO)))
MM=(MON-MOF)/payda; tot=(np.sum(son)-np.sum(sof))/payda
def pp(q):
    if q==0: return "art"
    return str(q)
print(f"[{veri}] bant τ={HED}  A²s2={tot:+.6f}  matris={MM.sum():+.6f}  (kul={bd['kul']})")
print("      "+"".join(f"{pp(q):>9}" for q in QS))
for a in range(n):
    print(f"{pp(QS[a]):>5} "+"".join((f"{MM[a,b]:+9.5f}" if b>=a else " "*9) for b in range(n)))
# 143 verteks etiketi
def ispp(x):
    return x>1 and len(factorint(int(x)))==1
print("\n en buyuk 12 hucre  (143 verteks: q_a*q_b ya da q_b/q_a asal-kuvvet mi?)")
idx=np.dstack(np.unravel_index(np.argsort(-np.abs(MM),axis=None),MM.shape))[0]
c=0
for a,b in idx:
    if b<a: continue
    qa,qb=QS[a],QS[b]
    et="artık" if (qa==0 or qb==0) else (
        ("Σ+Δ" if ispp(qa*qb) and (qb%qa==0 and ispp(qb//qa)) else
         ("Σ" if ispp(qa*qb) else ("Δ" if (qb%qa==0 and ispp(qb//qa)) or (qa%qb==0 and ispp(qa//qb)) else "—"))))
    print(f"  ({pp(qa)},{pp(qb)}): {MM[a,b]:+9.5f}   143-verteks: {et}")
    c+=1
    if c>=12: break
