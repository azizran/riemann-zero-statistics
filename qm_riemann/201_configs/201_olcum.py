"""KALEM 201 — faz kilidi ölçümü (YALNIZ sıfır konumları). Pencere başına, basamak b ∈ {0,1,2(ε)},
asal p ∈ {2,3,5,7,11,13}, harmonik m ∈ {1,2}: 64 t-bloğunda Σcos(mφ), Σsin(mφ), n. φ = t·log p (mod 2π),
t = taban(tamsayı) + ofset; taban·log p mod 2π mpmath ile. --sentetik: M6 için sentetik konumlar."""
import numpy as np, mpmath as mp, importlib.util, sys, json, argparse
from pathlib import Path
HERE=Path(__file__).resolve().parent; QM=HERE.parent
spec=importlib.util.spec_from_file_location('v',QM/'200_configs/200c_veri.py'); V=importlib.util.module_from_spec(spec); spec.loader.exec_module(V)
P=[2,3,5,7,11,13]; M=[1,2]; NB=64; EPS=[0.1,0.2,0.3]
def pencereler():
    z=np.load(QM/'128_odl_zeros6_2e6_zeros.npz')['zeros']
    for ad,(i0,i1) in {'W1':(20868,198239),'W2':(198239,598742),'W3':(598742,2001052)}.items(): yield ad,0,z[i0:i1]
    for ad,f in (('C1','zeros_8846000'),('C2','zeros_99146000'),('C3','zeros_997946000'),('C4','zeros_30599546000')):
        tb,o,_=V.oku(QM/f'veri_lmfdb/{f}.dat',n_max=1_500_000); yield ad,int(tb),o
def blok_toplam(base,off,t_a,t_b,ofs_blk):
    blk=np.clip(((ofs_blk-t_a)/(t_b-t_a)*NB).astype(int),0,NB-1)
    out=np.zeros((len(P),len(M),NB,3))
    for i,p in enumerate(P):
        faz0=float(mp.fmod(mp.mpf(base)*mp.log(p),2*mp.pi)); lp=float(mp.log(p))
        for j,m in enumerate(M):
            ph=m*(faz0+off*lp)
            out[i,j,:,0]=np.bincount(blk,np.cos(ph),NB); out[i,j,:,1]=np.bincount(blk,np.sin(ph),NB); out[i,j,:,2]=np.bincount(blk,minlength=NB)
    return out
def olc(ad,base,o,rng):
    t_a,t_b=o[0],o[-1]; L=lambda x: np.log((base+x)/(2*np.pi))
    R={'L1':float(L(o).mean())}
    R['b1']=blok_toplam(base,o,t_a,t_b,o)
    u=rng.uniform(t_a,t_b,1_000_000); R['b0']=blok_toplam(base,u,t_a,t_b,u)
    m=(o[1:]+o[:-1])/2; dt=np.diff(o)*L(m)/(2*np.pi)
    for e in EPS:
        s=dt<e; R[f'b2_{e}']=blok_toplam(base,m[s],t_a,t_b,m[s]); R[f'L2_{e}']=float(L(m[s]).mean())
    return R
if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--sentetik',action='store_true'); a=ap.parse_args()
    rng=np.random.default_rng(201); sonuc={}
    if a.sentetik:
        # M6: yoğunluğu LG ile modüle edilmiş sentetik süreç (inceltme); çiftler için orta noktalar n(t)^k ile
        L0=11.0; T0=5e5; span=6e5; n0=L0/(2*np.pi)
        def mod(t): return 1-(2/L0)*sum(np.log(p)/np.sqrt(p)*np.cos(t*np.log(p))+np.log(p)/p*np.cos(2*t*np.log(p)) for p in P)
        c=np.sort(rng.uniform(T0,T0+span,int(1.3*n0*span))); keep=rng.uniform(0,1.3,len(c))<mod(c); z=c[keep]
        S={}
        for k in (1,2,4):
            cm=rng.uniform(T0,T0+span,400000); w=np.clip(mod(cm),0,None)**k; kk=rng.uniform(0,w.max(),len(cm))<w; S[k]=cm[kk]
        sonuc['S1']={'L1':L0,'b1':blok_toplam(0,z,T0,T0+span,z)}
        for k,mm in S.items(): sonuc[f'S2_k{k}']={'L1':L0,'b1':blok_toplam(0,mm,T0,T0+span,mm)}
        np.savez('201_configs/201_sentetik.npz',**{f'{a}__{b}':v for a,d in sonuc.items() for b,v in d.items() if not isinstance(v,float)},**{f'{a}__L':d['L1'] for a,d in sonuc.items()})
    else:
        for ad,base,o in pencereler(): sonuc[ad]=olc(ad,base,o,rng); print(ad,'tamam',flush=True)
        np.savez('201_configs/201_olcum.npz',**{f'{a}__{b}':v for a,d in sonuc.items() for b,v in d.items()})
    print('kaydedildi')
