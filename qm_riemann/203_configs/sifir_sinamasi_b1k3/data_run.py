import sys, time, json, numpy as np
sys.dont_write_bytecode=True
sys.path.insert(0,'.')
import importlib.util, os
QM=os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','..'))   # qm_riemann/
from cumlib import *; from shifts import *
spec=importlib.util.spec_from_file_location('veri',os.path.join(QM,'200_configs','200c_veri.py'))
v=importlib.util.module_from_spec(spec); spec.loader.exec_module(v)
taban,off,bl=v.oku(os.path.join(QM,'veri_lmfdb','zeros_8846000.dat'),n_max=1_500_000)
N=len(off); KMAX=int(sys.argv[1]) if len(sys.argv)>1 else 800
SNAPS=[k for k in (50,100,200,400,800,1600) if k<=KMAX]
alphas=np.array(alpha_list()); na=len(alphas)
i0=KMAX; i1=N-KMAX
NB0=128
n_use=((i1-i0)//NB0)*NB0; i1=i0+n_use
Li=np.log((taban+off)/(2*np.pi))
Lbar=float(Li[i0:i1].mean()); a_raw=alphas/Lbar
print('N',N,'eval',i0,i1,n_use,'Lbar',Lbar,flush=True)
re=np.zeros((na,N)); im=np.zeros((na,N)); a2=a_raw**2
out={'Lbar':Lbar,'alphas':alphas.tolist(),'K':{}}
def build_X(K):
    """complex X arrays on eval slice for window +-K"""
    idx=np.arange(i0,i1)
    Wi=(off[idx+K]-off[idx-K])/2+np.pi/Li[idx]
    X=[]
    for k in range(na):
        x=(1.0/a_raw[k]+re[k,i0:i1])+1j*im[k,i0:i1]
        x=x-Li[i0:i1]/2 + (Li[i0:i1]/np.pi)*np.arctan(a_raw[k]/Wi)
        X.append(x)
    return X
def analyse(X,sl=None,NB=32):
    ai={round(a,10):k for k,a in enumerate(alphas)}
    g=lambda al: X[ai[round(al,10)]][sl] if sl is not None else X[ai[round(al,10)]]
    res={'triple':{}, 'pair_XXB':{}, 'pair_XX':{}, 'mean':{}}
    for lab,nom,act in TRIPLES:
        a,b,c=act
        r1=cum3_jk(g(a),g(b),np.conj(g(c)),NB)
        r2=cum3_jk(g(a),g(b),g(c),NB)
        res['triple'][lab]={'XXB':[r1[0].real,r1[1],r1[0].imag,r1[2]],'XXX':[r2[0].real,r2[1],r2[0].imag,r2[2]]}
    for (a,b) in PAIRS_XXB:
        r=cum2_jk(g(a),np.conj(g(b)),NB); res['pair_XXB'][f'{a},{b}']=[r[0].real,r[1],r[0].imag,r[2]]
    for (a,b) in PAIRS_XX:
        r=cum2_jk(g(a),g(b),NB); res['pair_XX'][f'{a},{b}']=[r[0].real,r[1],r[0].imag,r[2]]
    for al in alphas: res['mean'][f'{al}']=[complex(g(al).mean()).real,complex(g(al).mean()).imag]
    return res
t0=time.time()
for d in range(1,KMAX+1):
    u=off[d:]-off[:-d]; u2=u*u
    for k in range(na):
        inv=1.0/(a2[k]+u2); w=a_raw[k]*inv
        re[k,:N-d]+=w; re[k,d:]+=w
        w2=u*inv
        im[k,:N-d]+=w2; im[k,d:]-=w2
    if d%50==0: print(d,round(time.time()-t0),'s',flush=True)
    if d in SNAPS:
        X=build_X(d)
        out['K'][d]={'full':analyse(X),'full_NB64':analyse(X,NB=64),'full_NB128':analyse(X,NB=128)}
        n=len(X[0]); h=n//2
        out['K'][d]['half1']=analyse(X,slice(0,h)); out['K'][d]['half2']=analyse(X,slice(h,n))
        out['K'][d]['mid']=analyse(X,slice(n//4,n//4+n//2))
        if d==400: np.save('X400.npy',np.array(X))
        json.dump(out,open('data_out.json','w'))
        print('snapshot',d,'done',flush=True)
json.dump(out,open('data_out.json','w'))
print('ALL DONE',time.time()-t0)
