# CUE Monte Carlo check of the small-gap hump theorem (no zeta data touched)
import numpy as np, sys, json
from importlib import import_module; pred=import_module("200t_tahmin").pred
rng=np.random.default_rng(int(sys.argv[3]) if len(sys.argv)>3 else 1)
N=int(sys.argv[1]); nmat=int(sys.argv[2]); EPS=0.4; G=41
def haar(b):
    Z=(rng.standard_normal((b,N,N))+1j*rng.standard_normal((b,N,N)))/np.sqrt(2)
    Q,R=np.linalg.qr(Z); d=np.diagonal(R,axis1=1,axis2=2)
    return Q*(d/abs(d))[:,None,:]
rows=[]
for _ in range(nmat//50000):
    th=np.sort(np.angle(np.linalg.eigvals(haar(50000))),axis=1)
    nxt=np.concatenate([th[:,1:],th[:,:1]+2*np.pi],axis=1)
    s=nxt-th; su=s*N/(2*np.pi)
    b,i=np.nonzero(su<EPS)
    a=th[b,i]; ss=s[b,i]
    u=np.linspace(0,1,G)[1:-1]
    grid=a[:,None]+ss[:,None]*u[None,:]                        # (m,G-2)
    ll=np.log(np.abs(2*np.sin((grid[:,:,None]-th[b][:,None,:])/2))).sum(-1)
    k=np.clip(ll.argmax(1),1,G-4); r=np.arange(len(k))
    y0,y1,y2=ll[r,k-1],ll[r,k],ll[r,k+1]; den=y0-2*y1+y2
    logM=y1-0.125*(y2-y0)**2/np.where(den==0,-1e-300,den)       # parabolic refine
    mid=a+ss/2
    oth=np.ones(th[b].shape,bool); oth[r,i]=False; oth[r,(i+1)%N]=False
    lr=np.where(oth,np.log(np.abs(2*np.sin((mid[:,None]-th[b])/2))),0).sum(1)  # log|Lambda_{N-2}(mid)|
    rows.append(np.c_[ss,su[b,i],logM-2*np.log(ss),lr-np.log(4)])
D=np.concatenate(rows)
P=pred(N); res={'N':N,'nmat':nmat,'pred':P,'bands':{}}
for e in (0.4,0.3,0.2,0.1,0.05):
    m=D[:,1]<e; x=D[m,2]; X=np.exp(x); n=m.sum()
    res['bands'][e]=dict(n=int(n), rate=float(n/nmat/(2*np.pi*e/N)**3),
        m1=float(X.mean()), m1_se=float(X.std()/np.sqrt(n)),
        m2=float((X**2).mean()), m2_se=float((X**2).std()/np.sqrt(n)),
        Elog=float(x.mean()), Elog_se=float(x.std()/np.sqrt(n)), Vlog=float(x.var()), k3=float(((x-x.mean())**3).mean()),
        corr_s_hump=float(np.corrcoef(D[m,1],x)[0,1]),
        approx_err=float(np.abs(x-D[m,3]).max()))
print(json.dumps(res,indent=1))
