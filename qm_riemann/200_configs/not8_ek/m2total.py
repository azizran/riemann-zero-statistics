# CUE MC: E sum_n M_n^2 over all N gaps, vs N^2 (no zeta data)
import numpy as np, sys
rng=np.random.default_rng(11)
def haar(b,N):
    Z=(rng.standard_normal((b,N,N))+1j*rng.standard_normal((b,N,N)))/np.sqrt(2)
    Q,R=np.linalg.qr(Z); d=np.diagonal(R,axis1=1,axis2=2)
    return Q*(d/abs(d))[:,None,:]
G=65
for N in (6,10,16,24,32):
    nmat=40000 if N<=16 else 20000; tot=[]
    for _ in range(nmat//2000):
        th=np.sort(np.angle(np.linalg.eigvals(haar(2000,N))),axis=1)
        nxt=np.concatenate([th[:,1:],th[:,:1]+2*np.pi],axis=1); s=nxt-th
        u=np.linspace(0,1,G)[1:-1]
        grid=th[:,:,None]+s[:,:,None]*u[None,None,:]          # (b,N,G-2)
        ll=np.log(np.abs(2*np.sin((grid[:,:,:,None]-th[:,None,None,:])/2))).sum(-1)
        k=np.clip(ll.argmax(-1),1,G-4)
        y0=np.take_along_axis(ll,(k-1)[...,None],-1)[...,0]; y1=np.take_along_axis(ll,k[...,None],-1)[...,0]; y2=np.take_along_axis(ll,(k+1)[...,None],-1)[...,0]
        den=y0-2*y1+y2; logM=y1-0.125*(y2-y0)**2/np.where(den==0,-1e-300,den)
        tot.append(np.exp(2*logM).sum(1))
    t=np.concatenate(tot)
    print(N, nmat, 'E sum M^2 / N^2 =', t.mean()/N**2, '+-', t.std()/np.sqrt(len(t))/N**2, ' /(N(N+1))', t.mean()/(N*(N+1)), flush=True)
print('(e^2-5)/2 =', (np.e**2-5)/2)
