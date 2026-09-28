"""KALEM 201 — analiz (önceden kayıtlı kurallar)."""
import numpy as np, json, sys
P=np.array([2,3,5,7,11,13]); NB=64
def lg(m,L): return -np.log(P)/(P**(m/2)*L)
def jk(B):  # B: (6 asal, NB, 3) → ĉ (6,), kovaryans (6,6)
    tot=B.sum(1); c=tot[:,0]/tot[:,2]
    reps=np.array([ (tot[:,0]-B[:,i,0])/(tot[:,2]-B[:,i,2]) for i in range(NB)])
    d=reps-reps.mean(0); return c,(NB-1)/NB*d.T@d
def gls(c,C,t):
    Ci=np.linalg.inv(C); v=1/(t@Ci@t); return float(v*(t@Ci@c)), float(np.sqrt(v))
def pool(xs):
    w=np.array([1/s**2 for _,s in xs]); x=np.array([a for a,_ in xs]); return float((w*x).sum()/w.sum()), float(1/np.sqrt(w.sum()))
def analiz(f,pencereler,etiket_b='b1',Lkey='L1'):
    d=np.load(f); out={}
    for w in pencereler:
        for mi,m in enumerate((1,2)):
            B=d[f'{w}__{etiket_b}'][:,mi]; c,C=jk(B); L=float(d[f'{w}__{Lkey}'])
            out[(w,m)]=gls(c,C,lg(m,L))+(L,)
    return out
if __name__=='__main__':
    if sys.argv[1:]==['--sentetik']:
        d=np.load('201_configs/201_sentetik.npz')
        for k in ('S1','S2_k1','S2_k2','S2_k4'):
            for mi,m in enumerate((1,2)):
                c,C=jk(d[f'{k}__b1'][:,mi]); R,s=gls(c,C,lg(m,float(d[f'{k}__L']))); print(f'{k} m={m}: R = {R:.3f} ± {s:.3f}')
        sys.exit()
    W=['W1','W2','W3','C1','C2','C3','C4']; d=np.load('201_configs/201_olcum.npz'); H={}
    r1=analiz('201_configs/201_olcum.npz',W,'b1','L1')
    R1=pool([r1[(w,1)][:2] for w in W]); R1m2=pool([r1[(w,2)][:2] for w in W])
    Ls=np.array([r1[(w,1)][2] for w in W]); Rs=np.array([r1[(w,1)][0] for w in W]); Ss=np.array([r1[(w,1)][1] for w in W])
    A=np.vstack([np.ones_like(Ls),np.log(Ls)]).T; Wt=np.diag(1/Ss**2); beta=np.linalg.solve(A.T@Wt@A,A.T@Wt@Rs); egim=float(beta[1]); egim_se=float(np.sqrt(np.linalg.inv(A.T@Wt@A)[1,1]))
    H['H-201-1']=dict(R1=R1,R1_m2=R1m2,egim_logL=(egim,egim_se),pencere={w:r1[(w,1)] for w in W},
        karar='TUTAR' if abs(R1[0]-1)<=0.15 and abs(R1[0]-1)<=3*R1[1] else 'TUTMADI', olcekleme='TUTAR' if abs(egim)<=0.1 else 'TUTMADI')
    for e in (0.1,0.2,0.3):
        r2=analiz('201_configs/201_olcum.npz',W,f'b2_{e}',f'L2_{e}'); R2=pool([r2[(w,1)][:2] for w in W]); R2m2=pool([r2[(w,2)][:2] for w in W])
        H[f'b2_eps{e}']=dict(R2=R2,R2_m2=R2m2,pencere={w:r2[(w,1)] for w in W})
    R2,s2=H['b2_eps0.2']['R2']; siniflar=[1,2,4]; karar='BELİRSİZ'
    for c in siniflar:
        if abs(R2-c)<=0.25*c and all(min(abs(R2-o*0.75),abs(R2-o*1.25)) >= 3*s2 or not (o*0.75<=R2<=o*1.25) for o in siniflar if o!=c):
            if all((R2<o*0.75-3*s2) or (R2>o*1.25+3*s2) for o in siniflar if o!=c): karar=f'KESİN: {c}'
    if karar=='BELİRSİZ' and not any(abs(R2-c)<=0.25*c for c in siniflar): karar='ARA DEĞER'
    H['H-201-2']=dict(R2=(R2,s2),karar=karar)
    b0=analiz('201_configs/201_olcum.npz',W,'b0','L1'); H['M1_b0_kontrol']={w:b0[(w,1)] for w in W}
    json.dump({k:(v if not isinstance(v,dict) else {str(kk):vv for kk,vv in v.items()}) for k,v in H.items()},open('201_configs/HUKUM_201.json','w'),indent=1,default=str)
    print('H-201-1: R1 = %.3f ± %.3f (m=2: %.3f ± %.3f); 1/L ölçekleme eğimi %.3f ± %.3f → %s / %s'%(R1[0],R1[1],R1m2[0],R1m2[1],egim,egim_se,H['H-201-1']['karar'],H['H-201-1']['olcekleme']))
    for w in W: print('   ',w,'L=%.2f R1=%.3f±%.3f'%(r1[(w,1)][2],r1[(w,1)][0],r1[(w,1)][1]))
    for e in (0.1,0.2,0.3): print('b2 ε=%s: R2 = %.3f ± %.3f (m=2: %.3f ± %.3f)'%(e,*H[f'b2_eps{e}']['R2'],*H[f'b2_eps{e}']['R2_m2']))
    print('H-201-2:',karar); print('b0 denetim R (≈0 olmalı, LG ölçeğinde):',{w:round(v[0],3) for w,v in H['M1_b0_kontrol'].items()})
