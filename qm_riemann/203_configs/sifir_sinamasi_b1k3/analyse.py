import json, numpy as np, sys
sys.dont_write_bytecode=True
from shifts import *
d=json.load(open('data_out.json')); f=json.load(open('formula_out.json'))
Ks=[int(k) for k in d['K'].keys()]
def val(K,sub,kind,lab,typ=None,idx=0):
    r=d['K'][str(K)][sub][kind][lab]
    if typ: r=r[typ]
    return r[idx]
def extrap(sub,kind,lab,typ=None,idx=0,Kset=(200,400,800),order=2):
    K=np.array(Kset,float); y=np.array([val(int(k),sub,kind,lab,typ,idx) for k in Kset])
    A=np.vander(1/K,order+1,increasing=True)
    c=np.linalg.solve(A,y); return c[0]
def ext_all(sub,kind,lab,typ=None,idx=0):
    e3=extrap(sub,kind,lab,typ,idx,(200,400,800),2)
    e2=extrap(sub,kind,lab,typ,idx,(400,800),1)
    e3b=extrap(sub,kind,lab,typ,idx,(100,200,400),2)
    e4=extrap(sub,kind,lab,typ,idx,(100,200,400,800),3)
    return e3,e2,e3b,e4
L=[]
P=lambda s='': L.append(s)
Lbar=d['Lbar']
P(f'Lbar = {Lbar:.6f}; data: zeros_8846000.dat, first 1.5e6 zeros (window C1); evaluated zeros idx 800..1499168 (n=1498368 used); 32-block jackknife')
out={}
def row(name,dat,se,ft,ct,extra=''):
    return f'| {name} | {dat:+.4f} ± {se:.4f} | {ft:+.4f} | {(dat-ft)/se:+.1f} | {ct:+.4f} | {(dat-ct)/se:+.1f} |{extra}'
P('\n## W-stability (window = ±K zeros, K=50..800; mean spacing 2π/L=0.4426; W=K·0.4426) – cumulant vs K')
P('K-dependence ~ c/K (systematic). Richardson (200,400,800; quadratic in 1/K) used as K→∞ estimate.')
P('| cumulant | K=50 | K=100 | K=200 | K=400 | K=800 | extrap(200,400,800) | extrap(400,800; lin) | extrap(100,200,400) | extrap(all4; cubic) |')
P('|---|---|---|---|---|---|---|---|---|---|')
allrows=[]
for lab,nom,act in TRIPLES:
    for typ in ('XXB','XXX'):
        vals=[val(K,'full','triple',lab,typ,0) for K in (50,100,200,400,800)]
        e=ext_all('full','triple',lab,typ,0)
        P(f'| {lab} {typ} | '+' | '.join(f'{v:.4f}' for v in vals)+' | '+' | '.join(f'{x:.4f}' for x in e)+' |')
for kind,pairs,typ in (('pair_XXB',PAIRS_XXB,'XXB'),('pair_XX',PAIRS_XX,'XX')):
    for (a,b) in pairs:
        lab=f'{a},{b}'
        vals=[val(K,'full',kind,lab,None,0) for K in (50,100,200,400,800)]
        e=ext_all('full',kind,lab,None,0)
        P(f'| {typ}2 ({a},{b}) | '+' | '.join(f'{v:.4f}' for v in vals)+' | '+' | '.join(f'{x:.4f}' for x in e)+' |')

P('\n## Main table: third-order joint zero-cumulants, real part (data = K→∞ extrapolated, SE = jackknife at K=800)')
P('Actual shifts used (units 1/L) in brackets; coincident nominal shifts perturbed by ±2% (formula singular at exact coincidence), data and formula evaluated at the SAME shifts.')
P('| triple, type | data ± SE | formula | pull_F | CUE(N=L) | pull_CUE |')
P('|---|---|---|---|---|---|')
res3={}
for lab,nom,act in TRIPLES:
    for typ in ('XXB','XXX'):
        dat=extrap('full','triple',lab,typ,0); se=val(800,'full','triple',lab,typ,1)
        ft=f['triple'][lab]['zeta'][0 if typ=='XXB' else 1]; ct=f['triple'][lab]['cue'][0 if typ=='XXB' else 1]
        im=val(800,'full','triple',lab,typ,2); ime=val(800,'full','triple',lab,typ,3)
        nm=f'{lab} [{",".join(f"{x:g}" for x in act)}] {"X X X̄" if typ=="XXB" else "X X X"}'
        P(row(nm,dat,se,ft,ct)); res3[(lab,typ)]=(dat,se,ft,ct,im,ime)
P('\n## Imaginary parts of the same third cumulants (K=800; expected ≈ 0 by symmetry)')
P('| triple, type | Im part ± SE | Im/SE |'); P('|---|---|---|')
for (lab,typ),(dat,se,ft,ct,im,ime) in res3.items():
    P(f'| {lab} {typ} | {im:+.4f} ± {ime:.4f} | {im/ime:+.1f} |')
P('\n## Second-order warm-up (real part; data K→∞ extrapolated; SE at K=800)')
P('| pair, type | data ± SE | formula | pull_F | CUE(N=L) | pull_CUE |'); P('|---|---|---|---|---|---|')
res2={}
for kind,pairs,typ in (('pair_XXB',PAIRS_XXB,'XXB'),('pair_XX',PAIRS_XX,'XX')):
    for (a,b) in pairs:
        lab=f'{a},{b}'; dat=extrap('full',kind,lab,None,0); se=val(800,'full',kind,lab,None,1)
        ft=f[kind][lab]['zeta']; ct=f[kind][lab]['cue']
        P(row(f'{"X X̄" if typ=="XXB" else "X X"} ({a},{b})',dat,se,ft,ct)); res2[(kind,lab)]=(dat,se,ft,ct)
# subsets
P('\n## Subset consistency (L-mixing / statistical): halves and middle 50% of the window, extrapolated data')
P('| item | full | half1 | half2 | mid50 | formula |'); P('|---|---|---|---|---|---|')
for lab,nom,act in TRIPLES:
    for typ in ('XXB','XXX'):
        vs=[extrap(s,'triple',lab,typ,0) for s in ('full','half1','half2','mid')]
        ses=[val(800,s,'triple',lab,typ,1) for s in ('full','half1','half2','mid')]
        ft=f['triple'][lab]['zeta'][0 if typ=='XXB' else 1]
        P(f'| {lab} {typ} | '+' | '.join(f'{v:.3f}±{s:.3f}' for v,s in zip(vs,ses))+f' | {ft:.3f} |')
P('\n## Jackknife block-count check of SE (K=800): NB=32 / 64 / 128')
P('| item | SE32 | SE64 | SE128 |'); P('|---|---|---|---|')
for lab,nom,act in TRIPLES:
    for typ in ('XXB','XXX'):
        s=[val(800,k,'triple',lab,typ,1) for k in ('full','full_NB64','full_NB128')]
        P(f'| {lab} {typ} | '+' | '.join(f'{x:.4f}' for x in s)+' |')
P('\n## Formula near-coincidence stability (formula value vs perturbation δ of coincident shifts; XXB / XXX)')
for lab,sc in f['delta_scan'].items():
    P(lab+': '+'; '.join(f'δ={k}: zeta {v["zeta"][0]:.4f}/{v["zeta"][1]:.4f}' for k,v in sc.items()))
open('report_tables.md','w').write('\n'.join(L)); print('\n'.join(L))
