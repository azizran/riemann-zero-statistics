import sys, os, json, time, numpy as np
sys.dont_write_bytecode=True
REPO=os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),'..'))
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
os.chdir(REPO); sys.path.insert(0,REPO)
from importlib import import_module
from shifts import *
zeta=import_module('203_b1k3_zeta'); cue=import_module('203_b1k3_cue')
Lbar=json.load(open(HERE+'/data_out.json'))['Lbar'] if len(sys.argv)<2 else float(sys.argv[1])
t=time.time(); Z=zeta.Zeta(Lbar); print('Zeta built',time.time()-t,flush=True)
def moments(funcs, L):
    """funcs: dict G,T,Q,Q2 (array args, raw shifts). returns closures for E_z moments"""
    G,T,Q,Q2=funcs
    arr=lambda x: np.atleast_1d(np.asarray(x,float))
    def mz1(a): a=arr(a); return G(a)/L
    def mzXX(a,b): a=arr(a);b=arr(b);z=0*a; return T(a,b,z)/L
    def mzXXb(a,c): a=arr(a);c=arr(c);z=0*a; return G(a+c)+(T(a,z,c)+T(c,z,a))/L
    def mz3XXb(a,b,c):
        a=arr(a);b=arr(b);c=arr(c);z=0*a
        return T(a,b,c)+(Q([a,b,z],c)+Q2(a,b,c,z))/L
    def mz3XXX(a,b,c):
        a=arr(a);b=arr(b);c=arr(c);z=0*a
        return Q([a,b,c],z)/L
    return mz1,mzXX,mzXXb,mz3XXb,mz3XXX
def cumulants(funcs,L,a,b,c):
    m1,mXX,mXXb,m3b,m3x=moments(funcs,L)
    kb=m3b(a,b,c)-mXX(a,b)*m1(c)-mXXb(a,c)*m1(b)-mXXb(b,c)*m1(a)+2*m1(a)*m1(b)*m1(c)
    kx=m3x(a,b,c)-mXX(a,b)*m1(c)-mXX(a,c)*m1(b)-mXX(b,c)*m1(a)+2*m1(a)*m1(b)*m1(c)
    return float(kb[0]),float(kx[0])
def cum2(funcs,L,a,b):
    m1,mXX,mXXb,_,_=moments(funcs,L)
    return float((mXXb(a,b)-m1(a)*m1(b))[0]), float((mXX(a,b)-m1(a)*m1(b))[0])
zf=(Z.S.G,Z.T,Z.Q,Z.Q2); _G,_T,_Q,_Q2=cue.funcs(Lbar); cf=(_G,_T,(lambda A,b:_Q(A[0],A[1],A[2],b)),_Q2)
res={'Lbar':Lbar,'triple':{},'pair_XXB':{},'pair_XX':{},'delta_scan':{}}
for lab,nom,act in TRIPLES:
    a,b,c=[x/Lbar for x in act]
    res['triple'][lab]={'zeta':cumulants(zf,Lbar,a,b,c),'cue':cumulants(cf,Lbar,a,b,c)}
for (a,b) in PAIRS_XXB:
    A,B=a/Lbar,b/Lbar
    res['pair_XXB'][f'{a},{b}']={'zeta':cum2(zf,Lbar,A,B)[0],'cue':cum2(cf,Lbar,A,B)[0]}
for (a,b) in PAIRS_XX:
    A,B=a/Lbar,b/Lbar
    res['pair_XX'][f'{a},{b}']={'zeta':cum2(zf,Lbar,A,B)[1],'cue':cum2(cf,Lbar,A,B)[1]}
json.dump(res,open(HERE+'/formula_out.json','w'),indent=1)
# delta-scan for coincident triples (noise / smoothness of formula near coincident shifts)
for lab,nom,act in TRIPLES:
    if len(set(nom))<3:
        sc={}
        for dd in (0.005,0.01,0.02,0.04,0.08):
            al=list(nom)
            # perturb the coincident pair(s) as in shifts.py: first two symmetric
            if lab=='(1,1,1)': tr=(1-dd,1+dd,1.0)
            elif lab=='(3,3,2)': tr=(3*(1-dd),3*(1+dd),2.0)
            else: tr=(2*(1-dd),2*(1+dd),4.0)
            a,b,c=[x/Lbar for x in tr]
            sc[dd]={'zeta':cumulants(zf,Lbar,a,b,c),'cue':cumulants(cf,Lbar,a,b,c)}
        res['delta_scan'][lab]=sc
json.dump(res,open(HERE+'/formula_out.json','w'),indent=1)
print(json.dumps(res['triple'],indent=1)); print(json.dumps(res['delta_scan'],indent=1))
