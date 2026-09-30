import sys, os, json, time, numpy as np
sys.dont_write_bytecode=True
REPO=os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),'..'))
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
os.chdir(REPO); sys.path.insert(0,REPO)
from importlib import import_module
from shifts import *
zeta=import_module('203_b1k3_zeta')
d=json.load(open(HERE+'/data_out.json')); Lbar=d['Lbar']
exec(open(HERE+'/formula_run.py').read().split("zf=(Z.S.G")[0].split("t=time.time()")[1].split("def moments")[1].join(["def moments",""]) if False else "")
# re-import helper defs from formula_run without running it
src=open(HERE+'/formula_run.py').read()
start=src.index('def moments'); end=src.index('zf=(Z.S.G')
exec(src[start:end])
out={'mean_check':{},'Lslope':{}}
for L in (Lbar,14.1576,14.2300):
    Z=zeta.Zeta(L); zf=(Z.S.G,Z.T,Z.Q,Z.Q2)
    if L==Lbar:
        for al in alpha_list():
            a=np.array([al/L]); out['mean_check'][str(al)]=[float(Z.S.G(a)[0]/L), d['K']['800']['full']['mean'][str(al)][0]]
    tr={}
    for lab,nom,act in TRIPLES:
        a,b,c=[x/L for x in act]; tr[lab]=cumulants(zf,L,a,b,c)
    out['Lslope'][str(L)]=tr
    print(L,tr,flush=True)
json.dump(out,open(HERE+'/extra_out.json','w'),indent=1)
print(json.dumps(out['mean_check'],indent=1))
