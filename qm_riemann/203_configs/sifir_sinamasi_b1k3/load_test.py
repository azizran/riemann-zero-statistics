import sys, time, numpy as np
sys.dont_write_bytecode=True
import importlib.util, os
QM=os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','..'))   # qm_riemann/
spec=importlib.util.spec_from_file_location('veri',os.path.join(QM,'200_configs','200c_veri.py'))
v=importlib.util.module_from_spec(spec); spec.loader.exec_module(v)
t=time.time()
taban,off,bl=v.oku(os.path.join(QM,'veri_lmfdb','zeros_8846000.dat'),n_max=1_500_000)
print(time.time()-t,taban,len(off),off[:3],off[-3:])
g=np.diff(off); print(g.min(),g.mean(),g.max(), np.all(g>0))
L=np.log((taban+off)/(2*np.pi)); print(L[0],L[-1],L.mean())
