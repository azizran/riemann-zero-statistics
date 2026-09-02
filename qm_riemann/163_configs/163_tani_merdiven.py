"""163 probe2 — sentetik gazin GERCEKLEsen merdiven genligi nominalin ne kadari?"""
import importlib, sys
from pathlib import Path
import numpy as np
QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
sys.path.insert(0, str(QM/"160_configs")); sys.path.insert(0, str(QM/"159_configs"))
C160 = importlib.import_module("160_cekirdek")
TWO_PI = 2*np.pi
for veri in ("keskin","A4","son"):
    z = C160.KOS155.veri_yukle(veri)
    Y = C160.Yerel160(z, veri, 0.40, 4000)
    m = Y.m0; L = Y.L
    ds_full = Y.ds                      # ds_n at mid_n
    dsm = ds_full - ds_full.mean()
    print(f"=== {veri}  L={L:.4f} Var(ds)={np.var(ds_full):.5f}")
    qm = C160.C154.pk_m(int(np.exp(1.0*L)))
    print("   q   tau      b_nom     |ds-line|   oran    arg")
    for q in (2,3,4,5,7,8,9,11,16,101,1009):
        if q not in qm: continue
        w = np.log(q); tau = w/L
        a = 1.0/(np.pi*qm[q]*np.sqrt(q)); b = 2*a*np.sin(np.pi*tau)
        c = 2*np.mean(dsm*np.exp(-1j*w*Y.mid))
        print(f"{q:5d} {tau:.4f}  {b:9.5f}  {abs(c):9.5f}  {abs(c)/b:6.3f}  {np.angle(c):+7.4f}")
