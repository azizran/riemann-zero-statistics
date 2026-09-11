# -*- coding: utf-8 -*-
"""
(1) K-FAKTORU: acik formul genligi 1/(pi k sqrt(q)); benim normum sqrt(q) idi.
    u x k birime donuyor mu?
(2) TARAK KOPRUSU (188): unfold kafeste omega=2*pi*k Bragg yansimasi.
    |G(2pi k)| olc; Gauss/DW tahmini exp(-2 pi^2 sigma_phi^2) ile oranla.
"""
import numpy as np, importlib.util, json
def yukle(ad,yol):
    s=importlib.util.spec_from_file_location(ad,yol); m=importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
o=yukle("o","02_olcum.py"); s7=yukle("s7","07_wv.py"); TWO_PI=2*np.pi
zeros=np.sort(np.asarray(np.load("../qm_riemann/128_odl_zeros6_2e6_zeros.npz")["zeros"],float))

print("=== (1) K-FAKTORU: kosulsuz icerik u(q) x k ===")
print(f"{'q':>4} {'p':>3} {'k':>2} {'u':>8} {'u*k':>8}")
tum=[]
for mc in [1.2e5,1.0e6]:
    d,M,tmid=s7.pencere_verisi(zeros,mc,40000)
    L0=float(np.log(np.mean(tmid)/o.TWO_PI)); Q=s7.taban_olustur(L0)
    _,u,_,_,_=s7.regresyon(d,M,tmid,Q,gap_kontrol=False)
    print(f"  -- L={L0:.3f}")
    for (q,p) in Q:
        if q>30: continue
        k=int(round(np.log(q)/np.log(p)))
        print(f"{q:4d} {p:3d} {k:2d} {u[q]:8.4f} {u[q]*k:8.4f}")
    tum.append((L0,u))
kuvvet=[(q,u[q]*int(round(np.log(q)/np.log(p)))) for (L0,u) in tum for (q,p) in s7.taban_olustur(L0) if q>1]
print(f"\n  tüm asal-kuvvetlerinde u*k ortalaması = {np.mean([x[1] for x in kuvvet]):.4f} ± {np.std([x[1] for x in kuvvet]):.4f}  (birim beklenir)")

print("\n=== (2) TARAK KOPRUSU: Bragg yansıması omega=2pi k (unfold) ===")
print(f"{'L':>7} {'|G(2pi)|':>9} {'sigma_phi':>10} {'Gauss e^-c':>11} {'oran':>7} | {'|G(4pi)|':>9} {'|G(6pi)|':>9}")
tarak=[]
for mc in [1.2e5,2.0e5,3.5e5,6.0e5,1.0e6]:
    d,M,tmid=s7.pencere_verisi(zeros,mc,40000)
    L=float(np.log(np.mean(tmid)/TWO_PI))
    x=tmid*L/TWO_PI                      # unfold: ortalama aralik 1
    phi=x-np.floor(x); phi=phi-phi.mean()
    s_phi=float(phi.std())
    def G(w): return abs(np.exp(1j*w*x).mean())
    g1,g2,g3=G(TWO_PI),G(2*TWO_PI),G(3*TWO_PI)
    gauss=float(np.exp(-0.5*(TWO_PI*s_phi)**2))     # DW: c=(2pi sigma)^2/2
    print(f"{L:7.3f} {g1:9.5f} {s_phi:10.5f} {gauss:11.5f} {g1/gauss:7.3f} | {g2:9.5f} {g3:9.5f}")
    tarak.append(dict(L=L, G1=g1, G2=g2, G3=g3, sigma_phi=s_phi, gauss=gauss, oran=g1/gauss))
print(f"\n  oran ortalaması = {np.mean([t['oran'] for t in tarak]):.4f} ± {np.std([t['oran'] for t in tarak]):.4f}   (makale: 0.88-0.89 sabit)")
print(f"  ikinci tarak |G(4pi)| ort = {np.mean([t['G2'] for t in tarak]):.5f} (gürültü tabanı?)")
json.dump(dict(tarak=tarak), open("15_tarak.json","w"), indent=1)
