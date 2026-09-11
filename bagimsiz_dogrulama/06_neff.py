"""Bagimsiz N_eff: zeta'nin (Pearson/Spearman) korelasyonunu CUE egrisine ters cevir."""
import numpy as np, json, importlib.util
spec=importlib.util.spec_from_file_location("o","02_olcum.py"); o=importlib.util.module_from_spec(spec); spec.loader.exec_module(o)
TWO_PI=o.TWO_PI; A_CG=o.A_CG
cue=json.load(open("05_cue_sonuc_20000.json")); Ns=np.array(sorted(int(k) for k in cue)); 
rP=np.array([cue[str(n)]["r_pearson"] for n in Ns]); rS=np.array([cue[str(n)]["r_spearman"] for n in Ns])

dz=np.load("../qm_riemann/128_odl_zeros6_2e6_zeros.npz"); zeros=np.sort(np.asarray(dz["zeros"],float))

# makalenin tablosu (arxiv_gap_amplitude.tex)
tablo=[(5.60,1029,0.8917,0.75,1.70),(6.30,2322,0.8735,0.76,1.80),(6.99,5174,0.8566,0.77,1.81),
       (7.69,11418,0.8394,0.87,1.93),(8.39,24988,0.8238,0.95,1.98),(9.08,54299,0.8096,0.99,2.00),
       (9.86,39992,0.7947,1.04,2.00),(10.37,39997,0.7849,1.13,2.01),(10.93,39998,0.7754,1.20,2.01),
       (11.47,39998,0.7663,1.31,2.14),(11.98,39999,0.7589,1.31,2.07),(12.45,40000,0.7521,1.30,2.10)]
def ters(r, r_egri):
    i=np.argsort(r_egri); rr=r_egri[i]; nn=Ns[i]
    if r>rr.max() or r<rr.min(): return np.nan
    return float(np.interp(r, rr, nn))
print(f"{'L':>6} {'n':>7} | {'r_P ben':>8} {'r_P mak':>8} | {'NeffP-L ben':>11} {'mak':>6} | {'NeffS-L ben':>11} {'mak':>6}")
for (L,n,r_mak,dP,dS) in tablo:
    tc=2*np.pi*np.exp(L); s=o.olc(zeros,tc,n)
    g=o.pencere(zeros,tc,n); mid=0.5*(g[:-1]+g[1:]); Lm=np.log(np.mean(mid)/TWO_PI)
    d=np.diff(g); M=np.array([o.M_aralik(g[i],g[i+1]) for i in range(len(d))])
    dt=d*Lm/TWO_PI; Mt=M/np.sqrt(A_CG*Lm+2.758-0.054/Lm)
    rPz=float(np.corrcoef(dt,Mt)[0,1])
    rg=np.argsort(np.argsort(dt)); rm=np.argsort(np.argsort(Mt)); rSz=float(np.corrcoef(rg,rm)[0,1])
    nP=ters(rPz,rP); nS=ters(rSz,rS)
    print(f"{Lm:6.2f} {len(d):7d} | {rPz:8.4f} {r_mak:8.4f} | {nP-Lm:11.2f} {dP:6.2f} | {nS-Lm:11.2f} {dS:6.2f}")
print("\n(r_P mak = makaledeki Pearson; NeffP/NeffS = CUE etkin boyut kaymasi)")
