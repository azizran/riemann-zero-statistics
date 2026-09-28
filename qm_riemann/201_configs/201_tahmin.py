# KALEM 201 — Landau–Gonek faz kilidi öngörüleri (yalnız formül; veri yok)
# Sıfırlarda φ_p = γ log p (mod 2π): E cos(mφ) = −log p / (p^{m/2} L) (b=1, önde gelen mertebe); çift için H_ÇİFT: 2×.
# Kilidin asal-izine (κ₃ kayması) etkisi: W_p = −log|1 − p^{-1/2} e^{-iφ}|, eğilmiş yoğunlukla kümülantlar (asallar bağımsız varsayımı).
import numpy as np, json
from sympy import primerange
phi=np.linspace(0,2*np.pi,4096,endpoint=False)
def cum(W,w):
    w=w/w.sum(); m=(W*w).sum(); d=W-m; return m,(d*d*w).sum(),(d**3*w).sum()
def kilit_etkisi(L,kat,P=2000):
    d1=d2=d3=0.0
    for p in primerange(2,P):
        W=-np.log(np.abs(1-p**-0.5*np.exp(-1j*phi)))
        f=np.ones_like(phi)
        for m in range(1,8): f+=2*kat*(-np.log(p)/(p**(m/2)*L))*np.cos(m*phi)
        a=cum(W,np.ones_like(phi)); b=cum(W,f)
        d1+=b[0]-a[0]; d2+=b[1]-a[1]; d3+=b[2]-a[2]
    return d1,d2,d3
out={}
print('Öngörülen faz harmonikleri (b=1), c1=E cos φ, c2=E cos 2φ:')
for L in (9.35,10.59,11.66,14.19,16.58,18.88,22.31):
    row={p:(-np.log(p)/(np.sqrt(p)*L), -np.log(p)/(p*L)) for p in (2,3,5,7,11,13)}
    out[str(L)]={'harmonik_b1':{str(p):v for p,v in row.items()}}
    print(f' L={L}: '+' '.join(f'p{p}: {c1:+.4f}/{c2:+.4f}' for p,(c1,c2) in row.items()))
print('Kilidin kümülant kaymaları (asallar ≤ 2000; tek kilit b=1 / çift kilit b=2):')
for L in (10.59,22.31):
    e1=kilit_etkisi(L,1); e2=kilit_etkisi(L,2)
    out[str(L)]['kilit_b1']=e1; out[str(L)]['kilit_b2']=e2
    print(f' L={L}: b1 Δκ1 {e1[0]:+.4f} Δκ2 {e1[1]:+.4f} Δκ3 {e1[2]:+.4f} | b2 Δκ1 {e2[0]:+.4f} Δκ2 {e2[1]:+.4f} Δκ3 {e2[2]:+.4f}')
json.dump(out,open('201_configs/201_tahmin.json','w'),indent=1)
