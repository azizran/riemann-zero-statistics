# Finite-L model spread: pure tilted CUE(L) | a_k x CUE(L) | hybrid(X) = random Euler product p<=X  (x)  tilted CUE(N_X), N_X = L/(e^gamma log X)
import mpmath as mp, importlib
from sympy import isprime
mp.mp.dps=25
_k=importlib.import_module("200t_merdiven").kap
kap=lambda r,N,b: float(_k(r,N,b))
def S12(x):
    s=mp.mpf(0); H=mp.mpf(0); m=1; t=x
    while True:
        m+=1; H+=mp.mpf(1)/(m-1); t*=x; term=H*t/m**2; s+=term
        if term<mp.mpf(10)**-25: return s
def euler(X):  # cumulants of sum_{p<=X} -log|1-p^{-1/2}e^{i a_p}|
    ps=[p for p in range(2,X+1) if isprime(p)]
    return float(mp.fsum(mp.polylog(2,mp.mpf(1)/p)/2 for p in ps)), float(mp.fsum(1.5*S12(mp.mpf(1)/p) for p in ps))
A2,A3=-0.088124,0.233653
for b in (0,1,2):
  for L in (10,12):
    row=[f'b={b} L={L}: CUE k2 {kap(2,L,b):.3f} k3 {kap(3,L,b):+.3f} | a_k k2 {kap(2,L,b)+A2:.3f} k3 {kap(3,L,b)+A3:+.3f}']
    for X in (2,3,5,7):
        NX=float(L/(mp.e**mp.euler*mp.log(X))); e2,e3=euler(X)
        if NX-b<=0.05: row.append(f'X={X}: N_X={NX:.2f} (yok)'); continue
        row.append(f'X={X} (N_X={float(NX):.2f}): k2 {kap(2,NX,b)+e2:.3f} k3 {kap(3,NX,b)+e3:+.3f}')
    print(' | '.join(row))
