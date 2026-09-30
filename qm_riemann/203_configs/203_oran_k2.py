"""203 (ısınma): oran sanısından (CFKRS/Conrey–Snaith 2007) log|ζ(½+it)| varyansı, sonlu yükseklikte.
Çift formülü (türetim bu dosyada, kendimiz; A'nın birinci türevleri P'de sıfır):
  E[(ζ'/ζ)(s+α)(ζ'/ζ)(1−s+β)] = G(α+β),
  G(x) = H(x) + e^{−Lx} ζ(1+x) ζ(1−x) A(x),   H(x) = (ζ'/ζ)'(1+x) − Σ_p log²p/(p^{1+x}−1)² = Σ_p log²p/(p^{1+x}−1),
  A(x) = Π_p (1 − p^{−1−x})(1 − 2/p + p^{−1−x}) / (1 − 1/p)².
log ζ(s) = −∫_0^∞ (ζ'/ζ)(s+a) da ve E[X X] = 0 (eşleniksiz ortalama) ⇒ κ₂(Re log ζ) = ½ ∫∫ G(a+b) da db = ½ ∫_0^∞ x G(x) dx."""
import numpy as np, mpmath as mp, sys
from scipy.integrate import quad
mp.mp.dps=30
def primes(n):
    s=np.ones(n+1,bool); s[:2]=False
    for i in range(2,int(n**.5)+1):
        if s[i]: s[i*i::i]=False
    return np.nonzero(s)[0].astype(float)
P=primes(2_000_000); lP=np.log(P); PM=P[-1]
def B(x):  # Σ log²p/(p^{1+x}−1)² + kuyruk
    return float(np.sum(lP**2/(P**(1+x)-1)**2) + np.log(PM)/((1+2*x)*PM**(1+2*x)))
def logA(x):
    y=P**(-1-x); return float(np.sum(np.log1p(-y)+np.log1p(-2/P+y)-2*np.log1p(-1/P)))
def G(x,L):
    s=mp.mpf(1)+x; z=mp.zeta(s); z1=mp.zeta(s,derivative=1); z2=mp.zeta(s,derivative=2)
    dlog2=z2/z-(z1/z)**2
    sw=mp.exp(-L*x)*z*mp.zeta(1-mp.mpf(x))*mp.exp(logA(x))
    return float(dlog2+sw)-B(x)
def k2(L):
    f=lambda x: x*G(x,L)
    return 0.5*sum(quad(f,a,b,limit=200,epsabs=1e-10)[0] for a,b in ((1e-12,1/L),(1/L,1.0),(1.0,4.0),(4.0,14.0)))
if __name__=='__main__':
    for L in map(float,sys.argv[1:]): print(f'L={L:.4f}  κ2_oran={k2(L):.5f}')
