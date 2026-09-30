"""203: oran sanısından κ₃(log|ζ(½+it)|), sonlu yükseklikte (türetim: 203_TURETIM notu).
κ₃ = c_P − (3/2) ∫_0^∞ dx₁ x₁ ∫_0^∞ dw Φ(x₁, x₁+w),   Φ(x₁,x₂) = F(x₁)K(x₂−x₁, x₁) + F(x₂)K(x₁−x₂, x₂)
  F(x) = e^{−Lx} ζ(1+x) ζ(1−x) A(x),  A(x) = Π_p (1−p^{−1−x})(1−2/p+p^{−1−x})/(1−1/p)²
  K(u,x) = −Σ_p log p · q/(1−q) · φ_p(x),  q = p^{−1−u},  φ_p = (1−1/p)(1−p^{−x})/(1−2/p+p^{−1−x})
         = (ζ'/ζ)(1+u) + P(1+u+x) + Q₂ − C   (analitik devam; w < 1'de yakınsak)
Kimlik terimi tam olarak c_P verir (FG'nin sabiti). CUE'de aynı çekirdek kesin κ₃(N)'i veriyor (N=1..10)."""
import numpy as np, mpmath as mp, sys
from scipy.integrate import quad
mp.mp.dps=20
def primes(n):
    s=np.ones(n+1,bool); s[:2]=False
    for i in range(2,int(n**.5)+1):
        if s[i]: s[i*i::i]=False
    return np.nonzero(s)[0].astype(float)
Pr=primes(1_000_000); lp=np.log(Pr); PM=Pr[-1]; lPM=np.log(PM)
cP=0.233653
def dz(s): s=mp.mpf(s); return float(mp.zeta(s,derivative=1)/mp.zeta(s))
def Pdir(s):   # s büyükse az asal yeter (kuyruk P^{1−s}/(s−1) ile)
    n=78498 if s<2.2 else (9592 if s<3.5 else 168); P0=Pr[n-1]
    return float(np.sum(lp[:n]*np.exp(-s*lp[:n])))+P0**(1-s)/(s-1)
def Pfun(s):
    if s>=1.6: return Pdir(s)
    return -dz(s)-sum(Pdir(m*s) for m in range(2,int(60/s)+2))
def logA(x):
    y=np.exp(-(1+x)*lp); return float(np.sum(np.log1p(-y)+np.log1p(-2/Pr+y)-2*np.log1p(-1/Pr)))
_F={}
def F(x,L):
    k=(x,L)
    if k not in _F:
        _F[k]=float(mp.exp(-L*x)*mp.zeta(1+mp.mpf(x))*mp.zeta(1-mp.mpf(x)))*np.exp(logA(x))
    return _F[k]
def K(u,x):
    q=np.exp(-(1+u)*lp); px=np.exp(-x*lp)
    s2=2+2*u+x; Q2=float(np.sum(lp*px*q*q/(1-q)))+PM**(1-s2)/(s2-1)
    psi1=(1/Pr-np.exp(-(1+x)*lp))/(1-2/Pr+np.exp(-(1+x)*lp))
    C=float(np.sum(lp*q/(1-q)*(1-px)*psi1))+PM**(-1-u)/(1+u)-2*PM**(-1-u-x)/(1+u+x)+PM**(-1-u-2*x)/(1+u+2*x)
    return dz(1+u)+Pfun(1+u+x)+Q2-C
def Phi(x1,w,L,Wc):
    t=F(x1,L)*K(w,x1)
    if w<Wc: t+=F(x1+w,L)*K(-w,x1+w)
    return t
def gl(bps,n):
    X,W=np.polynomial.legendre.leggauss(n); xs=[];ws=[]
    for a,b in zip(bps,bps[1:]):
        if b>a: xs+=list((b-a)/2*X+(a+b)/2); ws+=list((b-a)/2*W)
    return np.array(xs),np.array(ws)
def k3_gl(L,Wc=0.75,n=16):
    xb=sorted(set([1e-6,1e-4,1e-3,0.01,0.3/L,1/L,3/L,8/L,1.0,3.0])); X1,W1=gl(xb,n); tot=0.0
    for x1,a1 in zip(X1,W1):
        wb=sorted(set([0.0,x1/2,x1,2*x1,1/L,3/L,Wc,2.0,30.0])); Wv,Ww=gl(wb,n)
        tot+=a1*x1*sum(b*Phi(x1,w,L,Wc) for w,b in zip(Wv,Ww))
    return cP-1.5*tot
def k3(L,Wc=0.75):
    inner=lambda x1: x1*quad(lambda w: Phi(x1,w,L,Wc),1e-9,30,limit=200,points=[min(x1,Wc*0.99),Wc],epsabs=1e-7)[0]
    return cP-1.5*quad(inner,1e-7,3.0,limit=200,points=[0.3/L,1/L,4/L],epsabs=1e-7)[0]
if __name__=='__main__':
    Wc=float(sys.argv[1]); n=int(sys.argv[2]); L=float(sys.argv[3]); print(f'L={L:.4f} Wc={Wc} n={n} κ3_oran={k3_gl(L,Wc,n):.5f}',flush=True)
