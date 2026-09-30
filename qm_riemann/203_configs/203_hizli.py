"""203-b1κ₃ adım 4: hızlı (vektörlü) ζ yapı taşları — tablolar + kesik asal toplamları; yavaş/kesin sürümlerle sınanır.
K(u,x) = [r(u) − 1/u] + P(1+u+x) + R(u,x),  r(u) = (ζ'/ζ)(1+u) + 1/u,  P(s) = Σ log p p^{−s},  R = Q₂ − C (203_oran_k3)."""
import numpy as np, mpmath as mp, sys
from scipy.interpolate import CubicSpline, RectBivariateSpline
sys.path.insert(0,'.'); from importlib import import_module
k3=import_module('203_oran_k3'); k2=import_module('203_oran_k2')
Pr=k3.Pr; lp=k3.lp
Ps=Pr[Pr<1e4]; lps=np.log(Ps)
EG=float(mp.euler)
def _logx_grid(lo,hi,n): return np.exp(np.linspace(np.log(lo),np.log(hi),n))
# --- 1B tablolar (L'den bağımsız) ---
_u=np.concatenate([np.linspace(-0.78,-1e-3,500),np.linspace(1e-3,3,900),_logx_grid(3.02,130,250)])
_rs=CubicSpline(_u,np.array([float(mp.zeta(1+mp.mpf(u),derivative=1)/mp.zeta(1+mp.mpf(u)))+1/u for u in _u]))
def r_u(u): return _rs(u)          # r(u) = (ζ'/ζ)(1+u) + 1/u, analitik (|u|<1e-3 aralığı spline ile köprülenir)
_s=_logx_grid(1e-5,200,1600)       # s−1
_Ps=CubicSpline(np.log(_s),np.array([k3.Pfun(1+x)-1/x for x in _s]))
def Pfun(s): sm=s-1; return _Ps(np.log(sm))+1/sm
# R(u,x) = Q₂ − C; geçerli bölge x ≥ max(0,−u). İki tablo: u ≥ 0 için (u,x); u < 0 için (w=−u, y=x−w ≥ 0).
def _Rterms(u,x,P,l):   # asal kümesi P için Q₂ − C (kuyruksuz), vektörlü
    q=np.exp(-(1+u)*l); px=np.exp(-x*l)
    Q2=l*px*q*q/(1-q)
    psi1=(1/P-np.exp(-(1+x)*l))/(1-2/P+np.exp(-(1+x)*l))
    return Q2-l*q/(1-q)*(1-px)*psi1
def _R_direct(u,x):     # p ≥ 3 (p = 2 ayrıca, kesin: 1−2/p+p^{−1−x} → 0 ⇒ 2^x büyümesi tabloya girmez)
    s2=2+2*u+x
    tail=k3.PM**(1-s2)/(s2-1)-(k3.PM**(-1-u)/(1+u)-2*k3.PM**(-1-u-x)/(1+u+x)+k3.PM**(-1-u-2*x)/(1+u+2*x))
    return float(np.sum(_Rterms(u,x,Pr[1:],lp[1:])))+tail
def R2(u,x): return _Rterms(np.asarray(u,float),np.asarray(x,float),2.0,np.log(2.0))
_gx=np.concatenate([[0.0],_logx_grid(1e-5,0.1,30),np.linspace(0.12,2.0,60),_logx_grid(2.1,60,30)])
_gu=np.concatenate([np.linspace(0,2.0,70),_logx_grid(2.1,60,30)])
_gw=np.linspace(0,0.78,50)
_Rp=RectBivariateSpline(_gu,_gx,np.array([[_R_direct(u,x) for x in _gx] for u in _gu]),kx=3,ky=3)
_Rm=RectBivariateSpline(_gw,_gx,np.array([[_R_direct(-w,y+w) for y in _gx] for w in _gw]),kx=3,ky=3)
def R_ux(u,x):
    u=np.asarray(u,float); x=np.asarray(x,float)
    return np.where(u>=0,_Rp.ev(np.maximum(u,0),x),_Rm.ev(np.minimum(-u,0.78),np.maximum(x+u,0)))
def K(u,x):
    u=np.clip(u,-0.78,120.0); x=np.clip(x,0.0,60.0)   # tablo aralıkları (dışarıda çarpan F zaten ~0)
    return r_u(u)-1/u+Pfun(1+u+x)+R_ux(u,x)+R2(u,x)
class Lset:
    """L'ye bağlı tablolar: F(x), H(x), G(x)."""
    def __init__(s,L):
        s.L=L; X=_logx_grid(1e-7,130,1400); s.X=X
        Fv=np.array([k3.F(x,L) for x in X]); Hv=np.array([k2.G(x,L)-k3.F(x,L) for x in X])
        s._F=CubicSpline(np.log(X),X**2*Fv); s._H=CubicSpline(np.log(X),X**2*Hv)
    def F(s,x): return np.where(x<3.0,s._F(np.log(np.minimum(x,3.0)))/np.minimum(x,3.0)**2,0.0)   # x ≥ 3: |F| ≤ e^{−3L}·O(1) ⇒ 0
    def H(s,x): return s._H(np.log(x))/x**2
    def G(s,x): return s.F(x)+s.H(x)
    def Aexp(s,x):   # A(x) = Π_p (1−p^{−1−x})(1−2/p+p^{−1−x})/(1−1/p)², tablo
        if not hasattr(s,'_A'):
            Xa=np.concatenate([[0.0],_logx_grid(1e-6,60,600)]); s._A=CubicSpline(Xa,np.array([np.exp(k3.logA(x)) for x in Xa]))
        return s._A(np.clip(x,0,60))
if __name__=='__main__':
    rng=np.random.default_rng(5)
    print('K hızlı vs kesin:')
    err=[]
    for _ in range(80):
        if rng.random()<0.5: u=rng.uniform(-0.75,0); x=-u+10**rng.uniform(-4,0.5)
        else: u=10**rng.uniform(-3,1.2); x=10**rng.uniform(-4,1.2)
        f=float(K(np.array(u),np.array(x))); e=k3.K(u,x); err.append(abs(f-e)/max(1,abs(e)))
        if err[-1]>1e-6: print(f'  büyük hata u={u:.4g} x={x:.4g}: hızlı {f:+.9f} kesin {e:+.9f}')
    print('  K göreli hata: maks %.2e medyan %.2e (80 rastgele nokta)'%(max(err),np.median(err)))
    S=Lset(11.66)
    for x in (1e-5,3e-3,0.08,0.7,5.0):
        print(f'  x={x}: F {float(S.F(np.array(x))):+.9e} kesin {k3.F(x,11.66):+.9e} | G {float(S.G(np.array(x))):+.9e} kesin {k2.G(x,11.66):+.9e}')
