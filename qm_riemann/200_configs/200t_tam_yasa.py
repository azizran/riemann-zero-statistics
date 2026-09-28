# Full predicted law of W = log X (small-gap hump, theta-units part) under H_R (tilted CUE) and H_A (x a_k):
# phi_R(w) = M_n(2+iw/2)/M_n(2) (Barnes G, analytic in n); phi_A = phi_R * a_{iw/2}; density by Fourier inversion.
import mpmath as mp, numpy as np, sys
from sympy import primerange
mp.mp.dps=20
N=float(sys.argv[1]) if len(sys.argv)>1 else 11.0; n=N-2
G=mp.barnesg
def logM(n,k): return 2*mp.log(G(1+k))+mp.log(G(n+1))+mp.log(G(n+1+2*k))-mp.log(G(1+2*k))-2*mp.log(G(n+1+k))
PR=list(primerange(2,20000))
def loga(k): return mp.fsum(k*k*mp.log(1-mp.mpf(1)/p)+mp.log(mp.hyp2f1(k,k,1,mp.mpf(1)/p)) for p in PR)
W=np.linspace(-30,30,1201); dw=W[1]-W[0]
L0=logM(n,2)
phR=np.array([complex(mp.e**(logM(n,2+1j*w/2)-L0)) for w in W])
phA=phR*np.array([complex(mp.e**loga(1j*w/2)) for w in W])
x=np.linspace(-1.5,4.5,1201)
def dens(ph): return np.real(np.exp(-1j*np.outer(x,W))@ph)*dw/(2*np.pi)
for name,ph in (('H_R',phR),('H_A',phA)):
    f=dens(ph); F=np.cumsum(f)*(x[1]-x[0])
    q=lambda a: np.interp(a,F,x)
    m=np.trapz(x*f,x); v=np.trapz((x-m)**2*f,x); k3=np.trapz((x-m)**3*f,x)
    kelly=(q(.9)+q(.1)-2*q(.5))/(q(.9)-q(.1)); bowley=(q(.75)+q(.25)-2*q(.5))/(q(.75)-q(.25))
    print(f'N={N} {name}: mass {F[-1]:.5f} min f {f.min():.2e} | mean {m:.4f} var {v:.4f} k3 {k3:+.4f} | Kelly {kelly:+.4f} Bowley {bowley:+.4f} | q10 {q(.1):.3f} q50 {q(.5):.3f} q90 {q(.9):.3f}')
