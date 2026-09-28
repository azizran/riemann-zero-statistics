# Tilt ladder: cumulants of log|Lambda| under |Lambda|^{2b}-tilted CUE(N-b), b=0,1,2 (N may be non-integer)
import mpmath as mp
mp.mp.dps=25
def Sr(r,n,a):  # sum_{j=1}^n psi^{(r-1)}(j+a), analytic continuation in n
    x=a+1
    if r==1: return (x+n-1)*mp.digamma(x+n)-(x-1)*mp.digamma(x)-n
    if r==2: return mp.digamma(x+n)+(x+n-1)*mp.psi(1,x+n)-mp.digamma(x)-(x-1)*mp.psi(1,x)
    if r==3: return 2*mp.psi(1,x+n)+(x+n-1)*mp.psi(2,x+n)-2*mp.psi(1,x)-(x-1)*mp.psi(2,x)
def kap(r,N,b):
    n=N-b; return Sr(r,n,2*b)-mp.mpf(2)**(1-r)*Sr(r,n,b)
A2,A3=-0.088124,0.233653  # converged (200t_aritmetik_kesin.py)
print('rung b | L | CUE k2  k3 | zeta-conj k2  k3 | skew CUE -> conj')
for b in (0,1,2):
  for L in (10,12,24.47):
    k2,k3=kap(2,L,b),kap(3,L,b); z2,z3=k2+A2,k3+A3
    print(b,'%6.2f | %.4f %+.4f | %.4f %+.4f | %+.3f -> %+.3f'%(L,k2,k3,z2,z3,k3/k2**1.5,z3/z2**1.5))
print('N->inf tilted(b=2) k3:',mp.nstr(kap(3,10**7,2),6),' b=1:',mp.nstr(kap(3,10**7,1),6))
print('Gaussian chi3 rival: k2 %.4f k3 %+.4f skew %+.3f'%(mp.psi(1,1.5)/4,mp.psi(2,1.5)/8,(mp.psi(2,1.5)/8)/(mp.psi(1,1.5)/4)**1.5))
