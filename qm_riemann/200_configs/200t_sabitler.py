# Large-N constants + prediction table for log Y (Y = small-gap limit of M/s~^2, unfolded)
import mpmath as mp
mp.mp.dps=30
def S1(n,a):  # sum_{j=1}^n psi(j+a), analytic in n
    x=a+1; return (x+n-1)*mp.digamma(x+n)-(x-1)*mp.digamma(x)-n
def S2(n,a):  # sum_{j=1}^n psi'(j+a)
    x=a+1; return mp.digamma(x+n)+(x+n-1)*mp.psi(1,x+n)-mp.digamma(x)-(x-1)*mp.psi(1,x)
def EY(N):  n=N-2; return S1(n,4)-S1(n,2)-mp.log(4)+2*mp.log(2*mp.pi/N)
def VY(N):  n=N-2; return S2(n,4)-S2(n,2)/2
# check vs integer sums
n=8; print('check', EY(10), mp.fsum(mp.digamma(j+4)-mp.digamma(j+2) for j in range(1,n+1))-mp.log(4)+2*mp.log(2*mp.pi/10))
print('E log Y (N->inf) = 2logpi+2gamma-10/3 =', mp.nstr(2*mp.log(mp.pi)+2*mp.euler-mp.mpf(10)/3,8), ' N=1e6:', mp.nstr(EY(10**6),8))
for N in (10**4,10**6,10**8): print(' c2(N)=VY-0.5logN', N, mp.nstr(VY(N)-mp.log(N)/2,10))
S=-0.17624781; A=S/2
chi3=mp.psi(1,1.5)/4
print('Gaussian rival Var log = psi1(3/2)/4 =',mp.nstr(chi3,6))
print(' L    | E_tilt(N=L)  E(L+1)  E(L+2) | V_tilt(L)  V(L+1)  V(L+2) | V_zeta_conj=V(L)+%.4f'%A)
for L in (6,8,10,11,12,16.58,18.89,21.19,23.49):
    print('%6.2f | %.4f %.4f %.4f | %.4f %.4f %.4f | %.4f'%(L,EY(L),EY(L+1),EY(L+2),VY(L),VY(L+1),VY(L+2),VY(L)+A))
