# Arithmetic cumulants of log|zeta| from a_k (Keating-Snaith): kappa_r = (d/dk)^r log a_k |_0 / 2^r
import mpmath as mp
from sympy import primerange
mp.mp.dps=25
def loga(k,P=200000):
    k=mp.mpf(k); s=mp.mpf(0)
    for p in primerange(2,P):
        x=mp.mpf(1)/p; s+=k*k*mp.log(1-x)+mp.log(mp.hyp2f1(k,k,1,x))
    return s
# numerical derivatives (small-prime-dominated; tail beyond P is O(k^2/P))
h=mp.mpf('0.02')
f={i:loga(i*h,20000) for i in (-3,-2,-1,0,1,2,3)}
d2=(f[1]-2*f[0]+f[-1])/h**2
d3=(f[2]-2*f[1]+2*f[-1]-f[-2])/(2*h**3)
print('kappa2_arith = %.5f  kappa3_arith = %.5f'%(d2/4,d3/8))
# direct check: independent Euler factors, -log|1-p^{-1/2}e^{ia}| cumulants (+ (1/2)log(1-1/p) in kappa2)
k2=k3=mp.mpf(0)
for p in primerange(2,20000):
    r=1/mp.sqrt(p); g=lambda a: -mp.log(abs(1-r*mp.expj(a)))
    m2=mp.quad(lambda a:g(a)**2,[0,mp.pi,2*mp.pi])/(2*mp.pi) if p<200 else mp.polylog(2,mp.mpf(1)/p)/2
    m3=mp.quad(lambda a:g(a)**3,[0,mp.pi,2*mp.pi])/(2*mp.pi) if p<200 else 0
    k2+=m2+mp.log(1-mp.mpf(1)/p)/2; k3+=m3
print('direct: kappa2 = %.5f  kappa3(p<200) = %.5f'%(k2,k3))
