# Converged arithmetic cumulants (audit fix): kappa2 = 1/2 sum_p [log(1-1/p)+Li2(1/p)], kappa3 = 3/2 sum_p S_{1,2}(1/p)
# S_{1,2}(x) = sum_{m>=2} H_{m-1} x^m / m^2  (Nielsen); tail p>P estimated by ~sum 1/(4p^2) (kappa3: 3/2*x^2/4 per prime)
import mpmath as mp
from sympy import primerange
mp.mp.dps=30
def S12(x):
    s=mp.mpf(0); H=mp.mpf(0); m=1; t=x
    while True:
        m+=1; H+=mp.mpf(1)/(m-1); t*=x; term=H*t/m**2; s+=term
        if term<mp.mpf(10)**-28: return s
P=3_000_000; k2=mp.mpf(0); k3=mp.mpf(0)
for p in primerange(2,P):
    x=mp.mpf(1)/p; k2+=mp.log(1-x)+mp.polylog(2,x)
    if p<200000: k3+=S12(x)
    else: k3+=x*x/4+x**3/2*mp.mpf(1)/3*mp.mpf(3)/2  # S12 ~ x^2/4 + x^3/6 ... (negligible)
print('kappa2_arith =',mp.nstr(k2/2,8),' kappa3_arith =',mp.nstr(1.5*k3,8))
