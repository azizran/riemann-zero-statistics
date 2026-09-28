# Theorem 2 check: tilted |Lambda_n(1)| = prod |1-xi_j|, xi_j indep. (BHNY/BNR-type), + arithmetic constant
import numpy as np, mpmath as mp
rng=np.random.default_rng(7)
def sample_xi(j,m):
    out=np.empty(m,complex); k=0
    while k<m:
        B=4*(m-k)+100
        if j==1: z=np.exp(1j*rng.uniform(0,2*np.pi,B))
        else:
            r2=rng.beta(1,j-1,B); z=np.sqrt(r2)*np.exp(1j*rng.uniform(0,2*np.pi,B))  # density ∝ (1-|z|^2)^{j-2}
        acc=rng.uniform(0,16,B)<np.abs(1-z)**4; z=z[acc][:m-k]; out[k:k+len(z)]=z; k+=len(z)
    return out
for N in (10,24):
    n=N-2; m=400000
    L=sum(np.log(np.abs(1-sample_xi(j,m))) for j in range(1,n+1))-np.log(4)
    Elog=float(mp.fsum(mp.digamma(j+4)-mp.digamma(j+2) for j in range(1,n+1))-mp.log(4))
    Vlog=float(mp.fsum(mp.psi(1,j+4)-mp.psi(1,j+2)/2 for j in range(1,n+1)))
    k3=float(mp.fsum(mp.psi(2,j+4)-mp.psi(2,j+2)/4 for j in range(1,n+1)))
    X=np.exp(L); lm=lambda n,k: mp.fsum(mp.loggamma(j)+mp.loggamma(j+2*k)-2*mp.loggamma(j+k) for j in range(1,n+1))
    print(N,'Elog %.4f vs %.4f | Vlog %.4f vs %.4f | k3 %.4f vs %.4f | E[X] %.3f vs %.3f'%(L.mean(),Elog,L.var(),Vlog,((L-L.mean())**3).mean(),k3,X.mean(),float(mp.e**(lm(n,2.5)-lm(n,2))/4)))
# arithmetic variance constant: (1/2) sum_p [log(1-1/p)+Li2(1/p)]
ps=[p for p in range(2,2_000_000) if mp.isprime(p)] if False else None
from sympy import primerange
S=mp.mpf(0)
for p in primerange(2,3_000_000): S+=mp.log(1-mp.mpf(1)/p)+mp.polylog(2,mp.mpf(1)/p)
print('sum_p[log(1-1/p)+Li2(1/p)] =',mp.nstr(S,8),' -> Var corr (1/2)*S =',mp.nstr(S/2,6))
