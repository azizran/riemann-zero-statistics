"""203-b1κ₃ adım 1: CUE dört-nokta ortalamaları (oran reçetesinden) ve küçük N'de kesin sınama.
X(α) = (Λ'/Λ)(α) = Σ_n Tr(U^n) e^{−nα};  x_ij = α_i + β_j.
 Q(α₁,α₂,α₃;β)   = E[X X X X̄]  = Σ_i F(x_i) K(α_j−α_i,x_i) K(α_k−α_i,x_i)
 Q2(α₁,α₂;β₁,β₂) = E[X X X̄ X̄]  = H₁₁H₂₂ + H₁₂H₂₁ + Σ_{i,j} F(x_ij)[H(x_{i'j'}) + K(α_{i'}−α_i,x_ij)K(β_{j'}−β_j,x_ij)] + (çift takas)
 çift takas = e^{−N(α₁+α₂+β₁+β₂)} Π_{k,l} z(−α_l−β_k)z(α_l+β_k) / [z(β₂−β₁)z(β₁−β₂)z(α₂−α₁)z(α₁−α₂)]."""
import numpy as np
z=lambda x: 1/(1-np.exp(-x)); zp=lambda u: -1/np.expm1(u)
H=lambda x: np.exp(x)/np.expm1(x)**2
def mk(N):
    F=lambda x: -np.exp(-N*x)*np.exp(x)/np.expm1(x)**2
    K=lambda u,x: zp(u)-zp(u+x)
    def Q(a,b):  # a: 3 kaydırma, b: 1
        tot=0.0
        for i in range(3):
            j,k=[m for m in range(3) if m!=i]; x=a[i]+b
            tot+=F(x)*K(a[j]-a[i],x)*K(a[k]-a[i],x)
        return tot
    def Q2(a,b):
        x=lambda i,j: a[i]+b[j]
        tot=H(x(0,0))*H(x(1,1))+H(x(0,1))*H(x(1,0))
        for i in range(2):
            for j in range(2):
                ip,jp=1-i,1-j
                tot+=F(x(i,j))*(H(x(ip,jp))+K(a[ip]-a[i],x(i,j))*K(b[jp]-b[j],x(i,j)))
        num=np.prod([z(-a[l]-b[k])*z(a[l]+b[k]) for k in range(2) for l in range(2)])
        den=z(b[1]-b[0])*z(b[0]-b[1])*z(a[1]-a[0])*z(a[0]-a[1])
        return tot+np.exp(-N*(sum(a)+sum(b)))*num/den
    return Q,Q2
def brute(N,a,b,conj,M=48):
    """Weyl yoğunluğuyla kesin ortalama (periyodik ızgara); conj: kaç X̄ var (b'nin uzunluğu)."""
    th=2*np.pi*np.arange(M)/M; grids=np.meshgrid(*([th]*N),indexing='ij'); T=[g.ravel() for g in grids]
    V=np.ones_like(T[0])
    for p in range(N):
        for q in range(p+1,N): V*=np.abs(np.exp(1j*T[p])-np.exp(1j*T[q]))**2
    V/=V.sum()
    X=lambda al: sum(np.exp(1j*t-al)/(1-np.exp(1j*t-al)) for t in T)
    Xb=lambda be: sum(np.exp(-1j*t-be)/(1-np.exp(-1j*t-be)) for t in T)
    f=np.ones_like(T[0],dtype=complex)
    for al in a: f*=X(al)
    for be in b: f*=Xb(be)
    return (V*f).sum()
if __name__=='__main__':
    rng=np.random.default_rng(2)
    for N,M in ((1,256),(2,96),(3,40)):
        Q,Q2=mk(N)
        for _ in range(2):
            a=rng.uniform(0.2,0.9,3); b=rng.uniform(0.2,0.9,1)
            print(f'N={N} Q : formül {Q(a,b[0]):+.8f}  kesin {brute(N,a,b,1,M).real:+.8f}')
            a2=rng.uniform(0.2,0.9,2); b2=rng.uniform(0.2,0.9,2)
            print(f'N={N} Q2: formül {Q2(a2,b2):+.8f}  kesin {brute(N,a2,b2,2,M).real:+.8f}')
