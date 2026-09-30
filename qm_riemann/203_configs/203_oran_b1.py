"""203-b1: log|ζ′(ρ)|'nun sıfırlar üzerindeki kümülantları, oran sanısından, sonlu yükseklikte.
Sıfır yoğunluğu ağırlığı: Σ_γ f(γ) = ∫ f(t)(1/2π)[L + 2 Re (ζ′/ζ)(½+it)] dt  ⇒
  E_sıfır[M] = E_t[M] + (1/L)(E_t[M X₀] + E_t[M X̄₀]),  X₀ = X(0).
log ζ′(ρ) = −∫_0^∞ [X(a) − 1_{a<1}/a] da ⇒ κ₂ = ½[Cov(Y,Ȳ) + Cov(Y,Y)]:
  E_z[X(a)X̄(b)] = G(a+b) + (T3(a,0;b) + T3(b,0;a))/L,  E_z[X(a)X(b)] = T3(a,b;0)/L,  E_z[X(a)] = G(a)/L.
T3(α₁,α₂;β) = D + F(x₁)K(α₂−α₁,x₁) + F(x₂)K(α₁−α₂,x₂) (203_oran_k3 ile aynı üç-nokta formülü).
CUE modu (--cue): aynı yol, z-fonksiyonlarıyla; kesin κ₂(N, b=1) = Σ_{j<N}[ψ′(j+2) − ½ψ′(j+1)] ile sınanır."""
import numpy as np, sys
from importlib import import_module
sys.path.insert(0,'.')
k3=import_module('203_oran_k3'); k2=import_module('203_oran_k2')
lp=k3.lp
def make(L,cue=False,Wc=0.75):
    if cue:
        N=L; zp=lambda u: -1/np.expm1(u)
        G=lambda x: np.exp(x)*(1-np.exp(-N*x))/np.expm1(x)**2
        F=lambda x: -np.exp(-N*x)*np.exp(x)/np.expm1(x)**2
        K=lambda u,x: zp(u)-zp(u+x); D=lambda x1,x2: 0.0
    else:
        cache={}
        def G(x):
            k=('G',x)
            if k not in cache: cache[k]=k2.G(x,L)
            return cache[k]
        F=lambda x: k3.F(x,L); K=k3.K
        def D(x1,x2):
            y1=np.exp(-(1+x1)*lp); y2=np.exp(-(1+x2)*lp)
            return -float(np.sum(lp**3*y1*y2/((1-y1)*(1-y2))))
    # Wc: K'nın Euler temsili u > −1'de geçerli; u < −Wc olan takas terimi e^{−L·x}, x ≥ |u| taşır ⇒ atılır (b0 ile aynı)
    def S(a1,a2,b):
        x1,x2=a1+b,a2+b; u=a2-a1; t=0.0
        if u>-Wc: t+=F(x1)*K(u,x1)
        if -u>-Wc: t+=F(x2)*K(-u,x2)
        return t
    def T3(a1,a2,b):
        u=a2-a1; h=1e-3*min(a1+b,a2+b,1.0)   # ölçek: x = α+β (köşegende S, u'da çift ⇒ S(±h) = S(0) + O(h²))
        if abs(u)<h:
            m=(a1+a2)/2; s=S(m-h/2,m+h/2,b)
        else: s=S(a1,a2,b)
        return D(a1+b,a2+b)+s
    def I(a,b):
        gg=G(a)*G(b)/L**2
        return 0.5*(G(a+b)+(T3(a,0,b)+T3(b,0,a))/L-gg+T3(a,b,0)/L-gg)
    return I
def gl(bps,n):
    X,W=np.polynomial.legendre.leggauss(n); xs=[];ws=[]
    for a,b in zip(bps,bps[1:]):
        if b>a: xs+=list((b-a)/2*X+(a+b)/2); ws+=list((b-a)/2*W)
    return np.array(xs),np.array(ws)
def k2_b1(L,cue=False,n=12,amax=40.0,Wc=0.75):
    # a → 0'da terimler ~1/a ile birbirini götürür; asal toplamı kesme hataları orada büyür. Gerçek iç integral
    # a → 0'da doğrusal sıfıra gider ⇒ [0, a_min] şeridi doğrusal ara değerle eklenir (a_min = 0.01/L; katkı ~1e-5).
    I=make(L,cue,Wc); amin=0.01/L
    def inner(a):
        B,WB=gl(sorted(set([a,a+a/2,2*a,a+1/L,a+3/L,a+1.0,a+3.0,a+8.0,a+16.0,a+amax])),n)
        return sum(wb*I(a,b) for b,wb in zip(B,WB))
    A,WA=gl(sorted(set([amin,0.03/L,0.1/L,0.3/L,1/L,3/L,8/L,1.0,3.0,8.0,16.0,amax])),n)
    tot=sum(wa*inner(a) for a,wa in zip(A,WA))+inner(amin)*amin/2
    return 2*tot
if __name__=='__main__':
    mode=sys.argv[1]; n=int(sys.argv[2]); L=float(sys.argv[3]); Wc=float(sys.argv[4]) if len(sys.argv)>4 else 0.75
    print(f'{mode} L={L} n={n} Wc={Wc} κ2_b1={k2_b1(L,mode=="cue",n,Wc=Wc):.6f}',flush=True)
