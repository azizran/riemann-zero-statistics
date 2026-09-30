"""203-b1κ₃ adım 3b: takas terimlerinin türevli aritmetik parçaları — kapalı biçim vs sayısal türev (asal başına).
log A_p = log E_θ,p − log Y_p (takaslı argümanlarda). Köşegen: γ=α, δ=β; türevler γ,δ sabitken alınır."""
import numpy as np, sys
sys.path.insert(0,'.'); from importlib import import_module; Y4=import_module('203_yerel4')
def logYp(p,Un_s,Cn_s,Ud_s,Cd_s):
    """Y_p: Π_{k,l} ζ(1+α'_k+β'_l) Π ζ(1+γ+δ) / (Π ζ(1+α'_k+δ) Π ζ(1+β'_l+γ)), ζ(1+s) → (1−p^{−1−s})^{−1}."""
    z=lambda s: -np.log(1-p**(-1.0-s)); t=0.0
    for a in Un_s:
        for b in Cn_s: t+=z(a+b)
    for g in Ud_s:
        for d in Cd_s: t+=z(g+d)
    for a in Un_s:
        for d in Cd_s: t-=z(a+d)
    for b in Cn_s:
        for g in Ud_s: t-=z(b+g)
    return t
def logA(p,Un_s,Cn_s,Ud_s,Cd_s):   # kaydırmalar: ζ(s+α') ↔ e = p^{−1/2−α'}, ζ(1−s+β') ↔ b = p^{−1/2−β'}
    f=lambda s: p**(-.5-s)
    E=Y4.Egen([f(x) for x in Un_s],[f(x) for x in Ud_s],[f(x) for x in Cn_s],[f(x) for x in Cd_s])
    return np.log(E)-logYp(p,Un_s,Cn_s,Ud_s,Cd_s)
h=1e-5
lp_=lambda p: np.log(p)
# --- Q (3+1), takas i=0: α'=(−β, α1, α2), β'=(−α0); γ=(α0,α1,α2) sabit, δ=β sabit ---
def Q_fjk_num(p,al,be):
    def L(a1,a2): return logA(p,[-be,a1,a2],[-al[0]],[al[0],al[1],al[2]],[be])
    return (L(al[1]+h,al[2]+h)-L(al[1]+h,al[2]-h)-L(al[1]-h,al[2]+h)+L(al[1]-h,al[2]-h))/(4*h*h)
def Q_fjk_cf(p,al,be):
    x=al[0]+be; qj=p**(-1-(al[1]-al[0])); qk=p**(-1-(al[2]-al[0]))
    om=(1-p**-x)*(1-1/p)*p**-x*(1-p**(-1+x))/(1-2/p+p**(-1-x))**2
    return lp_(p)**2*om*qj*qk/((1-qj)*(1-qk))
# --- Q2 (2+2), tek takas (i=0,j=0): α'=(−β0, α1), β'=(−α0, β1); γ=(α0,α1), δ=(β0,β1) sabit ---
def Q2_parts_num(p,al,be):
    def L(a1,b1): return logA(p,[-be[0],a1],[-al[0],b1],[al[0],al[1]],[be[0],be[1]])
    fb=(L(al[1],be[1]+h)-L(al[1],be[1]-h))/(2*h)
    fab=(L(al[1]+h,be[1]+h)-L(al[1]+h,be[1]-h)-L(al[1]-h,be[1]+h)+L(al[1]-h,be[1]-h))/(4*h*h)
    return fb,fab
def Q2_parts_cf(p,al,be):
    lg=lp_(p); x=al[0]+be[0]; u=be[1]-be[0]; x2=al[1]+be[1]
    # f_{j'} aritmetik kısmı (K'nın asal-yerel aritmetik kalıntısı ile aynı olmalı: K_p − yerel ζ'/ζ farkı)
    q=p**(-1-u); phi=(1-1/p)*(1-p**-x)/(1-2/p+p**(-1-x))
    Kp=-lg*q/(1-q)*phi; zeta_part=(-lg*q/(1-q))-(-lg*q*p**-x/(1-q*p**-x))   # yerel (ζ'/ζ)(1+u) − (ζ'/ζ)(1+u+x)
    fb=Kp-zeta_part
    # f_{i'j'} = ∂α₁∂β₁ log E_p (kapalı): E'lerden
    b1=p**(-.5+al[0]); b2=p**(-.5-be[1]); dj=p**(-.5-be[0]); e1=p**(-.5+be[0]); e2=p**(-.5-al[1]); ci=p**(-.5-al[0])
    hv=lambda v: (1-v*ci)/(1-v*e1)
    C1=1-dj/b1; E=1+C1*(hv(b1)-1)
    dC1=(1-dj/b1)*(-lg*b2/b1)/(1-b2/b1); dC2=-lg*(1-dj/b2)/(1-b1/b2)
    q1=b1*e2; q2=b2*e2; l1=-lg*q1/(1-q1); l2=-lg*q2/(1-q2)
    dE_b=dC1*(hv(b1)-1)+dC2*(hv(b2)-1); dE_a=C1*hv(b1)*l1; dE_ab=dC1*hv(b1)*l1+dC2*hv(b2)*l2
    fab=dE_ab/E-dE_a*dE_b/E**2
    return fb,fab
if __name__=='__main__':
    rng=np.random.default_rng(11)
    for p in (2,3,11,101):
        al=rng.uniform(0.05,0.35,3); be=rng.uniform(0.05,0.35)
        print(f'p={p} Q  f_jk: sayısal {Q_fjk_num(p,al,be):+.9f}  kapalı {Q_fjk_cf(p,al,be):+.9f}')
        al2=rng.uniform(0.05,0.35,2); be2=rng.uniform(0.05,0.35,2)
        (nb,nab),(cb,cab)=Q2_parts_num(p,al2,be2),Q2_parts_cf(p,al2,be2)
        # sayısal f_{j'} = ∂β₁ log A_p toplamı; ζ-kısmı ∂β₁ log Y_sw ile birlikte K'yı vermeli: burada yalnız A_p türevini karşılaştır
        print(f'p={p} Q2 ∂β log A_p: sayısal {nb:+.9f}  (K yerel aritmetik) {cb:+.9f} | ∂α∂β log A_p: sayısal {nab:+.9f}')
