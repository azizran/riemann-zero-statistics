"""203-b1κ₃ adım 3: dört-nokta aritmetik parçaların asal başına kapalı biçimleri + θ-integraliyle sınama.
Genel yerel ortalama: E_θ[Π_{e∈Un}(1−w e)^{-1} Π_{c∈Ud}(1−w c) Π_{b∈Cn}(1−w̄ b)^{-1} Π_{d∈Cd}(1−w̄ d)]
  = 1 + Σ_l C_l (h(b_l) − 1),  h(v) = Π(1−v c)/Π(1−v e),  C_l = Π_r(1−d_r/b_l)/Π_{l'≠l}(1−b_{l'}/b_l)   (|Cd| = |Cn|)."""
import numpy as np
from scipy.integrate import quad
def Egen(Un,Ud,Cn,Cd):
    h=lambda v: np.prod([1-v*c for c in Ud],axis=0)/np.prod([1-v*e for e in Un],axis=0)
    tot=1.0
    for l,b in enumerate(Cn):
        C=np.prod([1-d/b for d in Cd],axis=0)
        for l2,b2 in enumerate(Cn):
            if l2!=l: C=C/(1-b2/b)
        tot=tot+C*(h(b)-1)
    return tot
def Enum(Un,Ud,Cn,Cd):
    def f(th,part):
        w=np.exp(2j*np.pi*th); wb=np.conj(w)
        v=np.prod([1-w*c for c in Ud])*np.prod([1-wb*d for d in Cd])/(np.prod([1-w*e for e in Un])*np.prod([1-wb*b for b in Cn]))
        return v.real if part==0 else v.imag
    return quad(lambda t:f(t,0),0,1,limit=400,epsabs=1e-14)[0]
def E4_cf(p,a1,a2,b1,b2):   # Σ_{m1+m2=m3+m4} u1^m1 u2^m2 v1^m3 v2^m4 (log⁴p hariç)
    u1,u2,v1,v2=p**(-.5-a1),p**(-.5-a2),p**(-.5-b1),p**(-.5-b2); s=lambda u,v: u*v/(1-u*v)
    return u1*u2*v1*v2/((u1-u2)*(v1-v2))*(s(u1,v1)-s(u1,v2)-s(u2,v1)+s(u2,v2))
def E4_num(p,a1,a2,b1,b2,M=200):
    u1,u2,v1,v2=p**(-.5-a1),p**(-.5-a2),p**(-.5-b1),p**(-.5-b2); tot=0.0
    for m1 in range(1,M):
        for m2 in range(1,M):
            n=m1+m2
            for m3 in range(1,n):
                tot+=u1**m1*u2**m2*v1**m3*v2**(n-m3)
            if u1**m1*u2**m2<1e-30: break
    return tot
if __name__=='__main__':
    rng=np.random.default_rng(7); a=lambda p,s: p**(-.5-s)
    print('— genel E_θ (2+2, çift takas yapısı) —')
    for p in (2,5,13):
        s=rng.uniform(-0.15,0.15,8)
        Un=[a(p,-s[0]),a(p,-s[1])]; Ud=[a(p,s[2]),a(p,s[3])]; Cn=[a(p,-s[4]),a(p,-s[5])]; Cd=[a(p,s[6]),a(p,s[7])]
        print(p,'kapalı',round(float(Egen(Un,Ud,Cn,Cd)),12),'θ-integral',round(Enum(Un,Ud,Cn,Cd),12))
    print('— E4 (dört-nokta kimlik terimi, asal başına) —')
    for p in (3,7):
        s=rng.uniform(0.05,0.4,4); print(p,'kapalı',E4_cf(p,*s),'doğrudan toplam',E4_num(p,*s,M=60))
