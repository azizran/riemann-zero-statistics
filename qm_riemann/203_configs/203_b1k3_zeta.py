"""203-b1κ₃ adım 5: ζ için log|ζ′(ρ)| üçüncü kümülantı (oran sanısı + sıfır yoğunluğu ağırlığı), vektörlü.
Yapı 203_b1k3_cue.py ile aynı (CUE'de kesin); aritmetik parçalar 203_yerel4*.py'de asal başına θ-integrali/sayısal
türevle doğrulandı. Asal toplamları p < 10⁴ + kuyruk integralleri; p = 2 dahil hepsi kesin formda."""
import numpy as np, mpmath as mp, sys
np.seterr(all='ignore')
sys.path.insert(0,'.'); sys.path.insert(0,'../200_configs'); from importlib import import_module
hz=import_module('203_hizli')
P=hz.Ps[:,None]; l=np.log(P); PM=float(hz.Ps[-1]); lPM=np.log(PM)
Psm=hz.Ps[hz.Ps<1000][:,None]; lsm=np.log(Psm)
XC=0.45; DSVAR=1
XMAX=3.0   # e^{−L·x}, x > XMAX ⇒ ≤ e^{−3L}: takas terimleri sıfırlanır
def tail_log(k,s):   # Σ_{p>PM} log^k p · p^{−s} ≈ ∫_PM^∞ log^{k−1}t t^{−s} dt  (k = 1..4), s > 1
    m=s-1.0; L1=lPM; e=PM**(-m)
    if k==1: return e/m
    if k==2: return e*(L1/m+1/m**2)
    if k==3: return e*(L1**2/m+2*L1/m**2+2/m**3)
    return e*(L1**3/m+3*L1**2/m**2+6*L1/m**3+6/m**4)
def psum(expr): return np.sum(expr,axis=0)
Wc=0.75
# ζ(1+s)·s tablosu (çift takas için), s ∈ [−6, 6]
_sg=np.concatenate([np.linspace(-6,-1e-4,700),np.linspace(1e-4,6,700)])
from scipy.interpolate import CubicSpline
_zs=CubicSpline(_sg,np.array([float(mp.zeta(1+mp.mpf(s)))*s for s in _sg]))
def zeta1(s): return _zs(s)/s
class Zeta:
    def __init__(s,L):
        s.L=L; s.S=hz.Lset(L)
    # --- üç nokta ---
    def D(s,x1,x2):
        y1=P**(-1-x1); y2=P**(-1-x2)
        return -psum(l**3*y1*y2/((1-y1)*(1-y2)))-tail_log(3,2+x1+x2)
    def T(s,a1,a2,b):
        F=s.S.F; x1,x2=a1+b,a2+b; u=a2-a1
        t=s.D(x1,x2)
        v1=np.nan_to_num(F(x1)*hz.K(u,x1),nan=0.0,posinf=0.0,neginf=0.0)
        v2=np.nan_to_num(F(x2)*hz.K(-u,x2),nan=0.0,posinf=0.0,neginf=0.0)
        t=t+np.where((u>-Wc)&(x1<XMAX),v1,0.0)+np.where((-u>-Wc)&(x2<XMAX),v2,0.0)
        return t
    # --- dört nokta: 3+1 ---
    def fjk(s,x,uj,uk):
        qj=P**(-1-uj); qk=P**(-1-uk); px=P**(-x)
        om=(1-px)*(1-1/P)*px*(1-P**(-1+x))/(1-2/P+P**(-1-x))**2
        sj=2+uj+uk
        return psum(l**2*om*qj*qk/((1-qj)*(1-qk)))+tail_log(2,sj+x)-tail_log(2,sj+2*x)   # TEFTİŞ DÜZELTMESİ (29 Eyl): log²p ⇒ k=2 (önce yanlışlıkla 3)
    def Q(s,A,b):
        F=s.S.F; y=[P**(-1-(a+b)) for a in A]
        tot=psum(l**4*y[0]/(1-y[0])*y[1]/(1-y[1])*y[2]/(1-y[2]))
        for i in range(3):
            j,k=[m for m in range(3) if m!=i]; x=A[i]+b; uj=A[j]-A[i]; uk=A[k]-A[i]
            ok=(uj>-Wc)&(uk>-Wc)
            ok=ok&(x<XMAX)
            val=F(x)*(s.fjk(x,uj,uk)+hz.K(uj,x)*hz.K(uk,x))
            tot=tot+np.where(ok,np.nan_to_num(val,nan=0.0,posinf=0.0,neginf=0.0),0.0)
        return tot
    # --- dört nokta: 2+2 ---
    def Q2(s,a1,a2,b1,b2):
        F=s.S.F; H=s.S.H; A=[a1,a2]; B=[b1,b2]; X=lambda i,j: A[i]+B[j]
        # kimlik: H11H22 + H12H21 + Σ c4
        Pp,lp_=Psm,lsm
        u1,u2=Pp**(-.5-a1),Pp**(-.5-a2); v1,v2=Pp**(-.5-b1),Pp**(-.5-b2)
        sf=lambda u,v: u*v/(1-u*v)
        E4=u1*u2*v1*v2/((u1-u2)*(v1-v2))*(sf(u1,v1)-sf(u1,v2)-sf(u2,v1)+sf(u2,v2))
        e=lambda y: lp_**2*y/(1-y)
        c4=lp_**4*E4-e(Pp**(-1-X(0,0)))*e(Pp**(-1-X(1,1)))-e(Pp**(-1-X(0,1)))*e(Pp**(-1-X(1,0)))
        PMs=float(Psm[-1,0]); m4=1+a1+a2+b1+b2
        tail4=PMs**(-m4)*(np.log(PMs)**3/m4+3*np.log(PMs)**2/m4**2+6*np.log(PMs)/m4**3+6/m4**4)
        tot=H(X(0,0))*H(X(1,1))+H(X(0,1))*H(X(1,0))+psum(c4)-tail4
        # tek takaslar
        for i in range(2):
            for j in range(2):
                ip,jp=1-i,1-j; x=X(i,j); ua=A[ip]-A[i]; ub=B[jp]-B[j]
                ok=(ua>-Wc)&(ub>-Wc)&(x<XMAX)
                val=F(x)*(s.fij(A[i],A[ip],B[j],B[jp])+hz.K(ua,x)*hz.K(ub,x))
                tot=tot+np.where(ok,np.nan_to_num(val,nan=0.0,posinf=0.0,neginf=0.0),0.0)
        # çift takas
        tot=tot+s.DS(a1,a2,b1,b2)
        return tot
    def fij(s,ai,aip,bj,bjp):
        """f_{i'j'} = Σ_p ∂α_{i'}∂β_{j'} log E_p (kapalı biçim; 203_yerel4_turev.Q2_parts_cf ile aynı)."""
        H=s.S.H; lg=l
        b1=P**(-.5+ai); b2=P**(-.5-bjp); dj=P**(-.5-bj); e1=P**(-.5+bj); e2=P**(-.5-aip); ci=P**(-.5-ai)
        hv=lambda v: (1-v*ci)/(1-v*e1)
        C1=1-dj/b1; E=1+C1*(hv(b1)-1)
        dC1=(1-dj/b1)*(-lg*b2/b1)/(1-b2/b1); dC2=-lg*(1-dj/b2)/(1-b1/b2)
        q1=b1*e2; q2=b2*e2; l1=-lg*q1/(1-q1); l2=-lg*q2/(1-q2)
        dEb=dC1*(hv(b1)-1)+dC2*(hv(b2)-1); dEa=C1*hv(b1)*l1; dEab=dC1*hv(b1)*l1+dC2*hv(b2)*l2
        term=dEab/E-dEa*dEb/E**2
        x2=aip+bjp; ey=lambda x: lg**2*P**(-1-x)/(1-P**(-1-x))
        # öncü davranış terim ≈ e(x2) = log²p p^{−1−x2}/(1−p^{−1−x2}) (fark ~p^{−2}; sayısal doğrulandı p ≤ 1e6)
        # (İLK SÜRÜM HATALIYDI: −H(x1+x2) tam eklenip kesik toplamla çıkarılıyordu ⇒ ~−22 sahte kayma)
        return H(x2)+psum(term-ey(x2))
    def DS(s,a1,a2,b1,b2):
        tot_sh=a1+a2+b1+b2; big=tot_sh>XMAX
        A=[a1,a2]; B=[b1,b2]; num=1.0
        for k in range(2):
            for m in range(2):
                x=np.minimum(A[m]+B[k],5.9); num=num*zeta1(x)*zeta1(-x)
        den=zeta1(b2-b1)*zeta1(b1-b2)*zeta1(a2-a1)*zeta1(a1-a2)
        Pp=Psm; f=lambda sh: Pp**(-.5-sh)
        Un=[f(-b1),f(-b2)]; Ud=[f(a1),f(a2)]; Cn=[f(-a1),f(-a2)]; Cd=[f(b1),f(b2)]
        h=lambda v: (1-v*Ud[0])*(1-v*Ud[1])/((1-v*Un[0])*(1-v*Un[1]))
        E=1.0
        for li,bv in enumerate(Cn):
            C=(1-Cd[0]/bv)*(1-Cd[1]/bv)/(1-Cn[1-li]/bv); E=E+C*(h(bv)-1)
        z=lambda sh: 1/(1-Pp**(-1.0-sh))
        Yp=1.0
        for m in range(2):
            for k in range(2): Yp=Yp*z(-A[m]-B[k])*z(A[m]+B[k])
        Yp=Yp/(z(b2-b1)*z(b1-b2)*z(a2-a1)*z(a1-a2))
        Yp=Yp*(1-1/Pp)**4          # sıfırlanan ζ(1−β_k+δ_k), ζ(1−α_l+γ_l) paydalarının yerel çarpanları (tam Y_p)
        r=E/Yp; Aval=np.prod(r,axis=0)          # işaret korunur
        # Euler çarpımı yalnız bütün x_lk < XC'de yakınsar (p^{−2+2x} terimleri). Dışında: tekil limitleri koruyan
        # ara değer — a₁→0 (x₁₂→0): A_DS → A(x₂₁); a₂→0 (x₂₂→0): A_DS → A(x₁₁)  [β₂ = 0, sayısal doğrulandı].
        Ax=lambda x: s.S.Aexp(x)
        x11,x12,x21,x22=a1+b1,a1+b2,a2+b1,a2+b2
        if DSVAR==1: w=x22/(x12+x22)
        else: w=x22**2/(x12**2+x22**2)
        interp=w*Ax(x12)*Ax(x21)+(1-w)*Ax(x11)*Ax(x22)
        xm=np.maximum(np.maximum(x11,x12),np.maximum(x21,x22))
        Aval=np.where(xm<XC,Aval,interp)
        val=np.exp(-s.L*tot_sh)*num/den*Aval
        return np.where(big,0.0,np.nan_to_num(val,nan=0.0,posinf=0.0,neginf=0.0))
    # --- b1 κ₃ integrandı (203_b1k3_cue.integrand ile aynı birleşim) ---
    def integrand(s,a,b,c):
        G=s.S.G; L=s.L; z0=0.0*a
        m=lambda u: G(u)/L; U=lambda u,v: s.T(u,v,z0)/L; V=lambda u,w: G(u+w)+(s.T(u,z0,w)+s.T(w,z0,u))/L
        M3=s.T(a,b,c)+(s.Q([a,b,z0],c)+s.Q2(a,b,c,z0))/L
        M3u=s.Q([a,b,c],z0)/L
        kxxb=M3-U(a,b)*m(c)-V(a,c)*m(b)-V(b,c)*m(a)+2*m(a)*m(b)*m(c)
        kxxx=M3u-U(a,b)*m(c)-U(a,c)*m(b)-U(b,c)*m(a)+2*m(a)*m(b)*m(c)
        return -(kxxx+3*kxxb)/4

def k3_b1(L,n=6,amin_f=0.01,Z=None,scale=(1.0,1.0003,0.9997)):
    import importlib; cue=importlib.import_module('203_b1k3_cue')
    Z=Z or Zeta(L)
    bps=lambda s: sorted(set([amin_f/L*s,0.03/L*s,0.1/L*s,0.3/L*s,1/L*s,3/L*s,8/L*s,1.0*s,3.0*s,8.0*s,16.0*s,40.0*s]))
    (A,WA),(B,WB),(C,WC)=[cue.gl(bps(s),n) for s in scale]
    aa,bb=np.meshgrid(A,B,indexing='ij'); aa=aa.ravel(); bb=bb.ravel(); wab=np.outer(WA,WB).ravel(); tot=0.0
    for c,wc in zip(C,WC):
        v=Z.integrand(aa,bb,c+0*aa); tot+=wc*np.nansum(wab*v)
    return tot
def k3_b1_delta(L,n=6,s0=0.01,Z=None,scale=(1.0,1.0003,0.9997)):
    """κ₃ = κ₃^CUE(L) (kesin Palm, analitik) + ∭ Δ,  Δ = I_ζ − I_CUE(N=L).
    Δ yalnız a,b,c ≥ s0 bölgesinde (hızlı ζ güvenilir; ışın sınaması) integre edilir; 0 ≤ (bir değişken) < s0 levhaları
    Δ(s0,·,·)·s0/2 ile (doğrusal uzatma) eklenir; kenar/köşe (s0², s0³) ihmal. Levha katkısı ayrıca raporlanır."""
    import importlib; cue=importlib.import_module('203_b1k3_cue'); o=importlib.import_module('200a_ortak')
    Z=Z or Zeta(L)
    D=lambda a,b,c: Z.integrand(a,b,c)-cue.integrand(a,b,c,L)
    bps=lambda s: sorted(set([s0*s]+[x*s for x in (2*s0,0.1/L,0.3/L,1/L,3/L,8/L,1.0,3.0,8.0,16.0,40.0) if x>s0*1.5]))   # ızgara s0'da başlar (üst üste binme yok)
    (A,WA),(B,WB),(C,WC)=[cue.gl(bps(s),n) for s in scale]
    aa,bb=np.meshgrid(A,B,indexing='ij'); aa=aa.ravel(); bb=bb.ravel(); wab=np.outer(WA,WB).ravel(); tot=0.0
    for c,wc in zip(C,WC):
        tot+=wc*np.nansum(wab*D(aa,bb,c+0*aa))
    # levhalar: a<s0, b<s0, c<s0 (her biri s0 kalınlığında, Δ sınırdaki değerle)
    slab=0.0   # levha: Δ a→0'da doğrusal sıfıra gider (ışın sınaması) ⇒ ∫_0^{s0} Δ ≈ s0·Δ(s0)/2
    bb2,cc2=np.meshgrid(B,C,indexing='ij'); wbc=np.outer(WB,WC).ravel(); slab+=s0*np.nansum(wbc*D(np.full(bb2.size,s0),bb2.ravel(),cc2.ravel()))
    aa2,cc2=np.meshgrid(A,C,indexing='ij'); wac=np.outer(WA,WC).ravel(); slab+=s0*np.nansum(wac*D(aa2.ravel(),np.full(aa2.size,s0),cc2.ravel()))
    slab+=s0*np.nansum(wab*D(aa,bb,np.full(aa.size,s0)))
    slab=slab/2
    kC=float(o.kap_vec(3,L,1))
    return kC+tot+slab, kC, tot, slab

if __name__=='__main__':
    L=float(sys.argv[1]); n=int(sys.argv[2]); DSVAR=int(sys.argv[3]) if len(sys.argv)>3 else 1
    af=float(sys.argv[4]) if len(sys.argv)>4 else 0.01
    if af<0: r=k3_b1_delta(L,n,s0=-af); print(f"L={L} n={n} s0={-af} κ3_b1={r[0]:.6f} (CUE {r[1]:.6f} Δ-iç {r[2]:+.6f} levha {r[3]:+.6f})",flush=True)
    else: print(f"L={L} n={n} DSVAR={DSVAR} amin={af}/L κ3_b1={k3_b1(L,n,amin_f=af):.6f}",flush=True)
