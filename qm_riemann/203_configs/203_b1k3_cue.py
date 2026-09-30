"""203-b1κ₃ adım 2 (CUE): sıfırlar (özdeğerler) üzerinde log|Λ′| üçüncü kümülantı, oran reçetesi + yoğunluk ağırlığı.
κ₃ = ¼[κ(Y,Y,Y) + 3κ(Y,Y,Ȳ)],  κ(Y,Y,Ȳ) = −∭ κ_z(X(a),X(b),X̄(c)),  κ(Y,Y,Y) = −∭ κ_z(X(a),X(b),X(c)).
E_z[M] = E[M] + (E[M X₀] + E[M X̄₀])/N. Kesin: κ₃(N,b=1) = Σ_{j<N}[ψ″(j+2) − ¼ψ″(j+1)]."""
import numpy as np, sys
np.seterr(all='ignore')
z=lambda x: 1/(1-np.exp(-x)); zp=lambda u: -1/np.expm1(u); H=lambda x: np.exp(x)/np.expm1(x)**2
def funcs(N):
    F=lambda x: -np.exp(-N*x)*np.exp(x)/np.expm1(x)**2
    K=lambda u,x: zp(u)-zp(u+x)
    G=lambda x: H(x)+F(x)
    def T(a1,a2,b):
        x1,x2=a1+b,a2+b; return F(x1)*K(a2-a1,x1)+F(x2)*K(a1-a2,x2)
    def Q(a1,a2,a3,b):
        A=[a1,a2,a3]; tot=0.0
        for i in range(3):
            j,k=[m for m in range(3) if m!=i]; x=A[i]+b
            tot=tot+F(x)*K(A[j]-A[i],x)*K(A[k]-A[i],x)
        return tot
    def Q2(a1,a2,b1,b2):
        A=[a1,a2]; B=[b1,b2]; x=lambda i,j: A[i]+B[j]
        tot=H(x(0,0))*H(x(1,1))+H(x(0,1))*H(x(1,0))
        for i in range(2):
            for j in range(2):
                ip,jp=1-i,1-j; tot=tot+F(x(i,j))*(H(x(ip,jp))+K(A[ip]-A[i],x(i,j))*K(B[jp]-B[j],x(i,j)))
        num=1.0
        for k in range(2):
            for l in range(2): num=num*z(-A[l]-B[k])*z(A[l]+B[k])
        den=z(b2-b1)*z(b1-b2)*z(a2-a1)*z(a1-a2)
        return tot+np.exp(-N*(a1+a2+b1+b2))*num/den
    return G,T,Q,Q2
def integrand(a,b,c,N):
    G,T,Q,Q2=funcs(N); z0=0.0*a
    m=lambda u: G(u)/N; U=lambda u,v: T(u,v,z0)/N; V=lambda u,w: G(u+w)+(T(u,z0,w)+T(w,z0,u))/N
    M3=T(a,b,c)+(Q(a,b,z0,c)+Q2(a,b,c,z0))/N
    M3u=Q(a,b,c,z0)/N
    kxxb=M3-U(a,b)*m(c)-V(a,c)*m(b)-V(b,c)*m(a)+2*m(a)*m(b)*m(c)
    kxxx=M3u-U(a,b)*m(c)-U(a,c)*m(b)-U(b,c)*m(a)+2*m(a)*m(b)*m(c)
    return -(kxxx+3*kxxb)/4        # κ₃ integrandı (−∭ ile birlikte)
def gl(bps,n):
    X,W=np.polynomial.legendre.leggauss(n); xs=[];ws=[]
    for a,b in zip(bps,bps[1:]):
        if b>a: xs+=list((b-a)/2*X+(a+b)/2); ws+=list((b-a)/2*W)
    return np.array(xs),np.array(ws)
def k3_b1_cue(N,n=8,amin_f=0.01,scale=(1.0,1.0003,0.9997)):
    bps=lambda s: sorted(set([amin_f/N*s,0.03/N*s,0.1/N*s,0.3/N*s,1/N*s,3/N*s,8/N*s,1.0*s,3.0*s,8.0*s,16.0*s,40.0*s]))
    (A,WA),(B,WB),(C,WC)=[gl(bps(s),n) for s in scale]   # hafif farklı ızgaralar: tam köşegen çakışmasından kaçınır
    tot=0.0
    for c,wc in zip(C,WC):
        aa,bb=np.meshgrid(A,B,indexing='ij'); w=np.outer(WA,WB)*wc
        v=integrand(aa,bb,c+0*aa,N); tot+=np.nansum(w*v)
    return tot
if __name__=='__main__':
    from scipy.special import polygamma
    for N in (2,5,10.0):
        ex=sum(polygamma(2,j+2)-0.25*polygamma(2,j+1) for j in range(1,int(N)))
        for n in (8,12): print(N,n,'formül',round(k3_b1_cue(N,n),6),'kesin',round(float(ex),6))
