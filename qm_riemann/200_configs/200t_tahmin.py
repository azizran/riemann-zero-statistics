# Exact tilted-CUE predictions for the small-gap hump law (theorem check)
import mpmath as mp
def logM(n,k):  # log Keating-Snaith moment M_n(k)=E|Lambda_n|^{2k}
    return mp.fsum(mp.loggamma(j)+mp.loggamma(j+2*k)-2*mp.loggamma(j+k) for j in range(1,n+1))
def pred(N):
    n=N-2
    out={}
    for k in (0.5,1.0):   # E[(M/s^2)^{2k}] (theta units)
        out[f'm{2*k:g}']=float(mp.e**(logM(n,k+2)-logM(n,2))/4**(2*k))
    out['Elog']=float(mp.fsum(mp.digamma(j+4)-mp.digamma(j+2) for j in range(1,n+1))-mp.log(4))
    out['Vlog']=float(mp.fsum(mp.psi(1,j+4)-mp.psi(1,j+2)/2 for j in range(1,n+1)))
    # unfolded: M/s~^2 = (2pi/N)^2 M/s^2
    c=(2*mp.pi/N)**2
    out['u_m1']=float(out['m1']*c); out['u_m2']=float(out['m2']*c**2)
    out['u_Elog']=float(out['Elog']+mp.log(c))
    out['rate']=float(N**2*(N**2-1)/(72*mp.pi))  # E#{gaps<eps}/eps^3, theta units
    return out
if __name__=='__main__':
    for N in (4,6,10,16,24):
        print(N,{k:round(v,5) for k,v in pred(N).items()})
