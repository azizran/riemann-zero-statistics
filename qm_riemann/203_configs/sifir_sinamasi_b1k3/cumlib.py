import numpy as np
def block_sums(x, NB):
    n=len(x)//NB
    return x[:n*NB].reshape(NB,n).sum(axis=1), n
def cum3_jk(A,B,C,NB=32):
    """joint 3rd cumulant of complex arrays A,B,C with block jackknife. returns (theta, se_re, se_im)"""
    M=len(A)//NB*NB; A=A[:M];B=B[:M];C=C[:M]
    arrs=[A,B,C,A*B,A*C,B*C,A*B*C]
    S=[]; 
    for x in arrs:
        s,n=block_sums(x,NB); S.append(s)
    S=np.array(S)          # (7,NB)
    tot=S.sum(axis=1); N=n*NB
    def f(m):
        a,b,c,ab,ac,bc,abc=m
        return abc-ab*c-ac*b-bc*a+2*a*b*c
    th=f(tot/N)
    lo=(tot[:,None]-S)/(N-n)
    thb=f(lo)
    var=lambda v:(NB-1)/NB*np.sum((v-v.mean())**2)
    return th, np.sqrt(var(thb.real)), np.sqrt(var(thb.imag))
def cum2_jk(A,B,NB=32):
    M=len(A)//NB*NB; A=A[:M];B=B[:M]
    arrs=[A,B,A*B]; S=[]
    for x in arrs:
        s,n=block_sums(x,NB); S.append(s)
    S=np.array(S); tot=S.sum(axis=1); N=n*NB
    f=lambda m: m[2]-m[0]*m[1]
    th=f(tot/N); lo=(tot[:,None]-S)/(N-n); thb=f(lo)
    var=lambda v:(NB-1)/NB*np.sum((v-v.mean())**2)
    return th, np.sqrt(var(thb.real)), np.sqrt(var(thb.imag))
