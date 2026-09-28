# KALEM 200-C predictions: height dependence of the prime-imprint ladder (kappa3 and kappa2 shifts), from sealed 200-A/B per-window values
import numpy as np, json
from scipy.optimize import least_squares
A={'k3':0.233653,'k2':-0.088124}
# (L, shift, se) per window, sealed 200-A (b0,b1) and 200-B (b2, eps=0.2)
D={'k3':{0:[(9.313,0.2983,0.0181),(10.582,0.2729,0.0195),(11.650,0.2791,0.0130)],
         1:[(9.343,0.1439,0.0027),(10.589,0.1457,0.0021),(11.658,0.1518,0.0010)],
         2:[(9.355,0.0764,0.0059),(10.594,0.0844,0.0052),(11.662,0.0928,0.0046)]},
   'k2':{0:[(9.313,-0.0857,0.0045),(10.582,-0.0787,0.0044),(11.650,-0.0818,0.0039)],
         1:[(9.343,-0.0733,0.0022),(10.589,-0.0716,0.0014),(11.658,-0.0731,0.0011)],
         2:[(9.355,-0.1001,0.0078),(10.594,-0.1047,0.0061),(11.662,-0.0938,0.0052)]}}
Lnew=[14.194,16.577,18.884,22.306]  # kesin L̄ (yalnız konumlardan, 27 Eyl pilot)
out={}
for r in ('k3','k2'):
    a=A[r]; res={}
    # H_C: constant per rung (weighted mean)
    HC={b:float(np.average([s for _,s,_ in v],weights=[1/e**2 for *_,e in v])) for b,v in D[r].items()}
    # H_S1: a + c_b/L (gamma=1), weighted LS for c_b
    HS1={}
    for b,v in D[r].items():
        L=np.array([x[0] for x in v]); s=np.array([x[1] for x in v]); e=np.array([x[2] for x in v])
        X=1/L; c=np.sum((s-a)*X/e**2)/np.sum(X**2/e**2); HS1[b]=float(c)
    # H_Sg: a + c_b L^-g, common g (fit over all rungs)
    def resid(p):
        g=p[0]; out=[]
        for i,(b,v) in enumerate(D[r].items()):
            for L,s,e in v: out.append((s-(a+p[1+i]*L**(-g)))/e)
        return out
    fit=least_squares(resid,[0.5,-1,-1,-1])
    J=fit.jac; cov=np.linalg.inv(J.T@J)*max(1,np.sum(np.array(fit.fun)**2)/(9-4))
    g,sg=fit.x[0],np.sqrt(cov[0,0]); chi2=float(np.sum(np.array(fit.fun)**2))
    res['H_C']=HC; res['H_S1_c']=HS1; res['H_Sg']={'gamma':float(g),'gamma_se':float(sg),'c':fit.x[1:].tolist(),'chi2_9_4':chi2}
    HC_se={b:float(1/np.sqrt(sum(1/e**2 for *_,e in v))) for b,v in D[r].items()}
    HS1_se={}
    for b,v in D[r].items():
        L=np.array([x[0] for x in v]); e=np.array([x[2] for x in v]); HS1_se[b]=float(1/np.sqrt(np.sum((1/L)**2/e**2)))
    tab={}
    for b in (0,1,2):
        tab[b]={}
        for Ln in Lnew:
            Jg=np.zeros(4); Jg[0]=-fit.x[1+b]*np.log(Ln)*Ln**(-g); Jg[1+b]=Ln**(-g)
            tab[b][f'{Ln}']={'H_C':[HC[b],HC_se[b]],'H_S1':[a+HS1[b]/Ln,HS1_se[b]/Ln],'H_Sg':[a+fit.x[1+b]*Ln**(-g),float(np.sqrt(Jg@cov@Jg))]}
    res['tablo']=tab; out[r]=res
    print(f'== {r}: H_Sg ortak gamma = {g:.3f} ± {sg:.3f} (chi2 {chi2:.2f}/5)  c_b = {np.round(fit.x[1:],3)}')
    for b in (0,1,2):
        print(f'  b={b}: H_C {HC[b]:+.4f}±{HC_se[b]:.4f} | ' + ' | '.join(f'L={Ln}: S1 {tab[b][str(Ln)]["H_S1"][0]:+.4f}±{tab[b][str(Ln)]["H_S1"][1]:.4f} Sg {tab[b][str(Ln)]["H_Sg"][0]:+.4f}±{tab[b][str(Ln)]["H_Sg"][1]:.4f}' for Ln in Lnew))
json.dump(out,open('200_configs/200c_tahmin.json','w'),indent=1)
