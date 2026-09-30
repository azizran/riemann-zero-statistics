"""Not 8 v3.2 Şekil: merdivenin ilk iki basamağı — mühürlü veri, oran sanısı (commit'li tahminler), eski a_k modeli."""
import json,sys,numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0,'200_configs'); from importlib import import_module; o=import_module('200a_ortak')
d=json.load(open('200_configs/HUKUM_200C.json'))['per_window']; C=('C1','C2','C3','C4')
L=np.array([9.35,10.59,11.66]+[d[w]['Lbar0'] for w in C])
obs={'b0k2':([1.8191,1.8899,1.9348]+[d[w]['point_k2k3_b0b1b2'][0] for w in C],[0.0045,0.0044,0.0039]+[d[w]['SE_k2_calib'][0] for w in C]),
     'b0k3':([-2.0887,-2.1237,-2.1240]+[d[w]['point_k2k3_b0b1b2'][1] for w in C],[0.0181,0.0195,0.0130]+[d[w]['SE_k3_calib'][0] for w in C]),
     'b1k2':([0.4672,0.5200,0.5586]+[d[w]['point_k2k3_b0b1b2'][2] for w in C],[0.0022,0.0014,0.0011]+[d[w]['SE_k2_calib'][1] for w in C]),
     'b1k3':([-0.0465,-0.0519,-0.0509]+[d[w]['point_k2k3_b0b1b2'][3] for w in C],[0.0027,0.0021,0.0010]+[d[w]['SE_k3_calib'][1] for w in C])}
pred={'b0k2':[1.82488,1.88578,1.93303,2.02998,2.10677,2.17138,2.25415],'b0k3':[-2.09511,-2.10863,-2.11837,-2.13636,-2.14879,-2.15812,-2.16867],
      'b1k2':[0.469063,0.519006,0.558432,0.640990,0.707828,0.764964,0.839246],'b1k3':[-0.048774,-0.049155,-0.049007,-0.047964,-0.046687,-0.045444,-0.043741]}
Lf=np.linspace(9,23,60)
ak={'b0k2':[float(o.kap_vec(2,x,0))-0.088124 for x in Lf],'b0k3':[float(o.kap_vec(3,x,0))+0.233653 for x in Lf],
    'b1k2':[float(o.kap_vec(2,x,1))-0.088124 for x in Lf],'b1k3':[float(o.kap_vec(3,x,1))+0.233653 for x in Lf]}
tit={'b0k2':r'$b=0$: $\kappa_2(\log|\zeta|)$','b0k3':r'$b=0$: $\kappa_3(\log|\zeta|)$','b1k2':r"$b=1$: $\kappa_2(\log|\zeta'(\rho)|)$",'b1k3':r"$b=1$: $\kappa_3(\log|\zeta'(\rho)|)$"}
chi={'b0k2':'6.4','b0k3':'2.5','b1k2':'10.9','b1k3':'11.6'}; chia={'b0k2':'','b0k3':'58','b1k2':'1780','b1k3':'16 826'}
pos={'b0k2':(0.03,0.93,'left'),'b0k3':(0.97,0.93,'right'),'b1k2':(0.03,0.93,'left'),'b1k3':(0.97,0.40,'right')}
plt.rcParams.update({'font.size':9,'font.family':'serif'})
fig,ax=plt.subplots(2,2,figsize=(7.2,5.4))
for a,k in zip(ax.ravel(),['b0k2','b0k3','b1k2','b1k3']):
    y,e=obs[k]; a.errorbar(L,y,yerr=e,fmt='o',ms=4,color='k',capsize=2,label='data (sealed)',zorder=3)
    a.plot(L,pred[k],'s-',ms=4,color='#c97a12',lw=1.3,label='ratios conjecture',zorder=2)
    a.plot(Lf,ak[k],'--',color='#3a5a9a',lw=1.1,label=r'tilted CUE ($N=L$) + $a_k$',zorder=1)
    a.set_title(tit[k]); a.set_xlabel(r'$L=\log(t/2\pi)$')
    x0,y0,ha=pos[k]; s=r'$\chi^2_{7}$: ratios '+chi[k]+(r', $a_k$ '+chia[k] if chia[k] else '')
    a.text(x0,y0,s,transform=a.transAxes,va='top',ha=ha,fontsize=8); a.grid(alpha=.25)
ax[0,0].legend(loc='lower right',fontsize=7.5,frameon=False)
fig.tight_layout(); fig.savefig('figures/ladder_ratios.pdf'); fig.savefig('figures/ladder_ratios.png',dpi=150)
