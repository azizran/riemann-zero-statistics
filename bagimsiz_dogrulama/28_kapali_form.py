# -*- coding: utf-8 -*-
"""
KAPALI FORM ARAYIŞI: |ζ|(τ) = taban + kesime-yakınlık terimi
Aday: f(x) = A + B·x^{−p},  x = τ_c − τ     (p = 1/2, 1, 3/2, 2 serbest p)
      f(x) = A + B·e^{−x/x0}
      f(x) = A + B·ln(1/x)
İki kesimde (0.86, 1.00) AYNI (A,B,p) tutuyor mu?  ← asıl test
"""
import json, numpy as np
from scipy.optimize import curve_fit

D = json.load(open("24_T3_sonuc.json"))
def veri(tc):
    sat = D[tc if tc in D else str(tc)]
    tau = np.array([s[0] for s in sat]); mod = np.array([s[1] for s in sat])
    m = tau <= float(tc) - 0.02            # kenar artefaktını dışla (ilk/son yarım pencere)
    return tau[m], mod[m]

formlar = {
 "A + B·x^(−1/2)":      (lambda x,A,B: A + B*x**-0.5, [0.3,0.01]),
 "A + B·x^(−1)":        (lambda x,A,B: A + B/x,       [0.3,0.005]),
 "A + B·x^(−2)":        (lambda x,A,B: A + B/x**2,    [0.3,0.0002]),
 "A + B·x^(−p) (serb.)":(lambda x,A,B,p: A + B*x**(-p), [0.3,0.005,1.0]),
 "A + B·e^(−x/x0)":     (lambda x,A,B,x0: A + B*np.exp(-x/x0), [0.3,0.2,0.1]),
 "A + B·ln(1/x)":       (lambda x,A,B: A + B*np.log(1/x), [0.3,0.03]),
}
print(f"{'form':>22} | {'τc=0.86 RMS':>12} {'τc=1.00 RMS':>12} | {'parametreler (0.86)':>34}")
en_iyi=None
for ad,(f,p0) in formlar.items():
    rms={}; par={}
    for tc in [0.86, 1.00]:
        tau,mod = veri(tc); x = tc - tau
        try:
            p,_ = curve_fit(f, x, mod, p0=p0, maxfev=20000)
            rms[tc]=float(np.sqrt(np.mean((mod-f(x,*p))**2))); par[tc]=p
        except Exception as e:
            rms[tc]=np.nan; par[tc]=None
    print(f"{ad:>22} | {rms[0.86]:12.5f} {rms[1.00]:12.5f} | "
          + (", ".join(f"{v:.4f}" for v in par[0.86]) if par[0.86] is not None else "—"))
    if en_iyi is None or rms[0.86] < en_iyi[1]: en_iyi=(ad, rms[0.86], par[0.86])
print(f"\nen iyi (τc=0.86): {en_iyi[0]}  RMS={en_iyi[1]:.5f}  parametreler={[f'{v:.4f}' for v in en_iyi[2]]}")
# serbest p'li formun iki kesimdeki parametreleri
tau,mod=veri(0.86); x=0.86-tau
f=lambda x,A,B,p: A+B*x**(-p)
p1,_=curve_fit(f,x,mod,p0=[0.3,0.005,1.0],maxfev=20000)
tau2,mod2=veri(1.00); x2=1.00-tau2
p2,_=curve_fit(f,x2,mod2,p0=[0.3,0.005,1.0],maxfev=20000)
print(f"\nserbest-p formu:  τc=0.86 → A={p1[0]:.4f} B={p1[1]:.5f} p={p1[2]:.4f}")
print(f"                  τc=1.00 → A={p2[0]:.4f} B={p2[1]:.5f} p={p2[2]:.4f}")
print("\nAynı (A,B,p) iki kesime de uyuyor mu? → ortak fit:")
tum_x=np.r_[x,x2]; tum_y=np.r_[mod,mod2]
pc,_=curve_fit(f,tum_x,tum_y,p0=[0.3,0.005,1.0],maxfev=20000)
print(f"  ortak: A={pc[0]:.4f} B={pc[1]:.5f} p={pc[2]:.4f}   "
      f"RMS(0.86)={np.sqrt(np.mean((mod-f(x,*pc))**2)):.5f}  RMS(1.00)={np.sqrt(np.mean((mod2-f(x2,*pc))**2)):.5f}")
json.dump(dict(en_iyi=en_iyi[0], rms=en_iyi[1], par=list(map(float,en_iyi[2])),
               ortak=dict(A=float(pc[0]),B=float(pc[1]),p=float(pc[2])),
               ayri=dict(c086=list(map(float,p1)), c100=list(map(float,p2)))),
          open("28_kapali_form.json","w"), indent=1)
