# -*- coding: utf-8 -*-
"""
KONVANSİYON-1: makalenin kısıtındaki w ve v TAM OLARAK nedir?
185b:  wg = Σ_b|â_q|/Σ_b ae ;  wog = Σ_b|c^öz|/Σ_b ae ;  mg = wg − wog
186c:  v_b = Σ_b|Ĝ_q| / (ḡ · Σ_b ae)
Kısıt iddiası: sqrt(w_b) = 1.017 − 0.884 v_b
"""
import json, importlib.util
import numpy as np
from pathlib import Path

QM = Path("/Users/ugur/Desktop/Deney/qm_riemann"); SCR = QM / "scratchpad"
S155, S184, S185, S186 = (SCR/d for d in ["155","184","185","186"])
DISI_QM="/Users/ugursezen/Desktop/arin/deney/qm_riemann"
DISI_SCR=("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
          "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
def yukle(ad,yol,ek=()):
    k=Path(yol).read_text(encoding="utf-8").replace(DISI_QM,str(QM)).replace(DISI_SCR,str(SCR))
    for a,b in ek: k=k.replace(a,b)
    sp=importlib.util.spec_from_loader(ad,loader=None); m=importlib.util.module_from_spec(sp)
    m.__file__=str(yol); __import__("sys").modules[ad]=m
    exec(compile(k,str(yol),"exec"),m.__dict__); return m
b185=yukle("b185", QM/"185_configs"/"185b_oz_muhasebe.py")

kenar=json.load(open(S186/"ONKAYIT_186.json"))["kenar"]
print("bantlar:", kenar)
for gaz in ["gercek","Hkeskin"]:
    D=np.load(S184/f"K1_{gaz}.npz"); OZ=np.load(S185/f"OZ_{gaz}.npz")
    tau=np.asarray(D["tau"],float); ae=np.asarray(D["aq_eff"],float); aq=np.asarray(D["aq"],float)
    ahat=np.asarray(D["ahat"],float)
    mask=[(tau>=kenar[b])&(tau<kenar[b+1]) for b in range(8)]
    # c^öz (185b) ve Ĝ (186c tarifi)
    aq_full=aq
    coz=b185.c_oz(OZ, aq, -1)
    N=int(OZ["N"]); gbar=float(OZ["gbar"])
    Gre=OZ["Gre"].sum(0); Gim=OZ["Gim"].sum(0)
    Ghat=(Gre+1j*Gim)/N
    sat=[]
    for m in mask:
        wg=ahat[m].sum()/ae[m].sum()
        wog=np.abs(coz[m]).sum()/ae[m].sum()
        v= np.abs(Ghat[m]).sum()/(gbar*ae[m].sum())
        sat.append((wg,wog,wg-wog,v))
    sat=np.array(sat)
    print(f"\n[{gaz}]  {'bant':>10} {'wg':>7} {'wog':>7} {'mg':>7} {'v_b':>7} | {'√wg':>6} {'√mg':>6} {'1.017−0.884v':>13}")
    for k in range(8):
        wg,wog,mg,v=sat[k]
        print(f"        {kenar[k]:.2f}-{kenar[k+1]:.2f} {wg:7.4f} {wog:7.4f} {mg:7.4f} {v:7.4f} | "
              f"{np.sqrt(max(wg,0)):6.4f} {np.sqrt(max(mg,0)):6.4f} {1.017-0.884*v:13.4f}")
    for ad,idx in [("wg",0),("mg",2)]:
        y=np.sqrt(np.maximum(sat[:,idx],0)); v=sat[:,3]
        c=np.polyfit(v,y,1); r=np.corrcoef(v,y)[0,1]; rms=np.sqrt(np.mean((y-np.polyval(c,v))**2))
        print(f"    sqrt({ad}) = {c[1]:.3f} + {c[0]:.3f}·v   r={r:.3f}  RMS={rms:.4f}   (makale: 1.017 − 0.884 v)")
