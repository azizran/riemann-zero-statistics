# -*- coding: utf-8 -*-
"""
187h — FİGÜR: 187_faz_iptal.png
================================
Sol üst: yakınsama defteri m^ya(τ_c)/m_ölç (bantlar + HAVUZ) → 1.
Sol alt: katman artımları / iptal-yoğunluğu ι/Δτ' (HAVUZ) + biçim kayıtları.
Sağ üst: ζ açı defteri (gerçek + ikiz; 180°±15° bandı) + |ζ| profili.
Sağ alt: M_HA4 öngörü-vs-ölçüm (K4).
"""
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
S187 = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
            "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/187")

K1 = json.load(open(S187 / "K1_katman_defteri.json"))
K2 = json.load(open(S187 / "K2_zeta.json"))
K3 = json.load(open(S187 / "K3_yapi.json"))
K4 = json.load(open(S187 / "K4_HA4.json"))
KK = json.load(open(S187 / "K1k_ikiz_sagirlik.json"))

izg = np.array(K1["izgara"])
mya = np.array(K1["m_ya"])          # (9, 8)
m_olc = np.array(K1["m_olc"])
tau_bar = np.array(K1["tau_bar"])
et = K1["etiket"]

fig, ax = plt.subplots(2, 2, figsize=(13.5, 9.5))
fig.suptitle("187 — Pencere-ötesi faz-iptal defteri "
             "(katmanlar, ζ açısı, ikiz-sağırlığı, HA4)", fontsize=13)

# --- sol üst: yakınsama defteri ---
a = ax[0, 0]
cmap = plt.cm.viridis(np.linspace(0.15, 0.9, 8))
for k in range(8):
    a.plot(izg, mya[k] / m_olc[k], "-o", ms=3, lw=1, color=cmap[k],
           alpha=0.7, label=f"τ̄={tau_bar[k]:.2f}")
a.plot(izg, mya[8] / m_olc[8], "-o", ms=5, lw=2.5, color="crimson",
       label="HAVUZ")
a.axhline(1.0, color="k", lw=0.8, ls="--")
a.set_xlabel(r"katman derinliği $\tau_c$")
a.set_ylabel(r"$m^{ya}(\tau_c)\,/\,m_{ölç}$")
a.set_title("K1 — yakınsama defteri (H-F1)")
a.legend(fontsize=7, ncol=2)
kal = K1["HF1"]["kalan_acik_pay"]
a.text(0.03, 0.05, f"HAVUZ: kapanan pay = {K1['HF1']['kapanan_pay']:.2f}, "
       f"kalan açık = {kal:.2f}\nkorr(ds$^{{ya}}$,ds): "
       f"{K1['korr'][0]:.3f}$\\to${K1['korr'][-1]:.3f}",
       transform=a.transAxes, fontsize=8,
       bbox=dict(fc="white", alpha=0.8, ec="0.7"))

# --- sol alt: iptal-yoğunluğu ---
a = ax[1, 0]
orta = 0.5 * (izg[:-1] + izg[1:])
dtau = np.diff(izg)
iota_h = np.array(K3["iota"])[8]
se_iota_h = np.array(K3["se_iota"])[8]
a.bar(orta, iota_h / dtau, width=dtau * 0.9, color="steelblue", alpha=0.7,
      yerr=se_iota_h / dtau, capsize=3, label=r"$\iota/\Delta\tau'$ (HAVUZ)")
bu = K3["bicim_iota_ustel"]
xx = np.linspace(0.86, 1.21, 200)
a.plot(xx, bu["A"] * np.exp(-(xx - 0.86) / bu["lambda"]) /
       np.interp(xx, orta, dtau), "r--", lw=1.5,
       label=(f"üstel kayıt: λ={bu['lambda']:.3f} "
              f"(χ²/dof={bu['chi2dof']:.1f})"))
a.set_xlabel(r"katman $\tau'$")
a.set_ylabel(r"iptal yoğunluğu $\iota/\Delta\tau'$")
a.set_title("K3 — iptal-yoğunluk profili ι(τ') (KAYIT)")
a.legend(fontsize=8)

# --- sağ üst: ζ açı defteri ---
a = ax[0, 1]
tb = tau_bar[:8]
for gaz, renk, mk in [("gercek", "crimson", "o"), ("Hkeskin", "gray", "s")]:
    aci = np.array(K2[gaz]["zeta_aci"])[:8]
    se_aci = np.array(K2[gaz]["se_aci"])[:8]
    aci_g = np.where(aci < 0, aci + 360, aci)   # 180 çevresine getir
    a.errorbar(tb, aci_g, yerr=se_aci, fmt=mk, color=renk, ms=5,
               capsize=3, label=f"{gaz} açı")
a.axhspan(165, 195, color="green", alpha=0.12, label="H-F2 penceresi 180°±15°")
a.axhline(180, color="k", lw=0.8, ls="--")
a.set_ylim(150, 210)
a.set_xlabel(r"bant $\bar\tau$")
a.set_ylabel("ζ açısı (derece)")
a.set_title("K2 — iptal açısı (H-F2)")
a2 = a.twinx()
for gaz, renk in [("gercek", "crimson"), ("Hkeskin", "gray")]:
    a2.plot(tb, np.array(K2[gaz]["zeta_mod"])[:8], "-", color=renk,
            alpha=0.4, lw=1.2)
a2.set_ylabel(r"$|\zeta|$ (soluk çizgiler)", fontsize=8)
a2.set_ylim(0, 0.45)
a.legend(fontsize=8, loc="lower left")

# --- sağ alt: M_HA4 ---
a = ax[1, 1]
M_olc = np.array(K4["M_olc"])[:8]
sigM = np.array(K4["sigM"])[:8]
Mp = np.array(K4["M_pred_086"])[:8]
seMp = np.array(K4["se_M_pred"])[:8]
Mpt = np.array(K4["M_pred_110tam"])[:8]
a.errorbar(tb, M_olc, yerr=sigM, fmt="ko", ms=5, capsize=3,
           label="M_HA4 ölçülü (186e)")
a.errorbar(tb, Mp, yerr=seMp, fmt="^", color="crimson", ms=6, capsize=3,
           label="öngörü: iptalsiz pencere-içi erfc (0.86)")
a.plot(tb, Mpt, "x", color="steelblue", ms=6,
       label="yan: tam-env (1.10)")
a.axhline(1.0, color="k", lw=0.8, ls="--")
a.set_xlabel(r"bant $\bar\tau$")
a.set_ylabel(r"$M_{HA4}$")
a.set_title(f"K4 — HA4 kapanışı (H-F4): χ²/dof = "
            f"{K4['chi2dof_birincil']:.2f}")
a.legend(fontsize=8)

plt.tight_layout(rect=[0, 0, 1, 0.965])
plt.savefig(QM / "187_faz_iptal.png", dpi=150)
print(f"-> {QM/'187_faz_iptal.png'}")
