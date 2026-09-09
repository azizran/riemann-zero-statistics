# -*- coding: utf-8 -*-
"""
186g — FİGÜR: 186_M_turetim.png
===============================
Sol: ölçülü M(τ) vs G1/G2 öngörüleri (+ G3 güç fiti).
Sağ: zincir r_pred vs 184 yasası (+ sabit-nokta r*).
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S186 = SCR / "186"

G1S = json.load(open(S186 / "G1_sonuc.json"))
G2S = json.load(open(S186 / "G2_sonuc.json"))
K3S = json.load(open(S186 / "K3_sonuc.json"))
G3S = json.load(open(S186 / "G3_bicim.json"))

tb = np.array(G1S["tau_bar"])
M = np.array(G1S["M_olc"])
sigM = np.array(G1S["sigM_defter"])
tt = np.linspace(0.45, 0.86, 200)

fig, ax = plt.subplots(1, 2, figsize=(12.5, 5.2))
fig.suptitle("186 — M(τ) karışım-iletiminin türetim sınavı: G1 öz-tutarlılık / "
             "G2 v-köprüsü İKİSİ DE ÖLDÜ (ön-kayıt 23a8130a)", fontsize=11)

a = ax[0]
a.axhline(1.0, color="0.75", lw=0.8, zorder=0)
a.errorbar(tb, M, yerr=sigM, fmt="k^", ms=7, lw=1.5, capsize=3,
           label="M ölçülü (185: m_g/m_Hk)", zorder=6)
k, be = G3S["guc"]["k"], G3S["guc"]["beta"]
a.plot(tt, 1 - k * tt ** be, "k--", lw=1,
       label=f"G3 güç fiti 1−{k:.3f}·τ^{be:.2f} (χ²/dof={G3S['guc']['chi2dof']:.2f})")
g1 = np.array(G1S["g_yasa"]["M_pred"])
g1se = np.array(G1S["g_yasa"]["se_M_pred"])
a.errorbar(tb, g1, yerr=g1se, fmt="o--", color="#c0392b", ms=5, capsize=2,
           label=f"G1 M^pred (ρ=184 yasası; χ²/dof={G1S['g_yasa']['chi2dof_M']:.0f}) ÖLDÜ")
a.plot(tb, G1S["g_bir"]["M_pred"], ":", color="#7f8c8d", lw=1.5,
       label="G1 kontrol ρ≡1 (makine mührü: >1 ✓)")
a.plot(tb, G1S["h_yasa"]["M_pred"], "d:", color="#e67e22", ms=4,
       label=f"G1 yan-lehçe Hk-çekirdek (χ²/dof={G1S['h_yasa']['chi2dof_M']:.0f})")
g2 = np.array(G2S["M_pred"])
a.errorbar(tb, g2, yerr=np.array(G2S["M_pred_se"]), fmt="s--", color="#2980b9",
           ms=5, capsize=2,
           label=f"G2 v-köprüsü (χ²/dof={G2S['chi2dof_M']:.0f}) ÖLDÜ")
a.set_xlabel("τ")
a.set_ylabel("M(τ) = m_g/m_Hk")
a.set_title("Karışım-iletimi: öngörüler 1'in ÜSTÜNDE, ölçüm 0.80→0.72")
a.legend(fontsize=7.5, loc="center left")
a.grid(alpha=0.25)

b = ax[1]
b.axhline(1.0, color="0.75", lw=0.8, zorder=0)
r = np.array(K3S["r_olc"])
sr = np.array(K3S["sig_r"])
b.errorbar(tb, r, yerr=sr, fmt="k^", ms=7, lw=1.5, capsize=3,
           label="r ölçülü (184 defteri)", zorder=6)
b.plot(tt, 1 - 0.149 * tt ** 1.30, "k-", lw=1.2,
       label="184: 1−0.149·τ^1.30 (χ²/dof=0.34)")
rp = np.array(K3S["r_pred"])
b.errorbar(tb, rp, yerr=np.array(K3S["se_r_pred"]), fmt="o--", color="#c0392b",
           ms=5, capsize=2,
           label=f"zincir r_pred (G1 teşhis; χ²/dof={K3S['zincir_chi2dof']:.0f}) "
                 f"KAPANMADI")
c_p, a_p, _ = K3S["fit"]["zincir"]
b.plot(tt, 1 - c_p * tt ** a_p, "--", color="#c0392b", alpha=0.4, lw=1,
       label=f"fit: 1−({c_p:+.3f})·τ^{a_p:.2f}")
b.plot(tb, K3S["sabit_nokta"]["rho_star"], "x", color="#8e44ad", ms=7, mew=2,
       label="sabit nokta r* (doğrusallaştırma) — anti-zarf")
b.set_xlabel("τ")
b.set_ylabel("r(τ)")
b.set_title("Zincir: öz-tutarlılık sabit noktası r*>1 — zarf öz-üretmiyor")
b.legend(fontsize=7.5, loc="center left")
b.grid(alpha=0.25)

fig.tight_layout(rect=[0, 0, 1, 0.95])
yol = QM / "186_M_turetim.png"
fig.savefig(yol, dpi=150)
print(f"-> {yol}  BİTTİ")
