# -*- coding: utf-8 -*-
"""
185f — FİGÜR: 185_turetim.png
=============================
Sol: üç faktörün τ-profili (+ karışım-oranı M(τ), K2-yorum türevi).
Sağ: r ölçülü vs r_pred^A/B + 184 güç-yasası + türetilmiş güç-yasaları.
Ayrıca M(τ) defterini (se yayılımıyla) ekrana basar (rapor için).
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
S185 = SCR / "185"

K1 = np.load(S185 / "K1_faktorler.npz", allow_pickle=True)
K3 = json.load(open(S185 / "K3_kapanis.json"))
tb = K1["tau_bar"][:8]
r, sig = K1["r"][:8], K1["sig_r"][:8]
F1, F2, F3i = K1["F1"][:8], K1["F2"][:8], K1["F3i"][:8]
wg, wh, wog, woh = K1["wg"][:8], K1["wh"][:8], K1["wog"][:8], K1["woh"][:8]
se = K1["se_tablo"]  # (9,11): sütun 7..10 = wg, wh, wog, woh

# --- karışım-oranı M(τ) = (w_g − w_g^öz)/(w_Hk − w_Hk^öz)  (K2-yorum türevi) ---
mg = wg - wog
mh = wh - woh
M = mg / mh
seM = M * np.sqrt((np.sqrt(se[:8, 7] ** 2 + se[:8, 9] ** 2) / mg) ** 2
                  + (np.sqrt(se[:8, 8] ** 2 + se[:8, 10] ** 2) / mh) ** 2)
print("KARIŞIM-ORANI DEFTERİ (K2-yorum türevi; hizalı-faz yaklaşımı):")
print(f"{'τ̄':>6} {'m_g=w_g−w_g^öz':>14} {'m_Hk':>8} {'M=m_g/m_Hk':>12}")
for k in range(8):
    print(f"{tb[k]:6.3f} {mg[k]:14.4f} {mh[k]:8.4f} {M[k]:8.4f}±{seM[k]:.4f}")

A = np.array(K3["r_pred_A"])
B = np.array(K3["r_pred_B"])
cA, aA, _ = K3["fit"]["A"]
cB, aB, _ = K3["fit"]["B"]
cO, aO, _ = K3["fit"]["olculu"]
tt = np.linspace(0.45, 0.86, 200)

fig, ax = plt.subplots(1, 2, figsize=(12.5, 5.2))
fig.suptitle("185 — Zarfın türetim muhasebesi: r = F_ANOM × F_KİN × F_KOMŞU⁻¹"
             "  (ön-kayıt d49a52a4)", fontsize=11)

a = ax[0]
a.axhline(1.0, color="0.75", lw=0.8, zorder=0)
a.plot(tb, F1, "o-", color="#c0392b", label="F_ANOMALİ = â_g/â_g^öz")
a.plot(tb, F2, "s-", color="#27ae60", label="F_KİNEMATİK = â_g^öz/â_Hk^öz")
a.plot(tb, F3i, "d-", color="#2980b9", label="F_KOMŞU⁻¹ = â_Hk^öz/â_Hk")
a.plot(tb, r, "k^-", lw=2, label="r = çarpım (ölçülü zarf)")
a.plot(tb, M, ":", color="#8e44ad", lw=2,
       label="M(τ) = karışım-oranı (yorum)")
a.set_xlabel("τ")
a.set_ylabel("faktör")
a.set_title("Üç faktörün τ-profili — zarfı ANOM×KOMŞU⁻¹ çifti taşıyor")
a.legend(fontsize=8, loc="center left")
a.grid(alpha=0.25)

b = ax[1]
b.errorbar(tb, r, yerr=sig, fmt="k^", ms=7, lw=1.5, capsize=3,
           label="r ölçülü (184 defteri)", zorder=5)
b.plot(tt, 1 - 0.149 * tt ** 1.30, "k-", lw=1.2,
       label="184: 1−0.149·τ^1.30 (χ²/dof=0.34)")
b.plot(tb, A, "s--", color="#27ae60",
       label=f"r_pred^A kinematik (χ²/dof={K3['chi2dof_A']:.0f}) ÖLDÜ")
b.plot(tt, 1 - cA * tt ** aA, "--", color="#27ae60", alpha=0.4, lw=1)
b.plot(tb, B, "o--", color="#e67e22",
       label=f"r_pred^B ortak-karışım (χ²/dof={K3['chi2dof_B']:.0f}) ÖLDÜ")
b.plot(tt, 1 - cB * tt ** aB, "--", color="#e67e22", alpha=0.4, lw=1)
b.axhline(1.0, color="0.75", lw=0.8, zorder=0)
b.set_xlabel("τ")
b.set_ylabel("r(τ)")
b.set_title("İleri-hesaplar zarfın TERS yönünde → zarf = anomali/karışım kanalı")
b.legend(fontsize=8, loc="center left")
b.grid(alpha=0.25)

fig.tight_layout(rect=[0, 0, 1, 0.95])
yol = QM / "185_turetim.png"
fig.savefig(yol, dpi=150)
print(f"\n-> {yol}  BİTTİ")
