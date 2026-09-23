# -*- coding: utf-8 -*-
"""
188e — FİGÜR: 188_zeta_kimligi.png
  sol üst : |K|(τ_b, τ') ısı haritası (41 ince bant × 88 dilim) + uydu/Bragg
            dikey çizgileri + köşegenler τ'=2−τ_b, τ'=τ_b
  sağ üst : |v(τ')| ile S(τ'), W(τ') (her biri kendi maksimumuna normalize)
  alt     : 8-bant |ζ| ölçülen (187c, tam derinlik) vs ζ^{≤1.30} yeniden
            kurulumu (z = −Re ζ: toplam, rang-1 bileşeni, artık)
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
S188 = SCR / "188"

C1, C2, C3 = "#2a78d6", "#eb6834", "#1baf7a"   # mavi / turuncu / su yeşili
GRI, MUREKKEP, IKINCIL = "#8a8985", "#0b0b0b", "#52514e"

H = np.load(S188 / "harita_K_gercek.npz")
A = json.load(open(S188 / "K3_analiz.json"))
ONK = json.load(open(S188 / "ONKAYIT_188.json"))
K = H["K"][:41]
kenar = H["kenar"]
ik = H["ince_kenar"]
orta = np.array(A["orta"])
tb = np.array(A["tb"])
tp = ONK["uydular"]["tau_p"]

plt.rcParams.update({"font.size": 9, "axes.edgecolor": GRI,
                     "axes.labelcolor": MUREKKEP, "xtick.color": IKINCIL,
                     "ytick.color": IKINCIL, "axes.spines.top": False,
                     "axes.spines.right": False})
fig = plt.figure(figsize=(12.5, 9.0), facecolor="#fcfcfb")
gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1.0], hspace=0.32,
                      wspace=0.22)

# --- sol üst: ısı haritası ---
ax = fig.add_subplot(gs[0, 0])
im = ax.pcolormesh(kenar, ik, np.abs(K), cmap="Blues", shading="flat")
cb = fig.colorbar(im, ax=ax, pad=0.01)
cb.set_label("|K(τ_b, τ')|")
for p in ("2", "3", "5", "7"):
    for sg in (1, -1):
        t = 1 + sg * tp[p]
        if kenar[0] < t < kenar[-1]:
            ax.axvline(t, color=C2, lw=0.8, ls="--", alpha=0.9)
            ax.text(t, ik[-1] + 0.003, f"1{'+' if sg > 0 else '−'}τ{p}",
                    rotation=90, fontsize=7, color=IKINCIL, ha="center",
                    va="bottom")
ax.axvline(1.0, color=MUREKKEP, lw=0.9, ls=":")
ax.text(1.0, ik[-1] + 0.003, "Bragg", rotation=90, fontsize=7,
        color=IKINCIL, ha="center", va="bottom")
yy = np.linspace(ik[0], ik[-1], 50)
ax.plot(2 - yy, yy, color=C3, lw=1.2, label="τ' = 2 − τ_b")
ax.set_xlim(kenar[0], kenar[-1])
ax.set_ylim(ik[0], ik[-1])
ax.set_xlabel("pencere-ötesi dilim τ'")
ax.set_ylabel("ince bant τ_b")
ax.set_title("İptal çekirdeği haritası |K| (gerçek; 188b)", loc="left",
             fontsize=10, color=MUREKKEP, pad=26)
ax.legend(loc="lower left", fontsize=7, frameon=False)

# --- sağ üst: |v| vs S, W ---
ax = fig.add_subplot(gs[0, 1])
v = np.array(A["v_mod"])
S = np.array(A["S"])
W = np.array(A["W"])
ax.plot(orta, v / v.max(), color=C1, lw=2, label="|v(τ')| (rang-1)")
ax.plot(orta, S / S.max(), color=C2, lw=2, label="S = Σ a|Ĝ_mid|")
ax.plot(orta, W / W.max(), color=C3, lw=2, label="W = Σ a")
SRe = -np.array(A["SRe"])
ax.plot(orta, SRe / SRe.max(), color=MUREKKEP, lw=1.2, ls="--",
        label="İKİNCİL −S_Re = −Σ a πτ' cos(πτ') Re Ĝ")
for p in ("2", "3", "5", "7"):
    for sg in (1, -1):
        t = 1 + sg * tp[p]
        if kenar[0] < t < kenar[-1]:
            ax.axvline(t, color=GRI, lw=0.6, ls="--")
ax.axvline(1.0, color=GRI, lw=0.8, ls=":")
ax.set_xlim(kenar[0], kenar[-1])
ax.set_ylim(0, 1.05)
ax.set_xlabel("τ'")
ax.set_ylabel("kendi maksimumuna normalize")
ax.grid(axis="y", color="#e6e5e0", lw=0.6)
ax.set_title(f"corr(|v|,S) = {A['ii']['corr_vS']:.3f}   "
             f"corr(|v|,W) = {A['ii']['corr_vW']:.3f}\n"
             f"[İKİNCİL corr(Re v, S_Re) = {A['ii']['ikincil_corr_RevSRe']:.3f}]"
             "  (188d)", loc="left", fontsize=10, color=MUREKKEP, pad=26)
ax.legend(loc="upper left", fontsize=8, frameon=False)

# --- alt: 8-bant yeniden kurulum ---
ax = fig.add_subplot(gs[1, :])
k3 = A["vi"]
x = np.arange(8)
et = [r["bant"] for r in k3]
zf = np.array([r["zeta_tam_187c"] for r in k3])
zfs = np.array([r["se_zeta_tam"] for r in k3])
z = np.array([r["z"] for r in k3])
zs = np.array([r["se_z"] for r in k3])
zr = np.array([r["z_r1"] for r in k3])
zrs = np.array([r["se_z_r1"] for r in k3])
za = np.array([r["z_art"] for r in k3])
zas = np.array([r["se_z_art"] for r in k3])
ax.errorbar(x, zf, zfs, color=MUREKKEP, marker="o", ms=6, lw=2,
            capsize=3, label="|ζ| ölçülen, tam derinlik (187c)")
ax.errorbar(x, z, zs, color=C1, marker="s", ms=6, lw=2, capsize=3,
            label="−Re ζ^{≤1.30} harita toplamı")
ax.errorbar(x, zr, zrs, color=C2, marker="^", ms=6, lw=2, capsize=3,
            label="rang-1 bileşeni (u·Σv)")
ax.errorbar(x, za, zas, color=C3, marker="v", ms=6, lw=2, capsize=3,
            label="rang-1 artığı")
ax.axhline(0, color=GRI, lw=0.6)
ax.set_xticks(x)
ax.set_xticklabels(et)
ax.set_xlabel("bant τ_b (184'ün 8 bandı)")
ax.set_ylabel("iptal payı")
ax.grid(axis="y", color="#e6e5e0", lw=0.6)
ax.set_title("8-bant ζ profili: ölçülen vs harita yeniden kurulumu "
             "(±jk se; 188d)", loc="left", fontsize=10, color=MUREKKEP)
ax.legend(loc="upper left", fontsize=8, frameon=False, ncol=2)

fig.suptitle("188 — ζ(τ)'nin kimliği: iptal çekirdeği haritası",
             x=0.01, ha="left", fontsize=12, color=MUREKKEP)
out = QM / "188_zeta_kimligi.png"
fig.savefig(out, dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
print(f"-> {out}")
