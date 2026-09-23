# -*- coding: utf-8 -*-
"""
189e — FİGÜR: 189_derin_ikiz.png
=================================
Sol  : r_D(τ̄) üç derinlik (1.00, 1.10, 1.20; ± σ_r 184 konv.) + ön-mühür
       öngörüleri (kesikli, açık işaret) + 1 − c·τ^α fitleri (ince) + r = 1.
Orta : üst σ_ε(D), alt |ζ_HD|(HAVUZ) ve harita z_D; gerçek değerler yatay çizgi.
Sağ  : üst c(D), alt α(D) (± jk; 184c makinesi).
Tek eksen kuralı: farklı ölçekli nicelikler ayrı panellerde.
"""
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
S189 = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
            "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/189")
HU = json.load(open(S189 / "HUKUM_189.json"))
ONK = json.load(open(S189 / "ONKAYIT_189.json"))
DER = {"1.00": "Hkeskin", "1.10": "Hderin110", "1.20": "Hderin120"}
RENK = {"1.00": "#2a78d6", "1.10": "#eb6834", "1.20": "#1baf7a"}
GERCEK = "#52514e"
Z = {D: json.load(open(S189 / f"zincir_{ad}.json")) for D, ad in DER.items()}
ong = ONK["ongoru"]["tablo_kalem_AYNEN"]
BANT = list(ong)
Dx = np.array([1.00, 1.10, 1.20])

plt.rcParams.update({"axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": "0.9", "grid.linewidth": 0.6})
fig = plt.figure(figsize=(15.5, 7.6))
gs = fig.add_gridspec(2, 3, width_ratios=[1.45, 1, 1], hspace=0.35, wspace=0.28)
fig.suptitle("189 — Derin ikiz: zarfın kimliği sınavı "
             "(Hkeskin τ≤1.00 → Hderin110 → Hderin120)", fontsize=13)

# ---------------- sol: r_D(τ) ----------------
a = fig.add_subplot(gs[:, 0])
tt = np.linspace(0.44, 0.87, 200)
for D in DER:
    j = Z[D]
    tb = np.array(j["tau_bar"][:8])
    r = np.array(j["r"][:8])
    s = np.array(j["sig_r"][:8])
    a.errorbar(tb, r, yerr=s, fmt="-o", ms=6, lw=2, capsize=2,
               color=RENK[D], label=f"ölçülen D={D} ({DER[D]})")
    f = HU["K3"]["fit"][D]
    a.plot(tt, 1 - f["c"] * tt ** f["alfa"], lw=0.9, color=RENK[D], alpha=0.6)
    if D != "1.00":
        col = 4 if D == "1.10" else 5
        a.plot(tb, [ong[b][col] for b in BANT], "--", marker="o", mfc="white",
               ms=6, lw=1.5, color=RENK[D], label=f"ön-mühür öngörü D={D}")
a.axhline(1.0, color="0.2", lw=1.0, ls=":")
a.text(0.555, 0.9965, "r = 1 (zarf kapalı)", fontsize=8, color="0.3", va="top")
a.set_ylim(0.868, 1.038)
a.set_xlabel(r"bant $\bar\tau$ (ae-ağırlıklı)")
a.set_ylabel(r"$r_D = w_g / w_{HD}$")
a.set_title("r_D(τ): ölçülen (dolu) vs ön-mühür (kesikli); ince: 1 − c·τ^α")
a.legend(fontsize=8, loc="lower left")
ha = HU["H-189a"]
if "r_havuz" in ha:
    a.text(0.02, 0.985,
           "HAVUZ r: " + " → ".join(f"{ha['r_havuz'][D]:.4f}" for D in DER)
           + f"\nf$_{{1.10}}$ = {ha['f110_havuz']:.3f}, "
             f"f$_{{1.20}}$ = {ha['f120_havuz']:.3f}   [H-189a {ha['hukum']}]",
           transform=a.transAxes, ha="left", va="top", fontsize=8.5,
           bbox=dict(fc="white", ec="0.7", alpha=0.9))

# ---------------- orta üst: σ_ε ----------------
hc = HU["H-189c"]
a = fig.add_subplot(gs[0, 1])
sg = [hc["sigma_eps"][D] for D in DER]
se = [hc["se"][D] for D in DER]
a.errorbar(Dx, sg, yerr=se, fmt="-o", ms=7, lw=2, color="0.3",
           capsize=3, label="ikiz σ_ε(D) (± jk)")
for D, x, y in zip(DER, Dx, sg):
    a.plot(x, y, "o", ms=8, color=RENK[D], zorder=5)
a.axhline(hc["gercek"][0], color=GERCEK, lw=1.5, ls="--",
          label=f"gerçek {hc['gercek'][0]:.4f}")
a.set_xticks(Dx)
a.set_xlabel("merdiven derinliği D (τ_üst)")
a.set_ylabel(r"$\sigma_\varepsilon = \mathrm{std}(g)/\bar g$")
a.set_title(f"H-189c geri-besleme kinematiği [{hc['hukum']}]")
a.legend(fontsize=8)

# ---------------- orta alt: ζ_HD ----------------
hd = HU["H-189d"]
a = fig.add_subplot(gs[1, 1])
zm = [Z[D]["zeta_mod"][8] for D in DER]
zs = [Z[D]["se_mod"][8] for D in DER]
a.errorbar(Dx, zm, yerr=zs, fmt="-o", ms=7, lw=2, color="0.3",
           capsize=3, label="|ζ_HD| ölçülen (187c, HAVUZ)")
for D, x, y in zip(DER, Dx, zm):
    a.plot(x, y, "o", ms=8, color=RENK[D], zorder=5)
zD = [ONK["ongoru"]["HAVUZ_bilgi"][f"z_{D}"] for D in DER]
a.plot(Dx, zD, "--", marker="s", mfc="white", ms=7, color="0.35",
       label="harita z_D (188, öngörü)")
for x, z in zip(Dx[1:], zD[1:]):
    a.fill_between([x - 0.02, x + 0.02], z * 0.9, z * 1.1, color="0.85",
                   alpha=0.7, lw=0, zorder=0)
zg = float(np.load(S189 / "zincir_gercek.npz")["z_tam"][8, 0])
a.axhline(zg, color=GERCEK, lw=1.5, ls="--",
          label=f"gerçek |ζ_g| = {zg:.4f} (tam derinlik)")
a.set_xticks(Dx)
a.set_xlabel("merdiven derinliği D (τ_üst)")
a.set_ylabel("iptal payı |ζ| (HAVUZ)")
a.set_title(f"H-189d harita geçerliliği [{hd['hukum']}] (gri: ±%10)")
a.set_ylim(0.04, 0.36)
a.legend(fontsize=8, loc="upper left", bbox_to_anchor=(0.0, 0.87))

# ---------------- sağ: c, α ----------------
k3 = HU["K3"]["fit"]
for row, (key, idx, lab) in enumerate((("c", 0, "c"), ("alfa", 1, "α"))):
    a = fig.add_subplot(gs[row, 2])
    v = [k3[D][key] for D in DER]
    e = [k3[D]["jk_se"][idx] for D in DER]
    a.errorbar(Dx, v, yerr=e, fmt="-o", ms=7, lw=2, color="0.3", capsize=3)
    for D, x, y in zip(DER, Dx, v):
        a.plot(x, y, "o", ms=8, color=RENK[D], zorder=5)
        if k3[D]["sinirda"]:
            a.annotate("sınırda", (x, y), xytext=(4, 6),
                       textcoords="offset points", fontsize=8)
    a.set_xticks(Dx)
    a.set_xlabel("merdiven derinliği D (τ_üst)")
    a.set_ylabel(lab)
    a.set_title(f"zarf fiti 1 − c·τ^α: {lab}(D) (± jk)")
    if key == "c":
        a.axhline(0, color="0.2", lw=0.8, ls=":")

out = QM / "189_derin_ikiz.png"
fig.savefig(out, dpi=140, bbox_inches="tight")
print(f"-> {out}")
