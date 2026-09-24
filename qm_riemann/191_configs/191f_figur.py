# -*- coding: utf-8 -*-
"""
191f — FİGÜR: 191_derin_ikiz_130.png
======================================
Sol   : r_D(τ̄) dört derinlik (1.00,1.10,1.20,1.30; ± σ_r) + 1 − c·τ^α fitleri.
Orta  : üst f_D, alt g_D — Δz_D'ye karşı dört nokta (D=1.00 orijin) +
        orijinden-geçen yasa çizgisi + ∞ imaları (Δz_∞ = 0.2529).
Sağ   : üst σ_ε(D) (+ gerçek çizgi), alt c(D) ve α(D) (twin eksen).
"""
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S189, S191 = SCR / "189", SCR / "191"
HU = json.load(open(S191 / "HUKUM_191.json"))
ONK = json.load(open(S191 / "ONKAYIT_191.json"))
DER = {"1.00": "Hkeskin", "1.10": "Hderin110", "1.20": "Hderin120",
      "1.30": "Hderin130"}
RENK = {"1.00": "#2a78d6", "1.10": "#eb6834", "1.20": "#1baf7a",
       "1.30": "#c0392b"}
GERCEK = "#52514e"
Z = {D: json.load(open(S189 / f"zincir_{ad}.json")) for D, ad in DER.items()}
Dx = np.array([1.00, 1.10, 1.20, 1.30])

plt.rcParams.update({"axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": "0.9",
                     "grid.linewidth": 0.6})
fig = plt.figure(figsize=(16.5, 7.8))
gs = fig.add_gridspec(2, 3, width_ratios=[1.45, 1, 1], hspace=0.38, wspace=0.30)
fig.suptitle("191 — Derin ikiz τ≤1.30: zarf tam mı kapanır, gerçek kalıntı "
             "mı kalır? (Hkeskin → Hderin110 → Hderin120 → Hderin130)",
             fontsize=13)

# ---------------- sol: r_D(τ) ----------------
a = fig.add_subplot(gs[:, 0])
tt = np.linspace(0.44, 0.87, 200)
for D in DER:
    j = Z[D]
    tb = np.array(j["tau_bar"][:8])
    r = np.array(j["r"][:8])
    s = np.array(j["sig_r"][:8])
    a.errorbar(tb, r, yerr=s, fmt="-o", ms=6, lw=2, capsize=2,
               color=RENK[D], label=f"D={D} ({DER[D]})")
    f = HU["K3"]["fit"][D]
    a.plot(tt, 1 - f["c"] * tt ** f["alfa"], lw=0.9, color=RENK[D], alpha=0.6)
a.axhline(1.0, color="0.2", lw=1.0, ls=":")
a.text(0.555, 0.9975, "r = 1 (zarf kapalı)", fontsize=8, color="0.3", va="top")
a.set_ylim(0.865, 1.005)
a.set_xlabel(r"bant $\bar\tau$ (ae-ağırlıklı)")
a.set_ylabel(r"$r_D = w_g / w_{HD}$")
a.set_title("r_D(τ): dört derinlik ölçülen; ince çizgi: 1 − c·τ^α fiti")
a.legend(fontsize=8, loc="lower left")
ha = HU["H-191a"]
a.text(0.02, 0.985,
       "HAVUZ r: " + " → ".join(f"{ha['r_havuz'][D]:.4f}" for D in DER)
       + f"\nf$_{{1.10}}$={ha['f110_havuz']:.3f} f$_{{1.20}}$="
         f"{ha['f120_havuz']:.3f} f$_{{1.30}}$={ha['f130_havuz']:.3f}"
         f"±{ha['f130_se']:.3f}"
       + f"\n[H-191a {ha['hukum']}]",
       transform=a.transAxes, ha="left", va="top", fontsize=8,
       bbox=dict(fc="white", ec="0.7", alpha=0.92))

# ---------------- orta üst: f_D vs Δz_D ----------------
dl = HU["K3"]["dogrusallik"]
dz_pts = np.array([0.0, dl["Delta_z"]["1.10"], dl["Delta_z"]["1.20"],
                   dl["Delta_z"]["1.30"]])
dz_inf = dl["Delta_z"]["inf"]
f_pts = np.array([0.0, dl["f"]["1.10"], dl["f"]["1.20"], dl["f"]["1.30"]])
f_se_pts = np.array([0.0, dl["f_se"]["1.10"], dl["f_se"]["1.20"],
                     dl["f_se"]["1.30"]])
g_pts = np.array([0.0, dl["g"]["1.10"], dl["g"]["1.20"], dl["g"]["1.30"]])
g_se_pts = np.array([0.0, dl["g_se"]["1.10"], dl["g_se"]["1.20"],
                     dl["g_se"]["1.30"]])
xx = np.linspace(0, dz_inf * 1.06, 100)

a = fig.add_subplot(gs[0, 1])
for D, x, y, e in zip(DER, dz_pts, f_pts, f_se_pts):
    a.errorbar([x], [y], yerr=[e], fmt="o", ms=8, color=RENK[D], capsize=3,
               zorder=5)
a.plot(xx, dl["f_egim_orijin"] * xx, "-", lw=1.3, color="0.35",
       label=f"yasa: eğim={dl['f_egim_orijin']:.3f}±{dl['f_egim_orijin_se']:.3f}")
a.plot([dz_inf], [dl["f_inf"]], "*", ms=16, color="#8e44ad", zorder=6,
       label=f"f$_\\infty$={dl['f_inf']:.3f}±{dl['f_inf_se']:.3f}")
a.axvline(dz_inf, color="#8e44ad", lw=0.8, ls=":")
a.axhline(1.0, color="0.6", lw=0.8, ls=":")
a.set_xlabel(r"$\Delta z_D$ (gerçek-kinematik harita)")
a.set_ylabel(r"$f_D = (r_D-r_{1.00})/(1-r_{1.00})$")
a.set_title("f_D vs Δz_D — harita-doğrusal yasa + f$_\\infty$ ima")
a.legend(fontsize=7.5, loc="lower right")

# ---------------- orta alt: g_D vs Δz_D ----------------
a = fig.add_subplot(gs[1, 1])
for D, x, y, e in zip(DER, dz_pts, g_pts, g_se_pts):
    a.errorbar([x], [y], yerr=[e], fmt="o", ms=8, color=RENK[D], capsize=3,
               zorder=5)
a.plot(xx, dl["g_egim_orijin"] * xx, "-", lw=1.3, color="0.35",
       label=f"yasa: eğim={dl['g_egim_orijin']:.3f}±{dl['g_egim_orijin_se']:.3f}")
a.plot([dz_inf], [dl["g_inf"]], "*", ms=16, color="#8e44ad", zorder=6,
       label=f"g$_\\infty$={dl['g_inf']:.3f}±{dl['g_inf_se']:.3f}")
a.axvline(dz_inf, color="#8e44ad", lw=0.8, ls=":")
a.axhline(1.0, color="0.6", lw=0.8, ls=":")
a.set_xlabel(r"$\Delta z_D$ (gerçek-kinematik harita)")
a.set_ylabel(r"$g_D$ = σ-kapanış payı")
a.set_title("g_D vs Δz_D — σ_ε kapanış payı + g$_\\infty$ ima")
a.legend(fontsize=7.5, loc="lower right")

# ---------------- sağ üst: σ_ε(D) ----------------
hc = HU["H-191c"]
a = fig.add_subplot(gs[0, 2])
sg = [hc["sigma_eps"][D] for D in DER]
se = [hc["se"][D] for D in DER]
a.errorbar(Dx, sg, yerr=se, fmt="-o", ms=7, lw=2, color="0.3", capsize=3,
           label="ikiz σ_ε(D)")
for D, x, y in zip(DER, Dx, sg):
    a.plot(x, y, "o", ms=8, color=RENK[D], zorder=5)
a.axhline(hc["gercek"][0], color=GERCEK, lw=1.5, ls="--",
          label=f"gerçek {hc['gercek'][0]:.4f}")
a.axhspan(hc["bant"][0], hc["bant"][1], color="0.85", alpha=0.5, lw=0,
         label="ön-mühür bandı (1.30)")
a.set_xticks(Dx)
a.set_xlabel("merdiven derinliği D")
a.set_ylabel(r"$\sigma_\varepsilon$")
a.set_title(f"H-191c σ_ε yakınsaması [{hc['hukum']}]")
a.legend(fontsize=7.5)

# ---------------- sağ alt: c(D), α(D) ----------------
k3 = HU["K3"]["fit"]
a = fig.add_subplot(gs[1, 2])
c_v = [k3[D]["c"] for D in DER]
c_e = [k3[D]["jk_se"][0] for D in DER]
al_v = [k3[D]["alfa"] for D in DER]
al_e = [k3[D]["jk_se"][1] for D in DER]
a.errorbar(Dx, c_v, yerr=c_e, fmt="-o", ms=7, lw=2, color="#2a78d6",
           capsize=3, label="c(D)")
a.set_xlabel("merdiven derinliği D")
a.set_ylabel("c", color="#2a78d6")
a.tick_params(axis="y", labelcolor="#2a78d6")
a.set_xticks(Dx)
b = a.twinx()
b.errorbar(Dx, al_v, yerr=al_e, fmt="-s", ms=7, lw=2, color="#c0392b",
           capsize=3, label="α(D)")
b.set_ylabel("α", color="#c0392b")
b.tick_params(axis="y", labelcolor="#c0392b")
b.grid(False)
a.set_title("zarf fiti 1 − c·τ^α: c(D) ve α(D) (± jk)")
lines1, labels1 = a.get_legend_handles_labels()
lines2, labels2 = b.get_legend_handles_labels()
a.legend(lines1 + lines2, labels1 + labels2, fontsize=8, loc="center left")

out = QM / "191_derin_ikiz_130.png"
fig.savefig(out, dpi=140, bbox_inches="tight")
print(f"-> {out}")
