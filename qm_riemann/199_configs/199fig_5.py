# -*- coding: utf-8 -*-
"""
199fig_5 — NOT 7 FİGÜRÜ 5: fig_n7_height_en.{png,pdf}
=======================================================
Mühürlü 198 sonucunu (H-198a H_S TUTAR) gösterir; YENİ ÖLÇÜM YOK. Yalnız json okunur.
κ_p(W) ± 8-grup jk se, üç pencere (W_alt L = 9.343 KÖR; W_düşük 10.484 referans; W_son 12.030
yarı-kör) — HUKUM_198.json KAYIT.pencereler.
Eğriler (ön-kayıtlı model 1 − κ_p(L) = B_p (L/L_d)^{−γ}, L_d = 10.484):
  düz      γ̂ uyumu: γ̂ = HUKUM_198 gamma.deger, B_p = gamma.B_p (profillenmiş)
  kesikli  H_S (γ = 1), B_p = 1 − κ_p^{197}   (ön-kayıtlı öngörü tablosunun çapası)
  noktalı  H_C (γ = 0): κ_p = κ_p^{197}
κ_p^{197} = A(+log p)/(s·c_BK(p)), HUKUM_197 KAYIT İKİNCİL (197 raporu 0.795/0.720/0.608/0.530);
ön-kayıttaki öngörü tablosu (H_S W_alt 0.770/0.686/0.560/0.473, W_son 0.821/0.756/0.659/0.591)
burada yeniden üretilip assert edilir.
Girdi: scratchpad/198/HUKUM_198.json, scratchpad/197/HUKUM_197.json
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
S197, S198 = QM / "scratchpad" / "197", QM / "scratchpad" / "198"
CIKTI = QM / "fig_n7_height_en"

MAVI, KIZIL, YESIL, MOR = "#0072B2", "#D55E00", "#009E73", "#AA3377"
INK, GRI = "#1a1a1a", "#8c8c8c"
plt.rcParams.update({"font.size": 9, "axes.labelsize": 10, "xtick.labelsize": 9,
                     "ytick.labelsize": 9, "legend.fontsize": 9, "axes.linewidth": 0.8,
                     "mathtext.fontset": "dejavusans", "savefig.dpi": 300,
                     "pdf.fonttype": 42, "ps.fonttype": 42})

H = json.load(open(S198 / "HUKUM_198.json"))
H7 = json.load(open(S197 / "HUKUM_197.json"))
PRIMES = (2, 3, 5, 7)
PENC = ("alt", "dusuk", "son")
L_D = float(H["KAYIT"]["pencereler"]["dusuk"]["L_W"])
LW = np.array([H["KAYIT"]["pencereler"][W]["L_W"] for W in PENC])
K = {p: np.array([H["KAYIT"]["pencereler"][W]["kappa"][str(p)]["deger"] for W in PENC])
     for p in PRIMES}
E = {p: np.array([H["KAYIT"]["pencereler"][W]["kappa"][str(p)]["se"] for W in PENC])
     for p in PRIMES}
g0 = H["gamma"]["deger"]
sg = H["gamma"]["sigma_eff"]
Bp = {p: H["gamma"]["B_p"][str(p)] for p in PRIMES}


def c_bk(p):
    return p / (p - 1) ** 2


ik = H7["KAYIT"]["IKINCIL_arti"]["konumlar_192"]
s7 = H7["nicelikler"]["s"]["deger"]
K197 = {p: ik[f"+log{p}"]["A"]["deger"] / (s7 * c_bk(p)) for p in PRIMES}

# ---- denetim ----
HS_alt = {p: 1 - (1 - K197[p]) * L_D / LW[0] for p in PRIMES}
HS_son = {p: 1 - (1 - K197[p]) * L_D / LW[2] for p in PRIMES}
assert np.allclose([K197[p] for p in PRIMES], [0.795, 0.720, 0.608, 0.530], atol=5e-4)
assert np.allclose([HS_alt[p] for p in PRIMES], [0.770, 0.686, 0.560, 0.473], atol=1.5e-3)  # kalem tablosu yuvarlak κ ile
assert np.allclose([HS_son[p] for p in PRIMES], [0.821, 0.756, 0.659, 0.591], atol=1.5e-3)
print(f"L_W = {np.round(LW, 3).tolist()}, L_d = {L_D:.3f}")
for p in PRIMES:
    print(f"κ_{p}: " + "  ".join(f"{k:.3f}±{e:.3f}" for k, e in zip(K[p], E[p])) +
          f"   | κ^197 = {K197[p]:.3f}; H_S son {HS_son[p]:.3f}; B_p = {Bp[p]:.4f}")
print(f"γ̂ = {g0:.3f} ± {sg:.3f} (σ_eff); |γ̂−1|/σ = {abs(g0 - 1) / sg:.2f}; "
      f"|γ̂|/σ = {abs(g0) / sg:.2f}")
print("hüküm:", H["hukum"])

# ---- figür ----
fig, axs = plt.subplots(2, 2, figsize=(6.5, 4.9), sharex=True, layout="constrained")
xx = np.linspace(9.0, 12.35, 200)
YAR = 0.115
for n, (p, ax) in enumerate(zip(PRIMES, axs.flat)):
    ax.plot(xx, np.full_like(xx, K197[p]), color=YESIL, lw=1.5, ls=(0, (1, 1.4)), zorder=1)
    ax.plot(xx, 1 - (1 - K197[p]) * (xx / L_D) ** (-1.0), color=KIZIL, lw=1.4,
            ls=(0, (4, 1.8)), zorder=2)
    ax.plot(xx, 1 - Bp[p] * (xx / L_D) ** (-g0), color=MAVI, lw=1.6, zorder=3)
    ax.errorbar(LW, K[p], yerr=E[p], fmt="o", ms=5.2, color=INK, mec="white", mew=0.6,
                elinewidth=1.2, capsize=2.5, zorder=4)
    mid = 0.5 * (min(K[p].min(), HS_alt[p]) + max(K[p].max(), HS_son[p]))
    ax.set_ylim(mid - YAR, mid + YAR)
    ax.set_xlim(9.0, 12.35)
    ax.set_xticks(np.round(LW, 2))
    ax.text(0.03, 0.95, f"({'abcd'[n]})  $p={p}$", transform=ax.transAxes, va="top",
            ha="left", fontsize=10, fontweight="bold",
            bbox=dict(fc="white", ec="none", pad=1.2))
    ax.grid(alpha=0.25, lw=0.5)
fig.supxlabel(r"window height $L_W=\langle\log(t/2\pi)\rangle$", fontsize=10)
fig.supylabel(r"$\eta_p=A(+\log p)\,/\,(s\,c_{\rm BK}(p))$", fontsize=10)
tut = [plt.Line2D([], [], marker="o", ls="none", color=INK, mec="white", ms=5.2,
                  label=r"measured $\pm$ jk s.e."),
       plt.Line2D([], [], color=MAVI, lw=1.6,
                  label=rf"fit, $\hat\gamma={g0:.2f}\pm{sg:.2f}$"),
       plt.Line2D([], [], color=KIZIL, lw=1.4, ls=(0, (4, 1.8)), label=r"$H_S$ ($\gamma=1$)"),
       plt.Line2D([], [], color=YESIL, lw=1.5, ls=(0, (1, 1.4)), label=r"$H_C$ ($\gamma=0$)")]
fig.legend(handles=tut, loc="outside upper center", ncol=4, frameon=False,
           handlelength=2.2, columnspacing=1.3)

fig.savefig(CIKTI.with_suffix(".png"), dpi=300)
fig.savefig(CIKTI.with_suffix(".pdf"))
print(f"-> {CIKTI}.png / .pdf")
