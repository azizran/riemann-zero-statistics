# -*- coding: utf-8 -*-
"""
194d — FİGÜR: 194_bos_pencere.png
Sol panel: Δω ekseninde boş (gri) ve μ=0 (kırmızı) pencerelerin κ ± jk se'si,
κ̄_boş bandı (±se, yatay). Sağ panel: 193'ün +log10/+log7/+log5 payları —
ham (açık) vs κ̄_boş ile düzeltilmiş (koyu) vs cos yasası öngörüsü (çizgi).
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
S194 = SCR / "194"
ONK = json.load(open(S194 / "ONKAYIT_194.json"))
K1 = json.load(open(S194 / "K1_bos_194.json"))
H = json.load(open(S194 / "HUKUM_194.json"))

INK, INK2, MUTED = "#0b0b0b", "#52514e", "#a3a29c"
GRI, KIRMIZI, YESIL, TURUNCU = "#8a8a86", "#e34948", "#008300", "#e08214"
AMBER = "#b58900"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.titlesize": 9.5})

fig, (axL, axR) = plt.subplots(1, 2, figsize=(14.5, 5.4))

# ---------------- sol: boş + mu0 pencereler ----------------
BOS = K1["bos_pencereler"]
MU0 = K1["mu0_pencereler"]
xb = [w["merkez"] for w in BOS.values()]
yb = [w["kappa"] for w in BOS.values()]
eb = [w["se_kappa"] for w in BOS.values()]
xm = [w["merkez"] for w in MU0.values()]
ym = [w["kappa"] for w in MU0.values()]
em = [w["se_kappa"] for w in MU0.values()]

axL.errorbar(xb, yb, yerr=eb, fmt="o", color=GRI, ecolor=GRI, capsize=3, ms=6,
            label=f"boş pencere (N={len(BOS)})", zorder=4)
axL.errorbar(xm, ym, yerr=em, fmt="s", color=KIRMIZI, ecolor=KIRMIZI, capsize=3, ms=6,
            label=f"μ=0 pencere (N={len(MU0)})", zorder=4)

kbar = H["kappa_bar_bos"]
sebar = H["se_bos"]
axL.axhspan(kbar - sebar, kbar + sebar, color=AMBER, alpha=0.18, zorder=1,
           label=f"κ̄_boş = {kbar:+.4f}±{sebar:.4f}")
axL.axhline(kbar, color=AMBER, lw=1.3, zorder=2)
axL.axhline(0, color=MUTED, lw=0.8, zorder=1)

kbar_mu0 = H["kappa_bar_mu0"]
sebar_mu0 = H["se_mu0"]
axL.axhline(kbar_mu0, color=KIRMIZI, lw=1.0, ls=(0, (4, 2)), zorder=2,
           label=f"κ̄_μ0 = {kbar_mu0:+.4f}±{sebar_mu0:.4f}")

axL.set_xlabel("Δω (blok-yerel pencere merkezi)")
axL.set_ylabel("κ (HAVUZ, −Re K)")
axL.legend(fontsize=7.6, frameon=False, loc="best")
axL.set_title(f"H-194a: {H['hukum_ozet']['H-194a']}   |   H-194d: {H['hukum_ozet']['H-194d']}\n"
             f"H-194b (eğim, KAYIT): {H['H_194b']['egim']:+.5f}±{H['H_194b']['se_egim']:.5f}",
             fontsize=8.6)

# ---------------- sağ: 193 payları ham vs düzeltilmiş vs öngörü ----------------
E = H["H_194e"]["hedefler"]
SIRA = ["+log10", "+log7", "+log5"]
etiketler = []
x_ham, x_duz, x_pred = [], [], []
i = 0
xt = []
for ad in SIRA:
    rs = sorted(E[ad]["siniflar"].keys(), key=lambda k: int(k[1:]))
    for r in rs:
        rv = E[ad]["siniflar"][r]
        etiketler.append(f"{ad}\n{r}")
        x_ham.append(rv["s_ham"])
        x_duz.append(rv["s_duzeltilmis"])
        x_pred.append(rv["s_pred"])
        xt.append(i)
        i += 1
    i += 0.6  # uydular arası boşluk

xt = np.array(xt)
axR.scatter(xt - 0.15, x_ham, marker="o", s=34, color=GRI, label="ham (193)", zorder=3)
renk_duz = []
for ad in SIRA:
    for r in sorted(E[ad]["siniflar"].keys(), key=lambda k: int(k[1:])):
        rv = E[ad]["siniflar"][r]
        renk_duz.append(YESIL if rv["bantta_010"] else KIRMIZI)
axR.scatter(xt + 0.15, x_duz, marker="D", s=42, color=renk_duz, zorder=4,
           label="düzeltilmiş (κ̄_boş ile)")
axR.scatter(xt, x_pred, marker="_", s=380, color=INK, linewidths=2.0, zorder=5,
           label="cos yasası öngörüsü")
for xi, lo, hi in zip(xt, np.array(x_pred) - 0.10, np.array(x_pred) + 0.10):
    axR.plot([xi, xi], [lo, hi], color=INK, lw=0.6, alpha=0.35, zorder=1)

axR.axhline(0, color=MUTED, lw=0.8, zorder=1)
axR.set_xticks(xt)
axR.set_xticklabels(etiketler, fontsize=7.2)
axR.set_ylabel("s_r = κ_r / κ_top")
axR.legend(fontsize=7.6, frameon=False, loc="best")
axR.set_title(f"H-194e: {H['hukum_ozet']['H-194e']}", fontsize=8.6)

fig.suptitle("194 — Boş pencere tabanı: κ̄_boş (sol, δ=0.03 blok-yerel Δω, HAVUZ) ve "
            "193 sınıf paylarının düzeltilmesi (sağ)  [ön-kayıt sha "
            f"{ONK['sha256'][:10]}]", fontsize=11, color=INK, x=0.02, ha="left")
fig.subplots_adjust(left=0.06, right=0.985, top=0.82, bottom=0.14, wspace=0.22)
out = QM / "194_bos_pencere.png"
fig.savefig(out, dpi=150)
print(f"-> {out}")
