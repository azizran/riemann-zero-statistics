# -*- coding: utf-8 -*-
"""
199fig_2 — NOT 7 FİGÜRÜ 2: fig_n7_classes_en.{png,pdf}
========================================================
Mühürlü 193 sonucunu (H-193a MÜHÜR, H-193b KAYIT, kontroller) gösterir; YENİ ÖLÇÜM YOK.
Sınıf payı s_r = κ_r/κ_top (bölünen sınıf payda dışı), ± 8-blok jk se; öngörü
s_r = cos(2π r b/a)/μ(a) (193a öngörü tablosu; K1_sinif.json 's_pred').
(a) +log10 (a = 10, b = 1; μ(10) = +1): r = 1, 3, 7, 9 — kör birincil hedef.
(b) +log7 (r = 1…6; kör ikincil; r2/r5 hükme girmeyen yan-KAYIT) ve kontroller
    +log3 (r = 1, 2; öngörü ½/½) ile −log3 (sıfır-yapı, keyfi mod-3 bölmesi; ½/½).
Öngörüler burada cos(2πrb/a)/μ(a)'dan YENİDEN hesaplanır ve s_pred ile assert edilir.
Girdi: scratchpad/193/{K1_sinif.json, HUKUM_193.json}
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
S193 = QM / "scratchpad" / "193"
CIKTI = QM / "fig_n7_classes_en"

MAVI, KIZIL, YESIL, MOR = "#0072B2", "#D55E00", "#009E73", "#AA3377"
INK, GRI = "#1a1a1a", "#8c8c8c"
plt.rcParams.update({"font.size": 9, "axes.labelsize": 10, "xtick.labelsize": 9,
                     "ytick.labelsize": 9, "legend.fontsize": 9, "axes.linewidth": 0.8,
                     "mathtext.fontset": "dejavusans", "savefig.dpi": 300,
                     "pdf.fonttype": 42, "ps.fonttype": 42})


def harf(ax, s, x=0.02, y=0.975):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=10, fontweight="bold",
            va="top", ha="left", bbox=dict(fc="white", ec="none", pad=1.2), zorder=6)


K1 = json.load(open(S193 / "K1_sinif.json"))["hedefler"]
H = json.load(open(S193 / "HUKUM_193.json"))


def mobius(n):
    k, p = 0, 2
    while p * p <= n:
        if n % p == 0:
            n //= p
            if n % p == 0:
                return 0
            k += 1
        p += 1
    if n > 1:
        k += 1
    return (-1) ** k


def seri(ad, a, b, siniflar):
    h = K1[ad]["siniflar"]
    s = np.array([h[f"r{r}"]["s"] for r in siniflar])
    e = np.array([h[f"r{r}"]["se_s"] for r in siniflar])
    p = np.array([h[f"r{r}"]["s_pred"] for r in siniflar])
    ong = np.array([np.cos(2 * np.pi * r * b / a) / mobius(a) for r in siniflar])
    if ad == "-log3":                        # sıfır-yapı kontrolü: a_teori = 1 ⇒ ½/½
        ong = np.full(len(siniflar), 0.5)
    assert np.allclose(ong, p, atol=5e-6), (ad, ong, p)
    return s, e, ong


s10, e10, p10 = seri("+log10", 10, 1, (1, 3, 7, 9))
s7, e7, p7 = seri("+log7", 7, 1, (1, 2, 3, 4, 5, 6))
s3, e3, p3 = seri("+log3", 3, 1, (1, 2))
m3, f3, q3 = seri("-log3", 3, 1, (1, 2))
for ad, s, e, p in (("+log10", s10, e10, p10), ("+log7", s7, e7, p7), ("+log3", s3, e3, p3),
                    ("-log3", m3, f3, q3)):
    print(f"{ad:>7}: s = {np.round(s, 4).tolist()} ± {np.round(e, 4).tolist()}  "
          f"öngörü {np.round(p, 4).tolist()}  oran {np.round(s / p, 3).tolist()}")
print("hüküm:", H.get("hukum"))

fig, (ax, bx) = plt.subplots(1, 2, figsize=(6.5, 2.95), layout="constrained",
                             gridspec_kw={"width_ratios": [4.6, 11.0]})
BAR = 0.30


def ciz(axx, x, s, e, p, etiket=True):
    for xi, pi in zip(x, p):
        axx.plot([xi - BAR, xi + BAR], [pi, pi], color=INK, lw=1.8, solid_capstyle="butt",
                 zorder=3)
    axx.errorbar(x, s, yerr=e, fmt="o", ms=5.2, color=MAVI, mfc=MAVI, mec="white", mew=0.6,
                 ecolor=MAVI, elinewidth=1.2, capsize=2.5, zorder=4)


# (a) +log10
x = np.arange(4)
ax.axhline(0, color=GRI, lw=0.6)
ciz(ax, x, s10, e10, p10)
ax.set_xticks(x)
ax.set_xticklabels(["1", "3", "7", "9"])
ax.set_xlim(-0.6, 3.6)
ax.set_ylim(-0.62, 1.12)
ax.set_xlabel(r"class $r$ ($q'\equiv r$ mod 10)")
ax.set_ylabel(r"class share $s_r=\kappa_r/\kappa_{\rm tot}$")
ax.text(0.5, 0.60, r"$\Delta\omega={+}\log 10$", transform=ax.transAxes, ha="center",
        va="center", fontsize=9)
harf(ax, "(a)")

# (b) +log7 | +log3 | −log3
x7 = np.arange(6)
x3 = np.array([7.2, 8.2])
xm = np.array([9.9, 10.9])
bx.axhline(0, color=GRI, lw=0.6)
for xs in (6.1, 8.8):
    bx.axvline(xs, color=GRI, lw=0.6, ls=(0, (2, 2)))
ciz(bx, x7, s7, e7, p7)
ciz(bx, x3, s3, e3, p3)
ciz(bx, xm, m3, f3, q3)
bx.set_xticks(np.concatenate([x7, x3, xm]))
bx.set_xticklabels(["1", "2", "3", "4", "5", "6", "1", "2", "1", "2"])
bx.set_xlim(-0.6, 11.5)
bx.set_ylim(-2.05, 2.35)
bx.set_xlabel(r"class $r$ ($q'\equiv r$ mod $a$)")
for xc, t in ((2.5, r"$+\log 7$"), (7.7, r"$+\log 3$"), (10.4, r"$-\log 3$")):
    bx.text(xc, 2.22, t, ha="center", va="top", fontsize=9)
harf(bx, "(b)", x=0.012)

h1 = plt.Line2D([], [], marker="o", ls="none", color=MAVI, mec="white", ms=5.2,
                label=r"measured $s_r\pm$ jk s.e.")
h2 = plt.Line2D([], [], color=INK, lw=1.8, label=r"$\cos(2\pi rb/a)/\mu(a)$")
bx.legend(handles=[h1, h2], loc="lower right", frameon=True, framealpha=0.95,
          edgecolor="none", borderaxespad=0.3, handlelength=1.6)

fig.savefig(CIKTI.with_suffix(".png"), dpi=300)
fig.savefig(CIKTI.with_suffix(".pdf"))
print(f"-> {CIKTI}.png / .pdf")
