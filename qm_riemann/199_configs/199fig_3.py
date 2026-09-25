# -*- coding: utf-8 -*-
"""
199fig_3 — NOT 7 FİGÜRÜ 3: fig_n7_characters_en.{png,pdf}
===========================================================
Mühürlü 195 sonucunu (B — çekirdek K̃, üst bölge; H-195B1 KAYIT, H-195B2 MÜHÜR) gösterir;
YENİ ÖLÇÜM YOK.
(a) Re K̃ ± jk se (symlog y) — beş L-fonksiyonu adası (β = χ₄, χ₃, χ₅ₑ, χ₈ₑ, χ₈ₒ; üst bölge
    L_χ ≥ 8.5) ve ζ referansı (K0c, düşük pencere, adalarla AYNI ΔL ızgarası). Renk: izinli
    (n'nin hiçbir asal çarpanı k'yı bölmez) mavi dolu / yasak (p | k) kızıl boş. −log5 B'den
    ön-kayıtla dışlandı (pencere-içi), çizilmez.
(b) ∠K̃: izinliler 180°, yasaklar 0° çevresinde (kırık eksen).
Rapordaki özet iddialar (∠180° ± 1.3°, yasaklar ∠0° ± 7.2°, +0.0014…+0.0063, z −55…−82,
ρ₈ = 1.009) burada veriden yeniden hesaplanıp ekrana basılır.
Girdi: scratchpad/195/{B_cekirdek.json, K0c.json, HUKUM_195.json}
"""
import json
from fractions import Fraction as F
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
S195 = QM / "scratchpad" / "195"
CIKTI = QM / "fig_n7_characters_en"

MAVI, KIZIL, YESIL, MOR = "#0072B2", "#D55E00", "#009E73", "#AA3377"
INK, GRI = "#1a1a1a", "#8c8c8c"
plt.rcParams.update({"font.size": 9, "axes.labelsize": 10, "xtick.labelsize": 9,
                     "ytick.labelsize": 9, "legend.fontsize": 9, "axes.linewidth": 0.8,
                     "mathtext.fontset": "dejavusans", "savefig.dpi": 300,
                     "pdf.fonttype": 42, "ps.fonttype": 42})


def harf(ax, s, x=0.012, y=0.975):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=10, fontweight="bold",
            va="top", ha="left", bbox=dict(fc="white", ec="none", pad=1.2), zorder=6)


B = json.load(open(S195 / "B_cekirdek.json"))["adalar"]
Z = json.load(open(S195 / "K0c.json"))["Y_yeni_Ktilde"]["dusuk_izgara_tum"]
HK = json.load(open(S195 / "HUKUM_195.json"))

ADALAR = [("beta", r"$\chi_4$", "o"), ("chi3", r"$\chi_3$", "s"), ("chi5e", r"$\chi_{5e}$", "^"),
          ("chi8e", r"$\chi_{8e}$", "D"), ("chi8o", r"$\chi_{8o}$", "v")]
UYDU = [("+log2", F(2), r"${+}\log 2$"), ("-log2", F(1, 2), r"$-\log 2$"),
        ("+log3", F(3), r"${+}\log 3$"), ("-log3", F(1, 3), r"$-\log 3$"),
        ("+log4", F(4), r"${+}\log 4$"), ("+log5", F(5), r"${+}\log 5$"),
        ("+log6", F(6), r"${+}\log 6$"), ("+log10", F(10), r"${+}\log 10$"),
        ("+log(3/2)", F(3, 2), r"${+}\log\frac{3}{2}$")]


def asallar(n):
    out, p = set(), 2
    while p * p <= n:
        while n % p == 0:
            out.add(p)
            n //= p
        p += 1
    if n > 1:
        out.add(n)
    return out


def yasak(r, k):
    return any(k % p == 0 for p in asallar(r.numerator) | asallar(r.denominator))


def mu0(r):   # μ(a) = 0 (a = pay) — +log4
    a = r.numerator
    return any(a % (p * p) == 0 for p in asallar(a))


# ---- denetim: rapordaki özet iddialar ----
izin_aci, yasak_aci, yasak_re, izin_imz = [], [], [], []
for ad, _, _ in ADALAR:
    k = B[ad]["k"]
    for u, r, _ in UYDU:
        d = B[ad]["ust"][u]
        assert d["n_pencere_ici"] == 0, (ad, u)
        if yasak(r, k):
            yasak_aci.append(d["aci"])
            yasak_re.append(d["re"])
        elif not mu0(r):
            izin_aci.append(abs(abs(d["aci"]) - 180.0))
            izin_imz.append(abs(d["z_im"]))
print(f"izinli (μ≠0) ∠: maks |∠ − 180°| = {max(izin_aci):.2f}°  (rapor ±1.3°); "
      f"maks |Im z| = {max(izin_imz):.2f} (rapor ≤ 2.1)")
print(f"yasak ∠: maks |∠| = {max(abs(a) for a in yasak_aci):.2f}° (rapor ±7.2°); "
      f"Re K̃ ∈ [{min(yasak_re):+.5f}, {max(yasak_re):+.5f}] (rapor +0.0014…+0.0063); "
      f"n = {len(yasak_re)}")
bir = [("chi3", "+log2"), ("chi5e", "+log2"), ("chi5e", "+log3"), ("beta", "+log3"),
       ("chi8e", "+log3"), ("chi8o", "+log3")]
print("birincil izinli z_Re:", {f"{a}:{u}": round(B[a]["ust"][u]["z_re"], 1) for a, u in bir})
r8 = B["chi8o"]["ust"]["+log3"]["re"] / B["chi8e"]["ust"]["+log3"]["re"]
assert abs(r8 - HK["H_195B2"]["rho8"]) < 1e-12
print(f"ρ₈(+log3) = {r8:.3f} ± {HK['H_195B2']['se_rho8']:.3f} (rapor 1.009 ± 0.024)")
print("hüküm:", {h: HK[h]["hukum"] for h in ("H_195A", "H_195B1", "H_195B2", "H_195B3")})

# ---- figür ----
fig, (ax, bx) = plt.subplots(2, 1, figsize=(6.5, 5.6), layout="constrained",
                             gridspec_kw={"height_ratios": [1.75, 1.0]})
OFF = np.linspace(-0.34, 0.34, 6)
ax.axhline(0, color=GRI, lw=0.6, zorder=1)
for i in range(len(UYDU) - 1):
    ax.axvline(i + 0.5, color="0.9", lw=0.6, zorder=0)
for j, (ad, et, mk) in enumerate(ADALAR):
    k = B[ad]["k"]
    for i, (u, r, _) in enumerate(UYDU):
        d = B[ad]["ust"][u]
        y = yasak(r, k)
        c = KIZIL if y else MAVI
        ax.errorbar(i + OFF[j], d["re"], yerr=d["se_re"], fmt=mk, ms=4.6, color=c,
                    mfc="white" if y else c, mec=c, mew=1.0, elinewidth=0.9, capsize=0,
                    zorder=3)
for i, (u, r, _) in enumerate(UYDU):
    d = Z[u]
    ax.errorbar(i + OFF[5], d["re"], yerr=d["se_re"], fmt="*", ms=7.0, color=INK,
                mfc=INK, mec=INK, elinewidth=0.9, capsize=0, zorder=3)
ax.set_yscale("symlog", linthresh=1e-3, linscale=0.45)
ax.set_ylim(-0.32, 0.028)
ax.set_yticks([-0.1, -0.01, -0.001, 0, 0.001, 0.01])
ax.set_xticks(range(len(UYDU)))
ax.set_xticklabels([s for _, _, s in UYDU])
ax.set_xlim(-0.55, len(UYDU) - 0.45)
ax.set_ylabel(r"Re $\tilde K$  (symlog)")
ax.set_xlabel(r"satellite $\Delta\omega=\log(a/b)$")
harf(ax, "(a)")
tut = [plt.Line2D([], [], marker=mk, ls="none", color=GRI, mfc=GRI, mec=GRI, ms=5,
                  label=et) for _, et, mk in ADALAR]
tut.append(plt.Line2D([], [], marker="*", ls="none", color=INK, ms=7, label=r"$\zeta$"))
tut += [plt.Line2D([], [], marker="o", ls="none", color=MAVI, mfc=MAVI, ms=5,
                   label=r"allowed ($p\nmid k$)"),
        plt.Line2D([], [], marker="o", ls="none", color=KIZIL, mfc="white", mew=1.0, ms=5,
                   label=r"forbidden ($p\mid k$)")]
ax.legend(handles=tut, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=8, frameon=False,
          borderaxespad=0.2, handletextpad=0.15, columnspacing=0.9)

# (b) faz: kırık eksen içinde iki pencere — sol ∠ ≈ 0°, sağ ∠ ≈ 180°
bx.set_axis_off()
gs = bx.get_subplotspec().subgridspec(1, 2, wspace=0.07)
b0 = fig.add_subplot(gs[0])
b1 = fig.add_subplot(gs[1], sharey=b0)
SATIR = [et for _, et, _ in ADALAR] + [r"$\zeta$"]
dy = np.linspace(-0.28, 0.28, len(UYDU))


def aci180(a):
    return a % 360.0          # −179.9 → 180.1


for j, (ad, et, mk) in enumerate(ADALAR):
    k = B[ad]["k"]
    for i, (u, r, _) in enumerate(UYDU):
        d = B[ad]["ust"][u]
        y = yasak(r, k)
        c = KIZIL if y else MAVI
        a = d["aci"]
        tgt = b0 if abs(a) < 90 else b1
        xx = a if abs(a) < 90 else aci180(a)
        tgt.plot(xx, j + dy[i], mk, ms=4.2, color=c, mfc="white" if y else c, mec=c, mew=1.0)
for i, (u, r, _) in enumerate(UYDU):
    a = Z[u]["aci"]
    tgt = b0 if abs(a) < 90 else b1
    tgt.plot(a if abs(a) < 90 else aci180(a), 5 + dy[i], "*", ms=6.0, color=INK)
b0.set_xlim(-12, 12)
b1.set_xlim(168, 192)
b0.axvline(0, color=GRI, lw=0.7, zorder=0)
b1.axvline(180, color=GRI, lw=0.7, zorder=0)
b0.set_yticks(range(6))
b0.set_yticklabels(SATIR)
b0.set_ylim(5.55, -0.55)
b0.set_xticks([-10, -5, 0, 5])
b1.set_xticks([175, 180, 185, 190])
plt.setp(b1.get_yticklabels(), visible=False)
b1.tick_params(axis="y", length=0)
b0.spines["right"].set_visible(False)
b1.spines["left"].set_visible(False)
for axx, xk in ((b0, 1.0), (b1, 0.0)):
    kw = dict(transform=axx.transAxes, color=INK, clip_on=False, lw=0.8)
    axx.plot([xk - 0.012, xk + 0.012], [-0.03, 0.03], **kw)
    axx.plot([xk - 0.012, xk + 0.012], [0.97, 1.03], **kw)
fig.supxlabel(r"phase $\angle\tilde K$ (degrees)", fontsize=10, x=0.54)
b0.text(0.02, 0.97, "(b)", transform=b0.transAxes, fontsize=10, fontweight="bold",
        va="top", ha="left", bbox=dict(fc="white", ec="none", pad=1.2))
b0.text(0.97, 0.03, "forbidden", transform=b0.transAxes, ha="right", va="bottom", fontsize=9,
        color=KIZIL)
b1.text(0.03, 0.03, "allowed", transform=b1.transAxes, ha="left", va="bottom", fontsize=9,
        color=MAVI)

fig.savefig(CIKTI.with_suffix(".png"), dpi=300)
fig.savefig(CIKTI.with_suffix(".pdf"))
print(f"-> {CIKTI}.png / .pdf")
