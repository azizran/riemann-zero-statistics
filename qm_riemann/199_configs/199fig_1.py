# -*- coding: utf-8 -*-
"""
199fig_1 — NOT 7 FİGÜRÜ 1: fig_n7_positions_en.{png,pdf}
==========================================================
Mühürlü 190 sonucunu (H-190a MÜHÜR + 190e ön-kayıtsız keşif (a)) gösterir; YENİ ÖLÇÜM YOK.
Profil kurgusu 190d_figur.py / 190e_kesif.py ile BİREBİR aynıdır:
  düşük : ω-dilim (0.025) blok profilleri, n_b ağırlıklı birleşik yoğunluk / 0.025
  son   : 188 τ'-dilim blok profilleri, ortak ızgaraya np.interp (left=right=0)
  normalizasyon: her profil kendi maks'ına, Δω ∈ [−1.3, 2.3] içinde (190d AYNEN)
(a) Δω = ω' − L_loc ekseninde iki profil; dikey çizgiler uydu konumları ±log(a/b):
    düz = 190'ın ön-kayıtlı hedefleri (0, ±log2, ±log3, +log6); noktalı = 192 kataloğunda
    iki pencerede de yanan diğer uydular (±log 3/2, +log 5/4, +log 5, +log 10).
(b) aynı profiller τ'-ölçekli eksende (Δω / L_W): 190e(a)'nın τ'-sabit okuması.
Korelasyonlar 190e ile aynı yoldan yeniden hesaplanır ve K_kesif_190.json ile assert edilir.
Girdi: scratchpad/190/{profiller_190.npz, ONKAYIT_190.json, HUKUM_190.json, K_kesif_190.json}
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
S190 = QM / "scratchpad" / "190"
CIKTI = QM / "fig_n7_positions_en"

# ---- ortak stil (Not 7) ----
MAVI, KIZIL, YESIL, MOR = "#0072B2", "#D55E00", "#009E73", "#AA3377"
INK, GRI = "#1a1a1a", "#8c8c8c"
plt.rcParams.update({"font.size": 9, "axes.labelsize": 10, "xtick.labelsize": 9,
                     "ytick.labelsize": 9, "legend.fontsize": 9, "axes.linewidth": 0.8,
                     "mathtext.fontset": "dejavusans", "savefig.dpi": 300,
                     "pdf.fonttype": 42, "ps.fonttype": 42})


def harf(ax, s, x=0.012, y=0.965):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=10, fontweight="bold",
            va="top", ha="left", bbox=dict(fc="white", ec="none", pad=1.2), zorder=6)


# ---- veri (190d / 190e AYNEN) ----
P = np.load(S190 / "profiller_190.npz")
ONK = json.load(open(S190 / "ONKAYIT_190.json"))
H = json.load(open(S190 / "HUKUM_190.json"))
KES = json.load(open(S190 / "K_kesif_190.json"))
L_D = float(ONK["pencere"]["L"])
L_S = float(ONK["L_son"])

m = P["om_merkez"]
nb = P["om_nb"].astype(float)
yog_d = (P["om_kap"] * (nb / nb.sum())[:, None]).sum(0) / 0.025
nbs = P["nb_s"].astype(float)


def son_yog(x):
    y = np.zeros_like(x)
    for b in range(8):
        y += (nbs[b] / nbs.sum()) * np.interp(x, P["ts_dw"][b], P["ts_kap"][b] / (0.005 * L_S),
                                              left=0.0, right=0.0)
    return y


ARALIK = (-1.3, 2.3)
msk = (m >= ARALIK[0]) & (m <= ARALIK[1])
c_dw = float(np.corrcoef(yog_d[msk], son_yog(m[msk]))[0, 1])
c_tau = float(np.corrcoef(yog_d[msk], son_yog(m[msk] * L_S / L_D))[0, 1])
assert abs(c_dw - KES["a_profil_corr_Δω"]) < 1e-12, c_dw
assert abs(c_tau - KES["a_profil_corr_tau_olcekli"]) < 1e-12, c_tau

izg = np.arange(-1.6, 3.0, 0.005)          # 190d ızgarası AYNEN
yog_s = son_yog(izg)
nd = yog_d / yog_d[msk].max()
izm = (izg >= ARALIK[0]) & (izg <= ARALIK[1])
ns = yog_s / yog_s[izm].max()

# düşük pencerede tüm 8 bloğun ω-kapsaması: Δω ≥ 0.86·L_D − min L_b (190a tanımı)
Lb = np.array(ONK["bloklar"]["L_b"])
KAPSAMA_ALT = 0.86 * L_D - Lb.min()

# blok-medyan tepeler (H-190a tablosu; yalnız konsol denetimi)
tepe = {k: v["medyan"] for k, v in H["H190a_tablo"].items()}
print(f"L_düşük = {L_D:.4f}, L_son = {L_S:.4f}")
print(f"corr Δω = {c_dw:.3f} | corr τ'-ölçekli = {c_tau:.3f}  (190e: 0.980 / 0.002)")
print(f"blok-medyan tepeler (190c): {tepe}")
print(f"düşük pencere tam blok kapsaması Δω ≥ {KAPSAMA_ALT:.3f}")

# ---- uydu konumları ----
ONKAYITLI = [(0.0, r"$0$"), (np.log(2), r"$\log 2$"), (-np.log(2), r"$-\log 2$"),
             (np.log(3), r"$\log 3$"), (-np.log(3), r"$-\log 3$"), (np.log(6), r"$\log 6$")]
DIGER = [(np.log(1.5), r"$\log\,3/2$"), (-np.log(1.5), r"$-\log\,3/2$"),
         (np.log(1.25), r"$\log\,5/4$"), (np.log(5), r"$\log 5$"), (np.log(10), r"$\log 10$")]

# ---- figür ----
fig, (ax, bx) = plt.subplots(2, 1, figsize=(6.5, 5.3), layout="constrained",
                             gridspec_kw={"height_ratios": [1.55, 1.0]})
XL = (-1.3, 2.36)

ax.axvspan(XL[0], KAPSAMA_ALT, color="0.93", lw=0, zorder=0)
for x, _ in ONKAYITLI:
    ax.axvline(x, color=INK, lw=0.7, alpha=0.55, zorder=1)
for x, _ in DIGER:
    ax.axvline(x, color=INK, lw=0.8, ls=(0, (1, 2)), alpha=0.6, zorder=1)
ax.axhline(0, color=GRI, lw=0.6, zorder=1)
ax.plot(izg, ns, color=KIZIL, lw=1.2, ls=(0, (4, 1.6)), zorder=3,
        label=rf"window $L={L_S:.2f}$")
ax.step(m, nd, where="mid", color=MAVI, lw=1.1, zorder=4,
        label=rf"window $L={L_D:.2f}$")
ax.set_xlim(*XL)
ax.set_ylim(-0.16, 1.12)
ax.set_xlabel(r"$\Delta\omega=\omega'-L_{\rm loc}$")
ax.set_ylabel(r"$\kappa$ density (norm.)")
hh, ll = ax.get_legend_handles_labels()
ax.legend(hh[::-1], ll[::-1], loc="upper right", frameon=True, framealpha=0.95,
          edgecolor="none", handlelength=2.6, borderaxespad=0.3)
ax.text(0.075, 0.955, f"corr = {c_dw:.3f}", transform=ax.transAxes, va="top", ha="left",
        fontsize=9, bbox=dict(fc="white", ec="none", pad=1.5), zorder=6)
top = ax.twiny()
top.set_xlim(*XL)
pos = ONKAYITLI + DIGER
top.set_xticks([x for x, _ in pos])
top.set_xticklabels([s for _, s in pos], rotation=90, fontsize=9)
top.tick_params(axis="x", length=2.5, pad=1.5)
harf(ax, "(a)")

# (b) τ'-ölçekli eksen: her profil kendi penceresinin L'siyle
bx.axhline(0, color=GRI, lw=0.6)
bx.plot(izg / L_S, ns, color=KIZIL, lw=1.2, ls=(0, (4, 1.6)), zorder=3,
        label=rf"window $L={L_S:.2f}$")
bx.step(m / L_D, nd, where="mid", color=MAVI, lw=1.1, zorder=4,
        label=rf"window $L={L_D:.2f}$")
for n in (2, 3, 6):
    for L, c in ((L_D, MAVI), (L_S, KIZIL)):
        bx.plot([np.log(n) / L] * 2, [1.09, 1.18], color=c, lw=1.4, solid_capstyle="butt")
bx.set_xlim(XL[0] / L_D, XL[1] / L_D)
bx.set_ylim(-0.16, 1.21)
bx.set_xlabel(r"$\Delta\omega/L$  ($\tau'$-scaled axis)")
bx.set_ylabel(r"$\kappa$ density (norm.)")
bx.text(0.075, 0.955, f"corr = {c_tau:.3f}", transform=bx.transAxes, va="top", ha="left",
        fontsize=9, bbox=dict(fc="white", ec="none", pad=1.5), zorder=6)
harf(bx, "(b)")

fig.savefig(CIKTI.with_suffix(".png"), dpi=300)
fig.savefig(CIKTI.with_suffix(".pdf"))
print(f"-> {CIKTI}.png / .pdf")
