# -*- coding: utf-8 -*-
"""
190d — FİGÜR: 190_tarak_evrensellik.png
  sol : blok-birleşik HAVUZ κ yoğunluğu Δω = ω' − L_b ekseninde — düşük pencere
        (ω-dilim 0.025, 190b) vs SON (188 τ'-dilim blok profilleri, ortak ızgaraya
        ara değerlenmiş); her biri kendi maksimumuna normalize. Hedefler
        (0, ±log2, ±log3, +log6) DÜZ, τ'-sabit rakipler KESİKLİ; düşük pencerenin
        blok-medyan tepeleri üçgenle.
  sağ : düşük pencere Re v(τ') ile −S_Re ve W (her biri max|·|'a normalize);
        corr değerleri başlıkta.
Girdi: 190/profiller_190.npz, 190/HUKUM_190.json.
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
S190 = SCR / "190"

C1, C2, C3 = "#2a78d6", "#eb6834", "#1baf7a"   # 188e paleti (mavi/turuncu/su yeşili)
GRI, MUREKKEP, IKINCIL = "#8a8985", "#0b0b0b", "#52514e"

P = np.load(S190 / "profiller_190.npz")
H = json.load(open(S190 / "HUKUM_190.json"))
ONK = json.load(open(S190 / "ONKAYIT_190.json"))
L_SON = ONK["L_son"]

# --- düşük: ω-dilim, bloklar nb-ağırlıklı birleşik yoğunluk ---
m = P["om_merkez"]
nb = P["om_nb"].astype(float)
dw = float(ONK["omega_dilim"]["genislik"])
yog_d = (P["om_kap"] * (nb / nb.sum())[:, None]).sum(0) / dw

# --- son: τ'-dilim blok profilleri → ortak ızgara ---
izg = np.arange(-1.6, 3.0, 0.005)
nbs = P["nb_s"].astype(float)
w_s = 0.005 * L_SON
yog_s = np.zeros_like(izg)
for b in range(8):
    x, y = P["ts_dw"][b], P["ts_kap"][b] / w_s
    yog_s += (nbs[b] / nbs.sum()) * np.interp(izg, x, y, left=0.0, right=0.0)

ARALIK = (-1.3, 2.3)


def norm(x, y):
    msk = (x >= ARALIK[0]) & (x <= ARALIK[1])
    return y / np.max(y[msk])


plt.rcParams.update({"font.size": 9, "axes.edgecolor": GRI,
                     "axes.labelcolor": MUREKKEP, "xtick.color": IKINCIL,
                     "ytick.color": IKINCIL, "axes.spines.top": False,
                     "axes.spines.right": False})
fig, (ax, ax2) = plt.subplots(1, 2, figsize=(13.0, 4.9), facecolor="#fcfcfb",
                              gridspec_kw={"width_ratios": [1.45, 1.0],
                                           "wspace": 0.18})

# ---------------- sol ----------------
ax.axhline(0, color=GRI, lw=0.6)
ax.plot(izg, norm(izg, yog_s), color=C2, lw=1.4, ls="--",
        label=f"son (L={L_SON:.2f}; 188 τ'-dilim, 0.060 ω)")
ax.step(m, norm(m, yog_d), where="mid", color=C1, lw=1.5,
        label=f"düşük (L={ONK['pencere']['L']:.2f}; ω-dilim 0.025)")
hed = {"0": 0.0, "+log2": np.log(2), "+log3": np.log(3), "+log6": np.log(6),
       "−log2": -np.log(2), "−log3": -np.log(3)}
for ad, h in hed.items():
    ax.axvline(h, color=MUREKKEP, lw=0.8, alpha=0.55)
    ax.text(h + 0.02, 1.13, ad, ha="left", va="bottom", fontsize=7.5, color=IKINCIL)
for ad, r in ONK["rakip_tau_sabit"].items():
    if r is not None:
        ax.axvline(r, color=C2, lw=0.9, ls=(0, (3, 2)), alpha=0.8)
tab = H["H190a_tablo"]
for ad, t in tab.items():
    ax.plot([t["medyan"]], [1.04], marker="v", ms=7, color=C1,
            markeredgecolor="white", markeredgewidth=0.8, zorder=5)
ax.set_xlim(*ARALIK)
ax.set_ylim(-0.45, 1.2)
ax.set_xlabel("Δω = ω' − L_b   (blok-yerel)")
ax.set_ylabel("HAVUZ κ yoğunluğu (normalize)")
hA = H["hukum"]["H-190a"]
ax.set_title(f"Δω profili — düz: hedef (0, ±log2, ±log3, +log6); kesikli turuncu: "
             f"τ'-sabit rakip\n▼ düşük pencere blok-medyan tepeleri   ·   H-190a: {hA}",
             fontsize=9, color=MUREKKEP, loc="left")
ax.legend(loc="lower left", frameon=False, fontsize=8)

# ---------------- sağ ----------------
o = P["orta"]
v = P["v_d"]
c = H["H190b_c"]["dusuk"]
cs = H["H190b_c"]["son_188_yeniden"]


def n1(y):
    return y / np.max(np.abs(y))


ax2.axhline(0, color=GRI, lw=0.6)
ax2.plot(o, n1(P["W_d"]), color=GRI, lw=1.0, label="W = Σa' (yapısız)")
ax2.plot(o, n1(-P["SRe_d"]), color=C2, lw=1.5, ls="--", label="−S_Re")
ax2.plot(o, n1(v.real), color=C1, lw=1.6, label="Re v(τ')")
ax2.axvline(1.0, color=MUREKKEP, lw=0.7, ls=":")
ax2.set_xlabel("τ' = ω'/L_düşük")
ax2.set_ylabel("normalize (max|·| = 1)")
ax2.set_title(f"düşük pencere: corr(Re v, S_Re) = {c['corr_RevSRe']:+.3f}±"
              f"{c['se_RevSRe']:.3f};  corr(Re v, W) = {c['corr_RevW']:+.3f}±"
              f"{c['se_RevW']:.3f}\n(son, 188: {cs['corr_RevSRe']:+.3f} / "
              f"{cs['corr_RevW']:+.3f})   H-190b: {H['hukum']['H-190b']}  ·  "
              f"H-190c: {H['hukum']['H-190c']}", fontsize=9, color=MUREKKEP,
              loc="left")
ax2.legend(loc="lower left", frameon=False, fontsize=8)

fig.suptitle("190 — Tarağın evrenselliği: düşük pencerede örneklem-dışı sınav",
             x=0.01, y=1.04, ha="left", fontsize=11, color=MUREKKEP)
yol = QM / "190_tarak_evrensellik.png"
fig.savefig(yol, dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
print(f"-> {yol}")
