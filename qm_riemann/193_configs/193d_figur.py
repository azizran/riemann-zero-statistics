# -*- coding: utf-8 -*-
"""
193d — FİGÜR: 193_sinif_yasasi.png
Her uydu (+log10, +log7, +log5, +log3, −log3) için ölçülen sınıf payları
s_r = κ_r/κ_top (çubuk, ± jk se) ve KALEM pay-yasası öngörüsü (kesikli çizgi);
renk: H-193a/b'ye giren sınıflarda yön+nicel ikisi de OK ise yeşil, yalnız yön
OK ise turuncu, ikisi de değilse kırmızı; KAYIT-yalnız sınıflarda (+log5, r2/r5
yan-KAYIT, +log3/−log3 kontrol) mavi. Başlıkta κ_top ± se ve hüküm.
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
S193 = SCR / "193"
ONK = json.load(open(S193 / "ONKAYIT_193.json"))
K1 = json.load(open(S193 / "K1_sinif.json"))
H = json.load(open(S193 / "HUKUM_193.json"))

INK, INK2, MUTED = "#0b0b0b", "#52514e", "#a3a29c"
YESIL, KIRMIZI, TURUNCU, MAVI = "#008300", "#e34948", "#e08214", "#2a78d6"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.titlesize": 8.7})

SIRA = ["+log10", "+log7", "+log5", "+log3", "-log3"]
fig, axes = plt.subplots(1, 5, figsize=(18, 4.6))
for ax, ad in zip(axes, SIRA):
    h = K1["hedefler"][ad]
    S_ = h["siniflar"]
    rs = sorted(S_.keys(), key=lambda x: int(x[1:]))
    s_vals = np.array([S_[r]["s"] for r in rs])
    se_vals = np.array([S_[r]["se_s"] for r in rs])
    preds = np.array([S_[r]["s_pred"] for r in rs])
    x = np.arange(len(rs))

    hd = H["hedefler"].get(ad) or (H["KAYIT"].get(ad) if ad in H["KAYIT"] else None)
    hukum_sinif = (hd or {}).get("siniflar", {}) if isinstance(hd, dict) else {}
    renkler = []
    for r in rs:
        info = hukum_sinif.get(r)
        if info is not None and "yon_ok" in info:
            ok_y, ok_n = info["yon_ok"], info["nicel_ok"]
            renkler.append(YESIL if (ok_y and ok_n) else (TURUNCU if ok_y else KIRMIZI))
        elif info is not None and "bantta_0.5pm0.1" in info:
            renkler.append(YESIL if info["bantta_0.5pm0.1"] else KIRMIZI)
        else:
            renkler.append(MAVI)
    ax.bar(x, s_vals, color=renkler, width=0.62, zorder=3, alpha=0.88)
    ax.errorbar(x, s_vals, yerr=se_vals, fmt="none", ecolor=INK, elinewidth=1.1,
               capsize=3, zorder=4)
    for xi, p in zip(x, preds):
        ax.plot([xi - 0.33, xi + 0.33], [p, p], color=INK, lw=2.0, ls=(0, (3, 2)),
                zorder=5)
    ax.axhline(0, color=MUTED, lw=0.8, zorder=1)
    ax.set_xticks(x)
    ax.set_xticklabels(rs)
    hukum_ad = ""
    if ad == "+log10":
        hukum_ad = H["hukum"]["H-193a"]
    elif ad == "+log7":
        hukum_ad = H["hukum"]["H-193b"]
    else:
        hepsi = H["KAYIT"].get(ad, {}).get("hepsi_bantta")
        hukum_ad = "KAYIT" if hepsi is None else ("KAYIT (bantta)" if hepsi else "KAYIT (bant dışı)")
    ax.set_title(f"{ad}\nκ_top={h['kappa_top']:+.5f}±{h['se_kappa_top']:.5f}\n{hukum_ad}",
                fontsize=7.6, color=INK)
    ax.set_ylim(min(-0.2, (s_vals - se_vals).min() * 1.15, preds.min() * 1.15),
               max(0.2, (s_vals + se_vals).max() * 1.15, preds.max() * 1.15))
axes[0].set_ylabel("s_r = κ_r / κ_top")

from matplotlib.lines import Line2D
leg = [Line2D([], [], color=INK, lw=2, ls=(0, (3, 2)), label="öngörü (pay yasası)"),
       plt.Rectangle((0, 0), 1, 1, color=YESIL, alpha=0.88, label="yön+nicel OK"),
       plt.Rectangle((0, 0), 1, 1, color=TURUNCU, alpha=0.88, label="yalnız yön OK"),
       plt.Rectangle((0, 0), 1, 1, color=KIRMIZI, alpha=0.88, label="yön HAYIR / bant dışı"),
       plt.Rectangle((0, 0), 1, 1, color=MAVI, alpha=0.88, label="yan-KAYIT (hükme girmez)")]
fig.legend(handles=leg, loc="lower center", ncol=5, fontsize=8, frameon=False,
          bbox_to_anchor=(0.5, -0.02))
fig.suptitle("193 — Kalıntı-sınıfı yasası: s_r = κ_r/κ_top ölçülen (çubuk ± jk se) vs "
            "cos(2πrb/a)/μ(a) öngörü (kesikli)  [δ=0.03, blok-yerel Δω, kör sınav]",
            fontsize=11.5, color=INK, x=0.02, ha="left")
fig.subplots_adjust(left=0.045, right=0.985, top=0.80, bottom=0.20, wspace=0.32)
out = QM / "193_sinif_yasasi.png"
fig.savefig(out, dpi=150)
print(f"-> {out}")
