"""
101g — DONMUŞ KOYLAR FİGÜRÜ (20 Ağustos)
==========================================================================
Veri: 101d çıktısı (sertifikalı-düz segmentler; ζ referans + 4 ada).
Mesaj: çözülme sınırı adanın İLK SAĞ KALAN çizgisidir — β'nın 2-çizgisi
ölü (χ₄(2)=0) → donma log3/L'e uzar; ötekiler log2/L'de çözülür.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

TAU = [0.040, 0.050, 0.060, 0.068, 0.075, 0.085, 0.095, 0.105, 0.113,
       0.125, 0.140, 0.160, 0.200, 0.300]
D = {
    "zeta": [0.041, 0.181, 0.240, 0.341, 0.482, 0.530, 0.642, 0.676,
             0.728, 0.788, 0.855, 0.903, 0.964, 1.052],
    "chi3": [0.069, 0.130, 0.314, 0.394, 0.533, 0.594, 0.605, 0.545,
             0.489, 0.583, 0.727, 0.817, 0.939, 1.007],
    "beta": [0.010, 0.016, 0.026, 0.040, 0.054, 0.116, 0.460, 0.520,
             0.705, 0.817, 0.878, 0.884, 0.909, 1.034],
    "chi5": [0.030, 0.112, 0.242, 0.350, 0.429, 0.517, 0.594, 0.653,
             0.711, 0.779, 0.781, 0.647, 0.920, 1.034],
    "chi7": [0.078, 0.220, 0.366, 0.465, 0.529, 0.619, 0.686, 0.737,
             0.783, 0.794, 0.863, 0.903, 0.927, 1.055],
}
VEKIL = [0.010, 0.015, 0.021, 0.028, 0.028, 0.036, 0.036, 0.031, 0.040,
         0.055, 0.049, 0.061, 0.056, 0.088]

fig, ax = plt.subplots(figsize=(8.4, 5.6))
ax.fill_between(TAU, 0, VEKIL, color="0.85", zorder=0,
                label="vekil taban (gürültü)")
stil = {"zeta": ("k", "-", "o", r"$\zeta$"),
        "chi3": ("#2a7", "-", "s", r"$\chi_3$"),
        "chi5": ("#27b", "-", "^", r"$\chi_5$"),
        "chi7": ("#a6c", "-", "v", r"$\chi_7$"),
        "beta": ("#d33", "-", "D", r"$\beta\ (\chi_4)$")}
for k in ["zeta", "chi3", "chi5", "chi7", "beta"]:
    c, ls, mk, lb = stil[k]
    lw = 2.6 if k == "beta" else 1.4
    ax.plot(TAU, D[k], ls, color=c, marker=mk, ms=4.5, lw=lw, label=lb)

ax.axvline(0.070, color="0.4", ls="--", lw=1)
ax.axvline(0.113, color="#d33", ls="--", lw=1)
ax.text(0.070, 1.12, r"$\log 2/L$" + "\n(ilk çizgi: ζ, χ₃, χ₅, χ₇)",
        ha="center", fontsize=8.5, color="0.3")
ax.text(0.113, 0.28, r"$\log 3/L$" + "\n(β'nın ilk\nsağ kalan çizgisi)",
        ha="left", fontsize=8.5, color="#d33")
ax.annotate("β burada hâlâ donuk:\n2-çizgisi ölü (χ₄(2)=0)",
            xy=(0.075, 0.054), xytext=(0.045, 0.62),
            arrowprops=dict(arrowstyle="->", color="#d33"),
            fontsize=9, color="#d33")
ax.set_xlabel(r"$\tau = \omega/L$")
ax.set_ylabel(r"$D_{\rm spont}(\tau)$  (spontane gap-yanıtı, CUE$\approx$1)")
ax.set_title("Donmuş koylar: çözülme sınırı adanın ilk sağ kalan çizgisidir\n"
             "(101d — sertifikalı-düz segmentler, 5 aile)")
ax.legend(loc="lower right", fontsize=9)
ax.set_xlim(0.03, 0.31); ax.set_ylim(0, 1.22)
ax.grid(alpha=0.25)
fig.tight_layout()
fig.savefig("101_donmus_koylar.png", dpi=150)
print("kaydedildi: 101_donmus_koylar.png")
