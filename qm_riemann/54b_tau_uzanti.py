"""
54b — τ YASASININ 160.000 KAT SIÇRAMASI (17 Ağustos 2026)
===========================================================

45'in w(p,L) noktaları (L=9.86..12.45, t≤1.6×10⁶) ile 54'ün Odlyzko
noktaları (L=24.48, t=2.7×10¹¹) aynı τ ekseninde. Kılavuz çizgi 45b'de
KÜÇÜK T'ye fit edilmişti — Odlyzko noktaları out-of-sample test.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

HERE = Path(__file__).resolve().parent
PRIMES = [2, 3, 5, 7, 11, 13]
LS_OLD = [9.86, 10.37, 10.93, 11.47, 11.98, 12.45]
W_OLD = np.array([
    [0.742, 0.616, 0.459, 0.363, 0.242, 0.198],
    [0.755, 0.633, 0.478, 0.384, 0.264, 0.223],
    [0.762, 0.644, 0.497, 0.410, 0.289, 0.253],
    [0.771, 0.657, 0.514, 0.427, 0.318, 0.275],
    [0.778, 0.665, 0.526, 0.441, 0.334, 0.293],
    [0.781, 0.675, 0.544, 0.454, 0.345, 0.311],
])
L_ODL = 24.4751
W_ODL = [0.871, 0.803, 0.718, 0.674, 0.609, 0.558]
SE_ODL = [0.007, 0.009, 0.011, 0.013, 0.017, 0.018]

fig, ax = plt.subplots(figsize=(8.5, 5.6))
for i, L in enumerate(LS_OLD):
    xs = [np.log(p) / L for p in PRIMES]
    ax.plot(xs, W_OLD[i], "o", ms=4, c="steelblue", alpha=0.55,
            label="t ≤ 1.6×10⁶ (36 nokta)" if i == 0 else None)
xs_o = [np.log(p) / L_ODL for p in PRIMES]
ax.errorbar(xs_o, W_ODL, yerr=SE_ODL, fmt="s", ms=7, c="firebrick", zorder=5,
            label="t = 2.7×10¹¹ (Odlyzko, out-of-sample)")
xx = np.linspace(0, 0.28, 50)
ax.plot(xx, 0.940 - 2.967 * xx, "k--", lw=1.2,
        label="45b kılavuzu (yalnız küçük-T'ye fit)")
ax.set_xlabel("τ = log p / L")
ax.set_ylabel("w — genlik kanalı iletimi")
ax.set_title("τ yasası, yükseklikte 160.000 kat sıçramayı atlıyor")
ax.legend(fontsize=9)
ax.grid(alpha=0.3)
plt.tight_layout()
out = HERE / "54b_tau_uzanti.png"
plt.savefig(out, dpi=110)
print(f"Grafik: {out.name}")
resid = np.array(W_ODL) - (0.940 - 2.967 * np.array(xs_o))
print("Odlyzko artıkları (kılavuza göre):", np.round(resid, 3))
