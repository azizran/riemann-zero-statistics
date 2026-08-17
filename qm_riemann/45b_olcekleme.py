"""
45b — w(p,L) ÖLÇEKLEME ÇÖKMESİ: TEK DEĞİŞKEN τ = log p / L (17 Ağustos 2026)
=============================================================================

45'in ham tablosu: w hem p ile düşüyor hem L ile yükseliyor.
Hipotez: ikisi tek değişkende birleşir — τ = log p / L
(Berry 1988 form-faktör değişkeni; sıfır gazının asal dalgayı
"perdeleme" etkinliği bu orana bağlı olmalı).

Veri: 45'in çıktı tablosu (yeniden üretilebilir: 45_wp_kalibrasyon.py).
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

HERE = Path(__file__).resolve().parent
PRIMES = [2, 3, 5, 7, 11, 13]
LS = [9.86, 10.37, 10.93, 11.47, 11.98, 12.45]
# 45 çıktısı: satır=L, sütun=p
W = np.array([
    [0.742, 0.616, 0.459, 0.363, 0.242, 0.198],
    [0.755, 0.633, 0.478, 0.384, 0.264, 0.223],
    [0.762, 0.644, 0.497, 0.410, 0.289, 0.253],
    [0.771, 0.657, 0.514, 0.427, 0.318, 0.275],
    [0.778, 0.665, 0.526, 0.441, 0.334, 0.293],
    [0.781, 0.675, 0.544, 0.454, 0.345, 0.311],
])

fig, axes = plt.subplots(1, 2, figsize=(12.5, 5))

# sol: ham eğriler
ax = axes[0]
for j, p in enumerate(PRIMES):
    ax.plot(LS, W[:, j], "o-", ms=4, label=f"p={p}")
ax.set_xlabel("L"); ax.set_ylabel("w(p, L)")
ax.set_title("Ham: iki değişken"); ax.legend(fontsize=8, ncol=2); ax.grid(alpha=0.3)

# sağ: çökme
ax = axes[1]
xs_all, ws_all = [], []
for i, L in enumerate(LS):
    xs = [np.log(p) / L for p in PRIMES]
    ax.plot(xs, W[i], "o", ms=5, alpha=0.75, label=f"L={L}")
    xs_all += xs; ws_all += list(W[i])
xs_all = np.array(xs_all); ws_all = np.array(ws_all)

# ampirik düz fit (kapalı form İDDİASI DEĞİL — sadece göz kılavuzu)
c = np.polyfit(xs_all, ws_all, 1)
xx = np.linspace(0, 0.28, 50)
ax.plot(xx, np.polyval(c, xx), "k--", lw=1,
        label=f"kılavuz: {c[1]:.2f} {c[0]:+.2f}·τ")
resid = ws_all - np.polyval(c, xs_all)
ax.set_xlabel("τ = log p / L"); ax.set_ylabel("w")
ax.set_title(f"Çökme: tek değişken (lineer artık RMS {np.sqrt((resid**2).mean()):.3f})")
ax.legend(fontsize=7, ncol=2); ax.grid(alpha=0.3)

plt.tight_layout()
out = HERE / "45b_olcekleme.png"
plt.savefig(out, dpi=110)

print(f"Lineer kılavuz: w ≈ {c[1]:.3f} {c[0]:+.3f}·τ, artık RMS = {np.sqrt((resid**2).mean()):.4f}")
print(f"τ→0 ekstrapolasyonu: w(0) ≈ {c[1]:.3f}  (1'e gitmesi 'tam iletim' demek olurdu)")
print("NOT: lineer form kapalı-form iddiası DEĞİL; sadece çökmenin kalitesini ölçüyor.")
print(f"Grafik: {out.name}")
