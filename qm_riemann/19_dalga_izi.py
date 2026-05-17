"""
19 — Riemann'ın Dalga İzi: Z(t)
==================================

Uğur'un sorusu: "milyarlarca sıfır var ve bu dalga izine bakılabilir"

Riemann-Siegel Z fonksiyonu:
    Z(t) = e^{iθ(t)} · ζ(1/2 + it)   (GERÇEK değerli)

Bu, ζ'nın σ=1/2 üzerindeki "döndürülmüş gerçel hali" — yani bir dalga.
Z(t) sıfırları ↔ ζ sıfırları (bire bir aynı yerde).

Çiziyoruz:
  - Z(t) için t ∈ [0, 100]: dalga görünür
  - Sıfırlar = sıfır geçişleri (kırmızı işaret)
  - Zirveler = lokal maks/min (yeşil)
  - Zirvelerin büyüklük dağılımı: bazı zirveler büyük (kuvvetli sıfır komşu),
    bazıları küçük (zayıf) → senin "dalga izi" sezgin
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import mpmath as mp
from scipy.signal import find_peaks

mp.mp.dps = 25

ROOT = Path(__file__).resolve().parent.parent
gamma = np.loadtxt(ROOT / "zeros_100k.txt")

# ============================================================
# Z(t) dalgası
# ============================================================
t_max = 100.0
N = 5000
t_grid = np.linspace(0.5, t_max, N)

print(f"Z(t) hesaplanıyor (mpmath siegelz, ~30 saniye)...")
Z = np.array([float(mp.siegelz(tt)) for tt in t_grid])
print(f"Hesaplandı. Z aralığı: [{Z.min():.2f}, {Z.max():.2f}]")

# Sıfır geçişleri (Z=0 olduğu yerler)
sign_changes = np.where(np.diff(np.sign(Z)) != 0)[0]
zero_crossings = t_grid[sign_changes]
print(f"Z'nin sıfır geçişleri: {len(zero_crossings)} (beklenen ≈ Riemann sıfır sayısı)")

# Gerçek Riemann sıfırları (karşılaştırma için)
gamma_window = gamma[gamma < t_max]
print(f"Bu pencerede gerçek Riemann sıfırı: {len(gamma_window)}")

# Zirveler
peaks_pos, _ = find_peaks(Z, distance=3)
peaks_neg, _ = find_peaks(-Z, distance=3)

# ============================================================
# Plot — üç katman
# ============================================================
fig = plt.figure(figsize=(14, 11))
gs = fig.add_gridspec(3, 2, height_ratios=[1.4, 1, 1], hspace=0.32, wspace=0.25)

# (1) Z(t) DALGASI — geniş
ax = fig.add_subplot(gs[0, :])
ax.plot(t_grid, Z, color="navy", lw=0.9, alpha=0.9)
ax.fill_between(t_grid, Z, 0, where=(Z > 0), alpha=0.18, color="blue")
ax.fill_between(t_grid, Z, 0, where=(Z < 0), alpha=0.18, color="red")
ax.axhline(0, color="black", lw=0.5)
# Sıfırlar
for g in gamma_window:
    ax.plot(g, 0, "o", ms=5, color="red", mfc="white", mec="red", mew=1)
# Zirveler
ax.plot(t_grid[peaks_pos], Z[peaks_pos], "^", ms=4, color="darkgreen", alpha=0.6)
ax.plot(t_grid[peaks_neg], Z[peaks_neg], "v", ms=4, color="darkgreen", alpha=0.6)

ax.set_xlim(0, t_max)
ax.set_xlabel("t")
ax.set_ylabel("Z(t)")
ax.set_title("Riemann-Siegel Z(t) — ζ'nın σ=1/2 üzerindeki GERÇEL dalgası\n"
             f"kırmızı halkalar = {len(gamma_window)} sıfır (Z=0 geçişleri), "
             "yeşil üçgenler = zirveler")
ax.grid(alpha=0.3)

# (2) Zirvelerin büyüklük dağılımı
ax = fig.add_subplot(gs[1, 0])
peak_heights = np.concatenate([Z[peaks_pos], -Z[peaks_neg]])
ax.hist(peak_heights, bins=40, color="seagreen", alpha=0.7, edgecolor="darkgreen")
ax.set_xlabel("|Z(zirve)|")
ax.set_ylabel("zirve sayısı")
ax.set_title(f"Zirvelerin büyüklük dağılımı (n={len(peak_heights)})\n"
             f"ortalama={peak_heights.mean():.2f}, "
             f"max={peak_heights.max():.2f}")
ax.grid(alpha=0.3)

# (3) Sıfır geçişleri vs gerçek Riemann sıfırları
ax = fig.add_subplot(gs[1, 1])
ax.plot(zero_crossings, np.zeros_like(zero_crossings), "o", ms=8, color="red",
        mfc="white", mec="red", mew=1.5, label=f"Z'nin sıfır geçişleri ({len(zero_crossings)})")
ax.plot(gamma_window, np.zeros_like(gamma_window) + 0.05, "s", ms=6,
        color="black", alpha=0.6, label=f"Gerçek γ_n ({len(gamma_window)})")
ax.set_xlim(0, t_max)
ax.set_ylim(-0.3, 0.5)
ax.set_xlabel("t")
ax.set_yticks([])
ax.set_title("Z'nin sıfır geçişleri ↔ gerçek Riemann sıfırları\n"
             "(bire bir aynı yerde olmalı)")
ax.legend(fontsize=9, loc="upper right")
ax.grid(alpha=0.3)

# (4) Komşu zirveler arası mesafe
ax = fig.add_subplot(gs[2, 0])
peak_positions = np.sort(np.concatenate([t_grid[peaks_pos], t_grid[peaks_neg]]))
peak_gaps = np.diff(peak_positions)
ax.hist(peak_gaps, bins=40, color="orange", alpha=0.7, edgecolor="darkorange")
ax.set_xlabel("zirveler arası mesafe (t birimi)")
ax.set_ylabel("sayı")
ax.set_title(f"Komşu zirveler arası mesafe (n={len(peak_gaps)})\n"
             f"ortalama={peak_gaps.mean():.3f}")
ax.grid(alpha=0.3)

# (5) Sıfır geçişleri arası mesafe (gerçek sıfır spacing)
ax = fig.add_subplot(gs[2, 1])
zero_gaps = np.diff(gamma_window)
ax.hist(zero_gaps, bins=40, color="crimson", alpha=0.7, edgecolor="darkred")
ax.set_xlabel("ardışık sıfırlar arası mesafe")
ax.set_ylabel("sayı")
ax.set_title(f"Ardışık Riemann sıfırları arası mesafe (n={len(zero_gaps)})\n"
             f"ortalama={zero_gaps.mean():.3f}")
ax.grid(alpha=0.3)

plt.suptitle("19 — Riemann'ın Dalga İzi: Z(t) ve Sıfırların Düşüş Deseni",
             fontsize=13, fontweight="bold", y=0.995)

out = Path(__file__).parent / "19_dalga_izi.png"
plt.savefig(out, dpi=130, bbox_inches="tight")
print(f"\nKaydedildi: {out}")

# ============================================================
# Sayısal — Uğur'un sezgisi: yoğunluk dağılımı
# ============================================================
print(f"\n=== DALGA İSTATİSTİĞİ ===")
print(f"Z dalga genliği:    [{Z.min():+.2f}, {Z.max():+.2f}]")
print(f"Zirve büyüklükleri: ort={peak_heights.mean():.2f}, std={peak_heights.std():.2f}")
print(f"  En büyük zirve t = {t_grid[peaks_pos[np.argmax(Z[peaks_pos])]]:.2f}")
print(f"  En derin dip   t = {t_grid[peaks_neg[np.argmax(-Z[peaks_neg])]]:.2f}")
print(f"Komşu zirveler arası mesafe: ort={peak_gaps.mean():.3f}, std={peak_gaps.std():.3f}")
print(f"Komşu sıfırlar arası mesafe: ort={zero_gaps.mean():.3f}, std={zero_gaps.std():.3f}")
print(f"\nGözlem: dalga her t'de farklı genlikte salınıyor.")
print(f"Bazı yerlerde 'küçük dalgalar' (zayıf zirve), bazı yerlerde 'büyük dalgalar' (güçlü zirve).")
print(f"Bu küçük/büyük örüntüsü = senin 'milyarlarca sıfırın dalga izi' sezgin.")
