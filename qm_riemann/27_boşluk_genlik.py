"""
27 — Boşluk-Genlik İlişkisi
==============================

Uğur'un gözlemi: "Riemann plakasında büyük boşlukların olduğu anlar var.
O boşluklarda arkadaki dalgalar özellikle peak yapıyor. Bilinçli tercih gibi."

Test: Ardışık sıfır çifti için
  - Aralık = γ_{n+1} - γ_n
  - Aralıktaki Z(t)'nin max |genliği| = M_n

İki nicelik arasındaki ilişki: pozitif korelasyon var mı?
Eğer evet → Uğur'un sezgisi sayısal doğrulanır.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import mpmath as mp
from scipy.stats import pearsonr, spearmanr

mp.mp.dps = 20

ROOT = Path(__file__).resolve().parent.parent
gamma = np.loadtxt(ROOT / "zeros_100k.txt")

# ============================================================
# Çalışma penceresi
# ============================================================
t_max = 1000.0
gamma_w = gamma[gamma < t_max]
print(f"Pencere içindeki sıfır: {len(gamma_w)}")
print(f"İlk: {gamma_w[0]:.3f}, son: {gamma_w[-1]:.3f}")

# ============================================================
# Her ardışık çift için: aralık + max |Z| ölçü
# ============================================================
print("\nHer aralıkta max |Z(t)| hesaplanıyor...")
intervals = np.diff(gamma_w)        # γ_{n+1} - γ_n
max_amplitudes = np.zeros(len(intervals))

for n in range(len(intervals)):
    g_lo, g_hi = gamma_w[n], gamma_w[n+1]
    # 30 nokta arada
    t_sub = np.linspace(g_lo + 0.01, g_hi - 0.01, 30)
    Z_sub = np.array([float(mp.siegelz(tt)) for tt in t_sub])
    max_amplitudes[n] = np.max(np.abs(Z_sub))
    if n % 100 == 0:
        print(f"  {n}/{len(intervals)}: aralık={intervals[n]:.2f}, max|Z|={max_amplitudes[n]:.2f}")

print(f"\nToplam çift: {len(intervals)}")
print(f"Aralık istatistik: ort={intervals.mean():.3f}, std={intervals.std():.3f}, min={intervals.min():.3f}, max={intervals.max():.3f}")
print(f"max|Z|  istatistik: ort={max_amplitudes.mean():.3f}, std={max_amplitudes.std():.3f}, min={max_amplitudes.min():.3f}, max={max_amplitudes.max():.3f}")

# ============================================================
# Korelasyon testi
# ============================================================
r_pearson, p_pearson = pearsonr(intervals, max_amplitudes)
r_spearman, p_spearman = spearmanr(intervals, max_amplitudes)

print(f"\n=== KORELASYON ===")
print(f"Pearson  r = {r_pearson:+.4f}, p = {p_pearson:.4e}")
print(f"Spearman r = {r_spearman:+.4f}, p = {p_spearman:.4e}")

# Lokal yoğunluğa göre normalize: t ile log etkisini çıkar
# Normalize aralık: aralık × ρ(t) burada ρ(t) ≈ log(t/2π)/(2π)
t_mid = 0.5 * (gamma_w[:-1] + gamma_w[1:])
rho_local = np.log(t_mid / (2*np.pi)) / (2*np.pi)
intervals_normalized = intervals * rho_local

r_pearson_norm, _ = pearsonr(intervals_normalized, max_amplitudes)
print(f"\nNormalize aralık (yerel yoğunluk hariç):")
print(f"Pearson r = {r_pearson_norm:+.4f}")

# ============================================================
# Plot
# ============================================================
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# (1) Saçılım grafiği: aralık vs max |Z|
ax = axes[0, 0]
ax.scatter(intervals, max_amplitudes, s=12, alpha=0.5, color="steelblue")
# Trend çizgisi (linear fit)
z = np.polyfit(intervals, max_amplitudes, 1)
x_fit = np.linspace(intervals.min(), intervals.max(), 100)
ax.plot(x_fit, np.polyval(z, x_fit), "r-", lw=2,
        label=f"linear fit: y = {z[0]:.2f}·x + {z[1]:.2f}")
ax.set_xlabel("aralık γ_{n+1} − γ_n")
ax.set_ylabel("aralıkta max |Z(t)|")
ax.set_title(f"AĞRALIK vs GENLİK — pearson r = {r_pearson:.3f}, p = {p_pearson:.2e}")
ax.legend(fontsize=10)
ax.grid(alpha=0.3)

# (2) Normalize aralık vs genlik
ax = axes[0, 1]
ax.scatter(intervals_normalized, max_amplitudes, s=12, alpha=0.5, color="seagreen")
z2 = np.polyfit(intervals_normalized, max_amplitudes, 1)
x_fit2 = np.linspace(intervals_normalized.min(), intervals_normalized.max(), 100)
ax.plot(x_fit2, np.polyval(z2, x_fit2), "r-", lw=2,
        label=f"r_pearson = {r_pearson_norm:.3f}")
ax.set_xlabel("normalize aralık (lokal yoğunluk dahil)")
ax.set_ylabel("max |Z(t)|")
ax.set_title("NORMALİZE AĞRALIK vs GENLİK — t etkisi çıkarıldı")
ax.legend(fontsize=10)
ax.grid(alpha=0.3)

# (3) Z(t) + sıfırlar + en yüksek genlik noktaları
ax = axes[1, 0]
t_show_max = 200
t_show = np.linspace(0.5, t_show_max, 3000)
Z_show = np.array([float(mp.siegelz(tt)) for tt in t_show])
ax.plot(t_show, Z_show, color="navy", lw=0.7, alpha=0.85)
ax.axhline(0, color="black", lw=0.5)
# Sıfırlar
for g in gamma_w[gamma_w < t_show_max]:
    ax.axvline(g, color="red", lw=0.3, alpha=0.5)
# En büyük 5 aralığın merkezi
top_5_idx = np.argsort(intervals)[-5:]
for idx in top_5_idx:
    if gamma_w[idx] < t_show_max and gamma_w[idx+1] < t_show_max:
        mid = (gamma_w[idx] + gamma_w[idx+1])/2
        ax.axvspan(gamma_w[idx], gamma_w[idx+1], alpha=0.2, color="gold",
                   label="büyük boşluk" if idx == top_5_idx[-1] else "")
ax.set_xlim(0, t_show_max)
ax.set_xlabel("t")
ax.set_ylabel("Z(t)")
ax.set_title("Z(t) — kırmızı çizgi: sıfır, altın bant: en büyük boşluklar\n"
             "büyük boşluklarda dalga gerçekten yüksek mi?")
ax.legend(fontsize=9)
ax.grid(alpha=0.3)

# (4) Dağılım: en büyük aralıklarda max|Z| dağılımı vs ortalama
ax = axes[1, 1]
sorted_intervals = np.argsort(intervals)
# En küçük %25 vs en büyük %25 vs ortadaki
n_q = len(intervals) // 4
small_intervals_amps = max_amplitudes[sorted_intervals[:n_q]]
large_intervals_amps = max_amplitudes[sorted_intervals[-n_q:]]
middle_intervals_amps = max_amplitudes[sorted_intervals[n_q:-n_q]]
bins = np.linspace(0, max_amplitudes.max(), 30)
ax.hist(small_intervals_amps, bins=bins, alpha=0.5, label="küçük aralık (alt %25)", color="green")
ax.hist(middle_intervals_amps, bins=bins, alpha=0.5, label="orta aralık", color="grey")
ax.hist(large_intervals_amps, bins=bins, alpha=0.5, label="büyük aralık (üst %25)", color="crimson")
ax.set_xlabel("max |Z| dağılımı")
ax.set_ylabel("sayı")
ax.set_title("Aralık büyüklüğüne göre max|Z| dağılımı\n"
             f"büyük aralık ort = {large_intervals_amps.mean():.2f}, küçük ort = {small_intervals_amps.mean():.2f}")
ax.legend(fontsize=9)
ax.grid(alpha=0.3)

plt.suptitle(f"27 — Uğur'un Sezgisi Testi: Aralık ↔ Genlik\n"
             f"Pearson r = {r_pearson:.3f}, n = {len(intervals)} çift",
             fontsize=13, fontweight="bold", y=0.998)
plt.tight_layout()

out = Path(__file__).parent / "27_bosluk_genlik.png"
plt.savefig(out, dpi=130, bbox_inches="tight")
print(f"\nKaydedildi: {out}")
