"""
20 — Pozitif Kayma Evrensel mi? (S1)
========================================

Uğur'un sezgisi: "+1 kısmına kayış her yerde var. her sayı biçiminde
doğanın bir yönü var."

Test: Z(t) dalgasında pozitif zirveler ile negatif dipler
  - Aynı sayıda mı?
  - Aynı büyüklükte mi?
  - Dağılımları özdeş mi?

Beklenti: fonksiyonel denklem ξ(s)=ξ(1-s) simetri zorlar.
Eğer doğanın yönü gerçekten varsa, küçük ama sistemik bir asimetri olmalı.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import mpmath as mp
from scipy.signal import find_peaks
from scipy.stats import ks_2samp

mp.mp.dps = 25

# ============================================================
# Z(t) için geniş hesap
# ============================================================
t_lo, t_hi = 10.0, 300.0   # geniş pencere
N = 12000                  # çözünürlük ~0.024
t_grid = np.linspace(t_lo, t_hi, N)

print(f"Z(t) hesaplanıyor: t ∈ [{t_lo}, {t_hi}], {N} nokta (~2 dakika)...")
Z = np.zeros(N)
for i in range(N):
    Z[i] = float(mp.siegelz(t_grid[i]))
    if i % 1000 == 0:
        print(f"  i={i}/{N}, t={t_grid[i]:.1f}")
print(f"Z aralığı: [{Z.min():+.3f}, {Z.max():+.3f}]")

# ============================================================
# Pozitif zirveler ve negatif dipler
# ============================================================
peaks_pos, _ = find_peaks(Z, distance=3, height=0)        # lokal max ve Z>0
peaks_neg, _ = find_peaks(-Z, distance=3, height=0)       # lokal min ve Z<0

pos_heights = Z[peaks_pos]
neg_depths  = -Z[peaks_neg]    # pozitif değerler (büyüklük)

print(f"\nPozitif zirve sayısı: {len(pos_heights)}")
print(f"Negatif dip sayısı:   {len(neg_depths)}")
print(f"\nPozitif zirve istatistik:  ort={pos_heights.mean():.4f}, std={pos_heights.std():.4f}, max={pos_heights.max():.4f}")
print(f"Negatif dip istatistik:    ort={neg_depths.mean():.4f}, std={neg_depths.std():.4f}, max={neg_depths.max():.4f}")
print(f"\nFARK (pozitif - negatif): ort={pos_heights.mean()-neg_depths.mean():+.4f}")

# Kolmogorov-Smirnov: iki dağılım özdeş mi?
ks_stat, p_value = ks_2samp(pos_heights, neg_depths)
print(f"\nKS testi: stat={ks_stat:.4f}, p-değer={p_value:.4f}")
if p_value < 0.05:
    print("  → DAĞILIMLAR FARKLI (anlamlı asimetri)")
else:
    print("  → DAĞILIMLAR ÖZDEŞ (asimetri yok)")

# Toplam pozitif alan vs negatif alan
pos_area = np.trapz(np.maximum(Z, 0), t_grid)
neg_area = np.trapz(np.maximum(-Z, 0), t_grid)
print(f"\nToplam pozitif alan (∫Z⁺ dt): {pos_area:.3f}")
print(f"Toplam negatif alan (∫Z⁻ dt): {neg_area:.3f}")
print(f"Fark / toplam: {(pos_area - neg_area)/(pos_area + neg_area)*100:+.3f}%")

# Pencerelere böl, lokal asimetri var mı bak
n_windows = 10
window_size = (t_hi - t_lo) / n_windows
window_asymmetries = []
for i in range(n_windows):
    t_start = t_lo + i * window_size
    t_end = t_start + window_size
    mask = (t_grid >= t_start) & (t_grid < t_end)
    Zw = Z[mask]
    pa = np.trapz(np.maximum(Zw, 0), t_grid[mask])
    na = np.trapz(np.maximum(-Zw, 0), t_grid[mask])
    asym = (pa - na) / (pa + na) if (pa + na) > 0 else 0
    window_asymmetries.append((t_start, asym))
    print(f"  Pencere [{t_start:.0f}, {t_end:.0f}]: asimetri = {asym*100:+.2f}%")

# ============================================================
# Plot
# ============================================================
fig = plt.figure(figsize=(14, 11))
gs = fig.add_gridspec(3, 2, height_ratios=[1.3, 1, 1], hspace=0.32, wspace=0.25)

# (1) Z(t) parça — gösterim
ax = fig.add_subplot(gs[0, :])
mask_show = (t_grid >= 10) & (t_grid <= 100)
ax.plot(t_grid[mask_show], Z[mask_show], color="navy", lw=0.8)
ax.fill_between(t_grid[mask_show], Z[mask_show], 0,
                where=(Z[mask_show] > 0), alpha=0.25, color="blue")
ax.fill_between(t_grid[mask_show], Z[mask_show], 0,
                where=(Z[mask_show] < 0), alpha=0.25, color="red")
ax.axhline(0, color="black", lw=0.5)
# zirveleri işaretle
mask_peaks_pos = t_grid[peaks_pos] < 100
mask_peaks_neg = t_grid[peaks_neg] < 100
ax.plot(t_grid[peaks_pos][mask_peaks_pos], pos_heights[mask_peaks_pos],
        "^", ms=6, color="darkgreen", label=f"pozitif zirve ({len(pos_heights)} toplam)")
ax.plot(t_grid[peaks_neg][mask_peaks_neg], -neg_depths[mask_peaks_neg],
        "v", ms=6, color="darkred", label=f"negatif dip ({len(neg_depths)} toplam)")
ax.set_xlim(10, 100)
ax.set_xlabel("t")
ax.set_ylabel("Z(t)")
ax.set_title(f"Z(t) dalgası — pozitif (mavi) ile negatif (kırmızı) doldurulmuş\n"
             f"yeşil üçgenler = pozitif zirveler, kırmızı üçgenler = negatif dipler")
ax.legend(fontsize=9)
ax.grid(alpha=0.3)

# (2) Pozitif vs negatif büyüklük dağılımları
ax = fig.add_subplot(gs[1, 0])
bins = np.linspace(0, max(pos_heights.max(), neg_depths.max())*1.1, 50)
ax.hist(pos_heights, bins=bins, alpha=0.55, color="seagreen",
        label=f"pozitif zirveler (ort={pos_heights.mean():.3f})")
ax.hist(neg_depths, bins=bins, alpha=0.55, color="indianred",
        label=f"negatif dipler  (ort={neg_depths.mean():.3f})")
ax.set_xlabel("zirve/dip büyüklüğü |Z|")
ax.set_ylabel("sayı")
ax.set_title(f"Büyüklük dağılımları yan yana\nKS test p-değer: {p_value:.4f}")
ax.legend(fontsize=9)
ax.grid(alpha=0.3)

# (3) Pencere-pencere asimetri
ax = fig.add_subplot(gs[1, 1])
ws, asyms = zip(*window_asymmetries)
ax.bar(np.array(ws) + window_size/2, [a*100 for a in asyms],
       width=window_size*0.8, color="purple", alpha=0.7, edgecolor="darkmagenta")
ax.axhline(0, color="black", lw=0.5)
ax.set_xlabel("t (pencere başlangıcı)")
ax.set_ylabel("asimetri  (∫Z⁺ - ∫Z⁻)/(toplam)  %")
ax.set_title(f"Pencere-pencere asimetri — küçük ama sistemik bir kayma var mı?")
ax.grid(alpha=0.3)

# (4) Empirik CDF karşılaştırma
ax = fig.add_subplot(gs[2, 0])
pos_sorted = np.sort(pos_heights)
neg_sorted = np.sort(neg_depths)
ax.plot(pos_sorted, np.linspace(0, 1, len(pos_sorted)), color="seagreen", lw=1.5,
        label="pozitif CDF")
ax.plot(neg_sorted, np.linspace(0, 1, len(neg_sorted)), color="indianred", lw=1.5,
        label="negatif CDF")
ax.set_xlabel("|Z|")
ax.set_ylabel("kümülatif olasılık")
ax.set_title("Kümülatif dağılımlar — KS testi farkları görsel")
ax.legend(fontsize=10)
ax.grid(alpha=0.3)

# (5) Pozitif minus negatif — sıralı çiftleme
ax = fig.add_subplot(gs[2, 1])
N_compare = min(len(pos_heights), len(neg_depths))
pos_subset = np.sort(pos_heights)[:N_compare]
neg_subset = np.sort(neg_depths)[:N_compare]
diff = pos_subset - neg_subset
ax.plot(diff, color="purple", lw=1)
ax.axhline(0, color="black", lw=0.5)
ax.axhline(diff.mean(), color="red", lw=1.2, ls="--",
           label=f"ortalama fark = {diff.mean():+.4f}")
ax.set_xlabel("sıralı çift indeksi")
ax.set_ylabel("|pozitif zirve| − |negatif dip|  (sıralı eşleştirilmiş)")
ax.set_title(f"Sıralı eşleştirme: pozitif tarafta sistemik avantaj var mı?")
ax.legend(fontsize=10)
ax.grid(alpha=0.3)

plt.suptitle("20 — S1: Pozitif Kayma Evrensel mi?\n"
             f"Z(t) dalgasında pozitif ↔ negatif simetrisi  (t ∈ [{t_lo:.0f}, {t_hi:.0f}])",
             fontsize=12, fontweight="bold", y=0.998)

out = Path(__file__).parent / "20_pozitif_kayma.png"
plt.savefig(out, dpi=130, bbox_inches="tight")
print(f"\nKaydedildi: {out}")

# ============================================================
# Yorum
# ============================================================
print(f"\n=== YORUM ===")
print(f"Pozitif zirve ortalaması:  {pos_heights.mean():.4f}")
print(f"Negatif dip ortalaması:    {neg_depths.mean():.4f}")
print(f"Fark / ortalama büyüklük:  {(pos_heights.mean()-neg_depths.mean())/((pos_heights.mean()+neg_depths.mean())/2)*100:+.3f}%")
print(f"")
print(f"Toplam alan farkı:         {(pos_area-neg_area)/(pos_area+neg_area)*100:+.3f}%")
print(f"")
if abs(pos_heights.mean() - neg_depths.mean()) > 2 * pos_heights.std() / np.sqrt(len(pos_heights)):
    print("ANLAMLI ASİMETRİ VAR — Uğur'un sezgisi sayısal olarak desteklendi")
else:
    print("ASİMETRİ YOK ya da örneklem-içi gürültü — büyük pencere lazım")
