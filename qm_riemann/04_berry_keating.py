"""
QM × Riemann — Berry-Keating Testi
====================================

Berry-Keating (1999): Hamiltonian H = (xp + px)/2 önerdiler. Klasik faz uzayı
hiperbolik (H = xp = sabit eğrileri), sınırlar: x ≥ ℓ_x = √(2πℏ), p ≥ ℓ_p = √(2πℏ).

Yarı-klasik Weyl yasası → bağlı durum sayısı:

    N_BK(E) = (E/2π) · log(E/(2π·e)) + 7/8

Bu Riemann sıfır sayım fonksiyonu N(T)'nin **TAM düzgün kısmıdır** (Riemann-Siegel).
Yani Berry-Keating'in iddiası: salınım = asallar, ortalama = QM.

Test: Odlyzko verisinden gerçek N(T) - N_BK(T)'ye bak.
  - Düzgün uyum varsa: BK'nin yarı-klasik tarafı doğru
  - Kalan salınım = ψ(x) tarzı asal toplamları olmalı
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
t = np.loadtxt(ROOT / "zeros_100k.txt")
N_emp = np.arange(1, len(t) + 1)   # N(t_n) = n

# ---------- Berry-Keating düzgün sayma ----------
def N_BK(T):
    """Berry-Keating yarı-klasik = Riemann-Siegel düzgün kısmı"""
    return (T / (2*np.pi)) * np.log(T / (2*np.pi*np.e)) + 7/8

N_pred = N_BK(t)
residual = N_emp - N_pred   # gerçek - öngörü

# ---------- Salınımın asal kaynağı: Riemann explicit formula ----------
# N(T) - N_BK(T) = (1/π) Σ_n Im[Li(T^ρ_n)] + ...  (salınım)
# Bizim test için: salınımın yapısına bak (rastgele değil, log-doğrusal yapı)
# Pratikte: residual'in (−1/π)·arg ζ(½+iT) ile özdeşliği
# Bunu doğrulamak için S(T) = (1/π) Im log ζ(½+iT) eşitliği bilinen.

# ---------- Plot ----------
fig, axes = plt.subplots(2, 1, figsize=(13, 8),
                         gridspec_kw={"height_ratios": [2, 1.3]})

# Üst: N(T) vs N_BK(T)
ax = axes[0]
# Tüm sıfırları çizmek yoğunluk yaratır → log eksende seyrelt
idx = np.unique(np.geomspace(1, len(t)-1, 800).astype(int))
ax.plot(t[idx], N_emp[idx], "o", ms=2.5, color="steelblue",
        label=f"ζ sıfır sayımı (gerçek, N={len(t)})")
T_grid = np.geomspace(t[0], t[-1], 500)
ax.plot(T_grid, N_BK(T_grid), "r-", lw=2,
        label=r"Berry-Keating: $\frac{T}{2\pi}\log\frac{T}{2\pi e}+\frac{7}{8}$")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("T")
ax.set_ylabel("N(T)")
ax.set_title("Berry-Keating yarı-klasik Weyl yasası vs gerçek ζ sıfır sayısı")
ax.legend(loc="lower right", fontsize=10)
ax.grid(alpha=0.3, which="both")

# Alt: residual = oscillating part = S(T) ⋅ ?
ax = axes[1]
# tüm noktaları göstermek mantıksız (100k), pencere içinde seçim:
T_window_max = 1000
mask = t < T_window_max
ax.plot(t[mask], residual[mask], color="steelblue", lw=0.6, alpha=0.85)
ax.axhline(0, color="r", lw=1.2, label="BK öngörüsü (salınım yok)")
ax.set_xlim(0, T_window_max)
ax.set_xlabel("T")
ax.set_ylabel(r"$N(T) - N_{BK}(T) = S(T) + \frac{1}{2}$")
ax.set_title("Kalan salınım: bu sıçramalar her sıfırda +1 — asalların izi")
ax.legend(loc="upper right", fontsize=9)
ax.grid(alpha=0.3)

plt.tight_layout()
out = Path(__file__).parent / "04_berry_keating.png"
plt.savefig(out, dpi=130, bbox_inches="tight")
print(f"Kaydedildi: {out}")

# ---------- Sayısal: ne kadar iyi uyuyor? ----------
# Yüksek T'de mutlak hata orantısal olarak küçülür mü?
print(f"\nBerry-Keating düzgün öngörünün performansı:")
print(f"{'T aralığı':>20} | {'|residual| mean':>17} | {'std':>10} | {'/(log T)':>10}")
for lo, hi in [(10, 100), (100, 1000), (1000, 10000), (10000, 75000)]:
    m = (t > lo) & (t < hi)
    r = residual[m]
    print(f"  [{lo:>5}, {hi:>6}]   |  {np.mean(np.abs(r)):>15.4f}  | {r.std():>10.4f} | {r.std()/np.sqrt(np.log(0.5*(lo+hi))):>10.4f}")

# Selberg: <S(T)²> ~ (1/π²) log log T   →  std(residual) ~ √(log log T) / π
# Çok yavaş büyüme. Test et.
print(f"\nSelberg öngörüsü: std(S(T)) ~ √(log log T)/π")
for lo, hi in [(100, 1000), (1000, 10000), (10000, 75000)]:
    m = (t > lo) & (t < hi)
    Tm = 0.5*(lo+hi)
    selberg_std = np.sqrt(np.log(np.log(Tm))) / np.pi
    # residual = S(T) + 1/2  →  std olarak aynı
    actual_std = residual[m].std()
    print(f"  T~{Tm:>6.0f}:  öngörü {selberg_std:.4f}  vs  gerçek {actual_std:.4f}")
