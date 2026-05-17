"""
28 — Çift Yarık × Riemann Paralel + Selberg Karşılaştırma
==============================================================

Önceki bulgu (27): aralık ↔ max|Z| Pearson r = 0.92 (normalize)

(a) Klasik çift yarık şiddet eğrisi ile Riemann Z(t)² karşılaştırma
(b) Bu korelasyon Selberg'in öngörüsünden mi türetiliyor, daha mı güçlü?
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import mpmath as mp
from scipy.stats import pearsonr

mp.mp.dps = 20

ROOT = Path(__file__).resolve().parent.parent
gamma = np.loadtxt(ROOT / "zeros_100k.txt")

# ============================================================
# (a) Klasik çift yarık: I(y) = cos²(πd·y/λL)·sinc²(πa·y/λL)
# ============================================================
def double_slit(y, d=1, L=10, lam=1, a=0.2):
    cos_t = np.cos(np.pi*d*y/(lam*L))**2
    sinc_arg = np.pi*a*y/(lam*L)
    sinc_t = np.where(np.abs(sinc_arg)<1e-9, 1.0, (np.sin(sinc_arg)/sinc_arg)**2)
    return cos_t * sinc_t

y = np.linspace(-15, 15, 2000)
I_classic = double_slit(y)

# Klasik plakada parlak şeritler (lokal maks) ve aralıkları
from scipy.signal import find_peaks
peaks_cl, _ = find_peaks(I_classic, height=0.02, distance=20)
y_peaks_cl = y[peaks_cl]
heights_cl = I_classic[peaks_cl]
intervals_cl = np.diff(y_peaks_cl)
print(f"Klasik çift yarık: {len(y_peaks_cl)} parlak şerit")
print(f"  aralık ort: {intervals_cl.mean():.3f}, std: {intervals_cl.std():.3f} (sabit beklenir!)")

# ============================================================
# Riemann Z(t)² hesapla — şiddet
# ============================================================
print("\nRiemann Z(t) hesaplanıyor (geniş pencere)...")
t_lo, t_hi = 0.5, 300
N_grid = 5000
t_grid = np.linspace(t_lo, t_hi, N_grid)
Z = np.array([float(mp.siegelz(tt)) for tt in t_grid])
I_riemann = Z**2

# ============================================================
# (b) Riemann verisi: aralık vs max|Z| (27'den)
# ============================================================
gamma_w = gamma[gamma < 1000]
print(f"\nRiemann analizi: ilk {len(gamma_w)} sıfır")
intervals_r = np.diff(gamma_w)
max_Z = np.zeros(len(intervals_r))
print("max|Z| ölçülüyor...")
for n in range(len(intervals_r)):
    t_sub = np.linspace(gamma_w[n]+0.01, gamma_w[n+1]-0.01, 25)
    Z_sub = np.array([float(mp.siegelz(tt)) for tt in t_sub])
    max_Z[n] = np.max(np.abs(Z_sub))
    if n % 200 == 0: print(f"  {n}/{len(intervals_r)}")

# Normalize aralık
t_mid = 0.5 * (gamma_w[:-1] + gamma_w[1:])
rho = np.log(t_mid/(2*np.pi))/(2*np.pi)
intervals_normalized = intervals_r * rho

r_emp, p_emp = pearsonr(intervals_normalized, max_Z)
print(f"\nBizim sonucumuz: r = {r_emp:.4f}")

# ============================================================
# SELBERG NAIVE ÖNGÖRÜ:
# Z(t) ≈ R(t) sin(θ(t)), R(t) ortalama genlik ~ √(log(t/2π))
# Eğer aralık λ ise: yarım periyot uzunlukta sin'in maksimum = R(t) (sabit)
# Yani naif: max|Z| ≈ R(t), aralığa bağlı değil → r naif ≈ 0
#
# Daha doğru: Z = lokal sin × yavaş-değişen R, max ≈ R × |sin maks|
# Aralık büyükse Z daha çok faz alır → daha çok salınım → max büyüyebilir
# Riemann-Siegel: max|Z|_n ≈ R(t) × g(spacing) burada g(λ) yarı-açık
# ============================================================

# Selberg null modeli: aralık etkisiz, sadece R(t) etkisi
# R(t) için ortalama proxy: yerel sıfır yoğunluğu × bir sabit
# R_avg(t) = sqrt(log(t/2π)) (Selberg ekstremal teoreminden)
R_local = np.sqrt(np.log(t_mid/(2*np.pi)))

r_null_R, _ = pearsonr(R_local, max_Z)
print(f"Selberg-naif null: r(R_local, max|Z|) = {r_null_R:.4f}")

# Eğer aralık katkısı varsa, residual = max_Z - α*R_local
# Bu residual aralık ile korele olmalı
slope_R, intercept_R = np.polyfit(R_local, max_Z, 1)
residual_Z = max_Z - (slope_R * R_local + intercept_R)
r_residual, _ = pearsonr(intervals_normalized, residual_Z)
print(f"R_local etkisi çıkarıldıktan sonra: r(aralık, kalan max|Z|) = {r_residual:.4f}")

# Bizim normalize edilmiş ile karşılaştır:
# Eğer r_residual hâlâ büyük (>0.5), bizim sonuç Selberg'i AŞAR
# Eğer r_residual küçükse (~0), Selberg açıklar

# ============================================================
# Plot — büyük birleşik görsel
# ============================================================
fig = plt.figure(figsize=(16, 12))
gs = fig.add_gridspec(3, 2, hspace=0.35, wspace=0.25)

# (1a) Klasik çift yarık şiddet
ax = fig.add_subplot(gs[0, 0])
ax.plot(y, I_classic, "gold", lw=2)
ax.fill_between(y, I_classic, 0, color="gold", alpha=0.3)
# Lokal şeritleri vurgula
ax.plot(y_peaks_cl, heights_cl, "ro", ms=6, label="parlak şerit")
ax.set_xlabel("plaka konumu y")
ax.set_ylabel("şiddet I(y)")
ax.set_title("KLASİK ÇİFT YARIK ŞİDDETİ\n"
             "şeritler EŞİT ARALIKLI (λL/d), yükseklik SINC² zarfı")
ax.legend(fontsize=10)
ax.grid(alpha=0.3)

# (1b) Riemann Z(t)² şiddet (aynı pencere)
ax = fig.add_subplot(gs[0, 1])
mask_show = t_grid < 100
ax.plot(t_grid[mask_show], I_riemann[mask_show], "cyan", lw=1)
ax.fill_between(t_grid[mask_show], I_riemann[mask_show], 0, color="cyan", alpha=0.3)
gamma_show = gamma_w[gamma_w<100]
ax.plot(gamma_show, np.zeros_like(gamma_show), "ro", ms=5, label=f"sıfır ({len(gamma_show)})")
ax.set_xlabel("t (plaka konumu)")
ax.set_ylabel("Z(t)²  (şiddet)")
ax.set_title("RIEMANN PLAKASI ŞİDDETİ — Z(t)²\n"
             "şeritler DEĞİŞKEN ARALIK + DEĞİŞKEN YÜKSEKLİK")
ax.legend(fontsize=10)
ax.grid(alpha=0.3)

# (2) Aralık histogramları yan yana
ax = fig.add_subplot(gs[1, 0])
ax.hist(intervals_cl, bins=15, alpha=0.6, color="gold", edgecolor="darkorange",
        label=f"klasik (ort={intervals_cl.mean():.2f}, std={intervals_cl.std():.3f})")
ax.set_xlabel("aralık (klasik birim)")
ax.set_ylabel("şerit sayısı")
ax.set_title("Klasik çift yarık: aralık DAĞILIMI dar (sabit)")
ax.legend(fontsize=10)
ax.grid(alpha=0.3)

ax = fig.add_subplot(gs[1, 1])
ax.hist(intervals_r, bins=40, alpha=0.6, color="cyan", edgecolor="navy",
        label=f"Riemann (ort={intervals_r.mean():.2f}, std={intervals_r.std():.3f})")
ax.set_xlabel("aralık")
ax.set_ylabel("şerit sayısı")
ax.set_title("Riemann: aralık DAĞILIMI geniş (logaritmik sıkışma + GUE)")
ax.legend(fontsize=10)
ax.grid(alpha=0.3)

# (3) Selberg karşılaştırma
ax = fig.add_subplot(gs[2, :])
ax.scatter(intervals_normalized, max_Z, s=15, alpha=0.5, color="steelblue",
           label=f"empirik (r={r_emp:.3f})")
# Selberg null: sadece R(t) açıklıyor — yatay
x_line = np.linspace(intervals_normalized.min(), intervals_normalized.max(), 100)
y_const = max_Z.mean() * np.ones_like(x_line)
ax.plot(x_line, y_const, "k--", lw=1, alpha=0.5, label=f"Selberg-naif null (r=0)")
# Linear fit
slope, intercept = np.polyfit(intervals_normalized, max_Z, 1)
ax.plot(x_line, slope*x_line + intercept, "r-", lw=2,
        label=f"linear fit: y = {slope:.2f}x + {intercept:.2f}")
ax.set_xlabel("normalize aralık")
ax.set_ylabel("max |Z|")
ax.set_title(f"(b) SELBERG KARŞILAŞTIRMA\n"
             f"r(empirik) = {r_emp:.3f}  ||  r(R(t) null) = {r_null_R:.3f}  ||  "
             f"R(t) çıkarıldıktan sonra r = {r_residual:.3f}\n"
             f"→ {r_residual:.3f} hâlâ büyükse: bizim sonuç R(t) etkisinin ÖTESİNDE bir şey")
ax.legend(fontsize=10)
ax.grid(alpha=0.3)

plt.suptitle("28 — (a) Çift Yarık × Riemann Paralel  &  (b) Selberg Karşılaştırma",
             fontsize=13, fontweight="bold", y=0.998)

out = Path(__file__).parent / "28_paralel_test.png"
plt.savefig(out, dpi=130, bbox_inches="tight")
print(f"\nKaydedildi: {out}")

# ============================================================
# Yorum
# ============================================================
print(f"\n=== YORUM ===")
print(f"(a) Yapısal karşılaştırma:")
print(f"  Klasik şeritler: ortalama aralık {intervals_cl.mean():.3f}, std {intervals_cl.std():.3f}")
print(f"  Riemann şeritler: ortalama {intervals_r.mean():.3f}, std {intervals_r.std():.3f}")
print(f"  → Klasikte SABİT aralık, Riemann'da DAĞINIK aralık")
print(f"")
print(f"(b) Selberg açıklıyor mu?")
print(f"  Bizim r(aralık, max|Z|) = {r_emp:.3f}")
print(f"  Selberg-naif r(R_local, max|Z|) = {r_null_R:.3f}")
print(f"  R(t) etkisi çıkarıldıktan sonra residual r = {r_residual:.3f}")
if abs(r_residual) > 0.5:
    print(f"  → r_residual = {r_residual:.3f} BÜYÜK → bizim sonuç Selberg'i AŞAN bir şey içeriyor")
else:
    print(f"  → r_residual = {r_residual:.3f} KÜÇÜK → Selberg + R(t) açıklıyor, yeni değil")
