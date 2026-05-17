"""
29 — Aralık-Genlik Korelasyonu: Gaussian-PSD Null Testi
=============================================================

Önceki testte (28) Selberg-naif null çok kaba idi. Şimdi doğru null:
  Z(t) için PSD hesapla → AYNI PSD'li Gaussian rastgele süreç üret
  Onun sıfır geçişlerini bul → aralık vs max|Z_null| hesapla
  Bizim r=0.92 ile karşılaştır

Eğer r_null ≈ 0.92: bizim sonuç Gaussian-süreç doğası (yeni değil)
Eğer r_null << 0.92: bizim sonuç ζ-özgün (gerçek pattern)
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import mpmath as mp
from scipy.stats import pearsonr
from scipy.signal import find_peaks

mp.mp.dps = 20

# ============================================================
# Z(t) hesapla — yüksek çözünürlük
# ============================================================
t_lo, t_hi = 50.0, 500.0
N_grid = 20000
t_grid = np.linspace(t_lo, t_hi, N_grid)
dt = t_grid[1] - t_grid[0]

print(f"Z(t) hesaplanıyor (N={N_grid}, t ∈ [{t_lo}, {t_hi}])...")
Z = np.zeros(N_grid)
for i in range(N_grid):
    Z[i] = float(mp.siegelz(t_grid[i]))
    if i % 2500 == 0:
        print(f"  {i}/{N_grid}")
print(f"Z aralığı: [{Z.min():+.2f}, {Z.max():+.2f}]")

# ============================================================
# ZERO CROSSING / GAP-AMPLITUDE FONKSİYONLARI
# ============================================================
def find_zero_crossings(t_arr, Z_arr):
    """Lineer interpolasyon ile tam sıfır geçişi pozisyonları"""
    sign_chg = np.where(np.diff(np.sign(Z_arr)) != 0)[0]
    zeros = []
    for i in sign_chg:
        # Linear interp
        t0, t1 = t_arr[i], t_arr[i+1]
        z0, z1 = Z_arr[i], Z_arr[i+1]
        if z1 != z0:
            tz = t0 - z0 * (t1 - t0) / (z1 - z0)
            zeros.append(tz)
    return np.array(zeros)

def gap_amplitude_corr(t_arr, Z_arr):
    """Sıfır geçişleri arası aralık ile aralıktaki max|Z| korelasyonu"""
    zeros = find_zero_crossings(t_arr, Z_arr)
    if len(zeros) < 5:
        return np.nan, len(zeros)
    gaps = np.diff(zeros)
    max_amps = np.zeros(len(gaps))
    for n in range(len(gaps)):
        mask = (t_arr >= zeros[n]) & (t_arr <= zeros[n+1])
        if mask.sum() < 2:
            max_amps[n] = 0
        else:
            max_amps[n] = np.max(np.abs(Z_arr[mask]))
    r, p = pearsonr(gaps, max_amps)
    return r, len(zeros)

# ============================================================
# Gerçek Z(t) için korelasyon
# ============================================================
r_real, n_real = gap_amplitude_corr(t_grid, Z)
print(f"\n=== GERÇEK ζ ===")
print(f"r(aralık, max|Z|) = {r_real:.4f},  {n_real} sıfır geçişi")

# ============================================================
# Gaussian-PSD null
# ============================================================
print(f"\n=== GAUSSIAN-PSD NULL ÜRETİLİYOR ===")
Z_centered = Z - Z.mean()
fft_Z = np.fft.rfft(Z_centered)
psd = np.abs(fft_Z)**2
freqs = np.fft.rfftfreq(N_grid, d=dt)

# Çok küçük PSD'leri (gürültü) filtrele
# Çok büyük frekansları (Nyquist civarı) kes
psd_smooth = psd.copy()

rng = np.random.default_rng(11)
N_realizations = 30
r_nulls = []

for trial in range(N_realizations):
    # Aynı PSD, rastgele faz
    phases = rng.uniform(0, 2*np.pi, size=len(psd))
    fft_null = np.sqrt(psd) * np.exp(1j * phases)
    fft_null[0] = 0   # DC = 0
    Z_null = np.fft.irfft(fft_null, n=N_grid)
    # Aynı norm
    Z_null = Z_null * (np.std(Z_centered) / np.std(Z_null))

    r_n, n_n = gap_amplitude_corr(t_grid, Z_null)
    r_nulls.append(r_n)
    if trial % 5 == 0:
        print(f"  trial {trial}: r_null = {r_n:.4f},  {n_n} geçiş")

r_nulls = np.array(r_nulls)
print(f"\nGaussian null: ortalama r = {r_nulls.mean():.4f} ± {r_nulls.std():.4f}")
print(f"Min: {r_nulls.min():.4f}, Max: {r_nulls.max():.4f}")

# ============================================================
# Z-skoru
# ============================================================
z_score = (r_real - r_nulls.mean()) / r_nulls.std()
print(f"\nBizim r = {r_real:.4f}")
print(f"Gaussian null = {r_nulls.mean():.4f} ± {r_nulls.std():.4f}")
print(f"Z-skoru: {z_score:+.2f}σ")
print(f"  (>2σ: ζ-özgün; ~0σ: Gaussian açıklar)")

# ============================================================
# Karşılaştırma örnek pencerede
# ============================================================
# Bir Gaussian null örneği daha üret görsel için
phases = rng.uniform(0, 2*np.pi, size=len(psd))
fft_null = np.sqrt(psd) * np.exp(1j * phases)
fft_null[0] = 0
Z_null_example = np.fft.irfft(fft_null, n=N_grid)
Z_null_example *= np.std(Z_centered) / np.std(Z_null_example)

# Her ikisi için saçılım
zeros_real = find_zero_crossings(t_grid, Z)
gaps_real = np.diff(zeros_real)
max_real = np.zeros(len(gaps_real))
for n in range(len(gaps_real)):
    mask = (t_grid >= zeros_real[n]) & (t_grid <= zeros_real[n+1])
    if mask.sum() > 0:
        max_real[n] = np.max(np.abs(Z[mask]))

zeros_null = find_zero_crossings(t_grid, Z_null_example)
gaps_null = np.diff(zeros_null)
max_null = np.zeros(len(gaps_null))
for n in range(len(gaps_null)):
    mask = (t_grid >= zeros_null[n]) & (t_grid <= zeros_null[n+1])
    if mask.sum() > 0:
        max_null[n] = np.max(np.abs(Z_null_example[mask]))

# ============================================================
# Plot
# ============================================================
fig = plt.figure(figsize=(14, 11))
gs = fig.add_gridspec(3, 2, hspace=0.35, wspace=0.25)

# (1) Z(t) örnek pencere
ax = fig.add_subplot(gs[0, 0])
mask_show = (t_grid >= 50) & (t_grid <= 100)
ax.plot(t_grid[mask_show], Z[mask_show], color="navy", lw=0.9)
ax.fill_between(t_grid[mask_show], Z[mask_show], 0, where=(Z[mask_show]>0), alpha=0.2, color="blue")
ax.fill_between(t_grid[mask_show], Z[mask_show], 0, where=(Z[mask_show]<0), alpha=0.2, color="red")
ax.axhline(0, color="black", lw=0.5)
ax.set_xlim(50, 100)
ax.set_title(f"GERÇEK Z(t) — ζ(½+it)\nr(aralık, max|Z|) = {r_real:.4f}")
ax.grid(alpha=0.3)
ax.set_xlabel("t"); ax.set_ylabel("Z(t)")

# (2) Z_null örnek pencere
ax = fig.add_subplot(gs[0, 1])
ax.plot(t_grid[mask_show], Z_null_example[mask_show], color="purple", lw=0.9)
ax.fill_between(t_grid[mask_show], Z_null_example[mask_show], 0,
                where=(Z_null_example[mask_show]>0), alpha=0.2, color="purple")
ax.fill_between(t_grid[mask_show], Z_null_example[mask_show], 0,
                where=(Z_null_example[mask_show]<0), alpha=0.2, color="orange")
ax.axhline(0, color="black", lw=0.5)
ax.set_xlim(50, 100)
ax.set_title(f"GAUSSIAN NULL (aynı PSD) — ζ değil\nr(aralık, max|Z|) = {r_nulls[-1]:.4f}")
ax.grid(alpha=0.3)
ax.set_xlabel("t"); ax.set_ylabel("Z_null(t)")

# (3) Gerçek saçılım
ax = fig.add_subplot(gs[1, 0])
ax.scatter(gaps_real, max_real, s=10, alpha=0.5, color="steelblue")
ax.set_xlabel("aralık")
ax.set_ylabel("max |Z|")
ax.set_title(f"Gerçek Z: aralık vs max|Z| (r={r_real:.4f})")
ax.grid(alpha=0.3)

# (4) Null saçılım
ax = fig.add_subplot(gs[1, 1])
ax.scatter(gaps_null, max_null, s=10, alpha=0.5, color="purple")
ax.set_xlabel("aralık")
ax.set_ylabel("max |Z_null|")
ax.set_title(f"Gaussian null: aralık vs max|Z| (r={r_nulls[-1]:.4f})")
ax.grid(alpha=0.3)

# (5) Tüm null r'lerin dağılımı vs gerçek r
ax = fig.add_subplot(gs[2, :])
ax.hist(r_nulls, bins=15, alpha=0.7, color="purple", edgecolor="darkmagenta",
        label=f"Gaussian null ({len(r_nulls)} trial): {r_nulls.mean():.3f} ± {r_nulls.std():.3f}")
ax.axvline(r_real, color="red", lw=2.5,
           label=f"GERÇEK ζ: r = {r_real:.3f}")
ax.set_xlabel("r(aralık, max|Z|)")
ax.set_ylabel("sayı")
ax.set_title(f"Gaussian null'lar vs Gerçek ζ — Z-skoru = {z_score:+.2f}σ\n"
             "Eğer kırmızı çizgi mor histogram içindeyse: Gaussian açıklar; "
             "uzaktaysa: ζ-özgün")
ax.legend(fontsize=10)
ax.grid(alpha=0.3)

plt.suptitle("29 — DOĞRU NULL TESTİ: aralık-genlik korelasyonu ζ-özgün mü?",
             fontsize=13, fontweight="bold", y=0.998)

out = Path(__file__).parent / "29_gauss_null.png"
plt.savefig(out, dpi=130, bbox_inches="tight")
print(f"\nKaydedildi: {out}")

print(f"\n=== KARAR ===")
if z_score > 2:
    print(f"✓ Bizim r = {r_real:.4f}, Gaussian null = {r_nulls.mean():.4f}")
    print(f"  Z-skoru +{z_score:.2f}σ → BİZİM SONUÇ ζ-ÖZGÜN")
elif z_score < -2:
    print(f"✓ Gaussian null daha güçlü, gerçek r daha zayıf")
    print(f"  Z-skoru {z_score:.2f}σ → BİZİM SONUÇ Gaussian'dan ZAYIF")
else:
    print(f"△ Z-skoru {z_score:+.2f}σ — sınırda, Gaussian açıklayabilir")
    print(f"  Gaussian r = {r_nulls.mean():.4f}, Gerçek r = {r_real:.4f}")
