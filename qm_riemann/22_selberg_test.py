"""
22 — Selberg Testi: Bizim Sonucumuz Yeni mi?
================================================

Dün bulduk: lokal RMS asimetri 58%, rastgele-işaret null'undan 5x güçlü.

Soru: bu **gerçekten** Z'ye özel bir özellik mi, yoksa
  (a) Selberg log-normallik teoreminin doğal sonucu mu (1946)
  (b) sadece "Z bir osilasyon, dolayısıyla yerel işaret korelasyonu var" mı?

İki katmanlı test:
  T1. log|Z(t)| dağılımı Selberg öngörüsüne (normal, std=√(½ log log T)) uyuyor mu?
  T2. Z'nin osilasyon yapısını KORUYAN daha iyi null kullansak hâlâ 5x mı?
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import mpmath as mp
from scipy.stats import norm, kstest
from scipy.signal import find_peaks

mp.mp.dps = 25

# ============================================================
# Z(t) — geniş hesap
# ============================================================
t_lo, t_hi = 10.0, 1000.0      # Selberg için büyük T lazım
N = 20000
t_grid = np.linspace(t_lo, t_hi, N)

print(f"Z(t) hesaplanıyor: t ∈ [{t_lo}, {t_hi}], {N} nokta...")
print(f"(yaklaşık 3-4 dakika)")
Z = np.zeros(N)
for i in range(N):
    Z[i] = float(mp.siegelz(t_grid[i]))
    if i % 2000 == 0:
        print(f"  i={i}/{N}, t={t_grid[i]:.1f}")
print(f"Tamamlandı. Z aralığı: [{Z.min():+.2f}, {Z.max():+.2f}]")

# ============================================================
# TEST 1: Selberg log-normallik
# ============================================================
print(f"\n=== TEST 1: SELBERG LOG-NORMALLİK ===")
# Z'nin sıfır olmayan değerleri için log|Z|
log_absZ = np.log(np.abs(Z) + 1e-12)
# Çok büyük negatif değerleri (sıfıra yakın) filtrele
mask = log_absZ > -10
log_absZ = log_absZ[mask]
print(f"Veri sayısı (sıfırlara çok yakın olanlar hariç): {len(log_absZ)}")

# Selberg öngörüsü:
# <log|ζ(1/2+it)|>_{t∈[T,2T]} ≈ 0
# Var ≈ (1/2) log log T (Selberg 1946)
t_mean = t_grid.mean()
selberg_std = np.sqrt(0.5 * np.log(np.log(t_mean)))
print(f"Selberg öngörü (T={t_mean:.0f}): std = √(½ log log T) = {selberg_std:.4f}")
print(f"Empirik mean: {log_absZ.mean():+.4f}")
print(f"Empirik std:  {log_absZ.std():.4f}")
print(f"Oran (empirik/Selberg): {log_absZ.std() / selberg_std:.3f}")

# Selberg "doğru" std ile karşılaştır
from scipy.stats import shapiro, normaltest
# Çok büyük örneklem için Shapiro yerine normaltest
stat, p_val = normaltest(log_absZ)
print(f"\nNormalité testi (D'Agostino-Pearson): stat={stat:.2f}, p={p_val:.4e}")
print(f"  p>>0 olsaydı normal denirdi; p<<0.05 ise normal değil")

# Çarpıklık (skewness) ve basıklık (kurtosis)
from scipy.stats import skew, kurtosis
print(f"Skewness: {skew(log_absZ):+.4f}  (normal=0)")
print(f"Kurtosis: {kurtosis(log_absZ):+.4f}  (normal=0)")

# ============================================================
# TEST 2: Akıllı null — osilasyon yapısını koruyan
# ============================================================
print(f"\n=== TEST 2: AKILLI NULL ===")

# Dün null'unu hatırla: Z*rastgele(±1) — bu doğal işaret korelasyonunu yok eder
# Akıllı null seçenekleri:
#   a) Random phase shift: Z(t + θ) farklı θ'larla → işaret korelasyonu korunur
#   b) Bootstrap (replikalı örnekleme)
#   c) Fit edilmiş Gaussian process simülasyonu

# (a) Phase shift null: Z'yi farklı yerlerden kes, ortala
def asymmetry(t_arr, Z_arr, n_windows):
    t_lo, t_hi = t_arr[0], t_arr[-1]
    ws = (t_hi - t_lo) / n_windows
    asyms = []
    for i in range(n_windows):
        ts, te = t_lo + i*ws, t_lo + (i+1)*ws
        mask = (t_arr >= ts) & (t_arr < te)
        Zw = Z_arr[mask]
        if len(Zw) < 5: continue
        pa = np.trapezoid(np.maximum(Zw, 0), t_arr[mask])
        na = np.trapezoid(np.maximum(-Zw, 0), t_arr[mask])
        if pa + na > 0:
            asyms.append((pa - na)/(pa + na))
    return np.array(asyms)

# Gerçek RMS asimetri (100 pencere)
n_w = 100
a_real = asymmetry(t_grid, Z, n_w)
rms_real = np.sqrt(np.mean(a_real**2)) * 100
print(f"Gerçek Z(t) (100 pencere): RMS asimetri = {rms_real:.2f}%")

# Null 1 (dün): rastgele işaret
rng = np.random.default_rng(42)
n_trial = 50
rms_null_sign = []
for _ in range(n_trial):
    Z_null = Z * rng.choice([-1, 1], size=len(Z))
    a_n = asymmetry(t_grid, Z_null, n_w)
    rms_null_sign.append(np.sqrt(np.mean(a_n**2)) * 100)
print(f"Null 1 (rastgele işaret): RMS = {np.mean(rms_null_sign):.2f} ± {np.std(rms_null_sign):.2f}%")
print(f"  Oran: gerçek/null1 = {rms_real / np.mean(rms_null_sign):.2f}")

# Null 2: Z'yi rastgele permüte et (parça parça)
# Yani yarı-pencereleri rastgele yerlerden al
rms_null_perm = []
for _ in range(n_trial):
    # 1000 birimlik dilime böl, dilimleri karıştır
    n_chunks = 20
    chunk_size = N // n_chunks
    perm = rng.permutation(n_chunks)
    Z_null = np.concatenate([Z[p*chunk_size:(p+1)*chunk_size] for p in perm])
    a_n = asymmetry(t_grid[:len(Z_null)], Z_null, n_w)
    rms_null_perm.append(np.sqrt(np.mean(a_n**2)) * 100)
print(f"Null 2 (parça permütasyonu): RMS = {np.mean(rms_null_perm):.2f} ± {np.std(rms_null_perm):.2f}%")
print(f"  Oran: gerçek/null2 = {rms_real / np.mean(rms_null_perm):.2f}")

# Null 3: AYNI Z dilimlerini RASTGELE İŞARETLİ olarak birleştir
# Bu "yerel osilasyon korunur ama bağlamı bozulur"
rms_null_chunk = []
for _ in range(n_trial):
    n_chunks = 100
    chunk_size = N // n_chunks
    Z_null = np.copy(Z)
    for k in range(n_chunks):
        s = rng.choice([-1, 1])
        Z_null[k*chunk_size:(k+1)*chunk_size] *= s
    a_n = asymmetry(t_grid, Z_null, n_w)
    rms_null_chunk.append(np.sqrt(np.mean(a_n**2)) * 100)
print(f"Null 3 (dilim işaret çevirme): RMS = {np.mean(rms_null_chunk):.2f} ± {np.std(rms_null_chunk):.2f}%")
print(f"  Oran: gerçek/null3 = {rms_real / np.mean(rms_null_chunk):.2f}")

# Null 4: Gaussian random process with same autocorrelation
# autocorr'u tahmin et, ona göre rastgele üret
from scipy.signal import correlate
Z_centered = Z - Z.mean()
acf = correlate(Z_centered, Z_centered, mode='full')[len(Z)-1:len(Z)+100]
acf /= acf[0]
# Power spectral density'den çekiş yapalım
fft_Z = np.fft.rfft(Z_centered)
psd = np.abs(fft_Z)**2
# Aynı PSD'ye sahip rastgele Gaussian üret
rms_null_gauss = []
for _ in range(n_trial):
    phases = rng.uniform(0, 2*np.pi, size=len(fft_Z))
    fft_null = np.sqrt(psd) * np.exp(1j * phases)
    Z_null = np.fft.irfft(fft_null, n=N)
    a_n = asymmetry(t_grid, Z_null, n_w)
    rms_null_gauss.append(np.sqrt(np.mean(a_n**2)) * 100)
print(f"Null 4 (aynı PSD, Gaussian): RMS = {np.mean(rms_null_gauss):.2f} ± {np.std(rms_null_gauss):.2f}%")
print(f"  Oran: gerçek/null4 = {rms_real / np.mean(rms_null_gauss):.2f}")

# ============================================================
# Plot
# ============================================================
fig = plt.figure(figsize=(14, 11))
gs = fig.add_gridspec(2, 2, hspace=0.32, wspace=0.25)

# (1) Selberg log-normallik testi
ax = fig.add_subplot(gs[0, 0])
ax.hist(log_absZ, bins=80, density=True, alpha=0.65, color="steelblue",
        edgecolor="navy", label=f"empirik log|Z(t)|")
x = np.linspace(log_absZ.min(), log_absZ.max(), 200)
ax.plot(x, norm.pdf(x, log_absZ.mean(), log_absZ.std()), "r-", lw=2,
        label=f"Normal fit (μ={log_absZ.mean():.2f}, σ={log_absZ.std():.2f})")
ax.plot(x, norm.pdf(x, 0, selberg_std), "g--", lw=2,
        label=f"Selberg öngörü (μ=0, σ={selberg_std:.2f})")
ax.set_xlabel("log |Z(t)|")
ax.set_ylabel("yoğunluk")
ax.set_title(f"Selberg log-normallik testi (T̄={t_mean:.0f})\n"
             f"normaltest p={p_val:.2e}, skewness={skew(log_absZ):+.3f}")
ax.legend(fontsize=9)
ax.grid(alpha=0.3)

# (2) Null karşılaştırma
ax = fig.add_subplot(gs[0, 1])
nulls = [("Gerçek Z", rms_real, 0),
         ("Null 1\n(rastgele işaret)", np.mean(rms_null_sign), np.std(rms_null_sign)),
         ("Null 2\n(parça permüt.)", np.mean(rms_null_perm), np.std(rms_null_perm)),
         ("Null 3\n(dilim işaret çev.)", np.mean(rms_null_chunk), np.std(rms_null_chunk)),
         ("Null 4\n(aynı PSD Gauss.)", np.mean(rms_null_gauss), np.std(rms_null_gauss))]
labels = [n[0] for n in nulls]
values = [n[1] for n in nulls]
errors = [n[2] for n in nulls]
colors = ["crimson"] + ["grey"]*4
ax.bar(range(len(nulls)), values, yerr=errors, color=colors, alpha=0.75,
       edgecolor="black", capsize=5)
ax.set_xticks(range(len(nulls)))
ax.set_xticklabels(labels, fontsize=9)
ax.set_ylabel("RMS lokal asimetri (%)")
ax.set_title("Gerçek Z vs çeşitli null'lar — hangi null en yakın?")
ax.grid(alpha=0.3, axis='y')

# (3) Z(t) örnek + Null 4 (aynı PSD)
ax = fig.add_subplot(gs[1, 0])
mask_show = (t_grid >= 100) & (t_grid <= 200)
ax.plot(t_grid[mask_show], Z[mask_show], color="navy", lw=0.8,
        label="Gerçek Z(t)")
# Bir tane Null 4 örneği
phases = rng.uniform(0, 2*np.pi, size=len(fft_Z))
fft_null = np.sqrt(psd) * np.exp(1j * phases)
Z_n4 = np.fft.irfft(fft_null, n=N)
ax.plot(t_grid[mask_show], Z_n4[mask_show] - 12, color="purple", lw=0.8, alpha=0.7,
        label="Null 4 örnek (aynı PSD)")
ax.set_xlim(100, 200)
ax.set_xlabel("t")
ax.set_ylabel("Z(t)  /  Null 4")
ax.set_title("Z(t) vs aynı PSD'li rastgele Gaussian — gözle benzer mi?")
ax.legend(fontsize=10)
ax.grid(alpha=0.3)

# (4) Özet karar
ax = fig.add_subplot(gs[1, 1])
ax.axis("off")
ratio_best = rms_real / np.mean(rms_null_gauss)
selberg_match = log_absZ.std() / selberg_std

if selberg_match > 0.95 and selberg_match < 1.05:
    selberg_verdict = "✓ Selberg öngörüsü mükemmel uyuyor"
elif selberg_match > 0.85 and selberg_match < 1.15:
    selberg_verdict = "△ Selberg yakın, küçük sapma"
else:
    selberg_verdict = "✗ Selberg'den sapma"

if ratio_best > 2.0:
    null_verdict = f"✓ Hâlâ {ratio_best:.1f}x güçlü — yapı gerçek"
elif ratio_best > 1.3:
    null_verdict = f"△ {ratio_best:.1f}x — sınırda, daha derin test lazım"
else:
    null_verdict = f"✗ Sadece {ratio_best:.1f}x — Gaussian process açıklıyor"

ax.text(0.05, 0.85, "ÖZET — Bizim sonucumuz YENİ mi?",
        fontsize=14, fontweight="bold", transform=ax.transAxes)
ax.text(0.05, 0.70, f"Selberg log-normallik:\n  {selberg_verdict}\n"
        f"  empirik σ = {log_absZ.std():.3f}\n  Selberg σ = {selberg_std:.3f}\n"
        f"  oran = {selberg_match:.3f}",
        fontsize=11, transform=ax.transAxes, family='monospace')
ax.text(0.05, 0.35, f"En akıllı null (aynı PSD):\n  {null_verdict}\n"
        f"  Gerçek RMS = {rms_real:.1f}%\n"
        f"  Null 4 RMS = {np.mean(rms_null_gauss):.1f}%\n"
        f"  oran = {ratio_best:.2f}",
        fontsize=11, transform=ax.transAxes, family='monospace')

plt.suptitle("22 — Selberg Testi: Bizim Sonucumuz Yeni mi yoksa Selberg'in Görüntüsü mü?",
             fontsize=13, fontweight="bold", y=0.998)

out = Path(__file__).parent / "22_selberg_test.png"
plt.savefig(out, dpi=130, bbox_inches="tight")
print(f"\nKaydedildi: {out}")

print(f"\n=== KARAR ===")
print(f"Selberg log-normallik: {selberg_verdict}")
print(f"Akıllı null karşı: {null_verdict}")
