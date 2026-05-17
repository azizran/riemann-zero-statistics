"""
30 — Aralık-Genlik Eğiminin Kapalı Formu
==============================================

Önceki: linear fit slope ≈ 4.41 (genel)
Soru: bu sabit mi yoksa t-bağımlı mı?

Test: t-pencerelerine böl, her birinde slope hesapla.
Eğer slope = const → kapalı form sabit ara (π√2, 2π/log, ...)
Eğer slope ∝ √log(t) → Selberg-tipi
Eğer slope ∝ log(t) → başka bir yapı
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
# Geniş aralıkta tüm sıfırlar için aralık + max|Z|
# Tek hesap, sonra t-pencerelerine böl
# ============================================================
T_MAX = 2000.0
gamma_w = gamma[gamma < T_MAX]
print(f"İşlenecek sıfır: {len(gamma_w)}, t ∈ [{gamma_w[0]:.2f}, {gamma_w[-1]:.2f}]")

intervals = np.diff(gamma_w)
max_amps = np.zeros(len(intervals))
print("max|Z| hesaplanıyor...")
for n in range(len(intervals)):
    t_sub = np.linspace(gamma_w[n]+0.01, gamma_w[n+1]-0.01, 20)
    Z_sub = np.array([float(mp.siegelz(tt)) for tt in t_sub])
    max_amps[n] = np.max(np.abs(Z_sub))
    if n % 200 == 0:
        print(f"  {n}/{len(intervals)}: t≈{gamma_w[n]:.0f}")

t_mid = 0.5 * (gamma_w[:-1] + gamma_w[1:])

# Normalize aralık (yerel yoğunluk × aralık)
rho = np.log(t_mid / (2*np.pi)) / (2*np.pi)
intervals_norm = intervals * rho

# ============================================================
# T-pencereleri içinde slope
# ============================================================
windows = [(20, 80), (80, 150), (150, 250), (250, 400),
           (400, 600), (600, 900), (900, 1300), (1300, 2000)]

print(f"\n{'t penceresi':>15} | {'n çift':>7} | {'slope':>8} | {'intercept':>9} | {'r':>7} | {'√log(t/2π)':>11}")
print("-"*75)

results = []
for t_lo, t_hi in windows:
    mask = (t_mid >= t_lo) & (t_mid < t_hi)
    if mask.sum() < 20:
        continue
    x = intervals_norm[mask]
    y = max_amps[mask]
    slope, intercept = np.polyfit(x, y, 1)
    r, _ = pearsonr(x, y)
    t_center = 0.5 * (t_lo + t_hi)
    sqrt_log = np.sqrt(np.log(t_center / (2*np.pi)))
    print(f"[{t_lo:>5}, {t_hi:>5}] | {mask.sum():>7} | {slope:>8.3f} | {intercept:>9.3f} | {r:>7.3f} | {sqrt_log:>11.3f}")
    results.append((t_center, slope, intercept, r, sqrt_log, mask.sum()))

results = np.array(results)
t_centers = results[:, 0]
slopes = results[:, 1]
intercepts = results[:, 2]
rs = results[:, 3]
sqrt_logs = results[:, 4]

# ============================================================
# Slope'un t-bağımlılığını ara
# ============================================================
print(f"\n=== SLOPE / √log(t/2π) ORANI (eğer Selberg-tipi: sabit) ===")
for i in range(len(t_centers)):
    print(f"  t={t_centers[i]:>6.0f}:  slope={slopes[i]:>6.3f},  /√log={slopes[i]/sqrt_logs[i]:>7.4f}")

print(f"\n=== SLOPE'UN KENDİSİ (eğer sabit-form: aynı) ===")
print(f"Ortalama slope: {slopes.mean():.4f}, std: {slopes.std():.4f}")
print(f"  π√2     = {np.pi*np.sqrt(2):.4f}")
print(f"  2π/√2   = {2*np.pi/np.sqrt(2):.4f}  (= π√2)")
print(f"  √(2π²)  = {np.sqrt(2*np.pi**2):.4f}")
print(f"  2π/log(2π)/sqrt(2) = {2*np.pi/np.log(2*np.pi)/np.sqrt(2):.4f}")

# ============================================================
# Slope vs √log(t/2π) — doğrusal mı?
# ============================================================
fit_slope_vs_sqrtlog = np.polyfit(sqrt_logs, slopes, 1)
print(f"\nLineer fit: slope = {fit_slope_vs_sqrtlog[0]:.4f} × √log(t/2π) + {fit_slope_vs_sqrtlog[1]:.4f}")

# Daha güçlü model: slope = c × √log(t/2π)  (intercept yok)
c_fit = np.sum(slopes * sqrt_logs) / np.sum(sqrt_logs**2)
print(f"Saf model: slope ≈ {c_fit:.4f} × √log(t/2π)")
print(f"  e^(1/π)  = {np.exp(1/np.pi):.4f}")
print(f"  √e       = {np.sqrt(np.e):.4f}")
print(f"  4/√π     = {4/np.sqrt(np.pi):.4f}")
print(f"  π/√2     = {np.pi/np.sqrt(2):.4f}")
print(f"  2 × (something)?")

# ============================================================
# Plot
# ============================================================
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# (1) Slope vs t
ax = axes[0, 0]
ax.plot(t_centers, slopes, "o-", lw=1.5, ms=8, color="purple")
ax.axhline(slopes.mean(), color="red", lw=1, ls="--",
           label=f"ortalama = {slopes.mean():.3f}")
ax.set_xlabel("t (pencere merkezi)")
ax.set_ylabel("slope (aralık → max|Z|)")
ax.set_title("Slope t ile değişiyor mu? Sabit mi yoksa monoton mu?")
ax.legend(fontsize=10)
ax.grid(alpha=0.3)

# (2) Slope vs √log(t/2π) — Selberg testi
ax = axes[0, 1]
ax.plot(sqrt_logs, slopes, "o-", lw=1.5, ms=8, color="green")
# Lineer fit çizgisi
x_fit = np.linspace(sqrt_logs.min(), sqrt_logs.max(), 100)
ax.plot(x_fit, c_fit * x_fit, "r--", lw=1.5,
        label=f"slope = {c_fit:.3f} × √log(t/2π)")
ax.plot(x_fit, fit_slope_vs_sqrtlog[0] * x_fit + fit_slope_vs_sqrtlog[1],
        "b:", lw=1.5,
        label=f"slope = {fit_slope_vs_sqrtlog[0]:.2f}×√log + {fit_slope_vs_sqrtlog[1]:.2f}")
ax.set_xlabel("√log(t/2π)")
ax.set_ylabel("slope")
ax.set_title("Slope ∝ √log(t/2π) mı? Eğer DOĞRUSAL → Selberg formu")
ax.legend(fontsize=9)
ax.grid(alpha=0.3)

# (3) Tüm verinin saçılımı renklere göre (t penceresi)
ax = axes[1, 0]
cmap = plt.cm.viridis
for i, (t_lo, t_hi) in enumerate(windows):
    mask = (t_mid >= t_lo) & (t_mid < t_hi)
    if mask.sum() < 20: continue
    color = cmap(i / len(windows))
    ax.scatter(intervals_norm[mask], max_amps[mask], s=8, alpha=0.5,
               color=color, label=f"t∈[{t_lo},{t_hi}]")
ax.set_xlabel("normalize aralık")
ax.set_ylabel("max |Z|")
ax.set_title("Tüm veri renklere göre t penceresi")
ax.legend(fontsize=7, loc="upper left")
ax.grid(alpha=0.3)

# (4) intercept vs t — sabit mi, artıyor mu?
ax = axes[1, 1]
ax.plot(t_centers, intercepts, "o-", lw=1.5, ms=8, color="orange")
ax.set_xlabel("t (pencere merkezi)")
ax.set_ylabel("intercept (aralık=0 limit)")
ax.set_title("Intercept t ile değişiyor mu?")
ax.grid(alpha=0.3)

plt.suptitle("30 — Aralık-Genlik Eğiminin Şekli: sabit mi, Selberg-tipi mi?",
             fontsize=13, fontweight="bold", y=0.998)
plt.tight_layout()

out = Path(__file__).parent / "30_slope_formul.png"
plt.savefig(out, dpi=130, bbox_inches="tight")
print(f"\nKaydedildi: {out}")
