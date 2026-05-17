"""
31 — slope² = π · log(t/2π) + b hipotezi
=============================================

Önceki test (30) → slope artıyor.
Yeni gözlem: slope² ile log(t/2π) DOĞRUSAL ilişki, eğim ≈ π!

Test:
  slope²(t) = a · log(t/2π) + b
  a ≈ π ?
  b ≈ ?
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
# Daha geniş t, daha çok pencere
# ============================================================
T_MAX = 3000.0
gamma_w = gamma[gamma < T_MAX]
print(f"İşlenecek sıfır: {len(gamma_w)}")

intervals = np.diff(gamma_w)
max_amps = np.zeros(len(intervals))
print("max|Z| hesaplanıyor (bu sefer çok sayıda nokta)...")
for n in range(len(intervals)):
    t_sub = np.linspace(gamma_w[n]+0.005, gamma_w[n+1]-0.005, 30)
    Z_sub = np.array([float(mp.siegelz(tt)) for tt in t_sub])
    max_amps[n] = np.max(np.abs(Z_sub))
    if n % 300 == 0:
        print(f"  {n}/{len(intervals)}: t≈{gamma_w[n]:.0f}")

t_mid = 0.5 * (gamma_w[:-1] + gamma_w[1:])
rho = np.log(t_mid / (2*np.pi)) / (2*np.pi)
intervals_norm = intervals * rho

# ============================================================
# Daha fazla pencere — slope² vs log(t/2π) doğrusallık testi
# ============================================================
n_windows = 15
edges = np.geomspace(20, T_MAX, n_windows+1)

results = []
print(f"\n{'pencere':>15} | {'n':>5} | {'slope':>7} | {'slope²':>7} | {'log(t/2π)':>10} | {'r':>6}")
print("-"*70)
for i in range(n_windows):
    t_lo, t_hi = edges[i], edges[i+1]
    mask = (t_mid >= t_lo) & (t_mid < t_hi)
    if mask.sum() < 15: continue
    x = intervals_norm[mask]
    y = max_amps[mask]
    slope, intercept = np.polyfit(x, y, 1)
    r, _ = pearsonr(x, y)
    t_center = np.exp(0.5 * (np.log(t_lo) + np.log(t_hi)))  # geometric center
    log_t = np.log(t_center / (2*np.pi))
    print(f"[{t_lo:>5.0f},{t_hi:>5.0f}] | {mask.sum():>5} | {slope:>7.3f} | {slope**2:>7.2f} | {log_t:>10.3f} | {r:>6.3f}")
    results.append((t_center, slope, slope**2, log_t, mask.sum()))

results = np.array(results)

# ============================================================
# slope² = a · log(t/2π) + b fit
# ============================================================
log_ts = results[:, 3]
slope_sqs = results[:, 2]

a_fit, b_fit = np.polyfit(log_ts, slope_sqs, 1)
print(f"\nFit: slope² = {a_fit:.4f} · log(t/2π) + {b_fit:.4f}")
print(f"  π = {np.pi:.4f}")
print(f"  a/π oranı: {a_fit/np.pi:.4f}")
print(f"  fark: a - π = {a_fit - np.pi:+.4f}")
print(f"  Eğer a = π ise b ≈ {(slope_sqs - np.pi*log_ts).mean():.4f}")

# a = π SABİT TUTULMUŞ fit
b_with_pi = (slope_sqs - np.pi * log_ts).mean()
b_std = (slope_sqs - np.pi * log_ts).std()
print(f"\na = π SABİT TUTULURSA: b = {b_with_pi:.4f} ± {b_std:.4f}")

# b adayları
print(f"\nb için aday kapalı formlar:")
adaylar_b = {
    "2π": 2*np.pi,
    "π+e": np.pi + np.e,
    "2π·√2/π": 2*np.sqrt(2),
    "π²/log(2π)": np.pi**2/np.log(2*np.pi),
    "3γ+π": 3*0.5772 + np.pi,
    "log(2π)·π": np.pi * np.log(2*np.pi),
    "e^(π/2)·...": np.exp(np.pi/2)/1.5,
}
for label, val in adaylar_b.items():
    print(f"  {label:>15} = {val:.4f}  (sapma: {abs(val-b_with_pi):.4f})")

# Residual analiz
residual_pi = slope_sqs - (np.pi * log_ts + b_with_pi)
print(f"\nResidual (a=π, b={b_with_pi:.3f} sabit):")
print(f"  ortalama: {residual_pi.mean():+.5f}")
print(f"  std:      {residual_pi.std():.4f}")
print(f"  max abs:  {np.max(np.abs(residual_pi)):.4f}")

# ============================================================
# Plot
# ============================================================
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# (1) slope² vs log(t/2π) — düz çizgi mi?
ax = axes[0, 0]
ax.plot(log_ts, slope_sqs, "o-", lw=1.5, ms=10, color="purple",
        label=f"empirik ({len(results)} pencere)")
xfit = np.linspace(log_ts.min(), log_ts.max(), 100)
ax.plot(xfit, a_fit * xfit + b_fit, "r--", lw=2,
        label=f"fit: y = {a_fit:.3f}·x + {b_fit:.3f}")
ax.plot(xfit, np.pi * xfit + b_with_pi, "g:", lw=2,
        label=f"a=π sabit: y = π·x + {b_with_pi:.3f}")
ax.set_xlabel("log(t/2π)")
ax.set_ylabel("slope²")
ax.set_title(f"slope² ↔ log(t/2π) → DOĞRUSAL!\nFit eğimi {a_fit:.4f} vs π = {np.pi:.4f}")
ax.legend(fontsize=10)
ax.grid(alpha=0.3)

# (2) Residual
ax = axes[0, 1]
ax.plot(log_ts, residual_pi, "o-", lw=1.2, ms=8, color="darkred")
ax.axhline(0, color="black", lw=0.5)
ax.set_xlabel("log(t/2π)")
ax.set_ylabel("residual (a=π, b sabit)")
ax.set_title(f"a=π hipotezi residuali — sıfıra ne kadar yakın?\n"
             f"std = {residual_pi.std():.4f}")
ax.grid(alpha=0.3)

# (3) Slope kendisi
ax = axes[1, 0]
ax.plot(np.exp(log_ts) * 2*np.pi, results[:, 1], "o-", lw=1.5, ms=8, color="green")
# Tahmin: √(π·log(t/2π) + b_with_pi)
t_curve = np.geomspace(20, T_MAX, 200)
slope_pred = np.sqrt(np.pi * np.log(t_curve/(2*np.pi)) + b_with_pi)
ax.plot(t_curve, slope_pred, "r--", lw=2,
        label=f"√(π·log(t/2π) + {b_with_pi:.2f})")
ax.set_xscale("log")
ax.set_xlabel("t")
ax.set_ylabel("slope")
ax.set_title("slope(t) öngörü: √(π·log(t/2π) + b)")
ax.legend(fontsize=10)
ax.grid(alpha=0.3, which="both")

# (4) a parametresi tek değil — fit kararlılığı
ax = axes[1, 1]
# Pencereler kümülatif fit
cumul_a = []
for k in range(3, len(results)+1):
    a_k, _ = np.polyfit(log_ts[:k], slope_sqs[:k], 1)
    cumul_a.append(a_k)
ax.plot(range(3, len(results)+1), cumul_a, "o-", lw=1.5, ms=8, color="orange")
ax.axhline(np.pi, color="red", lw=1.5, ls="--", label=f"π = {np.pi:.4f}")
ax.set_xlabel("kullanılan pencere sayısı")
ax.set_ylabel("fit edilmiş eğim a")
ax.set_title("a'nın kararlılığı — π'ye yakınsıyor mu?")
ax.legend(fontsize=10)
ax.grid(alpha=0.3)

plt.suptitle("31 — HİPOTEZ: slope² = π · log(t/2π) + b",
             fontsize=13, fontweight="bold", y=0.998)
plt.tight_layout()

out = Path(__file__).parent / "31_pi_formul.png"
plt.savefig(out, dpi=130, bbox_inches="tight")
print(f"\nKaydedildi: {out}")
