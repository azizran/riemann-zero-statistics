"""
32 — a Katsayısının T → ∞ Yakınsaması
=========================================

Önceki test (31): slope² ≈ 3.03·log(t/2π) + 5.75
  - b ≈ π·log(2π) = 5.774 (sapma %0.4) — çok güzel
  - a ≈ 3.03 vs π = 3.14 (sapma %3.6) — yakın ama tam değil

Soru: Daha büyük T (T=10000) için a katsayısı π'ye yakınsıyor mu?

Eğer evet → slope² = π · log(t) kapalı formu doğrulanır
Eğer hayır → a farklı bir sabit
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
# Geniş T
# ============================================================
T_MAX = 10000.0
gamma_w = gamma[gamma < T_MAX]
print(f"İşlenecek sıfır: {len(gamma_w)}, t ∈ [{gamma_w[0]:.2f}, {gamma_w[-1]:.2f}]")

intervals = np.diff(gamma_w)
max_amps = np.zeros(len(intervals))

print("max|Z| hesaplanıyor (10000'lik veri)...")
for n in range(len(intervals)):
    t_sub = np.linspace(gamma_w[n]+0.005, gamma_w[n+1]-0.005, 12)
    Z_sub = np.array([float(mp.siegelz(tt)) for tt in t_sub])
    max_amps[n] = np.max(np.abs(Z_sub))
    if n % 1000 == 0:
        print(f"  {n}/{len(intervals)}: t≈{gamma_w[n]:.0f}")

t_mid = 0.5 * (gamma_w[:-1] + gamma_w[1:])
rho = np.log(t_mid / (2*np.pi)) / (2*np.pi)
intervals_norm = intervals * rho

# ============================================================
# Çok sayıda pencere — a vs t kararlılık testi
# ============================================================
n_windows = 25
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
    t_center = np.exp(0.5 * (np.log(t_lo) + np.log(t_hi)))
    log_t = np.log(t_center / (2*np.pi))
    print(f"[{t_lo:>5.0f},{t_hi:>5.0f}] | {mask.sum():>5} | {slope:>7.3f} | {slope**2:>7.2f} | {log_t:>10.3f} | {r:>6.3f}")
    results.append((t_center, slope, slope**2, log_t, mask.sum()))

results = np.array(results)
t_centers = results[:, 0]
slopes = results[:, 1]
slope_sqs = results[:, 2]
log_ts = results[:, 3]

# ============================================================
# Kümülatif fit — a'nın T büyüdükçe yakınsama
# ============================================================
print(f"\n=== a'nın T → ∞ YAKINSAMASI ===")
print(f"{'T üst sınır':>12} | {'kullanılan pencere':>18} | {'a fit':>8} | {'b fit':>10} | {'b/π':>6}")
print("-"*65)

T_uppers = [200, 500, 1000, 2000, 4000, 6000, 8000, 10000]
for T_up in T_uppers:
    mask_T = t_centers < T_up
    if mask_T.sum() < 3: continue
    a_T, b_T = np.polyfit(log_ts[mask_T], slope_sqs[mask_T], 1)
    print(f"{T_up:>12.0f} | {mask_T.sum():>18} | {a_T:>8.4f} | {b_T:>10.4f} | {b_T/np.pi:>6.3f}")

# ============================================================
# Yakınsama plot
# ============================================================
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# (1) slope² vs log(t/2π) — tüm 25 pencere
ax = axes[0, 0]
ax.plot(log_ts, slope_sqs, "o", lw=1.2, ms=8, color="purple", label="empirik")
# fit
xfit = np.linspace(log_ts.min(), log_ts.max(), 100)
a_full, b_full = np.polyfit(log_ts, slope_sqs, 1)
ax.plot(xfit, a_full*xfit + b_full, "r--", lw=2,
        label=f"fit: {a_full:.3f}·log + {b_full:.3f}")
ax.plot(xfit, np.pi*xfit + np.pi*np.log(2*np.pi), "g:", lw=2,
        label=f"hipotez: π·log + π·log(2π)")
ax.set_xlabel("log(t/2π)")
ax.set_ylabel("slope²")
ax.set_title(f"slope² vs log(t/2π) — T={T_MAX}, {n_windows} pencere\n"
             f"Fit a={a_full:.4f} vs π={np.pi:.4f}")
ax.legend(fontsize=10)
ax.grid(alpha=0.3)

# (2) Kümülatif a(T) — π'ye yakınsıyor mu?
ax = axes[0, 1]
T_grid = np.geomspace(50, T_MAX, 50)
a_kum = []
for T_up in T_grid:
    mask_T = t_centers < T_up
    if mask_T.sum() < 3:
        a_kum.append(np.nan)
        continue
    a_T, _ = np.polyfit(log_ts[mask_T], slope_sqs[mask_T], 1)
    a_kum.append(a_T)
ax.semilogx(T_grid, a_kum, "o-", lw=1.5, ms=6, color="orange")
ax.axhline(np.pi, color="red", lw=1.5, ls="--", label=f"π = {np.pi:.4f}")
ax.set_xlabel("T üst sınırı")
ax.set_ylabel("kümülatif a fit")
ax.set_title("a → π yaklaşıyor mu?")
ax.legend(fontsize=10)
ax.grid(alpha=0.3, which="both")

# (3) b/π vs T — log(2π) sabit mi?
ax = axes[1, 0]
b_kum = []
for T_up in T_grid:
    mask_T = t_centers < T_up
    if mask_T.sum() < 3:
        b_kum.append(np.nan)
        continue
    _, b_T = np.polyfit(log_ts[mask_T], slope_sqs[mask_T], 1)
    b_kum.append(b_T)
ax.semilogx(T_grid, np.array(b_kum)/np.pi, "o-", lw=1.5, ms=6, color="green")
ax.axhline(np.log(2*np.pi), color="red", lw=1.5, ls="--",
           label=f"log(2π) = {np.log(2*np.pi):.4f}")
ax.set_xlabel("T üst sınırı")
ax.set_ylabel("b/π")
ax.set_title("b/π → log(2π) mı?")
ax.legend(fontsize=10)
ax.grid(alpha=0.3, which="both")

# (4) Hipotez: slope² = π·log(t) — direkt test
ax = axes[1, 1]
ax.plot(log_ts, slope_sqs / (np.pi * (log_ts + np.log(2*np.pi))), "o-",
        lw=1.5, ms=8, color="darkviolet")
ax.axhline(1, color="red", lw=1.5, ls="--", label="hipotez = 1")
ax.set_xlabel("log(t/2π)")
ax.set_ylabel("slope² / [π·log(t)]")
ax.set_title("Hipotez direkt test: slope² / π·log(t) → 1 olmalı")
ax.legend(fontsize=10)
ax.grid(alpha=0.3)

plt.suptitle(f"32 — a YAKINSAMASI: slope² = π·log(t) hipotezi T={T_MAX}'a kadar",
             fontsize=13, fontweight="bold", y=0.998)
plt.tight_layout()

out = Path(__file__).parent / "32_a_yakinsama.png"
plt.savefig(out, dpi=130, bbox_inches="tight")
print(f"\nKaydedildi: {out}")

print(f"\n=== KARAR ===")
print(f"En son (T=10000) fit: a = {a_full:.4f}, b = {b_full:.4f}")
print(f"Hipotez:                a = π = {np.pi:.4f}, b = π·log(2π) = {np.pi*np.log(2*np.pi):.4f}")
print(f"a sapma: {(np.pi - a_full)/np.pi*100:+.2f}%")
print(f"b sapma: {(np.pi*np.log(2*np.pi) - b_full)/(np.pi*np.log(2*np.pi))*100:+.2f}%")
