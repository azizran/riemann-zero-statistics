"""
33 — T = 30000 BÜYÜK TEST
==============================

slope² ≈ 3·log(t/2π) + 2π aday formülünün gerçek mi yoksa küçük-N artefakt mı
olduğunu T büyüterek test ediyoruz.

Beklenti:
  Eğer formül GERÇEK ise → a katsayısı 3'te kalır, b 2π'de kalır
  Eğer N-artefakt ise → a azalmaya devam eder (FHK log log T ölçeğine yakınsar)

Yöntem:
  - T = 30000'e kadar tüm sıfırlar
  - Her aralık için max|Z| (5 nokta + parabolik fit)
  - 25 logaritmik pencere
  - Kümülatif fit a(T), b(T)
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import mpmath as mp
from scipy.stats import pearsonr
from scipy.optimize import minimize_scalar
import time

mp.mp.dps = 12  # daha düşük precision, daha hızlı

ROOT = Path(__file__).resolve().parent.parent
gamma = np.loadtxt(ROOT / "zeros_100k.txt")

T_MAX = 30000.0
gamma_w = gamma[gamma < T_MAX]
print(f"İşlenecek sıfır: {len(gamma_w)}", flush=True)
print(f"Tahmini süre: {len(gamma_w) * 0.03 / 60:.1f} dakika", flush=True)

intervals = np.diff(gamma_w)
max_amps = np.zeros(len(intervals))

t_start = time.time()
print("Hesap başlıyor...", flush=True)

for n in range(len(intervals)):
    a, b = gamma_w[n] + 0.005, gamma_w[n+1] - 0.005
    # minimize_scalar bounded — hızlı
    res = minimize_scalar(
        lambda t: -abs(float(mp.siegelz(t))),
        bounds=(a, b),
        method='bounded',
        options={'xatol': 0.02, 'maxiter': 15}
    )
    max_amps[n] = -res.fun

    if (n+1) % 2000 == 0:
        elapsed = time.time() - t_start
        pct = (n+1) / len(intervals) * 100
        eta = elapsed / (n+1) * (len(intervals) - n - 1) / 60
        print(f"  {n+1}/{len(intervals)} ({pct:.1f}%), t≈{gamma_w[n]:.0f}, "
              f"geçen {elapsed/60:.1f}dk, ETA {eta:.1f}dk", flush=True)

elapsed_total = time.time() - t_start
print(f"\nToplam süre: {elapsed_total/60:.1f} dakika", flush=True)

# Veriyi sakla — bir daha hesaplamayı asla istemeyiz
np.savez(Path(__file__).parent / "33_max_amps_T30000.npz",
         intervals=intervals, max_amps=max_amps, t_mid=0.5*(gamma_w[:-1]+gamma_w[1:]))
print("Veri kaydedildi: 33_max_amps_T30000.npz", flush=True)

t_mid = 0.5 * (gamma_w[:-1] + gamma_w[1:])
rho = np.log(t_mid / (2*np.pi)) / (2*np.pi)
intervals_norm = intervals * rho

# 30 pencere
n_windows = 30
edges = np.geomspace(20, T_MAX, n_windows+1)

results = []
print(f"\n{'pencere':>14} | {'n':>6} | {'slope':>7} | {'slope²':>7} | {'log(t/2π)':>10}", flush=True)
print("-"*65, flush=True)
for i in range(n_windows):
    t_lo, t_hi = edges[i], edges[i+1]
    mask = (t_mid >= t_lo) & (t_mid < t_hi)
    if mask.sum() < 20: continue
    x = intervals_norm[mask]; y = max_amps[mask]
    slope, intercept = np.polyfit(x, y, 1)
    t_center = np.exp(0.5*(np.log(t_lo)+np.log(t_hi)))
    log_t = np.log(t_center/(2*np.pi))
    results.append((t_center, slope, slope**2, log_t, mask.sum()))
    print(f"[{t_lo:>5.0f},{t_hi:>5.0f}] | {mask.sum():>6} | {slope:>7.3f} | {slope**2:>7.2f} | {log_t:>10.3f}", flush=True)

results = np.array(results)
t_centers = results[:, 0]
slopes = results[:, 1]
slope_sqs = results[:, 2]
log_ts = results[:, 3]

# Kümülatif fit
print(f"\n=== KÜMÜLATİF FIT — a YAKINSAMASI ===", flush=True)
print(f"{'T_upper':>10} | {'pencere':>8} | {'a':>8} | {'b':>10} | {'a/π':>6} | {'b/2π':>6}", flush=True)
T_uppers = [500, 1000, 2000, 4000, 8000, 12000, 18000, 25000, 30000]
for T_up in T_uppers:
    mask_T = t_centers < T_up
    if mask_T.sum() < 4: continue
    a_T, b_T = np.polyfit(log_ts[mask_T], slope_sqs[mask_T], 1)
    print(f"{T_up:>10.0f} | {mask_T.sum():>8} | {a_T:>8.4f} | {b_T:>10.4f} | "
          f"{a_T/np.pi:>6.3f} | {b_T/(2*np.pi):>6.3f}", flush=True)

# Tam fit
a_full, b_full = np.polyfit(log_ts, slope_sqs, 1)
print(f"\nTAM FIT (T={T_MAX}): a = {a_full:.4f}, b = {b_full:.4f}", flush=True)
print(f"  Hipotez: a=3, b=2π={2*np.pi:.4f}", flush=True)
print(f"  Sapma a: {(a_full-3)/3*100:+.2f}%", flush=True)
print(f"  Sapma b: {(b_full-2*np.pi)/(2*np.pi)*100:+.2f}%", flush=True)

# log log T testi — eğer a azalmaya devam ediyorsa, log log skalasına yakınsayabilir
loglog_t = np.log(np.log(t_centers))
a_vs_loglog, b_vs_loglog = np.polyfit(loglog_t, slope_sqs, 1)
print(f"\nAlternatif fit: slope² = {a_vs_loglog:.4f}·log log(t) + {b_vs_loglog:.4f}", flush=True)

# Plot
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

ax = axes[0, 0]
ax.plot(log_ts, slope_sqs, "o", ms=7, color="purple", label="empirik")
xfit = np.linspace(log_ts.min(), log_ts.max(), 100)
ax.plot(xfit, a_full*xfit + b_full, "r--", lw=2,
        label=f"fit: {a_full:.3f}·log + {b_full:.3f}")
ax.plot(xfit, 3*xfit + 2*np.pi, "g:", lw=2,
        label="hipotez: 3·log + 2π")
ax.plot(xfit, np.pi*xfit + np.pi*np.log(2*np.pi), "b:", lw=1.5, alpha=0.7,
        label="eski hipotez: π·log + π·log(2π)")
ax.set_xlabel("log(t/2π)"); ax.set_ylabel("slope²")
ax.set_title(f"T={T_MAX}, {len(results)} pencere")
ax.legend(fontsize=9); ax.grid(alpha=0.3)

ax = axes[0, 1]
T_grid = np.geomspace(50, T_MAX, 60)
a_kum = []; b_kum = []
for T_up in T_grid:
    mask_T = t_centers < T_up
    if mask_T.sum() < 4:
        a_kum.append(np.nan); b_kum.append(np.nan); continue
    a_T, b_T = np.polyfit(log_ts[mask_T], slope_sqs[mask_T], 1)
    a_kum.append(a_T); b_kum.append(b_T)
ax.semilogx(T_grid, a_kum, "o-", color="orange", label="a(T)")
ax.axhline(3, color="green", lw=1.5, ls="--", label="hipotez = 3")
ax.axhline(np.pi, color="red", lw=1, ls=":", alpha=0.7, label=f"π")
ax.set_xlabel("T üst sınırı"); ax.set_ylabel("a fit")
ax.set_title("a(T) — 3'te kalıyor mu, aşağı mı iniyor?")
ax.legend(); ax.grid(alpha=0.3, which="both")

ax = axes[1, 0]
ax.semilogx(T_grid, np.array(b_kum)/(2*np.pi), "s-", color="teal")
ax.axhline(1, color="red", lw=1.5, ls="--", label="b/2π = 1")
ax.set_xlabel("T üst sınırı"); ax.set_ylabel("b / 2π")
ax.set_title("b(T)/2π → 1 mi?")
ax.legend(); ax.grid(alpha=0.3, which="both")

# Hipotez direkt test: residual = slope² − (3·log(t/2π) + 2π)
ax = axes[1, 1]
predicted_3_2pi = 3*log_ts + 2*np.pi
residual = slope_sqs - predicted_3_2pi
ax.plot(log_ts, residual, "o-", lw=1.2, ms=7, color="darkred")
ax.axhline(0, color="black", lw=0.5)
ax.set_xlabel("log(t/2π)"); ax.set_ylabel("slope² - (3·log(t/2π) + 2π)")
ax.set_title(f"Hipotez residuali — std = {residual.std():.3f}")
ax.grid(alpha=0.3)

plt.suptitle(f"33 — T=30000 BÜYÜK TEST: slope² → 3·log(t/2π) + 2π asimptotik kalıyor mu?",
             fontsize=13, fontweight="bold", y=0.998)
plt.tight_layout()
out = Path(__file__).parent / "33_T_buyuk.png"
plt.savefig(out, dpi=130, bbox_inches="tight")
print(f"\nKaydedildi: {out}", flush=True)
