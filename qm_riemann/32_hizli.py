"""
32 — a Katsayısının T → ∞ Yakınsaması (HIZLI VERSİYON)
========================================================

Önceki yavaş: 12 nokta × her aralık. Yeni: 5 nokta + scipy.optimize.minimize_scalar
ile gerçek max bul.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import mpmath as mp
from scipy.stats import pearsonr
from scipy.optimize import minimize_scalar

mp.mp.dps = 15   # daha az hassas, daha hızlı

ROOT = Path(__file__).resolve().parent.parent
gamma = np.loadtxt(ROOT / "zeros_100k.txt")

T_MAX = 6000.0
gamma_w = gamma[gamma < T_MAX]
print(f"İşlenecek sıfır: {len(gamma_w)}", flush=True)

intervals = np.diff(gamma_w)
max_amps = np.zeros(len(intervals))

print("max|Z| hesaplanıyor (hızlı)...", flush=True)
for n in range(len(intervals)):
    a, b = gamma_w[n]+0.01, gamma_w[n+1]-0.01
    # Z burada tek bir extremum yapacak — minimize -|Z|
    res = minimize_scalar(lambda t: -abs(float(mp.siegelz(t))),
                          bounds=(a, b), method='bounded',
                          options={'xatol': 0.01})
    max_amps[n] = -res.fun
    if n % 500 == 0:
        print(f"  {n}/{len(intervals)}: t≈{gamma_w[n]:.0f}", flush=True)

t_mid = 0.5 * (gamma_w[:-1] + gamma_w[1:])
rho = np.log(t_mid / (2*np.pi)) / (2*np.pi)
intervals_norm = intervals * rho

# Çok pencere
n_windows = 20
edges = np.geomspace(20, T_MAX, n_windows+1)

results = []
for i in range(n_windows):
    t_lo, t_hi = edges[i], edges[i+1]
    mask = (t_mid >= t_lo) & (t_mid < t_hi)
    if mask.sum() < 15: continue
    x = intervals_norm[mask]; y = max_amps[mask]
    slope, intercept = np.polyfit(x, y, 1)
    r, _ = pearsonr(x, y)
    t_center = np.exp(0.5*(np.log(t_lo)+np.log(t_hi)))
    log_t = np.log(t_center/(2*np.pi))
    results.append((t_center, slope, slope**2, log_t, mask.sum()))

results = np.array(results)
t_centers = results[:, 0]
slopes = results[:, 1]
slope_sqs = results[:, 2]
log_ts = results[:, 3]

print(f"\n{'pencere t merkez':>18} | {'n':>5} | {'slope':>7} | {'slope²':>7} | {'log(t/2π)':>10}", flush=True)
for r in results:
    print(f"{r[0]:>18.0f} | {r[4]:>5.0f} | {r[1]:>7.3f} | {r[2]:>7.2f} | {r[3]:>10.3f}", flush=True)

# Kümülatif a vs T_upper
print(f"\n=== KÜMÜLATİF FIT — a YAKINSAMASI ===", flush=True)
print(f"{'T_upper':>10} | {'pencere':>8} | {'a fit':>7} | {'b fit':>8} | {'a/π':>6} | {'b/π':>6}", flush=True)
T_uppers = [200, 500, 1000, 2000, 3000, 4000, 5000, 6000]
for T_up in T_uppers:
    mask_T = t_centers < T_up
    if mask_T.sum() < 3: continue
    a_T, b_T = np.polyfit(log_ts[mask_T], slope_sqs[mask_T], 1)
    print(f"{T_up:>10.0f} | {mask_T.sum():>8} | {a_T:>7.4f} | {b_T:>8.4f} | {a_T/np.pi:>6.3f} | {b_T/np.pi:>6.3f}", flush=True)

# Tam fit
a_full, b_full = np.polyfit(log_ts, slope_sqs, 1)
print(f"\nTAM FIT (T=6000): a = {a_full:.4f}, b = {b_full:.4f}", flush=True)
print(f"  π = {np.pi:.4f}, π·log(2π) = {np.pi*np.log(2*np.pi):.4f}", flush=True)
print(f"  a sapma: {(np.pi - a_full)/np.pi*100:+.2f}%", flush=True)
print(f"  b sapma: {(np.pi*np.log(2*np.pi) - b_full)/(np.pi*np.log(2*np.pi))*100:+.2f}%", flush=True)

# Plot
fig, axes = plt.subplots(1, 2, figsize=(13, 5))

ax = axes[0]
ax.plot(log_ts, slope_sqs, "o", ms=8, color="purple", label="empirik")
xfit = np.linspace(log_ts.min(), log_ts.max(), 100)
ax.plot(xfit, a_full*xfit+b_full, "r--", lw=2, label=f"fit: {a_full:.3f}·log+{b_full:.3f}")
ax.plot(xfit, np.pi*xfit + np.pi*np.log(2*np.pi), "g:", lw=2,
        label=f"hipotez: π·log+π·log(2π)")
ax.set_xlabel("log(t/2π)"); ax.set_ylabel("slope²")
ax.set_title(f"slope² vs log(t/2π), T={T_MAX}")
ax.legend(); ax.grid(alpha=0.3)

ax = axes[1]
T_grid = np.geomspace(50, T_MAX, 50)
a_kum = []; b_kum = []
for T_up in T_grid:
    mask_T = t_centers < T_up
    if mask_T.sum() < 3:
        a_kum.append(np.nan); b_kum.append(np.nan); continue
    a_T, b_T = np.polyfit(log_ts[mask_T], slope_sqs[mask_T], 1)
    a_kum.append(a_T); b_kum.append(b_T)
ax.semilogx(T_grid, a_kum, "o-", color="orange", label="a fit")
ax.axhline(np.pi, color="red", lw=1.5, ls="--", label=f"π")
ax.set_xlabel("T üst sınırı"); ax.set_ylabel("a")
ax.set_title("a → π yakınsama")
ax.legend(); ax.grid(alpha=0.3, which="both")

plt.tight_layout()
out = Path(__file__).parent / "32_a_yakinsama.png"
plt.savefig(out, dpi=130, bbox_inches="tight")
print(f"\nKaydedildi: {out}", flush=True)
