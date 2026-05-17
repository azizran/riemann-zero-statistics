"""
36 — T = 100000 GECE TESTİ
============================

Hipotez doğrulama: slope² ≈ (π³/4)·log(t/2π)^(2/π)
T=30000'de β=2/π neredeyse tam yapıştı (sapma %0.3).
T=100000 ile çürür mü kalır mı?

Gece çalışsın, sabah cevap.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import mpmath as mp
from scipy.optimize import minimize_scalar, curve_fit
import time

mp.mp.dps = 10  # düşük precision (hızlı)

ROOT = Path(__file__).resolve().parent.parent
gamma = np.loadtxt(ROOT / "zeros_100k.txt")
print(f"Tüm 100k sıfır kullanılıyor. T_max = {gamma[-1]:.0f}", flush=True)

gamma_w = gamma  # hepsi
print(f"İşlenecek: {len(gamma_w)} sıfır", flush=True)
print(f"Tahmini süre: 2-4 saat", flush=True)

intervals = np.diff(gamma_w)
max_amps = np.zeros(len(intervals))

t_start = time.time()
for n in range(len(intervals)):
    a, b = gamma_w[n] + 0.005, gamma_w[n+1] - 0.005
    res = minimize_scalar(
        lambda t: -abs(float(mp.siegelz(t))),
        bounds=(a, b), method='bounded',
        options={'xatol': 0.05, 'maxiter': 10}
    )
    max_amps[n] = -res.fun
    if (n+1) % 5000 == 0:
        el = time.time() - t_start
        pct = (n+1)/len(intervals)*100
        eta = el/(n+1)*(len(intervals)-n-1)/60
        print(f"  {n+1}/{len(intervals)} ({pct:.1f}%), t≈{gamma_w[n]:.0f}, "
              f"el {el/60:.1f}dk, ETA {eta:.1f}dk", flush=True)

print(f"\nToplam: {(time.time()-t_start)/60:.1f} dakika", flush=True)
np.savez(Path(__file__).parent / "36_T100k.npz",
         intervals=intervals, max_amps=max_amps,
         t_mid=0.5*(gamma_w[:-1]+gamma_w[1:]))
print("Veri kaydedildi: 36_T100k.npz", flush=True)

t_mid = 0.5*(gamma_w[:-1]+gamma_w[1:])
rho = np.log(t_mid/(2*np.pi))/(2*np.pi)
intervals_norm = intervals * rho

# 40 pencere
n_w = 40
edges = np.geomspace(20, gamma_w[-1]+10, n_w+1)
results = []
for i in range(n_w):
    t_lo, t_hi = edges[i], edges[i+1]
    mask = (t_mid >= t_lo) & (t_mid < t_hi)
    if mask.sum() < 30: continue
    slope, _ = np.polyfit(intervals_norm[mask], max_amps[mask], 1)
    t_c = np.exp(0.5*(np.log(t_lo)+np.log(t_hi)))
    results.append((t_c, slope**2, np.log(t_c/(2*np.pi))))

results = np.array(results)
ts, slope_sqs, log_ts = results[:,0], results[:,1], results[:,2]

# β serbest fit
M = lambda x, a, b: a * x**b
p_free, _ = curve_fit(M, log_ts, slope_sqs, p0=[7.7, 0.64])
print(f"\nβ SERBEST: a = {p_free[0]:.5f}, β = {p_free[1]:.5f}", flush=True)
print(f"  2/π    = {2/np.pi:.5f}, sapma β-2/π = {p_free[1]-2/np.pi:+.5f}", flush=True)
print(f"  π³/4   = {np.pi**3/4:.5f}, sapma a-π³/4 = {p_free[0]-np.pi**3/4:+.5f}", flush=True)

# β = 2/π sabit
M_fix = lambda x, a: a * x**(2/np.pi)
p_fix, _ = curve_fit(M_fix, log_ts, slope_sqs, p0=[7.7])
print(f"\nβ = 2/π SABİT: a = {p_fix[0]:.5f}", flush=True)
print(f"  π³/4 sapma: {(p_fix[0]-np.pi**3/4)/(np.pi**3/4)*100:+.3f}%", flush=True)

# Hipotez sıfır-parametre
pred_zero = (np.pi**3/4) * log_ts**(2/np.pi)
rms_zero = np.sqrt(np.mean((slope_sqs - pred_zero)**2))
rms_free = np.sqrt(np.mean((slope_sqs - M(log_ts, *p_free))**2))
print(f"\nHipotez (sıfır parametre): slope² = (π³/4)·log^(2/π)", flush=True)
print(f"  RMS sıfır-param: {rms_zero:.4f}", flush=True)
print(f"  RMS serbest:     {rms_free:.4f}", flush=True)
print(f"  Oran: {rms_zero/rms_free:.4f}", flush=True)

# Kümülatif β trendi
print(f"\n=== KÜMÜLATİF β(T) ===", flush=True)
for T_up in [1000, 5000, 10000, 30000, 50000, 75000, gamma_w[-1]]:
    mask_T = ts < T_up
    if mask_T.sum() < 5: continue
    try:
        p, _ = curve_fit(M, log_ts[mask_T], slope_sqs[mask_T], p0=[7.7, 0.64])
        print(f"  T≤{T_up:>6.0f}: a={p[0]:.4f}, β={p[1]:.4f}, β-2/π={p[1]-2/np.pi:+.4f}", flush=True)
    except: pass

# Plot
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

ax = axes[0]
ax.plot(log_ts, slope_sqs, "ko", ms=6, label="empirik")
xfit = np.linspace(log_ts.min(), log_ts.max(), 100)
ax.plot(xfit, M(xfit, *p_free), "r--", lw=2,
        label=f"serbest: a={p_free[0]:.3f}, β={p_free[1]:.4f}")
ax.plot(xfit, (np.pi**3/4) * xfit**(2/np.pi), "g-", lw=2,
        label="hipotez: π³/4·log^(2/π)")
ax.set_xlabel("log(t/2π)"); ax.set_ylabel("slope²")
ax.set_title(f"T=100k, β={p_free[1]:.4f} vs 2/π={2/np.pi:.4f}")
ax.legend(); ax.grid(alpha=0.3)

ax = axes[1]
T_grid = np.geomspace(200, ts[-1], 60)
beta_kum = []
for T_up in T_grid:
    mask_T = ts < T_up
    if mask_T.sum() < 5: beta_kum.append(np.nan); continue
    try:
        p, _ = curve_fit(M, log_ts[mask_T], slope_sqs[mask_T], p0=[7.7, 0.64])
        beta_kum.append(p[1])
    except: beta_kum.append(np.nan)
ax.semilogx(T_grid, beta_kum, "o-", color="purple", label="β fit kümülatif")
ax.axhline(2/np.pi, color="green", lw=1.5, ls="--", label=f"2/π = {2/np.pi:.4f}")
ax.set_xlabel("T"); ax.set_ylabel("β")
ax.set_title("β(T) — 2/π'ye gerçekten yakınsıyor mu?")
ax.legend(); ax.grid(alpha=0.3, which="both")

plt.tight_layout()
plt.savefig(Path(__file__).parent / "36_T100k.png", dpi=130, bbox_inches="tight")
print(f"\nGörsel kaydedildi: 36_T100k.png", flush=True)
