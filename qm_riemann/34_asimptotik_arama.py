"""
34 — Asimptotik Form Arayışı
================================

T=30000'in kayıtlı verisinden, slope²(t) için en iyi modeli ara.

Aday modeller:
  M1: a·log(t/2π) + b               (lineer — yetersiz, biliyoruz)
  M2: a·log(t/2π) + b·log log(t)    (iki ölçek)
  M3: a·log(t/2π)·(1 − c/log(t/2π)) (FHK düzeltmesi)
  M4: a·(log(t/2π))^β               (kuvvet üs)
  M5: a·log(t/2π) − b·log log(t)    (Selberg-Najnudel benzeri)
  M6: a·log(t)                      (saf log t, b=0)
  M7: a·log(t/2π) + b·log²log(t)    (FHK alt-mertebe)
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.optimize import curve_fit
from scipy.stats import pearsonr

# Veriyi yükle
data = np.load(Path("/Users/ugur/Desktop/Deney/qm_riemann/33_max_amps_T30000.npz"))
intervals = data['intervals']
max_amps = data['max_amps']
t_mid = data['t_mid']
rho = np.log(t_mid / (2*np.pi)) / (2*np.pi)
intervals_norm = intervals * rho

print(f"Veri: {len(intervals)} çift, t ∈ [{t_mid[0]:.0f}, {t_mid[-1]:.0f}]")

# 30 logaritmik pencere
n_windows = 30
edges = np.geomspace(20, t_mid[-1]+10, n_windows+1)
results = []
for i in range(n_windows):
    t_lo, t_hi = edges[i], edges[i+1]
    mask = (t_mid >= t_lo) & (t_mid < t_hi)
    if mask.sum() < 20: continue
    x = intervals_norm[mask]; y = max_amps[mask]
    slope, intercept = np.polyfit(x, y, 1)
    t_center = np.exp(0.5*(np.log(t_lo)+np.log(t_hi)))
    log_t = np.log(t_center/(2*np.pi))
    loglog_t = np.log(np.log(t_center))
    results.append((t_center, slope, slope**2, log_t, loglog_t, mask.sum()))

results = np.array(results)
ts = results[:, 0]
slopes = results[:, 1]
slope_sqs = results[:, 2]
log_ts = results[:, 3]
loglog_ts = results[:, 4]

print(f"\n{n_windows} pencere fit edildi. T aralığı: [{ts[0]:.0f}, {ts[-1]:.0f}]")

# ============================================================
# Model fitleri
# ============================================================
print(f"\n{'Model':>40} | {'parametreler':>30} | {'RMS':>7} | {'R²':>6}")
print("-"*95)

models = []

# M1: a·log + b
M1 = lambda x, a, b: a*x + b
p1, _ = curve_fit(M1, log_ts, slope_sqs)
res1 = slope_sqs - M1(log_ts, *p1)
rms1 = np.sqrt(np.mean(res1**2))
r2_1 = 1 - np.var(res1)/np.var(slope_sqs)
print(f"{'M1: a·log(t/2π) + b':>40} | a={p1[0]:.3f}, b={p1[1]:.3f}      | {rms1:>7.4f} | {r2_1:>6.4f}")
models.append(('M1', M1, p1, res1, rms1))

# M2: a·log + b·log log t
M2 = lambda x_pair, a, b: a*x_pair[0] + b*x_pair[1]
p2, _ = curve_fit(M2, (log_ts, loglog_ts), slope_sqs)
res2 = slope_sqs - M2((log_ts, loglog_ts), *p2)
rms2 = np.sqrt(np.mean(res2**2))
r2_2 = 1 - np.var(res2)/np.var(slope_sqs)
print(f"{'M2: a·log(t/2π) + b·log log(t)':>40} | a={p2[0]:.3f}, b={p2[1]:.3f}     | {rms2:>7.4f} | {r2_2:>6.4f}")
models.append(('M2', M2, p2, res2, rms2))

# M3: a·log·(1 − c/log) = a·log − a·c → a·log − constant — denenen ekstra fit
M3 = lambda x, a, c: a*x - a*c
p3, _ = curve_fit(M3, log_ts, slope_sqs)
res3 = slope_sqs - M3(log_ts, *p3)
rms3 = np.sqrt(np.mean(res3**2))
r2_3 = 1 - np.var(res3)/np.var(slope_sqs)
print(f"{'M3: a·log·(1-c/log) ~ M1':>40} | a={p3[0]:.3f}, c={p3[1]:.3f}     | {rms3:>7.4f} | {r2_3:>6.4f}")
models.append(('M3', M3, p3, res3, rms3))

# M4: a · log^β
def M4(x, a, b):
    return a * x**b
p4, _ = curve_fit(M4, log_ts, slope_sqs, p0=[3, 1])
res4 = slope_sqs - M4(log_ts, *p4)
rms4 = np.sqrt(np.mean(res4**2))
r2_4 = 1 - np.var(res4)/np.var(slope_sqs)
print(f"{'M4: a·log(t/2π)^β':>40} | a={p4[0]:.3f}, β={p4[1]:.4f}    | {rms4:>7.4f} | {r2_4:>6.4f}")
models.append(('M4', M4, p4, res4, rms4))

# M5: a·log − b·log log
M5 = lambda x_pair, a, b: a*x_pair[0] - b*x_pair[1]
p5, _ = curve_fit(M5, (log_ts, loglog_ts), slope_sqs)
res5 = slope_sqs - M5((log_ts, loglog_ts), *p5)
rms5 = np.sqrt(np.mean(res5**2))
r2_5 = 1 - np.var(res5)/np.var(slope_sqs)
print(f"{'M5: a·log(t/2π) − b·log log(t)':>40} | a={p5[0]:.3f}, b={p5[1]:.3f}     | {rms5:>7.4f} | {r2_5:>6.4f}")
models.append(('M5', M5, p5, res5, rms5))

# M6: saf log t (b=0 zorlanmış)
M6 = lambda x, a: a*x
p6, _ = curve_fit(M6, np.log(ts), slope_sqs)
res6 = slope_sqs - M6(np.log(ts), *p6)
rms6 = np.sqrt(np.mean(res6**2))
r2_6 = 1 - np.var(res6)/np.var(slope_sqs)
print(f"{'M6: a·log(t) (intercept=0)':>40} | a={p6[0]:.4f}                  | {rms6:>7.4f} | {r2_6:>6.4f}")
models.append(('M6', M6, p6, res6, rms6))

# M7: a·log + b·(log log)²
M7 = lambda x_pair, a, b: a*x_pair[0] + b*x_pair[1]**2
p7, _ = curve_fit(M7, (log_ts, loglog_ts), slope_sqs)
res7 = slope_sqs - M7((log_ts, loglog_ts), *p7)
rms7 = np.sqrt(np.mean(res7**2))
r2_7 = 1 - np.var(res7)/np.var(slope_sqs)
print(f"{'M7: a·log + b·(log log)²':>40} | a={p7[0]:.3f}, b={p7[1]:.3f}     | {rms7:>7.4f} | {r2_7:>6.4f}")
models.append(('M7', M7, p7, res7, rms7))

# ============================================================
# Bonus: Selberg-Najnudel formu
# max|Z| ~ log T - (3/4)log log T + sabit
# slope ~ (max|Z|)/normalize_gap, normalize_gap ~ 1
# Yani slope² ~ (log T - (3/4)log log T + sabit)²
# Test:
# ============================================================
M8 = lambda x_pair, c0, c1, c2: (c0 + c1*x_pair[0] + c2*x_pair[1])**2
try:
    p8, _ = curve_fit(M8, (log_ts, loglog_ts), slope_sqs, p0=[0, 1, -0.5])
    res8 = slope_sqs - M8((log_ts, loglog_ts), *p8)
    rms8 = np.sqrt(np.mean(res8**2))
    r2_8 = 1 - np.var(res8)/np.var(slope_sqs)
    print(f"{'M8: (c0 + c1·log + c2·loglog)²':>40} | c0={p8[0]:.2f},c1={p8[1]:.2f},c2={p8[2]:.2f} | {rms8:>7.4f} | {r2_8:>6.4f}")
    models.append(('M8', M8, p8, res8, rms8))
except Exception as e:
    print(f"M8 fit hatası: {e}")

# ============================================================
# En iyi modeli seç + ekstrapolasyon
# ============================================================
print(f"\n=== EN İYİ FIT ===")
best = min(models, key=lambda m: m[4])
print(f"En küçük RMS: {best[0]}")

# β = M4'ün üs değeri özellikle ilginç
print(f"\nM4 üs değeri: β = {p4[1]:.4f}")
print(f"  β = 1   → saf log")
print(f"  β = 1/2 → √log (Selberg-tipi)")
print(f"  β = 2   → log²")
print(f"  Yorum: β ≈ {p4[1]:.3f} → log'un ne kadar üzerinde/altında")

# Modelin geleceğe ekstrapolasyonu — T=10^5, 10^6'da ne öngörür
print(f"\nM4 ile ekstrapolasyon:")
for T_test in [1e4, 1e5, 1e6, 1e9]:
    log_T = np.log(T_test/(2*np.pi))
    pred = M4(log_T, *p4)
    print(f"  T={T_test:.0e}: slope² ≈ {pred:.2f}, slope ≈ {np.sqrt(pred):.3f}")

# ============================================================
# Plot
# ============================================================
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# (1) Tüm modeller üst üste
ax = axes[0, 0]
ax.plot(log_ts, slope_sqs, "ko", ms=8, label="empirik")
xfit = np.linspace(log_ts.min(), log_ts.max(), 100)
ax.plot(xfit, M1(xfit, *p1), "r--", lw=1.5, label=f"M1: lineer")
ax.plot(xfit, M4(xfit, *p4), "g-", lw=2, label=f"M4: a·log^β (β={p4[1]:.3f})")
ax.set_xlabel("log(t/2π)"); ax.set_ylabel("slope²")
ax.set_title("Empirik vs model fitleri")
ax.legend(fontsize=10); ax.grid(alpha=0.3)

# (2) M4 üs vs lineerlik
ax = axes[0, 1]
ax.plot(log_ts, slope_sqs - M1(log_ts, *p1), "ro-", label=f"M1 residual (RMS={rms1:.3f})")
ax.plot(log_ts, slope_sqs - M4(log_ts, *p4), "go-", label=f"M4 residual (RMS={rms4:.3f})")
if 'M8' in [m[0] for m in models]:
    ax.plot(log_ts, slope_sqs - M8((log_ts, loglog_ts), *p8), "bo-",
            label=f"M8 residual (RMS={rms8:.3f})")
ax.axhline(0, color="black", lw=0.5)
ax.set_xlabel("log(t/2π)"); ax.set_ylabel("residual")
ax.set_title("Hangi model residuali sıfıra en yakın?")
ax.legend(fontsize=9); ax.grid(alpha=0.3)

# (3) slope² / log(t/2π) — eğer M1 doğruysa sabit, değilse trend
ax = axes[1, 0]
ax.semilogx(ts, slope_sqs / log_ts, "o-", color="purple", lw=1.5)
ax.set_xlabel("T (log ölçek)"); ax.set_ylabel("slope² / log(t/2π)")
ax.set_title("Bu oran sabit mi? (M1 doğruysa olur)")
ax.grid(alpha=0.3, which="both")

# (4) Model tahminleri ileride
ax = axes[1, 1]
ax.plot(log_ts, slope_sqs, "ko", ms=8, label="empirik (T ≤ 30000)")
xfit_long = np.linspace(log_ts.min(), 15, 200)  # T~10^6'a kadar
ax.plot(xfit_long, M1(xfit_long, *p1), "r--", lw=1.5, alpha=0.7,
        label=f"M1 ekstrapolasyon")
ax.plot(xfit_long, M4(xfit_long, *p4), "g-", lw=1.5,
        label=f"M4 ekstrapolasyon (β={p4[1]:.2f})")
ax.set_xlabel("log(t/2π)"); ax.set_ylabel("slope²")
ax.set_title("İleriye ekstrapolasyon: T büyüdükçe modeller ne öngörür?")
ax.legend(fontsize=10); ax.grid(alpha=0.3)

plt.suptitle("34 — Asimptotik Form Arayışı: hangi model verimi en iyi açıklar?",
             fontsize=13, fontweight="bold", y=0.998)
plt.tight_layout()

out = Path(__file__).parent / "34_asimptotik_arama.png"
plt.savefig(out, dpi=130, bbox_inches="tight")
print(f"\nKaydedildi: {out}")
