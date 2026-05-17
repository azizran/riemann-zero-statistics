"""
35 — β = 2/π Hipotez Testi
============================

slope² = a · log(t/2π)^(2/π) hipotezi.
β'yı 2/π = 0.6366'ya SABİT tutarak fit kalitesi düşüyor mu bak.
Eğer düşmüyorsa: hipotez güçlü. Düşüyorsa: koincidans.

Ek: T-pencereleri içinde β fit değeri sabit mi yoksa T ile değişiyor mu?
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.optimize import curve_fit

# Veriyi yükle
data = np.load("/Users/ugur/Desktop/Deney/qm_riemann/33_max_amps_T30000.npz")
intervals = data['intervals']
max_amps = data['max_amps']
t_mid = data['t_mid']
rho = np.log(t_mid / (2*np.pi)) / (2*np.pi)
intervals_norm = intervals * rho

# 30 pencere
n_windows = 30
edges = np.geomspace(20, t_mid[-1]+10, n_windows+1)
results = []
for i in range(n_windows):
    t_lo, t_hi = edges[i], edges[i+1]
    mask = (t_mid >= t_lo) & (t_mid < t_hi)
    if mask.sum() < 20: continue
    x = intervals_norm[mask]; y = max_amps[mask]
    slope, _ = np.polyfit(x, y, 1)
    t_center = np.exp(0.5*(np.log(t_lo)+np.log(t_hi)))
    log_t = np.log(t_center/(2*np.pi))
    results.append((t_center, slope, slope**2, log_t, mask.sum()))

results = np.array(results)
ts = results[:, 0]
slopes = results[:, 1]
slope_sqs = results[:, 2]
log_ts = results[:, 3]

# ============================================================
# Karşılaştır: β serbest vs β=2/π sabit
# ============================================================
print(f"=== β SERBEST vs β=2/π SABİT FITS ===\n")

# β serbest
M_free = lambda x, a, b: a * x**b
p_free, _ = curve_fit(M_free, log_ts, slope_sqs, p0=[7, 0.6])
pred_free = M_free(log_ts, *p_free)
rms_free = np.sqrt(np.mean((slope_sqs - pred_free)**2))
print(f"β serbest:      a = {p_free[0]:.4f}, β = {p_free[1]:.4f}")
print(f"                RMS = {rms_free:.4f}")

# β = 2/π sabit
beta_fixed = 2/np.pi
print(f"\nβ = 2/π = {beta_fixed:.6f} (sabit)")
M_fixed = lambda x, a: a * x**beta_fixed
p_fixed, _ = curve_fit(M_fixed, log_ts, slope_sqs, p0=[7])
pred_fixed = M_fixed(log_ts, *p_fixed)
rms_fixed = np.sqrt(np.mean((slope_sqs - pred_fixed)**2))
print(f"                a = {p_fixed[0]:.4f}")
print(f"                RMS = {rms_fixed:.4f}")
print(f"  RMS oranı (sabit/serbest): {rms_fixed/rms_free:.4f}")
print(f"  (1.00'a yakınsa: sabit β kalitesi düşürmüyor → hipotez güçlü)")

# Aday a değeri için kapalı formlar
print(f"\na = {p_fixed[0]:.4f}'a aday kapalı formlar:")
adaylar = {
    "2π": 2*np.pi,
    "e²": np.e**2,
    "π²/log(2π)": np.pi**2 / np.log(2*np.pi),
    "8/(√π)": 8/np.sqrt(np.pi),
    "4·log(2π)": 4*np.log(2*np.pi),
    "π·e": np.pi*np.e,
    "5·log(2π)": 5*np.log(2*np.pi),
}
for label, val in adaylar.items():
    print(f"  {label:>20} = {val:.4f}  (sapma: {abs(val-p_fixed[0])/p_fixed[0]*100:+.1f}%)")

# ============================================================
# T-pencereleri içinde β'nın sabit olup olmadığı testi
# ============================================================
print(f"\n=== T-pencerelerinde β nasıl davranıyor? ===\n")
T_uppers = [500, 1000, 2000, 5000, 10000, 20000, 30000]
for T_up in T_uppers:
    mask_T = ts < T_up
    if mask_T.sum() < 5: continue
    p, _ = curve_fit(M_free, log_ts[mask_T], slope_sqs[mask_T], p0=[7, 0.6])
    print(f"  T ≤ {T_up:>6.0f}: a = {p[0]:>7.4f}, β = {p[1]:>7.4f}, β-2/π = {p[1]-2/np.pi:+.4f}")

# ============================================================
# Daha hassas: lokal β (kayan pencere)
# ============================================================
print(f"\n=== KAYAN PENCERE (her 5 pencereyi fit et) β trendi ===")
print(f"{'T merkez':>10} {'β lokal':>10}")
for i in range(0, len(ts)-4, 3):
    sub_log = log_ts[i:i+5]
    sub_sq = slope_sqs[i:i+5]
    try:
        p, _ = curve_fit(M_free, sub_log, sub_sq, p0=[7, 0.6])
        t_c = np.exp(np.mean(np.log(ts[i:i+5])))
        print(f"  {t_c:>10.0f} {p[1]:>10.4f}")
    except: pass

# ============================================================
# Plot
# ============================================================
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# (1) Empirik vs iki model
ax = axes[0, 0]
ax.plot(log_ts, slope_sqs, "ko", ms=8, label="empirik")
xfit = np.linspace(log_ts.min(), log_ts.max()+1, 100)
ax.plot(xfit, M_free(xfit, *p_free), "r--", lw=2,
        label=f"β serbest = {p_free[1]:.4f} (RMS={rms_free:.3f})")
ax.plot(xfit, M_fixed(xfit, *p_fixed), "g-", lw=2,
        label=f"β=2/π sabit  (RMS={rms_fixed:.3f})")
ax.set_xlabel("log(t/2π)"); ax.set_ylabel("slope²")
ax.set_title(f"β serbest vs β=2/π=0.6366 sabit\nRMS oranı: {rms_fixed/rms_free:.3f}")
ax.legend(fontsize=10); ax.grid(alpha=0.3)

# (2) Residuallar
ax = axes[0, 1]
ax.plot(log_ts, slope_sqs - M_free(log_ts, *p_free), "ro-", ms=6,
        label=f"β serbest residual (RMS={rms_free:.3f})")
ax.plot(log_ts, slope_sqs - M_fixed(log_ts, *p_fixed), "go-", ms=6,
        label=f"β=2/π residual (RMS={rms_fixed:.3f})")
ax.axhline(0, color="black", lw=0.5)
ax.set_xlabel("log(t/2π)"); ax.set_ylabel("residual")
ax.set_title("Residual karşılaştırma")
ax.legend(fontsize=10); ax.grid(alpha=0.3)

# (3) β kümülatif fit T ile
ax = axes[1, 0]
T_grid = np.geomspace(200, ts[-1], 30)
beta_kum = []
a_kum = []
for T_up in T_grid:
    mask_T = ts < T_up
    if mask_T.sum() < 5:
        beta_kum.append(np.nan); a_kum.append(np.nan); continue
    try:
        p, _ = curve_fit(M_free, log_ts[mask_T], slope_sqs[mask_T], p0=[7, 0.6])
        beta_kum.append(p[1]); a_kum.append(p[0])
    except:
        beta_kum.append(np.nan); a_kum.append(np.nan)
ax.semilogx(T_grid, beta_kum, "o-", color="purple", label="β fit (kümülatif)")
ax.axhline(2/np.pi, color="green", lw=1.5, ls="--", label=f"2/π = {2/np.pi:.4f}")
ax.axhline(0.5, color="red", lw=1, ls=":", alpha=0.6, label="√log")
ax.axhline(1.0, color="blue", lw=1, ls=":", alpha=0.6, label="saf log")
ax.set_xlabel("T üst sınırı"); ax.set_ylabel("β")
ax.set_title("β(T) — 2/π'ye yakınsıyor mu?")
ax.legend(fontsize=9); ax.grid(alpha=0.3, which="both")

# (4) Ekstrapolasyon
ax = axes[1, 1]
ax.plot(log_ts, slope_sqs, "ko", ms=8, label="empirik (T ≤ 30000)")
xfit_long = np.linspace(log_ts.min(), 18, 200)
ax.plot(xfit_long, M_free(xfit_long, *p_free), "r--", lw=1.5, alpha=0.7,
        label=f"β serbest ekstrapolasyon")
ax.plot(xfit_long, M_fixed(xfit_long, *p_fixed), "g-", lw=1.5,
        label=f"β=2/π ekstrapolasyon")
ax.set_xlabel("log(t/2π)"); ax.set_ylabel("slope²")
ax.set_title("İleri ekstrapolasyon (T → 10⁸'e kadar)")
ax.legend(fontsize=10); ax.grid(alpha=0.3)

plt.suptitle("35 — β = 2/π Hipotezi: gerçek mi koincidans mi?",
             fontsize=13, fontweight="bold", y=0.998)
plt.tight_layout()

out = Path(__file__).parent / "35_iki_pi_hipotez.png"
plt.savefig(out, dpi=130, bbox_inches="tight")
print(f"\nKaydedildi: {out}")
