"""
14 — Artığın Autocorrelation Testi
===================================

Soru: 13_gonek_test'teki artık (empirik δ − explicit formula öngörüsü)
ardışık değerlerinde GUE level repulsion imzası taşıyor mu?

Test:
  - C(k) = ⟨artık(n) · artık(n+k)⟩ / ⟨artık²⟩
  - Eğer rastgele (sayısal artık): C(k) ≈ 0 her k > 0 için
  - Eğer GUE level repulsion: C(1) < 0, sonra salınım
  - Baseline: artığı KARIŞTIR (shuffle), autocorrelation hesapla → null hypothesis

GUE öngörüsü: ardışık seviyeler birbirini iter → komşu artıkları ters yönde.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.fft import rfft, rfftfreq

ROOT = Path(__file__).resolve().parent.parent
gamma = np.loadtxt(ROOT / "zeros_100k.txt")

# ============================================================
# 13_gonek_test'teki HASSAS öngörüyü yeniden üret
# ============================================================
def N_smooth_hassas(t):
    base = (t/(2*np.pi)) * np.log(t/(2*np.pi*np.e)) + 7/8
    corr1 = 1.0 / (48 * np.pi * t)
    corr2 = 7.0 / (5760 * np.pi * t**3)
    return base + corr1 + corr2

def t_of_n_hassas(n_arr):
    target = n_arr - 0.5
    t = np.maximum(2*np.pi * n_arr / np.maximum(np.log(n_arr+2), 1), 1.0)
    for _ in range(80):
        F = N_smooth_hassas(t) - target
        Fp = np.log(t/(2*np.pi))/(2*np.pi) \
             - 1.0/(48*np.pi*t**2) \
             - 21.0/(5760*np.pi*t**4)
        t = t - F/Fp
        t = np.maximum(t, 1.0)
    return t

def primes_upto(N):
    s = np.ones(N+1, bool); s[:2] = False
    for i in range(2, int(N**0.5)+1):
        if s[i]: s[i*i::i] = False
    return np.flatnonzero(s)

N_use = 4000
n_arr = np.arange(1, N_use+1)
gamma_use = gamma[:N_use]
t_pred = t_of_n_hassas(n_arr)
delta_emp = gamma_use - t_pred

# Explicit formula öngörü
P = primes_upto(500)
delta_pred = np.zeros(N_use)
for i, gn in enumerate(gamma_use):
    log_factor = np.log(gn / (2*np.pi))
    s = 0.0
    pk_lim = (2*np.pi*gn) ** 0.5
    for p in P:
        if p > pk_lim: break
        log_p = np.log(p)
        for k in [1, 2, 3]:
            pk = p**k
            if pk > pk_lim: break
            s += np.sin(gn * k * log_p) / (k * np.sqrt(pk))
    delta_pred[i] = (2.0 / log_factor) * s

artik = delta_emp - delta_pred
artik_centered = artik - artik.mean()

print(f"Artık istatistik:")
print(f"  mean = {artik.mean():+.5f}")
print(f"  std  = {artik.std():.5f}")
print(f"  N    = {len(artik)}")

# ============================================================
# Autocorrelation
# ============================================================
def autocorr(x, max_lag):
    n = len(x)
    var = np.var(x)
    return np.array([np.mean(x[:n-k] * x[k:]) / var for k in range(max_lag+1)])

max_lag = 50
C_emp = autocorr(artik_centered, max_lag)

# Baseline: artığı KARIŞTIR (null hypothesis)
rng = np.random.default_rng(7)
N_shuffle = 100
C_null_all = []
for _ in range(N_shuffle):
    shuffled = rng.permutation(artik_centered)
    C_null_all.append(autocorr(shuffled, max_lag))
C_null = np.array(C_null_all)
C_null_mean = C_null.mean(axis=0)
C_null_std = C_null.std(axis=0)

# Sapma normalize: kaç sigma uzakta?
sigma = (C_emp - C_null_mean) / C_null_std

print(f"\nAutocorrelation (ilk 10 gecikme):")
print(f"{'k':>3} {'C_emp':>10} {'C_null (μ±σ)':>20} {'sigma':>8}")
for k in range(11):
    print(f"{k:>3} {C_emp[k]:>+10.4f} {C_null_mean[k]:>+10.4f}±{C_null_std[k]:.4f}  {sigma[k]:>+8.2f}σ")

# ============================================================
# Plot
# ============================================================
fig, axes = plt.subplots(2, 1, figsize=(14, 9))

# (1) Autocorrelation
ax = axes[0]
ax.plot(range(max_lag+1), C_emp, "o-", color="crimson", lw=1.5, ms=5,
        label="EMPİRİK artık autocorr")
# Null bandı (±2σ)
ax.fill_between(range(max_lag+1),
                C_null_mean - 2*C_null_std,
                C_null_mean + 2*C_null_std,
                color="grey", alpha=0.25, label="null (shuffle) ±2σ")
ax.plot(range(max_lag+1), C_null_mean, "k--", lw=0.7, alpha=0.5, label="null ortalaması")
ax.axhline(0, color="black", lw=0.5)
ax.set_xlabel("gecikme k (sıfır)")
ax.set_ylabel("C(k) = ⟨artık(n) · artık(n+k)⟩")
ax.set_title(f"Artığın AUTOCORRELATION'ı — sıfırların kuantum imzası var mı?\n"
             f"C(1) = {C_emp[1]:+.4f}, null beklenti = {C_null_mean[1]:+.4f} ± {C_null_std[1]:.4f}, "
             f"yani {sigma[1]:+.1f}σ")
ax.legend(fontsize=10)
ax.grid(alpha=0.3)
ax.set_xlim(0, max_lag)

# (2) Sigma (kaç sigma anlamlı)
ax = axes[1]
ax.bar(range(max_lag+1), sigma, color="purple", alpha=0.7, edgecolor="darkmagenta")
ax.axhline(0, color="black", lw=0.5)
ax.axhline(2, color="red", lw=1, ls="--", alpha=0.5, label="anlamlılık eşiği (2σ)")
ax.axhline(-2, color="red", lw=1, ls="--", alpha=0.5)
ax.set_xlabel("gecikme k (sıfır)")
ax.set_ylabel("(C_emp − C_null) / σ_null")
ax.set_title("Anlamlılık testi — hangi gecikmelerde gerçek korelasyon var?")
ax.legend(fontsize=10)
ax.grid(alpha=0.3)
ax.set_xlim(-0.5, max_lag+0.5)

plt.suptitle("14 — Artık Autocorrelation: GUE Level Repulsion İmzası mı, Sayısal Gürültü mü?",
             fontsize=12, fontweight="bold", y=0.998)
plt.tight_layout()

out = Path(__file__).parent / "14_autocorrelation.png"
plt.savefig(out, dpi=130, bbox_inches="tight")
print(f"\nKaydedildi: {out}")

# ============================================================
# Yorumlama
# ============================================================
significant_neg = np.where((sigma < -2) & (np.arange(max_lag+1) > 0))[0]
significant_pos = np.where((sigma > 2) & (np.arange(max_lag+1) > 0))[0]

print(f"\n=== YORUM ===")
print(f"Anlamlı NEGATİF korelasyonlu gecikmeler (k>0): {len(significant_neg)}")
if len(significant_neg) > 0:
    print(f"  Liste: {list(significant_neg)}")
print(f"Anlamlı POZİTİF korelasyonlu gecikmeler (k>0): {len(significant_pos)}")
if len(significant_pos) > 0:
    print(f"  Liste: {list(significant_pos)}")

if abs(sigma[1]) < 2:
    print(f"\nC(1) anlamsız → ardışık artıklar bağımsız → SAYISAL ARTIK")
elif sigma[1] < -2:
    print(f"\nC(1) güçlü negatif → ardışık artıklar TERS YÖNDE → LEVEL REPULSION imzası")
elif sigma[1] > 2:
    print(f"\nC(1) güçlü pozitif → ardışık artıklar AYNI YÖNDE → trend/clustering var")
