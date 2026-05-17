"""
13 — Gonek-Hassas Test: Konvansiyon + Higher-Order Düzeltmeler
==================================================================

İki düzeltme:
  (1) N(γ_n) = n − 1/2  (konvansiyon hatası düzelt)
  (2) Stirling açılımının yüksek mertebe terimleri:
      N_hassas(t) = (t/2π)log(t/2πe) + 7/8 + 1/(48πt) + 7/(5760π t³) + ...

Sonra: empirik δ − bu hassas öngörü = artık → orada gerçekten yapı var mı?
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.fft import rfft, rfftfreq
from scipy.signal import find_peaks

ROOT = Path(__file__).resolve().parent.parent
gamma = np.loadtxt(ROOT / "zeros_100k.txt")

# ============================================================
# Hassas N(t) — yüksek mertebe asimptotik (S(t) hariç, sadece düzgün kısım)
# ============================================================
def N_smooth_hassas(t):
    """θ(t) Stirling açılımıyla N(t) düzgün kısmı"""
    base = (t/(2*np.pi)) * np.log(t/(2*np.pi*np.e)) + 7/8
    corr1 = 1.0 / (48 * np.pi * t)
    corr2 = 7.0 / (5760 * np.pi * t**3)
    return base + corr1 + corr2

def t_of_n_hassas(n_arr):
    """N_hassas(t̃) = n − 1/2 çöz (Newton)"""
    target = n_arr - 0.5    # DÜZELTİLEN konvansiyon
    t = np.maximum(2*np.pi * n_arr / np.maximum(np.log(n_arr+2), 1), 1.0)
    for _ in range(80):
        F = N_smooth_hassas(t) - target
        Fp = np.log(t/(2*np.pi))/(2*np.pi) \
             - 1.0/(48*np.pi*t**2) \
             - 21.0/(5760*np.pi*t**4)
        t = t - F/Fp
        t = np.maximum(t, 1.0)
    return t

# Eski (kaba) Berry-Keating
def t_of_n_BK_eski(n_arr):
    t = np.maximum(2*np.pi * n_arr / np.maximum(np.log(n_arr+2), 1), 1.0)
    for _ in range(60):
        F = (t/(2*np.pi))*np.log(t/(2*np.pi*np.e)) + 7/8 - n_arr   # n − 1/2 değil
        Fp = np.log(t/(2*np.pi*np.e))/(2*np.pi) + 1/(2*np.pi)
        t = t - F/Fp
        t = np.maximum(t, 1.0)
    return t

# ============================================================
# Empirik sapma — iki versiyonla karşılaştır
# ============================================================
N_use = 4000
n_arr = np.arange(1, N_use+1)
gamma_use = gamma[:N_use]

t_eski = t_of_n_BK_eski(n_arr)
t_yeni = t_of_n_hassas(n_arr)

delta_eski = gamma_use - t_eski
delta_yeni = gamma_use - t_yeni

print(f"İki öngörünün karşılaştırması (N={N_use} sıfır)")
print(f"")
print(f"ESKİ (N(t)=n, sadece 7/8):")
print(f"  mean = {delta_eski.mean():+.5f}")
print(f"  std  = {delta_eski.std():.5f}")
print(f"")
print(f"YENİ (N(t)=n−1/2, + higher-order):")
print(f"  mean = {delta_yeni.mean():+.5f}")
print(f"  std  = {delta_yeni.std():.5f}")
print(f"")
print(f"İki t̃ arasındaki ortalama fark: {(t_yeni-t_eski).mean():+.5f}")

# ============================================================
# Şimdi YENİ artık vs explicit formula öngörüsü
# ============================================================
def primes_upto(N):
    s = np.ones(N+1, bool); s[:2] = False
    for i in range(2, int(N**0.5)+1):
        if s[i]: s[i*i::i] = False
    return np.flatnonzero(s)

P = primes_upto(500)
print(f"\nExplicit formula için {len(P)} asal kullanılıyor")

# δ_pred = (2/log(γ/2π)) Σ sin(γ × k log p) / (k √(p^k))
delta_pred_yeni = np.zeros(N_use)
for i, gn in enumerate(gamma_use):
    log_factor = np.log(gn / (2*np.pi))
    s = 0.0
    for p in P:
        log_p = np.log(p)
        pk_lim = (2*np.pi*gn) ** 0.5     # GUE truncation seviyesi
        if p > pk_lim: break
        for k in [1, 2, 3]:
            pk = p**k
            if pk > pk_lim: break
            s += np.sin(gn * k * log_p) / (k * np.sqrt(pk))
    delta_pred_yeni[i] = (2.0 / log_factor) * s

# YENİ konvansiyonla artık
artik_yeni = delta_yeni - delta_pred_yeni

print(f"\nDüzeltilmiş empirik vs explicit formula:")
print(f"  Korelasyon: {np.corrcoef(delta_yeni, delta_pred_yeni)[0,1]:.6f}")
print(f"  Artık mean: {artik_yeni.mean():+.5f}")
print(f"  Artık std:  {artik_yeni.std():.5f}")
print(f"  Artık/δ_yeni: {artik_yeni.std()/delta_yeni.std():.2%}")

# Artığın spektral yapısı
artik_centered = artik_yeni - artik_yeni.mean()
fft_artik = np.abs(rfft(artik_centered))
freqs = rfftfreq(N_use, d=1.0)
fft_artik[0] = 0

# En güçlü frekanslar
peaks, _ = find_peaks(fft_artik, height=fft_artik.max()*0.3, distance=5)
print(f"\nARTIK sinyalin önemli pikler:")
for k in peaks[:15]:
    if freqs[k] > 0.005:
        # Asal eşleşmesi var mı?
        log_factor_bar = np.log(gamma_use.mean() / (2*np.pi))
        # Bu frekansa en yakın asal
        f_emp = freqs[k]
        best_match = None
        best_err = 1e9
        for p in P:
            for kk in [1, 2, 3]:
                f_pred = kk * np.log(p) / log_factor_bar
                if abs(f_pred - f_emp) < best_err:
                    best_err = abs(f_pred - f_emp)
                    best_match = (p, kk)
        if best_match:
            p, kk = best_match
            label = f"{p}^{kk}" if kk > 1 else f"{p}"
            print(f"  freq = {freqs[k]:.4f}, periyot = {1/freqs[k]:.2f}, amp = {fft_artik[k]:.3f}  → asal {label} (sapma {best_err:.4f})")

# ============================================================
# Plot
# ============================================================
fig, axes = plt.subplots(3, 1, figsize=(14, 11))

# (1) Sapma karşılaştırması: eski vs yeni öngörü
ax = axes[0]
ax.plot(n_arr, delta_eski, color="grey", lw=0.5, alpha=0.7, label=f"eski (mean={delta_eski.mean():+.3f})")
ax.plot(n_arr, delta_yeni, color="crimson", lw=0.5, alpha=0.85, label=f"yeni hassas (mean={delta_yeni.mean():+.3f})")
ax.axhline(0, color="black", lw=0.5)
ax.set_xlabel("n")
ax.set_ylabel("δ_n = γ_n − t̃_n")
ax.set_title("ESKİ vs YENİ öngörü — konvansiyon ve higher-order düzeltmesinin etkisi")
ax.legend(fontsize=10)
ax.grid(alpha=0.3)

# (2) Yeni artık: empirik − explicit formula
ax = axes[1]
ax.plot(n_arr, artik_yeni, color="purple", lw=0.6, alpha=0.85)
ax.axhline(0, color="black", lw=0.5)
ax.set_xlabel("n")
ax.set_ylabel("artık")
ax.set_title(f"GERÇEK ARTIK (yeni öngörü − explicit formula)\n"
             f"std/δ = {artik_yeni.std()/delta_yeni.std():.2%}, "
             f"korelasyon empirik-öngörü = {np.corrcoef(delta_yeni, delta_pred_yeni)[0,1]:.3f}")
ax.grid(alpha=0.3)

# (3) Artığın spektrumu
ax = axes[2]
ax.semilogy(freqs[1:], fft_artik[1:] / fft_artik[1:].max(),
            color="purple", lw=0.7, alpha=0.85)
# Asal pozisyonları (gri)
log_factor_bar = np.log(gamma_use.mean() / (2*np.pi))
for p in [2,3,5,7,11,13,17,19,23,29,31]:
    f = np.log(p) / log_factor_bar
    ax.axvline(f, color="grey", lw=0.5, ls=":", alpha=0.5)
    ax.text(f, 1.2, str(p), ha="center", fontsize=7, color="grey")
ax.set_xlim(0, 0.5)
ax.set_ylim(1e-3, 2)
ax.set_xlabel("frekans (1/sıfır)")
ax.set_ylabel("|FFT(artık)| (normalize)")
ax.set_title("Artığın spektrumu — eğer asal pozisyonlarında pik kalmadıysa, "
             "explicit formula RESMİ TAMAMLAMIŞ demek")
ax.grid(alpha=0.3, which="both")

plt.suptitle("13 — Gonek-Hassas Test: gerçekten 'fazla' bir şey var mı?",
             fontsize=12, fontweight="bold", y=0.998)
plt.tight_layout()

out = Path(__file__).parent / "13_gonek_test.png"
plt.savefig(out, dpi=130, bbox_inches="tight")
print(f"\nKaydedildi: {out}")
