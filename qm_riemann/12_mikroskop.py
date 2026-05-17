"""
12 — Tek Asal Tepesini Mikroskopla İncele
=============================================

Soru: Asal 2 tepesinin etrafındaki ondalik dalga deseni
  (a) sadece sonlu-pencere artefakı (sinc/Dirichlet zarfı) mi?
  (b) yoksa fazla bilgi (asal komşulukları, ikiz asallar) mı taşıyor?

Test: explicit formula öngörüsünü TAM olarak hesapla, empirik sapma ile karşılaştır.
Eğer ikisi mükemmel uyuşuyorsa → sadece sinc + girişim.
Eğer fark varsa → asal-üstü bir yapı (örn. asal korelasyonları).
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.fft import rfft, rfftfreq

ROOT = Path(__file__).resolve().parent.parent
gamma = np.loadtxt(ROOT / "zeros_100k.txt")

# Berry-Keating öngörü
def t_of_n_BK_vec(n_arr):
    t = np.maximum(2*np.pi * n_arr / np.maximum(np.log(n_arr+2), 1), 1.0)
    for _ in range(60):
        F = (t/(2*np.pi))*np.log(t/(2*np.pi*np.e)) + 7/8 - n_arr
        Fp = np.log(t/(2*np.pi*np.e))/(2*np.pi) + 1/(2*np.pi)
        t = t - F/Fp
        t = np.maximum(t, 1.0)
    return t

# ============================================================
# Empirik sapma
# ============================================================
N_use = 4000
n_arr = np.arange(1, N_use+1)
t_pred_BK = t_of_n_BK_vec(n_arr)
gamma_use = gamma[:N_use]
delta_emp = gamma_use - t_pred_BK

# ============================================================
# Explicit formula öngörüsü — δ_n^pred
# ============================================================
def primes_upto(N):
    s = np.ones(N+1, bool); s[:2] = False
    for i in range(2, int(N**0.5)+1):
        if s[i]: s[i*i::i] = False
    return np.flatnonzero(s)

P = primes_upto(200)   # bütün ilgili asallar
print(f"Asal sayısı: {len(P)}")

delta_pred = np.zeros(N_use)
for i, gn in enumerate(gamma_use):
    log_factor = np.log(gn / (2*np.pi))
    s = 0.0
    for p in P:
        log_p = np.log(p)
        for k in [1, 2, 3]:
            pk = p**k
            if pk > 1000: break
            s += np.sin(gn * k * log_p) / (k * np.sqrt(pk))
    delta_pred[i] = (2.0 / log_factor) * s

# Karşılaştırma
print(f"\nEmpirik δ:  mean = {delta_emp.mean():+.4f}, std = {delta_emp.std():.4f}")
print(f"Öngörü δ:    mean = {delta_pred.mean():+.4f}, std = {delta_pred.std():.4f}")
print(f"Korelasyon:  {np.corrcoef(delta_emp, delta_pred)[0,1]:.6f}")
print(f"Fark        std: {(delta_emp - delta_pred).std():.4f}")
print(f"Fark/orijinal: {(delta_emp - delta_pred).std() / delta_emp.std():.3%}")

# ============================================================
# FFT — empirik ve öngörü
# ============================================================
fft_emp = np.abs(rfft(delta_emp - delta_emp.mean()))
fft_pred = np.abs(rfft(delta_pred - delta_pred.mean()))
freqs = rfftfreq(N_use, d=1.0)

# Asal 2'nin tahmini frekansı
t_bar = gamma_use.mean()
log_factor_bar = np.log(t_bar / (2*np.pi))
f_prime2 = np.log(2) / log_factor_bar
print(f"\nAsal 2'nin merkez frekansı: {f_prime2:.5f}")

# ============================================================
# Plot — 3 panel
# ============================================================
fig, axes = plt.subplots(3, 1, figsize=(14, 11))

# (1) Tam spektrum karşılaştırma
ax = axes[0]
ax.semilogy(freqs[1:], fft_emp[1:] / fft_emp[1:].max(),
            color="navy", lw=0.7, label="empirik (gerçek sıfırlar)")
ax.semilogy(freqs[1:], fft_pred[1:] / fft_pred[1:].max(),
            color="crimson", lw=0.7, alpha=0.7, label="explicit formula öngörü")
for p in [2,3,5,7,11,13]:
    f = np.log(p) / log_factor_bar
    ax.axvline(f, color="grey", lw=0.5, ls=":", alpha=0.5)
    ax.text(f, 1.3, str(p), ha="center", fontsize=8, color="grey")
ax.set_xlim(0, 0.5)
ax.set_ylim(1e-4, 2)
ax.set_xlabel("frekans (1/sıfır)")
ax.set_ylabel("|FFT(δ)| (normalize)")
ax.set_title("Tam spektrum: empirik vs explicit formula öngörü")
ax.legend(fontsize=9)
ax.grid(alpha=0.3, which="both")

# (2) Asal 2 mikroskobu
ax = axes[1]
f_lo, f_hi = f_prime2 - 0.02, f_prime2 + 0.02
mask = (freqs > f_lo) & (freqs < f_hi)
ax.plot(freqs[mask], fft_emp[mask], "o-", color="navy", ms=3, lw=0.8,
        label="empirik")
ax.plot(freqs[mask], fft_pred[mask], "s-", color="crimson", ms=3, lw=0.8,
        alpha=0.7, label="öngörü")
ax.axvline(f_prime2, color="red", lw=1.5, ls="--", alpha=0.7, label=f"asal 2 (f={f_prime2:.4f})")
# Diğer asallar bu pencerede var mı
for p in [3,5,7,11,13]:
    f = np.log(p) / log_factor_bar
    if f_lo < f < f_hi:
        ax.axvline(f, color="orange", lw=1, ls=":", alpha=0.7)
        ax.text(f, ax.get_ylim()[1]*0.9, str(p), fontsize=8, color="orange")
ax.set_xlabel("frekans (1/sıfır)")
ax.set_ylabel("|FFT(δ)|")
ax.set_title(f"Asal 2 mikroskobu (zoom: f ∈ [{f_lo:.3f}, {f_hi:.3f}])\n"
             "yan loblar = sinc artefakı + asal komşuluk girişimi")
ax.legend(fontsize=9)
ax.grid(alpha=0.3)

# (3) ARTIK: empirik - öngörü
ax = axes[2]
residual = (delta_emp - delta_emp.mean()) - (delta_pred - delta_pred.mean())
ax.plot(n_arr, residual, color="purple", lw=0.5, alpha=0.85)
ax.axhline(0, color="black", lw=0.5)
ax.set_xlabel("n")
ax.set_ylabel(r"δ_emp − δ_öngörü")
ax.set_title(f"ARTIK — explicit formula açıklayamadığı kısım\n"
             f"std artık / std orijinal = {residual.std()/delta_emp.std():.2%}")
ax.grid(alpha=0.3)

plt.suptitle("12 — Asal 2 Tepesi Mikroskop Altında\n"
             "ondalik dalgalar: sadece sinc mi, yoksa fazla yapı mı?",
             fontsize=12, fontweight="bold", y=0.998)
plt.tight_layout()

out = Path(__file__).parent / "12_mikroskop.png"
plt.savefig(out, dpi=130, bbox_inches="tight")
print(f"\nKaydedildi: {out}")

# ============================================================
# Artığın yapısı — Fourier'ini de al
# ============================================================
res_fft = np.abs(rfft(residual - residual.mean()))
top_k = np.argsort(res_fft)[-10:][::-1]
print(f"\nARTIK sinyalin en güçlü Fourier frekansları (varsa yapı belirtir):")
for k in top_k:
    if freqs[k] < 1e-3: continue
    print(f"  freq = {freqs[k]:.4f}, periyot = {1/freqs[k]:.2f} sıfır, "
          f"amplitüd = {res_fft[k]:.3f}")
