"""
10 — Riemann Sıfırları Bir Dalga Deseni mi?
=============================================

Uğur'un sorusu: γ_n'leri parçacık gibi düşün — plakada nereye düşüyorlar?
Düzgün eşit aralıkta mı (bir merkez), yoksa dalga gibi salınıyorlar mı
(çok merkez)?

Test:
  1) Berry-Keating düzgün öngörü: t̃_n = N_BK⁻¹(n)
  2) Gerçek sıfır γ_n
  3) Sapma δ_n = γ_n - t̃_n  →  dalga gibi mi?

Eğer δ_n bir dalga deseni gösteriyorsa, sıfırlar gerçekten "girişim
çıktısı" gibi davranıyor demektir.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
gamma = np.loadtxt(ROOT / "zeros_100k.txt")

# ============================================================
# Berry-Keating düzgün öngörü: N_BK(t) = (t/2π)log(t/(2πe)) + 7/8
# t̃_n: N_BK(t̃_n) = n  (Newton ile çöz)
# ============================================================
def t_of_n_BK(n):
    """n. sıfırın Berry-Keating yarı-klasik öngörüsü"""
    # ilk tahmin: t ~ 2πn / log(n)
    t = 2*np.pi * n / max(np.log(n+2), 1)
    for _ in range(60):
        F = (t/(2*np.pi))*np.log(t/(2*np.pi*np.e)) + 7/8 - n
        Fp = np.log(t/(2*np.pi*np.e))/(2*np.pi) + 1/(2*np.pi)
        t = t - F/Fp
        if t < 1: t = 1.0
    return t

N_use = 1000
t_pred = np.array([t_of_n_BK(n) for n in range(1, N_use+1)])
gamma_use = gamma[:N_use]

# Sapma (parçacığın "yer kayması")
delta = gamma_use - t_pred

print(f"Kullanılan sıfır sayısı: {N_use}")
print(f"Sapma istatistik: mean={delta.mean():.4f}, std={delta.std():.4f}")
print(f"Min sapma: {delta.min():.4f}, Max sapma: {delta.max():.4f}")

# ============================================================
# Plot — üç katman
# ============================================================
fig, axes = plt.subplots(3, 1, figsize=(14, 11))

# (1) Sıfırlar "plakada" düşmüş parçacıklar
ax = axes[0]
ax.scatter(gamma_use, np.zeros(N_use), s=12, c="red",
           alpha=0.6, edgecolors="darkred", linewidths=0.4)
ax.set_yticks([])
ax.set_xlim(0, gamma_use[-1])
ax.set_xlabel("t — plaka ekseni (sayma yönü)", fontsize=10)
ax.set_title(f"İlk {N_use} ζ sıfırı, 'plakada' düşmüş parçacıklar gibi",
             fontsize=11, fontweight="bold")
# Yer yer eşit aralık kontrol çizgileri
for tt in np.arange(0, gamma_use[-1], 100):
    ax.axvline(tt, color="grey", lw=0.3, alpha=0.3)
ax.text(gamma_use[-1]*0.5, 0.3,
        "Gözle: noktalar düzgün dağılmış gibi ama yer yer sıkışık/seyrek",
        ha="center", fontsize=9, style="italic", color="dimgrey")
ax.set_ylim(-0.5, 0.5)

# (2) Beklenen (Berry-Keating düzgün) vs gerçek
ax = axes[1]
n_arr = np.arange(1, N_use+1)
ax.plot(n_arr, gamma_use, "ro-", ms=2.5, lw=0.5, label="gerçek γ_n", alpha=0.7)
ax.plot(n_arr, t_pred, "b-", lw=1.0, alpha=0.6, label="Berry-Keating düzgün öngörü")
ax.set_xlabel("n", fontsize=10)
ax.set_ylabel("t değeri", fontsize=10)
ax.set_title("Gerçek sıfırlar düzgün öngörüye yapışık görünüyor — büyük ölçekte tek hat",
             fontsize=11, fontweight="bold")
ax.legend(loc="lower right", fontsize=10)
ax.grid(alpha=0.3)

# (3) SAPMA — dalga deseni var mı?
ax = axes[2]
ax.plot(n_arr, delta, color="crimson", lw=0.8, alpha=0.85)
ax.axhline(0, color="black", lw=0.6)
ax.fill_between(n_arr, delta, 0, where=(delta > 0), alpha=0.2, color="red")
ax.fill_between(n_arr, delta, 0, where=(delta < 0), alpha=0.2, color="blue")
ax.set_xlabel("n", fontsize=10)
ax.set_ylabel("γ_n − t̃_n  (sapma)", fontsize=10)
ax.set_title("SAPMA: gerçek sıfırlar düzgün öngörünün etrafında DALGA gibi salınıyor\n"
             "(0 etrafında: ne sadece bir yönde toplanma var ne de rastgele)",
             fontsize=11, fontweight="bold")
ax.grid(alpha=0.3)
ax.set_xlim(1, N_use)

# Sayısal: kaç salınım?
zero_crossings = np.sum(np.diff(np.sign(delta)) != 0)
print(f"Sıfır geçişi sayısı (delta'da): {zero_crossings}")
print(f"Yani yaklaşık {zero_crossings/2:.0f} 'dalga' var ilk {N_use} sıfırda")
print(f"Ortalama 'dalga genişliği': {N_use / max(zero_crossings/2, 1):.1f} sıfır")

plt.suptitle("Riemann Sıfırları PARÇACIK olarak — bir merkez mi, çok merkez mi?",
             fontsize=13, fontweight="bold", y=0.998)
plt.tight_layout()

out = Path(__file__).parent / "10_dalga_deseni.png"
plt.savefig(out, dpi=130, bbox_inches="tight")
print(f"\nKaydedildi: {out}")

# ============================================================
# Bonus: sapmanın kendisinin FOURIER'i
# (eğer dalga deseniyse, periyodik bileşenleri olmalı)
# ============================================================
from scipy.fft import rfft, rfftfreq
# n eksen mesafe = 1, sample rate = 1
fft_vals = np.abs(rfft(delta - delta.mean()))
freqs = rfftfreq(len(delta), d=1.0)

print(f"\nSapma sinyalinin en güçlü Fourier frekansları:")
top_k = np.argsort(fft_vals)[-5:][::-1]
for k in top_k:
    if freqs[k] > 0:
        period = 1/freqs[k]
        print(f"  frekans = {freqs[k]:.4f}, periyot = {period:.1f} sıfır, amplitüd = {fft_vals[k]:.2f}")
