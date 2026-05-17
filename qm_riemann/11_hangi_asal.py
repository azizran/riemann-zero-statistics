"""
11 — Sapmadaki dalgalar hangi asallardan geliyor?
====================================================

Riemann explicit formula:
  δ_n ≈ -S(γ_n)/ρ(γ_n) = (2/log(t̄/2π)) · Σ_{p,k} sin(γ_n · k·log p)/(k·√(p^k))

Yani her asal p (ve onun her k. kuvveti):
  - n-periyot:  T_n(p,k) = log(t̄/2π) / (k·log p)
  - n-genlik:   A_n(p,k) ∝ 1 / (k·√(p^k))

Test: sapma δ_n'in Fourier spektrumu (n bazında) tam bu noktalarda pik vermeli.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.fft import rfft, rfftfreq

ROOT = Path(__file__).resolve().parent.parent
gamma = np.loadtxt(ROOT / "zeros_100k.txt")

# Berry-Keating düzgün öngörü (Newton inversion)
def t_of_n_BK(n):
    t = max(2*np.pi * n / max(np.log(n+2), 1), 1.0)
    for _ in range(60):
        F = (t/(2*np.pi))*np.log(t/(2*np.pi*np.e)) + 7/8 - n
        Fp = np.log(t/(2*np.pi*np.e))/(2*np.pi) + 1/(2*np.pi)
        t = t - F/Fp
        if t < 1: t = 1.0
    return t

# ============================================================
# Sapma serisi
# ============================================================
N_use = 4000           # daha çok sıfır → daha iyi frekans çözünürlüğü
t_pred = np.array([t_of_n_BK(n) for n in range(1, N_use+1)])
gamma_use = gamma[:N_use]
delta = gamma_use - t_pred
delta_centered = delta - delta.mean()

t_bar = gamma_use.mean()
log_factor = np.log(t_bar / (2*np.pi))
print(f"N={N_use} sıfır, t̄={t_bar:.1f}, log(t̄/2π)={log_factor:.4f}")

# ============================================================
# Fourier (n bazında)
# ============================================================
fft_vals = np.abs(rfft(delta_centered))
freqs = rfftfreq(len(delta_centered), d=1.0)     # 1/sıfır birimi
# Çok düşük frekansları (mean drift) bastır
fft_vals[0] = 0

# ============================================================
# Asalların öngörü pozisyonları
# ============================================================
def primes_upto(N):
    s = np.ones(N+1, bool); s[:2] = False
    for i in range(2, int(N**0.5)+1):
        if s[i]: s[i*i::i] = False
    return np.flatnonzero(s)

P = primes_upto(100)

# Beklenen pik listesi: (p, k, frekans = 1/periyot, genlik öngörüsü)
expected = []
for p in P:
    for k in [1, 2, 3]:
        period = log_factor / (k * np.log(p))
        if period < 1.8: break    # çözünürlük altı
        freq  = 1.0 / period
        amp   = 1.0 / (k * np.sqrt(p**k))
        expected.append((p, k, freq, period, amp))

print(f"\nBeklenen pik tablosu (asal + kuvvet):")
print(f"{'p':>5} {'k':>3} {'periyot_n':>10} {'frekans':>10} {'genlik':>10}")
for p, k, f, T, A in expected[:20]:
    print(f"{p:>5} {k:>3} {T:>10.3f} {f:>10.4f} {A:>10.4f}")

# ============================================================
# Empirik pikleri bul
# ============================================================
from scipy.signal import find_peaks
# normalize amplitüdü
fft_normalized = fft_vals / fft_vals.max()
peaks, props = find_peaks(fft_normalized, height=0.15, distance=4)
peak_freqs = freqs[peaks]
peak_amps = fft_normalized[peaks]

print(f"\nEmpirik tepe sayısı (yükseklik > %15): {len(peaks)}")
print(f"{'empirik freq':>14} {'empirik T':>12} {'en yakın asal':>18} {'sapma':>10}")
matches = []
for f_emp, a_emp in zip(peak_freqs, peak_amps):
    if f_emp < 1e-3: continue
    T_emp = 1/f_emp
    # en yakın beklenen pik
    best = min(expected, key=lambda e: abs(e[2] - f_emp))
    p, k, f_pred, T_pred, A_pred = best
    label = str(p) if k == 1 else f"{p}^{k}"
    print(f"{f_emp:>14.4f} {T_emp:>12.3f} {label:>18} {f_emp-f_pred:>+10.4f}")
    matches.append((f_emp, T_emp, p, k, T_pred))

# ============================================================
# Plot — empirik spektrum + asal öngörüleri
# ============================================================
fig, axes = plt.subplots(2, 1, figsize=(14, 9))

# Üst: sapmanın kendisi
ax = axes[0]
ax.plot(np.arange(1, N_use+1), delta, color="crimson", lw=0.5, alpha=0.8)
ax.axhline(0, color="black", lw=0.6)
ax.set_xlabel("n")
ax.set_ylabel("δ_n = γ_n − t̃_n")
ax.set_title(f"Sıfırların sapması (N={N_use} sıfır) — bu dalgaların kaynağı asallar mı?")
ax.grid(alpha=0.3)

# Alt: spektrum + asal işaretleri
ax = axes[1]
ax.semilogy(freqs[1:], fft_normalized[1:], color="navy", lw=0.7, alpha=0.85,
            label="empirik FFT |Re δ_n|")
# Asal pozisyonları
for p, k, f_pred, T_pred, A_pred in expected:
    if f_pred > freqs.max(): continue
    color = "crimson" if k == 1 else "orange"
    ax.axvline(f_pred, color=color, lw=1.0, ls="--", alpha=0.55)
    label = str(p) if k == 1 else f"{p}^{k}"
    if A_pred > 0.15:
        ax.text(f_pred, 1.2, label, ha="center", fontsize=8.5,
                color=color, fontweight="bold")

ax.set_xlabel("frekans (1/sıfır)")
ax.set_ylabel("|FFT(δ)|  (normalize)")
ax.set_title("Sapmanın Fourier spektrumu — kırmızı kesikli: asal öngörüleri; "
             "turuncu: asal-üs")
ax.set_xlim(0, 0.6)
ax.set_ylim(1e-3, 2)
ax.grid(alpha=0.3, which="both")
ax.legend(loc="upper right", fontsize=9)

plt.suptitle("11 — DALGALAR HANGİ ASALLARDAN GELİYOR?\n"
             "Explicit formula öngörüsü vs empirik FFT",
             fontsize=13, fontweight="bold", y=1.00)
plt.tight_layout()

out = Path(__file__).parent / "11_hangi_asal.png"
plt.savefig(out, dpi=130, bbox_inches="tight")
print(f"\nKaydedildi: {out}")

# ============================================================
# Sayısal toplam: kaç pikin en yakın "asal" mesafesi nedir
# ============================================================
if matches:
    deltas_pred = [abs(m[0] - 1/m[4]) for m in matches]
    print(f"\nTüm eşleşmelerin sapma istatistiği:")
    print(f"  Medyan sapma:     {np.median(deltas_pred):.5f}")
    print(f"  Maksimum sapma:   {max(deltas_pred):.5f}")
    print(f"  Frekans çözünürlüğü Δf = 1/N = {1/N_use:.5f}")
