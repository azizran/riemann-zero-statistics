"""
QM × Riemann — Asalların Sesi
==============================

Form factor'un GUE'den sapması rastgele değil. Riemann'ın explicit formula'sı:

    Σ_γ f(γ)  ↔  Σ_p (log p)/√p · f̂(log p) + ...

Yani sıfırlar (γ) ve asallar (p) Fourier eşleridir. Sıfırların Fourier
dönüşümünde **asalların logaritmalarında** pikler görmeliyiz.

Bu test: Σ_n cos(γ_n · u) çiz, u-ekseni boyunca log(p), log(p²),... noktalarına bak.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
t_all = np.loadtxt(ROOT / "zeros_100k.txt")

# ---------- Asalları üret (Eratosthenes) ----------
def primes_upto(N):
    sieve = np.ones(N+1, bool); sieve[:2] = False
    for i in range(2, int(N**0.5)+1):
        if sieve[i]: sieve[i*i::i] = False
    return np.flatnonzero(sieve)

P = primes_upto(200)   # u ekseninde log p, log p², log p³ vs

# Prime power'ları üret: n = p^k, ağırlık Λ(n)/√n = log(p)/√(p^k)
pp = []   # (log(n), weight)
for p in P:
    pk = p
    while pk <= 200:
        pp.append((np.log(pk), np.log(p) / np.sqrt(pk), int(pk), p))
        pk *= p
pp.sort()

# ---------- ζ sıfırlarının Fourier dönüşümü ----------
# F(u) = Σ_n cos(γ_n · u)   ,  γ_n ∈ zero seti, sınırlı pencerede
# pencere: orta yükseklik civarı, Gaussian ağırlık ile (kenar etkisini bastır)
GAMMA_LOW, GAMMA_HIGH = 100, 1000       # γ ∈ [100, 1000]
mask = (t_all > GAMMA_LOW) & (t_all < GAMMA_HIGH)
g = t_all[mask]
print(f"Kullanılan sıfır sayısı: {len(g)}, [{g[0]:.2f}, {g[-1]:.2f}]")

# Gaussian pencere (orta=550, σ=300)
center, sigma = 0.5*(GAMMA_LOW+GAMMA_HIGH), 0.4*(GAMMA_HIGH-GAMMA_LOW)
w = np.exp(-(g - center)**2 / (2*sigma**2))

# u-ekseni: log p civarında çözünürlük lazım, log 200 ≈ 5.3
u = np.linspace(0.3, 5.5, 4000)

# F(u) = Σ_n w_n · cos(γ_n · u)
F = (w[:, None] * np.cos(np.outer(g, u))).sum(axis=0)
F /= w.sum()    # normalize

# ---------- Plot ----------
fig, axes = plt.subplots(2, 1, figsize=(13, 8),
                         gridspec_kw={"height_ratios": [3, 1]})

# Üst: Fourier dönüşümü + asal işaretleri
ax = axes[0]
ax.plot(u, F, color="steelblue", lw=0.8, alpha=0.8)
ax.axhline(0, color="k", lw=0.4)

# Asal yerleri işaretle
for log_n, weight, n, p in pp:
    if log_n < u[0] or log_n > u[-1]: continue
    # ağırlık ile orantılı dikey çizgi
    h = weight * 0.4
    color = "crimson" if n == p else "orange"      # asal kırmızı, asal-üs turuncu
    ax.axvline(log_n, ymin=0.5, ymax=0.5+h*0.5, color=color, lw=1.5, alpha=0.85)
    if weight > 0.15:
        label = f"{p}" if n == p else f"{p}^{int(np.log(n)/np.log(p))}"
        ax.text(log_n, 0.02 + h*1.05, label, ha="center", fontsize=8, color=color)

ax.set_xlabel("u  (Fourier ekseni)")
ax.set_ylabel(r"$F(u) = \sum_n w_n \cos(\gamma_n u)$")
ax.set_title("ζ sıfırlarının Fourier dönüşümü — asalların log-pozisyonlarında pikler\n"
             "(900 sıfır, γ ∈ [100, 1000])")
ax.set_xlim(u[0], u[-1])
ax.grid(alpha=0.25)

# Alt: gerçek asalların konumlarını ayrı bantta göster (legend gibi)
ax2 = axes[1]
for log_n, weight, n, p in pp:
    if log_n < u[0] or log_n > u[-1]: continue
    color = "crimson" if n == p else "orange"
    ax2.vlines(log_n, 0, weight, color=color, lw=2)
ax2.set_xlim(u[0], u[-1])
ax2.set_xlabel("u  =  log(asal veya asal-üs)")
ax2.set_ylabel(r"ağırlık $\Lambda(n)/\sqrt{n}$")
ax2.set_title("Beklenen pik yerleri (von Mangoldt ağırlığı)")
ax2.grid(alpha=0.25)

# legend
from matplotlib.lines import Line2D
ax2.legend(handles=[
    Line2D([0],[0], color="crimson", lw=2, label="asal: p"),
    Line2D([0],[0], color="orange", lw=2, label="asal-üs: p^k"),
], loc="upper right", fontsize=9)

plt.tight_layout()
out = Path(__file__).parent / "03_asallarin_sesi.png"
plt.savefig(out, dpi=130, bbox_inches="tight")
print(f"Kaydedildi: {out}")

# ---------- Sayısal: pik yerleri ζ ve beklenti ne kadar uyuşuyor? ----------
# Empirik tepe noktaları bul (basit: lokal max + threshold)
from scipy.signal import find_peaks
peaks, _ = find_peaks(np.abs(F), height=0.15*np.max(np.abs(F)), distance=20)
u_peaks = u[peaks]
F_peaks = F[peaks]

print(f"\nEmpirik tepe sayısı: {len(u_peaks)}")
print(f"{'u_peak':>8} | {'|F|':>8} | en yakın log(p^k) | sapma")
for up, fp in zip(u_peaks, F_peaks):
    # en yakın beklenen pik
    nearest = min(pp, key=lambda x: abs(x[0] - up))
    log_n, weight, n, p = nearest
    delta = up - log_n
    label = str(p) if n == p else f"{p}^{int(np.log(n)/np.log(p))}"
    print(f"{up:8.3f} | {fp:+8.3f} | log({label})={log_n:.3f}  | Δ={delta:+.3f}")
