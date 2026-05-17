"""
16 — Asansör Büyüdükçe
=========================

N=3 → N=5 → N=20 → N=100 → N=1000 GUE asansörleri.
Her birinde özdeğerleri (oturma pozisyonlarını) topla, dağılımına bak.

Beklenti: N büyüdükçe dağılım **YARIM DAİRE** şeklini alır (Wigner semicircle).

Bu Wigner'in 1955'te keşfettiği yasadır:
  ρ(λ) = (1/(2π)) √(4N − λ²),   λ ∈ [−2√N, +2√N]

Asansörün doğal şekli yarı-daire. Sonsuz büyük asansörde herkesin oturduğu
yer rastgele görünür, ama büyük resim mükemmel bir yarı çember çizer.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

rng = np.random.default_rng(11)

def gue_matrix(N):
    A = rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))
    return (A + A.conj().T) / 2

def all_eigs(N, N_trials):
    """N_trials deneme, hepsinin özdeğerlerini tek bir dizide topla"""
    arr = np.zeros((N_trials, N))
    for i in range(N_trials):
        H = gue_matrix(N)
        arr[i] = np.sort(np.linalg.eigvalsh(H))
    return arr

def semicircle(x, N):
    """Wigner yarı-daire yoğunluğu (1/(2πN))·√(4N - x²)
       ∫ρ dx = 1 olacak şekilde normalize"""
    inside = np.maximum(4*N - x**2, 0)
    return np.sqrt(inside) / (2 * np.pi * N)

# ============================================================
# Verileri üret
# ============================================================
Ns = [3, 5, 20, 100, 1000]
trials_per_N = [50000, 30000, 5000, 500, 50]

datasets = {}
for N, T in zip(Ns, trials_per_N):
    print(f"N={N:>4}, {T:>5} deneme üretiliyor...")
    eigs = all_eigs(N, T)
    datasets[N] = eigs.flatten()
    print(f"   sınırlar: [{eigs.min():.2f}, {eigs.max():.2f}],  beklenen yarı-daire sınırı: ±{2*np.sqrt(N):.2f}")

# ============================================================
# Plot — 5 panel yan yana, her biri farklı N
# ============================================================
fig, axes = plt.subplots(1, 5, figsize=(20, 5))

for ax, N in zip(axes, Ns):
    data = datasets[N]
    x_lim = 2 * np.sqrt(N) * 1.2

    # Histogram
    bins = 70
    counts, edges = np.histogram(data, bins=bins, range=(-x_lim, x_lim), density=True)
    centers = 0.5 * (edges[1:] + edges[:-1])
    ax.bar(centers, counts, width=(edges[1]-edges[0])*0.95,
           color="steelblue", alpha=0.65, edgecolor="navy", linewidth=0.3)

    # Yarı-daire eğrisi
    x_grid = np.linspace(-2*np.sqrt(N), 2*np.sqrt(N), 400)
    ax.plot(x_grid, semicircle(x_grid, N), "r-", lw=2.2,
            label="Wigner yarı-daire")

    ax.set_title(f"N = {N} kişilik asansör", fontsize=12, fontweight="bold")
    ax.set_xlabel("oturma pozisyonu")
    if N == 3:
        ax.set_ylabel("yoğunluk")
    ax.set_xlim(-x_lim, x_lim)
    ax.legend(fontsize=8, loc="lower center")
    ax.grid(alpha=0.3)

plt.suptitle("16 — Asansör Büyüdükçe: özdeğer dağılımı YARI-DAİREYE yakınsıyor",
             fontsize=13, fontweight="bold", y=1.03)
plt.tight_layout()

out = Path(__file__).parent / "16_asansor_buyuyor.png"
plt.savefig(out, dpi=130, bbox_inches="tight")
print(f"\nKaydedildi: {out}")

# ============================================================
# Bonus: aynı eksende karşılaştırma (normalize)
# ============================================================
fig2, ax = plt.subplots(figsize=(11, 6))
colors = ["crimson", "darkorange", "gold", "seagreen", "steelblue"]
for N, c in zip(Ns, colors):
    data = datasets[N]
    # Normalize: özdeğerleri 2√N'ye böl → hepsi [-1, +1]'e gelsin
    x_norm = data / (2 * np.sqrt(N))
    ax.hist(x_norm, bins=80, density=True, alpha=0.45, color=c,
            label=f"N={N}", histtype="stepfilled")

# Teorik yarı-daire normalize: ρ(x) = (2/π)√(1 - x²)
x_grid = np.linspace(-1, 1, 400)
y_grid = (2/np.pi) * np.sqrt(np.maximum(1 - x_grid**2, 0))
ax.plot(x_grid, y_grid, "k-", lw=2.5, label="Wigner yarı-daire (teori)")

ax.set_xlabel("normalize pozisyon  (özdeğer / 2√N)")
ax.set_ylabel("yoğunluk")
ax.set_title("Aynı eksende: tüm asansör boyutları YARIM DAİREYE yakınsıyor")
ax.legend(fontsize=10)
ax.grid(alpha=0.3)
ax.set_xlim(-1.3, 1.3)

plt.tight_layout()
out2 = Path(__file__).parent / "16_asansor_yakinsama.png"
plt.savefig(out2, dpi=130, bbox_inches="tight")
print(f"Kaydedildi: {out2}")

# ============================================================
# Sayısal: yarı-daireye yakınlık
# ============================================================
print(f"\n=== YARI-DAİRE YAKINLIĞI (KS-tipi metrik) ===")
for N in Ns:
    data = datasets[N] / (2 * np.sqrt(N))
    # data ∈ [-1, 1] içinde olsun
    data = data[np.abs(data) <= 1]
    # Teorik CDF: ∫_{-1}^{x} (2/π)√(1-t²) dt = (1/2) + (1/π)(x√(1-x²) + arcsin(x))
    sorted_d = np.sort(data)
    empirical_cdf = np.arange(1, len(sorted_d)+1) / len(sorted_d)
    theoretical_cdf = 0.5 + (sorted_d * np.sqrt(1 - sorted_d**2) + np.arcsin(sorted_d)) / np.pi
    ks_dist = np.max(np.abs(empirical_cdf - theoretical_cdf))
    print(f"  N={N:>4}:  KS mesafesi = {ks_dist:.4f}  (küçük = yarı-daire yakın)")
