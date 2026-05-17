"""
QM × Riemann — ilk resim
========================

Soru: Riemann ζ sıfırlarının istatistiği gerçekten GUE (kuantum kaotik
sistemlerin) istatistiği ile örtüşüyor mu? Gözle göster.

İki test:
  (1) Nearest-neighbor spacing dağılımı vs GUE Wigner surmise
  (2) Pair correlation R_2(r) vs 1 - sinc^2(πr)

Veri: zeros_100k.txt (Odlyzko, ilk 100000 ζ sıfırının imajiner kısmı)
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
zeros_path = ROOT / "zeros_100k.txt"

# ---------- 1. Veriyi yükle ----------
t = np.loadtxt(zeros_path)
print(f"Yüklenen sıfır sayısı: {len(t)}")
print(f"İlk: {t[0]:.6f}, son: {t[-1]:.6f}")

# ---------- 2. Lokal ortalama yoğunlukla normalize et ----------
# Riemann-von Mangoldt: N(T) ~ (T/2π) log(T/2π) - T/2π
# Yerel sıfır yoğunluğu:  ρ(t) = (1/2π) log(t/2π)
# Normalize edilmiş sıfırlar:  ξ_n = t_n * ρ(t_n) — ardışık fark ~1
rho = np.log(t / (2 * np.pi)) / (2 * np.pi)
xi = t * rho  # unfolded
spacings = np.diff(xi)
spacings = spacings / np.mean(spacings)  # ortalama 1'e ayarla

print(f"Spacing istatistik: mean={spacings.mean():.4f}, std={spacings.std():.4f}")

# ---------- 3. GUE Wigner surmise ----------
def gue_wigner(s):
    return (32 / np.pi**2) * s**2 * np.exp(-4 * s**2 / np.pi)

def goe_wigner(s):
    return (np.pi / 2) * s * np.exp(-np.pi * s**2 / 4)

def poisson(s):
    return np.exp(-s)

# ---------- 4. Pair correlation (Montgomery) ----------
# R_2(r) ≈ 1 - (sin(πr)/(πr))^2 (GUE prediction)
def gue_pair_corr(r):
    pi_r = np.pi * r
    return 1 - np.where(np.abs(pi_r) < 1e-9, 1.0, (np.sin(pi_r) / pi_r) ** 2)

# Pair correlation'ı sınırlı pencerede ölç (N=20000 sıfır yeterli, hızlı)
N_pc = 20000
xi_sub = xi[:N_pc]
# tüm çiftlerin farkları (mod yerel ortalama) — sadece pozitif farkları al
all_diffs = []
window = 50  # ardışık ±50 sıfır içinde
for i in range(N_pc):
    j_max = min(i + window, N_pc)
    diffs = xi_sub[i+1:j_max] - xi_sub[i]
    all_diffs.append(diffs)
all_diffs = np.concatenate(all_diffs)
# normalize: ortalama spacing = 1 olacak şekilde zaten unfolded

# ---------- 5. Plot ----------
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Sol: nearest-neighbor spacing
ax = axes[0]
ax.hist(spacings, bins=80, density=True, alpha=0.55,
        color="steelblue", edgecolor="navy", label=f"ζ sıfırları (N={len(spacings)})")
s_grid = np.linspace(0, 3.5, 400)
ax.plot(s_grid, gue_wigner(s_grid), "r-", lw=2.2, label="GUE (Wigner)")
ax.plot(s_grid, goe_wigner(s_grid), "g--", lw=1.4, alpha=0.7, label="GOE")
ax.plot(s_grid, poisson(s_grid), "k:", lw=1.4, alpha=0.7, label="Poisson (bağımsız)")
ax.set_xlabel("normalize edilmiş ardışık fark  s")
ax.set_ylabel("P(s)")
ax.set_title("Nearest-Neighbor Spacing")
ax.set_xlim(0, 3.5)
ax.legend(loc="upper right", fontsize=9)
ax.grid(alpha=0.3)

# Sağ: pair correlation
ax = axes[1]
counts, edges = np.histogram(all_diffs, bins=200, range=(0, 3.5), density=True)
centers = 0.5 * (edges[1:] + edges[:-1])
# normalize: uzun mesafede R_2 → 1
# basit normalizasyon: histogramı tail ortalamasına böl
tail = counts[(centers > 2.5) & (centers < 3.5)].mean()
R2_emp = counts / tail
ax.plot(centers, R2_emp, "o", ms=3, color="steelblue", alpha=0.7,
        label="ζ sıfırları (empirik)")
r_grid = np.linspace(0.001, 3.5, 600)
ax.plot(r_grid, gue_pair_corr(r_grid), "r-", lw=2.2,
        label=r"GUE: $1-\mathrm{sinc}^2(\pi r)$")
ax.axhline(1, color="k", lw=0.7, alpha=0.5)
ax.set_xlabel("normalize edilmiş mesafe  r")
ax.set_ylabel(r"$R_2(r)$")
ax.set_title("Pair Correlation")
ax.set_xlim(0, 3.5)
ax.set_ylim(-0.05, 1.5)
ax.legend(loc="lower right", fontsize=9)
ax.grid(alpha=0.3)

plt.suptitle("Riemann ζ sıfırları vs GUE — kuantum kaotik istatistiğin imzası",
             fontsize=12, y=1.02)
plt.tight_layout()
out_path = Path(__file__).parent / "01_resim.png"
plt.savefig(out_path, dpi=130, bbox_inches="tight")
print(f"\nKaydedildi: {out_path}")

# ---------- 6. Sayısal karşılaştırma ----------
# Kolmogorov-Smirnov-tipi basit metrik: histogram mesafesi
hist, edges = np.histogram(spacings, bins=80, range=(0, 3.5), density=True)
centers = 0.5 * (edges[1:] + edges[:-1])
gue_pred = gue_wigner(centers)
goe_pred = goe_wigner(centers)
poi_pred = poisson(centers)
def L1(a, b): return np.mean(np.abs(a - b))
print(f"\nL1 mesafesi (küçük = iyi uyum):")
print(f"  ζ vs GUE:     {L1(hist, gue_pred):.4f}")
print(f"  ζ vs GOE:     {L1(hist, goe_pred):.4f}")
print(f"  ζ vs Poisson: {L1(hist, poi_pred):.4f}")
