"""
18 — Neden σ = 1/2? — ζ'nın σ-t Manzarası
=============================================

8. grafikte bir tek t değeri için σ taraması yapmıştık (V şekli).
Şimdi her şeyi birlikte görelim:

  - σ ∈ [0.05, 0.95], t ∈ [10, 50] grid
  - ζ(σ + it) hesapla, |ζ| yüzeyini renkle çiz
  - Sıfırlar (γ₁..γ₆) σ=1/2 çizgisinde kırmızı noktalar
  - Yan panel: σ-yönünde KUVVET = -∂|ζ|/∂σ  (sıfırı σ=1/2'ye iten kuvvet)

Bu görsel, "σ=1/2 niye?" sorusunun yarı-fiziksel cevabını verir:
σ ekseni üzerinde ζ'nın doğal kuyusu var, sıfır kuyunun dibine düşüyor.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import mpmath as mp

mp.mp.dps = 25

ROOT = Path(__file__).resolve().parent.parent
gamma = np.loadtxt(ROOT / "zeros_100k.txt")

# ============================================================
# σ-t gridinde ζ hesapla
# ============================================================
sigma_grid = np.linspace(0.05, 0.95, 80)
t_grid = np.linspace(10, 50, 200)

print(f"Grid: {len(sigma_grid)} × {len(t_grid)} = {len(sigma_grid)*len(t_grid)} nokta")
print(f"ζ hesaplanıyor (mpmath, ~30 saniye)...")

Z = np.zeros((len(sigma_grid), len(t_grid)))
for i, s in enumerate(sigma_grid):
    for j, tt in enumerate(t_grid):
        Z[i, j] = float(abs(mp.zeta(mp.mpc(s, tt))))
    if i % 20 == 0: print(f"  σ={s:.2f}")

# ============================================================
# σ yönünde KUVVET = -∂|ζ|/∂σ
# (sıfırın σ=1/2'ye doğru kayması için bu kuvvete sıfır gerekir)
# ============================================================
dZ_dsigma = -np.gradient(Z, sigma_grid, axis=0)

# Sıfırların bulunduğu t değerleri (bu pencerede)
gamma_window = gamma[(gamma > t_grid[0]) & (gamma < t_grid[-1])]

# ============================================================
# Plot
# ============================================================
fig = plt.figure(figsize=(15, 11))
gs = fig.add_gridspec(2, 2, height_ratios=[1.3, 1], width_ratios=[1.5, 1],
                     hspace=0.3, wspace=0.25)

# Üst sol: |ζ| yüzey
ax = fig.add_subplot(gs[0, 0])
extent = [t_grid[0], t_grid[-1], sigma_grid[0], sigma_grid[-1]]
im = ax.imshow(np.log10(Z + 1e-12), origin='lower', aspect='auto',
               extent=extent, cmap='magma_r')
plt.colorbar(im, ax=ax, label="log₁₀ |ζ(σ+it)|", fraction=0.04)
# Sıfırlar
for g in gamma_window:
    ax.plot(g, 0.5, "o", ms=10, color="cyan", mec="white", mew=1.5)
ax.axhline(0.5, color="cyan", lw=1, ls="--", alpha=0.6)
ax.set_xlabel("t")
ax.set_ylabel("σ")
ax.set_title(f"|ζ(σ+it)| manzarası — sıfırlar (cyan) σ=1/2'de\n"
             f"karanlık bölgeler = ζ küçük (kuyu), sıfırlar dipte oturuyor")

# Üst sağ: σ KESİTLERİ (her sıfır için ayrı V kuyusu)
ax = fig.add_subplot(gs[0, 1])
for g in gamma_window:
    j = np.argmin(np.abs(t_grid - g))
    ax.plot(sigma_grid, Z[:, j], lw=1.4, alpha=0.85,
            label=f"t = {g:.2f}")
ax.axvline(0.5, color="red", lw=1.5, ls="--", alpha=0.6, label="σ=1/2")
ax.set_xlabel("σ")
ax.set_ylabel("|ζ(σ + i·γ_n)|")
ax.set_title("Her sıfır γ_n için σ kesiti — hepsi σ=1/2'de V-kuyusu")
ax.legend(fontsize=8, loc="upper right")
ax.grid(alpha=0.3)
ax.set_ylim(-0.2, 4)

# Alt sol: σ yönünde KUVVET ALANI
ax = fig.add_subplot(gs[1, 0])
v = np.max(np.abs(dZ_dsigma)) * 0.5
im2 = ax.imshow(dZ_dsigma, origin='lower', aspect='auto',
                extent=extent, cmap='RdBu_r', vmin=-v, vmax=v)
plt.colorbar(im2, ax=ax, label="σ-yönünde kuvvet (-∂|ζ|/∂σ)", fraction=0.04)
for g in gamma_window:
    ax.plot(g, 0.5, "o", ms=8, color="black", mec="white", mew=1)
ax.axhline(0.5, color="black", lw=1, ls="--", alpha=0.5)
ax.set_xlabel("t")
ax.set_ylabel("σ")
ax.set_title("σ-yönünde KUVVET — kırmızı: σ büyütme isteği, mavi: σ küçültme isteği\n"
             "σ=1/2'de tam denge (sıfır)")

# Alt sağ: σ=1/2 çizgisi boyunca kuvvet
ax = fig.add_subplot(gs[1, 1])
sigma_idx_half = np.argmin(np.abs(sigma_grid - 0.5))
ax.plot(t_grid, dZ_dsigma[sigma_idx_half], color="purple", lw=1.5,
        label=f"σ=1/2 üzerinde kuvvet")
ax.axhline(0, color="black", lw=0.5)
for g in gamma_window:
    ax.axvline(g, color="cyan", lw=0.7, alpha=0.6)
ax.set_xlabel("t")
ax.set_ylabel("σ-yön kuvvet")
ax.set_title("σ=1/2 üzerinde kuvvet — sıfırlarda (cyan) sıfırlanıyor")
ax.legend(fontsize=9)
ax.grid(alpha=0.3)

plt.suptitle("18 — σ=1/2 NEDEN? — ζ'nın σ-t Manzarası ve Kuvvet Alanı",
             fontsize=13, fontweight="bold", y=0.995)

out = Path(__file__).parent / "18_yari_kuyu.png"
plt.savefig(out, dpi=130, bbox_inches="tight")
print(f"\nKaydedildi: {out}")

# ============================================================
# Sayısal: her sıfırda kuvvet ne kadar?
# ============================================================
print(f"\n=== σ=1/2'de KUVVET ÖLÇÜMÜ ===")
print(f"Her sıfır γ_n için, σ yönündeki etkili kuvvet (kuyu dibi olduğu için ≈ 0 olmalı):")
for g in gamma_window:
    j = np.argmin(np.abs(t_grid - g))
    F = dZ_dsigma[sigma_idx_half, j]
    Z_val = Z[sigma_idx_half, j]
    print(f"  t = {g:6.2f}:  |ζ| = {Z_val:.5f},  F_σ = {F:+.5f}")
