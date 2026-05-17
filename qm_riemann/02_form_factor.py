"""
QM × Riemann — Form Factor K(τ)
================================

Üçüncü yüz: pair correlation'ın Fourier eşi.

GUE öngörüsü (Dyson):
    K(τ) = τ           for 0 ≤ τ ≤ 1   ← DOĞRUSAL YÜKSELİM
    K(τ) = 1           for τ ≥ 1       ← doyum

Doğrusal yükselim = "kuantum kaotik" imzası (level repulsion'ın Fourier hali).
İntegre sistemlerde K(τ) = sabit (Poisson).
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
t = np.loadtxt(ROOT / "zeros_100k.txt")

# Unfold
rho = np.log(t / (2 * np.pi)) / (2 * np.pi)
xi = t * rho
print(f"N = {len(xi)}, mean spacing (unfold) = {np.mean(np.diff(xi)):.5f}")

# ---------- Form factor (window-averaged) ----------
# K(τ) = <(1/M) |Σ_n e^{2πi τ ξ_n}|²>_windows
# pratik: penceredeki ξ'leri merkez etrafında kaydır ki düşük τ'da sayısal sorun olmasın
WINDOW = 1000              # her pencerede 1000 sıfır
N_WIN  = 80                # 80 farklı pencere ortala
tau    = np.linspace(0.005, 2.0, 220)

rng = np.random.default_rng(0)
starts = rng.integers(0, len(xi) - WINDOW - 1, size=N_WIN)

K_emp = np.zeros_like(tau)
for s in starts:
    seg = xi[s:s+WINDOW]
    seg = seg - seg.mean()                       # merkezle
    # K(τ) için kullanılan sum
    # |Σ e^{2πi τ ξ}|² / M
    phase = np.exp(2j * np.pi * np.outer(tau, seg))  # (Ntau, M)
    contrib = np.abs(phase.sum(axis=1))**2 / WINDOW
    K_emp += contrib
K_emp /= N_WIN

# GUE öngörüsü
K_gue = np.where(tau < 1, tau, 1.0)
# Poisson (integre / bağımsız) için K = sabit = 1 (her zaman)
K_poi = np.ones_like(tau)

# ---------- Plot ----------
fig, ax = plt.subplots(figsize=(8.5, 5.5))

ax.plot(tau, K_emp, "o-", ms=3, lw=1.0, color="steelblue", alpha=0.85,
        label="ζ sıfırları (empirik)")
ax.plot(tau, K_gue, "r-", lw=2.4, label="GUE: K(τ)=min(τ,1)")
ax.plot(tau, K_poi, "k:", lw=1.4, alpha=0.7, label="Poisson: K(τ)=1")

ax.axvline(1, color="grey", lw=0.7, alpha=0.5)
ax.text(1.02, 0.05, "τ=1\n(Heisenberg time)", fontsize=8, color="grey")

ax.set_xlabel("τ  (zaman / spacing⁻¹ ölçeği)")
ax.set_ylabel("K(τ)")
ax.set_title("Form Factor — ζ sıfırları vs GUE\n"
             "doğrusal yükselim = kuantum kaotik imza")
ax.set_xlim(0, 2)
ax.set_ylim(-0.05, 1.4)
ax.grid(alpha=0.3)
ax.legend(loc="lower right", fontsize=10)

out = Path(__file__).parent / "02_form_factor.png"
plt.tight_layout()
plt.savefig(out, dpi=130, bbox_inches="tight")
print(f"Kaydedildi: {out}")

# ---------- Sayısal rapor ----------
# Küçük τ bölgesinde doğru eğim?
mask_low = (tau > 0.05) & (tau < 0.5)
slope, intercept = np.polyfit(tau[mask_low], K_emp[mask_low], 1)
print(f"\nKüçük τ (0.05<τ<0.5) doğrusal fit:")
print(f"  eğim     = {slope:.4f}   (GUE: 1.0000)")
print(f"  intercept = {intercept:.4f}   (GUE: 0.0000)")

mask_high = tau > 1.2
print(f"\nBüyük τ (τ>1.2) ortalama:")
print(f"  K_emp    = {K_emp[mask_high].mean():.4f}   (GUE: 1.0000)")
