"""
QM × Riemann — Sierra Genişletmesi Sayısal Diagonalizasyon
=============================================================

İki Hamiltonian'ı yan yana diagonalize et:

  (A) NAİF Berry-Keating:  H = (xp + px)/2  → log koord. -i(∂_u + 1/2)
      Periyodik BC ile özdeğer: E_n = 2πn/L + 1/2  ← EŞIT ARALIKLI (yanlış)

  (B) SIERRA fix:          H = x(p + ℓ²/p)
      Sınırlı yörüngeler. log-yoğunluklu spektrum çıkmalı.
      Yarı-klasik N(E) ≈ (E/2π) log(E/(2πe)) + 7/8 ← Riemann'a yapışık

Özdeğerleri ilk 30 Riemann γ_n ile karşılaştır.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.linalg import eigh

ROOT = Path(__file__).resolve().parent.parent
gamma_riemann = np.loadtxt(ROOT / "zeros_100k.txt")

# ============================================================
# (A) NAİF Berry-Keating — log latis, periyodik BC
# ============================================================
# u ∈ [0, U],  ψ periyodik → -i ∂_u'nin özdeğerleri 2πn/U
# +1/2 ofseti dahil

U_box = 30.0       # box uzunluğu
N = 800            # latis nokta sayısı
du = U_box / N
u = np.arange(N) * du

# -i ∂_u: anti-Hermitian first-order — periyodik BC, merkezi fark
# (D ψ)_k = (ψ_{k+1} - ψ_{k-1}) / (2 du)
D = np.zeros((N, N), complex)
for k in range(N):
    D[k, (k+1) % N] += 1.0 / (2*du)
    D[k, (k-1) % N] -= 1.0 / (2*du)
H_naive = -1j * D + 0.5 * np.eye(N)
H_naive = 0.5 * (H_naive + H_naive.conj().T)   # Hermitize (sayısal yuvarlama)

eigs_naive = np.sort(np.linalg.eigvalsh(H_naive))
# pozitif olanlar
eigs_naive_pos = eigs_naive[eigs_naive > 0][:30]

# ============================================================
# (B) SIERRA  H = x(p + ℓ²/p)
# ============================================================
# x ∈ [1, X_max] log-latis: x_k = exp(u_k), u_k = k*du, u ∈ [0, log(X_max)]
# Konum bazında p = -i d/dx,  1/p formal — Fourier üzerinden çalışmak gerekir
#
# Pratik: doğrudan u-temsilinde çalışalım.
#   xp + px = -i (1 + 2 x ∂_x) = -i (1 + 2 ∂_u)
#   x · (ℓ²/p) terimi ise: -i ℓ² · (1/p) · x → bunu spectral metodla yapacağız
#
# Daha temiz yol: x-temsilinde matris kur, p = -i d/dx central difference,
# 1/p'yi pseudoinverse ile yaklaş. SVD ile çözer.

X_min, X_max = 1.0, 200.0
N2 = 600
x = np.linspace(X_min, X_max, N2)
dx = x[1] - x[0]
ell2 = 1.0   # ℏ=1, ℓ² parametresi

# p = -i d/dx (Dirichlet BC: ψ=0 dışta), merkezi fark
P = np.zeros((N2, N2), complex)
for k in range(N2):
    if k+1 < N2: P[k, k+1] += 1.0/(2*dx)
    if k-1 >= 0: P[k, k-1] -= 1.0/(2*dx)
P = -1j * P
P = 0.5 * (P + P.conj().T)   # Hermitize

X = np.diag(x.astype(complex))

# 1/p: spektral hesapla (P'nin sıfır olmayan özdeğerlerinin tersi)
eigP, UP = np.linalg.eigh(P)
# sıfırdan uzak özdeğerler için 1/E, çok küçük olanları kes (regularize)
threshold = 0.01 * np.max(np.abs(eigP))
inv_eig = np.where(np.abs(eigP) > threshold, 1.0/eigP, 0.0)
P_inv = UP @ np.diag(inv_eig) @ UP.conj().T

# H_sierra = x · (p + ℓ²/p)
H_sierra = X @ (P + ell2 * P_inv)
# Hermitize (simetrik ordering): (xp + px)/2 + ℓ²(x/p + 1/p · x)/2
xp_sym = 0.5 * (X @ P + P @ X)
inv_sym = 0.5 * (X @ P_inv + P_inv @ X)
H_sierra = xp_sym + ell2 * inv_sym
H_sierra = 0.5 * (H_sierra + H_sierra.conj().T)

eigs_sierra = np.sort(np.linalg.eigvalsh(H_sierra))
eigs_sierra_pos = eigs_sierra[eigs_sierra > 1.0][:30]   # ilk 30 pozitif

# ============================================================
# Riemann ilk 30 sıfır
# ============================================================
gamma30 = gamma_riemann[:30]

# ============================================================
# Plot
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(13, 6))

# Sol: 1. naif BK vs Riemann
ax = axes[0]
ax.plot(np.arange(1, 31), gamma30, "ro-", lw=1.5, ms=6, label="Riemann γ_n (ilk 30)")
ax.plot(np.arange(1, 31), eigs_naive_pos, "bs--", lw=1.0, ms=5, alpha=0.7,
        label="NAİF BK: -i(∂_u+½), periyodik")
ax.set_xlabel("n")
ax.set_ylabel("E_n  /  γ_n")
ax.set_title("Naif Berry-Keating: EŞIT ARALIKLI spektrum\n(Riemann ile yapısal uyuşmazlık)")
ax.legend(loc="upper left", fontsize=9)
ax.grid(alpha=0.3)

# Sağ: 2. Sierra fix vs Riemann
ax = axes[1]
ax.plot(np.arange(1, len(gamma30)+1), gamma30, "ro-", lw=1.5, ms=6, label="Riemann γ_n")
ax.plot(np.arange(1, len(eigs_sierra_pos)+1), eigs_sierra_pos, "g^-",
        lw=1.0, ms=5, alpha=0.85, label="SIERRA H=x(p+ℓ²/p) sayısal")
# Berry-Keating yarı-klasik öngörü (smooth)
n_grid = np.arange(1, 31)
# E_n smooth: N(E) = n ⇒ inverse. (E/2π) log(E/2πe) + 7/8 ≈ n  → Newton
def E_of_n(n):
    """Yarı-klasik N(E)^{-1}"""
    E = 2*np.pi * n / np.log(n+1) * 0.5   # ilk tahmin
    for _ in range(40):
        F = (E/(2*np.pi))*np.log(E/(2*np.pi*np.e)) + 7/8 - n
        Fp = np.log(E/(2*np.pi*np.e))/(2*np.pi) + 1/(2*np.pi)
        E = E - F/Fp
    return E
E_semi = np.array([E_of_n(n) for n in n_grid])
ax.plot(n_grid, E_semi, "k:", lw=1.4, alpha=0.7, label="Yarı-klasik N⁻¹ (BK öngörü)")

ax.set_xlabel("n")
ax.set_ylabel("E_n  /  γ_n")
ax.set_title("Sierra genişlemesi vs Riemann\nyarı-klasik kıvrım yakın, tek tek sıfırlar yok")
ax.legend(loc="upper left", fontsize=9)
ax.grid(alpha=0.3)

plt.tight_layout()
out = Path(__file__).parent / "05_sierra_diagonalize.png"
plt.savefig(out, dpi=130, bbox_inches="tight")
print(f"Kaydedildi: {out}")

# ============================================================
# Sayısal rapor
# ============================================================
print(f"\n{'n':>3} | {'γ_n (Riemann)':>15} | {'Naif BK':>10} | {'Sierra':>10} | {'Yarı-klasik':>12}")
print("-"*65)
for i in range(15):
    print(f"{i+1:>3} | {gamma30[i]:>15.4f} | {eigs_naive_pos[i]:>10.4f} | "
          f"{eigs_sierra_pos[i]:>10.4f} | {E_semi[i]:>12.4f}")

print(f"\nNaif BK aralık (E_2 - E_1) = {eigs_naive_pos[1]-eigs_naive_pos[0]:.4f}  (sabit, yanlış)")
print(f"Riemann (γ_2 - γ_1)       = {gamma30[1]-gamma30[0]:.4f}  (değişken, gerçek)")
print(f"Sierra (E_2 - E_1)        = {eigs_sierra_pos[1]-eigs_sierra_pos[0]:.4f}")
