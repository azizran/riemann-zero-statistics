"""
QM × Riemann — Üç Dünya Yan Yana
==================================

Bohigas-Giannoni-Schmit konjektürü: kuantum kaotik bir sistemin enerji
seviyeleri büyük ölçüde rastgele Hermitian matrislerin özdeğerleri gibidir.

Bu dosya üç gerçek nesneyi yan yana koyar:

  (1) ζ Riemann sıfırları              — sayı teorisi
  (2) GUE matris özdeğerleri (zaman tersi yok)   — kuantum kaotik
  (3) GOE matris özdeğerleri (zaman tersi simetrik) — uranyum çekirdeği

(2) ile ζ üst üste binmeli (Montgomery-Dyson). (3) farklı durmalı.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
t = np.loadtxt(ROOT / "zeros_100k.txt")

# ============================================================
# (1) ζ unfolded spacings (mevcut yöntem)
# ============================================================
rho = np.log(t / (2 * np.pi)) / (2 * np.pi)
xi = t * rho
s_zeta = np.diff(xi)
s_zeta = s_zeta / s_zeta.mean()

# ============================================================
# (2) GUE matris — Gaussian Unitary Ensemble
# ============================================================
# H = (A + A†)/√2  with A_ij ~ CN(0, 1/N)
# özdeğerler Wigner semicircle dağılımı
N_mat = 2500
N_real = 20    # 20 farklı realizasyon, istatistik için
rng = np.random.default_rng(7)

s_gue_all = []
for _ in range(N_real):
    A = (rng.standard_normal((N_mat, N_mat)) +
         1j * rng.standard_normal((N_mat, N_mat))) / np.sqrt(2*N_mat)
    H = (A + A.conj().T) / np.sqrt(2)
    eigs = np.sort(np.linalg.eigvalsh(H))
    # bulk'a kalsın (kenarlardaki seviye yoğunluğu eşit değil)
    eigs = eigs[int(0.2*N_mat):int(0.8*N_mat)]
    # unfold: lokal yoğunluk = semicircle ρ(E) = (1/π)√(2-E²)
    rho_sc = lambda E: np.sqrt(np.maximum(0, 2 - E**2)) / np.pi
    # kümülatif sayım: N(E) = ∫_{-√2}^E ρ ≈ analytic form
    def N_semicircle(E, N):
        # bulk normalization: integral over [-√2,√2] = 1
        # local count = N · ∫_{-√2}^E ρ(E')dE'
        from numpy import arcsin, sqrt
        return N * (E*sqrt(2-E**2)/(2*np.pi) + arcsin(E/np.sqrt(2))/np.pi + 0.5)
    xi_gue = N_semicircle(eigs, N_mat)
    sp = np.diff(xi_gue)
    sp = sp[sp > 0]
    s_gue_all.append(sp / sp.mean())
s_gue = np.concatenate(s_gue_all)

# ============================================================
# (3) GOE matris — Gaussian Orthogonal Ensemble (gerçek simetrik)
# ============================================================
s_goe_all = []
for _ in range(N_real):
    A = rng.standard_normal((N_mat, N_mat)) / np.sqrt(N_mat)
    H = (A + A.T) / np.sqrt(2)
    eigs = np.sort(np.linalg.eigvalsh(H))
    eigs = eigs[int(0.2*N_mat):int(0.8*N_mat)]
    xi_goe = N_semicircle(eigs, N_mat)
    sp = np.diff(xi_goe)
    sp = sp[sp > 0]
    s_goe_all.append(sp / sp.mean())
s_goe = np.concatenate(s_goe_all)

print(f"ζ spacing örneklem: {len(s_zeta)}")
print(f"GUE spacing örneklem: {len(s_gue)}")
print(f"GOE spacing örneklem: {len(s_goe)}")

# ============================================================
# Teorik öngörüler
# ============================================================
def wigner_gue(s):
    return (32/np.pi**2) * s**2 * np.exp(-4*s**2/np.pi)

def wigner_goe(s):
    return (np.pi/2) * s * np.exp(-np.pi*s**2/4)

# ============================================================
# Plot — 3 panel + birleşik son panel
# ============================================================
fig, axes = plt.subplots(2, 2, figsize=(13, 10))
s_grid = np.linspace(0, 3.5, 400)

# Üst-sol: ζ sıfırları
ax = axes[0, 0]
ax.hist(s_zeta, bins=80, density=True, alpha=0.6, color="steelblue",
        edgecolor="navy", label=f"ζ sıfırları (N={len(s_zeta)})")
ax.plot(s_grid, wigner_gue(s_grid), "r-", lw=2.2, label="GUE teorik")
ax.plot(s_grid, wigner_goe(s_grid), "g--", lw=1.2, alpha=0.6, label="GOE teorik")
ax.set_xlim(0, 3.5)
ax.set_xlabel("s")
ax.set_ylabel("P(s)")
ax.set_title("(1) ζ sıfırları — sayı teorisi nesnesi")
ax.legend(fontsize=9)
ax.grid(alpha=0.3)

# Üst-sağ: GUE matris
ax = axes[0, 1]
ax.hist(s_gue, bins=80, density=True, alpha=0.6, color="crimson",
        edgecolor="darkred", label=f"GUE matris (N={N_mat}, {N_real} realizasyon)")
ax.plot(s_grid, wigner_gue(s_grid), "r-", lw=2.2, label="GUE teorik")
ax.plot(s_grid, wigner_goe(s_grid), "g--", lw=1.2, alpha=0.6, label="GOE teorik")
ax.set_xlim(0, 3.5)
ax.set_xlabel("s")
ax.set_ylabel("P(s)")
ax.set_title("(2) GUE rastgele matris\n(zaman tersi kırık QM kaotik sistem)")
ax.legend(fontsize=9)
ax.grid(alpha=0.3)

# Alt-sol: GOE matris (uranyum çekirdeği gibi)
ax = axes[1, 0]
ax.hist(s_goe, bins=80, density=True, alpha=0.6, color="seagreen",
        edgecolor="darkgreen", label=f"GOE matris (N={N_mat})")
ax.plot(s_grid, wigner_gue(s_grid), "r--", lw=1.2, alpha=0.6, label="GUE teorik")
ax.plot(s_grid, wigner_goe(s_grid), "g-", lw=2.2, label="GOE teorik")
ax.set_xlim(0, 3.5)
ax.set_xlabel("s")
ax.set_ylabel("P(s)")
ax.set_title("(3) GOE rastgele matris\n(uranyum çekirdeği, kaotik biliardo)")
ax.legend(fontsize=9)
ax.grid(alpha=0.3)

# Alt-sağ: ÜÇÜ ÜST ÜSTE
ax = axes[1, 1]
ax.hist(s_zeta, bins=80, density=True, alpha=0.4, color="steelblue",
        label=f"ζ sıfırları")
ax.hist(s_gue, bins=80, density=True, alpha=0.4, color="crimson",
        label=f"GUE matris")
ax.hist(s_goe, bins=80, density=True, alpha=0.4, color="seagreen",
        label=f"GOE matris")
ax.plot(s_grid, wigner_gue(s_grid), "r-", lw=2, label="GUE teorik")
ax.plot(s_grid, wigner_goe(s_grid), "g--", lw=1.5, label="GOE teorik")
ax.set_xlim(0, 3.5)
ax.set_xlabel("s")
ax.set_ylabel("P(s)")
ax.set_title("Üç dünya üst üste\nζ + GUE matris çakışıyor, GOE ayrı duruyor")
ax.legend(fontsize=8, loc="upper right")
ax.grid(alpha=0.3)

plt.suptitle("Bohigas-Giannoni-Schmit: ζ ≈ GUE rastgele matris ≠ GOE",
             fontsize=13, y=1.00)
plt.tight_layout()
out = Path(__file__).parent / "07_uc_dunya.png"
plt.savefig(out, dpi=130, bbox_inches="tight")
print(f"Kaydedildi: {out}")

# ============================================================
# Sayısal: hangi çiftler ne kadar yakın?
# ============================================================
def hist_dist(s1, s2, bins=100):
    h1, edges = np.histogram(s1, bins=bins, range=(0, 4), density=True)
    h2, _ = np.histogram(s2, bins=bins, range=(0, 4), density=True)
    return np.mean(np.abs(h1 - h2))

print(f"\nİki örneklem arası L1 mesafesi (küçük = aynı dağılım):")
print(f"  ζ      vs GUE matris : {hist_dist(s_zeta, s_gue):.4f}")
print(f"  ζ      vs GOE matris : {hist_dist(s_zeta, s_goe):.4f}")
print(f"  GUE    vs GOE matris : {hist_dist(s_gue, s_goe):.4f}")

# Teorik öngörüye uyum
def fit_dist(samp, pred_fn, bins=100):
    h, edges = np.histogram(samp, bins=bins, range=(0, 4), density=True)
    centers = 0.5*(edges[1:]+edges[:-1])
    return np.mean(np.abs(h - pred_fn(centers)))

print(f"\nL1 hata teorik öngörüye karşı:")
print(f"  ζ    vs GUE Wigner : {fit_dist(s_zeta, wigner_gue):.4f}")
print(f"  ζ    vs GOE Wigner : {fit_dist(s_zeta, wigner_goe):.4f}")
print(f"  GUE  vs GUE Wigner : {fit_dist(s_gue,  wigner_gue):.4f}")
print(f"  GOE  vs GOE Wigner : {fit_dist(s_goe,  wigner_goe):.4f}")
