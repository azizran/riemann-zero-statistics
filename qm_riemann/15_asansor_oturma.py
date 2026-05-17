"""
15 — GUE Asansöründe Oturma Düzeni
=====================================

3 kişilik ve 5 kişilik asansörler (N=3 ve N=5 GUE matrisleri).
10000 deneme yap, her seferinde özdeğerleri sırala, kim nereye oturdu bak.

Sol panel:  her özdeğerin (1., 2., 3., ...) ayrı ayrı dağılımı
Sağ panel:  100 ayrı denemenin "oturma fotoğrafı" (her satır bir deneme)

Bonus: ardışık özdeğerler arası mesafe (gap) dağılımı — Wigner surmise.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

rng = np.random.default_rng(42)

def gue_matrix(N):
    """N×N GUE matris üret"""
    A = rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))
    return (A + A.conj().T) / 2

def sorted_eigs(N, N_trials):
    """N_trials deneme, her birinin N özdeğerini sıralı şekilde döndür"""
    arr = np.zeros((N_trials, N))
    for i in range(N_trials):
        H = gue_matrix(N)
        eigs = np.sort(np.linalg.eigvalsh(H))
        arr[i] = eigs
    return arr

# ============================================================
# Verileri üret
# ============================================================
N_trials = 10000
eigs_3 = sorted_eigs(3, N_trials)
eigs_5 = sorted_eigs(5, N_trials)

print(f"N=3 asansör, {N_trials} deneme:")
for k in range(3):
    mean = eigs_3[:, k].mean()
    std = eigs_3[:, k].std()
    print(f"  Kişi {k+1}: ortalama pozisyon = {mean:+.3f},  yayılım = ±{std:.3f}")

print(f"\nN=5 asansör, {N_trials} deneme:")
for k in range(5):
    mean = eigs_5[:, k].mean()
    std = eigs_5[:, k].std()
    print(f"  Kişi {k+1}: ortalama pozisyon = {mean:+.3f},  yayılım = ±{std:.3f}")

# ============================================================
# Plot
# ============================================================
fig, axes = plt.subplots(3, 2, figsize=(15, 13))

# Renkler
colors_3 = ["crimson", "darkorange", "steelblue"]
colors_5 = ["crimson", "darkorange", "gold", "seagreen", "steelblue"]

# === ÜST SOL — N=3 oturma dağılımı ===
ax = axes[0, 0]
x = np.linspace(-4, 4, 200)
for k in range(3):
    ax.hist(eigs_3[:, k], bins=80, density=True, alpha=0.55,
            color=colors_3[k], edgecolor=colors_3[k],
            label=f"{k+1}. kişi (sıralı)")
ax.set_xlabel("pozisyon (özdeğer)")
ax.set_ylabel("yoğunluk")
ax.set_title("3 KİŞİLİK ASANSÖR — her kişinin oturma dağılımı")
ax.legend(fontsize=10)
ax.grid(alpha=0.3)
ax.set_xlim(-4, 4)

# === ÜST SAĞ — 100 deneme oturma fotoğrafı ===
ax = axes[0, 1]
N_show = 80
for i in range(N_show):
    for k in range(3):
        ax.plot(eigs_3[i, k], i, "o", ms=8, color=colors_3[k], alpha=0.8)
ax.set_xlim(-4, 4)
ax.set_ylim(-1, N_show)
ax.set_xlabel("pozisyon")
ax.set_ylabel("deneme numarası")
ax.set_title(f"N=3: {N_show} farklı denemede insanların oturma yerleri\n"
             "(her satır bir asansör fotoğrafı)")
ax.grid(alpha=0.3)
for k, c in enumerate(colors_3):
    ax.axvline(eigs_3[:, k].mean(), color=c, lw=0.5, ls="--", alpha=0.6)

# === ORTA SOL — N=5 oturma dağılımı ===
ax = axes[1, 0]
for k in range(5):
    ax.hist(eigs_5[:, k], bins=80, density=True, alpha=0.5,
            color=colors_5[k], edgecolor=colors_5[k],
            label=f"{k+1}. kişi")
ax.set_xlabel("pozisyon")
ax.set_ylabel("yoğunluk")
ax.set_title("5 KİŞİLİK ASANSÖR — her kişinin oturma dağılımı")
ax.legend(fontsize=9, ncol=2)
ax.grid(alpha=0.3)
ax.set_xlim(-5, 5)

# === ORTA SAĞ — N=5 fotoğraf ===
ax = axes[1, 1]
for i in range(N_show):
    for k in range(5):
        ax.plot(eigs_5[i, k], i, "o", ms=7, color=colors_5[k], alpha=0.8)
ax.set_xlim(-5, 5)
ax.set_ylim(-1, N_show)
ax.set_xlabel("pozisyon")
ax.set_ylabel("deneme numarası")
ax.set_title(f"N=5: {N_show} farklı denemede 5 kişinin oturma yerleri")
ax.grid(alpha=0.3)
for k, c in enumerate(colors_5):
    ax.axvline(eigs_5[:, k].mean(), color=c, lw=0.5, ls="--", alpha=0.6)

# === ALT SOL — Komşu kişiler arası mesafe (N=3) ===
ax = axes[2, 0]
gaps_3_12 = eigs_3[:, 1] - eigs_3[:, 0]    # 1 ile 2 arası
gaps_3_23 = eigs_3[:, 2] - eigs_3[:, 1]    # 2 ile 3 arası
ax.hist(gaps_3_12, bins=80, density=True, alpha=0.55,
        color="darkmagenta", label="kişi 1-2 arası mesafe")
ax.hist(gaps_3_23, bins=80, density=True, alpha=0.55,
        color="teal", label="kişi 2-3 arası mesafe")
ax.set_xlabel("mesafe (boşluk)")
ax.set_ylabel("yoğunluk")
ax.set_title("N=3: komşular arası mesafe — 0'dan uzak duruyorlar (level repulsion)")
ax.legend(fontsize=10)
ax.grid(alpha=0.3)
ax.set_xlim(0, 6)

# === ALT SAĞ — N=5 komşular arası mesafe ===
ax = axes[2, 1]
for k in range(4):
    g = eigs_5[:, k+1] - eigs_5[:, k]
    ax.hist(g, bins=80, density=True, alpha=0.4, label=f"kişi {k+1}-{k+2}")
ax.set_xlabel("mesafe")
ax.set_ylabel("yoğunluk")
ax.set_title("N=5: 4 farklı komşu çifti arası mesafe — hiçbiri 0'da pik yapmıyor")
ax.legend(fontsize=9, ncol=2)
ax.grid(alpha=0.3)
ax.set_xlim(0, 6)

plt.suptitle("15 — GUE Asansörünün İçi: kim nerede oturuyor, hangi sırayla, hangi mesafelerle",
             fontsize=13, fontweight="bold", y=0.998)
plt.tight_layout()

from pathlib import Path
out = Path(__file__).parent / "15_asansor_oturma.png"
plt.savefig(out, dpi=130, bbox_inches="tight")
print(f"\nKaydedildi: {out}")

# ============================================================
# Sayısal: en küçük mesafe ne kadar?
# ============================================================
print(f"\n=== ARALIK ANALİZİ ===")
print(f"N=3: kişi 1-2 arası mesafe   ortalama={gaps_3_12.mean():.3f},  min={gaps_3_12.min():.3f}")
print(f"N=3: kişi 2-3 arası mesafe   ortalama={gaps_3_23.mean():.3f},  min={gaps_3_23.min():.3f}")
print(f"\n%5 olasılıkla en yakın mesafe (N=3, 1-2): {np.percentile(gaps_3_12, 5):.4f}")
print(f"Yani: 100 denemenin 5'inde komşular bu mesafeden DAHA yakın oturuyor")
print(f"Asla yapışmıyorlar — 'level repulsion' canlı veride")
