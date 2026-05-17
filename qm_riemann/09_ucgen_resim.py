"""
09 — 0 / 1/2 / 1 Üçgeni
=========================

Üç sistemde aynı yapı:
  - 0       : yasak / yokluk
  - 1       : ateşleyen / birim / başlangıç
  - 1/2     : ortanın dengesi, yığılma yeri

Üst sıra: üç sistem ayrı ayrı (çift yarık, Riemann, Sezen voting)
Alt sıra: ortak yapının tek-bakışta haritası

Bu Uğur'un dakikalar içinde indiği resmin görsel ifadesi.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle, FancyArrowPatch
from matplotlib.lines import Line2D
from pathlib import Path
import mpmath as mp

ROOT = Path(__file__).resolve().parent.parent
gamma = np.loadtxt(ROOT / "zeros_100k.txt")

fig = plt.figure(figsize=(15, 11))
gs = fig.add_gridspec(2, 3, height_ratios=[1.5, 1], hspace=0.35, wspace=0.25)

# ============================================================
# ÜST SOL — Çift Yarık
# ============================================================
ax = fig.add_subplot(gs[0, 0])
ax.set_xlim(-0.5, 6)
ax.set_ylim(-3, 3)
ax.set_aspect("equal")
ax.axis("off")
ax.set_title("Çift Yarık Deneyi", fontsize=13, fontweight="bold", pad=10)

# Tabanca (kaynak)
ax.add_patch(Rectangle((-0.3, -0.4), 0.7, 0.8, color="dimgrey"))
ax.text(0, -1.0, "kaynak\n(parçacık)", ha="center", fontsize=9, color="dimgrey")

# Yarık duvarı
ax.plot([2.5, 2.5], [-2.5, -0.6], color="black", lw=4)
ax.plot([2.5, 2.5], [-0.4, 0.4], color="black", lw=4)
ax.plot([2.5, 2.5], [0.6, 2.5], color="black", lw=4)
ax.text(2.5, -2.9, "iki yarık", ha="center", fontsize=9, color="black")

# Plaka
ax.plot([5.2, 5.2], [-2.5, 2.5], color="navy", lw=2)
# Girişim deseni: 5 parlak
for y, intensity in [(-2.0, 0.3), (-1.2, 0.6), (-0.4, 0.9),
                     (0.4, 0.9), (1.2, 0.6), (2.0, 0.3)]:
    ax.plot(5.2, y, "o", ms=8+15*intensity, color="gold", alpha=0.6+0.4*intensity)
ax.text(5.2, -2.9, "plaka\n(girişim deseni)", ha="center", fontsize=9, color="navy")

# Ok: kaynaktan yarığa
ax.annotate("", xy=(2.3, 0), xytext=(0.5, 0),
            arrowprops=dict(arrowstyle="->", color="dimgrey", lw=1.5))
# Yarıktan plakaya — iki ışın
for slit_y in [-0.5, 0.5]:
    for fringe_y in [-1.5, -0.5, 0.5, 1.5]:
        ax.plot([2.5, 5.2], [slit_y, fringe_y], color="gold", lw=0.5, alpha=0.25)

# Etiket: 0 ve 1
ax.text(2.5, 3.1, "Yarık ekseni →", ha="center", fontsize=8,
        style="italic", color="grey")

# ============================================================
# ÜST ORTA — Riemann
# ============================================================
ax = fig.add_subplot(gs[0, 1])
ax.set_xlim(-0.2, 1.2)
ax.set_ylim(0, 55)
ax.set_title("Riemann ζ Sıfırları", fontsize=13, fontweight="bold", pad=10)
ax.set_xlabel("σ (gerçek kısım)", fontsize=10)
ax.set_ylabel("t (yükseklik = sayma yönü)", fontsize=10)

# Kritik şerit (izinli bölge)
ax.add_patch(Rectangle((0, 0), 1, 55, facecolor="lightyellow", alpha=0.5))
# Yasak: σ ≤ 0 ve σ ≥ 1
ax.add_patch(Rectangle((-0.2, 0), 0.2, 55, facecolor="lightgrey", alpha=0.6))
ax.add_patch(Rectangle((1, 0), 0.3, 55, facecolor="lightgrey", alpha=0.6))

# Kritik çizgi σ=1/2
ax.axvline(0.5, color="crimson", lw=2.5, alpha=0.8, label="σ = 1/2 (kritik çizgi)")
# Sıfırlar
for g in gamma[gamma < 55]:
    ax.plot(0.5, g, "o", ms=10, color="crimson",
            mfc="red", mec="darkred", mew=1.2)
    ax.text(0.55, g, f"γ={g:.2f}", fontsize=7, va="center")

# Yasak duvarlar
ax.axvline(0, color="black", lw=2)
ax.axvline(1, color="black", lw=2)
ax.text(-0.13, 28, "σ=0\nyasak\n(ζ≠0)", ha="center", fontsize=8.5, color="black")
ax.text(1.13, 28, "σ=1\nateşleyen\n(birim)", ha="center", fontsize=8.5, color="darkblue")
ax.text(0.5, -3, "σ=1/2\nyığılma", ha="center", fontsize=9,
        fontweight="bold", color="crimson")

ax.set_xticks([0, 0.5, 1])
ax.set_xticklabels(["0", "½", "1"], fontsize=10)
ax.grid(alpha=0.2)

# ============================================================
# ÜST SAĞ — Sezen Voting
# ============================================================
ax = fig.add_subplot(gs[0, 2])
ax.set_xlim(-6, 6)
ax.set_ylim(-0.5, 5)
ax.set_title("Sezen Voting", fontsize=13, fontweight="bold", pad=10)
ax.set_xlabel("oy değeri", fontsize=10)
ax.set_ylabel("yığılma yoğunluğu", fontsize=10)

# Yasak bölge (-1, +1) gri
ax.add_patch(Rectangle((-1, 0), 2, 5, facecolor="lightgrey", alpha=0.6))
ax.text(0, 4.5, "YASAK\n(0 ve etrafı)", ha="center", fontsize=9, color="dimgrey")

# Yığılma çubukları (Sezen voting karakteristik dağılım)
rng = np.random.default_rng(7)
yiglma_oranlari = [4.2, 2.5, 1.4, 0.8, 0.4,    # +1 -> +5
                   0.4, 0.8, 1.4, 2.5, 4.2]    # -5 -> -1 (ters)
values_pos = [1, 2, 3, 4, 5]
values_neg = [-1, -2, -3, -4, -5]
for v, h in zip(values_pos, [4.2, 2.5, 1.4, 0.8, 0.4]):
    ax.bar(v, h, width=0.6, color="steelblue", edgecolor="navy", alpha=0.8)
for v, h in zip(values_neg, [4.2, 2.5, 1.4, 0.8, 0.4]):
    ax.bar(v, h, width=0.6, color="indianred", edgecolor="darkred", alpha=0.8)

# ±1 vurgu
ax.text(1, -0.3, "+1", ha="center", fontsize=11, fontweight="bold", color="navy")
ax.text(-1, -0.3, "−1", ha="center", fontsize=11, fontweight="bold", color="darkred")
ax.text(0, -0.3, "0\nyasak", ha="center", fontsize=9, color="black")

ax.set_xticks([-5,-4,-3,-2,-1,1,2,3,4,5])
ax.grid(alpha=0.2)

# ============================================================
# ALT — Ortak yapı haritası (kavramsal)
# ============================================================
ax = fig.add_subplot(gs[1, :])
ax.set_xlim(-0.5, 11)
ax.set_ylim(-2, 2.5)
ax.axis("off")

# Üç sistem için aynı eksen
def draw_axis(y, system_name, color_zero, color_half, color_one,
              label_zero, label_half, label_one):
    # ana çizgi
    ax.plot([1, 9], [y, y], color="black", lw=2.5)
    # 0 noktası
    ax.plot(1, y, "s", ms=18, color=color_zero, mec="black", mew=1.5)
    ax.text(1, y-0.45, label_zero, ha="center", fontsize=9, color=color_zero,
            fontweight="bold")
    # 1/2 noktası
    ax.plot(5, y, "o", ms=24, color=color_half, mec="black", mew=1.5)
    ax.text(5, y-0.45, label_half, ha="center", fontsize=9, color=color_half,
            fontweight="bold")
    # 1 noktası
    ax.plot(9, y, "^", ms=18, color=color_one, mec="black", mew=1.5)
    ax.text(9, y-0.45, label_one, ha="center", fontsize=9, color=color_one,
            fontweight="bold")
    # sistem adı
    ax.text(-0.2, y, system_name, ha="left", va="center", fontsize=10.5,
            fontweight="bold")

# 1) Riemann
draw_axis(1.5, "Riemann ζ",
          "lightgrey", "crimson", "lightblue",
          "σ=0\nyasak", "σ=1/2\nsıfırlar (γ_n)", "σ=1\nbirim")

# 2) Sezen voting
draw_axis(0.0, "Sezen voting",
          "lightgrey", "steelblue", "lightgreen",
          "0\nyasak", "±1\nyığılma", "±∞\nuç yön")

# 3) Çift yarık
draw_axis(-1.5, "Çift yarık",
          "lightgrey", "gold", "lightblue",
          "yarık\nyok", "girişim\nmerkezi", "tabanca\n(kaynak)")

# Üst başlık
ax.text(5, 2.2, "ORTAK YAPI: 0 (yasak/yok)  ──  1/2 (yığılma/denge)  ──  1 (birim/kaynak)",
        ha="center", fontsize=12, fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.5", facecolor="lightyellow",
                  edgecolor="black"))

# Açıklama notu
ax.text(5, -2.0,
        "Üç sistem de aynı üçgeni çiziyor: yasak boşluk ile birim/kaynak arasında, denge tam ortada (1/2).\n"
        "Sezen voting'de denge ±1'de (kenar) kristalleşir; Riemann'da denge tam ortada (σ=1/2); çift yarıkta plakanın ortasında.",
        ha="center", fontsize=9.5, style="italic", color="dimgrey")

plt.suptitle("0 ─ 1/2 ─ 1 ÜÇGENİ  —  Üç sistemde tek yapı\n"
             "(Uğur'un sezgisi, 16 Mayıs 2026)",
             fontsize=14, fontweight="bold", y=0.995)

out = Path(__file__).parent / "09_ucgen_resim.png"
plt.savefig(out, dpi=130, bbox_inches="tight")
print(f"Kaydedildi: {out}")
