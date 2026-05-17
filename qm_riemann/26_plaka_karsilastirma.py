"""
26 — Plaka Karşılaştırma: Çift Yarık vs Riemann
====================================================

Uğur'un sezgisi: 1D ile 2D plaka aslında gözlemcinin perspektifinden
aynı düzlemde. Boyut görsel değil yapısal.

Şimdi iki plakayı yan yana çiziyoruz:
  - Klasik çift yarık plakası (Young's experiment)
  - Riemann plakası: σ=1/2 çizgisi üzerinde sıfırlar parlak şeritler

Soru: gözle bakınca aynı tür desen mi görünüyor?
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
gamma = np.loadtxt(ROOT / "zeros_100k.txt")

# ============================================================
# Klasik çift yarık: girişim deseni
# ============================================================
def classical_pattern(y, d, L, lam):
    """Young's çift yarık: I(y) = cos²(πd y / (λL)) × sinc²(πa y / (λL))"""
    a = d / 5    # tek yarık genişliği
    cos_term = np.cos(np.pi * d * y / (lam * L))**2
    sinc_arg = np.pi * a * y / (lam * L)
    sinc_term = np.where(np.abs(sinc_arg) < 1e-9, 1.0,
                          (np.sin(sinc_arg)/sinc_arg)**2)
    return cos_term * sinc_term

# Boyutsuz parametreler
y_classical = np.linspace(-10, 10, 2000)
I_classical = classical_pattern(y_classical, d=1, L=10, lam=1)

# ============================================================
# Riemann plakası: σ=1/2 boyunca γ_n şeritleri
# ============================================================
# Plaka uzunluğu: t ∈ [0, 200] gösterelim
t_lo, t_hi = 0, 200
gamma_window = gamma[gamma < t_hi]
print(f"Pencere içindeki sıfır: {len(gamma_window)}")

# Plaka görselleştirme — siyah arka plan, parlak şeritler
plaka_resolution = 4000
plaka_t = np.linspace(t_lo, t_hi, plaka_resolution)
# Her sıfırın etrafında parlak nokta (Gaussian smear)
plaka_intensity = np.zeros(plaka_resolution)
for g in gamma_window:
    plaka_intensity += np.exp(-((plaka_t - g) / 0.3)**2)

# ============================================================
# Plot — iki plaka yan yana
# ============================================================
fig = plt.figure(figsize=(15, 12), facecolor="black")
gs = fig.add_gridspec(4, 1, height_ratios=[1.2, 2, 2, 1.5], hspace=0.4)

# (1) Klasik çift yarık plakası — geleneksel görüntü
ax = fig.add_subplot(gs[0])
ax.set_facecolor("black")
y_plaka = np.linspace(-1, 1, 100)
Y, Y_int = np.meshgrid(y_classical, y_plaka)
intensity_2d = np.tile(I_classical, (100, 1))
ax.imshow(intensity_2d, aspect='auto', cmap='inferno',
          extent=[y_classical[0], y_classical[-1], -1, 1],
          interpolation='nearest', vmin=0, vmax=0.7)
ax.set_xticks([])
ax.set_yticks([])
ax.text(0, 0, "ÇİFT YARIK PLAKASI (Young 1801)\n"
        "parçacık dalga gibi davranır,\n"
        "plakada parlak/karanlık şeritler",
        ha="center", va="center", fontsize=12, color="white",
        bbox=dict(boxstyle="round,pad=0.5", facecolor="black", alpha=0.7,
                  edgecolor="white"))

# (2) Şiddet eğrisi klasik
ax = fig.add_subplot(gs[1])
ax.set_facecolor("black")
ax.plot(y_classical, I_classical, color="gold", lw=1.5)
ax.fill_between(y_classical, I_classical, 0, color="gold", alpha=0.3)
ax.set_xlabel("plaka konumu y", color="white")
ax.set_ylabel("şiddet I(y)", color="white")
ax.set_title("Klasik çift yarık şiddeti — cos²(πd y/λL)·sinc²",
             color="white", fontsize=11)
ax.tick_params(colors="white")
ax.spines["bottom"].set_color("white")
ax.spines["left"].set_color("white")
ax.spines["top"].set_color("white")
ax.spines["right"].set_color("white")
ax.grid(alpha=0.2)

# (3) Riemann plakası — aynı görsel dilde
ax = fig.add_subplot(gs[2])
ax.set_facecolor("black")
plaka_2d = np.tile(plaka_intensity, (100, 1))
ax.imshow(plaka_2d, aspect='auto', cmap='inferno',
          extent=[t_lo, t_hi, -1, 1], interpolation='bilinear',
          vmin=0, vmax=0.7)
ax.set_xticks(np.arange(0, t_hi+1, 25))
ax.set_yticks([])
ax.set_xlabel("plaka konumu t (σ=1/2 boyunca)", color="white")
ax.tick_params(colors="white")
ax.text((t_lo+t_hi)/2, 0, "RIEMANN PLAKASI (σ=1/2 çizgisi)\n"
        f"{len(gamma_window)} sıfır = parlak şeritler\n"
        "asal sayılar her biri kuantum kaynak",
        ha="center", va="center", fontsize=12, color="white",
        bbox=dict(boxstyle="round,pad=0.5", facecolor="black", alpha=0.7,
                  edgecolor="white"))

# (4) Şiddet eğrisi Riemann — Z(t)² gösterelim
ax = fig.add_subplot(gs[3])
ax.set_facecolor("black")
# Yoğunluk: zirve yükseklik proxy = 1/spacing or kümülatif yoğunluk
# Basitçe: her sıfırı stem olarak çiz
ax.stem(gamma_window, np.ones_like(gamma_window),
        linefmt="gold", markerfmt="o", basefmt="white")
ax.set_xlim(t_lo, t_hi)
ax.set_ylim(0, 1.2)
ax.set_xlabel("plaka konumu t", color="white")
ax.set_ylabel("sıfırlar", color="white")
ax.set_title(f"Riemann sıfır pozisyonları (parlak şerit olarak)",
             color="white", fontsize=11)
ax.tick_params(colors="white")
ax.spines["bottom"].set_color("white")
ax.spines["left"].set_color("white")
ax.grid(alpha=0.2)

plt.suptitle("26 — İKİ PLAKA, AYNI MANTIK\n"
             "klasik çift yarık (üst) — Riemann (alt)",
             fontsize=14, fontweight="bold", color="white", y=0.998)

out = Path(__file__).parent / "26_plaka_karsilastirma.png"
plt.savefig(out, dpi=130, bbox_inches="tight", facecolor="black")
print(f"Kaydedildi: {out}")
plt.close()

# ============================================================
# 2. resim: SADECE RIEMANN PLAKASI — geniş görüntü, gözle pattern ara
# ============================================================
fig = plt.figure(figsize=(20, 6), facecolor="black")
ax = fig.add_subplot(1, 1, 1)
ax.set_facecolor("black")

# Daha geniş — 500'e kadar
t_hi2 = 500
gamma_w2 = gamma[gamma < t_hi2]
plaka_t2 = np.linspace(0, t_hi2, 8000)
plaka_int2 = np.zeros(8000)
for g in gamma_w2:
    plaka_int2 += np.exp(-((plaka_t2 - g) / 0.25)**2)

plaka_2d = np.tile(plaka_int2, (300, 1))
ax.imshow(plaka_2d, aspect='auto', cmap='inferno',
          extent=[0, t_hi2, 0, 1], interpolation='bilinear',
          vmin=0, vmax=0.8)
ax.set_xticks(np.arange(0, t_hi2+1, 50))
ax.set_yticks([])
ax.set_xlabel("σ = 1/2 ÇİZGİSİ ÜZERİNDE PLAKA  →  t  →", color="white", fontsize=12)
ax.tick_params(colors="white")
ax.set_title(f"RIEMANN PLAKASI: {len(gamma_w2)} sıfır, σ=1/2 boyunca\n"
             "gözle pattern var mı?",
             color="white", fontsize=14, fontweight="bold")

out2 = Path(__file__).parent / "26_riemann_plaka_genis.png"
plt.savefig(out2, dpi=140, bbox_inches="tight", facecolor="black")
print(f"Kaydedildi: {out2}")
plt.close()

# ============================================================
# 3. resim: PLAKA + ALTTA Z(t) dalgası
# ============================================================
import mpmath as mp
mp.mp.dps = 20

print("Z(t) hesaplanıyor (geniş pencere)...")
t_for_Z = np.linspace(0.5, t_hi2, 6000)
Z = np.array([float(mp.siegelz(tt)) for tt in t_for_Z])

fig, axes = plt.subplots(2, 1, figsize=(20, 8), facecolor="black",
                         gridspec_kw={"height_ratios": [1, 1.5], "hspace": 0.1})

# Üst: plaka
ax = axes[0]
ax.set_facecolor("black")
plaka_2d = np.tile(plaka_int2, (100, 1))
ax.imshow(plaka_2d, aspect='auto', cmap='inferno',
          extent=[0, t_hi2, 0, 1], interpolation='bilinear', vmin=0, vmax=0.8)
ax.set_xticks([])
ax.set_yticks([])
ax.set_title("PLAKA — σ=1/2 çizgisi üzerinde sıfırlar (parlak şerit)",
             color="white", fontsize=12)

# Alt: Z(t) dalgası — KAYNAK DALGAYI göster
ax = axes[1]
ax.set_facecolor("black")
ax.plot(t_for_Z, Z, color="cyan", lw=0.8)
ax.fill_between(t_for_Z, Z, 0, where=(Z>0), alpha=0.25, color="cyan")
ax.fill_between(t_for_Z, Z, 0, where=(Z<0), alpha=0.25, color="orange")
ax.axhline(0, color="white", lw=0.5)
# Sıfırlar
for g in gamma_w2:
    ax.axvline(g, color="gold", lw=0.3, alpha=0.6)
ax.set_xlim(0, t_hi2)
ax.set_xlabel("t", color="white", fontsize=12)
ax.set_ylabel("Z(t) (dalga genliği)", color="white", fontsize=11)
ax.tick_params(colors="white")
ax.spines["bottom"].set_color("white")
ax.spines["left"].set_color("white")
ax.set_title("ARKADAKİ DALGA — Z(t) salınıyor; her sıfır geçişi plakada parlak şerit",
             color="white", fontsize=11)
ax.grid(alpha=0.15)

out3 = Path(__file__).parent / "26_plaka_dalga_birlikte.png"
plt.savefig(out3, dpi=140, bbox_inches="tight", facecolor="black")
print(f"Kaydedildi: {out3}")
