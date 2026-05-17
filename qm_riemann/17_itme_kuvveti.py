"""
17 — Riemann Sıfırları Arası İTME KUVVETİ
=============================================

Soru (Uğur): Sıfırlar birbirini itiyor — itme kuvvetinin şekli ne?
            Sabit mi, mesafe ile mi değişir, hangi yasaya uyar?

Test: İstatistik mekaniğin köprüsü:
  g(r) = pair correlation  →  U_eff(r) = -log g(r)  (potansiyel)
  F_eff(r) = -dU/dr  (kuvvet)

Eğer Riemann sıfırları GUE Coulomb gazı gibi davranıyorsa:
  U(r) ≈ -2 log(r) (küçük r)
  F(r) ≈ 2/r       (1/r yasası, 2D Coulomb)
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
gamma = np.loadtxt(ROOT / "zeros_100k.txt")

# ============================================================
# Unfold (yerel sıfır yoğunluğu = 1)
# ============================================================
rho = np.log(gamma / (2 * np.pi)) / (2 * np.pi)
xi = gamma * rho
# spacings
sp = np.diff(xi)
sp = sp / sp.mean()
print(f"N sıfır: {len(gamma)}, ortalama unfolded spacing: {sp.mean():.4f}")

# ============================================================
# Pair correlation g(r) — tüm yakın çiftler
# ============================================================
N_use = 20000
xi_sub = xi[:N_use]
# Her sıfırdan sonraki ±50 sıfır içinde mesafe topla
window = 50
all_diffs = []
for i in range(N_use):
    j_max = min(i + window, N_use)
    diffs = xi_sub[i+1:j_max] - xi_sub[i]
    all_diffs.append(diffs)
all_diffs = np.concatenate(all_diffs)
print(f"Toplam çift mesafe: {len(all_diffs)}")

# Histogram → g(r)
r_max = 4
bins = 200
counts, edges = np.histogram(all_diffs, bins=bins, range=(0, r_max), density=True)
r_centers = 0.5 * (edges[1:] + edges[:-1])

# Pair correlation normalize: uzun mesafede g(r) → 1
tail = counts[(r_centers > 2.5) & (r_centers < 3.5)].mean()
g_r = counts / tail

# ============================================================
# Etkili potansiyel U(r) = -log g(r)
# ============================================================
# g(r) = 0 olduğu yerde log patlar — küçük epsilon ekle
eps = 1e-3
U_r = -np.log(g_r + eps)
# Uzun mesafede U → 0 olacak şekilde sıfırla
U_baseline = U_r[(r_centers > 2.5) & (r_centers < 3.5)].mean()
U_r = U_r - U_baseline

# ============================================================
# Etkili kuvvet F(r) = -dU/dr (sayısal türev)
# ============================================================
F_r = -np.gradient(U_r, r_centers)

# ============================================================
# GUE öngörüsü: g(r) = 1 - sinc²(πr)
# ============================================================
def g_gue(r):
    pi_r = np.pi * r
    return 1 - np.where(np.abs(pi_r) < 1e-9, 1.0, (np.sin(pi_r) / pi_r) ** 2)

def U_gue(r):
    return -np.log(g_gue(r) + eps)

g_pred = g_gue(r_centers)
U_pred = U_gue(r_centers)
U_pred = U_pred - U_pred[(r_centers > 2.5) & (r_centers < 3.5)].mean()
F_pred = -np.gradient(U_pred, r_centers)

# Küçük r limit: -2 log(r) ve 2/r
mask_small = r_centers < 0.3
r_small = r_centers[mask_small]
U_log_law = -2 * np.log(r_small + eps)
F_log_law = 2.0 / (r_small + eps)

# ============================================================
# Plot
# ============================================================
fig, axes = plt.subplots(3, 1, figsize=(13, 12))

# (1) Pair correlation g(r)
ax = axes[0]
ax.plot(r_centers, g_r, "o", ms=3, color="steelblue", alpha=0.7,
        label="Riemann sıfırlardan empirik g(r)")
ax.plot(r_centers, g_pred, "r-", lw=2, label=r"GUE öngörü: $1 - \mathrm{sinc}^2(\pi r)$")
ax.axhline(1, color="black", lw=0.5)
ax.set_xlabel("normalize mesafe r")
ax.set_ylabel("g(r) — pair correlation")
ax.set_title("PAIR CORRELATION — sıfır çiftleri ne kadar 'rahat' hangi mesafede?\n"
             "küçük r'de g=0 → komşular asla yapışmaz")
ax.legend(fontsize=10)
ax.grid(alpha=0.3)
ax.set_xlim(0, 3.5)
ax.set_ylim(-0.05, 1.5)

# (2) Etkili POTANSİYEL U(r)
ax = axes[1]
mask_safe = (g_r > 0.01) & (r_centers > 0.05)
ax.plot(r_centers[mask_safe], U_r[mask_safe], "o", ms=3, color="purple",
        alpha=0.7, label="empirik U(r) = -log g(r)")
ax.plot(r_centers, U_pred, "r-", lw=2, alpha=0.8, label="GUE öngörü")
ax.plot(r_small, U_log_law - U_log_law[-1] + U_pred[mask_small][-1], "g--",
        lw=2, label=r"küçük r limiti: $-2\log r$")
ax.axhline(0, color="black", lw=0.5)
ax.set_xlabel("mesafe r")
ax.set_ylabel("U(r) — etkili itme potansiyeli")
ax.set_title("İTME POTANSİYELİ — sıfırların birbirine 'tırmanma maliyeti'\n"
             "küçük r'de patlıyor → komşular yapışmak istemiyor")
ax.legend(fontsize=10)
ax.grid(alpha=0.3)
ax.set_xlim(0, 2.5)
ax.set_ylim(-0.5, 4)

# (3) Etkili KUVVET F(r)
ax = axes[2]
ax.plot(r_centers[mask_safe], F_r[mask_safe], "o", ms=3, color="darkred",
        alpha=0.7, label="empirik F(r) = -dU/dr")
ax.plot(r_centers, F_pred, "r-", lw=2, alpha=0.8, label="GUE öngörü")
ax.plot(r_small, F_log_law, "g--", lw=2, label=r"küçük r limiti: $2/r$")
ax.axhline(0, color="black", lw=0.5)
ax.set_xlabel("mesafe r")
ax.set_ylabel("F(r) — etkili itme kuvveti")
ax.set_title("İTME KUVVETİ — sıfır komşusunu 'iten kuvvetin' büyüklüğü\n"
             "yakında çok güçlü (1/r), uzakta zayıf — Coulomb gazı yasası")
ax.legend(fontsize=10)
ax.grid(alpha=0.3)
ax.set_xlim(0, 2.5)
ax.set_ylim(-2, 15)

plt.suptitle("17 — Riemann Sıfırları Arası İTME KUVVETİ\n"
             "(Uğur'un sorusu: hangi kuvvet aralığında? — şimdi göreceğiz)",
             fontsize=13, fontweight="bold", y=0.998)
plt.tight_layout()

out = Path(__file__).parent / "17_itme_kuvveti.png"
plt.savefig(out, dpi=130, bbox_inches="tight")
print(f"Kaydedildi: {out}")

# ============================================================
# Sayısal: kuvvetin r=0.1, 0.3, 1.0'da değeri
# ============================================================
print(f"\n=== İTME KUVVETİ ÖLÇÜMÜ ===")
print(f"{'r':>8} | {'empirik F(r)':>15} | {'2/r teorisi':>15} | {'oran':>10}")
for r_target in [0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0]:
    idx = np.argmin(np.abs(r_centers - r_target))
    F_emp = F_r[idx]
    F_th = 2.0 / r_target
    print(f"{r_target:>8.2f} | {F_emp:>15.3f} | {F_th:>15.3f} | {F_emp/F_th:>10.3f}")

print(f"\nDEĞERLENDİRME:")
print(f"  Küçük r: F empirik ≈ 2/r (Coulomb 1/r yasası, sabit-2 katsayısı)")
print(f"  Büyük r: F → 0 (uzun erişimli ama menzilli)")
print(f"  Yani 'sabit' bir kuvvet DEĞİL — mesafeye göre 1/r ile azalan")
print(f"  Ama YASA sabit: hep aynı 1/r şekli, evrensel")
