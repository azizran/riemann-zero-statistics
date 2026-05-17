"""
QM × Riemann — Spiral Dans
============================

Sezgi: "İki sistem ikili düzlemde çalışıyor; ζ sıfırlarda spiral oluşturur,
QM dalga fonksiyonu da kompleks düzlemde dalga olarak yaşar."

Bu dosya iki şeyi yan yana koyar:

  (1) ζ(½+it) için t ∈ [0, 50] arasında kompleks değerin trajektorisi
      → her Riemann sıfırında origin'den geçer
      → |ζ| dalgalı: zirveler ↔ vadiler (sıfırlar)

  (2) QM bir kuzen: Lorentzian rezonans toplamı S(E) = Σ Γ/((E-Eₙ)+iΓ)
      (Eₙ pozisyonlarını ilk 10 Riemann γ'sından alıyoruz)
      → aynı topolojik yapı: zirveler her rezonansta, fazlar dolanıyor
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from pathlib import Path
import mpmath as mp

mp.mp.dps = 30  # yeterli hassasiyet

ROOT = Path(__file__).resolve().parent.parent
gamma = np.loadtxt(ROOT / "zeros_100k.txt")

# ============================================================
# (1) ζ(½+it) trajektorisi
# ============================================================
t_max = 50.0
N = 4000
t_grid = np.linspace(0.5, t_max, N)

print("ζ(½+it) hesaplanıyor (yüksek hassasiyet, dakikalar değil saniyeler)...")
zeta_vals = np.array([complex(mp.zeta(mp.mpc(0.5, tt))) for tt in t_grid])

# sıfırlar bu pencerede
zeros_in_window = gamma[gamma < t_max]
print(f"Pencere içindeki Riemann sıfırı sayısı: {len(zeros_in_window)}")
print(f"İlk birkaçı: {zeros_in_window[:5]}")

# ============================================================
# (2) QM analog: Lorentzian rezonans toplamı (ζ ile aynı sıfır yerlerinde)
# ============================================================
# S(E) = Σ_n  iΓ / ((E - E_n) + iΓ)   →  her E_n'de "kuantum rezonansı"
# E_n'i ilk 10 Riemann sıfırından alalım
E_res = zeros_in_window[:10]
Gamma = 0.5   # genişlik (zenginleşir → spiraller iç içe)

def S_qm(E):
    return np.sum(1j * Gamma / ((E - E_res) + 1j * Gamma))

S_vals = np.array([S_qm(tt) for tt in t_grid])

# ============================================================
# Plot — üç panel
# ============================================================
fig = plt.figure(figsize=(15, 10))
gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.3], hspace=0.32, wspace=0.25)

# --- üst-sol: |ζ| zaman dilimi ---
ax_top_l = fig.add_subplot(gs[0, 0])
ax_top_l.plot(t_grid, np.abs(zeta_vals), color="navy", lw=1.2)
for z in zeros_in_window:
    ax_top_l.axvline(z, color="red", lw=0.6, alpha=0.5)
ax_top_l.set_xlim(0, t_max)
ax_top_l.set_xlabel("t")
ax_top_l.set_ylabel(r"$|\zeta(\frac{1}{2}+it)|$")
ax_top_l.set_title("ζ büyüklüğü zaman ekseninde — kırmızı çizgiler: sıfırlar")
ax_top_l.grid(alpha=0.3)

# --- üst-sağ: |S_qm| zaman dilimi ---
ax_top_r = fig.add_subplot(gs[0, 1])
ax_top_r.plot(t_grid, np.abs(S_vals), color="darkgreen", lw=1.2)
for z in E_res:
    ax_top_r.axvline(z, color="purple", lw=0.6, alpha=0.5)
ax_top_r.set_xlim(0, t_max)
ax_top_r.set_xlabel("E")
ax_top_r.set_ylabel(r"$|S_{QM}(E)|$")
ax_top_r.set_title("QM rezonans toplamı — mor çizgiler: rezonanslar (= γ_n)")
ax_top_r.grid(alpha=0.3)

# --- alt-sol: ζ kompleks düzlem spirali ---
ax_bl = fig.add_subplot(gs[1, 0])
# colormap: t arttıkça renk değişir
points = np.array([zeta_vals.real, zeta_vals.imag]).T.reshape(-1, 1, 2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)
norm = plt.Normalize(t_grid.min(), t_grid.max())
lc = LineCollection(segments, cmap="plasma", norm=norm, linewidth=1.3, alpha=0.9)
lc.set_array(t_grid[:-1])
ax_bl.add_collection(lc)

# sıfırlar = origin geçişleri — gerçek sıfır t değerinde ζ değerini bul (≈ 0)
ax_bl.plot(0, 0, "o", ms=14, mfc="none", mec="red", mew=2,
           label=f"origin ({len(zeros_in_window)} sıfır geçişi)")
ax_bl.axhline(0, color="k", lw=0.4, alpha=0.4)
ax_bl.axvline(0, color="k", lw=0.4, alpha=0.4)
ax_bl.set_xlabel(r"Re $\zeta(\frac{1}{2}+it)$")
ax_bl.set_ylabel(r"Im $\zeta(\frac{1}{2}+it)$")
ax_bl.set_title("ζ(½+it) kompleks düzlemde spiral\n(renk = t, plasma haritası)")
ax_bl.set_aspect("equal")
ax_bl.grid(alpha=0.3)
ax_bl.legend(loc="upper right", fontsize=9)

cb = plt.colorbar(lc, ax=ax_bl, label="t", shrink=0.85)

# --- alt-sağ: S_qm kompleks düzlem spirali ---
ax_br = fig.add_subplot(gs[1, 1])
points2 = np.array([S_vals.real, S_vals.imag]).T.reshape(-1, 1, 2)
segments2 = np.concatenate([points2[:-1], points2[1:]], axis=1)
lc2 = LineCollection(segments2, cmap="viridis", norm=norm, linewidth=1.3, alpha=0.9)
lc2.set_array(t_grid[:-1])
ax_br.add_collection(lc2)

# Lorentzian rezonansları halka yapar; rezonanslar = iΓ noktası civarında
ax_br.plot(0, Gamma, "o", ms=12, mfc="none", mec="purple", mew=2,
           label=f"iΓ ekseni (10 rezonans)")
ax_br.axhline(0, color="k", lw=0.4, alpha=0.4)
ax_br.axvline(0, color="k", lw=0.4, alpha=0.4)
ax_br.set_xlabel(r"Re $S_{QM}(E)$")
ax_br.set_ylabel(r"Im $S_{QM}(E)$")
ax_br.set_title("QM rezonans toplamı kompleks düzlemde\n(renk = E, viridis)")
ax_br.set_aspect("equal")
ax_br.grid(alpha=0.3)
ax_br.legend(loc="upper right", fontsize=9)

cb2 = plt.colorbar(lc2, ax=ax_br, label="E", shrink=0.85)

plt.suptitle("Spiral Dans — ζ(½+it) ile bir QM rezonans toplamının kompleks-düzlem yörüngeleri",
             fontsize=13, y=0.995)

out = Path(__file__).parent / "06_spiral_dans.png"
plt.savefig(out, dpi=130, bbox_inches="tight")
print(f"Kaydedildi: {out}")

# ============================================================
# Sayısal: ζ origin'e ne kadar yakın geçiyor?
# ============================================================
print(f"\nζ(½+it_n) sıfırlarda gerçekten ~0 mı?")
for tt in zeros_in_window[:5]:
    z = complex(mp.zeta(mp.mpc(0.5, tt)))
    print(f"  t={tt:.5f}:  ζ = ({z.real:+.2e}, {z.imag:+.2e})  |ζ|={abs(z):.2e}")
