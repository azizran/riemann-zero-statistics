"""
QM × Riemann — Bizim Açımızdan Hilbert-Pólya Saldırısı
=========================================================

POLYGON-π SEZGİSİ:
  yasak bölge   →   yığılma kenarı
  sıfır yasak   →   ±1 yığılma   (Sezen voting)
  sıfıra varamama → fire (c, c_iç, ε)  (polygon)

RIEMANN PARALELİ:
  fonksiyonel simetri ξ(s)=ξ(1-s) → σ=1/2 yığılma çizgisi
  σ ≠ 1/2 "yasak" (RH varsayar)   → tüm önemsiz sıfırlar σ=1/2'de

YAPISAL HİPOTEZ:
  Polygon fire sabitlerimiz (c, c_iç, ε), ζ'nin GUE istatistiğinde
  doğal işaretler taşıyor mu? Yoksa bizim "fire" ve Riemann'ın "kritik
  çizgi" mekanizmaları sadece kavramsal akraba mı?

TEST:
  K(τ) ve R_2(r)'yi c, c_iç, ε noktalarında değerlendir.
  Integraller, oranlar — bir şey "denk geliyor" mu?
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import mpmath as mp

# ============================================================
# Polygon-π fire sabitleri (MEMORY'den)
# ============================================================
c     = 0.85393950843     # OEIS A394396
c_ic  = 1.18055935988     # OEIS A394496 (= 1/c + γ_geo)
eps   = 0.79443           # ≈ π/4 + γ_geo
gamma_geo = 300 / (31800 - np.pi)
kappa = c * c_ic - 1      # ≈ 0.00806

print("Polygon-π fire sabitleri:")
print(f"  c        = {c:.8f}    (~ 4-π aileden)")
print(f"  c_iç     = {c_ic:.8f}  (= 1/c + γ_geo)")
print(f"  ε        = {eps:.8f}    (~ π/4 + γ_geo)")
print(f"  κ        = {kappa:.8f}  (c·c_iç - 1)")
print(f"  γ_geo    = {gamma_geo:.10f}")
print()

# ============================================================
# (1) Sezgi: K(τ) bizim sabitlerimizde özel mi?
# ============================================================
ROOT = Path(__file__).resolve().parent.parent
t_zeros = np.loadtxt(ROOT / "zeros_100k.txt")

# Unfold
rho = np.log(t_zeros / (2 * np.pi)) / (2 * np.pi)
xi = t_zeros * rho

# K(τ) ince ızgarada hesapla
WINDOW = 1500
N_WIN  = 100
tau = np.linspace(0.001, 1.5, 800)

rng = np.random.default_rng(42)
starts = rng.integers(0, len(xi) - WINDOW - 1, size=N_WIN)
K_emp = np.zeros_like(tau)
for s in starts:
    seg = xi[s:s+WINDOW]
    seg = seg - seg.mean()
    phase = np.exp(2j * np.pi * np.outer(tau, seg))
    K_emp += np.abs(phase.sum(axis=1))**2 / WINDOW
K_emp /= N_WIN

# K(c), K(c_iç), K(ε), K(κ), K(γ_geo) değerleri
def Ktau_at(t):
    """Yakın τ değerinden interpolasyon"""
    return np.interp(t, tau, K_emp)

print("K(τ) değerleri özel noktalarda:")
print(f"  K(γ_geo  ={gamma_geo:.5f}) = {Ktau_at(gamma_geo):.5f}   GUE: {min(gamma_geo,1):.5f}")
print(f"  K(κ      ={kappa:.5f}) = {Ktau_at(kappa):.5f}   GUE: {min(kappa,1):.5f}")
print(f"  K(ε      ={eps:.5f}) = {Ktau_at(eps):.5f}   GUE: {min(eps,1):.5f}")
print(f"  K(c      ={c:.5f}) = {Ktau_at(c):.5f}   GUE: {min(c,1):.5f}")
print(f"  K(c_iç   ={c_ic:.5f}) = {Ktau_at(c_ic):.5f}   GUE: {min(c_ic,1):.5f}")
print()

# ============================================================
# (2) "Yasak bölge → yığılma" yapısal hipotez:
#     ζ sıfır olmayan değerleri, σ ≠ 1/2'de "yasaklanmış" gibi
#     davranıyor mu? Sayısal kontrol: ζ(σ + it)'nin σ = 1/2'den
#     ne hızla "kaçtığını" ölç. Bizim Sezen ±1 yığılma analoğu var mı?
# ============================================================
mp.mp.dps = 30
# Bir sıfır komşuluğu (t ≈ γ_1 = 14.1347)
t0 = 14.13472514
sigma_grid = np.linspace(0.05, 0.95, 100)
zeta_abs = np.array([abs(complex(mp.zeta(mp.mpc(s, t0)))) for s in sigma_grid])

# ============================================================
# (3) Form factor küçük τ doğrusal eğimi
#     Sezen voting'de "yığılma kenarındaki yoğunluk"a paralel
# ============================================================
mask = (tau > 0.01) & (tau < c)   # bizim c'ye kadar
slope_to_c, _ = np.polyfit(tau[mask], K_emp[mask], 1)
mask2 = (tau > c) & (tau < 1.0)
slope_after_c, _ = np.polyfit(tau[mask2], K_emp[mask2], 1)
print(f"K(τ) eğimi 0→c bölgesinde: {slope_to_c:.4f}")
print(f"K(τ) eğimi c→1 bölgesinde: {slope_after_c:.4f}")
print(f"  (GUE: ikisi de 1.0; sapma yapısal mı yoksa örneklem gürültüsü mü?)")

# ============================================================
# Plot
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

# Sol: K(τ) + bizim sabitlerimiz işaretli
ax = axes[0]
ax.plot(tau, K_emp, color="steelblue", lw=0.8, alpha=0.8, label="ζ empirik K(τ)")
ax.plot(tau, np.minimum(tau, 1), "r-", lw=1.8, alpha=0.7, label="GUE: min(τ,1)")
for label, val, color in [
    ("γ_geo", gamma_geo, "olive"),
    ("κ", kappa, "purple"),
    ("ε", eps, "orange"),
    ("c", c, "crimson"),
    ("c_iç", c_ic, "darkblue"),
]:
    ax.axvline(val, color=color, lw=1.2, ls="--", alpha=0.7,
               label=f"{label} = {val:.4f}")
ax.set_xlim(0, 1.5)
ax.set_ylim(-0.05, 1.4)
ax.set_xlabel("τ")
ax.set_ylabel("K(τ)")
ax.set_title("Bizim polygon-π sabitleri ζ form factor üzerinde\n"
             "(c ve c_iç sınır geçişin iki yanında)")
ax.legend(loc="upper left", fontsize=8, ncol=2)
ax.grid(alpha=0.3)

# Sağ: ζ(σ+it₀) σ ile değişim — "yasak bölge" sezgisi gerçek mi?
ax = axes[1]
ax.plot(sigma_grid, zeta_abs, color="darkred", lw=1.8)
ax.axvline(0.5, color="red", lw=1.2, ls="-", alpha=0.7, label="σ=1/2 kritik")
ax.axhline(0, color="k", lw=0.4)
ax.set_xlabel("σ")
ax.set_ylabel(r"$|\zeta(\sigma + i \cdot 14.135)|$")
ax.set_title("ζ(½+iγ₁) sıfır. σ'yu kaydırınca yokoluyor mu?\n"
             "(sıfır σ=1/2'ye 'sıkışmış' mı, kayınca açılıyor mu?)")
ax.legend(fontsize=10)
ax.grid(alpha=0.3)

plt.tight_layout()
out = Path(__file__).parent / "08_bizim_acidan.png"
plt.savefig(out, dpi=130, bbox_inches="tight")
print(f"\nKaydedildi: {out}")

# ============================================================
# Yapısal numerikler
# ============================================================
print(f"\n--- 'Yasak bölge → yığılma' yapısal test ---")
print(f"σ=1/2'den uzaklaştıkça |ζ| büyüklüğü (yığılma 'noktasal' mı?):")
for s in [0.3, 0.4, 0.45, 0.49, 0.5, 0.51, 0.55, 0.6, 0.7]:
    val = abs(complex(mp.zeta(mp.mpc(s, t0))))
    print(f"  σ={s:.2f}:  |ζ| = {val:.5f}")
# Beklenti: σ=0.5'te 0, başka her yerde > 0
# Sezen voting'de: x=0'da yasak, |x|=1'de yığılma; burada σ=0.5'te yığılma (ters)
print(f"\nÇıkarım: Sezen voting analojisi TERS yönlü — orada YASAK noktada hiç yok,")
print(f"burada KRİTİK çizgide tam sıfırlar var. Yapısal paralel kavramsal, mekanik DEĞIL.")
