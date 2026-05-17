"""
21 — Lokal Kayma Sistematik mi? (S1 düzeltmesi)
====================================================

Uğur'un düzeltmesi: "yön değişiyor ama kayma var"

Test: küçük pencerelerde Z(t)'nin yerel asimetrisi nasıl dağılıyor?
  - Mean = ?  (küresel kayma)
  - RMS  = ?  (lokal kayma sıfırdan farklı mı)
  - Dağılım normal mi?

Eğer mean=0 ama RMS>0 ise: küresel simetri + lokal kayma
  = kuantum "süperpozisyon + collapse" imzası
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import mpmath as mp

mp.mp.dps = 25

# ============================================================
# Z(t) için ÇOK daha geniş hesap
# ============================================================
t_lo, t_hi = 10.0, 500.0   # daha geniş pencere (12000 nokta yetiyor)
N = 15000
t_grid = np.linspace(t_lo, t_hi, N)

print(f"Z(t) hesaplanıyor: t ∈ [{t_lo}, {t_hi}], {N} nokta...")
Z = np.zeros(N)
for i in range(N):
    Z[i] = float(mp.siegelz(t_grid[i]))
    if i % 2000 == 0:
        print(f"  i={i}/{N}, t={t_grid[i]:.1f}")
print(f"Z aralığı: [{Z.min():+.3f}, {Z.max():+.3f}]")

# ============================================================
# Pencere asimetrileri — birçok pencere boyu için
# ============================================================
def window_asymmetries(t, Z, n_windows):
    t_lo, t_hi = t[0], t[-1]
    ws = (t_hi - t_lo) / n_windows
    asyms = []
    for i in range(n_windows):
        ts, te = t_lo + i*ws, t_lo + (i+1)*ws
        mask = (t >= ts) & (t < te)
        Zw = Z[mask]
        if len(Zw) < 5: continue
        pa = np.trapezoid(np.maximum(Zw, 0), t[mask])
        na = np.trapezoid(np.maximum(-Zw, 0), t[mask])
        if pa + na > 0:
            asyms.append((pa - na)/(pa + na))
    return np.array(asyms)

# Farklı pencere sayıları için test
print(f"\n{'#pencere':>10} {'pencere boyu':>14} {'mean':>10} {'RMS':>10} {'std':>10}")
all_asyms = {}
for n_w in [5, 10, 20, 50, 100, 200, 500]:
    a = window_asymmetries(t_grid, Z, n_w)
    all_asyms[n_w] = a
    rms = np.sqrt(np.mean(a**2))
    print(f"{n_w:>10} {(t_hi-t_lo)/n_w:>14.2f} {a.mean()*100:>+9.2f}% {rms*100:>+9.2f}% {a.std()*100:>+9.2f}%")

# ============================================================
# Plot
# ============================================================
fig = plt.figure(figsize=(14, 10))
gs = fig.add_gridspec(2, 2, hspace=0.32, wspace=0.25)

# (1) RMS vs pencere boyu
ax = fig.add_subplot(gs[0, 0])
n_ws = sorted(all_asyms.keys())
window_lens = [(t_hi-t_lo)/n for n in n_ws]
rmss = [np.sqrt(np.mean(all_asyms[n]**2))*100 for n in n_ws]
means = [all_asyms[n].mean()*100 for n in n_ws]
ax.semilogx(window_lens, rmss, "o-", color="purple", lw=2, ms=8, label="RMS lokal kayma")
ax.semilogx(window_lens, [abs(m) for m in means], "s--", color="orange",
            lw=1.5, ms=6, label="|küresel kayma|")
ax.axhline(0, color="black", lw=0.5)
ax.set_xlabel("pencere boyu (t-birimi)")
ax.set_ylabel("kayma (%)")
ax.set_title("Lokal kayma vs pencere boyu\n"
             "küresel kayma → 0; lokal RMS → büyük (sistemik)")
ax.legend(fontsize=10)
ax.grid(alpha=0.3, which="both")

# (2) Asimetri dağılımı (n_w=100 için)
ax = fig.add_subplot(gs[0, 1])
a100 = all_asyms[100] * 100
ax.hist(a100, bins=30, color="steelblue", alpha=0.7, edgecolor="navy", density=True)
# Normal dağılım üst çiz
from scipy.stats import norm
x = np.linspace(a100.min(), a100.max(), 200)
ax.plot(x, norm.pdf(x, a100.mean(), a100.std()), 'r-', lw=2,
        label=f"Normal(ort={a100.mean():.2f}, std={a100.std():.2f})")
ax.axvline(0, color="black", lw=1, ls="--", alpha=0.6, label="simetri çizgisi")
ax.axvline(a100.mean(), color="red", lw=1.5, label=f"ortalama={a100.mean():.2f}%")
ax.set_xlabel("pencere asimetrisi (%)")
ax.set_ylabel("yoğunluk")
ax.set_title(f"100 pencerenin asimetri dağılımı\n"
             "küresel ortalama 0'a yakın, ama yayılım büyük")
ax.legend(fontsize=9)
ax.grid(alpha=0.3)

# (3) Z(t) bir kısım göster
ax = fig.add_subplot(gs[1, 0])
mask = (t_grid >= 100) & (t_grid <= 200)
ax.plot(t_grid[mask], Z[mask], color="navy", lw=0.7)
ax.fill_between(t_grid[mask], Z[mask], 0, where=(Z[mask] > 0), alpha=0.25, color="blue")
ax.fill_between(t_grid[mask], Z[mask], 0, where=(Z[mask] < 0), alpha=0.25, color="red")
ax.axhline(0, color="black", lw=0.5)
ax.set_xlim(100, 200)
ax.set_xlabel("t")
ax.set_ylabel("Z(t)")
ax.set_title("Z(t) örnek pencere (t ∈ [100, 200])\n"
             "her küçük pencerede +/- alanlar AYNI değil")
ax.grid(alpha=0.3)

# (4) Asimetrinin t boyunca değişimi
ax = fig.add_subplot(gs[1, 1])
n_w = 100
window_size = (t_hi - t_lo) / n_w
window_centers = t_lo + (np.arange(n_w) + 0.5) * window_size
ax.bar(window_centers, all_asyms[100]*100, width=window_size*0.85,
       color=["red" if a > 0 else "blue" for a in all_asyms[100]], alpha=0.65)
ax.axhline(0, color="black", lw=0.7)
ax.axhline(rmss[n_ws.index(100)], color="purple", lw=1.2, ls="--",
           label=f"+RMS = {rmss[n_ws.index(100)]:.1f}%")
ax.axhline(-rmss[n_ws.index(100)], color="purple", lw=1.2, ls="--")
ax.set_xlabel("t (pencere merkezi)")
ax.set_ylabel("asimetri (%)")
ax.set_title(f"100 pencerede asimetri — yön rastgele ama büyüklük sıfırdan farklı")
ax.legend(fontsize=10)
ax.grid(alpha=0.3)

plt.suptitle("21 — Lokal Kayma Sistematik mi?\n"
             "Uğur'un düzeltmesi: 'yön değişir ama kayma her zaman var'",
             fontsize=13, fontweight="bold", y=0.998)

out = Path(__file__).parent / "21_lokal_kayma.png"
plt.savefig(out, dpi=130, bbox_inches="tight")
print(f"\nKaydedildi: {out}")

# ============================================================
# Kuantum collapse karşılaştırma — null hypothesis
# ============================================================
print(f"\n=== KUANTUM COLLAPSE İMZASI TESTİ ===")
print(f"")
print(f"Eğer kayma 'gerçek bir sistemik özellik' ise:")
print(f"  - küresel mean → 0  (uzun ortalama simetrik)")
print(f"  - lokal RMS    > 0  (her küçük pencerede kayma var)")
print(f"")
print(f"Sayılar:")
for n_w in [20, 100, 500]:
    a = all_asyms[n_w]
    rms = np.sqrt(np.mean(a**2)) * 100
    se_mean = a.std()/np.sqrt(len(a)) * 100   # ortalamanın standart hatası
    z_score = abs(a.mean()*100) / se_mean if se_mean > 0 else 0
    print(f"  {n_w} pencere:  mean = {a.mean()*100:+.2f}% (z={z_score:.1f}σ),  RMS = {rms:.2f}%")

# Null test: aynı Z'nin işaretini rastgele karıştırırsak ne çıkar?
rng = np.random.default_rng(42)
Z_shuffled = Z * rng.choice([-1, 1], size=len(Z))
a_null = window_asymmetries(t_grid, Z_shuffled, 100)
print(f"\nKarıştırılmış null hipotez (100 pencere):")
print(f"  mean = {a_null.mean()*100:+.2f}%,  RMS = {np.sqrt(np.mean(a_null**2))*100:.2f}%")
print(f"\nGerçek RMS / null RMS oranı: {rmss[n_ws.index(100)] / (np.sqrt(np.mean(a_null**2))*100):.2f}")
