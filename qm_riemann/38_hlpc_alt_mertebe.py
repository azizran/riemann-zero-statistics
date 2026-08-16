"""
38 — HLP-C ALT-MERTEBE KONTROLÜ (16 Ağustos 2026)
===================================================

Hughes–Lugmayer–Pearce-Crump (arXiv:2411.05573, Teorem 1, RH altında):

    Σ_{0<γ≤T} max² = (e²−5)/2 · (T/2π)L² + α₋₁(T/2π)L + (T/2π)Σ αₙ/Lⁿ

    L = log(T/2π)
    α₋₁ = 5 − e² − 10γ₀ + 2e²γ₀            (γ₀ = Euler–Mascheroni)
    α₀  = 12γ₁ − 4e²γ₁ − 5 + e² + 10γ₀ − 2e²γ₀ − 4γ₀²   (γ₁ = 1. Stieltjes)

Bizim ölçüm pencere-bazlı YEREL ortalama. Toplamın türevinden yerel öngörü:
    S(t) = (t/2π)f(L),  f(L) = AL² + α₋₁L + α₀ + ...
    dS/dt = (1/2π)(f + f'),   dN/dt = L/2π
    ⟹ mean(max² | yükseklik t) = (f+f')/L
                                = A·L + (α₋₁+2A) + (α₀+α₋₁)/L + O(1/L²)

Yani 37'de ölçtüğümüz b sabit terimi teorik olarak: b_teori ≈ α₋₁ + 2A
(+ küçük 1/L düzeltmesi). 37 ölçümü: b = 2.776. Tutuyor mu?
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import mpmath as mp

mp.mp.dps = 30
HERE = Path(__file__).resolve().parent

# --- Teorik katsayılar ---
e2 = float(mp.e**2)
g0 = float(mp.euler)          # γ₀
g1 = float(mp.stieltjes(1))   # γ₁ = -0.0728158...
A = (e2 - 5) / 2

alpha_m1 = 5 - e2 - 10*g0 + 2*e2*g0
alpha_0 = 12*g1 - 4*e2*g1 - 5 + e2 + 10*g0 - 2*e2*g0 - 4*g0**2

b_lead = alpha_m1 + 2*A          # yerel sabit terim (1/L'siz)
c_1L = alpha_0 + alpha_m1        # 1/L katsayısı

print("HLP-C Teorem 1 katsayıları:")
print(f"  A = ½(e²−5)      = {A:.6f}")
print(f"  α₋₁              = {alpha_m1:.6f}")
print(f"  α₀               = {alpha_0:.6f}")
print(f"\nYerel ortalama öngörüsü: mean = A·L + {b_lead:.4f} + {c_1L:.4f}/L\n")

# --- Verimiz (36) ile karşılaştırma ---
d = np.load(HERE / "36_T100k.npz")
max2 = d["max_amps"]**2
logt = np.log(d["t_mid"] / (2*np.pi))

n_w = 40
edges = np.geomspace(d["t_mid"][0], d["t_mid"][-1]*1.0001, n_w+1)
rows = []
for i in range(n_w):
    m = (d["t_mid"] >= edges[i]) & (d["t_mid"] < edges[i+1])
    if m.sum() < 50:
        continue
    L = logt[m].mean()
    # pencere içi standart hata (bağımsızlık varsayımıyla, yaklaşık)
    rows.append((L, max2[m].mean(), max2[m].std()/np.sqrt(m.sum()), m.sum()))
w_L, w_mean, w_se, w_n = map(np.array, zip(*rows))

pred = A*w_L + b_lead + c_1L/w_L
resid = w_mean - pred
z = resid / w_se

print("Pencere karşılaştırması (teori öngörüsü sıfır-parametreli!):")
print(f"  ortalama artık        : {resid.mean():+.4f}")
print(f"  artık RMS             : {np.sqrt((resid**2).mean()):.4f}")
print(f"  |z| > 2 pencere sayısı: {(abs(z) > 2).sum()}/{len(z)}")
print(f"  ortalama |z|          : {abs(z).mean():.2f}")

# 37'deki serbest fit ile kıyas
Amat = np.vstack([w_L, np.ones_like(w_L)]).T
coef, *_ = np.linalg.lstsq(Amat*np.sqrt(w_n)[:, None], w_mean*np.sqrt(w_n), rcond=None)
a_fit, b_fit = coef
# teori b'sini aynı L aralığına indirge: b_eff = b_lead + c_1L·mean(1/L)
b_eff = b_lead + c_1L*np.mean(1/w_L)
print(f"\nSabit terim kıyası:")
print(f"  ölçülen b (37 serbest fit) = {b_fit:.4f}")
print(f"  teorik  b_eff              = {b_eff:.4f}")
print(f"  sapma                      = {100*(b_fit/b_eff-1):+.2f}%")

# --- Grafik ---
fig, axes = plt.subplots(1, 2, figsize=(12.5, 5))
ax = axes[0]
ax.errorbar(w_L, w_mean, yerr=w_se, fmt="o", ms=4, c="steelblue",
            label="veri (40 pencere, 100k aralık)")
xs = np.linspace(w_L.min(), w_L.max(), 200)
ax.plot(xs, A*xs + b_lead + c_1L/xs, "g-", lw=1.6,
        label="HLP-C öngörüsü (0 serbest parametre)")
ax.plot(xs, A*xs + b_lead, "g:", lw=1,
        label="… 1/L terimi olmadan")
ax.set_xlabel("L = log(t/2π)")
ax.set_ylabel("mean(max|Z|²)")
ax.set_title("HLP-C Teorem 1 vs bizim veri")
ax.legend(fontsize=9)
ax.grid(alpha=0.3)

ax = axes[1]
ax.errorbar(w_L, resid, yerr=w_se, fmt="o", ms=4, c="firebrick")
ax.axhline(0, color="g", lw=1.4)
ax.set_xlabel("L = log(t/2π)")
ax.set_ylabel("veri − teori")
ax.set_title(f"Artıklar (RMS = {np.sqrt((resid**2).mean()):.3f})")
ax.grid(alpha=0.3)

plt.tight_layout()
out = HERE / "38_hlpc_alt_mertebe.png"
plt.savefig(out, dpi=110)
print(f"\nGrafik: {out.name}")
