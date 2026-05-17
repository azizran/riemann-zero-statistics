"""
23 — Mertens Katsayısı → Riemann ζ Köprüsü
================================================

Polygon-π paper'ındaki açık problem (Bölüm 4.3 + 5.3):
  σ_S(s) = σ_S^cont(s) + σ_S^Mertens(s, N)
  σ_S^cont(s) = -2k·s²/(1-s²)
  σ_S^Mertens(s, N) = b(s)/log N  (b(s)'nin kapalı formu açık)

Tek bilinen ipucu: s=0.5'te b ≈ -19.99 ≈ -2π² (paper).

BUGÜN TEST: b(s)'i birden çok s değeri için hesapla, kapalı form ara,
           Riemann ζ değerleriyle bağlantı kur.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import mpmath as mp
from math import log

mp.mp.dps = 25

# ============================================================
# Asalları üret (sieve)
# ============================================================
def primes_upto(N):
    s = np.ones(N+1, bool); s[:2] = False
    for i in range(2, int(N**0.5)+1):
        if s[i]: s[i*i::i] = False
    return np.flatnonzero(s).astype(np.float64)

print("Asallar üretiliyor (N = 10^6, yaklaşık 78498 asal)...")
N_max = 10**6
primes = primes_upto(N_max)
print(f"Asal sayısı: {len(primes)}")

# ============================================================
# Hızlı M_+(s) ve M_-(s) hesabı (paper'dan O(N))
# ============================================================
def fast_M_plus(s, seq, Npairs):
    """M_+(s) = (1/Npairs) Σ_{a<b} (a/b)^s"""
    cum = 0.0; S = 0.0
    for x in seq:
        xs = x**s
        if cum > 0: S += cum / xs
        cum += xs
    return S / Npairs

def fast_M_minus(s, seq, Npairs):
    """M_-(s) = (1/Npairs) Σ_{a<b} (b/a)^s"""
    cum_inv = 0.0; S = 0.0
    for q in seq:
        qs = q**s
        if cum_inv > 0: S += qs * cum_inv
        cum_inv += 1.0 / qs
    return S / Npairs

def sigma_emp(seq, s, N_for_log):
    n = len(seq); Npairs = n*(n-1)//2
    Mp = fast_M_plus(s, seq, Npairs)
    Mm = fast_M_minus(s, seq, Npairs)
    return (1/Mp + 1/Mm - 2) * log(N_for_log)**2

def sigma_cont(s, k=1):
    """Paper Teorem 2"""
    return -2*k*s**2 / (1 - s**2)

# ============================================================
# s grid'inde σ_emp hesapla
# ============================================================
s_grid = np.array([0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5,
                   0.55, 0.6, 0.65, 0.7, 0.75, 0.8])

print(f"\nσ_emp(s) hesaplanıyor, {len(s_grid)} s değeri için...")
sigma_e = []
sigma_c = []
b_values = []
log_N = log(N_max)

for i, s in enumerate(s_grid):
    se = sigma_emp(primes, s, N_max)
    sc = sigma_cont(s)
    b = (se - sc) * log_N
    sigma_e.append(se)
    sigma_c.append(sc)
    b_values.append(b)
    print(f"  s={s:.2f}:  σ_emp={se:+.4f}  σ_cont={sc:+.4f}  b={b:+.3f}")

sigma_e = np.array(sigma_e)
sigma_c = np.array(sigma_c)
b_values = np.array(b_values)

# ============================================================
# b(s) için kapalı form adayları
# ============================================================
print(f"\n=== b(s) için aday kapalı formlar ===")
print(f"{'s':>5} {'b_emp':>10} {'-2π²s²/(1-s²)':>15} {'-2γs²/(1-s²)':>15} {'-2π²s²':>10}")
gamma = 0.5772156649
adaylar = {
    "-2π²·s²/(1-s²)":  -2*np.pi**2 * s_grid**2 / (1 - s_grid**2),
    "-2γ·s²/(1-s²)":   -2*gamma * s_grid**2 / (1 - s_grid**2),
    "-2π²·s²":         -2*np.pi**2 * s_grid**2,
    "-2·s²·ζ(2)":      -2 * s_grid**2 * np.pi**2/6,
    "ratio·σ_cont":    None,   # b/σ_cont — sabit mi?
}

# Tablo: bilgileri görsel olarak
print(f"\n{'s':>5} | {'b_emp':>10} | {'-2π²s²/(1-s²)':>15} | {'ratio':>10}")
for i, s in enumerate(s_grid):
    pred = -2*np.pi**2 * s**2 / (1 - s**2)
    ratio = b_values[i]/pred if abs(pred) > 1e-10 else np.nan
    print(f"{s:>5.2f} | {b_values[i]:>+10.3f} | {pred:>+15.3f} | {ratio:>10.3f}")

# Ana karşılaştırma: b(s) / σ_cont(s) sabit mi?
ratio_b_sc = b_values / sigma_c
print(f"\n=== b(s) / σ_cont(s) oranı (eğer sabitse → b = c · σ_cont) ===")
for i, s in enumerate(s_grid):
    print(f"  s={s:.2f}:  b/σ_cont = {ratio_b_sc[i]:+.3f}")
print(f"\nOrtalama oran: {ratio_b_sc.mean():.3f}, std: {ratio_b_sc.std():.3f}")
print(f"Eğer std küçükse → b(s) = c · σ_cont(s), basit formül!")

# Bu sabit ne olabilir? π², 2π², 2π²/log...
const = ratio_b_sc.mean()
print(f"\nSabit ≈ {const:.3f}")
print(f"  π²    = {np.pi**2:.3f}")
print(f"  2π²   = {2*np.pi**2:.3f}")
print(f"  e²    = {np.e**2:.3f}")
print(f"  3π    = {3*np.pi:.3f}")

# ============================================================
# Plot
# ============================================================
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# (1) σ_emp vs σ_cont
ax = axes[0, 0]
ax.plot(s_grid, sigma_e, "ro-", lw=1.5, ms=7, label="σ_emp (N=10⁶)")
ax.plot(s_grid, sigma_c, "b-", lw=2, label="σ_cont = -2s²/(1-s²) (teorem)")
ax.set_xlabel("s")
ax.set_ylabel("σ_prime(s)")
ax.set_title("σ_prime: empirik (N=10⁶) vs teorem leading")
ax.legend(fontsize=10)
ax.grid(alpha=0.3)

# (2) Mertens deviation b(s)
ax = axes[0, 1]
ax.plot(s_grid, b_values, "ko-", lw=1.5, ms=7, label="b_emp(s)")
# Aday kapalı formlar
for label, vals in adaylar.items():
    if vals is None: continue
    ax.plot(s_grid, vals, "--", lw=1.2, alpha=0.7, label=label)
ax.set_xlabel("s")
ax.set_ylabel("b(s) = (σ_emp - σ_cont)·log N")
ax.set_title(f"Mertens deviation b(s) — kapalı formu var mı?")
ax.legend(fontsize=8)
ax.grid(alpha=0.3)

# (3) b(s) / σ_cont(s) oranı — sabit mi?
ax = axes[1, 0]
ax.plot(s_grid, ratio_b_sc, "go-", lw=1.5, ms=7)
ax.axhline(ratio_b_sc.mean(), color="red", lw=1.2, ls="--",
           label=f"ortalama = {ratio_b_sc.mean():.3f}")
ax.set_xlabel("s")
ax.set_ylabel("b(s) / σ_cont(s)")
ax.set_title(f"b/σ_cont oranı — sabitse: b(s) = c · σ_cont(s)\n"
             f"std = {ratio_b_sc.std():.3f}")
ax.legend(fontsize=10)
ax.grid(alpha=0.3)

# (4) Riemann ζ ile karşılaştırma
ax = axes[1, 1]
# log|ζ(1+s)·ζ(1-s)| — Riemann tarafı çiftli ifade
zeta_combo = []
for s in s_grid:
    z1 = abs(complex(mp.zeta(mp.mpc(1+s, 0.0001))))
    z2 = abs(complex(mp.zeta(mp.mpc(1-s, 0.0001))))
    zeta_combo.append(log(z1 * z2))
zeta_combo = np.array(zeta_combo)
ax.plot(s_grid, b_values, "ko-", lw=1.5, ms=6, label="b_emp(s)")
ax.plot(s_grid, zeta_combo, "r--", lw=1.2, label="log|ζ(1+s)·ζ(1-s)|")
# Normalize bir kıvrım
ax.plot(s_grid, -2*np.pi**2 * s_grid**2/(1-s_grid**2), "b:", lw=1.5,
        label="-2π²·s²/(1-s²)")
ax.set_xlabel("s")
ax.set_ylabel("değer")
ax.set_title("Riemann ζ ile karşılaştırma")
ax.legend(fontsize=9)
ax.grid(alpha=0.3)

plt.suptitle("23 — Mertens Katsayısı b(s): kapalı form ve Riemann ζ köprüsü arayışı",
             fontsize=13, fontweight="bold", y=0.998)
plt.tight_layout()

out = Path(__file__).parent / "23_mertens_kopru.png"
plt.savefig(out, dpi=130, bbox_inches="tight")
print(f"\nKaydedildi: {out}")
