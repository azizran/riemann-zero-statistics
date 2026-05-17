"""
25 — Open Problem 2: Polynomial-Density Diziler
=================================================

Paper Open Problem 2: kareler, küpler için σ_S(s) → 0 ama "non-trivial
subleading structure" var. Henüz karakterize edilmemiş.

Test:
  - Kareler {n² : n ≤ M=√N}
  - Küpler {n³ : n ≤ M=N^(1/3)}
  - σ_S(s, N) için N tarama
  - Yakınsama hızı: σ · N^p → c(s) hangi p için?
  - c(s) için kapalı form var mı?
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from math import log

def fast_M_plus(s, seq, Npairs):
    cum = 0.0; S = 0.0
    for x in seq:
        xs = x**s
        if cum > 0: S += cum / xs
        cum += xs
    return S / Npairs

def fast_M_minus(s, seq, Npairs):
    cum_inv = 0.0; S = 0.0
    for q in seq:
        qs = q**s
        if cum_inv > 0: S += qs * cum_inv
        cum_inv += 1.0 / qs
    return S / Npairs

def sigma_emp(seq, s, N_for_log=None):
    """σ_S(s) = (1/M_+ + 1/M_- - 2) · log²N (paper tanımı)
       AMA polynomial için log²N normalize yanlış olabilir.
       Burada log faktörü olmadan ham deviation döndür."""
    n = len(seq); Npairs = n*(n-1)//2
    Mp = fast_M_plus(s, seq, Npairs)
    Mm = fast_M_minus(s, seq, Npairs)
    raw = 1/Mp + 1/Mm - 2     # ham fark
    return raw, Mp, Mm

# ============================================================
# Kareler ve küpler için seqler üret
# ============================================================
print("Diziler hazırlanıyor...")
N_max = 10**8

squares_dict = {}
cubes_dict = {}
for N in [10**5, 10**6, 10**7, 10**8]:
    M_sq = int(np.sqrt(N))
    M_cu = int(N**(1/3))
    squares_dict[N] = np.arange(2, M_sq+1).astype(np.float64)**2
    cubes_dict[N]   = np.arange(2, M_cu+1).astype(np.float64)**3
    print(f"  N=10^{int(np.log10(N))}: kareler {len(squares_dict[N])}, küpler {len(cubes_dict[N])}")

# ============================================================
# σ_raw hesapla, s tarama
# ============================================================
s_grid = np.array([0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45])

print(f"\n=== KARELER ===")
print(f"{'s':>5} | {'N=10⁵':>15} | {'N=10⁶':>15} | {'N=10⁷':>15} | {'N=10⁸':>15}")
sq_results = {}
for s in s_grid:
    row = [s]
    for N in [10**5, 10**6, 10**7, 10**8]:
        raw, Mp, Mm = sigma_emp(squares_dict[N], s)
        row.append(raw)
        sq_results[(s, N)] = raw
    print(f"{s:>5.2f} | {row[1]:>+15.8f} | {row[2]:>+15.8f} | {row[3]:>+15.8f} | {row[4]:>+15.8f}")

print(f"\n=== KÜPLER ===")
print(f"{'s':>5} | {'N=10⁵':>15} | {'N=10⁶':>15} | {'N=10⁷':>15} | {'N=10⁸':>15}")
cu_results = {}
for s in s_grid:
    row = [s]
    for N in [10**5, 10**6, 10**7, 10**8]:
        raw, Mp, Mm = sigma_emp(cubes_dict[N], s)
        row.append(raw)
        cu_results[(s, N)] = raw
    print(f"{s:>5.2f} | {row[1]:>+15.8f} | {row[2]:>+15.8f} | {row[3]:>+15.8f} | {row[4]:>+15.8f}")

# ============================================================
# Yakınsama hızı: σ_raw · M^p → c(s) hangi p?
# ============================================================
# M kareler için √N, küpler için N^(1/3)
# Log-log fit: log|σ_raw| vs log M
print(f"\n=== KARELER YAKINSAMA HIZI (log|σ| vs log M) ===")
print(f"{'s':>5} | {'eğim p (lineer fit)':>22} | {'tahmin formu':>30}")
Ns = [10**5, 10**6, 10**7, 10**8]
Ms_sq = [int(np.sqrt(N)) for N in Ns]
for s in s_grid:
    vals = [abs(sq_results[(s, N)]) for N in Ns]
    log_M = np.log(Ms_sq)
    log_v = np.log(vals)
    slope, _ = np.polyfit(log_M, log_v, 1)
    print(f"{s:>5.2f} | {slope:>22.4f} | M^{slope:.2f}")

print(f"\n=== KÜPLER YAKINSAMA HIZI ===")
Ms_cu = [int(N**(1/3)) for N in Ns]
for s in s_grid:
    vals = [abs(cu_results[(s, N)]) for N in Ns]
    log_M = np.log(Ms_cu)
    log_v = np.log(vals)
    slope, _ = np.polyfit(log_M, log_v, 1)
    print(f"{s:>5.2f} | {slope:>22.4f} | M^{slope:.2f}")

# ============================================================
# Eğer yakınsama M^-α ise, α(s) ne? Sabit mi yoksa s-bağımlı?
# Sabit ise: σ_raw · M^α → c(s), bu c(s)'in kapalı formunu ara
# ============================================================
print(f"\n=== ÖZGÜN YENİ SORU: yakınsama üssü α(s) ===")
print(f"Kareler için, M=√N: eğer α=1 (σ ~ 1/M), σ · M ~ c(s)")
print(f"Bu c(s)'ye bakalım:\n")
print(f"{'s':>5} | {'σ_raw · M':>15} {'(N=10⁵)':>10} {'(N=10⁶)':>10} {'(N=10⁷)':>10} {'(N=10⁸)':>10}")
for s in s_grid:
    vals = [sq_results[(s, N)] * Ms_sq[i] for i, N in enumerate(Ns)]
    print(f"{s:>5.2f} | {vals[0]:>+15.4f} {'':>10} {vals[1]:>+10.4f} {vals[2]:>+10.4f} {vals[3]:>+10.4f}")

# Ya M^2 normalize ise?
print(f"\n{'s':>5} | {'σ_raw · M²':>15} {'(N=10⁵)':>10} {'(N=10⁶)':>10} {'(N=10⁷)':>10} {'(N=10⁸)':>10}")
for s in s_grid:
    vals = [sq_results[(s, N)] * Ms_sq[i]**2 for i, N in enumerate(Ns)]
    print(f"{s:>5.2f} | {vals[0]:>+15.4f} {'':>10} {vals[1]:>+10.2f} {vals[2]:>+10.2f} {vals[3]:>+10.2f}")

# ============================================================
# Plot
# ============================================================
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# (1) Kareler σ_raw vs s, farklı N
ax = axes[0, 0]
colors = ["lightblue", "steelblue", "navy", "midnightblue"]
for N, c in zip(Ns, colors):
    vals = [sq_results[(s, N)] for s in s_grid]
    ax.plot(s_grid, vals, "o-", lw=1.5, ms=6, color=c,
            label=f"N=10^{int(np.log10(N))}")
ax.axhline(0, color="red", lw=1, ls="--", label="σ→0 (k=0 sınıfı)")
ax.set_xlabel("s")
ax.set_ylabel("σ_raw = 1/M_+ + 1/M_- - 2")
ax.set_title("KARELER: σ_raw vs s, farklı N\nN büyüdükçe sıfıra yakınsıyor (k=0 doğrulanır)")
ax.legend(fontsize=9)
ax.grid(alpha=0.3)

# (2) log-log yakınsama
ax = axes[0, 1]
for s in [0.1, 0.2, 0.3, 0.4]:
    vals = [abs(sq_results[(s, N)]) for N in Ns]
    ax.loglog(Ms_sq, vals, "o-", lw=1.2, ms=6, label=f"s={s}")
ax.set_xlabel("M = √N")
ax.set_ylabel("|σ_raw|")
ax.set_title("Kareler — log-log yakınsama (eğim = -α)")
ax.legend(fontsize=10)
ax.grid(alpha=0.3, which="both")

# (3) σ · M (eğer α=1 ise sabit olmalı)
ax = axes[1, 0]
for s in [0.1, 0.2, 0.3, 0.4]:
    vals = [sq_results[(s, N)] * Ms_sq[i] for i, N in enumerate(Ns)]
    ax.semilogx(Ms_sq, vals, "o-", lw=1.5, ms=6, label=f"s={s}")
ax.axhline(0, color="black", lw=0.5)
ax.set_xlabel("M")
ax.set_ylabel("σ_raw · M")
ax.set_title("Eğer α=1: σ·M → c(s) sabit olmalı")
ax.legend(fontsize=10)
ax.grid(alpha=0.3)

# (4) σ · M² (eğer α=2 ise)
ax = axes[1, 1]
for s in [0.1, 0.2, 0.3, 0.4]:
    vals = [sq_results[(s, N)] * Ms_sq[i]**2 for i, N in enumerate(Ns)]
    ax.semilogx(Ms_sq, vals, "o-", lw=1.5, ms=6, label=f"s={s}")
ax.axhline(0, color="black", lw=0.5)
ax.set_xlabel("M")
ax.set_ylabel("σ_raw · M²")
ax.set_title("Eğer α=2: σ·M² → c(s) sabit olmalı")
ax.legend(fontsize=10)
ax.grid(alpha=0.3)

plt.suptitle("25 — Polynomial Density (Kareler) — Subleading Structure Aranıyor",
             fontsize=13, fontweight="bold", y=0.998)
plt.tight_layout()

out = Path(__file__).parent / "25_polynomial_density.png"
plt.savefig(out, dpi=130, bbox_inches="tight")
print(f"\nKaydedildi: {out}")
