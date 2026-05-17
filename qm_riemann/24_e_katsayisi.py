"""
24 — e Katsayısı Doğrulama Testi
====================================

Önceki gözlem (23_mertens_kopru): düşük s (s ≤ 0.4) için
  b(s) / [-2π²·s²/(1-s²)] ≈ 2.72 ≈ e

Test: aynı şey daha büyük N (10⁷) ve daha geniş s aralığında doğru mu?
Eğer evet → asıl hipotez: b(s) ≈ -2π²·e · s²/(1-s²) küçük s için.

Eğer doğrulanırsa: bu **YENİ BİR KAPALI FORM** Mertens correction için
ve paper'daki açık problem 1'in kısmi cevabı.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from math import log, pi, e

def primes_upto(N):
    s = np.ones(N+1, bool); s[:2] = False
    for i in range(2, int(N**0.5)+1):
        if s[i]: s[i*i::i] = False
    return np.flatnonzero(s).astype(np.float64)

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

def sigma_emp(seq, s, N_for_log):
    n = len(seq); Npairs = n*(n-1)//2
    Mp = fast_M_plus(s, seq, Npairs)
    Mm = fast_M_minus(s, seq, Npairs)
    return (1/Mp + 1/Mm - 2) * log(N_for_log)**2

def sigma_cont(s, k=1):
    return -2*k*s**2 / (1 - s**2)

# ============================================================
# N tarama: 10⁵, 10⁶, 10⁷
# ============================================================
Ns = [10**5, 10**6, 10**7]
s_grid = np.array([0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50,
                   0.55, 0.60])

results = {}
for N in Ns:
    print(f"\n=== N = {N} ===")
    primes = primes_upto(N)
    print(f"  {len(primes)} asal üretildi")
    rows = []
    for s in s_grid:
        se = sigma_emp(primes, s, N)
        sc = sigma_cont(s)
        b  = (se - sc) * log(N)
        target = -2*pi**2 * s**2 / (1 - s**2)
        ratio = b / target if abs(target) > 1e-10 else np.nan
        rows.append((s, se, sc, b, ratio))
        print(f"  s={s:.2f}: σ_emp={se:+.4f} σ_cont={sc:+.4f} b={b:+.3f} oran={ratio:.4f}")
    results[N] = rows

# ============================================================
# Karşılaştırma: oran s ve N ile nasıl değişiyor?
# ============================================================
print(f"\n=== ORAN b(s,N) / [-2π²s²/(1-s²)] ===")
print(f"{'s':>5} | {'N=10^5':>10} | {'N=10^6':>10} | {'N=10^7':>10} | {'e':>8}")
for i, s in enumerate(s_grid):
    r5 = results[10**5][i][4]
    r6 = results[10**6][i][4]
    r7 = results[10**7][i][4]
    print(f"{s:>5.2f} | {r5:>10.4f} | {r6:>10.4f} | {r7:>10.4f} | {e:>8.5f}")

# ============================================================
# Düşük s ortalaması — N → ∞ extrapolasyonu
# ============================================================
print(f"\n=== DÜŞÜK s ORTALAMASI (s ≤ 0.3, hipotez sabit ≈ e) ===")
for N in Ns:
    ratios_low = [r[4] for r in results[N] if r[0] <= 0.3]
    print(f"  N={N}: ratio_mean = {np.mean(ratios_low):.5f},  std = {np.std(ratios_low):.5f}")
print(f"  e = {e:.5f}")

# ============================================================
# Plot
# ============================================================
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# (1) Oran vs s, farklı N
ax = axes[0, 0]
colors = ["lightblue", "steelblue", "darkblue"]
for N, c in zip(Ns, colors):
    ratios = [r[4] for r in results[N]]
    ax.plot(s_grid, ratios, "o-", lw=1.5, ms=6, color=c, label=f"N=10^{int(np.log10(N))}")
ax.axhline(e, color="red", lw=1.5, ls="--", label=f"e = {e:.5f}")
ax.axhline(2*pi**2/(2*pi**2), color="green", lw=1, ls=":", alpha=0.6)
ax.set_xlabel("s")
ax.set_ylabel("b(s, N) / [-2π²·s²/(1-s²)]")
ax.set_title("Oran s ile değişimi — düşük s'de e'ye yapışıyor mu?")
ax.legend(fontsize=10)
ax.grid(alpha=0.3)
ax.set_ylim(-2, 5)

# (2) σ_emp / σ_cont — Mertens factor M_1(s, N)
ax = axes[0, 1]
for N, c in zip(Ns, colors):
    ratios = [r[1]/r[2] for r in results[N]]
    ax.plot(s_grid, ratios, "o-", lw=1.5, ms=6, color=c, label=f"N=10^{int(np.log10(N))}")
ax.axhline(e, color="red", lw=1.5, ls="--", label=f"e")
ax.set_xlabel("s")
ax.set_ylabel("σ_emp(s, N) / σ_cont(s)  =  M₁(s, N)")
ax.set_title("Mertens factor M₁(s, N) — sabit mi, s-bağımlı mı?")
ax.legend(fontsize=10)
ax.grid(alpha=0.3)
ax.set_ylim(0, 4)

# (3) Düşük-s yakınlığı (s=0.1, 0.2, 0.3) ortalama oran vs N
ax = axes[1, 0]
mean_ratios = []
for N in Ns:
    rs = [r[4] for r in results[N] if r[0] <= 0.3]
    mean_ratios.append(np.mean(rs))
ax.semilogx(Ns, mean_ratios, "o-", lw=2, ms=10, color="purple")
ax.axhline(e, color="red", lw=1.5, ls="--", label=f"e = {e:.5f}")
ax.set_xlabel("N")
ax.set_ylabel("ortalama oran (s ≤ 0.3)")
ax.set_title("N artarken oran → e mi?")
ax.legend(fontsize=10)
ax.grid(alpha=0.3, which="both")

# (4) Eğer hipotez doğruysa: b(s) - (-2π²·e · s²/(1-s²)) artığı
ax = axes[1, 1]
for N, c in zip(Ns, colors):
    b_vals = [r[3] for r in results[N]]
    predict = -2*pi**2 * e * s_grid**2 / (1 - s_grid**2)
    residual = np.array(b_vals) - predict
    ax.plot(s_grid, residual, "o-", lw=1.5, ms=6, color=c, label=f"N=10^{int(np.log10(N))}")
ax.axhline(0, color="black", lw=0.7)
ax.set_xlabel("s")
ax.set_ylabel("b(s) − (-2π²·e·s²/(1-s²))")
ax.set_title("Hipotez artığı — eğer 0'a yakınsa hipotez doğru")
ax.legend(fontsize=10)
ax.grid(alpha=0.3)

plt.suptitle("24 — Hipotez: b(s) = -2π²·e · s²/(1-s²) küçük s için\n"
             "(önceki gözlem: oran ≈ e = 2.71828)",
             fontsize=13, fontweight="bold", y=0.998)
plt.tight_layout()

out = Path(__file__).parent / "24_e_katsayisi.png"
plt.savefig(out, dpi=130, bbox_inches="tight")
print(f"\nKaydedildi: {out}")
