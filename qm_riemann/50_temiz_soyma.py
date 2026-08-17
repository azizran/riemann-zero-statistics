"""
50 — TEMİZ SOYMA: ENDOJENLİKSİZ, PLASEBOLU, BİRLEŞİK (17 Ağustos 2026)
=======================================================================

49'un şüphesi: t_pk fazları endojen (tepe, dalganın sırtını arar → katsayı
şişer). Düzeltme: TÜM fazlar tmid'de (sıfır konumlarından belirlenir,
genlik ölçümünden bağımsız → egzojen). Z hesabı gerekmez, saklı veri yeter.

Tutarlı soyma (49-2 matematiği, temiz fazlarla):
  log g̃ → asal terimleri çıkar (v kanalı)          → g̃_res
  log ã → KONTROLSÜZ asal terimleri çıkar (u+f'·v) → ã_res
  Model doğruysa (g̃_res, ã_res) ≈ (g̃_C, m̃_C) = saf matris çifti.

PLASEBO: aynı makine, 6 asal-olmayan frekansla → r KIPIRDAMAMALI.
Kıpırdarsa soyma mekanik olarak r şişiriyor demektir ve asal sonucu da çöp.

KARAR KRİTERİ:
  soyulmuş r CUE aralığına iner VE Pearson/Spearman N_eff'leri buluşursa
  → HALKA KAPANDI: ζ = CUE(N*) + ölçülmüş asal dalgaları
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import pearsonr, spearmanr
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
A = (np.e**2 - 5) / 2
B0, B1 = 2.7580, -0.0543
PRIMES = [2, 3, 5, 7, 11, 13]
PLACEBO = [2.31, 3.85, 5.55, 7.77, 11.3, 13.7]  # asal değil, asal kuvveti değil

FIT_P = (0.5238, 3.703, -8.486)
FIT_S = (0.7704, 1.800, -4.099)

def inv_fit(fit, r):
    f = lambda N: fit[0] + fit[1] / N + fit[2] / N**2 - r
    try:
        return brentq(f, 5.0, 80)
    except ValueError:
        return np.nan

def wave_cols(tvals, freqs):
    cols = []
    for q in freqs:
        cols += [np.cos(tvals * np.log(q)), np.sin(tvals * np.log(q))]
    return cols

def strip(y, tvals, freqs):
    """y'den freq bileşenlerini regresyonla çıkar; (artık, çıkarılan varyans)."""
    C = np.vstack([np.ones_like(y)] + wave_cols(tvals, freqs)).T
    b = np.linalg.lstsq(C, y, rcond=None)[0]
    part = C[:, 1:] @ b[1:]
    return y - part, float(part.var()), b[1:]

d41 = np.load(HERE / "41_bigT_windows.npz")
keys = sorted({k.split("_")[1] for k in d41.files}, key=lambda s: int(s[:-1]))

rows = []
print(f"{'L':>6} | {'r_P ham':>8} {'soyul':>7} {'plasebo':>8} | "
      f"{'r_S ham':>8} {'soyul':>7} {'plasebo':>8} | {'ΔVar_a':>7} {'ΔVar_g':>7}")
for wkey in keys:
    gaps = d41[f"gaps_{wkey}"]; tmid = d41[f"tmid_{wkey}"]; amps = d41[f"amps_{wkey}"]
    Lw = np.log(tmid / TWO_PI); L = float(Lw.mean())
    g_u = gaps * Lw / TWO_PI
    a_u = amps / np.sqrt(A * Lw + B0 + B1 / Lw)
    a_u /= np.sqrt((a_u**2).mean())
    y_a, y_g = np.log(a_u), np.log(g_u)

    rP0, _ = pearsonr(g_u, a_u); rS0, _ = spearmanr(g_u, a_u)

    # ASAL soyma (tmid fazları — egzojen)
    ya_r, va_removed, coef_a = strip(y_a, tmid, PRIMES)
    yg_r, vg_removed, coef_g = strip(y_g, tmid, PRIMES)
    a_res = np.exp(ya_r); a_res /= np.sqrt((a_res**2).mean())
    g_res = np.exp(yg_r); g_res *= g_u.mean() / g_res.mean()
    rP1, _ = pearsonr(g_res, a_res); rS1, _ = spearmanr(g_res, a_res)

    # PLASEBO soyma
    ya_p, _, _ = strip(y_a, tmid, PLACEBO)
    yg_p, _, _ = strip(y_g, tmid, PLACEBO)
    a_pl = np.exp(ya_p); a_pl /= np.sqrt((a_pl**2).mean())
    g_pl = np.exp(yg_p); g_pl *= g_u.mean() / g_pl.mean()
    rP2, _ = pearsonr(g_pl, a_pl); rS2, _ = spearmanr(g_pl, a_pl)

    rows.append((L, rP0, rP1, rP2, rS0, rS1, rS2, va_removed, vg_removed,
                 coef_a, coef_g))
    print(f"{L:>6.2f} | {rP0:>8.4f} {rP1:>7.4f} {rP2:>8.4f} | "
          f"{rS0:>8.4f} {rS1:>7.4f} {rS2:>8.4f} | {va_removed:>7.4f} {vg_removed:>7.4f}")

# kanal katsayıları (tmid bazlı; 45/47 ile kıyas için)
print("\ntmid-bazlı kanal katsayıları (L=12.45; katsayı/p^-1/2):")
La, ca, cg = rows[-1][0], rows[-1][9], rows[-1][10]
for i, p in enumerate(PRIMES):
    u_tot = np.hypot(ca[2 * i], ca[2 * i + 1]) / p**-0.5
    v_ch = np.hypot(cg[2 * i], cg[2 * i + 1]) / p**-0.5
    print(f"  p={p:>2}: |u_toplam|={u_tot:.3f}  |v|={v_ch:.3f}  "
          f"(45'in t_pk w'si: {[0.782,0.674,0.544,0.454,0.345,0.311][i]:.3f})")

# N_eff buluşma testi
arr = np.array([r[:9] for r in rows])
print(f"\n{'L':>6} | {'N_eff−L (soyulmuş)':^24}")
print(f"{'':>6} | {'Pearson':>10} {'Spearman':>10}")
NP_list, NS_list = [], []
for L, rP0, rP1, rP2, rS0, rS1, rS2, va, vg in arr:
    NP = inv_fit(FIT_P, rP1); NS = inv_fit(FIT_S, rS1)
    NP_list.append(NP - L if not np.isnan(NP) else np.nan)
    NS_list.append(NS - L if not np.isnan(NS) else np.nan)
    fp = f"{NP-L:+10.2f}" if not np.isnan(NP) else f"{'aralık dışı':>10}"
    fs = f"{NS-L:+10.2f}" if not np.isnan(NS) else f"{'aralık dışı':>10}"
    print(f"{L:>6.2f} | {fp} {fs}")

NP_a, NS_a = np.nanmean(NP_list), np.nanmean(NS_list)
print(f"\nSOYULMUŞ ortalama: N_eff−L  Pearson {NP_a:+.3f}  Spearman {NS_a:+.3f}"
      f"  → fark {NS_a-NP_a:+.3f}")
print("HAM (42b): Pearson +0.93, Spearman +1.86 → fark +0.93 idi")
print(f"PLASEBO kontrol: r değişimi P {np.abs(arr[:,3]-arr[:,1]).max():.4f}, "
      f"S {np.abs(arr[:,6]-arr[:,4]).max():.4f} (≈0 olmalı)")

# grafik
fig, ax = plt.subplots(figsize=(8.5, 5.4))
ax.plot(arr[:, 0], [inv_fit(FIT_P, r) - L for (L, r) in zip(arr[:, 0], arr[:, 1])],
        "o-", c="firebrick", label="ham Pearson")
ax.plot(arr[:, 0], [inv_fit(FIT_S, r) - L for (L, r) in zip(arr[:, 0], arr[:, 4])],
        "s-", c="darkorange", label="ham Spearman")
ax.plot(arr[:, 0], NP_list, "o--", c="steelblue", label="soyulmuş Pearson")
ax.plot(arr[:, 0], NS_list, "s--", c="teal", label="soyulmuş Spearman")
ax.axhline(0, color="k", lw=0.8)
ax.set_xlabel("L"); ax.set_ylabel("N_eff − L")
ax.set_title("Temiz soyma: iki metrik buluşuyor mu?")
ax.legend(fontsize=9); ax.grid(alpha=0.3)
plt.tight_layout()
out = HERE / "50_temiz_soyma.png"
plt.savefig(out, dpi=110)
print(f"\nGrafik: {out.name}")
