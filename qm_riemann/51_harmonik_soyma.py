"""
51 — HARMONİK SOYMA: ÇEKİRDEĞİN GERÇEK SIKILIĞI (17 Ağustos 2026)
===================================================================

50'de yalnız p≤13, k=1 soyulmuştu → çekirdek r≈0.94/0.97 (alt sınır).
Şimdi explicit formulanın tüm anlamlı terimleri:

  katman 1: p ≤ 13, k=1              (50 ile aynı — replikasyon)
  katman 2: + asal kuvvetleri p^k     (4,8,9,16,25,27,32,49; ağırlık p^{-k/2}/k)
  katman 3: + kuyruk asalları 17..47  (τ-yasasına göre hâlâ iletimde)

Plasebo: 23 asal-olmayan frekans (katman 3 ile aynı sütun sayısı) —
overfit tabanını ölçer.

Bonus: kuvvet katsayıları / (p^{-k/2}/k) → toplam kuralı harmoniklerde
de u ≈ 1 mi?
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import pearsonr, spearmanr

rng = np.random.default_rng(51)
HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
A = (np.e**2 - 5) / 2
B0, B1 = 2.7580, -0.0543

L1 = [2, 3, 5, 7, 11, 13]
POWERS = [4, 8, 9, 16, 25, 27, 32, 49]           # p^k
TAIL = [17, 19, 23, 29, 31, 37, 41, 43, 47]
L2 = L1 + POWERS
L3 = L2 + TAIL

# explicit formula ağırlığı: q = p^k için p^{-k/2}/k
def ef_weight(q):
    for p in [2, 3, 5, 7]:
        k = round(np.log(q) / np.log(p))
        if k >= 2 and abs(p**k - q) < 0.5:
            return p**(-k / 2) / k
    return q**-0.5  # k=1

# plasebo: kullanılan tüm frekanslardan uzak 23 sahte frekans
used = np.log(np.array(L3, dtype=float))
plc = []
while len(plc) < len(L3):
    cand = np.exp(rng.uniform(np.log(np.log(2)), np.log(np.log(50))))
    if np.abs(used - cand).min() > 0.02 and all(abs(c - cand) > 0.02 for c in plc):
        plc.append(cand)
PLACEBO_FREQS = np.array(plc)  # doğrudan frekans (log q eşdeğeri)

def wave_cols_freq(tvals, freqs):
    cols = []
    for om in freqs:
        cols += [np.cos(tvals * om), np.sin(tvals * om)]
    return cols

def strip_freqs(y, tvals, freqs):
    C = np.vstack([np.ones_like(y)] + wave_cols_freq(tvals, freqs)).T
    b = np.linalg.lstsq(C, y, rcond=None)[0]
    part = C[:, 1:] @ b[1:]
    return y - part, float(part.var()), b[1:]

d41 = np.load(HERE / "41_bigT_windows.npz")
keys = sorted({k.split("_")[1] for k in d41.files}, key=lambda s: int(s[:-1]))

LEVELS = [("katman1 (p≤13)", np.log(np.array(L1, float))),
          ("katman2 (+p^k)", np.log(np.array(L2, float))),
          ("katman3 (+kuyruk)", np.log(np.array(L3, float))),
          ("plasebo (23 frek.)", PLACEBO_FREQS)]

print(f"{'L':>6} | {'ham r_P':>8} " + " ".join(f"{n.split()[0]:>10}" for n, _ in LEVELS)
      + f" | {'ham r_S':>8} " + " ".join(f"{n.split()[0]:>10}" for n, _ in LEVELS))
results = []
for wkey in keys:
    gaps = d41[f"gaps_{wkey}"]; tmid = d41[f"tmid_{wkey}"]; amps = d41[f"amps_{wkey}"]
    Lw = np.log(tmid / TWO_PI); L = float(Lw.mean())
    g_u = gaps * Lw / TWO_PI
    a_u = amps / np.sqrt(A * Lw + B0 + B1 / Lw)
    a_u /= np.sqrt((a_u**2).mean())
    y_a, y_g = np.log(a_u), np.log(g_u)
    rP0, _ = pearsonr(g_u, a_u); rS0, _ = spearmanr(g_u, a_u)

    rowP, rowS, va_list, coefs = [], [], [], None
    for name, freqs in LEVELS:
        ya, va, ca = strip_freqs(y_a, tmid, freqs)
        yg, vg, _ = strip_freqs(y_g, tmid, freqs)
        ar = np.exp(ya); ar /= np.sqrt((ar**2).mean())
        gr = np.exp(yg); gr *= g_u.mean() / gr.mean()
        rP, _ = pearsonr(gr, ar); rS, _ = spearmanr(gr, ar)
        rowP.append(rP); rowS.append(rS); va_list.append(va)
        if name.startswith("katman3"):
            coefs = ca
    results.append((L, rP0, rowP, rS0, rowS, va_list, coefs))
    print(f"{L:>6.2f} | {rP0:>8.4f} " + " ".join(f"{r:>10.4f}" for r in rowP)
          + f" | {rS0:>8.4f} " + " ".join(f"{r:>10.4f}" for r in rowS))

# toplam kuralı: katman-3 katsayıları / explicit-formula ağırlığı (L=12.45)
print("\nToplam kuralı kontrolü (L=12.45, katsayı/EF-ağırlığı):")
L_, _, _, _, _, _, coefs = results[-1]
for i, q in enumerate(L3):
    u = np.hypot(coefs[2 * i], coefs[2 * i + 1]) / ef_weight(q)
    tag = "kuvvet" if q in POWERS else ("kuyruk" if q in TAIL else "temel")
    print(f"  q={q:>3} ({tag:>6}): |u| = {u:.3f}")

# varyans bütçesi
print("\nGenlikten soyulan log-varyans (L=12.45):")
for (name, _), va in zip(LEVELS, results[-1][5]):
    print(f"  {name:>20}: {va:.4f}")

# grafik: r_P ve r_S katman ilerlemesi (pencere ort.)
arr_P = np.array([[r0] + row for _, r0, row, _, _, _, _ in results])
arr_S = np.array([[r0] + row for _, _, _, r0, row, _, _ in results])
labels = ["ham"] + [n for n, _ in LEVELS]
fig, ax = plt.subplots(figsize=(9, 5.4))
x = np.arange(len(labels))
ax.errorbar(x - 0.05, arr_P.mean(axis=0), yerr=arr_P.std(axis=0), fmt="o-",
            c="firebrick", label="Pearson (6 pencere ort.)")
ax.errorbar(x + 0.05, arr_S.mean(axis=0), yerr=arr_S.std(axis=0), fmt="s-",
            c="teal", label="Spearman")
ax.set_xticks(x); ax.set_xticklabels(labels, rotation=15, fontsize=9)
ax.set_ylabel("soyulmuş r")
ax.set_title("Çekirdeğin sıkılığı: soyma katmanları")
ax.legend(); ax.grid(alpha=0.3)
plt.tight_layout()
out = HERE / "51_harmonik_soyma.png"
plt.savefig(out, dpi=110)
print(f"\nGrafik: {out.name}")
