"""
49 — ASAL SOYMA: KALAN SAF CUE MU? (17 Ağustos 2026)
======================================================

48'in dersi: generatif yön (CUE+asal→ζ) kırılgan. Tersi sağlam:
ζ verisinden ÖLÇÜLMÜŞ asal bileşenlerini regresyonla soy
(genlik kanalı 45, boşluk kanalı 47 — sadece asal terimleri çıkar),
kalıntı çiftin (g̃_res, ã_res) istatistiğine bak.

BÜYÜK ÖNGÖRÜ: iki-sabit anomalisi (N_eff−L: Pearson +0.93 ≠ Spearman
+1.86) asalların eseriyse, soyma sonrası İKİ METRİK AYNI N_eff'TE
BULUŞMALI. Buluşursa: "ζ = tek bir sonlu-CUE + ölçülmüş asal dalgaları"
cümlesi kanıtlanır ve halka kapanır.

CUE r(N) referansı: 42b pürüzsüz fitleri (r∞ + c₁/N + c₂/N²).
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

# 42b pürüzsüz CUE fitleri
FIT_P = (0.5238, 3.703, -8.486)   # Pearson : r∞, c1, c2
FIT_S = (0.7704, 1.800, -4.099)   # Spearman

def inv_fit(fit, r):
    f = lambda N: fit[0] + fit[1] / N + fit[2] / N**2 - r
    try:
        return brentq(f, 5.0, 80)
    except ValueError:
        return np.nan  # fit aralığı dışı — r tablodaki tüm CUE değerlerinin dışında

def theta(t):
    return t / 2 * np.log(t / TWO_PI) - t / 2 - np.pi / 8 + 1 / (48 * t) + 7 / (5760 * t**3)

def Z_rs(t, chunk=20000):
    t = np.asarray(t, dtype=np.float64)
    out = np.empty_like(t)
    for s in range(0, len(t), chunk):
        tt = t[s:s + chunk]
        a = np.sqrt(tt / TWO_PI)
        N = a.astype(np.int64)
        th = theta(tt)
        z = np.zeros_like(tt)
        for Nv in np.unique(N):
            m = N == Nv
            n = np.arange(1, Nv + 1)
            ph = th[m, None] - tt[m, None] * np.log(n)[None, :]
            z[m] = 2 * (np.cos(ph) @ (n**-0.5))
        p = a - N
        cp = np.cos(TWO_PI * p)
        cp = np.where(np.abs(cp) < 1e-8, 1e-8 * np.sign(cp + 1e-300), cp)
        psi = np.cos(TWO_PI * (p**2 - p - 1 / 16)) / cp
        z += (-1) ** (N - 1) * (tt / TWO_PI) ** -0.25 * psi
        out[s:s + chunk] = z
    return out

def prime_cols(tvals):
    cols = []
    for p in PRIMES:
        cols += [np.cos(tvals * np.log(p)), np.sin(tvals * np.log(p))]
    return cols

d41 = np.load(HERE / "41_bigT_windows.npz")
keys = sorted({k.split("_")[1] for k in d41.files}, key=lambda s: int(s[:-1]))

rows = []
for wkey in keys:
    gaps = d41[f"gaps_{wkey}"]; tmid = d41[f"tmid_{wkey}"]
    Lw = np.log(tmid / TWO_PI); L = float(Lw.mean())
    g_lo = tmid - gaps / 2

    # tepe konumu (45 ile aynı grid)
    u = np.arange(1, 25) / 25
    tt = g_lo[:, None] + gaps[:, None] * u[None, :]
    Zg = np.abs(Z_rs(tt.ravel())).reshape(tt.shape)
    j = np.argmax(Zg, axis=1)
    t_pk = tt[np.arange(len(j)), j]
    amps = Zg[np.arange(len(j)), j]

    g_u = gaps * Lw / TWO_PI
    a_u = amps / np.sqrt(A * Lw + B0 + B1 / Lw)
    a_u /= np.sqrt((a_u**2).mean())
    y_a, y_g = np.log(a_u), np.log(g_u)

    # HAM r'ler
    rP_raw, _ = pearsonr(g_u, a_u)
    rS_raw, _ = spearmanr(g_u, a_u)

    # GENLİK soyması — TUTARLILIK KRİTİK: g̃ kontrolü KOYMA!
    # (Kontrollü fit yalnız w'yi ölçer; boşluk-aracılı β·S_v genlikte kalır
    #  ve g̃ soyulunca çift uyumsuzlaşır — ilk denemede Spearman bu yüzden
    #  çöktü. Kontrolsüz fit w + β·v toplam asal payını yakalar.)
    pc_a = prime_cols(t_pk)
    Xa = np.vstack([np.ones_like(y_a)] + pc_a).T
    ba = np.linalg.lstsq(Xa, y_a, rcond=None)[0]
    prime_part_a = np.vstack(pc_a).T @ ba[1:]
    a_res = np.exp(y_a - prime_part_a)
    a_res /= np.sqrt((a_res**2).mean())

    # BOŞLUK soyması
    pc_g = prime_cols(tmid)
    Xg = np.vstack([np.ones_like(y_g)] + pc_g).T
    bg = np.linalg.lstsq(Xg, y_g, rcond=None)[0]
    prime_part_g = np.vstack(pc_g).T @ bg[1:]
    g_res = np.exp(y_g - prime_part_g)
    g_res *= g_u.mean() / g_res.mean()

    rP_res, _ = pearsonr(g_res, a_res)
    rS_res, _ = spearmanr(g_res, a_res)

    NP_raw = inv_fit(FIT_P, rP_raw); NS_raw = inv_fit(FIT_S, rS_raw)
    NP_res = inv_fit(FIT_P, rP_res); NS_res = inv_fit(FIT_S, rS_res)
    rows.append((L, NP_raw - L, NS_raw - L, NP_res - L, NS_res - L,
                 rP_raw, rP_res, rS_raw, rS_res))
    print(f"L={L:5.2f} | HAM: N_eff−L  P {NP_raw-L:+.2f}  S {NS_raw-L:+.2f}"
          f" | SOYULMUŞ: P {NP_res-L:+.2f}  S {NS_res-L:+.2f}"
          f" | r_P {rP_raw:.4f}→{rP_res:.4f}  r_S {rS_raw:.4f}→{rS_res:.4f}")

arr = np.array(rows)
print("\nÖZET (6 pencere ortalaması):")
print(f"  HAM      : N_eff−L  Pearson {arr[:,1].mean():+.3f}  Spearman {arr[:,2].mean():+.3f}"
      f"  → AYRIK (fark {arr[:,2].mean()-arr[:,1].mean():+.3f})")
print(f"  SOYULMUŞ : N_eff−L  Pearson {arr[:,3].mean():+.3f}  Spearman {arr[:,4].mean():+.3f}"
      f"  → fark {arr[:,4].mean()-arr[:,3].mean():+.3f}")

fig, ax = plt.subplots(figsize=(8.5, 5.4))
ax.plot(arr[:, 0], arr[:, 1], "o-", c="firebrick", label="ham Pearson")
ax.plot(arr[:, 0], arr[:, 2], "s-", c="darkorange", label="ham Spearman")
ax.plot(arr[:, 0], arr[:, 3], "o--", c="steelblue", label="soyulmuş Pearson")
ax.plot(arr[:, 0], arr[:, 4], "s--", c="teal", label="soyulmuş Spearman")
ax.axhline(0, color="k", lw=0.8)
ax.set_xlabel("L"); ax.set_ylabel("N_eff − L")
ax.set_title("Asal soyma: iki metrik buluşuyor mu?")
ax.legend(fontsize=9); ax.grid(alpha=0.3)
plt.tight_layout()
out = HERE / "49_asal_soyma.png"
plt.savefig(out, dpi=110)
print(f"\nGrafik: {out.name}")
