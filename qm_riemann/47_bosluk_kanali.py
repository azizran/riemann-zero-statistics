"""
47 — BOŞLUK KANALI: ASAL DALGASI BOŞLUKLARI DA SÜRÜYOR MU? (17 Ağustos 2026)
=============================================================================

46'nın negatif sonucunun teşhisi: bağımsız-fazlı asal gürültüsü r'yi gerçekte
olduğundan 7 kat fazla bozuyor. Eksik kanal adayı: asal dalgası BOŞLUKLARI da
aynı fazla modüle ediyorsa (explicit formula: sıfır konumları asallara cevap
verir), ortak faz korelasyonu yukarı geri iter.

Test (anlık, Z hesabı yok): log g̃_n ~ Σ_p [c_p cos(t_mid log p) + d_p sin(...)]
→ v_p = katsayı / p^{-1/2}. Sahte frekans kontrolleri dahil.

Beklenti: v_p ≠ 0 ve işaretlerin w_p ile ilişkisi ortak-faz mekanizmasını
belirler. (Dikkat: sin bileşeni burada ≠0 olabilir — boşluk yanıtı türev
kanalıdır, faz kayması DOĞALDIR. Genlik |v| = √(c²+d²) raporlanır.)
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
PRIMES = [2, 3, 5, 7, 11, 13]
FAKES = [2.5, 6.0]

d41 = np.load(HERE / "41_bigT_windows.npz")
keys = sorted({k.split("_")[1] for k in d41.files}, key=lambda s: int(s[:-1]))

print(f"{'L':>6} | " + " ".join(f"|v|({p})".rjust(8) for p in PRIMES)
      + " |  faz(2)  faz(3) | sahte RMS")
res = {p: [] for p in PRIMES}
for wkey in keys:
    gaps = d41[f"gaps_{wkey}"]; tmid = d41[f"tmid_{wkey}"]
    Lw = np.log(tmid / TWO_PI); L = float(Lw.mean())
    g_u = gaps * Lw / TWO_PI
    y = np.log(g_u)

    cols = [np.ones_like(y)]
    for p in PRIMES:
        cols += [np.cos(tmid * np.log(p)), np.sin(tmid * np.log(p))]
    for f in FAKES:
        cols += [np.cos(tmid * np.log(f)), np.sin(tmid * np.log(f))]
    X = np.vstack(cols).T
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    cov = resid.var() * np.linalg.inv(X.T @ X)
    se = np.sqrt(np.diag(cov))

    vabs, phases = [], []
    for i, p in enumerate(PRIMES):
        c_p, d_p = beta[1 + 2 * i], beta[2 + 2 * i]
        v = np.hypot(c_p, d_p) / p**-0.5
        se_v = se[1 + 2 * i] / p**-0.5
        res[p].append((v, se_v))
        vabs.append(v)
        phases.append(np.arctan2(d_p, c_p))
    foff = 1 + 2 * len(PRIMES)
    fake_rms = np.sqrt(np.mean(beta[foff:foff + 4]**2))
    print(f"{L:>6.2f} | " + " ".join(f"{v:8.4f}" for v in vabs)
          + f" | {phases[0]:+7.2f} {phases[1]:+7.2f} | {fake_rms:9.5f}")

print("\nAğırlıklı ortalama |v(p)| ve gürültü tabanıyla kıyas:")
for p in PRIMES:
    vs = np.array([r[0] for r in res[p]])
    es = np.array([r[1] for r in res[p]])
    vm = np.sum(vs / es**2) / np.sum(1 / es**2)
    em = 1 / np.sqrt(np.sum(1 / es**2))
    # |v| pozitif yanlı (Rayleigh): gürültü tabanı ~ se·√(π/2)
    floor = em * np.sqrt(np.pi / 2) * np.sqrt(len(vs))
    sig = "SİNYAL" if vm > 3 * floor else "gürültü seviyesinde"
    print(f"  |v({p:>2})| = {vm:.4f} ± {em:.4f}   ({sig})")

# w ile kıyas grafiği (45 tablosundan w ortalamaları)
W_MEAN = {2: 0.764, 3: 0.647, 5: 0.501, 7: 0.410, 11: 0.295, 13: 0.255}
fig, ax = plt.subplots(figsize=(7.5, 5))
pm = np.array(PRIMES, dtype=float)
vm_arr = [np.sum(np.array([r[0] for r in res[p]]) / np.array([r[1] for r in res[p]])**2)
          / np.sum(1 / np.array([r[1] for r in res[p]])**2) for p in PRIMES]
ax.plot(pm, [W_MEAN[p] for p in PRIMES], "s-", c="firebrick", label="w(p) — genlik kanalı (45)")
ax.plot(pm, vm_arr, "o-", c="steelblue", label="|v(p)| — boşluk kanalı (47)")
ax.set_xlabel("p"); ax.set_ylabel("iletim katsayısı")
ax.set_title("İki kanal: asal dalgası genliği vs boşlukları nasıl sürüyor")
ax.legend(); ax.grid(alpha=0.3)
plt.tight_layout()
out = HERE / "47_bosluk_kanali.png"
plt.savefig(out, dpi=110)
print(f"\nGrafik: {out.name}")
