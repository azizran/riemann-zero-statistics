"""
58 — v-KANALI L≈45-47'DE: YASANIN EN UZAK TESTİ (17 Ağustos 2026, gece)
=========================================================================

zeros4 (10²¹'inci sıfır civarı, L=44.58) ve zeros5 (10²²'inci, L=46.83):
genlik hesabı bu yüksekliklerde imkânsız ama v-KANALI yalnız sıfır
konumlarından ölçülür — tablolar yeter.

Test: v(τ) yasası (küçük-τ'da v ≈ a·τ, önceki verilerden a fit edilir)
τ₂ = 0.0155'e kadar tutuyor mu? Bu, t'de 5×10⁸ kat daha uzak bir
out-of-sample sıçrama (toplam kaldıraç: 18 büyüklük mertebesi).

Fazlar tamsayı tabanlardan çapalı (mpmath, kesin mod 2π).
Kontroller: sahte frekanslar + sin/cos faz yapısı.
"""

import numpy as np
import mpmath as mp
import matplotlib.pyplot as plt
from pathlib import Path

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
PRIMES = [2, 3, 5, 7, 11, 13]
FAKES = [2.31, 6.7]

FILES = [
    ("odlyzko_zeros4.txt", 144176897509546973000, "10^21"),
    ("odlyzko_zeros5.txt", 1370919909931995300000, "10^22"),
]

# önceki v ölçümleri (47 pencereleri + 54 Odlyzko) — çökme grafiği için
PREV = []
V47 = {9.86: [0.1405, 0.2245, 0.3254, 0.3877, 0.4663, 0.4928],
       10.37: [0.1358, 0.2159, 0.3108, 0.3713, 0.4473, 0.4744],
       10.93: [0.1279, 0.2034, 0.2962, 0.3576, 0.4299, 0.4581],
       11.47: [0.1222, 0.1950, 0.2823, 0.3403, 0.4172, 0.4404],
       11.98: [0.1174, 0.1847, 0.2701, 0.3258, 0.3988, 0.4220],
       12.45: [0.1114, 0.1787, 0.2628, 0.3123, 0.3807, 0.4082]}
for L, vs in V47.items():
    for p, v in zip(PRIMES, vs):
        PREV.append((np.log(p) / L, v))
for p, v in zip(PRIMES, [0.054, 0.091, 0.126, 0.163, 0.202, 0.202]):
    PREV.append((np.log(p) / 24.475, v))
PREV = np.array(PREV)

# küçük-τ kılavuzu: önceki verinin τ<0.11 kısmına v = a·τ fiti
m = PREV[:, 0] < 0.11
a_guide = float(np.sum(PREV[m, 0] * PREV[m, 1]) / np.sum(PREV[m, 0] ** 2))
print(f"Küçük-τ kılavuzu (önceki veriden): v ≈ {a_guide:.3f}·τ\n")

mp.mp.dps = 50
results = []
for fname, base, label in FILES:
    offs = []
    with open(HERE / fname) as f:
        for line in f:
            try:
                offs.append(float(line.strip()))
            except ValueError:
                continue
    offs = np.array(offs)
    Lh = float(mp.log(mp.mpf(base) / (2 * mp.pi)))
    gaps = np.diff(offs)
    tmid = 0.5 * (offs[:-1] + offs[1:])
    g_u = gaps * Lh / TWO_PI
    y = np.log(g_u)
    n = len(y)
    print(f"=== {label} (L={Lh:.3f}, {n} aralık, mean g̃={g_u.mean():.4f}) ===")

    cols = [np.ones_like(y)]
    for q in PRIMES + FAKES:
        ph0 = float(mp.fmod(mp.mpf(base) * mp.log(q), 2 * mp.pi))
        arg = ph0 + tmid * np.log(q)
        cols += [np.cos(arg), np.sin(arg)]
    X = np.vstack(cols).T
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ b
    se = np.sqrt(res.var() * np.diag(np.linalg.inv(X.T @ X)))

    for i, p in enumerate(PRIMES):
        c_p, d_p = b[1 + 2 * i], b[2 + 2 * i]
        v = np.hypot(c_p, d_p) / p**-0.5
        sev = se[1 + 2 * i] / p**-0.5
        tau = np.log(p) / Lh
        pred = a_guide * tau
        sig = np.hypot(c_p, d_p) / se[1 + 2 * i]
        results.append((tau, v, sev, Lh, p))
        print(f"  p={p:>2} (τ={tau:.4f}): |v| = {v:.4f} ± {sev:.4f}"
              f"   kılavuz {pred:.4f}   (sinyal {sig:.1f}σ)")
    foff = 1 + 2 * len(PRIMES)
    fr = np.hypot(b[foff], b[foff + 1]), np.hypot(b[foff + 2], b[foff + 3])
    print(f"  sahte frekans genlikleri: {fr[0]:.4f}, {fr[1]:.4f} (≈0 olmalı)\n")

# grafik: tam çökme, 18 büyüklük mertebesi
res = np.array([(t, v, s) for t, v, s, _, _ in results])
fig, ax = plt.subplots(figsize=(8.5, 5.6))
ax.plot(PREV[:, 0], PREV[:, 1], "o", ms=4, c="steelblue", alpha=0.6,
        label="t ≤ 2.7×10¹¹ (42 nokta)")
ax.errorbar(res[:, 0], res[:, 1], yerr=res[:, 2], fmt="s", ms=7,
            c="firebrick", zorder=5, label="t = 1.4×10²⁰ ve 1.4×10²¹ (YENİ)")
xx = np.linspace(0, 0.28, 50)
ax.plot(xx, a_guide * xx, "k--", lw=1.1, label=f"kılavuz v = {a_guide:.2f}τ (eski veriye fit)")
ax.set_xlabel("τ = log p / L")
ax.set_ylabel("|v| — boşluk kanalı emilimi")
ax.set_title("v-yasası: 18 büyüklük mertebesinde")
ax.legend(fontsize=9)
ax.grid(alpha=0.3)
# küçük-τ zoom
axins = ax.inset_axes([0.55, 0.08, 0.42, 0.38])
axins.plot(PREV[:, 0], PREV[:, 1], "o", ms=3, c="steelblue", alpha=0.6)
axins.errorbar(res[:, 0], res[:, 1], yerr=res[:, 2], fmt="s", ms=5, c="firebrick")
axins.plot(xx, a_guide * xx, "k--", lw=1)
axins.set_xlim(0, 0.07); axins.set_ylim(0, 0.16)
axins.grid(alpha=0.3)
plt.tight_layout()
out = HERE / "58_v_kanali_L45.png"
plt.savefig(out, dpi=110)
print(f"Grafik: {out.name}")
