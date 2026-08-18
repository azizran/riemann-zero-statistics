"""
62 — v-KANALI 36 PENCERELERİNDE: GERÇEK "EIGHTEEN ORDERS" (18 Ağustos 2026)
=============================================================================

Denetim K3: v-yasası "18 mertebe" iddiası için alçak pencerelerde (36 verisi,
t=1.2×10³'ten itibaren) v hiç ölçülmemişti. Ölçüyoruz — Z hesabı gerekmez,
36_T100k.npz'deki boşluklar yeter.

Ek değer: 36 (L=5.6-9.1) ile 41 (L=9.9-12.5) pencereleri τ'da örtüşür →
motorlar-arası v tutarlılık kontrolü. Ve büyük-τ ucu (p=13@L=5.6: τ=0.46,
τ*'ın ötesi!) v-eğrisinin yeni bölgesi.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
PRIMES = [2, 3, 5, 7, 11, 13]
FAKES = [2.31, 6.7]

d36 = np.load(HERE / "36_T100k.npz")
t_all = d36["t_mid"]; g_all = d36["intervals"]
edges = np.geomspace(t_all[0], t_all[-1] * 1.0001, 13)

results = []
print(f"{'L':>6} {'t_ort':>9} {'n':>6} | " +
      " ".join(f"|v|({p})".rjust(8) for p in PRIMES) + " | sahte")
for i in range(12):
    m = (t_all >= edges[i]) & (t_all < edges[i + 1])
    if m.sum() < 500:
        continue
    tmid = t_all[m]; gaps = g_all[m]
    Lw = np.log(tmid / TWO_PI); L = float(Lw.mean())
    y = np.log(gaps * Lw / TWO_PI)
    cols = [np.ones_like(y)]
    for q in PRIMES + FAKES:
        arg = tmid * np.log(q)
        cols += [np.cos(arg), np.sin(arg)]
    X = np.vstack(cols).T
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ b
    se = np.sqrt(res.var() * np.diag(np.linalg.inv(X.T @ X)))
    vs, ss = [], []
    for j, p in enumerate(PRIMES):
        v = np.hypot(b[1 + 2 * j], b[2 + 2 * j]) / p**-0.5
        sv = se[1 + 2 * j] / p**-0.5
        vs.append(v); ss.append(sv)
        results.append((np.log(p) / L, v, sv, L, p))
    foff = 1 + 2 * len(PRIMES)
    fr = max(np.hypot(b[foff], b[foff + 1]), np.hypot(b[foff + 2], b[foff + 3]))
    print(f"{L:>6.2f} {tmid.mean():>9.0f} {m.sum():>6} | " +
          " ".join(f"{v:8.4f}" for v in vs) + f" | {fr:.4f}")

# 41/47 verileriyle örtüşme kontrolü ve tam çökme grafiği
V47 = {9.86: [0.1405, 0.2245, 0.3254, 0.3877, 0.4663, 0.4928],
       12.45: [0.1114, 0.1787, 0.2628, 0.3123, 0.3807, 0.4082]}
res = np.array([(t, v, s) for t, v, s, _, _ in results])

fig, ax = plt.subplots(figsize=(8.5, 5.4))
ax.errorbar(res[:, 0], res[:, 1], yerr=res[:, 2], fmt="o", ms=4,
            c="darkorange", alpha=0.8, label="36 pencereleri (YENİ, t ≥ 1.2×10³)")
for L, vs in V47.items():
    xs = [np.log(p) / L for p in PRIMES]
    ax.plot(xs, vs, "s", ms=4, c="steelblue", alpha=0.6,
            label="41 pencereleri (47)" if L == 9.86 else None)
xx = np.linspace(0, 0.30, 50)
ax.plot(xx, 2.014 * xx, "k--", lw=1, label="küçük-τ kılavuzu 2.01τ")
ax.axvline(0.40, color="gray", ls=":", lw=1, label="τ* ≈ 0.40")
ax.set_xlabel("τ = log p / L"); ax.set_ylabel("|v|")
ax.set_title("v-eğrisi: alçak pencereler eklendi (t=1.2×10³'ten)")
ax.legend(fontsize=9); ax.grid(alpha=0.3)
plt.tight_layout()
out = HERE / "62_v_36_pencereleri.png"
plt.savefig(out, dpi=110)
print(f"\nGrafik: {out.name}")
print(f"\nEn alçak pencere t ≈ {t_all[ (t_all >= edges[6]) ][0]:.0f} → "
      f"L=44.6-46.8 tablolarıyla menzil: "
      f"log10(1.37e21 / {t_all[(t_all >= edges[6])][0]:.0f}) = "
      f"{np.log10(1.37e21 / t_all[(t_all >= edges[6])][0]):.2f} mertebe")
