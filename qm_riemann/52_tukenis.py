"""
52 — TÜKENİŞ EĞRİSİ: KALAN %32 NE? (17 Ağustos 2026, gece)
============================================================

51'de q≤49 soyunca r 0.974'e çıktı. Şimdi soyma kesim noktası Q'yu
büyütüyoruz: tüm asal kuvvetleri q = p^k ≤ Q, Q = 13 → 300.

  r*(Q) doyarsa  → doyma değeri = çekirdeğin gerçek sıkılığı;
                   kalan varyans V_res = sıfır gazının ÖZ-rastgeleliği
                   (yeni ölçülebilir nicelik: V_res(L) eğrisi)
  doymazsa       → tek-aralık yasası tamamen asal-deterministik'e gidiyor

Her Q'da plasebo AYNI SAYIDA sahte frekansla koşuyor → overfit tabanı
her noktada görünür. (200 sütun / 40k satır'da bile taban ~%0.5 olmalı.)
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import pearsonr, spearmanr
from sympy import primerange

rng = np.random.default_rng(52)
HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
A = (np.e**2 - 5) / 2
B0, B1 = 2.7580, -0.0543

def prime_powers_upto(Q):
    """log-frekans listesi: tüm q = p^k ≤ Q."""
    out = []
    for p in primerange(2, Q + 1):
        pk = p
        while pk <= Q:
            out.append(np.log(float(pk)))
            pk *= p
    return np.array(sorted(out))

def placebo_freqs(n, lo, hi, avoid, seed_rng):
    out = []
    while len(out) < n:
        c = seed_rng.uniform(lo, hi)
        if np.abs(avoid - c).min() > 0.015 and all(abs(o - c) > 0.015 for o in out):
            out.append(c)
    return np.array(out)

def strip_freqs(y, tvals, freqs):
    cols = [np.ones_like(y)]
    for om in freqs:
        cols += [np.cos(tvals * om), np.sin(tvals * om)]
    C = np.vstack(cols).T
    b = np.linalg.lstsq(C, y, rcond=None)[0]
    part = C[:, 1:] @ b[1:]
    return y - part, float(part.var())

QS = [13, 25, 50, 100, 200, 300]
d41 = np.load(HERE / "41_bigT_windows.npz")
keys = sorted({k.split("_")[1] for k in d41.files}, key=lambda s: int(s[:-1]))

curves = {}   # wkey -> (L, [r_P(Q)], [r_S(Q)], [V_res(Q)], [r_P placebo], V_toplam)
for wkey in keys:
    gaps = d41[f"gaps_{wkey}"]; tmid = d41[f"tmid_{wkey}"]; amps = d41[f"amps_{wkey}"]
    Lw = np.log(tmid / TWO_PI); L = float(Lw.mean())
    g_u = gaps * Lw / TWO_PI
    a_u = amps / np.sqrt(A * Lw + B0 + B1 / Lw)
    a_u /= np.sqrt((a_u**2).mean())
    y_a, y_g = np.log(a_u), np.log(g_u)
    V_tot = float(y_a.var())

    rP_list, rS_list, vres_list, plc_list = [], [], [], []
    for Q in QS:
        freqs = prime_powers_upto(Q)
        ya, va = strip_freqs(y_a, tmid, freqs)
        yg, _ = strip_freqs(y_g, tmid, freqs)
        ar = np.exp(ya); ar /= np.sqrt((ar**2).mean())
        gr = np.exp(yg); gr *= g_u.mean() / gr.mean()
        rP, _ = pearsonr(gr, ar); rS, _ = spearmanr(gr, ar)
        rP_list.append(rP); rS_list.append(rS)
        vres_list.append(V_tot - va)

        # K2 düzeltmesi: çekim aralığındaki TÜM asal kuvvetlerinden kaçın
        avoid_all = np.log(np.array(prime_powers_upto(360)))
        pf = placebo_freqs(len(freqs), np.log(2) * 0.9, np.log(300) * 1.02,
                           avoid_all, rng)
        yap, _ = strip_freqs(y_a, tmid, pf)
        ygp, _ = strip_freqs(y_g, tmid, pf)
        arp = np.exp(yap); arp /= np.sqrt((arp**2).mean())
        grp = np.exp(ygp); grp *= g_u.mean() / grp.mean()
        rPp, _ = pearsonr(grp, arp)
        plc_list.append(rPp)
    curves[wkey] = (L, rP_list, rS_list, vres_list, plc_list, V_tot)
    n_last = len(prime_powers_upto(QS[-1]))
    print(f"L={L:5.2f}: r_P " + " ".join(f"{r:.4f}" for r in rP_list)
          + f" | plasebo(son) {plc_list[-1]:.4f} | V_res(son) {vres_list[-1]:.4f}"
          f"/{V_tot:.3f} | son frekans sayısı {n_last}")

# özet
print(f"\nQ değerleri: {QS}")
print("\n=== DOYMA ANALİZİ ===")
for wkey, (L, rP, rS, vres, plc, V_tot) in curves.items():
    d_last = rP[-1] - rP[-2]
    print(f"  L={L:5.2f}: r_P son iki adım {rP[-2]:.4f}→{rP[-1]:.4f} (Δ={d_last:+.4f})"
          f" | V_res/V_tot = {vres[-1]/V_tot:.3f}"
          f" | plasebo kayması {plc[-1] - rP[0] + (rP[0]-plc[0]):+.4f}")

print("\nV_res(L) — sıfır gazının öz-rastgeleliği (aday yeni nicelik):")
for wkey, (L, rP, rS, vres, plc, V_tot) in curves.items():
    print(f"  L={L:5.2f}: V_res = {vres[-1]:.4f}  (toplamın %{100*vres[-1]/V_tot:.1f}'i)"
          f"  r*_P = {rP[-1]:.4f}  r*_S = {rS[-1]:.4f}")

# grafik
fig, axes = plt.subplots(1, 2, figsize=(13, 5.2))
ax = axes[0]
for wkey, (L, rP, rS, vres, plc, V_tot) in curves.items():
    ax.plot(QS, rP, "o-", ms=4, label=f"L={L:.1f}")
    ax.plot(QS, plc, ":", c="gray", lw=0.8)
ax.plot([], [], ":", c="gray", label="plasebolar")
ax.set_xscale("log"); ax.set_xlabel("soyma kesimi Q")
ax.set_ylabel("soyulmuş r (Pearson)")
ax.set_title("Tükeniş eğrisi: r*(Q)")
ax.legend(fontsize=8, ncol=2); ax.grid(alpha=0.3)

ax = axes[1]
Ls = [c[0] for c in curves.values()]
Vr = [c[3][-1] for c in curves.values()]
ax.plot(Ls, Vr, "o-", c="firebrick")
ax.set_xlabel("L"); ax.set_ylabel("V_res (Q=300 sonrası)")
ax.set_title("Öz-rastgelelik V_res(L)")
ax.grid(alpha=0.3)
plt.tight_layout()
out = HERE / "52_tukenis.png"
plt.savefig(out, dpi=110)
print(f"\nGrafik: {out.name}")
