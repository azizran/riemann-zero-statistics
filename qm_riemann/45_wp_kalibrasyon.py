"""
45 — w(p) KALİBRASYONU: ζ ASAL DALGASININ YÜZDE KAÇINI TAŞIYOR? (17 Ağustos 2026)
==================================================================================

Explicit formula: log|Z(t)| ⊃ Σ_p w_p · p^{-1/2} cos(t log p) + ...
(tam Euler çarpanı = w_p = 1; 44'teki aşırı-çıkarma w<1'i VEYA istatistik
artefaktını gösteriyor olabilirdi — bu ölçüm ikisini AYIRIR.)

Yöntem (CUE'suz, bölmesiz, doğrudan):
  Her aralıkta tepe konumu t_pk ve unfold tepe ã ölçülü.
  Regresyon: log ã ~ 1 + g̃ + g̃² + Σ_p [a_p cos(t_pk log p) + b_p sin(t_pk log p)]
  → w_p = a_p / p^{-1/2}

Yerleşik kontroller:
  - b_p ≈ 0 olmalı (teori saf kosinüs; faz kayması = sorun işareti)
  - Sahte frekanslar (ω = log 2.5, log 6.0) → katsayı ≈ 0 olmalı
  - w_p(L) pencereler arası tutarlılık

p büyüdükçe uyarı: dalga periyodu 2π/log p aralık boyutuna yaklaşır (p=13'te
~5 aralık) → "aralık içinde sabit" varsayımı zayıflar, hafif yanlılık olabilir.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
A = (np.e**2 - 5) / 2
B0, B1 = 2.7580, -0.0543

# --- Riemann-Siegel (41 ile aynı) ---
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

PRIMES = [2, 3, 5, 7, 11, 13]
FAKES = [2.5, 6.0]  # kontrol frekansları (asal değil, log'ları asal log'u değil)

d41 = np.load(HERE / "41_bigT_windows.npz")
keys = sorted({k.split("_")[1] for k in d41.files}, key=lambda s: int(s[:-1]))

results = {p: [] for p in PRIMES}
fakes_res = {f: [] for f in FAKES}
sins_res = {p: [] for p in PRIMES}
Ls = []

print(f"{'L':>6} | " + " ".join(f"w({p})".rjust(7) for p in PRIMES)
      + " | " + " ".join(f"sahte{f}".rjust(8) for f in FAKES))
for wkey in keys:
    gaps = d41[f"gaps_{wkey}"]; tmid = d41[f"tmid_{wkey}"]
    Lw = np.log(tmid / TWO_PI)
    L = float(Lw.mean())
    g_lo = tmid - gaps / 2

    # tepe konumu + genlik (24 nokta grid — 44 ile aynı)
    u = np.arange(1, 25) / 25
    tt = g_lo[:, None] + gaps[:, None] * u[None, :]
    Zg = np.abs(Z_rs(tt.ravel())).reshape(tt.shape)
    j = np.argmax(Zg, axis=1)
    t_pk = tt[np.arange(len(j)), j]
    amps = Zg[np.arange(len(j)), j]

    g_u = gaps * Lw / TWO_PI
    a_u = amps / np.sqrt(A * Lw + B0 + B1 / Lw)
    y = np.log(a_u)

    # tasarım matrisi
    cols = [np.ones_like(y), g_u, g_u**2]
    for p in PRIMES:
        cols += [np.cos(t_pk * np.log(p)), np.sin(t_pk * np.log(p))]
    for f in FAKES:
        cols += [np.cos(t_pk * np.log(f))]
    X = np.vstack(cols).T
    beta, res_, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    sigma2 = resid.var()
    cov = sigma2 * np.linalg.inv(X.T @ X)
    se = np.sqrt(np.diag(cov))

    Ls.append(L)
    off = 3
    wvals = []
    for i, p in enumerate(PRIMES):
        a_p, b_p = beta[off + 2 * i], beta[off + 2 * i + 1]
        w_p = a_p / p**-0.5
        se_w = se[off + 2 * i] / p**-0.5
        results[p].append((w_p, se_w))
        sins_res[p].append(b_p / p**-0.5)
        wvals.append(w_p)
    foff = off + 2 * len(PRIMES)
    fvals = []
    for i, f in enumerate(FAKES):
        fv = beta[foff + i]
        fakes_res[f].append((fv, se[foff + i]))
        fvals.append(fv)
    print(f"{L:>6.2f} | " + " ".join(f"{w:7.3f}" for w in wvals)
          + " | " + " ".join(f"{v:8.4f}" for v in fvals))

# özetler
print("\nAğırlıklı ortalama w(p) (6 pencere):")
for p in PRIMES:
    ws = np.array([r[0] for r in results[p]])
    es = np.array([r[1] for r in results[p]])
    wm = np.sum(ws / es**2) / np.sum(1 / es**2)
    em = 1 / np.sqrt(np.sum(1 / es**2))
    sin_rms = np.sqrt(np.mean(np.array(sins_res[p])**2))
    print(f"  w({p:>2}) = {wm:+.4f} ± {em:.4f}   (sin bileşeni RMS {sin_rms:.3f} — 0 olmalı)")
print("Sahte frekans kontrolleri (0 olmalı):")
for f in FAKES:
    vs = np.array([r[0] for r in fakes_res[f]])
    print(f"  ω=log{f}: ortalama katsayı {vs.mean():+.4f} (RMS {np.sqrt((vs**2).mean()):.4f})")

# grafik
fig, axes = plt.subplots(1, 2, figsize=(12.5, 5))
ax = axes[0]
for p in PRIMES:
    ws = [r[0] for r in results[p]]
    es = [r[1] for r in results[p]]
    ax.errorbar(Ls, ws, yerr=es, fmt="o-", ms=4, label=f"p={p}")
ax.axhline(1, color="k", ls="--", lw=1, label="tam içerik (w=1)")
ax.axhline(0, color="gray", lw=0.7)
ax.set_xlabel("L = log(t/2π)"); ax.set_ylabel("w(p)")
ax.set_title("ζ'nın taşıdığı asal-dalga oranı")
ax.legend(fontsize=8, ncol=2); ax.grid(alpha=0.3)

ax = axes[1]
pm = np.array(PRIMES, dtype=float)
wm_arr, em_arr = [], []
for p in PRIMES:
    ws = np.array([r[0] for r in results[p]])
    es = np.array([r[1] for r in results[p]])
    wm_arr.append(np.sum(ws / es**2) / np.sum(1 / es**2))
    em_arr.append(1 / np.sqrt(np.sum(1 / es**2)))
ax.errorbar(pm, wm_arr, yerr=em_arr, fmt="s", ms=6, c="firebrick",
            label="ölçülen w(p)")
ax.axhline(1, color="k", ls="--", lw=1)
for f in FAKES:
    vs = np.array([r[0] for r in fakes_res[f]])
    ax.plot(f, vs.mean(), "x", ms=9, c="gray")
ax.plot([], [], "x", c="gray", label="sahte frekans kontrolleri")
ax.set_xlabel("p"); ax.set_ylabel("w(p) (pencere ortalaması)")
ax.set_title("Asal içerik spektrumu")
ax.legend(fontsize=9); ax.grid(alpha=0.3)

plt.tight_layout()
out = HERE / "45_wp_kalibrasyon.png"
plt.savefig(out, dpi=110)
print(f"\nGrafik: {out.name}")
