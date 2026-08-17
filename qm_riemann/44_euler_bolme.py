"""
44 — EULER BÖLME TESTİ: AÇIĞI ASALLAR MI YARATIYOR? (17 Ağustos 2026)
======================================================================

43'ün bulgusu: ζ vs CUE açığı tamamen koşullu varyansta, büyük boşluklarda.
Hipotez: bunu küçük asalların uzun-periyotlu modülasyonu yaratıyor
(p=2: periyot 9.1, p=3: 5.7, p=5: 3.9 — hepsi ortalama aralıktan kat kat uzun).

Test (GHK hibrit çarpım çerçevesi, Gonek-Hughes-Keating 2007):
    ζ ≈ P_X(asal kısmı) × Z_X(sıfır kısmı ≈ CUE)
|Z_removed(t)| = |Z(t)| · Π_{p∈P} |1 − p^{-1/2-it}|   (p-Euler çarpanı bölünür;
sıfırlar AYNI kalır, sadece genlik zarfı değişir).

Öngörü: P büyüdükçe
  (a) büyük-boşluk varyans fazlası ÇÖKMELİ,
  (b) r, CUE(N=L) değerine YÜKSELMELİ,
  (c) N_eff − L kayması 0'a inmeli.
Çökmezse → mekanizma asal dalgaları değil, daha derin bir şey.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import pearsonr

rng = np.random.default_rng(44)
HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
A = (np.e**2 - 5) / 2

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

def euler_factor(t, primes):
    """Π_p |1 − p^{-1/2} e^{-it log p}|  (bölme = bu çarpanla ÇARPMA)."""
    E = np.ones_like(t)
    for p in primes:
        E *= np.abs(1 - p**-0.5 * np.exp(-1j * t * np.log(p)))
    return E

# --- 43 ile aynı binleme ---
BINS = np.concatenate([[0.05], np.linspace(0.25, 2.2, 14), [3.5]])
BC = 0.5 * (BINS[:-1] + BINS[1:])

def bin_stats(g, a, min_n=200):
    mean = np.full(len(BC), np.nan)
    var = np.full(len(BC), np.nan)
    idx = np.digitize(g, BINS) - 1
    for i in range(len(BC)):
        m = idx == i
        if m.sum() >= min_n:
            mean[i] = a[m].mean(); var[i] = a[m].var()
    return mean, var

# --- CUE referansı (43 yöntemi) ---
def haar_unitary_batch(M, N):
    G = (rng.standard_normal((M, N, N)) + 1j * rng.standard_normal((M, N, N))) / np.sqrt(2)
    Q, R = np.linalg.qr(G)
    diag = np.einsum("mii->mi", R)
    return Q * (diag / np.abs(diag))[:, None, :]

def cue_ref(N, n_gaps=300000, grid=24, chunk=1200):
    M = int(np.ceil(n_gaps / N))
    gs, ms = [], []
    for s in range(0, M, chunk):
        m = min(chunk, M - s)
        U = haar_unitary_batch(m, N)
        ph = np.sort(np.angle(np.linalg.eigvals(U)), axis=1)
        gaps = np.diff(np.concatenate([ph, ph[:, :1] + TWO_PI], axis=1), axis=1)
        u = np.arange(1, grid + 1) / (grid + 1)
        th = ph[:, :, None] + gaps[:, :, None] * u[None, None, :]
        df = th[:, :, :, None] - ph[:, None, None, :]
        amp = np.prod(2 * np.abs(np.sin(df / 2)), axis=-1)
        gs.append(gaps.ravel()); ms.append(amp.max(axis=-1).ravel())
    g = np.concatenate(gs) * N / TWO_PI
    a = np.concatenate(ms); a /= np.sqrt((a**2).mean())
    r, _ = pearsonr(g, a)
    return bin_stats(g, a), r

# ---------------------------------------------------------------
# ANA DÖNGÜ: pencere × asal kümesi
# ---------------------------------------------------------------
d41 = np.load(HERE / "41_bigT_windows.npz")
PSETS = [((), "ham"), ((2,), "p=2"), ((2, 3), "p≤3"),
         ((2, 3, 5, 7), "p≤7"), ((2, 3, 5, 7, 11, 13), "p≤13")]
WINDOWS = ["120k", "1600k"]

results = {}
for wkey in WINDOWS:
    gaps = d41[f"gaps_{wkey}"]; tmid = d41[f"tmid_{wkey}"]
    L = float(np.log(tmid / TWO_PI).mean())
    g_lo = tmid - gaps / 2; g_hi = tmid + gaps / 2
    # 24 iç nokta gridi (ham ve bölünmüş aynı gridle — elma-elma)
    u = np.arange(1, 25) / 25
    tt = g_lo[:, None] + gaps[:, None] * u[None, :]
    print(f"\n=== Pencere {wkey} (L={L:.2f}): Z değerlendiriliyor... ===")
    Zg = np.abs(Z_rs(tt.ravel())).reshape(tt.shape)

    # CUE referansı N=⌊L⌋, ⌊L⌋+1 → N=L'ye interp
    N0 = int(np.floor(L)); w = L - N0
    (cM0, cV0), rC0 = cue_ref(N0)
    (cM1, cV1), rC1 = cue_ref(N0 + 1)
    cV = (1 - w) * cV0 + w * cV1
    r_cue = (1 - w) * rC0 + w * rC1

    g_u = gaps * np.log(tmid / TWO_PI) / TWO_PI
    rows = []
    for primes, label in PSETS:
        E = euler_factor(tt.ravel(), primes).reshape(tt.shape) if primes else 1.0
        amps = (Zg * E).max(axis=1)
        a_u = amps / np.sqrt((amps**2).mean())
        r, _ = pearsonr(g_u, a_u)
        zM, zV = bin_stats(g_u, a_u)
        ok = ~np.isnan(zV) & ~np.isnan(cV) & (BC < 2.2)  # uç bin hariç (43 uyarısı)
        big = ok & (BC > 1.3)
        dV_big = np.nanmean((zV - cV)[big])
        rows.append((label, r, dV_big))
        print(f"  {label:>5}: r = {r:.4f}   büyük-boşluk ΔVar = {dV_big:+.4f}")
    print(f"  CUE(N={L:.2f}) referans r = {r_cue:.4f}")
    results[wkey] = (L, rows, r_cue)

# ---------------------------------------------------------------
# GRAFİK
# ---------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(12.5, 5))
for ax, wkey in zip(axes, WINDOWS):
    L, rows, r_cue = results[wkey]
    labels = [r[0] for r in rows]
    rs = [r[1] for r in rows]
    dvs = [r[2] for r in rows]
    x = np.arange(len(rows))
    ax2 = ax.twinx()
    ax.bar(x - 0.18, rs, 0.36, color="firebrick", label="r (sol eksen)")
    ax.axhline(r_cue, color="steelblue", ls="--", lw=1.5,
               label=f"CUE hedefi r={r_cue:.3f}")
    ax2.bar(x + 0.18, dvs, 0.36, color="purple", alpha=0.7,
            label="büyük-boşluk ΔVar (sağ)")
    ax2.axhline(0, color="k", lw=0.7)
    ax.set_xticks(x); ax.set_xticklabels(labels)
    ax.set_ylim(min(rs) - 0.01, max(max(rs), r_cue) + 0.01)
    ax.set_ylabel("Pearson r"); ax2.set_ylabel("ΔVar (g̃>1.3)")
    ax.set_title(f"L = {L:.2f} — Euler bölme dizisi")
    ax.legend(fontsize=8, loc="upper left"); ax2.legend(fontsize=8, loc="upper right")
plt.tight_layout()
out = HERE / "44_euler_bolme.png"
plt.savefig(out, dpi=110)
print(f"\nGrafik: {out.name}")
