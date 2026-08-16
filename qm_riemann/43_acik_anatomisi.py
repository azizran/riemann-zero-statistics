"""
43 — AÇIĞIN ANATOMİSİ: −0.019 NEREDE YAŞIYOR? (16 Ağustos 2026, gece)
======================================================================

r_ζ − r_CUE(N=L) = −0.019. İki aday mekanizma:

(a) KOŞULLU ORTALAMA: E[ã | g̃] eğrisi ζ'da daha yatık
    → boşluk-tepe bağının kendisi zayıf
(b) KOŞULLU VARYANS: eğri aynı, ama ζ'da boşluktan BAĞIMSIZ ekstra
    genlik varyansı var → asal-salınım mekanizması (explicit formula
    dalgaları tek aralıktan uzun periyotlu: p=2 için periyot 2π/log2 ≈ 9.1,
    L=10'da ~15 aralık — aralığa göre "sabit arka plan" gibi davranır)

Öngörü (b) ise: Var_ζ(ã|g̃) − Var_CUE(ã|g̃) ≈ SABİT > 0 tüm g̃ binlerinde.

Test: 3 ζ penceresi (L ≈ 7.0, 9.9, 12.4) vs CUE (N=L'ye interp edilmiş
koşullu eğriler). Her iki taraf da unfold edilmiş: mean(g̃)=1, mean(ã²)=1.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

rng = np.random.default_rng(43)
HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
A = (np.e**2 - 5) / 2
B0, B1 = 2.7580, -0.0543

BINS = np.concatenate([[0.05], np.linspace(0.25, 2.2, 14), [3.5]])  # g̃ bin kenarları
BC = 0.5 * (BINS[:-1] + BINS[1:])

def bin_stats(g, a):
    """g̃ binlerinde koşullu ortalama ve varyans."""
    mean = np.full(len(BC), np.nan)
    var = np.full(len(BC), np.nan)
    cnt = np.zeros(len(BC), dtype=int)
    idx = np.digitize(g, BINS) - 1
    for i in range(len(BC)):
        m = idx == i
        if m.sum() >= 200:
            mean[i] = a[m].mean()
            var[i] = a[m].var()
            cnt[i] = m.sum()
    return mean, var, cnt

# ---------------------------------------------------------------
# 1) ZETA pencereleri
# ---------------------------------------------------------------
d36 = np.load(HERE / "36_T100k.npz")
d41 = np.load(HERE / "41_bigT_windows.npz")

def zeta_window(gaps, amps, tmid):
    L = np.log(tmid / TWO_PI)
    g = gaps * L / TWO_PI
    a = amps / np.sqrt(A * L + B0 + B1 / L)
    a /= np.sqrt((a**2).mean())  # tam normalize (şekil kıyası için)
    return L.mean(), g, a

# L≈7.0 penceresi (36'dan), L≈9.9 ve 12.4 (41'den)
edges = np.geomspace(d36["t_mid"][0], d36["t_mid"][-1] * 1.0001, 13)
m7 = (d36["t_mid"] >= edges[8]) & (d36["t_mid"] < edges[9])
zw = [
    zeta_window(d36["intervals"][m7], d36["max_amps"][m7], d36["t_mid"][m7]),
    zeta_window(d41["gaps_120k"], d41["amps_120k"], d41["tmid_120k"]),
    zeta_window(d41["gaps_1600k"], d41["amps_1600k"], d41["tmid_1600k"]),
]
print("ζ pencereleri:", [f"L={w[0]:.2f} (n={len(w[1])})" for w in zw])

# ---------------------------------------------------------------
# 2) CUE koşullu eğrileri (N=6..14, 300k aralık/N)
# ---------------------------------------------------------------
def haar_unitary_batch(M, N):
    G = (rng.standard_normal((M, N, N)) + 1j * rng.standard_normal((M, N, N))) / np.sqrt(2)
    Q, R = np.linalg.qr(G)
    diag = np.einsum("mii->mi", R)
    return Q * (diag / np.abs(diag))[:, None, :]

def cue_curves(N, n_gaps=300000, grid=48, chunk=1200):
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
    g = np.concatenate(gs) * N / TWO_PI          # unfold
    a = np.concatenate(ms)
    a /= np.sqrt((a**2).mean())
    return bin_stats(g, a)

print("CUE eğrileri hesaplanıyor (N=6..14)...")
Ns = np.arange(6, 15)
cueM, cueV = {}, {}
for N in Ns:
    mN, vN, _ = cue_curves(N)
    cueM[N], cueV[N] = mN, vN
    print(f"  N={N} tamam")

def cue_interp(L, table):
    """N=L'ye lineer interpolasyon (komşu iki tam N'den)."""
    N0 = int(np.floor(L)); w = L - N0
    N0 = max(min(N0, Ns[-1] - 1), Ns[0])
    return (1 - w) * table[N0] + w * table[N0 + 1]

# ---------------------------------------------------------------
# 3) KARŞILAŞTIRMA
# ---------------------------------------------------------------
fig, axes = plt.subplots(2, 3, figsize=(15.5, 8.5))
print("\nSONUÇLAR:")
for j, (Lv, g, a) in enumerate(zw):
    zM, zV, zc = bin_stats(g, a)
    cM = cue_interp(Lv, cueM)
    cV = cue_interp(Lv, cueV)
    ok = ~np.isnan(zM) & ~np.isnan(cM)

    dM = zM - cM
    dV = zV - cV
    # ağırlıklı (bin nüfusu) özet
    w = zc[ok]
    print(f"\n  L={Lv:.2f}:")
    print(f"    ortalama eğri farkı  <ζ−CUE> : {np.average(dM[ok], weights=w):+.4f} "
          f"(bin RMS {np.sqrt(np.average(dM[ok]**2, weights=w)):.4f})")
    print(f"    varyans farkı  <Varζ−VarCUE> : {np.average(dV[ok], weights=w):+.4f}")
    print(f"    varyans farkı bin-bin        : "
          + " ".join(f"{x:+.3f}" for x in dV[ok]))

    ax = axes[0, j]
    ax.plot(BC[ok], zM[ok], "o-", ms=4, c="firebrick", label=f"ζ (L={Lv:.2f})")
    ax.plot(BC[ok], cM[ok], "s--", ms=3, c="steelblue", label=f"CUE (N={Lv:.2f})")
    ax.set_xlabel("g̃ (unfold aralık)"); ax.set_ylabel("E[ã | g̃]")
    ax.set_title("Koşullu ortalama"); ax.legend(fontsize=8); ax.grid(alpha=0.3)

    ax = axes[1, j]
    ax.plot(BC[ok], zV[ok], "o-", ms=4, c="firebrick", label="ζ")
    ax.plot(BC[ok], cV[ok], "s--", ms=3, c="steelblue", label="CUE")
    ax.plot(BC[ok], dV[ok], "^-", ms=4, c="purple", label="fark (ζ−CUE)")
    ax.axhline(0, color="k", lw=0.7)
    ax.set_xlabel("g̃"); ax.set_ylabel("Var[ã | g̃]")
    ax.set_title("Koşullu varyans"); ax.legend(fontsize=8); ax.grid(alpha=0.3)

plt.tight_layout()
out = HERE / "43_acik_anatomisi.png"
plt.savefig(out, dpi=110)
print(f"\nGrafik: {out.name}")
