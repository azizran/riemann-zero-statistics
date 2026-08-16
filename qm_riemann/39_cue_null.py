"""
39 — CUE NULL: GAP-GENLİK KORELASYONU RMT-EVRENSEL Mİ? (16 Ağustos 2026)
=========================================================================

Soru (Mayıs'ın 2. açık adımı): r=0.83'ün ne kadarı "beklenen"?

29'daki Gaussian-PSD null spektrumu tutuyordu ama sıfır-yapısını tutmuyordu
(r_null = 0.505, gerçek 0.825, +6.67σ). Daha keskin null: CUE karakteristik
polinomu — Keating–Snaith'ten beri ζ'nın kritik doğru modeli. GUE itmesini
VE genlik yapısını birlikte içerir.

Model: N×N Haar-random üniter U, özfaz φ_j; |Λ(θ)| = Π_k 2|sin((θ-φ_k)/2)|.
Yükseklik t'deki ζ ↔ N ≈ L = log(t/2π) boyutlu CUE (yoğunluk eşleme).

Test: her N için (özfaz aralığı, aralıktaki max|Λ|) Pearson r'si.
ζ tarafı: 36_T100k verisi, log-pencerelerde aynı r (pencere içinde L≈sabit,
normalizasyon Pearson'ı değiştirmez).

  r_CUE(N≈L) ≈ r_ζ(L)  → korelasyon RMT-evrensel (asal gerekmez)
  r_CUE ≠ r_ζ           → fark ζ-özgü (aritmetik) yapı
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import pearsonr, spearmanr

rng = np.random.default_rng(39)
HERE = Path(__file__).resolve().parent

# ---------------------------------------------------------------
# 1) ZETA TARAFI — pencere bazlı r (36 verisi)
# ---------------------------------------------------------------
d = np.load(HERE / "36_T100k.npz")
gaps_z, amps_z, t_mid = d["intervals"], d["max_amps"], d["t_mid"]
Lz = np.log(t_mid / (2 * np.pi))

n_w = 12
edges = np.geomspace(t_mid[0], t_mid[-1] * 1.0001, n_w + 1)
zw = []
for i in range(n_w):
    m = (t_mid >= edges[i]) & (t_mid < edges[i + 1])
    if m.sum() < 500:
        continue
    r, _ = pearsonr(gaps_z[m], amps_z[m])
    rs, _ = spearmanr(gaps_z[m], amps_z[m])
    se = (1 - r**2) / np.sqrt(m.sum())
    zw.append((Lz[m].mean(), r, rs, se, m.sum()))
zw_L, zw_r, zw_rs, zw_se, zw_n = map(np.array, zip(*zw))

# tüm veri, normalize aralıkla (27 konvansiyonu)
rho = Lz / (2 * np.pi)
r_all, _ = pearsonr(gaps_z * rho, amps_z)
r_raw, _ = pearsonr(gaps_z, amps_z)
print("ZETA (100k aralık, T≤75k):")
print(f"  tüm veri: r_raw = {r_raw:.4f}, r_norm = {r_all:.4f}")
print(f"  pencere r aralığı: {zw_r.min():.3f} … {zw_r.max():.3f}  (L: {zw_L.min():.1f}…{zw_L.max():.1f})\n")

# ---------------------------------------------------------------
# 2) CUE TARAFI
# ---------------------------------------------------------------
def haar_unitary_batch(M, N):
    G = (rng.standard_normal((M, N, N)) + 1j * rng.standard_normal((M, N, N))) / np.sqrt(2)
    Q, R = np.linalg.qr(G)
    diag = np.einsum("mii->mi", R)
    return Q * (diag / np.abs(diag))[:, None, :]

def cue_gap_max(N, n_gaps_target=60000, grid=48, chunk=1500):
    """Her matriste N (döngüsel) aralık: (aralık, max|Λ|) çiftleri."""
    M = int(np.ceil(n_gaps_target / N))
    all_gaps, all_maxs = [], []
    for start in range(0, M, chunk):
        m = min(chunk, M - start)
        U = haar_unitary_batch(m, N)
        phases = np.sort(np.angle(np.linalg.eigvals(U)), axis=1)  # (m,N) ∈ (-π,π]
        # döngüsel aralıklar
        gaps = np.diff(np.concatenate([phases, phases[:, :1] + 2 * np.pi], axis=1), axis=1)  # (m,N)
        # aralık içi grid (uçlar hariç)
        u = (np.arange(1, grid + 1) / (grid + 1))  # (grid,)
        theta = phases[:, :, None] + gaps[:, :, None] * u[None, None, :]  # (m,N,grid)
        # |Λ(θ)| = Π_k 2|sin((θ-φ_k)/2)|
        diff = theta[:, :, :, None] - phases[:, None, None, :]  # (m,N,grid,N)
        amp = np.prod(2 * np.abs(np.sin(diff / 2)), axis=-1)  # (m,N,grid)
        all_gaps.append(gaps.ravel())
        all_maxs.append(amp.max(axis=-1).ravel())
    return np.concatenate(all_gaps), np.concatenate(all_maxs)

print("CUE (her N için ~60k aralık):")
Ns = np.arange(5, 15)
cue_r, cue_rs, cue_se = [], [], []
for N in Ns:
    g, a = cue_gap_max(N)
    r, _ = pearsonr(g, a)
    rs, _ = spearmanr(g, a)
    se = (1 - r**2) / np.sqrt(len(g))
    cue_r.append(r); cue_rs.append(rs); cue_se.append(se)
    print(f"  N={N:>2}: Pearson r = {r:.4f} ± {se:.4f}   Spearman = {rs:.4f}   ({len(g)} aralık)")
cue_r, cue_rs, cue_se = map(np.array, (cue_r, cue_rs, cue_se))

# ---------------------------------------------------------------
# 3) EŞLEŞMİŞ KARŞILAŞTIRMA: r_ζ(L) vs r_CUE(N=L) (lineer interp.)
# ---------------------------------------------------------------
cue_at_L = np.interp(zw_L, Ns, cue_r)
diff = zw_r - cue_at_L
zsc = diff / zw_se
print("\nEŞLEŞMİŞ KARŞILAŞTIRMA (pencere bazında):")
print(f"  {'L':>5} {'r_ζ':>8} {'r_CUE':>8} {'fark':>8} {'z':>6}")
for i in range(len(zw_L)):
    print(f"  {zw_L[i]:>5.2f} {zw_r[i]:>8.4f} {cue_at_L[i]:>8.4f} {diff[i]:>+8.4f} {zsc[i]:>6.1f}")
print(f"\n  ortalama fark = {diff.mean():+.4f},  ortalama |z| = {abs(zsc).mean():.1f}")

# ---------------------------------------------------------------
# 4) GRAFİK
# ---------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(12.5, 5))
ax = axes[0]
ax.errorbar(zw_L, zw_r, yerr=2 * zw_se, fmt="o", ms=5, c="firebrick",
            label="ζ verisi (pencereler, ±2σ)")
ax.errorbar(Ns, cue_r, yerr=2 * cue_se, fmt="s-", ms=4, c="steelblue",
            label="CUE null (N ↔ L)")
ax.set_xlabel("L = log(t/2π)  ↔  N (CUE boyutu)")
ax.set_ylabel("Pearson r (aralık, max genlik)")
ax.set_title("Gap-genlik korelasyonu: ζ vs CUE")
ax.legend(fontsize=9)
ax.grid(alpha=0.3)

ax = axes[1]
ax.errorbar(zw_L, diff, yerr=2 * zw_se, fmt="o", ms=5, c="purple")
ax.axhline(0, color="k", lw=1)
ax.set_xlabel("L = log(t/2π)")
ax.set_ylabel("r_ζ − r_CUE")
ax.set_title("Fark (0 = RMT-evrensel)")
ax.grid(alpha=0.3)

plt.tight_layout()
out = HERE / "39_cue_null.png"
plt.savefig(out, dpi=110)
print(f"\nGrafik: {out.name}")

# ---------------------------------------------------------------
# 5) KARAR
# ---------------------------------------------------------------
print("\n" + "=" * 62)
if abs(zsc).mean() < 2:
    print("SONUÇ: r_ζ ≈ r_CUE → korelasyon RMT-EVRENSEL görünüyor.")
    print("Yeni katkı çerçevesi: 'ζ, gap-genlik kuplajında RMT öngörüsünü")
    print("X% hassasiyetle izliyor' (yine değerli ama farklı cümle).")
else:
    yon = "ÜSTÜNDE" if diff.mean() > 0 else "ALTINDA"
    print(f"SONUÇ: r_ζ, CUE öngörüsünün sistematik {yon} → ζ-ÖZGÜ fark var!")
    print("Aritmetik (asal) katkı adayı — makalenin ana cümlesi bu olabilir.")
