"""
40 — UNFOLD + N_eff: 39'DAKİ FARK GERÇEK Mİ? (16 Ağustos 2026)
================================================================

39'da r_ζ − r_CUE ≈ −0.022 çıktı. İki şüphe:

(a) Pencere-içi kayma: pencerede t ~2× değişiyor; ortalama aralık ∝ 2π/L
    DÜŞERKEN genlik ölçeği ∝ √(AL+b) YÜKSELİYOR → zıt kayma r_ζ'yi yapay
    düşürür. Çözüm: UNFOLD — g̃ = gap·L/2π, ã = amp/√(AL + 2.758 − 0.054/L)
    (38'in doğrulanmış HLP-C ölçeği!). Kayma ölür, saf dalgalanma kalır.

(b) N↔L eşleme belirsizliği: dr_CUE/dN ≈ −0.022 → ΔN=1 konvansiyon farkı
    tüm farkı açıklayabilir. Çözüm: N_eff(L) tanımla — r_CUE(N_eff) = r_ζ(L)
    olacak şekilde. N_eff − L SABIT çıkarsa "kalibrasyon farkı" (yapısal
    iddia yok); L ile DEĞİŞİRSE gerçek şekil farkı var demektir.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import pearsonr

rng = np.random.default_rng(40)
HERE = Path(__file__).resolve().parent

A = (np.e**2 - 5) / 2
B0, B1 = 2.7580, -0.0543  # 38'den (HLP-C, doğrulanmış)

# ---------------------------------------------------------------
# 1) ZETA — UNFOLD EDİLMİŞ pencere r'leri
# ---------------------------------------------------------------
d = np.load(HERE / "36_T100k.npz")
t_mid = d["t_mid"]
L = np.log(t_mid / (2 * np.pi))
g_u = d["intervals"] * L / (2 * np.pi)                  # ortalama ≈ 1
a_u = d["max_amps"] / np.sqrt(A * L + B0 + B1 / L)      # ortalama² ≈ 1

print("Unfold kontrol: mean(g̃) = %.4f, mean(ã²) = %.4f" % (g_u.mean(), (a_u**2).mean()))

def win_r(gv, av, tv, n_w):
    edges = np.geomspace(tv[0], tv[-1] * 1.0001, n_w + 1)
    out = []
    for i in range(n_w):
        m = (tv >= edges[i]) & (tv < edges[i + 1])
        if m.sum() < 500:
            continue
        r, _ = pearsonr(gv[m], av[m])
        out.append((np.log(tv[m].mean() / (2 * np.pi)), r,
                    (1 - r**2) / np.sqrt(m.sum()), m.sum()))
    return map(np.array, zip(*out))

zw_L, zw_r, zw_se, zw_n = win_r(g_u, a_u, t_mid, 12)
zw_L24, zw_r24, zw_se24, _ = win_r(g_u, a_u, t_mid, 24)

r_pool, _ = pearsonr(g_u, a_u)
print(f"Unfold sonrası havuz r = {r_pool:.4f} (39'daki norm r: 0.8196)")
print(f"Pencere r (12 pencere): {zw_r.round(4)}\n")

# ---------------------------------------------------------------
# 2) CUE — daha sıkı istatistik (100k aralık/N), N=5..16
# ---------------------------------------------------------------
def haar_unitary_batch(M, N):
    G = (rng.standard_normal((M, N, N)) + 1j * rng.standard_normal((M, N, N))) / np.sqrt(2)
    Q, R = np.linalg.qr(G)
    diag = np.einsum("mii->mi", R)
    return Q * (diag / np.abs(diag))[:, None, :]

def cue_r(N, n_gaps=100000, grid=48, chunk=1500):
    M = int(np.ceil(n_gaps / N))
    gs, ms = [], []
    for s in range(0, M, chunk):
        m = min(chunk, M - s)
        U = haar_unitary_batch(m, N)
        ph = np.sort(np.angle(np.linalg.eigvals(U)), axis=1)
        gaps = np.diff(np.concatenate([ph, ph[:, :1] + 2 * np.pi], axis=1), axis=1)
        u = np.arange(1, grid + 1) / (grid + 1)
        th = ph[:, :, None] + gaps[:, :, None] * u[None, None, :]
        df = th[:, :, :, None] - ph[:, None, None, :]
        amp = np.prod(2 * np.abs(np.sin(df / 2)), axis=-1)
        gs.append(gaps.ravel()); ms.append(amp.max(axis=-1).ravel())
    g = np.concatenate(gs); a = np.concatenate(ms)
    r, _ = pearsonr(g, a)
    return r, (1 - r**2) / np.sqrt(len(g))

print("CUE r(N):")
Ns = np.arange(5, 17)
cue = np.array([cue_r(N) for N in Ns])
for N, (r, se) in zip(Ns, cue):
    print(f"  N={N:>2}: r = {r:.4f} ± {se:.4f}")
cue_r_arr, cue_se_arr = cue[:, 0], cue[:, 1]

# ---------------------------------------------------------------
# 3) KARŞILAŞTIRMA + N_eff
# ---------------------------------------------------------------
cue_at = np.interp(zw_L, Ns, cue_r_arr)
diff = zw_r - cue_at
zsc = diff / zw_se

# N_eff: r_CUE(N_eff) = r_ζ  (monoton azalan → ters interp)
N_eff = np.interp(-zw_r, -cue_r_arr, Ns.astype(float))
shift = N_eff - zw_L

print("\nUNFOLD SONRASI KARŞILAŞTIRMA:")
print(f"  {'L':>5} {'r_ζ':>8} {'r_CUE(L)':>9} {'fark':>8} {'z':>6} {'N_eff':>6} {'N_eff−L':>8}")
for i in range(len(zw_L)):
    print(f"  {zw_L[i]:>5.2f} {zw_r[i]:>8.4f} {cue_at[i]:>9.4f} {diff[i]:>+8.4f} "
          f"{zsc[i]:>6.1f} {N_eff[i]:>6.2f} {shift[i]:>+8.2f}")
print(f"\n  ortalama fark = {diff.mean():+.4f} (39'da −0.0222 idi)")
print(f"  N_eff − L: ort = {shift.mean():+.3f}, std = {shift.std():.3f}")

# N_eff−L trendi (sabit mi?)
sl = np.polyfit(zw_L, shift, 1)
print(f"  N_eff−L eğimi (L'ye karşı): {sl[0]:+.3f} / birim L")

# ---------------------------------------------------------------
# 4) GRAFİK
# ---------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.6))
ax = axes[0]
ax.errorbar(zw_L, zw_r, yerr=2 * zw_se, fmt="o", ms=5, c="firebrick",
            label="ζ unfold (12 pencere)")
ax.errorbar(zw_L24, zw_r24, yerr=2 * zw_se24, fmt=".", ms=3, c="salmon",
            alpha=0.6, label="ζ unfold (24 pencere)")
ax.errorbar(Ns, cue_r_arr, yerr=2 * cue_se_arr, fmt="s-", ms=4, c="steelblue",
            label="CUE")
ax.set_xlabel("L  ↔  N"); ax.set_ylabel("Pearson r")
ax.set_title("Unfold sonrası: ζ vs CUE"); ax.legend(fontsize=8); ax.grid(alpha=0.3)

ax = axes[1]
ax.errorbar(zw_L, diff, yerr=2 * zw_se, fmt="o", ms=5, c="purple")
ax.axhline(0, color="k", lw=1)
ax.set_xlabel("L"); ax.set_ylabel("r_ζ − r_CUE(N=L)")
ax.set_title(f"Fark (ort {diff.mean():+.4f})"); ax.grid(alpha=0.3)

ax = axes[2]
ax.plot(zw_L, shift, "o-", c="darkgreen")
ax.axhline(shift.mean(), color="gray", ls="--",
           label=f"ort = {shift.mean():+.2f}")
ax.set_xlabel("L"); ax.set_ylabel("N_eff − L")
ax.set_title("Eşleme kayması: sabit mi?"); ax.legend(fontsize=9); ax.grid(alpha=0.3)

plt.tight_layout()
out = HERE / "40_unfold_recheck.png"
plt.savefig(out, dpi=110)
print(f"\nGrafik: {out.name}")
