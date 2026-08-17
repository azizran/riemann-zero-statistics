"""
60 — N_eff NİHAİ KARARI (K1 çözümü) (17-18 Ağustos gecesi)
============================================================

59'un hassas CUE tablosu (1.5M aralık/N, ampirik hatalar) + PCHIP
(form varsayımı YOK) ile N_eff(L) yeniden türetilir.

KARAR SORUSU: N_eff − L sabit mi (eski iddia), L ile yükseliyor mu
(denetçinin bulgusu)? Ve iki-metrik yarılması (P vs S farkı ~0.9)
interpolant değişince hayatta kalıyor mu?

Hatalar: ζ r se'leri + CUE interpolasyon hatası (se/|eğim|) birlikte.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import pearsonr, spearmanr
from scipy.interpolate import PchipInterpolator

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
A = (np.e**2 - 5) / 2
B0, B1 = 2.7580, -0.0543

# ---------------------------------------------------------------
# 1) ζ pencereleri (42 ile aynı: unfold + r)
# ---------------------------------------------------------------
def unfolded_r(gaps, amps, tmid):
    L = np.log(tmid / TWO_PI)
    g = gaps * L / TWO_PI
    a = amps / np.sqrt(A * L + B0 + B1 / L)
    rP, _ = pearsonr(g, a)
    rS, _ = spearmanr(g, a)
    n = len(g)
    return (float(L.mean()), rP, (1 - rP**2) / np.sqrt(n),
            rS, (1 - rS**2) / np.sqrt(n))

rows = []
d36 = np.load(HERE / "36_T100k.npz")
edges = np.geomspace(d36["t_mid"][0], d36["t_mid"][-1] * 1.0001, 13)
for i in range(12):
    m = (d36["t_mid"] >= edges[i]) & (d36["t_mid"] < edges[i + 1])
    if m.sum() < 500:
        continue
    rows.append(unfolded_r(d36["intervals"][m], d36["max_amps"][m], d36["t_mid"][m]))
d41 = np.load(HERE / "41_bigT_windows.npz")
for k in sorted({x.split("_")[1] for x in d41.files}, key=lambda s: int(s[:-1])):
    rows.append(unfolded_r(d41[f"gaps_{k}"], d41[f"amps_{k}"], d41[f"tmid_{k}"]))
zw = np.array(rows)  # L, rP, seP, rS, seS

# ---------------------------------------------------------------
# 2) Hassas CUE + PCHIP
# ---------------------------------------------------------------
c = np.load(HERE / "59_cue_hassas.npz")
Ns, rP_c, seP_c, rS_c, seS_c = c["N"], c["rP"], c["seP"], c["rS"], c["seS"]
print("CUE tablosu:", " ".join(f"{n}:{r:.4f}" for n, r in zip(Ns, rP_c)))
fP = PchipInterpolator(Ns, rP_c)
fS = PchipInterpolator(Ns, rS_c)

def invert(f, r, lo=5.0, hi=20.0):
    xs = np.linspace(lo, hi, 4000)
    ys = f(xs)
    i = np.argmin(np.abs(ys - r))
    return xs[i]

def slope(f, x, h=0.05):
    return (f(x + h) - f(x - h)) / (2 * h)

# ---------------------------------------------------------------
# 3) N_eff ve trend analizi
# ---------------------------------------------------------------
print(f"\n{'L':>6} {'N_effP−L':>9} {'±':>5} {'N_effS−L':>9} {'±':>5} {'S−P farkı':>9}")
res = []
for L, rP, seP, rS, seS in zw:
    NP = invert(fP, rP); NS_ = invert(fS, rS)
    sl_P = abs(slope(fP, NP)); sl_S = abs(slope(fS, NS_))
    # hata: ζ se + CUE interp hatası (yerel se'lerin ortalaması)
    cue_seP = np.interp(NP, Ns, seP_c); cue_seS = np.interp(NS_, Ns, seS_c)
    eP = np.hypot(seP, cue_seP) / sl_P
    eS = np.hypot(seS, cue_seS) / sl_S
    res.append((L, NP - L, eP, NS_ - L, eS))
    print(f"{L:>6.2f} {NP-L:>+9.3f} {eP:>5.2f} {NS_-L:>+9.3f} {eS:>5.2f} "
          f"{NS_-NP:>+9.3f}")
res = np.array(res)

def trend(shift, err, label):
    w = 1 / err**2
    cw = np.sum(w * shift) / np.sum(w)
    chi_c = np.sum(w * (shift - cw)**2)
    p = np.polyfit(res[:, 0], shift, 1, w=np.sqrt(w))
    chi_l = np.sum(w * (shift - np.polyval(p, res[:, 0]))**2)
    n = len(shift)
    # eğim hatası
    X = np.vstack([res[:, 0], np.ones(n)]).T
    cov = np.linalg.inv((X * w[:, None]).T @ X)
    se_slope = np.sqrt(cov[0, 0])
    print(f"\n{label}:")
    print(f"  SABİT : {cw:+.3f}  χ²/dof = {chi_c/(n-1):.2f}")
    print(f"  LİNEER: eğim {p[0]:+.4f} ± {se_slope:.4f} /L  χ²/dof = {chi_l/(n-2):.2f}")
    print(f"  Δχ² = {chi_c-chi_l:.1f} (1 dof) → eğim anlamlılığı {abs(p[0])/se_slope:.1f}σ")
    return cw, p

cwP, pP = trend(res[:, 1], res[:, 2], "PEARSON N_eff − L")
cwS, pS = trend(res[:, 3], res[:, 4], "SPEARMAN N_eff − L")

diff = res[:, 3] - res[:, 1]
print(f"\nİKİ-METRİK YARILMASI (S − P): ort {diff.mean():+.3f}, std {diff.std():.3f}")
print("(Bu fark interpolant-bağımsız yapısal bulgu — hayatta kalması beklenir)")

# ---------------------------------------------------------------
# 4) Grafik
# ---------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8.8, 5.6))
ax.errorbar(res[:, 0], res[:, 1], yerr=res[:, 2], fmt="o", ms=5, c="firebrick",
            label="Pearson")
ax.errorbar(res[:, 0], res[:, 3], yerr=res[:, 4], fmt="s", ms=5, c="teal",
            label="Spearman")
xs = np.linspace(res[:, 0].min(), res[:, 0].max(), 50)
ax.plot(xs, np.polyval(pP, xs), "--", c="firebrick", lw=1)
ax.plot(xs, np.polyval(pS, xs), "--", c="teal", lw=1)
ax.set_xlabel("L = log(t/2π)"); ax.set_ylabel("N_eff − L")
ax.set_title("N_eff nihai: hassas CUE + PCHIP (form varsayımsız)")
ax.legend(); ax.grid(alpha=0.3)
plt.tight_layout()
out = HERE / "60_neff_nihai.png"
plt.savefig(out, dpi=110)
print(f"\nGrafik: {out.name}")
