"""
67 — İKİ-AKIŞKAN MODELİ: U = DC PERDELEME + w'NİN NEDENSELLİK-EŞİ
===================================================================
(18 Ağustos 2026 — Sıçrama 1, 2. tur)

66b'nin reti "tek skaler dielektrik" içindi. Parite çelişkisinin (U(0)≠0)
fizikteki standart çözümü: İKİ-AKIŞKAN — iletken/süperakışkan bileşen
(DC, ω→0'da bile yanıt verir) + normal bileşen (nedensel gevşeme).

TEST A: U(τ) ≈ c₀ + μ·Im χ_w(τ)
  Im χ_w tamamen w VERİSİNDEN inşa edilir (pozitif Debye fiti → nedensellik-eşi).
  Sadece 2 serbest parametre (c₀, μ). Tutarsa: konum kanalı =
  sabit perdeleme + genlik kanalının KK-partneri → yapı kurtulur.

TEST B: kuplaj artığının şekli, sum-rule sapması d(τ) = w + βv − 1 ile
  aynı mı? (44 eşleşmiş (w,v) çifti; 55 pencerelerinin v'si burada ölçülür.)
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.optimize import lsq_linear

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
_src = open(HERE / "66_kramers_kronig.py").read()
exec(_src[:_src.index("def chi(")])

tw, w, sw = W_DATA.T
tv, v, sv = V_DATA.T
U = v / tv
sU = sv / tv

# ---- w'nin pozitif Debye spektrumu (66b ile aynı) ----
S = np.geomspace(0.02, 2.0, 25)
def design(tau, kind):
    x = tau[:, None] / S[None, :]
    return 1.0 / (1.0 + x**2) if kind == "re" else x / (1.0 + x**2)

def solve_nn(A, y, ridge=1e-3):
    n = A.shape[1]
    Ar = np.vstack([A, ridge * np.eye(n)])
    yr = np.concatenate([y, np.zeros(n)])
    return lsq_linear(Ar, yr, bounds=(0, np.inf), method="bvls", max_iter=300).x

Aw = np.hstack([design(tw, "re"), np.ones((len(tw), 1)), -np.ones((len(tw), 1))])
rho_w = solve_nn(Aw / sw[:, None], w / sw)
rms_w = np.sqrt(np.mean((Aw @ rho_w - w) ** 2))
print(f"w-spektrumu fiti: RMS = {rms_w:.4f}")

Im_w = design(tv, "im") @ rho_w[:25]   # nedensellik-eşi, v-grid üstünde

# ---- TEST A: U = c0 + mu * Im_w ----
X = np.vstack([np.ones_like(tv), Im_w]).T
Wt = 1 / sU**2
beta_ls, *_ = np.linalg.lstsq(X * np.sqrt(Wt)[:, None], U * np.sqrt(Wt), rcond=None)
c0, mu = beta_ls
res = U - X @ beta_ls
chi2 = np.sum(Wt * res**2)
# kıyas tavanı: U'nun kendi serbest Debye fiti (66b: chi2/dof ~ 1.1)
Au = np.hstack([design(tv, "im"), np.ones((len(tv), 1)), -np.ones((len(tv), 1))])
ru = solve_nn(Au / sU[:, None], U / sU)
chi2_free = np.sum(((Au @ ru - U) / sU) ** 2)
print(f"\nTEST A — iki-akışkan: U = c₀ + μ·Imχ_w")
print(f"  c₀ = {c0:.3f}  (DC perdeleme; U(0)≈2.01 ile kıyasla)")
print(f"  μ  = {mu:.3f}")
print(f"  χ²/dof = {chi2/(len(tv)-2):.2f}   (serbest 27-parametreli tavan: "
      f"{chi2_free/(len(tv)-27):.2f})")
print(f"  U-RMS = {np.sqrt(np.mean(res**2)):.4f} (serbest fit RMS: "
      f"{np.sqrt(np.mean((Au@ru - U)**2)):.4f})")

# ---- TEST B: eşleşmiş sum-rule sapması d(τ) ----
PRIMES4 = [2, 3, 5, 7]
V47 = {9.86:[0.1405,0.2245,0.3254,0.3877],10.37:[0.1358,0.2159,0.3108,0.3713],
       10.93:[0.1279,0.2034,0.2962,0.3576],11.47:[0.1222,0.1950,0.2823,0.3403],
       11.98:[0.1174,0.1847,0.2701,0.3258],12.45:[0.1114,0.1787,0.2628,0.3123],
       24.48:[0.054,0.091,0.126,0.163]}
# 55 pencerelerinin v'si (anlık regresyon)
for f, key in [("55_win_1e+08.npz","16.58"),("55_win_1e+09.npz","18.89"),
               ("55_win_1e+10.npz","21.19"),("55_win_1e+11.npz","23.49")]:
    d = np.load(HERE / f)
    gaps, tmid = d["gaps"], d["tmid"]
    Lw = np.log(tmid / TWO_PI)
    y = np.log(gaps * Lw / TWO_PI)
    cols = [np.ones_like(y)]
    for p in PRIMES4:
        arg = tmid * np.log(p)
        cols += [np.cos(arg), np.sin(arg)]
    b, *_ = np.linalg.lstsq(np.vstack(cols).T, y, rcond=None)
    V47[float(key)] = [float(np.hypot(b[1+2*i], b[2+2*i]) / p**-0.5)
                       for i, p in enumerate(PRIMES4)]

# w tablosundan eşleşmiş (tau, w, v) üçlüleri
WMAP = {}  # (p, L) -> (tau, w)
WL = [(2,24.48),(2,23.49),(2,21.19),(2,18.89),(2,16.58),(3,24.48),(3,23.49),
      (3,21.19),(2,12.45),(2,11.98),(3,18.89),(2,11.47),(2,10.93),(5,24.48),
      (3,16.58),(2,10.37),(5,23.49),(2,9.86),(5,21.19),(7,24.48),(7,23.49),
      (5,18.89),(3,12.45),(3,11.98),(7,21.19),(3,11.47),(5,16.58),(3,10.93),
      (7,18.89),(3,10.37),(3,9.86),(7,16.58),(5,12.45),(5,11.98),(5,11.47),
      (5,10.93),(5,10.37),(7,12.45),(7,11.98),(5,9.86),(7,11.47),(7,10.93),
      (7,10.37),(7,9.86)]
pairs = []
for (p, L), (t_, w_, s_) in zip(WL, W_DATA):
    if L in V47:
        vv = V47[L][PRIMES4.index(p)]
        pairs.append((t_, w_, vv))
pairs = np.array(sorted(pairs))
tp, wp, vp = pairs.T
BETA = 2.14  # faz muhasebesinden (U(0)/w(0))
d_sum = wp + BETA * vp - 1.0
# kuplaj artığı: iki-akışkan modelin aynı noktalardaki artığı
Im_w_p = design(tp, "im") @ rho_w[:25]
res_p = vp / tp - (c0 + mu * Im_w_p)
print(f"\nTEST B — {len(pairs)} eşleşmiş çift, β = {BETA} (faz muhasebesi):")
print(f"  d(τ) küçük-τ ort: {d_sum[tp<0.06].mean():+.3f}  büyük-τ ort: "
      f"{d_sum[tp>0.12].mean():+.3f}")
if np.std(res_p) > 0:
    cc = np.corrcoef(d_sum, res_p)[0, 1]
    print(f"  corr(d, iki-akışkan artığı) = {cc:+.3f}")

# ---- grafik ----
fig, axes = plt.subplots(1, 2, figsize=(13, 5.2))
tt = np.linspace(1e-3, 0.5, 400)
Im_tt = design(tt, "im") @ rho_w[:25]
ax = axes[0]
ax.errorbar(tv, U, yerr=sU, fmt="s", ms=3.5, c="teal", label="U = v/τ (90)")
ax.plot(tt, c0 + mu * Im_tt, "k-", lw=1.4,
        label=f"iki-akışkan: {c0:.2f} + {mu:.2f}·Imχ_w")
ax.plot(tt, np.full_like(tt, c0), "--", c="gray", lw=1, label=f"DC bileşen c₀={c0:.2f}")
ax.set_xlabel("τ"); ax.set_ylabel("U"); ax.legend(fontsize=9); ax.grid(alpha=0.3)
ax.set_title("Test A: konum kanalı = DC + w'nin KK-eşi")
ax = axes[1]
ax.plot(tp, d_sum, "o", ms=4, c="firebrick", label=f"d = w + {BETA}·v − 1")
ax.plot(tp, res_p, "s", ms=4, c="purple", label="iki-akışkan artığı (U-modeli)")
ax.axhline(0, color="k", lw=0.7)
ax.set_xlabel("τ"); ax.set_ylabel("sapma")
ax.set_title("Test B: iki anomali aynı şekil mi?")
ax.legend(fontsize=9); ax.grid(alpha=0.3)
plt.tight_layout()
out = HERE / "67_iki_akiskan.png"
plt.savefig(out, dpi=110)
print(f"\nGrafik: {out.name}")
