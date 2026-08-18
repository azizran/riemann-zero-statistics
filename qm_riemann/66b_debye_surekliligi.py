"""
66b — DEBYE SÜREKLİLİĞİ: TEK POZİTİF SPEKTRUM İKİ KANALI BESLİYOR MU?
=======================================================================
(18 Ağustos 2026 — Sıçrama 1, rafine)

66'nın dersleri: (a) rezonans ailesi yanlış — aşırı-sönümlü (Debye) rejim;
(b) faz muhasebesi: boşluk = sayımın türevi → w EŞ-FAZLI, U := v/τ
ÇEYREK-FAZLI → doğal Re/Im çifti (w, U). Bonus: U(0)/w(0) ≈ 2.14 ≈ β.

Model (nedensellik = pozitiflik):
  χ(τ) = ε + Σ_j ρ_j · 1/(1 − i τ/s_j),  ρ_j ≥ 0
  Re χ = ε + Σ ρ_j/(1+x_j²),  Im χ = Σ ρ_j x_j/(1+x_j²),  x_j = τ/s_j

Test:  w ≈ Re χ  ve  U ≈ μ·Im χ  — AYNI ρ ile.
  - ayrı fitler (w kendi ρ_w'su, U kendi ρ_U'su) vs ortak ρ
  - "bağ maliyeti": ortak χ² − ayrı χ²'ler toplamı. Küçükse → tek analitik
    yanıt fonksiyonu İKİ KANALI BİRDEN taşıyor (KK-tipi sonuç, gevşeme formunda)
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.optimize import lsq_linear

HERE = Path(__file__).resolve().parent
_src = open(HERE / "66_kramers_kronig.py").read()
exec(_src[:_src.index("def chi(")])  # 66'nın veri bloklarını yükle (fit kısmı hariç)

tw, w, sw = W_DATA.T
tv, v, sv = V_DATA.T
U = v / tv
sU = sv / tv

S = np.geomspace(0.02, 2.0, 25)  # gevşeme ölçekleri gridi

def design(tau, kind):
    x = tau[:, None] / S[None, :]
    if kind == "re":
        return 1.0 / (1.0 + x**2)
    return x / (1.0 + x**2)

RIDGE = 1e-3

def solve_nn(A, y):
    n = A.shape[1]
    Ar = np.vstack([A, RIDGE * np.eye(n)])
    yr = np.concatenate([y, np.zeros(n)])
    r = lsq_linear(Ar, yr, bounds=(0, np.inf), method="bvls", max_iter=300)
    return r.x

def fit_sep():
    Aw = np.hstack([design(tw, "re"), np.ones((len(tw), 1)), -np.ones((len(tw), 1))])
    rw = solve_nn(Aw / sw[:, None], w / sw)
    res_w = Aw @ rw - w
    Au = np.hstack([design(tv, "im"), np.ones((len(tv), 1)), -np.ones((len(tv), 1))])
    ru = solve_nn(Au / sU[:, None], U / sU)
    res_u = Au @ ru - U
    return (np.sum((res_w / sw) ** 2), np.sqrt(np.mean(res_w**2)),
            np.sum((res_u / sU) ** 2), np.sqrt(np.mean(res_u**2)))

def fit_joint(mu):
    top = np.hstack([design(tw, "re"), np.ones((len(tw), 1)), -np.ones((len(tw), 1))])
    bot = np.hstack([mu * design(tv, "im"), np.zeros((len(tv), 2))])
    A = np.vstack([top / sw[:, None], bot / sU[:, None]])
    y = np.concatenate([w / sw, U / sU])
    r = solve_nn(A, y)
    resid = A @ r - y
    return np.sum(resid**2), r

chi_w, rms_w, chi_u, rms_u = fit_sep()
print(f"AYRI fitler: w-RMS = {rms_w:.4f} (χ²={chi_w:.0f}), "
      f"U-RMS = {rms_u:.4f} (χ²={chi_u:.0f})")
print("  (kıyas tabanı: τ-çökmesinin per-asal yarılma tabanı ~0.010-0.015)")

best = None
for mu in np.linspace(0.5, 4.0, 71):
    c, r = fit_joint(mu)
    if best is None or c < best[0]:
        best = (c, mu, r)
chi_j, mu_b, rho = best
top = np.hstack([design(tw, "re"), np.ones((len(tw), 1)), -np.ones((len(tw), 1))])
res_wj = top @ rho - w
res_uj = mu_b * (design(tv, "im") @ rho[:25]) - U
print(f"ORTAK fit  : w-RMS = {np.sqrt(np.mean(res_wj**2)):.4f}, "
      f"U-RMS = {np.sqrt(np.mean(res_uj**2)):.4f}  (μ = {mu_b:.2f}, χ²={chi_j:.0f})")
print(f"BAĞ MALİYETİ: Δχ² = {chi_j - (chi_w+chi_u):.0f}; "
      f"RMS artışları: w {np.sqrt(np.mean(res_wj**2))-rms_w:+.4f}, "
      f"U {np.sqrt(np.mean(res_uj**2))-rms_u:+.4f}")
print(f"\nε (ortak fit) = {rho[25] - rho[26]:+.3f}")
print("Aktif gevşeme ölçekleri (ρ_j > 0.01):")
for s, rj in zip(S, rho[:25]):
    if rj > 0.01:
        print(f"  s = {s:.3f}: ρ = {rj:.3f}")

# grafik
tt = np.linspace(1e-3, 0.5, 400)
Xre = np.hstack([design(tt, "re"), np.ones((len(tt), 1)), -np.ones((len(tt), 1))])
Xim = design(tt, "im")
w_fit = Xre @ rho
U_fit = mu_b * (Xim @ rho[:25])
fig, axes = plt.subplots(1, 2, figsize=(13, 5.2))
ax = axes[0]
ax.errorbar(tw, w, yerr=sw, fmt="o", ms=3.5, c="firebrick", label="w (44 nokta)")
ax.plot(tt, w_fit, "k-", lw=1.3, label="ortak-ρ Debye fiti")
ax.set_xlabel("τ"); ax.set_ylabel("w = Re χ"); ax.legend(fontsize=9); ax.grid(alpha=0.3)
ax.set_title("Eş-fazlı kanal")
ax = axes[1]
ax.errorbar(tv, U, yerr=sU, fmt="s", ms=3.5, c="teal", label="U = v/τ (90 nokta)")
ax.plot(tt, U_fit, "k-", lw=1.3, label=f"AYNI ρ → μ·Im χ (μ={mu_b:.2f})")
ax.set_xlabel("τ"); ax.set_ylabel("U = v/τ"); ax.legend(fontsize=9); ax.grid(alpha=0.3)
ax.set_title("Çeyrek-fazlı kanal — aynı spektrumdan")
plt.tight_layout()
out = HERE / "66b_debye.png"
plt.savefig(out, dpi=110)
print(f"\nGrafik: {out.name}")
