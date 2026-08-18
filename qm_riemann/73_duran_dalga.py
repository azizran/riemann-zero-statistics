"""
73 — s≈3.3'ÜN TEORİSİ: DURAN DALGA SÖNÜMÜ (18 Ağustos 2026)
=============================================================

Hipotez: Bragg yansımasında kristal içi alan 2k periyotlu DURAN DALGADIR;
koherent sönüm |⟨e^{i·2k·u}⟩| ile olur (naif |⟨e^{iku}⟩| değil).
Gauss + bağımsız: s=4; gerçek dağılım/korelasyonla s = ölçülür.

Test — hiçbir Gauss varsayımı yok, AMPİRİK karakteristik fonksiyon:
  φ₁(k) = |⟨e^{i k u}⟩|   (naif DW)
  φ₂(k) = |⟨e^{i 2k u}⟩|  (duran dalga)
  u: pencere yerdeğiştirme serisi (tmid − RvM-akışı)

Model yarışı (yalnız A serbest, B=plato sabit):
  w(τ) = A·(1−2τ)·Φ(k=τL) + B,  Φ ∈ {φ₁, φ₂}  → hangisi taban RMS'e iner?
Ek ölçümler: u kurtosisi, komşu korelasyon ρ₁, s_eff = ln φ₂ / ln φ₁ oranı.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
A_CG = (np.e**2 - 5) / 2
B0, B1 = 2.7580, -0.0543
PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]

def rvm_N(t):
    x = t / TWO_PI
    return x * np.log(x / np.e) + 7 / 8

def u_series(gaps, tmid):
    t0 = tmid[0] - gaps[0] / 2
    k = np.arange(len(tmid)) + 0.5
    t = t0 + k * TWO_PI / np.log(t0 / TWO_PI)
    for _ in range(6):
        f = rvm_N(t) - rvm_N(t0) - k
        t = t - f / (np.log(t / TWO_PI) / TWO_PI)
    return tmid - t

def w_channel(gaps, amps, tmid):
    Lw = np.log(tmid / TWO_PI)
    g_u = gaps * Lw / TWO_PI
    a_u = amps / np.sqrt(A_CG * Lw + B0 + B1 / Lw)
    a_u /= np.sqrt((a_u**2).mean())
    ya = np.log(a_u)
    cw = [np.ones_like(ya), g_u, g_u**2]
    for q in PRIMES:
        arg = tmid * np.log(q)
        cw += [np.cos(arg), np.sin(arg)]
    X = np.vstack(cw).T
    b, *_ = np.linalg.lstsq(X, ya, rcond=None)
    se = np.sqrt((ya - X @ b).var() * np.diag(np.linalg.inv(X.T @ X)))
    L = float(Lw.mean())
    return [(np.log(p) / L, b[3+2*i] / p**-0.5, se[3+2*i] / p**-0.5, np.log(p))
            for i, p in enumerate(PRIMES)]

WNDS = []
d36 = np.load(HERE / "36_T100k.npz")
edges = np.geomspace(d36["t_mid"][0], d36["t_mid"][-1] * 1.0001, 13)
for i in range(12):
    m = (d36["t_mid"] >= edges[i]) & (d36["t_mid"] < edges[i + 1])
    if m.sum() >= 500:
        WNDS.append((d36["intervals"][m], d36["max_amps"][m], d36["t_mid"][m]))
d41 = np.load(HERE / "41_bigT_windows.npz")
for k in sorted({x.split("_")[1] for x in d41.files}, key=lambda s: int(s[:-1])):
    WNDS.append((d41[f"gaps_{k}"], d41[f"amps_{k}"], d41[f"tmid_{k}"]))

# ---- u istatistikleri + ampirik karakteristik fonksiyonlar ----
print("U-SERİSİ İSTATİSTİKLERİ:")
tau_l, w_l, sw_l, p1_l, p2_l = [], [], [], [], []
s_effs = []
for gaps, amps, tmid in WNDS:
    L = float(np.log(tmid / TWO_PI).mean())
    u = u_series(gaps, tmid)
    u = u - u.mean()
    kurt = float(np.mean(u**4) / np.mean(u**2)**2 - 3)
    rho1 = float(np.corrcoef(u[:-1], u[1:])[0, 1])
    row = w_channel(gaps, amps, tmid)
    for tau, w_, s_, kk in row:
        p1 = float(np.abs(np.mean(np.exp(1j * kk * u))))
        p2 = float(np.abs(np.mean(np.exp(2j * kk * u))))
        tau_l.append(tau); w_l.append(w_); sw_l.append(s_)
        p1_l.append(p1); p2_l.append(p2)
        if 0.25 < tau < 0.45 and p1 > 1e-3 and p2 > 1e-3:
            s_effs.append(np.log(p2) / np.log(p1))
    print(f"  L={L:5.2f}: σ_u={u.std():.4f}  kurtosis={kurt:+.3f}  ρ₁(komşu)={rho1:+.3f}")

tau = np.array(tau_l); w = np.array(w_l); sw = np.array(sw_l)
phi1 = np.array(p1_l); phi2 = np.array(p2_l)
B_PLATO = float(np.mean(w[tau > 0.5]))
print(f"\ns_eff = ln φ₂ / ln φ₁ (τ∈0.25-0.45 bölgesi): "
      f"{np.mean(s_effs):.2f} ± {np.std(s_effs):.2f}  (duran-dalga Gauss'u 4 der; 72'nin fiti 3.3 istedi)")

# ---- model yarışı ----
def fit_with(Phi, name):
    f = least_squares(lambda p: ((p[0]*(1-2*tau)*Phi*(tau<0.5) + B_PLATO) - w) / sw, [1.05])
    res = (f.x[0]*(1-2*tau)*Phi*(tau<0.5) + B_PLATO) - w
    rms = np.sqrt(np.mean(res**2))
    xr = np.linspace(0.0, 0.499, 2000)
    Phi_i = np.interp(xr, np.sort(tau), Phi[np.argsort(tau)])
    cross = xr[np.argmin(np.abs(f.x[0]*(1-2*xr)*Phi_i + B_PLATO))]
    print(f"  {name:<22} A={f.x[0]:.3f}  RMS={rms:.4f}  kesiş τ={cross:.3f}")
    return rms

print("\nMODEL YARIŞI (yalnız A serbest, B=plato):")
fit_with(np.ones_like(tau), "sönümsüz (kontrol)")
fit_with(phi1, "naif DW  φ₁(k)")
fit_with(phi2, "DURAN DALGA φ₂(2k)")

# grafik
fig, ax = plt.subplots(figsize=(9, 5.6))
o = np.argsort(tau)
ax.errorbar(tau, w, yerr=sw, fmt="o", ms=3, c="firebrick", alpha=0.55, label="w (132)")
f2 = least_squares(lambda p: ((p[0]*(1-2*tau)*phi2*(tau<0.5) + B_PLATO) - w) / sw, [1.05])
ax.plot(tau[o], (f2.x[0]*(1-2*tau)*phi2*(tau<0.5) + B_PLATO)[o], "k-", lw=1.4,
        label=f"duran-dalga modeli (A={f2.x[0]:.2f})")
f1 = least_squares(lambda p: ((p[0]*(1-2*tau)*phi1*(tau<0.5) + B_PLATO) - w) / sw, [1.05])
ax.plot(tau[o], (f1.x[0]*(1-2*tau)*phi1*(tau<0.5) + B_PLATO)[o], "--", c="gray", lw=1,
        label="naif DW")
ax.axhline(0, color="gray", lw=0.6); ax.axvline(0.5, color="gray", ls=":", lw=1)
ax.set_xlabel("τ"); ax.set_ylabel("w")
ax.set_title("s'nin teorisi: duran dalga (2k) sönümü — ampirik φ ile")
ax.legend(fontsize=9); ax.grid(alpha=0.3)
plt.tight_layout()
out = HERE / "73_duran_dalga.png"
plt.savefig(out, dpi=110)
print(f"\nGrafik: {out.name}")
