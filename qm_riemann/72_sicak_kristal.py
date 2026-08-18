"""
72 — SICAK KRİSTAL: DEBYE-WALLER TESTİ (18 Ağustos 2026)
==========================================================

Bragg resmi: sıfır örgüsü kristal, asal dalgası k = log p ile vuran ışın.
Isıl titreşim (sıfır jitter'ı u) koherent yanıtı e^{−k²σ_u²/2} söndürür.
σ_u ∝ 1/L (unfold sabit) olduğundan DW çarpanı SAF τ-yasası:
  DW(τ) = exp(−c·τ²),  c = (L·σ_u)²/2  ← JİTTERDAN HESAPLANIR, FİT DEĞİL

Model:  w(τ) = A·(1−2τ)·DW(τ)·[τ<1/2]  +  B
  (1−2τ): ideal Bragg yapısı (kesiş τ=1/2 = RS ufku = "2 salıncak")
  B: sıfır-aracılı dolaylı kanal (ufuk-ötesi platodan ölçülür)
  → serbest parametre: yalnız A (+c'nin serbest versiyonu çapraz-kilit için)

Testler:
  T1 jitter ölçümü 12 pencerede → unfold-sabitlik → c_jitter
  T2 model fiti: c SABİT (jitterdan) → RMS; kesiş konumu ≈ 0.40 çıkıyor mu?
  T3 c SERBEST → c_fit ≈ c_jitter mi? (sıcaklıkların buluşması)
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

def jitter_rms(gaps, tmid):
    t0 = tmid[0] - gaps[0] / 2
    k = np.arange(len(tmid)) + 0.5
    t = t0 + k * TWO_PI / np.log(t0 / TWO_PI)
    for _ in range(6):
        f = rvm_N(t) - rvm_N(t0) - k
        t = t - f / (np.log(t / TWO_PI) / TWO_PI)
    return float((tmid - t).std())

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
    return [(np.log(p) / L, b[3 + 2*i] / p**-0.5, se[3 + 2*i] / p**-0.5)
            for i, p in enumerate(PRIMES)]

# ---- pencereler ----
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

# ---- T1: jitter → c ----
print("T1 — JİTTER (kristal sıcaklığı):")
cs = []
for gaps, amps, tmid in WNDS:
    L = float(np.log(tmid / TWO_PI).mean())
    s = jitter_rms(gaps, tmid)
    su = s * L / TWO_PI   # unfold
    c = (L * s) ** 2 / 2
    cs.append(c)
    print(f"  L={L:5.2f}: σ_u = {s:.4f}  (unfold {su:.3f})  c = {c:.3f}")
c_jit = float(np.mean(cs))
print(f"  → c_jitter ortalaması = {c_jit:.3f}  (unfold jitter ~sabit mi? yukarı bak)\n")

# ---- veri: 132 w-noktası ----
pts = []
for gaps, amps, tmid in WNDS:
    pts += w_channel(gaps, amps, tmid)
pts = np.array(sorted(pts))
tau, w, sw = pts.T
B_PLATO = float(np.mean(w[tau > 0.5]))
print(f"B (ufuk-ötesi plato, ölçülen) = {B_PLATO:+.4f}  "
      f"({int((tau>0.5).sum())} nokta)\n")

def model(par, c_val):
    A_, = par[:1]
    coh = A_ * (1 - 2 * tau) * np.exp(-c_val * tau**2) * (tau < 0.5)
    return coh + B_PLATO

# ---- T2: c sabit (jitterdan), yalnız A serbest ----
f2 = least_squares(lambda p: (model(p, c_jit) - w) / sw, [1.1])
A2 = f2.x[0]
res2 = model(f2.x, c_jit) - w
rms2 = np.sqrt(np.mean(res2**2))
xr = np.linspace(0, 0.499, 2000)
cross2 = xr[np.argmin(np.abs(A2 * (1 - 2*xr) * np.exp(-c_jit * xr**2) + B_PLATO))]
print(f"T2 — c SABİT ({c_jit:.2f}): A = {A2:.3f}")
print(f"  RMS = {rms2:.4f}  (per-asal taban ~0.012-0.016)")
print(f"  modelin sıfır-geçişi: τ = {cross2:.3f}  (ölçülen τ* ≈ 0.395-0.40)\n")

# ---- T3: c serbest → sıcaklıklar buluşuyor mu? ----
f3 = least_squares(lambda p: ((p[0]*(1-2*tau)*np.exp(-p[1]*tau**2)*(tau<0.5)
                               + B_PLATO) - w) / sw, [1.1, c_jit],
                   bounds=([0.3, 0.0], [3.0, 10.0]))
A3, c3 = f3.x
res3 = (A3*(1-2*tau)*np.exp(-c3*tau**2)*(tau<0.5) + B_PLATO) - w
# bootstrap c hatası
rng = np.random.default_rng(72)
c_bs = []
for _ in range(200):
    ii = rng.integers(0, len(tau), len(tau))
    try:
        r = least_squares(lambda p: ((p[0]*(1-2*tau[ii])*np.exp(-p[1]*tau[ii]**2)
                                      *(tau[ii]<0.5) + B_PLATO) - w[ii]) / sw[ii],
                          [A3, c3], bounds=([0.3, 0.0], [3.0, 10.0]))
        c_bs.append(r.x[1])
    except Exception:
        pass
print(f"T3 — c SERBEST: A = {A3:.3f},  c_fit = {c3:.3f} ± {np.std(c_bs):.3f}")
print(f"  c_jitter = {c_jit:.3f}  →  fark {abs(c3-c_jit)/np.std(c_bs):.1f}σ")
print(f"  RMS = {np.sqrt(np.mean(res3**2)):.4f}")

# ---- grafik ----
fig, ax = plt.subplots(figsize=(9, 5.6))
ax.errorbar(tau, w, yerr=sw, fmt="o", ms=3.5, c="firebrick", alpha=0.65,
            label="w ölçümleri (132)")
tt = np.linspace(0.001, 0.62, 500)
mm = A2 * (1 - 2*tt) * np.exp(-c_jit * tt**2) * (tt < 0.5) + B_PLATO
ax.plot(tt, mm, "k-", lw=1.5,
        label=f"sıcak kristal: {A2:.2f}·(1−2τ)·e^(−{c_jit:.2f}τ²) {B_PLATO:+.2f}")
ax.axvline(0.5, color="gray", ls="--", lw=1, label="ideal Bragg / RS ufku (τ=1/2)")
ax.axvline(cross2, color="teal", ls=":", lw=1.2, label=f"modelin kesişi τ={cross2:.3f}")
ax.axhline(0, color="gray", lw=0.6)
ax.set_xlabel("τ = log p / L"); ax.set_ylabel("w")
ax.set_title("Sıcak kristal (Debye-Waller) modeli — c jitterdan, tek serbest A")
ax.legend(fontsize=9); ax.grid(alpha=0.3)
plt.tight_layout()
out = HERE / "72_sicak_kristal.png"
plt.savefig(out, dpi=110)
print(f"\nGrafik: {out.name}")
