"""
114 — KALEM OTURUMU HAKEMİ: −(3/4)σκ² DOĞRULAMASI + LAG-MOMENTLER (26 Ağu)
==========================================================================
TÜRETİLDİ (Isserlis, sabah): Gaussian kinematik çarpıklık yasası
   R₃(σ³) = −(3/4)·σ_ds·κ²
Mekanizma: dilatasyon eğriliğinin gürültü-karesini doğrultması; kritik
incelik: fit ⟨N⟩'yi merkezler → katsayı 12σ⁴(P−Q) değil 6σ⁴(P−Q).
Dünkü üç yüksek-κ hücresi %2-7 içinde vurdu (−0.60/−0.885/−1.20 vs
−0.584/−0.949/−1.148).

BU SCRIPT:
  T1  Katsayı doğrulaması: ince (σ,τ) ızgarası, yüksek-κ bölgesi;
      ölçüm/öngörü oranı tablosu (hedef: 1.00±0.1; düşük-τ hücreleri
      gürültü-tabanı nedeniyle dışarıda, dürüstçe işaretli).
  T2  GERÇEK η'NIN LAG-MOMENTLERİ (ζ, iki pencere) — hedef-2'nin
      (ζ-boyalı düz −2.8) korelasyonlu-Isserlis türetiminin girdileri:
      ⟨η_n²η_{n+1}⟩, ⟨η_nη_{n+1}²⟩, ⟨η_n²η_{n+2}⟩, ⟨η_nη_{n+1}η_{n+2}⟩
      (σ³ birimli) + karşılaştırma: Gaussian-iid'de hepsi 0.
"""

import numpy as np
from pathlib import Path
from sympy import primerange

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
rng = np.random.default_rng(114)
L0 = 10.37
NZ = 200000
gbar = TWO_PI / L0
A_DS = 0.10

def pk_list(lim):
    out = []
    for p in primerange(2, int(lim) + 1):
        q = p
        while q <= lim:
            out.append(q); q *= p
    return sorted(set(out))

def rvm_N(t):
    x = t / TWO_PI
    return x * np.log(x / np.e) + 7 / 8

t0 = TWO_PI * np.exp(L0)
idx = np.arange(NZ + 1, dtype=float)
xg = t0 + idx * gbar
for _ in range(8):
    xg = xg - (rvm_N(xg) - rvm_N(t0) - idx) / (np.log(xg / TWO_PI) / TWO_PI)

def r3_olc(z, om):
    g = np.diff(z)
    m = 0.5 * (z[:-1] + z[1:])
    ds = g / gbar - 1
    ds = ds - ds.mean()
    X = np.vstack([np.ones_like(m), np.cos(om * m), np.sin(om * m),
                   np.cos(2 * om * m), np.sin(2 * om * m)]).T
    b, *_ = np.linalg.lstsq(X, ds, rcond=None)
    A1 = np.hypot(b[1], b[2]); ph = np.arctan2(b[2], b[1])
    eta = ds - X @ b
    s2 = float((eta**2).mean()); m3 = float((eta**3).mean())
    b3, *_ = np.linalg.lstsq(X, eta**3 - m3, rcond=None)
    P3 = b3[1] * np.cos(ph) + b3[2] * np.sin(ph)
    return P3 / (s2**1.5 * A1)

print("T1 — katsayı doğrulaması (öngörü R₃ = −0.75·σ·κ²):")
print(f"{'σ_ds':>6} {'τ':>5} {'ölçüm':>8} {'öngörü':>8} {'oran':>6}")
for sig in [0.08, 0.12, 0.148, 0.20]:
    for tau in [0.30, 0.38, 0.45, 0.50]:
        om = tau * L0
        kap = 2 * np.pi * tau
        UA = A_DS / (2 * np.sin(np.pi * tau)) * gbar
        w = rng.normal(0, sig / np.sqrt(2) * gbar, NZ + 1)
        y = np.sort(xg + w)
        z = np.sort(y + UA * np.cos(om * y))
        r3 = r3_olc(z, om)
        ong = -0.75 * sig * kap**2
        print(f"{sig:>6.3f} {tau:>5.2f} {r3:>+8.4f} {ong:>+8.4f} "
              f"{r3/ong:>6.2f}", flush=True)

print("\nT2 — GERÇEK η lag-momentleri (σ³ birimli; Gaussian-iid: 0):")
d41 = np.load(HERE / "41_bigT_windows.npz")

def chunked_fit(y, tmid, freqs, chunk=40000):
    tt = (tmid - tmid.mean()) / (tmid[-1] - tmid[0])
    C = 3 + 2 * len(freqs)
    XtX = np.zeros((C, C)); Xty = np.zeros(C)
    n = len(y)
    def cols(sl):
        c = [np.ones(sl.stop - sl.start), tt[sl], tt[sl]**2]
        for f in freqs:
            arg = f * tmid[sl]
            c += [np.cos(arg), np.sin(arg)]
        return np.vstack(c).T
    for s0 in range(0, n, chunk):
        sl = slice(s0, min(s0 + chunk, n))
        Xc = cols(sl)
        XtX += Xc.T @ Xc; Xty += Xc.T @ y[sl]
    b = np.linalg.solve(XtX, Xty)
    fit = np.empty(n)
    for s0 in range(0, n, chunk):
        sl = slice(s0, min(s0 + chunk, n))
        fit[sl] = cols(sl) @ b
    return fit

for k in ["120k", "200k"]:
    gz, tm = d41[f"gaps_{k}"], d41[f"tmid_{k}"]
    zz = np.empty(len(gz) + 1)
    zz[0] = tm[0] - gz[0] / 2
    zz[1:] = zz[0] + np.cumsum(gz)
    g = np.diff(zz)
    m = 0.5 * (zz[:-1] + zz[1:])
    Lw = np.log(m / TWO_PI)
    L = float(Lw.mean())
    ds = g * Lw / TWO_PI - 1
    qs = pk_list(min(np.exp(0.52 * L), 720))
    eta = ds - chunked_fit(ds, m, [np.log(q) for q in qs])
    s3 = float((eta**2).mean())**1.5
    e0, e1, e2 = eta[:-2], eta[1:-1], eta[2:]
    print(f"[{k}] skew={float((eta**3).mean())/s3:+.3f}  "
          f"⟨η²η₊⟩/σ³={float((e0**2*e1).mean())/s3:+.4f}  "
          f"⟨ηη₊²⟩/σ³={float((e0*e1**2).mean())/s3:+.4f}  "
          f"⟨η²η₊₊⟩/σ³={float((e0**2*e2).mean())/s3:+.4f}  "
          f"⟨ηη₊η₊₊⟩/σ³={float((e0*e1*e2).mean())/s3:+.4f}", flush=True)
