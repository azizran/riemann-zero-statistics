"""
129 — ZEROS6 SINAVI: YASA 10⁶ İSTATİSTİKLE + ÇÖZÜNÜRLÜK EĞRİSİ (27 Ağu)
==========================================================================
128'in kazanımı (zeros6: ilk 2.001.052 sıfır, L→12.10) ile 127'nin
çözülür-rejim testi büyük istatistikle; ve pencere-küçültme merdiveniyle
alet sınırının İLK ÖLÇÜMÜ.

ÖN-MÜHÜRLER:
  Y1  Büyük pencerelerde (n ≥ 10⁵) R_nn/(−2cosπτ) → 1 ± 0.05
      (127'de ±0.10 istatistik-sınırlıydı).
  Y2  Sabit L≈12'de n-merdiveni (10⁴ → 3·10⁴ → 10⁵ → 3·10⁵): derin-blok
      bozulması (σ_η şişmesi, c₁ kayması, R işaret bozulması) BURADA
      GÖRÜLMEMELİ — 127 tanısının testi: bozulma L'nin (orman
      yoğunluğu) işiydi, n'in değil. n=10⁴'te yalnız İSTATİSTİK
      gürültüsü büyümeli, tanı bayrakları temiz kalmalı.
  Y3  R_p düzlüğü dar hatalarla (içsel-β anlatısının büyük-n hali).
Pencereler: zeros6'nın son 1M (L 11.47-12.10), son 300k (11.95-12.10),
son 100k (12.06-12.10); + orta 300k (L≈10.9-11.3 civarı).

SONUÇ (27 Ağustos sabahı) — Y1 ✓ Y2 ✓ Y3 ✓:
  Y1: 1M-pencerede ilk yedi asal (τ≤0.24) oran 0.97-1.04 — YASA
  ±%3-4 DÜZEYİNDE DOĞRULANDI (üç pencere tutarlı). YENİ İNCE YAPI:
  τ≳0.26'da sistematik %8-12 eksik (oran 0.88-0.93; 126-D3 lab'ının
  %5-12'siyle aynı desen) — yasanın bir-üst-mertebe düzeltmesi;
  kalem-cilasının hedef eğrisi artık hassas: −2cosπτ·(1+δ(τ)),
  δ(0.26-0.30) ≈ −0.08..−0.12.
  Y2: TANI ÖLÇÜLDÜ — sabit L≈12'de n 3·10⁵→10⁴: σ_η 0.151-0.152,
  c₁/σ² −0.51 SABİT, oranlar ~1 (yalnız istatistik saçılması büyür,
  n=10⁴'te ±0.09). Derin-blok patolojisi YOK → 127'nin yüksek-L
  bozulması L'nin (orman yoğunluğu) işiydi, n'in değil — tanı artık
  ölçülmüş gerçek.
  Y3: R_p = −0.67..−0.90, dar hatalarla; bilinen kinematik-bulaşıklı
  yumuşak τ-eğilimi.
"""

import numpy as np
from pathlib import Path
from sympy import primerange

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi

def pk_list(lim):
    out = []
    for p in primerange(2, int(lim) + 1):
        q = p
        while q <= lim:
            out.append(q); q *= p
    return sorted(set(out))

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
    return b, fit

def coef(b, freqs, f):
    i = freqs.index(f)
    return b[3 + 2 * i], b[3 + 2 * i + 1]

def olc(zz, primler):
    g = np.diff(zz)
    m = 0.5 * (zz[:-1] + zz[1:])
    Lw = np.log(m / TWO_PI)
    L = float(Lw.mean())
    ds = g * Lw / TWO_PI - 1
    qs = pk_list(min(int(np.exp(0.52 * L)), 720))
    freqs = [np.log(q) for q in qs]
    b1, f1 = chunked_fit(ds, m, freqs)
    eta = ds - f1
    s2 = float((eta**2).mean())
    b2, _ = chunked_fit(eta**2 - s2, m, freqs)
    ee = eta[:-1] * eta[1:]
    c1 = float(ee.mean())
    mm = 0.5 * (m[:-1] + m[1:])
    b3, _ = chunked_fit(ee - c1, mm, freqs)
    rows = []
    for p in primler:
        f = np.log(p)
        tau = f / L
        cg, sg = coef(b1, freqs, f)
        A1 = np.hypot(cg, sg); ph = np.arctan2(sg, cg)
        Rp = ((coef(b2, freqs, f)[0] * np.cos(ph) +
               coef(b2, freqs, f)[1] * np.sin(ph)) / (2 * s2 * A1))
        Rn = ((coef(b3, freqs, f)[0] * np.cos(ph) +
               coef(b3, freqs, f)[1] * np.sin(ph)) / (2 * c1 * A1))
        rows.append((p, tau, Rp, Rn, -2 * np.cos(np.pi * tau)))
    return rows, np.sqrt(s2), c1 / s2, L

d = np.load(HERE / "128_odl_zeros6_2e6_zeros.npz")
Z = np.sort(np.asarray(d["zeros"], dtype=float))
N = len(Z)
PR = [2, 3, 5, 7, 11, 13, 17, 23, 31]

def bas(ad, rows, se, c1r, L):
    print(f"\n[{ad}] L={L:.2f}  σ_η={se:.3f}  c₁/σ²={c1r:+.2f}")
    print(f"{'p':>4} {'τ':>6} {'R_p':>7} {'R_nn':>7} {'−2cosπτ':>8} {'oran':>6}")
    for p, tau, Rp, Rn, hed in rows:
        print(f"{p:>4} {tau:>6.3f} {Rp:>+7.3f} {Rn:>+7.3f} {hed:>+8.3f} "
              f"{Rn/hed:>6.2f}", flush=True)

print("Y1/Y3 — büyük pencereler:")
for ad, sl in [("son-1M", slice(N - 1000000, N)),
               ("son-300k", slice(N - 300000, N)),
               ("orta-300k", slice(700000, 1000000))]:
    zz = Z[sl]
    pr = [p for p in PR if np.log(p) / np.log(zz.mean() / TWO_PI) < 0.31]
    rows, se, c1r, L = olc(zz, pr)
    bas(ad, rows, se, c1r, L)

print("\nY2 — n-merdiveni (hepsi zeros6'nın SONUNDAN; sabit L≈12):")
for n_w in [300000, 100000, 30000, 10000]:
    zz = Z[N - n_w:N]
    rows, se, c1r, L = olc(zz, [2, 3, 5, 7, 13])
    oranlar = [f"{r[3]/r[4]:.2f}" for r in rows]
    print(f"  n={n_w:>6}: L={L:.2f} σ_η={se:.3f} c₁/σ²={c1r:+.2f}  "
          f"R_nn/hedef oranları: {' '.join(oranlar)}", flush=True)
