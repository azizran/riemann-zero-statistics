"""
120 — Ö1: ADALARDA NEFES R_p — SICAKLIK AYRIŞTIRICISI (26 Ağustos)
==========================================================================
Kalem defteri Ö1 (mühürlü): yedi adada R_p ilk kez.
  H-sabitlenme: aşırı-termal sabitlenme çerçevesinde soğuk ölü-2
     üçlüsünün (σ² −%20) R'si canlılardan ayrışabilir.
  H-evrensel: D gibi (108b) R de adadan bağımsızdır → sabit dalgalanma
     bütçesi YALNIZ evrensel büyüklüklerden kurulu olmalı (türetimi
     daraltan güçlü ipucu).
Makine: 109'un zinciri ada başına — sertifikalı sıfırlar (101f/105b;
118 taper-dökümü gelince çapraz kontrol edilir), üst %75, tam-taban
fit (iletken yoğunluğu), η²-dalgası gap-dalgası yönüne izdüşüm,
sağ kalan asal çizgiler τ∈[0.06,0.30]. R_nn de basılır. ζ referansı:
R_p ≈ −0.85 düz (109).

SONUÇ (26 Ağustos) — H-EVRENSEL KAZANDI + YENİ KAPALI FORM:
  Ada ortalamaları ⟨R_p⟩ = −0.880..−0.941; ölü-2 üçlüsü (−0.890/
  −0.898/−0.914) canlılardan AYRIŞMIYOR (çizgi-içi saçılma ±0.05
  içinde); hepsi ζ'nın −0.85'iyle uyumlu. NEFES DE EVRENSEL (D gibi).
  → Sabit dalgalanma bütçesi YALNIZ evrensel büyüklüklerden kurulmalı;
  sıcaklık/çizgi-envanteri nefes katsayısına girmiyor.
  YENİ KAPALI FORM (tablodan): R_nn ≈ −2·cos(πτ) — yedi adada VE
  ζ'da (110c modelinin β+kin_nn yapısında bond kanalı kinematik-temiz;
  içsel β = −2.0 EVRENSEL SABİT olarak en temiz bond'dan okunuyor;
  düşük-orta τ'da ±%3, yüksek τ'da ~%10).
  YAN BULGU: σ_η = 0.149-0.157 ve c₁/σ² = −0.58..−0.60 YEDİ ADADA
  ÖZDEŞ — artık-gürültü istatistiği evrensel; 105'in −%20 sıcaklık
  farkı tamamen DETERMİNİSTİK DALGA payında yaşıyor (tutarlılık ✓).
  Türetim hedefi keskinleşti: "β = −2 neden evrensel?" (iki kanal/
  iki kuadratür bütçesinin −1+−1'i mi?)
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

SETS = [
    ("chi3",  3, "101f_chi3_zeros.npz",  {3}),
    ("beta",  4, "101f_beta_zeros.npz",  {2}),
    ("chi5",  5, "101f_chi5_zeros.npz",  {5}),
    ("chi7",  7, "101f_chi7_zeros.npz",  {7}),
    ("chi5e", 5, "105b_chi5e_zeros.npz", {5}),
    ("chi8e", 8, "105b_chi8e_zeros.npz", {2}),
    ("chi8o", 8, "105b_chi8o_zeros.npz", {2}),
]
OLU2 = {"beta", "chi8e", "chi8o"}
PRL = [2, 3, 5, 7, 11, 13, 17]

SON = {}
for ad, q, dosya, olu in SETS:
    zc = np.load(HERE / dosya)["zeros"]
    lo = np.exp(np.log(zc[0] + 1) + 0.25 * (np.log(zc[-1]) - np.log(zc[0] + 1)))
    zc = zc[zc >= lo]
    g = np.diff(zc)
    m = 0.5 * (zc[:-1] + zc[1:])
    Lw = np.log(q * m / TWO_PI)
    L = float(Lw.mean())
    ds = g * Lw / TWO_PI - 1
    qs = pk_list(min(np.exp(0.52 * L), 720))
    freqs = [np.log(qq) for qq in qs] + [np.log(p) + 0.037 for p in PRL]
    b1, f1 = chunked_fit(ds, m, freqs)
    eta = ds - f1
    s2 = float((eta**2).mean())
    b2, _ = chunked_fit(eta**2 - s2, m, freqs)
    ee = eta[:-1] * eta[1:]
    c1 = float(ee.mean())
    mm = 0.5 * (m[:-1] + m[1:])
    b3, _ = chunked_fit(ee - c1, mm, freqs)
    SON[ad] = {}
    for p in PRL:
        if p in olu:
            continue
        f = np.log(p)
        tau = f / L
        if not (0.055 < tau < 0.31):
            continue
        cg, sg = coef(b1, freqs, f)
        A1 = np.hypot(cg, sg); ph = np.arctan2(sg, cg)
        cv, sv = coef(b2, freqs, f)
        Rp = (cv * np.cos(ph) + sv * np.sin(ph)) / (2 * s2 * A1)
        cn, sn = coef(b3, freqs, f)
        Rnn = (cn * np.cos(ph) + sn * np.sin(ph)) / (2 * c1 * A1)
        fpl = np.log(p) + 0.037
        cp, sp = coef(b2, freqs, fpl)
        pl = np.hypot(cp, sp) / (2 * s2 * A1)
        SON[ad][p] = (tau, Rp, Rnn, pl)
    print(f"[{ad}] L={L:.2f} σ_η={np.sqrt(s2):.4f} c₁/σ²={c1/s2:+.2f}",
          flush=True)

print(f"\nR_p tablosu (satır=p; ζ referans −0.85 düz):")
print(f"{'p':>3}", end="")
for ad, *_ in SETS:
    print(f" {ad:>7}", end="")
print()
for p in PRL:
    print(f"{p:>3}", end="")
    for ad, *_ in SETS:
        v = SON[ad].get(p)
        print(f" {v[1]:>+7.3f}" if v else f" {'—':>7}", end="")
    print()
print(f"\nR_nn tablosu:")
for p in PRL:
    print(f"{p:>3}", end="")
    for ad, *_ in SETS:
        v = SON[ad].get(p)
        print(f" {v[2]:>+7.3f}" if v else f" {'—':>7}", end="")
    print()
print("\nAda ortalamaları (τ∈[0.06,0.31] asalları):")
for ad, *_ in SETS:
    vals = [v[1] for v in SON[ad].values()]
    pls = [v[3] for v in SON[ad].values()]
    et = "ÖLÜ-2" if ad in OLU2 else "canlı"
    print(f"  {ad:>6} [{et}]: ⟨R_p⟩ = {np.mean(vals):+.3f} "
          f"(n={len(vals)}; plasebo {np.mean(pls):.3f})")
