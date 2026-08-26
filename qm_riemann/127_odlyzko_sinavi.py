"""
127 — ODLYZKO SINAVI: MEKANİZMA YASASI BÜYÜK T'DE (26 Ağustos gecesi)
==========================================================================
126'nın türetimi (anti-nefes = katlanmış kuyruğun çarpımsal çift-
girişimi) geleneğin en eski kuralından geçmeli: BÜYÜK T'DE TEST.
Yasa: R_nn = −2cos(πτ_p), τ_p = log p / L — L büyüdükçe her asalın
τ'su küçülür, R_nn → −2'ye yaklaşmalı; biçim DEĞİŞMEMELİ.

VERİ MERDİVENİ (L: 9.86 → 44.6; yükseklikte ~10¹⁰ kat):
  A) 41'in ALTI penceresi (L 9.86/10.37/10.93/11.47/11.98/12.45)
  B) 53 Odlyzko penceresi (L=24.475, n=9999; çapalı ofsetler —
     R faz-izdüşümlü olduğundan mutlak-t kayması ETKİSİZ)
  C) odlyzko_zeros3/4/5.txt (10^k'ıncı sıfır blokları; başlıktan taban
     okunur; n=10⁴'er)
Taban notu (dürüst): L≥24'te tam-taban patlar; taban τ≤~0.30 (B) /
q≤720 (C) ile kesilir — dışarıda kalan çizgiler zaten mekanizmanın
kuyruğudur, η'da kalırlar; σ_η büyür, R-tanımı korunur (kayıtlı).

ÖN-MÜHÜR:
  O1  R_nn/(−2cosπτ) = 1 ± 0.15, TÜM yüksekliklerde (out-of-sample).
  O2  R_p ≈ −0.75..−0.95 civarı kalır (kinematik-bulaşıklı kanal).
  O3  En derin blokta (L≈44.6) tüm küçük asallar τ→0 → R_nn ≈ −2.0
      öbeklenmeli — yasanın en çıplak görünümü.
"""

import numpy as np
import re
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

def olc(zz, L, qcap, primler):
    g = np.diff(zz)
    m = 0.5 * (zz[:-1] + zz[1:])
    ds = g * L / TWO_PI - 1
    qs = pk_list(qcap)
    freqs = [np.log(q) for q in qs]
    b1, f1 = chunked_fit(ds, m, freqs)
    eta = ds - f1
    s2 = float((eta**2).mean())
    b2, _ = chunked_fit(eta**2 - s2, m, freqs)
    ee = eta[:-1] * eta[1:]
    c1 = float(ee.mean())
    mm = 0.5 * (m[:-1] + m[1:])
    b3, _ = chunked_fit(ee - c1, mm, freqs)
    out = []
    for p in primler:
        f = np.log(p)
        tau = f / L
        cg, sg = coef(b1, freqs, f)
        A1 = np.hypot(cg, sg); ph = np.arctan2(sg, cg)
        Rp = ((coef(b2, freqs, f)[0] * np.cos(ph) +
               coef(b2, freqs, f)[1] * np.sin(ph)) / (2 * s2 * A1))
        Rn = ((coef(b3, freqs, f)[0] * np.cos(ph) +
               coef(b3, freqs, f)[1] * np.sin(ph)) / (2 * c1 * A1))
        out.append((p, tau, Rp, Rn, -2 * np.cos(np.pi * tau)))
    return out, np.sqrt(s2), c1 / s2

def bas(ad, L, rows, se, c1r):
    print(f"\n[{ad}] L={L:.2f}  σ_η={se:.3f}  c₁/σ²={c1r:+.2f}")
    print(f"{'p':>4} {'τ':>6} {'R_p':>7} {'R_nn':>7} {'−2cosπτ':>8} {'oran':>6}")
    for p, tau, Rp, Rn, hed in rows:
        print(f"{p:>4} {tau:>6.3f} {Rp:>+7.3f} {Rn:>+7.3f} {hed:>+8.3f} "
              f"{Rn/hed:>6.2f}", flush=True)

# A) 41 pencereleri
d41 = np.load(HERE / "41_bigT_windows.npz")
K41 = sorted({x.split("_")[1] for x in d41.files}, key=lambda s: int(s[:-1]))
for k in K41:
    gz, tm = d41[f"gaps_{k}"], d41[f"tmid_{k}"]
    zz = np.empty(len(gz) + 1)
    zz[0] = tm[0] - gz[0] / 2
    zz[1:] = zz[0] + np.cumsum(gz)
    L = float(np.log(tm / TWO_PI).mean())
    pr = [p for p in [2, 3, 5, 7, 11, 13, 17, 23, 31] if np.log(p) / L < 0.31]
    rows, se, c1r = olc(zz, L, min(int(np.exp(0.52 * L)), 720), pr)
    bas(f"41-{k}", L, rows, se, c1r)

# B) 53 Odlyzko (L=24.5; çapalı ofsetler — faz-izdüşüm bağışıklığı)
d53 = np.load(HERE / "53_odlyzko_amps.npz")
g53, t53 = d53["gaps"], d53["t_mid"]
L53 = float(d53["L"])
zz = np.empty(len(g53) + 1)
zz[0] = t53[0] - g53[0] / 2
zz[1:] = zz[0] + np.cumsum(g53)
rows, se, c1r = olc(zz, L53, 1550, [2, 3, 5, 7, 11, 13])
bas("53-Odlyzko", L53, rows, se, c1r)

# C) Odlyzko metin blokları
for f in ["odlyzko_zeros3.txt", "odlyzko_zeros4.txt", "odlyzko_zeros5.txt"]:
    yol = HERE / f
    if not yol.exists():
        continue
    metin = yol.read_text().split("\n")
    taban = None
    vals = []
    for satir in metin:
        mm_ = re.search(r"gamma\s*-\s*([0-9]+)", satir)
        if mm_ and taban is None:
            taban = float(mm_.group(1))
        s = satir.strip()
        if re.fullmatch(r"[0-9]+\.[0-9]+", s):
            vals.append(float(s))
    if taban is None or len(vals) < 1000:
        print(f"\n[{f}] ATLANDI (taban={taban}, n={len(vals)})")
        continue
    zz = np.sort(np.array(vals))
    L = float(np.log((taban + zz.mean()) / TWO_PI))
    rows, se, c1r = olc(zz, L, 720, [2, 3, 5, 7, 11, 13])
    bas(f, L, rows, se, c1r)
