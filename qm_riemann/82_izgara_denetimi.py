"""
82 — ÖRNEKLEME-IZGARASI DENETİMİ: ÜÇLEMENİN KAPI TESTİ (20 Ağustos 2026)
==========================================================================
81'in mekanizması (ızgara asal dalgası taşır → atlanmış-değişken
transferi × √p büyütmesi) kat dibinde okumanın %88'ini hayalet yaptı.
SORU: üçlemenin manşet sayıları (küçük-orta τ) bağışık mı?

Öngörü: bağışık — büyütme √p, küçük τ'da √2..√31 ≤ 5.6; kat dibinde ~20.
Ama öngörü yetmez, ölçeriz:

  T1  KORUNUM YASASI: 77'nin manşet fiti (44 çift, 11 pencere) eski taban
      (yalnız P4) vs TAM taban (tüm p^k ≤ min(e^{0.45L}, 720)) —
      (a, b, RMS) ve a≡1 eğimi kaç σ oynuyor? (77 blok-bootstrap ±0.012)
  T2  τ-YASALARI: 132 nokta (P11 × 12 pencere) eski taban (P11) vs tam
      taban (p^k ≤ e^{0.55L}) — Δw/σ dağılımı, τ-binli profil.
  T3  ÖNYARGI BÜTÇESİ: hedef kolonun geri-kalan-taban R²'si × √p — τ
      profili (neden küçük τ güvenli, kat dibi değil).

Bellek: kaya-normal-denklem (parça parça XtX) — Odlyzko 1.5M satır dahil.

SONUÇ (20 Ağustos): KAPI KAPANDI — ÜÇLEME BEKLEMEYE ALINDI.
  T1: yasa fiti P4→TAM: a 0.990→1.017, b 1.025→0.884 (12σ!), RMS
      0.0096→0.0075 (İYİLEŞİYOR). Merdiven TAM/2→TAM'da yakınsıyor.
  T2: 132 noktanın HEPSİ oynuyor (ort 14.5σ, küçük τ'da 20σ); Δv/v %13-20.
  T3: bütçe küçük τ'da da büyük — R² yeterince düşmüyor (öngörüm yanlıştı).
  T4 (yer-gerçeği): gerçek ızgara + bilinen sinyal → P4 tabanı w'yi
      %15-25 DÜŞÜK ölçer (21σ); TAM taban ~%1 içinde doğru bulur.
      → Kısıtlı-taban kanal sayıları sistematik taşıyor; fiziksel okuma
      TAM-taban okuması. Not 2-3 sayıları YENİDEN ÖLÇÜLMELİ (Not 1 muaf:
      kanal regresyonu içermiyor). Yasa formu sağlam (RMS iyileşiyor),
      katsayıları revize: a≈1.02, b≈0.88.
"""

import numpy as np
import mpmath as mp
from pathlib import Path
from sympy import primerange

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
A_CG = (np.e**2 - 5) / 2
B0, B1 = 2.7580, -0.0543
T0_ODL = 267653395647
P4 = [2, 3, 5, 7]
P11 = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
mp.mp.dps = 30

def pk_list(lim, kmin=1):
    out = []
    for p in primerange(2, int(lim) + 1):
        q, k = p, 1
        while q <= lim:
            if k >= kmin:
                out.append(q)
            q *= p; k += 1
    return sorted(set(out))

def chunked_reg(y, tmid, qs, anchored, g_u=None, chunk=40000):
    C = (3 if g_u is not None else 1) + 2 * len(qs)
    XtX = np.zeros((C, C)); Xty = np.zeros(C)
    ph = {q: (float(mp.fmod(T0_ODL * mp.log(q), 2 * mp.pi)) if anchored else 0.0)
          for q in qs}
    n = len(y); ss_y = 0.0
    for s0 in range(0, n, chunk):
        sl = slice(s0, min(s0 + chunk, n))
        cols = [np.ones(sl.stop - sl.start)]
        if g_u is not None:
            cols += [g_u[sl], g_u[sl]**2]
        for q in qs:
            arg = ph[q] + tmid[sl] * np.log(q)
            cols += [np.cos(arg), np.sin(arg)]
        Xc = np.vstack(cols).T
        XtX += Xc.T @ Xc; Xty += Xc.T @ y[sl]
    b = np.linalg.solve(XtX, Xty)
    # artık varyansı için ikinci geçiş
    for s0 in range(0, n, chunk):
        sl = slice(s0, min(s0 + chunk, n))
        cols = [np.ones(sl.stop - sl.start)]
        if g_u is not None:
            cols += [g_u[sl], g_u[sl]**2]
        for q in qs:
            arg = ph[q] + tmid[sl] * np.log(q)
            cols += [np.cos(arg), np.sin(arg)]
        Xc = np.vstack(cols).T
        r = y[sl] - Xc @ b
        ss_y += (r**2).sum()
    sig2 = ss_y / n
    se = np.sqrt(sig2 * np.diag(np.linalg.inv(XtX)))
    return b, se

def unfold(gaps, amps, tmid, anchored=False):
    Lw = np.log(((T0_ODL + tmid) if anchored else tmid) / TWO_PI)
    g_u = gaps * Lw / TWO_PI
    a_u = amps / np.sqrt(A_CG * Lw + B0 + B1 / Lw)
    a_u /= np.sqrt((a_u**2).mean())
    return np.log(a_u), np.log(g_u), g_u, float(Lw.mean())

def kanal(gaps, amps, tmid, targets, qs, anchored):
    ya, yg, g_u, L = unfold(gaps, amps, tmid, anchored)
    bw, sew = chunked_reg(ya, tmid, qs, anchored, g_u=g_u)
    bv, sev = chunked_reg(yg, tmid, qs, anchored, g_u=None)
    out = []
    for p in targets:
        i = qs.index(p)
        iw, iv = 3 + 2 * i, 1 + 2 * i
        out.append((np.log(p) / L,
                    bw[iw] / p**-0.5, sew[iw] / p**-0.5,
                    float(np.hypot(bv[iv], bv[iv + 1]) / p**-0.5)))
    return out

# ---- pencereler
LAW = []
d41 = np.load(HERE / "41_bigT_windows.npz")
K41 = sorted({x.split("_")[1] for x in d41.files}, key=lambda s: int(s[:-1]))
for k in K41[:6]:
    LAW.append((d41[f"gaps_{k}"], d41[f"amps_{k}"], d41[f"tmid_{k}"], False))
for f in ["55_win_1e+08.npz", "55_win_1e+09.npz", "55_win_1e+10.npz", "55_win_1e+11.npz"]:
    d = np.load(HERE / f)
    LAW.append((d["gaps"], d["amps"], d["tmid"], False))
d53 = np.load(HERE / "53_odlyzko_amps.npz")
LAW.append((d53["gaps"], d53["max_amps"], d53["t_mid"], True))

EXT = []
d36 = np.load(HERE / "36_T100k.npz")
edges = np.geomspace(d36["t_mid"][0], d36["t_mid"][-1] * 1.0001, 13)
for i in range(12):
    m = (d36["t_mid"] >= edges[i]) & (d36["t_mid"] < edges[i + 1])
    if m.sum() >= 3000:
        EXT.append((d36["intervals"][m], d36["max_amps"][m], d36["t_mid"][m], False))
for k in K41:
    EXT.append((d41[f"gaps_{k}"], d41[f"amps_{k}"], d41[f"tmid_{k}"], False))

# ============ T1: KORUNUM YASASI, ESKİ vs TAM TABAN ============
print("T1 — KORUNUM YASASI (44 çift):")
for isim, taban_al in [("eski (yalnız P4)", lambda L: P4),
                       ("TAM (p^k ≤ min(e^{0.45L},720))",
                        lambda L: sorted(set(P4) | set(pk_list(min(np.exp(0.45*L), 720)))))]:
    pairs = []
    for gaps, amps, tmid, anch in LAW:
        L = float(np.log(((T0_ODL + tmid) if anch else tmid) / TWO_PI).mean())
        qs = taban_al(L)
        for tau, w_, s_, v_ in kanal(gaps, amps, tmid, P4, qs, anch):
            if w_ > 0:
                pairs.append((v_, np.sqrt(w_), w_))
    arr = np.array(pairs)
    vv, sq, wv = arr.T
    X = np.vstack([np.ones_like(vv), vv]).T
    c, *_ = np.linalg.lstsq(X, sq, rcond=None)
    rms = np.sqrt(np.mean((sq - X @ c)**2))
    b1 = np.sum((1 - sq) * vv) / np.sum(vv * vv)
    print(f"  {isim}: √w = {c[0]:.4f} − {-c[1]:.4f}·v  RMS {rms:.4f}  |  a≡1: b = {b1:.4f}"
          f"  (n={len(arr)})")

# ============ T2: τ-YASALARI, 132 NOKTA ============
print("\nT2 — τ-YASALARI (P11 × 12 pencere): Δw/σ dağılımı")
rows = []
for gaps, amps, tmid, anch in EXT:
    L = float(np.log(tmid / TWO_PI).mean())
    eski = kanal(gaps, amps, tmid, P11, list(P11), anch)
    tamq = sorted(set(P11) | set(pk_list(min(np.exp(0.55 * L), 1000))))
    yeni = kanal(gaps, amps, tmid, P11, tamq, anch)
    for (t, w0, s0, v0), (_, w1, s1, v1) in zip(eski, yeni):
        rows.append((t, w0, s0, w1, s1, v0, v1))
R = np.array(rows)
dz = (R[:, 3] - R[:, 1]) / R[:, 2]
dv = np.abs(R[:, 6] - R[:, 5]) / np.maximum(R[:, 5], 1e-9)
print(f"  Δw/σ: ort |Δ/σ| = {np.abs(dz).mean():.2f}   maks = {np.abs(dz).max():.2f}"
      f"   |Δ/σ|>1 oranı = {(np.abs(dz) > 1).mean()*100:.0f}%")
print(f"  Δv/v: ort = {dv.mean()*100:.1f}%   maks = {dv.max()*100:.1f}%")
print("  τ-binli |Δw/σ| profili:")
for lo, hi in [(0.0, 0.15), (0.15, 0.25), (0.25, 0.35), (0.35, 0.45), (0.45, 0.62)]:
    m = (R[:, 0] >= lo) & (R[:, 0] < hi)
    if m.sum():
        print(f"    τ∈[{lo:.2f},{hi:.2f}): n={int(m.sum()):>3}  ⟨|Δw/σ|⟩ = "
              f"{np.abs(dz[m]).mean():.2f}  maks {np.abs(dz[m]).max():.2f}")

# ============ T3: ÖNYARGI BÜTÇESİ ============
print("\nT3 — ÖNYARGI BÜTÇESİ (hedef kolon R² × √p; L=10.37 penceresi):")
gaps, amps, tmid, anch = EXT[7]
ya, yg, g_u, L = unfold(gaps, amps, tmid)
tamq = sorted(set(pk_list(min(np.exp(0.55 * L), 1000))))
hedefler = [2, 5, 13, 31] + [p for p in primerange(int(np.exp(0.5*L)),
                                                   int(np.exp(0.55*L)))][:2]
for p in hedefler:
    c = np.cos(tmid * np.log(p)); c = c - c.mean()
    digq = [q for q in tamq if q != p][:150]
    Xo = np.vstack([f(tmid * np.log(q)) for q in digq
                    for f in (np.cos, np.sin)]).T
    bb, *_ = np.linalg.lstsq(Xo, c, rcond=None)
    r2 = 1 - ((c - Xo @ bb)**2).sum() / (c**2).sum()
    print(f"  p={p:>4} (τ={np.log(p)/L:.3f}): R² = {r2:.2e}  →  ölçek √p·√R² = "
          f"{np.sqrt(p * max(r2, 0)):.4f}")

# ============ T4: MERDİVEN + SENTETİK YER-GERÇEĞİ ============
print("\nT4a — T1 merdiveni:")
for isim, taban_al in [
    ("P4", lambda L: list(P4)),
    ("P11", lambda L: list(P11)),
    ("TAM/2 (q<=min(e^.45L,360))", lambda L: sorted(set(P4) | set(pk_list(min(np.exp(0.45*L), 360))))),
    ("TAM (q<=min(e^.45L,720))", lambda L: sorted(set(P4) | set(pk_list(min(np.exp(0.45*L), 720))))),
]:
    pairs = []
    for gaps, amps, tmid, anch in LAW:
        L = float(np.log(((T0_ODL + tmid) if anch else tmid) / TWO_PI).mean())
        for tau, w_, s_, v_ in kanal(gaps, amps, tmid, P4, taban_al(L), anch):
            if w_ > 0:
                pairs.append((v_, np.sqrt(w_)))
    arr = np.array(pairs); vv, sq = arr.T
    X = np.vstack([np.ones_like(vv), vv]).T
    c, *_ = np.linalg.lstsq(X, sq, rcond=None)
    rms = np.sqrt(np.mean((sq - X @ c)**2))
    b1 = np.sum((1 - sq) * vv) / np.sum(vv * vv)
    print(f"  {isim:>28}: a={c[0]:.4f} b={-c[1]:.4f} RMS={rms:.4f} | a≡1 b={b1:.4f}")

print("\nT4b — sentetik yer-gerçeği (gerçek ızgara, bilinen sinyal):")
gaps = d41[f"gaps_{K41[1]}"]; amps = d41[f"amps_{K41[1]}"]; tmid = d41[f"tmid_{K41[1]}"]
Lw = np.log(tmid / TWO_PI); L = float(Lw.mean())
g_u = gaps * Lw / TWO_PI
rng = np.random.default_rng(82)
qs_all = pk_list(min(np.exp(0.55 * L), 1000))
def w_true(tau):
    return (1 - 2*tau) * np.exp(-1.2 * tau**2) if tau < 0.5 else -0.05
y = 0.35 * (g_u - g_u.mean()) - 0.06 * (g_u**2 - (g_u**2).mean())
for q in qs_all:
    y = y + w_true(np.log(q) / L) * q**-0.5 * np.cos(tmid * np.log(q))
y = y + rng.normal(0, 0.55, len(y))
b4, se4 = chunked_reg(y, tmid, list(P4), False, g_u=g_u)
bt, _ = chunked_reg(y, tmid, qs_all, False, g_u=g_u)
for p in P4:
    i4 = P4.index(p); it = qs_all.index(p)
    wt = w_true(np.log(p) / L)
    w4 = b4[3 + 2*i4] / p**-0.5; wT = bt[3 + 2*it] / p**-0.5
    s4 = se4[3 + 2*i4] / p**-0.5
    print(f"  p={p}: gerçek {wt:.4f} | P4-taban {w4:.4f} ({abs(w4-wt)/s4:.0f}σ sapma) | TAM {wT:.4f}")
