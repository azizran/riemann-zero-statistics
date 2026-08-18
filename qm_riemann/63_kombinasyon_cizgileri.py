"""
63 — KOMBİNASYON ÇİZGİLERİ: İKİNCİ-MERTEBE ASAL İÇERİĞİ (18 Ağustos 2026)
===========================================================================

Denetim S2 bulgusu: genlik serisinde log(p·q) ve log(p/q) frekanslarında
11-13σ, FAZ-KİLİTLİ (cos negatif, sin≈0) içerik; boşluk serisi temiz.

Sistematik ölçüm + mekanizma testi:
  Eğer log M'nin asal dalgasına yanıtı kuadratik terim taşıyorsa
  (yanıt = w·S + κ·S² + ..., S = Σ a_p cos θ_p), S² çapraz terimleri
  a_p·a_q·[cos(θ_p−θ_q) + cos(θ_p+θ_q)] üretir:
    → toplam ve fark frekanslarında EŞİT genlik,
    → genlik ∝ a_p·a_q (tek katsayı κ her çizgiyi açıklamalı),
    → işaret κ'nın işareti (cos-kilitli).

12 kombinasyon çizgisi × 3 pencere; kontroller: g̃, g̃², 6 asal çizgisi,
sahte frekanslar. Test: amp(pq) / (a_p·a_q) sabit mi?
"""

import numpy as np
from pathlib import Path

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
A = (np.e**2 - 5) / 2
B0, B1 = 2.7580, -0.0543
PRIMES = [2, 3, 5, 7, 11, 13]

# kombinasyonlar: (etiket, frekans-değeri q_eff, p, q, tip)
COMBOS = [
    ("2*3", 6.0, 2, 3, "+"), ("2*5", 10.0, 2, 5, "+"), ("2*7", 14.0, 2, 7, "+"),
    ("3*5", 15.0, 3, 5, "+"), ("3*7", 21.0, 3, 7, "+"), ("5*7", 35.0, 5, 7, "+"),
    ("3/2", 1.5, 3, 2, "-"), ("5/2", 2.5, 5, 2, "-"), ("5/3", 5/3, 5, 3, "-"),
    ("7/2", 3.5, 7, 2, "-"), ("7/3", 7/3, 7, 3, "-"), ("7/5", 1.4, 7, 5, "-"),
]
FAKES = [2.31, 6.7, 11.9]

d41 = np.load(HERE / "41_bigT_windows.npz")
all_ratio = []
for wkey in ["1600k", "350k", "120k"]:
    gaps = d41[f"gaps_{wkey}"]; tmid = d41[f"tmid_{wkey}"]; amps = d41[f"amps_{wkey}"]
    Lw = np.log(tmid / TWO_PI); L = float(Lw.mean())
    g_u = gaps * Lw / TWO_PI
    a_u = amps / np.sqrt(A * Lw + B0 + B1 / Lw)
    a_u /= np.sqrt((a_u**2).mean())
    y = np.log(a_u)
    n = len(y)

    # tasarım: kontroller + asal çizgileri + kombinasyonlar + sahteler
    cols = [np.ones_like(y), g_u, g_u**2]
    for p in PRIMES:
        arg = tmid * np.log(p)
        cols += [np.cos(arg), np.sin(arg)]
    for _, qv, _, _, _ in COMBOS:
        arg = tmid * np.log(qv)
        cols += [np.cos(arg), np.sin(arg)]
    for f in FAKES:
        arg = tmid * np.log(f)
        cols += [np.cos(arg), np.sin(arg)]
    X = np.vstack(cols).T
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    resd = y - X @ b
    se = np.sqrt(resd.var() * np.diag(np.linalg.inv(X.T @ X)))

    # birinci-mertebe iletilmiş genlikler a_p = w_p·p^{-1/2} (bu pencerede)
    a1 = {}
    for i, p in enumerate(PRIMES):
        a1[p] = np.hypot(b[3 + 2 * i], b[4 + 2 * i])

    print(f"\n=== L={L:.2f} ({wkey}) ===")
    print(f"{'çizgi':>5} {'genlik':>8} {'cos':>8} {'sin':>8} {'σ':>6} "
          f"{'a_p·a_q':>8} {'oran κ':>7}")
    coff = 3 + 2 * len(PRIMES)
    for i, (lab, qv, p, q, typ) in enumerate(COMBOS):
        c_, s_ = b[coff + 2 * i], b[coff + 2 * i + 1]
        amp = np.hypot(c_, s_)
        sig = amp / se[coff + 2 * i]
        prod = a1[p] * a1[q]
        ratio = amp / prod if prod > 0 else np.nan
        all_ratio.append((L, lab, typ, amp, sig, ratio, np.sign(c_)))
        print(f"{lab:>5} {amp:>8.4f} {c_:>+8.4f} {s_:>+8.4f} {sig:>6.1f} "
              f"{prod:>8.4f} {ratio:>7.3f}")
    foff = coff + 2 * len(COMBOS)
    fr = [np.hypot(b[foff + 2 * i], b[foff + 2 * i + 1]) for i in range(len(FAKES))]
    print(f"  sahteler: " + " ".join(f"{x:.4f}" for x in fr))

    # boşluk serisi kontrolü (temiz olmalı) — yalnız L=12.45'te
    if wkey == "1600k":
        yg = np.log(g_u)
        bg, *_ = np.linalg.lstsq(X, yg, rcond=None)
        rg = yg - X @ bg
        seg = np.sqrt(rg.var() * np.diag(np.linalg.inv(X.T @ X)))
        sigs = [np.hypot(bg[coff + 2 * i], bg[coff + 2 * i + 1]) / seg[coff + 2 * i]
                for i in range(len(COMBOS))]
        print(f"  BOŞLUK serisi kombinasyon σ'ları: " +
              " ".join(f"{x:.1f}" for x in sigs) + "  (temiz olmalı)")

# mekanizma özeti
ar = [r for r in all_ratio if r[4] > 3]
print(f"\n=== MEKANİZMA TESTİ ===")
print(f"σ>3 çizgi sayısı: {len(ar)}/{len(all_ratio)}")
rat = np.array([r[5] for r in ar])
sgn = np.array([r[6] for r in ar])
print(f"κ oranları: medyan {np.median(rat):.3f}, aralık [{rat.min():.3f}, {rat.max():.3f}]")
print(f"cos işareti: {int((sgn<0).sum())}/{len(sgn)} negatif")
print("(tek-κ modeli tutarsa: oranlar ~sabit, işaretler tek-tip)")
