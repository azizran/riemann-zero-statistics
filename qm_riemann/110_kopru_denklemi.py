"""
110 — KÖPRÜ DENKLEMİ MERDİVENİ: NEFES → PERDE SENTETİK TESTİ (24 Ağu)
==========================================================================
Soru: tam-ters-nefes (R≈−1) tek başına ölçülen perdeyi (D≈1−0.36τ)
ÜRETİR Mİ? Kalem-türetme (tek-gap teleskop, 2. mertebe):
   D = 1/(1 + (κ/2)tan(κ/2)·σ_η²·|R|)  → τ=0.3'te yalnız 0.97 (ölçülen
   0.884) ve biçim κ² — YETERSİZ görünüyor; eksik aday: komşu-bağı
   (R_nn) ve kümülatif-faz terimleri. Sentetik merdivenle sınanır.

MERDİVEN (hepsi aynı 108b/109 ölçüm zinciriyle okunur):
  M-a  boyalı kontrol: pürüzsüz örgü + dalga + iid konum-jitteri
       (gap-gürültüsü c₁/σ² = −0.5 ✓ gerçeğe yakın). Beklenti: D=1, R=+1.
  M-b  ANTİ-NEFES: jitter genliği dalga fazıyla modüle —
       w_n → w_n(1 + R_b·W(x_n)/A·...), R_b = −1 reçeteli.
       ÖN-MÜHÜR: R ölçümü ≈ −1 çıkmalı (reçete sağlaması); D ise
       teleskop formülü kadar (≈0.97 @ τ0.3) düşmeli. Ölçülen 0.884'e
       İNMİYORSA → nefes tek başına köprüyü kurmuyor; M-c'ye geç.
  M-c  anti-nefes + KOMŞU-BAĞI modülasyonu (jitter'i ardışık-korele
       yapan bileşen de fazla modüle) — R_nn reçeteli ~−1.5.
Karar: hangi merdiven basamağı D(τ)'nun hem düzeyini hem lineer-τ
biçimini üretiyorsa, köprünün mekanik içeriği odur; kapalı form o
basamağın analitiğinden yazılır.
"""

import numpy as np
from pathlib import Path

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
rng = np.random.default_rng(110)
L0 = 10.37
NZ = 200000
gbar = TWO_PI / L0

def rvm_N(t):
    x = t / TWO_PI
    return x * np.log(x / np.e) + 7 / 8

t0 = TWO_PI * np.exp(L0)
idx = np.arange(NZ + 1, dtype=float)
x = t0 + idx * gbar
for _ in range(8):
    x = x - (rvm_N(x) - rvm_N(t0) - idx) / (np.log(x / TWO_PI) / TWO_PI)

SIG_DS = 0.148                       # hedef artık-gürültü (109: σ_η)
sw = SIG_DS / np.sqrt(2) * gbar      # iid konum-jitteri (gap c₁/σ²=−0.5)
A_DS = 0.10                          # dalga genliği (ds birimi)

def olc(z, om, tau):
    g = np.diff(z)
    m = 0.5 * (z[:-1] + z[1:])
    Lw = np.log(m / TWO_PI)
    ds = g * Lw / TWO_PI - 1
    tt = (m - m.mean()) / (m[-1] - m[0])
    X = np.vstack([np.ones_like(m), tt, tt**2,
                   np.cos(om * m), np.sin(om * m),
                   np.cos(2 * om * m), np.sin(2 * om * m)]).T
    b, *_ = np.linalg.lstsq(X, ds, rcond=None)
    A1 = np.hypot(b[3], b[4]); ph = np.arctan2(b[4], b[3])
    eta = ds - X @ b
    s2 = float((eta**2).mean())
    b2, *_ = np.linalg.lstsq(X, eta**2 - s2, rcond=None)
    Rp = (b2[3] * np.cos(ph) + b2[4] * np.sin(ph)) / (2 * s2 * A1)
    ee = eta[:-1] * eta[1:]
    c1 = float(ee.mean())
    mm = 0.5 * (m[:-1] + m[1:])
    tt2 = (mm - mm.mean()) / (mm[-1] - mm[0])
    X2 = np.vstack([np.ones_like(mm), tt2, tt2**2,
                    np.cos(om * mm), np.sin(om * mm),
                    np.cos(2 * om * mm), np.sin(2 * om * mm)]).T
    b3, *_ = np.linalg.lstsq(X2, ee - c1, rcond=None)
    Rnn = (b3[3] * np.cos(ph) + b3[4] * np.sin(ph)) / (2 * c1 * A1)
    return A1, Rp, Rnn, s2, c1 / s2

print(f"{'hücre':>6} {'τ':>5} {'D=A1/A':>7} {'R_p':>7} {'R_nn':>7} "
      f"{'σ_η':>6} {'c₁/σ²':>6}   [teleskop-öngörü D]")
for tau in [0.15, 0.30, 0.45]:
    om = tau * L0
    kap = 2 * np.pi * tau
    UA = A_DS / (2 * np.sin(np.pi * tau)) * gbar   # konum-dalga genliği
    Wx = np.cos(om * x)
    tel = 1.0 / (1 + (kap / 2) * np.tan(kap / 2) * SIG_DS**2)
    for ad, mod, modnn in [("M-a", 0.0, 0.0), ("M-b", -1.0, 0.0),
                           ("M-c", -1.0, -1.5)]:
        w = rng.normal(0, sw, NZ + 1)
        gen = 1 + mod * A_DS * Wx
        wj = w * np.clip(gen, 0.2, None)
        if modnn != 0.0:
            ort = 0.5 * (w[1:] + w[:-1])           # komşu-ortak bileşen
            wj = wj + modnn * A_DS * np.concatenate(
                [[0], ort * np.cos(om * 0.5 * (x[1:] + x[:-1]))])
        z = np.sort(x + UA * Wx + wj)
        A1, Rp, Rnn, s2, c1r = olc(z, om, tau)
        print(f"{ad:>6} {tau:>5.2f} {A1/A_DS:>7.3f} {Rp:>+7.3f} "
              f"{Rnn:>+7.3f} {np.sqrt(s2):>6.3f} {c1r:>+6.2f}   [{tel:.3f}]",
              flush=True)
print("\nGERÇEK HEDEF: D_fiz ≈ 0.913/0.884/0.845 @ τ=0.2/0.3/0.4 (~1−0.36τ);"
      "\n              R_p(ham ζ) ≈ −0.85, R_nn ≈ −1.5, c₁/σ² = −0.52")
