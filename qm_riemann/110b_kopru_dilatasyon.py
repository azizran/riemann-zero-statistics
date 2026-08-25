"""
110b — KÖPRÜ MERDİVENİ, DİLATASYON MİMARİSİYLE (24 Ağustos)
==========================================================================
110'un dersi: eklemeli mimari (z = x + dalga + jitter) nefes üretmez
(R≈0) — boyalı-dilatasyon teoremi gürültünün dalganın İÇİNDEN geçmesini
ister. Doğru mimari: önce gürültülü örgü y = x + jitter, sonra
z = y + U·cos(ω y) (dilatasyon). Anti-nefes reçetesi: jitter genliği
β·A·cos(ωx) ile modüle; net R ≈ 1 + β·cos(κ/2) → β=−2 ~ R≈−1.

ÖN-MÜHÜR:
  N1  β=0 (saf dilatasyon): R=+1, D=1 (kinematik teorem sağlaması).
  N2  β=−2: R≈−1; D'nin düşüşü kalem-öngörüsü ΔA1/A ≈ β·σ_η²·cos(κ/2)
      ≈ −0.04 civarı (küçük!). Ölçülen perde (0.88 @ τ0.3) BURADAN
      ÇIKMIYORSA hüküm: nefes tek başına köprüyü kurmaz — perdenin
      gövdesi kaynak-renormalizasyonudur (gerçek çok-cisim verteksi);
      nefes onun faz-uzayı gölgesidir, nedeni değil.

SONUÇ (24 Ağustos) — KÖPRÜNÜN AYAKLARI DİKİLDİ:
  FAZ DERSİ (K4-K6 teşhisi): nefes GAP-dalgasının (türev, −sin)
  fazında yaşar; yerdeğiştirme (cos) fazına boyanan modülasyon ortalama-
  dalgaya DİK kalır ve R'de görünmez — 110'un ilk nullünün açıklaması,
  ve başlı başına mekanizma bilgisi.
  DÜZELTİLMİŞ MERDİVEN: N1 ✓ (β=0: R +0.93/+0.69/+0.26, D≈1).
  β=−2: R −0.86/−0.49/−0.11, D 0.982/0.961/0.913.
  β=−3: R −1.78/−1.19/−0.35, D 0.971/0.917/0.872.
  HÜKÜM: anti-nefes perdeyi MEKANİK üretir ve biçim LİNEER-τ (kalem
  formülünün κ²'si değil — gerçek 1−0.36τ biçimi ✓). Büyüklük: orta/
  yüksek τ'da gerçek açığın ~%70-75'i (R orta-τ'da eşlenince);
  düşük-τ ucu eksik. KALAN AYAK: komşu-bağı dalgası kanalı (gerçekte
  R_nn≈−1.5, burada reçetesiz) + gerçek gürültünün iid-dışı yapısı —
  M-c reçeteli basamak ve kapalı form SIRADAKİ oturuma.
"""

import numpy as np
from pathlib import Path

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
rng = np.random.default_rng(1100)
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

SIG_DS = 0.148
sw = SIG_DS / np.sqrt(2) * gbar
A_DS = 0.10

def olc(z, om):
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
    return A1, Rp, s2

print(f"{'β':>4} {'τ':>5} {'D=A1/A':>7} {'R_p':>7} {'σ_η':>6}"
      f"   [kalem ΔD≈β·σ²·cos(κ/2)]")
for tau in [0.15, 0.30, 0.45]:
    om = tau * L0
    kap = 2 * np.pi * tau
    UA = A_DS / (2 * np.sin(np.pi * tau)) * gbar
    for beta in [0.0, -1.0, -2.0, -3.0]:
        w = rng.normal(0, sw, NZ + 1)
        # DÜZELTME (K4-K6 teşhisi): modülasyon GAP-dalgası fazına
        # (türev, −sin) boyanmalı — cos fazı ortalama-dalgaya DİK kalır
        wj = w * np.clip(1 - beta * A_DS * np.sin(om * x), 0.15, None)
        y = np.sort(x + wj)
        z = np.sort(y + UA * np.cos(om * y))
        A1, Rp, s2 = olc(z, om)
        kalem = beta * SIG_DS**2 * np.cos(kap / 2)
        print(f"{beta:>4.0f} {tau:>5.2f} {A1/A_DS:>7.3f} {Rp:>+7.3f} "
              f"{np.sqrt(s2):>6.3f}   [{1+kalem:.3f}]", flush=True)
print("\nGERÇEK: D_fiz ≈ 0.93/0.88/0.83 @ τ=0.15/0.30/0.45; R_p ≈ −0.85")
