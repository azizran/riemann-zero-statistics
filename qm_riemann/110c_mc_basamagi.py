"""
110c — M-c BASAMAĞI: ÜÇLÜNÜN (R_p, R_nn, D) ORTAK YASASI (25 Ağustos)
==========================================================================
110b'nin kapalı formu (kalemle türetildi, sayılarla üç ondalık doğrulandı):
   R_p(β,τ) = πτ·cot(πτ) + β·cos(πτ)
   — kinematik terim = 88'İN BENEK-GEOMETRİ FONKSİYONU πτcot(πτ)!
Burada: (1) ince (β,τ) ızgarasında R_p formülü + R_nn ve D ölçümü;
konum-jitter sınıfında bond KİLİTLİDİR: beklenti R_nn ≈ β + kin_nn(τ)
(tam güç, cos'suz) — düşük τ'da R_nn/R_p ≈ gerçeğin 2.7'sine yaklaşmalı.
(2) D için ampirik yasa fiti: D ≈ 1 − c·|β|·τ ?
(3) GERÇEKLE HİZALAMA: R_p=−0.85'i veren β_eff(τ) = (−0.85−πτcotπτ)/cosπτ;
o β'da modelin D'si vs gerçek 1−0.36τ → nefesin perdedeki payı.

ÖN-MÜHÜR: β_eff ≈ −2..−2.5; model-D açığı gerçek açığın ~%40-55'i
(110b kestirimi) — kalan pay bond-ağır/daha zengin gürültü çekirdeği.

SONUÇ (25 Ağustos) — KÖPRÜ %60-80 KURULDU, İKİ KAPALI FORM:
  (1) R_p(β,τ) = πτ·cot(πτ) + β·cos(πτ) — 25 hücrede ±0.05 içinde ✓.
      πτcot(πτ) ÜÇÜNCÜ kez sahnede (benek sönümü 88, bin-rampa 102c,
      şimdi nefes kinematiği) — örnekleme geometrisinin ana fonksiyonu.
  (2) D(β,τ) = 1 − 0.111·|β|·τ (ampirik, artık std 0.008; c=0.111
      σ_ds=0.148 içindir, σ-ölçeklemesi AÇIK).
  HİZALAMA (R_p=−0.85 hedefiyle, gerçek tarafta sıfır serbest param):
      τ=0.10/0.15/0.20/0.30 → model açığı gerçeğin %59/62/66/81'i.
      τ=0.45 geçersiz (cosπτ→0, β_eff ıraksar; gerçek R o bölgede
      ölçülmedi). HÜKÜM: NEFES KÖPRÜSÜ ölçülü bölgede perdenin
      %60-80'ini taşıyor; kalan %20-40 (düşük-τ ağırlıklı) daha zengin
      gürültü çekirdeği ister (bond kilidi gerçeğin düşük-τ oranını
      [2.75] vermiyor — konum-jitter sınıfının sınırı).
  AÇIK: c'nin analitik türetimi + σ-ölçeklemesi; yüksek-τ R ölçümü
  (kompozit asallarla?); bond-zengin çekirdek sınıfı (Δ² bileşenleri).
"""

import numpy as np
from pathlib import Path

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
rng = np.random.default_rng(1102)
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
    ds = g / gbar - 1
    ds = ds - ds.mean()
    X = np.vstack([np.ones_like(m), np.cos(om * m), np.sin(om * m),
                   np.cos(2 * om * m), np.sin(2 * om * m)]).T
    b, *_ = np.linalg.lstsq(X, ds, rcond=None)
    A1 = np.hypot(b[1], b[2]); ph = np.arctan2(b[2], b[1])
    eta = ds - X @ b
    s2 = float((eta**2).mean())
    b2, *_ = np.linalg.lstsq(X, eta**2 - s2, rcond=None)
    Rp = (b2[1] * np.cos(ph) + b2[2] * np.sin(ph)) / (2 * s2 * A1)
    ee = eta[:-1] * eta[1:]
    c1 = float(ee.mean())
    mm = 0.5 * (m[:-1] + m[1:])
    X2 = np.vstack([np.ones_like(mm), np.cos(om * mm), np.sin(om * mm),
                    np.cos(2 * om * mm), np.sin(2 * om * mm)]).T
    b3, *_ = np.linalg.lstsq(X2, ee - c1, rcond=None)
    Rnn = (b3[1] * np.cos(ph) + b3[2] * np.sin(ph)) / (2 * c1 * A1)
    return A1, Rp, Rnn

TAUS = [0.10, 0.15, 0.20, 0.30, 0.45]
BETAS = [0.0, -1.0, -2.0, -2.5, -3.0]
print(f"{'β':>5} {'τ':>5} {'D':>6} {'R_p':>7} {'form':>7} {'R_nn':>7}")
SON = {}
for tau in TAUS:
    om = tau * L0
    UA = A_DS / (2 * np.sin(np.pi * tau)) * gbar
    for beta in BETAS:
        w = rng.normal(0, sw, NZ + 1)
        wj = w * np.clip(1 - beta * A_DS * np.sin(om * x), 0.15, None)
        y = np.sort(x + wj)
        z = np.sort(y + UA * np.cos(om * y))
        A1, Rp, Rnn = olc(z, om)
        form = np.pi * tau / np.tan(np.pi * tau) + beta * np.cos(np.pi * tau)
        SON[(beta, tau)] = (A1 / A_DS, Rp, Rnn)
        print(f"{beta:>5.1f} {tau:>5.2f} {A1/A_DS:>6.3f} {Rp:>+7.3f} "
              f"{form:>+7.3f} {Rnn:>+7.3f}", flush=True)

print("\nAmpirik D-yasası fiti (D = 1 − c·|β|·τ):")
xs, ys = [], []
for (b, t), (D, Rp, Rnn) in SON.items():
    if b < 0:
        xs.append(abs(b) * t); ys.append(1 - D)
c = float(np.polyfit(xs, ys, 1)[0])
res = np.std(np.array(ys) - c * np.array(xs))
print(f"  c = {c:.4f}  (artık std {res:.4f})")

print("\nGERÇEKLE HİZALAMA (R_p hedefi −0.85):")
print(f"{'τ':>5} {'β_eff':>6} {'model-D':>8} {'gerçek-D':>8} {'pay':>5}")
for tau in TAUS:
    kin = np.pi * tau / np.tan(np.pi * tau)
    beff = (-0.85 - kin) / np.cos(np.pi * tau)
    Dm = 1 - c * abs(beff) * tau
    Dg = 1 - 0.36 * tau
    print(f"{tau:>5.2f} {beff:>6.2f} {Dm:>8.3f} {Dg:>8.3f} "
          f"{(1-Dm)/(1-Dg):>5.2f}")
