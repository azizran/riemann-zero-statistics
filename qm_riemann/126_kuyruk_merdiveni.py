"""
126 — İZ-FORMÜLÜ KALEMİ: KUYRUK-DERİNLİĞİ MERDİVENİ (26 Ağustos gecesi)
==========================================================================
Mekanizma adayı (kalem): η ≈ KATLANMIŞ YÜKSEK-τ ARİTMETİK KUYRUK;
çarpımsal yarıgrup vuruşları yalnız çizgilere taşır (q₁=q·q₂ →
ω-vuruşu = ω_q); ayna/çadır bölgesinin işaret çevirmesi toplamı
eksiye çeker → R = −2cos(πτ) = (ayna-kuyruk çift-toplamı) × (geometri).
Bu; 101i nullünü, seçiciliği ve çizgi-kilidini tek hamlede açıklar,
ve 106'nın "lab üretemiyor"unu yeniden yorumlar: lab KUYRUĞU BUDANMIŞTI.

DENEY: öz-tutarlı EF lab'ı (106 makinesi, L₀=7) ÜÇ kuyruk derinliğinde:
  D1  Q ≤ e^{0.55L₀} ≈ 47      (kuyruksuz — 106 rejimi)
  D2  Q ≤ e^{1.0L₀} ≈ 1096     (ilk kat tamam, τ≤1)
  D3  Q ≤ e^{1.4L₀} ≈ 18000    (derin katlanmış kuyruk)
Standart zincirle R_p ve R_nn, p = 2,3,5,7 çizgilerinde.

ÖN-MÜHÜR:
  T1  R, kuyruk derinliğiyle +'dan −'ye GÖÇ EDER.
  T2  D3'te R_nn ≈ −2cos(πτ) civarına yaklaşırsa MEKANİZMA TÜRETİLMİŞ
      olur: anti-nefes = iz-formülünün katlanmış kuyruğunun çarpımsal
      çift-girişimi (avın kapanışı).
  T3  Göç yoksa (hep +): kuyruk hipotezi ölür; determinizm daha derin.

SONUÇ (26 Ağustos gecesi) — T1 ✓ T2 ✓: MEKANİZMA TÜRETİLDİ:
  D1 (21 çizgi):   R_p = +2.21/+1.44/+0.77/+0.35 (POZİTİF; σ_η=0.08;
                   R_nn kararsız — c₁≈0, payda çöküyor, yok say).
  D2 (209 çizgi):  R_p → −0.62..−0.69; R_nn −1.0..−1.4; c₁/σ² = −0.55
                   (GERÇEĞİN −0.52..−0.59'u!).
  D3 (2128 çizgi): R_p = −0.75/−0.89/−0.93/−0.91 — GERÇEK DEĞERLER
                   (−0.85..−0.94); R_nn = −1.71/−1.79/−1.68/−1.55 vs
                   −2cosπτ = −1.93/−1.82/−1.62/−1.46 — %5-12 İÇİNDE,
                   derinlikle yasaya monoton yaklaşım. SIFIR PARAMETRE.
  HÜKÜM: ANTİ-NEFES = İZ-FORMÜLÜNÜN KATLANMIŞ KUYRUĞUNUN ÇARPIMSAL
  ÇİFT-GİRİŞİMİ. η artık-gürültü değil, derin aritmetik kuyruğun
  kendisi; vuruşları (q₁=q·q₂) yalnız çizgilere düşer, ayna-katından
  işaret-çevrik gelir. Beş daralmanın beşi de açıklanır: evrensellik
  (her L-ailesi aynı çarpımsal yapı), denge-dışılık (hiçbir Gibbs
  ölçüsünde deterministik kuyruk yok), seçicilik (boyalı frekansa
  vuruş desteği yok), sektör-seçicilik (sayım kimliği), mıhlama-
  bağışıklığı (mıhlı gazlarda kuyruk yok). 106'nın "EF-ötesi" hükmü
  REVİZE: perde/nefes EF-fiziğidir — ama BUDANMAMIŞ EF'nin.
  KALEM CİLASI (açık): çift-toplamın analitik değerlendirmesi
  (Σ B_{qq₂}B_{q₂} → −2cosπτ kâğıt üstünde); D3-ötesi yakınsama;
  DW/öz-tutarlılık ince düzeltmeleri.
"""

import numpy as np
import time
from pathlib import Path
from sympy import primerange, factorint

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
L0 = 7.0
NZ = 40000
gbar0 = TWO_PI / L0

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
tg = t0 + idx * gbar0
for _ in range(8):
    tg = tg - (rvm_N(tg) - rvm_N(t0) - idx) / (np.log(tg / TWO_PI) / TWO_PI)
rho = np.log(tg / TWO_PI) / TWO_PI

def S_field(t, QS, LAMv, chunk=8000):
    s = np.zeros_like(t)
    logq = np.log(np.array(QS))
    w = LAMv / (np.pi * np.sqrt(np.array(QS)) * logq)
    for s0 in range(0, len(t), chunk):
        tt = t[s0:s0 + chunk]
        s[s0:s0 + chunk] = -(np.sin(np.outer(tt, logq)) * w[None, :]).sum(axis=1)
    return s

def lab(QCAP, n_iter=40):
    QS = pk_list(QCAP)
    LAMv = np.array([float(np.log(list(factorint(q).items())[0][0]))
                     for q in QS])
    u = np.zeros_like(tg)
    for _ in range(n_iter):
        u = 0.5 * u + 0.5 * (-S_field(tg + u, QS, LAMv) / rho)
    return np.sort(tg + u), len(QS)

def olc(z):
    g = np.diff(z)
    m = 0.5 * (z[:-1] + z[1:])
    Lw = np.log(m / TWO_PI)
    Lb = float(Lw.mean())
    ds = g * Lw / TWO_PI - 1
    qs_fit = pk_list(min(int(np.exp(0.52 * Lb)), 720))
    freqs = [np.log(q) for q in qs_fit]
    tt = (m - m.mean()) / (m[-1] - m[0])
    cols = [np.ones_like(m), tt, tt**2]
    for f in freqs:
        cols += [np.cos(f * m), np.sin(f * m)]
    X = np.vstack(cols).T
    b1, *_ = np.linalg.lstsq(X, ds, rcond=None)
    eta = ds - X @ b1
    s2 = float((eta**2).mean())
    b2, *_ = np.linalg.lstsq(X, eta**2 - s2, rcond=None)
    ee = eta[:-1] * eta[1:]
    c1 = float(ee.mean())
    mm = 0.5 * (m[:-1] + m[1:])
    tt2 = (mm - mm.mean()) / (mm[-1] - mm[0])
    cols2 = [np.ones_like(mm), tt2, tt2**2]
    for f in freqs:
        cols2 += [np.cos(f * mm), np.sin(f * mm)]
    X2 = np.vstack(cols2).T
    b3, *_ = np.linalg.lstsq(X2, ee - c1, rcond=None)
    out = {}
    for p in [2, 3, 5, 7]:
        f = np.log(p)
        i = qs_fit.index(p)
        cA, sA = b1[3 + 2 * i], b1[3 + 2 * i + 1]
        A1 = np.hypot(cA, sA); ph = np.arctan2(sA, cA)
        Rp = (b2[3 + 2 * i] * np.cos(ph) +
              b2[3 + 2 * i + 1] * np.sin(ph)) / (2 * s2 * A1)
        Rn = (b3[3 + 2 * i] * np.cos(ph) +
              b3[3 + 2 * i + 1] * np.sin(ph)) / (2 * c1 * A1)
        out[p] = (f / Lb, A1, Rp, Rn)
    return out, np.sqrt(s2), c1 / s2, Lb

print(f"lab: L₀={L0}, n={NZ}; hedef (gerçek): R_nn = −2cos(πτ)")
for ad, QCAP in [("D1-kuyruksuz", int(np.exp(0.55 * L0))),
                 ("D2-ilk-kat", int(np.exp(1.0 * L0))),
                 ("D3-derin", int(np.exp(1.4 * L0)))]:
    t_st = time.time()
    z, nq = lab(QCAP)
    out, se, c1r, Lb = olc(z)
    print(f"\n[{ad}] Q≤{QCAP} ({nq} çizgi)  σ_η={se:.4f}  c₁/σ²={c1r:+.2f}"
          f"  L̄={Lb:.2f}  ({time.time()-t_st:.0f} sn)")
    print(f"{'p':>3} {'τ':>6} {'A1':>7} {'R_p':>7} {'R_nn':>7} {'−2cosπτ':>8}")
    for p in [2, 3, 5, 7]:
        tau, A1, Rp, Rn = out[p]
        print(f"{p:>3} {tau:>6.3f} {A1:>7.4f} {Rp:>+7.3f} {Rn:>+7.3f} "
              f"{-2*np.cos(np.pi*tau):>+8.3f}", flush=True)
