"""
130 — KALEM CİLASI ÇATALI: KULE-RETİ + t₀/FAZ TESTLERİ (27 Ağustos)
==========================================================================
Kalem düzeltmesi: KESKİN rezonans çiftleri yalnız aynı-taban kuleleri
(Q₁=pQ₂ ⟹ ikisi de p-kuvveti) ve sayısal olarak İHMAL edilebilir
(R_kule(2) ~ 0.007 « 0.9). Yakın-rezonans genişbandı büyük-t₀'da
dekohere. → 126'nın göçünü taşıyan yapı t₀-bağımsız YEREL 2.-mertebe
terimler olmalı (öz-tutarlılık + orta-nokta örnekleme çarpımları).

ÇATAL TESTLERİ (D3-derinlikli lab, L₀=7):
  F1  t₀-TESTİ: aynı kuyruk, pencere e^{0.3} kaydırılmış → R DEĞİŞMEMELİ
      (gerçek verinin pencere-tutarlılığı lab'da; donmuş-gürültü ölür).
  F2  FAZ-KARIŞTIRMA: kuyruk çizgilerinin (τ>0.52) fazları rastgele →
      R ÖLÜRSE koherans esaslı (çift-toplam yolu yaşar, dekoherans
      muhasebesi yeniden düşünülür); YAŞARSA zarf-istatistiği esaslı
      (yerel yol; kalem oraya döner). İki sonuç da kalemi yönlendirir.
  F3  (kontrol) taban çizgilerinin fazları sabit tutulur — mean-wave
      referansı bozulmasın.

SONUÇ (27 Ağustos) — ÇATAL KESİN KONUŞTU:
  F1 ✓ t₀-DEĞİŞMEZ: kaydırılmış pencerede R'ler ±%5 içinde aynı
     (−0.75/−0.89/−0.93/−0.91 → −0.78/−0.87/−0.92/−0.91) — donmuş-
     gürültü ÖLDÜ; gerçek verinin pencere-tutarlılığı lab'da da var.
  F2 — KOHERANS ESASLI: kuyruk fazları karıştırılınca ETKİ ÖLÜYOR
     (R_p +0.3..+0.6'ya döner, R_nn ≈ 0) ve σ_η 0.163→0.243 fırlar
     (cos-kilitli kuyruk örgü tarafından kısmen EMİLİYORMUŞ; karıştırma
     onu serbest bırakıyor — yan bulgu). Zarf-istatistiği yolu ÖLDÜ.
  KALEMİN YENİ SINIR KOŞULLARI (üçü birden): mekanizma (i) koherans-
  esaslı (cos-kilit şart), (ii) t₀-değişmez, (iii) kule-rezonansları
  ihmal (R_kule ~ 0.007). Üçünü birden sağlayan tek yapı: KUYRUĞUN
  p-DALGASININ KENDİSİYLE öz-tutarlı 2.-mertebe etkileşimi — kuyruk-
  kuyruk vuruşları değil, her kuyruk çizgisinin p-dalgasınca yerinden
  oynatılmış örneklenişi (u_p(m) ω_p-kilitli → çarpımlar t₀-değişmez;
  işaret, çizgilerin p-dalgasına göre fazından → karıştırma öldürür ✓).
  DÜNKÜ ANLATI REVİZE: "vuruşlar çizgilere düşer" sezgisi doğru yönde
  ama yanlış defterdi; koheran kanal kuyruk×dalga etkileşimi. Kapanış
  cümlesi güncel: "gaz kendi aritmetiğini tanır — kuyruğu, taşıdığı
  her dalgayla cos-kilitli konuştuğu için." KALAN HESAP (sınırlı ve
  iyi-tanımlı): tek kuyruk çizgisinin p-dalgalı öz-tutarlı 2.-mertebe
  pertürbasyonu, ayna ağırlıklarıyla kuyruk üstünden toplanacak →
  hedef: −2cosπτ·(1+δ(τ)), δ(0.26-0.30) ≈ −0.10 (129'un hassas eğrisi).
"""

import numpy as np
import time
from pathlib import Path
from sympy import primerange, factorint

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
rng = np.random.default_rng(130)
L0 = 7.0
NZ = 40000

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

def grid(t0):
    idx = np.arange(NZ + 1, dtype=float)
    tg = t0 + idx * TWO_PI / L0
    for _ in range(8):
        tg = tg - (rvm_N(tg) - rvm_N(t0) - idx) / (np.log(tg / TWO_PI) / TWO_PI)
    return tg, np.log(tg / TWO_PI) / TWO_PI

QCAP = int(np.exp(1.4 * L0))
QS = pk_list(QCAP)
LOGQ = np.log(np.array(QS))
LAMv = np.array([float(np.log(list(factorint(q).items())[0][0])) for q in QS])
W0 = LAMv / (np.pi * np.sqrt(np.array(QS)) * LOGQ)
TAU0 = LOGQ / L0
KUYRUK = TAU0 > 0.52

def S_field(t, fazlar, chunk=8000):
    s = np.zeros_like(t)
    for s0 in range(0, len(t), chunk):
        tt = t[s0:s0 + chunk]
        s[s0:s0 + chunk] = -(np.sin(np.outer(tt, LOGQ) + fazlar[None, :])
                             * W0[None, :]).sum(axis=1)
    return s

def lab(t0, fazlar, n_iter=40):
    tg, rho = grid(t0)
    u = np.zeros_like(tg)
    for _ in range(n_iter):
        u = 0.5 * u + 0.5 * (-S_field(tg + u, fazlar) / rho)
    return np.sort(tg + u)

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
        out[p] = (f / Lb, Rp, Rn)
    return out, np.sqrt(s2), c1 / s2

T0A = TWO_PI * np.exp(L0)
faz0 = np.zeros(len(QS))
faz_k = faz0.copy()
faz_k[KUYRUK] = rng.uniform(0, TWO_PI, KUYRUK.sum())

print(f"lab: {len(QS)} çizgi ({KUYRUK.sum()} kuyruk); hücreler:", flush=True)
for ad, t0, fz in [("F0-taban", T0A, faz0),
                   ("F1-t0-kaymış", T0A * np.exp(0.3), faz0),
                   ("F2-faz-karışık", T0A, faz_k)]:
    ts = time.time()
    z = lab(t0, fz)
    out, se, c1r = olc(z)
    print(f"\n[{ad}] σ_η={se:.4f} c₁/σ²={c1r:+.2f} ({time.time()-ts:.0f} sn)")
    print(f"{'p':>3} {'τ':>6} {'R_p':>7} {'R_nn':>7} {'−2cosπτ':>8}")
    for p in [2, 3, 5, 7]:
        tau, Rp, Rn = out[p]
        print(f"{p:>3} {tau:>6.3f} {Rp:>+7.3f} {Rn:>+7.3f} "
              f"{-2*np.cos(np.pi*tau):>+8.3f}", flush=True)
