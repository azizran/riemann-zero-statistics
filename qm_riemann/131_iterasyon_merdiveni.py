"""
131 — İTERASYON MERDİVENİ: KOHERANS KANALININ DİSEKSİYONU (27 Ağustos)
==========================================================================
130'un sınır koşulları (koherans-esaslı + t₀-değişmez + kule-ihmal)
kalemin üç kanalını tek tek sıkıştırdı: (A) ebeveyn×öz-yanbant F2'ye
bağışık olmalıydı (ama F2 öldürdü), (B) genişband naifçe t₀-rastgele,
(C) kule ihmal. Ayrıştırıcı: ÖZ-TUTARLILIK DERİNLİĞİ.

DENEY: D3-lab (2128 çizgi), n_iter = 1 / 2 / 5 / 40.
ÖN-MÜHÜR:
  I1  R naif tek-adımda (n_iter=1: u = −S(t)/ρ̄, yanbant üretimi yok)
      TAM GÜÇTEYSE → kanal ölçüm-zinciri çarpımları + genişband (B);
      t₀-değişmezliğin kaynağı öz-tutarlılık değil örnekleme cebiri.
  I2  R iterasyonla BÜYÜYORSA → yanbant üretimi (öz-tutarlılık) esas;
      kanal (A)-türevi ama F2-uyumlu bir varyantı (ör. yanbant×yanbant
      veya karışık-faz-duyarlı üçlü çarpımlar) — kalem oraya döner.
Not: naif lab'da gap-kesişmeleri olabilir (106'da 671 idi; burada tail
derin) — kesişme sayısı raporlanır; sort sonrası ds tanımlıdır.

SONUÇ (27 Ağustos) — I2: ÖZ-TUTARLILIK ESAS, ~5 İTERASYONDA YAKINSAR:
  n_iter=1 (naif): R zayıf (R_nn −0.33..−0.84), σ_η 0.297 şişkin,
    916 kesişme — naif süperpozisyon yasayı ÜRETMİYOR.
  n_iter=2: kesişme 1'e düşer, σ_η 0.155; R yarı-yolda (−0.99..−1.22).
  n_iter=5: YAKINSADI — R_nn −1.81/−1.78/−1.66/−1.54 ≈ n_iter=40 ≈
    hedef −2cosπτ. Etki İLK BİRKAÇ öz-tutarlı düzeltmede doğuyor
    (pertürbatif rejim — kalem için iyi haber).
  KALEMİN NİHAİ NESNESİ KESİNLEŞTİ: 2. iterasyonun çapraz terimi
    δu = −S'(t)·u₁/ρ̄ = S'S/ρ̄² — ARİTMETİK ALANIN SAYIM-ÖZDEŞLİĞİ
    NONLİNEERLİĞİ ÜZERİNDEN KENDİSİYLE ETKİLEŞİMİ. Çift-çarpım yapısı
    faz-duyarlı (F2 ✓: φ_Q−φ_Q' taşır), yakınsak (I2 ✓), t₀-kararlılığı
    sayısal olarak sabit (F1 ✓; cebirsel muhasebesi hesabın parçası).
  KALAN (sınırlı, tek nesne): S'S çift-toplamının gap-fark +
    η²-regresyon zincirinden okunması, ayna ağırlıklı kuyruk toplamı →
    hedef −2cosπτ·(1+δ(τ)), δ≈−0.10 @ 0.26-0.30. Şiirsel hali:
    "gaz, aritmetik alanının kendi türeviyle çarpımını taşır —
    tanıma S'·S'tir."
"""

import numpy as np
import time
from pathlib import Path
from sympy import primerange, factorint

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
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

t0 = TWO_PI * np.exp(L0)
idx = np.arange(NZ + 1, dtype=float)
tg = t0 + idx * TWO_PI / L0
for _ in range(8):
    tg = tg - (rvm_N(tg) - rvm_N(t0) - idx) / (np.log(tg / TWO_PI) / TWO_PI)
rho = np.log(tg / TWO_PI) / TWO_PI

QCAP = int(np.exp(1.4 * L0))
QS = pk_list(QCAP)
LOGQ = np.log(np.array(QS))
LAMv = np.array([float(np.log(list(factorint(q).items())[0][0])) for q in QS])
W0 = LAMv / (np.pi * np.sqrt(np.array(QS)) * LOGQ)

def S_field(t, chunk=8000):
    s = np.zeros_like(t)
    for s0 in range(0, len(t), chunk):
        tt = t[s0:s0 + chunk]
        s[s0:s0 + chunk] = -(np.sin(np.outer(tt, LOGQ)) * W0[None, :]).sum(axis=1)
    return s

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

print(f"lab: {len(QS)} çizgi; iterasyon merdiveni:", flush=True)
for n_iter in [1, 2, 5, 40]:
    ts = time.time()
    u = np.zeros_like(tg)
    for _ in range(n_iter):
        u = 0.5 * u + 0.5 * (-S_field(tg + u) / rho) if n_iter > 1 else \
            (-S_field(tg) / rho)
    zraw = tg + u
    kes = int((np.diff(zraw) <= 0).sum())
    z = np.sort(zraw)
    out, se, c1r = olc(z)
    print(f"\n[n_iter={n_iter}] kesişme={kes} σ_η={se:.4f} c₁/σ²={c1r:+.2f} "
          f"({time.time()-ts:.0f} sn)")
    print(f"{'p':>3} {'τ':>6} {'R_p':>7} {'R_nn':>7} {'−2cosπτ':>8}")
    for p in [2, 3, 5, 7]:
        tau, Rp, Rn = out[p]
        print(f"{p:>3} {tau:>6.3f} {Rp:>+7.3f} {Rn:>+7.3f} "
              f"{-2*np.cos(np.pi*tau):>+8.3f}", flush=True)
