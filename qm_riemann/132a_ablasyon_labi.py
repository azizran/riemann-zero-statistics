"""
132a — ABLASYON LAB'I: MERTEBE AYRIŞTIRMASI (27 Ağustos)
==========================================================================
131 kalemin nihai nesnesini "2. iterasyonun çapraz terimi"
δu₂ = −S'(t)·u₁(t)/ρ̄ olarak işaretledi, ama bunu DOĞRUDAN ölçmedi:
merdivende yalnız n_iter = 1/2/5/40 vardı. Bu lab, örgüyü mertebe
mertebe kurup R'nin hangi mertebede doğduğunu ayrıştırır.

Cebir: S(t) = −Σ_Q a_Q sin(ω_Q t), ω_Q = log Q, a_Q = Λ(Q)/(π√Q logQ);
u₁ = −S(t)/ρ̄ = (1/ρ̄) Σ a_Q sin(ω_Q t);
u₂ = −S(t+u₁)/ρ̄ ≈ u₁ − S'(t)u₁/ρ̄ ⟹ δu₂ = −S'(t)u₁/ρ̄ = S'S/ρ̄².

ÖRGÜLER (aynı ızgara, aynı olc ölçüm zinciri — 131 kalıbı):
  A0  z = t + u₁                      (naif; 131'de zayıf — referans)
  A1  z = t + u₁ + δu₂                (KESİN 2. mertebe, iterasyonsuz)
  A2  z = t + u₂,  u₂ = −S(t+u₁)/ρ̄    (A1 + tüm yüksek terimler)
  A3  n_iter=5 sönümlü (131 referansı)
  A4  n_iter=40 sönümlü (yakınsak referans)

ÖN-MÜHÜR:
  M1  A1, A3'ün R_nn'inin ≥%70'ini veriyorsa 2.-mertebe terim kalemin
      nesnesi olarak DOĞRULANIR; kalan pay yüksek mertebelerin cilası.
  M2  A1 zayıf ama A2 güçlüyse etki 3.+ mertebede doğuyor demektir —
      kalem Taylor'ı bir basamak daha açmalı (rapor edilir).
  M3  A1 ≈ A2 ≈ A3 ise seri 2. mertebede pratikte kapanıyor.

SONUÇ (27 Ağustos) — M1 ✗ M2 ✗ M3 ✗: SERİ YAKINSAMIYOR, KESİLEMEZ:
  A0 (u₁ naif):   kesişme 916, σ_η 0.297, c₁/σ² −0.47
                  R_nn = −0.33/−0.66/−0.82/−0.84  (131'in n_iter=1'i ✓)
  A1 (u₁+δu₂):    kesişme 3429(!), σ_η 0.437, c₁/σ² −0.33
                  R_nn = −0.12/−0.29/−0.39/−0.54  — A0'DAN DA KÖTÜ.
  A2 (tam u₂):    kesişme 136, σ_η 0.270, c₁/σ² −0.37
                  R_nn = +0.28/−0.68/−0.84/−0.98  — A0 seviyesinde.
  A3 (n=5):       kesişme 6, σ_η 0.163, c₁/σ² −0.57
                  R_nn = −1.81/−1.78/−1.66/−1.54 (hedefin %94/98/102/105)
  A4 (n=40):      A3'ün %99.3'ü — yakınsama n=5'te ✓ (131 doğrulandı).
  PAY (R_nn/A3): A0 %40, A1 %20, A2 %34, A4 %99.
  SEBEP ÖLÇÜLDÜ — AÇILIM PARAMETRESİ BİRDEN BÜYÜK:
    rms(S'/ρ̄) = 1.172  (kuyruk 1.113 / taban 0.419), max = 4.24
    rms(δu₂)/rms(u₁) = 1.07   [ḡ = 0.898; rms(δu₂) = 0.333]
  Yani δu₂ "düzeltme" değil, öncü terimden BÜYÜK; Lagrange serisi
  z = t − S/ρ̄ + S'S/ρ̄² − … yakınsama yarıçapının DIŞINDA. λ=1'de örgü
  yıkılıyor (kesişme patlaması) ve ölçüm zinciri anlamsızlaşıyor.
  HÜKÜM: 131'in "kalemin nihai nesnesi δu₂ = S'S/ρ̄²" işareti YÖN olarak
  doğru (nonlineerlik oradan), DEFTER olarak yanlış: yasayı üreten şey
  SONLU BİR MERTEBE DEĞİL, öz-tutarlılığın RESUMLADIĞI seri. Kalem
  Taylor kesmez — yerdeğiştirmeyi argümanın içinde tutmalı:
      ds = −x/(1+x),  x = S'(z)/ρ̄,  z = t + u,  u = −S(z)/ρ̄.
  Diseksiyon bu yüzden iki katta yürütülür (→ 132b): Taylor blokları
  yalnız DOĞRUSAL-TEPKİ rejiminde (λ«1, eğim ölçümü) anlamlı; asıl
  atıf resumlanmış yerdeğiştirme sınıflarıyla (C-hücreleri) yapılır.
"""

import numpy as np
import time
from sympy import primerange, factorint

TWO_PI = 2 * np.pi
L0 = 7.0
NZ = 40000
PS = [2, 3, 5, 7]


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
    """S(t) = −Σ_Q a_Q sin(ω_Q t)"""
    s = np.zeros_like(t)
    for s0 in range(0, len(t), chunk):
        tt = t[s0:s0 + chunk]
        s[s0:s0 + chunk] = -(np.sin(np.outer(tt, LOGQ)) * W0[None, :]).sum(axis=1)
    return s


def Sp_field(t, chunk=8000):
    """S'(t) = −Σ_Q a_Q ω_Q cos(ω_Q t)"""
    w = W0 * LOGQ
    s = np.zeros_like(t)
    for s0 in range(0, len(t), chunk):
        tt = t[s0:s0 + chunk]
        s[s0:s0 + chunk] = -(np.cos(np.outer(tt, LOGQ)) * w[None, :]).sum(axis=1)
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
    for p in PS:
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


# ---- örgüler --------------------------------------------------------
u1 = -S_field(tg) / rho


def orgu(ad):
    if ad == "A0":
        return tg + u1
    if ad == "A1":
        return tg + u1 + (-Sp_field(tg) * u1 / rho)
    if ad == "A2":
        return tg + (-S_field(tg + u1) / rho)
    n_iter = 5 if ad == "A3" else 40
    u = np.zeros_like(tg)
    for _ in range(n_iter):
        u = 0.5 * u + 0.5 * (-S_field(tg + u) / rho)
    return tg + u


ETIKET = {"A0": "u₁ (naif)",
          "A1": "u₁+δu₂ (kesin 2. mertebe)",
          "A2": "u₂ = −S(t+u₁)/ρ̄",
          "A3": "n_iter=5 (131 ref)",
          "A4": "n_iter=40"}

print(f"lab: {len(QS)} çizgi (L₀={L0}, NZ={NZ}); ablasyon merdiveni", flush=True)
res = {}
for ad in ["A0", "A1", "A2", "A3", "A4"]:
    ts = time.time()
    zraw = orgu(ad)
    kes = int((np.diff(zraw) <= 0).sum())
    out, se, c1r = olc(np.sort(zraw))
    res[ad] = out
    print(f"\n[{ad}] {ETIKET[ad]}  kesişme={kes} σ_η={se:.4f} "
          f"c₁/σ²={c1r:+.2f} ({time.time()-ts:.0f} sn)")
    print(f"{'p':>3} {'τ':>6} {'R_p':>7} {'R_nn':>7} {'−2cosπτ':>8} {'R_nn/hedef':>11}")
    for p in PS:
        tau, Rp, Rn = out[p]
        hed = -2 * np.cos(np.pi * tau)
        print(f"{p:>3} {tau:>6.3f} {Rp:>+7.3f} {Rn:>+7.3f} {hed:>+8.3f} "
              f"{Rn/hed:>10.2%}", flush=True)

print("\n=== PAY TABLOSU: R_nn(hücre) / R_nn(A3) ===", flush=True)
print(f"{'p':>3} " + " ".join(f"{a:>8}" for a in ["A0", "A1", "A2", "A4"]))
for p in PS:
    ref = res["A3"][p][2]
    print(f"{p:>3} " + " ".join(f"{res[a][p][2]/ref:>7.1%}" for a in ["A0", "A1", "A2", "A4"]),
          flush=True)
for a in ["A0", "A1", "A2", "A4"]:
    pay = np.mean([res[a][p][2] / res["A3"][p][2] for p in PS])
    print(f"ORT {a}: {pay:.1%}", flush=True)
