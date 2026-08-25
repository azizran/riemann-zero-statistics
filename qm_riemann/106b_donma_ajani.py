"""
106b — DONMA AJANI: LAB-GERÇEK FARKI DÜŞÜK-ω'DA (24 Ağustos)
==========================================================================
106'nın işaretleri: (1) v_sc ≈ katı (perde öz-tutarlılık değil; ama
naif kanal nan — burada lineer kanalla temiz oran), (2) lab spontane
eğrisi τ≥0.15'te gerçeği %2-6 içinde üretir (CUE-uyumu = örnekleme
optiği) AMA lab τ=0.04'te DONUK DEĞİL (0.66 vs 0.04) → donma EF-ötesi.

ÖN-MÜHÜRLÜ ÖNGÖRÜLER:
  P1' lineer kanalda D_lab = v_sc/v_naif ≈ 0.9-1.05 tüm τ'da, düşüş
      yok → perde (0.91→0.48) gerçek gazın EF-ötesi yanıtı (92'nin
      "etkileşim kanıtı" güçlenir: lab kinematiği perdeyi ÜRETEMEZ).
  P2' lab orta-nokta I(ω) düşük-ω'da (ω<1, çizgi-dışı, Hann) gerçeğin
      ÇOK üstünde (gerçek 102a: 4.12e-6; lab için ≥100× bekliyorum) —
      gerçek gaz düşük-ω yan bantlarını İPTAL ediyor (donma ajanı).
  P3' lab sabit-ω taramasında vuruş tepeleri VAR (log(3/2)=0.4055'te
      101i-tarzı ±0.004 bantla yerel tepe); gerçekte YOKTU (101i) →
      iptalin imzası: gerçek gaz çizgi-vuruşlarını söndürür.

SONUÇ (24 Ağustos) — İKİ-DAL GİZEMİ ÇÖZÜLDÜ, ÜÇ HÜKÜM:
  P1' oran kanalı kesişme-kirli (671 kesişme; yükselen oran 0.95→2.2
    ARTEFAKT — v_naif yüksek ω'da kesişme gürültüsüyle çöküyor). Temiz
    hüküm 106'dan: v_sc ≈ katı-analitik (%±5, düz) → PERDE D(τ)
    GAZIN EF-ÖTESİ ÇOK-CİSİM YANITI; öz-tutarlılık hipotezi RET.
    Asıl açık soru ("perdeyi ne üretir") artık temiz izole.
  P2' lab orta-nokta düşük-ω'da gerçeğin 109 KATI parlak — ve lab
    SIFIR-ızgarası da parlak (4.7e-4 vs gerçek ~7e-16!). YORUM:
    gerçek sıfırların boşlukları TAM boştur çünkü açık formül KESİN
    özdeşliktir (yoğunluk dalgalanması = tam olarak çizgiler); lab'ın
    yerdeğiştirme inşası (t+u) yaklaşık ters-çevirimdir ve u²'nin
    vuruş terimlerini SAHTE olarak üretir. → DONMANIN AÇIKLAMASI:
    boşluğun boşluğu aritmetik özdeşlik; lab boşluk modeli DEĞİLDİR
    (orta-ω örnekleme yan bantları için geçerliliği [102b %6] durur).
  P3' lab'da vuruş-tepe kontrastı net değil (0.43-0.81 geniş salınım;
    kontroller vuruşlardan ayrışmıyor) — BELİRSİZ, iddia edilmedi.
  SENTEZ: (1) donma = atomik tayfın kesinliği; (2) spontane CUE-uyumu
    = örnekleme optiği (92-T4 yeniden yorumlanmalı, Not 4 §branches'e
    kayıt); (3) perde = gerçek çok-cisim fiziği. İki dal tek resimde:
    biri özdeşlik, biri optik; fizik yalnız perdede.
"""

import numpy as np
from pathlib import Path
from sympy import primerange, factorint

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
L0 = 10.37
NZ = 60000

def pk_list(lim):
    out = []
    for p in primerange(2, int(lim) + 1):
        q = p
        while q <= lim:
            out.append(q); q *= p
    return sorted(set(out))

QS = pk_list(720)
LAM = {}
for q in QS:
    (pp, kk), = factorint(q).items()
    LAM[q] = float(np.log(pp))

def S_field(t):
    s = np.zeros_like(t)
    for q in QS:
        s -= (LAM[q] / (np.pi * np.sqrt(q) * np.log(q))) * np.sin(t * np.log(q))
    return s

def rvm_N(t):
    x = t / TWO_PI
    return x * np.log(x / np.e) + 7 / 8

t0 = TWO_PI * np.exp(L0)
idx = np.arange(NZ + 1, dtype=float)
t = t0 + idx * TWO_PI / L0
for _ in range(8):
    t = t - (rvm_N(t) - rvm_N(t0) - idx) / (np.log(t / TWO_PI) / TWO_PI)
rho = np.log(t / TWO_PI) / TWO_PI
u_na = -S_field(t) / rho
z_na = t + u_na
u_sc = np.zeros_like(t)
for _ in range(60):
    u_sc = 0.5 * u_sc + 0.5 * (-S_field(t + u_sc) / rho)
z_sc = t + u_sc
print(f"lab hazır; naif kesişme sayısı: {(np.diff(z_na) <= 0).sum()}", flush=True)

# ---- P1': lineer kanal oranı
def v_lin(z, om):
    g = np.diff(z)
    m = 0.5 * (z[:-1] + z[1:])
    ds = g * np.log(m / TWO_PI) / TWO_PI - 1
    X = np.vstack([np.ones_like(m), np.cos(om * m), np.sin(om * m)]).T
    b, *_ = np.linalg.lstsq(X, ds, rcond=None)
    return np.hypot(b[1], b[2])

PRC = [2, 3, 5, 7, 11, 13, 17, 23, 31, 41, 53, 79, 107, 149]
GERCEK = {0.05: 0.908, 0.15: 0.863, 0.25: 0.797, 0.35: 0.726,
          0.45: 0.634, 0.55: 0.484}
print("\nP1' — lineer kanal: D_lab = v_sc/v_naif")
print(f"{'p':>4} {'τ':>6} {'v_naif':>8} {'v_sc':>8} {'D_lab':>6}")
DL = []
for p in PRC:
    om = np.log(p)
    tau = om / L0
    vn, vs = v_lin(z_na, om), v_lin(z_sc, om)
    DL.append((tau, vs / vn))
    print(f"{p:>4} {tau:>6.3f} {vn:>8.5f} {vs:>8.5f} {vs/vn:>6.3f}", flush=True)
DL = np.array(DL)
for c, dger in GERCEK.items():
    m = (DL[:, 0] >= c - 0.05) & (DL[:, 0] < c + 0.05)
    if m.sum():
        print(f"  τ≈{c:.2f}: D_lab = {DL[m, 1].mean():.3f}  vs gerçek {dger:.3f}")

# ---- P2': lab orta-nokta I(ω), düşük-ω çizgi-dışı, Hann
def I_hann(pts, oms):
    w = np.hanning(len(pts))
    sw2 = (w**2).sum()
    out = []
    for s0 in range(0, len(oms), 40):
        ob = oms[s0:s0 + 40]
        F = (np.exp(1j * np.outer(ob, pts)) * w[None, :]).sum(axis=1)
        out.extend((np.abs(F)**2 / sw2).tolist())
    return np.array(out)

LINES = [np.log(q) for q in QS if q <= 200]
mids_sc = 0.5 * (z_sc[:-1] + z_sc[1:])
rng = np.random.default_rng(106)
oml = rng.uniform(0.12, 0.66, 500)
oml = np.array([o for o in oml if min(abs(o - l) for l in LINES) > 0.01])
I_lab = I_hann(mids_sc, oml)
I_zero = I_hann(z_sc, oml)
print(f"\nP2' — düşük-ω (0.12-0.66, çizgi-dışı, Hann) medyan I:")
print(f"  lab orta-nokta: {np.median(I_lab):.3e}")
print(f"  lab sıfır-ızgarası: {np.median(I_zero):.3e}")
print(f"  GERÇEK orta-nokta (102a, ω<1): 4.12e-06 → lab/gerçek ≈ "
      f"{np.median(I_lab)/4.12e-6:.0f}×", flush=True)

# ---- P3': lab vuruş taraması (101i usulü, ±0.004 çözünürlük)
def D_omega(z, om0, halfw=0.004, nf=60):
    gaps = np.diff(z)
    mids = 0.5 * (z[:-1] + z[1:])
    ds = gaps * np.log(mids / TWO_PI) / TWO_PI - 1
    tt = (mids - mids.mean()) / (mids[-1] - mids[0])
    P = np.vstack([np.ones_like(tt), tt, tt**2, tt**3]).T
    ds = ds - P @ np.linalg.lstsq(P, ds, rcond=None)[0]
    oms = om0 + np.linspace(-halfw, halfw, nf)
    num = 0.0 + 0j; den = 0.0
    for s0 in range(0, len(oms), 30):
        ob = oms[s0:s0 + 30]
        rr = np.exp(1j * np.outer(ob, z)).sum(axis=1)
        GG = (np.exp(1j * np.outer(ob, mids)) * ds[None, :]).sum(axis=1)
        num += (GG * np.conj(rr)).sum()
        den += (np.abs(rr)**2).sum()
    tau0 = om0 / L0
    kap = 2 * np.pi * tau0
    return abs(num) / ((2 * np.sin(kap / 2) / kap) * den)

print("\nP3' — lab sabit-ω taraması (±0.004):")
NOKTA = [("vuruş 3-2", 0.4055), ("kontrol", 0.3700), ("kontrol", 0.4400),
         ("vuruş 7-5", 0.3365), ("vuruş 5-3", 0.5108), ("kontrol", 0.4750)]
for ad, om0 in NOKTA:
    print(f"  {ad:>10} ω={om0:.4f}: D = {D_omega(z_sc, om0):.3f}", flush=True)
