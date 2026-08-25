"""
106 — İKİ DALIN LABORATUVARI: PERDE = ÖZ-TUTARLILIK MI? (24 Ağustos)
==========================================================================
92'nin gizemi: asal dalı perdelenir (D 0.91→0.48), spontane dal CUE
gibi geçer (1.08-1.18). 102 sonrası yeni resim ("her şey çizgidir";
çizgi-dışı içerik = örnekleme yan bantları) iki dalın TEK mekanizmaya
inebileceğini söylüyor. Belirleyici deney: gürültüsüz açık-formül lab'ı.

KİLİT HİLE: öz-tutarlı örgü (ρ̄u = −S(t+u); 100'ün makinesi) ile NAİF
örgünün (u = −S(t)/ρ̄; katı fark alma) v-genlik ORANI boru-hattı
çarpanını iptal eder → kalibrasyonsuz perde ölçümü.

ÖN-MÜHÜRLÜ ÖNGÖRÜLER (koşudan önce):
  P0  naif lab v-genlikleri katı eğriye ((2/π)sin(πτ) biçimi) oturur.
  P1  ÖZ-TUTARLILIK = PERDE: D_lab(q) = v_sc(q)/v_naif(q) gerçek
      D(τ)'yu üretir (0.85-0.95 @ τ≈0.1 → 0.45-0.65 @ τ≈0.5;
      92-T1: 0.908/0.863/0.797/0.726/0.634/0.484, binler 0.05→0.55).
      Mekanizma adayı: öteki modların varyansı çizginin gap-kanalını
      DW-gibi perdeler (Kepler-tipi dalga dikleşmesi).
  P2  lab'ın spontane eğrisi (çizgi-dışı bantlar, 101d estimatörü)
      gerçeğinkini ~%15 içinde üretir: donuk @0.04, ~CUE @0.3-0.5
      (gerçek: 0.041 / 0.66 @0.1 / 0.956 @0.2 / 1.075 / 1.088 / 1.177).
      ÜRETİRSE: "ζ spontane modlara CUE gibi yanıt verir" örnekleme
      optiğidir (lab'da CUE fiziği YOK). ÜRETEMEZSE: spontane dal
      gerçek gaz fiziğidir. İkisi de büyük sonuç.

SONUÇ (24 Ağustos):
  P0/P1: log kanalı naif örgüde nan (671 gap-kesişmesi); temiz hüküm
    analitik-katı karşılaştırmasından: v_sc/katı = 0.92-1.01, τ'da DÜZ
    → ÖZ-TUTARLILIK PERDE ÜRETMİYOR (gerçek 0.91→0.48'e inmiyor).
  P2: τ≥0.15'te lab-öz spontane eğrisi GERÇEĞİ %2-6 İÇİNDE üretir
    (1.041/1.112/1.192 vs 1.075/1.088/1.177 @ 0.3/0.4/0.5) — CUE
    fiziği olmayan örgüden → CUE-uyumu ÖRNEKLEME OPTİĞİ.
    τ=0.04'te lab DONUK DEĞİL (0.66 vs gerçek 0.041) → bkz. 106b:
    bu fark lab'ın kusuru (sahte vuruş parlaklığı), donmanın açıklaması
    açık formülün kesinliği çıktı.
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

# pürüzsüz örgü (100'ün makinesi)
t0 = TWO_PI * np.exp(L0)
idx = np.arange(NZ + 1, dtype=float)
t = t0 + idx * TWO_PI / L0
for _ in range(8):
    t = t - (rvm_N(t) - rvm_N(t0) - idx) / (np.log(t / TWO_PI) / TWO_PI)
rho = np.log(t / TWO_PI) / TWO_PI

# NAİF örgü: u = −S(t)/ρ̄  (katı fark alma)
u_na = -S_field(t) / rho
z_na = t + u_na
# ÖZ-TUTARLI örgü: ρ̄u = −S(t+u)
u_sc = np.zeros_like(t)
for _ in range(60):
    u_sc = 0.5 * u_sc + 0.5 * (-S_field(t + u_sc) / rho)
z_sc = t + u_sc
print(f"lab: L={L0}, n={NZ}; σ_u naif {u_na.std():.4f}, öz-tut {u_sc.std():.4f}",
      flush=True)

def v_genlik(z, om):
    g = np.diff(z)
    m = 0.5 * (z[:-1] + z[1:])
    yg = np.log(g * np.log(m / TWO_PI) / TWO_PI)
    X = np.vstack([np.ones_like(m), np.cos(om * m), np.sin(om * m)]).T
    b, *_ = np.linalg.lstsq(X, yg, rcond=None)
    return np.hypot(b[1], b[2])

# ---- P0/P1: asal dalı — çizgi başına v, iki örgüde
PRC = [2, 3, 5, 7, 11, 13, 17, 23, 31, 41, 53, 79, 107, 149]
GERCEK = {0.05: 0.908, 0.15: 0.863, 0.25: 0.797, 0.35: 0.726,
          0.45: 0.634, 0.55: 0.484}   # 92-T1 binleri
print(f"\nP0/P1 — asal dalı (D_lab = v_sc/v_naif; katı = (2/π)sin(πτ)·Λ/(π√q lnq)·ω/ρ̄):")
print(f"{'p':>4} {'τ':>6} {'v_naif':>8} {'katı-ör':>8} {'v_sc':>8} {'D_lab':>6}")
DL = []
for p in PRC:
    om = np.log(p)
    tau = om / L0
    vn = v_genlik(z_na, om)
    vs = v_genlik(z_sc, om)
    A = LAM[p] / (np.pi * np.sqrt(p) * np.log(p)) / rho.mean()
    katix = A * om * (2 * np.sin(np.pi * tau) / (2 * np.pi * tau)) * (2 * np.pi * tau) / (2 * np.pi / L0) / L0
    # katı differencing genliği (göreli-gap): A·(2/ḡ)·sin(κ/2), κ=2πτ
    kati = A * (2 / (TWO_PI / L0)) * np.sin(np.pi * tau)
    DL.append((tau, vn, vs, vs / vn))
    print(f"{p:>4} {tau:>6.3f} {vn:>8.5f} {kati:>8.5f} {vs:>8.5f} {vs/vn:>6.3f}",
          flush=True)
print("\nbin karşılaştırması (gerçek 92-T1):")
DL = np.array(DL)
for c, dger in GERCEK.items():
    m = (DL[:, 0] >= c - 0.05) & (DL[:, 0] < c + 0.05)
    if m.sum():
        print(f"  τ≈{c:.2f}: D_lab = {DL[m, 3].mean():.3f}  vs gerçek {dger:.3f}")

# ---- P2: spontane dal — 101d estimatörü, iki örgüde
LINES_ALL = [np.log(q) for q in QS if q <= 200]
def D_spont(z, taus):
    gaps = np.diff(z)
    mids = 0.5 * (z[:-1] + z[1:])
    ds = gaps * np.log(mids / TWO_PI) / TWO_PI - 1
    tt = (mids - mids.mean()) / (mids[-1] - mids[0])
    P = np.vstack([np.ones_like(tt), tt, tt**2, tt**3]).T
    ds = ds - P @ np.linalg.lstsq(P, ds, rcond=None)[0]
    out = []
    for tau0 in taus:
        oms = tau0 * L0 + np.linspace(-0.02 * L0, 0.02 * L0, 160)
        oms = np.array([o for o in oms
                        if min(abs(o - l) for l in LINES_ALL) > 0.01])
        num = 0.0 + 0j; den = 0.0
        for s0 in range(0, len(oms), 40):
            ob = oms[s0:s0 + 40]
            rr = np.exp(1j * np.outer(ob, z)).sum(axis=1)
            GG = (np.exp(1j * np.outer(ob, mids)) * ds[None, :]).sum(axis=1)
            num += (GG * np.conj(rr)).sum()
            den += (np.abs(rr)**2).sum()
        kap = 2 * np.pi * tau0
        out.append(abs(num) / ((2 * np.sin(kap / 2) / kap) * den))
    return out

TAUS = [0.04, 0.068, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50]
GER_SP = [0.041, 0.34, 0.66, 0.93, 0.956, 1.075, 1.088, 1.177]
print("\nP2 — spontane dal (çizgi-dışı bantlar):", flush=True)
Dsc = D_spont(z_sc, TAUS)
Dna = D_spont(z_na, TAUS)
print(f"{'τ':>6} {'lab-öz':>7} {'lab-naif':>8} {'gerçek':>7}")
for i, tt0 in enumerate(TAUS):
    print(f"{tt0:>6.3f} {Dsc[i]:>7.3f} {Dna[i]:>8.3f} {GER_SP[i]:>7.3f}")
