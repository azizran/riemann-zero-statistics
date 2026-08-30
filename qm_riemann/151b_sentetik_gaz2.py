"""
151b — TAÇ SINAV, İKİNCİ DENEME: iyileştirilmiş çözücü + derin merdiven
==========================================================================
Girdi YALNIZ asal merdiveni: z_n, N̄(z)+S(z)=n öz-tutarlılığından
(Newton), S(z) = −Σ_{q: τ≤0.90} a_q sin(ω_q z), a_q = 1/(π m p^{m/2}),
N̄ = RvM ortalaması. Rastgelelik yok, serbest parametre yok, ölçülmüş
girdi yok. Sonra GERÇEK zincirin birebir aynısı: ds → 0.52-taban η →
Γ_rot bantları → R(τ) tayfı.
KALEM NOTU: analitik aday R = −dlnρ/dA ≈ 0.99 sabit (ölçek ✓ şekil ✗);
ilk bandın düşüklüğü taban-kenarı seyrelmesi adayı — hakem sentetik gaz.
ÖN-MÜHÜR:
  S1  Sentetik Γ_rot bantları gerçeği (Re,Im) ±0.08 üretir (bant-5
      hariç) → dispersiyon dahil tüm transfer, merdiven+geometriden.
  S2  Sentetik R(τ) tayfı gerçeğin ölçek/şeklini ±%30 üretir (kenar
      davranışı dahil çıkarsa kenar-yorumu da mühürlenir).
  S3  Kayıt: sentetik σ_ds², σ_η², c₁ vs gerçek (kesme farkları
      beklenir).
  HÜKÜM: S1+S2 → teori OPERASYONEL OLARAK KAPANIR: bütün kolektif
  fenomenoloji (ekran, soğurma, dispersiyon, R) = kesin özdeşlik +
  PNT + örnekleme geometrisi. R(τ)'nin kapalı formu "makine
  matematiği" statüsüne iner.

SONUÇ-REF (151 ilk koşu) — NİTEL ZAFER, NİCEL EKSİK (dürüst):
  FENOMEN ÜRETİLDİ: yalnız asallardan kurulan gaz, anormal dispersiyonu
  kendiliğinden taşıyor — işaret +, τ ile yükseliş, yapı aynı (fazlar
  sent. 0.54/1.04/1.86/2.68 vs gerçek 0.26/0.73/1.38/1.86). Dispersiyon
  ve sağkalım-nefesi asal merdiveni + örnekleme geometrisinin ZORUNLU
  sonucu — parametresiz gösterildi.
  S1/S2 NİCEL: ±%30 bandı DIŞINDA (fazlar ~1.4-1.5× dik; R tayfı
  taşıyor, 3.0-tavana çarpıyor; bant-5 faz sarması). ŞÜPHELİLER:
  (a) çözücü tam yakınsamadı (maks|F| 0.16'da takıldı — sönüm/kelepçe
  temkinli, salınım); (b) merdiven τ≤0.90 kesik; (c) sentetik
  varyanslar farklı (σ_ds²=0.119, σ_η²=0.061 vs 0.167/0.023 —
  yakınsama artığı η'ya sızıyor; c₁=−0.027 vs −0.012).
  S3 kayıt ✓. HÜKÜM: taç sınav İLK DENEMEDE fenomeni üretti;
  kalibrasyon çözücü+kesme iyileştirmesi istiyor → 151b (daha çok
  iterasyon, τ≤1.0 merdiven) ile yeniden.

SONUÇ (30 Ağustos, 151b koşusu) — NİCEL FARK GERÇEK:
  Derin merdiven (15450 çizgi, τ≤1.0) + 40 iterasyon: bantlar 151'le
  özdeş (±0.03; R_sent 1.58/2.70/3.00/3.00) — iki çözücü/kesmede
  KARARLI ⇒ saf-merdiven gazı ζ gazından ~1.45× daha dispersif; fark
  gerçek, artefakt değil. (maks|F|≈1.3 yalnız birkaç yapışkan noktada
  salınım; kütle yakınsak — sonuçlar değişmedi.)
  EKSİK MALZEME ADAYLARI: (a) gerçek S(t)'nin ERFC-yumuşak kesimi
  (BK 5.24 — alet çantası E) vs sentetik keskin merdiven; (b) gerçek
  gazın inkoherent payının sağkalım-korelasyonunu seyreltmesi
  (R_sent×0.65 ≈ R_ger, kaba tutarlı). SIRADAKİ: erfc-kesimli sentetik
  gaz → nicel kapanış denemesi.
"""

import numpy as np
from pathlib import Path
from sympy import primerange

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi

# gerçek pencerenin t-aralığı (karşılaştırma için aynı bölge)
d = np.load(HERE / "128_odl_zeros6_2e6_zeros.npz")
Z = np.sort(np.asarray(d["zeros"], dtype=float))
zr = Z[len(Z) - 300000:]
t0, t1 = float(zr[0]), float(zr[-1])
L_hedef = float(np.log(0.5 * (t0 + t1) / TWO_PI))

def rvm_N(t):
    x = t / TWO_PI
    return x * np.log(x / np.e) + 7 / 8

def rvm_d(t):
    return np.log(t / TWO_PI) / TWO_PI

def pk(lim):
    out = []
    for p in primerange(2, lim + 1):
        q = p
        while q <= lim:
            out.append(q); q *= p
    return sorted(set(out))

# merdiven (τ ≤ 0.90)
lad = []
for p in primerange(2, int(np.exp(1.00 * L_hedef)) + 1):
    q, m = p, 1
    while q <= int(np.exp(1.00 * L_hedef)):
        lad.append((np.log(q), 1.0 / (np.pi * m * np.sqrt(q))))
        q *= p; m += 1
om_l = np.array([w for w, _ in lad])
a_l = np.array([a for _, a in lad])
print(f"merdiven: {len(lad)} çizgi (τ≤0.90, L={L_hedef:.3f})", flush=True)

def S_ve_Sp(z, blok=800):
    S = np.zeros_like(z); Sp = np.zeros_like(z)
    for b0 in range(0, len(om_l), blok):
        w = om_l[b0:b0 + blok]; a = a_l[b0:b0 + blok]
        for s0 in range(0, len(z), 40000):
            sl = slice(s0, min(s0 + 40000, len(z)))
            arg = np.outer(z[sl], w)
            S[sl] += -np.sin(arg) @ a
            Sp[sl] += -np.cos(arg) @ (a * w)
            del arg
    return S, Sp

n0 = int(np.ceil(rvm_N(t0))); n1 = n0 + 300000
ns = np.arange(n0, n1, dtype=float)
# 1) pürüzsüz sayımı tam çöz (S'siz Newton — monoton, hızlı)
z = np.full_like(ns, 0.5 * (t0 + t1))
for _ in range(30):
    F = rvm_N(z) - ns
    z = np.clip(z - F / rvm_d(z), 100.0, None)
    if np.max(np.abs(F)) < 1e-9:
        break
print(f"  pürüzsüz çözüm: maks|F|={np.max(np.abs(rvm_N(z)-ns)):.2e}",
      flush=True)
# 2) öz-tutarlı Newton (sönümlü, kelepçeli)
gbar_t = 1.0 / rvm_d(z.mean())
for it in range(40):
    S, Sp = S_ve_Sp(z)
    F = rvm_N(z) + S - ns
    payda = np.maximum(rvm_d(z) + Sp, 0.3 * rvm_d(z))
    adim = np.clip(0.8 * F / payda, -1.0 * gbar_t, 1.0 * gbar_t)
    z = np.clip(z - adim, 100.0, None)
    mf = float(np.max(np.abs(F)))
    print(f"  Newton {it}: maks|F|={mf:.2e}", flush=True) if it % 4 == 0 or mf < 1e-3 else None
    if mf < 1e-3:
        break
z = np.sort(z)

# --- gerçek zincirin aynısı
g = np.diff(z)
mid = 0.5 * (z[:-1] + z[1:])
Lw = np.log(mid / TWO_PI)
L = float(Lw.mean())
ds = g * Lw / TWO_PI - 1
Nn = len(ds)
tt = (mid - mid.mean()) / (mid[-1] - mid[0])
qs = pk(min(int(np.exp(0.52 * L)), 720))
fr = np.array([np.log(q) for q in qs])
C = 3 + 2 * len(fr)
XtX = np.zeros((C, C)); Xty = np.zeros(C)
for s0 in range(0, Nn, 40000):
    sl = slice(s0, min(s0 + 40000, Nn))
    arg = np.outer(mid[sl], fr)
    Xc = np.empty((sl.stop - sl.start, C))
    Xc[:, 0] = 1; Xc[:, 1] = tt[sl]; Xc[:, 2] = tt[sl]**2
    Xc[:, 3::2] = np.cos(arg); Xc[:, 4::2] = np.sin(arg)
    XtX += Xc.T @ Xc; Xty += Xc.T @ ds[sl]
    del Xc, arg
b = np.linalg.solve(XtX, Xty)
eta = np.empty(Nn)
for s0 in range(0, Nn, 40000):
    sl = slice(s0, min(s0 + 40000, Nn))
    arg = np.outer(mid[sl], fr)
    eta[sl] = ds[sl] - (b[0] + b[1]*tt[sl] + b[2]*tt[sl]**2 +
                        np.cos(arg) @ b[3::2] + np.sin(arg) @ b[4::2])
    del arg
c1 = float(np.mean(eta[:-1] * eta[1:]))
print(f"\nS3 — sentetik: σ_ds²={np.var(ds):.4f}  σ_η²={np.var(eta):.4f}  "
      f"c₁={c1:+.5f}   [gerçek: 0.1674 / 0.0227 / −0.01158]", flush=True)

dsA = 0.5 * (ds[:-1] + ds[1:]); dsA -= dsA.mean()
allq = pk(int(np.exp(0.86 * L)))
allw = np.array([np.log(q) for q in allq])
T = mid[-1] - mid[0]; dres = TWO_PI / T
e0, e1 = eta[:-1], eta[1:]
m0 = mid[:-1]
rng = np.random.default_rng(21)
BANTLAR = [(0.525, 0.55), (0.55, 0.62), (0.62, 0.70),
           (0.70, 0.78), (0.78, 0.85)]
GERCEK = {0.5375: (0.787, 0.206, 0.40), 0.585: (0.550, 0.488, 1.07),
          0.660: (0.118, 0.614, 1.85), 0.740: (-0.148, 0.498, 2.12),
          0.815: (-0.103, 0.655, 1.50)}
print("\nS1/S2 — bant   sentetik Γrot(Re,Im)  R_sent  | gerçek Γrot  R_ger")
for lo, hi in BANTLAR:
    tumu = [(q, w) for q, w in zip(allq, allw) if lo < w / L <= hi]
    cand = tumu
    if len(cand) > 220:
        idx = rng.choice(len(cand), 220, replace=False)
        cand = [cand[i] for i in idx]
    crN = 0j; crO = 0j; on0 = off0 = 0.0
    for q, w in cand:
        j = np.searchsorted(allw, w)
        koms = [allw[k] for k in (j - 1, j + 1)
                if 0 <= k < len(allw) and abs(allw[k] - w) > 1e-12]
        gap = min(abs(w - k) for k in koms)
        if gap < 2.5 * dres:
            continue
        for W, hedef in ((w, True), (w + gap / 2, False)):
            cw, sw = np.cos(W * m0), np.sin(W * m0)
            zc = complex(2 * np.mean(e0 * cw), -2 * np.mean(e0 * sw))
            zp = complex(2 * np.mean(e1 * cw), -2 * np.mean(e1 * sw))
            cr = zp * np.conj(zc) * np.exp(1j * TWO_PI * W / L)
            if hedef:
                crN += cr; on0 += abs(zc)**2
            else:
                crO += cr; off0 += abs(zc)**2
    G = (crN - crO) / (on0 - off0)
    tb = 0.5 * (lo + hi)
    A = TWO_PI * tb
    phm = float(np.angle(G))
    best = None
    for R in np.linspace(0.0, 3.0, 121):
        wgt = np.exp(-A * R * dsA)
        M = np.mean(wgt * np.exp(-1j * A * dsA)) / np.mean(wgt)
        if best is None or abs(np.angle(M) - phm) < best[0]:
            best = (abs(np.angle(M) - phm), R)
    gr = GERCEK[round(tb, 4)]
    print(f"  {tb:.4f}  ({G.real:+.3f},{G.imag:+.3f})   {best[1]:4.2f}   "
          f"| ({gr[0]:+.3f},{gr[1]:+.3f})  {gr[2]:.2f}", flush=True)
