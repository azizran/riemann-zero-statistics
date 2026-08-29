"""
144 — KOLEKTİF K(τ): ÇEKİRDEK, ÜÇ-DALGA MAKİNESİNDEN (29 Ağu)
==========================================================================
Zincir: (1) standart 0.52-taban η'sını kur (K ölçümlerinin kendi
konvansiyonu); (2) η'nın kuyruk-soğurma profilini ρ(τ) bantlar halinde
ölç (çizgi-üstü/ara-nokta spektroskopisi, çözünür bölge); (3) kolektif
çekirdek oranını kur:
   κ_B = Σ ρ(τ)·w(τ)·πτcos(3πτ)/sin(πτ) / Σ ρ(τ)·w(τ)·cos(2πτ),
   w = a²sin²(πτ), toplam τ∈(0.52, 1] tam merdiven üstünde
   (135'in oranı, kuyruk artık ÇIPLAK değil — kolektif soğurmalı);
(4) paydayı ölçülen c₁=⟨η_nη_{n+1}⟩ ile çapraz-sağla; (5) kalıntı
C(τ_p) = K_öl(τ_p) − κ_B'nin şeklini üç-dalga C-kanalı imzası
q^{-1/2} ile karşılaştır (K_öl: 138 son-300k sütunu).
ÖN-MÜHÜR:
  Ö1  ρ-ağırlıklı payda iptal-vari değil ve ölçülen c₁'i ±%25 vurur.
  Ö2  κ_B, K bandına düşer (−2.8..−3.4) — 132a/136 değer sorununun
      kolektif şifası.
  Ö3  Kalıntı C(τ_p) küçük; sıfır değilse q^{-1/2} şekliyle korele
      (üç-dalga doğrudan kanalı imzası).

SONUÇ (29 Ağustos, koşu + çapraz-ayrım testi):
  ρ(τ) HARİTASI (yeni): dört bantta 0.381/0.266/0.159/0.080 →
  ρ(τ) = exp(2.62 − 6.20τ) [L=12.03'te üs = −u/2: yine kritik ½!
  τ-yasası mı u-yasası mı ikinci pencerede ayrışır — AÇIK].
  Ö1 KISMİ ✓: payda ŞİFALI (−7.5e-3, iptal yok — 135'in kesim-
  duyarlılığı kolektif ağırlıkla çözüldü); 2·Σρw·cos2πτ = −0.0150 vs
  c₁ = −0.0116 (%30 fazla — ikinci mertebe kaydı).
  Ö2 DÜŞTÜ: κ_B = −1.845 — K bandında değil; pozisyon-modülasyon
  kanalı tek başına çekirdeği vermiyor. (Not: 136'nın κ_DW = −1.850'si
  ile çakışması kaza değil — DW-Gauss ve ρ-üstel ağırlıklar etkin
  bantta benzer.)
  Ö3 KISMİ, SİRENİYLE: pencere-içi C(τ_p), q^{-1/2}'ye r=0.986 oturdu
  AMA çapraz-pencere tanıkları (eş-τ farklı-q üç çift) q-modelini ÜÇ
  kez yanlış işaretle reddetti; τ-modeli üçünü de vurdu (rms 0.041 vs
  0.046). HÜKÜM: K yalnız τ'nun fonksiyonu — 138 L-değişmezliğiyle
  tutarlı; q^{-1/2} tek-pencere taklidiydi (siren protokolü çapraz
  testle yakaladı).
  AMPİRİK YASALAR: K(τ) = −2.786 − 1.897τ (21 nokta, rms 0.041);
  C(τ) = K − κ_B = −0.94 − 1.90τ (türetilecek kanal); tutarlılık:
  κ_ad = −2−K = +0.79+1.90τ → ölçülen κ_ad eğrisini %0-8 ile vuruyor.
  AÇIK: ρ'nun τ/u kimliği (orta pencerede bant testi); C(τ) kanalının
  üç-dalga türetimi (τ-yalnız normalizasyonla); c₁'in %30'u.
"""

import numpy as np
from pathlib import Path
from sympy import primerange, factorint

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi

d = np.load(HERE / "128_odl_zeros6_2e6_zeros.npz")
Z = np.sort(np.asarray(d["zeros"], dtype=float))
zz = Z[len(Z) - 300000:]
g = np.diff(zz)
mid = 0.5 * (zz[:-1] + zz[1:])
Lw = np.log(mid / TWO_PI)
L = float(Lw.mean())
ds = g * Lw / TWO_PI - 1
Nn = len(ds)
tt = (mid - mid.mean()) / (mid[-1] - mid[0])

def pk(lim):
    out = []
    for p in primerange(2, lim + 1):
        q = p
        while q <= lim:
            out.append(q); q *= p
    return sorted(set(out))

# (1) standart 0.52-taban zinciri (K ölçümlerinin konvansiyonu)
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
print(f"L={L:.3f}  taban {len(qs)} çizgi (τ≤0.52)  "
      f"σ_η²={np.var(eta):.4f}  c₁=⟨ηη⟩={c1:+.5f}", flush=True)

# (2) η kuyruk-soğurma profili ρ(τ) — bant spektroskopisi
allq = pk(int(np.exp(0.85 * L)))
allw = np.array([np.log(q) for q in allq])
T = mid[-1] - mid[0]
dres = TWO_PI / T
BANTLAR = [(0.525, 0.60), (0.60, 0.70), (0.70, 0.78), (0.78, 0.84)]
rng = np.random.default_rng(11)
RHO = []
print("\nρ(τ) bantları (çizgi-üstü − ara-nokta / çıplak):")
for lo, hi in BANTLAR:
    cand = [(q, w) for q, w in zip(allq, allw) if lo < w / L <= hi]
    if len(cand) > 260:
        idx = rng.choice(len(cand), 260, replace=False)
        cand = [cand[i] for i in idx]
    on = off = onp = 0.0; kul = 0
    for q, w in cand:
        j = np.searchsorted(allw, w)
        koms = [allw[k] for k in (j - 1, j + 1)
                if 0 <= k < len(allw) and abs(allw[k] - w) > 1e-12]
        gap = min(abs(w - k) for k in koms)
        if gap < 2.5 * dres:
            continue
        cg = 2 * np.mean(eta * np.cos(w * mid))
        sg = 2 * np.mean(eta * np.sin(w * mid))
        (pp, kk), = factorint(q).items()
        a = 1 / (np.pi * kk * np.sqrt(q))
        on += cg * cg + sg * sg
        onp += (2 * a * np.sin(np.pi * w / L))**2
        wm = w + gap / 2
        cg = 2 * np.mean(eta * np.cos(wm * mid))
        sg = 2 * np.mean(eta * np.sin(wm * mid))
        off += cg * cg + sg * sg
        kul += 1
    rho = (on - off) / onp
    tmid_b = 0.5 * (lo + hi)
    RHO.append((tmid_b, rho))
    print(f"  τ∈({lo:.3f},{hi:.2f}]  çizgi={kul:3d}  ρ={rho:.3f}", flush=True)

# ρ(τ) modeli: ölçülen bantlara log-doğrusal fit
tb = np.array([t for t, _ in RHO]); rb = np.array([r for _, r in RHO])
mfit = rb > 0
cf = np.polyfit(tb[mfit], np.log(rb[mfit]), 1)
def rho_f(t):
    return np.exp(cf[1] + cf[0] * t)
print(f"  model: ρ(τ) = exp({cf[1]:.2f} {cf[0]:+.2f}·τ)")

# (3) kolektif çekirdek oranı, tam merdiven τ∈(0.52,1]
num = den = 0.0
for p in primerange(2, int(np.exp(L)) + 1):
    q, mm = p, 1
    while q <= int(np.exp(L)):
        t = np.log(q) / L
        if t > 0.52:
            a2 = 1 / (np.pi**2 * mm**2 * q)
            w = a2 * np.sin(np.pi * t)**2 * rho_f(t)
            num += w * np.pi * t * np.cos(3 * np.pi * t) / np.sin(np.pi * t)
            den += w * np.cos(2 * np.pi * t)
        q *= p; mm += 1
kB = num / den
# payda çapraz-sağlaması: c₁ =? 2·den (per-line ⟨ηη⟩ katkısı 2a²sin²cos2πτ·ρ)
print(f"\n(3) κ_B(kolektif) = {kB:+.3f}   [payda={den:+.2e}]")
print(f"(4) payda sağlaması: 2·Σρw·cos2πτ = {2*den:+.5f}  vs  "
      f"c₁(ölçülen) = {c1:+.5f}  (oran {2*den/c1:.2f})")

# (5) kalıntı C(τ_p) ve q^{-1/2} şekli (K_öl: 138 son-300k)
K138 = [(2, 0.058, -2.867), (3, 0.091, -2.939), (5, 0.134, -3.017),
        (7, 0.162, -3.107), (11, 0.199, -3.151), (13, 0.213, -3.181),
        (17, 0.236, -3.247)]
print(f"\n(5) kalıntı C(τ_p) = K_öl − κ_B  ve  q^(-1/2) şekli:")
Cs = []; xs = []
for p, tp, Km in K138:
    Cp = Km - kB
    Cs.append(Cp); xs.append(p**-0.5)
    print(f"  p={p:>2} τ={tp:.3f}  K_öl={Km:+.3f}  C={Cp:+.3f}  "
          f"q^-½={p**-0.5:.3f}")
Cs = np.array(Cs); xs = np.array(xs)
A = np.polyfit(xs, Cs, 1)
r = np.corrcoef(xs, Cs)[0, 1]
print(f"  fit: C ≈ {A[1]:+.3f} {A[0]:+.3f}·q^(-1/2)   korelasyon r={r:.3f}")
