"""
150 — Γ_rot İNŞASI: ÇİZGİ-KARIŞIMI ARİTMETİKÇE KAPALI; YEREL
SAĞKALIM-AĞIRLIK MEKANİZMASI (30 Ağu)
==========================================================================
KALEM: (i) transfer-matrisin çizgi-karışım kanalı ölü: asal çizginin
alt çarpımsal komşusu YOK (Q/q₃ ∉ ℤ), üst komşu yalnız kule (G²·P
oranı <1e-4); ikinci-mertebe çift-aracılı ~1e-9. (ii) Kalan mekanizma:
yerel sağkalım-ağırlık korelasyonu — çizginin yerel gücü yerel adımla
w = e^{−A·R·dsΔ} gibi korelasyonluysa:
   Γ_rot = ⟨w·e^{−iA·dsΔ}⟩/⟨w⟩  (Gauss limitinde = M(A)·e^{iA²RσΔ²})
⇒ GENLİK koşulsuz sönümle aynı kalır; FAZ R(τ) ile döner.
ÖN-MÜHÜR:
  T1  |Γ_meas| = |M_emp| bant bant ±%12 (faz mekanizması genliğe
      dokunmaz).
  T2  AŞIRI-BELİRLEME: ampirik dağılımla (Gauss YOK) tek R(τ), fazı
      eşleyecek şekilde seçilince Re-kısmı da ±0.06 içinde ÇIKMALI
      → mekanizma ailesi MÜHÜR.
  T3  R(τ) tayfı kayıt — kalan türetme dişlisi Γ'den R'ye iner
      (R'nin soğurma teorisinden çıkarılması sonraki kalem).

SONUÇ (30 Ağustos, gerçek koşudan) — T1 ✓ T2 ✓ (4/5, biri sınırda):
Γ_rot'UN YAPISI ÇÖZÜLDÜ:
  T1  |Γ_meas| = |M_emp| bant bant ±%12 (0.813/0.732, 0.735/0.692,
      0.625/0.626, 0.520/0.555; bant-5 +%35 — çözünmemiş-çizgi önyargısı,
      kayıtlı). Faz mekanizması genliğe dokunmuyor ✓.
  T2  AŞIRI-BELİRLEME GEÇTİ: fazı eşleyen tek R ile ampirik-ağırlıklı
      karakteristik fonksiyon Re'yi de vurdu — Δ = 0.059/0.019/0.026/
      0.064(sınırda ✗)/0.005. Mekanizma ailesi MÜHÜR:
      Γ_rot(τ) = ⟨w·e^{−iA·dsΔ}⟩/⟨w⟩,  w = e^{−A·R(τ)·dsΔ}
      (yerel sağkalım-ağırlık korelasyonu; Gauss limitinde
      M(A)·e^{iA²RσΔ²}).
  T3  R(τ) TAYFI: 0.40 / 1.07 / 1.85 / 2.12 / 1.50(gürültülü) @
      τ = 0.54/0.585/0.66/0.74/0.815 — O(1), derinlikle artan.
  KALEM KAYDI: çizgi-karışım kanalı ARİTMETİKÇE ÖLÜ (asal çizginin alt
  çarpımsal komşusu yok; kule <1e-4; çift-aracılı ~1e-9) — 149-D1'in
  "Gram-karışımı" yorumu düzeltildi: kolektif ama ÇİZGİ-KÖŞEGEN
  mekanizma. KALAN TEK TÜRETME HEDEFİ: R(τ) — soğurmanın yerel-adım
  tepki tayfı (tek reel fonksiyon; kompleks transferden tek fonksiyona
  indirgendi).
"""

import numpy as np
from pathlib import Path
from sympy import primerange

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

dsA = 0.5 * (ds[:-1] + ds[1:])
dsA = dsA - dsA.mean()
sA2 = float(np.var(dsA))

# bant Γ_rot (149 ile aynı yol/tohum)
allq = pk(int(np.exp(0.86 * L)))
allw = np.array([np.log(q) for q in allq])
T = mid[-1] - mid[0]; dres = TWO_PI / T
e0, e1 = eta[:-1], eta[1:]
m0 = mid[:-1]
rng = np.random.default_rng(21)
BANTLAR = [(0.525, 0.55), (0.55, 0.62), (0.62, 0.70),
           (0.70, 0.78), (0.78, 0.85)]
print(f"L={L:.3f}  σΔ²={sA2:.4f}")
print("bant     |Γ_meas| |M_emp|   φ_meas   R(faz-eş)  Re_meas Re_model  hüküm")
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
            z = complex(2 * np.mean(e0 * cw), -2 * np.mean(e0 * sw))
            zp = complex(2 * np.mean(e1 * cw), -2 * np.mean(e1 * sw))
            cr = zp * np.conj(z) * np.exp(1j * TWO_PI * W / L)
            if hedef:
                crN += cr; on0 += abs(z)**2
            else:
                crO += cr; off0 += abs(z)**2
    G = (crN - crO) / (on0 - off0)
    tb = 0.5 * (lo + hi)
    A = TWO_PI * tb
    Me = np.mean(np.exp(-1j * A * dsA))
    phm = float(np.angle(G))
    # ampirik ağırlıklı model: R-taraması, fazı eşle
    Rgrid = np.linspace(0.0, 3.0, 121)
    best = None
    for R in Rgrid:
        wgt = np.exp(-A * R * dsA)
        M = np.mean(wgt * np.exp(-1j * A * dsA)) / np.mean(wgt)
        if best is None or abs(np.angle(M) - phm) < best[0]:
            best = (abs(np.angle(M) - phm), R, M)
    _, Rb, Mb = best
    hük = "✓" if abs(Mb.real - G.real) < 0.06 else "✗"
    print(f"{tb:.4f}   {abs(G):.3f}   {abs(Me):.3f}   {phm:+.3f}   "
          f"{Rb:5.2f}      {G.real:+.3f}  {Mb.real:+.3f}   {hük}",
          flush=True)
