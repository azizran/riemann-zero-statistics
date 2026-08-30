"""
149 — DİSPERSİYON-NEFESİ: KAYNAK AYRIMI + BOND DEFTERİ + K₀ KAPANIŞ DENEMESİ
==========================================================================
(30 Ağu) 148'in kompleks transferi Γ_rot(τ) üstüne üç hakem:
D1  KAYNAK: M_emp(A) = ⟨e^{−iA·dsΔ}⟩ doğrudan (dsΔ=(ds_n+ds_{n+1})/2,
    HAM ds). Im M_emp ≪ 148-Im (≈0.08 vs 0.6) ⇒ koşulsuz-kinematik RED,
    Gram-karışım kaynağı destek. (Tersi çıkarsa: adım çarpıklığı — basit
    kaynak.)
D2  BOND DEFTERİ: B(τ) = Re[e^{−iA}Γ_rot] (bant başına, aynı örneklem)
    ve S₀(bant) = ölçülü lag-0 içerik → c₁_rec = ½ΣS₀·B; ölçülen
    c₁'i ±%15 vurursa 144'ün 0.77 gizemi ÇÖZÜLDÜ (cos2πτ → B).
D3  K₀ KAPANIŞI (tamamen ölçülmüş girdilerle):
    κ_bond = Σ S₀·B·2πτcot(πτ) + Σ S₀·τ·dB/dτ  /  Σ S₀·B
    (genlik-mod + adım-nefesi; ρ nefes almaz [147]); −2.87'yi ±%15
    vurursa C₀ ÇÖZÜLDÜ. (dB/dτ bant-farklarından — kaba; dürüst hata.)

SONUÇ (30 Ağustos, gerçek koşudan) — D1 ✓ D2 ✓ D3 ✓: C₀ ÇÖZÜLDÜ.
  D1  KAYNAK MÜHRÜ: M_emp'in Im'i 0.017-0.042 (koşulsuz adım-çarpıklığı
      φ≈0.08'lik bile değil) vs 148'in 0.49-0.87'si — 10-20×.
      Koşulsuz kinematik RED ⇒ dispersiyon KOLEKTİF (Gram-karışım:
      ρ-gradyanlı asimetrik çizgi karışımı; işaret+büyüklük+trend uyumlu).
  D2  BOND DEFTERİ KAPANDI: c₁_rec = ΣS₀·B = −0.01276 vs ölçülen
      −0.01158 (oran 1.10; 144'ün 1.30'u çözüldü — bond içeriği
      cos2πτ değil B(τ) taşıyor).
  D3  ÇEKİRDEK KAPANIŞI: tamamen ölçülmüş girdilerle
      κ_bond = (genlik −1.50) + (dispersiyon-nefesi −1.25) = −2.755 —
      hedef K₀ = −2.87±0.03'e %4 (dB/dτ beş banttan, kaba; %4 bu
      hatanın içinde). C₀ = −1.03'ün kimliği: DİSPERSİYON-NEFESİ.
  ⇒ −2'NİN TAM ZİNCİRİ (fenomenolojik düzeyde) KAPANDI:
  −2 = κ_ad + K₀·(ekran-gölgesi); K₀ = genlik-mod + dispersiyon-nefesi
  (ölçülü bant yapısından, %4); varyans kanalı %3 (147); dispersiyon
  kolektif (D1); L-değişmezlik PNT+τ-regülatör (138); mikro-verteks
  üç-dalga yasası (143). KALAN TEORİK ADIM: Γ_rot(τ)'nin G-yasasından
  (transfer-matris) inşası — son analitik dişli.
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

# standart 0.52-taban η
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
print(f"L={L:.3f}  σ_η²={np.var(eta):.4f}  c₁={c1:+.5f}", flush=True)

# D1 — adımın karakteristik fonksiyonu (HAM ds)
dsA = 0.5 * (ds[:-1] + ds[1:])
dsA = dsA - dsA.mean()
sA2 = float(np.var(dsA))
k3 = float(np.mean(dsA**3))
print(f"\nD1 — adım: σΔ²={sA2:.4f}  κ₃Δ={k3:+.5f} "
      f"(çarpıklık {k3/sA2**1.5:+.3f})")
print("   τ̄      A     M_emp(Re,Im)      148-Γrot(Re,Im)")
G148 = {0.5375: None, 0.585: (0.550, 0.490), 0.660: (0.118, 0.628),
        0.740: (-0.148, 0.564), 0.815: (-0.103, 0.867)}
for tb, gr in G148.items():
    A = TWO_PI * tb
    M = np.mean(np.exp(-1j * A * dsA))
    s = f"  {tb:.4f}  {A:.3f}  ({M.real:+.3f},{M.imag:+.3f})"
    if gr:
        s += f"     ({gr[0]:+.3f},{gr[1]:+.3f})"
    print(s, flush=True)

# D2/D3 — bant spektroskopisi: S₀ (toplam ölçekli) + Γ_rot + B
allq = pk(int(np.exp(0.86 * L)))
allw = np.array([np.log(q) for q in allq])
T = mid[-1] - mid[0]; dres = TWO_PI / T
e0, e1 = eta[:-1], eta[1:]
m0 = mid[:-1]
rng = np.random.default_rng(21)
BANTLAR = [(0.525, 0.55), (0.55, 0.62), (0.62, 0.70),
           (0.70, 0.78), (0.78, 0.85)]
SATIR = []
print("\nD2 — bantlar (S₀ toplam-ölçekli; B = Re[e^{-iA}Γ_rot]):")
for lo, hi in BANTLAR:
    tumu = [(q, w) for q, w in zip(allq, allw) if lo < w / L <= hi]
    cand = tumu
    if len(cand) > 220:
        idx = rng.choice(len(cand), 220, replace=False)
        cand = [cand[i] for i in idx]
    on0 = off0 = 0.0
    crN = 0j; crO = 0j
    kul = 0
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
                on0 += abs(z)**2; crN += cr
            else:
                off0 += abs(z)**2; crO += cr
        kul += 1
    olcek = len(tumu) / kul
    S0 = 0.5 * (on0 - off0) * olcek          # bant lag-0 varyans içeriği
    Grot = (crN - crO) / (on0 - off0)
    tb = 0.5 * (lo + hi)
    A = TWO_PI * tb
    B = np.cos(A) * Grot.real + np.sin(A) * Grot.imag
    SATIR.append((tb, S0, B))
    print(f"  τ∈({lo},{hi}] N={len(tumu):4d} kul={kul:3d}  S₀={S0:.5f}  "
          f"Γrot=({Grot.real:+.3f},{Grot.imag:+.3f})  B={B:+.3f}",
          flush=True)

c1_rec = sum(S0 * B for _, S0, B in SATIR)
print(f"\nD2 hüküm: c₁_rec = ΣS₀·B = {c1_rec:+.5f}  vs  "
      f"c₁ = {c1:+.5f}  (oran {c1_rec/c1:.2f})")

# D3 — çekirdek kapanışı (tamamen ölçülmüş girdiler)
tbs = np.array([t for t, _, _ in SATIR])
S0s = np.array([s for _, s, _ in SATIR])
Bs = np.array([bb for _, _, bb in SATIR])
dB = np.gradient(Bs, tbs)
num_amp = float(np.sum(S0s * Bs * TWO_PI * tbs / np.tan(np.pi * tbs)))
num_nef = float(np.sum(S0s * tbs * dB))
den = float(np.sum(S0s * Bs))
print(f"\nD3 — κ_bond = (genlik {num_amp:+.5f} + nefes {num_nef:+.5f})"
      f" / {den:+.5f} = {(num_amp+num_nef)/den:+.3f}   [hedef −2.87]")
print(f"    ayrıştırma: genlik/den = {num_amp/den:+.3f}   "
      f"nefes/den = {num_nef/den:+.3f}")
