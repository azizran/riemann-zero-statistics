"""
148 — BANT-ÇÖZÜNÜRLÜKLÜ BOND SPEKTROSKOPİSİ: γ₁(τ) HARİTASI (30 Ağu)
==========================================================================
Gram-tam tasarım: aynı η, aynı çizgi, iki izdüşüm — ζ (η_n, m_n'de) ve
ζ₊ (η_{n+1}, m_n'e göre); bilinen faz ilerlemesi 2πτ ile döndür:
  γ₁(Q) = Re[ζ₊ · conj(ζ) · e^{i2πτ_Q}] / |ζ|²   (bant-toplu, ara-nokta
  referans çıkarımlı). Pay/payda aynı konvansiyon → atıf kancası oranda
  sadeleşmeli (Γ2 bunu sınar).
ÖN-MÜHÜR:
  Γ1  γ₁(τ) haritası; toplulaştırılmış γ̄₁ ≈ 0.77±0.10 (144'ün c₁
      açığı kapanır).
  Γ2  GRAM-TAMLIK: γ₁ oranı taban-bağımsız (0.40 vs 0.52 tabanı,
      örtüşen bantlarda ±%10) — ilk konvansiyonsuz değişmez adayı.
  Γ3  C₀ ZİNCİRİ: γ₁ = e^{−λ(τ)} ise nefes kanalı λ'(τ) = 2π·R_bond
      = 1.80±0.4 öngörür → ölçülen log-eğim bunu vurursa C₀ mekanizması
      MÜHÜR; γ₁ düz ise C₀ yeniden açık.
  Bonus: Im-kısmı ≈ 0 (faz-doğruluk bond adımında).

SONUÇ (30 Ağustos, gerçek koşudan) — SKALER γ₁ DÜŞTÜ; YENİ NESNE:
BOND-ADIMI DİSPERSİYONU:
  Transfer KOMPLEKS çıktı: |γ| = 0.74/0.64/0.58/0.87 (yüksek koherans)
  ama 2πτ çıkarıldıktan sonra ANORMAL FAZ kalıyor: φ_ekstra =
  0.73/1.39/1.83(/1.69 gürültülü) rad → fit φ(τ) ≈ 7.1·τ − 3.4.
  Gerçek kısmın 0.55→−0.15 düşüşü koherans kaybı değil, fazın π/2
  geçişi — Γ1 (0.77) ve Γ3 (eğim 1.8) BU HALİYLE HÜKÜMSÜZ/DÜŞTÜ.
  Γ2 KISMİ: kompleks değerler iki tabanda ilk üç bantta ±%10-25
  kararlı — oran gölgelerden sağlam ama tam değişmez değil; doğru dil
  TRANSFER-MATRİSİ (çizgi-uzayı, Gram-tam).
  YENİ YASA ADAYI: kuyruk gürültüsünün gap-geçiş dispersiyonu
  Θ(τ) = 2πτ + φ(τ); c₁ defteri cos(Θ) ile yeniden kurulmalı (144'ün
  0.77'sinin asıl adresi bu olabilir). C₀ YENİ ADAY: dispersiyon-nefesi
  δΘ = τ(2π+φ′)·ds ≈ 13.4τ·ds — kaba tahmin O(−1) katkı: SIRADAKİ
  KALEM. Im≠0: faz-doğruluk teoremi yalnız EŞZAMANLI okuma içindi;
  bond-adımı (zaman-öteleme) dispersif — çelişki yok, yeni fizik.
"""

import numpy as np
from pathlib import Path
from sympy import primerange, factorint

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
d = np.load(HERE / "128_odl_zeros6_2e6_zeros.npz")
Z = np.sort(np.asarray(d["zeros"], dtype=float))
zz = Z[len(Z) - 300000:]

def pk(lim):
    out = []
    for p in primerange(2, lim + 1):
        q = p
        while q <= lim:
            out.append(q); q *= p
    return sorted(set(out))

def fit_eta(taumax):
    g = np.diff(zz); mid = 0.5 * (zz[:-1] + zz[1:])
    Lw = np.log(mid / TWO_PI); L = float(Lw.mean())
    ds = g * Lw / TWO_PI - 1
    Nn = len(ds)
    tt = (mid - mid.mean()) / (mid[-1] - mid[0])
    qs = [q for q in pk(min(int(np.exp(taumax * L)), 100000))
          if np.log(q) / L <= taumax]
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
    return L, mid, eta

def gamma_bant(L, mid, eta, lo, hi, rng):
    e0, e1 = eta[:-1], eta[1:]
    m0 = mid[:-1]
    allq = pk(int(np.exp(min(hi + 0.02, 0.86) * L)))
    allw = np.array([np.log(q) for q in allq])
    T = mid[-1] - mid[0]; dres = TWO_PI / T
    cand = [(q, w) for q, w in zip(allq, allw) if lo < w / L <= hi]
    if len(cand) > 220:
        idx = rng.choice(len(cand), 220, replace=False)
        cand = [cand[i] for i in idx]
    NUM = 0.0; IMM = 0.0; DEN = 0.0
    NUMo = 0.0; DENo = 0.0; kul = 0
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
            cross = zp * np.conj(z) * np.exp(1j * TWO_PI * W / L)
            if hedef:
                NUM += cross.real; IMM += cross.imag
                DEN += abs(z)**2
            else:
                NUMo += cross.real; DENo += abs(z)**2
        kul += 1
    g1 = (NUM - NUMo) / (DEN - DENo)
    im = (IMM) / (DEN - DENo)
    return g1, im, kul, DEN - DENo

rng = np.random.default_rng(21)
BANTLAR = [(0.55, 0.62), (0.62, 0.70), (0.70, 0.78), (0.78, 0.85)]
SON = {}
for tb in (0.52, 0.40):
    L, mid, eta = fit_eta(tb)
    print(f"taban τ≤{tb}  (L={L:.3f}):", flush=True)
    for lo, hi in BANTLAR:
        g1, im, kul, den = gamma_bant(L, mid, eta, lo, hi, rng)
        SON.setdefault(tb, []).append((0.5 * (lo + hi), g1))
        print(f"  τ∈({lo},{hi}]  çizgi={kul:3d}  γ₁={g1:+.3f}  "
              f"Im={im:+.3f}", flush=True)

print("\nΓ2 — taban-bağımsızlık (γ₁ oranları, örtüşen bantlar):")
for i, (t, _) in enumerate(SON[0.52]):
    print(f"  τ̄={t:.3f}:  0.52-taban {SON[0.52][i][1]:.3f}  "
          f"0.40-taban {SON[0.40][i][1]:.3f}")

print("\nΓ3 — λ(τ) = −ln γ₁ eğimi:")
t_ = np.array([t for t, _ in SON[0.52]])
g_ = np.array([g for _, g in SON[0.52]])
m = g_ > 0
lam = -np.log(g_[m])
A = np.polyfit(t_[m], lam, 1)
print(f"  λ(τ) ≈ {A[1]:+.3f} {A[0]:+.3f}·τ   [Γ3 öngörü eğim ≈ 1.80±0.4]")
print(f"  γ̄₁ (bant ort) = {np.mean(g_[m]):.3f}   [Γ1 hedef ≈ 0.77]")
