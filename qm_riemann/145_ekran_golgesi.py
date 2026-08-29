"""
145 — ÇEKİRDEK SABİT Mİ? "KAYMA = EKRANIN GÖLGESİ" SINAVI (29 Ağu)
==========================================================================
Kalem: R_nn tahmincisi paydada ÖLÇÜLEN A1'i (= r(τ)·çıplak; Gram-
atfetme, 142) kullanır; çekirdeği süren konum alanı kesin özdeşlikçe
ÇIPLAKTIR. ⇒ K_meas(τ) = K₀ / r(τ): kayma sabit çekirdeğin ekran
gölgesi olabilir. O zaman C = K₀ − κ_B SABİT (−1.0 civ.) — türetme
hedefi basitleşir; "δ(τ) ince yapısı" çözülür.
A) Üç 138 penceresinde 0.52-taban fit → çizgi başına r(τ_p) →
   K_öl(138)·r(kendi penceresi) 21 noktada sabitlik testi.
B) ρ kimliği: orta pencerede bant spektroskopisi; τ-yasası ⇒ son ile
   aynı τ'da çakışır; u-yasası (q^-½) ⇒ orta ~%16 yüksek.
ÖN-MÜHÜR:
  H1  K·r'nin τ-eğimi, K'nin eğiminin ≤⅓'üne düşer ve rms ≤0.03 →
      çekirdek SABİT K₀ (değer kaydedilir); kayma ekran gölgesi.
  H2  ρ bantları: çakışma → τ-yasası; +%10-20 sistematik → u-yasası.

SONUÇ (29 Ağustos, gerçek koşudan) — H1 ✓ H2 ✓, İKİ BÜYÜK MÜHÜR:
  H1  ÇEKİRDEK SABİT: K·r eğimi −1.896 → +0.320 (6×↓), sabit-etraf
      rms 0.131 → 0.048. 21 noktadan 20'si K₀ = −2.87 ± 0.03
      (tek aykırı: erken-p17, 138'de zaten şüpheliydi). "δ(τ) ince
      yapısı" ÇÖZÜLDÜ: kayma, sabit çekirdeğin Gram-ekran gölgesiydi
      (payda ölçülen A1=r·çıplak; süren alan çıplak). Türetme hedefi
      tek sayıya indi: K₀ = κ_B(−1.845) + C₀(−1.03).
  H2  ρ τ-YASASI: orta bantları 0.379/0.264/0.156/0.078 ↔ son
      0.381/0.266/0.159/0.080 — binde-birkaç çakışık; u-yasasının
      +%16'sı RED. ρ(τ) L-değişmez τ-regülatörü (138 mekanizmasının
      istediği sınıf) artık ölçülü. KAPALI-FORM ADAYI (mühürsüz):
      eğim 6.20±0.15 ≈ 2π ⇒ ρ ∝ e^{−2πτ} = e^{−ωḡ} — sağkalım,
      dalganın TEK ORTALAMA ARALIKTAKİ faz ilerlemesiyle cezalı.
  AÇIK: kalan muhasebe gerilimi — boyalı κ_ad yüksek τ'da öz-adyabatik
  defterden ~%15 hızlı yükseliyor (boyalı≠öz olabilir); C₀=−1.03'ün
  türetimi; ρ=Ae^{−2πτ}'nin A'sı ve mühürlü testi.
"""

import numpy as np
from pathlib import Path
from sympy import primerange, factorint

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
d = np.load(HERE / "128_odl_zeros6_2e6_zeros.npz")
Z = np.sort(np.asarray(d["zeros"], dtype=float))
N = len(Z)
PEN = {"erken": Z[200000:500000], "orta": Z[850000:1150000],
       "son": Z[N - 300000:]}
K138 = {
    "erken": [(2, -2.929), (3, -3.035), (5, -3.143), (7, -3.179),
              (11, -3.250), (13, -3.289), (17, -3.169)],
    "orta":  [(2, -2.867), (3, -2.957), (5, -3.034), (7, -3.061),
              (11, -3.196), (13, -3.245), (17, -3.265)],
    "son":   [(2, -2.867), (3, -2.939), (5, -3.017), (7, -3.107),
              (11, -3.151), (13, -3.181), (17, -3.247)],
}
PR = [2, 3, 5, 7, 11, 13, 17]

def pk(lim):
    out = []
    for p in primerange(2, lim + 1):
        q = p
        while q <= lim:
            out.append(q); q *= p
    return sorted(set(out))

def fit052(zz):
    g = np.diff(zz); mid = 0.5 * (zz[:-1] + zz[1:])
    Lw = np.log(mid / TWO_PI); L = float(Lw.mean())
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
    return qs, fr, b, L, mid, eta

# A) K·r sabitlik testi
tumK, tumKr, tumTau = [], [], []
ETA = {}
for ad, zz in PEN.items():
    qs, fr, b, L, mid, eta = fit052(zz)
    ETA[ad] = (L, mid, eta)
    print(f"{ad}: L={L:.3f}")
    for (p, Km) in K138[ad]:
        i = qs.index(p)
        A1 = float(np.hypot(b[3 + 2*i], b[4 + 2*i]))
        a = 1 / (np.pi * np.sqrt(p))
        tau = np.log(p) / L
        r = A1 / (2 * a * np.sin(np.pi * tau))
        tumK.append(Km); tumKr.append(Km * r); tumTau.append(tau)
        print(f"  p={p:>2} τ={tau:.3f}  r={r:.3f}  K={Km:+.3f}  "
              f"K·r={Km*r:+.3f}", flush=True)
K_, Kr_, T_ = map(np.array, (tumK, tumKr, tumTau))
for etik, Y in (("K", K_), ("K·r", Kr_)):
    A = np.polyfit(T_, Y, 1)
    rms = float(np.sqrt(np.mean((Y - np.polyval(A, T_))**2)))
    rms0 = float(np.std(Y))
    print(f"A) {etik}: eğim={A[0]:+.3f}  sabit-etraf rms={rms0:.3f}  "
          f"(doğru-etraf {rms:.3f})  ort={Y.mean():+.3f}")

# B) ρ kimliği — orta pencerede bantlar
print("\nB) ρ bantları (orta; son'un değerleri: 0.381/0.266/0.159/0.080):")
L, mid, eta = ETA["orta"]
allq = pk(int(np.exp(0.85 * L)))
allw = np.array([np.log(q) for q in allq])
T = mid[-1] - mid[0]; dres = TWO_PI / T
rng = np.random.default_rng(11)
for lo, hi in [(0.525, 0.60), (0.60, 0.70), (0.70, 0.78), (0.78, 0.84)]:
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
    print(f"  τ∈({lo:.3f},{hi:.2f}]  çizgi={kul:3d}  ρ_orta={(on-off)/onp:.3f}",
          flush=True)
