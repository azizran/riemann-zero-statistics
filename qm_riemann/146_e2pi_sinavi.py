"""
146 — ρ = A·e^{−2πτ} SINAVI (29 Ağu)
==========================================================================
145: ρ(τ) L-değişmez, eğim 6.20±0.15 ≈ 2π adayı. Sınav: bant menzilini
iki tabanla genişlet (τ≤0.40 tabanı → (0.42,0.52] bantları da açılır),
iki pencerede eğimi ±0.08'e sıkıştır.
ÖN-MÜHÜR:
  E1  (0.52,0.84] eğimi = 2π ± 0.10 → e^{−2πτ} MÜHÜR (ρ = A·e^{−ωḡ}:
      sağkalım, tek ortalama aralıktaki faz ilerlemesiyle cezalı).
  E2  Taban-bağımsızlık: (0.55,0.70] bantları iki tabanda ±%10 aynı
      (ρ, r'nin aksine, konvansiyonsuz bir gözlemlenebilir mi?).
  E1b Model uzatması kayda: (0.42,0.52] bandında A·e^{−2πτ} ≈ 0.75-0.8
      öngörür — tutar mı kırılır mı (yalnız kayıt).
  E3  A katsayısı kaydı (ln A ≈ 2.62; aday yorum sonra).

SONUÇ (29 Ağustos, gerçek koşudan) — E1 MÜHÜRSÜZ, E2 DÜŞTÜ (ÖĞRETİCİ):
  E1  Eğim = 6.03 ± 0.43 — 2π ile tutarlı (0.6σ) ama MÜHÜR DEĞİL.
  E2  DÜŞTÜ — ρ TABAN-BAĞIMLI: (0.55,0.70] bantları 0.40-tabanda
      0.56/0.43, 0.52-tabanda 0.38/0.27 (~1.5×). Küçük tabanda kalan
      güçlü çizgiler (0.40-0.52) çarpımsal Gram-ortaklıklarıyla üst
      bantlara güç taşıyor. ρ da r gibi Gram-atfetme gölgesi;
      "mutlak soğurma eğrisi" diye konvansiyonsuz nesne YOK.
      Eğim de konvansiyona bağlı (0.52-tabanda ~6.0; 0.40-tabanda
      ~2.6) ⇒ e^{−2πτ} KAPALI-FORM ADAYI BU HALİYLE DÜŞTÜ.
  Geçerli kalanlar: 145'in L-değişmezlik mührü (aynı konvansiyon içinde
  pencereler çakışık — o karşılaştırma meşru); 144/145'in ρ-ağırlıklı
  hesapları (K ölçümleriyle AYNI 0.52-konvansiyonunda — tutarlı).
  DERS: tek-çizgi VE bant atıfları gölge; fizik ancak Gram-tam
  (konvansiyonsuz) nesnelerde: C(n) çukurları, toplam yasalar, kanal
  toplamları. Nihai teori Gram-tam dille yazılmalı.
"""

import numpy as np
from pathlib import Path
from sympy import primerange, factorint

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
d = np.load(HERE / "128_odl_zeros6_2e6_zeros.npz")
Z = np.sort(np.asarray(d["zeros"], dtype=float))
N = len(Z)
PEN = {"son": Z[N - 300000:], "orta": Z[850000:1150000]}

def pk(lim):
    out = []
    for p in primerange(2, lim + 1):
        q = p
        while q <= lim:
            out.append(q); q *= p
    return sorted(set(out))

def fit_eta(zz, taumax):
    g = np.diff(zz); mid = 0.5 * (zz[:-1] + zz[1:])
    Lw = np.log(mid / TWO_PI); L = float(Lw.mean())
    ds = g * Lw / TWO_PI - 1
    Nn = len(ds)
    tt = (mid - mid.mean()) / (mid[-1] - mid[0])
    qs = [q for q in pk(int(np.exp(taumax * L)) + 1)
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

def bant_rho(L, mid, eta, lo, hi, rng):
    allq = pk(int(np.exp(min(hi + 0.02, 0.86) * L)))
    allw = np.array([np.log(q) for q in allq])
    T = mid[-1] - mid[0]; dres = TWO_PI / T
    cand = [(q, w) for q, w in zip(allq, allw) if lo < w / L <= hi]
    if len(cand) > 240:
        idx = rng.choice(len(cand), 240, replace=False)
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
    return (on - off) / onp, kul

rng = np.random.default_rng(11)
SONS = {}
for ad, zz in PEN.items():
    print(f"=== {ad} ===", flush=True)
    for tb, bantlar in ((0.40, [(0.42, 0.47), (0.47, 0.52),
                                (0.55, 0.62), (0.62, 0.70)]),
                        (0.52, [(0.525, 0.60), (0.60, 0.70),
                                (0.70, 0.78), (0.78, 0.84)])):
        L, mid, eta = fit_eta(zz, tb)
        for lo, hi in bantlar:
            rho, kul = bant_rho(L, mid, eta, lo, hi, rng)
            tmid = 0.5 * (lo + hi)
            SONS.setdefault(ad, []).append((tb, tmid, rho))
            print(f"  taban≤{tb}: τ∈({lo:.3f},{hi:.2f}] çizgi={kul:3d} "
                  f"ρ={rho:.3f}", flush=True)

print("\nE1 — eğim fiti (0.52-0.84, iki pencere birlikte, 0.52-taban):")
ts, rs = [], []
for ad in SONS:
    for tb, t, r in SONS[ad]:
        if tb == 0.52 and r > 0:
            ts.append(t); rs.append(np.log(r))
A = np.polyfit(ts, rs, 1)
res = np.array(rs) - np.polyval(A, ts)
serr = float(np.sqrt(np.sum(res**2) / (len(ts) - 2) /
                     np.sum((np.array(ts) - np.mean(ts))**2)))
print(f"  eğim = {-A[0]:.3f} ± {serr:.3f}   [2π = 6.283]   lnA = {A[1]:.3f}")
