"""
137 — SON İNCE SORU: ADYABATİK ARTIK, ZEROS6 HASSASİYETİYLE (28 Ağu)
==========================================================================
"−2 tam mı, ≈−2.15 mi?" Üç ayaklı muhasebe (zeros6 son-300k penceresi;
129: 300k ile 1M özdeş, hata ±0.03):
  A1  Gerçek çizgilerde κ_arit(p) = R_nn/cos(πτ_p), p=2..17.
  A2  BOYALI çizgi-dışı dalgalar (A_ds=0.1) → κ_ad(τ*) = R_nn/cos(πτ*)
      — gerçek gazın dış sese bond-yanıtı (adyabatik artık, doğrudan).
  A3  Kapalı-form κ_DW (135/136) bu pencerenin kendi σ_u'suyla
      (tarak-termometre, indekssiz).
ÖN-MÜHÜR:
  S1  κ_arit − κ_ad ≈ κ_DW (±%5) ise: "−2.00 = κ_DW + κ_ad" kapanır;
      −2 sihirli tamsayı DEĞİL, DW-toplam + küçük artıktır.
  S2  κ_ad ≈ 0 ve κ_arit ≈ −2.00 ≠ κ_DW(−2.14) ise: %7 gerilim —
      DW-modeli/kuvvet-katkıları inceltmesi gerekir (sonraki av).

SONUÇ (28 Ağustos) — ÜÇÜNCÜ KAPI; İNCE SORU CEVAPLANDI:
  İki mühürlü şık da düşmedi — gerçek kapı üçüncüsü:
  A1  κ_arit = −1.975 ± 0.04 (yedi çizgi, düz) → TAM −2 İLE UYUMLU
      (%1.3); "≈−2.15" kapısı KAPANDI.
  A2  κ_ad = +0.89/+1.03/+1.02/+1.34 (τ*=0.08-0.24) — adyabatik artık
      KÜÇÜK DEĞİL, TAM BOY (+~1.0, τ ile yükseliyor; lab 109b eğrisiyle
      uyumlu). GERÇEK VERİDE İLK DOĞRUDAN κ_ad ÖLÇÜMÜ.
  Defter: fark = κ_arit − κ_ad = −2.84 → −3.24 (ort −3.05) — lab'ın
  keskin-kesim çekirdeği −3.18'e %4; DW'nin −1.85'ine DEĞİL.
  ⇒ Gerçek gazın ayrışımı LAB'IN ORİJİNALİ: (−3.05) + (+1.05) = −2.00.
  136'nın değer-okuması (κ ≈ κ_DW + küçük artık) REDDEDİLDİ; DW
  Gauss-sönümü fazla agresif — gazın kendi titreşimi büyük ölçüde
  KOHERENT (dalgaların kendisi), kendi dalgasını dekorele edemez.
  Model uyarısı (ek tarama): kapalı-form oranın paydası ΣA²cos2πτ derin
  kuyrukta (τ→1.4, L=12) neredeyse sıfırlanıyor → oran kesime aşırı
  duyarlı (σ_u taraması −0.6..−2.65 arası savruluyor). Çekirdek modelinin
  L=12 inceltmesi AYRI AV; ölçüm sonucu bundan bağımsız.
  Fark'ın τ-kayması (−2.84→−3.24): çekirdek tam cos(πτ)-çarpanlı düz
  değil — δ(τ) ince yapısıyla tutarlı (132c notu).
"""

import numpy as np
from pathlib import Path
from sympy import primerange

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi

def pk_list(lim):
    out = []
    for p in primerange(2, int(lim) + 1):
        q = p
        while q <= lim:
            out.append(q); q *= p
    return sorted(set(out))

def chunked_fit(y, tmid, freqs, chunk=40000):
    tt = (tmid - tmid.mean()) / (tmid[-1] - tmid[0])
    C = 3 + 2 * len(freqs)
    XtX = np.zeros((C, C)); Xty = np.zeros(C)
    n = len(y)
    def cols(sl):
        c = [np.ones(sl.stop - sl.start), tt[sl], tt[sl]**2]
        for f in freqs:
            arg = f * tmid[sl]
            c += [np.cos(arg), np.sin(arg)]
        return np.vstack(c).T
    for s0 in range(0, n, chunk):
        sl = slice(s0, min(s0 + chunk, n))
        Xc = cols(sl)
        XtX += Xc.T @ Xc; Xty += Xc.T @ y[sl]
    b = np.linalg.solve(XtX, Xty)
    fit = np.empty(n)
    for s0 in range(0, n, chunk):
        sl = slice(s0, min(s0 + chunk, n))
        fit[sl] = cols(sl) @ b
    return b, fit

def coef(b, freqs, f):
    i = freqs.index(f)
    return b[3 + 2 * i], b[3 + 2 * i + 1]

def rnn_at(z, hedefler, ekstra=()):
    g = np.diff(z)
    m = 0.5 * (z[:-1] + z[1:])
    Lw = np.log(m / TWO_PI)
    L = float(Lw.mean())
    ds = g * Lw / TWO_PI - 1
    qs = pk_list(min(int(np.exp(0.52 * L)), 720))
    freqs = [np.log(q) for q in qs] + list(ekstra)
    b1, f1 = chunked_fit(ds, m, freqs)
    eta = ds - f1
    ee = eta[:-1] * eta[1:]
    c1 = float(ee.mean())
    mm = 0.5 * (m[:-1] + m[1:])
    b3, _ = chunked_fit(ee - c1, mm, freqs)
    out = {}
    for f in hedefler:
        tau = f / L
        cg, sg = coef(b1, freqs, f)
        A1 = np.hypot(cg, sg); ph = np.arctan2(sg, cg)
        Rn = (coef(b3, freqs, f)[0] * np.cos(ph) +
              coef(b3, freqs, f)[1] * np.sin(ph)) / (2 * c1 * A1)
        out[f] = (tau, Rn, Rn / np.cos(np.pi * tau))
    return out, L

d = np.load(HERE / "128_odl_zeros6_2e6_zeros.npz")
Z = np.sort(np.asarray(d["zeros"], dtype=float))
zz = Z[len(Z) - 300000:]
L0 = float(np.log(zz.mean() / TWO_PI))
gbar = TWO_PI / L0

def rvm_N(t):
    x = t / TWO_PI
    return x * np.log(x / np.e) + 7 / 8

# A3 hazırlığı: pencerenin kendi σ_u'su (tarak)
xw = rvm_N(zz)
comb = abs(np.exp(2j * np.pi * xw).mean())
su = np.sqrt(-2 * np.log(comb)) / TWO_PI
print(f"pencere: L={L0:.2f}  σ_u(tarak)={su:.4f}", flush=True)

# A1: gerçek çizgiler
PR = [2, 3, 5, 7, 11, 13, 17]
out, L = rnn_at(zz, [np.log(p) for p in PR])
print(f"\nA1 — gerçek çizgiler (κ_arit = R_nn/cosπτ):")
for p in PR:
    tau, Rn, kap = out[np.log(p)]
    print(f"  p={p:>2} τ={tau:.3f}  R_nn={Rn:+.3f}  κ_arit={kap:+.3f}",
          flush=True)

# A2: boyalı çizgi-dışı dalgalar
LINES = [np.log(q) for q in pk_list(720)]
print(f"\nA2 — boyalı (κ_ad = R_nn/cosπτ*):")
KADS = []
for ts in [0.08, 0.13, 0.18, 0.24]:
    om = ts * L
    while min(abs(om - l) for l in LINES) < 0.012:
        om += 0.013
    U = 0.1 / (2 * np.sin(np.pi * om / L))
    zp = np.sort(zz + U * gbar * np.cos(om * zz))
    outp, _ = rnn_at(zp, [om], ekstra=[om])
    tau, Rn, kap = outp[om]
    KADS.append((tau, kap))
    print(f"  τ*={tau:.3f}  R_nn={Rn:+.3f}  κ_ad={kap:+.3f}", flush=True)

# A3: kapalı-form κ_DW (bu pencerenin σ_u'suyla)
from sympy import factorint
sig_t = su * TWO_PI / L
QALL = pk_list(200000)
num = den = 0.0
for q in QALL:
    om = np.log(q)
    tau = om / L
    if tau <= 0.52:
        continue
    w = np.exp(-om**2 * sig_t**2)
    if w < 1e-8:
        break
    (pp, kk), = factorint(q).items()
    a = np.log(pp) / (np.pi * np.sqrt(q) * om)
    A2c = (2 * a * np.sin(np.pi * tau))**2 * w
    num += A2c * np.pi * tau * np.cos(3 * np.pi * tau) / np.sin(np.pi * tau)
    den += A2c * np.cos(2 * np.pi * tau)
kdw = num / den
print(f"\nA3 — kapalı form: κ_DW = {kdw:+.3f}")

print(f"\nDEFTER (S1 testi): κ_arit − κ_ad =? κ_DW")
ka_t = [t for t, k in KADS]; ka_v = [k for t, k in KADS]
for p in PR:
    tau, Rn, kap = out[np.log(p)]
    kad = float(np.interp(tau, ka_t, ka_v))
    print(f"  p={p:>2}: κ_arit={kap:+.3f}  κ_ad={kad:+.3f}  "
          f"fark={kap-kad:+.3f}  [κ_DW={kdw:+.3f}]", flush=True)
