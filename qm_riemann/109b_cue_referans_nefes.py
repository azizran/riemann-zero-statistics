"""
109b — NEFESİN HAKİKAT MATRİSİ: CUE / BOYALI / LAB REFERANSLARI (24 Ağu)
==========================================================================
109 bulgusu: ζ kendi aritmetik dalgalarına ANTİ-ADYABATİK nefesle yanıt
verir (R ≈ −0.85). KİNEMATİK TEOREM (ön-mühürden önce türetildi):
boyanmış yerdeğiştirme dalgası u(t) her süreçte gap'leri (1+u′) ile
çarpar → η'nın varyans dalgası = 2σ²·(ort dalga) → R = +1 ZORUNLU.
Dört hücreli hakikat matrisi:
  S1  CUE + boyalı dalga → R = +1 beklenir (alet + termal kontrol).
  S2  ζ + boyalı dalga (çizgi-dışı ω*) → R = +1 beklenir; öyleyse
      ters-nefes YALNIZ gazın kendi dalgalarına özgü (seçiciliğin
      ikinci-moment sureti — 108c'nin ikizi).
  S3  ζ + kendi dalgaları: R = −0.85 (109'da ölçüldü).
  S4  EF-lab (öz-tutarlı) + kendi dalgaları → lab boyalıdır → R = +1
      beklenir; öyleyse ters-nefes EF-ÖTESİ dinamik işaretidir (106
      ailesine yeni üye).
Hücreler böyle çıkarsa: gaz, kendi aritmetiğinin fazında gürültüsünü
AKTİF yeniden düzenler — termal/kinematik hiçbir referans bunu yapmaz.

SONUÇ (24 Ağustos) — MATRİS DOLDU, SEÇİCİLİK 2. MOMENTTE KANITLANDI:
  S1 CUE+boyalı:  R_p = +0.94 / +0.65 / +0.15  (τ = 0.10/0.30/0.45)
     → kinematik +1 küçük τ'da tutuyor; 2.-moment kanalının KENDİ
     τ-transferi var (1. momentin a_v≈1'inden farklı) — referans eğrisi.
  S2 ζ+boyalı:    R_p = +0.95 / +0.77 / +0.27 — CUE ile örtüşüyor:
     boyalı dalgaya ζ da termal/kinematik davranır.
  S3 ζ kendi:     R_p ≈ −0.85 DÜZ (109) — işaret TERS, ayrım devasa.
  S4 lab kendi:   R_p = +0.59/+0.28/+0.22 — POZİTİF: ters nefes lab'da
     YOK → EF-ötesi dinamik işareti (106 ailesine yeni üye).
  TRANSFER-DÜZELTMELİ OKUMA: R_true = R_meas/transfer(S1-S2 eğrisi) ≈
     −0.90..−1.24 → R_true ≈ −1: gaz kendi dalgasına ADYABATİĞİN TAM
     TERSİ, EŞİT GÜÇTE nefesle yanıt veriyor ("tam ters nefes").
  YORUM (hipotez): anti-faz gürültü dalgası, dinamik perdeleyen ortamın
     imzası — sıkıştırılan bölgede artan gürültü, koheran okumadan
     (1−D)'yi yiyen karşı-yanıt olabilir; D≈1−0.36τ ile nicel köprü
     sıradaki teori ödevi.
"""

import numpy as np
from pathlib import Path
from sympy import primerange, factorint

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
rng = np.random.default_rng(1099)

def pk_list(lim):
    out = []
    for p in primerange(2, int(lim) + 1):
        q = p
        while q <= lim:
            out.append(q); q *= p
    return sorted(set(out))

def fit_ve_R(ds, mids, freqs, hedef, tt_norm=None, chunk=40000):
    """Tam-taban fit → A1, sonra η²/η-komşu dalgalarının R'leri (hedefte)."""
    if tt_norm is None:
        tt = (mids - mids.mean()) / (mids[-1] - mids[0])
    else:
        tt = tt_norm
    def cols(x, tloc, sl):
        c = [np.ones(sl.stop - sl.start), tloc[sl], tloc[sl]**2]
        for f in freqs:
            arg = f * x[sl]
            c += [np.cos(arg), np.sin(arg)]
        return np.vstack(c).T
    n = len(ds)
    C = 3 + 2 * len(freqs)
    XtX = np.zeros((C, C)); Xty = np.zeros(C)
    for s0 in range(0, n, chunk):
        sl = slice(s0, min(s0 + chunk, n))
        Xc = cols(mids, tt, sl)
        XtX += Xc.T @ Xc; Xty += Xc.T @ ds[sl]
    b1 = np.linalg.solve(XtX, Xty)
    eta = np.empty(n)
    for s0 in range(0, n, chunk):
        sl = slice(s0, min(s0 + chunk, n))
        eta[sl] = ds[sl] - cols(mids, tt, sl) @ b1
    s2 = float((eta**2).mean())
    def kanal(y, x, tloc):
        XtX2 = np.zeros((C, C)); Xty2 = np.zeros(C)
        for s0 in range(0, len(y), chunk):
            sl = slice(s0, min(s0 + chunk, len(y)))
            Xc = cols(x, tloc, sl)
            XtX2 += Xc.T @ Xc; Xty2 += Xc.T @ y[sl]
        return np.linalg.solve(XtX2, Xty2)
    b2 = kanal(eta**2 - s2, mids, tt)
    ee = eta[:-1] * eta[1:]
    c1 = float(ee.mean())
    mm = 0.5 * (mids[:-1] + mids[1:])
    tt2 = (mm - mm.mean()) / (mm[-1] - mm[0]) if tt_norm is None else tt[:-1]
    b3 = kanal(ee - c1, mm, tt2)
    out = {}
    for f in hedef:
        i = freqs.index(f)
        cA, sA = b1[3 + 2 * i], b1[3 + 2 * i + 1]
        A1 = np.hypot(cA, sA); ph = np.arctan2(sA, cA)
        P2 = b2[3 + 2 * i] * np.cos(ph) + b2[3 + 2 * i + 1] * np.sin(ph)
        P3 = b3[3 + 2 * i] * np.cos(ph) + b3[3 + 2 * i + 1] * np.sin(ph)
        out[f] = (A1, P2 / (2 * s2 * A1), P3 / (2 * c1 * A1))
    return out, s2, c1

# ---- S1: CUE + boyalı dalga
print("S1 — CUE + boyalı (N=256, M=600):", flush=True)
N, M = 256, 600
for mmode in [26, 77, 115]:
    tau = mmode / N
    U = 0.1 / (2 * np.sin(np.pi * tau))          # ds-genliği 0.1 hedefi
    gbar = TWO_PI / N
    DS, MD = [], []
    for s in range(M):
        A = (rng.normal(size=(N, N)) + 1j * rng.normal(size=(N, N))) / np.sqrt(2)
        Q, Rq = np.linalg.qr(A)
        Q = Q * (np.diagonal(Rq) / np.abs(np.diagonal(Rq)))
        th = np.sort(np.angle(np.linalg.eigvals(Q)))
        th = th + U * gbar * np.cos(mmode * th)
        th = np.sort(th)
        dth = np.diff(np.concatenate([th, [th[0] + TWO_PI]]))
        mid = (th + dth / 2)
        DS.append(dth * N / TWO_PI - 1)
        MD.append(mid)
    ds = np.concatenate(DS); mid = np.concatenate(MD)
    tt0 = np.linspace(-0.5, 0.5, len(mid))        # zararsız nüisans rampa
    out, s2, c1 = fit_ve_R(ds, mid, [float(mmode), float(2 * mmode),
                                     float(mmode + 41)], [float(mmode)],
                           tt_norm=tt0)
    A1, Rp, Rnn = out[float(mmode)]
    print(f"  τ={tau:.3f}: A1={A1:.4f} (hedef ~0.1)  R_p={Rp:+.3f}  "
          f"R_nn={Rnn:+.3f}  [σ²={s2:.3f}, c₁/σ²={c1/s2:+.3f}]", flush=True)

# ---- S2: ζ + boyalı dalga (çizgi-dışı)
print("\nS2 — ζ-120k + boyalı (çizgi-dışı ω*):", flush=True)
d41 = np.load(HERE / "41_bigT_windows.npz")
gz, tm = d41["gaps_120k"], d41["tmid_120k"]
zz = np.empty(len(gz) + 1)
zz[0] = tm[0] - gz[0] / 2
zz[1:] = zz[0] + np.cumsum(gz)
L = float(np.log(tm / TWO_PI).mean())
gbar = TWO_PI / L
qs = pk_list(min(np.exp(0.52 * L), 720))
LINES = [np.log(q) for q in qs]
for taus in [0.15, 0.30, 0.45]:
    om = taus * L
    while min(abs(om - l) for l in LINES) < 0.012:
        om += 0.013
    U = 0.1 / (2 * np.sin(np.pi * om / L))
    zp = zz + U * gbar * np.cos(om * zz)
    gp = np.diff(zp)
    mp_ = 0.5 * (zp[:-1] + zp[1:])
    ds = gp * np.log(mp_ / TWO_PI) / TWO_PI - 1
    freqs = [np.log(q) for q in qs] + [om, om + 0.037]
    out, s2, c1 = fit_ve_R(ds, mp_, freqs, [om])
    A1, Rp, Rnn = out[om]
    print(f"  τ*={taus:.2f}: A1={A1:.4f}  R_p={Rp:+.3f}  R_nn={Rnn:+.3f}",
          flush=True)

# ---- S4: EF lab (öz-tutarlı) kendi dalgaları
print("\nS4 — EF lab (öz-tutarlı), kendi çizgileri:", flush=True)
L0, NZ = 10.37, 60000
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
u = np.zeros_like(t)
for _ in range(60):
    u = 0.5 * u + 0.5 * (-S_field(t + u) / rho)
zl = t + u
gl = np.diff(zl)
ml = 0.5 * (zl[:-1] + zl[1:])
dsl = gl * np.log(ml / TWO_PI) / TWO_PI - 1
qsl = pk_list(min(np.exp(0.52 * L0), 720))
freqs = [np.log(q) for q in qsl]
hedef = [np.log(p) for p in [3, 7, 13]]
out, s2, c1 = fit_ve_R(dsl, ml, freqs, hedef)
print(f"  [σ_η²={s2:.4f}, c₁/σ²={c1/s2:+.3f}]")
for p, f in zip([3, 7, 13], hedef):
    A1, Rp, Rnn = out[f]
    print(f"  p={p}: A1={A1:.4f}  R_p={Rp:+.3f}  R_nn={Rnn:+.3f}", flush=True)
print("\nS3 (109'dan): ζ kendi dalgaları R_p ≈ −0.85, R_nn ≈ −1.2..−1.9")
