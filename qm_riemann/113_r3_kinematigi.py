"""
113 — R₃ KİNEMATİĞİ: HAKEM EĞRİLERİ (25 Ağustos, gece)
==========================================================================
112 sürprizi: boyalı referans R₃ = −4.4 (naif +1 değil). Kalem dökümü
üç aday terim veriyor; burada her birini ifşa edecek hakem eğrileri:

ÖN-MÜHÜRLÜ KALEM ÖNGÖRÜLERİ:
  K1  GAUSSIAN jitter + boyalı dalga: saf dilatasyonda ⟨η³⟩ ≡ 0
      (Wick) → σ³-normlu R₃ ≈ 0. Sapma = katlanma (g≥0 duvarı) +
      komşu-lag terimi (⟨η_n³·δm_{n±1}⟩ = 3σ⁴ ≠ 0, Isserlis) —
      ikisi de ωσ ile büyümeli.
  K2  ÇARPIK jitter (gamma, skew≈0.24) + boyalı: saf dilatasyon
      R₃(m₃-norm) = (κ/2)/sin(κ/2) ≈ +1.04..+1.16. Ölçülen bundan
      saparsa fark = kübik karışım terimleri.
  K3  ζ+boyalı τ-eğrisi (112'nin 2 noktası → 5 nokta): biçimi K1/K2
      tabanlarıyla kıyaslanınca baskın terim seçilir.
Not: gamma ile skew 0.24 seçilince kurt 0.087 kalır (gerçek η kurt
0.86 — tek-parametreli tür; sınırlama dürüstçe kayıtlı).

SONUÇ (25 Ağustos gecesi) — ÜÇ AMPİRİK YASA + BİR YAPISAL KEŞİF:
  K1  Gaussian kinematiği SIFIR DEĞİL: R₃(σ³) ≈ −0.3·κ³·σ_ds
      (σ'da lineer ✓, τ'da ~kübik; skew(η)≈0 kalıyor → katlanma değil,
      komşu-lag mekanizması [⟨η³δm⟩=3σ⁴] aday). Kalem hedefi 1.
  K2  YAPISAL KEŞİF — konum-jitter FARKI gap-çarpıklığını ÖLDÜRÜR
      (⟨(Δw)³⟩ ≡ 0, her dağılımda): m₃-norm patladı (+46, anlamsız),
      skew(η)≈0. GERÇEK η'nın +0.24 çarpıklığı GAP-DÜZEYİ asimetri
      (itme fiziği) — jitter sınıfının tamamen DIŞINDA. Bu, kalan-%25
      için en somut iz: gerçek gürültünün çift-fonksiyon (repulsion)
      asimetrisi hiçbir konum-jitter modelinde yok.
  K3  ζ+boyalı eğrisi: R₃(σ³) = −2.71/−2.88/−3.18/−2.61/−2.33
      (τ* 0.10-0.40) — yaklaşık DÜZ; K1-Gaussian öngörüsünün (−0.02..
      −0.95) çok üstünde → baskın terim gerçek gürültünün kendi
      yapısıyla dilatasyonun etkileşimi. İşaret ANTİ-ölçek-sürükleme:
      naif "m₃ ∝ σ³ sürüklenir" +2.1 verirdi (σ³-norm), ölçüm −2.8 —
      büyüklük ~1.3×, işaret ters. Kalem hedefi 2.
  Gerçek kendi-dalga R₃ (σ³-norm −1.0..−1.45) / boyalı (−2.7..−3.2)
  ≈ 0.4 — 112'nin "kinematik yanıtın ~0.4'e bastırılması" okuması
  eğri düzeyinde doğrulandı.
  KALEM OTURUMU HEDEFLERİ (dinç kafayla): (1) −0.3κ³σ'nın Isserlis
  türetimi; (2) anti-ölçek-sürüklemenin kübik muhasebesi (gerçek
  gap-dağılımının biçim-değişimi); (3) gap-düzeyi asimetrili gürültü
  sınıfıyla D-merdiveni (kalan-%25'in yeni adayı).
"""

import numpy as np
from pathlib import Path
from sympy import primerange

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
rng = np.random.default_rng(113)
L0 = 10.37
NZ = 200000
gbar = TWO_PI / L0

def pk_list(lim):
    out = []
    for p in primerange(2, int(lim) + 1):
        q = p
        while q <= lim:
            out.append(q); q *= p
    return sorted(set(out))

def rvm_N(t):
    x = t / TWO_PI
    return x * np.log(x / np.e) + 7 / 8

t0 = TWO_PI * np.exp(L0)
idx = np.arange(NZ + 1, dtype=float)
xg = t0 + idx * gbar
for _ in range(8):
    xg = xg - (rvm_N(xg) - rvm_N(t0) - idx) / (np.log(xg / TWO_PI) / TWO_PI)

A_DS = 0.10

def r3_olc(z, om):
    g = np.diff(z)
    m = 0.5 * (z[:-1] + z[1:])
    ds = g / gbar - 1
    ds = ds - ds.mean()
    X = np.vstack([np.ones_like(m), np.cos(om * m), np.sin(om * m),
                   np.cos(2 * om * m), np.sin(2 * om * m)]).T
    b, *_ = np.linalg.lstsq(X, ds, rcond=None)
    A1 = np.hypot(b[1], b[2]); ph = np.arctan2(b[2], b[1])
    eta = ds - X @ b
    s2 = float((eta**2).mean()); m3 = float((eta**3).mean())
    y3 = eta**3 - m3
    b3, *_ = np.linalg.lstsq(X, y3, rcond=None)
    P3 = b3[1] * np.cos(ph) + b3[2] * np.sin(ph)
    r3_m3 = P3 / (3 * m3 * A1) if abs(m3) > 1e-8 else np.nan
    r3_s3 = P3 / (s2**1.5 * A1)
    return A1, r3_m3, r3_s3, m3 / s2**1.5

def jitter(dagilim, sig_ds, n):
    sw = sig_ds / np.sqrt(2) * gbar
    if dagilim == "gauss":
        return rng.normal(0, sw, n)
    k = (2 / 0.24)**2                       # gamma: skew 0.24
    raw = rng.gamma(k, 1.0, n)
    return (raw - k) / np.sqrt(k) * sw

print("K1 — GAUSSIAN + boyalı (σ³-norm; kalem: ≈0, sapma ∝ ωσ):")
print(f"{'σ_ds':>6} {'τ':>5} {'R₃(σ³)':>8} {'skew(η)':>8}")
for sig in [0.10, 0.148, 0.20]:
    for tau in [0.15, 0.30, 0.45]:
        om = tau * L0
        UA = A_DS / (2 * np.sin(np.pi * tau)) * gbar
        w = jitter("gauss", sig, NZ + 1)
        y = np.sort(xg + w)
        z = np.sort(y + UA * np.cos(om * y))
        A1, r3m, r3s, sk = r3_olc(z, om)
        print(f"{sig:>6.3f} {tau:>5.2f} {r3s:>+8.4f} {sk:>+8.3f}", flush=True)

print("\nK2 — ÇARPIK (gamma, skew 0.24) + boyalı "
      "(kalem: R₃ = (κ/2)/sin(κ/2) ≈ +1.0..+1.2):")
print(f"{'τ':>5} {'R₃(m₃)':>8} {'kalem':>6} {'R₃(σ³)':>8} {'skew(η)':>8}")
for tau in [0.10, 0.15, 0.30, 0.45]:
    om = tau * L0
    UA = A_DS / (2 * np.sin(np.pi * tau)) * gbar
    w = jitter("gamma", 0.148, NZ + 1)
    y = np.sort(xg + w)
    z = np.sort(y + UA * np.cos(om * y))
    A1, r3m, r3s, sk = r3_olc(z, om)
    kalem = (np.pi * tau) / np.sin(np.pi * tau)
    print(f"{tau:>5.2f} {r3m:>+8.3f} {kalem:>6.3f} {r3s:>+8.4f} {sk:>+8.3f}",
          flush=True)

print("\nK3 — ζ + boyalı τ-eğrisi (112'nin 2 noktası → 5):")
d41 = np.load(HERE / "41_bigT_windows.npz")
gz, tm = d41["gaps_120k"], d41["tmid_120k"]
zz = np.empty(len(gz) + 1)
zz[0] = tm[0] - gz[0] / 2
zz[1:] = zz[0] + np.cumsum(gz)
Lz = float(np.log(tm / TWO_PI).mean())
gz_bar = TWO_PI / Lz
qs = pk_list(min(np.exp(0.52 * Lz), 720))
LINES = [np.log(q) for q in qs]

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

print(f"{'τ*':>5} {'R₃(m₃)':>8} {'R₃(σ³)':>8}")
for taus in [0.10, 0.15, 0.22, 0.30, 0.40]:
    om = taus * Lz
    while min(abs(om - l) for l in LINES) < 0.012:
        om += 0.013
    U = 0.1 / (2 * np.sin(np.pi * om / Lz))
    zp = zz + U * gz_bar * np.cos(om * zz)
    gp = np.diff(zp)
    mp_ = 0.5 * (zp[:-1] + zp[1:])
    ds = gp * np.log(mp_ / TWO_PI) / TWO_PI - 1
    allf = [np.log(q) for q in qs] + [om]
    b1, fit1 = chunked_fit(ds, mp_, allf)
    eta = ds - fit1
    s2 = float((eta**2).mean()); m3 = float((eta**3).mean())
    y3 = eta**3 - m3
    b3, _ = chunked_fit(y3, mp_, allf)
    i = allf.index(om)
    cA, sA = b1[3 + 2 * i], b1[3 + 2 * i + 1]
    A1 = np.hypot(cA, sA); ph = np.arctan2(sA, cA)
    c3, s3 = b3[3 + 2 * i], b3[3 + 2 * i + 1]
    P3 = c3 * np.cos(ph) + s3 * np.sin(ph)
    print(f"{taus:>5.2f} {P3/(3*m3*A1):>+8.3f} {P3/(s2**1.5*A1):>+8.4f}",
          flush=True)
