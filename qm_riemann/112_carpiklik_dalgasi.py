"""
112 — ÇARPIKLIK DALGASI: ÜÇÜNCÜ MOMENTİN FAZ-ÇÖZÜMÜ (25 Ağustos)
==========================================================================
111 hükmü: Gaussian sınıflar perdenin ~%70 tavanında; kalan pay
Gaussian-ötesi. İlk bakılacak yer: ⟨η³⟩(faz) — hiç ölçülmemiş kanal.
Adyabatik yerel-ölçekleme referansı: η = η⁰(1+ε) → ⟨η³⟩ dalga genliği
= 3⟨η³⟩·ε → R₃ ≡ P₃/(3⟨η³⟩A1) = +1. Boyalı kontrol aynı scriptte.

ÖN-MÜHÜR:
  S1  η'nın global skew/kurt'u ilk kez kayda girer (alt-Gauss merdiveni
      için gerekli; u için 73: kurt −0.75 idi — η için bilinmiyor).
  S2  Çarpıklık dalgası ÖLÇÜLEBİLİRSE (plasebo üstü): işaret ve büyüklük
      yeni gözlemlenebilir. Tanıma 3. momentte de sürüyorsa R₃ < 0
      (ters nefesin ailesi); boyalı kontrol +1 civarı beklenir.
  S3  Boyalı kontrol R₃: kinematik referans (dilatasyon → +1 küçük τ).
Veri: ζ 120k+200k; asallar 2..17; plasebo +0.037.

SONUÇ (25 Ağustos) — İKİ SÜRPRİZ, YORUM AÇIK:
  S1  η AĞIR-KUYRUKLU: skew(η)=+0.24, kurt(η)=+0.86 (iki pencerede
      tutarlı). 111'in "alt-Gauss yapı" adayı η DÜZEYİNDE RET (73'ün
      −0.75 kurt'u u içindi; η'ya taşınmıyor). Yeni defter sayısı.
  S2  ÇARPIKLIK DALGASI VAR ve anti-fazlı: R₃ = −1.4..−2.0 (σ³-normlu
      −1.0..−1.45), plasebo 0.1-0.67 (3.-moment kanalı gürültülü;
      sinyal tabanın ~3-6 katı).
  S3  SÜRPRİZ — BOYALI REFERANS DA NEGATİF ve DAHA GÜÇLÜ: R₃ = −4.5/
      −4.1. Naif adyabatik referans (+1) bu kanal için YANLIŞMIŞ
      (ön-mühürün çerçevesi kısmen geçersiz — dürüst kayıt). Kinematik
      3.-moment yanıtının kendi teorisi gerekiyor (dilatasyonun kübik
      muhasebesi; sort/asimetri terimleri aday).
  OKUMA (geçici): gerçek gazın çarpıklık dalgası kinematik referansın
      ~0.35-0.45'i — tanıma bu kanalda işaret-çevirme değil GÜÇLÜ
      BASTIRMA olarak görünüyor. Köprüye sayılması, R₃-kinematiğinin
      teorisi yazılana dek BEKLEMEDE. Kalan-%25 aday listesi güncel:
      (1) alt-Gauss RET; (2) 3.-moment yapısı ölçüldü, yorum açık;
      (3) kaynak-verteksi hâlâ masada.
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

d41 = np.load(HERE / "41_bigT_windows.npz")
PR = [2, 3, 5, 7, 11, 13, 17]
SONUC = {}
for k in ["120k", "200k"]:
    gz, tm = d41[f"gaps_{k}"], d41[f"tmid_{k}"]
    zz = np.empty(len(gz) + 1)
    zz[0] = tm[0] - gz[0] / 2
    zz[1:] = zz[0] + np.cumsum(gz)
    g = np.diff(zz)
    m = 0.5 * (zz[:-1] + zz[1:])
    Lw = np.log(m / TWO_PI)
    L = float(Lw.mean())
    ds = g * Lw / TWO_PI - 1
    qs = pk_list(min(np.exp(0.52 * L), 720))
    freqs = [np.log(q) for q in qs]
    plas = [np.log(p) + 0.037 for p in PR]
    allf = freqs + plas
    b1, fit1 = chunked_fit(ds, m, allf)
    eta = ds - fit1
    s2 = float((eta**2).mean())
    m3 = float((eta**3).mean())
    skew = m3 / s2**1.5
    kurt = float((eta**4).mean()) / s2**2 - 3
    print(f"[{k}] L={L:.2f}  σ_η={np.sqrt(s2):.4f}  skew(η)={skew:+.3f}  "
          f"kurt(η)={kurt:+.3f}  ⟨η³⟩={m3:+.3e}", flush=True)
    y3 = eta**3 - m3
    b3, _ = chunked_fit(y3, m, allf)
    for p in PR:
        f = np.log(p)
        cA, sA = coef(b1, allf, f)
        A1 = np.hypot(cA, sA); ph = np.arctan2(sA, cA)
        c3, s3 = coef(b3, allf, f)
        P3 = c3 * np.cos(ph) + s3 * np.sin(ph)
        R3 = P3 / (3 * m3 * A1)
        R3s = P3 / (s2**1.5 * A1)          # σ³-normlu (m3-bağımsız)
        fp = np.log(p) + 0.037
        c3p, s3p = coef(b3, allf, fp)
        plv = np.hypot(c3p, s3p) / (abs(3 * m3) * A1)
        SONUC.setdefault(p, []).append((A1, R3, R3s, plv))

print(f"\n{'p':>3} {'⟨A1⟩':>7} {'R₃':>7} {'R₃(σ³)':>8} {'plasebo':>8}"
      f"   (adyabatik referans R₃=+1)")
for p in PR:
    arr = np.array(SONUC[p])
    a = arr.mean(axis=0)
    print(f"{p:>3} {a[0]:>7.4f} {a[1]:>+7.3f} {a[2]:>+8.4f} {a[3]:>8.3f}")

# ---- S3: boyalı kontrol (ζ-120k, çizgi-dışı ω*)
print("\nS3 — boyalı kontrol (ζ-120k):", flush=True)
gz, tm = d41["gaps_120k"], d41["tmid_120k"]
zz = np.empty(len(gz) + 1)
zz[0] = tm[0] - gz[0] / 2
zz[1:] = zz[0] + np.cumsum(gz)
L = float(np.log(tm / TWO_PI).mean())
gbar = TWO_PI / L
qs = pk_list(min(np.exp(0.52 * L), 720))
LINES = [np.log(q) for q in qs]
for taus in [0.15, 0.30]:
    om = taus * L
    while min(abs(om - l) for l in LINES) < 0.012:
        om += 0.013
    U = 0.1 / (2 * np.sin(np.pi * om / L))
    zp = zz + U * gbar * np.cos(om * zz)
    gp = np.diff(zp)
    mp_ = 0.5 * (zp[:-1] + zp[1:])
    ds = gp * np.log(mp_ / TWO_PI) / TWO_PI - 1
    allf = [np.log(q) for q in qs] + [om]
    b1, fit1 = chunked_fit(ds, mp_, allf)
    eta = ds - fit1
    s2 = float((eta**2).mean()); m3 = float((eta**3).mean())
    y3 = eta**3 - m3
    b3, _ = chunked_fit(y3, mp_, allf)
    cA, sA = coef(b1, allf, om)
    A1 = np.hypot(cA, sA); ph = np.arctan2(sA, cA)
    c3, s3 = coef(b3, allf, om)
    P3 = c3 * np.cos(ph) + s3 * np.sin(ph)
    print(f"  τ*={taus:.2f}: A1={A1:.4f}  R₃={P3/(3*m3*A1):+.3f}", flush=True)
