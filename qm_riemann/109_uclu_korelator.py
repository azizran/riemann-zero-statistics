"""
109 — ÜÇLÜ-KORELATÖR: GAZIN KENDİ DALGASINA FAZ-ÇÖZÜMLÜ YANITI (24 Ağu)
==========================================================================
108'in bilançosu: perde evrensel, sıcaklık-bağımsız, alet-temiz ve
SEÇİCİ (yalnız gazın kendi dalgaları). Kalan mekanizma sınıfı: gap-
dalgalanmalarının aritmetik dalganın FAZIYLA korelasyonu. Burada bu
korelasyon DOĞRUDAN ölçülür: dalga-ortalaması söküldükten sonra kalan
gürültü η'nın İKİNCİ momentleri çizgi fazına kilitli mi?

Kanallar (hepsi tam-taban ortak regresyonla — tek-çizgi tuzağı yok):
  M1  ⟨ds⟩(faz): birinci harmonik (bilinen v; çapa)
  M3  ⟨η²⟩(faz): VARYANS DALGASI — gazın "nefesi"
  M4  ⟨η_n η_{n+1}⟩(faz): komşu-bağı dalgası
Adyabatik yerel-ölçekleme öngörüsü (parametresiz): yerel yoğunluk
dalgayla ölçeklenirse ds|faz ~ ε(faz) + (1+ε)η₀ → varyans dalgası
genliği = 2σ_η²·ε_amp, yani R_p = b_var/(2σ_η²·b_ort) = 1; komşu-bağı
için R_nn = b_nn/(2c₁·b_ort) = 1 (c₁ = ⟨η_nη_{n+1}⟩ < 0).

ÖN-MÜHÜR:
  Q1  Varyans dalgası VAR (plasebo tabanının belirgin üstünde).
  Q2  R_p ~ 1 ise saf adyabatik; R_p < 1 sistematikse az-nefes
      (katılık) — perdenin faz-uzayı izi; işaret/faz raporlanır.
  Q3  R_nn, R_p ile aynı ölçekte.
Plasebo: aynı okuma çizgi+0.037 kaydırılmış sahte frekanslarda.
Veri: ζ 120k+200k (ayrı ayrı, ortalama). Asallar 2..17.

SONUÇ (24 Ağustos) — İKİ MANŞET + BİR TUZAK-KONTROLÜ:
  M0 (yan manşet): var(ds)=0.164 (GUE-yakın ✓ veri sağlam; 109b-teşhis)
    ama σ_η²=0.022 → GAP VARYANSININ ~%87'Sİ DETERMİNİSTİK açık-formül
    dalgaları (Σ A²/2 ≈ 0.14, ortogonallik muhasebesi tutuyor; GUE-biçim
    dağılım ~60 deterministik dalganın süperpozisyonundan doğuyor —
    Bogomolny-Keating kokusu). Artık c₁/σ² = −0.52 (GUE tam-gap −0.27'den
    güçlü anti-korelasyon).
  Q1 ✓ NEFES DALGASI VAR: plasebo 0.03-0.28, sinyal 0.69-0.92.
  Q2 — TERS NEFES (ANTİ-ADYABATİK): R_p ≈ −0.85 (−0.69..−0.92),
    7 asalda DÜZ. Adyabatik yerel-ölçekleme +1 derdi: gerilen bölge
    SESSİZLEŞİYOR, sıkışan bölge GÜRÜLTÜLENİYOR — işaret ters, %85 güç.
  Q3 — komşu-bağı dalgası R_nn = −1.9→−1.2 (p ile iniyor), plasebo temiz.
  TUZAK-KONTROLÜ (dalga-çifti özdeşliği): η²'nin ω=log p'deki sinyaline
    kaldırılan dalgaların çapraz çarpımları (q1/q2=p çiftleri) karışabilir;
    böyle çift p=2,3 için güçlü, 5,7 için zayıf, 11,13,17 İÇİN YOK —
    ama R_p yedi asalda DÜZ → beat-özdeşliği kaynak OLAMAZ (açık hesapla
    doğrulama yine de yapılacak). SIRADAKİ: CUE-referans R (dış dalga
    enjekte edilmiş CUE'da aynı ölçüm); adalar evrenselliği; perde
    D≈1−0.36τ ile nicel köprü.
"""

import numpy as np
from pathlib import Path
from sympy import primerange

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
Q_CAP = 720

def pk_list(lim):
    out = []
    for p in primerange(2, int(lim) + 1):
        q = p
        while q <= lim:
            out.append(q); q *= p
    return sorted(set(out))

def chunked_fit(y, tmid, freqs, chunk=40000):
    """y'yi sabit+2 sürüklenme+cos/sin(freqs) üstüne düş; (b, fit) döner."""
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
    qs = pk_list(min(np.exp(0.52 * L), Q_CAP))
    freqs = [np.log(q) for q in qs]
    plas = [np.log(p) + 0.037 for p in PR]          # sahte çizgiler
    allf = freqs + plas
    b1, fit1 = chunked_fit(ds, m, allf)
    eta = ds - fit1
    s2 = float((eta**2).mean())
    y2 = eta**2 - s2
    b2, _ = chunked_fit(y2, m, allf)
    ee = eta[:-1] * eta[1:]
    c1 = float(ee.mean())
    mm = 0.5 * (m[:-1] + m[1:])
    y3 = ee - c1
    b3, _ = chunked_fit(y3, mm, allf)
    print(f"[{k}] L={L:.2f}  σ_η²={s2:.4f}  c₁={c1:.4f}  "
          f"(c₁/σ² = {c1/s2:.3f})", flush=True)
    for p in PR:
        f = np.log(p)
        c1v, s1v = coef(b1, allf, f)
        A1 = np.hypot(c1v, s1v); ph1 = np.arctan2(s1v, c1v)
        c2v, s2v = coef(b2, allf, f)
        # ort-dalga yönüne işaretli izdüşüm
        P2 = (c2v * np.cos(ph1) + s2v * np.sin(ph1))
        c3v, s3v = coef(b3, allf, f)
        P3 = (c3v * np.cos(ph1) + s3v * np.sin(ph1))
        fp = np.log(p) + 0.037
        c2p, s2p = coef(b2, allf, fp)
        c3p, s3p = coef(b3, allf, fp)
        SONUC.setdefault(p, []).append(
            (A1, P2 / (2 * s2 * A1), np.hypot(c2p, s2p) / (2 * s2 * A1),
             P3 / (2 * c1 * A1), np.hypot(c3p, s3p) / (2 * abs(c1) * A1)))

print(f"\n{'p':>3} {'⟨A1⟩':>7} {'R_p':>7} {'plas_R':>7} {'R_nn':>7} "
      f"{'plas_nn':>8}   (iki pencere ort.)")
for p in PR:
    arr = np.array(SONUC[p])
    a = arr.mean(axis=0)
    print(f"{p:>3} {a[0]:>7.4f} {a[1]:>7.3f} {a[2]:>7.3f} {a[3]:>7.3f} "
          f"{a[4]:>8.3f}")
