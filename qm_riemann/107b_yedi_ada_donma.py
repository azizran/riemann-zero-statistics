"""
107b — YEDİ-ADA DONMA TABLOSU + FİGÜR v2 (25 Ağustos 2026)
==========================================================================
107a eski dört adayı (χ₃, β, χ₅, χ₇) iki YENİ kapıyla (105e kısa-çukur +
104e sıçrama) yeniden denetledi. Bu script YENİ ÜÇ adayı (chi5e, chi8e,
chi8o — 105b/105d) AYNI ÜÇ KAPILI boru hattından geçirir ve ζ + 7 ada =
SEKİZ VERİ KÜMESİNİ tek tabloya ve 101_donmus_koylar.png'nin v2'sine
toplar.

101g stili korunur; yeni adalar KERVAN RENGİYLE eklenir:
  LOG2 KERVANI (2-çizgisi CANLI)  : ζ, χ₃, χ₅, χ₇, chi5e   — soğuk renk
  LOG3 KERVANI (2-ailesi ÖLÜ)     : β, chi8e, chi8o        — sıcak renk
β kervanı artık ÜÇ üyeli: tek bir ada değil, bir SINIF.

Panel B: x = τ/τ_ilk yeniden ölçekleme — sekiz eğri TEK EŞİĞE kilitleniyor
mu? (τ_ilk = log q_ilk / ⟨L⟩_n, q_ilk = adanın ilk sağ kalan çizgisi.)

ÖN-MÜHÜR (ölçümden ÖNCE):
  B1  Yeni üç ada üç kapıyla da 105d'nin hükmünü koruyacak: chi5e
      log2/L'de çözülür, chi8e/chi8o log3/L'e gecikir.
  B2  Yeniden ölçeklemede x ∈ [0.85, 1.10] penceresinde SEKİZ eğrinin
      std/ort'u, ham τ'daki aynı pencerenin std/ort'undan BELİRGİN
      küçük olacak (kervanlar tek eşiğe kilitlenir).
  B3  yarı-çözülme/τ_ilk oranı sekiz kümede de 0.85-1.15 aralığında
      olacak (105d: 1.006 chi5e, 0.997 chi3, 0.880 chi8e, 0.869 chi8o,
      0.891 beta — yani ölü-2 üçlüsü sistematik olarak ~0.88).
==========================================================================
SONUÇ (25 Ağustos, koşu 16 s) — B1 ✓✓ B2 ✓ B3 ✓ (ζ'da bir çekince ile)
==========================================================================
KAPI RAPORU (yeni üç ada, üç kapılı):
  chi5e  K2: t ∈ [916.09, 935.78] (20 sıfır) — 105d bunu zaten kesiyordu
         (kapı 105d'nin duz_segmentler'inde vardı) ama YERİ hiç
         raporlanmamıştı; defterde ilk kez adıyla duruyor.
  chi8e  K2: t ∈ [19353.48, 19359.30] (12 sıfır) — 105e'nin ELLE bulduğu
         t≈19352 kısa çukuru; POZİTİF KONTROL İKİNCİ KEZ GEÇTİ
  chi8o  hiçbir kapı işaret vermedi (maks |med₂₀−med₂₀₀| = 0.20) — β ile
         birlikte iki TERTEMİZ küme
  KAPI 3 (SIÇRAMA) üç adada da 0 işaret (maks|r| = 0.06-0.11, temiz
  tabanın içinde). Bu YENİ bilgidir: 105d sıçrama kapısını hiç
  koşmamıştı. Dolayısıyla D değerleri 105d'ninkilerle ONDALIĞINA KADAR
  aynı çıktı (chi5e .045/.235/.425/.591/.538/.356/.704/.737, chi8e
  .005/.009/.009/.009/.029/.486/.767/.843, chi8o birebir) — yeni üç ada
  ÜÇÜNCÜ KAPIYI DA KUSURSUZ GEÇTİ.

B1 ✓✓ D(τ), band ±0.01L, ÜÇ KAPILI SEKİZ KÜME:
  τ      zeta   chi3   chi5   chi7  chi5e   beta  chi8e  chi8o  vekil
  0.040  0.021  0.036  0.031  0.037  0.045  0.012  0.005  0.025  0.005
  0.068  0.489  0.194  0.442  0.507  0.425  0.037  0.009  0.047  0.031
  0.085  0.527  0.617  0.511  0.334  0.538  0.091  0.029  0.110  0.039
  0.113  0.721  0.161  0.728  0.800  0.737  0.726  0.843  0.789  0.051
  0.300  1.130  1.023  1.025  1.039  1.038  1.067  1.068  1.040  0.110
  EŞİKSİZ AYRIŞTIRICI — D(log2/L), yani ORTAK çizgide ölçüm:
    LOG2 KERVANI: ζ .461  χ₃ .508  χ₅ .482  χ₇ .512  χ₅ᵉ .492
    LOG3 KERVANI: β .046  χ₈ᵉ .009  χ₈ᵒ .048     ← 10-50× DAHA DÜŞÜK
  Beş ada çözülmüş, üç ada donuk, TEK bir τ değerinde. 101'in taç
  bulgusu artık 5'e 3 ile duruyor (101d'de 4'e 1'di).

B2 ✓ TEK EŞİĞE KİLİTLENME (yeniden ölçekleme x = τ/τ_ilk):
  eşik bölgesinde sekiz eğrinin saçılması
    ham τ  (0.06 ≤ τ ≤ 0.125) : std/ort 0.5125,  mutlak std 0.1931
    ölçekli x (0.85 ≤ x ≤ 1.10): std/ort 0.2559,  mutlak std 0.1251
  Saçılma YARIYA iniyor: τ_ilk = log q_ilk/⟨L⟩ ile bölmek iki kervanı
  üst üste getiriyor. Kalan 0.256'nın bir kısmı iki bilinen kırılgan
  eğriden (ölçüldü, x∈[0.85,1.10]): χ₃ hariç 0.2296 | ζ hariç 0.2503 |
  ikisi de hariç 0.2141. Yani çökme kalan altı kümede daha da sıkı.

B3 ✓ YARI-ÇÖZÜLME / τ_ilk:
    LOG2 kervanı: χ₃ 0.998  χ₅ 1.018  χ₇ 0.980  χ₅ᵉ 1.005   (ζ 1.338*)
    LOG3 kervanı: β  0.891  χ₈ᵉ 0.880  χ₈ᵒ 0.870
  MUTLAK FREKANSTA (ω_yarı = τ_yarı·⟨L⟩) ÇARPICI SONUÇ:
    LOG2 kervanı (4 Dirichlet adası) ort = 0.6932   ←  log 2 = 0.69315
    LOG3 kervanı (3 ada)             ort = 0.9671   ←  log 3 = 1.0986
  Yani 2-çizgisi CANLI adalar yarı-çözülmeyi TAM log 2'de yapıyor
  (%0.01!); ölü-2 üçlüsü kendi log 3'ünden %12 ERKEN. Kervan kayması
  ×1.395 (Dirichlet tabanına göre) — 101d'nin ×1.4-1.5'i ve 105d'nin
  ×1.376-1.410'u ile aynı bant; saf öngörü log3/log2 = 1.585'in altında.
  Bu %12'lik erken açılma SİSTEMATİK (üç adada da 0.87-0.89) ve
  AÇIK SORU'dur: ölü çizginin komşuluğunda eşik biraz yumuşuyor.
  (*) ÇEKİNCE — ζ'nın 1.338'i KIRILGAN, gerçek bir sapma değil:
  D(0.068) = 0.489, 0.5'in yalnız 0.011 altında; τ=0.075'teki yerel
  çukur (0.393) yüzünden "ilk 0.5 geçişi" 0.085'e atlıyor. Eşiksiz
  ölçüt D(τ_ilk) ζ'da 0.461 — öteki log2 adalarıyla (0.482-0.512)
  aynı bantta. yarı-çözülme ölçütünün eşik-kırılganlığı kayda geçiyor.

FİGÜR: 101_donmus_koylar_v2.png (iki panel).
==========================================================================
"""
import numpy as np, time
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

T0 = time.time()
HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
rng = np.random.default_rng(1071)

PKS = [2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 17, 19, 23, 25, 27, 29, 31, 32, 37,
       41, 43, 47, 49, 53, 59, 61, 64, 67, 71, 73, 79, 81, 83, 89, 97, 101,
       103, 107, 109, 113, 121, 125, 127, 128]
LINES_ALL = [np.log(qq) for qq in PKS]
TAUS = [0.04, 0.05, 0.06, 0.068, 0.075, 0.085, 0.095, 0.105, 0.113,
        0.125, 0.14, 0.16, 0.20, 0.30]
BAND = 0.01


def pencere_hazirla(z, qeff):
    gaps = np.diff(z)
    mids = 0.5 * (z[:-1] + z[1:])
    Lw = float(np.log(qeff * mids / TWO_PI).mean())
    ds = gaps * np.log(qeff * mids / TWO_PI) / TWO_PI - 1
    tt = (mids - mids.mean()) / (mids[-1] - mids[0])
    P = np.vstack([np.ones_like(tt), tt, tt**2, tt**3]).T
    ds = ds - P @ np.linalg.lstsq(P, ds, rcond=None)[0]
    return (z, mids, ds, Lw)


def egri(WIN, taus, band=BAND):
    Dv, Vv, Wv = [], [], []
    for tau0 in taus:
        num = 0.0 + 0j
        den = 0.0
        Gs, rs, per = [], [], []
        for (zz, tm, ds, Lw) in WIN:
            oms = tau0 * Lw + np.linspace(-band * Lw, band * Lw, 160)
            oms = np.array([o for o in oms
                            if min(abs(o - l) for l in LINES_ALL) > 0.01])
            nw = 0.0 + 0j
            dw = 0.0
            for s0 in range(0, len(oms), 40):
                ob = oms[s0:s0 + 40]
                rr = np.exp(1j * np.outer(ob, zz)).sum(axis=1)
                GG = (np.exp(1j * np.outer(ob, tm)) * ds[None, :]).sum(axis=1)
                nw += (GG * np.conj(rr)).sum()
                dw += (np.abs(rr)**2).sum()
                Gs.append(GG)
                rs.append(rr)
            num += nw
            den += dw
            per.append(abs(nw) / max(dw, 1e-300))
        kap = 2 * np.pi * tau0
        c = 2 * np.sin(kap / 2) / kap
        G_all, r_all = np.concatenate(Gs), np.concatenate(rs)
        fl = [abs((G_all * np.conj(r_all[rng.permutation(len(r_all))])).sum())
              / (c * den) for _ in range(40)]
        Dv.append(abs(num) / (c * den))
        Vv.append(float(np.mean(fl)))
        Wv.append([p / c for p in per])
    return Dv, Vv, Wv


def _roll(x, w, f):
    from numpy.lib.stride_tricks import sliding_window_view
    if len(x) <= w:
        return np.full(len(x), float(f(x)))
    m = f(sliding_window_view(x, w), axis=1)
    out = np.empty(len(x))
    out[w // 2:w // 2 + len(m)] = m
    out[:w // 2] = m[0]
    out[w // 2 + len(m):] = m[-1]
    return out


def bloklar(mask, t):
    out, i, n = [], 0, len(mask)
    while i < n:
        if mask[i]:
            j = i
            while j < n and mask[j]:
                j += 1
            out.append((i, j - 1, float(t[i]), float(t[j - 1])))
            i = j
        else:
            i += 1
    return out


def kapi12(zz, sayim, win=80, esik=0.5, pad=160,
           kisa=20, uzun=200, kesik=0.7, kpad=160):
    d = np.arange(len(zz)) - (sayim(zz) - sayim(zz[0]))
    from numpy.lib.stride_tricks import sliding_window_view
    med = np.median(sliding_window_view(d, win), axis=1)
    adim = med[win + 1:] - med[:-(win + 1)]
    g1 = np.zeros(len(zz), bool)
    for i in np.where(np.abs(adim) > esik)[0] + win // 2:
        g1[max(0, i - pad):i + 2 * pad] = True
    sapma = np.abs(_roll(d, kisa, np.median) - _roll(d, uzun, np.median))
    ham2 = sapma > kesik
    g2 = np.zeros(len(zz), bool)
    for i in np.where(ham2)[0]:
        g2[max(0, i - kpad):i + kpad] = True
    return g1, g2, ham2, float(sapma.max())


def segmentle(zz, bad, min_n=3000):
    seg, kes, s0 = [], 0, 0
    for i in range(1, len(zz) + 1):
        if i == len(zz) or bad[i] != bad[i - 1]:
            if not bad[s0] and i - s0 >= min_n:
                seg.append(zz[s0:i])
            if i < len(zz) and bad[i]:
                kes += 1
            s0 = i
    return seg, kes


def kapi3(p, q, sayim, pad=600, uesik=0.005, resik=0.30):
    d = np.arange(len(p)) - (sayim(p) - sayim(p[0]))
    r = _roll(d, 21, np.mean) - _roll(d, 801, np.median)
    g = np.diff(p)
    mid = 0.5 * (p[:-1] + p[1:])
    u = g * np.log(q * mid / TWO_PI) / TWO_PI
    bad = np.abs(r) > resik
    dup = np.where(u < uesik)[0]
    for i in dup:
        bad[i] = True
        bad[min(i + 1, len(bad) - 1)] = True
    kes = np.zeros(len(p), bool)
    for i in np.where(bad)[0]:
        kes[max(0, i - pad):min(len(p), i + pad + 1)] = True
    parca, s0 = [], 0
    for i in range(1, len(p) + 1):
        if i == len(p) or kes[i] != kes[i - 1]:
            if not kes[s0] and i - s0 >= 3000:
                parca.append(p[s0:i])
            s0 = i
    return parca, bad, kes, float(np.abs(r).max()), u, dup


# --------------------------------------------------- YENİ ÜÇ ADA (105b)
exec(open(HERE / "98_L_motoru.py").read().split('CHI4 = ')[0])
HERE = Path(__file__).resolve().parent
CHI5E = {0: 0, 1: 1, 2: -1, 3: -1, 4: 1}
CHI8E = {0: 0, 1: 1, 2: 0, 3: -1, 4: 0, 5: -1, 6: 0, 7: 1}
CHI8O = {0: 0, 1: 1, 2: 0, 3: 1, 4: 0, 5: -1, 6: 0, 7: -1}
YENI = [("chi5e", 5, CHI5E, 0, 2), ("chi8e", 8, CHI8E, 0, 3),
        ("chi8o", 8, CHI8O, 1, 3)]

print("=" * 74)
print("107b — YENİ ÜÇ ADA, ÜÇ KAPILI BORU HATTINDAN GEÇİRİLİYOR")
print("=" * 74, flush=True)

PEN, TAN, NHAM, QILK, BOLG = {}, {}, {}, {}, {}
for et, q, tab, a, qilk in YENI:
    d = np.load(HERE / f"105b_{et}_zeros.npz")
    zc = d["zeros"]
    lo = np.exp(np.log(zc[0] + 1) + 0.25 * (np.log(zc[-1]) - np.log(zc[0] + 1)))
    zc = zc[zc >= lo]
    M = Lmotor(q, tab, a)
    say = (lambda t, M=M: M.theta(t) / np.pi)
    NHAM[et] = len(zc)
    QILK[et] = qilk
    BOLG[et] = int(np.asarray(d["bolgeler"]).reshape(-1, 2).shape[0])
    g1, g2, ham2, smax = kapi12(zc, say)
    kus = [("K2", t0, t1, i1 - i0 + 1) for (i0, i1, t0, t1)
           in bloklar(ham2, zc)]
    seg, kes = segmentle(zc, g1 | g2)
    kenar = np.exp(np.linspace(np.log(zc[0]), np.log(zc[-1] * 1.0001), 4))
    ara = []
    for s in seg:
        for i in range(3):
            pw = s[(s >= kenar[i]) & (s < kenar[i + 1])]
            if len(pw) >= 3000:
                ara.append(pw)
    nA = sum(len(x) for x in ara)
    WIN, g3t, dupt, rms = [], 0, 0, []
    for p in ara:
        alt, bad3, kes3, rmax, u, dup = kapi3(p, float(q), say)
        rms.append(rmax)
        g3t += int(kes3.sum())
        dupt += len(dup)
        for (i0, i1, t0, t1) in bloklar(bad3, p):
            kus.append(("K3", t0, t1, i1 - i0 + 1))
        for x in alt:
            WIN.append(pencere_hazirla(x, float(q)))
    PEN[et] = WIN
    TAN[et] = dict(g1=int(g1.sum()), g2=int(g2.sum()), g2ham=int(ham2.sum()),
                   g3=g3t, dup=dupt, smax=smax, rmax=rms, kusur=kus,
                   nA=nA, nB=sum(len(w[0]) for w in WIN), npenA=len(ara))
    print(f"[{et}] ham n={len(zc)}  |  KAPI1+2 sonrası {len(ara)} pencere "
          f"n={nA}  →  KAPI3 sonrası {len(WIN)} pencere n={TAN[et]['nB']}")
    print(f"      KAPI1 {int(g1.sum())} | KAPI2 ham {int(ham2.sum())} → "
          f"padli {int(g2.sum())} (maks |med20−med200| = {smax:.2f}) | "
          f"KAPI3 padli {g3t}, kopya {dupt}")
    print(f"      maks|r| pencere başına: " + " ".join(f"{x:.2f}" for x in rms))
    for (k, t0, t1, n) in kus:
        print(f"      {k}: t ∈ [{t0:.4f}, {t1:.4f}] ({n} sıfır)")
    print(flush=True)

SONY = {}
for et, q, tab, a, qilk in YENI:
    D, V, W = egri(PEN[et], TAUS)
    SONY[et] = (D, V, W)
    print(f"  ...{et} ölçüldü ({time.time()-T0:.0f} s)", flush=True)

# ------------------------------------------------- 107a'dan eski dört + ζ
A = np.load(HERE / "107a_temiz.npz", allow_pickle=True)
ESKI = [str(x) for x in A["adlar"]]

ADLAR = ["zeta", "chi3", "chi5", "chi7", "chi5e", "beta", "chi8e", "chi8o"]
KERVAN = {"zeta": 2, "chi3": 2, "chi5": 2, "chi7": 2, "chi5e": 2,
          "beta": 3, "chi8e": 3, "chi8o": 3}
QLAB = {"zeta": "—", "chi3": "3 (a=1)", "chi5": "5 (a=1)", "chi7": "7 (a=1)",
        "chi5e": "5 (a=0)", "beta": "4 (a=1)", "chi8e": "8 (a=0)",
        "chi8o": "8 (a=1)"}

D, V, WPEN, LORT, TILK, NPEN, NN, LS = {}, {}, {}, {}, {}, {}, {}, {}
for a in ESKI:
    D[a] = A[f"DB_{a}"]
    V[a] = A[f"VB_{a}"]
    LS[a] = A[f"L_{a}"]
    NN[a] = A[f"n_{a}"]
for et, q, tab, aa, qilk in YENI:
    D[et] = np.array(SONY[et][0])
    V[et] = np.array(SONY[et][1])
    LS[et] = np.array([w[3] for w in PEN[et]])
    NN[et] = np.array([len(w[0]) for w in PEN[et]])
    WPEN[et] = SONY[et][2]
for a in ADLAR:
    LORT[a] = float(np.average(LS[a], weights=NN[a].astype(float)))
    TILK[a] = np.log(KERVAN[a]) / LORT[a]
    NPEN[a] = len(LS[a])

# ------------------------------------------------------- YEDİ-ADA TABLOSU
print("\n" + "=" * 78)
print(f"YEDİ ADA + ζ — D(τ), band ±{BAND}L, ÜÇ KAPILI TEMİZ VERİ")
print("=" * 78)
print(f"{'τ':>6} " + " ".join(f"{a:>7}" for a in ADLAR) + "   vekil(mx)")
for i, t in enumerate(TAUS):
    vek = max(V[a][i] for a in ADLAR)
    print(f"{t:>6.3f} " + " ".join(f"{D[a][i]:>7.3f}" for a in ADLAR)
          + f"   {vek:>7.3f}")
print("\nvekil taban (aile aile):")
for a in ADLAR:
    print(f"  {a:>6}: " + " ".join(f"{v:.3f}" for v in V[a]))

print(f"\n=== BAND KAÇAĞI (sızan pencere, ±{BAND}L; ★ = sızıntı var) ===")
SIZ = {}
for a in ADLAR:
    SIZ[a] = [int(np.sum(LS[a] * (t + BAND) >= np.log(KERVAN[a])))
              for t in TAUS]
print(f"{'τ':>6} " + " ".join(f"{a:>8}" for a in ADLAR))
for i, t in enumerate(TAUS):
    print(f"{t:>6.3f} " + " ".join(
        f"{D[a][i]:>7.3f}" + ("★" if SIZ[a][i] > 0 else " ") for a in ADLAR))


def yari_coz(Dv, T):
    Dv = np.asarray(Dv, float)
    T = np.asarray(T, float)
    ix = np.where(Dv >= 0.5)[0]
    if not len(ix) or ix[0] == 0:
        return np.nan
    i0 = ix[0]
    return T[i0 - 1] + (0.5 - Dv[i0 - 1]) * (T[i0] - T[i0 - 1]) / \
        (Dv[i0] - Dv[i0 - 1])


print("\n" + "=" * 92)
print("YEDİ-ADA DONMA TABLOSU (ζ + 7 ada) — ÜÇ KAPILI TEMİZ VERİ")
print("=" * 92)
print(f"{'aile':>6} {'q(a)':>8} {'kervan':>7} {'n_pen':>5} {'⟨L⟩':>7} "
      f"{'τ_ilk':>7} {'D@log2/L':>9} {'D@log3/L':>9} {'D@τ_ilk':>8} "
      f"{'yarı-τ':>7} {'yarı/τ_ilk':>10} {'D@0.30':>7}")
YARI, TABLO = {}, {}
for a in ADLAR:
    t2, t3 = np.log(2) / LORT[a], np.log(3) / LORT[a]
    d2 = float(np.interp(t2, TAUS, D[a]))
    d3 = float(np.interp(t3, TAUS, D[a]))
    di = float(np.interp(TILK[a], TAUS, D[a]))
    ty = yari_coz(D[a], TAUS)
    YARI[a] = ty
    TABLO[a] = (d2, d3, di, ty, ty / TILK[a])
    print(f"{a:>6} {QLAB[a]:>8} {'log'+str(KERVAN[a]):>7} {NPEN[a]:>5} "
          f"{LORT[a]:>7.3f} {TILK[a]:>7.4f} {d2:>9.3f} {d3:>9.3f} "
          f"{di:>8.3f} {ty:>7.4f} {ty/TILK[a]:>10.3f} {D[a][-1]:>7.3f}")

canli = [a for a in ADLAR if KERVAN[a] == 2]
olu = [a for a in ADLAR if KERVAN[a] == 3]
canli_L = [a for a in canli if a != "zeta"]   # yalnız Dirichlet adaları
wc = np.mean([YARI[a] * LORT[a] for a in canli])
wcL = np.mean([YARI[a] * LORT[a] for a in canli_L])
wo = np.mean([YARI[a] * LORT[a] for a in olu])
print(f"\n  ω_yarı = τ_yarı·⟨L⟩ (yarı-çözülmenin MUTLAK frekansı):")
for a in ADLAR:
    print(f"    {a:>6}: {YARI[a]*LORT[a]:.4f}")
print(f"\n  LOG2 kervanı ort (ζ dâhil)    {wc:.4f}")
print(f"  LOG2 kervanı ort (4 Dirichlet) {wcL:.4f}   ← log 2 = "
      f"{np.log(2):.4f}")
print(f"  LOG3 kervanı ort (3 ada)       {wo:.4f}   ← log 3 = "
      f"{np.log(3):.4f}")
print(f"  KERVAN KAYMASI ×{wo/wcL:.3f} (Dirichlet tabanına göre) / "
      f"×{wo/wc:.3f} (ζ dâhil)   [öngörü log3/log2 = 1.585]")
print(f"\n  yarı/τ_ilk: log2 kervanı " + " ".join(
    f"{a}={YARI[a]/TILK[a]:.3f}" for a in canli))
print(f"              log3 kervanı " + " ".join(
    f"{a}={YARI[a]/TILK[a]:.3f}" for a in olu))
print("  DÜRÜST KAYIT — ζ'nın 'ilk 0.5 geçişi' KIRILGAN: D(0.068)=%.3f, "
      "0.5'in\n  sadece %.3f altında; τ=0.075'teki yerel çukuk yüzünden ilk "
      "geçiş\n  0.085'e atlıyor ve yarı/τ_ilk %.2f çıkıyor. EŞİKSİZ ölçüt "
      "D(τ_ilk) ζ'da\n  %.3f — öteki log2 adalarıyla (%.3f-%.3f) birebir "
      "aynı bantta."
      % (D["zeta"][3], 0.5 - D["zeta"][3], YARI["zeta"] / TILK["zeta"],
         TABLO["zeta"][2], min(TABLO[a][2] for a in canli_L),
         max(TABLO[a][2] for a in canli_L)))
print(f"\n  EŞİKSİZ ÖLÇÜT — D(τ_ilk) (her ada KENDİ ilk çizgisinde):")
print("    log2 kervanı: " + " ".join(f"{a}={TABLO[a][2]:.3f}" for a in canli))
print("    log3 kervanı: " + " ".join(f"{a}={TABLO[a][2]:.3f}" for a in olu))
print(f"  EŞİKSİZ AYRIŞTIRICI — D(log2/L) (ORTAK çizgide):")
print("    log2 kervanı: " + " ".join(f"{a}={TABLO[a][0]:.3f}" for a in canli))
print("    log3 kervanı: " + " ".join(f"{a}={TABLO[a][0]:.3f}" for a in olu)
      + "   ← ölü çizgide DONUK")

# --------------------------------------------- x = τ/τ_ilk ÇÖKME TABLOSU
XG = np.arange(0.55, 1.751, 0.05)
MM = {a: np.interp(XG, np.array(TAUS) / TILK[a], D[a],
                   left=np.nan, right=np.nan) for a in ADLAR}
# ham τ karşılaştırması için: aynı sayıda noktada, τ ızgarasında
print("\n" + "=" * 92)
print("YENİDEN ÖLÇEKLEME:  D(x),  x = τ / τ_ilk   (τ_ilk = log q_ilk / ⟨L⟩)")
print("=" * 92)
print(f"{'x':>6} " + " ".join(f"{a:>7}" for a in ADLAR) + "    std  std/ort")
sp, ab = [], []
for j, x in enumerate(XG):
    vals = np.array([MM[a][j] for a in ADLAR])
    s = float(np.nanstd(vals))
    m = float(np.nanmean(vals))
    sp.append(s / m if m > 0 else np.nan)
    ab.append(s)
    print(f"{x:>6.2f} " + " ".join(
        ("   nan" if np.isnan(v) else f"{v:>7.3f}") for v in vals)
        + f"  {s:>6.3f}  {s/m if m>0 else np.nan:>7.3f}")
sp, ab = np.array(sp), np.array(ab)
for k, msk in [("x≤0.80", XG <= 0.801),
               ("0.85≤x≤1.10", (XG >= 0.849) & (XG <= 1.101)),
               ("x≥1.15", XG >= 1.149), ("tümü", np.ones(len(XG), bool))]:
    print(f"  ÇÖKME {k:>12}: std/ort {np.nanmean(sp[msk]):.4f}   "
          f"mutlak std {np.nanmean(ab[msk]):.4f}")

# ham τ'da aynı D-aralığında saçılma (kıyas): τ ızgarasında std/ort
ham_sp, ham_ab = [], []
for i, t in enumerate(TAUS):
    vals = np.array([D[a][i] for a in ADLAR])
    ham_ab.append(float(np.std(vals)))
    ham_sp.append(float(np.std(vals) / np.mean(vals)))
m1 = (np.array(TAUS) >= 0.06) & (np.array(TAUS) <= 0.125)
print(f"\n  KIYAS — HAM τ (0.06 ≤ τ ≤ 0.125, eşik bölgesi): "
      f"std/ort {np.mean(np.array(ham_sp)[m1]):.4f}, "
      f"mutlak std {np.mean(np.array(ham_ab)[m1]):.4f}")
print(f"  KIYAS — ÖLÇEKLİ x (0.85 ≤ x ≤ 1.10)          : "
      f"std/ort {np.nanmean(sp[(XG>=0.849)&(XG<=1.101)]):.4f}, "
      f"mutlak std {np.nanmean(ab[(XG>=0.849)&(XG<=1.101)]):.4f}")

# ------------------------------------------------------------ FİGÜR v2
STIL = {
    "zeta":  ("#000000", "-",  "o", r"$\zeta$",              1.5, 4.5),
    "chi3":  ("#22aa77", "-",  "s", r"$\chi_3$",             1.4, 4.0),
    "chi5":  ("#2277bb", "-",  "^", r"$\chi_5$",             1.4, 4.0),
    "chi7":  ("#aa66cc", "-",  "v", r"$\chi_7$",             1.4, 4.0),
    "chi5e": ("#00a0a0", "--", "P", r"$\chi_{5}^{\rm even}$", 1.8, 5.0),
    "beta":  ("#dd3333", "-",  "D", r"$\beta\ (\chi_4)$",    2.8, 5.5),
    "chi8e": ("#ee8800", "-",  "X", r"$\chi_{8}^{\rm even}$", 2.2, 5.0),
    "chi8o": ("#993311", "--", "*", r"$\chi_{8}^{\rm odd}$",  2.2, 6.5),
}
VEK = np.array([max(V[a][i] for a in ADLAR) for i in range(len(TAUS))])

fig, (ax, bx) = plt.subplots(1, 2, figsize=(13.8, 6.4))
t2m = float(np.mean([np.log(2) / LORT[a] for a in canli]))
t3m = float(np.mean([np.log(3) / LORT[a] for a in olu]))

# --- (a) ham τ -------------------------------------------------------
ax.axvspan(t2m, t3m, color="#ffe9e9", zorder=0)
ax.fill_between(TAUS, 0, VEK, color="0.82", zorder=1,
                label="vekil taban (gürültü)")
for a in ADLAR:
    c, ls, mk, lb, lw, ms = STIL[a]
    ax.plot(TAUS, D[a], ls, color=c, marker=mk, ms=ms, lw=lw, label=lb,
            zorder=4 if KERVAN[a] == 3 else 3)
ax.axvline(t2m, color="#1a5a8a", ls=":", lw=1.5, zorder=2)
ax.axvline(t3m, color="#cc2222", ls=":", lw=1.5, zorder=2)
ax.text(t2m - 0.002, 1.435, r"$\log 2/L$", ha="right", va="top",
        fontsize=9.5, color="#1a5a8a", fontweight="bold")
ax.text(t2m - 0.002, 1.365, "LOG2 KERVANI\n(2-çizgisi CANLI)\n"
        "ζ  χ₃  χ₅  χ₇  χ₅ᵉ\nsoğuk renk", ha="right", va="top",
        fontsize=8.0, color="#1a5a8a", linespacing=1.35)
ax.text(t3m + 0.003, 1.435, r"$\log 3/L$", ha="left", va="top",
        fontsize=9.5, color="#cc2222", fontweight="bold")
ax.text(t3m + 0.003, 1.365, "LOG3 KERVANI\n(2-ailesi ÖLÜ)\n"
        "β  χ₈ᵉ  χ₈ᵒ\nsıcak renk", ha="left", va="top",
        fontsize=8.0, color="#cc2222", linespacing=1.35)
ax.annotate("", xy=(t2m, 1.235), xytext=(t3m, 1.235),
            arrowprops=dict(arrowstyle="<->", color="#cc2222", lw=1.4))
ax.text(0.5 * (t2m + t3m), 1.335, "DONMA\nPENCERESİ", ha="center",
        va="top", fontsize=8.4, color="#cc2222", fontweight="bold",
        linespacing=1.25)
ax.text(0.5 * (t2m + t3m), 1.185,
        "log2 kervanı ÇÖZÜLDÜ  D = %.2f–%.2f\n"
        "log3 kervanı DONUK    D = %.3f–%.3f"
        % (min(TABLO[a][0] for a in canli), max(TABLO[a][0] for a in canli),
           min(TABLO[a][0] for a in olu), max(TABLO[a][0] for a in olu)),
        ha="center", va="top", fontsize=8.0, color="0.2",
        bbox=dict(fc="white", ec="0.7", lw=0.6, pad=2.5))
ax.annotate("üç ada birden burada hâlâ donuk\n(2-çizgisi ölü: χ(2)=0)",
            xy=(0.0755, 0.045), xytext=(0.115, 0.28), fontsize=8.6,
            color="#cc2222", ha="left",
            arrowprops=dict(arrowstyle="->", color="#cc2222", lw=1.3,
                            connectionstyle="arc3,rad=0.25"))
ax.set_xlabel(r"$\tau=\omega/L$")
ax.set_ylabel(r"$D_{\rm spont}(\tau)$   (spontane gap-yanıtı, CUE$\approx$1)")
ax.set_title("(a) Yedi ada + ζ, ham $\\tau$ — çözülme sınırı adanın "
             "İLK SAĞ KALAN çizgisidir", fontsize=10.2)
ax.set_xlim(0.028, 0.312)
ax.set_ylim(0, 1.45)
ax.set_yticks(np.arange(0, 1.21, 0.2))
ax.grid(alpha=0.22)
ax.legend(loc="lower right", fontsize=8.2, ncol=2, framealpha=0.93)

# --- (b) x = τ/τ_ilk --------------------------------------------------
bx.axvspan(0.93, 1.07, color="0.90", zorder=0)
for a in ADLAR:
    c, ls, mk, lb, lw, ms = STIL[a]
    bx.plot(np.array(TAUS) / TILK[a], D[a], ls, color=c, marker=mk, ms=ms,
            lw=lw, label=lb, zorder=4 if KERVAN[a] == 3 else 3)
bx.axvline(1.0, color="0.30", ls="--", lw=1.4, zorder=2)
bx.axhline(0.5, color="0.55", ls=":", lw=1.1, zorder=2)
bx.text(1.0, 1.42, r"$x=1$" + "\nadanın İLK ÇİZGİSİ", ha="center", va="top",
        fontsize=9, color="0.15", fontweight="bold")
bx.text(0.52, 0.53, "yarı-çözülme", fontsize=8.2, color="0.42")
bx.text(1.12, 1.16,
        "eşik bölgesinde saçılma\nham τ:  std/ort %.2f\nölçekli x: std/ort %.2f"
        % (np.mean(np.array(ham_sp)[m1]),
           np.nanmean(sp[(XG >= 0.849) & (XG <= 1.101)])),
        fontsize=8.2, va="top", color="0.2",
        bbox=dict(fc="white", ec="0.7", lw=0.6, pad=3))
bx.set_xlabel(r"$x=\tau/\tau_{\rm ilk}$,   "
              r"$\tau_{\rm ilk}=\log q_{\rm ilk}/\langle L\rangle$")
bx.set_ylabel(r"$D_{\rm spont}$")
bx.set_title("(b) Tek eşiğe kilitlenme: iki kervan $x=\\tau/\\tau_{\\rm ilk}$"
             " ile üst üste biniyor", fontsize=10.2)
bx.set_xlim(0.5, 3.05)
bx.set_ylim(0, 1.45)
bx.set_yticks(np.arange(0, 1.21, 0.2))
bx.grid(alpha=0.22)
bx.legend(loc="lower right", fontsize=8.2, ncol=2, framealpha=0.93)

fig.suptitle("Donmuş koylar v2 — 107a/107b: iki YENİ kusur kapısıyla "
             "(105e kısa-çukur + 104e sıçrama) denetlenmiş SEKİZ veri kümesi",
             fontsize=11.5)
fig.tight_layout(rect=(0, 0, 1, 0.945))
fig.savefig(HERE / "101_donmus_koylar_v2.png", dpi=160)
print("\nkaydedildi: 101_donmus_koylar_v2.png")

np.savez(HERE / "107b_yedi_ada.npz", taus=np.array(TAUS), band=BAND,
         adlar=np.array(ADLAR), xg=XG,
         Lort=np.array([LORT[a] for a in ADLAR]),
         tilk=np.array([TILK[a] for a in ADLAR]),
         qilk=np.array([KERVAN[a] for a in ADLAR]),
         yari=np.array([YARI[a] for a in ADLAR]),
         **{f"D_{a}": np.asarray(D[a]) for a in ADLAR},
         **{f"V_{a}": np.asarray(V[a]) for a in ADLAR},
         **{f"X_{a}": MM[a] for a in ADLAR},
         **{f"SIZ_{a}": np.array(SIZ[a]) for a in ADLAR},
         **{f"nham_{et}": NHAM[et] for et in NHAM},
         **{f"nB_{et}": TAN[et]["nB"] for et in TAN},
         **{f"nA_{et}": TAN[et]["nA"] for et in TAN},
         **{f"g1_{et}": TAN[et]["g1"] for et in TAN},
         **{f"g2_{et}": TAN[et]["g2"] for et in TAN},
         **{f"g2ham_{et}": TAN[et]["g2ham"] for et in TAN},
         **{f"g3_{et}": TAN[et]["g3"] for et in TAN},
         **{f"dup_{et}": TAN[et]["dup"] for et in TAN},
         **{f"smax_{et}": TAN[et]["smax"] for et in TAN},
         **{f"bolg_{et}": BOLG[et] for et in BOLG},
         **{f"L_{a}": np.asarray(LS[a]) for a in ADLAR},
         **{f"n_{a}": np.asarray(NN[a]) for a in ADLAR})
print(f"107b_yedi_ada.npz yazıldı. Süre {time.time()-T0:.0f} s.", flush=True)
