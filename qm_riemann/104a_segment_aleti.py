"""
104a — H-seg: χ₇ τ=0.04 ANOMALİSİ SEGMENTASYON ALETİ Mİ? (25 Ağustos)
==========================================================================
AÇIK KAYIT (101h, 102 "dürüst kayıtlar"): χ₇'nin D_spont(τ=0.04) değeri
band daraldıkça YÜKSELİYOR — ±0.02L 0.078 → ±0.01L 0.106 → ±0.005L
0.125 — ve orada HİÇBİR pencere ilk çizgiyi kapsamıyor (band kaçağı
kapısı 0/6, üç bantta da). Öteki dört ailede aynı τ'da düşüş/durağanlık
var (ζ .041→.020, χ₃ .069→.018, χ₅ .030→.034, β .010→.013).

χ₇ ADASININ YAPISAL FARKI (bu koşudan önce ölçüldü, ön-mühür verisi):
  aile   pencere   uzunluklar (sıfır sayısı)                 L aralığı
  χ₃      3        6720, 6102, 63549                         7.94-9.65
  χ₅      4        3491, 8490, 49235, 7244                   8.16-10.50
  χ₇      6        11539, 5589, 4447, 25033, 3610, 9431      8.91-10.68
χ₇ EN ÇOK PENCERELİ ve EN UZUN PENCERESİ EN KÜÇÜK aile (25033 vs
63549/49235); ayrıca τ=0.04'te pencere-başına ω merkezi (τ·L_w)
0.357-0.427 ile ailelerin EN YÜKSEĞİ (χ₃ 0.317-0.386).

HİPOTEZ H-seg: kısa/çok parçalı segmentasyonda estimatörün tabanı band
daraldıkça yükselir; anomali χ₇'ye değil, χ₇'nin PARÇA YAPISINA aittir.

ÖN-MÜHÜR (ölçümden ÖNCE yazıldı) — H-seg doğruysa:
  S1  χ₃ ve χ₅'in sertifikalı sıfırları YAPAY OLARAK χ₇'nin segment-
      uzunluk profiline ([11539,5589,4447,25033,3610,9431]) bölünürse,
      onlar da τ=0.04'te band daraldıkça YÜKSELECEK (D(±0.005L) >
      D(±0.02L), en az ×1.3). (χ₅'in havuzu profili birebir almazsa
      profil ölçeklenir; ölçek raporlanır.)
  S1' EŞİT-PARÇA MERDİVENİ (asıl ayrıştırıcı, ortak taban): her ailenin
      EN UZUN penceresi ortak n* = 24000'e kırpılır (χ₃ 63549, χ₅ 49235,
      χ₇ 25033 — üçü de kaldırır) ve m = 1,2,4,6,8 EŞİT parçaya
      bölünür. Toplam n, aile, t-aralığı sabit; yalnız PARÇA SAYISI
      değişir. H-seg doğruysa D(±0.005L)/D(±0.02L) oranı m ile
      BÜYÜYECEK — üç ailede de. m=1'de (hiç bölünmemiş tek parça)
      χ₇ yükselişi HÂLÂ varsa H-seg ORADA ÖLÜR.
  S2  Aynı toplam n'i AZ ve UZUN parçaya bölen kontrol YÜKSELMEYECEK
      (düşecek ya da düz kalacak) — yani etki "veri azlığı" değil
      "parça kısalığı/çokluğu".
  S3  χ₇'nin PENCERE-BAŞINA ayrıştırımında yükseliş KISA pencerelerde
      olacak, en uzun penceresinde (25033) olmayacak.
  S4  Vekil taban V(n, band) haritası: V hem n küçüldükçe hem band
      daraldıkça YÜKSELECEK; donuk-D'nin anlamlılık sınırı D/V oranıdır.
  S5  AYRIŞTIRICI (H-seg'in alternatifi "L-yayılımı/hizalama"): tüm
      pencereleri ORTAK ω₀ = 0.04·⟨L⟩ merkezinde, ORTAK ω-bandıyla
      ölçmek (τ-bandı yerine) yükselişi ÖLDÜRÜRSE suç segmentasyonda
      değil pencere-L yayılımındadır. Bu ayrı bir ALET hükmüdür.
  Eş zamanlı 6-eşit-parça kontrolü (aynı n, aynı parça sayısı, eşit
  uzunluk) profil-biçiminin mi yoksa yalnız parça sayısının mı
  önemli olduğunu ayırır.

KOD-EŞDEĞERLİK KAPISI G0: χ₇ doğal pencerelerle τ=0.04, üç bant
101h'nin 0.078 / 0.106 / 0.125'ini BİREBİR vermeli.

Makine 101d/101h/101i'den BİREBİR kopya (pencere_hazirla, egri,
egri_omega, duz_segmentler); yalnız "parçalama" yardımcıları eklendi.

==========================================================================
SONUÇ (25 Ağustos, koşu 38 s) — S1 RET, S1' RET, S2 boş, S3 RET, S4 ✓,
S5 RET  →  H-seg REDDEDİLDİ. Ve ANOMALİNİN TANIMI DEĞİŞTİ (aşağıda).
==========================================================================
G0 ✓ BİREBİR: χ₇ doğal 0.0783 / 0.1062 / 0.1253 (101h: .078/.106/.125).

S1 RET — ama beklenmedik biçimde: χ₃'ü χ₇ profiline bölmek YÜKSELİŞ
ÜRETİYOR, fakat χ₃'ü BAŞKA HERHANGİ BİR biçimde bölmek de üretiyor:
  düzenek            npar  n_top   ±0.02L ±0.01L ±0.005L  oran
  chi3/doğal           3   76371   0.0685 0.0375 0.0178   0.26×
  chi3/χ₇profil        6   59649   0.0398 0.0466 0.0482   1.21×
  chi3/2uzun           2   59649   0.0544 0.0613 0.0622   1.14×
  chi3/6eşit           6   59649   0.0399 0.0439 0.0453   1.13×
  chi5/doğal           4   68460   0.0297 0.0327 0.0340   1.14×
  chi5/χ₇profil        6   56667   0.0394 0.0440 0.0453   1.15×
  chi5/2uzun           2   47719   0.0439 0.0504 0.0538   1.23×
  chi5/6eşit           6   50702   0.0289 0.0345 0.0357   1.23×
  chi7/doğal           6   59649   0.0783 0.1062 0.1253   1.60×
  χ₇ profili (1.21×) ile 2-uzun-parça (1.14×) ve 6-eşit (1.13×)
  AYIRT EDİLEMEZ. Parça profili bilgi taşımıyor → S1 RET.
  (χ₅'in havuzu profilleri birebir almadı; ölçekler 0.95/0.80/0.85
  olarak raporlandı — nicel karşılaştırma bu yüzden aile-içi yapıldı.)

S1' RET — EŞİT-PARÇA MERDİVENİ (ortak n*=24000, tek uzun pencereden):
  aile  m=1    m=2    m=4    m=6    m=8   (oran dar/geniş)
  χ₃   1.28×  1.27×  1.21×  1.05×  1.21×
  χ₅   1.23×  1.21×  1.21×  0.97×  1.11×
  χ₇   1.34×  1.07×  1.05×  0.68×  1.06×
  PARÇA SAYISIYLA HİÇBİR EĞİLİM YOK; m=1'de (hiç bölünmemiş) oran en
  büyük. Segmentasyon yükselişi ÜRETMİYOR, hatta hafifçe SÖNDÜRÜYOR.

S2 boş — eş-n uzun kontrol de yükseldiği için ayrıştırıcı olmadı.

S3 RET + ASIL BULGU: χ₇ PENCERE-BAŞINA (τ=0.04):
  pencere      n       L     ±0.02L ±0.01L ±0.005L  oran   V(±0.005L)
  w1       11539   8.912     0.3036 0.3956 0.4451   1.47×   0.0437
  w2        5589   9.570     0.0578 0.0648 0.0669   1.16×   0.0095
  w3        4447   9.820     0.0200 0.0206 0.0212   1.06×   0.0046
  w4       25033  10.261     0.0202 0.0260 0.0289   1.43×   0.0077
  w5        3610  10.572     0.1158 0.1335 0.1377   1.19×   0.0119
  w6        9431  10.684     0.0344 0.0296 0.0301   0.87×   0.0026
  w1 (t ∈ [2915,11063], EN DÜŞÜK L, EN DÜŞÜK t) τ=0.04'te D ≈ 0.30-0.45
  — yani PRATİKTE HİÇ DONMAMIŞ; öteki beş pencere 0.020-0.138 ile
  donuk. Havuzlanmış 0.078-0.125'in tamamı bu tek pencerenin
  seyreltilmiş hâlidir. Kısa pencereler (3610/4447/5589) düz.
  → Yükseliş "kısa parça" imzası DEĞİL; TEK PENCERE imzası.
  Karşılaştırma: χ₃-w1 (n=6720, L=7.935) tersine 0.160→0.032 (0.20×)
  — χ₃'ün doğal 0.26×'i de tek pencereden geliyor. Yani hem χ₇'nin
  "yükselişi" hem χ₃'ün "düşüşü" AYNI CİNSTEN: en düşük-L, en düşük-t
  penceresinin tekil davranışı.

S4 ✓ VEKİL TABAN HARİTASI (tek pencere, τ=0.04, n taraması):
  kaynak    n      L      V(±0.02/0.01/0.005L)   D(±0.02/0.01/0.005L)
   χ₃     3000   8.850   0.0021 0.0017 0.0019   0.0171 0.0211 0.0222
   χ₃     6000   8.918   0.0030 0.0028 0.0033   0.0255 0.0297 0.0307
   χ₃    12000   9.036   0.0028 0.0025 0.0027   0.0263 0.0318 0.0332
   χ₃    25000   9.242   0.0030 0.0024 0.0024   0.0108 0.0119 0.0119
   χ₃    50000   9.528   0.0035 0.0028 0.0028   0.0216 0.0261 0.0273
   χ₇     3000   9.973   0.0057 0.0049 0.0055   0.0551 0.0647 0.0668
   χ₇     6000  10.020   0.0067 0.0049 0.0057   0.0671 0.0729 0.0742
   χ₇    12000  10.105   0.0068 0.0058 0.0053   0.0434 0.0452 0.0460
   χ₇    25000  10.260   0.0050 0.0047 0.0041   0.0379 0.0444 0.0466
  KURAL 1: vekil taban V BAND'DAN BAĞIMSIZ (±%15 dalgalanma, eğilim
    yok) ve n'DEN de pratikte BAĞIMSIZ (χ₃ 0.002-0.0035; χ₇ 0.004-0.007
    tüm n'lerde) — V ~ 1/√n DEĞİL, çünkü havuz tek pencere ve
    permütasyon frekanslar arası yapılıyor. V AİLEYE bağlı.
  KURAL 2: temiz veride bile D band daraldıkça 1.1-1.3× YÜKSELİR
    (9 satırın 9'unda). YANİ "BAND DARALINCA YÜKSELME" ESTİMATÖRÜN
    NORMAL DAVRANIŞIDIR, anomali değildir. Anomali olan χ₇'nin
    1.60'ı değil, χ₃'ün DOĞAL 0.26'sıdır (aşağıya bak).

S5 RET — ORTAK ω₀ + ORTAK ω-BANDI (τ↔ω hizalama ayrıştırıcısı):
  χ₇ ω₀=0.3995: D 0.0809 → 0.1050 → 0.1121   oran 1.39×
  χ₃ ω₀=0.3766: D 0.0832 → 0.0598 → 0.0517   oran 0.62×
  χ₅ ω₀=0.3917: D 0.0320 → 0.0349 → 0.0354   oran 1.11×
  Pencere-L yayılımı / τ↔ω çevrimi SORUMLU DEĞİL; sıralama korunuyor.

HÜKÜM (H-seg): KESİN RET. Anomali segmentasyonun eseri değil — parça
uzunluğu, parça sayısı, parça profili, eş-n veri miktarı ve τ↔ω
hizalamasının hiçbiri onu üretmiyor ya da söndürmüyor.
ANOMALİNİN YENİDEN TANIMI (104a'nın asıl kazancı): "band daraldıkça
yükselme" 20+ düzenekte GENEL DAVRANIŞ (1.1-1.3×). χ₇ bu ölçeğin
üstünde (1.60×) ama tek bir pencereden (w1, D≈0.30-0.45 = donmamış)
geliyor; χ₃ ise ters yönde (0.26×) ve o da tek bir pencereden
(χ₃-w1, 0.20×). Sıradaki iki soru 104b ve 104c'ye devredildi:
w1'de KUSUR mu var (104b: pencere-içi sayım denetimi + w1-dışlama),
yoksa χ₇'nin donuk bölgesinde FREKANS yapısı mı var (104c)?
"""

import numpy as np
import time
from pathlib import Path

T0 = time.time()
HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
rng = np.random.default_rng(1040)

PKS = [2,3,4,5,7,8,9,11,13,16,17,19,23,25,27,29,31,32,37,41,43,47,49,53,
       59,61,64,67,71,73,79,81,83,89,97,101,103,107,109,113,121,125,127,128]
LINES_ALL = [np.log(qq) for qq in PKS]

# ---------------- 101d/101h/101i makinesi (birebir kopya) ----------------

def pencere_hazirla(z, qeff):
    gaps = np.diff(z)
    mids = 0.5 * (z[:-1] + z[1:])
    Lw = float(np.log(qeff * mids / TWO_PI).mean())
    ds = gaps * np.log(qeff * mids / TWO_PI) / TWO_PI - 1
    tt = (mids - mids.mean()) / (mids[-1] - mids[0])
    P = np.vstack([np.ones_like(tt), tt, tt**2, tt**3]).T
    ds = ds - P @ np.linalg.lstsq(P, ds, rcond=None)[0]
    return (z, mids, ds, Lw)

def egri(WIN, taus, band=0.02):
    """101h'nin egri'si: τ-bandı ±band·L_w, pencereler-arası havuz."""
    Dv, Vv, Nf = [], [], []
    for tau0 in taus:
        num = 0.0 + 0j; den = 0.0
        Gs, rs = [], []
        nf_tot, nf_ham = 0, 0
        for (zz, tm, ds, Lw) in WIN:
            oms = tau0 * Lw + np.linspace(-band * Lw, band * Lw, 160)
            nf_ham += len(oms)
            oms = np.array([o for o in oms
                            if min(abs(o - l) for l in LINES_ALL) > 0.01])
            nf_tot += len(oms)
            for s0 in range(0, len(oms), 40):
                ob = oms[s0:s0 + 40]
                rr = np.exp(1j * np.outer(ob, zz)).sum(axis=1)
                GG = (np.exp(1j * np.outer(ob, tm)) * ds[None, :]).sum(axis=1)
                num += (GG * np.conj(rr)).sum()
                den += (np.abs(rr)**2).sum()
                Gs.append(GG); rs.append(rr)
        kap = 2 * np.pi * tau0
        c = 2 * np.sin(kap / 2) / kap
        G_all, r_all = np.concatenate(Gs), np.concatenate(rs)
        fl = []
        for _ in range(40):
            perm = rng.permutation(len(r_all))
            fl.append(abs((G_all * np.conj(r_all[perm])).sum()) / (c * den))
        Dv.append(abs(num) / (c * den)); Vv.append(float(np.mean(fl)))
        Nf.append(nf_tot / max(nf_ham, 1))
    return Dv, Vv, Nf

def egri_omega(WIN, omegas, band=0.02, nfrek=160, cizgi_filtre=True):
    """101i: sabit-ω, ORTAK ω-bandı; c-faktörü pencere başına."""
    Dv, Vv = [], []
    for om in omegas:
        num = 0.0 + 0j; den = 0.0
        Gs, rs = [], []
        for (zz, tm, ds, Lw) in WIN:
            tau_w = om / Lw
            kap = 2 * np.pi * tau_w
            c = 2 * np.sin(kap / 2) / kap
            oms = om + np.linspace(-band, band, nfrek)
            if cizgi_filtre:
                oms = np.array([o for o in oms
                                if min(abs(o - l) for l in LINES_ALL) > 0.01])
            for s0 in range(0, len(oms), 40):
                ob = oms[s0:s0 + 40]
                rr = np.exp(1j * np.outer(ob, zz)).sum(axis=1)
                GG = (np.exp(1j * np.outer(ob, tm)) * ds[None, :]).sum(axis=1)
                GG = GG / c
                num += (GG * np.conj(rr)).sum()
                den += (np.abs(rr)**2).sum()
                Gs.append(GG); rs.append(rr)
        G_all, r_all = np.concatenate(Gs), np.concatenate(rs)
        fl = []
        for _ in range(40):
            perm = rng.permutation(len(r_all))
            fl.append(abs((G_all * np.conj(r_all[perm])).sum()) / den)
        Dv.append(abs(num) / den); Vv.append(float(np.mean(fl)))
    return Dv, Vv

def duz_segmentler(zz, sayim, min_n=3000, win=80, esik=0.5, pad=160):
    d = np.arange(len(zz)) - (sayim(zz) - sayim(zz[0]))
    from numpy.lib.stride_tricks import sliding_window_view
    med = np.median(sliding_window_view(d, win), axis=1)
    adim = med[win + 1:] - med[:-(win + 1)]
    bad = np.zeros(len(zz), bool)
    for i in np.where(np.abs(adim) > esik)[0] + win // 2:
        bad[max(0, i - pad):i + 2 * pad] = True
    seg, kes = [], 0
    s0 = 0
    for i in range(1, len(zz) + 1):
        if i == len(zz) or bad[i] != bad[i - 1]:
            if not bad[s0] and i - s0 >= min_n:
                seg.append(zz[s0:i])
            if i < len(zz) and bad[i]:
                kes += 1
            s0 = i
    return seg, kes

# ------------------------- doğal pencereler -------------------------

exec(open(HERE / "98_L_motoru.py").read().split('CHI4 = ')[0])
HERE = Path(__file__).resolve().parent
CHI4 = {0: 0, 1: 1, 2: 0, 3: -1}
CHI3 = {0: 0, 1: 1, 2: -1}
CHI5 = {0: 0, 1: 1, 2: 1j, 3: -1j, 4: -1}
z6c = np.exp(1j * np.pi / 3)
CHI7 = {0: 0, 1: 1, 3: z6c, 2: z6c**2, 6: z6c**3, 4: z6c**4, 5: z6c**5}
AILE = {"chi3": (3, CHI3), "chi5": (5, CHI5), "chi7": (7, CHI7)}

HAM, PEN, SEGHAM = {}, {}, {}
for etiket, (q, tab) in AILE.items():
    zc = np.load(HERE / f"101f_{etiket}_zeros.npz")["zeros"]
    lo = np.exp(np.log(zc[0] + 1) + 0.25 * (np.log(zc[-1]) - np.log(zc[0] + 1)))
    zc = zc[zc >= lo]
    M = Lmotor(q, tab, 1)
    seg, kes = duz_segmentler(zc, lambda t: M.theta(t) / np.pi)
    SEGHAM[etiket] = seg
    kenar = np.exp(np.linspace(np.log(zc[0]), np.log(zc[-1] * 1.0001), 4))
    WIN, PARCA = [], []
    for s in seg:
        for i in range(3):
            p = s[(s >= kenar[i]) & (s < kenar[i + 1])]
            if len(p) >= 3000:
                WIN.append(pencere_hazirla(p, float(q)))
                PARCA.append(p)
    PEN[etiket] = WIN
    HAM[etiket] = PARCA
    ns = [len(w[0]) for w in WIN]
    print(f"[{etiket}] {len(WIN)} pencere n={ns} toplam={sum(ns)} "
          f"L={[round(w[3],3) for w in WIN]}", flush=True)

PROFIL = [len(w[0]) for w in PEN["chi7"]]
NTOP = sum(PROFIL)
print(f"\nχ₇ PROFİLİ: {PROFIL}  toplam={NTOP}", flush=True)

# ------------------------- yapay parçalama -------------------------

def paketle(havuz, uzunluklar):
    """FFD: her havuz parçasında sığan EN UZUN bloğu ardışık yerleştir.
    Bir blok tek bir havuz parçası içinden çıkar (kusur sınırını aşmaz).
    Sığmazsa None döner."""
    out, kalan = [], sorted(uzunluklar, reverse=True)
    for h in sorted(havuz, key=len, reverse=True):
        off = 0
        while kalan:
            se = next((k for k in kalan if off + k <= len(h)), None)
            if se is None:
                break
            kalan.remove(se)
            out.append(h[off:off + se])
            off += se
    return None if kalan else out

def parcala(etiket, uzunluklar, q, ad=""):
    """Profili sığdır; sığmazsa ölçekleyerek küçült (ölçek raporlanır)."""
    for s in np.arange(1.0, 0.39, -0.05):
        uz = [max(3000, int(round(u * s))) for u in uzunluklar]
        out = paketle(HAM[etiket], uz)
        if out is not None:
            if s < 0.999:
                print(f"    [{etiket}/{ad}] profil ölçeği {s:.2f} "
                      f"(uzunluklar {sorted([len(p) for p in out], reverse=True)})",
                      flush=True)
            return [pencere_hazirla(p, float(q)) for p in out]
    raise RuntimeError(f"{etiket}/{ad}: profil hiçbir ölçekte sığmadı")

TAUS = [0.04, 0.05, 0.06]
BANDS = [0.02, 0.01, 0.005]

DUZENEK = {}
DUZENEK[("chi3", "doğal")] = PEN["chi3"]
DUZENEK[("chi5", "doğal")] = PEN["chi5"]
DUZENEK[("chi7", "doğal")] = PEN["chi7"]
for et, q in [("chi3", 3), ("chi5", 5)]:
    DUZENEK[(et, "χ₇profil")] = parcala(et, PROFIL, q, "χ₇profil")
    DUZENEK[(et, "2uzun")] = parcala(et, [NTOP // 2, NTOP - NTOP // 2], q, "2uzun")
    e6 = [NTOP // 6] * 5 + [NTOP - 5 * (NTOP // 6)]
    DUZENEK[(et, "6eşit")] = parcala(et, e6, q, "6eşit")

print("\n=== S1/S2: YAPAY PARÇALAMA, D(τ=0.04) ÜÇ BANTTA ===")
print(f"{'düzenek':>22} {'npar':>5} {'n_top':>7} "
      + " ".join(f"±{b}L" for b in BANDS) + "   oran(dar/geniş)   L aralığı")
SON = {}
SIRA = [("chi3", "doğal"), ("chi3", "χ₇profil"), ("chi3", "2uzun"),
        ("chi3", "6eşit"), ("chi5", "doğal"), ("chi5", "χ₇profil"),
        ("chi5", "2uzun"), ("chi5", "6eşit"), ("chi7", "doğal")]
for key in SIRA:
    W = DUZENEK[key]
    vals, vek = [], []
    for b in BANDS:
        D, V, _ = egri(W, TAUS, band=b)
        SON[(key, b)] = (D, V)
        vals.append(D[0]); vek.append(V[0])
    Ls = [w[3] for w in W]
    print(f"{key[0]+'/'+key[1]:>22} {len(W):>5} {sum(len(w[0]) for w in W):>7} "
          + " ".join(f"{v:>6.4f}" for v in vals)
          + f"   {vals[-1]/vals[0]:>6.2f}×    "
          + f"[{min(Ls):.2f},{max(Ls):.2f}]", flush=True)
print("  (vekil tabanlar τ=0.04)")
for key in SIRA:
    print(f"{key[0]+'/'+key[1]:>22}       " +
          " ".join(f"{SON[(key,b)][1][0]:>6.4f}" for b in BANDS))

print("\n  aynı düzenekler τ=0.05 ve τ=0.06:")
for j, t in enumerate(TAUS[1:], start=1):
    print(f"  --- τ={t} ---")
    for key in SIRA:
        print(f"{key[0]+'/'+key[1]:>22}       " +
              " ".join(f"{SON[(key,b)][0][j]:>6.4f}" for b in BANDS))

# ---------------- S1': EŞİT-PARÇA MERDİVENİ (ortak n*) ----------------
NYILDIZ = 24000
print(f"\n=== S1': EŞİT-PARÇA MERDİVENİ — ortak n*={NYILDIZ}, "
      f"her ailenin EN UZUN penceresinden ===")
print(f"{'aile':>6} {'m':>3} {'parça n':>8} " + " ".join(f"±{b}L" for b in BANDS)
      + "   oran   " + " ".join(f"V±{b}L" for b in BANDS))
MERD = {}
for et, q in [("chi3", 3), ("chi5", 5), ("chi7", 7)]:
    kaynak = max(HAM[et], key=len)[:NYILDIZ]
    for m in [1, 2, 4, 6, 8]:
        kesim = np.array_split(kaynak, m)
        W = [pencere_hazirla(p, float(q)) for p in kesim]
        Ds, Vs = [], []
        for b in BANDS:
            D, V, _ = egri(W, [0.04], band=b)
            Ds.append(D[0]); Vs.append(V[0])
        MERD[(et, m)] = (Ds, Vs)
        print(f"{et:>6} {m:>3} {len(kesim[0]):>8} "
              + " ".join(f"{d:>6.4f}" for d in Ds)
              + f"  {Ds[-1]/Ds[0]:>5.2f}× "
              + " ".join(f"{v:>6.4f}" for v in Vs), flush=True)

# ------------------------- S3: pencere-başına -------------------------
print("\n=== S3: χ₇ PENCERE-BAŞINA D(τ=0.04) ===")
print(f"{'pencere':>8} {'n':>7} {'L':>7} " + " ".join(f"±{b}L" for b in BANDS)
      + "   oran   " + " ".join(f"V±{b}L" for b in BANDS))
for i, w in enumerate(PEN["chi7"]):
    vals, veks = [], []
    for b in BANDS:
        D, V, _ = egri([w], [0.04], band=b)
        vals.append(D[0]); veks.append(V[0])
    print(f"{'w'+str(i+1):>8} {len(w[0]):>7} {w[3]:>7.3f} "
          + " ".join(f"{v:>6.4f}" for v in vals)
          + f"  {vals[-1]/vals[0]:>5.2f}×  "
          + " ".join(f"{v:>6.4f}" for v in veks), flush=True)
print("  karşılaştırma — χ₃ ve χ₅ pencere-başına:")
for et in ["chi3", "chi5"]:
    for i, w in enumerate(PEN[et]):
        vals = [egri([w], [0.04], band=b)[0][0] for b in BANDS]
        print(f"{et+'-w'+str(i+1):>8} {len(w[0]):>7} {w[3]:>7.3f} "
              + " ".join(f"{v:>6.4f}" for v in vals)
              + f"  {vals[-1]/vals[0]:>5.2f}×", flush=True)

# ------------------------- S4: taban haritası -------------------------
print("\n=== S4: VEKİL TABAN HARİTASI V(n, band) ve D(n, band), τ=0.04 ===")
print("  (kaynak: χ₃'ün en uzun sertifikalı parçasının baş kısmı)")
uzun3 = max(HAM["chi3"], key=len)
uzun7 = max(HAM["chi7"], key=len)
print(f"{'kaynak':>8} {'n':>7} {'L':>7} " + " ".join(f"V±{b}L" for b in BANDS)
      + "  |  " + " ".join(f"D±{b}L" for b in BANDS) + "   V·√n")
HARITA = []
for kaynak, arr, q in [("χ₃", uzun3, 3), ("χ₇", uzun7, 7)]:
    for n in [3000, 6000, 12000, 25000, 50000]:
        if n > len(arr):
            continue
        w = pencere_hazirla(arr[:n], float(q))
        Ds, Vs = [], []
        for b in BANDS:
            D, V, _ = egri([w], [0.04], band=b)
            Ds.append(D[0]); Vs.append(V[0])
        HARITA.append((kaynak, n, w[3], Vs, Ds))
        print(f"{kaynak:>8} {n:>7} {w[3]:>7.3f} "
              + " ".join(f"{v:>6.4f}" for v in Vs) + "  |  "
              + " ".join(f"{d:>6.4f}" for d in Ds)
              + f"   {np.mean(Vs)*np.sqrt(n):>6.3f}", flush=True)

# ------------------------- S5: ortak ω-bandı -------------------------
print("\n=== S5: ORTAK ω₀ ve ORTAK ω-BANDI (τ↔ω hizalama ayrıştırıcısı) ===")
for et in ["chi7", "chi3", "chi5"]:
    W = PEN[et]
    ns = np.array([len(w[0]) for w in W], float)
    Lort = float(np.average([w[3] for w in W], weights=ns))
    om0 = 0.04 * Lort
    print(f"  {et}: ⟨L⟩_n={Lort:.4f}  ω₀={om0:.4f}")
    vals = []
    for b in BANDS:
        D, V = egri_omega(W, [om0], band=b * Lort)
        vals.append(D[0])
        print(f"     Δω={b*Lort:.4f} (≡{b}L)   D={D[0]:.4f}   V={V[0]:.4f}",
              flush=True)
    print(f"     oran(dar/geniş) = {vals[-1]/vals[0]:.2f}×")

np.savez(HERE / "104a_segment.npz",
         profil=np.array(PROFIL), bands=np.array(BANDS), taus=np.array(TAUS),
         **{f"D_{k[0]}_{k[1]}_{b}": np.array(SON[(k, b)][0])
            for k in SIRA for b in BANDS},
         **{f"V_{k[0]}_{k[1]}_{b}": np.array(SON[(k, b)][1])
            for k in SIRA for b in BANDS})
print(f"\n104a_segment.npz yazıldı. Süre {time.time()-T0:.0f} s.")
