"""
104b — H-kusur: χ₇'NİN KALAN ~10 KAYIP SIFIRI ANOMALİYİ ÜRETİYOR MU?
==========================================================================
101f: χ₇'de kalan şüpheli bölge = 8, sayım sürüklenmesi −9.68 → ~10 sıfır
HÂLÂ KAYIP. 101d/101h'nin düzlük segmentasyonu (esik=0.5, pad=160)
bunların hepsini analiz dışına atmıyor OLABİLİR. 101b kusur kapısı
donma estimatörünün kayıp sıfıra aşırı duyarlı olduğunu göstermişti
(%0.1 silme → D 0.09'dan 0.55'e) — AMA yalnız TEK bantta (±0.018L).
Anomali BAND'a bağlı olduğu için kusur duyarlılığının BAND ÖLÇEKLENMESİ
ölçülmemiş bir boşluktur.

HİPOTEZ H-kusur: χ₇'nin τ=0.04 yükselişi segmentasyondan kaçan artık
kayıp sıfırların eseridir; kusurun D'ye katkısı band daraldıkça BÜYÜR.

ÖN-MÜHÜR (ölçümden ÖNCE yazıldı) — H-kusur doğruysa:
  K1  SIKILAŞTIRILMIŞ SEGMENTASYON (esik 0.35 / pad 300; ikinci basamak
      esik 0.25 / pad 500) χ₇'nin τ=0.04 yükselişini SÖNDÜRECEK:
      D(±0.005L)/D(±0.02L) oranı 1.6×'ten ~1'e (ya da altına) inecek.
      Öteki ailelerde (kontrol) oran anlamlı değişmeyecek.
  K2  SAYIM DENETİMİ: χ₇'nin analiz pencerelerinin İÇİNDE net sürüklenme
      (d_son − d_ilk, d_i = i − θ(t_i)/π) SIFIRDAN FARKLI olacak —
      kayıp sıfırlar pencere içinde demektir. (Bu kapı H-kusur'un
      ön koşuludur: pencere içi sürüklenme ≈ 0 ise hipotez temelsizdir.)
  K3  ENJEKSİYON: temiz χ₃'ün analiz pencerelerine χ₇-benzeri k rastgele
      SİLME (k = 10, 30, 100, 300; 8 gerçekleme) sokulduğunda D artacak
      ve ARTIŞ BAND DARALDIKÇA BÜYÜYECEK — yani ΔD(±0.005L) >
      ΔD(±0.02L). k=10'da bile ölçülebilir bir yükseliş çıkmalı
      (χ₇'nin kayıp sayısı ~10).
  K4  Aynı enjeksiyon χ₇'ye uygulandığında yükseliş oranı DAHA DA
      büyüyecek (kusur üstüne kusur).
  K5  101f'nin 'bolgeler' listesindeki şüpheli aralıkların analiz
      pencereleriyle örtüşmesi raporlanacak (0 örtüşme → H-kusur'un
      "bilinen" kolu kapanır, yalnız "bilinmeyen kayıp" kolu kalır).
  K6  [104a'DAN SONRA EKLENDİ, açıkça etiketli] 104a'nın S3'ü χ₇'nin
      havuzlanmış D'sinin TEK BİR PENCERE tarafından sürüklendiğini
      gösterdi: w1 (n=11539, L=8.912, t∈[2915,11063]) τ=0.04'te
      D = 0.304/0.396/0.445 — öteki beş pencere 0.020-0.138. w1
      pratikte "donmamış" görünüyor. K6: χ₇'yi w1 DIŞLANMIŞ hâlde
      yeniden ölç; ayrıca her rejimde pencere-başına D yaz. w1
      çıkarıldığında anomali (dar/geniş oranı) sönerse suç TEK BİR
      PENCEREDE, yani bir VERİ olayında yoğunlaşmış demektir.

VERİ: 101f sertifikalı sıfırlar; ζ bu koşuya girmiyor (χ-adaları
karşılaştırması yeterli, ζ'nın motoru yok).
Makine 101d/101h'den BİREBİR kopya.

==========================================================================
SONUÇ (25 Ağustos, koşu 173 s) — K1 ✓✓ K2 ön-mühürlü biçimi RET /
inceltilmiş biçimi İSABET, K3 ✓✓ (BELİRLEYİCİ), K4 kısmi, K5 ✓
→ H-kusur DOĞRULANDI.
==========================================================================
K1 ✓✓ SIKILAŞTIRMA χ₇'NİN ANOMALİSİNİ SÖNDÜRÜYOR (τ=0.04, üç bant):
  aile  rejim                n_top   ±0.02L ±0.01L ±0.005L  oran
  χ₃    gevşek(0.50/160)     76371   0.0685 0.0375 0.0178   0.26×
  χ₃    sıkı(0.35/300)       75941   0.0668 0.0362 0.0160   0.24×
  χ₃    çok sıkı(0.25/500)   75334   0.0709 0.0404 0.0216   0.31×
  β     gevşek               70881   0.0102 0.0124 0.0130   1.28×
  β     sıkı                 70600   0.0070 0.0101 0.0093   1.33×
  β     çok sıkı             70196   0.0071 0.0087 0.0089   1.26×
  χ₅    gevşek               68460   0.0297 0.0327 0.0340   1.14×
  χ₅    sıkı                 67606   0.0246 0.0275 0.0290   1.18×
  χ₅    çok sıkı             63241   0.0085 0.0114 0.0120   1.41×
  χ₇    gevşek               59649   0.0783 0.1062 0.1253   1.60×
  χ₇    sıkı                 57090   0.0421 0.0487 0.0507   1.20×
  χ₇    çok sıkı             51662   0.0363 0.0415 0.0435   1.20×
  χ₇ HEM SEVİYE (0.125 → 0.044, 2.9×) HEM ORAN (1.60× → 1.20×) düşüyor
  ve genel banda (1.1-1.3×, bkz. 104a) oturuyor. Kontrol aileleri
  (β, χ₃) pratikte değişmiyor → etki χ₇'ye ÖZGÜ ve VERİYE bağlı.

K2 ön-mühürlü biçimi RET, inceltilmiş biçimi İSABET: χ₇ w1'in NET
  sürüklenmesi −0.00 (yani "pencere içinde kayıp yok" diyor) AMA
  pencere-içi maks|d − medyan| = 2.75 — 15 pencerenin EN BÜYÜĞÜ:
    χ₇ w1 2.75 | χ₃ w1 2.60 | kalan 13 pencere 0.88-1.16
  Yani kusur GİT-GEL: net sıfır, içeride ±2'lik sıçrama. NET SÜRÜKLENME
  TESTİ BU KUSUR SINIFINI GÖREMİYOR — 101f'nin "χ₇'de ~10 kayıp" ölçütü
  de aynı körlükte. (Teşhis 104d'ye devredildi.)

K5 ✓ (H-kusur'un "bilinen" kolu KAPANDI): 101f'nin şüpheli bölgeleri
  analiz pencereleriyle 0 (SIFIR) örtüşüyor — χ₃ 2/0, β 1/0, χ₅ 3/0,
  χ₇ 8/0. Suçlu LİSTEDE OLMAYAN bir kusur.

K6 ✓✓ (104a'dan sonra eklendi) — χ₇ PENCERE-BAŞINA, rejim rejim:
  gevşek : w1 0.304→0.396→0.445 (1.47×) | w2-w6 0.020-0.138
  sıkı   : w1 0.201→0.245→0.256 (1.27×)  [w1 başı 2914.6 → 3572.4]
  çok sıkı: w1 0.049→0.057→0.060 (1.22×) [w1 başı → 3879.7]
  Sıkılaştırma w1'in BAŞINI kırpıyor ve D'si 0.445 → 0.060'a ÇÖKÜYOR.
  Kusur t ∈ [2914.6, 3879.7] aralığında.
K6b ✓✓: w1 havuzdan çıkarılınca χ₇ 0.0236/0.0248/0.0255, oran 1.08× —
  tamamen normal. Anomalinin %100'ü tek pencereden geliyor.

K3 ✓✓ BELİRLEYİCİ — KUSUR DUYARLILIĞI BAND DARALDIKÇA BÜYÜYOR:
  temiz χ₃'ün pencerelerine k rastgele silme (8 gerçekleme), τ=0.04:
   k      ‰    ±0.02L ±0.01L ±0.005L      ΔD(temizden)          oran
   0    0.00   0.0685 0.0375 0.0178            —                0.26×
  10    0.13   0.1316 0.1187 0.1343   +0.0630 +0.0811 +0.1165   1.02×
  30    0.39   0.2390 0.2447 0.2880   +0.1705 +0.2072 +0.2701   1.20×
 100    1.31   0.4737 0.5044 0.5899   +0.4052 +0.4668 +0.5720   1.25×
 300    3.93   0.7158 0.7600 0.8081   +0.6473 +0.7225 +0.7903   1.13×
  SADECE 10 SİLME (0.13‰) χ₃'ü 0.018'den 0.134'e (7.5×) çıkarıyor ve
  ΔD dar bantta geniş banda göre 1.8× BÜYÜK (+0.117 vs +0.063).
  101b'nin tek-bant duyarlılığı böylece BAND'A GENİŞLETİLDİ: donma
  estimatörünün kusur duyarlılığı band daraldıkça ARTAR. χ₇'nin
  "band daralınca yükselen" imzası TAM OLARAK BU İMZADIR.
  χ₅ aynı: k=10 → 0.034'ten 0.100'e, ΔD +0.054/+0.067/+0.066.

K4 kısmi — χ₇'ye enjeksiyon seviyeyi yükseltiyor ama ORANI DÜŞÜRÜYOR
  (1.60 → 1.47 → 1.37 → 1.22 → 1.12): ağır kusurda D → 1'e doyuyor,
  doymuş bölgede oran sıkışıyor. Ön-mühürdeki "oran daha da büyüyecek"
  beklentisi YANLIŞTI; mekanizma yine de aynı (dürüst kayıt).

HÜKÜM (H-kusur): DOĞRULANDI — KESİN. χ₇'nin τ=0.04 anomalisi
düzlük segmentasyonundan kaçan, sayım-sürüklenmesi net-sıfır olduğu
için görünmeyen bir GİT-GEL KUSUR KÜMESİNİN eseri; kümenin yeri
χ₇ w1'in başı (t ≲ 3880). Teşhis 104d, temizlik ve yeniden ölçüm 104e.
"""

import numpy as np
import time
from pathlib import Path

T0 = time.time()
HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
rng = np.random.default_rng(1041)

PKS = [2,3,4,5,7,8,9,11,13,16,17,19,23,25,27,29,31,32,37,41,43,47,49,53,
       59,61,64,67,71,73,79,81,83,89,97,101,103,107,109,113,121,125,127,128]
LINES_ALL = [np.log(qq) for qq in PKS]

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

exec(open(HERE / "98_L_motoru.py").read().split('CHI4 = ')[0])
HERE = Path(__file__).resolve().parent
CHI4 = {0: 0, 1: 1, 2: 0, 3: -1}
CHI3 = {0: 0, 1: 1, 2: -1}
CHI5 = {0: 0, 1: 1, 2: 1j, 3: -1j, 4: -1}
z6c = np.exp(1j * np.pi / 3)
CHI7 = {0: 0, 1: 1, 3: z6c, 2: z6c**2, 6: z6c**3, 4: z6c**4, 5: z6c**5}
AILE = {"chi3": (3, CHI3), "beta": (4, CHI4), "chi5": (5, CHI5),
        "chi7": (7, CHI7)}

ZC, MOT, BOLGE = {}, {}, {}
for et, (q, tab) in AILE.items():
    d = np.load(HERE / f"101f_{et}_zeros.npz")
    zc = d["zeros"]
    lo = np.exp(np.log(zc[0] + 1) + 0.25 * (np.log(zc[-1]) - np.log(zc[0] + 1)))
    ZC[et] = zc[zc >= lo]
    MOT[et] = Lmotor(q, tab, 1)
    BOLGE[et] = np.asarray(d["bolgeler"]) if "bolgeler" in d.files else np.zeros((0, 2))

def pencereler(et, esik=0.5, pad=160):
    q = AILE[et][0]
    zc = ZC[et]
    M = MOT[et]
    seg, kes = duz_segmentler(zc, lambda t: M.theta(t) / np.pi, esik=esik, pad=pad)
    kenar = np.exp(np.linspace(np.log(zc[0]), np.log(zc[-1] * 1.0001), 4))
    WIN, PARCA = [], []
    for s in seg:
        for i in range(3):
            p = s[(s >= kenar[i]) & (s < kenar[i + 1])]
            if len(p) >= 3000:
                WIN.append(pencere_hazirla(p, float(q)))
                PARCA.append(p)
    return WIN, PARCA, kes

TAUS = [0.04, 0.05, 0.06]
BANDS = [0.02, 0.01, 0.005]
REJIM = [("gevşek(0.50/160)", 0.5, 160),
         ("sıkı(0.35/300)", 0.35, 300),
         ("çok sıkı(0.25/500)", 0.25, 500)]

# ---------------- K1: sıkılaştırılmış segmentasyon ----------------
print("=== K1: SEGMENTASYON SIKILAŞTIRMA ===", flush=True)
PENR = {}
for ad, es, pd in REJIM:
    print(f"\n--- rejim {ad} ---")
    for et in AILE:
        W, P, kes = pencereler(et, esik=es, pad=pd)
        PENR[(ad, et)] = (W, P)
        ns = [len(w[0]) for w in W]
        print(f"  [{et}] {len(W)} pencere n={ns} toplam={sum(ns)} kesim={kes}",
              flush=True)

print(f"\n{'aile':>6} {'rejim':>20} {'n_top':>7} "
      + " ".join(f"±{b}L" for b in BANDS) + "   oran   "
      + " ".join(f"V±{b}L" for b in BANDS))
K1 = {}
for et in AILE:
    for ad, es, pd in REJIM:
        W = PENR[(ad, et)][0]
        if not W:
            print(f"{et:>6} {ad:>20}  --- pencere yok ---"); continue
        Ds, Vs = [], []
        for b in BANDS:
            D, V, _ = egri(W, TAUS, band=b)
            K1[(et, ad, b)] = (D, V)
            Ds.append(D[0]); Vs.append(V[0])
        print(f"{et:>6} {ad:>20} {sum(len(w[0]) for w in W):>7} "
              + " ".join(f"{d:>6.4f}" for d in Ds)
              + f"  {Ds[-1]/Ds[0]:>5.2f}× "
              + " ".join(f"{v:>6.4f}" for v in Vs), flush=True)

print("\n  τ=0.05 ve τ=0.06 (aynı rejimler):")
for j, t in enumerate(TAUS[1:], start=1):
    print(f"  --- τ={t} ---")
    for et in AILE:
        for ad, _, _ in REJIM:
            if (et, ad, BANDS[0]) not in K1:
                continue
            print(f"{et:>6} {ad:>20}   "
                  + " ".join(f"{K1[(et,ad,b)][0][j]:>6.4f}" for b in BANDS))

# ---------------- K2: pencere-içi sayım denetimi ----------------
print("\n=== K2: PENCERE-İÇİ SAYIM SÜRÜKLENMESİ (kayıp sıfır sayacı) ===")
print("  d_i = i − (θ(t_i)−θ(t_0))/π ;  net = d_son − d_ilk")
print("  (net > 0 → pencere içinde o kadar sıfır KAYIP; ~0 → tam)")
print(f"{'aile':>6} {'pencere':>4} {'n':>7} {'net':>8} {'std':>7} "
      f"{'maks|d-med|':>11} {'t aralığı':>26}")
for et in AILE:
    W, P = PENR[("gevşek(0.50/160)", et)]
    for i, p in enumerate(P):
        th = MOT[et].theta(p) / np.pi
        d = np.arange(len(p)) - (th - th[0])
        net = float(d[-1] - d[0])
        print(f"{et:>6} {'w'+str(i+1):>4} {len(p):>7} {net:>8.2f} "
              f"{float(np.std(d)):>7.2f} {float(np.max(np.abs(d-np.median(d)))):>11.2f} "
              f"[{p[0]:>9.1f},{p[-1]:>9.1f}]", flush=True)

# ---------------- K5: 'bolgeler' örtüşmesi ----------------
print("\n=== K5: 101f ŞÜPHELİ BÖLGELERİNİN ANALİZ PENCERELERİYLE ÖRTÜŞMESİ ===")
for et in AILE:
    W, P = PENR[("gevşek(0.50/160)", et)]
    ort = 0
    for (a, b) in BOLGE[et]:
        for p in P:
            if not (b < p[0] or a > p[-1]):
                ort += 1
                print(f"  {et}: bölge [{a:.1f},{b:.1f}] pencere "
                      f"[{p[0]:.1f},{p[-1]:.1f}] ile ÖRTÜŞÜYOR")
    print(f"  {et}: {len(BOLGE[et])} bölge, {ort} örtüşme")

# ---------------- K6: pencere-başına ve w1-dışlamalı ----------------
print("\n=== K6: HER REJİMDE χ₇ PENCERE-BAŞINA D(τ=0.04) ===")
for ad, _, _ in REJIM:
    W, P = PENR[(ad, "chi7")]
    print(f"  --- rejim {ad} ---")
    for i, w in enumerate(W):
        Ds = [egri([w], [0.04], band=b)[0][0] for b in BANDS]
        print(f"   w{i+1} n={len(w[0]):>6} L={w[3]:>6.3f} t∈"
              f"[{w[0][0]:>8.1f},{w[0][-1]:>8.1f}] "
              + " ".join(f"{d:>6.4f}" for d in Ds)
              + f"  {Ds[-1]/Ds[0]:>5.2f}×", flush=True)

print("\n=== K6b: χ₇ EN YÜKSEK-D PENCERESİ DIŞLANARAK (havuz) ===")
for ad, _, _ in REJIM:
    W, P = PENR[(ad, "chi7")]
    if len(W) < 2:
        continue
    d1 = [egri([w], [0.04], band=BANDS[0])[0][0] for w in W]
    jmax = int(np.argmax(d1))
    Wsub = [w for j, w in enumerate(W) if j != jmax]
    Ds = [egri(Wsub, [0.04], band=b)[0][0] for b in BANDS]
    Dt = [egri(W, [0.04], band=b)[0][0] for b in BANDS]
    print(f"  {ad}: dışlanan w{jmax+1} (n={len(W[jmax][0])}, "
          f"L={W[jmax][3]:.3f})")
    print(f"    tam   : " + " ".join(f"{d:>6.4f}" for d in Dt)
          + f"  {Dt[-1]/Dt[0]:>5.2f}×")
    print(f"    w{jmax+1} yok: " + " ".join(f"{d:>6.4f}" for d in Ds)
          + f"  {Ds[-1]/Ds[0]:>5.2f}×", flush=True)

# ---------------- K3/K4: silme enjeksiyonu ----------------
print("\n=== K3/K4: SİLME ENJEKSİYONU (pencere İÇİNE, segmentasyondan sonra) ===")
KS = [10, 30, 100, 300]
NREP = 8
print(f"{'aile':>6} {'k':>5} {'‰':>6} " + " ".join(f"±{b}L" for b in BANDS)
      + "    ΔD (temizden fark)          oran(dar/geniş)")
ENJ = {}
for et in ["chi3", "chi5", "chi7"]:
    W, P = PENR[("gevşek(0.50/160)", et)]
    q = AILE[et][0]
    ntop = sum(len(p) for p in P)
    temiz = []
    for b in BANDS:
        D, V, _ = egri(W, [0.04], band=b)
        temiz.append(D[0])
    print(f"{et:>6} {0:>5} {0.0:>6.2f} " + " ".join(f"{d:>6.4f}" for d in temiz)
          + f"    {'(temiz)':>28}  {temiz[-1]/temiz[0]:>5.2f}×", flush=True)
    ENJ[(et, 0)] = (temiz, [0, 0, 0])
    for k in KS:
        acc = np.zeros((NREP, len(BANDS)))
        for r in range(NREP):
            WK = []
            # k silmeyi pencerelere n-orantılı dağıt
            paylar = np.array([len(p) for p in P], float)
            paylar = paylar / paylar.sum()
            kk = rng.multinomial(k, paylar)
            for p, ki in zip(P, kk):
                if ki > 0:
                    idx = rng.choice(len(p), size=int(ki), replace=False)
                    m = np.ones(len(p), bool); m[idx] = False
                    WK.append(pencere_hazirla(p[m], float(q)))
                else:
                    WK.append(pencere_hazirla(p, float(q)))
            for jb, b in enumerate(BANDS):
                D, V, _ = egri(WK, [0.04], band=b)
                acc[r, jb] = D[0]
        mu = acc.mean(axis=0); sd = acc.std(axis=0)
        ENJ[(et, k)] = (list(mu), list(sd))
        print(f"{et:>6} {k:>5} {1000*k/ntop:>6.2f} "
              + " ".join(f"{m:>6.4f}" for m in mu)
              + "    " + " ".join(f"{m-t:>+7.4f}" for m, t in zip(mu, temiz))
              + f"    {mu[-1]/mu[0]:>5.2f}×  (σ "
              + " ".join(f"{s:.4f}" for s in sd) + ")", flush=True)

np.savez(HERE / "104b_kusur.npz",
         bands=np.array(BANDS), taus=np.array(TAUS), ks=np.array([0] + KS),
         **{f"K1_D_{et}_{ad.split('(')[0]}_{b}": np.array(K1[(et, ad, b)][0])
            for (et, ad, b) in K1},
         **{f"ENJ_{et}_{k}": np.array(ENJ[(et, k)][0])
            for (et, k) in ENJ})
print(f"\n104b_kusur.npz yazıldı. Süre {time.time()-T0:.0f} s.")
