"""
101h — DAR BANT + YENİDEN ÖLÇEKLEME: ÇÖZÜLME YAMACININ KESKİNLİĞİ (20 Ağu)
==========================================================================
101d ölçümü ω-bandı τ₀·L ± 0.02·L ile yapıldı. Band ne kadar genişse
eşik-öncesi τ'larda banda eşik-üstü frekanslar sızar → "yumuşak öncü"
(donuk bölgede sıfır olmayan artık D) kısmen band-kaçağı olabilir.
Test: bandı ikiye ve dörde daralt (±0.01L, ±0.005L), τ ≤ 0.16 ızgarası,
beş aile (ζ + 4 ada). Makine 101d ile BİREBİR aynı (pencere_hazirla /
egri / duz_segmentler kopyalandı; egri'ye yalnız band parametresi
eklendi — band=0.02 ile 101d sayıları birebir çıkmalı, kod-eşdeğerlik
kapısı).

VERİ DİSİPLİNİ: 101f sertifikalı sıfırlar + düzlük segmentasyonu (101b
kusur kapısının gereği), ζ için 41_bigT_windows.npz ilk 6 penceresi.
Her ölçümde 40-permütasyonlu vekil taban.

YENİDEN ÖLÇEKLEME: her ailenin ilk sağ kalan çizgisi τ_ilk =
log(q_ilk)/⟨L⟩ (ζ,χ₃,χ₅,χ₇: q_ilk=2; β: q_ilk=3, çünkü χ₄(2)=0 ölü).
⟨L⟩ = pencere L'lerinin sıfır-sayısı ağırlıklı ortalaması (estimatörün
havuzlaması da ~n ağırlıklı: den = Σ|r|² ~ n). x = τ/τ_ilk ekseninde
beş eğri tek eğriye çökmeli.

ÖN-MÜHÜR (ölçümden ÖNCE yazıldı):
  H1  Dar bantta eşik-öncesi "yumuşak öncü" KESKİNLEŞECEK: eşik-altı
      (x ≲ 0.9) D değerleri düşecek.
  H2  Çözülme yamacı DİKLEŞECEK (x≈1 civarı dD/dx artacak).
  H3  β'nın anahtar-biçimi (tabanda sürünme + dik sıçrama) KALACAK.
  H4  x = τ/τ_ilk çökmesi dar bantta İYİLEŞECEK mi? — açık soru,
      raporlanacak (çökme ölçütü: x-ızgarasında aileler-arası std).

==========================================================================
SONUÇ (20 Ağustos, koşu 95 s) — H1 ✓ H2 ✓✓ H3 ✓✓ H4 kısmen
==========================================================================
KOD-EŞDEĞERLİK KAPISI ✓: band=0.02L sütunu 101d ile BİREBİR
(τ=0.040: 0.041/0.069/0.010/0.030/0.078; τ=0.068: 0.341/0.394/0.040/
0.350/0.465; β 0.085→0.095→0.113: 0.116/0.460/0.705).

⟨L⟩_n (sıfır-ağırlıklı) ve τ_ilk:
  ζ 11.174 → 0.06203 | χ₃ 9.414 → 0.07363 | β 9.611 → 0.11431 (log3!)
  χ₅ 9.791 → 0.07079 | χ₇ 9.988 → 0.06940

D(τ), band ±0.02L → ±0.01L → ±0.005L (vekil taban her yerde ≤0.07):
  τ      ζ                χ₃               β                χ₅               χ₇
  0.040  .041 .021 .020   .069 .038 .018   .010 .012 .013   .030 .033 .034   .078 .106 .125
  0.050  .181 .134 .092   .130 .115 .092   .016 .018 .018   .112 .069 .065   .220 .166 .155
  0.060  .240 .451 .479   .314 .286 .327   .026 .026 .026   .242 .248 .266   .366 .430 .398
  0.068  .341 .489 .480   .394 .420 .392   .040 .037 .036   .350 .443 .501   .465 .585 .675
  0.075  .482 .393 .334   .533 .606 .718   .054 .053 .051   .429 .543 .547   .529 .618 .603
  0.085  .530 .527 .479   .594 .636 .651   .116 .091 .091   .517 .517 .558   .619 .488 .429
  0.095  .642 .649 .692   .605 .590 .554   .460 .178 .156   .594 .468 .293   .686 .638 .559
  0.105  .676 .753 .726   .545 .565 .560   .520 .647 .522   .653 .675 .760   .737 .802 .829
  0.113  .728 .721 .772   .489 .512 .479   .705 .726 .847   .711 .730 .698   .783 .832 .848
  0.125  .788 .747 .689   .583 .320 .232   .817 .872 .788   .779 .787 .778   .794 .761 .764
  0.140  .855 .909 .868   .727 .766 .837   .878 .954 .728   .781 .819 .818   .863 .818 .760
  0.160  .903 .926 .937   .817 .767 .841   .884 .948 1.057  .647 .645 .612   .903 .897 .843

H1 ✓ (İSABET, ama MEKANİZMASI BEKLENENDEN BAŞKA): eşik-öncesi D
  düşüyor — τ=0.050'de ζ .181→.092, χ₃ .130→.092, χ₅ .112→.065,
  χ₇ .220→.155 (β zaten tabanda: .016→.018). AMA yeni BAND KAÇAĞI
  KAPISI gösterdi ki bu düşüşün büyük kısmı ARTEFAKT TEMİZLİĞİ:
  pencere bandı [L(τ−b), L(τ+b)] ilk çizgiyi İÇERİYORSA (çizgi-dışı
  filtre yalnız ±0.01'i atıyor, Bragg omzu bantta kalıyor) eşik-altı
  ölçüme çizgi SIZAR. τ=0.050'de sızan pencere sayısı: ζ 5/6 → 2/6 →
  0/6; χ₅ 2/4 → 0/4; χ₇ 3/6 → 0/6. Yani 101d'nin "yumuşak öncüsü"
  KISMEN kendi bandının kaçağıydı. Kalan (kaçaksız) taban gerçek:
  τ=0.050, ±0.005L'de ζ .092, χ₃ .092, χ₅ .065, χ₇ .155 — vekilin
  (.009-.013) 7-15 katı. β .018 — pratikte taban.
  ANOMALİ (gizlenmiyor, AÇIKLANAMADI): χ₇ τ=0.040'ta band daralınca
  YÜKSELİYOR (.078→.106→.125) ve orada HİÇBİR pencere sızmıyor (0/6,
  üç bantta da). İlk akla gelen açıklama, χ₇ pencerelerinin ω=τ·L
  aralığının (0.356-0.427) (3,2) vuruşunu (ω=0.4055) kapsamasıydı —
  AMA 101i vuruş hipotezini duyarlılığı kanıtlanmış bir aletle
  REDDETTİ. Dolayısıyla bu anomali şu an açıklamasız duruyor;
  χ₇'nin en çok segmentli (6) ve en kısa parçalı aile olması bir
  segmentasyon/örnekleme etkisini akla getiriyor, ama test edilmedi.

H2 ✓✓ (İSABET, beş ailede de tekdüze): eşik geçiş dikliği
  Δ D/Δx, x: 0.80→1.00, band ±0.02L → ±0.01L → ±0.005L
    ζ  0.45 → 1.66 → 1.95     χ₃ 1.06 → 1.51 → 1.77
    β  1.89 → 2.97 → 3.54     χ₅ 0.91 → 1.48 → 1.60
    χ₇ 0.89 → 1.40 → 1.86
  Band dörde bölününce yamaç 1.8-4.3 kat dikleşiyor. Çözülme bir
  GEÇİŞ, geniş bantta yayvanlaşan bir sızıntı değil.

H3 ✓✓ (İSABET, üstelik KESKİNLEŞEREK): β dar bantta daha uzun
  sürünüp daha dik sıçrıyor. ±0.02L'de τ=0.095'te 0.460 görünüyordu;
  o noktada 1/2 pencere sızıyordu (log3 banda giriyor). ±0.005L'de
  hiç sızıntı yok ve β τ=0.095'te 0.156'da SÜRÜNMEYE DEVAM, sonra
  0.105'te 0.522, 0.113'te 0.847 (τ_ilk = log3/⟨L⟩ = 0.11431 —
  sıçrama tam çizginin üstünde). β'nın erken çözülüyor görünmesi
  bandın kendi kaçağıydı; temizlenince çizgiye kilit DAHA da net.

H4 KISMEN: çökme ölçütü (aileler-arası std/ort), band ±0.02→0.01→0.005
    x ≤ 0.80        0.2846 → 0.3704 → 0.4231   (mutlak std .052 .055 .061)
    0.85 ≤ x ≤ 1.10 0.2972 → 0.1992 → 0.1838   (mutlak std .130 .100 .094)
    x ≥ 1.15        0.1963 → 0.2536 → 0.2517   (mutlak std .122 .154 .139)
    tümü            0.2543 → 0.2552 → 0.2588   (mutlak std .111 .112 .105)
  EŞİK KOMŞULUĞUNDA ÇÖKME BELİRGİN İYİLEŞİYOR (%38 daha sıkı, mutlak
  std'de de). Donuk kuyrukta göreli std kötüleşiyor ama MUTLAK std
  neredeyse sabit — kötüleşme paydanın (ortalamanın) H1 gereği sıfıra
  çökmesinden; bilgi değil, normalizasyon eseri. x≥1.15'te dar bant
  gerçekten daha gürültülü (bağımsız frekans sayısı düşüyor).
  Küresel ölçüt düz: dar bant çökmeyi topyekûn İYİLEŞTİRMİYOR, ama
  ÖNEMLİ OLDUĞU YERDE — geçişin kendisinde — iyileştiriyor.

HÜKÜM: çözülme sınırı band daraltıldıkça KESKİNLEŞİYOR, yayvanlaşmıyor
— yani gerçek bir eşik. 101d'nin eşik-altı artığı iki parçaymış:
(i) bandın ilk çizgiden kaçağı (artefakt, dar bantta yok oluyor),
(ii) kalan gerçek taban (vekilin ~10 katı) — 101i bunun VURUŞ
moiré'si olup olmadığını soruyor.
"""

import numpy as np
from pathlib import Path

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
rng = np.random.default_rng(1010)

PKS = [2,3,4,5,7,8,9,11,13,16,17,19,23,25,27,29,31,32,37,41,43,47,49,53,
       59,61,64,67,71,73,79,81,83,89,97,101,103,107,109,113,121,125,127,128]
LINES_ALL = [np.log(qq) for qq in PKS]

# ---- 101d makinesi (birebir kopya; egri'ye band parametresi eklendi) ----

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
    """93-usulü: pencereler-arası havuzlanmış D(τ) + vekil taban."""
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
        Dv.append(abs(num) / (c * den)); Vv.append(np.mean(fl))
        Nf.append(nf_tot / max(nf_ham, 1))
    return Dv, Vv, Nf

def duz_segmentler(zz, sayim, min_n=3000, win=80, esik=0.5, pad=160):
    """Sayım-sürüklenmesi düz parçalara böl; basamak çevresini at."""
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

# ------------------------- pencereler -------------------------

def rvm_sayim(t):
    x = t / TWO_PI
    return x * np.log(x / np.e)

d41 = np.load(HERE / "41_bigT_windows.npz")
K41 = sorted({x.split("_")[1] for x in d41.files}, key=lambda s: int(s[:-1]))
WINZ, nk_z = [], 0
for k in K41[:6]:
    gz, tm = d41[f"gaps_{k}"], d41[f"tmid_{k}"]
    zz = np.empty(len(gz) + 1)
    zz[0] = tm[0] - gz[0] / 2
    zz[1:] = zz[0] + np.cumsum(gz)
    seg, kes = duz_segmentler(zz, rvm_sayim)
    nk_z += kes
    WINZ += [pencere_hazirla(s, 1.0) for s in seg]

exec(open(HERE / "98_L_motoru.py").read().split('CHI4 = ')[0])
HERE = Path(__file__).resolve().parent
CHI4 = {0: 0, 1: 1, 2: 0, 3: -1}
CHI3 = {0: 0, 1: 1, 2: -1}
CHI5 = {0: 0, 1: 1, 2: 1j, 3: -1j, 4: -1}
z6c = np.exp(1j * np.pi / 3)
CHI7 = {0: 0, 1: 1, 3: z6c, 2: z6c**2, 6: z6c**3, 4: z6c**4, 5: z6c**5}

PEN = {"zeta": WINZ}
for etiket, q, tab in [("chi3", 3, CHI3), ("beta", 4, CHI4),
                       ("chi5", 5, CHI5), ("chi7", 7, CHI7)]:
    zc = np.load(HERE / f"101f_{etiket}_zeros.npz")["zeros"]
    lo = np.exp(np.log(zc[0] + 1) + 0.25 * (np.log(zc[-1]) - np.log(zc[0] + 1)))
    zc = zc[zc >= lo]
    M = Lmotor(q, tab, 1)
    seg, kes = duz_segmentler(zc, lambda t: M.theta(t) / np.pi)
    kenar = np.exp(np.linspace(np.log(zc[0]), np.log(zc[-1] * 1.0001), 4))
    WIN = []
    for s in seg:
        for i in range(3):
            p = s[(s >= kenar[i]) & (s < kenar[i + 1])]
            if len(p) >= 3000:
                WIN.append(pencere_hazirla(p, float(q)))
    PEN[etiket] = WIN

ADLAR = ["zeta", "chi3", "beta", "chi5", "chi7"]
QILK = {"zeta": 2, "chi3": 2, "beta": 3, "chi5": 2, "chi7": 2}
LORT, TILK = {}, {}
for a in ADLAR:
    Ls = np.array([w[3] for w in PEN[a]])
    ns = np.array([len(w[0]) for w in PEN[a]], float)
    LORT[a] = float(np.average(Ls, weights=ns))
    TILK[a] = np.log(QILK[a]) / LORT[a]
    print(f"[{a}] {len(PEN[a])} segment  ⟨L⟩_n={LORT[a]:.3f} "
          f"(düz ort {Ls.mean():.3f})  τ_ilk={TILK[a]:.5f}  "
          f"L: " + " ".join(f"{x:.2f}" for x in Ls), flush=True)

TAUS = [0.04, 0.05, 0.06, 0.068, 0.075, 0.085, 0.095, 0.105, 0.113,
        0.125, 0.14, 0.16]
BANDS = [0.02, 0.01, 0.005]

SONUC = {}
for band in BANDS:
    for a in ADLAR:
        SONUC[(band, a)] = egri(PEN[a], TAUS, band=band)
        print(f"  ...band ±{band}L / {a} bitti", flush=True)

# ------------------------- tablolar -------------------------
for band in BANDS:
    print(f"\n=== BAND ±{band}L ===")
    print(f"{'τ':>6} " + " ".join(f"{a:>7}" for a in ADLAR)
          + "   vekil(maks)  frek-kalan")
    for i, t in enumerate(TAUS):
        vek = max(SONUC[(band, a)][1][i] for a in ADLAR)
        fk = min(SONUC[(band, a)][2][i] for a in ADLAR)
        print(f"{t:>6.3f} "
              + " ".join(f"{SONUC[(band, a)][0][i]:>7.3f}" for a in ADLAR)
              + f"   {vek:>7.3f}   {fk:>5.2f}")
    print("  vekil taban (aile aile):")
    for a in ADLAR:
        print(f"    {a:>5}: " + " ".join(f"{v:.3f}"
              for v in SONUC[(band, a)][1]))

# ------------------- yeniden ölçekleme (x = τ/τ_ilk) -------------------
XG = np.arange(0.70, 1.401, 0.05)
print("\n=== YENİDEN ÖLÇEKLEME: D(x), x = τ/τ_ilk ===")
COLL = {}
for band in BANDS:
    print(f"\n--- band ±{band}L ---")
    print(f"{'x':>6} " + " ".join(f"{a:>7}" for a in ADLAR) + "    std   std/ort")
    M = {}
    for a in ADLAR:
        xs = np.array(TAUS) / TILK[a]
        M[a] = np.interp(XG, xs, SONUC[(band, a)][0],
                         left=np.nan, right=np.nan)
    spreads, absstd = [], []
    for j, x in enumerate(XG):
        vals = np.array([M[a][j] for a in ADLAR])
        s = float(np.nanstd(vals)); m = float(np.nanmean(vals))
        spreads.append(s / m if m > 0 else np.nan)
        absstd.append(s)
        print(f"{x:>6.2f} " + " ".join(f"{v:>7.3f}" for v in vals)
              + f"  {s:>6.3f}  {s/m if m>0 else np.nan:>7.3f}")
    spreads = np.array(spreads); absstd = np.array(absstd)
    # x'in üç dilimi: donuk kuyruk / eşik komşuluğu / eşik üstü
    dilim = {"x≤0.80": XG <= 0.801,
             "0.85≤x≤1.10": (XG >= 0.849) & (XG <= 1.101),
             "x≥1.15": XG >= 1.149}
    COLL[band] = {k: (float(np.nanmean(spreads[m_])),
                      float(np.nanmean(absstd[m_]))) for k, m_ in dilim.items()}
    COLL[band]["tümü"] = (float(np.nanmean(spreads)), float(np.nanmean(absstd)))
    for k, (r, s) in COLL[band].items():
        print(f"  ÇÖKME {k:>12}: std/ort {r:.4f}   mutlak std {s:.4f}")

print("\n=== ÇÖKME KARŞILAŞTIRMASI (std/ort | mutlak std) ===")
print(f"{'dilim':>12} " + " ".join(f"{'±'+str(b)+'L':>18}" for b in BANDS))
for k in ["x≤0.80", "0.85≤x≤1.10", "x≥1.15", "tümü"]:
    print(f"{k:>12} " + " ".join(f"{COLL[b][k][0]:>8.4f} |{COLL[b][k][1]:>8.4f}"
                                 for b in BANDS))

print("\n=== EŞİK-ÖNCESİ D ve EŞİK GEÇİŞ DİKLİĞİ (x: 0.80 → 1.00) ===")
print(f"{'aile':>6} " + "".join(f"   ±{b}L: D(0.8)  D(1.0)  Δ/Δx" for b in BANDS))
for a in ADLAR:
    satir = f"{a:>6} "
    for band in BANDS:
        xs = np.array(TAUS) / TILK[a]
        d08 = np.interp(0.80, xs, SONUC[(band, a)][0])
        d10 = np.interp(1.00, xs, SONUC[(band, a)][0])
        satir += f"     {d08:>6.3f}  {d10:>6.3f} {(d10-d08)/0.20:>6.2f}"
    print(satir)

print("\n=== HAM τ IZGARASINDA EŞİK-ÖNCESİ (τ=0.040, 0.050) ===")
print(f"{'aile':>6}  " + "  ".join(f"τ={t}: " + " ".join(f"±{b}L" for b in BANDS)
                                   for t in [0.040, 0.050]))
for a in ADLAR:
    satir = f"{a:>6}  "
    for i in [0, 1]:
        satir += "  " + " ".join(f"{SONUC[(b, a)][0][i]:>6.3f}" for b in BANDS)
    print(satir)

print("\n=== BAND KAÇAĞI KAPISI (kritik artefakt denetimi) ===")
print("  Pencere w'nin bandı [L_w(τ−band), L_w(τ+band)]. İlk çizgi log q_ilk")
print("  banda GİRERSE eşik-altı ölçüme çizgi sızar → sahte 'yumuşak öncü'.")
print("  Sızma eşiği: τ_sız(w) = log q_ilk / L_w − band.")
print(f"{'aile':>6} {'logq':>7} " +
      "".join(f"  ±{b}L: τ_sız(ilk) kaçpencere@τ_ilk" for b in BANDS))
for a in ADLAR:
    Ls = np.array([w[3] for w in PEN[a]])
    ln = np.log(QILK[a])
    satir = f"{a:>6} {ln:>7.4f} "
    for b in BANDS:
        tsz = ln / Ls - b          # her pencere için sızma eşiği
        n_sz = int(np.sum(tsz <= TILK[a]))   # τ_ilk'te sızan pencere sayısı
        satir += f"      {tsz.min():>7.4f}        {n_sz}/{len(Ls)}    "
    print(satir)
print("\n  τ ızgarasında SIZAN PENCERE SAYISI (n_sız / n_pencere):")
for b in BANDS:
    print(f"  --- band ±{b}L ---")
    print(f"{'aile':>6} " + " ".join(f"{t:>6.3f}" for t in TAUS))
    for a in ADLAR:
        Ls = np.array([w[3] for w in PEN[a]])
        ln = np.log(QILK[a])
        print(f"{a:>6} " + " ".join(
            f"{int(np.sum(Ls*(t+b) >= ln)):>6d}" for t in TAUS))

np.savez(HERE / "101h_darbant.npz",
         taus=np.array(TAUS), bands=np.array(BANDS),
         adlar=np.array(ADLAR), Lort=np.array([LORT[a] for a in ADLAR]),
         **{f"D_{b}_{a}": np.array(SONUC[(b, a)][0])
            for b in BANDS for a in ADLAR},
         **{f"V_{b}_{a}": np.array(SONUC[(b, a)][1])
            for b in BANDS for a in ADLAR})
print("\n101h_darbant.npz yazıldı.")
