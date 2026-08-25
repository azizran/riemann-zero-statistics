"""
107c — SERTİFİKA DEFTERİ: SEKİZ VERİ KÜMESİNİN NİHAİ DURUMU (25 Ağu 2026)
==========================================================================
Not 5'in reproducibility bölümüne hazır malzeme. Sekiz veri kümesi
(ζ + 101f'nin dört eski adası + 105b'nin üç yeni adası) için TEK TABLO:

  n_ham        kaynak dosyadaki sıfır sayısı
  n_analiz     lo-kesiminden (alt %25 log-aralık atılır) sonra kalan
  kalan bölge  101f/105b sertifikasının kapatamadığı şüpheli bölge sayısı
  K1/K2/K3     üç kusur kapısının ham ihlal sayısı (pad'siz işaret)
  kopya        normalize boşluk u < 0.005 olan çiftler
  maskelenen   analiz penceresine GİRMEYEN sıfır sayısı (pad + n<3000
               düşen parçalar dâhil)
  n_pencere    nihai analiz penceresi sayısı, ⟨L⟩_n, τ_ilk

Ayrıca ADLI KUSUR DENETİMİ: 104d/105e'nin elle bulduğu kusur
koordinatlarının nihai pencerelerin DIŞINDA kaldığı tek tek doğrulanır.

Ölçüm yok (egri çağrılmaz) — bu yalnız veri disiplininin defteridir.
==========================================================================
SONUÇ (25 Ağustos, koşu 9 s) — SEKİZ KÜMENİN NİHAİ SERTİFİKASI
==========================================================================
 küme   n_ham  n_analiz kalan   K1   K2   K3  kopya kopya maskelenen  %   n_pen  ⟨L⟩   τ_ilk
                (lo-kes) bölge  ham  ham  ham   ada  pen.    sıfır
 zeta  239990   239990     0     0    0    0     0     0         0   0.0    6  11.17 0.0620
 chi3   80211    79672     2    85   24    0     0     0      4117   5.2    3   9.43 0.0735
 chi5   72836    72279     3   254   23    0     0     0      3966   5.5    4   9.79 0.0708
 chi7   65055    64494     8   789  121   22     6     0      6447  10.0    6  10.04 0.0691
 chi5e  72829    72275     2   169   20    0     0     0      2685   3.7    2   9.79 0.0708
 beta   74430    73887     1    84    0    0     0     0      3006   4.1    2   9.61 0.1143
 chi8e  65948    65377     2    85   12    0     1     0      3518   5.4    4  10.13 0.1085
 chi8o  65943    65372     0     0    0    0     0     0      2629   4.0    2  10.13 0.1085

OKUMA:
 • ζ ve chi8o HİÇBİR kapıyı tetiklemiyor (ζ 239990 sıfırın tamamı analize
   giriyor: maskelenen 0). β yalnız K1'i tetikliyor, K2/K3 sıfır —
   105e'nin "β'da maks fark 0.37 → 0 kusur" yanlış-pozitif denetimi
   bağımsız olarak tekrarlandı (bu defterde 0.37 aynen çıktı).
 • KUSUR SIRALAMASI: χ₇ ≫ χ₅ ≈ χ₃ ≈ chi5e ≈ chi8e > β ≈ chi8o ≈ ζ.
   χ₇ tek başına K1'in %54'ünü, K2'nin %61'ini ve K3'ün TAMAMINI taşıyor;
   analizden düşen sıfırların oranı da onda en yüksek (%10.0).
 • KOPYA (u<0.005) ENVANTERİ: χ₇'de 6, chi8e'de 1, kalan altı kümede 0.
   YEDİSİ DE nihai analiz pencerelerinin DIŞINDA ("kopya pen." sütunu
   sekiz kümede de 0). Bu, KAPI 3'ün u-ölçütünün neden hiç
   tetiklenmediğini açıklar: K1/K2 ve n<3000 elemesi onları zaten
   dışarıda bırakıyor. Kopya ölçütü GEREKSİZ değil — SON EMNİYET
   KİLİDİ olarak duruyor.
 • KALAN BÖLGE = 101f/105b sertifikasının kapatamadığı şüpheli bölge:
   toplam 18 (χ₇ 8, χ₅ 3, χ₃ 2, chi5e 2, chi8e 2, β 1, chi8o 0, ζ 0).
   Bu sayı bir KUSUR ölçüsü DEĞİLDİR: 104 χ₇ için 101f'nin 8 şüpheli
   bölgesiyle analiz pencerelerinin örtüşmesini 0 bulmuştu (suçlu
   listede değildi). Defterde bilgi olarak durur; hüküm K1/K2/K3'ündür.

ADLI KUSUR DENETİMİ — HEPSİ ANALİZ DIŞINDA (8/8 ✓):
  χ₇ 3081.6688 / 3081.6695 / 3081.8885 / 3081.8891 (104d kopya çiftleri)
  χ₇ ~3095 (104d |d|>1 bloğunun ortası)
  χ₃ 3872.54 ve 3890.24 (104d χ₃ w1 kusuru)
  chi8e 19352 (105e kısa çukuru)
  Yani iki YENİ kapı, önceki iki seferde ELLE bulunan kusurların
  tamamını OTOMATİK olarak analiz dışında tutuyor.

NİHAİ ANALİZ PENCERELERİ (Not 5 için birebir alıntılanabilir):
  zeta  [107252,132748] [187880,212120] [338501,361499] [589042,610958]
        [989509,1010491] [1589905,1610095]
  chi3  [4026,8666] [9090,13553] [13554,55000]
  chi5  [3114,5692] [6222,12226] [12227,43330] [43668,48000]
  chi7  [4132,11063] [11064,14732] [15100,17945] [18304,33637]
        [33974,36119] [36454,42000]
  chi5e [3109,12215] [12215,47999]
  beta  [3445,12596] [12597,50000]
  chi8e [2907,11049] [11050,19252] [19459,27319] [27657,42000]
  chi8o [2913,11062] [11063,42000]
  (χ₇'nin ilk penceresi 2914 yerine 4132'de başlıyor: K3'ün t≈3672
   sıçraması + pad 600 orayı kesiyor — 104d'nin "i₀ ≥ 1000, t ≥ 3680"
   elle bulduğu sınır kapı tarafından bağımsız olarak yeniden bulundu.)
==========================================================================
"""
import numpy as np, time
from pathlib import Path

T0 = time.time()
HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi


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


def segmentle(zz, bad, min_n=3000):
    seg, kes, s0 = [], 0, 0
    for i in range(1, len(zz) + 1):
        if i == len(zz) or bad[i] != bad[i - 1]:
            if not bad[s0] and i - s0 >= min_n:
                seg.append((s0, i))
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
CHI5E = {0: 0, 1: 1, 2: -1, 3: -1, 4: 1}
CHI8E = {0: 0, 1: 1, 2: 0, 3: -1, 4: 0, 5: -1, 6: 0, 7: 1}
CHI8O = {0: 0, 1: 1, 2: 0, 3: 1, 4: 0, 5: -1, 6: 0, 7: -1}

# (etiket, dosya, q, tablo, a, q_ilk, log-bölme?)
KUME = [("zeta",  None,                   1, None,  None, 2, False),
        ("chi3",  "101f_chi3_zeros.npz",  3, CHI3,  1, 2, True),
        ("chi5",  "101f_chi5_zeros.npz",  5, CHI5,  1, 2, True),
        ("chi7",  "101f_chi7_zeros.npz",  7, CHI7,  1, 2, True),
        ("chi5e", "105b_chi5e_zeros.npz", 5, CHI5E, 0, 2, True),
        ("beta",  "101f_beta_zeros.npz",  4, CHI4,  1, 3, True),
        ("chi8e", "105b_chi8e_zeros.npz", 8, CHI8E, 0, 3, True),
        ("chi8o", "105b_chi8o_zeros.npz", 8, CHI8O, 1, 3, True)]

# 104d / 105e'nin ADLI kusurları — nihai pencerelerin dışında kalmalı
ADLI = [("chi7", 3081.6688, "104d kopya çifti #1"),
        ("chi7", 3081.6695, "104d kopya çifti #1"),
        ("chi7", 3081.8885, "104d kopya çifti #2"),
        ("chi7", 3081.8891, "104d kopya çifti #2"),
        ("chi7", 3095.0,    "104d |d|>1 bloğunun ortası"),
        ("chi3", 3872.54,   "104d χ₃ w1 kusuru (baş)"),
        ("chi3", 3890.24,   "104d χ₃ w1 kusuru (son)"),
        ("chi8e", 19352.0,  "105e kısa çukur (2 kayıp + telafi)")]


def rvm_sayim(t):
    x = t / TWO_PI
    return x * np.log(x / np.e)


def kapilari_kos(zc, sayim, q, logbol):
    n = len(zc)
    d = np.arange(n) - (sayim(zc) - sayim(zc[0]))
    from numpy.lib.stride_tricks import sliding_window_view
    med = np.median(sliding_window_view(d, 80), axis=1)
    adim = med[81:] - med[:-81]
    ham1 = np.zeros(n, bool)
    ham1[np.clip(np.where(np.abs(adim) > 0.5)[0] + 40, 0, n - 1)] = True
    g1 = np.zeros(n, bool)
    for i in np.where(ham1)[0]:
        g1[max(0, i - 160):i + 320] = True
    sap = np.abs(_roll(d, 20, np.median) - _roll(d, 200, np.median))
    ham2 = sap > 0.7
    g2 = np.zeros(n, bool)
    for i in np.where(ham2)[0]:
        g2[max(0, i - 160):i + 160] = True
    seg, kes12 = segmentle(zc, g1 | g2)
    parcalar = []
    if logbol:
        kenar = np.exp(np.linspace(np.log(zc[0]), np.log(zc[-1] * 1.0001), 4))
        for (a, b) in seg:
            s = zc[a:b]
            for i in range(3):
                m = (s >= kenar[i]) & (s < kenar[i + 1])
                if m.sum() >= 3000:
                    j = np.where(m)[0]
                    parcalar.append((a + j[0], a + j[-1] + 1))
    else:
        parcalar = seg
    # KAPI 3 (pencere içi)
    son, ham3, dup3, rmax = [], 0, [], []
    for (a, b) in parcalar:
        p = zc[a:b]
        dd = np.arange(len(p)) - (sayim(p) - sayim(p[0]))
        r = _roll(dd, 21, np.mean) - _roll(dd, 801, np.median)
        gp = np.diff(p)
        mid = 0.5 * (p[:-1] + p[1:])
        u = gp * np.log(q * mid / TWO_PI) / TWO_PI
        bad = np.abs(r) > 0.30
        for i in np.where(u < 0.005)[0]:
            bad[i] = True
            bad[min(i + 1, len(bad) - 1)] = True
            dup3.append(float(p[i]))
        ham3 += int(bad.sum())
        rmax.append(float(np.abs(r).max()))
        kes = np.zeros(len(p), bool)
        for i in np.where(bad)[0]:
            kes[max(0, i - 600):min(len(p), i + 601)] = True
        s0 = 0
        for i in range(1, len(p) + 1):
            if i == len(p) or kes[i] != kes[i - 1]:
                if not kes[s0] and i - s0 >= 3000:
                    son.append((a + s0, a + i))
                s0 = i
    # tüm ada boyunca kopya sayımı (pencere dışı dâhil, kayda geçsin)
    gp = np.diff(zc)
    mid = 0.5 * (zc[:-1] + zc[1:])
    u_all = gp * np.log(q * mid / TWO_PI) / TWO_PI
    kop_all = [float(zc[i]) for i in np.where(u_all < 0.005)[0]]
    return dict(n=n, ham1=int(ham1.sum()), ham2=int(ham2.sum()), ham3=ham3,
                g12=int((g1 | g2).sum()), kes12=kes12,
                sapmax=float(sap.max()), rmax=rmax,
                dup_pencere=dup3, dup_ada=kop_all, son=son,
                nson=sum(b - a for (a, b) in son), npen=len(son),
                ara=parcalar)


print("=" * 96)
print("107c — SERTİFİKA DEFTERİ: SEKİZ VERİ KÜMESİ")
print("=" * 96, flush=True)

DEF = {}
PENCERE_T = {}
for et, dosya, q, tab, a, qilk, logbol in KUME:
    if et == "zeta":
        d41 = np.load(HERE / "41_bigT_windows.npz")
        K41 = sorted({x.split("_")[1] for x in d41.files},
                     key=lambda s: int(s[:-1]))
        tot = dict(n=0, ham1=0, ham2=0, ham3=0, g12=0, kes12=0, sapmax=0.0,
                   rmax=[], dup_pencere=[], dup_ada=[], nson=0, npen=0)
        PENCERE_T[et] = []
        Ls, Ns = [], []
        for k in K41[:6]:
            gz, tm = d41[f"gaps_{k}"], d41[f"tmid_{k}"]
            zz = np.empty(len(gz) + 1)
            zz[0] = tm[0] - gz[0] / 2
            zz[1:] = zz[0] + np.cumsum(gz)
            r = kapilari_kos(zz, rvm_sayim, 1.0, False)
            for kk in ["n", "ham1", "ham2", "ham3", "g12", "kes12", "nson",
                       "npen"]:
                tot[kk] += r[kk]
            tot["sapmax"] = max(tot["sapmax"], r["sapmax"])
            tot["rmax"] += r["rmax"]
            tot["dup_pencere"] += r["dup_pencere"]
            tot["dup_ada"] += r["dup_ada"]
            for (x, y) in r["son"]:
                p = zz[x:y]
                PENCERE_T[et].append((float(p[0]), float(p[-1])))
                mids = 0.5 * (p[:-1] + p[1:])
                Ls.append(float(np.log(1.0 * mids / TWO_PI).mean()))
                Ns.append(len(p))
        tot["nham"] = tot["n"]
        tot["bolge"] = 0
        tot["Lort"] = float(np.average(Ls, weights=np.array(Ns, float)))
        DEF[et] = tot
        DEF[et]["qilk"] = qilk
        continue

    dd = np.load(HERE / dosya)
    zraw = dd["zeros"]
    bolge = int(np.asarray(dd["bolgeler"]).reshape(-1, 2).shape[0])
    lo = np.exp(np.log(zraw[0] + 1) + 0.25 *
                (np.log(zraw[-1]) - np.log(zraw[0] + 1)))
    zc = zraw[zraw >= lo]
    M = Lmotor(q, tab, a)
    say = (lambda t, M=M: M.theta(t) / np.pi)
    r = kapilari_kos(zc, say, float(q), logbol)
    r["nham"] = len(zraw)
    r["bolge"] = bolge
    r["qilk"] = qilk
    Ls, Ns = [], []
    PENCERE_T[et] = []
    for (x, y) in r["son"]:
        p = zc[x:y]
        PENCERE_T[et].append((float(p[0]), float(p[-1])))
        mids = 0.5 * (p[:-1] + p[1:])
        Ls.append(float(np.log(q * mids / TWO_PI).mean()))
        Ns.append(len(p))
    r["Lort"] = float(np.average(Ls, weights=np.array(Ns, float)))
    r["zc"] = zc
    DEF[et] = r
    print(f"  [{et}] işlendi ({time.time()-T0:.0f} s)", flush=True)

SIRA = ["zeta", "chi3", "chi5", "chi7", "chi5e", "beta", "chi8e", "chi8o"]

print("\n" + "=" * 118)
print("SERTİFİKA TABLOSU — SEKİZ VERİ KÜMESİ (Not 5 §reproducibility)")
print("=" * 118)
print(f"{'küme':>6} {'n_ham':>7} {'n_analiz':>8} {'kalan':>5} "
      f"{'K1':>5} {'K2':>5} {'K3':>5} {'kopya':>6} {'kopya':>6} "
      f"{'maskelenen':>10} {'%':>5} {'n_pen':>5} {'⟨L⟩':>6} {'τ_ilk':>7}")
print(f"{'':>6} {'':>7} {'(lo-kes)':>8} {'bölge':>5} "
      f"{'ham':>5} {'ham':>5} {'ham':>5} {'ada':>6} {'pen.':>6} "
      f"{'sıfır':>10} {'':>5} {'':>5} {'':>6} {'':>7}")
for a in SIRA:
    r = DEF[a]
    mask = r["n"] - r["nson"]
    tilk = np.log(r["qilk"]) / r["Lort"]
    print(f"{a:>6} {r['nham']:>7} {r['n']:>8} {r['bolge']:>5} "
          f"{r['ham1']:>5} {r['ham2']:>5} {r['ham3']:>5} "
          f"{len(r['dup_ada']):>6} {len(r['dup_pencere']):>6} "
          f"{mask:>10} {100*mask/r['n']:>5.1f} {r['npen']:>5} "
          f"{r['Lort']:>6.2f} {tilk:>7.4f}")

print("\n  K1 = düzlük/basamak (101d, |Δmed₈₀|>0.5)")
print("  K2 = kısa-çukur (105e-K4, |med₂₀−med₂₀₀|>0.7)")
print("  K3 = sıçrama (104e, |ort₂₁−med₈₀₁|>0.30), pencere-içi ölçülür")
print("  'kopya ada' = tüm adada u<0.005 çift sayısı (analiz dışı dâhil)")
print("  'kopya pen.' = nihai analiz penceresine kadar sağ kalan kopya")

print(f"\n{'küme':>6} {'maks|med₂₀−med₂₀₀|':>19} {'K3 maks|r| pencere başına':>34}")
for a in SIRA:
    r = DEF[a]
    print(f"{a:>6} {r['sapmax']:>19.2f}   "
          + " ".join(f"{x:.2f}" for x in r["rmax"]))

print("\n" + "=" * 96)
print("ADLI KUSUR DENETİMİ — 104d/105e'nin kusurları analiz DIŞINDA mı?")
print("=" * 96)
tum_temiz = True
for (et, t, ad) in ADLI:
    ic = any(x <= t <= y for (x, y) in PENCERE_T[et])
    if ic:
        tum_temiz = False
    print(f"  [{et}] t = {t:>12.4f}  {ad:<34} "
          + ("!!! PENCERE İÇİNDE — KAPI KAÇIRDI" if ic
             else "✓ analiz dışı"))
print(f"\n  HÜKÜM: " + ("TÜM ADLI KUSURLAR ANALİZ DIŞINDA"
                        if tum_temiz else "EN AZ BİR KUSUR İÇERİDE"))

print("\n" + "=" * 96)
print("NİHAİ ANALİZ PENCERELERİ (t-aralıkları)")
print("=" * 96)
for a in SIRA:
    print(f"  [{a}] {DEF[a]['npen']} pencere: " + "  ".join(
        f"[{x:.0f},{y:.0f}]" for (x, y) in PENCERE_T[a]))

print("\n" + "=" * 96)
print("TÜM ADALARDA KOPYA SIFIR ENVANTERİ (u < 0.005)")
print("=" * 96)
for a in SIRA:
    ka = DEF[a]["dup_ada"]
    if not ka:
        print(f"  [{a}] kopya YOK")
    else:
        print(f"  [{a}] {len(ka)} kopya: " + " ".join(f"{t:.4f}" for t in ka))
        for t in ka:
            ic = any(x <= t <= y for (x, y) in PENCERE_T[a])
            print(f"        t={t:.4f}  " + ("ANALİZ İÇİNDE (!)" if ic
                                            else "analiz dışı ✓"))

np.savez(HERE / "107c_sertifika.npz", adlar=np.array(SIRA),
         **{f"{k}_{a}": DEF[a][k] for a in SIRA
            for k in ["nham", "n", "bolge", "ham1", "ham2", "ham3", "g12",
                      "kes12", "nson", "npen", "sapmax", "Lort", "qilk"]},
         **{f"dupada_{a}": np.array(DEF[a]["dup_ada"]) for a in SIRA},
         **{f"pen_{a}": np.array(PENCERE_T[a]).reshape(-1, 2) for a in SIRA})
print(f"\n107c_sertifika.npz yazıldı. Süre {time.time()-T0:.0f} s.", flush=True)
