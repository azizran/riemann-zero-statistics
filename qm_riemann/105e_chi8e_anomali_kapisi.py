"""
105e — chi8e EŞİK-ALTI ANOMALİSİNİN KAPISI (25 Ağustos 2026)
==========================================================================
105d'nin PARİTE KONTROLÜ eşik ÜSTÜNDE kusursuz (τ≥0.095'te |Δ|≤0.10,
yarı-çözülme 0.0944 vs 0.0943) AMA eşik ALTINDA ayrışıyor:
  τ      0.040 0.050 0.060 0.068 0.075 0.085
  chi8e  0.122 0.192 0.227 0.218 0.168 0.127
  chi8o  0.025 0.031 0.039 0.047 0.063 0.110
İki ada AYNI iletken (q=8), AYNI ölü aile (2'nin kuvvetleri), AYNI
çizgi envanteri — yalnız parite farklı. P1/P3 paritenin hiçbir yasaya
girmediğini söylüyor; öyleyse bu fark ya (i) gerçek bir parite etkisi
(P1'i kısmen çürütür), ya (ii) SEGMENTASYON/PENCERE artefaktı.

ŞÜPHE: 105d'de chi8e 3 segment-pencere (L = 9.04, 10.08, 10.69; 1
kesim) chi8o ise 2 (L = 9.04, 10.38; 0 kesim) kullanıyor — düzlük
segmentasyonu iki adayı FARKLI böldü. 101h'nin açıklanamayan
"χ₇ τ=0.04 anomalisi" de en çok segmentli adadaydı (6 segment) —
aynı sınıf şüphe.

KAPILAR:
  K1  PENCERE AYRIŞTIRMASI: her adanın her penceresi için ayrı D(τ).
      Anomali TEK pencereye lokalize ise → aletsel.
  K2  EŞLEŞTİRİLMİŞ PENCERE: iki adayı AYNI t-aralığında, AYNI log-
      pencere kenarlarıyla, segmentasyonu ORTAK maskeyle (iki adanın
      kesimlerinin BİRLEŞİMİ) yeniden ölç. Fark kalırsa gerçek.
  K3  BOŞLUK İSTATİSTİĞİ: sahte/kayıp sıfır imzası (s≈0 fazlalığı).
  K4  SEGMENTASYON DÜZELTMESİ (K1-K3'ten sonra eklendi).
HÜKÜM: |Δ| eşik-altı τ'larda ≲0.05'e inerse ANOMALİ ALETSELDİR
(P1 ayakta); inmezse P1 eşik-altı dinamikte KISMEN ihlal edilmiştir
ve açıkça öyle raporlanır.

SONUÇ: K4 ile eşik-altı maks |Δ| 0.188 → 0.045 (vekil taban 0.005-0.018)
— ANOMALİ ALETSEL, P1 AYAKTA. chi8e artık eşik altında chi8o'dan bile
temiz (0.002-0.065 vs 0.025-0.110); eşik üstü zaten aynıydı.

==========================================================================
K1-K3 SONUÇLARI ve K4'ÜN GEREKÇESİ (ölçüldükten sonra yazıldı)
==========================================================================
K1 ✓ LOKALİZE: anomali chi8e'nin TEK penceresinde (L=10.08, t∈[11062,
   27368]) yaşıyor: D = 0.305/0.425/0.485/0.468/0.403/0.267 —
   diğer iki penceresi (L=9.04 ve 10.69) 0.02-0.17 ile chi8o gibi.
   ⟹ ARİTMETİK OLAMAZ: karakter tüm t'de aynıdır; t-yerel bir fark
   ancak VERİ kusurudur.
K2 △ eşleştirilmiş pencere/maske farkı KAPATMADI (0.188 → 0.151):
   demek ki kusur, sayım sertifikasının BASAMAK dedektörünün
   göremediği bir yerde.
K3 ✓ boşluk istatistiği TEMİZ (chi8e min normalize boşluk 0.0045,
   s<0.05 sayısı 24 vs chi8o 22) — kaba kopya/kayıp imzası YOK.

TEŞHİS (iki adımda; ilk adım YANILTICIYDI, kayda geçiyor):
 (1) Kaba sürüklenme haritası chi8e p3'te (t 11062-42000) std 1.038 ve
     uç-uca +2.014 gösterdi → "yavaş birikim" sanıldı. Bu okumaya göre
     yazılan KÜMÜLATİF SÜRÜKLENME kapısı BAŞARISIZ oldu (|Δ| 0.188 →
     0.214); çünkü o +2.014 zaten bilinen bir basamaktan (t≈27368,
     sertifikanın kalan 2 bölgesinden biri) geliyordu.
 (2) BLOK AYRIŞTIRMASI kesin sonucu verdi: anomalinin TAMAMI
     t∈[15000,19500) bloğunda (D 0.295-0.486; komşu bloklar 0.008-0.086)
     ve orada ham d, t≈19352'de 2.5 basamak DÜŞÜP ~10 sıfır sonra
     KAPANIYOR — 2 sıfır kaybı + hemen telafi, NET BASAMAK YOK.
     Ölçek: yuvarlanan medyan aralığı win=20 → 1.60, win=40 → 0.25,
     win=80 → 0.16 (temiz taban 0.30/0.22/0.17). 101'in win=80
     dedektörü bu kusuru YIKIYOR; 101b ise donma estimatörünün tam da
     buna aşırı duyarlı olduğunu göstermişti (%0.1 → 0.09→0.55).

K4 — DÜZELTME (ÇALIŞAN): KISA-ÖLÇEK KUSUR KAPISI — kısa (20) ve uzun
(200) yuvarlanan medyanların farkı > 0.7 olan sıfırlar kusurludur,
±160 komşuluğuyla atılır. Yanlış-pozitif denetimi: β'da maks fark 0.37
→ 0 kusur; beş adada ada başına yalnız 20-38 sıfır işaretleniyor.
Bu kapı 105d'ye standart olarak taşındı.
"""

import numpy as np
from pathlib import Path

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
rng = np.random.default_rng(1010)

PKS = [2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 17, 19, 23, 25, 27, 29, 31, 32, 37,
       41, 43, 47, 49, 53, 59, 61, 64, 67, 71, 73, 79, 81, 83, 89, 97, 101,
       103, 107, 109, 113, 121, 125, 127, 128]
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


def egri(WIN, taus, band=0.01):
    Dv, Vv = [], []
    for tau0 in taus:
        num = 0.0 + 0j; den = 0.0
        Gs, rs = [], []
        for (zz, tm, ds, Lw) in WIN:
            oms = tau0 * Lw + np.linspace(-band * Lw, band * Lw, 160)
            oms = np.array([o for o in oms
                            if min(abs(o - l) for l in LINES_ALL) > 0.01])
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
        fl = [abs((G_all * np.conj(r_all[rng.permutation(len(r_all))])).sum())
              / (c * den) for _ in range(40)]
        Dv.append(abs(num) / (c * den)); Vv.append(np.mean(fl))
    return Dv, Vv


def bad_maskesi(zz, sayim, win=80, esik=0.5, pad=160):
    d = np.arange(len(zz)) - (sayim(zz) - sayim(zz[0]))
    from numpy.lib.stride_tricks import sliding_window_view
    med = np.median(sliding_window_view(d, win), axis=1)
    adim = med[win + 1:] - med[:-(win + 1)]
    bad = np.zeros(len(zz), bool)
    for i in np.where(np.abs(adim) > esik)[0] + win // 2:
        bad[max(0, i - pad):i + 2 * pad] = True
    return bad


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


exec(open(HERE / "98_L_motoru.py").read().split('CHI4 = ')[0])
HERE = Path(__file__).resolve().parent
CHI8E = {0: 0, 1: 1, 2: 0, 3: -1, 4: 0, 5: -1, 6: 0, 7: 1}
CHI8O = {0: 0, 1: 1, 2: 0, 3: 1, 4: 0, 5: -1, 6: 0, 7: -1}
TAUS = [0.040, 0.050, 0.060, 0.068, 0.075, 0.085, 0.095, 0.105, 0.113]

ADA = [("chi8e", CHI8E, 0), ("chi8o", CHI8O, 1)]
Z, BAD, M_ = {}, {}, {}
for et, tab, a in ADA:
    zc = np.load(HERE / f"105b_{et}_zeros.npz")["zeros"]
    lo = np.exp(np.log(zc[0] + 1) + 0.25 * (np.log(zc[-1]) - np.log(zc[0] + 1)))
    zc = zc[zc >= lo]
    M = Lmotor(8, tab, a)
    Z[et] = zc; M_[et] = M
    BAD[et] = bad_maskesi(zc, lambda t: M.theta(t) / np.pi)
    print(f"[{et}] n={len(zc)}, t∈[{zc[0]:.0f},{zc[-1]:.0f}], "
          f"kusurlu sıfır {int(BAD[et].sum())}", flush=True)

# ---------------- K1: pencere ayrıştırması ----------------
print("\n=== K1 — PENCERE AYRIŞTIRMASI (105d'nin kendi bölmesi) ===",
      flush=True)
for et, tab, a in ADA:
    zc = Z[et]
    seg, kes = segmentle(zc, BAD[et])
    kenar = np.exp(np.linspace(np.log(zc[0]), np.log(zc[-1] * 1.0001), 4))
    WIN = []
    for s in seg:
        for i in range(3):
            pw = s[(s >= kenar[i]) & (s < kenar[i + 1])]
            if len(pw) >= 3000:
                WIN.append(pencere_hazirla(pw, 8.0))
    print(f"  {et}: {len(WIN)} pencere ({kes} kesim)", flush=True)
    for j, w in enumerate(WIN):
        D, V = egri([w], TAUS)
        print(f"    pencere {j+1} (L={w[3]:.2f}, n={len(w[0])}): "
              + " ".join(f"{d:.3f}" for d in D), flush=True)
        print(f"      vekil:                        "
              + " ".join(f"{v:.3f}" for v in V), flush=True)

# ---------------- K2: eşleştirilmiş pencere ----------------
print("\n=== K2 — EŞLEŞTİRİLMİŞ PENCERE (ortak maske + ortak kenarlar) ===",
      flush=True)
# ortak t-aralığı
t0 = max(Z["chi8e"][0], Z["chi8o"][0])
t1 = min(Z["chi8e"][-1], Z["chi8o"][-1])
kenar = np.exp(np.linspace(np.log(t0), np.log(t1 * 1.0001), 4))
print(f"  ortak aralık [{t0:.0f}, {t1:.0f}]; kenarlar "
      + " ".join(f"{k:.0f}" for k in kenar), flush=True)

# ortak KESİM t-aralıkları: iki adanın kusurlu bölgelerinin BİRLEŞİMİ
kesim_araliklari = []
for et, tab, a in ADA:
    zc, bad = Z[et], BAD[et]
    i = 0
    while i < len(zc):
        if bad[i]:
            j = i
            while j < len(zc) and bad[j]:
                j += 1
            kesim_araliklari.append((zc[i], zc[j - 1]))
            i = j
        else:
            i += 1
print(f"  birleşik kusur aralığı sayısı: {len(kesim_araliklari)}", flush=True)

SON = {}
for et, tab, a in ADA:
    zc = Z[et]
    kotu = np.zeros(len(zc), bool)
    for (x, y) in kesim_araliklari:
        kotu |= (zc >= x) & (zc <= y)
    seg, kes = segmentle(zc, kotu)
    WIN = []
    for s in seg:
        for i in range(3):
            pw = s[(s >= kenar[i]) & (s < kenar[i + 1])]
            if len(pw) >= 3000:
                WIN.append(pencere_hazirla(pw, 8.0))
    D, V = egri(WIN, TAUS)
    SON[et] = (D, V, [w[3] for w in WIN], [len(w[0]) for w in WIN])
    print(f"  {et}: {len(WIN)} pencere, L = "
          + " ".join(f"{x:.2f}" for x in SON[et][2])
          + ", n = " + " ".join(str(x) for x in SON[et][3]), flush=True)

print(f"\n{'τ':>6} {'chi8e':>7} {'chi8o':>7} {'Δ':>7} {'vekil(mx)':>9} "
      f"{'105d Δ':>8}")
D105 = {"chi8e": [0.122, 0.192, 0.227, 0.218, 0.168, 0.127, 0.524, 0.770,
                  0.859],
        "chi8o": [0.025, 0.031, 0.039, 0.047, 0.063, 0.110, 0.528, 0.671,
                  0.789]}
for i, t in enumerate(TAUS):
    de, do = SON["chi8e"][0][i], SON["chi8o"][0][i]
    vk = max(SON["chi8e"][1][i], SON["chi8o"][1][i])
    print(f"{t:>6.3f} {de:>7.3f} {do:>7.3f} {de-do:>+7.3f} {vk:>9.3f} "
          f"{D105['chi8e'][i]-D105['chi8o'][i]:>+8.3f}")

alt = slice(0, 6)
d_yeni = max(abs(SON["chi8e"][0][i] - SON["chi8o"][0][i]) for i in range(6))
d_eski = max(abs(D105["chi8e"][i] - D105["chi8o"][i]) for i in range(6))
print(f"\n  EŞİK-ALTI (τ≤0.085) maks |Δ|: 105d {d_eski:.3f} → "
      f"eşleştirilmiş {d_yeni:.3f}")
print("  HÜKÜM: " + ("ANOMALİ ALETSEL (segmentasyon/pencere) — P1 ayakta"
                     if d_yeni <= 0.05 else
                     "ANOMALİ KISMEN AZALDI ama sürüyor — açık kayıt"
                     if d_yeni < d_eski * 0.6 else
                     "ANOMALİ SÜRÜYOR — P1 eşik-altı dinamikte sorgulanır"))


# ---------------- K4: KISA-ÖLÇEK KUSUR KAPISI ----------------
# (İlk denenen "kümülatif sürüklenme" varyantı BAŞARISIZ oldu — eşik-altı
#  |Δ| 0.188 → 0.214; kayda geçiyor. Teşhis oradan sonra keskinleşti:)
# BLOK AYRIŞTIRMASI (t 11062-27368'i dörde bölerek):
#   chi8e t[15000,19500): D = 0.295/0.431/0.481/0.486/0.419  ← anomalinin
#          TAMAMI burada; komşu bloklar 0.008-0.086 (chi8o gibi)
#   ve o blokta ham d 2.5 basamak DÜŞÜP ~10 sıfır sonra KAPANIYOR
#   (t≈19352): 2 sıfır kaybı + hemen telafi → NET basamak yok.
# ÖLÇEK TEŞHİSİ (yuvarlanan medyan aralığı, aynı blok):
#   win=20 → 1.60   win=40 → 0.25   win=80 → 0.16  (taban 0.30/0.22/0.17)
# ⟹ 101'in win=80 dedektörü bu kusuru YIKIYOR. Kapı: kısa (20) ve uzun
#   (200) medyanların farkı > 0.7 → kusurlu. Yanlış-pozitif denetimi:
#   β'da maks fark 0.37 → 0 kusur; chi5e/chi8e/chi8o/chi3'te ada başına
#   yalnız 20-38 sıfır işaretleniyor.


def kisa_olcek_segmentler(zz, sayim, min_n=3000, win=80, esik=0.5, pad=160,
                          kisa_win=20, uzun_win=200, kisa_esik=0.7):
    d = np.arange(len(zz)) - (sayim(zz) - sayim(zz[0]))
    from numpy.lib.stride_tricks import sliding_window_view

    def cmed(x, w):
        m = np.median(sliding_window_view(x, w), axis=1)
        out = np.empty(len(x))
        out[w // 2:w // 2 + len(m)] = m
        out[:w // 2] = m[0]
        out[w // 2 + len(m):] = m[-1]
        return out

    med = np.median(sliding_window_view(d, win), axis=1)
    adim = med[win + 1:] - med[:-(win + 1)]
    bad = np.zeros(len(zz), bool)
    for i in np.where(np.abs(adim) > esik)[0] + win // 2:
        bad[max(0, i - pad):i + 2 * pad] = True
    sapma = np.abs(cmed(d, kisa_win) - cmed(d, uzun_win))
    print(f"    [{'':>0}maks |med20−med200| = {sapma.max():.2f}; "
          f"eşik {kisa_esik}]", flush=True)
    for i in np.where(sapma > kisa_esik)[0]:
        bad[max(0, i - pad):i + pad] = True
    return segmentle(zz, bad, min_n)


print("\n=== K4 — KISA-ÖLÇEK KUSUR KAPISI ===", flush=True)
SON4 = {}
for et, tab, a in ADA:
    zc = Z[et]
    seg, kes = kisa_olcek_segmentler(zc, lambda t: M_[et].theta(t) / np.pi)
    WIN = []
    for s in seg:
        for i in range(3):
            pw = s[(s >= kenar[i]) & (s < kenar[i + 1])]
            if len(pw) >= 3000:
                WIN.append(pencere_hazirla(pw, 8.0))
    D, V = egri(WIN, TAUS)
    SON4[et] = (D, V)
    print(f"  {et}: {len(seg)} segment ({kes} kesim) → {len(WIN)} pencere, "
          f"L = " + " ".join(f"{w[3]:.2f}" for w in WIN)
          + ", n = " + " ".join(str(len(w[0])) for w in WIN)
          + f"  (atılan sıfır: {len(zc) - sum(len(w[0]) for w in WIN)})",
          flush=True)

print(f"\n{'τ':>6} {'chi8e':>7} {'chi8o':>7} {'Δ':>7} {'vekil(mx)':>9} "
      f"{'105d Δ':>8} {'K2 Δ':>7}")
for i, t in enumerate(TAUS):
    de, do = SON4["chi8e"][0][i], SON4["chi8o"][0][i]
    vk = max(SON4["chi8e"][1][i], SON4["chi8o"][1][i])
    print(f"{t:>6.3f} {de:>7.3f} {do:>7.3f} {de-do:>+7.3f} {vk:>9.3f} "
          f"{D105['chi8e'][i]-D105['chi8o'][i]:>+8.3f} "
          f"{SON['chi8e'][0][i]-SON['chi8o'][0][i]:>+7.3f}")

d_k4 = max(abs(SON4["chi8e"][0][i] - SON4["chi8o"][0][i]) for i in range(6))
print(f"\n  EŞİK-ALTI (τ≤0.085) maks |Δ|: 105d {d_eski:.3f} → K2 {d_yeni:.3f}"
      f" → K4 {d_k4:.3f}")
print("  K4 HÜKMÜ: " + ("ANOMALİ ALETSEL (kısa-ölçek kusur) — P1 AYAKTA"
                        if d_k4 <= 0.05 else
                        "büyük ölçüde aletsel, artık var — açık kayıt"
                        if d_k4 < 0.5 * d_eski else
                        "ANOMALİ SÜRÜYOR — P1 eşik-altı dinamikte sorgulanır"))

np.savez(HERE / "105e_anomali.npz", taus=np.array(TAUS),
         **{f"D_{et}": np.array(SON[et][0]) for et in SON},
         **{f"V_{et}": np.array(SON[et][1]) for et in SON},
         **{f"D4_{et}": np.array(SON4[et][0]) for et in SON4},
         **{f"V4_{et}": np.array(SON4[et][1]) for et in SON4})
print("\n105e_anomali.npz yazıldı.", flush=True)
