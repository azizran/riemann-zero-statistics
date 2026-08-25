"""
105d — YENİ ADALARDA DONMA EĞRİLERİ (25 Ağustos 2026)
==========================================================================
101'in TAÇ BULGUSU ("donma sınırı adanın İLK SAĞ KALAN ÇİZGİSİDİR")
İKİ BAĞIMSIZ ADADA replike edilir; ayrıca parite/gamma'nın hükme
girmediği sınanır.

MAKİNE: 101d estimatörü (yerel katlama + kübik detrend + havuzlanmış
D(τ) + 40-permütasyon vekil taban) + DÜZLÜK SEGMENTASYONU (101b kusur
kapısının gereği) + 105e'nin KISA-ÖLÇEK KUSUR KAPISI (bkz. duz_segmentler
docstring'i: 101'in 80'lik basamak dedektörü kısa çukur kusurları
göremiyor; ilk koşuda chi8e'nin eşik-altı D'sini 0.02→0.48 şişirdi)
+ 101h dersi: BAND ±0.01L (±0.02L artefakt taşıyor)
ve BAND KAÇAĞI KAPISI — bandı ilk çizgiyi kapsayan pencereler her
τ'da SAYILIR ve raporlanır (sızan ölçüm hükme sokulmaz).

AİLELER: chi5e, chi8e, chi8o (105b) + β, χ₃ (101f, kıyas).
τ ızgarası 101d'ninki.

ÖN-MÜHÜRLÜ ÖNGÖRÜ P3 (kampanya mühründen, en keskin):
  chi5e  — 2-çizgisi CANLI (χ(2)=−1≠0) → log2/L'de ÇÖZÜLÜR
           (χ₃/χ₅/χ₇/ζ gibi; τ_ilk = log2/⟨L⟩ ≈ 0.07)
  chi8e  — 2-ailesi ÖLÜ → log3/L'e GECİKİR (β gibi; τ_ilk ≈ 0.115)
  chi8o  — 2-ailesi ÖLÜ → log3/L'e GECİKİR (β gibi)
  TUTARSA: donma sınırı YALNIZ ÇİZGİ ENVANTERİNİN fonksiyonudur —
  paritenin (chi8e a=0 vs chi8o a=1 AYNI davranmalı), gamma faktörünün
  ve iletkenin (q=4 vs q=8, ikisi de 2-ailesiz) DEĞİL.
  RET KOŞULU: chi8e ile chi8o birbirinden ayrışırsa (parite hükme
  giriyor demektir) veya biri log2/L'de çözülürse.

==========================================================================
SONUÇ (25 Ağustos) — P3 İSABET; 101'İN TAÇ BULGUSU İKİ ADADA REPLİKE
==========================================================================
(İLK koşu chi8e'de eşik-altı bir anomali gösterdi; 105e onu KISA-ÖLÇEK
VERİ KUSURU olarak teşhis etti ve kapı buraya taşındı. Aşağıdaki
sayılar kapılı koşudan.)

D(τ) — band ±0.01L (vekil taban ≤0.10 her yerde):
  τ      chi5e   chi8e   chi8o    beta    chi3
  0.040  0.045   0.005   0.025   0.012   0.036
  0.060  0.235   0.009   0.039   0.026   0.124
  0.068  0.425   0.009   0.047   0.037   0.194
  0.075  0.591   0.009   0.063   0.053   0.597
  0.085  0.538   0.029   0.110   0.091   0.617
  0.095  0.356   0.486   0.528   0.178   0.426
  0.105  0.704   0.767   0.671   0.647   0.199
  0.113  0.737   0.843   0.789   0.726   0.161
  0.300  1.038   1.068   1.040   1.067   1.023

HÜKÜM TABLOSU:
  aile   q_ilk  τ_ilk    D@log2/L  D@log3/L  yarı-çöz τ   ω_yarı
  chi5e    2   0.0708      0.492     0.734     0.0712     0.697
  chi8e    3   0.1085      0.009     0.800     0.0955     0.967
  chi8o    3   0.1085      0.048     0.722     0.0943     0.955
  beta     3   0.1143      0.046     0.742     0.1019     0.979
  chi3     2   0.0735      0.508     0.153     0.0733     0.692

P3 ✓✓✓ ÜÇ MADDESİ DE TUTTU:
 (a) chi5e (2-çizgisi CANLI, a=0) log2/L'de ÇÖZÜLÜR: yarı-çözülme
     τ=0.0712 vs τ_ilk=0.0708 — ORAN 1.006. χ₃ ile birebir (0.0733
     vs 0.0735). Çift karakter, tek karakterler gibi.
 (b) chi8e VE chi8o İKİSİ DE log3/L'e GECİKİR: log2/L'de D = 0.009
     (chi8e — VEKİL TABANINDA, 0.003!) ve 0.048 (chi8o); log3/L'de
     0.800 / 0.722. Ölü çizgide çözülme YOK.
 (c) Kayma ×1.376-1.410 (β 1.410; öngörü log3/log2 = 1.585; 101d
     ölçümü ×1.4-1.5 — aynı bant).
 τ=0.30'da beş aile 1.02-1.07: tek CUE.

PARİTE KONTROLÜ ✓ chi8e (a=0) vs chi8o (a=1): yarı-çözülme 0.0955 vs
0.0943 (%1.3); eşik-altı maks |Δ| = 0.081, eşik-üstü ≤0.12 (vekil
gürültüsü mertebesinde). AYNI ÇİZGİ ENVANTERİ = AYNI DONMA SINIRI,
parite ve gamma faktörü hükme GİRMİYOR.

DÜRÜST KAYITLAR:
 • chi3'ün τ=0.105-0.140 çukuru (0.199/0.161/0.134/0.531) yeni kapıyla
   DERİNLEŞTİ (kapısız: 0.565/0.512/0.320/0.766). 101d/101h'de de
   vardı (101h ±0.005L'de τ=0.125 → 0.232). chi3'ün orta-τ değerleri
   pencere bölmesine kırılgan; EŞİK hükmü (log2/L) etkilenmiyor.
 • chi5e'nin τ=0.095'teki 0.356 düşüşü de aynı sınıf salınım.
 • chi8e 4 pencereye bölündü (2 kesim), chi8o 2 (0 kesim) — eşit
   olmayan bölme kalıntı bir sistematik bırakabilir (105e K2 bunu
   eşleştirilmiş bölmeyle sınadı: hüküm değişmiyor).
 • Band kaçağı kapısı: eşik-altı τ'ların HİÇBİRİNDE mod-8 adalarında
   sızıntı yok (log3 banda τ≥0.095'te giriyor); chi5e ve chi3'te
   τ≥0.060/0.068'de sızıntı var ve ★ ile işaretli.
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


def _cmed(x, w):
    """Ortalanmış yuvarlanan medyan (kenarlar uçtaki değerle doldurulur)."""
    from numpy.lib.stride_tricks import sliding_window_view
    m = np.median(sliding_window_view(x, w), axis=1)
    out = np.empty(len(x))
    out[w // 2:w // 2 + len(m)] = m
    out[:w // 2] = m[0]
    out[w // 2 + len(m):] = m[-1]
    return out


def duz_segmentler(zz, sayim, min_n=3000, win=80, esik=0.5, pad=160,
                   kisa_win=20, uzun_win=200, kisa_esik=0.7):
    """101d'nin düzlük segmentasyonu + 105e'nin KISA-ÖLÇEK KUSUR KAPISI.

    101'in dedektörü 80'lik medyan farkıyla BASAMAK arar; 105e, chi8e'de
    ~10 sıfır boyunca süren ve sonra kapanan bir ÇUKUR kusurunun (2 sıfır
    kaybı, t≈19352) 80'lik medyanda tamamen yıkandığını (0.16 vs taban
    0.17) ama donma estimatörünü 0.02'den 0.48'e fırlattığını gösterdi.
    Çare: kısa (20) ve uzun (200) yuvarlanan medyanların farkı > 0.7 olan
    yerler de KUSURLU sayılır. β'da 0 kusur (maks fark 0.37) verir —
    yani kapı yanlış-pozitif üretmiyor.
    """
    d = np.arange(len(zz)) - (sayim(zz) - sayim(zz[0]))
    from numpy.lib.stride_tricks import sliding_window_view
    med = np.median(sliding_window_view(d, win), axis=1)
    adim = med[win + 1:] - med[:-(win + 1)]
    bad = np.zeros(len(zz), bool)
    for i in np.where(np.abs(adim) > esik)[0] + win // 2:
        bad[max(0, i - pad):i + 2 * pad] = True
    # --- kısa-ölçek kusur kapısı (105e) ---
    sapma = np.abs(_cmed(d, kisa_win) - _cmed(d, uzun_win))
    for i in np.where(sapma > kisa_esik)[0]:
        bad[max(0, i - pad):i + pad] = True
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
CHI5E = {0: 0, 1: 1, 2: -1, 3: -1, 4: 1}
CHI8E = {0: 0, 1: 1, 2: 0, 3: -1, 4: 0, 5: -1, 6: 0, 7: 1}
CHI8O = {0: 0, 1: 1, 2: 0, 3: 1, 4: 0, 5: -1, 6: 0, 7: -1}

AILELER = [("chi5e", 5, CHI5E, 0, "105b_chi5e_zeros.npz", 2),
           ("chi8e", 8, CHI8E, 0, "105b_chi8e_zeros.npz", 3),
           ("chi8o", 8, CHI8O, 1, "105b_chi8o_zeros.npz", 3),
           ("beta",  4, CHI4,  1, "101f_beta_zeros.npz",  3),
           ("chi3",  3, CHI3,  1, "101f_chi3_zeros.npz",  2)]

PEN, QILK, KES = {}, {}, {}
ADLAR = []
for etiket, q, tab, a, dosya, qilk in AILELER:
    p = HERE / dosya
    if not p.exists():
        print(f"[{etiket}] {dosya} YOK — atlanıyor", flush=True); continue
    zc = np.load(p)["zeros"]
    lo = np.exp(np.log(zc[0] + 1) + 0.25 * (np.log(zc[-1]) - np.log(zc[0] + 1)))
    zc = zc[zc >= lo]
    M = Lmotor(q, tab, a)
    seg, kes = duz_segmentler(zc, lambda t: M.theta(t) / np.pi)
    kenar = np.exp(np.linspace(np.log(zc[0]), np.log(zc[-1] * 1.0001), 4))
    WIN = []
    for s in seg:
        for i in range(3):
            pw = s[(s >= kenar[i]) & (s < kenar[i + 1])]
            if len(pw) >= 3000:
                WIN.append(pencere_hazirla(pw, float(q)))
    PEN[etiket] = WIN; QILK[etiket] = qilk; KES[etiket] = kes
    ADLAR.append(etiket)
    print(f"[{etiket}] n={len(zc)}, {len(WIN)} segment-pencere ({kes} kesim), "
          f"L'ler: " + " ".join(f"{w[3]:.2f}" for w in WIN), flush=True)

LORT, TILK = {}, {}
for a in ADLAR:
    Ls = np.array([w[3] for w in PEN[a]])
    ns = np.array([len(w[0]) for w in PEN[a]], float)
    LORT[a] = float(np.average(Ls, weights=ns))
    TILK[a] = np.log(QILK[a]) / LORT[a]
    print(f"  {a}: ⟨L⟩_n = {LORT[a]:.3f}, q_ilk = {QILK[a]}, "
          f"τ_ilk = {TILK[a]:.5f}", flush=True)

TAUS = [0.04, 0.05, 0.06, 0.068, 0.075, 0.085, 0.095, 0.105, 0.113,
        0.125, 0.14, 0.16, 0.20, 0.30]
BAND = 0.01

SONUC = {}
for a in ADLAR:
    SONUC[a] = egri(PEN[a], TAUS, band=BAND)
    print(f"  ...{a} bitti", flush=True)

print(f"\n=== D(τ), BAND ±{BAND}L ===", flush=True)
print(f"{'τ':>6} " + " ".join(f"{a:>7}" for a in ADLAR) + "   vekil(maks)")
for i, t in enumerate(TAUS):
    vek = max(SONUC[a][1][i] for a in ADLAR)
    print(f"{t:>6.3f} " + " ".join(f"{SONUC[a][0][i]:>7.3f}" for a in ADLAR)
          + f"   {vek:>7.3f}")
print("\nvekil taban (aile aile):")
for a in ADLAR:
    print(f"  {a:>6}: " + " ".join(f"{v:.3f}" for v in SONUC[a][1]))

print(f"\n=== BAND KAÇAĞI KAPISI: SIZAN PENCERE SAYISI (band ±{BAND}L) ===")
print("  pencere w sızar ⇔ L_w·(τ+band) ≥ log q_ilk")
print(f"{'aile':>6} {'n_pen':>5} " + " ".join(f"{t:>6.3f}" for t in TAUS))
SIZ = {}
for a in ADLAR:
    Ls = np.array([w[3] for w in PEN[a]])
    ln = np.log(QILK[a])
    s = [int(np.sum(Ls * (t + BAND) >= ln)) for t in TAUS]
    SIZ[a] = s
    print(f"{a:>6} {len(Ls):>5} " + " ".join(f"{x:>6d}" for x in s))

print("\n=== KAÇAKSIZ TABLO (sızan τ'lar ★ ile işaretli) ===")
print(f"{'τ':>6} " + " ".join(f"{a:>8}" for a in ADLAR))
for i, t in enumerate(TAUS):
    huc = []
    for a in ADLAR:
        mark = "★" if SIZ[a][i] > 0 else " "
        huc.append(f"{SONUC[a][0][i]:>7.3f}{mark}")
    print(f"{t:>6.3f} " + " ".join(huc))

XG = np.arange(0.60, 1.601, 0.05)
print("\n=== YENİDEN ÖLÇEKLEME: D(x), x = τ/τ_ilk ===")
print(f"{'x':>6} " + " ".join(f"{a:>7}" for a in ADLAR) + "    std  std/ort")
MM = {}
for a in ADLAR:
    xs = np.array(TAUS) / TILK[a]
    MM[a] = np.interp(XG, xs, SONUC[a][0], left=np.nan, right=np.nan)
sp, ab = [], []
for j, x in enumerate(XG):
    vals = np.array([MM[a][j] for a in ADLAR])
    s = float(np.nanstd(vals)); m = float(np.nanmean(vals))
    sp.append(s / m if m > 0 else np.nan); ab.append(s)
    print(f"{x:>6.2f} " + " ".join(f"{v:>7.3f}" for v in vals)
          + f"  {s:>6.3f}  {s/m if m>0 else np.nan:>7.3f}")
sp, ab = np.array(sp), np.array(ab)
for k, msk in [("x≤0.80", XG <= 0.801), ("0.85≤x≤1.10",
               (XG >= 0.849) & (XG <= 1.101)), ("x≥1.15", XG >= 1.149),
               ("tümü", np.ones(len(XG), bool))]:
    print(f"  ÇÖKME {k:>12}: std/ort {np.nanmean(sp[msk]):.4f}   "
          f"mutlak std {np.nanmean(ab[msk]):.4f}")

print("\n=== P3 HÜKÜM TABLOSU ===")
print(f"{'aile':>6} {'q_ilk':>5} {'τ_ilk':>7} {'D@log2/L':>9} "
      f"{'D@log3/L':>9} {'D@τ_ilk':>8} {'yarı-çöz τ':>11} {'D@0.30':>7}")
YARI = {}
for a in ADLAR:
    Ls = LORT[a]
    t2, t3 = np.log(2) / Ls, np.log(3) / Ls
    D = np.array(SONUC[a][0]); T = np.array(TAUS)
    d2 = float(np.interp(t2, T, D)); d3 = float(np.interp(t3, T, D))
    di = float(np.interp(TILK[a], T, D))
    # yarı-çözülme: D'nin 0.5'i ilk geçtiği τ (lineer ara değer)
    ix = np.where(D >= 0.5)[0]
    if len(ix) and ix[0] > 0:
        i0 = ix[0]
        ty = T[i0 - 1] + (0.5 - D[i0 - 1]) * (T[i0] - T[i0 - 1]) / \
            (D[i0] - D[i0 - 1])
    else:
        ty = np.nan
    YARI[a] = ty
    print(f"{a:>6} {QILK[a]:>5} {TILK[a]:>7.4f} {d2:>9.3f} {d3:>9.3f} "
          f"{di:>8.3f} {ty:>11.4f} {D[-1]:>7.3f}")

print("\n  yarı-çözülme kayması (2-ailesiz / 2-canlı), öngörü log3/log2=1.585:")
canli = [a for a in ADLAR if QILK[a] == 2 and not np.isnan(YARI[a])]
olu = [a for a in ADLAR if QILK[a] == 3 and not np.isnan(YARI[a])]
if canli and olu:
    yc = np.mean([YARI[a] * LORT[a] for a in canli])   # ω cinsinden
    for a in olu:
        print(f"    {a:>6}: ω_yarı = {YARI[a]*LORT[a]:.4f} vs canlı ort "
              f"{yc:.4f} → ×{YARI[a]*LORT[a]/yc:.3f}", flush=True)

print("\n  PARİTE KONTROLÜ (chi8e a=0 vs chi8o a=1 — AYNI olmalı):")
if "chi8e" in SONUC and "chi8o" in SONUC:
    de = np.array(SONUC["chi8e"][0]); do = np.array(SONUC["chi8o"][0])
    print(f"{'τ':>6} {'chi8e':>7} {'chi8o':>7} {'Δ':>7}")
    for i, t in enumerate(TAUS):
        print(f"{t:>6.3f} {de[i]:>7.3f} {do[i]:>7.3f} {de[i]-do[i]:>+7.3f}")
    print(f"  maks |Δ| = {np.abs(de-do).max():.3f}; "
          f"eşik-altı (τ≤0.085) maks |Δ| = {np.abs(de-do)[:6].max():.3f}")

np.savez(HERE / "105d_donma.npz", taus=np.array(TAUS), band=BAND,
         adlar=np.array(ADLAR),
         Lort=np.array([LORT[a] for a in ADLAR]),
         tilk=np.array([TILK[a] for a in ADLAR]),
         qilk=np.array([QILK[a] for a in ADLAR]),
         **{f"D_{a}": np.array(SONUC[a][0]) for a in ADLAR},
         **{f"V_{a}": np.array(SONUC[a][1]) for a in ADLAR},
         **{f"SIZ_{a}": np.array(SIZ[a]) for a in ADLAR})
print("\n105d_donma.npz yazıldı.", flush=True)
