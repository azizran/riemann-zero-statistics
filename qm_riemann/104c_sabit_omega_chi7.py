"""
104c — H-fizik: χ₇'NİN DONUK BÖLGESİNDE (ω = 0.20-0.65) YAPI VAR MI?
==========================================================================
101h'nin açık anomalisi: χ₇, τ=0.04'te band daraldıkça YÜKSELİYOR
(0.078 → 0.106 → 0.125), kaçak yok. 101i vuruş/moiré hipotezini
duyarlılığı kanıtlanmış bir aletle REDDETTİ ama yalnız ζ ve β'da
baktı — χ₇'nin kendi donuk bölgesi SABİT-ω ile HİÇ TARANMADI.
τ=0.04'te χ₇ pencerelerinin ω merkezleri (τ·L_w) 0.357-0.427; üç bandın
kapsadığı birleşik ω aralığı ±0.02L'de ~[0.17,0.62], ±0.005L'de
~[0.31,0.48]. Yani anomali ω ∈ 0.2-0.65 penceresinde yaşıyor.

ÇİZGİ ENVANTERİ bu aralıkta: log2 = 0.6931 aralığın ÜSTÜNDE; en yakın
alt çizgi yok (log q ≥ log2 tüm asal kuvvetler için). χ₇'ye özgü ÖLÜ
çizgi log7 = 1.9459 — bu aralıkta DEĞİL. Yani 0.20-0.65 aritmetik
olarak BOŞ bir penceredir; oradaki her yerel yapı çizgi-dışıdır.

YÖNTEM: 101i'nin sabit-ω estimatörü (pencere başına τ_w = ω/L_w ve
c_w), band ±0.004 (101i'nin çözünürlük kapısını geçen ayar),
40-permütasyonlu vekil taban. Veri 101f sertifikalı sıfırlar +
101d düzlük segmentasyonu; ζ = 41_bigT_windows ilk 6 pencere.

ÖN-MÜHÜR (ölçümden ÖNCE yazıldı):
  F1  ÇÖZÜNÜRLÜK KAPISI — ALET χ₇'DE DE GÖRÜYOR MU? χ₇, ω = log2 =
      0.6931 (CANLI) çevresinde keskin tepe vermeli (kontrast ≥ 1.5×,
      101i'de ζ 1.93× vermişti). Bu kapı geçilmezse null anlamsızdır.
  F2  ÖLÜ/CANLI ÇAPRAZI (χ₇'ye özgü, iki yönlü):
        ω = log7 = 1.9459 → χ₇ DÜZ (χ₇(7)=0), χ₃ ve χ₅ TEPE.
        ω = log3 = 1.0986 → χ₃ DÜZ (χ₃(3)=0), χ₅ ve χ₇ TEPE.
      Bu çapraz hem aleti hem veri hattını doğrular.
  F3  TARAMA: ω ∈ [0.20, 0.65], adım 0.01, dört aile (ζ, χ₃, χ₅, χ₇).
      H-fizik doğruysa χ₇'nin D(ω)'sinde pürüzsüz eğilimin üstünde
      YEREL TEPE olacak: artık z-skoru ≥ 3 (eğilim = log D'ye 3.
      dereceden polinom), ve aynı ω'da öteki üç aile TEPE VERMEYECEK.
  F4  Tepe bulunursa ince tarama (±0.03, adım 0.0025) ile frekansı
      raporlanacak. YORUM YOK — yalnız frekans, genlik, z.
  F5  H-fizik YANLIŞSA: χ₇'nin D(ω)'si pürüzsüz ve tekdüze artan
      olacak (101i'nin ζ/β'da bulduğu "çizgi omzu" vadisi gibi),
      χ₇ yalnızca ÖTEKİLERDEN YÜKSEK bir seviyede duracak.

==========================================================================
SONUÇ (25 Ağustos, koşu 230 s) — F1 ✓, F2 yarım (bir kolu bilgisiz),
F3 RET, F4 RET, F5 ✓  →  H-fizik REDDEDİLDİ; ve tarama BAĞIMSIZ OLARAK
kusur teşhisini destekledi.
==========================================================================
F1 ✓ ÇÖZÜNÜRLÜK KAPISI GEÇİLDİ — alet χ₇'de de çizgiyi görüyor
  (ω = log2 = 0.6931, band ±0.004, çizgi filtresi kapalı):
   aile   0.6531  0.6731  0.6931  0.7131  0.7331   kontrast
   ζ      0.3250  0.7691  1.0045  0.6838  0.4422    1.81×  TEPE
   χ₃     0.4592  0.7758  1.0045  0.8466  0.6446    1.47×  (sınırda)
   χ₅     0.4561  0.7088  1.0037  0.8617  0.4960    1.59×  TEPE
   χ₇     0.4724  0.7725  1.0050  0.8043  0.5442    1.55×  TEPE
  → χ₇'de bulunacak null GERÇEK bir null'dur.

F2 YARIM:
  ω = log3 = 1.0986 → χ₃ 1.03× DÜZ (χ₃(3)=0, ölü) ✓; ζ 1.21×, χ₅ 1.18×,
    χ₇ 1.16× — zayıf ama var. Ölü/canlı ayrımı yönü DOĞRU: χ₃ çizgide
    tepe yapmıyor, tekdüze tırmanıyor (0.417→0.543→0.607→0.684→0.710).
  ω = log7 = 1.9459 → DÖRT AİLE DE DÜZ (1.02-1.07×). Ön-mühürdeki
    "χ₃ ve χ₅ tepe verecek" beklentisi ÇIKMADI. Nedeni ölçümden sonra
    anlaşıldı (dürüst kayıt, ön-mühürde yoktu): ω=1.946'da τ = ω/⟨L⟩
    ≈ 0.19-0.20 ve orada BEŞ AİLE DE ÇOKTAN ÇÖZÜLMÜŞ (D ≈ 0.93-1.10),
    yani çizgi DOYMUŞ bir zemine düşüyor — kontrast fiziksel olarak
    imkânsız. Bu kapı ORADA BİLGİSİZDİR, ölü-çizgi hükmünün nulli
    DEĞİLDİR. (Kayıt: χ₇ o civarda TEK aile olarak 1'in ALTINDA —
    0.854/0.776/0.888/0.834/0.982; ötekiler 0.93-1.10. Yorumlanmadı.)

F3 RET — χ₇'NİN DONUK BÖLGESİNDE YEREL YAPI YOK:
  ω ∈ [0.20, 0.65], adım 0.01, band ±0.004; eğilim = log D'ye 3.
  dereceden polinom; artık z-skorları:
   aile  artık std(log D)   maks z @ ω          min z @ ω
   ζ        0.0699          +3.50 @ 0.650      −1.52 @ 0.610
   χ₃       0.4638          +1.83 @ 0.250      −2.47 @ 0.290
   χ₅       0.0558          +1.60 @ 0.410      −2.27 @ 0.260
   χ₇       0.3030          +2.00 @ 0.380      −2.28 @ 0.330
  χ₇'nin en güçlü artığı z = +2.00 (@ ω=0.380) — ön-mühürdeki z ≥ 3
  ölçütünün ALTINDA. ζ'nın +3.50'si tarama UCUNDA (0.650, log2'ye
  yaklaşma; polinom eğiliminin uç kusuru), yapı değil.
F4 RET — ince tarama (0.350-0.410, adım 0.0025, band ±0.002): χ₇'de
  tekrarlanabilir tepe yok, saçılma var (maks 0.2328 @ 0.3800, min
  0.0221 @ 0.3650 — komşu ω'lar arasında 10× sıçramalar).

F5 ✓ ve ASIL KAZANÇ — SAÇILMA KUSURUN İMZASI:
  D(ω) profili ζ ve χ₅'te PÜRÜZSÜZ (artık std 0.070 ve 0.056, 101i'nin
  "çizgi omzu vadisi" resmiyle uyumlu), χ₃ ve χ₇'de 4-8 KAT GÜRÜLTÜLÜ
  (0.464 ve 0.303). Gürültülü iki aile, 104b/104d'nin kusur taşıyan
  iki ailesiyle BİREBİR AYNI (χ₃ w1 ve χ₇ w1). Yani sabit-ω taraması,
  H-fizik'i reddederken H-kusur'u BAĞIMSIZ BİR KANALDAN doğruladı.
  Seviye farkı (yapı değil) sürüyor: χ₇/ζ 3.6-13.5×, D/V(χ₇) 3.7-26.6.

HÜKÜM (H-fizik): RET. ω ∈ [0.20,0.65] penceresinde χ₇'ye özgü hiçbir
yerel frekans yapısı yok (en güçlü aday z=+2.0 @ ω=0.380, ölçütün
altında ve ince taramada tekrarlanmıyor). Bulunan tek "yapı" ARTIK
SAÇILMASIDIR ve kaynağı fizik değil veri kusurudur.
"""

import numpy as np
import time
from pathlib import Path

T0 = time.time()
HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
rng = np.random.default_rng(1042)

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

def egri_omega(WIN, omegas, band=0.004, nfrek=160, cizgi_filtre=True):
    """101i birebir: sabit-ω, ortak ω-bandı, c pencere başına."""
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
            if len(oms) == 0:
                continue
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
    return np.array(Dv), np.array(Vv)

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

def rvm_sayim(t):
    x = t / TWO_PI
    return x * np.log(x / np.e)

d41 = np.load(HERE / "41_bigT_windows.npz")
K41 = sorted({x.split("_")[1] for x in d41.files}, key=lambda s: int(s[:-1]))
WINZ = []
for k in K41[:6]:
    gz, tm = d41[f"gaps_{k}"], d41[f"tmid_{k}"]
    zz = np.empty(len(gz) + 1)
    zz[0] = tm[0] - gz[0] / 2
    zz[1:] = zz[0] + np.cumsum(gz)
    seg, _ = duz_segmentler(zz, rvm_sayim)
    WINZ += [pencere_hazirla(s, 1.0) for s in seg]

exec(open(HERE / "98_L_motoru.py").read().split('CHI4 = ')[0])
HERE = Path(__file__).resolve().parent
CHI3 = {0: 0, 1: 1, 2: -1}
CHI5 = {0: 0, 1: 1, 2: 1j, 3: -1j, 4: -1}
z6c = np.exp(1j * np.pi / 3)
CHI7 = {0: 0, 1: 1, 3: z6c, 2: z6c**2, 6: z6c**3, 4: z6c**4, 5: z6c**5}

PEN = {"zeta": WINZ}
for et, q, tab in [("chi3", 3, CHI3), ("chi5", 5, CHI5), ("chi7", 7, CHI7)]:
    zc = np.load(HERE / f"101f_{et}_zeros.npz")["zeros"]
    lo = np.exp(np.log(zc[0] + 1) + 0.25 * (np.log(zc[-1]) - np.log(zc[0] + 1)))
    zc = zc[zc >= lo]
    M = Lmotor(q, tab, 1)
    seg, _ = duz_segmentler(zc, lambda t: M.theta(t) / np.pi)
    kenar = np.exp(np.linspace(np.log(zc[0]), np.log(zc[-1] * 1.0001), 4))
    W = []
    for s in seg:
        for i in range(3):
            p = s[(s >= kenar[i]) & (s < kenar[i + 1])]
            if len(p) >= 3000:
                W.append(pencere_hazirla(p, float(q)))
    PEN[et] = W

ADLAR = ["zeta", "chi3", "chi5", "chi7"]
for a in ADLAR:
    Ls = np.array([w[3] for w in PEN[a]])
    ns = np.array([len(w[0]) for w in PEN[a]], float)
    print(f"[{a}] {len(PEN[a])} pencere  n={int(ns.sum())}  "
          f"⟨L⟩_n={np.average(Ls, weights=ns):.3f}  "
          f"L∈[{Ls.min():.2f},{Ls.max():.2f}]", flush=True)

BK = 0.004

# ---------------- F1/F2: çözünürlük kapısı + ölü/canlı çaprazı --------
print("\n=== F1/F2: ÇÖZÜNÜRLÜK KAPISI ve ÖLÜ/CANLI ÇAPRAZI ===")
print("  (band ±0.004, çizgi filtresi KAPALI; tepe = çizgide/yan ort.)")
KAPI = {}
for adk, cizgi in [("log2", np.log(2)), ("log3", np.log(3)),
                   ("log7", np.log(7))]:
    gg = cizgi + np.array([-0.040, -0.020, 0.0, 0.020, 0.040])
    print(f"\n  --- ω = {adk} = {cizgi:.4f} ---")
    print(f"{'aile':>6} " + " ".join(f"{g:>8.4f}" for g in gg)
          + "   kontrast(çizgi/yan-ort)  durum")
    for a in ADLAR:
        D, V = egri_omega(PEN[a], list(gg), band=BK, cizgi_filtre=False)
        KAPI[(adk, a)] = (D, V)
        yan = float(np.mean([D[0], D[1], D[3], D[4]]))
        kon = D[2] / yan if yan > 0 else np.nan
        durum = "TEPE" if kon >= 1.5 else ("düz" if kon <= 1.15 else "ara")
        print(f"{a:>6} " + " ".join(f"{d:>8.4f}" for d in D)
              + f"   {kon:>8.2f}×  {durum}", flush=True)
    print(f"{'  vekil':>6} " + "  (ör. χ₇) " +
          " ".join(f"{v:>7.4f}" for v in KAPI[(adk, 'chi7')][1]))

# ---------------- F3: donuk bölge taraması ----------------
SCAN = np.round(np.arange(0.20, 0.6501, 0.01), 4)
print(f"\n=== F3: TARAMA ω ∈ [0.20, 0.65], adım 0.01, band ±{BK} ===")
TAR = {}
for a in ADLAR:
    TAR[a] = egri_omega(PEN[a], list(SCAN), band=BK)
    print(f"  ...{a} tarama bitti ({time.time()-T0:.0f} s)", flush=True)

print(f"\n{'ω':>7} " + " ".join(f"{a:>8}" for a in ADLAR)
      + "  |  " + " ".join(f"V:{a[:4]:>6}" for a in ADLAR))
for i, om in enumerate(SCAN):
    print(f"{om:>7.3f} " + " ".join(f"{TAR[a][0][i]:>8.4f}" for a in ADLAR)
          + "  |  " + " ".join(f"{TAR[a][1][i]:>8.4f}" for a in ADLAR))

print("\n=== ARTIK ANALİZİ (eğilim: log D'ye 3. derece polinom) ===")
ARTIK = {}
for a in ADLAR:
    y = np.log(np.maximum(TAR[a][0], 1e-12))
    P = np.polyfit(SCAN, y, 3)
    r = y - np.polyval(P, SCAN)
    ARTIK[a] = r / r.std()
    print(f"  {a}: artık std(log D) = {r.std():.4f}  "
          f"maks z = {ARTIK[a].max():+.2f} @ ω={SCAN[ARTIK[a].argmax()]:.3f}  "
          f"min z = {ARTIK[a].min():+.2f} @ ω={SCAN[ARTIK[a].argmin()]:.3f}")
print(f"\n{'ω':>7} " + " ".join(f"z:{a:>7}" for a in ADLAR))
for i, om in enumerate(SCAN):
    print(f"{om:>7.3f} " + " ".join(f"{ARTIK[a][i]:>+9.2f}" for a in ADLAR))

print("\n  χ₇/öteki oranları (seviye farkı, yapı değil):")
for i in range(0, len(SCAN), 5):
    om = SCAN[i]
    print(f"   ω={om:.2f}: χ₇/ζ {TAR['chi7'][0][i]/TAR['zeta'][0][i]:.2f}× "
          f" χ₇/χ₃ {TAR['chi7'][0][i]/TAR['chi3'][0][i]:.2f}× "
          f" χ₇/χ₅ {TAR['chi7'][0][i]/TAR['chi5'][0][i]:.2f}×"
          f"  D/V(χ₇) {TAR['chi7'][0][i]/TAR['chi7'][1][i]:.1f}")

# ---------------- F4: ince tarama (en güçlü χ₇ artığı) ----------------
om_pik = float(SCAN[ARTIK["chi7"].argmax()])
z_pik = float(ARTIK["chi7"].max())
print(f"\n=== F4: İNCE TARAMA χ₇ en güçlü artık çevresi ω*={om_pik:.3f} "
      f"(z={z_pik:+.2f}) ===")
INCE = np.round(np.arange(om_pik - 0.03, om_pik + 0.0301, 0.0025), 5)
FIN = {}
for a in ADLAR:
    FIN[a] = egri_omega(PEN[a], list(INCE), band=0.002)
    print(f"  ...{a} ince tarama bitti ({time.time()-T0:.0f} s)", flush=True)
print(f"{'ω':>8} " + " ".join(f"{a:>8}" for a in ADLAR)
      + "  |  V(χ₇)")
for i, om in enumerate(INCE):
    print(f"{om:>8.4f} " + " ".join(f"{FIN[a][0][i]:>8.4f}" for a in ADLAR)
          + f"  |  {FIN['chi7'][1][i]:>7.4f}")
for a in ADLAR:
    y = FIN[a][0]
    print(f"  {a}: ince tarama maks {y.max():.4f} @ ω={INCE[y.argmax()]:.4f}, "
          f"min {y.min():.4f}, tepe/ort {y.max()/y.mean():.2f}×")

np.savez(HERE / "104c_omega.npz", scan=SCAN, ince=INCE, band=BK,
         **{f"D_{a}": TAR[a][0] for a in ADLAR},
         **{f"V_{a}": TAR[a][1] for a in ADLAR},
         **{f"z_{a}": ARTIK[a] for a in ADLAR},
         **{f"F_{a}": FIN[a][0] for a in ADLAR},
         **{f"kapi_{k[0]}_{k[1]}": KAPI[k][0] for k in KAPI})
print(f"\n104c_omega.npz yazıldı. Süre {time.time()-T0:.0f} s.")
