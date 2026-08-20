"""
102a — KARANLIK ALANIN TAPERLİ YENİDEN ÖLÇÜMÜ (86-T4 DENETİMİ) (20 Ağustos)
==========================================================================
DENETİM işi. 101j, ζ'nın log2 çevresindeki "karanlık taban"ının TAPERSİZ
ölçümde tamamen pencere sızıntısı olduğunu gösterdi (kusursuz örgüden
ayırt edilemiyor). 86-T4'ün hiperuniformluk sayısı ("çizgi-dışı taban
2e-4 = 1/√n'in 25 KAT altı, karıştırılmış-gap vekilinin de altı") aynı
tapersiz estimatörle kurulmuştu. Bu script o ölçümü 101j'nin taper
makinesi + kalibrasyon kapılarıyla yeniden yapar.

BİRİM: I(ω) = |Σ_n w_n e^{iω t_n}|² / Σ_n w_n².  Bu, taper'dan BAĞIMSIZ
doğru normalizasyondur:
  • w=1 iken I = n|Ĝ|² (86/88'in F'i birebir),
  • Poisson süreci → I ≈ 1 (atım gürültüsü birimi),
  • bin-ortalaması Montgomery rampasına (F=α) oturur (türetim: çizgi
    toplamı Σ_q|A_q|²·L/W = ω/L; herhangi bir w için aynı).
86'nın |Ĝ| birimine dönüş: |Ĝ| = √(I/n).

ÖLÇÜMLER
  G0 SAYISAL TABAN KAPISI: faz büyüklüğü t~2e5 olduğu için exp(iωt)
     mutlak faz hatası ~1e-10 rad. Her toplamda t−t_merkez kullanılır
     (|Ĝ| değişmez, faz hatası ~16× küçülür). Kusursuz örgü + rastgele
     ω ile sayısal taban ölçülür; ölçülen tabanlar bunun ÜSTÜNDE mi?
  G1 KALİBRASYON: Poisson (I≈1) ve kusursuz örgü (I≈0), taperli/tapersiz.
  A1 86-T4 TEKRARI: 86'nın tam kurulumu (200k penceresi, L=10.37,
     ω∈[0.15,4.2] 1200 nokta, çizgilere >0.015), beş ızgara (gerçek,
     karıştırılmış-gap, RvM pürüzsüz, kusursuz örgü, Poisson), iki ayar
     (tapersiz / Hann). 86'nın bantları: hepsi, ω<1, ω>3.
  A2 DUYARLILIK MERDİVENİ (kritik kapı): gerçek ızgaraya ε genliğinde
     BAĞIMSIZ (iid) Gauss titreşimi enjekte edilir. Einstein/DW
     öngörüsü: I_difüz = 1 − e^{−ω²ε²}. Estimatörün kör OLMADIĞI ve
     hangi düzeyden itibaren gördüğü böyle kalibre edilir → "taban
     ≥X× rampanın altında" alt sınırı ALET-SINIRLI değil ÖLÇÜLMÜŞ olur.
  A3 ALTI PENCERE: aynı ölçüm 6 pencerede (L=9.86…12.45) özet tablo.
  A4 RAMPA KIYASI: bant ortasında α=ω/L rampası ile ölçülen taban oranı.
  A5 [EK, ön-mühür SONRASI eklendi — sebebi aşağıda dürüstçe kayıtlı]
     IZGARA MERDİVENİ: 101j sıfırların KENDİSİNDE (z_n) taperli tabanı
     ~1e-15 ölçmüştü; A1 ise ORTA NOKTALARDA (t_n = z_n + g_n/2, yani
     86/87/88'in ızgarası) ~1e-6…1e-2 buluyor. Aradaki 9-10 kademe
     ölçüm hatası değil, İKİ FARKLI NESNE. Kimlik:
        Ĝ_orta(ω) = (1/n) Σ_n e^{iωz_n} · e^{iωg_n/2}
     yani orta-nokta ızgarası, sıfır örgüsünü DALGALANAN yarım-gap
     kaymasıyla örnekler. ⟨e^{iωg/2}⟩ = φ(ω) bir Debye-Waller çarpanı
     (GAP dağılımının), dalgalanma ise difüz saçılma verir. Merdiven:
       z (sıfırlar) | z+ḡ/2 (sabit kayma, kontrol) | z+g_n/2 (86'nın
       ızgarası) | z+c(g_n−ḡ)/2, c=0.25/0.5/2 (c² kapısı) |
       z+g_karıştırılmış/2 (gap korelasyonları silinir → beyaz sınır).
     Beyaz üst sınır: 1 − |φ(ω)|².

ÖN-MÜHÜR (ölçümden ÖNCE yazıldı):
  S1 86'nın tapersiz tabanı (I≈1.6e-3, |Ĝ|≈2e-4) kusursuz-örgü
     tapersiz tabanıyla AYNI MERTEBEDE çıkacak → 86'nın "25×" sayısı
     fizik değil alet (101j'nin log2 civarındaki bulgusunun geniş-bant
     doğrulaması). Yön (hiperuniform = karanlık) sağlam kalacak.
  S2 Taperli tabanda gerçek ızgara ile kusursuz örgü yine ayırt
     edilemeyecek (I ≲ 1e-12), yani ÖLÇÜLEBİLİR fiziksel difüz bileşen
     ÇIKMAYACAK — açık formülün atomik tayfıyla tutarlı.
  S3 Duyarlılık merdiveni: ε ≥ 0.001 seviyesindeki iid titreşim rahatça
     görülecek; dolayısıyla gerçek ızgaranın çizgi-dışı iid-titreşim
     içeriği için ε ≲ 1e-3 üst sınırı verilebilecek.
  RİSK: taperli ölçümde ω>3 bandında (DW zayıflaması büyür) gerçek
     ızgaranın örgü üstünde KALMASI mümkün — o zaman S2 kısmen düşer ve
     ilk gerçek difüz bileşen orada bulunur. Bu, işin ilginç sonucu olur.

==========================================================================
SONUÇ (20 Ağustos, koşu 43 s) — S1 KISMEN ✓, S2 RET, S3 ✓, RİSK GERÇEKLEŞTİ
KARANLIK ALANIN SAHİBİ SIFIR ÖRGÜSÜ DEĞİL, ORTA-NOKTA IZGARASIDIR.
==========================================================================
KALİBRASYON ✓: Poisson I = 0.68–0.75 (≈1 ✓); kusursuz örgü tapersiz
  3.3e-5, Hann 5.4e-21 (sayısal taban); RvM pürüzsüz ≡ kusursuz örgü.
  86'nın sayıları BİREBİR yeniden üretildi: gerçek |Ĝ| = 1.53e-4
  (86: 2e-4), karıştırılmış-gap 1.45e-3 (86: 1.4e-3), 1/√n = 5.0e-3.

A1 — 86-T4 TEKRARI (200k penceresi, I birimi; |Ĝ| = √(I/n)):
  ızgara / ayar              hepsi        ω<1        ω>3
  gerçek ORTA, tapersiz    9.37e-04   4.19e-04   5.35e-03
  gerçek ORTA, Hann        4.21e-04   4.12e-06   4.38e-03
  SIFIRLAR z_n, tapersiz   2.82e-04   4.42e-04   3.33e-04
  SIFIRLAR z_n, Hann       2.38e-14   7.29e-16   6.11e-13
  karıştırılmış-gap, Hann  7.99e-02   1.24e-01   4.05e-02
  kusursuz örgü, Hann      5.44e-21   8.54e-18   1.41e-22
  ★ (1) 86'NIN TABANI SIZINTI DÜZEYİNDE MİYDİ? KISMEN. ω<1'de EVET
    (gerçek/kusursuz-örgü = 1.20×, ayırt edilemiyor — 101j'nin log2
    bulgusunun geniş-bant doğrulaması); ama bütün bantta 28×, ω>3'te
    293× sızıntının ÜSTÜNDE. 86'nın "25×" sayısı yine de yanlış, ama
    BEKLENENİN TERSİ yönde: taperli tabanda karanlık DAHA DERİN —
    tüm bantta 1/√n'in 48 katı, ω<1'de 490 katı altında.
  ★ (2) KUSURSUZ ÖRGÜ ÜSTÜNDE FİZİKSEL BİLEŞEN VAR MI? VAR, DEVASA:
    taperli gerçek/örgü = 4.8e11 (ω<1), 7.7e16 (hepsi), 3.1e19 (ω>3).
    AMA SIFIRLARIN KENDİSİNDE taban 7.3e-16 … 6.1e-13 — orta noktalarla
    arasında 8–10 KADEME fark (A5).

A2 — DUYARLILIK MERDİVENİ ✓ (estimatör kör değil): iid ε enjeksiyonunda
  ölçülen artık / ω²ε² oranı ε ≥ 3e-3'te 0.55–1.00 (medyan/ortalama
  farkı ln2 = 0.69 çarpanını açıklar). Orta-nokta ızgarasında saptama
  eşiği ε ≈ 3e-3 (kendi tabanı yüzünden); SIFIR ızgarasında (A2b)
  ε = 1e-7 bile oran 0.84–1.10 ile görülüyor.

A2c — SIFIR TABANI VERİ-SINIRLI (dürüstlük kapısı): 41_bigT_scan.py'nin
  18 bisection'ı sıfırları ±5.78e-8 ile veriyor (σ_kuant = 3.34e-8).
  Bu kuantizasyon ENJEKTE edilince sıfır tabanı ω<1'de 1.52×, tüm
  bantta 1.28× artıyor → ölçülen 7.3e-16 fiziksel değil VERİ tabanıdır;
  gerçek sıfır-örgüsü karanlık alanı bundan da koyu olabilir. ÜST SINIR.

A5 — ★ SEFERİN EN BÜYÜK BULGUSU: KARANLIK ALAN KİMİN?
  Kimlik: Ĝ_orta(ω) = (1/n) Σ_n e^{iωz_n}·e^{iωg_n/2}. Orta noktalar,
  sıfır örgüsünü DALGALANAN yarım-gap kaymasıyla örnekler; ⟨e^{iωg/2}⟩
  bir Debye-Waller çarpanı, dalgalanması ise difüz saçılmadır.
  ÖRNEKLEME-FAZI EĞRİSİ  t(c) = z_n + ḡ/2 + c·(g_n−ḡ)/2, I(hepsi):
    c=0.00  2.41e-14   (= z_n + sabit kayma → SIFIR ÖRGÜSÜ)
    c=0.25  7.78e-05
    c=0.50  2.32e-04
    c=1.00  4.21e-04   ← 86/87/88'in ızgarası (orta noktalar)
    c=1.50  2.42e-04
    c=2.00  2.41e-14   (= z_{n+1} − ḡ/2 → YİNE SIFIR ÖRGÜSÜ)
  Karanlık alan, İKİ ÖZDEŞ ve 10 KADEME DAHA KOYU örgü arasında salınan
  bir ÖRNEKLEME FAZI etkisidir; c=1 (orta nokta) tam tepede.
  Gap korelasyonları silinince (karıştırılmış g): 2.65e-02 (63× yukarı).
  Beyaz üst sınır 1 − |⟨e^{iωg/2}⟩|²: ω<1 4.72e-3, hepsi 5.23e-2,
  ω>3 1.73e-1 → gerçek gap dizisi bu sınırın 1146× / 124× / 39× ALTINDA.
  ★ HİPERUNİFORMLUĞUN DÜRÜST ÖLÇÜSÜ BUDUR: "taban 1/√n'in 25 katı
  altında" değil, "gap-difüz kanalı bağımsız-gap vekilinin 40–1150 katı
  altında" (+ sıfır örgüsünün kendisi ölçüm sınırına kadar TAM KARANLIK).

A4 — RAMPA KIYASI (α = ω/L; alet-sınırlı ALT SINIRLAR):
  bant   ⟨α⟩   orta-Hann I    α/I      SIFIR-Hann I     α/I
  hepsi 0.190   4.21e-04    4.5e+02      2.38e-14    8.0e+12
  ω<1   0.055   4.12e-06    1.3e+04      7.29e-16    7.5e+13
  ω>3   0.342   4.38e-03    7.8e+01      6.11e-13    5.6e+11
  DÜRÜST İFADE: "orta-nokta ızgarasının çizgi-dışı tabanı rampanın en
  az 78× (ω>3) – 1.3e4× (ω<1) altındadır; SIFIR örgüsünün çizgi-dışı
  tabanı rampanın en az 5.6e11× – 7.5e13× altındadır (veri-sınırlı)."

A3 — ALTI PENCEREDE AYNI: orta-Hann 3.1e-4 … 4.2e-4 (L ile hafif
  düşüş), kusursuz örgü 3.9e-21 … 1.7e-20. Sonuç L'den bağımsız.

ARTEFAKT ŞÜPHELERİM (gizlemiyorum):
  • Sıfır tabanı veri kuantizasyonuna 1.3–1.5× yakın → ÜST SINIR.
  • Medyan kullanıldı (I üstel dağılımlı; ortalama ≈ medyan/ln2).
    Tüm oranlar aynı istatistikle alındı.
  • 86'nın "vekilin de altında" hükmü YÖN olarak ayakta (vekil 7.99e-2
    vs gerçek 4.21e-4, 190×), sayısı değil.
  • Çizgi maskesi 86'dan sıkı (q ≤ 200 tam asal-kuvvet listesi): 956
    çizgi-dışı nokta (86: 978). Bu, tabanı hafifçe DÜŞÜRÜR.
"""

import numpy as np
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sympy import primerange

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi

# ---------------- çizgi listesi (tam asal-kuvvet) ----------------
def prime_powers(qmax):
    out = []
    for p in primerange(2, qmax + 1):
        q = p
        while q <= qmax:
            out.append(q)
            q *= p
    return sorted(out)

# ---------------- taper makinesi (101j'den) ----------------
def hann(t):
    return 0.5 * (1 - np.cos(TWO_PI * (t - t[0]) / (t[-1] - t[0])))

def I_omega(t, omegas, taper=False, chunk=256, center=True):
    """I(ω) = |Σ w e^{iω(t−c)}|² / Σ w².  Poisson→1, w=1 iken n|Ĝ|²."""
    omegas = np.asarray(omegas, float)
    w = hann(t) if taper else np.ones_like(t)
    nrm = float((w**2).sum())
    tc = t - 0.5 * (t[0] + t[-1]) if center else t   # faz küçültme
    out = np.empty(len(omegas))
    for s0 in range(0, len(omegas), chunk):
        ob = omegas[s0:s0 + chunk]
        S = (np.exp(1j * np.outer(ob, tc)) * w[None, :]).sum(axis=1)
        out[s0:s0 + chunk] = np.abs(S)**2 / nrm
    return out

def rvm_N(t):
    x = t / TWO_PI
    return x * np.log(x / np.e) + 7 / 8

def smooth_grid(gaps, tmid):
    """86'nın RvM pürüzsüz ızgarası (aynı sayım, aynı başlangıç)."""
    t0 = tmid[0] - gaps[0] / 2
    kk = np.arange(len(tmid)) + 0.5
    ts = t0 + kk * TWO_PI / np.log(t0 / TWO_PI)
    for _ in range(6):
        fdel = rvm_N(ts) - rvm_N(t0) - kk
        ts = ts - fdel / (np.log(ts / TWO_PI) / TWO_PI)
    return ts

# ---------------- veri ----------------
d41 = np.load(HERE / "41_bigT_windows.npz")
K41 = sorted({x.split("_")[1] for x in d41.files}, key=lambda s: int(s[:-1]))
WNDS = [(d41[f"gaps_{k}"], d41[f"amps_{k}"], d41[f"tmid_{k}"]) for k in K41[:6]]
gaps, amps, tmid = WNDS[1]                 # 86-T4'ün penceresi
L = float(np.log(tmid / TWO_PI).mean())
n = len(tmid)
print(f"[102a] 86-T4 penceresi: L={L:.4f}  n={n}  "
      f"span={tmid[-1]-tmid[0]:.4g}  2π/T={TWO_PI/(tmid[-1]-tmid[0]):.3e}",
      flush=True)

rng = np.random.default_rng(1021)

# ---------------- ω taraması (86'nın birebir kurulumu) ----------------
LINES = [np.log(q) for q in prime_powers(200)]
om_scan = np.linspace(0.15, 4.2, 1200)
om_scan = np.array([o for o in om_scan
                    if min(abs(o - l) for l in LINES) > 0.015])
B_ALL = np.ones(len(om_scan), bool)
B_LO = om_scan < 1.0
B_HI = om_scan > 3.0
BANDS = [("hepsi", B_ALL), ("ω<1", B_LO), ("ω>3", B_HI)]
print(f"[102a] çizgi-dışı ω sayısı: {len(om_scan)} / 1200 "
      f"(86: 978/1200; çizgi listesi q≤200 → daha sıkı)", flush=True)

# ---------------- ızgaralar ----------------
ts_sm = smooth_grid(gaps, tmid)
g_sh = rng.permutation(gaps)
t_sh = tmid[0] + np.cumsum(g_sh) - g_sh / 2
t_lat = np.linspace(tmid[0], tmid[-1], n)
t_poi = np.sort(rng.uniform(tmid[0], tmid[-1], n))

zz = np.empty(n + 1)                        # SIFIRLARIN kendisi
zz[0] = tmid[0] - gaps[0] / 2
zz[1:] = zz[0] + np.cumsum(gaps)

GRIDS = [("gerçek ızgara", tmid), ("SIFIRLAR z_n", zz),
         ("karıştırılmış-gap", t_sh),
         ("RvM pürüzsüz", ts_sm), ("kusursuz örgü", t_lat),
         ("Poisson", t_poi)]

# ---------------- G0: sayısal taban kapısı ----------------
print("\n=== G0) SAYISAL TABAN KAPISI ===", flush=True)
om_rnd = rng.uniform(0.15, 4.2, 200)
for isim, taper in [("tapersiz", False), ("Hann", True)]:
    a = I_omega(t_lat, om_rnd, taper=taper)
    b = I_omega(t_lat, om_rnd, taper=taper, center=False)
    print(f"  kusursuz örgü, {isim:>8}: medyan I = {np.median(a):.3e}  "
          f"(merkezleme yokken {np.median(b):.3e})", flush=True)
print("  → merkezlenmiş sürüm kullanılıyor; ölçülen tabanlar bu sayısal")
print("    tabanın üstündeyse anlamlıdır.")

# ---------------- G1 + A1 ----------------
print("\n=== G1/A1) 86-T4 TEKRARI: beş ızgara × iki ayar ===", flush=True)
print(f"  atım gürültüsü referansı: I=1 ⇔ |Ĝ|=1/√n={1/np.sqrt(n):.5f}")
res = {}
for isim, t in GRIDS:
    for taper in [False, True]:
        Iv = I_omega(t, om_scan, taper=taper)
        res[(isim, taper)] = Iv
        print(f"  {isim:>18} {'Hann' if taper else 'ham ':>4}: ", end="")
        for bn, bm in BANDS:
            med = np.median(Iv[bm])
            print(f"{bn}: I={med:.3e} (|Ĝ|={np.sqrt(med/n):.2e})  ", end="")
        print("", flush=True)

print("\n  86 ile doğrudan kıyas (|Ĝ| medyanı, tapersiz):")
print(f"    86 dedi: gerçek 0.0002, karıştırılmış 0.0014, 1/√n = 0.005")
for isim in ["gerçek ızgara", "karıştırılmış-gap", "RvM pürüzsüz",
             "kusursuz örgü"]:
    m = np.median(res[(isim, False)])
    print(f"    {isim:>18}: |Ĝ| = {np.sqrt(m/n):.5f}")

print("\n  ★ ALET KAPISI: gerçek/kusursuz-örgü oranı (I birimi)")
for taper in [False, True]:
    print(f"    {'Hann' if taper else 'tapersiz'}: ", end="")
    for bn, bm in BANDS:
        r = np.median(res[("gerçek ızgara", taper)][bm]) / \
            np.median(res[("kusursuz örgü", taper)][bm])
        print(f"{bn} {r:.2f}×  ", end="")
    print("", flush=True)

# ---------------- A2: duyarlılık merdiveni ----------------
print("\n=== A2) DUYARLILIK MERDİVENİ (iid Gauss titreşim enjeksiyonu) ===",
      flush=True)
print("  öngörü: I_difüz = 1 − e^{−ω²ε²} ≈ ω²ε²  (Einstein/DW limiti)")
print(f"{'ε':>8} {'bant':>6} {'öngörü':>11} {'ölçülen(Hann)':>14} "
      f"{'ölçülen−taban':>14} {'oran':>7}")
TAB = {bn: np.median(res[("gerçek ızgara", True)][bm]) for bn, bm in BANDS}
LAD = []
for eps in [0.0, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2]:
    t_inj = tmid + rng.normal(0, eps, n) if eps > 0 else tmid
    Iv = I_omega(t_inj, om_scan, taper=True)
    row = [eps]
    for bn, bm in BANDS:
        om_m = np.median(om_scan[bm])
        pred = 1 - np.exp(-om_m**2 * eps**2)
        meas = np.median(Iv[bm])
        exc = meas - TAB[bn]
        row += [pred, meas]
        print(f"{eps:>8.1e} {bn:>6} {pred:>11.3e} {meas:>14.3e} "
              f"{exc:>14.3e} {(exc/pred if pred > 0 else np.nan):>7.2f}")
    LAD.append(row)
LAD = np.array(LAD)

# ---------------- A4: rampa kıyası ----------------
print("\n=== A4) RAMPA KIYASI (F=α; taban kaç kat altında?) ===", flush=True)
print(f"{'bant':>6} {'⟨α⟩':>8} {'I_ham':>11} {'ham: α/I':>10} "
      f"{'I_Hann':>11} {'Hann: α/I':>12}")
for bn, bm in BANDS:
    al = float(np.mean(om_scan[bm]) / L)
    ih = np.median(res[("gerçek ızgara", False)][bm])
    it = np.median(res[("gerçek ızgara", True)][bm])
    print(f"{bn:>6} {al:>8.4f} {ih:>11.3e} {al/ih:>10.1f} "
          f"{it:>11.3e} {al/it:>12.3e}")

# ---------------- A5: IZGARA MERDİVENİ ----------------
print("\n=== A5) IZGARA MERDİVENİ: karanlık alan kimin? ===", flush=True)
print("  Ĝ_orta(ω) = (1/n)Σ e^{iωz_n}·e^{iωg_n/2}  → φ(ω)=⟨e^{iωg/2}⟩ (DW)")
gbar = float(gaps.mean())
LAD2 = [("z sıfırlar", zz),
        ("z + ḡ/2 (kontrol)", zz[:-1] + gbar / 2),
        ("z + 0.25·(g−ḡ)/2", zz[:-1] + gbar / 2 + 0.25 * (gaps - gbar) / 2),
        ("z + 0.50·(g−ḡ)/2", zz[:-1] + gbar / 2 + 0.50 * (gaps - gbar) / 2),
        ("z + g/2 = 86 ızgarası", tmid),
        ("z + 1.50·(g−ḡ)/2", zz[:-1] + gbar / 2 + 1.50 * (gaps - gbar) / 2),
        ("z + 2.0·(g−ḡ)/2 ≡ z", zz[:-1] + gbar / 2 + 2.0 * (gaps - gbar) / 2),
        ("z + karıştırılmış g/2", zz[:-1] + rng.permutation(gaps) / 2)]
print(f"{'ızgara':>24} " + " ".join(f"{bn:>13}" for bn, _ in BANDS))
A5 = {}
for isim, tt in LAD2:
    tt = np.sort(tt)
    Iv = I_omega(tt, om_scan, taper=True)
    A5[isim] = Iv
    print(f"{isim:>24} " + " ".join(f"{np.median(Iv[bm]):>13.3e}"
                                    for _, bm in BANDS), flush=True)
print("  BEYAZ ÜST SINIR 1−|φ(ω)|² (gap'ler bağımsız olsaydı):")
for bn, bm in BANDS:
    om_m = float(np.median(om_scan[bm]))
    phi = np.abs(np.mean(np.exp(1j * om_m * gaps / 2)))
    print(f"    {bn:>6} ω≈{om_m:.2f}: |φ|={phi:.4f}  1−|φ|²={1-phi**2:.3e}  "
          f"ölçülen(86 ızgarası)={np.median(A5['z + g/2 = 86 ızgarası'][bm]):.3e}"
          f"  bastırma={(1-phi**2)/np.median(A5['z + g/2 = 86 ızgarası'][bm]):.1f}×")
print("  ÖRNEKLEME-FAZI EĞRİSİ  t(c) = z_n + ḡ/2 + c(g_n−ḡ)/2:")
print("    c=0 → z_n+ḡ/2 (sıfır örgüsü),  c=2 → z_{n+1}−ḡ/2 (YİNE sıfır")
print("    örgüsü!), c=1 → 86'nın orta noktaları. Yani karanlık alan, iki")
print("    ÖZDEŞ karanlık örgü arasında salınan bir ÖRNEKLEME FAZI etkisidir.")
for c, k in [(0.0, "z + ḡ/2 (kontrol)"), (0.25, "z + 0.25·(g−ḡ)/2"),
             (0.5, "z + 0.50·(g−ḡ)/2"), (1.0, "z + g/2 = 86 ızgarası"),
             (1.5, "z + 1.50·(g−ḡ)/2"), (2.0, "z + 2.0·(g−ḡ)/2 ≡ z")]:
    v = np.median(A5[k][B_ALL])
    print(f"    c={c:>4}: I(hepsi)={v:.3e}")

# ---------------- A2b: SIFIR IZGARASINDA DUYARLILIK ----------------
print("\n=== A2b) SIFIR IZGARASINDA DUYARLILIK (taban ne kadar gerçek?) ===",
      flush=True)
print(f"{'ε':>9} {'bant':>6} {'öngörü ω²ε²':>13} {'ölçülen':>12} "
      f"{'ölçülen−taban':>14} {'oran':>7}")
TABZ = {bn: np.median(res[("SIFIRLAR z_n", True)][bm]) for bn, bm in BANDS}
for eps in [1e-7, 3e-7, 1e-6, 3e-6, 1e-5]:
    Iv = I_omega(np.sort(zz + rng.normal(0, eps, len(zz))), om_scan, taper=True)
    for bn, bm in BANDS:
        om_m = float(np.median(om_scan[bm]))
        pred = 1 - np.exp(-om_m**2 * eps**2)
        meas = np.median(Iv[bm]); exc = meas - TABZ[bn]
        print(f"{eps:>9.1e} {bn:>6} {pred:>13.3e} {meas:>12.3e} "
              f"{exc:>14.3e} {exc/pred:>7.2f}")

# ---------------- A2c: VERİ SETİNİN KENDİ KUANTİZASYONU ----------------
print("\n=== A2c) SIFIR-BULUCUNUN KUANTİZASYON TABANI (dürüstlük kapısı) ===",
      flush=True)
print("  41_bigT_scan.py: adım (2π/L)/20 üzerinden 18 bisection iterasyonu")
Lstep = (TWO_PI / L) / 20
half = Lstep / 2**18 / 2
print(f"  → sıfır konumu belirsizliği ±{half:.3e} (σ_kuant={half/np.sqrt(3):.3e})")
Iq = I_omega(np.sort(zz + rng.uniform(-half, half, len(zz))), om_scan, taper=True)
for bn, bm in BANDS:
    print(f"    {bn:>6}: sıfır tabanı {np.median(res[('SIFIRLAR z_n', True)][bm]):.3e}"
          f"   +kuantizasyon {np.median(Iq[bm]):.3e}"
          f"   artış {np.median(Iq[bm])/np.median(res[('SIFIRLAR z_n', True)][bm]):.2f}×")
print("  Artış ~2× ise ölçülen sıfır tabanı VERİ-SINIRLIDIR (üst sınır).")

# ---------------- A3: altı pencere ----------------
print("\n=== A3) ALTI PENCERE (çizgi-dışı medyan I) ===", flush=True)
print(f"{'L':>7} {'n':>7} {'ham gerçek':>12} {'ham örgü':>11} "
      f"{'Hann gerçek':>13} {'Hann örgü':>11} {'oran(Hann)':>11}")
A3 = []
for gg, aa, tt in WNDS:
    Lw = float(np.log(tt / TWO_PI).mean())
    tl = np.linspace(tt[0], tt[-1], len(tt))
    r1 = np.median(I_omega(tt, om_scan, taper=False))
    r2 = np.median(I_omega(tl, om_scan, taper=False))
    r3 = np.median(I_omega(tt, om_scan, taper=True))
    r4 = np.median(I_omega(tl, om_scan, taper=True))
    A3.append((Lw, len(tt), r1, r2, r3, r4))
    print(f"{Lw:>7.3f} {len(tt):>7} {r1:>12.3e} {r2:>11.3e} "
          f"{r3:>13.3e} {r4:>11.3e} {r3/r4:>11.2f}", flush=True)
A3 = np.array(A3)

# ---------------- figür ----------------
fig, ax = plt.subplots(1, 2, figsize=(13.5, 5.0))
for isim, c in [("gerçek ızgara", "firebrick"), ("karıştırılmış-gap", "steelblue"),
                ("kusursuz örgü", "0.55")]:
    ax[0].semilogy(om_scan, np.sqrt(res[(isim, False)] / n), ".", ms=2,
                   alpha=0.55, color=c, label=isim)
ax[0].axhline(1 / np.sqrt(n), color="k", lw=0.9, ls="--", label="1/√n")
ax[0].set_xlabel("ω (çizgi-dışı)"); ax[0].set_ylabel("|Ĝ(ω)|")
ax[0].set_title("TAPERSİZ (86-T4): gerçek ızgara = kusursuz örgü")
ax[0].legend(fontsize=8); ax[0].grid(alpha=0.25)
for isim, c in [("gerçek ızgara", "firebrick"), ("SIFIRLAR z_n", "seagreen"),
                ("karıştırılmış-gap", "steelblue"), ("kusursuz örgü", "0.55")]:
    ax[1].semilogy(om_scan, res[(isim, True)], ".", ms=2, alpha=0.55,
                   color=c, label=isim)
ax[1].axhline(1.0, color="k", lw=0.9, ls="--", label="atım gürültüsü I=1")
ax[1].plot(om_scan, om_scan / L, "k:", lw=1.2, label="Montgomery rampası α")
ax[1].set_yscale("log")
ax[1].set_xlabel("ω (çizgi-dışı)"); ax[1].set_ylabel("I(ω) = n|Ĝ|²")
ax[1].set_title("HANN TAPERLİ: karanlık alan gerçekten karanlık")
ax[1].legend(fontsize=8); ax[1].grid(alpha=0.25)
fig.tight_layout()
fig.savefig(HERE / "102a_karanlik_alan.png", dpi=125)

np.savez(HERE / "102a_karanlik.npz", om=om_scan, L=L, n=n,
         **{f"I_{i}_{int(tp)}": res[(k, tp)]
            for i, (k, _) in enumerate(GRIDS) for tp in [False, True]},
         ladder=LAD, A3=A3)
print("\n102a_karanlik_alan.png + 102a_karanlik.npz yazıldı.", flush=True)
