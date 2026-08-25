"""
107a — ESKİ DÖRT ADANIN İKİ YENİ KUSUR KAPISIYLA YENİDEN DENETİMİ
       (25 Ağustos 2026)
==========================================================================
105'in kapanış cümlesi bir iş bırakmıştı: "χ₃'ün orta-τ çukuru yeni
kapıyla DERİNLEŞTİ → 101f verilerinin (4 eski ada) kısa-çukur kapısıyla
yeniden denetimi SIRADAKİ İŞ." Bu script o işi yapar ve iki kapıyı
BİRLİKTE, tek boru hattında çalıştırır:

  KAPI 1 (101d, eski)  DÜZLÜK/BASAMAK: 80'lik kayan medyanın adımı
                       |Δmed| > 0.5 → ±160 sıfır at.
  KAPI 2 (105e-K4, YENİ) KISA-ÇUKUR: |med₂₀ − med₂₀₀| > 0.7 → ±160 at.
                       (2 kayıp + anında telafi; 80'lik dedektör YIKAR.)
  KAPI 3 (104e, YENİ)  SIÇRAMA: r = ort₂₁(d) − med₈₀₁(d), |r| > 0.30
                       → ±600 at; ek olarak normalize boşluk u < 0.005
                       (kopya sıfır). Pencere-içi uygulanır (104e'nin
                       kalibrasyonu pencere başınaydı).

Ölçüm: 101d estimatörü (yerel katlama + kübik detrend + havuzlanmış
D(τ) + 40-permütasyon vekil taban), τ ızgarası 101d'ninki, BAND ±0.01L
(101h/105d dersi: ±0.02L artefakt taşıyor), band-kaçağı kapısı açık,
ve 104'ün dersi gereği havuz DAİMA PENCERE-AYRIŞTIRMALI raporlanır.

ÖN-MÜHÜR (ölçümden ÖNCE yazıldı):
  A1  POZİTİF KONTROL: kapılar 104d'nin elle bulduğu iki kusuru
      yeniden bulacak — χ₇ w1'de t ≈ 3081.67 ve 3081.89 kopya çiftleri
      (u ≈ 0.0008) ve χ₃ w1'de t ∈ [3872.5, 3890.2] bölgesi. Bulamazsa
      kapı tarifi yanlış demektir.
  A2  β ve χ₅ TEMİZ kalacak: sıçrama kapısı 104e'de ikisinde de 0
      işaret vermişti; kısa-çukur kapısı da ada başına ≲40 sıfırdan
      fazlasını işaretlemeyecek (105e: β'da maks fark 0.37 → 0 kusur).
  A3  RİSKLİ: χ₃'ün orta-τ çukuru (101d τ=0.105-0.140'ta 0.545/0.489/
      0.583) TEMİZLİKLE DÜZELMEYECEK. Çünkü 104 onu VERİ KUSURUNA
      değil, HAVUZLAMA GİRİŞİMİNE bağladı (D = |Σ_w num_w| / Σ_w den_w
      koherent toplamdır). Pencere-ayrıştırmalı raporda tek tek
      pencereler ≳0.4-0.6 iken HAVUZ onların altında kalacak. Çukur
      temizlikle DEĞİL, ancak pencere-ayrıştırmasıyla "açıklanacak".
  A4  Eşik hükmü değişmeyecek: ζ/χ₃/χ₅/χ₇ log2/⟨L⟩ ≈ 0.068-0.074'te
      çözülür (D ≳ 0.3), β log3/⟨L⟩ = 0.114'e dek donuk (≲0.10).

==========================================================================
SONUÇ (25 Ağustos, koşu 81 s) — A1 ✓✓ A2 ✓ A3 ✓/✗ (İKİYE BÖLÜNDÜ) A4 ✓✓
==========================================================================
A1 ✓✓ POZİTİF KONTROL TAM İSABET. Kapılar 104d'nin ELLE bulduğu iki
  kusuru otomatik olarak yeniden buldu:
    χ₃  K2: t ∈ [3872.5365, 3892.7874]  (24 sıfır) — 104d [3872.54, 3890.24]
    χ₇  K2: t ∈ [3081.8885, 3109.5478]  (35 sıfır) — 104d [3081.67, 3107.56]
  Ayrıca 104d'nin "yumuşak kuyruğu" da yakalandı:
    χ₇  K3: t ∈ [3664.5970, 3680.8290]  (22 sıfır) — 104d'nin "i₀ ≥ 1000,
        t ≥ 3680" gözlemiyle birebir. χ₇ w1'in maks|r|'si kuyruk hâlâ
        içerideyken 0.72 (104e'de 2.01'di; K2 sıçramanın gövdesini
        önceden aldığı için düşük).
  ÜÇ YENİ KUSUR (daha önce hiç raporlanmamış), hepsi K2'den:
    χ₅  t ∈ [5810.6964, 5828.4765]  (23 sıfır, sapma 2.01) ← χ₅ "tertemiz"
        sanılıyordu (104e sıçrama kapısı 0 işaret vermişti); kısa-çukur
        kapısı onda da bir kusur buldu.
    χ₇  t ∈ [ 767.6227,  831.4470]  (70 sıfır) ve
    χ₇  t ∈ [1903.7250, 1915.3388]  (16 sıfır) — ikisi de χ₇'nin en alt
        log-penceresinde; o pencere zaten n=2578 < 3000 olduğu için
        analize GİRMİYORDU (ölçümü etkilemez, ama sertifika kaydına
        girer). Not: [767.6, 831.4] bloğu dizinin sol UCUNDA (i ∈ [0,69])
        ve ortalanmış kayan medyanın kenar doldurması orada yanlış-pozitif
        üretebilir; DENETLENDİ — ham d gerçekten i=60→80 arası −2.3
        basamak atıyor, yani GERÇEK kusur (K1 de aynı yeri işaretliyor).
  KOPYA ENVANTERİ (107c): χ₇'de u<0.005 olan ALTI çift var
  (t = 935.4419, 1531.0633, 1913.7446, 3081.6688, 3081.8885, 36263.8528)
  — 104d yalnız ikisini bulmuştu. Hepsi nihai analiz pencerelerinin
  DIŞINDA (kimi n<3000 düşen pencerede, kimi K1/K2 padinde). Bu yüzden
  KAPI 3'ün u-ölçütü hiçbir adada tetiklenmedi: kopyalar ona sıra
  gelmeden zaten temizlenmiş oluyor.

A2 ✓ β TERTEMİZ (K1 dışında hiçbir kapı işaret vermedi; maks
  |med₂₀−med₂₀₀| = 0.37, eşiğin altında — 105e'nin yanlış-pozitif
  denetimi bağımsız olarak tekrarlandı). χ₅ ise ≤40 sınırında ama
  SIFIR DEĞİL (23 işaret): "χ₅ temiz" hükmü ARTIK GEÇERSİZ.

A3 — ÖZEL SORU: χ₃'ün orta-τ çukuru TEMİZLİKLE DÜZELMİYOR, DERİNLEŞİYOR.
  τ       0.095  0.105  0.113  0.125  0.140
  ÖNCE    0.590  0.565  0.512  0.320  0.766     (yalnız KAPI 1)
  SONRA   0.426  0.199  0.161  0.134  0.531     (üç kapı)
  Yalnız 816 sıfır (%1.1) maskelendi ve D(0.105) 0.565 → 0.199'a düştü —
  101b'nin "estimatör kusura aşırı duyarlı" uyarısının ters yönde
  doğrulanışı. (Bu satır 105d'nin χ₃ satırını ONDALIKLARINA KADAR
  yeniden üretiyor: 0.036/0.124/0.194/0.597/0.617/0.426/0.199/0.161 —
  boru hattı bağımsız olarak yeniden kuruldu, aynı sayılar çıktı.
  β satırı da birebir: 0.012/0.026/0.037/0.053/0.091/0.178/0.647/0.726.)
  A3'ÜN İKİNCİ YARISI ÇÜRÜDÜ (dürüst kayıt): çukuru HAVUZLAMA GİRİŞİMİNE
  bağlamayı öngörmüştüm; PENCERE AYRIŞTIRMASI bunu ÇÜRÜTTÜ.
    τ=0.105: pencereler 0.282 / 0.134 / 0.137, havuz 0.199 → havuz/maks
             0.71; havuz pencerelerin ORTASINDA, ALTINDA değil.
    τ=0.113: 0.215 / 0.110 / 0.156, havuz 0.161 (oran 0.75)
    τ=0.125: 0.174 / 0.088 / 0.143, havuz 0.134 (oran 0.77)
  ÜÇ PENCERENİN ÜÇÜ DE ÇUKURDA. Yani χ₃'ün orta-τ çukuru ALETSEL DEĞİL
  (ne kusur ne havuzlama) — n=63549'luk BÜYÜK pencere de dâhil, χ₃'ün
  kendi eğrisinin bir özelliği. Buna karşılık τ=0.068'de yıkıcı girişim
  GERÇEKTEN var: pencereler 0.320/0.145/0.544, havuz 0.194 (oran 0.36) —
  104'ün teşhisi O τ NOKTASINDA doğru, orta-τ çukurunda değil.
  Dört ada + ζ'nın havuz/maks oranları τ=0.105-0.140'ta 0.59-0.92 (χ₃
  0.67-0.77 ile bandın içinde) → orta-τ'da genel bir havuzlama patolojisi
  YOK. AÇIK SORU olarak kayda geçiyor: χ₃'ün (q=3, en küçük iletken,
  en küçük ⟨L⟩ = 9.44) orta-τ çukuru fiziksel mi, yoksa hâlâ görülmemiş
  bir üçüncü kusur sınıfı mı?

A4 ✓✓ ANA HÜKÜM SAĞLAM (band ±0.01L, temizlenmiş):
  τ=0.068'de ζ 0.489, χ₃ 0.194(*), χ₅ 0.442, χ₇ 0.507 — β 0.037
  (vekil taban 0.004). β log3/⟨L⟩ = 0.1143'e dek donuk (0.012-0.091)
  ve tam orada 0.726'ya sıçrıyor. τ=0.30'da beşi 1.02-1.13 (tek CUE).
  (*) χ₃'ün 0.194'ü yukarıdaki τ=0.068 yıkıcı girişiminin sonucu;
  pencere-başına maks 0.544 ve D(τ_ilk = 0.0735) = 0.508.

YAN KAZANÇ — χ₇'NİN τ=0.04 ANOMALİSİ BAĞIMSIZ OLARAK KAPANDI:
  ÖNCE 0.106 → SONRA 0.037 (104e'de 0.125 → 0.038, farklı bant ve
  farklı kapı sırasıyla). χ₇ ayrıca τ=0.085'te 0.488 → 0.334.
==========================================================================
"""
import numpy as np, time, json
from pathlib import Path

T0 = time.time()
HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
rng = np.random.default_rng(1070)

PKS = [2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 17, 19, 23, 25, 27, 29, 31, 32, 37,
       41, 43, 47, 49, 53, 59, 61, 64, 67, 71, 73, 79, 81, 83, 89, 97, 101,
       103, 107, 109, 113, 121, 125, 127, 128]
LINES_ALL = [np.log(qq) for qq in PKS]

TAUS = [0.04, 0.05, 0.06, 0.068, 0.075, 0.085, 0.095, 0.105, 0.113,
        0.125, 0.14, 0.16, 0.20, 0.30]
BAND = 0.01


# ----------------------------------------------------------------- alet
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
    """101d estimatörü. DÖNÜŞ: (havuz D, vekil V, pencere-başına D).

    104'ün PRATİK KURALI (4): havuzlanmış D daima pencere-başına
    ayrıştırılarak raporlanır — D = |Σ_w num_w| / (c Σ_w den_w)
    KOHERENT bir toplamdır ve pencereler ters fazda gelirse havuz
    her bir pencerenin ALTINA düşebilir (yıkıcı girişim).
    """
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
    """Ortalanmış kayan istatistik (kenarlar uçtaki değerle doldurulur)."""
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
    """Maskedeki bitişik True bloklarını (i0, i1, t0, t1) olarak döndür."""
    out, i = [], 0
    n = len(mask)
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
    """KAPI 1 (basamak, 101d) + KAPI 2 (kısa-çukur, 105e-K4).

    İkisi de HAM diziye uygulanır; işaretli bölgeler ±pad ile atılır.
    DÖNÜŞ: bad maskesi, kapı-1 ham işaretleri, kapı-2 ham işaretleri.
    """
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
    return g1, g2, ham2, float(sapma.max()), d


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
    """KAPI 3 — SIÇRAMA (104e). Pencere-İÇİ uygulanır (kalibrasyon öyle).

    r = 21'lik kayan ORTALAMA − 801'lik kayan MEDYAN. 104e: temiz 17
    pencerede maks|r| = 0.06-0.12; kusurlu ikisinde 1.99/2.01. Eşik
    0.30, pad 600 (yumuşak kuyruk ~1000 sıfır sürüyor).
    Ek ölçüt: normalize boşluk u < 0.005 → kopya sıfır.
    """
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


# ------------------------------------------------------------ veri hazır
exec(open(HERE / "98_L_motoru.py").read().split('CHI4 = ')[0])
HERE = Path(__file__).resolve().parent

CHI4 = {0: 0, 1: 1, 2: 0, 3: -1}
CHI3 = {0: 0, 1: 1, 2: -1}
CHI5 = {0: 0, 1: 1, 2: 1j, 3: -1j, 4: -1}
z6c = np.exp(1j * np.pi / 3)
CHI7 = {0: 0, 1: 1, 3: z6c, 2: z6c**2, 6: z6c**3, 4: z6c**4, 5: z6c**5}

ADA = [("chi3", 3, CHI3, 1, 2), ("beta", 4, CHI4, 1, 3),
       ("chi5", 5, CHI5, 1, 2), ("chi7", 7, CHI7, 1, 2)]


def rvm_sayim(t):
    x = t / TWO_PI
    return x * np.log(x / np.e)


ZZ, SAY, QQ, QILK, NHAM, BOLG = {}, {}, {}, {}, {}, {}

# --- ζ referansı: 41'in ilk 6 penceresi, tek diziye eklenmez (ayrı ayrı)
d41 = np.load(HERE / "41_bigT_windows.npz")
K41 = sorted({x.split("_")[1] for x in d41.files}, key=lambda s: int(s[:-1]))
ZETA_PARCA = []
for k in K41[:6]:
    gz, tm = d41[f"gaps_{k}"], d41[f"tmid_{k}"]
    zz = np.empty(len(gz) + 1)
    zz[0] = tm[0] - gz[0] / 2
    zz[1:] = zz[0] + np.cumsum(gz)
    ZETA_PARCA.append(zz)
SAY["zeta"] = rvm_sayim
QQ["zeta"] = 1.0
QILK["zeta"] = 2
NHAM["zeta"] = sum(len(z) for z in ZETA_PARCA)
BOLG["zeta"] = 0

for et, q, tab, a, qilk in ADA:
    d = np.load(HERE / f"101f_{et}_zeros.npz")
    zc = d["zeros"]
    lo = np.exp(np.log(zc[0] + 1) + 0.25 * (np.log(zc[-1]) - np.log(zc[0] + 1)))
    zc = zc[zc >= lo]
    M = Lmotor(q, tab, a)
    ZZ[et] = zc
    SAY[et] = (lambda t, M=M: M.theta(t) / np.pi)
    QQ[et] = float(q)
    QILK[et] = qilk
    NHAM[et] = len(zc)
    BOLG[et] = int(np.asarray(d["bolgeler"]).reshape(-1, 2).shape[0])

ADLAR = ["zeta", "chi3", "beta", "chi5", "chi7"]


# ------------------------------------------------- BORU HATTI: ÖNCE/SONRA
def boru(ad, yeni_kapilar):
    """yeni_kapilar=False → 101d'nin hâli (yalnız KAPI 1).
       yeni_kapilar=True  → KAPI 1 + KAPI 2 + KAPI 3."""
    tan = {"g1": 0, "g2": 0, "g2_ham": 0, "g3": 0, "dup": 0,
           "sapma_max": 0.0, "rmax": [], "kusur": [], "kopya": []}
    if ad == "zeta":
        kaynak = ZETA_PARCA
    else:
        kaynak = [ZZ[ad]]
    WIN = []
    for zc in kaynak:
        g1, g2, ham2, smax, _ = kapi12(zc, SAY[ad])
        tan["sapma_max"] = max(tan["sapma_max"], smax)
        tan["g1"] += int(g1.sum())
        tan["g2"] += int(g2.sum())
        tan["g2_ham"] += int(ham2.sum())
        if yeni_kapilar:
            for (i0, i1, t0, t1) in bloklar(ham2, zc):
                tan["kusur"].append(("K2", ad, t0, t1, i1 - i0 + 1))
        bad = g1 | (g2 if yeni_kapilar else np.zeros(len(zc), bool))
        seg, _ = segmentle(zc, bad)
        if ad == "zeta":
            parcalar = seg
        else:
            kenar = np.exp(np.linspace(np.log(zc[0]),
                                       np.log(zc[-1] * 1.0001), 4))
            parcalar = []
            for s in seg:
                for i in range(3):
                    pw = s[(s >= kenar[i]) & (s < kenar[i + 1])]
                    if len(pw) >= 3000:
                        parcalar.append(pw)
        for p in parcalar:
            if not yeni_kapilar:
                WIN.append(pencere_hazirla(p, QQ[ad]))
                continue
            alt, bad3, kes3, rmax, u, dup = kapi3(p, QQ[ad], SAY[ad])
            tan["rmax"].append(rmax)
            tan["g3"] += int(kes3.sum())
            tan["dup"] += len(dup)
            for (i0, i1, t0, t1) in bloklar(bad3, p):
                tan["kusur"].append(("K3", ad, t0, t1, i1 - i0 + 1))
            for i in dup:
                tan["kopya"].append((ad, float(p[i]), float(p[i + 1]),
                                     float(p[i + 1] - p[i]), float(u[i])))
            for x in alt:
                WIN.append(pencere_hazirla(x, QQ[ad]))
    return WIN, tan


print("=" * 74)
print("107a — ESKİ DÖRT ADA + ζ : İKİ YENİ KAPIYLA YENİDEN DENETİM")
print("=" * 74, flush=True)

PEN_A, PEN_B, TAN = {}, {}, {}
for a in ADLAR:
    PEN_A[a], _ = boru(a, False)
    PEN_B[a], TAN[a] = boru(a, True)
    nA = sum(len(w[0]) for w in PEN_A[a])
    nB = sum(len(w[0]) for w in PEN_B[a])
    t = TAN[a]
    print(f"[{a}] ham n={NHAM[a]}  |  ÖNCE {len(PEN_A[a])} pencere n={nA}"
          f"  →  SONRA {len(PEN_B[a])} pencere n={nB}"
          f"   (maskelenen {nA - nB})")
    print(f"      KAPI1 basamak {t['g1']:>6} sıfır | "
          f"KAPI2 kısa-çukur ham {t['g2_ham']:>4} → padli {t['g2']:>6} "
          f"(maks |med20−med200| = {t['sapma_max']:.2f}) | "
          f"KAPI3 sıçrama padli {t['g3']:>6}, kopya {t['dup']}")
    print(f"      KAPI3 pencere-başına maks|r|: "
          + " ".join(f"{x:.2f}" for x in t["rmax"]), flush=True)

# ------------------------------------------------ A1: POZİTİF KONTROL
print("\n" + "=" * 74)
print("A1 — POZİTİF KONTROL: 104d'NİN İKİ KUSURU YENİDEN BULUNDU MU?")
print("=" * 74)
print("  (104d elle bulmuştu: χ₇ t≈3081.6688/3081.6695 ve 3081.8885/"
      "3081.8891\n   kopya çiftleri; χ₃ t∈[3872.54, 3890.24] bölgesi)")
print("\n  --- İŞARETLENEN KUSUR BÖLGELERİ (kapı, ada, t-aralığı, sıfır) ---")
for a in ADLAR:
    ks = TAN[a]["kusur"]
    if not ks:
        print(f"  [{a}] TEMİZ — hiçbir kapı işaret vermedi")
        continue
    for (kap, _, t0, t1, n) in ks:
        print(f"  [{a}] {kap}: t ∈ [{t0:.4f}, {t1:.4f}]  ({n} sıfır, "
              f"genişlik {t1 - t0:.3f})")
print("\n  --- KOPYA SIFIRLAR (u < 0.005) ---")
kop_var = False
for a in ADLAR:
    for (_, t1, t2, dt, uu) in TAN[a]["kopya"]:
        kop_var = True
        print(f"  [{a}] t = {t1:.4f} ve {t2:.4f}   Δt = {dt:.2e}   u = {uu:.5f}")
if not kop_var:
    print("  (hiç kopya yok)")

# ------------------------------------------------------- EĞRİLER
print("\n" + "=" * 74)
print(f"D(τ) — ÖNCE (yalnız KAPI 1) ve SONRA (üç kapı), band ±{BAND}L")
print("=" * 74, flush=True)
SON = {}
for a in ADLAR:
    DA, VA, WA = egri(PEN_A[a], TAUS)
    DB, VB, WB = egri(PEN_B[a], TAUS)
    SON[a] = dict(DA=DA, VA=VA, WA=WA, DB=DB, VB=VB, WB=WB)
    print(f"  ...{a} bitti ({time.time()-T0:.0f} s)", flush=True)

print(f"\n{'τ':>6} " + "  ".join(f"{a:>15}" for a in ADLAR) + "   vekil(mx)")
print(f"{'':>6} " + "  ".join(f"{'önce':>7}{'sonra':>8}" for a in ADLAR))
for i, t in enumerate(TAUS):
    vek = max(SON[a]["VB"][i] for a in ADLAR)
    hu = []
    for a in ADLAR:
        hu.append(f"{SON[a]['DA'][i]:>7.3f}{SON[a]['DB'][i]:>8.3f}")
    print(f"{t:>6.3f} " + "  ".join(hu) + f"   {vek:>7.3f}")

print("\nvekil taban SONRA (aile aile):")
for a in ADLAR:
    print(f"  {a:>6}: " + " ".join(f"{v:.3f}" for v in SON[a]["VB"]))

# ------------------------------------ ⟨L⟩, τ_ilk, band kaçağı
LORT, TILK, SIZ = {}, {}, {}
for a in ADLAR:
    Ls = np.array([w[3] for w in PEN_B[a]])
    ns = np.array([len(w[0]) for w in PEN_B[a]], float)
    LORT[a] = float(np.average(Ls, weights=ns))
    TILK[a] = np.log(QILK[a]) / LORT[a]
    SIZ[a] = [int(np.sum(Ls * (t + BAND) >= np.log(QILK[a]))) for t in TAUS]

print(f"\n{'aile':>6} {'⟨L⟩_n':>8} {'q_ilk':>5} {'τ_ilk':>8}")
for a in ADLAR:
    print(f"{a:>6} {LORT[a]:>8.3f} {QILK[a]:>5} {TILK[a]:>8.5f}")

print(f"\n=== BAND KAÇAĞI KAPISI (sızan pencere sayısı, ±{BAND}L) ===")
print(f"{'aile':>6} {'n_pen':>5} " + " ".join(f"{t:>6.3f}" for t in TAUS))
for a in ADLAR:
    print(f"{a:>6} {len(PEN_B[a]):>5} " + " ".join(f"{x:>6d}" for x in SIZ[a]))
print("\n=== KAÇAKSIZ TABLO — SONRA (sızan τ'lar ★) ===")
print(f"{'τ':>6} " + " ".join(f"{a:>9}" for a in ADLAR))
for i, t in enumerate(TAUS):
    hu = [f"{SON[a]['DB'][i]:>8.3f}" + ("★" if SIZ[a][i] > 0 else " ")
          for a in ADLAR]
    print(f"{t:>6.3f} " + " ".join(hu))

# ------------------------------- A3: χ₃ ORTA-τ ÇUKURU, PENCERE AYRIŞTIRMASI
print("\n" + "=" * 74)
print("A3 — ÖZEL SORU: χ₃'ÜN ORTA-τ ÇUKURU TEMİZLİKLE DÜZELİYOR MU?")
print("=" * 74)
print("  (101d ±0.02L: τ=0.095/0.105/0.113/0.125 → 0.605/0.545/0.489/0.583;"
      "\n   105d ±0.01L kapılı: 0.426/0.199/0.161/0.134)")
print(f"\n  χ₃ HAVUZ vs PENCERE-BAŞINA (SONRA, band ±{BAND}L, "
      f"{len(PEN_B['chi3'])} pencere)")
print(f"  L'ler: " + " ".join(f"{w[3]:.2f}" for w in PEN_B["chi3"])
      + "   n: " + " ".join(str(len(w[0])) for w in PEN_B["chi3"]))
print(f"{'τ':>6} {'ÖNCE hav':>9} {'SONRA hav':>10} {'pencere-başına (SONRA)':>28}"
      f" {'maks_p':>7} {'hav/maks':>9}")
for i, t in enumerate(TAUS):
    w = SON["chi3"]["WB"][i]
    mx = max(w)
    print(f"{t:>6.3f} {SON['chi3']['DA'][i]:>9.3f} "
          f"{SON['chi3']['DB'][i]:>10.3f}  "
          + " ".join(f"{x:>6.3f}" for x in w)
          + f"   {mx:>6.3f} {SON['chi3']['DB'][i]/max(mx,1e-9):>8.2f}")

print("\n  --- AYNI AYRIŞTIRMA, DÖRT ADA + ζ (yıkıcı girişim denetimi) ---")
print(f"{'aile':>6} {'τ':>6} {'havuz':>7} {'maks_pencere':>13} {'oran':>6}")
for a in ADLAR:
    for i, t in enumerate(TAUS):
        if t not in (0.068, 0.105, 0.113, 0.125, 0.14):
            continue
        w = SON[a]["WB"][i]
        mx = max(w)
        print(f"{a:>6} {t:>6.3f} {SON[a]['DB'][i]:>7.3f} {mx:>13.3f} "
              f"{SON[a]['DB'][i]/max(mx,1e-9):>6.2f}")

# ------------------------------------------------------------- kayıt
np.savez(HERE / "107a_temiz.npz",
         taus=np.array(TAUS), band=BAND, adlar=np.array(ADLAR),
         Lort=np.array([LORT[a] for a in ADLAR]),
         tilk=np.array([TILK[a] for a in ADLAR]),
         qilk=np.array([QILK[a] for a in ADLAR]),
         nham=np.array([NHAM[a] for a in ADLAR]),
         nA=np.array([sum(len(w[0]) for w in PEN_A[a]) for a in ADLAR]),
         nB=np.array([sum(len(w[0]) for w in PEN_B[a]) for a in ADLAR]),
         npenA=np.array([len(PEN_A[a]) for a in ADLAR]),
         npenB=np.array([len(PEN_B[a]) for a in ADLAR]),
         bolgeler=np.array([BOLG[a] for a in ADLAR]),
         **{f"DA_{a}": np.array(SON[a]["DA"]) for a in ADLAR},
         **{f"DB_{a}": np.array(SON[a]["DB"]) for a in ADLAR},
         **{f"VA_{a}": np.array(SON[a]["VA"]) for a in ADLAR},
         **{f"VB_{a}": np.array(SON[a]["VB"]) for a in ADLAR},
         **{f"SIZ_{a}": np.array(SIZ[a]) for a in ADLAR},
         **{f"L_{a}": np.array([w[3] for w in PEN_B[a]]) for a in ADLAR},
         **{f"n_{a}": np.array([len(w[0]) for w in PEN_B[a]]) for a in ADLAR},
         **{f"g1_{a}": TAN[a]["g1"] for a in ADLAR},
         **{f"g2_{a}": TAN[a]["g2"] for a in ADLAR},
         **{f"g2ham_{a}": TAN[a]["g2_ham"] for a in ADLAR},
         **{f"g3_{a}": TAN[a]["g3"] for a in ADLAR},
         **{f"dup_{a}": TAN[a]["dup"] for a in ADLAR},
         **{f"smax_{a}": TAN[a]["sapma_max"] for a in ADLAR},
         **{f"rmax_{a}": np.array(TAN[a]["rmax"]) for a in ADLAR})
with open(HERE / "107a_kusurlar.json", "w") as f:
    json.dump({a: {"kusur": TAN[a]["kusur"], "kopya": TAN[a]["kopya"],
                   "g1": TAN[a]["g1"], "g2": TAN[a]["g2"],
                   "g2_ham": TAN[a]["g2_ham"], "g3": TAN[a]["g3"],
                   "dup": TAN[a]["dup"], "sapma_max": TAN[a]["sapma_max"],
                   "rmax": TAN[a]["rmax"]} for a in ADLAR}, f,
              ensure_ascii=False, indent=1)
print(f"\n107a_temiz.npz + 107a_kusurlar.json yazıldı. "
      f"Süre {time.time()-T0:.0f} s.", flush=True)
