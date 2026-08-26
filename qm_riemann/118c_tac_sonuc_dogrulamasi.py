"""
118c — TAÇ SONUÇ DOĞRULAMASI: TAPER'LI VERİDE YEDİ-ADA DONMASI
       (26 Ağustos 2026)
==========================================================================
ÖN-MÜHÜR (ölçümden ÖNCE yazıldı)
==========================================================================
118b yedi GL(1) adasını taper'lı motorla YENİDEN DÖKTÜ. 117'nin dersi
gereği yeni liste eskisinin "düzeltilmiş kopyası" değil BAĞIMSIZ bir
dökümdür (sıfır konumları ~1e-2⟨g⟩ oynar). Bu script 107b'nin TAÇ
SONUÇLARINI yeni veriyle yeniden ölçer.

BORU HATTI 107a/107b ile BİREBİR AYNI (tek değişen: sıfır kaynağı):
  K1 düzlük/basamak + K2 kısa-çukur (ham dizide, ±160 pad)
  → düzlük segmentasyonu (n ≥ 3000) → log-üçleme
  → K3 sıçrama + kopya (pencere içi, ±600 pad)
  → 101d estimatörü, τ ızgarası 101d'ninki, BAND ±0.01L,
    band-kaçağı kapısı açık, havuz DAİMA pencere-ayrıştırmalı (104 dersi).
ζ DEĞİŞMEZ (41_bigT_windows, RVM üretimi — L-motoruyla ilgisi yok);
referans olarak aynı boru hattından geçirilir.

DÖRT SORU ve ÖN-MÜHÜR:
  C1 ★ (AYAKTA KALMALI) LOG2/LOG3 KERVAN AYRIMI. 107b: D(log2/L) →
     LOG2 ζ .461 χ₃ .508 χ₅ .482 χ₇ .512 χ₅ᵉ .492 | LOG3 β .046
     χ₈ᵉ .009 χ₈ᵒ .048 (10-50× kontrast). Veri kalitesi arttığı için
     kontrast KORUNMALI, hatta ARTMALI. Çökerse taç sonuç aletselmiş
     demektir — o zaman da dürüstçe yazılır.
  C2 ★ (AYAKTA KALMALI) ω_yarı = τ_yarı·⟨L⟩ log2 kervanının dört
     Dirichlet adasında log 2 = 0.69315'e oturuyordu (%0.01). Bu
     rakam ARİTMETİK bir sabite kilit olduğu için veri kalitesine
     duyarsız olmalı; %1'den iyi kalmalı.
  C3 ? (AÇIK) LOG3 kervanının %12 ERKEN açılması (ω_yarı 0.9671 vs
     log 3 = 1.0986; yarı/τ_ilk 0.870-0.891, üç adada da sistematik).
     Bu sistematik kayma FİZİKSEL mi, yoksa keskin motorun kaçırdığı
     sıfırların kalıntısı mı? Temiz veride ölçülecek. Kayma KALIRSA
     fiziktir; KAPANIRSA 107'nin açık sorusu bir ALET KUSURU olarak
     kapanır. İki sonuç da bilgi.
  C4 ? (AÇIK) χ₃'ün orta-τ ÇUKURU (107a: τ=0.105'te havuz 0.199,
     üç pencerenin ÜÇÜ DE çukurda: 0.282/0.134/0.137). 107 bunu
     "ne kusur ne havuzlama; ya fiziksel ya da GÖRÜLMEMİŞ BİR ÜÇÜNCÜ
     KUSUR SINIFI" diye bıraktı. Kaldırılmış dip TAM DA böyle bir
     üçüncü sınıftır (sertifikanın basamak dedektörü onu ancak net
     basamak bıraktığında görür). ⟹ Bu seferin en olası sürprizi:
     çukur temiz veride DÜZELİR. Düzelmezse fiziksel olma ihtimali
     ciddi biçimde artar.
  EK: x = τ/τ_ilk yeniden ölçekleme çökmesi (107b: std/ort 0.5125 →
     0.2559) de yeniden ölçülür.

DÜRÜST KAYIT (önceden): D(τ) estimatörü kusura AŞIRI duyarlıdır
(101b: %0.1 silme 0.09→0.55). Bu yüzden sayıların ondalığına kadar
aynı çıkması BEKLENMİYOR; beklenen HÜKÜMLERİN aynı kalmasıdır.
==========================================================================
SONUÇ (26 Ağustos, koşu 96 s) — C1 ✓✓ C2 ✓/DÜZELTME C3 ✓ (AÇIK KALDI,
       ARTIK ALETSEL DEĞİL) C4 ✓ (AÇIK KALDI, ÜÇÜNCÜ KUSUR SINIFI ELENDİ)
==========================================================================
VERİ DİSİPLİNİ: yedi adanın YEDİSİNDE de K1 = K2 = 0, K3 yalnız χ₅'te
(bir kopya çifti), kalan şüpheli bölge 0. maks |med₂₀−med₂₀₀| 0.20-0.29
(eşik 0.7). ζ değişmedi (referans).

C1 ✓✓ KERVAN AYRIMI AYAKTA — ve ölçüt bazında GÜÇLENDİ.
  D(log2/L), ORTAK çizgide (eşiksiz ayrıştırıcı):
    LOG2  ζ .461  χ₃ .508→.533  χ₅ .482→.516  χ₇ .512→.512  χ₅ᵉ .492→.492
    LOG3  β .046→.032  χ₈ᵉ .009→.025  χ₈ᵒ .048→.061
  Beş ada çözülmüş (0.461-0.533), üç ada donuk (0.025-0.061), TEK bir
  τ'da. Kervan ortalamaları: LOG2 .4910→.5028, LOG3 .0343→.0393 ⟹ oran
  14.3× → 12.8× (en kötü-hâl oranı 9.6× → 7.5×). Kontrast rakamı biraz
  düştü çünkü LOG3'ün taban gürültüsü de yükseldi (vekil taban o
  τ'larda .005-.012); AYRIM tartışmasız duruyor.
  D(τ_ilk) — her ada KENDİ ilk çizgisinde: LOG3 kervanı .742/.800/.722
  → .756/.757/.722. Üç donuk ada kendi log 3'ünde ÇÖZÜLÜYOR ve dağılım
  DARALDI. "Donma sınırı adanın İLK SAĞ KALAN ÇİZGİSİDİR" hükmü
  temiz veride 5'e 3 ile aynen duruyor.
  YENİDEN ÖLÇEKLEME ÇÖKMESİ İYİLEŞTİ: x = τ/τ_ilk penceresinde
  std/ort 0.2559 → 0.2184 (ham τ 0.5125 → 0.5089 sabit). Yani veri
  temizlendikçe iki kervan tek eşiğe DAHA SIKI kilitleniyor — bu,
  çökmenin gerçek olduğunun bağımsız delili.

C2 ✓/DÜZELTME — ω_yarı ≈ log 2 SAĞ, ama "%0.01" BİR TESADÜFMÜŞ.
    ada     ESKİ      YENİ
    χ₃    0.6918    0.6857
    χ₅    0.7054    0.6834
    χ₇    0.6792    0.6840
    χ₅ᵉ   0.6966    0.6966
    ort   0.6932    0.6874     log 2 = 0.69315
    sapma  %0.01     %0.83
  DÜRÜST OKUMA: dört adanın SAÇILMASI YARIYA İNDİ (örneklem std
  0.0110 → 0.0062; sem 0.0055 → 0.0031). Eski ortalamanın log 2'ye
  %0.01 oturması bu ±%1.6'lık saçılmanın içinde bir RASTLANTIYDI.
  Yeni ortalama log 2'nin %0.83 altında ve log 2 hâlâ 1.9 sem
  mesafede — yani UYUMLU. Doğru ifade: "ω_yarı log 2'ye ~%1
  doğrulukla oturur", "%0.01" DEĞİL. Not 5'te bu şekilde yazılmalı.

C3 ✓ LOG3 KERVANININ ERKEN AÇILMASI GERÇEK — ALETSEL DEĞİL.
    β 0.9790→0.9783, χ₈ᵉ 0.9668→0.9982, χ₈ᵒ 0.9555→0.9561
    ort 0.9671 → 0.9775;  log 3 = 1.0986
    ERKENLİK %12.0 → %11.0  (log 3'ten 9.9 sem uzakta — anlamlı)
    yarı/τ_ilk: 0.891/0.880/0.870 → 0.890/0.909/0.870 (üç adada da
    sistematik ~0.89, hiç bozulmadı)
    KERVAN KAYMASI ×1.395 → ×1.422 (saf öngörü log3/log2 = 1.585)
  ⟹ "keskin motorun kaçırdığı sıfırların kalıntısı" AÇIKLAMASI ÖLDÜ.
  %11'lik erken açılma FİZİKSELDİR ve mekanizması AÇIK SORUDUR.
  (Kervan kayması 1.585 yerine ~1.42 çıkması da aynı olgunun başka
  yüzü ve aynı ölçüde sağlam.)

C4 ✓ χ₃'ÜN ORTA-τ ÇUKURU TEMİZ VERİDE DURUYOR — "ÜÇÜNCÜ KUSUR
    SINIFI" HİPOTEZİ ELENDİ.
    τ      0.095  0.105  0.113  0.125  0.140
    ESKİ   0.426  0.199  0.161  0.134  0.531   (107a, üç kapı)
    YENİ   0.506  0.237  0.210  0.191  0.634   (taper, sıfır ihlal)
    kıyas: χ₅ 0.517/0.734/0.754/0.741, χ₇ 0.036*/0.716/0.827/0.767
  Çukur biraz sığlaştı (+0.04…+0.06) ama komşu adaların 3-4 katı
  altında KALDI. 107'nin bıraktığı üç şıktan biri artık KAPALI:
  bu bir veri kusuru DEĞİL (yeni veride hiçbir kapı tetiklenmiyor,
  kalan bölge 0, kaldırılmış dip nüfusu 60 kat azaldı) ve bir
  havuzlama patolojisi de DEĞİL.
  PENCERE + BAND AYRIŞTIRMASI (C5) çukuru YERİNE OTURTTU:
    τ=0.105, χ₃  band ±0.005L: havuz 0.217 [0.31 / 0.16]
                  band ±0.010L: havuz 0.237 [0.33 / 0.13]
                  band ±0.020L: havuz 0.412 [0.58 / 0.16]
  İki pencerede de düşük, üç bandın üçünde de düşük, ve BÜYÜK
  pencere (L=9.65, n=63549, t ∈ [13554,55000]) her koşulda 0.11-0.20.
  ⟹ ÇUKUR χ₃'ÜN YÜKSEK-t YARISINDA YAŞIYOR ve band'a sağlam.
  Geriye FİZİKSEL şık kalıyor. AÇIK SORU olarak devam, ama artık
  "alet" savunması tükendi.

(*) ARTEFAKT AVI — TEK YENİ ANOMALİ, GİZLENMEDİ:
  χ₇'nin D(0.095)'i 0.527 → 0.036'ya düştü (tek nokta çukuru).
  C5 teşhis etti: band ±0.01L'de havuz 0.036 iken PENCERELER
  0.42 / 0.42 — havuz İKİ PENCERENİN DE ÇOK ALTINDA. Bu 104'ün
  YIKICI GİRİŞİM patolojisinin ders kitabı örneği (D koherent bir
  toplamdır). Band değişince yok oluyor: ±0.005L'de 0.361
  [0.40/0.52], ±0.02L'de 0.602 [0.59/0.61]. ⟹ ALETSEL, fiziksel
  değil. Taç sonuçları ETKİLEMEZ: χ₇'nin τ_ilk'i 0.0693, yarı-τ'sı
  0.0684 — ikisi de bu noktadan uzak.
  İKİNCİ (küçük) kayıt: χ₈ᵉ'nin D(0.140) 0.638→0.474 ve D(0.200)
  0.953→0.778 oynadı; pencere oranı 0.90 (girişim yok), estimatörün
  bilinen duyarlılık bandı içinde.
==========================================================================
"""
import numpy as np, time
from pathlib import Path
from numpy.lib.stride_tricks import sliding_window_view

T0 = time.time()
HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
rng = np.random.default_rng(1070)          # 107a ile AYNI tohum

exec(open(HERE / "118a_taperli_L_motoru.py").read().split(
    'def kesin_Z')[0].split('if __name__')[0])
HERE = Path(__file__).resolve().parent

PKS = [2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 17, 19, 23, 25, 27, 29, 31, 32, 37,
       41, 43, 47, 49, 53, 59, 61, 64, 67, 71, 73, 79, 81, 83, 89, 97, 101,
       103, 107, 109, 113, 121, 125, 127, 128]
LINES_ALL = [np.log(qq) for qq in PKS]
TAUS = [0.04, 0.05, 0.06, 0.068, 0.075, 0.085, 0.095, 0.105, 0.113,
        0.125, 0.14, 0.16, 0.20, 0.30]
BAND = 0.01


# --------------------------------------------------------------- alet
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
    """101d estimatörü. DÖNÜŞ: (havuz D, vekil V, pencere-başına D)."""
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


def rvm_sayim(t):
    x = t / TWO_PI
    return x * np.log(x / np.e)


# ------------------------------------------------------------ veri hazır
KERVAN = {"zeta": 2, "chi3": 2, "chi5": 2, "chi7": 2, "chi5e": 2,
          "beta": 3, "chi8e": 3, "chi8o": 3}
QLAB = {"zeta": "—", "chi3": "3 (a=1)", "chi5": "5 (a=1)", "chi7": "7 (a=1)",
        "chi5e": "5 (a=0)", "beta": "4 (a=1)", "chi8e": "8 (a=0)",
        "chi8o": "8 (a=1)"}
ADLAR = ["zeta", "chi3", "chi5", "chi7", "chi5e", "beta", "chi8e", "chi8o"]

print("=" * 78, flush=True)
print("118c — TAÇ SONUÇ DOĞRULAMASI (taper'lı yedi ada + ζ)", flush=True)
print("=" * 78, flush=True)

KAYNAK, SAY, QQ, NHAM, BOLG = {}, {}, {}, {}, {}
d41 = np.load(HERE / "41_bigT_windows.npz")
K41 = sorted({x.split("_")[1] for x in d41.files}, key=lambda s: int(s[:-1]))
ZETA_PARCA = []
for k in K41[:6]:
    gz, tm = d41[f"gaps_{k}"], d41[f"tmid_{k}"]
    zz = np.empty(len(gz) + 1)
    zz[0] = tm[0] - gz[0] / 2
    zz[1:] = zz[0] + np.cumsum(gz)
    ZETA_PARCA.append(zz)
KAYNAK["zeta"] = ZETA_PARCA
SAY["zeta"] = rvm_sayim
QQ["zeta"] = 1.0
NHAM["zeta"] = sum(len(z) for z in ZETA_PARCA)
BOLG["zeta"] = 0

for et, q, tab, a, T1, _, _ in globals()["ADALAR"]:
    d = np.load(HERE / f"118b_{et}_zeros.npz")
    zc = d["zeros"]
    lo = np.exp(np.log(zc[0] + 1) + 0.25 * (np.log(zc[-1]) - np.log(zc[0] + 1)))
    zc = zc[zc >= lo]
    M = LmotorT(q, tab, a, 0.5)
    KAYNAK[et] = [zc]
    SAY[et] = (lambda t, M=M: M.theta(t) / np.pi)
    QQ[et] = float(q)
    NHAM[et] = len(zc)
    BOLG[et] = int(np.asarray(d["bolgeler"]).reshape(-1, 2).shape[0])


def boru(ad):
    """107a/107b'nin ÜÇ KAPILI boru hattı."""
    tan = {"g1": 0, "g2": 0, "g2ham": 0, "g3": 0, "dup": 0,
           "smax": 0.0, "rmax": [], "kusur": []}
    WIN = []
    for zc in KAYNAK[ad]:
        g1, g2, ham2, smax = kapi12(zc, SAY[ad])
        tan["smax"] = max(tan["smax"], smax)
        tan["g1"] += int(g1.sum())
        tan["g2"] += int(g2.sum())
        tan["g2ham"] += int(ham2.sum())
        for (i0, i1, t0, t1) in bloklar(ham2, zc):
            tan["kusur"].append(("K2", t0, t1, i1 - i0 + 1))
        seg, _ = segmentle(zc, g1 | g2)
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
            alt, bad3, kes3, rmax, u, dup = kapi3(p, QQ[ad], SAY[ad])
            tan["rmax"].append(rmax)
            tan["g3"] += int(kes3.sum())
            tan["dup"] += len(dup)
            for (i0, i1, t0, t1) in bloklar(bad3, p):
                tan["kusur"].append(("K3", t0, t1, i1 - i0 + 1))
            for x in alt:
                WIN.append(pencere_hazirla(x, QQ[ad]))
    return WIN, tan


PEN, TAN = {}, {}
for a in ADLAR:
    PEN[a], TAN[a] = boru(a)
    nB = sum(len(w[0]) for w in PEN[a])
    t = TAN[a]
    print(f"[{a}] ham n={NHAM[a]}  →  {len(PEN[a])} pencere n={nB}"
          f"   (maskelenen {NHAM[a]-nB}, %{100*(NHAM[a]-nB)/NHAM[a]:.1f}; "
          f"kalan bölge {BOLG[a]})", flush=True)
    print(f"      K1 {t['g1']:>6} | K2 ham {t['g2ham']:>4} → padli "
          f"{t['g2']:>6} (maks |med₂₀−med₂₀₀| = {t['smax']:.2f}) | "
          f"K3 padli {t['g3']:>6}, kopya {t['dup']}", flush=True)
    print(f"      K3 maks|r| pencere başına: "
          + " ".join(f"{x:.2f}" for x in t["rmax"]), flush=True)
    for (k, t0, t1, n) in t["kusur"]:
        print(f"      {k}: t ∈ [{t0:.4f}, {t1:.4f}] ({n} sıfır)", flush=True)

# ------------------------------------------------------------ ÖLÇÜM
print("\n--- D(τ) ölçülüyor ---", flush=True)
D, V, W = {}, {}, {}
for a in ADLAR:
    Dv, Vv, Wv = egri(PEN[a], TAUS)
    D[a], V[a], W[a] = np.array(Dv), np.array(Vv), Wv
    print(f"  ...{a} ({time.time()-T0:.0f} s)", flush=True)

LORT, TILK, NPEN, SIZ = {}, {}, {}, {}
for a in ADLAR:
    Ls = np.array([w[3] for w in PEN[a]])
    ns = np.array([len(w[0]) for w in PEN[a]], float)
    LORT[a] = float(np.average(Ls, weights=ns))
    TILK[a] = np.log(KERVAN[a]) / LORT[a]
    NPEN[a] = len(Ls)
    SIZ[a] = [int(np.sum(Ls * (t + BAND) >= np.log(KERVAN[a]))) for t in TAUS]

# 107b'nin ESKİ tablosu (karşılaştırma için)
E = np.load(HERE / "107b_yedi_ada.npz", allow_pickle=True)
DE = {a: E[f"D_{a}"] for a in ADLAR}
LE = {a: float(E["Lort"][list(E["adlar"]).index(a)]) for a in ADLAR}
YE = {a: float(E["yari"][list(E["adlar"]).index(a)]) for a in ADLAR}

print("\n" + "=" * 100, flush=True)
print(f"D(τ) — ESKİ (107b, keskin) vs YENİ (118b, taper), band ±{BAND}L",
      flush=True)
print("=" * 100, flush=True)
print(f"{'τ':>6} " + "  ".join(f"{a:>15}" for a in ADLAR), flush=True)
print(f"{'':>6} " + "  ".join(f"{'eski':>7}{'yeni':>8}" for a in ADLAR),
      flush=True)
for i, t in enumerate(TAUS):
    print(f"{t:>6.3f} " + "  ".join(
        f"{DE[a][i]:>7.3f}{D[a][i]:>8.3f}" for a in ADLAR), flush=True)
print("\nvekil taban (yeni): ", flush=True)
for a in ADLAR:
    print(f"  {a:>6}: " + " ".join(f"{v:.3f}" for v in V[a]), flush=True)


def yari_coz(Dv, T):
    Dv = np.asarray(Dv, float); T = np.asarray(T, float)
    ix = np.where(Dv >= 0.5)[0]
    if not len(ix) or ix[0] == 0:
        return np.nan
    i0 = ix[0]
    return T[i0 - 1] + (0.5 - Dv[i0 - 1]) * (T[i0] - T[i0 - 1]) / \
        (Dv[i0] - Dv[i0 - 1])


print("\n" + "=" * 100, flush=True)
print("YEDİ-ADA DONMA TABLOSU — YENİ (taper) / ESKİ (107b)", flush=True)
print("=" * 100, flush=True)
print(f"{'aile':>6} {'q(a)':>8} {'kervan':>7} {'n_pen':>5} {'⟨L⟩':>7} "
      f"{'τ_ilk':>7} {'D@log2/L':>17} {'D@τ_ilk':>17} {'yarı-τ':>16} "
      f"{'yarı/τ_ilk':>10}", flush=True)
YARI, TAB = {}, {}
for a in ADLAR:
    t2 = np.log(2) / LORT[a]
    d2 = float(np.interp(t2, TAUS, D[a]))
    d2e = float(np.interp(np.log(2) / LE[a], TAUS, DE[a]))
    di = float(np.interp(TILK[a], TAUS, D[a]))
    die = float(np.interp(np.log(KERVAN[a]) / LE[a], TAUS, DE[a]))
    ty = yari_coz(D[a], TAUS)
    YARI[a] = ty
    TAB[a] = (d2, di, ty)
    print(f"{a:>6} {QLAB[a]:>8} {'log'+str(KERVAN[a]):>7} {NPEN[a]:>5} "
          f"{LORT[a]:>7.3f} {TILK[a]:>7.4f} "
          f"{d2e:>7.3f}→{d2:>7.3f}  {die:>7.3f}→{di:>7.3f}  "
          f"{YE[a]:>7.4f}→{ty:>7.4f} {ty/TILK[a]:>10.3f}", flush=True)

canli = [a for a in ADLAR if KERVAN[a] == 2]
olu = [a for a in ADLAR if KERVAN[a] == 3]
canli_L = [a for a in canli if a != "zeta"]

print("\n" + "-" * 78, flush=True)
print("C1 — EŞİKSİZ AYRIŞTIRICI D(log2/L) (ORTAK çizgide)", flush=True)
print("-" * 78, flush=True)
for isim, grup in (("LOG2 KERVANI", canli), ("LOG3 KERVANI", olu)):
    print(f"  {isim}  YENİ: " + "  ".join(
        f"{a} {TAB[a][0]:.3f}" for a in grup), flush=True)
    print(f"  {'':<12}  ESKİ: " + "  ".join(
        f"{a} {float(np.interp(np.log(2)/LE[a], TAUS, DE[a])):.3f}"
        for a in grup), flush=True)
k_yeni = min(TAB[a][0] for a in canli) / max(TAB[a][0] for a in olu)
k_eski = (min(float(np.interp(np.log(2)/LE[a], TAUS, DE[a])) for a in canli)
          / max(float(np.interp(np.log(2)/LE[a], TAUS, DE[a])) for a in olu))
print(f"  KONTRAST (en düşük log2 / en yüksek log3): "
      f"ESKİ {k_eski:.1f}×  →  YENİ {k_yeni:.1f}×", flush=True)

print("\n" + "-" * 78, flush=True)
print("C2/C3 — ω_yarı = τ_yarı·⟨L⟩ (yarı-çözülmenin MUTLAK frekansı)",
      flush=True)
print("-" * 78, flush=True)
for a in ADLAR:
    print(f"    {a:>6}: ESKİ {YE[a]*LE[a]:.4f}  →  YENİ {YARI[a]*LORT[a]:.4f}",
          flush=True)
wcL = np.mean([YARI[a] * LORT[a] for a in canli_L])
wo = np.mean([YARI[a] * LORT[a] for a in olu])
wcL_e = np.mean([YE[a] * LE[a] for a in canli_L])
wo_e = np.mean([YE[a] * LE[a] for a in olu])
print(f"\n  LOG2 (4 Dirichlet) ort: ESKİ {wcL_e:.4f} → YENİ {wcL:.4f}   "
      f"log 2 = {np.log(2):.5f}   sapma ESKİ %{100*abs(wcL_e/np.log(2)-1):.2f}"
      f" → YENİ %{100*abs(wcL/np.log(2)-1):.2f}", flush=True)
print(f"  LOG3 (3 ada) ort      : ESKİ {wo_e:.4f} → YENİ {wo:.4f}   "
      f"log 3 = {np.log(3):.5f}   ERKENLİK ESKİ %{100*(1-wo_e/np.log(3)):.1f}"
      f" → YENİ %{100*(1-wo/np.log(3)):.1f}", flush=True)
print(f"  KERVAN KAYMASI: ESKİ ×{wo_e/wcL_e:.3f} → YENİ ×{wo/wcL:.3f}   "
      f"[saf öngörü log3/log2 = 1.585]", flush=True)
print(f"  yarı/τ_ilk  log2: " + " ".join(
    f"{a}={YARI[a]/TILK[a]:.3f}" for a in canli), flush=True)
print(f"              log3: " + " ".join(
    f"{a}={YARI[a]/TILK[a]:.3f}" for a in olu), flush=True)

print("\n" + "-" * 78, flush=True)
print("C4 — χ₃'ÜN ORTA-τ ÇUKURU: TEMİZ VERİDE NE OLDU?", flush=True)
print("-" * 78, flush=True)
print(f"  (107a ±0.01L üç kapılı: 0.426/0.199/0.161/0.134/0.531 "
      f"@ τ=0.095…0.140)", flush=True)
print(f"  χ₃ pencereleri (yeni): L = "
      + " ".join(f"{w[3]:.2f}" for w in PEN["chi3"]) + "   n = "
      + " ".join(str(len(w[0])) for w in PEN["chi3"]), flush=True)
print(f"{'τ':>6} {'ESKİ hav':>9} {'YENİ hav':>9} "
      f"{'pencere-başına (YENİ)':>30} {'maks_p':>7} {'hav/maks':>9}",
      flush=True)
for i, t in enumerate(TAUS):
    w = W["chi3"][i]
    mx = max(w)
    print(f"{t:>6.3f} {DE['chi3'][i]:>9.3f} {D['chi3'][i]:>9.3f}  "
          + " ".join(f"{x:>6.3f}" for x in w)
          + f"   {mx:>6.3f} {D['chi3'][i]/max(mx,1e-9):>8.2f}", flush=True)

print("\n  --- AYNI AYRIŞTIRMA, SEKİZ KÜME (yıkıcı girişim denetimi) ---",
      flush=True)
print(f"{'aile':>6} {'τ':>6} {'havuz':>7} {'maks_pen':>9} {'oran':>6}",
      flush=True)
for a in ADLAR:
    for i, t in enumerate(TAUS):
        if t not in (0.068, 0.105, 0.113, 0.125, 0.14):
            continue
        mx = max(W[a][i])
        print(f"{a:>6} {t:>6.3f} {D[a][i]:>7.3f} {mx:>9.3f} "
              f"{D[a][i]/max(mx,1e-9):>6.2f}", flush=True)

# ---------------------------------- C5: BAND DUYARLILIĞI / ANOMALİ DENETİMİ
# İlk koşuda χ₇'nin D(0.095)'i 0.527 → 0.036'ya düştü (tek nokta çukuru).
# 104'ün pratik kuralı 4 gereği önce PENCERE AYRIŞTIRMASI, sonra BAND
# duyarlılığı sorulur: havuz iki pencerenin de ALTINA düşüyorsa yıkıcı
# girişimdir, fiziksel değildir.
print("\n" + "-" * 78, flush=True)
print("C5 — BAND DUYARLILIĞI (χ₇'nin τ=0.095 çukuru + χ₃'ün orta-τ çukuru)",
      flush=True)
print("-" * 78, flush=True)
TT5 = [0.085, 0.090, 0.093, 0.095, 0.097, 0.100, 0.105]
for a in ("chi7", "chi3"):
    print(f"  [{a}] pencereler L = " + " ".join(f"{w[3]:.2f}" for w in PEN[a]),
          flush=True)
    for bnd in (0.005, 0.01, 0.02):
        Db, Vb, Wb = egri(PEN[a], TT5, band=bnd)
        print(f"    band ±{bnd:.3f}L: " + "  ".join(
            f"{t:.3f}:{Db[i]:.3f}[" + "/".join(f"{x:.2f}" for x in Wb[i]) + "]"
            for i, t in enumerate(TT5)), flush=True)
print("  OKUMA: köşeli parantez = pencere-başına D. Havuz İKİ pencerenin de"
      "\n  ALTINDAYSA yıkıcı girişim (aletsel); İKİSİ DE ÇUKURDAYSA gerçek.",
      flush=True)

# ----------------------------------------- x = τ/τ_ilk ÇÖKME
XG = np.arange(0.55, 1.751, 0.05)
MM = {a: np.interp(XG, np.array(TAUS) / TILK[a], D[a],
                   left=np.nan, right=np.nan) for a in ADLAR}
sp, ab = [], []
for j, x in enumerate(XG):
    vals = np.array([MM[a][j] for a in ADLAR])
    s, m = float(np.nanstd(vals)), float(np.nanmean(vals))
    sp.append(s / m if m > 0 else np.nan)
    ab.append(s)
sp, ab = np.array(sp), np.array(ab)
ham_sp = [float(np.std([D[a][i] for a in ADLAR])
                / np.mean([D[a][i] for a in ADLAR])) for i in range(len(TAUS))]
m1 = (np.array(TAUS) >= 0.06) & (np.array(TAUS) <= 0.125)
mx_ = (XG >= 0.849) & (XG <= 1.101)
print("\n" + "-" * 78, flush=True)
print("EK — YENİDEN ÖLÇEKLEME ÇÖKMESİ (x = τ/τ_ilk)", flush=True)
print("-" * 78, flush=True)
print(f"  HAM τ  (0.06-0.125): std/ort {np.mean(np.array(ham_sp)[m1]):.4f}   "
      f"[107b: 0.5125]", flush=True)
print(f"  ÖLÇEKLİ x (0.85-1.10): std/ort {np.nanmean(sp[mx_]):.4f}   "
      f"[107b: 0.2559]", flush=True)

print("\n=== BAND KAÇAĞI (★ = sızıntı) ===", flush=True)
print(f"{'τ':>6} " + " ".join(f"{a:>8}" for a in ADLAR), flush=True)
for i, t in enumerate(TAUS):
    print(f"{t:>6.3f} " + " ".join(
        f"{D[a][i]:>7.3f}" + ("★" if SIZ[a][i] > 0 else " ")
        for a in ADLAR), flush=True)

np.savez(HERE / "118c_yedi_ada.npz", taus=np.array(TAUS), band=BAND,
         adlar=np.array(ADLAR), xg=XG,
         Lort=np.array([LORT[a] for a in ADLAR]),
         tilk=np.array([TILK[a] for a in ADLAR]),
         qilk=np.array([KERVAN[a] for a in ADLAR]),
         yari=np.array([YARI[a] for a in ADLAR]),
         **{f"D_{a}": D[a] for a in ADLAR},
         **{f"V_{a}": V[a] for a in ADLAR},
         **{f"X_{a}": MM[a] for a in ADLAR},
         **{f"SIZ_{a}": np.array(SIZ[a]) for a in ADLAR},
         **{f"L_{a}": np.array([w[3] for w in PEN[a]]) for a in ADLAR},
         **{f"n_{a}": np.array([len(w[0]) for w in PEN[a]]) for a in ADLAR})
print(f"\n118c_yedi_ada.npz yazıldı. Süre {time.time()-T0:.0f} s.", flush=True)
