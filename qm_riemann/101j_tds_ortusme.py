"""
101j — TDS ÖRTÜŞMESİ: ÇÖZÜLME YAMACI, ÇİZGİNİN OMZU MU? (20 Ağustos)
==========================================================================
KEŞİFSEL. Hipotez: çözülme yamacının BİÇİMİ, ilk sağ kalan çizginin
termal-difüz-saçılma (TDS) omzudur. Kırınımda bir Bragg beneğinin
çevresinde fonon TDS'i I ∝ 1/Δk² ile düşer; burada "benek" ω=log2
çizgisi, "TDS" ise çizgi çevresindeki difüz omuz. Eğer donma→çözülme
geçişi çizginin omzuna binmekse, D(τ)'nun yamacı |Ĝ(ω)|²'nin omuz
profilini İZLEMELİ.

101h ve 101i bu hipoteze iki bağımsız ipucu vermişti:
  - 101h: eşik-altı D'nin bir kısmı bandın ilk çizgiden KAÇAĞIYDI.
  - 101i: vuruş (moiré) hipotezi REDDEDİLDİ; yerine D(ω) taraması
    çizgiler arasında pürüzsüz bir VADİ çizdi.

ÖLÇÜMLER
  A) OMUZ PROFİLİ: I(ω) = |Σ_n w_n e^{iωt_n}|² / Σ_n w_n², ζ'nın
     41-pencere düz segmentlerinde, pencereler arası ŞİDDET ortalaması.
     İKİ AYARDA: (i) w=1 (101d/87 geleneği, TAPERSİZ), (ii) w=Hann.
     Taper şart, çünkü tapersiz sonlu pencere ω=0'daki dev tepeyi
     1/ω² kuyruğuyla her yere sızdırır — "karanlık taban" sanılan şey
     bu olabilir. Kontrol: Poisson → I≈1, kusursuz örgü → I≈0.
  B) ÇÖZÜLME YAMACI (τ): 101h'nin ±0.005L dar bandı, ζ, τ_ilk
     çevresinde ince ızgara; x = τ·⟨L⟩ − log2.
  C) ÇÖZÜLME YAMACI (sabit-ω): 101i'nin egri_omega'sı, band ±0.02,
     ω ∈ log2±0.15. x = ω − log2 DOĞRUDAN, ⟨L⟩ çevrimi yok.
  E) [EK, ön-mühür sonrası] TAPER KAPISI: D(ω) ölçümü r ve G
     toplamlarını TAPERSİZ alıyor. Aynı ölçüm Hann taperiyle
     tekrarlanır. Eşik-altı D taperle çökerse, "yumuşak öncü"
     fiziksel değil PENCERE KAÇAĞIDIR.
  D) ÖRTÜŞME: profiller kendi tabanlarına normalize; serbest parametre
     yok. log(P−taban) ~ p·log|Δω| eğimleri karşılaştırılır.

ÖN-MÜHÜR (ölçümden ÖNCE yazıldı):
  T1  I(ω) log2 çevresinde karanlık tabanın belirgin üstünde bir omuz
      gösterecek, Δω büyüdükçe düşecek.
  T2  D profilinin (B ve C) biçimi omuz profilini izleyecek.
  T3  B ve C birbirini tutacak (hizalama kapısı).
  KEŞİFSEL iş: nitel/yarı-nicel hüküm yeterli.

==========================================================================
SONUÇ (20 Ağustos, koşu 114 s) — T1 RET, T2 RET, T3 ✓, E sonuçsuz.
TDS ÖRTÜŞMESİ YOK; ÇÜNKÜ ÖLÇÜLEBİLİR BİR TDS OMZU DA YOK.
==========================================================================
KALİBRASYON KAPISI ✓ (sekil_faktoru doğru mu?):
  Poisson süreci → I = 0.9924 / 1.1104 / 1.0761 (ω≈0.35/0.55/0.95) ✓
  KUSURSUZ ÖRGÜ  → I = 1.01e-3 / 4.13e-4 / 1.41e-4 (tapersiz)
  Dikkat: kusursuz örgünün "tabanı" sıfır değil ~1e-3 — bu tamamen
  pencere kenarı artefaktıdır. Ve ζ'nın tapersiz tabanı da AYNI:
  9.39e-4 (ω=0.35), 3.79e-4 (0.45), 5.83e-4 (0.55), 3.29e-4 (0.90),
  1.93e-4 (0.95). ζ'nın karanlık alanı KUSURSUZ ÖRGÜDEN AYIRT
  EDİLEMİYOR → ölçülen taban fizik değil, alet.

A) OMUZ PROFİLİ — HÜKÜM: ÖLÇÜLEN "OMUZ" TAMAMEN ALETİN KENDİSİ.
  Hann taperli I(ω), log2 çevresinde 15 KADEME düşüyor; düşüş Hann
  çekirdeğinin yan-lobu olan Δω^{−6} yasasını birebir izliyor
  (uydurulan üs Δω∈[0.0005,0.05] için p = −6.18; yerel eğimler
  −6.9, −6.2, −6.1, −5.9, −6.6, −5.0, −6.9, −5.8) ve iki yaka
  4 haneye kadar SİMETRİK:
    Δω     I(log2−Δω)  I(log2+Δω)   tapersiz(sağ)  p_yerel(tapersiz)
    0      5.23e+01 (çizgi)          7.84e+01
    0.0005 2.075e-01   2.075e-01     1.011e+00        —
    0.001  1.757e-03   1.757e-03     3.180e-01      −1.67
    0.002  2.414e-05   2.414e-05     8.007e-02      −1.99
    0.005  9.162e-08   9.162e-08     1.665e-02      −1.71
    0.010  1.515e-09   1.516e-09     5.110e-03      −1.70
    0.020  1.610e-11   1.611e-11     1.054e-03      −2.28
    0.030  2.106e-12   2.110e-12     7.712e-04      −0.77
    0.050  6.295e-14   6.149e-14     4.410e-04      −1.09
    0.080  3.944e-15   3.959e-15     3.224e-04      −0.67
    0.120  6.838e-16   1.021e-15     4.241e-04      +0.68
    0.150  4.588e-16   1.642e-15     3.383e-04      −1.01
  Δω≳0.08'de taperli I, uzak kontrol düzeyine (2.6e-16 … 2.6e-15)
  oturuyor. FİZİKSEL DİFÜZ BİLEŞEN YOK: karanlık alan ölçebildiğimiz
  kadarıyla GERÇEKTEN KARANLIK (I ≲ 3e-16, yani tapersiz "taban"ın
  ~10^12 katı altında). Açık formülle tutarlı: sıfır örgüsünün
  dalgalanma tayfı yalnız ω=log q çizgilerinde destekli, arada
  sürekli bileşen yok.
  ★ EN TEHLİKELİ BULGU: TAPERSİZ ölçümün yerel üssü çizgi civarında
  p ≈ −1.7 … −2.3, yani TAM OLARAK TDS'nin beklediği 1/Δω² YASASI.
  Ama bu, dikdörtgen pencerenin kendi çekirdeğidir (sinc² ~ Δω^−2),
  fizik değil. Tapersiz bir ölçüm bu hipotezi SAHTE OLARAK
  DOĞRULARDI. Taper olmadan TDS iddiası kurulamaz.
  ARTEFAKT ŞÜPHESİ (gizlemiyorum): 87_kapali_devre.py'nin T3 hükmü
  ("karanlık alan tabanı 2e-4/5e-4; taşıyıcılar bunun 4-15 katı →
  TDS") da TAPERSİZ Ĝ ile kurulmuştu ve o taban tam bu sızıntı
  düzeyinde. 87'nin yeniden denetlenmesi gerekir (bu işin kapsamı
  dışında; ayrı bir iş olarak not edildi).

B/C) ÇÖZÜLME YAMACI (x = Δω = ω − log2; D(τ) dar bant ±0.005L,
  D(ω) sabit-ω ±0.02; vekil taban 0.005-0.040):
    x      D(τ)  D(τ)/tab   D(ω)  D(ω)/tab
  −0.15   0.044    1.00    0.040    1.00
  −0.10   0.205    4.64    0.072    1.82
  −0.05   0.437    9.90    0.255    6.46
  −0.02   0.484   10.97    0.738   18.66
  −0.01   0.507   11.48    0.823   20.82
  +0.01   0.484   10.95    0.734   18.57
  +0.03   0.491   11.11    0.646   16.34
  +0.05   0.491   11.11    0.398   10.06
  +0.10   0.424    9.61    0.226    5.73
  +0.15   0.329    7.44    0.192    4.84
  T1 RET — taperli omuz yok; tapersiz "omuz" Δω≈0.03'te zaten tabana
    inip 0.15'e dek DÜZ kalıyor (p_yerel ≈ 0).
  T2 RET — örtüşme yok, ölçek uyuşmazlığı ezici: aynı Δω aralığında
    (0.01 → 0.15) taperli omuz 13 KADEME, tapersiz omuz ~1 kademe,
    D(ω) ise yalnızca 0.63 kademe (20.8× → 4.8×) düşüyor. Üsler:
    D(τ) p(sol) −1.31 / p(sağ) −0.21; D(ω) p(sol) −2.53 / p(sağ)
    −0.82; taperli omuz p ≈ −6. D(ω)'nun sol yakasındaki −2.53
    TDS'ye benziyor ama sağ yakası −0.82 ve asimetri büyük; omuz ise
    kusursuz simetrik. Biçim örtüşmesi YOK.
  T3 ✓ — B ile C uyuşuyor: r(D(τ), D(ω)) = +0.744 (n=30). ⟨L⟩ ile
    hizalama kendi başına artefakt üretmiyor; kapı geçildi.

E) TAPER KAPISI — TEKNİK OLARAK BAŞARISIZ, SONUÇSUZ (dürüst kayıt):
  D estimatörünü Hann taperiyle koşma denemesi ÇÖKTÜ. Sebep: payda
  den = Σ|r|², taperli hâlde karanlık bölgede ~1e-15'e iniyor, yani
  D = |pay|/payda bir 0/0 oranına dönüşüyor. Sayılar bunu açıkça
  itiraf ediyor: D_hann 0.87'den 4322'ye kadar savruluyor ve VEKİL
  TABAN da onunla birlikte savruluyor (3.3 → 2365). D_hann/vekil
  oranı taranan HER x'te 1'in ALTINDA (0.22-0.85) — yani taperli
  ölçümde hiçbir yerde sinyal yok, sadece gürültü.
  HÜKÜM: bu kapı ne doğruluyor ne yalanlıyor; estimatör taper'a
  taşınabilir değil, yeniden tasarım ister (payda için ayrı,
  taperlenmemiş bir normalizasyon ya da çizgi-modeli çıkarma).
  Bu yüzden "eşik-altı D artığı sızıntıdır" iddiası KANITLANMADI —
  yalnızca I(ω) düzeyinde tabanın sızıntı olduğu kesinleşti.
  (r(D(τ), D_hann) = −0.176; r(D(ω)ham, D_hann) = −0.324 — beklendiği
  gibi bağlantısız.)

GENEL HÜKÜM: "Çözülme yamacı ilk çizginin TDS omzudur" hipotezi
REDDEDİLDİ, ama beklenenden ilginç bir sebeple: ÖLÇÜLEBİLİR BİR TDS
OMZU YOK. Karanlık alan gerçekten karanlık (I ≲ 3e-16, kusursuz
örgüden ayırt edilemez); çizgi, pencere çözünürlüğünde saf bir delta.
Dolayısıyla çözülme yamacının GENİŞLİĞİ fiziksel bir difüz omuzla
açıklanamaz. Neyle açıklandığı AÇIK KALDI (E kapısı çöktüğü için).
101d'nin ana bulgusu bundan ETKİLENMEZ: o, aynı aletle ölçülen
aileler ARASI karşılaştırmadır ve 101i'nin ölü-çizgi kontrastı
(ζ 1.93× vs β 0.97×, aynı ω, aynı band) doğrudan fizikseldir.
Etkilenen tek şey, yamacın genişliğinin fiziksel yorumu.
"""

import numpy as np
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
rng = np.random.default_rng(1012)
LOG2 = float(np.log(2))

PKS = [2,3,4,5,7,8,9,11,13,16,17,19,23,25,27,29,31,32,37,41,43,47,49,53,
       59,61,64,67,71,73,79,81,83,89,97,101,103,107,109,113,121,125,127,128]
LINES_ALL = [np.log(qq) for qq in PKS]

# ---- 101d/101h/101i makinesi (kopya) ----

def pencere_hazirla(z, qeff):
    gaps = np.diff(z)
    mids = 0.5 * (z[:-1] + z[1:])
    Lw = float(np.log(qeff * mids / TWO_PI).mean())
    ds = gaps * np.log(qeff * mids / TWO_PI) / TWO_PI - 1
    tt = (mids - mids.mean()) / (mids[-1] - mids[0])
    P = np.vstack([np.ones_like(tt), tt, tt**2, tt**3]).T
    ds = ds - P @ np.linalg.lstsq(P, ds, rcond=None)[0]
    return (z, mids, ds, Lw)

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

def hann(t):
    return 0.5 * (1 - np.cos(TWO_PI * (t - t[0]) / (t[-1] - t[0])))

def egri(WIN, taus, band=0.005):
    """101h'nin dar-bant D(τ)'su (tapersiz — özgün makine)."""
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
        Dv.append(abs(num) / (c * den)); Vv.append(float(np.mean(fl)))
    return Dv, Vv

def egri_omega(WIN, omegas, band=0.02, nfrek=160, taper=False):
    """101i'nin sabit-ω D(ω)'sı; taper=True ise Hann ile."""
    Dv, Vv = [], []
    for om in omegas:
        num = 0.0 + 0j; den = 0.0
        Gs, rs = [], []
        for (zz, tm, ds, Lw) in WIN:
            wz = hann(zz) if taper else 1.0
            wm = hann(tm) if taper else 1.0
            tau_w = om / Lw
            kap = 2 * np.pi * tau_w
            c = 2 * np.sin(kap / 2) / kap
            oms = om + np.linspace(-band, band, nfrek)
            oms = np.array([o for o in oms
                            if min(abs(o - l) for l in LINES_ALL) > 0.01])
            for s0 in range(0, len(oms), 40):
                ob = oms[s0:s0 + 40]
                rr = (np.exp(1j * np.outer(ob, zz)) * wz).sum(axis=1)
                GG = (np.exp(1j * np.outer(ob, tm))
                      * (ds * wm)[None, :]).sum(axis=1) / c
                num += (GG * np.conj(rr)).sum()
                den += (np.abs(rr)**2).sum()
                Gs.append(GG); rs.append(rr)
        G_all, r_all = np.concatenate(Gs), np.concatenate(rs)
        fl = [abs((G_all * np.conj(r_all[rng.permutation(len(r_all))])).sum())
              / den for _ in range(40)]
        Dv.append(abs(num) / den); Vv.append(float(np.mean(fl)))
    return Dv, Vv

def sekil_faktoru(WIN, omegas, taper=True, chunk=64):
    """I(ω) = |Σ w_n e^{iωt_n}|² / Σ w_n², pencereler arası ortalama."""
    omegas = np.asarray(omegas)
    out = np.zeros(len(omegas))
    for (zz, tm, ds, Lw) in WIN:
        w = hann(zz) if taper else np.ones_like(zz)
        nrm = (w**2).sum()
        acc = np.empty(len(omegas))
        for s0 in range(0, len(omegas), chunk):
            ob = omegas[s0:s0 + chunk]
            S = (np.exp(1j * np.outer(ob, zz)) * w[None, :]).sum(axis=1)
            acc[s0:s0 + chunk] = np.abs(S)**2 / nrm
        out += acc
    return out / len(WIN)

# ------------------------- ζ pencereleri -------------------------

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
Ls = np.array([w[3] for w in WINZ])
ns = np.array([len(w[0]) for w in WINZ], float)
LORT = float(np.average(Ls, weights=ns))
print(f"[zeta] {len(WINZ)} pencere  ⟨L⟩_n={LORT:.4f}  "
      f"τ_ilk={LOG2/LORT:.5f}", flush=True)

# ---- KALİBRASYON: sekil_faktoru doğru mu? (Poisson ≈ 1, örgü ≈ 0) ----
z0 = WINZ[0][0]
rg = np.random.default_rng(7)
poi = [(np.sort(rg.uniform(z0[0], z0[-1], len(z0))), None, None, WINZ[0][3])]
lat = [(np.linspace(z0[0], z0[-1], len(z0)), None, None, WINZ[0][3])]
# tek bir ω'da I üstel dağılımlı (bağıl std %100) — bandda ortalama şart
for om0 in [0.35, 0.55, 0.95]:
    kal = om0 + np.linspace(-0.01, 0.01, 201)
    print(f"[KALİBRASYON ω≈{om0}] Poisson I = "
          f"{sekil_faktoru(poi, kal, taper=False).mean():.4f} (≈1 olmalı) | "
          f"kusursuz örgü I = "
          f"{sekil_faktoru(lat, kal, taper=False).mean():.3e} (≈0 olmalı)",
          flush=True)

# ---------------- A) OMUZ PROFİLİ ----------------
DW = np.array([0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.03, 0.05,
               0.08, 0.12, 0.15])
OM_SIM = np.concatenate([LOG2 - DW[::-1], [LOG2], LOG2 + DW])
I_TAP = sekil_faktoru(WINZ, OM_SIM, taper=True)
I_HAM = sekil_faktoru(WINZ, OM_SIM, taper=False)
UZAK = np.array([0.30, 0.35, 0.45, 0.55, 0.90, 0.95, 1.00])
IU_TAP = sekil_faktoru(WINZ, UZAK, taper=True)
IU_HAM = sekil_faktoru(WINZ, UZAK, taper=False)

print("\n=== A) OMUZ PROFİLİ: log2 çevresi, Hann taperli vs TAPERSİZ ===")
print(f"  çizgide (Δω=0): taperli {I_TAP[len(DW)]:.4e}  "
      f"tapersiz {I_HAM[len(DW)]:.4e}")
print(f"{'Δω':>8} {'tap(sol)':>11} {'tap(sağ)':>11} {'p_yerel':>8} "
      f"{'ham(sol)':>11} {'ham(sağ)':>11} {'p_yerel':>8}")
for j, d in enumerate(DW):
    a, b = I_TAP[len(DW)-1-j], I_TAP[len(DW)+1+j]
    ah, bh = I_HAM[len(DW)-1-j], I_HAM[len(DW)+1+j]
    if j:
        pl = np.log(b / I_TAP[len(DW)+j]) / np.log(d / DW[j-1])
        ph = np.log(bh / I_HAM[len(DW)+j]) / np.log(d / DW[j-1])
    else:
        pl = ph = np.nan
    print(f"{d:>8.4f} {a:>11.3e} {b:>11.3e} {pl:>8.2f} "
          f"{ah:>11.3e} {bh:>11.3e} {ph:>8.2f}")
print("  (p_yerel = ardışık iki nokta arası log-log eğimi, sağ yaka;")
print("   Hann yan-lobu ⇒ p=−6, fiziksel TDS omzu ⇒ p≈−2, düz taban ⇒ p≈0)")
print(f"  UZAK KONTROL (çizgilerden uzak ω):")
for o, a, b in zip(UZAK, IU_TAP, IU_HAM):
    print(f"    ω={o:.2f}  taperli {a:.3e}   tapersiz {b:.3e}")
lg = np.log(DW[DW <= 0.051]); li = np.log(I_TAP[len(DW)+1:][DW <= 0.051])
print(f"  TAPERLİ OMUZ ÜSSÜ (Δω∈[0.0005,0.05], sağ): "
      f"p = {np.polyfit(lg, li, 1)[0]:+.2f}  (Hann yan-lobu ⇒ −6)")

# ---------------- B) ÇÖZÜLME YAMACI (τ) ----------------
XT = np.round(np.arange(-0.15, 0.15001, 0.01), 4)
TAUS = (LOG2 + XT) / LORT
DB, VB = egri(WINZ, list(TAUS), band=0.005)
print("\n[B] dar-bant D(τ) bitti", flush=True)

# ---------------- C) ÇÖZÜLME YAMACI (sabit-ω, tapersiz) ----------------
XW = XT[np.abs(XT) > 0.005]
DC, VC = egri_omega(WINZ, list(LOG2 + XW), band=0.02, taper=False)
print("[C] sabit-ω D(ω) (tapersiz) bitti", flush=True)

# ---------------- E) TAPER KAPISI ----------------
DE, VE = egri_omega(WINZ, list(LOG2 + XW), band=0.02, taper=True)
print("[E] sabit-ω D(ω) (Hann) bitti", flush=True)

# ---------------- D) ÖRTÜŞME ----------------
TABAN_B, TABAN_C, TABAN_E = min(DB), min(DC), min(DE)
print(f"\n=== D) ÖRTÜŞME (x = Δω) ===")
print(f"  tabanlar: D(τ)_min={TABAN_B:.4f}  D(ω)_min={TABAN_C:.4f}  "
      f"D_hann(ω)_min={TABAN_E:.4f}")
print(f"{'x=Δω':>7} {'D(τ)':>7} {'/tab':>6} {'D(ω)':>7} {'/tab':>6} "
      f"{'vekil':>6} | {'D_hann':>7} {'/tab':>6} {'vekil_h':>7} {'ham/hann':>8}")
for i, x in enumerate(XT):
    j = np.where(np.abs(XW - x) < 1e-6)[0]
    if len(j):
        j = j[0]
        print(f"{x:>7.2f} {DB[i]:>7.3f} {DB[i]/TABAN_B:>6.2f} "
              f"{DC[j]:>7.3f} {DC[j]/TABAN_C:>6.2f} {VC[j]:>6.3f} | "
              f"{DE[j]:>7.3f} {DE[j]/TABAN_E:>6.2f} {VE[j]:>7.3f} "
              f"{DC[j]/DE[j]:>8.2f}")
    else:
        print(f"{x:>7.2f} {DB[i]:>7.3f} {DB[i]/TABAN_B:>6.2f} "
              f"{'—':>7} {'—':>6} {'—':>6} | {'—':>7} {'—':>6} "
              f"{'—':>7} {'—':>8}")

print("\n=== E2) TAPER KAPISI HÜKMÜ ===")
uz = np.abs(XW) >= 0.09      # eşikten uzak (donuk taraf ve ötesi)
yk = np.abs(XW) <= 0.03      # çizgi komşuluğu
DCa, DEa = np.array(DC), np.array(DE)
print(f"  çizgi komşuluğu |Δω|≤0.03 : tapersiz {DCa[yk].mean():.3f}  "
      f"Hann {DEa[yk].mean():.3f}  oran {DCa[yk].mean()/DEa[yk].mean():.2f}")
print(f"  uzak       |Δω|≥0.09     : tapersiz {DCa[uz].mean():.3f}  "
      f"Hann {DEa[uz].mean():.3f}  oran {DCa[uz].mean()/DEa[uz].mean():.2f}")
print(f"  KONTRAST (yakın/uzak): tapersiz "
      f"{DCa[yk].mean()/DCa[uz].mean():.2f}×  "
      f"Hann {DEa[yk].mean()/DEa[uz].mean():.2f}×")
print(f"  vekil taban: tapersiz {np.mean(VC):.3f}  Hann {np.mean(VE):.3f}")

def us_uydur(x, y, taban, isaret):
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = (np.sign(x) == isaret) & (np.abs(x) >= 0.019) & (np.abs(x) <= 0.151)
    yy = y[m] - taban
    m2 = yy > 0
    if m2.sum() < 4:
        return np.nan
    return float(np.polyfit(np.log(np.abs(x[m][m2])), np.log(yy[m2]), 1)[0])

print("\n=== D2) ÜS KARŞILAŞTIRMASI (log(P−taban) ~ p·log|Δω|) ===")
print(f"  D(τ)      : p(sol) {us_uydur(XT, DB, TABAN_B, -1):+.2f}   "
      f"p(sağ) {us_uydur(XT, DB, TABAN_B, +1):+.2f}")
print(f"  D(ω) ham  : p(sol) {us_uydur(XW, DC, TABAN_C, -1):+.2f}   "
      f"p(sağ) {us_uydur(XW, DC, TABAN_C, +1):+.2f}")
print(f"  D(ω) Hann : p(sol) {us_uydur(XW, DE, TABAN_E, -1):+.2f}   "
      f"p(sağ) {us_uydur(XW, DE, TABAN_E, +1):+.2f}")
print(f"  I(ω) taperli omuz: p ≈ −6 (pencere yan-lobu, fizik değil)")

print("\n=== D3) HİZALAMA KAPISI (Pearson r, ortak x) ===")
ortk = [i for i, x in enumerate(XT) if np.any(np.abs(XW - x) < 1e-6)]
b = np.array([DB[i] for i in ortk])
print(f"  r(D(τ), D(ω) ham)  = {np.corrcoef(b, DCa)[0,1]:+.3f}  (n={len(b)})")
print(f"  r(D(τ), D(ω) Hann) = {np.corrcoef(b, DEa)[0,1]:+.3f}")
print(f"  r(D(ω) ham, Hann)  = {np.corrcoef(DCa, DEa)[0,1]:+.3f}")

# ---------------- figür ----------------
fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
sag = OM_SIM[len(DW)+1:] - LOG2
ax[0].loglog(sag, I_TAP[len(DW)+1:], "o-", ms=4, color="C0",
             label="I(ω) Hann taperli")
ax[0].loglog(sag, I_HAM[len(DW)+1:], "s-", ms=4, color="C1",
             label="I(ω) TAPERSİZ (sızıntı)")
ax[0].loglog(sag, I_TAP[len(DW)+1] * (sag / sag[0])**-6.0, "k:",
             lw=1, label="Δω⁻⁶ (Hann yan-lobu)")
ax[0].set_xlabel("Δω = ω − log2"); ax[0].set_ylabel("I(ω)")
ax[0].set_title("ζ: log2 çizgisinin 'omzu' = pencere yan-lobu")
ax[0].legend(fontsize=8); ax[0].grid(alpha=0.25, which="both")
ax[1].plot(XT, DB, "o-", ms=4, color="C3", label="D(τ) dar bant ±0.005L")
ax[1].plot(XW, DC, "s-", ms=4, color="C0", label="D(ω) sabit-ω ±0.02")
ax[1].plot(XW, VC, ":", color="0.5", label="vekil taban")
ax[1].axvline(0, color="k", lw=0.8, ls=":")
ax[1].set_xlabel("Δω = ω − log2"); ax[1].set_ylabel("D")
ax[1].set_title("çözülme yamacı (özgün, tapersiz estimatör)")
ax[1].legend(fontsize=8); ax[1].grid(alpha=0.25)
ax[2].semilogy(XW, DE, "^-", ms=4, color="C2", label="D(ω) Hann")
ax[2].semilogy(XW, VE, ":", color="0.5", label="vekil taban (Hann)")
ax[2].axvline(0, color="k", lw=0.8, ls=":")
ax[2].set_xlabel("Δω = ω − log2"); ax[2].set_ylabel("D (log)")
ax[2].set_title("E) TAPER KAPISI ÇÖKTÜ: D ≲ vekil, her yerde")
ax[2].legend(fontsize=8); ax[2].grid(alpha=0.25, which="both")
fig.tight_layout()
fig.savefig(HERE / "101j_tds_ortusme.png", dpi=130)

np.savez(HERE / "101j_tds.npz", om_sim=OM_SIM, I_tap=I_TAP, I_ham=I_HAM,
         uzak=UZAK, IU_tap=IU_TAP, IU_ham=IU_HAM, xt=XT,
         DB=np.array(DB), VB=np.array(VB), xw=XW, DC=DCa, VC=np.array(VC),
         DE=DEa, VE=np.array(VE), Lort=LORT)
print("\n101j_tds_ortusme.png + 101j_tds.npz yazıldı.")
