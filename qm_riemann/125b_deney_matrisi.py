"""
125b — DENEY MATRİSİ: ÇOKLU-MOD × FAZ YAPISI × GENLİK REJİMİ (26 Ağustos)
==========================================================================
BİLİM SORUSU (kalem defteri 3'ün son iki adayından biri): gerçek sıfır
gazı ~46 aritmetik çizgiyi ORTAK-ORİJİNLİ sinüs fazlarıyla birden taşır
ve nefesi R_p ≈ −0.9 (içsel β = −2). Tek-mod mıhlama (123) 10× termal
sabitlemede bile +1'de çakılı kaldı. ÇOKLU-MOD + FAZ YAPISI −2'yi
üretir mi? 124'ün ilk geçişi sonuçsuzdu (kabul %11, dengelenme şüpheli);
bu, 125a'nın ayarlı/kapılı aletiyle yapılan karar koşusudur.

MATRİS: {A-aritmetik-faz, B-karışık-faz} × {katı-genlik (124),
PERDELİ-genlik Ã_q·(1−0.36τ_q) (108'in ölçülen fiziksel perdesi)}
= 4 mühürlü hücre; her hücrede 3 bağımsız tohum.
+ ARTEFAKT KONTROL ÇİFTİ (mühürlü, 125a'nın ön-teşhisinden doğdu):
{A,B} × KESİK tayf (τ ≤ 0.25) — bu tayfta hedef profil ρ_hedef(θ) HER
YERDE POZİTİF, yani gaz onu gerçekten gerçekleyebilir. Tam tayfta (46
çizgi) A hedefi ρ_min ≈ −2.4'e düşer: fiziksel olmayan bir hedef.
Kontrol çifti A/B karşılaştırmasını bu konfundan arındırır.

ÖLÇÜM: her hücrede 5 hedef çizgide (m=17,27,40,48,63 ↔ p=2,3,5,7,13)
R_p (varyans kanalı) ve R_nn (bond kanalı — gerçekte en temiz gösterge,
120: R_nn ≈ −2cos(πτ)); zincirler-arası hata çubuğu (3 tohum, t₉₅).
Ayrıca patolojik bölgeyi (ρ_hedef < 0.25) dışlayan MASKELİ R_p.

ÖN-MÜHÜR (124'ten devralınan ayrıştırıcı, kapıya bağlı):
  Ç1  A negatife kayar, B +1'de kalır → ANOMALİNİN AJANI FAZ
      TUTARLILIĞI (tanımanın mekanik karşılığı — büyük kapanış).
  Ç2  İkisi de negatif → çoklu-mod/etkileşim yeter, faz gerekmez.
  Ç3  İkisi de +1 → H-ÇOKLU DA ÖLÜR; tek aday determinizm/iz-formülü.
  HÜKÜM KURALI: Ç1/Ç2/Ç3 ancak (a) R̂ < 1.1 kapısı geçilmişse ve
  (b) hata çubukları ayrımı destekliyorsa verilir. Desteklemiyorsa
  "SONUÇSUZ + neyin gerektiği" dürüstçe yazılır.

SONUÇ (26 Ağustos) — Ç3: H-ÇOKLU DA ÖLDÜ. BEŞİNCİ DARALMA.
  ALET SAĞLAMASI ✓✓ 123'ün tek-mod kurulumu bu aletle yeniden ölçüldü:
    A1 = 0.2915±0.0001, R_p = +0.962±0.031 (123: 0.289, +0.95), R̂=1.001.
    R zinciri sağlam — aşağıdaki sayılar ölçüm hatası değil.
  KAPI: 6 hücrenin 5'i geçti. TEK DÜŞEN A-katı (R̂[A1] max 1.297,
    n_eff 8) — tam da hedef profilin en derin negatife (−2.90) düştüğü,
    yani fiziksel olarak gerçeklenemez olan hücre. Kapı kusuru bulmakta
    beklendiği yerde çalıştı.
  KARAR TABLOSU R_p (hücre ort. ± t₉₅, 3 tohum):
    A-katı    +0.884±0.255 | B-katı    +0.670±0.120   (kapı A'da ✗)
    A-perdeli +0.842±0.160 | B-perdeli +0.631±0.131
    A-kesik   +0.808±0.029 | B-kesik   +0.805±0.032   ← FİZİKSEL HEDEF
  R_nn (bond): A-katı +1.77 | B-katı +0.34 | A-perdeli +1.43 |
    B-perdeli +0.15 | A-kesik +0.675±0.188 | B-kesik +0.654±0.082.
  HİÇBİR HÜCRE NEGATİF DEĞİL. Gerçek gaz R_p ≈ −0.88, R_nn ≈ −1.9..−1.4.
  Ç1 RET (kesin): hedefi fiziksel olan tek çiftte (kesik tayf, R̂≈1.00,
    n_eff 206-1500, hata ±0.03) A − B = +0.003 ± 0.043. ARİTMETİK FAZ
    TUTARLILIĞININ NEFESE ETKİSİ SIFIR. Üstelik 46-çizgili hücrelerde
    işaret Ç1'in BEKLEDİĞİNİN TERSİ: A hep B'den YÜKSEK (daha adyabatik).
  Ç2 RET: dört hücrede de (kapı geçenler) R_p pozitif; alet dalından
    ζ'ya kat edilen yol en fazla %15 (B-katı %12.7, B-perdeli %14.9) ve
    o da KARIŞIK fazlı hücrelerde; A hücrelerinde %0.6-3.0.
  Ç3 KABUL: çoklu-mod mıhlama + faz yapısı −2'yi ÜRETMİYOR; mıhlanmış
    gaz 1 modda da 9 modda da 46 modda da ADYABATİK DALDA kalıyor.
  ARTEFAKT TEŞHİSİ (125b_profiller.png): A hücrelerinde gaz koherans
    sivrisini GERÇEKTEN kuruyor (⟨ρ⟩ tepesi 4.9'a çıkıyor) ve hedefin
    negatif bölgesini BOŞLUKLA karşılıyor (ρ→0'da kırpılma). Yani
    sivrilme "patolojik ama gerçekleşen". Maskeli R_p (ρ_hedef<0.25
    bölgesi dışlanmış) hükümde HİÇBİR ŞEYİ değiştirmiyor (A-katı +0.906,
    B-katı +0.885, A-kesik +0.814) → sivrilme R'yi sürüklemiyor; hüküm
    patolojiye dayanmıyor. Kesik tayfta hedef ile gerçekleşen ⟨ρ⟩
    neredeyse tam örtüşüyor (mıh sadık).
  POST-HOC EK (mühürsüz, aşağıdaki tarama): MIHLANMIŞ DENGE GAZININ
    ADYABATİK DALI DÜZ +1 DEĞİL. Tek-mod, sabit Ã=0.281, τ∈[0.066,0.43]:
    R_p = +0.991 → +0.413; en iyi uyum R_p ≈ 1.136 − 1.604·τ (artık RMS
    0.024) — "sabit +1" RMS 0.195, "+1.06·cos(πτ)" RMS 0.085. Bond
    kanalı R_nn ≈ 1.20 − 4.09·τ, τ≈0.29'da SIFIRI GEÇİYOR ve τ=0.43'te
    −0.56'ya iniyor; ama gerçek gazın R_nn = −2cos(πτ)'si τ=0.066'da
    zaten −1.96 — iki eğri hiçbir τ'da örtüşmüyor. (İlk geçişte hücre
    hükümleri düz +1 referansına göre verilmiş ve "SONUÇSUZ" çıkmıştı;
    referans aletin kendi ölçtüğü dala düzeltildi — hüküm Ç3.)
  DÜRÜST SINIR (yapısal, 125a'nın ön-teşhisi): dairesel modelde
    m_q = round(N log q/L₀) TAM SAYI olmak zorunda; bu, her periyotta
    gerçek gazda bulunmayan bir eşzamanlı-koherans noktası yaratır ve
    rijit ötelemeyle kaldırılamaz (φ_q = c·m_q sadece profili kaydırır).
    Yani "ortak-orijinli aritmetik faz" bu model sınıfında ancak
    bozulmuş biçimde temsil edilebiliyor. Kesik-tayf çifti bu kusurdan
    arınmıştır ve orada etki TAM SIFIR ölçüldü; ama gerçek gazın
    ölçülemez-oranlı (incommensurate) faz örgüsünü sınamak periyodik
    olmayan bir gaz gerektirir — bu, açık iş olarak kalır.
"""

import importlib.util
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location(
    "m125a", str(HERE / "125a_ayarli_mcmc.py"))
M = importlib.util.module_from_spec(_spec)
sys.modules["m125a"] = M
_spec.loader.exec_module(M)

TWO_PI = M.TWO_PI
N = M.N
HEDEF = M.HEDEF

# ---------------------------- koşu bütçesi: 125a'nın G1/G2/G3 kapılarından
# G2: λ=4 tek değer hem mıh kapısını (σ_C/σ_termal max 0.121 ≤ 0.15) hem
#     kabul bandını (%34) geçiyor (λ=2 → 0.165, kapıda kalıyor).
# G3: 12000+20000 süpürmede, aşırı-dağıtılmış üç başlangıçtan R̂ ≤ 1.06
#     ve enerji izi düz (eğim ≤ 0.05 σ_E/1k süpürme) → bu bütçe yeterli.
LAM = 4.0
N_BURN = 12000
N_SAMP = 20000
HER = 40
TOHUMLAR = [1251, 1252, 1253]
JITTER = {1251: 1.0, 1252: 0.5, 1253: 0.0}   # aşırı-dağıtılmış başlangıç:
# hedef profilden / yarı yoldan / düzgün dağılımdan — R̂ gerçek bir sınav
# kalsın diye (aynı noktadan başlayan zincirlerin R̂'si iyimserdir).
T95_3 = 4.303                      # 3 tohum, 2 serbestlik derecesi

HUCRELER = [
    ("A-katı",     "A", dict(perde=False, taumax=0.5)),
    ("B-karışık-katı", "B", dict(perde=False, taumax=0.5)),
    ("A-perdeli",  "A", dict(perde=True, taumax=0.5)),
    ("B-karışık-perdeli", "B", dict(perde=True, taumax=0.5)),
    ("A-kesik(τ≤.25)", "A", dict(perde=True, taumax=0.25)),
    ("B-kesik(τ≤.25)", "B", dict(perde=True, taumax=0.25)),
]


def _hucre_tayf(kw):
    return M.tayf(kw["perde"], kw["taumax"])


# POZİTİF KONTROL: 123'ün tek-mod deneyi (m=26, C₀=36 ⇒ Ã=2C₀/N=0.281,
# C₀/σ_termal = 10) bu ALETLE yeniden ölçülür. 123 sonucu: A1=0.289,
# R_p=+0.95. Alet doğruysa aynı sayıyı vermeli — R zincirinin sağlaması.
KONTROL = dict(ad="kontrol-tekmod(123)", MODS=np.array([26]),
               AMPS=np.array([2 * 36.0 / N]), fit=np.array([26, 52]),
               hedefler=[26])


def is_zincir(arg):
    """Tek zincir: koş → mıh kapısı → R_p/R_nn (tam ve maskeli)."""
    hidx, tohum = arg
    ad, tip, kw = HUCRELER[hidx]
    MODS, AMPS = _hucre_tayf(kw)
    faz = M.fazlar_uret(tip, MODS, tohum)
    r = M.kos(MODS, AMPS, faz, tohum=tohum, LAM=LAM,
              n_burn=N_BURN, n_samp=N_SAMP, her=HER,
              jitter=JITTER.get(tohum, 1.0))
    oran, kay = M.mih_kapisi(r["ps_C"], r["ps_S"], MODS, AMPS, faz)
    out, se, c1s = M.R_olc(r["ornekler"], MODS)

    def maske(mid):
        rr = 1 + (AMPS[:, None] *
                  np.cos(MODS[:, None] * mid[None, :] - faz[:, None])).sum(0)
        return rr > 0.25
    outm, sem_, _ = M.R_olc(r["ornekler"], MODS, agirlik_fn=maske)

    _, rho_h = M.hedef_profil(MODS, AMPS, faz, 20001)
    return dict(
        hidx=hidx, tohum=tohum, ad=ad,
        kabul=r["kabul"], adim=r["adim"], sure=r["sure"], n_orn=r["n_orn"],
        ps_var=r["ps_var"], ps_A1=r["ps_A1"], rho=r["rho"],
        iz=M.durgunluk(r["iz_E"], r["iz_sw"], N_BURN, 25),
        oran_max=float(oran.max()), oran_med=float(np.median(oran)),
        kay_max=float(kay.max()),
        se=se, c1s=c1s, sem=sem_,
        A1=np.array([out[m]["A1"] for m in HEDEF]),
        R_p=np.array([out[m]["R_p"] for m in HEDEF]),
        R_nn=np.array([out[m]["R_nn"] for m in HEDEF]),
        R_p_msk=np.array([outm[m]["R_p"] for m in HEDEF]),
        rho_hedef_min=float(rho_h.min()),
        rho_hedef_neg=float((rho_h < 0).mean()),
        hedef_A=np.array([AMPS[list(MODS).index(m)] for m in HEDEF]),
    )


def is_kontrol(tohum):
    """123'ün tek-mod kurulumu, 125a aletiyle (pozitif kontrol)."""
    MODS = KONTROL["MODS"]; AMPS = KONTROL["AMPS"]
    faz = np.full(1, -np.pi / 2)
    r = M.kos(MODS, AMPS, faz, tohum=tohum, LAM=LAM,
              n_burn=N_BURN, n_samp=N_SAMP, her=HER,
              jitter=JITTER.get(tohum, 1.0))
    out, se, c1s = M.R_olc(r["ornekler"], KONTROL["fit"],
                           hedefler=KONTROL["hedefler"])
    oran, kay = M.mih_kapisi(r["ps_C"], r["ps_S"], MODS, AMPS, faz)
    return dict(tohum=tohum, kabul=r["kabul"], se=se, c1s=c1s,
                oran_max=float(oran.max()),
                A1=out[26]["A1"], R_p=out[26]["R_p"], R_nn=out[26]["R_nn"],
                ps_var=r["ps_var"], sure=r["sure"])


# ------------------------------------------------------------- POST-HOC EK
# MÜHÜRSÜZ (matrisin sonucundan DOĞDU, önceden kayıtlı değil): matrisin en
# temiz hücrelerinde R_p per-mod düşüşü R_p/cos(πτ) = 0.928±0.006 (A-kesik)
# ve 0.927±0.023 (B-kesik) gibi ŞAŞIRTICI DÜZ bir orana oturdu. Bu, gerçek
# gazın R_nn ≈ −2cos(πτ) kapalı formuyla (120) AYNI KİNEMATİK ÇARPAN.
# Estimator artefaktı mı, mıhlanmış gazın yasası mı? SINAV: 123'ün tek-mod
# kurulumunu SABİT genlikle geniş bir τ bandında tara. Yasa gerçekse
# R_p/cos(πτ) τ'dan bağımsız kalmalı ve τ→0.5'te R_p→0'a gitmeli.
POSTHOC_M = [17, 27, 40, 48, 63, 78, 96, 110]
POSTHOC_A = 2 * 36.0 / N          # 123 ile aynı genlik (Ã=0.281)


def is_posthoc(arg):
    m, tohum = arg
    MODS = np.array([m]); AMPS = np.array([POSTHOC_A])
    faz = np.full(1, -np.pi / 2)
    r = M.kos(MODS, AMPS, faz, tohum=tohum, LAM=LAM, n_burn=N_BURN,
              n_samp=N_SAMP, her=HER, jitter=JITTER.get(tohum, 1.0))
    out, se, c1s = M.R_olc(r["ornekler"], np.array([m, 2 * m]), hedefler=[m])
    oran, kay = M.mih_kapisi(r["ps_C"], r["ps_S"], MODS, AMPS, faz)
    return dict(m=m, tohum=tohum, kabul=r["kabul"], se=se, c1s=c1s,
                oran_max=float(oran.max()), A1=out[m]["A1"],
                R_p=out[m]["R_p"], R_nn=out[m]["R_nn"], ps_var=r["ps_var"])


def ort_hata(v):
    """3 tohumluk diziden (ortalama, t₉₅ yarı-genişlik)."""
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    if len(v) < 2:
        return (float(v[0]) if len(v) else np.nan), np.nan
    return float(v.mean()), float(T95_3 * v.std(ddof=1) / np.sqrt(len(v)))


# ---------------------------------------------------------------- referanslar
# POST-HOC taramasından (bu SCRIPT'in kendi aleti, tek-mod, Ã=0.281, λ=4):
# mıhlanmış DENGE gazının adyabatik dalı τ'da DÜZ DEĞİL, düşüyor —
# R_p ≈ 1.136 − 1.604·τ (artık RMS 0.024; "sabit +1" varsayımı RMS 0.195).
# Bu yüzden "denge" sınavı düz +1'e değil, ALETİN KENDİ ÖLÇTÜĞÜ tek-mod
# dalına karşı yapılır. (İlk geçişte düz +1 kullanılmış ve hücreler
# +0.81'de oldukları için "SONUÇSUZ" etiketi almışlardı — o etiket yanlış
# referanstan geliyordu, düzeltildi.)
REF_TEK_P = {17: 0.991, 27: 0.958, 40: 0.899, 48: 0.857, 63: 0.770}
REF_TEK_N = {17: 0.931, 27: 0.823, 40: 0.628, 48: 0.510, 63: 0.232}
REF_TEK_P_ORT = float(np.mean(list(REF_TEK_P.values())))     # +0.895
REF_TEK_N_ORT = float(np.mean(list(REF_TEK_N.values())))     # +0.625
REF_ZETA_P = -0.88                                            # 109/120
REF_ZETA_N = float(np.mean([-2 * np.cos(np.pi * m / N) for m in HEDEF]))
REF_TOL = 0.10          # tek-mod referansının kendi belirsizlik payı


def sinif(mu, h, ref):
    """CI'ye göre etiket; 'denge' sınavı ALETİN tek-mod dalına karşı."""
    if not np.isfinite(h):
        return "?"
    if mu + h < 0:
        return "NEGATİF"
    if abs(mu - ref) < h + REF_TOL:
        return "ADYABATİK"
    if mu - h > 0:
        return "POZİTİF-ama-sapmış"
    return "AYRIŞMIYOR"


def mesafe(mu, ref, hedef):
    """Adyabatik daldan gerçek gaza kat edilen yolun oranı (%)."""
    if abs(ref - hedef) < 1e-9:
        return np.nan
    return 100.0 * (ref - mu) / (ref - hedef)


def posthoc_kos():
    """POST-HOC: tek-mod mıhlanmış gazda R_p(τ) yasası taraması."""
    import multiprocessing as mp
    t0 = time.time()
    print("=" * 78)
    print("POST-HOC EK (MÜHÜRSÜZ) — TEK-MOD MIHLANMIŞ GAZDA R_p(τ) YASASI")
    print(f"Ã={POSTHOC_A:.3f} sabit, λ={LAM}, {N_BURN}+{N_SAMP} süpürme, "
          f"3 tohum; m = {POSTHOC_M}")
    print("=" * 78, flush=True)
    isler = [(m, t) for m in POSTHOC_M for t in TOHUMLAR]
    with mp.get_context("fork").Pool(8) as pool:
        S = pool.map(is_posthoc, isler)
    K = {}
    for s in S:
        K.setdefault(s["m"], []).append(s)
    print(f"{'m':>4} {'τ':>6} {'cos(πτ)':>8} {'kabul':>6} {'σC/σt':>6} "
          f"{'R̂var':>6} {'A1':>6} {'R_p':>16} {'R_p/cos':>15} "
          f"{'R_nn':>16} {'R_nn/cos':>10}")
    orn = []
    for m in POSTHOC_M:
        Z = K[m]
        tau = m / N; c = float(np.cos(np.pi * tau))
        rp, hp = ort_hata([z["R_p"] for z in Z])
        rn, hn = ort_hata([z["R_nn"] for z in Z])
        rv, _ = M.rhat_neff(np.vstack([z["ps_var"] for z in Z]))
        print(f"{m:>4} {tau:>6.3f} {c:>8.4f} "
              f"{np.mean([z['kabul'] for z in Z]):>6.3f} "
              f"{max(z['oran_max'] for z in Z):>6.3f} {rv:>6.3f} "
              f"{np.mean([z['A1'] for z in Z]):>6.3f} "
              f"{rp:>+9.3f}±{hp:>5.3f} {rp/c:>+9.3f}±{hp/c:>4.3f} "
              f"{rn:>+9.3f}±{hn:>5.3f} {rn/c:>+10.3f}", flush=True)
        orn.append((tau, c, rp, hp, rn))
    tau = np.array([o[0] for o in orn]); c = np.array([o[1] for o in orn])
    rp = np.array([o[2] for o in orn]); rn = np.array([o[4] for o in orn])
    kp = float((rp * c).sum() / (c * c).sum())
    kn = float((rn * c).sum() / (c * c).sum())
    res_p = float(np.sqrt(np.mean((rp - kp * c) ** 2)))
    # rakip model: R_p sabit (τ'dan bağımsız)
    sab = float(rp.mean()); res_s = float(np.sqrt(np.mean((rp - sab) ** 2)))
    # rakip model: lineer 1−aτ
    Ax = np.vstack([np.ones_like(tau), tau]).T
    bl = np.linalg.lstsq(Ax, rp, rcond=None)[0]
    res_l = float(np.sqrt(np.mean((rp - Ax @ bl) ** 2)))
    print(f"\n  EN İYİ UYUM  R_p = ({kp:+.3f})·cos(πτ)   artık RMS {res_p:.4f}")
    print(f"  rakip: R_p = sabit {sab:+.3f}            artık RMS {res_s:.4f}")
    print(f"  rakip: R_p = {bl[0]:+.3f} {bl[1]:+.3f}·τ        "
          f"artık RMS {res_l:.4f}")
    print(f"  bond:  R_nn = ({kn:+.3f})·cos(πτ)")
    print(f"\n  GERÇEK GAZ (120): R_nn = −2·cos(πτ);  bu gaz: "
          f"R_p = {kp:+.2f}·cos(πτ)  → ORAN {(-2.0)/kp:+.2f}")
    print(f"\nsüre {time.time()-t0:.0f} sn")


if __name__ == "__main__":
    import multiprocessing as mp

    if "posthoc" in sys.argv:
        posthoc_kos()
        sys.exit(0)

    t00 = time.time()
    print("=" * 78)
    print("125b — DENEY MATRİSİ (çoklu-mod × faz × genlik)")
    print(f"λ={LAM}, {N_BURN} ısınma + {N_SAMP} örnekleme (her {HER}), "
          f"{len(TOHUMLAR)} tohum × {len(HUCRELER)} hücre")
    print("=" * 78, flush=True)

    print("\n[HÜCRE TAYFLARI ve MIHLANAN HEDEF PROFİL]")
    print(f"{'hücre':>20} {'n_mod':>6} {'ΣA':>6} {'dalga payı':>11} "
          f"{'min ρ_hedef':>12} {'ρ<0':>7}")
    for ad, tip, kw in HUCRELER:
        MODS, AMPS = _hucre_tayf(kw)
        faz = M.fazlar_uret(tip, MODS, TOHUMLAR[0])
        _, rho = M.hedef_profil(MODS, AMPS, faz, 20001)
        print(f"{ad:>20} {len(MODS):>6} {AMPS.sum():>6.2f} "
              f"{0.5*(AMPS**2).sum():>11.4f} {rho.min():>+12.3f} "
              f"{(rho<0).mean():>7.4f}", flush=True)
    print("  (ζ dalga payı ≈ 0.142; ρ_hedef<0 ⇒ hedef fiziksel değil)")

    print("\n[POZİTİF KONTROL — 123'ün tek-mod kurulumu bu aletle] "
          "(beklenen: A1≈0.289, R_p≈+0.95)", flush=True)
    with mp.get_context("fork").Pool(3) as pool:
        KON = pool.map(is_kontrol, TOHUMLAR)
    ka, kh = ort_hata([k["A1"] for k in KON])
    ra, rh = ort_hata([k["R_p"] for k in KON])
    na, nh = ort_hata([k["R_nn"] for k in KON])
    rk_v, ne_k = M.rhat_neff(np.vstack([k["ps_var"] for k in KON]))
    print(f"  kabul={np.mean([k['kabul'] for k in KON]):.3f}  "
          f"σC/σt={max(k['oran_max'] for k in KON):.3f}  "
          f"σ_η={np.mean([k['se'] for k in KON]):.3f}  "
          f"c₁/σ²={np.mean([k['c1s'] for k in KON]):+.3f}  "
          f"R̂[Var]={rk_v:.3f}")
    print(f"  A1 = {ka:.4f}±{kh:.4f}   R_p = {ra:+.3f}±{rh:.3f}   "
          f"R_nn = {na:+.3f}±{nh:.3f}")
    kon_ok = np.isfinite(rh) and abs(ra - 0.95) < max(rh, 0.15)
    _ks = "✓ 123 yeniden üretildi" if kon_ok else \
        "✗ 123 YENİDEN ÜRETİLEMEDİ — R zinciri şüpheli"
    print(f"  ALET SAĞLAMASI: {_ks}", flush=True)

    isler = [(h, t) for h in range(len(HUCRELER)) for t in TOHUMLAR]
    print(f"\n{len(isler)} zincir koşuyor (6 paralel)...", flush=True)
    with mp.get_context("fork").Pool(6) as pool:
        SON = pool.map(is_zincir, isler)
    print(f"  bitti ({time.time()-t00:.0f} sn)", flush=True)

    KUT = {}
    for s in SON:
        KUT.setdefault(s["hidx"], []).append(s)

    # ---------------- yakınsama kapısı ----------------
    print("\n" + "=" * 78)
    print("YAKINSAMA KAPISI (R̂ < 1.1; 3 zincir, bölünmüş)")
    print("=" * 78)
    print(f"{'hücre':>20} {'kabul':>6} {'adım/ḡ':>7} {'σC/σt':>6} "
          f"{'|kay|':>6} {'σ_η':>6} {'c₁/σ²':>7} {'R̂var':>6} "
          f"{'R̂A1max':>8} {'n_eff':>7} {'iz eğim':>8} {'kapı':>5}")
    KAPI = {}
    for h, ad_tip in enumerate(HUCRELER):
        ad = ad_tip[0]
        Z = sorted(KUT[h], key=lambda s: s["tohum"])
        rh_v, ne_v = M.rhat_neff(np.vstack([z["ps_var"] for z in Z]))
        rhs = []; nes = []
        for k in range(len(HEDEF)):
            rr, nn = M.rhat_neff(np.vstack([z["ps_A1"][:, k] for z in Z]))
            rhs.append(rr); nes.append(nn)
        rh_max = float(np.nanmax([rh_v] + rhs))
        ne_min = float(np.nanmin([ne_v] + nes))
        ok = (rh_max < M.RHAT_KAPI) and \
             all(0.25 <= z["kabul"] <= 0.55 for z in Z) and \
             max(z["oran_max"] for z in Z) <= M.MIH_KAPI
        KAPI[h] = dict(ok=ok, rh_var=rh_v, rh_max=rh_max, ne_min=ne_min)
        print(f"{ad:>20} {np.mean([z['kabul'] for z in Z]):>6.3f} "
              f"{np.mean([z['adim'] for z in Z]):>7.4f} "
              f"{max(z['oran_max'] for z in Z):>6.3f} "
              f"{max(z['kay_max'] for z in Z):>6.3f} "
              f"{np.mean([z['se'] for z in Z]):>6.3f} "
              f"{np.mean([z['c1s'] for z in Z]):>+7.3f} "
              f"{rh_v:>6.3f} {float(np.nanmax(rhs)):>8.3f} "
              f"{ne_min:>7.0f} "
              f"{np.mean([z['iz'] for z in Z]):>+8.3f} "
              f"{'✓' if ok else '✗':>5}", flush=True)

    # ---------------- karar tablosu ----------------
    def tablo(alan, baslik, ref_mod, ref_ort, zeta):
        print("\n" + "=" * 78)
        print(baslik)
        print("=" * 78)
        bas = f"{'hücre':>20}"
        for m in HEDEF:
            bas += f" {('m='+str(m)):>14}"
        bas += f" {'HÜCRE ORT.':>16} {'hüküm':>20} {'ζ yolu%':>8}"
        print(bas)
        sat = f"{'· ALET tek-mod dalı':>20}"
        for m in HEDEF:
            sat += f" {ref_mod[m]:>+8.2f}      "
        sat += f" {ref_ort:>+9.3f}       {'(referans)':>20} {0.0:>8.1f}"
        print(sat)
        sat = f"{'· GERÇEK GAZ (ζ)':>20}"
        for m in HEDEF:
            sat += f" {zeta(m):>+8.2f}      "
        zo = float(np.mean([zeta(m) for m in HEDEF]))
        sat += f" {zo:>+9.3f}       {'(hedef)':>20} {100.0:>8.1f}"
        print(sat)
        hucre_ort = {}
        for h, ad_tip in enumerate(HUCRELER):
            ad = ad_tip[0]
            Z = sorted(KUT[h], key=lambda s: s["tohum"])
            sat = f"{ad:>20}"
            for k in range(len(HEDEF)):
                mu, hh = ort_hata([z[alan][k] for z in Z])
                sat += f" {mu:>+7.2f}±{hh:>5.2f}" if np.isfinite(hh) \
                    else f" {mu:>+7.2f}±  ?  "
            mu, hh = ort_hata([np.mean(z[alan]) for z in Z])
            hucre_ort[h] = (mu, hh)
            sat += (f" {mu:>+9.3f}±{hh:>5.3f} {sinif(mu, hh, ref_ort):>20}"
                    f" {mesafe(mu, ref_ort, zo):>8.1f}")
            print(sat, flush=True)
        return hucre_ort

    ORT_p = tablo("R_p", "KARAR TABLOSU — R_p (VARYANS KANALI)",
                  REF_TEK_P, REF_TEK_P_ORT, lambda m: REF_ZETA_P)
    ORT_n = tablo("R_nn", "KARAR TABLOSU — R_nn (BOND KANALI)",
                  REF_TEK_N, REF_TEK_N_ORT,
                  lambda m: -2 * np.cos(np.pi * m / N))
    ORT_m = tablo("R_p_msk", "SAĞLAMLIK — MASKELİ R_p (ρ_hedef<0.25 "
                             "bölgesi dışlanmış)",
                  REF_TEK_P, REF_TEK_P_ORT, lambda m: REF_ZETA_P)
    print("\n  'ζ yolu%' = adyabatik daldan gerçek gaza kat edilen yolun "
          "yüzdesi (0 = alet dalında, 100 = ζ'da).")

    print("\n[GERÇEKLEŞEN vs HEDEF DALGA GENLİĞİ] (mıh ne kadar tutuyor?)")
    print(f"{'hücre':>20}" + "".join(f" {('m='+str(m)):>13}" for m in HEDEF))
    for h, ad_tip in enumerate(HUCRELER):
        Z = sorted(KUT[h], key=lambda s: s["tohum"])
        sat = f"{ad_tip[0]:>20}"
        for k in range(len(HEDEF)):
            a = np.mean([z["A1"][k] for z in Z])
            t = Z[0]["hedef_A"][k]
            sat += f" {a:>6.3f}/{t:<6.3f}"
        print(sat)
    print("  (gerçekleşen A1 / mıhlanan Ã; R_p bu A1'e normalize edilir)")

    # ---------------- hüküm ----------------
    print("\n" + "=" * 78)
    print("HÜKÜM")
    print("=" * 78)
    AILE = [("KATI GENLİK (124 rejimi)", 0, 1),
            ("PERDELİ GENLİK (gerçek dalga düzeyi)", 2, 3),
            ("KESİK TAYF τ≤0.25 (FİZİKSEL HEDEF — artefakt kontrolü)", 4, 5)]
    for ad, ha, hb in AILE:
        mua, ha_ = ORT_p[ha]; mub, hb_ = ORT_p[hb]
        na, nb = ORT_n[ha][0], ORT_n[hb][0]
        kapi = KAPI[ha]["ok"] and KAPI[hb]["ok"]
        ayr = np.isfinite(ha_) and np.isfinite(hb_) and \
            abs(mua - mub) > (ha_ + hb_)
        sa = sinif(mua, ha_, REF_TEK_P_ORT)
        sb = sinif(mub, hb_, REF_TEK_P_ORT)
        poz = ("ADYABATİK", "POZİTİF-ama-sapmış")
        if not kapi:
            hk = "HÜKÜM ASKIDA — yakınsama kapısı geçilmedi"
        elif sa == "NEGATİF" and sb in poz and ayr:
            hk = "Ç1 — AJAN FAZ TUTARLILIĞI"
        elif sa == "NEGATİF" and sb == "NEGATİF":
            hk = "Ç2 — çoklu-mod/etkileşim yeter, faz gerekmez"
        elif sa in poz and sb in poz:
            hk = ("Ç3 — H-ÇOKLU ÖLÜR (her iki faz yapısı da adyabatik "
                  "dalda); tek aday determinizm/iz-formülü")
        else:
            hk = "SONUÇSUZ (hata çubukları ayrımı desteklemiyor)"
        c1a = np.mean([z["c1s"] for z in KUT[ha]])
        c1b = np.mean([z["c1s"] for z in KUT[hb]])
        nn_ok = min(abs(c1a), abs(c1b)) > 0.05
        print(f"\n  {ad}")
        print(f"    A: R_p = {mua:+.3f}±{ha_:.3f} [{sa}], R_nn = {na:+.3f}")
        print(f"    B: R_p = {mub:+.3f}±{hb_:.3f} [{sb}], R_nn = {nb:+.3f}")
        print(f"    bond kanalı c₁/σ² = {c1a:+.3f}/{c1b:+.3f} → R_nn "
              f"{'güvenilir' if nn_ok else 'GÜVENİLMEZ (c₁≈0, payda patlar)'}")
        dab, hab = mua - mub, np.hypot(ha_, hb_)
        print(f"    A−B = {dab:+.3f}±{hab:.3f} → ayrım "
              f"{'VAR' if ayr else 'YOK'}   kapı: {'✓' if kapi else '✗'}")
        print(f"    adyabatik daldan (alet tek-mod {REF_TEK_P_ORT:+.3f}) "
              f"gerçek gaza ({REF_ZETA_P:+.2f}) kat edilen yol: "
              f"A %{mesafe(mua, REF_TEK_P_ORT, REF_ZETA_P):.1f}, "
              f"B %{mesafe(mub, REF_TEK_P_ORT, REF_ZETA_P):.1f}")
        print(f"    → {hk}", flush=True)

    # ---------------- profil figürü ----------------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axs = plt.subplots(3, 2, figsize=(13, 10), sharex=True)
        for h, (ad, tip, kw) in enumerate(HUCRELER):
            ax = axs[h // 2, h % 2]
            MODS, AMPS = _hucre_tayf(kw)
            Z = sorted(KUT[h], key=lambda s: s["tohum"])
            faz = M.fazlar_uret(tip, MODS, Z[0]["tohum"])
            th, rho_h = M.hedef_profil(MODS, AMPS, faz, 4001)
            nb = len(Z[0]["rho"])
            thb = (np.arange(nb) + 0.5) * TWO_PI / nb
            ax.axhline(0, color="k", lw=0.6)
            ax.plot(th, rho_h, lw=1.0, color="crimson", label="hedef ρ")
            ax.plot(thb, np.mean([z["rho"] for z in Z], axis=0), lw=1.2,
                    color="navy", label="gerçekleşen ⟨ρ⟩")
            ax.set_title(f"{ad}  (min ρ_hedef = {Z[0]['rho_hedef_min']:+.2f})",
                         fontsize=9)
            ax.set_ylabel("ρ/ρ̄")
            if h == 0:
                ax.legend(fontsize=8)
        for ax in axs[2]:
            ax.set_xlabel("θ")
        fig.suptitle("125b — mıhlanan hedef profil vs gerçekleşen ⟨ρ(θ)⟩ "
                     "(yoğunluk-sivrilmesi teşhisi)")
        fig.tight_layout()
        fig.savefig(HERE / "125b_profiller.png", dpi=110)
        print(f"\n  figür: 125b_profiller.png")
    except Exception as e:
        print(f"\n  [figür atlandı: {e}]")

    np.savez(HERE / "125b_matris.npz",
             hucreler=np.array([h[0] for h in HUCRELER]),
             tohumlar=np.array(TOHUMLAR), hedef=np.array(HEDEF),
             R_p=np.array([[z["R_p"] for z in sorted(KUT[h],
                            key=lambda s: s["tohum"])]
                           for h in range(len(HUCRELER))]),
             R_nn=np.array([[z["R_nn"] for z in sorted(KUT[h],
                             key=lambda s: s["tohum"])]
                            for h in range(len(HUCRELER))]),
             R_p_msk=np.array([[z["R_p_msk"] for z in sorted(KUT[h],
                                key=lambda s: s["tohum"])]
                               for h in range(len(HUCRELER))]),
             A1=np.array([[z["A1"] for z in sorted(KUT[h],
                           key=lambda s: s["tohum"])]
                          for h in range(len(HUCRELER))]),
             rho=np.array([[z["rho"] for z in sorted(KUT[h],
                            key=lambda s: s["tohum"])]
                           for h in range(len(HUCRELER))]))
    print("  önbellek: 125b_matris.npz")
    print(f"\nTOPLAM SÜRE: {time.time()-t00:.0f} sn")
