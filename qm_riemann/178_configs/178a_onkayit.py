# -*- coding: utf-8 -*-
"""
178a — ÖN-KAYIT: BANT-AĞIRLIKLI θ-KESTİRİCİSİ (θ_bw) ve H-F1b^ba
=================================================================
ÖN-MÜHÜR (koşudan önce yazıldı; bu docstring mühürlenen metnin
parçasıdır — betiğin sha256'sı bu satırları da kapsar).

══════════════════════════════════════════════════════════════════
0. AÇIK BEYAN — BU BİR KESTİRİCİ DEĞİŞİKLİĞİDİR
══════════════════════════════════════════════════════════════════
BU BİR VARYANS-KÜÇÜLTME KESTİRİCİ DEĞİŞİKLİĞİDİR; EŞİK GEVŞETME
DEĞİLDİR. Değişen tek şey θ'nın BÖLENİNDEKİ g_X'in okunuşudur:
tüm-ızgara tek-skaler EKK eğimi yerine, DONMUŞ hüküm penceresinin
beş bandında bant-başına eğim + ikizin jackknife hatasından gelen
ters-varyans ağırlıklı birleşim. Eşiklerin FORMÜLLERİ 176'nın F3
formülleriyle HARF HARF aynıdır; yalnız içine giren çapa sayıları
yeni para biriminde yeniden hesaplanır (Ş1.1, Ş1.2). Hiçbir eşik
sayısı elle girilmez, hiçbir eşik gevşetilmez.

177'nin "H-F1b KALICI HÜKÜMSÜZ" hükmü ESKİ kestirici için AYAKTA
KALIR ve geri alınmaz. 178'in sınadığı hipotez ayrı adlıdır:
    H-F1b^ba  —  "θ'nın kilit-duyarlılığı, bant-ağırlıklı g_X
                  para biriminde"
Bu, aynı dört tohuma İKİNCİ BAKIŞTIR: bakış sayacı k = 2.

Gerekçe (177e tanısı): Var(log θ)'nin %92.4'ü g_X'ten, %7.5'i M'den,
%0.08'i g_E'den geliyordu. Asgari cerrahi: yalnız g_X'e dokunulur.
M (= KALİB_u2, lo=0.60) ve g_E BİT DÜZEYİNDE DEĞİŞMEZ (Ş1.4).

══════════════════════════════════════════════════════════════════
1. KÖRLÜK — BU BETİK VEKİL DOSYASI AÇMAZ
══════════════════════════════════════════════════════════════════
VF1, VF2, VF3, VF4 adları bu betikte YALNIZCA yasak-listesi ve
dürüstlük beyanı olarak geçer. Çalışma zamanında bir bekçi
(`_ac`) her dosya açılışını süzer: yolunda vekil adı geçen hiçbir
dosya açılamaz (SystemExit). Kestirici kodu mühürden önce hiçbir
vekil üzerinde koşmaz — smoke test dahil (Ş2.1).

Çapalar (Hkeskin = ikiz, son = gerçek, HA4 = erfc ikizi) SERBESTTİR:
onlar sınanan nicelik değil, hüküm CETVELİDİR (Ş2.2). Cetvel
mühürden önce çizilir ki vekil ölçümü cetveli oynatamasın.

DÜRÜSTLÜK BEYANI (Ş2.3) — körlük İDDİASI YOKTUR. 176/177'nin bütün
vekil sayıları bu seferin tayfası tarafından okunmuştur; ön-kayıt
onları `durustluk` bloğunda listeler. Özellikle ESKİ vekil g_X'leri
(0.8134745 / 0.8333309 / 0.8154829 / 0.8068090) bilinmektedir —
değişen bileşen tam da g_X olduğu için bu, ağırlıkları sonuca doğru
eğme riskidir. Buna karşı tek güvence Ş3'ün KÖKEN kuralıdır:
ağırlıklar YALNIZ İKİZDEN gelir, vekillerin hiçbir sayısı ağırlığa
giremez; ve ağırlıklar bu ön-kayıtta SAYI olarak donar.

══════════════════════════════════════════════════════════════════
2. KESTİRİCİNİN KESİN TANIMI (tek birincil kestirici — Ş2.4)
══════════════════════════════════════════════════════════════════
ÖLÇÜM ZİNCİRİ DEĞİŞMEZ. Her gaz için 167/165 zinciri aynen koşar:
    Y  = 165_cekirdek.gaz(ad, taban=0.40, cap=4000)
    Mo = 165_cekirdek.Model165(ad, 0.40, 4000, 0.95, "olculen",
                               Y=Y, tau_c=0.95).sec("olculen", 0.95)
    Mo.alanlar(kmax=3)
Bu, 176c_olcum'un `O167.kos(ad, 0.40, 0.95, 0, 0)` çağrısının alan
katmanıyla BİREBİR aynı çağrıdır; ürettiği g_E ve global g_X,
diskteki C_<ad>.json["artik"] değerleriyle 1e−12 içinde eşleşmek
ZORUNDADIR (Ş1.7 kimlik kapısı).

BANT KÜMESİ (SABİT, tüm gazlarda AYNI — Ş3.3):
    IZGARA_T1 = izgara(0.44, 0.80, 0.04)'ün lo ∈ {0.52, 0.56, 0.60,
    0.64, 0.68} bantları  ⇒  τ ∈ (0.52, 0.72], beş bant.
    Üyelik kuralı  lo < τ_q ≤ hi  — 163_cekirdek.bant_adaylari'nin
    kuralıyla birebir.
GEREKÇE (yalnız DONMUŞ kararlara atıf — Ş3.5): (i) 167_olcum.LO_MIN
= 0.52, ölçüm makinesinin kendi donmuş "hüküm penceresi" tabanı;
(ii) 176c_olcum.LO_LO/LO_HI = 0.52/0.68, 176a'nın F0 hüküm
penceresi; (iii) M'nin okunduğu lo = 0.60 bandı bu pencerenin tam
merkezidir — pay ve payda AYNI donmuş pencerede yaşar. İkinci aday
(dokuz bandın tamamı) DENENMEZ; iki listeyi deneyip beğenileni
seçmek Ş3.3(ii) gereği yasaktır.

BANT-BAŞINA EĞİM (donmuş global g_X ile AYNI cebirsel biçim):
    msk_b   = (τ_q > lo) ∧ (τ_q ≤ hi)          [hepsi τ ≤ τ_c = 0.95]
    X_b     = sentez(s, w[msk_b], y[msk_b]) ,  X_b ← X_b − ⟨X_b⟩
    x1      = Y.Xtil0
    g_X,b  := ⟨X_b · x1⟩ / ⟨X_b · X_b⟩
(165_cekirdek.alanlar satır 287'nin gX = dot(X,x1)/dot(X,X) biçimi,
banda kısıtlanmış hâli. Toplam alan X = Σ_b X_b doğrusal olarak
ayrışır; ortalama çıkarma her banda ayrı uygulanır.)

BİRLEŞİM (doğrusal, ters-varyans):
    g_X,bw := Σ_b W_b · g_X,b ,   Σ_b W_b = 1 ,   W_b ∝ 1/s_b²
`s_b` YALNIZ İKİZDEN (Hkeskin) gelir (Ş3.2) ve bu ön-kayıtta SAYI
olarak donar; ölçüm sırasında ASLA yeniden hesaplanmaz. Vekillerin
hiçbir sayısı ve tohumlar-arası varyans ağırlığa GİREMEZ.

İKİZİN JACKKNIFE HATASI s_b (166_T1._jk biçimi, njack = 8):
    bandın çizgileri q'ya göre artan sıralanır; grup(i) = i mod 8
    (163_cekirdek.bant_adaylari'nin round-robin şeması, njack=8);
    j = 0..7 için g_X,b^(−j) = bandın j'inci grubu SİLİNMİŞ eğimi;
    s_b = sqrt( (n−1)/n · Σ_j (v_j − v̄)² ) ,  n = 8.
VARYANS TABANI (floor) YOKTUR; kırpma YOKTUR; log-uzay birleşim
DEĞİL doğrusal birleşim (tek seçim, burada donuyor — Ş3.5).

θ'NIN YENİ PARA BİRİMİ:
    θ_bw := M / ( g_E · g_X,bw² ) ,
    M    := KALİB_u2(lo = 0.60)  — C_<ad>.json'dan, DEĞİŞMEDEN,
    g_E  := artik.gE             — DEĞİŞMEDEN.
Özdeşlik θ_bw ≡ M/(g_E·g_X,bw²) yapı gereği kalıntısızdır (Ş1.4) ve
her gazda denetlenir.

İNDİRGEME BELGESİ (Ş3.1): W_b yerine GÜÇ ağırlıkları ⟨X_b,X_b⟩
konsaydı birleşim Σ_b⟨X_b,x1⟩ / Σ_b⟨X_b,X_b⟩ olurdu; bu, bant
desteği S üzerindeki "köşegen-Gram" EKK eğimidir ve Gram çaprazları
sıfır olsa ve S bütün ızgarayı kapsasaydı donmuş global g_X'i
AYNEN geri verirdi. Üç sayı (global g_X, güç-ağırlıklı g_X^S,
kısıtlı-destek EKK g_X^{S,EKK}) her çapa için raporlanır.

══════════════════════════════════════════════════════════════════
3. EŞİKLER — 176'NIN FORMÜLLERİ, YENİ ÇAPALAR (Ş1.2, Ş1.3)
══════════════════════════════════════════════════════════════════
    D        := log θ_bw(son) − log θ_bw(Hkeskin)
    eşik (a) := θ_bw(Hk)·e^{−D} ≡ θ_bw(Hk)² / θ_bw(son)
    eşik (b) := θ_bw(Hk)
    E        := log θ_bw(HA4) − log θ_bw(Hkeskin)
    f_θ      := D / E
    kesim_θ  := f_θ · D
    kilit_θ  := (1 − f_θ) · D
Sağlama: eski para biriminde bu formüller 0.8933338²/0.9219382 =
0.8656170 ve 0.8933338 verir — 176'nın F3 eşikleriyle özdeş.

GEÇERLİLİK ÖNKOŞULLARI (Ş4.0, mühür anında çapalardan görülür):
    P1: D > 0        değilse H-F1b^ba TANIMSIZ — SORU DÜŞTÜ.
    P2: kilit_θ > 0  değilse H2 TANIMSIZ (yalnız H1 koşulur).

══════════════════════════════════════════════════════════════════
4. HÜKÜM KURALI (176/177'den DEĞİŞMEDEN taşınır — Ş1.5)
══════════════════════════════════════════════════════════════════
y_i := θ_bw(VF_i), i = 1..4 (DÖRDÜ DE girer; eleme/ağırlıklama/
dışlama YASAK). ȳ = ort(y), s = std(y, ddof=1),
    SAÇ_n := 2·s/√n ,  n = 4  ⇒  SAÇ_4 = s      (n = 4 TAVANDIR;
beşinci tohum yok, merdiven yok, yeni gaz koşusu yok — Ş0.2.)
Katı okuma SAÇ^kat = max − min BAĞLAYICI DEĞİL, yalnız raporlanır.

H1 (176-F3'ün yeni-para karşılığı):
  (a) ȳ ≤ eşik_a  VE  |ȳ − eşik_a| ≥ SAÇ_4   ⇒  H-F1b^ba YAŞADI
  (b) ȳ ≥ eşik_b  VE  |ȳ − eşik_b| ≥ SAÇ_4   ⇒  H-F1b^ba ÖLDÜ
  aksi her durum                              ⇒  HÜKÜMSÜZ, ve n=4
      tavanı + Ş2.6 gereği DERHAL KALICI HÜKÜMSÜZ (2. para
      biriminde de). Bu dört tohumda ÜÇÜNCÜ bir kestirici DENENMEZ;
      soru ancak yeni veriyle (yeni tohum/yeni gaz, ayrı sefer,
      ayrı bütçe) açılabilir.

H2 (yalnız P2 sağlanırsa): Δ_i := log θ_bw(Hk) − log θ_bw(VF_i),
Δ̄ = ort(Δ), SAÇ_4(Δ) = std(Δ, ddof=1), ω_θ = kilit_θ / Δ̄.
  üst kenar: Δ̄ ≥ kilit_θ  (⟺ ω ≤ 1) ;  alt kenar: Δ̄ > 0 (⟺ ω > 0)
  her kenarda KESİNLİK ⟺ marj ≥ SAÇ_4(Δ).
  iki kenar da kesin ve sağlanıyor ⇒ FATURA ÖDENEBİLİR — KESİN
  ilgili kenar kesin ihlâlde       ⇒ ÖDENEMEZ — KESİN
  aksi                             ⇒ H2 HÜKÜMSÜZ
Tohum-başına ω dağılımı raporlanır; medyan-temelli hiçbir sayı
hükme giremez.

H1 yalnız başına H-F1b^ba'nın YAŞADI/ÖLDÜ hükmünü verir; H2 yalnız
F9'un θ satırını etiketler. Dokuz bileşimin sonuç cümlesi aşağıda
`bilesim_tablosu`'nda ÖNCEDEN yazılıdır — ölçümden sonra cümle
SEÇİLMEZ.

KAPI ÖLÜMÜ (Ş4.4): bir VEKİLDE bir bantta ⟨X_b,X_b⟩ ≤ 0 ya da
θ_bw tanımsız kalırsa o vekil ortalamadan DÜŞÜRÜLEMEZ; H1 ve H2 o
anda HÜKÜMSÜZ yazılır ve kalıcılık maddesi işler. Bant ikamesi /
onarımı YASAK.

══════════════════════════════════════════════════════════════════
5. MÜHÜR ÖNCESİ KAPILAR (bu betiğin kendisi koşar)
══════════════════════════════════════════════════════════════════
K-KİMLİK (Ş1.7): üç çapada da |g_E^yeni − g_E^disk| ≤ 1e−12 ve
    |g_X^yeni − g_X^disk| ≤ 1e−12; ayrıca θ^eski yeniden üretilir.
    Tutmazsa ölçüm makinesi aynı makine değildir ⇒ ENGELLENDİ.
ÖLÇÜT-T (Ş3.4, TEK ölçüt, SIFIR serbest parametre): ikizin beş
    g_X,b'sinin hepsi sonlu, paydası ⟨X_b,X_b⟩ > 0 ve işareti
    ikizin donmuş global g_X'inin işaretiyle AYNI olmalı. Değilse
    bant-ağırlıklama anlamsızdır ⇒ vekil verisine HİÇ BAKILMADAN,
    MÜHÜRLENMEDEN ENGELLENDİ. (Sayısal bir yayılım eşiği KASITLI
    OLARAK KONMAMIŞTIR: gerekçesi yalnız 178'in içinde doğacak bir
    sabit Ş3.5'i ihlâl ederdi. Ağırlıklı ki-kare ve max/min yayılımı
    TANI olarak raporlanır, kapı DEĞİLDİR.)
G-GÜÇ (Ş2.7): donmuş eşiklerle H1(a), H1(b), H2-üst, H2-alt için
    KESİN kümelerin BOŞ OLMADIĞI, ölçüm içermeyen cebirsel bir
    tanıkla mühürden önce gösterilir.

══════════════════════════════════════════════════════════════════
6. YENİDEN AÇILMAYANLAR (Ş0.3)
══════════════════════════════════════════════════════════════════
F2 (g_E ölümü), F4 (M ölümü), F6 (T-3), F8 (H-F2), R_η hükmü,
G1–G5 inşa kapıları — yeniden yargılanmaz. R_bant raporlanır,
HİÇBİR HÜKME DAYANAK YAPILMAZ (ne lehte ne aleyhte). Git'e
dokunulmaz. Ölümler kurtarmasız.

Kullanım: 178a_onkayit.py     (argümansız; ikinci koşuda yazmaz)
"""
import hashlib
import importlib
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
for _p in ("167_configs", "166_configs", "165_configs", "163_configs",
           "160_configs", "159_configs", "155_configs", "154_configs"):
    sys.path.insert(0, str(QM / _p))
K = importlib.import_module("165_cekirdek")
ORT = importlib.import_module("167_ortak")
K.PENCERE.update(ORT.pencere_dict())

SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S167 = SCR / "167"
S178 = SCR / "178"
OUT = S178 / "ONKAYIT.json"
BU = Path(__file__).resolve()

# --- DONMUŞ SABİTLER --------------------------------------------------
TABAN, CAP, TAU_MAX, TAU_C = 0.40, 4000, 0.95, 0.95
NJACK = 8
BANT_LO = (0.52, 0.56, 0.60, 0.64, 0.68)     # 176a/F0 hüküm penceresi
BANT_GEN = 0.04                              # IZGARA_T1 genişliği
LO_M = 0.60                                  # M'nin bandı (172b.oku)
IKIZ, GERCEK, ERFC = "Hkeskin", "son", "HA4"
CAPALAR = (IKIZ, GERCEK, ERFC)
VEKILLER = ("VF1", "VF2", "VF3", "VF4")      # YASAK LİSTE (körlük)
N_TAVAN = 4
KIMLIK_TOL = 1e-12

# --- ESKİ PARA BİRİMİ (dürüstlük beyanı; hükme GİREMEZ) ---------------
ESKI = dict(
    theta=dict(Hkeskin=0.8933338132966645, son=0.9219381611740145,
               HA4=0.7235729225877114, VF1=0.8666712252245586,
               VF2=0.8321358400382856, VF3=0.8797129149135897,
               VF4=0.8917115274642419),
    gX_vekil=dict(VF1=0.8134744638, VF2=0.8333309067,
                  VF3=0.8154829076, VF4=0.8068089541),
    esik_a=0.8656169530534596, esik_b=0.8933338132966645,
    kilit_theta=0.036231144031667745, f_theta=-0.14954441057282633,
    sac4=0.025733574899454686, ort_theta=0.867557876910169,
    omega_theta=1.2235332975434974,
    hukum_177="H-F1b KALICI HÜKÜMSÜZ (eski kestirici; geri alınmaz)")
ONCEKI_ONKAYIT = dict(
    onkayit_176="6173943f0f4ef0488dced150f9e5a767e13c6bb24af2837527d2"
                "bd3149690244",
    onkayit_177="03a9b9addedf738b4f4b03f0dedd0d54848fd11764cb4100e9490"
                "630cd12aebe")


# ---------------------------------------------------------------------
def _bekci(ad):
    """KÖRLÜK BEKÇİSİ — vekil adı geçen hiçbir şey açılamaz (Ş2.1)."""
    s = str(ad)
    for v in VEKILLER:
        if v in s:
            raise SystemExit(f"KÖRLÜK İHLÂLİ: vekil dosyası açılamaz -> {s}")
    return s


def _ac(p):
    return json.load(open(_bekci(p)))


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def _jk(v):
    """166_T1._jk — birebir."""
    v = np.asarray([x for x in v if np.isfinite(x)], float)
    if len(v) < 4:
        return float("nan")
    return float(np.sqrt((len(v) - 1) / len(v) * np.sum((v - v.mean()) ** 2)))


def _egim(Xb, x1):
    """⟨X_b·x1⟩/⟨X_b·X_b⟩ — 165.alanlar'ın gX biçimi, banda kısıtlı."""
    Xb = Xb - Xb.mean()
    pay, payda = float(np.dot(Xb, x1)), float(np.dot(Xb, Xb))
    return (pay / payda if payda > 0 else float("nan")), pay, payda


# ---------------------------------------------------------------------
def alanlar_ve_bantlar(ad):
    """167/165 zincirini AYNEN koşar; bant alanlarını ayrıca sentezler."""
    _bekci(ad)
    t0 = time.time()
    Y = K.gaz(ad, TABAN, CAP)
    Mo = K.Model165(ad, TABAN, CAP, TAU_MAX, "olculen", Y=Y, tau_c=TAU_C)
    Mo.sec("olculen", TAU_C).alanlar(kmax=3)
    A = Mo.artik
    x1 = Y.Xtil0
    tau, w, q = Mo.M["tau"], Mo.M["w"], Mo.M["q"]

    B = []
    for lo in BANT_LO:
        hi = lo + BANT_GEN
        mb = Mo.msk & (tau > lo) & (tau <= hi + 1e-12)
        idx = np.flatnonzero(mb)
        idx = idx[np.argsort(q[idx], kind="stable")]        # q artan
        grup = np.arange(len(idx)) % NJACK                  # round-robin
        Xg = [K.sentez(Mo.s, w[idx[grup == j]], Mo.y[idx[grup == j]])
              for j in range(NJACK)]
        Xb = np.sum(Xg, axis=0)
        g_b, pay, payda = _egim(Xb, x1)
        jk = [_egim(Xb - Xg[j], x1)[0] for j in range(NJACK)]
        B.append(dict(lo=float(lo), hi=float(hi), n=int(len(idx)),
                      gX_b=g_b, pay=pay, payda=payda,
                      jk=[float(v) for v in jk], s_b=_jk(jk)))
    d = _ac(S167 / f"C_{ad}.json")
    orta = [b for b in d["bant"] if abs(b["lo"] - LO_M) < 1e-9][0]
    return dict(ad=ad, gE=A["gE"], gX_global=A["gX"], gcal=A["gcal"],
                M=orta["KALIB_u2"], nline=A["nline"],
                gE_disk=d["artik"]["gE"], gX_disk=d["artik"]["gX"],
                th_eski=orta["KALIB_u2"] / d["artik"]["gcal"],
                bant=B, sure_s=time.time() - t0)


def birlestir(B, W):
    return float(sum(w * b["gX_b"] for w, b in zip(W, B)))


def theta_bw(R, W):
    g = birlestir(R["bant"], W)
    return R["M"] / (R["gE"] * g * g), g


# ---------------------------------------------------------------------
def main():
    S178.mkdir(parents=True, exist_ok=True)
    if OUT.exists():
        d = json.load(open(OUT))
        print("ÖN-KAYIT ZATEN VAR — ÜZERİNE YAZILMAZ.")
        print("  zaman : %s" % d["zaman"])
        print("  sha256: %s" % d["sha256"])
        return d

    print("=" * 74)
    print("178a — ÖN-KAYIT: BANT-AĞIRLIKLI θ (θ_bw), H-F1b^ba, k=2. bakış")
    print("=" * 74)
    print("BEYAN: bu bir VARYANS-KÜÇÜLTME KESTİRİCİ DEĞİŞİKLİĞİDİR;")
    print("       EŞİK GEVŞETME DEĞİLDİR. Eşik FORMÜLLERİ 176-F3 ile")
    print("       harf harf aynı; yalnız çapa sayıları yeni parada.")
    print("       177'nin eski-para hükmü (KALICI HÜKÜMSÜZ) ayakta kalır.")
    print("KÖRLÜK: vekil (VF1..VF4) dosyaları AÇILMAZ — bekçi etkin.\n")

    # ---- çapalar ----------------------------------------------------
    R = {}
    for g in CAPALAR:
        R[g] = alanlar_ve_bantlar(g)
        r = R[g]
        print("  %-8s g_E=%.10f  g_X(global)=%.10f  M=%.10f  "
              "çizgi=%d  (%.0f s)"
              % (g, r["gE"], r["gX_global"], r["M"], r["nline"],
                 r["sure_s"]))

    # ---- K-KİMLİK kapısı (Ş1.7) --------------------------------------
    print("\n  K-KİMLİK KAPISI (tol %.0e):" % KIMLIK_TOL)
    kim, kimlik_ok = {}, True
    for g in CAPALAR:
        r = R[g]
        dE = abs(r["gE"] - r["gE_disk"])
        dX = abs(r["gX_global"] - r["gX_disk"])
        dT = abs(r["th_eski"] - ESKI["theta"][g])
        ok = (dE <= KIMLIK_TOL) and (dX <= KIMLIK_TOL) and (dT <= 1e-12)
        kimlik_ok &= ok
        kim[g] = dict(d_gE=dE, d_gX=dX, d_theta_eski=dT, gecti=bool(ok))
        print("    %-8s |Δg_E|=%.3e  |Δg_X|=%.3e  |Δθ_eski|=%.3e   %s"
              % (g, dE, dX, dT, "✓" if ok else "✗ KAPI TUTMADI"))
    if not kimlik_ok:
        raise SystemExit("ENGELLENDİ — K-KİMLİK kapısı tutmadı; ölçüm "
                         "makinesi aynı makine değil. Mühür YOK.")

    # ---- ağırlıklar: YALNIZ İKİZDEN (Ş3.2) ---------------------------
    Bi = R[IKIZ]["bant"]
    inv = [1.0 / (b["s_b"] ** 2) for b in Bi]
    top = sum(inv)
    W = [v / top for v in inv]
    print("\n  AĞIRLIKLAR (kaynak: YALNIZ ikiz %s; njack=%d; floor YOK):"
          % (IKIZ, NJACK))
    for b, wv in zip(Bi, W):
        print("    bant τ∈(%.2f,%.2f]  n=%3d  g_X,b=%+.8f  s_b=%.3e  "
              "W_b=%.6f" % (b["lo"], b["hi"], b["n"], b["gX_b"],
                            b["s_b"], wv))
    print("    Σ W_b = %.15f" % sum(W))

    # ---- ÖLÇÜT-T kapısı (Ş3.4) ---------------------------------------
    isa = np.sign(R[IKIZ]["gX_global"])
    t_ok = all(np.isfinite(b["gX_b"]) and b["payda"] > 0
               and np.sign(b["gX_b"]) == isa for b in Bi)
    gb = np.array([b["gX_b"] for b in Bi], float)
    sb = np.array([b["s_b"] for b in Bi], float)
    gort = float(np.sum(gb / sb ** 2) / np.sum(1.0 / sb ** 2))
    ki2 = float(np.sum((gb - gort) ** 2 / sb ** 2))
    yay = float(gb.max() / gb.min()) if gb.min() > 0 else float("nan")
    print("\n  ÖLÇÜT-T (TEK kapı, sıfır serbest parametre): sonluluk + "
          "payda>0 + işaret(g_X global) ⇒ %s" % ("✓" if t_ok else "✗"))
    print("    TANI (kapı DEĞİL): ağırlıklı ort=%.8f  χ²=%.4f  "
          "χ²/sd=%.4f (sd=%d)  yayılım max/min=%.4f"
          % (gort, ki2, ki2 / (len(gb) - 1), len(gb) - 1, yay))
    if not t_ok:
        raise SystemExit("ENGELLENDİ — ÖLÇÜT-T tutmadı; bant-ağırlıklama "
                         "anlamsız. Vekil verisine BAKILMADI. Mühür YOK.")

    # ---- yeni para biriminde çapalar ---------------------------------
    th, gbw, ind = {}, {}, {}
    for g in CAPALAR:
        th[g], gbw[g] = theta_bw(R[g], W)
        Bg = R[g]["bant"]
        pay = sum(b["pay"] for b in Bg)
        payda = sum(b["payda"] for b in Bg)
        ind[g] = dict(gX_global_donmus=R[g]["gX_global"],
                      gX_S_guc=pay / payda, gX_bw=gbw[g],
                      ozdeslik_kalinti=abs(
                          th[g] - R[g]["M"] / (R[g]["gE"] * gbw[g] ** 2)))
    print("\n  ÇAPALAR — YENİ PARA BİRİMİ:")
    for g in CAPALAR:
        print("    %-8s g_X,bw=%.10f   θ_bw=%.10f   (eski θ=%.10f)"
              % (g, gbw[g], th[g], ESKI["theta"][g]))
    print("\n  İNDİRGEME BELGESİ (Ş3.1) — g_X üç okuma:")
    for g in CAPALAR:
        i = ind[g]
        print("    %-8s global=%.8f  güç-ağırlıklı(S)=%.8f  bw=%.8f  "
              "|özdeşlik kalıntı|=%.2e"
              % (g, i["gX_global_donmus"], i["gX_S_guc"], i["gX_bw"],
                 i["ozdeslik_kalinti"]))

    # ---- eşikler: 176'nın AYNI formülleri (Ş1.2, Ş1.3) ---------------
    D = float(np.log(th[GERCEK]) - np.log(th[IKIZ]))
    E = float(np.log(th[ERFC]) - np.log(th[IKIZ]))
    esik_a = th[IKIZ] * float(np.exp(-D))
    esik_a2 = th[IKIZ] ** 2 / th[GERCEK]
    esik_b = th[IKIZ]
    f_th = D / E
    kesim = f_th * D
    kilit = (1.0 - f_th) * D
    P1 = bool(D > 0)
    P2 = bool(kilit > 0 and abs(E) > 1e-12)
    print("\n  EŞİKLER (176-F3 formülleri, yeni çapalar):")
    print("    D  = Δlog θ_bw(son←Hk)  = %+.15f   ⇒ P1 %s"
          % (D, "✓" if P1 else "✗ SORU DÜŞTÜ"))
    print("    E  = Δlog θ_bw(erfc←Hk) = %+.15f" % E)
    print("    eşik (a) = θ_bw(Hk)·e^{−D} = %.15f" % esik_a)
    print("               (θ_bw(Hk)²/θ_bw(son) = %.15f ; fark %.2e)"
          % (esik_a2, abs(esik_a - esik_a2)))
    print("    eşik (b) = θ_bw(Hk)        = %.15f" % esik_b)
    print("    f_θ = D/E = %+.15f   kesim_θ = %+.15f   kilit_θ = %+.15f"
          % (f_th, kesim, kilit))
    print("    P2 (kilit_θ > 0) %s" % ("✓" if P2 else "✗ H2 TANIMSIZ"))
    print("    SAĞLAMA (eski para): %.10f²/%.10f = %.13f  (176-F3: %.13f)"
          % (ESKI["theta"]["Hkeskin"], ESKI["theta"]["son"],
             ESKI["theta"]["Hkeskin"] ** 2 / ESKI["theta"]["son"],
             ESKI["esik_a"]))

    # ---- G-GÜÇ: KESİN kümeler boş değil (Ş2.7) -----------------------
    def sac(v):
        return 2.0 * float(np.std(np.asarray(v, float), ddof=1)) / \
            np.sqrt(len(v))

    ta = [esik_a * (0.90 + 0.001 * k) for k in range(4)]  # dal (a) tanığı
    tb = [esik_b * (1.10 + 0.001 * k) for k in range(4)]  # dal (b) tanığı
    guc = dict(
        dal_a=dict(tanik=ta, ort=float(np.mean(ta)), sac4=sac(ta),
                   kesin=bool(np.mean(ta) <= esik_a
                              and abs(np.mean(ta) - esik_a) >= sac(ta))),
        dal_b=dict(tanik=tb, ort=float(np.mean(tb)), sac4=sac(tb),
                   kesin=bool(np.mean(tb) >= esik_b
                              and abs(np.mean(tb) - esik_b) >= sac(tb))))
    if P2:
        du = [kilit * (2.0 + 0.02 * k) for k in range(4)]  # H2 üst kenar
        dl = [kilit * (0.10 + 0.001 * k) for k in range(4)]  # H2 alt kenar
        guc["H2_ust"] = dict(tanik=du, ort=float(np.mean(du)),
                             sac4=sac(du),
                             kesin=bool(np.mean(du) >= kilit and
                                        abs(np.mean(du) - kilit) >= sac(du)))
        guc["H2_alt"] = dict(tanik=dl, ort=float(np.mean(dl)),
                             sac4=sac(dl),
                             kesin=bool(np.mean(dl) > 0 and
                                        abs(np.mean(dl)) >= sac(dl)))
    guc_ok = all(v["kesin"] for v in guc.values())
    print("\n  G-GÜÇ (Ş2.7) — KESİN kümeler boş değil (cebirsel tanık):")
    for k, v in guc.items():
        print("    %-7s ort=%+.8f  SAÇ_4=%.8f  KESİN %s"
              % (k, v["ort"], v["sac4"], "✓" if v["kesin"] else "✗"))
    if not guc_ok:
        raise SystemExit("ENGELLENDİ — KESİN küme boş; kurgulanmış "
                         "hükümsüzlük. Mühür YOK.")

    # ---- ön-kayıt kaydı ---------------------------------------------
    tarih = subprocess.run(["date"], capture_output=True,
                           text=True).stdout.strip()
    bilesim = {}
    for h1, c1 in (("YASADI", "H-F1b^ba YAŞADI"),
                   ("OLDU", "H-F1b^ba ÖLDÜ"),
                   ("HUKUMSUZ", "H-F1b^ba KALICI HÜKÜMSÜZ (2. para "
                                "biriminde de; 3. kestirici DENENMEZ)")):
        for h2, c2 in (("ODENEBILIR", "F9-θ satırı ✓ fatura ödenebilir "
                                      "— KESİN"),
                       ("ODENEMEZ", "F9-θ satırı ✗ fatura ödenemez "
                                    "— KESİN"),
                       ("HUKUMSUZ", "F9-θ satırı hükümsüz")):
            m = (h1 == "YASADI" and h2 == "ODENEBILIR")
            bilesim[f"{h1}|{h2}"] = dict(
                cumle=f"{c1}; {c2}.", muhur_adayi=bool(m))

    rec = dict(
        zaman=time.strftime("%Y-%m-%d %H:%M:%S %z"), date_cikti=tarih,
        betik=BU.name, betik_yolu=str(BU), sha256=sha(BU),
        hipotez="H-F1b^ba — θ'nın kilit-duyarlılığı, bant-ağırlıklı g_X "
                "para biriminde",
        bakis_sayaci=2,
        beyan="BU BİR VARYANS-KÜÇÜLTME KESTİRİCİ DEĞİŞİKLİĞİDİR; EŞİK "
              "GEVŞETME DEĞİLDİR. Değişen tek bileşen θ'nın bölenindeki "
              "g_X'in okunuşudur. M (KALİB_u2, lo=0.60) ve g_E "
              "değişmez. Eşik FORMÜLLERİ 176-F3 ile harf harf aynıdır; "
              "yalnız içlerine giren çapalar yeni para biriminde "
              "yeniden hesaplanmıştır. 177'nin eski-para hükmü "
              "(H-F1b KALICI HÜKÜMSÜZ) ayakta kalır, geri alınmaz.",
        onceki_onkayitlar=ONCEKI_ONKAYIT,
        kestirici=dict(
            ad="theta_bw",
            formul="theta_bw := M / (g_E * g_X_bw**2) ;  "
                   "g_X_bw := SUM_b W_b * g_X_b ;  "
                   "g_X_b := <X_b . x1> / <X_b . X_b> ;  "
                   "X_b := sentez(s, w[lo<tau<=hi], y[lo<tau<=hi]) - mean ;"
                   "  x1 := Y.Xtil0 ;  M := KALIB_u2(lo=0.60) "
                   "(DEĞİŞMEDİ) ;  g_E := artik.gE (DEĞİŞMEDİ)",
            olcum_zinciri="165_cekirdek.gaz(ad,0.40,4000) -> "
                          "Model165(ad,0.40,4000,0.95,'olculen',"
                          "tau_c=0.95).sec('olculen',0.95).alanlar(kmax=3)"
                          " — 176c/167 zinciriyle BİREBİR",
            bant_kumesi=dict(lo=list(BANT_LO), genislik=BANT_GEN,
                             tau_araligi="(0.52, 0.72]",
                             uyelik="lo < tau_q <= hi (163.bant_adaylari "
                                    "kuralı)",
                             kaynak="176a/F0 hüküm penceresi (176c "
                                    "LO_LO=0.52, LO_HI=0.68) + "
                                    "167_olcum.LO_MIN=0.52; M'nin bandı "
                                    "lo=0.60 pencerenin merkezi",
                             ikinci_aday_denenmedi="dokuz bandın tamamı "
                                                   "(lo 0.44..0.76) "
                                                   "DENENMEDİ — Ş3.3(ii)"),
            agirlik=dict(
                bicim="W_b ∝ 1/s_b², normalize Σ W_b = 1",
                kaynak=f"YALNIZ ikiz {IKIZ}; vekillerin hiçbir sayısı ve "
                       "tohumlar-arası varyans GİRMEZ",
                jackknife=f"bant içi silme-1-grup, njack={NJACK}, grup(i)="
                          "i mod 8 (q artan sıra; 163 round-robin şeması), "
                          "s_b = 166_T1._jk biçimi",
                floor="YOK", kirpma="YOK", uzay="doğrusal (log DEĞİL)",
                W=W, s_b=[b["s_b"] for b in Bi],
                gX_b_ikiz=[b["gX_b"] for b in Bi],
                n_cizgi=[b["n"] for b in Bi]),
            indirgeme="güç ağırlıkları ⟨X_b,X_b⟩ ile birleşim "
                      "Σ⟨X_b,x1⟩/Σ⟨X_b,X_b⟩ = bant desteği üzerindeki "
                      "köşegen-Gram EKK eğimi; Gram çaprazları sıfır ve "
                      "destek tüm ızgara olsaydı donmuş global g_X'i "
                      "aynen verirdi (sayılar 'indirgeme' bloğunda)"),
        capalar=dict(
            theta_bw={g: th[g] for g in CAPALAR},
            gX_bw={g: gbw[g] for g in CAPALAR},
            gE={g: R[g]["gE"] for g in CAPALAR},
            M={g: R[g]["M"] for g in CAPALAR},
            gX_global_donmus={g: R[g]["gX_global"] for g in CAPALAR},
            nline={g: R[g]["nline"] for g in CAPALAR},
            bant_defteri={g: R[g]["bant"] for g in CAPALAR}),
        indirgeme=ind,
        esikler=dict(
            D_dlog_theta_son_Hk=D, E_dlog_theta_erfc_Hk=E,
            esik_a=esik_a, esik_a_ikinci_yol=esik_a2, esik_b=esik_b,
            f_theta=f_th, kesim_theta=kesim, kilit_theta=kilit,
            formul_a="esik_a := theta_bw(Hk) * exp(-D) == "
                     "theta_bw(Hk)**2 / theta_bw(son)",
            formul_b="esik_b := theta_bw(Hk)",
            formul_kilit="f = D/E ; kesim = f*D ; kilit = (1-f)*D",
            saglama_eski_para=dict(
                esik_a_eski_formulden=ESKI["theta"]["Hkeskin"] ** 2
                / ESKI["theta"]["son"], esik_a_176=ESKI["esik_a"])),
        onkosullar=dict(P1_D_pozitif=P1, P2_kilit_pozitif=P2,
                        P1_aksi="H-F1b^ba TANIMSIZ — SORU DÜŞTÜ",
                        P2_aksi="H2 TANIMSIZ; yalnız H1 koşulur"),
        kural=dict(
            n=N_TAVAN, n_tavan="n=4 ÜST SINIRDIR — beşinci tohum yok, "
                               "merdiven yok, yeni gaz koşusu yok",
            sac="SAÇ_n := 2*std(y, ddof=1)/sqrt(n) ; n=4 ⇒ SAÇ_4 = "
                "std(y, ddof=1)",
            sac_kati="SAÇ^kat = max−min — BAĞLAYICI DEĞİL, raporlanır",
            tohum_kurali="VF1..VF4 DÖRDÜ DE ortalamaya girer; eleme / "
                         "ağırlıklama / dışlama YASAK",
            H1_a="ȳ ≤ esik_a VE |ȳ − esik_a| ≥ SAÇ_4 ⇒ H-F1b^ba YAŞADI",
            H1_b="ȳ ≥ esik_b VE |ȳ − esik_b| ≥ SAÇ_4 ⇒ H-F1b^ba ÖLDÜ",
            H1_aksi="HÜKÜMSÜZ ⇒ (n=4 tavanı + Ş2.6) KALICI HÜKÜMSÜZ",
            H2="Δ_i = log θ_bw(Hk) − log θ_bw(VF_i) ; ω_θ = kilit_θ/Δ̄ ; "
               "üst kenar Δ̄ ≥ kilit_θ (ω≤1), alt kenar Δ̄ > 0 (ω>0); "
               "her kenarda KESİNLİK ⟺ marj ≥ SAÇ_4(Δ)",
            kapi_olumu="bir vekilde payda ≤ 0 / θ_bw tanımsız ⇒ vekil "
                       "DÜŞÜRÜLEMEZ; H1 ve H2 HÜKÜMSÜZ + kalıcılık",
            ucuncu_deneme="Ş2.6 — bu dört tohumda ÜÇÜNCÜ kestirici "
                          "DENENMEZ; soru ancak yeni veriyle açılır",
            yeniden_acilmayanlar=["F2", "F4", "F6", "F8", "R_eta",
                                  "G1-G5"],
            R_bant="raporlanır; HİÇBİR HÜKME DAYANAK YAPILMAZ"),
        bilesim_tablosu=bilesim,
        kapilar=dict(
            K_KIMLIK=dict(tol=KIMLIK_TOL, gecti=bool(kimlik_ok),
                          ayrinti=kim),
            OLCUT_T=dict(
                tanim="ikizin beş g_X,b'si sonlu + payda>0 + işaret = "
                      "işaret(g_X global); sayısal yayılım eşiği "
                      "KASITLI OLARAK YOK (Ş3.5)",
                gecti=bool(t_ok), agirlikli_ort=gort, ki2=ki2,
                ki2_sd=ki2 / (len(gb) - 1), sd=len(gb) - 1,
                yayilim_max_min=yay),
            G_GUC=dict(gecti=bool(guc_ok), tanik=guc)),
        durustluk=dict(
            korluk_iddiasi="YOK — 176/177'nin bütün vekil sayıları "
                           "okunmuştur; aşağıda listelidir",
            eski_para=ESKI,
            guvence="ağırlıklar YALNIZ ikizden gelir ve bu ön-kayıtta "
                    "SAYI olarak donmuştur; vekil dosyaları bu betikte "
                    "AÇILMAMIŞTIR (bekçi etkin)",
            vekil_dosyasi_acilmadi=True,
            git="GİT'E DOKUNULMADI",
            bosa_giden_kosu="YOK — bu betik tek koşuda tamamlandı",
            sira="mühür ANINDA çapalar hesaplandı; vekil θ_bw'leri "
                 "mühürden SONRA hesaplanacak"),
        veri_yollari=dict(
            olcum_onbellegi=[str(S167 / f"C_{g}.json") for g in CAPALAR],
            tayf_onbellegi=[str(SCR / "165" /
                                f"tayf_{g}_t{TABAN}_tm{TAU_MAX}.npz")
                            for g in CAPALAR],
            kod=[str(QM / "165_configs" / "165_cekirdek.py"),
                 str(QM / "167_configs" / "167_olcum.py"),
                 str(QM / "172_configs" / "172b_gE_yasasi.py"),
                 str(QM / "176_configs" / "176c_olcum.py")]))

    OUT.write_text(json.dumps(rec, indent=1, ensure_ascii=False))
    print("\n" + "=" * 74)
    print("  MÜHÜR")
    print("=" * 74)
    print("  zaman : %s" % rec["zaman"])
    print("  date  : %s" % tarih)
    print("  sha256: %s" % rec["sha256"])
    print("  -> %s" % OUT)
    print("\n  ÖZET:")
    print("    θ_bw(ikiz  %s) = %.10f" % (IKIZ, th[IKIZ]))
    print("    θ_bw(gerçek %s) = %.10f" % (GERCEK, th[GERCEK]))
    print("    θ_bw(erfc  %s)  = %.10f" % (ERFC, th[ERFC]))
    print("    eşik (a) = %.10f   eşik (b) = %.10f" % (esik_a, esik_b))
    print("    kilit_θ  = %+.10f  f_θ = %+.10f" % (kilit, f_th))
    print("    W_b = [%s]" % ", ".join("%.6f" % v for v in W))
    print("    KURAL: SAÇ_4 = 2s/√4 = s ; n=4 TAVAN ; dört tohum da "
          "girer ; aksi ⇒ KALICI HÜKÜMSÜZ")
    return rec


if __name__ == "__main__":
    main()
