# -*- coding: utf-8 -*-
"""
182a — K0: DONMUŞ ÖN-KAYIT (GAUSS-ALTININ ANATOMİSİ + VF AİLESİ n=8)
=====================================================================
ÖN-MÜHÜR. Bu betik HİÇBİR yeni gaz kurulmadan, HİÇBİR yeni ölçüm
yapılmadan ÖNCE koşar; kendi sha256'sını ve `date` damgasını
`182/ONKAYIT_182.json`'a yazar; dosya bir daha DEĞİŞTİRİLMEZ.

NE DONDURULUR
-------------
(0) DEVRALINAN, GEVŞETİLMEYEN EŞİKLER — 180a ve 181a'dan, birebir:
      · GAUSS3 = (2/π)^{3/2}, GAUSS1 = √(2/π);  ρ₃(Hkeskin)=0.5284,
        ρ₃(son)=0.5356;
      · 180e'nin ÇÖZÜNÜRLÜK maddesi: saçılım := max ρ₃ − min ρ₃ ,
        HÜKÜMSÜZ eşiği 0.5·|ρ₃(Hk) − GAUSS3| = 0.0102254562630361 ;
      · 181a'nın okunurluk ölçütü O1: |ZARF_N| ≥ 2·se(ZARF_N), k=2.0 ;
      · 176a/180a inşa kapıları G1..G5 (sha(A) bit-bit, φ≡0 sağlaması
        0.0/0.0, ilk-kök 300000/300000, maks|F| ≤ 1e−8 aşan 0,
        sıralılık TAM).
(1) VF5..VF8'in tanımı (176b, tohum 5,6,7,8) ve inşa kapıları.
(2) n=8 istatistiğinin KARAR KURALLARI: Gauss-altı hükmü, kırpma
    hakeminin yeniden değerlendirmesi, O1'in n=8 hâli.
(3) GAUSS REFERANSININ TÜRETİM VARSAYIMLARI — açıkça yazılır (H-G1'in
    sınav nesnesi budur).
(4) H-G1 / H-G2 / H-G3'ün ÖLÜM ve YAŞAMA koşulları, sayısal.
(5) ÖN-KAYITLI ÖNGÖRÜLER (mevcut n=4 verisinden, yeni veri gelmeden).

BU BETİKTE HESAPLANMAYAN (ve hesaplanmaması ŞART olan): VF5..VF8'in
hiçbir niceliği; n=8 ⟨ρ₃⟩, z_G, saçılım; κ_E; vekil surrogate'ları;
n=8 ZARF/KİLİT payları. Hepsi bu dosya diske yazıldıktan SONRA ölçülür.

Kullanım: 182a_onkayit.py
"""
import hashlib
import json
import math
import subprocess
import time
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S182, S180, S181, S169, S176 = (SCR / "182", SCR / "180", SCR / "181",
                                SCR / "169", SCR / "176")
BU = Path(__file__).resolve()
KALEM = QM / "KALEM_GAUSS_ALTI_08EYL2026.md"

G1 = math.sqrt(2 / math.pi)                 # 0.7978845608028654
G3 = G1 ** 3                                # 0.5079490874739278
HAM = ("M", "Q_E", "rho_E", "Q_X", "rho_X")
ORAN = ("g_E", "g_X", "theta", "g_cal")


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def istat(v):
    v = np.asarray(v, float)
    n = len(v)
    sd = float(np.std(v, ddof=1)) if n > 1 else 0.0
    return float(v.mean()), sd, sd / math.sqrt(n), n


# =====================================================================
# GAUSS REFERANSININ TÜRETİMİ — H-G1'İN SINAV NESNESİ
# =====================================================================
GAUSS_TURETIM = {
    "olculen_nesne": (
        "169_k2b: G = E·X_a·X_b üçlü çarpımının, hüküm bantlarındaki "
        "çizgi frekanslarında ölçülen spektral yanıtı P[varyant]; "
        "ρ(varyant) := P[varyant]/P[000_tam]. clip(F) := σ_F·sgn(F)."),
    "V1_marjinal": (
        "Her bacak F ayrı ayrı GAUSS'tur ⇒ ⟨clip(F)·F⟩/⟨F²⟩ = "
        "E|F|/σ_F = √(2/π) = 0.7978845608028654. "
        "ÖLÇÜLEBİLİR ADI: κ_F (169_k2b bunu zaten basıyor)."),
    "V2_ortak_gaussluk": (
        "E, X_a, X_b BİRLİKTE Gauss'tur ⇒ Price teoremi: herhangi bir "
        "G ile ⟨sgn(F)·G⟩ = √(2/π)·⟨F·G⟩/σ_F. Bacak kırpması çarpımın "
        "HER çapraz momentini aynı √(2/π) çarpanıyla ölçekler."),
    "V3_carpimsallik": (
        "Bacaklar birbirinden BAĞIMSIZ ⇒ k bacak kırpıldığında çarpan "
        "(√(2/π))^k; 3 bacak için GAUSS3 = (2/π)^{3/2} = "
        "0.5079490874739278."),
    "V4_bant_yansizligi": (
        "P'nin bant ağırlıklandırması (pN·A_on² − pO·A_off²)/Δ ve "
        "s_den_J zinciri, kırpmanın yol açtığı spektral yeniden "
        "dağılıma KÖR değildir; referans bu ağırlığın ρ'yu ölçeklemediğini "
        "varsayar."),
    "not": ("Gauss referansı bir ÖLÇÜM DEĞİL bir TÜRETİMDİR: V1∧V2∧V3∧V4. "
            "H-G1 = 'referansın kendisi bu gazda geçersiz'. Bu dört "
            "varsayımın her biri ayrı ayrı sınanır (182g)."),
}

# =====================================================================
# HİPOTEZLERİN ÖLÜM / YAŞAMA KOŞULLARI (sayısal, veriden önce)
# =====================================================================
HIP = {
  "H-G1": {
    "iddia": ("Gauss referansı (√(2/π))^k bu gazda geçersizdir — "
              "genlik dağılımı / basıklık varsayımı (V1) veya referansın "
              "kendisi (V2–V4) kırılıyor."),
    "olculen": [
        "κ_E, κ_Xa, κ_Xb  (= ⟨clip(F)F⟩/⟨F²⟩, 169_k2b'nin bastığı sayı)",
        "basıklık kurt(E), kurt(Xa), kurt(Xb) (fazla basıklık, Gauss=0)",
        "ρ°(E): ÜÇ bacak da bağımsız Gauss surrogate (spektrum korunur, "
        "fazlar bağımsız U(0,2π)) ile aynı kestirimciden geçirilir",
    ],
    "OLUM": ("(i) 8 tohumun HEPSİNDE |κ_E/√(2/π) − 1| ≤ 0.05 iken "
             "(√(2/π) − ρ(E))/√(2/π) ≥ 0.10  VE  "
             "(ii) |ρ°(E)/√(2/π) − 1| ≤ 0.02  (8 tohumun ortalamasında). "
             "Yani hem marjinal Gaussluk sağlanıyor hem de referans, "
             "KENDİ varsayımları altında kestirimciden doğru çıkıyor ⇒ "
             "açığı referansın türetimi taşıyamaz. KURTARMA YOK."),
    "YASAR": ("|ρ°(E)/√(2/π) − 1| > 0.05, ya da κ_E'nin sapması ρ(E) "
              "açığının en az yarısını (≥ 0.5·(√(2/π) − ρ(E))) tek başına "
              "açıklıyor."),
    "ARADA": "0.02 < |ρ°(E)/√(2/π) − 1| ≤ 0.05 ⇒ HÜKÜMSÜZ yazılır.",
  },
  "H-G2": {
    "iddia": ("İnşa öz-tutarlılığı: sadakatli zincir (ızgara braketi + "
              "SIRALI ilk-kök + n−½) fazlar rastgele olsa bile aralık↔alan "
              "ilişkisine yapı sokar; karıştırma TAM bağımsızlık üretmez."),
    "olculen": [
        "korr(X_a, X_b) — bölünmüş merdivenin iki YARISI; tam bağımsızlıkta "
        "beklenen 0",
        "R_η = P/Var(η) ve ÖZDEŞLİK Kov(η_artık, η)/Var(η) = 1 − R_η "
        "(174-K1d); işaret dağılımı 8 tohumda",
        "ρ*(E): YALNIZ E Gauss surrogate'la değiştirilir (X_a, X_b aynen); "
        "geri kazanım Rec := (ρ*(E) − ρ(E))/(√(2/π) − ρ(E))",
    ],
    "OLUM": ("(i) ⟨korr(X_a,X_b)⟩_{n=8} sıfırla uyumlu (|⟨r⟩| ≤ 2·se) "
             "VEYA işareti 8 tohumda tekdüze değil;  VEYA "
             "(ii) Rec ≤ 0.30 (E'nin faz yapısını silmek açığı "
             "kapatmıyor). KURTARMA YOK."),
    "YASAR": ("|⟨korr(X_a,X_b)⟩| ≥ 5·se ve işaret 8/8 tekdüze,  VE  "
              "Rec ≥ 0.70 (en az 7/8 tohumda)."),
    "ARADA": "0.30 < Rec < 0.70 ⇒ HÜKÜMSÜZ yazılır.",
  },
  "H-G3": {
    "iddia": ("Kırpma-zinciri artefaktı: ρ'nun düşüşü 169_k2b'nin "
              "ölçüm-ağırlıklı kestirimcisinden (pN·A_on² − pO·A_off², "
              "s_den_J, bant birleştirmesi) doğuyor."),
    "olculen": [
        "BAĞIMSIZ KESTİRİMCİ ρ_J(varyant) := "
        "Σ_k Re[J_var(W_k)·conj(J_000(W_k))] / Σ_k |J_000(W_k)|² — "
        "yalnız MODEL tarafı; ölçülen gazın hiçbir niceliği (pN, A, h, "
        "ρ_ort) girmez, bant ağırlığı yoktur, s_den_J kullanılmaz.",
        "aynı W_k kümesi (Wall) kullanılır — karşılaştırma adil.",
    ],
    "OLUM": ("8 tohumun HEPSİNDE |ρ_J(E)/ρ(E) − 1| ≤ 0.05 "
             "VE ρ_J(E)'nin Gauss'tan açığı da ≥ 0.10 (yani bağımsız "
             "kestirimci aynı açığı görüyor) ⇒ artefakt değil. "
             "KURTARMA YOK."),
    "YASAR": "⟨|ρ_J(E)/ρ(E) − 1|⟩ > 0.10.",
    "ARADA": "0.05 < ⟨|ρ_J/ρ − 1|⟩ ≤ 0.10 ⇒ HÜKÜMSÜZ yazılır.",
  },
}

# =====================================================================
# ÇİZGİ / BANT AYRIŞIMI (E-bacağının anatomisi) — tanım donuyor
# =====================================================================
AYRISIM_E = {
    "carpimsal_tanim": "ρ(E) ≡ κ_E · β_E ,  β_E := ρ(E)/κ_E",
    "CIZGI_payi": ("κ_E — alanın MARJİNAL biçim çarpanı E|E|/σ_E; "
                   "kırpmanın tek-alan (genlik dağılımı) etkisi."),
    "BANT_payi": ("β_E — geri kalan her şey: kırpmanın ortak/spektral "
                  "yeniden dağılımı, bant ağırlığı, çapraz momentler."),
    "hangisi_tasiyor": (
        "Gauss'ta κ=√(2/π), β=1. Açık log dilinde bölünür: "
        "log ρ(E) − log √(2/π) = [log κ_E − log √(2/π)] + log β_E. "
        "PAY_ÇİZGİ := (log κ_E − log √(2/π)) / (log ρ(E) − log √(2/π))."),
    "tau_egimi": ("Her gazda 5 hüküm bandında ρ(E) vs τ_eff doğrusal eğim "
                  "s_E; aynısı ρ(Xa), ρ(Xb) için. ÖN-KAYIT: "
                  "|⟨s_E⟩| ≥ ⟨|s_X|⟩ + 3·se(s_E) ise düşüş BANT "
                  "karakterlidir (spektral), ÇİZGİ değil."),
}

# =====================================================================
# n = 8 KARAR KURALLARI
# =====================================================================
KARAR = {
  "A_gauss_alti": {
    "istatistik": ("ρ₃(g) := 169_k2b.main(g,0.40,0.95)['ortalama']"
                   "['111_hepsi'], g ∈ VF1..VF8;  ⟨ρ₃⟩, sd(ddof=1), "
                   "se = sd/√8;  z_G := (⟨ρ₃⟩ − GAUSS3)/se."),
    "HUKUM_GAUSS_ALTI": ("z_G ≤ −5 VE 8 tohumun HEPSİ tek tek "
                         "ρ₃ < GAUSS3 ⇒ 'karışık vekil gaz GAUSS'UN "
                         "ALTINDADIR' yazılır (n=8 ile)."),
    "HUKUM_YOK": ("z_G > −5 ya da bir tohum bile ρ₃ ≥ GAUSS3 ⇒ "
                  "Gauss-altı hükmü VERİLMEZ; 180'in n=4 gözlemi "
                  "sağlamlaşmadı diye yazılır."),
    "not": ("Bu hüküm, 180e'nin DAL A/DAL B hakeminden AYRI bir "
            "iddiadır; hakemin çözünürlük maddesini ne gevşetir ne "
            "yerine geçer."),
  },
  "B_kirpma_hakemi": {
    "kural": ("180e'nin dondurulmuş kuralı AYNEN, yalnız gaz kümesi "
              "VF1..VF8. Eşikler: DAL_A C≥0.70 ∧ |z_G|≤3; DAL_B C≤0.30; "
              "saçılım = max−min > 0.0102254562630361 ⇒ HÜKÜMSÜZ."),
    "UREME_KAPISI": ("182e, aynı kodla VF1..VF4 koşulduğunda "
                     "180/K2_HAKEM.json'un dal, ⟨ρ₃⟩, sd, se, C, z_G, "
                     "saçılım alanlarını TAM SIFIR farkla üretmeli; "
                     "üretmezse betik ÖLÜR."),
    "ONGORU_P1": ("Saçılım istatistiği RANGE'dir ve n'de MONOTON "
                  "ARTANDIR. VF1..VF4'te 0.016091 ölçüldü; VF5..VF8 "
                  "eklenince saçılım ≥ 0.016091 KESİNDİR ⇒ "
                  "0.0102254562630361 eşiği kaçınılmaz olarak aşılır ⇒ "
                  "HÜKÜMSÜZ hükmü n=8'de de KURTULAMAZ. Bu, ölçümden "
                  "önce yazılmış YAPISAL bir öngörüdür: 180'in "
                  "çözünürlük maddesi tohum EKLEYEREK aşılamaz; "
                  "aşılabilseydi ancak tohumların YAYILIMI küçülerek "
                  "aşılırdı, oysa range küçülemez. EŞİK GEVŞETİLMEZ."),
    "EK_KAYIT": ("Karar vermeyen, yalnız kaydedilen ölçek özetleri: "
                 "sd, se, IQR, ve 'saçılım/√n' — hiçbiri hükmü "
                 "değiştirmez."),
  },
  "C_okunurluk_O1": {
    "kural": "181a'nın O1'i AYNEN: |ZARF_N| ≥ 2·se(ZARF_N), k_sigma = 2.0",
    "makine": ("181c_defter_n4.ayrisim İTHAL edilir (kopyalanmaz); "
               "VS = VS1..VS4 (n=4, değişmez), VF = VF1..VF8 (n=8)."),
    "UREME_KAPISI": ("182h, VF=VF1..VF4 ile koşulduğunda "
                     "181/K3_DEFTER_n4.json'un her alanını TAM SIFIR "
                     "farkla üretmeli; üretmezse betik ÖLÜR."),
    "hukum": ("Bir satır n=8'de |Z|/se ≥ 2.0 ise OKUNUR, değilse "
              "HÜKÜMSÜZ. Yeni okunan satırlar 180 defterine işlenir; "
              "okunmayanlar süslenmeden HÜKÜMSÜZ kalır."),
  },
}


def main():
    t0 = time.time()
    S182.mkdir(parents=True, exist_ok=True)
    dmg = subprocess.run(["date"], capture_output=True, text=True).stdout.strip()
    print("=" * 78, flush=True)
    print("182a — K0: DONMUŞ ÖN-KAYIT (GAUSS-ALTININ ANATOMİSİ, n=8)",
          flush=True)
    print(f"  zaman (`date`): {dmg}", flush=True)
    print("=" * 78, flush=True)

    OK180 = json.load(open(S180 / "ONKAYIT_K0.json"))
    OK181 = json.load(open(S181 / "ONKAYIT_181_K0.json"))
    HK = json.load(open(S180 / "K2_HAKEM.json"))
    D4 = json.load(open(S181 / "K3_DEFTER_n4.json"))
    RHK = OK180["H180a"]["capalar"]["rho3_Hkeskin"]
    RSON = OK180["H180a"]["capalar"]["rho3_son"]
    # Çapa DEVRALINIR, yeniden türetilmez: 180a'nın dondurduğu bit-bit
    # değer kullanılır (G1**3 son bitte 1 ULP farklı çıkıyor; 180'in
    # sayısı esastır).
    G3f = OK180["H180a"]["capalar"]["GAUSS3"]
    G1f = OK180["H180a"]["capalar"]["GAUSS1"]
    ESIK = 0.5 * abs(RHK - G3f)

    print(f"\n  DEVRALINAN EŞİKLER (gevşetme YOK)")
    print(f"    GAUSS1 = {G1f:.16f}   GAUSS3 = {G3f:.16f}   "
          f"(180a'dan DEVRALINDI; G1**3 = {G1**3:.16f}, 1 ULP fark)")
    print(f"    ρ₃(Hkeskin) = {RHK}   ρ₃(son) = {RSON}")
    print(f"    180e çözünürlük eşiği = 0.5·|ρ₃(Hk)−GAUSS3| = "
          f"{ESIK:.16f}")
    print(f"    (180'in kaydettiği eşik = {HK['esik_sacilim']:.16f}  "
          f"fark {abs(ESIK-HK['esik_sacilim']):.1e})")
    print(f"    181a O1: |ZARF| ≥ {OK181['okunurluk']['k_sigma']}·se(ZARF)")
    assert abs(ESIK - HK["esik_sacilim"]) == 0.0, "EŞİK YENİDEN ÜRETİLEMEDİ"
    assert OK181["okunurluk"]["k_sigma"] == 2.0

    # ---- ÖN-KAYITLI ÖNGÖRÜ P1: saçılım monotonluğu -------------------
    r3_4 = HK["rho3"]
    sac4 = HK["sacilim"]
    print(f"\n  ÖNGÖRÜ P1 (yapısal, kaçınılmaz):")
    print(f"    n=4 saçılım = {sac4:.6f} > eşik {ESIK:.6f}")
    print(f"    range monotondur ⇒ n=8 saçılım ≥ {sac4:.6f} ⇒ "
          f"HÜKÜMSÜZ hükmü DEĞİŞMEZ (eşik gevşetilmez)")
    print(f"    n=4 tek tek ρ₃: " +
          " ".join(f"{k}={v:.4f}" for k, v in r3_4.items()))

    # ---- ÖN-KAYITLI ÖNGÖRÜ P2: O1'in n=8 hâli ------------------------
    print(f"\n  ÖNGÖRÜ P2 (O1, n_VF: 4 → 8; varsayım: sd_VF ve merkezler "
          f"sabit)")
    print("  satır      konv.     |Z|/se(n=4)  se_VF(4)  se_VF(8)ö  "
          "se_ZARF(8)ö  |Z|/se(8)ö  hüküm(4)  ÖNGÖRÜ(8)")
    P2 = {}
    for k in HAM + ORAN:
        for konv in ("dogrusal", "log"):
            r = D4[konv][k]
            seVF8 = r["VF_sd"] / math.sqrt(8.0)
            seZ8 = math.sqrt(r["VS_se"] ** 2 + seVF8 ** 2)
            t4 = abs(r["ZARF"]) / r["se_ZARF"]
            t8 = abs(r["ZARF"]) / seZ8
            h4 = "OKUNUR" if t4 >= 2.0 else "HÜKÜMSÜZ"
            h8 = "OKUNUR" if t8 >= 2.0 else "HÜKÜMSÜZ"
            # yapısal taban: se_ZARF ≥ se_VS her zaman ⇒ ulaşılabilirlik
            ulas = bool(abs(r["ZARF"]) >= 2.0 * r["VS_se"])
            print(f"  {k:10s} {konv:9s} {t4:10.3f} {r['VF_se']:10.6f} "
                  f"{seVF8:10.6f} {seZ8:11.6f} {t8:11.3f}  {h4:9s} {h8}"
                  f"{'' if ulas else '   [VF tohumuyla ULAŞILAMAZ: |Z|<2·se_VS]'}",
                  flush=True)
            P2[f"{konv}/{k}"] = dict(
                t_n4=t4, hukum_n4=h4, se_VF_n4=r["VF_se"], VF_sd_n4=r["VF_sd"],
                se_VF_n8_ongoru=seVF8, se_ZARF_n8_ongoru=seZ8,
                t_n8_ongoru=t8, ongorulen_hukum=h8, ZARF_n4=r["ZARF"],
                se_VS_n4=r["VS_se"], ulasilabilir_VF_ile=ulas)

    yeni = [k for k, v in P2.items()
            if v["hukum_n4"] == "HÜKÜMSÜZ" and v["ongorulen_hukum"] == "OKUNUR"]
    ulasilamaz = sorted({k.split("/")[1] for k, v in P2.items()
                         if not v["ulasilabilir_VF_ile"]})
    print(f"\n    ÖNGÖRÜLEN YENİ OKUNUR SATIRLAR: "
          f"{yeni if yeni else 'YOK'}")
    print(f"    VF tohumu EKLEYEREK ULAŞILAMAZ (|ZARF| < 2·se_VS): "
          f"{ulasilamaz if ulasilamaz else 'YOK'}")
    print("    (yapısal: se_ZARF ≥ se_VS her zaman; VF tarafı sıfıra "
          "gitse bile bu satırlar açılmaz)")

    # ---- kayıt --------------------------------------------------------
    rec = dict(
        gorev=182, ad="GAUSS-ALTININ ANATOMİSİ + VF AİLESİ BÜYÜTMESİ",
        zaman=dmg, zaman_unix=time.time(), betik=BU.name, sha256=sha(BU),
        kalem=KALEM.name, sha_kalem=sha(KALEM),
        sha_176b=sha(QM / "176_configs" / "176b_vekil_insa.py"),
        sha_176c=sha(QM / "176_configs" / "176c_olcum.py"),
        sha_180d=sha(QM / "180_configs" / "180d_kirpma_kos.py"),
        sha_180e=sha(QM / "180_configs" / "180e_hakem.py"),
        sha_181c=sha(QM / "181_configs" / "181c_defter_n4.py"),
        sha_169k2b=sha(QM / "169_configs" / "169_k2b.py"),
        onkayit_180_sha=OK180["sha256"], onkayit_181_sha=OK181["sha256"],
        capalar=dict(GAUSS1=G1f, GAUSS2=G1f ** 2, GAUSS3=G3f,
                     rho3_Hkeskin=RHK, rho3_son=RSON,
                     rho3_VF_n4=r3_4, sacilim_n4=sac4,
                     esik_sacilim=ESIK, k_sigma=2.0),
        yeni_gazlar=[dict(ad=f"VF{i}", tohum=i, tip="vekil_faz",
                          zarf="Hkeskin (a_q = 1/(πk√q), τ_q ≤ 1.00)",
                          insa="176b_vekil_insa.main aynen")
                     for i in (5, 6, 7, 8)],
        insa_kapilari=dict(
            G1="maks|A(vekil) − A(Hkeskin)| = 0.0 VE sha256(A) eşit",
            G2="φ≡0 sağlaması: maks|ΔS| = maks|ΔS'| = 0.0",
            G3="ilk-kök hücre benzersiz = 300000 / 300000",
            G4="maks|F| ≤ 1e−8, aşan tekne = 0",
            G5="sıralılık TAM (min Δz > 0)"),
        gauss_turetim=GAUSS_TURETIM, hipotezler=HIP, ayrisim_E=AYRISIM_E,
        karar=KARAR, P2_ongoru=P2, P2_yeni_okunur_ongoru=yeni,
        P2_ulasilamaz=ulasilamaz,
        kapsam_disi=["W_pos borcu (179 ŞART ③)",
                     "VS ailesi (n=4'te kalır, büyütülmez)",
                     "gerçek gazın (son) hiçbir yeniden ölçümü"],
        sure_s=time.time() - t0)
    p = S182 / "ONKAYIT_182.json"
    p.write_text(json.dumps(rec, indent=1, ensure_ascii=False, default=float))
    print(f"\n-> {p}")
    print(f"   betik sha256 = {rec['sha256']}")
    print(f"   KALEM  sha256 = {rec['sha_kalem']}")
    print(f"   ({time.time()-t0:.1f} s)", flush=True)


if __name__ == "__main__":
    main()
