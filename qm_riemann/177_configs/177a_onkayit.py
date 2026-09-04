# -*- coding: utf-8 -*-
"""
177a — ÖN-KAYIT (K2): ÜÇÜNCÜ TOHUM HAKEMLİĞİ
=============================================
Bu betik HİÇBİR 177 ÖLÇÜMÜ KOŞMADAN ÖNCE koşar ve `177/ONKAYIT_177.json`
dosyasını bir kez yazar (varlık kontrolü; üzerine yazmayı reddeder).

TEK SORU
--------
176'nın `F3`'ü (θ vekilde çöküyor mu?) **HÜKÜMSÜZ** kaldı: iki tohumun
θ'ları 0.8666712 / 0.8321358, ortalama 0.8494035, F7 saçılımı
|y1−y2| = 0.0345354; dal (a) eşiğine uzaklık 0.0162135 < saçılım.
177 yalnızca **tohum sayısını artırır** (2 → 3, gerekirse 4) ve
aynı hükmü yeniden sorar.

DEĞİŞMEYEN (176'dan bit-bit kopyalanır, bu betik hiçbirini yeniden
hesaplamaz):
  * F3'ün iki eşiği: dal (a) 0.8656169530534596, dal (b) 0.8933338132966645
  * θ çapaları: Hkeskin 0.8933338132966645, son 0.9219381611740145,
    HA4 0.7235729225877114
  * Δlog θ(son←Hk) = +0.031517828888066625,
    Δlog θ(erfc←Hk) = −0.2107589897030476
  * F9'un ödeme-oranı bandı: ω_j ∈ (0, 1]
  * F0'ın inşa kapıları G1–G5 (176b'nin kendi kodunda; dokunulmadı)
  * Vekil reçetesi: zarf birebir Hkeskin, φ_q ~ U(0,2π),
    `164_insa.coz_sadakatli`, h = 0.015, nz = 300000, c = −½
  * Ölçüm zinciri: `176c_olcum.py` (167_olcum.kos, 0.40, 0.95, kule 0,
    düz 0) — θ := KALİB_u2 / g_cal (lo = 0.60 bandı)

DEĞİŞEN: **YALNIZ TOHUM SAYISI.** Yeni tohum(lar) `tohum = 3` (VF3) ve
gerekirse `tohum = 4` (VF4) — 176'nın 1, 2 dizisinin doğal devamı;
tohum seçimi ölçümden önce dondurulur (tohum alışverişi YOK).

F7'NİN n-TOHUMLU BİÇİMİ (eşik gevşetme değil; ÖZDEŞ TALEP)
----------------------------------------------------------
176'nın F7'si iki tohumda `|y1 − y2| > |ȳ − eşik| ⇒ HÜKÜMSÜZ` der.
n = 2 için, s örnek standart sapması (ddof = 1) olmak üzere
    s = |y1 − y2| / √2  ve  SEM = s/√n = |y1 − y2| / 2
olduğundan **|y1 − y2| ≡ 2·SEM**. Yani F7'nin literal metni,
"ortalamanın eşiğe uzaklığı **iki standart hatadan** büyük olmalı"
talebinin ta kendisidir. n-tohumlu biçim bu talebi AYNEN korur:

    SAÇ_n := 2·s/√n ,   s = std(y, ddof=1)
    HÜKÜMSÜZ  ⟺  SAÇ_n > |ȳ − eşik|

n = 2'de SAÇ_2 ≡ |y1 − y2| (betik bunu 176'nın iki sayısıyla SAYISAL
olarak doğrular ve sonucu ön-kayda yazar). Eşik değişmedi; güven talebi
(2 standart hata) değişmedi; yalnız aynı talep daha çok tohumla
değerlendiriliyor.

Yanına, BAĞLAYICI OLMAYAN bir katı okuma da raporlanır:
    SAÇ^kat_n := max(y) − min(y)     (n = 2'de yine ≡ |y1 − y2|)
İkisi ayrışırsa rapor bunu açıkça yazar; hüküm SAÇ_n'e göre verilir.

HÜKÜM MERDİVENİ (dondurulmuş)
-----------------------------
  n = 3 ile (VF1, VF2, VF3):
    H1 (F3):  ȳ ≤ 0.8656169530534596 ve |ȳ − eşik| ≥ SAÇ_3  ⇒ H-F1b YAŞADI
              ȳ ≥ 0.8933338132966645 ve |ȳ − eşik| ≥ SAÇ_3  ⇒ H-F1b ÖLDÜ
              aksi hâlde: n ← 4
    H2 (ω_θ): Δ̄ := ort_i[ log θ(Hk) − log θ(VF_i) ] ;  ω_θ = kilit_θ / Δ̄
              ω_θ ∈ (0,1]  ⟺  Δ̄ ≥ kilit_θ  (kilit_θ > 0 olduğu için ÖZDEŞ)
              |Δ̄ − kilit_θ| ≥ SAÇ_n(Δ) ise KESİN, değilse n ← 4
  n = 4 ile (VF1…VF4) aynı iki hüküm.
  n = 4'te de kesinleşmezse: **KALICI HÜKÜMSÜZ** yazılır.
  EŞİK GEVŞETME YOKTUR; tohum sayısı 4'ü aşmaz.

kilit_θ (parametresiz, 176'nın F9 okumasıyla aynı — bu betikte
çapalardan yeniden türetilir, elle girilmez):
    f_θ    = Δlog θ(son←Hk) / Δlog θ(erfc←Hk)
    kesim_θ = f_θ · Δlog θ(son←Hk)
    kilit_θ = (1 − f_θ) · Δlog θ(son←Hk)

İNŞA KAPISI ÖLÜMÜ (ölçümden önce dondurulur)
--------------------------------------------
Bir tohum G1–G5'ten birini tutturamazsa: o gaz **ÖLÇÜLMEZ**, ortalamaya
**GİRMEZ**, rapora yazılır ve merdiven bir sonraki tohum indeksine geçer
(4, sonra 5). Ölçülen tohum sayısı yine 3'e (sonra 4'e) ulaşmak
zorundadır. En çok **2 kapı ölümüne** izin verilir; üçüncüsünde sefer
"KALICI HÜKÜMSÜZ" ile kapanır.

F0'IN `R_bant` ALT-KAPISI
-------------------------
176'da ıskaladı (VF1 0.9125, VF2 0.8925 < 0.98) ve kurtarılmadı; ıskanın
kendisi 176 §2c'de bir ÖLÇÜM BULGUSU olarak yazıldı. 177'de yeni tohum
için de raporlanır ve **hiçbir hükme dayanak yapılmaz** — ne lehte ne
aleyhte. Bu, ön-kayıt anında böyle yazılmıştır.

DÜRÜSTLÜK BEYANI
----------------
1. KÖRLÜK İDDİASI YOKTUR. 176'nın bütün sayıları (θ: 0.8666712 /
   0.8321358, ω_θ = +0.716, saçılım 0.0345354) bu tayfa tarafından
   okunmuştur. Ön-kayıt edilen şey KURALDIR.
2. ÖN-KAYIT ANINDA ÖLÇÜLMEMİŞ OLAN: `tohum = 3` (ve 4) ile kurulan
   vekil gazın HİÇBİR niceliği. Böyle bir gaz kurulmamıştır.
3. Ön-kayıt sonrası hiçbir eşik değiştirilmez. Ölümler kurtarmasız.
4. Git'e dokunulmaz.
"""
import hashlib
import json
import math
import time
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S176, S177 = SCR / "176", SCR / "177"
BEN = QM / "177_configs" / "177a_onkayit.py"

# ── 176'nın ön-kaydından BİT-BİT alınan çapalar ────────────────────────
ONK176 = json.load(open(S176 / "ONKAYIT_K2.json"))
F3_176 = ONK176["kural"]["F3"]
TH_A = F3_176["dal_a_yasar"]                 # 0.8656169530534596
TH_B = F3_176["dal_b_olur"]                  # 0.8933338132966645
TH_HK = F3_176["capa"]["Hkeskin"]
TH_SON = F3_176["capa"]["son"]
TH_HA4 = F3_176["capa"]["HA4"]
D_SON = F3_176["capa"]["dlog_th(son<-Hk)"]
D_ERFC = F3_176["capa"]["dlog_th(erfc<-Hk)"]

# 176'da ÖLÇÜLMÜŞ iki tohum (bu betik onları yeniden ölçmez, okur)
TH_VF1 = json.load(open(S176 / "G_VF1.json"))["th"]
TH_VF2 = json.load(open(S176 / "G_VF2.json"))["th"]


def sac_n(y):
    """F7'nin n-tohumlu biçimi: 2·s/√n  (n=2'de ≡ |y1−y2|)."""
    y = np.asarray(y, dtype=float)
    n = len(y)
    if n < 2:
        return float("nan")
    return float(2.0 * np.std(y, ddof=1) / math.sqrt(n))


def sac_kat(y):
    """Bağlayıcı olmayan katı okuma: aralık (n=2'de yine ≡ |y1−y2|)."""
    y = np.asarray(y, dtype=float)
    return float(y.max() - y.min()) if len(y) > 1 else float("nan")


def main():
    p = S177 / "ONKAYIT_177.json"
    S177.mkdir(parents=True, exist_ok=True)
    if p.exists():
        raise SystemExit(f"ÖN-KAYIT ZATEN VAR — üzerine YAZILMAZ: {p}")

    sha = hashlib.sha256(BEN.read_bytes()).hexdigest()
    zaman = time.strftime("%Y-%m-%d %H:%M:%S %z")

    print("=" * 78)
    print("177a — ÖN-KAYIT (K2): ÜÇÜNCÜ TOHUM HAKEMLİĞİ")
    print("=" * 78)
    print(f"  zaman  : {zaman}")
    print(f"  betik  : {BEN.name}")
    print(f"  sha256 : {sha}")

    # ── F7 ÖZDEŞLİK SAĞLAMASI (n=2'de yeni biçim = 176'nın literal metni)
    y2 = [TH_VF1, TH_VF2]
    literal = abs(y2[0] - y2[1])
    yeni = sac_n(y2)
    katı = sac_kat(y2)
    fark = abs(literal - yeni)
    print("\n" + "-" * 78)
    print("  F7 ÖZDEŞLİK SAĞLAMASI (n = 2, 176'nın iki tohumu)")
    print(f"    176'nın literal metni |y1−y2| = {literal:.16f}")
    print(f"    177'nin biçimi 2·s/√n         = {yeni:.16f}")
    print(f"    katı okuma  max−min           = {katı:.16f}")
    print(f"    |fark| = {fark:.3e}   ⇒  {'ÖZDEŞ ✓' if fark < 1e-15 else '✗ ÖZDEŞ DEĞİL'}")
    if fark >= 1e-15:
        raise SystemExit("F7 genellemesi n=2'de literal metne indirgenmiyor "
                         "— ön-kayıt YAZILMAZ.")

    # ── kilit_θ: 176'nın F9 okumasıyla, çapalardan türetilir ───────────
    f_th = D_SON / D_ERFC
    kesim_th = f_th * D_SON
    kilit_th = (1.0 - f_th) * D_SON
    print("\n" + "-" * 78)
    print("  kilit_θ (parametresiz, çapalardan türetildi)")
    print(f"    f_θ     = Δ(son←Hk)/Δ(erfc←Hk) = {f_th:+.7f}")
    print(f"    kesim_θ = f_θ·Δ(son←Hk)        = {kesim_th:+.7f}")
    print(f"    kilit_θ = (1−f_θ)·Δ(son←Hk)    = {kilit_th:+.7f}")
    print(f"    ω_θ ∈ (0,1]  ⟺  Δ̄(Hk←vekil) ≥ {kilit_th:.7f}")

    # ── 176'nın iki-tohumlu durumu (kayda geçer; hüküm DEĞİLDİR) ───────
    ort2 = float(np.mean(y2))
    print("\n" + "-" * 78)
    print("  176'nın bıraktığı yer (yeniden ölçüm YOK, yalnız okundu)")
    print(f"    θ(VF1) = {TH_VF1:.7f}   θ(VF2) = {TH_VF2:.7f}   "
          f"ort = {ort2:.7f}")
    print(f"    dal (a) eşiği {TH_A:.7f}  → uzaklık {abs(ort2-TH_A):.7f}"
          f"   SAÇ_2 = {yeni:.7f}  ⇒ HÜKÜMSÜZ")
    print(f"    dal (b) eşiği {TH_B:.7f}  → uzaklık {abs(ort2-TH_B):.7f}"
          f"   ⇒ (b) KESİN olarak dışlanmıştı")

    kural = dict(
        soru="H-F1b hakemliği: θ'nın kilit-duyarlılığı YAŞADI mı ÖLDÜ mü?",
        degisen="YALNIZ TOHUM SAYISI (2 → 3, gerekirse 4). "
                "176'nın hiçbir eşiği değiştirilmemiştir.",
        F3=dict(dal_a_yasar=TH_A, dal_b_olur=TH_B,
                capa=dict(Hkeskin=TH_HK, son=TH_SON, HA4=TH_HA4,
                          dlog_th_son_Hk=D_SON, dlog_th_erfc_Hk=D_ERFC),
                olculen176=dict(VF1=TH_VF1, VF2=TH_VF2, ort2=ort2,
                                sacilim2=literal)),
        F7n=dict(formul="SAÇ_n := 2·std(y, ddof=1)/√n",
                 ozdeslik="n=2'de SAÇ_2 ≡ |y1−y2| (176'nın literal metni); "
                          "sayısal sağlama bu ön-kayıtta yapıldı",
                 sagla_literal=literal, sagla_yeni=yeni, sagla_fark=fark,
                 kati_okuma="SAÇ^kat_n := max(y)−min(y) — BAĞLAYICI DEĞİL, "
                            "yalnız raporlanır",
                 kati2=katı,
                 kural="SAÇ_n > |ȳ − eşik|  ⇒  HÜKÜMSÜZ"),
        H1=dict(ad="F3 / H-F1b",
                yasar="ȳ ≤ %.16f  ve  |ȳ−eşik| ≥ SAÇ_n" % TH_A,
                olur="ȳ ≥ %.16f  ve  |ȳ−eşik| ≥ SAÇ_n" % TH_B,
                aksi="n ← n+1 (en çok 4); n=4'te de olmazsa KALICI HÜKÜMSÜZ"),
        H2=dict(ad="F9 / ω_θ bandı",
                delta="Δ_i := log θ(Hkeskin) − log θ(VF_i);  Δ̄ = ort_i Δ_i",
                omega="ω_θ = kilit_θ / Δ̄",
                kilit_th=kilit_th, kesim_th=kesim_th, f_th=f_th,
                band="ω_θ ∈ (0,1]  ⟺  Δ̄ ≥ kilit_θ",
                kesinlik="|Δ̄ − kilit_θ| ≥ SAÇ_n(Δ) ise KESİN"),
        merdiven=dict(n3=["VF1", "VF2", "VF3"], n4=["VF1", "VF2", "VF3", "VF4"],
                      tohumlar={"VF3": 3, "VF4": 4},
                      ust_sinir=4,
                      kapi_olumu="G1–G5'ten biri tutmazsa gaz ÖLÇÜLMEZ ve "
                                 "ortalamaya girmez; merdiven bir sonraki "
                                 "tohum indeksine geçer (en çok 2 kapı "
                                 "ölümü)"),
        F0=dict(insa_kapilari="G1–G5, 176b'nin kendi kodunda; değiştirilmedi",
                R_bant="176'da ıskaladı (0.9125 / 0.8925 < 0.98) ve "
                       "kurtarılmadı; 177'de raporlanır, HİÇBİR HÜKME "
                       "DAYANAK YAPILMAZ"),
        makine=dict(insa="176_configs/176b_vekil_insa.py (değiştirilmeden, "
                         "alt-süreç olarak çağrılır)",
                    olcum="176_configs/176c_olcum.py (değiştirilmeden, "
                          "alt-süreç olarak çağrılır)",
                    theta="θ := KALİB_u2 / g_cal (lo=0.60), G_<ad>.json'daki "
                          "'th' alanı — başka kestirici YOK"),
        kapsam="YALNIZ bu soru. R_η, g_E, M, DC kaçağı, H-F2 176'da "
               "hükme bağlandı; 177 onları yeniden açmaz.",
    )

    rec = dict(zaman=zaman, betik=BEN.name, sha256=sha,
               onkayit176=dict(zaman=ONK176["zaman"],
                               sha256=ONK176["sha256"]),
               kural=kural,
               durustluk=[
                   "korluk iddiasi yoktur; 176'nin butun sayilari okunmustur",
                   "on-kayit aninda olculmemis olan: tohum=3 (ve 4) vekilin "
                   "HICBIR niceligi",
                   "176'nin esikleri degistirilmedi; degisen yalniz tohum "
                   "sayisidir",
                   "F7'nin n-tohumlu bicimi n=2'de literal metne SAYISAL "
                   "olarak indirgenir (bu betikte dogrulandi)",
                   "olumler kurtarmasiz; git'e dokunulmaz",
               ])
    p.write_text(json.dumps(rec, indent=1, ensure_ascii=False))
    print("\n" + "=" * 78)
    print(f"  ÖN-KAYIT YAZILDI (bir kez, üzerine yazılmaz): {p}")
    print("=" * 78)
    return rec


if __name__ == "__main__":
    main()
