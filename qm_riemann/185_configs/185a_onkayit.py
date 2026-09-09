# -*- coding: utf-8 -*-
"""
185a — K0: KURAL-ÖNCE ÖN-KAYIT (zarf türetimi: neden 1 − 0.149·τ^1.30?)
=======================================================================
KALEM_W_TURETIMI_09EYL2026 uyarınca: â^öz/Ĝ/faktör tanımları, bant ızgarası
(184'ünkü AYNEN), jackknife düzeni, ölüm eşikleri ve H-W0..W3 formları
SAYILAR DOĞMADAN donar. Bu betik hiçbir zarf sayısı ölçmez; yalnız
tanımları yazar ve kendi sha256'sını damgalar.

Ön-kayıttan ÖNCE yapılmış tek şey MAKİNE-KONVANSİYONU doğrulamasıdır
(184'ün YAYIMLANMIŞ defter sayılarının hangi muhasebeyle çıktığının
tespiti — yeni sayı üretmez):
  - bant değeri  w̄_b = Σ_b â / Σ_b ae   (ae = aq_eff; genlik-ağırlık)
  - bant zarfı   r_b = Σ_b â_g / Σ_b â_Hk  (184 tablosunu 3-4 hane tutturur)
  - rapor se'si  = w_g ve w_Hk loo-jackknife se'lerinin BAĞIMSIZ yayılımı
    (184 defterindeki ±.002…±.013'ü hane hane tutturur)
  - gerçek gazın g_n'si eta önbelleğinden KESİN inversiyon:
    g = (ds+1)·2π/log(mid/2π)  (maks kalıntı 2.3e-10; 184b veri yolu
    yalnız mid+ds taşıdığından g bu yoldan gelir)

Çıktı: scratchpad/185/ONKAYIT_185.json
"""
import hashlib
import json
import time
from pathlib import Path

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S185 = SCR / "185"

ONKAYIT = {
    "gorev": "185 — ZARFIN TÜRETİMİ (K0 on-kayit; 184 zarfı 1-0.149·tau^1.30)",
    "zaman": None,      # damga aşağıda
    "sha256": None,     # bu betiğin sha'sı aşağıda
    "betik": str(QM / "185_configs" / "185a_onkayit.py"),

    # ---------------- donmuş tanımlar ----------------
    "donmus_tanim": {
        "ozdeslik": "ds_n = Σ_Q 2 a_Q sin(w_Q g_n/2) cos(w_Q m_n)  (132c, kesin)",
        "aq": "a_q = 1/(pi*m*sqrt(q)), q=p^m; w_q = log q; tau = w_q/L",
        "ahat_olculu": "184 K1 npz'lerinden AYNEN (re_blok,im_blok; yeniden ölçüm YOK)",
        "c_oz": ("c_q^oz = (4 a_q / N) Σ_n sin(w_q g_n/2) cos(w_q m_n) e^{-i w_q m_n}"
                 "  — KESİN beklenen-değer, Taylor YOK; a_q NOMİNAL"),
        "ahat_oz": "â_q^oz = |c_q^oz|",
        "oz_ayrisim": ("cos(P)e^{-iP} = 1/2 + e^{-2iP}/2  ⇒  c^oz = (2a_q/N)(S_dc + S_2w);"
                       " S_dc = Σ sin(w_q g/2), S_2w = Σ sin(w_q g/2) e^{-2i w_q m}"
                       " (anatomi için ayrı saklanır)"),
        "G_hat": "Ĝ_q = < (g_n - ḡ) e^{-i w_q m_n} >,  ḡ = tüm-pencere ortalama aralık",
        "kinematik_kaynak": {
            "gercek": ("mid,ds: 155/eta_son_t0.4_c4000.npz; g = (ds+1)·2π/log(mid/2π)"
                       " (kesin inversiyon, kalıntı 2.3e-10)"),
            "Hkeskin": "z: 155/z_Hkeskin.npy; g=diff(z), m=(z_i+z_{i+1})/2 (184b2 yolu)",
            "HA4": "z: 155/z_HA4.npy; aynı yol",
        },
        "uc_faktor": ("r = [â_g/â_g^oz] × [â_g^oz/â_Hk^oz] × [â_Hk^oz/â_Hk]"
                      " = F_ANOMALI × F_KINEMATIK × F_KOMSU^-1  (kesin, teleskopik)"),
        "faktor_bant_defteri": ("her faktör bant düzeyinde TOPLAM-ORANI ile:"
                                " F1_b=Σâ_g/Σâ_g^oz, F2_b=Σâ_g^oz/Σâ_Hk^oz,"
                                " F3inv_b=Σâ_Hk^oz/Σâ_Hk  (184'ün r_b lehçesi;"
                                " çarpım r_b'ye TAM teleskoplanır)"),
        "cizgi_lehcesi": ("H-W1b için ayrıca çizgi-düzeyi lehçe: F̄i_b ="
                          " Σ ae·Fi_q / Σ ae  (çizgi-başına faktörlerin ae-ağırlıklı"
                          " ortalaması; teleskoplanmaz, sınav budur)"),
    },

    # ---------------- bant/jackknife (184 AYNEN) ----------------
    "kenar": [0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.86],
    "kuyruk_tau": 0.70,
    "njack": 8,
    "jackknife": ("8 bitişik blok, kenarlar linspace(0,N,9) (N=299999; nb 184 ile"
                  " özdeş); leave-one-out; se = sqrt(7/8 Σ(θ_i-θ̄)²)"),
    "bant_agirlik": "ae = aq_eff(gercek); bant maskeleri gercek'in tau'suyla",
    "hata_konvansiyonu": {
        "birincil": ("184 DEFTER se'si: σ_r(b) = r_b·sqrt((se_wg/w̄g)²+(se_wHk/w̄Hk)²)"
                     " — bağımsız yayılım; tüm χ²/dof ve 2σ hükümleri BUNUNLA"),
        "yan_sutun": "ortak-(korelasyonlu)-jackknife se'si yalnız bilgi amaçlı raporlanır",
    },

    # ---------------- hedef (184 defteri; kopya) ----------------
    "hedef_184": {
        "kazanan_fit": "r = 1 - 0.149·tau^1.30, chi2/dof = 0.34",
        "r_bant": [0.944, 0.935, 0.925, 0.917, 0.910, 0.902, 0.898, 0.893],
        "r_se":   [0.002, 0.002, 0.003, 0.004, 0.004, 0.004, 0.008, 0.013],
        "not": "χ² hedefi bu defter; merkezler tam hassasiyetle 184 npz'den yeniden kurulur",
    },

    # ---------------- hipotezler + ölüm eşikleri (kurtarmasız) ----------------
    "hipotez": {
        "H-W0": {
            "ad": "L-akışı sıfır-adayı",
            "form": ("pencere içi L kayması: L_min/L_maks son penceresinin ilk/son"
                     " sıfırından; δτ = -τ·δL/L̄; â^oz'e birinci-mertebe etkisi"
                     " gerçek-ikiz FARKINDA sınırlanır (ortak-mod düşer)"),
            "esik": "zarfa katkı sınırı < %1 (1-r'nin yüzdesi olarak) → KAPANIR",
        },
        "H-W1": {
            "ad": "ÖZ-MUHASEBE",
            "a_yapisal": ("toplam-oranı lehçesinde çarpım r_b'ye teleskopiktir;"
                          " maks |Π F_b − r_b| makine-kapanışı olarak raporlanır"
                          " (beklenti ~1e-15; bu tek başına MÜHÜR SAYILMAZ)"),
            "b_cizgi": ("çizgi-lehçesi rekonstrüksiyonu r̂_b = F̄1·F̄2·F̄3inv,"
                        " 184 defterine karşı χ²/dof (birincil σ ile), dof=8"),
            "esik": "ÖLÜM: H-W1b bant-χ²/dof > 2  → zarf özdeşlik-makinesinin dışında",
        },
        "H-W2": {
            "ad": "BASKIN TERİM + KÖPRÜ",
            "baskin_tanim": ("kuyruk bandında (τ>0.70) log-pay: pay_i = log F̄i /"
                             " log r; en büyük |pay| baskındır (toplam-oranı lehçesi)"),
            "ileri_hesap_adaylari": {
                "A_carpimsal_kinematik": ("r_pred^A_b = F2_b = Σâ_g^oz/Σâ_Hk^oz —"
                                          " yalnız aralık+orta-nokta serilerinden;"
                                          " ölçülü â HİÇ girmez"),
                "B_toplamsal_ortak_karisim": ("c_g^pred = (c_Hk^olc − c_Hk^oz) + c_g^oz"
                                              " (kompleks, çizgi çizgi);"
                                              " r_pred^B_b = Σ|c_g^pred|/Σâ_Hk —"
                                              " ölçülü â_g HİÇ girmez; komşu-karışım"
                                              " ortak-mod varsayımı"),
            },
            "mekanizma_esleme": ("F2 baskınsa birincil = A; F1/F3 çifti baskınsa"
                                 " birincil = B (öz-terim değişimi). İkisi de raporlanır."),
            "esik": ("MÜHÜR: birincil ileri-hesap tüm 8 bantta |Δ|≤2σ_b VE χ²/dof≤2;"
                     " ÖLÜM: kuyruk bandı >2σ VEYA χ²/dof>2; arası KISMİ"),
        },
        "H-W3": {
            "ad": "kuvvet-açığı (script-100 mekanik yasası)",
            "form": ("f ≈ τ^3.3 biçimli zarf adayı 1 − c3·τ^3.3; ölçek bilgisi"
                     " okuma listesinde YOKSA c3 için yalnız EN-İYİ-ÖLÇEK sınırı"
                     " (tek-parametreli şekil sınavı) verilir; ölçeğin kendisi"
                     " 'erişilemedi' diye yazılır, uydurulmaz"),
            "esik": "dürüst kayıt: α=3.3 vs 1.30 kıyası + en-iyi-ölçek χ²/dof",
        },
    },

    # ---------------- K2 anatomi (yalnız yorum) ----------------
    "anatomi": {
        "T1": "(πτ)·cot(πτ)·Re[Ĝ_q]/ḡ  (uyumlu geri-besleme; YALNIZ YORUM)",
        "T2": "-(πτ)²·Δσ_ε²/2,  ε=(g-ḡ)/ḡ, Δ = gerçek − ikiz  (YALNIZ YORUM)",
        "not": "πτσ_ε ~ 1 ⇒ hesap DAİMA kesin beklenen-değerle; Taylor teşhis dili",
        "G_profil": "Ĝ(τ) bant profili ayrıca kaydedilir (H-W2 köprüsü)",
    },

    # ---------------- K3 kapanış ----------------
    "K3": {
        "r_pred_fit": ("birincil r_pred'e 1−c·τ^α fiti (birincil σ ile ağırlıklı,"
                       " bant τ̄'larında) → türetilmiş (c,α), 184'ün (0.149,1.30)'uyla"
                       " yan yana tablo"),
    },

    # ---------------- K4 HA4 sınavı (örneklem-dışı; sıfır yeni ayar) ----------------
    "K4": {
        "erfc_zarf": "env(τ) = 0.5·erfc((τ−0.68)/0.125)  (164/184 defterinden; donmuş)",
        "oz_HA4": ("c_HA4^oz,erfc = env(τ)·c_HA4^oz(nominal formül, HA4'ün kendi g,m'si)"
                   " — HA4'ün merdiveni erfc-sönümlü olduğundan öz-terimi env ile ölçeklenir"),
        "birincil_b": ("TOPLAMSAL ÖZ-DEĞİŞİM: c_HA4^pred = (c_Hk^olc − c_Hk^oz)"
                       " + c_HA4^oz,erfc; w_HA4^pred_b = Σ|c^pred|/Σae, ölçülü"
                       " w_HA4 defterine (K1_HA4.npz) karşı χ²/dof (HA4 jk se'yle,"
                       " bağımsız yayılım konvansiyonu). GEREKÇE (ölçüm-öncesi):"
                       " çarpımsal taşıma ölü-merdiven kuyruğunda sonlu â üretemez;"
                       " kalem 'yeniden üretmeli' der, tek aday toplamsal değişimdir."),
        "yan_a": ("ÇARPIMSAL: r_HA4^pred = Σ|c_HA4^oz,erfc|/Σ|c_Hk^oz| vs"
                  " r_HA4^olc = Σâ_HA4/Σâ_Hk — yalnız teşhis sütunu"),
        "esik": "bonus sınav: bant-χ²/dof raporlanır; ≤2 MÜHÜR, >2 dürüst ÇATLAK",
    },

    # ---------------- süreç ----------------
    "surec": ("TEK DALGA; ölümler kurtarmasız; sayı uydurmak/eşik gevşetmek YASAK;"
              " uzun koşular nohup+kısa yoklama; git'e dokunulmaz"),
    "girdiler": {
        "eta_son": "scratchpad/155/eta_son_t0.4_c4000.npz",
        "z_Hkeskin": "scratchpad/155/z_Hkeskin.npy",
        "z_HA4": "scratchpad/155/z_HA4.npy",
        "K1_olculu": "scratchpad/184/K1_{gercek,Hkeskin,HA4}.npz (AYNEN)",
        "onkayit_184": "scratchpad/184/ONKAYIT_184.json (sha 1e0de637…)",
    },
    "ciktilar": {
        "K1": "scratchpad/185/K1_faktorler.npz + OZ_<gaz>.npz",
        "K2": "scratchpad/185/K2_anatomi.json (+ G profili K1 npz içinde)",
        "K3": "scratchpad/185/K3_kapanis.json",
        "K4": "scratchpad/185/K4_HA4.json",
        "rapor": "qm_riemann/185_turetim_RAPOR.md + 185_turetim.png",
    },
}

if __name__ == "__main__":
    S185.mkdir(parents=True, exist_ok=True)
    ONKAYIT["zaman"] = time.strftime("%a %b %e %H:%M:%S %z %Y")
    sha = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    ONKAYIT["sha256"] = sha
    yol = S185 / "ONKAYIT_185.json"
    yol.write_text(json.dumps(ONKAYIT, indent=1, ensure_ascii=False))
    print(f"ONKAYIT_185.json yazıldı  sha256={sha}")
    print(f"damga: {ONKAYIT['zaman']}")
