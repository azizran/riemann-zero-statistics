# -*- coding: utf-8 -*-
"""
187a — K0: KURAL-ÖNCE ÖN-KAYIT (PENCERE-ÖTESİ FAZ-İPTAL DEFTERİ)
=================================================================
KALEM_FAZ_IPTAL_10EYL2026 AYNEN. SAYI ÖLÇÜLMEDEN donan tanımlar:
katman ızgaraları, kestirimciler (katman-artım izdüşümü, ζ ve açısı,
bant/8-blok jackknife 184-186 AYNEN), H-F1 çatalının iki okuması,
H-F2 açı penceresi, H-F3/H-F4 eşikleri, K3 tanımları.

Çıktı: scratchpad/187/ONKAYIT_187.json (sha256 = bu betiğin sha'sı + damga).
"""
import hashlib
import json
import time
from pathlib import Path

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S187 = SCR / "187"
S187.mkdir(exist_ok=True)

ONKAYIT = {
    "gorev": ("187 — PENCERE-ÖTESİ FAZ-İPTAL DEFTERİ (186 bonus teşhisinin "
              "defteri: katman defteri, ζ+açı, ikiz-sağırlığı, HA4 %26)"),
    "kalem": "KALEM_FAZ_IPTAL_10EYL2026.md AYNEN",
    "L": 12.029593241726252,
    "katman_izgara_gercek": [0.86, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20],
    "katman_izgara_kontrol": [0.86, 0.90, 1.00, 1.10],
    "katman_tanim": ("katman i = τ' ∈ (ızgara[i-1], ızgara[i]] dilimindeki "
                     "TÜM asal-kuvvetler q (Λ(q)>0), τ'=log(q)/L; taban = 184 "
                     "evreni (3425 çizgi, τ≤0.86). Nominal genlik a_q = "
                     "Λ(q)/(π·√q·log q) (184 AYNEN; makine kontrolüyle "
                     "doğrulandı, maks sapma 2.8e-17)."),
    "merdiven": ("katman-artım serisi dds_i(n) = Σ_{q∈katman i} 2·a_q·"
                 "sin(ω_q·g_n/2)·cos(ω_q·m_n), ω_q=log q, ρ≡1 (NOMİNAL — "
                 "132c özdeşliğinin katman dilimi); kinematik b185.kinematik "
                 "AYNEN (gercek: eta_son kesin inversiyon; ikizler: z-diff)."),
    "izdusum": ("c^{Δi}_q = 2⟨dds_i(n)·e^{−iω_q m_n}⟩ (184 konvansiyonu; "
                "8 jk-blok toplamlarıyla, loo destekli). "
                "c^kesik_q = 186 G1_proj_<gaz>.npz re_bir/im_bir'den AYNEN "
                "(yeniden hesap YOK; makine mührü: taban serisinin "
                "dyy/dxy momentleri G1_proj dyy_bir/dxy_bir ile <1e-6 "
                "bağıl farkla örtüşmeli). "
                "c^ya_q(τ_c) = c^kesik_q + Σ_{i: ızgara[i]≤τ_c} c^{Δi}_q."),
    "m_ya": ("m^ya_b(τ_c) = [Σ_b |c^ya_q(τ_c)| − Σ_b â^öz_q] / Σ_b ae_q; "
             "â^öz ρ≡1 öz-terim 185 OZ_<gaz>.npz'den AYNEN (b185.c_oz); "
             "ae = aq_eff nominal. Hedef m_ölç,b = w_b − w^öz_b (184 K1 + "
             "185 OZ, 186b defteri AYNEN). Δm_i,b = m^ya_b(ızgara[i]) − "
             "m^ya_b(ızgara[i-1]). se: 8-blok loo-jackknife (184-186 AYNEN, "
             "se=√(7/8·Σ(θ_i−θ̄)²)); korr(ds^ya(τ_c), ds) yan sütun "
             "(tam-örneklem Pearson; ds referansı eta_son ds AYNEN)."),
    "bant": {"kenar": [0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.86],
             "havuz": ("HAVUZ sütunu = birleşik maske τ∈[0.45,0.86) "
                       "(8 bandın birliği; H-F1 çatalı BU sütunda hükmedilir; "
                       "bant tablosu yapı içindir)"),
             "njack": 8, "not": "184/185/186 AYNEN"},
    "H_F1": {
        "hipotez": ("gerçek kinematikte m^ya(τ_c) m_ölç'e MONOTON yaklaşır; "
                    "τ_c=1.20'ye dek açığın ≥%50'si kapanır. "
                    "açık(τ_c) = m^ya(τ_c) − m_ölç (HAVUZ); "
                    "kapanan pay = 1 − açık(1.20)/açık(0.86)."),
        "catal_a": ("SIĞ İPTAL/yakınsadı: art arda İKİ katmanda "
                    "|Δm(HAVUZ)| < 2·se(Δm) → mühür"),
        "catal_b": ("DERİN KUYRUK: τ_c=1.20'de kalan açık > %50 "
                    "(açık(1.20)/açık(0.86) > 0.5) → bulgu (ölçümdür, "
                    "kurtarma değil)"),
        "not": ("iki okuma da ÖN-KAYITLI ve meşru; ikisi birden ateşlerse "
                "ikisi de yazılır. Bant-bazlı tablo yapı bilgisidir.")},
    "zeta": {
        "tanim": ("Δ_q = c^ölç_q − c^kesik_q (kompleks; c^ölç 184 K1 "
                  "bloklarından AYNEN). karışım^kesik_q = c^kesik_q − c^öz_q "
                  "(ρ≡1 öz-terim; 186b m_pred tanımıyla tutarlı: kesik "
                  "merdivenin izdüşümünden öz-payın çıkarılmışı). "
                  "ζ_b = Σ_{q∈b} Δ_q·conj(karışım^kesik_q) / "
                  "Σ_{q∈b} |karışım^kesik_q|²  (bant-içi kompleks "
                  "en-küçük-kare oranı; modül = iptal payı, açı = faz)."),
        "aci": ("açı derece, (−180,180]; jackknife açı-se: loo replika "
                "açıları tam-örneklem açısına (−180,180]'e sarılarak"),
        "yan_sutun": ("genlik-okuması 1 − m_ölç,b/m^kesik_b (bilgi; "
                      "beklenen ~0.31-0.37 bölgesi ζ modülü için)"),
        "ikiz": "aynı ölçüm Hkeskin'de (küçük/gürültülü beklenir; ölçülüp yazılır)"},
    "H_F2": {"pencere": "ζ açısı ∈ 180°±15° (bant başına) → 'saf yıkıcı' mührü",
             "disi": ("penceresi dışı → döndürülmüş faz-örgüsü KAYDI "
                      "(ölüm değil, karakter tayini)")},
    "H_F3": {
        "hipotez": ("ikiz (Hkeskin) kinematiğinde τ'>1.00 katman artımları "
                    "|Δm| ≤ 2·se — kontrol ızgarasında katman (1.00,1.10]. "
                    "Hüküm HAVUZ sütununda; bant tablosu bilgi."),
        "olum": ("kontrol ÖLÜRSE (|Δm(HAVUZ)| > 2·se) 186'nın 'iptal gerçeğe "
                 "özgü' okuması düşer — açıkça raporlanır"),
        "makine_muhru": ("ikizin kendi merdiveni τ≤1.00: τ'≤1.00 katmanları "
                         "(0.86,0.90]+(0.90,1.00] eklenince m^ya_Hk(1.00) "
                         "ölçülü m_Hk'yi yeniden vermeli (bant bant bağıl "
                         "fark raporlanır; makine mührü)")},
    "H_F4": {
        "birincil": ("M^pred_HA4,b = m^env_HA4,b(τ_c=0.86; HA4-KİNEMATİĞİ) / "
                     "m_Hk,b(ölçülü). m^env_HA4,b = [Σ_b |c^ya,env_q| − "
                     "Σ_b env_q·â^öz_HA4,q]/Σ_b ae_q; c^ya,env = pencere-içi "
                     "(3425 çizgi) erfc-ağırlıklı merdivenin (katsayı "
                     "2·a_q·env(τ_q), env(τ)=0.5·erfc((τ−0.68)/0.125) 185-K4 "
                     "DONMUŞ) HA4'ün KENDİ kinematiğindeki izdüşümü — "
                     "'erfc-ağırlıklı pencere-içi karışım, İPTALSİZ; payda = "
                     "Hkeskin'in tam karışımı' (KALEM AYNEN). 186e'de bu koşu "
                     "Hk/gerçek kinematiklerindeydi; HA4-kinematiği İLK KEZ."),
        "yan_sutun": ("aynı defter + HA4-kin kontrol-ızgarası katmanlarının "
                      "env-ağırlıklı artımları (τ_c=1.10'a tam-env; 'iptalsiz' "
                      "varsayımının kontrolü; ölçüm-öncesi beklenti: küçük). "
                      "Makine mührü: m^ya,env_HA4(1.10) vs ölçülü m_HA4 "
                      "(erfc-inşanın kapanışı)."),
        "hedef": ("ölçülü M_HA4,b = m_HA4,b/m_Hk,b (m_HA4 = w_HA4 − "
                  "w^öz,erfc; 186e AYNEN yeniden üretilir, kontrol); "
                  "z = (M^pred − M_ölç)/σ_defter(M_ölç); "
                  "bant-χ²/dof ≤ 2 → 186-K4 çatlağı KAPANIR; geçemezse ölüm "
                  "dürüst.")},
    "K3": {
        "iota": ("ι(τ'_i; b) = −Δm_i,b / m^kesik_b (BİRİNCİL; katman "
                 "defterinden iptal-payı; pozitif = iptal). Yan sütun: "
                 "koherent katman-ζ_i,b = Σ_b c^{Δi}_q·conj(karışım^kesik_q)/"
                 "Σ_b |karışım^kesik_q|² (Σ_i ζ_i = ζ teleskopik). "
                 "Yoğunluk gösterimi ι/Δτ' figürde."),
        "C_rekon": ("C(τ̄_b) = m^kesik_b/m_ölç_b − 1 (186'nın 0.46→0.58 "
                    "defteri, aynı makinede yeniden); C_rekon(τ̄_b) = "
                    "m^kesik_b/m^ya_b(1.20) − 1 (katman defterinden). "
                    "Tutarlılık mührü: |C_rekon − C| bant bant jk-se'yle "
                    "kıyas — türetim DEĞİL, KAYIT."),
        "bicim": ("C(τ̄)'ye güç-yasası k·τ^β (ağırlık 1/σ²; 185d fit makinesi "
                  "α-tarama [0.05,6.0) adım 0.001, kapalı-form k); ι(τ') "
                  "HAVUZ profiline üstel A·exp(−(τ'_orta−0.86)/λ) ve güç "
                  "A·(τ'_orta−0.86)^{−p} kayıtları. TÜRETİM DEĞİL, KAYIT — "
                  "raporda açık cümleyle.")},
    "ikincil": ("187g (isteğe bağlı): gerçek, (1.20,1.30] katmanı 100k "
                "alt-örneklem (ilk 100k nokta DEĞİL; her 3. nokta n%3==0 — "
                "pencere boyu korunur), ayrı se, İKİNCİL damgası; ana "
                "hükümlere GİRMEZ. Zaman kalmazsa 'koşulmadı' yazılır."),
    "girdiler": {
        "eta_son": str(SCR / "155/eta_son_t0.4_c4000.npz"),
        "z_Hkeskin": str(SCR / "155/z_Hkeskin.npy"),
        "z_HA4": str(SCR / "155/z_HA4.npy"),
        "K1_olculu": str(SCR / "184") + "/K1_{gercek,Hkeskin,HA4}.npz (AYNEN)",
        "OZ": str(SCR / "185") + "/OZ_{gercek,Hkeskin,HA4}.npz (AYNEN)",
        "G1_proj": str(SCR / "186") + "/G1_proj_{gercek,Hkeskin}.npz (AYNEN; "
                   "kesik izdüşüm önbelleği, yeniden hesap yok)",
        "onkayit_186": "sha 23a8130a… (m/M tanımlarının kaynağı, AYNEN)"},
    "ciktilar": {
        "K1": "scratchpad/187/katman_<i>_gercek.npy (kontrol noktası) + "
              "katmanproj_gercek.npz + K1_katman_defteri.json",
        "K2": "scratchpad/187/K2_zeta.json",
        "K1_kontrol": "scratchpad/187/katman_<i>_{Hkeskin,HA4}.npy + "
                      "katmanproj_{Hkeskin,HA4}.npz + K3_ikiz_sagirlik.json",
        "K4": "scratchpad/187/K4_HA4.json",
        "K3_yapi": "scratchpad/187/K3_yapi.json",
        "rapor": "qm_riemann/187_faz_iptal_RAPOR.md + 187_faz_iptal.png"},
    "surec": ("TEK DALGA; ölümler kurtarmasız; sayı uydurmak/eşik gevşetmek/"
              "ölü kurtarmak YASAK; ölçülemeyene 'erişilemedi'; ana koşu "
              "nohup + ≤30 sn yoklama, 5 dk ön-plan YASAK; git'e DOKUNULMAZ; "
              "sonuç ORTAK TEFTİŞE, commit kaptanda."),
}

if __name__ == "__main__":
    sha = hashlib.sha256(
        (QM / "187_configs" / "187a_onkayit.py").read_bytes()).hexdigest()
    ONKAYIT["sha256"] = sha
    ONKAYIT["betik"] = str(QM / "187_configs" / "187a_onkayit.py")
    ONKAYIT["zaman"] = time.strftime("%a %b %d %H:%M:%S %z %Y")
    yol = S187 / "ONKAYIT_187.json"
    json.dump(ONKAYIT, open(yol, "w"), indent=1, ensure_ascii=False)
    print(f"ON-KAYIT DONDU: {yol}")
    print(f"  sha256 = {sha}")
    print(f"  damga  = {ONKAYIT['zaman']}")
