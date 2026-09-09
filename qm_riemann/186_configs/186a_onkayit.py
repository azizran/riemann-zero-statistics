# -*- coding: utf-8 -*-
"""
186a — K0: KURAL-ÖNCE ÖN-KAYIT (M(τ) iletiminin türetimi; üç güzergâh)
======================================================================
KALEM_M_ILETIMI_09EYL2026 uyarınca; hiçbir sayı ölçülmeden yazılır ve
sha256+damga ile mühürlenir. Tanımlar 184/185 makinesinden AYNEN.
"""
import hashlib
import json
import subprocess
from pathlib import Path

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S186 = SCR / "186"
S186.mkdir(exist_ok=True)

BETIK = QM / "186_configs" / "186a_onkayit.py"

onkayit = {
    "gorev": "186 — M(τ) İLETİMİNİN TÜRETİMİ (185'in ölçtüğü 0.796→0.721; üç güzergâh)",
    "zaman": None,   # aşağıda damgalanır
    "sha256": None,  # bu dosyanın sha256'sı (185 usulü)
    "betik": str(BETIK),

    "donmus_tanim": {
        "m_M": ("w-lehçesi (185 AYNEN): w_b = Σ_b â_q/Σ_b ae_q, "
                "w_b^öz = Σ_b â_q^öz/Σ_b ae_q; m_b = w_b − w_b^öz; "
                "M_b = m_g,b/m_Hk,b. Ölçülü â'lar 184 K1 npz (AYNEN), "
                "â^öz'ler 185 OZ_*.npz bloklarından (AYNEN; yeniden hesap yok)."),
        "M_olculu_sigma": ("BİRİNCİL σ(M_ölç) = 185f defter konvansiyonu: "
                           "seM = M·√[(se_mg/mg)² + (se_mh/mh)²], "
                           "se_m = √(se_w² + se_wöz²) (bağımsız yayılım, jk-se'lerden). "
                           "Tüm G1/G2 χ²/dof hükümleri BU σ ile. "
                           "Ortak-jk se yalnız yan sütun."),

        "G1_oz_tutarlilik": {
            "merdiven": ("ds^ya_n = Σ_q' 2 a_q'·ρ(τ_q')·sin(ω_q' g_n/2)·cos(ω_q' m_n); "
                         "ρ(τ') = 1 − 0.149·τ'^1.30 (184 yasası, DONMUŞ). "
                         "BİRİNCİL kinematik = GERÇEK (KALEM AYNEN): "
                         "mid,ds eta_son önbelleğinden; g=(ds+1)·2π/log(mid/2π) "
                         "(kesin inversiyon), m_n = mid. Çizgi evreni = 184 K1 "
                         "(3425 çizgi, nominal a_q; τ≤0.86)."),
            "izdusum": "c^ya_q = 2⟨ds^ya_n e^{−iω_q m_n}⟩  (184 c_q konvansiyonu)",
            "oz_terim": ("yeniden-ağırlıklı öz-terim = ρ(τ_q)·c_q^öz(gerçek) — "
                         "KESİN (merdiven q-teriminin izdüşümü tam ρ_q·c^öz_q'dur); "
                         "OZ_gercek.npz bloklarından."),
            "m_pred": "m^pred_b = [Σ_b |c^ya_q| − Σ_b ρ_q·â_q^öz]/Σ_b ae_q",
            "M_pred": "M^pred_b = m^pred_b / m_Hk,b(ölçülü)",
            "kontrol": ("ρ≡1 AYNI makinede. MAKİNE MÜHRÜ: kontrol zinciri "
                        "r_pred(ρ≡1) = (w_g^öz + M_pred(ρ≡1)·m_Hk)/(w_Hk^öz + m_Hk) "
                        "8 bandın HEPSİNDE > 1 (185-B ters-yönü) OLMALI. "
                        "Vermezse: makine ŞÜPHELİ (hüküm değil, şüphe) — raporlanır; "
                        "teşhis: ds^ya(ρ≡1) vs gerçek ds sızıntı/kesme defteri "
                        "(korelasyon, varyans oranı, artık-varyans) zorunlu kayıt. "
                        "SEBEP (ölçüm-öncesi): 132c özdeşliği TÜM çizgilerle kesindir; "
                        "evren τ≤0.86'ya kesiktir — kontrol, kesik pencere-içi karışımın "
                        "ölçülü karışımın neresinde durduğunu da ölçer."),
            "yan_lehce": ("DONMUŞ YAN-LEHÇE (bilgi; hüküm birincille): aynı iki koşu "
                          "(ρ=yasa, ρ≡1) Hkeskin kinematiği üzerinde — 'ikiz-çekirdek' "
                          "okuması: m^pred = Σ ρ·X^Hk (karışım çekirdeği ikizden, "
                          "genlikler zarfla). ρ≡1 bu lehçede m_Hk'yi (kesme payı "
                          "kadarıyla) geri vermeli; 185-B ölümü bu uçta yapısaldır."),
        },

        "G2_v_kanali": {
            "v_tanim": ("v̂_q = |Ĝ_q| / (ḡ · ae_q);  Ĝ_q = ⟨(g_n−ḡ)e^{−iω_q m_n}⟩ "
                        "(185 OZ bloklarından AYNEN; dg tam-örneklem ḡ ile donmuş); "
                        "ae_q = 2a_q sin(πτ_q) NOMİNAL (iki denizde ÖZDEŞ dizi); "
                        "ḡ = o denizin kendi tam-pencere ortalama aralığı (OZ gbar). "
                        "İKİ DENİZDE ÖZDEŞ TANIM; HA4 dahil."),
            "bant": "v_b = Σ_b |Ĝ_q| / (ḡ · Σ_b ae_q)  (toplam-oranı lehçesi, 184 usulü)",
            "birincil_lehce": "ORAN V_b = v_g,b/v_Hk,b (konvansiyon riskine karşı)",
            "kopru": ("M^pred_b = [(1.017 − 0.884·v_g,b)/(1.017 − 0.884·v_Hk,b)]² "
                      "(KALEM AYNEN; kısıtın yüksek-τ'ya uzanması sınanan yasadır)"),
            "yan_kontrol": ("(a) her denizde ayrı ayrı √w_b vs 1.017−0.884·v_b "
                            "(kısıt-seviye kontrolü; konvansiyon farkı burada görünür); "
                            "(b) köprü² vs r_b yan-okuma (köprü cebirsel olarak w_g/w_Hk "
                            "öngörür; M'ye karşı sınav KALEM'in yazdığı gibi BİRİNCİL, "
                            "r'ye karşı yan sütun) — ikisi de bilgi."),
            "dusuk_tau_notu": ("v_g'nin en düşük banttaki değeri Not 2 doyumu "
                               "(~0.55-0.60 @ τ≈0.4) ile ORAN lehçesinde kıyaslanır; "
                               "konvansiyon farkı raporda şerh edilir."),
        },

        "G3_bicim_defteri": {
            "adaylar": ["M = a + b·τ  (doğrusal)",
                        "M = 1 − k·τ^β  (güç; β-tarama [0.05,6.0) adım 0.001, "
                        "kapalı-form k; 185d fit makinesi AYNEN)",
                        "M = a + b·V(τ)  (v-afin; V = ölçülü v_g/v_Hk, G2'den)"],
            "agirlik": "1/σ² (σ = birincil defter-seM); dof = 8−2 = 6",
            "jackknife": "parametreler 8 loo-replikada yeniden fit → jk-se",
        },
    },

    "kenar": [0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.86],
    "kuyruk_tau": 0.70,
    "njack": 8,
    "jackknife": ("8 bitişik blok, kenarlar linspace(0,N,9) (184/185 AYNEN); "
                  "leave-one-out; se = √(7/8 Σ(θ_i−θ̄)²)"),
    "bant_agirlik": "ae = aq_eff (nominal, evrenden); maskeler gerçek'in τ'suyla (AYNEN)",

    "olum_esikleri": {
        "G1": ("M^pred bant-χ²/dof ≤ 2 (8 bant; σ = birincil defter-seM) → MÜHÜR; "
               "> 2 → ÖLDÜ (kurtarmasız). Makine mührü ayrıca raporlanır."),
        "G2": "aynı eşik, aynı σ.",
        "K3_zincir": ("kazanan M^pred zincirde: r_pred_b = (w_g,b^öz + M^pred_b·m_Hk,b)"
                      "/(w_Hk,b^öz + m_Hk,b); 184 defterine (K1_faktorler r, sig_r) "
                      "χ²/dof ≤ 2 → zincir KAPANDI; değilse KAPANMADI."),
        "kazanan_kurali": ("ölüm eşiğini geçen güzergâh; ikisi geçerse χ²'si küçük olan. "
                           "HİÇBİRİ geçmezse: zincir, χ²'si küçük olana YALNIZ TEŞHİS "
                           "olarak takılır; hüküm ÖLDÜ yazılır, kurtarma yok."),
        "K4": "bonus: bant-χ²/dof ≤ 2 MÜHÜR; > 2 dürüst ÇATLAK.",
    },

    "K3": {
        "zincir": "r_pred_b = (w_g,b^öz + M^pred_b·m_Hk,b)/(w_Hk,b^öz + m_Hk,b)",
        "fit": ("1 − c·τ^α; α-tarama [0.05,6.0) adım 0.001, kapalı-form c, "
                "ağırlık 1/sig_r² (184 defter-σ); türetilmiş (c,α) vs (0.149,1.30)"),
        "sabit_nokta": ("G1 kazanırsa: ρ₁(τ) = 1 − c₁·τ^α₁ (K3 fitinden) ile G1 "
                        "makinesinde BİR yineleme daha → M^pred(2) → r_pred(2) → "
                        "(c₂,α₂); |c₂−c₁|,|α₂−α₁| yakınsama kaydı + güç-yasası "
                        "biçiminin öz-üretimi (fit-χ²/dof) rapora. Sayısal, kısa."),
    },

    "K4_HA4": {
        "hedef": ("ölçülü M_HA4,b = m_HA4,b/m_Hk,b; m_HA4 = w_HA4 − w_HA4^öz,erfc; "
                  "w_HA4^öz,erfc = Σ env(τ_q)·â_HA4,q^öz/Σae, "
                  "env(τ)=0.5·erfc((τ−0.68)/0.125) (185-K4 DONMUŞ). "
                  "185'in ~%8-10 karışım-eşitsizliği bu defterde."),
        "G1_kazanirsa": ("ds^ya(ρ=env) Hkeskin kinematiği üzerinde (HA4'ün genlik "
                         "yasası temiz denize katlanır) → m^pred_HA4; ayrıca gerçek-"
                         "kinematik birincil lehçeyle aynı koşu yan sütun."),
        "G2_kazanirsa": ("v-köprüsü: M^pred_HA4 = [(1.017−0.884·v_HA4)/"
                         "(1.017−0.884·v_Hk)]² (v_HA4 aynı donmuş kestirimciyle)."),
        "chi2": "χ²/dof, σ = defter konvansiyonuyla (M_HA4 bağımsız yayılım).",
    },

    "surec": ("TEK DALGA; ölümler kurtarmasız; sayı uydurmak/eşik gevşetmek YASAK; "
              ">2 dk koşular nohup+arka plan, ≤30 sn yoklama; git'e DOKUNULMAZ; "
              "sonuç KAPTAN+KULLANICI ortak teftişine sunulur, commit sonra."),

    "girdiler": {
        "eta_son": str(SCR / "155/eta_son_t0.4_c4000.npz"),
        "z_Hkeskin": str(SCR / "155/z_Hkeskin.npy"),
        "z_HA4": str(SCR / "155/z_HA4.npy"),
        "K1_olculu": str(SCR / "184") + "/K1_{gercek,Hkeskin,HA4}.npz (AYNEN)",
        "OZ": str(SCR / "185") + "/OZ_{gercek,Hkeskin,HA4}.npz (AYNEN)",
        "K1_faktorler": str(SCR / "185/K1_faktorler.npz"),
        "onkayit_185": "sha d49a52a4… (m/M tanımlarının kaynağı)",
    },
    "ciktilar": {
        "K1": "scratchpad/186/G1_proj.npz + G1_sonuc.json",
        "K2": "scratchpad/186/G2_sonuc.json",
        "K3": "scratchpad/186/K3_sonuc.json (+ sabit-nokta: G1_proj_iter2.npz)",
        "K4": "scratchpad/186/K4_HA4.json",
        "G3": "scratchpad/186/G3_bicim.json",
        "rapor": "qm_riemann/186_M_turetim_RAPOR.md + 186_M_turetim.png",
    },
}

if __name__ == "__main__":
    damga = subprocess.run(["date"], capture_output=True, text=True).stdout.strip()
    onkayit["zaman"] = damga
    # sha: bu dosyanın kendisi (185 usulü — zaman/sha JSON'da, dosyada değil)
    onkayit["sha256"] = hashlib.sha256(BETIK.read_bytes()).hexdigest()
    yol = S186 / "ONKAYIT_186.json"
    json.dump(onkayit, open(yol, "w"), indent=1, ensure_ascii=False)
    print(f"ON-KAYIT yazıldı: {yol}")
    print(f"  sha256 = {onkayit['sha256']}")
    print(f"  damga  = {damga}")
