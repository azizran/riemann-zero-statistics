# -*- coding: utf-8 -*-
"""
188a — K0b: KURAL-ÖNCE ÖN-KAYIT (ζ(τ)'NİN KİMLİĞİ — İPTAL ÇEKİRDEĞİ HARİTASI)
=============================================================================
KALEM_ZETA_KIMLIGI_23EYL2026 AYNEN. K0 yeniden-kurulum kapısı GEÇTİKTEN SONRA,
harita ÖLÇÜLMEDEN donan tanımlar: dilim/bant ızgaraları, K/Ĝ/S/W tanımları,
uydu konumları ve pencereleri, H-188a..d eşikleri + hüküm kuralları, K3
yeniden-kurulum kuralları, alt-örneklem ve jackknife.

Çıktı: scratchpad/188/ONKAYIT_188.json (sha256 = bu betiğin sha'sı + damga).
"""
import hashlib
import json
import time
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S188 = SCR / "188"
S188.mkdir(exist_ok=True)

L = 12.029593241726252
TAU_P = {str(p): float(np.log(p) / L) for p in (2, 3, 5, 7)}

ONKAYIT = {
    "gorev": ("188 — ζ(τ)'NİN KİMLİĞİ: İPTAL ÇEKİRDEĞİ HARİTASI K(τ_b, τ') "
              "(tarak köprüsü)"),
    "kalem": "KALEM_ZETA_KIMLIGI_23EYL2026.md AYNEN",
    "L": L,
    # ---------------- ızgaralar ----------------
    "dilim_izgara_gercek": [round(0.86 + 0.005 * k, 3) for k in range(89)],
    "dilim_tanim": ("dilim s = τ' ∈ (ızgara[s], ızgara[s+1]] içindeki TÜM "
                    "asal-kuvvetler q (Λ(q)>0), τ'=log q/L; 88 dilim, genişlik "
                    "0.005, (0.86, 1.30]. Çizgi listesi 187b.asal_kuvvetler "
                    "AYNEN (sympy primerange); a_q = Λ(q)/(π√q·log q) nominal."),
    "dilim_izgara_ikiz": [round(0.86 + 0.01 * k, 2) for k in range(35)],
    "ikiz_tanim": ("İKİZ (Hkeskin) haritası: τ' ∈ (0.86, 1.20], genişlik 0.01, "
                   "34 dilim; aynı tanımlar, Hkeskin kinematiği (b185.kinematik "
                   "AYNEN: z_Hkeskin → g=diff z, m=orta nokta), aynı L ile τ'."),
    "ince_bant_kenar": [round(0.45 + 0.01 * k, 2) for k in range(42)],
    "ince_bant_tanim": ("41 ince bant τ_b ∈ [lo, hi) genişlik 0.01, "
                        "[0.45, 0.86); pencere çizgileri = 184 evreni "
                        "(K1_gercek.npz; 187 maskeleri AYNEN: tau>=lo & tau<hi)."),
    "bant8_kenar": [0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.86],
    "havuz": "HAVUZ = τ ∈ [0.45, 0.86) birleşik maske (187 AYNEN)",
    # ---------------- seriler + izdüşüm ----------------
    "seri": ("dds_s(n) = Σ_{q'∈s} 2·a_{q'}·sin(ω_{q'} g_n/2)·cos(ω_{q'} m_n), "
             "ω=log q', ρ≡1 — 187b katman_serisi tanımı AYNEN (makine mührü: "
             "bir dilimde 187b.katman_serisi ile bit-bit eşitlik)."),
    "alt_orneklem": ("n = 0, 3, 6, … (her 3. nokta; gerçek N=299999 → 100000 "
                     "nokta). Seri ve izdüşüm AYNI alt-örneklemde. jk blokları: "
                     "tam-örneklem 8 bitişik blok kenarları (linspace(0,N,9)) — "
                     "alt-örneklem noktası orijinal indisinin bloğuna düşer; loo."),
    "tam_orneklem_kontrol": ("iki dilimde (1.050,1.055] ve (0.900,0.905] seri "
                             "TAM örneklemde (299999) de kurulur/izdüşürülür; "
                             "K_alt vs K_tam karmaşık uyumu raporlanır "
                             "(|ΔK|, |ΔK|/jk-se; hükme girmez, makine kontrolü)."),
    "izdusum": ("c^{(s)}_q = 2⟨dds_s(n)·e^{−iω_q m_n}⟩_alt (184 konvansiyonu; "
                "jk-blok toplamları, loo). Pencere çizgileri τ∈[0.45,0.86)."),
    "K_tanim": ("K(b,s) = Σ_{q∈b} c^{(s)}_q·conj(karışım_q) / Σ_{q∈b} "
                "|karışım_q|²; karışım_q = c^kesik_q − c^öz_q (187c AYNEN: "
                "c^kesik 186 G1_proj re_bir/im_bir, c^öz 185 OZ ρ≡1; TAM "
                "örneklem, aynı loo bloğu dışarıda). Σ_s K(b,s) = ζ'nin "
                "(0.86,1.30] payı (ζ^{≤1.30})."),
    "G_mid": ("Ĝ_mid(ω) = mean_n exp(+iω m_n) MUTLAK m_n ile (ortalama "
              "çıkarma YOK), her pencere-ötesi çizgide ω=ω_{q'}; alt-örneklemde."),
    "S_W": ("S(s) = Σ_{q'∈s} a_{q'}·|Ĝ_mid(ω_{q'})| ; W(s) = Σ_{q'∈s} a_{q'} "
            "(yapısız ağırlık)."),
    "uydular": {
        "tau_p": TAU_P,
        "konumlar": {f"1{sg}tau_{p}": (1.0 + s * TAU_P[str(p)])
                     for p in (2, 3, 5, 7) for sg, s in (("+", 1), ("-", -1))},
        "bragg": 1.0,
        "pencere": 0.006,
        "dilim_uyelik": ("dilim, merkezi pencerenin (konum ± 0.006) içinde "
                         "ise pencerededir (yarı-örtüşme kuralıyla aynı sonucu "
                         "verdiği kontrol edildi)."),
    },
    "iptal_yonu": ("κ(s) = −Re K_HAVUZ(s) (iptal-yönlü bileşen; 187 H-F2 "
                   "mührü: açı 180°). Tüm 'pay' hesapları κ ile (karmaşık "
                   "toplamlar — 187g şerhi: modül değil karmaşık toplam)."),
    # ---------------- hipotezler ----------------
    "H_188a": {
        "hipotez": ("AYRIK GERİ-BESLEME (tarak köprüsü): K ≈ u(τ_b)·v(τ') "
                    "rang-1; |v| S'yi izler; tepeler Bragg + asal uyduları."),
        "matris": ("41 ince bant × 88 dilim KARMAŞIK K; AĞIRLIKSIZ SVD; "
                   "rang-1 varyans payı = σ1²/Σσ_k²; u = σ1·u1, v = conj(V1) "
                   "(faz: Σ_s v(s) reel-pozitif)."),
        "kapi_i": "rang-1 varyans payı ≥ 0.80",
        "kapi_ii": ("corr(|v|, S) ≥ 0.80 VE corr(|v|,S) − corr(|v|,W) ≥ 0.15 "
                    "(Pearson, 88 dilim)"),
        "kapi_iii": ("(1.05,1.10] katkısının (Σ κ, 10 dilim) ≥ %60'ı 1+τ_2 ve "
                     "1+τ_3 uydu pencerelerindeki (±0.006) dilimlerde; sıfır "
                     "beklentisi = pencere dilim sayısı / 10 (raporlanır)."),
        "olum": "(i) < 0.60 VEYA corr(|v|,S) < 0.50",
        "hukum": ("üç kapı birden → MÜHÜR; ölüm koşulu → ÖLDÜ; arası → KAYIT. "
                  "H-188c tutarsa (sıfır hipotezi) H-188a kapılardan bağımsız "
                  "ÖLDÜ (tarak köprüsü ölür)."),
        "ikincil_kayit": ("gürültü-düzeltmeli rang-1 payı σ1²/(Σσ² − tr N), "
                          "N = hücre jk-varyansı toplamı; ve imzalı öngörücü "
                          "S_Re(s) = Σ a' πτ' cos(πτ') Re Ĝ_mid(ω') ile "
                          "corr(Re v, S_Re) — İKİNCİL KAYIT, hükme girmez."),
    },
    "H_188b": {
        "hipotez": ("KÖŞEGEN / KESİME-YAKINLIK: rang-1 artığı R = K − K_1 "
                    "köşegen sırtlar taşır (τ' = 2 − τ_b katlanması ve/veya "
                    "τ' − τ_b küçük) ve |ζ|'nin kesime doğru yükselişinin "
                    "≥ %60'ını taşır."),
        "serit": ("τ_b ≥ 0.70 ince bantları (16 bant, alt kenar ≥ 0.70); hücre "
                  "şeritte ⇔ |τ'_merkez − (2 − τ_b,merkez)| ≤ 0.02 VEYA "
                  "|τ'_merkez − τ_b,merkez| ≤ 0.02."),
        "kontrast": ("mean|R| (şerit hücreleri) / mean|R| (aynı bantlarda "
                     "şerit-dışı hücreler) ≥ 3."),
        "yukselis": ("8-bant ζ^{≤1.30}(B) = Σ_s K8(B,s); K8 = ince bantların "
                     "Σ|karışım|²-ağırlıklı toplamı (doğrusal, kesin). "
                     "Ayrışım ζ8 = ζ8_r1 + ζ8_art (K_1 ve R'den aynı "
                     "ağırlıkla). z = −Re ζ8. Yükseliş Δ = z(0.80-0.86) − "
                     "min_B z(B); artık payı = [z_art(0.80-0.86) − "
                     "z_art(B_min)]/Δ ≥ 0.60."),
        "hukum": ("kontrast ≥ 3 VE artık payı ≥ 0.60 → MÜHÜR; kontrast < 3 → "
                  "ÖLDÜ (köşegen sırt yok); kontrast ≥ 3 ama pay < 0.60 → "
                  "KAYIT (sırt var, yükselişi taşımıyor)."),
    },
    "H_188c": {
        "hipotez": ("YAYVAN TABAN (sıfır): yapı çarpanı hiçbir şey eklemiyor: "
                    "corr(|v|,W) ≥ corr(|v|,S)."),
        "hukum": ("koşul sağlanırsa H-188c MÜHÜR (yayvan taban) ve tarak "
                  "köprüsü (H-188a) ÖLÜR; sağlanmazsa H-188c ÖLDÜ."),
    },
    "H_188d": {
        "hipotez": ("KİNEMATİK-GENEL / İÇERİK AYRIMI: İKİZ haritası (τ'≤1.20, "
                    "dilim 0.01) şekilce gerçeğinkine benzer."),
        "profil": ("κ(s) = −Re K_HAVUZ(s) dilim-profili; gerçek harita "
                   "0.005 dilimlerinin ikili toplamıyla 0.01 dilimlere "
                   "(34 dilim, (0.86,1.20]) indirgenir (K doğrusal — kesin)."),
        "kapi": "Pearson corr(κ_gerçek, κ_ikiz) ≥ 0.80 (34 dilim)",
        "hukum": ("≥ 0.80 → MÜHÜR (sonda kinematik-genel, ζ farkı içerik); "
                  "< 0.80 → ÖLDÜ. İKİNCİL KAYIT: |K| profilleri korelasyonu."),
    },
    "K3": {
        "tanim": ("u(τ_b) ve R'den 8-bant |ζ^{≤1.30}(τ_b)| profilinin şekli "
                  "(orta taban + kesime yükseliş) yeniden kurulur; bant-bant "
                  "pay tablosu: z_r1(B), z_art(B), z(B) (z = −Re ζ8); ölçülen "
                  "tam-derinlik |ζ| (187c) ile oran."),
        "turetilen_olculen": ("rapor 'türetilen vs ölçülen' ayrımını açık "
                              "cümleyle yazar; bu kalemdeki tüm sayılar "
                              "ÖLÇÜLEN (harita); türetim yalnız kazanan "
                              "çekirdekle ileri-hesaplanabilirse kayıt."),
        "alfa_bag": ("α=1.30 bağı: ζ-makası → M → r zinciri kazanan çekirdekle "
                     "ileri-hesaplanabiliyorsa KAYIT; değilse dürüstçe 'açık'."),
        "hukum": "KAYIT (yapı bölümü; eşik yok)",
    },
    "kontroller": {
        "K1_tutarlilik": ("HAVUZ Σ_{s: τ'≤1.20} K(b,s) ≈ 187 K3 Σζ_i(≤1.20) = "
                          "0.228∠180° (%69 koherent kapsama) — jk se ile."),
        "makine": ("seri bit-bit eşitliği (187b.katman_serisi, bir dilim, "
                   "alt-örneklem); K0'da 187c ζ defteri hane hane."),
    },
    "jackknife": "8-blok loo (184-187 AYNEN), se = √(7/8·Σ(θ_i−θ̄)²)",
    "kurallar": ("TEK DALGA; ölümler kurtarmasız; eşik gevşetme yok; ölçülemeyen "
                 "'erişilemedi'; git'e dokunulmaz; sonuç ORTAK TEFTİŞE."),
    "zaman": None,
    "sha256": None,
}

if __name__ == "__main__":
    sha = hashlib.sha256(
        (QM / "188_configs" / "188a_onkayit.py").read_bytes()).hexdigest()
    ONKAYIT["sha256"] = sha
    ONKAYIT["zaman"] = time.strftime("%a %b %d %H:%M:%S %z %Y")
    yol = S188 / "ONKAYIT_188.json"
    json.dump(ONKAYIT, open(yol, "w"), indent=1, ensure_ascii=False)
    print(f"ONKAYIT_188.json yazıldı -> {yol}")
    print(f"  sha256 = {sha}")
    print(f"  damga  = {ONKAYIT['zaman']}")
