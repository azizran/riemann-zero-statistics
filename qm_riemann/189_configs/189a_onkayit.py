# -*- coding: utf-8 -*-
"""
189a — K0: ÖN-KAYIT (DERİN İKİZ — ZARFIN KİMLİĞİ SINAVI)
=========================================================
KALEM_DERIN_IKIZ_23EYL2026.md AYNEN. İNŞADAN ÖNCE donan: ön-mühür öngörü
tablosu (kalemden AYNEN kopya + kaynağından yeniden hesap ve eşitlik kontrolü),
H-189a..d eşikleri ve hüküm kuralları, ikiz adları (Hderin110, Hderin120),
ızgara sınavı kuralı, inşa kapıları, ölçüm zinciri, makine mührü, K3 kimlik
kaydı tanımları (fit makinesi + D→∞ modelleri).

Çıktı: scratchpad/189/ONKAYIT_189.json (sha256 = bu betiğin sha'sı + damga).
"""
import hashlib
import importlib.util
import json
import time
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S189 = SCR / "189"
S189.mkdir(exist_ok=True)

BANTLAR = ["0.45-0.50", "0.50-0.55", "0.55-0.60", "0.60-0.65",
           "0.65-0.70", "0.70-0.75", "0.75-0.80", "0.80-0.86"]

# ---- KALEMDEN AYNEN (ön-mühür tablosu; kaptan, 23 Eyl, 188 haritasından) ----
# sütunlar: r_1.00(184), z_1.00, z_1.10, z_1.20, r_1.10 öng, r_1.20 öng,
#           kapanan % 1.10, kapanan % 1.20
ONGORU_KALEM = {
    "0.45-0.50": [0.9445, 0.0920, 0.1884, 0.2480, 0.9802, 1.0037, 64, 107],
    "0.50-0.55": [0.9351, 0.0841, 0.1792, 0.2421, 0.9747, 1.0028, 61, 104],
    "0.55-0.60": [0.9253, 0.0821, 0.1808, 0.2418, 0.9708, 1.0013, 61, 102],
    "0.60-0.65": [0.9169, 0.0803, 0.1788, 0.2403, 0.9666, 1.0003, 60, 100],
    "0.65-0.70": [0.9096, 0.0705, 0.1694, 0.2314, 0.9624, 0.9987, 58, 99],
    "0.70-0.75": [0.9021, 0.0678, 0.1752, 0.2323, 0.9628, 0.9985, 62, 99],
    "0.75-0.80": [0.8978, 0.0710, 0.1886, 0.2619, 0.9677, 1.0171, 68, 117],
    "0.80-0.86": [0.8931, 0.0637, 0.1810, 0.2629, 0.9638, 1.0203, 66, 119],
}
SUTUNLAR = ["r_1.00", "z_1.00", "z_1.10", "z_1.20", "r_1.10_ong",
            "r_1.20_ong", "kapanan_pct_1.10", "kapanan_pct_1.20"]


def yeniden_hesap():
    """Öngörüyü KAYNAĞINDAN yeniden hesapla (185/K1_faktorler + 188 ikiz haritası)."""
    f = np.load(SCR / "185" / "K1_faktorler.npz", allow_pickle=True)
    h = np.load(SCR / "188" / "harita_K_Hkeskin.npz", allow_pickle=True)
    et = list(h["etiket"])
    orta = h["orta"]
    K = h["K"]
    out = {}
    for b, ad in enumerate(BANTLAR):
        assert str(f["etiket"][b]) == ad
        i = et.index("B" + ad)
        z = {D: float(-np.real(K[i, orta <= D + 1e-12].sum()))
             for D in (1.00, 1.10, 1.20)}
        wg, wh, wog, woh = (float(f[k][b]) for k in ("wg", "wh", "wog", "woh"))
        r0 = float(f["r"][b])
        mh = wh - woh
        mk = mh / (1.0 - z[1.00])
        rD = {D: wg / (woh + mk * (1.0 - z[D])) for D in z}
        out[ad] = [r0, z[1.00], z[1.10], z[1.20], rD[1.10], rD[1.20],
                   100 * (rD[1.10] - r0) / (1 - r0),
                   100 * (rD[1.20] - r0) / (1 - r0)]
    # HAVUZ (kalemde YOK; aynı aritmetik, bilgi amaçlı)
    spec = importlib.util.spec_from_file_location(
        "b185", QM / "185_configs" / "185b_oz_muhasebe.py")
    b185 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(b185)
    G = np.load(SCR / "184" / "K1_gercek.npz")
    H = np.load(SCR / "184" / "K1_Hkeskin.npz")
    OZh = np.load(SCR / "185" / "OZ_Hkeskin.npz")
    tau, ae, aq = G["tau"], G["aq_eff"], G["aq"]
    m = (tau >= 0.45) & (tau < 0.86)
    sae = ae[m].sum()
    wg = np.abs(b185.c_olculu(G))[m].sum() / sae
    wh = np.abs(b185.c_olculu(H))[m].sum() / sae
    woh = np.abs(b185.c_oz(OZh, aq))[m].sum() / sae
    i = et.index("HAVUZ")
    z = {D: float(-np.real(K[i, orta <= D + 1e-12].sum()))
         for D in (1.00, 1.10, 1.20)}
    r0 = wg / wh
    mk = (wh - woh) / (1.0 - z[1.00])
    rD = {D: wg / (woh + mk * (1.0 - z[D])) for D in z}
    havuz = {"r_1.00": float(r0), "z_1.00": z[1.00], "z_1.10": z[1.10],
             "z_1.20": z[1.20], "r_1.10_ong": float(rD[1.10]),
             "r_1.20_ong": float(rD[1.20]),
             "kapanan_pct_1.10": float(100 * (rD[1.10] - r0) / (1 - r0)),
             "kapanan_pct_1.20": float(100 * (rD[1.20] - r0) / (1 - r0)),
             "w_g": float(wg), "w_Hk": float(wh), "w_oz_Hk": float(woh)}
    return out, havuz


ONKAYIT = {
    "gorev": "189 — DERİN İKİZ: ZARFIN KİMLİĞİ SINAVI",
    "kalem": "KALEM_DERIN_IKIZ_23EYL2026.md AYNEN",
    "kalem_sha256": None,
    # ---------------- ikizler ----------------
    "ikizler": {
        "Hderin110": {"tau_ust": 1.10, "konfig": ["sadakatli", {"c": -0.5}]},
        "Hderin120": {"tau_ust": 1.20, "konfig": ["sadakatli", {"c": -0.5}]},
        "ortak": ("164_insa.py DÜZENLENMEDEN, sarmalayıcıyla (189b_insa.py): "
                  "modül importlib ile yüklenir, merdiven = partial(merdiven, "
                  "tau_ust=D) modül özniteliği olarak değiştirilir, KONFIG'e "
                  "Hderin110/Hderin120 = ('sadakatli', {'c': -0.5}) eklenir, "
                  "main(ad, h, 300000) çağrılır. Diğer her şey Hkeskin'in "
                  "aynısı: sadakatli çözücü (ızgara-braketi + sıralı ilk-kök + "
                  "korumalı Newton), seviye c = −½, nz = 300000, zeros6 son 300k "
                  "penceresi, L_hedef = log(½(t0+t1)/2π) (inşa-kesimi "
                  "lim = exp(D·L_hedef)), NWORK = 7. Çıktı main AYNEN: "
                  "scratchpad/164 ve scratchpad/155'e z_<ad>.npy + "
                  "164/insa_<ad>.json. İki inşa aynı anda KOŞULMAZ; önce 1.10."),
    },
    "izgara_sinavi": {
        "dilim": ("ilk 5000 seviye: n = n0 … n0+4999 (n0 = ceil(N̄(t0))), "
                  "c = −½; aynı derin merdiven; coz_sadakatli doğrudan "
                  "(main'in dosya yazmaması için)"),
        "kiyas": "h = 0.015 vs h = 0.0075",
        "fark_tanim": "ilk-kök farkı: |z(h=0.015) − z(h=0.0075)| > 1e-6 olan tekne",
        "kural": ("fark ≤ 5/5000 ⇒ h = 0.015 ile inşa; aksi hâlde h = 0.0075 "
                  "ile inşa (raporlanır). Her derinlik için AYRI sınav, "
                  "inşa h'si betik tarafından sınav JSON'undan otomatik seçilir."),
    },
    "insa_kapilari": {
        "maks_F": "maks|F| ≤ 1e-8 (tam S ile son değerlendirme, 164 AYNEN)",
        "siralilik": "TAM (min Δz > 0, sıra bozan çift = 0)",
        "L": "round(L, 8) == 12.02959324 (Hkeskin ile aynı pencere)",
        "kapi_kalirsa": "kapıdan kalan ikiz ZİNCİRE GİRMEZ; 'erişilemedi' yazılır",
    },
    # ---------------- zincir ----------------
    "zincir": {
        "adimlar": ("184b2.kos(ad, scratchpad/155/z_<ad>.npy, q, w, mq) → "
                    "185b.oz_hesapla(ad, w) [+ 185b defter AYNEN, HAVUZ satırı "
                    "eklenerek] → 186b.merdiven_izdusum(ad, w, aq, {'bir': ρ≡1}, "
                    "…) → 187c ζ defteri AYNEN (c_kesik = G1_proj re_bir/im_bir, "
                    "c_öz = OZ ρ≡1, c_ölç = K1)."),
        "surucu": ("189c_zincir.py: 184b2/185b/186b/187c modülleri importlib "
                   "ile yüklenir; yeni gaz adları için yalnız çalışma-anı "
                   "yaması: b185.kinematik (her iki b185 örneği: 185b ve "
                   "186b'nin içindeki) yeni adlarda z_<ad>.npy'den 185b'nin "
                   "Hkeskin kolu ile AYNI satırlarla (g=diff z, m=orta nokta, "
                   "L=mean log(m/2π)); çıktı dizinleri (S184 of 184b2, S185 of "
                   "185b) scratchpad/189'a yönlendirilir. __main__ blokları "
                   "(defter/ζ) sarmalayıcıda aynı işlem sırasıyla yeniden "
                   "kurulur. DOSYALAR DÜZENLENMEZ; mevcut 184-188 önbellekleri "
                   "okunur, ÜZERİNE YAZILMAZ."),
        "cikti": ("189/K1_<ad>.npz, 189/OZ_<ad>.npz, 189/G1_proj_<ad>.npz, "
                  "189/zincir_<ad>.json"),
        "makine_muhru": ("AYNI sarmalayıcı zinciri Hkeskin'e uygulanır (çıktılar "
                         "189/ altına): 189/K1_Hkeskin = 184/K1_Hkeskin, "
                         "189/OZ_Hkeskin = 185/OZ_Hkeskin, 189/G1_proj_Hkeskin "
                         "re_bir/im_bir = 186/G1_proj_Hkeskin, 185-defter 9 "
                         "satırı = 185/K1_faktorler (r, sig_r, F'ler, w'ler, "
                         "se_tablo), 187c Hkeskin satırları = 187/K2_zeta.json "
                         "— HEPSİ BİT-BİT (np.array_equal / ==). Tutmazsa "
                         "zincir ŞÜPHELİ, hüküm verilmez."),
    },
    "olcumler": {
        "bantlar": BANTLAR + ["HAVUZ"],
        "maskeler": ("185b/187c AYNEN: bant = (τ ≥ lo) & (τ < hi), τ = 184 "
                     "evreninin τ'su (K1_gercek); HAVUZ = τ ∈ [0.45, 0.86)"),
        "w": "w_X = Σ_b |c^ölç_X| / Σ_b ae (185b defter AYNEN)",
        "w_oz": "w^öz_X = Σ_b |c^öz_X| / Σ_b ae (185 OZ, ρ≡1)",
        "m": "m_X = w_X − w^öz_X",
        "r": "r_D = w_g / w_HD (toplam-oranı lehçesi)",
        "sig_r": ("184/185 konvansiyonu: σ_r = r·√((se_wg/wg)² + (se_wh/wh)²), "
                  "se'ler 8-blok loo jackknife √(7/8·Σ(θ_i−θ̄)²)"),
        "sigma_eps": ("185b oz_hesapla AYNEN: σ_ε = √(Σ(g−ḡ)²/N)/ḡ; jk se: loo "
                      "blok toplamları (sdg, sdg2) ile tam loo yeniden hesap"),
        "zeta": ("187c AYNEN: ζ_b = Σ_b Δ·conj(karışım)/Σ_b |karışım|², "
                 "Δ = c^ölç − c^kesik, karışım = c^kesik − c^öz; modül + açı "
                 "(jk; açı-se sarmalı); genlik-okuma 1 − m_ölç/m^kesik; "
                 "m_ölç, m^kesik"),
    },
    # ---------------- öngörü ----------------
    "ongoru": {
        "formul": ("m^kesik_Hk = m_Hk/(1 − z_1.00); m_HD = m^kesik_Hk·(1 − z_D); "
                   "r_D = w_g/(w^öz_Hk + m_HD); z_D = −Re Σ_{orta ≤ D} K_Hk "
                   "(188 ikiz haritası, satır 'B<bant>'); kaynak: "
                   "185/K1_faktorler.npz [wg, wh, wog, woh] + "
                   "188/harita_K_Hkeskin.npz"),
        "sutunlar": SUTUNLAR,
        "tablo_kalem_AYNEN": ONGORU_KALEM,
        "yeniden_hesap": None,
        "yeniden_hesap_esit_4hane": None,
        "HAVUZ_bilgi": None,
        "okuma": ("birincil öngörü YÖN ve KABA HIZ (sabit kinematikle 1.10'da "
                  "~%60, 1.20'de ~%100 kapanış); sayısal kıyas ikincil; derin "
                  "ikizin kinematiği de değişecek (F_KİN)."),
        "mekanizma": ("τ'≈1 çizgilerinin çekirdeği sin(πτ'(1+ds)) ≈ −π·ds: "
                      "aralık sapmasına NEGATİF GERİ-BESLEME ⇒ derin ikizin "
                      "σ_ε'si AZALMALI (gerçeğin 0.4092'sine doğru)."),
    },
    # ---------------- hipotezler ----------------
    "hipotezler": {
        "H-189a": {
            "ad": "ZARF = MERDİVEN DERİNLİĞİ (birincil)",
            "kosul": ("HAVUZ'da VE ≥7/8 bantta: r_1.00 < r_1.10 VE "
                      "r_1.10 ≤ r_1.20 + 2σ_Δ, σ_Δ = √(σ_r(1.10)² + σ_r(1.20)²) "
                      "(184 konvansiyonu, bağımsız yayılım); VE HAVUZ "
                      "f_1.10 = (r_1.10 − r_1.00)/(1 − r_1.00) ≥ 0.30. "
                      "r_1.00 = aynı sarmalayıcı zincirinin Hkeskin değeri."),
            "olum": ("HAVUZ f_1.10 < 0.15 YA DA r HAVUZ'da azalır: "
                     "r_1.10 < r_1.00 YA DA r_1.20 < r_1.10 − 2σ_Δ"),
            "hukum": "koşul → MÜHÜR; ölüm → ÖLDÜ; ikisi de değil → KAYIT (ara bölge)",
        },
        "H-189b": {
            "ad": "SABİT-KİNEMATİK ARİTMETİĞİ (ikincil, nicel)",
            "kosul": ("her D ∈ {1.10, 1.20} ve her bant: |r_D − r_D^öng| ≤ "
                      "max(2σ_r(D), 0.01); r_D^öng = tablo_kalem_AYNEN"),
            "hukum": ("D başına: 8/8 bant tutarsa TUTTU. İki derinlik de TUTTU → "
                      "MÜHÜR; en az biri tutmazsa ÖLDÜ (kurtarmasız). Her durumda "
                      "AYRIŞIM TABLOSU (KAYIT): Δw = w_HD^ölç − w_HD^öng = "
                      "(w^öz_HD − w^öz_Hk) [kinematik] + (m_HD − m_HD^öng) "
                      "[karışım]; kinematik payı = Δw^öz/Δw; ayrıca "
                      "r^kin = w_g/(w^öz_HD + m_HD^öng) (yalnız kinematik "
                      "değiştirilmiş öngörü)."),
        },
        "H-189c": {
            "ad": "GERİ-BESLEME KİNEMATİĞİ",
            "kosul": ("σ_ε(Hderin110) < σ_ε(Hkeskin) [aynı zincirde, 0.4313] VE "
                      "σ_ε(Hderin120) < σ_ε(Hderin110) (nokta değerler)"),
            "olum": "σ_ε derinlikle artar: σ_ε(110) > σ_ε(Hk) YA DA σ_ε(120) > σ_ε(110)",
            "kayit": ("gerçeğin 0.4092'sini aşarsa (altına inerse) KAYIT; her "
                      "farkın 2·jk-se'ye oranı ve 'gerçeğe doğru kapanan pay' "
                      "(σ_Hk − σ_D)/(σ_Hk − σ_g) raporlanır (hükmü değiştirmez)"),
            "hukum": "koşul → MÜHÜR; ölüm → ÖLDÜ; eşitlik → KAYIT",
        },
        "H-189d": {
            "ad": "HARİTA GEÇERLİLİĞİ (kapı-niteliğinde)",
            "kosul": ("HAVUZ δ_D = | |ζ_HD| − z_D | / z_D ≤ 0.10, her D ∈ "
                      "{1.10, 1.20}; |ζ_HD| = 187c modülü (187 raporunun ζ "
                      "sütunu), z_D = −Re Σ_{orta≤D} K_Hk(HAVUZ) (z_1.10 = "
                      "0.1794, z_1.20 = 0.2430). Açı 180°±15° ayrıca raporlanır; "
                      "İKİNCİL: −Re ζ_HD ile aynı oran; bant bazında δ tablosu."),
            "hukum": "iki derinlik de ≤ %10 → MÜHÜR; en az biri > %10 → ÖLDÜ",
        },
        "K3": {
            "ad": "ZARFIN KİMLİĞİ (KAYIT)",
            "fit": ("r_D(τ) 8 bant (D = 1.00, 1.10, 1.20) → 1 − c·τ^α; 184c "
                    "makinesi AYNEN (184c.bant_defteri: maske (lo, hi], τ̄ "
                    "ae-ağırlıklı, se_r; 184c.fit + H_Z2: curve_fit sigma=se_r "
                    "absolute_sigma, p0 = [0.3, 3.0], sınır ([0, 0.1], [10, 20])). "
                    "Makine mührü: D = 1.00 (Hkeskin) → c = 0.149, α = 1.30, "
                    "χ²/dof = 0.34 (184). (c, α) se: fit-perr + loo-jackknife "
                    "yeniden fit."),
            "D_sonsuz_A_birincil": (
                "KİNEMATİK 1.20'DE DONMUŞ + GERÇEK-DERİNLİK İPTALİ: "
                "r_∞ = w_g / (w^öz_H120 + m^kesik_H120·(1 − Γ_g)), m^kesik_H120 = "
                "187c m_kesik (derin ikizin KENDİ pencere-içi karışımı), Γ_g = "
                "gerçeğin tam-derinlik iptal payı 187c genlik-okuması "
                "(1 − m_g/m^kesik_g; HAVUZ 0.3335); gerekçe 188: derinlik-eşli "
                "ζ iki denizde aynı. Duyarlılık: Γ_g yerine |ζ_g| (0.3287). "
                "jk se: ortak 8-blok loo."),
            "D_sonsuz_B_ikincil": (
                "KİNEMATİK EĞİLİMİ SÜRER: w^öz_HD ve m^kesik_HD, üç derinlikte "
                "(1.00, 1.10, 1.20) iptal payı Γ_HD'ye (187c genlik-okuması) "
                "karşı doğrusal en-küçük-kare; Γ = Γ_g'de değerlendirilir; "
                "r_∞^B = w_g/(w^öz(Γ_g) + m^kesik(Γ_g)·(1 − Γ_g))."),
            "paylar": ("zarf (1 − r_1.00) içinde DERİNLİK payı P_der = "
                       "(r_∞ − r_1.00)/(1 − r_1.00); KİNEMATİK payı = 1 − P_der "
                       "= (1 − r_∞)/(1 − r_1.00) (r_∞ > 1 ⇒ negatif: ters yönlü "
                       "kinematik fark). Ölçülen (modelsiz) f_1.10, f_1.20 ayrıca. "
                       "Kinematik fark ayrışımı (A): w_H∞ − w_g = (w^öz_H120 − "
                       "w^öz_g) + (m^kesik_H120 − m^kesik_g)(1 − Γ_g)."),
            "turetilen_olculen": ("öngörü tablosu TÜRETİLEN (sabit-kinematik "
                                  "aritmetiği, 188 haritası); r_D, σ_ε, ζ_HD, "
                                  "(c, α) ÖLÇÜLEN; D→∞ MODEL (A/B açıkça)."),
            "hukum": "KAYIT (eşik yok)",
        },
    },
    "jackknife": "8-blok loo (184-187 AYNEN), se = √(7/8·Σ(θ_i−θ̄)²)",
    "kurallar": ("TEK DALGA; ölümler kurtarmasız; eşik gevşetme yok; ölçülemeyen "
                 "'erişilemedi'; git'e dokunulmaz; sonuç ORTAK TEFTİŞE."),
    "zaman": None,
    "sha256": None,
}

if __name__ == "__main__":
    yh, havuz = yeniden_hesap()
    esit = True
    fark_maks = 0.0
    for ad in BANTLAR:
        k = ONGORU_KALEM[ad]
        y = yh[ad]
        for j in range(6):
            fark_maks = max(fark_maks, abs(round(y[j], 4) - k[j]))
            if round(y[j], 4) != k[j]:
                esit = False
        for j in (6, 7):
            if round(y[j]) != k[j]:
                esit = False
    ONKAYIT["ongoru"]["yeniden_hesap"] = {ad: [float(x) for x in yh[ad]]
                                         for ad in BANTLAR}
    ONKAYIT["ongoru"]["yeniden_hesap_esit_4hane"] = bool(esit)
    ONKAYIT["ongoru"]["HAVUZ_bilgi"] = havuz
    ONKAYIT["kalem_sha256"] = hashlib.sha256(
        (QM / "KALEM_DERIN_IKIZ_23EYL2026.md").read_bytes()).hexdigest()
    sha = hashlib.sha256(
        (QM / "189_configs" / "189a_onkayit.py").read_bytes()).hexdigest()
    ONKAYIT["sha256"] = sha
    ONKAYIT["zaman"] = time.strftime("%a %b %d %H:%M:%S %z %Y")
    yol = S189 / "ONKAYIT_189.json"
    json.dump(ONKAYIT, open(yol, "w"), indent=1, ensure_ascii=False)
    print(f"öngörü tablosu kaynağından yeniden hesap = kalem (4 hane): "
          f"{'EŞİT' if esit else 'EŞİT DEĞİL'} (maks fark {fark_maks:.1e})")
    print("HAVUZ (bilgi): " + ", ".join(f"{k}={v:.4f}" for k, v in havuz.items()))
    print(f"ONKAYIT_189.json yazıldı -> {yol}")
    print(f"  sha256 = {sha}")
    print(f"  damga  = {ONKAYIT['zaman']}")
