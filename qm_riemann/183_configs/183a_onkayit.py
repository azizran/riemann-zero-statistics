# -*- coding: utf-8 -*-
"""
183a — K0: DONMUŞ ÖN-KAYIT (KURAL-ÖNCE, VERİDEN ÖNCE)
======================================================
Bu betik HİÇBİR YENİ VERİ ÜRETMEZ. Yalnız
  (1) 182'nin ölçülmüş ÇAPA değerlerini diskten okur ve dondurur,
  (2) üç-nokta ayrım formüllerini ve para birimini yazar,
  (3) çatal tablo şablonunu yazar,
  (4) H-İ1'in ÖLÜM/YAŞAMA koşullarını SAYIYLA yazar,
  (5) derinlik protokolünü ve hangi kapıların uygulanmayacağını yazar,
  (6) K2 surrogate protokolünü (SUR-A / SUR-B / SUR-C) yazar,
  (7) kendi sha256'sını + KALEM'in sha256'sını + damgayı basar.

EŞİK İCADI YOKTUR: her tolerans ya 182/176/180'den DEVRALINIR (sd/se
olarak) ya da 182'nin ölçülmüş değerine ORAN olarak tanımlanır.

Kullanım: 183a_onkayit.py
"""
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S183 = SCR / "183"
KALEM = QM / "KALEM_INSA_KILIDI_08EYL2026.md"
BEN = Path(__file__).resolve()


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    S183.mkdir(parents=True, exist_ok=True)
    damga = subprocess.run(["date"], capture_output=True,
                           text=True).stdout.strip()
    print("=" * 78)
    print("183a — K0 DONMUŞ ÖN-KAYIT (İNŞA-KİLİDİ)")
    print("=" * 78)
    print(f"  damga        : {damga}")
    print(f"  183a sha256  : {sha(BEN)}")
    print(f"  KALEM sha256 : {sha(KALEM)}")

    # ================================================================
    # 1. ÇAPALAR — 182/176/180/174'ten OKUNUR (yeni ölçüm yok)
    # ================================================================
    AN = [json.load(open(SCR / "182" / f"ANATOMI_VF{i}.json"))
          for i in range(1, 9)]
    K1 = {g: json.load(open(SCR / ("176" if g.startswith("VF") else "174")
                            / f"K1_{g}.json"))
          for g in [f"VF{i}" for i in range(1, 9)] + ["Hkeskin", "son"]}
    K2B = {g: json.load(open(SCR / "169" / f"K2b_{g}.json"))
           for g in [f"VF{i}" for i in range(1, 9)]}

    r_vf = np.array([a["korr_XaXb"] for a in AN])
    r_sur = np.array([a["korr_XaXb_sur"] for a in AN])
    kap_vf = np.array([K1[f"VF{i}"]["eta"]["kov_art_x"]
                       / K1[f"VF{i}"]["eta"]["Var"] for i in range(1, 9)])
    R_vf = np.array([K1[f"VF{i}"]["eta"]["R"] for i in range(1, 9)])
    rho3_vf = np.array([K2B[f"VF{i}"]["ortalama"]["111_hepsi"]
                        for i in range(1, 9)])

    CAPA = dict(
        # --- 182 §3.4 / §3.5 ---
        korr_XaXb_VF=float(r_vf.mean()), korr_XaXb_VF_sd=float(r_vf.std(ddof=1)),
        korr_XaXb_VF1=float(r_vf[0]),
        korr_XaXb_SUR=float(r_sur.mean()), korr_XaXb_SUR_sd=float(r_sur.std(ddof=1)),
        kovoran_VF=float(kap_vf.mean()), kovoran_VF_sd=float(kap_vf.std(ddof=1)),
        kovoran_VF1=float(kap_vf[0]),
        R_eta_VF=float(R_vf.mean()), R_eta_VF_sd=float(R_vf.std(ddof=1)),
        R_eta_VF1=float(R_vf[0]),
        rho3_VF=float(rho3_vf.mean()), rho3_VF_sd=float(rho3_vf.std(ddof=1)),
        rho3_VF1=float(rho3_vf[0]),
        rho3_sur_VF=float(np.mean([a["rho_sur"]["hepsi"] for a in AN])),
        # --- gerçek/ikiz (174/176/180) ---
        R_eta_Hkeskin=float(K1["Hkeskin"]["eta"]["R"]),
        R_eta_son=float(K1["son"]["eta"]["R"]),
        kovoran_Hkeskin=float(K1["Hkeskin"]["eta"]["kov_art_x"]
                              / K1["Hkeskin"]["eta"]["Var"]),
        kovoran_son=float(K1["son"]["eta"]["kov_art_x"]
                          / K1["son"]["eta"]["Var"]),
        kov_art_eta_Hkeskin=float(K1["Hkeskin"]["eta"]["kov_art_x"]),
        kov_art_eta_son=float(K1["son"]["eta"]["kov_art_x"]),
        kov_art_eta_VF=float(np.mean([K1[f"VF{i}"]["eta"]["kov_art_x"]
                                      for i in range(1, 9)])),
        Var_eta_Hkeskin=float(K1["Hkeskin"]["eta"]["Var"]),
        Var_eta_son=float(K1["son"]["eta"]["Var"]),
        Var_eta_VF=float(np.mean([K1[f"VF{i}"]["eta"]["Var"]
                                  for i in range(1, 9)])),
        # Kov(η_çizgi, η_artık) ≡ P − Var(η_çizgi)   (176 §2b'nin −0.0557'si)
        kov_ciz_art_Hkeskin=float(K1["Hkeskin"]["eta"]["P"]
                                  - K1["Hkeskin"]["eta"]["Var_cizgi"]),
        kov_ciz_art_son=float(K1["son"]["eta"]["P"]
                              - K1["son"]["eta"]["Var_cizgi"]),
        kov_ciz_art_VF=float(np.mean([K1[f"VF{i}"]["eta"]["P"]
                                      - K1[f"VF{i}"]["eta"]["Var_cizgi"]
                                      for i in range(1, 9)])),
        m3_cizgi_Hkeskin=float(K1["Hkeskin"]["m3"]["cizgi"]),
        m3_cizgi_son=float(K1["son"]["m3"]["cizgi"]),
        m3_cizgi_VF=float(np.mean([K1[f"VF{i}"]["m3"]["cizgi"]
                                   for i in range(1, 9)])),
        m3_cizgi_VF_sd=float(np.std([K1[f"VF{i}"]["m3"]["cizgi"]
                                     for i in range(1, 9)], ddof=1)),
        # 176 §2b tablosu (μ̂²_E) — 176 raporundan devralınan ölçülmüş satır
        mu2E_Hkeskin=0.056843, mu2E_son=0.040498,
        mu2E_VF=float(np.mean([json.load(open(SCR / "176" / f"G_VF{i}.json"))
                               ["E"]["mu2"] for i in range(1, 9)])),
        rho_var_E_VF=float(np.mean([json.load(open(SCR / "176" / f"G_VF{i}.json"))
                                    ["E"]["rho_var"] for i in range(1, 9)])),
        # 180a devralınan
        GAUSS3=0.5079490874739278, GAUSS1=0.7978845608028654,
        rho3_Hkeskin=0.5284, rho3_son=0.5356,
    )
    print("\n  --- ÇAPALAR (182/180/176/174'ten okundu; yeni ölçüm yok) ---")
    for k in sorted(CAPA):
        print(f"    {k:26s} = {CAPA[k]:+.6f}")

    # ================================================================
    # 2. ÜÇ-NOKTA AYRIM FORMÜLLERİ + PARA BİRİMİ
    # ================================================================
    AYRIM = {
        "para_birimi": (
            "Her gözlenebilir O için ayrım DOĞRUSALDIR ve iki para "
            "biriminde AYRI AYRI yazılır: (H) HAM = niceliğin kendi "
            "doğal birimi (kovaryans, moment, ρ) ; (O) ORAN = "
            "boyutsuz oran (Kov/Var, R). İkisi arasında dönüşüm "
            "gazın kendi Var(η)'sıyla yapılır ve BU DÖNÜŞÜM GAZDAN "
            "GAZA DEĞİŞİR — ayrımın para-birimi-değişmez olup "
            "olmadığı K3'te ÖLÇÜLÜR, varsayılmaz."),
        "alt_sinir_gercek_kilit": "GK_alt(O) = O(ikiz) − O(VF)",
        "insa_kilidi":            "IK(O)     = O(VF)   − O(SUR)",
        "ust_kol_gercek_kilit":   "GK_ust(O) = O(ikiz) − O(SUR)",
        "ozdeslik":               "GK_alt + IK ≡ GK_ust  (doğrusal, tanım gereği)",
        "catal":                  "ÇATAL(O) = [GK_alt(O), GK_ust(O)]",
        "ikiz": "Hkeskin (birincil). 'son' (gerçek gaz) ikinci satır olarak "
                "yazılır; hiçbir gerçek-gaz niceliği YENİDEN ölçülmez.",
        "VF":   "n=8 VF ailesinin ortalaması (182 §3.4/§2).",
        "SUR":  "alan-surrogate tabanı (K2; SUR-A/B/C protokolü).",
    }
    CATAL_SABLON = ["O", "O(ikiz)", "O(VF)±se", "O(SUR)±se",
                    "GK_alt = ikiz−VF", "IK = VF−SUR", "GK_ust = ikiz−SUR",
                    "ÇATAL [alt, üst]", "IK/GK_ust  (inşa payı)"]

    # ================================================================
    # 3. DERİNLİK PROTOKOLÜ (K1)
    # ================================================================
    DERINLIK = {
        "tohum": 1, "zarf": "Hkeskin (a_q = 1/(πk√q), τ_q ≤ 1.00, 15450 çizgi)",
        "cekirdek": "164_insa.coz_sadakatli AYNEN — TEK SATIR DEĞİŞMEZ. "
                    "Derinlik, çekirdeğin çağırdığı SSp_par'ı saran bir "
                    "sayaç-sarmalayıcıyla verilir (183b).",
        "tanim": (
            "derinlik d := ızgara-braketi başlangıcına uygulanan KORUMALI "
            "NEWTON GÜNCELLEMESİ sayısı. d=0: yalnız ızgara braketi + "
            "doğrusal ara değer (GERİ-BESLEME YOK). d=∞: çekirdeğin kendi "
            "yakınsaması (= mevcut VF1)."),
        "merdiven": [0, 1, 2, 4, "inf"],
        "tek_kosu": ("Bütün derinlikler TEK çözücü koşusundan alınır: "
                     "sarmalayıcı her Newton yinelemesinin GİRİŞİNDE z'yi "
                     "anlık-görüntüler; d yinelemesinin girişindeki z, tam "
                     "olarak d güncelleme almış z'dir. Böylece d={0,1,2,4} "
                     "ve d=∞ AYNI koşunun iç durumlarıdır — tohum, ızgara, "
                     "braket, faz dizisi bit-bit ORTAKTIR."),
        "R_KAPISI": ("d=∞ çıktısı 176/z_VF1.npy ile BİT-BİT aynı olmalı "
                     "(maks|Δz| = 0.0). Olmazsa 183b HÜKÜMSÜZDÜR."),
        "kapilar": {
            "uygulanan": ["G1 (sha(A) bit-bit)", "G2 (φ≡0 → 0.0/0.0)",
                          "G3 (ilk-kök hücre benzersizliği — braket ortak)"],
            "UYGULANMAYAN": ["G4 (maks|F| ≤ 1e−8)", "G5 (sıralılık TAM)"],
            "gerekce": ("G4 ve G5 sığ derinlikte YAPISAL OLARAK tutmaz "
                        "(yakınsama yoktur). Sığ-derinlik gazları bu yüzden "
                        "İNŞA GAZI DEĞİL, **TEŞHİS GAZI**dır; deftere "
                        "'TEŞHİS' etiketiyle girer ve 176/180/182'nin "
                        "hiçbir hükmünü değiştiremez."),
        },
        "olcum_zinciri": ("155_kos.veri_yukle diziyi SIRALAR (np.sort) — "
                          "sığ gazlar bataryaya SIRALANMIŞ girer; ham "
                          "sıralılık ihlali sayısı ayrıca kaydedilir."),
        "gozlenebilirler": ["korr(Xa,Xb) [183c]", "R_η [176d aynen]",
                            "Kov(η_artık,η)/Var(η) ≡ 1−R_η", "ρ₃ [180d aynen]"],
    }

    # ================================================================
    # 4. H-İ1 — ÖLÜM / YAŞAMA (SAYIYLA, VERİDEN ÖNCE)
    # ================================================================
    # Referans: aynı tohumun (VF1) yakınsamış değerleri — 182'de ölçülmüş.
    # Tolerans: 182'nin AİLE saçılımı (sd, n=8) — devralınır, icat edilmez.
    HI1 = {
        "onerme": ("İnşa-kilidi öz-tutarlılık DERİNLİĞİYLE birikir: "
                   "korr(Xa,Xb) ve Kov(η_artık,η)/Var(η) derinlik-0'da "
                   "~0'dan yakınsamış değere TEKDÜZE tırmanır."),
        "referans": {"korr_XaXb(inf)": CAPA["korr_XaXb_VF1"],
                     "kovoran(inf)": CAPA["kovoran_VF1"]},
        "pay": "f(d) := O(d) / O(∞)   (her iki gözlenebilir için ayrı)",
        "tol": {"korr_XaXb": CAPA["korr_XaXb_VF_sd"],
                "kovoran": CAPA["kovoran_VF_sd"],
                "kaynak": "182 §3.4 aile sd'si (n=8) — DEVRALINDI"},
        "YASAMA": ("f(0) ≤ 0.25 HER İKİ gözlenebilirde  VE  merdiven "
                   "TEKDÜZE (her ardışık çiftte O(d_{i+1}) ≥ O(d_i) − tol) "
                   "HER İKİ gözlenebilirde"),
        "OLUM": ("f(0) ≥ 0.75  gözlenebilirlerden EN AZ BİRİNDE "
                 "⇒ mekanizma yineleme derinliği DEĞİL, kök-seçimi/"
                 "sıralılıktır"),
        "HUKUMSUZ": "arada kalan her şey (kurtarma yok)",
        "OLUM_SONRASI": ("ÖLÜM hâlinde ikinci tur: kök-kuralı varyantı — "
                         "164_insa.coz_enyakin ('en yakın kök') ile aynı "
                         "tohum/zarf, d=∞; korr(Xa,Xb) ve Kov-oranı "
                         "yeniden ölçülür. Bu tur da 183'ün içinde kalır; "
                         "başka hiçbir eşik değişmez."),
    }

    # ================================================================
    # 5. K2 — SURROGATE PROTOKOLÜ (yeniden çözüm YOK)
    # ================================================================
    SUR = {
        "ortak_islem": ("FAZ-RASTGELELEŞTİRME: amp_q → |amp_q|·e^{iψ_q}, "
                        "ψ ~ U(0,2π) bağımsız. Güç tayfı BİT-BİT korunur; "
                        "siteler ve frekanslar AYNEN kalır; SIFIR YENİDEN "
                        "ÇÖZÜLMEZ."),
        "SUR_A": {
            "makine": "182g_anatomi.py (H-G1 makinesi) — AYNEN, diskten",
            "kapsam": ["korr(Xa°,Xb°)", "ρ°(E)", "ρ°₃ (üç-bacak)"],
            "n": 8, "kaynak": "182/ANATOMI_VF1..8.json (YENİDEN KOŞULMAZ)",
        },
        "SUR_B": {
            "makine": "162_cekirdek (174b'nin makinesi) — 183e",
            "tanim": ("η° := cizgi_kur(m, |c_η|·e^{iψ}, w) + η_artık(VF)  "
                      "→ c° := cizgi_cikar(m, η°, w)  → "
                      "R°_η := Σ|c°|²/2 / Var(η°) ; Kov-oranı° := 1 − R°_η. "
                      "Aynı işlem X̃ kanalında; m3°_çizgi := m3(η°_çiz, X̃°_çiz)."),
            "kapsam": ["R_η", "Kov(η_artık,η)/Var(η)", "Kov(η_çiz,η_art)",
                       "m3_çizgi"],
            "n": 4, "gaz": ["VF1", "VF2", "VF3", "VF4"],
            "yapisal_ongoru": ("m3°_çizgi ≡ 0 (rastgele fazlı çizgi alanının "
                               "üçüncü momenti özdeş sıfırdır — 176a/F5); "
                               "ölçüm yalnız SONLU-ÖRNEK saçılımını verir."),
        },
        "SUR_C": {
            "nicelik": "μ̂²_E",
            "TURETIM": ("Surrogate model alanı E° ölçülen e1'den BAĞIMSIZDIR "
                        "⇒ kor_E° = 0 + O(N^{-1/2}), V_R° = V_O + V_M ⇒ "
                        "Kv° = ½(V_O+V_M−V_R°) = 0 ⇒ S° = 0 ⇒ "
                        "μ̂²_E° = (0 − V_M)/V_O = −ρ_var(E)."),
            "DEJENERE": ("Bu bir eşik değil bir ÖZDEŞLİKTİR: surrogate, "
                         "μ̂²_E'nin ölçtüğü model↔ölçüm eşleşmesini TAMAMEN "
                         "yok eder. Dolayısıyla μ̂²_E'nin ÇATALI ön-kayıtta "
                         "**DEJENERE** etiketlidir ve üst kolu bir bilgi "
                         "taşımaz. Yine de sayıyla yazılır; gizlenmez."),
            "sinav": "183e, korr(E°, E) ≈ 0'ı doğrudan ölçüp basar.",
        },
    }

    # ================================================================
    # 6. K3 / K4 KURALLARI
    # ================================================================
    K3 = {
        "isaret_turnusolu": (
            "182 §3.5 son paragrafı Kov(η_artık,η)(Hkeskin)'i −0.0557 diye "
            "yazar. DENETİM (183a, veriden önce): 174/K1_Hkeskin.json'da "
            f"kov_art_x = {CAPA['kov_art_eta_Hkeskin']:.6f} ve "
            f"P − Var(η_çizgi) = {CAPA['kov_ciz_art_Hkeskin']:.6f}. "
            "Yani −0.0557 **Kov(η_çizgi, η_artık)**'tır, Kov(η_artık, η) "
            "DEĞİL (176 §2b'nin son sütunu). 183 turnusolü İKİ nesneyi "
            "AYRI AYRI ve DOĞRU ADLARIYLA yürütür; 182'nin etiket hatası "
            "raporda ERRATA olarak yazılır."),
        "turnusol_testi": (
            "ORAN ayrımı (κ = Kov/Var) ile HAM ayrımı (Kov) aynı bölüşümü "
            "veriyor mu? GK_alt^ham / GK_alt^oran bir Var(η) değeri verir; "
            "bu değer Var_η(ikiz) ya da Var_η(VF) ile uyuşuyor mu? "
            "Uyuşmuyorsa AYRIM PARA-BİRİMİ-DEĞİŞMEZ DEĞİLDİR ve bu aynen "
            "yazılır (doğrusal toplanmıyor olabilir; o da bulgudur)."),
    }
    K4 = {
        "soru": ("180/181/182 defterinin ÇAPASIZ ZARF/KİLİT payları "
                 "düzeltilmiş (surrogate) tabanla değişir mi?"),
        "yapisal_gerekce_sinavi": (
            "ZARF payı iki KARIŞIK ailenin (VS ve VF) FARKIdır: "
            "ZARF = ⟨N⟩_VS − ⟨N⟩_VF. İnşa-kilidi her iki ailede de "
            "AYNI çözücüyle üretildiğinden ORTAK olabilir ve farkta "
            "SADELEŞİR. Bu bir varsayım değil: 183g bunu ÖLÇER — "
            "VS ve VF ailelerinin inşa-kilidi imzaları (korr(Xa,Xb), "
            "Kov-oranı) karşılaştırılır. Ortaksa (|Δ| ≤ 2·birleşik se) "
            "ZARF payı DEĞİŞMEZ ve gerekçesi ölçülmüş olur; ortak "
            "değilse fark sayıyla yazılır."),
        "n_VS": 4, "n_VF": 8,
    }

    ON = dict(gorev=183, ad="INSA-KILIDI", damga=damga,
              sha_183a=sha(BEN), sha_KALEM=sha(KALEM),
              python=sys.version.split()[0], numpy=np.__version__,
              CAPA=CAPA, AYRIM=AYRIM, CATAL_SABLON=CATAL_SABLON,
              DERINLIK=DERINLIK, HI1=HI1, SUR=SUR, K3=K3, K4=K4)
    p = S183 / "ONKAYIT_183.json"
    p.write_text(json.dumps(ON, indent=1, ensure_ascii=False))

    print("\n  --- ÜÇ-NOKTA AYRIMI ---")
    for k, v in AYRIM.items():
        print(f"    {k}: {v}")
    print("\n  --- H-İ1 ---")
    for k, v in HI1.items():
        print(f"    {k}: {v}")
    print("\n  --- DERİNLİK PROTOKOLÜ ---")
    for k, v in DERINLIK.items():
        print(f"    {k}: {v}")
    print("\n  --- SURROGATE PROTOKOLÜ ---")
    for k, v in SUR.items():
        print(f"    {k}: {v}")
    print("\n  --- K3 / K4 ---")
    for k, v in list(K3.items()) + list(K4.items()):
        print(f"    {k}: {v}")
    print(f"\n  -> {p}")
    print(f"  ({time.strftime('%H:%M:%S')})")
    return ON


if __name__ == "__main__":
    main()
