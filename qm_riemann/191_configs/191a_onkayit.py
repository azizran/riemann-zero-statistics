# -*- coding: utf-8 -*-
"""
191a — K0: ÖN-KAYIT (DERİN İKİZ τ≤1.30 — ZARF TAM MI KAPANIR?)
================================================================
KALEM_DERIN_IKIZ_130_23EYL2026.md AYNEN (commit 4efc958, Wed Sep 23
18:31:28 +0300 2026 — İNŞA'dan (K1, 191b_insa.py, grid-sınavı logu
18:45, inşa logu bitiş 23:36) ÖNCE push'landı). Bu betik ZİNCİRDEN
(K2, 191c_zincir.py) ÖNCE koşar ve ölçümden önce donan her şeyi
scratchpad/191/ONKAYIT_191.json'a yazar: ön-mühür öngörü bandı
(f_1.30, σ_ε(1.30)), harita-doğrusal yasa sabitleri, H-191a/b/c eşikleri,
K3 kimlik-kaydı tanımları, K4 güç-analizi kuralı.

Kaynağından yeniden hesap (K1 DEĞİL — yalnız öngörünün girdisi olan
188/189 sayılarının tutarlılık denetimi): z_D (gerçek-kinematik harita,
188b harita_K_gercek.npz HAVUZ) ve HAVUZ f_1.10/f_1.20 (189/HUKUM_189.json)
kaynaktan okunur, KALEM'in yazdığı sayılarla (4 hane) karşılaştırılır.

Çıktı: scratchpad/191/ONKAYIT_191.json (sha256 = bu betiğin sha'sı + damga).
"""
import hashlib
import json
import subprocess
import time
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S188, S189, S191 = SCR / "188", SCR / "189", SCR / "191"
S191.mkdir(exist_ok=True)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def yeniden_hesap():
    """z_D (gerçek-kinematik harita, HAVUZ) + harita-doğrusal yasa
    sabitlerini (f_1.10/Δz_1.10, f_1.20/Δz_1.20) kaynağından yeniden
    hesapla; KALEM'in yazdığı sayılarla (0.0758/0.1717/0.2261/0.2566,
    tam 0.3287, 3.38/3.39) karşılaştır."""
    hg = np.load(S188 / "harita_K_gercek.npz", allow_pickle=True)
    et = list(hg["etiket"])
    i = et.index("HAVUZ")
    orta, K = hg["orta"], hg["K"]
    z = {D: float(-np.real(K[i, orta <= D + 1e-12].sum()))
         for D in (1.00, 1.10, 1.20, 1.30)}
    zg_full = np.load(S189 / "zincir_gercek.npz")
    z_inf = float(zg_full["z_tam"][8, 0])          # |ζ_g| HAVUZ (187c)

    h189 = json.load(open(S189 / "HUKUM_189.json"))
    f110 = h189["H-189a"]["f110_havuz"]
    f120 = h189["H-189a"]["f120_havuz"]
    dz110 = z[1.10] - z[1.00]
    dz120 = z[1.20] - z[1.00]
    dz130 = z[1.30] - z[1.00]
    dz_inf = z_inf - z[1.00]
    slope110 = f110 / dz110
    slope120 = f120 / dz120
    slope_c = 0.5 * (slope110 + slope120)
    f130_pred = slope_c * dz130
    f_inf_pred = slope_c * dz_inf

    # σ_ε kapanış payı g_D benzer yasa (bilgi amaçlı, aynı tutarlılık)
    sHk = h189["H-189c"]["sigma_eps"]["1.00"]
    s110 = h189["H-189c"]["sigma_eps"]["1.10"]
    s120 = h189["H-189c"]["sigma_eps"]["1.20"]
    sg = h189["H-189c"]["gercek"][0]
    g110 = (sHk - s110) / (sHk - sg)
    g120 = (sHk - s120) / (sHk - sg)
    gslope110 = g110 / dz110
    gslope120 = g120 / dz120
    g130_lo = min(gslope110, gslope120) * dz130
    g130_hi = max(gslope110, gslope120) * dz130
    sig130_lo = sHk - g130_hi * (sHk - sg)
    sig130_hi = sHk - g130_lo * (sHk - sg)

    return dict(
        z_D={f"{d:.2f}": v for d, v in z.items()}, z_inf=z_inf,
        Delta_z_D={"1.10": dz110, "1.20": dz120, "1.30": dz130,
                   "inf": dz_inf},
        f110_havuz=f110, f120_havuz=f120,
        yasa_egim={"1.10": slope110, "1.20": slope120, "merkez": slope_c},
        f130_pred_kaynak=f130_pred, f_inf_pred_kaynak=f_inf_pred,
        sigma_yasa_egim={"1.10": gslope110, "1.20": gslope120},
        sigma130_pred_kaynak_bant=[sig130_lo, sig130_hi],
    )


# ============================== ÖN-MÜHÜR (KALEM'DEN AYNEN) ================
ONGORU = {
    "z_D_gercek_kinematik": {"1.00": 0.0758, "1.10": 0.1717, "1.20": 0.2261,
                              "1.30": 0.2566, "tam_187c": 0.3287},
    "Delta_z_D": {"1.00": 0.0959, "1.10": 0.1503, "1.20": 0.1808,
                  "inf": 0.2529},
    "harita_dogrusal_yasa": {
        "tanim": "f_D / Δz_D  (HAVUZ; D=1.10, 1.20'de ölçülmüş f'lerle)",
        "1.10": 3.38, "1.20": 3.39, "kullanilan_merkez": 3.385,
    },
    "f_1.30": {"merkez": 0.612, "bant": [0.58, 0.64],
               "formul": "3.385 × 0.1808", "capraz_kontrol_geometrik": 0.615},
    "sigma_eps_1.30": {"merkez": 0.4150, "bant": [0.4140, 0.4160],
                       "formul": "g = (4.05…4.18) × 0.1808 = 0.732…0.756; "
                                 "σ_D = σ_Hk − g·(σ_Hk − σ_gercek)"},
    "sonsuz_derinlik_ima": {
        "f_inf": 0.86, "yorum": "zarfın ~%14'ü derinlikle KAPANMAZ "
                                "(derinlik-dışı kalıntı); 189 A/B modelleri "
                                "~%100 diyordu (P_der A %134, B %102)",
        "g_inf_bant": [1.02, 1.06],
        "g_inf_yorum": "aralık saçılımı derinlikle gerçeğe TAM yakınsar",
    },
}

HIPOTEZLER = {
    "H-191a": {
        "ad": "ZARF SÜRER (birincil, yön)",
        "kosul": ("HAVUZ f_1.30 > f_1.20 + 2·se(f) VE 8 bantta "
                  "r_1.20 ≤ r_1.30 + 2σ (σ = √(σ_r(1.20)²+σ_r(1.30)²), "
                  "184 konvansiyonu, bağımsız yayılım)"),
        "olum": "HAVUZ f_1.30 ≤ f_1.20 + 2·se(f)  (duraklama/geri dönüş)",
        "hukum": "koşul → MÜHÜR; ölüm → ÖLDÜ; ikisi de değil → KAYIT",
    },
    "H-191b": {
        "ad": "HARİTA-DOĞRUSAL YASA (ikincil, nicel)",
        "kosul": "f_1.30 ∈ [0.58, 0.64]",
        "yorum_tuttu": ("yasa üç derinlikte (1.10,1.20,1.30) geçerli ⇒ "
                        "f_∞ ≈ 0.86 okuması (zarfın ~%14'ü kalıntı) güç "
                        "kazanır"),
        "yorum_ustunde": "f_1.30 > 0.64 ⇒ yasa bozulur, tam kapanış "
                         "yolunda (189 A/B lehine)",
        "yorum_altinda": "f_1.30 < 0.58 ⇒ yasa bozulur, erken doyum",
        "hukum": "bant içi → MÜHÜR; dışı → ÖLDÜ (kurtarma yok)",
    },
    "H-191c": {
        "ad": "σ_ε YAKINSAMASI",
        "kosul": ("σ_ε(1.30) < σ_ε(1.20) − 2·se VE σ_ε(1.30) ≥ 0.40921 − "
                  "2·se (aşma yok); nicel bant [0.4140, 0.4160]"),
        "olum": ("σ_ε(1.30) ≥ σ_ε(1.20) − 2·se YA DA σ_ε(1.30) < "
                 "0.40921 − 2·se (gerçeği aşar) YA DA bant dışı"),
        "hukum": "koşul VE bant içi → MÜHÜR; ölüm → ÖLDÜ; ara → KAYIT",
    },
    "K3": {
        "ad": "KİMLİK KAYDI (hüküm yok, KAYIT)",
        "fit": ("dört derinlikte (1.00,1.10,1.20,1.30) r_D(τ) = 1−c·τ^α; "
                "184c makinesi AYNEN (189d.fit184c ile aynı çağrı: "
                "bant_defteri + H_Z2 curve_fit + loo-jk yeniden fit)."),
        "dogrusallik": ("f_D ve σ-kapanış payı g_D'nin Δz_D'ye karşı "
                        "doğrusallığı: dört nokta (D=1.00 orijin, f=g=0), "
                        "eğim ± se; hem orijinden-geçen (tek parametre) "
                        "hem serbest (iki parametre) fit."),
        "sonsuz_derinlik": ("yasa (orijinden geçen fit) tutarsa f_∞ = "
                            "eğim × Δz_∞(0.2529) ± se; 'zarfın derinlik "
                            "payı' = f_∞, 'kalıntı payı' = 1 − f_∞; "
                            "189'un D→∞ modelleri A (P_der %134) ve B "
                            "(P_der %102) ile yan yana konur."),
    },
}

K4_KURAL = {
    "ad": "KİLİT GÜÇ ANALİZİ (191e_guc.py; deep-twin KİLİT ölçümünden ÖNCE "
          "hesaplanır)",
    "girdi_180": {
        "kaynak": "180_nedensel_defter_RAPOR.md (T4 log-ayrışım + hüküm ii)",
        "Delta_log_M_son_Hk": 0.057542,
        "Delta_log_Q_E_son_Hk": -0.044475,
        "not_se_ham_fark": ("180a'nın kendi formülü: ΔN 'tohumsuz' (tek "
                            "gaz çifti, tohum yok) — ham farkın kendi "
                            "ayrı se'si YOKTUR (se(KİLİT_N)=se(ZARF_N) "
                            "180'in kuralı); scratchpad/180 (166-183) bu "
                            "oturumda TEMİZLENMİŞ, yeniden hesaplanamaz. "
                            "σ_Δ(1.00) = 0 (nokta değer) varsayılır; bu "
                            "AÇIKÇA yazılan bir sınırlamadır."),
        "s_kilit_M_log": {"merkez": 0.185, "se": 0.113},
        "s_kilit_QE_log": {"merkez": 0.206, "se": 0.070},
    },
    "iki_okuma": {
        "kilit_de_derinlik": "Δ(D) = Δ(1.00)·(1 − f_D)",
        "kilit_tabani": ("Δ(D) = Δ(1.00)·[(1 − s_kilit)·(1 − f_D) + "
                         "s_kilit]"),
        "f_D": "ölçülen f_1.30 (191c/191d), se = K2 jackknife se",
        "ayrisim_ozdeslik": ("Ayrım(D) := oku_ii − oku_i = Δ(1.00)·"
                             "s_kilit·f_D  (cebirsel özdeşlik)"),
        "sigma_birlesik": ("|Ayrım|·√((se_s_kilit/s_kilit)² + "
                           "(se_f_D/f_D)²)  (Δ(1.00) sabit/se'siz "
                           "varsayılır, yukarıdaki sınırlama)"),
    },
    "karar": ("Ayrım/σ_birleşik ≥ 3 (HEM M HEM Q_E'de) İSE VE 180'in para "
             "birimi hattı (167_olcum→172b) yeni bir gaz için 60 "
             "dakikadan kısa sürede koşulabiliyorsa, Hderin110/120/130 "
             "için ham para birimleri (M, Q_E) ölçülür ve hükmedilir; "
             "aksi hâlde ÖLÇME YAPILMAZ — 'güç yetersiz' (ayrım < 3σ) ya "
             "da 'erişilemedi' (hat 60 dk'da koşamıyor/altyapı yok) diye "
             "kaydedilir. Kurtarma yok, eşik gevşetme yok."),
}

KURALLAR = ("TEK DALGA; 2 dakikayı aşan koşular nohup+arka plan (≤30 sn "
           "yoklama), 5 dk'lık tek ön-plan komutu YASAK; sayı uydurmak / "
           "eşik gevşetmek / ölü kurtarmak YASAK; ölçülemeyen 'erişilemedi'; "
           "3-4 anlamlı hane + jk se; git'e DOKUNULMAZ; sonuç ORTAK "
           "TEFTİŞE (sabah).")

DER = {"1.00": "Hkeskin", "1.10": "Hderin110", "1.20": "Hderin120",
      "1.30": "Hderin130"}

if __name__ == "__main__":
    t0 = time.time()
    dmg = subprocess.run(["date"], capture_output=True, text=True).stdout.strip()
    print("=" * 78)
    print("191a — K0: ÖN-KAYIT (DERİN İKİZ τ≤1.30)")
    print("=" * 78, flush=True)
    print(f"  zaman (`date`): {dmg}", flush=True)

    # git commit doğrulaması (yalnız OKUMA; git'e dokunulmuyor)
    try:
        commit_log = subprocess.run(
            ["git", "-C", str(QM.parent), "log", "-1", "--format=%H %ad",
             "--date=iso", "--", "qm_riemann/KALEM_DERIN_IKIZ_130_23EYL2026.md"],
            capture_output=True, text=True).stdout.strip()
    except Exception as e:
        commit_log = f"OKUNAMADI: {e}"
    print(f"  KALEM commit (git log, salt okunur): {commit_log}", flush=True)

    yh = yeniden_hesap()
    print("\n  kaynaktan yeniden hesap (188b harita_K_gercek + 189 HUKUM):")
    for D in ("1.00", "1.10", "1.20", "1.30"):
        kalem_z = ONGORU["z_D_gercek_kinematik"][D]
        hes_z = yh["z_D"][D]
        print(f"    z_{D}: kalem={kalem_z:.4f}  yeniden-hesap={hes_z:.4f}  "
              f"fark={abs(kalem_z-hes_z):.1e}")
    print(f"    z_tam(187c): kalem={ONGORU['z_D_gercek_kinematik']['tam_187c']:.4f}"
          f"  yeniden-hesap={yh['z_inf']:.4f}")
    print(f"    yasa eğimi: 1.10={yh['yasa_egim']['1.10']:.4f} (kalem 3.38)  "
          f"1.20={yh['yasa_egim']['1.20']:.4f} (kalem 3.39)  "
          f"merkez={yh['yasa_egim']['merkez']:.4f} (kalem 3.385)")
    print(f"    f_1.30 öngörü (kaynaktan): {yh['f130_pred_kaynak']:.4f}  "
          f"(kalem 0.612)")
    print(f"    σ_ε(1.30) öngörü bandı (kaynaktan, dar): "
          f"[{yh['sigma130_pred_kaynak_bant'][0]:.4f}, "
          f"{yh['sigma130_pred_kaynak_bant'][1]:.4f}]  "
          f"(kalem [0.4140,0.4160], geniş/güvenli)")
    esit_z = all(abs(ONGORU["z_D_gercek_kinematik"][D] - yh["z_D"][D]) < 5e-5
                for D in ("1.00", "1.10", "1.20", "1.30"))
    esit_slope = (abs(yh["yasa_egim"]["1.10"] - 3.38) < 0.01
                 and abs(yh["yasa_egim"]["1.20"] - 3.39) < 0.01)
    print(f"\n  z_D dörtlüsü kalemle EŞİT (±5e-5): {esit_z}")
    print(f"  yasa eğimleri kalemle EŞİT (±0.01): {esit_slope}")

    ONKAYIT = {
        "gorev": "191 — DERİN İKİZ τ≤1.30: zarf tam mı kapanır, gerçek "
                "kalıntı mı kalır?",
        "kalem": "KALEM_DERIN_IKIZ_130_23EYL2026.md AYNEN",
        "kalem_sha256": sha(QM / "KALEM_DERIN_IKIZ_130_23EYL2026.md"),
        "kalem_git_commit": commit_log,
        "insa_K1_durumu": ("YAPILDI (191b_insa.py; kalem commitinden SONRA, "
                           "bu ön-kayıt JSON'undan önce) — 425863 çizgi, "
                           "kapılar GEÇTİ: maks|F|=1.863e-09, sıralılık "
                           "TAM, L=12.02959324 "
                           "(scratchpad/189/insa_kapi_Hderin130.json)"),
        "gaz_adlari": DER,
        "onceki_taban_189": {
            "aciklama": "189'un ölçtüğü 1.00/1.10/1.20 sayıları (191'in "
                        "karşılaştırma tabanı; 191 bunları YENİDEN "
                        "ÖLÇMEZ, 189/HUKUM_189.json'dan okunur)",
            "f110_havuz": yh["f110_havuz"], "f120_havuz": yh["f120_havuz"],
            "sigma_eps": {"1.00": 0.43135, "1.10": 0.42248, "1.20": 0.41788,
                          "gercek": 0.40921},
            "K3_fit_c_alfa": {
                "1.00": [0.14921, 1.30407], "1.10": [0.11553, 1.68064],
                "1.20": [0.09971, 2.22458]},
        },
        "ongoru": ONGORU,
        "yeniden_hesap_kaynaktan": yh,
        "yeniden_hesap_tutarli": {"z_D": bool(esit_z),
                                  "yasa_egim": bool(esit_slope)},
        "hipotezler": HIPOTEZLER,
        "K4_guc_analizi": K4_KURAL,
        "zincir": {
            "makine": ("189c_zincir.py + 189d_hukum.py AYNEN (importlib "
                       "ile yüklenir, DOSYALAR DEĞİŞTİRİLMEZ); 191c ince "
                       "sarmalayıcı yalnız 'Hderin130' gaz adını ekler ve "
                       "189c.zincir('Hderin130') çağırır. Hkeskin YENİDEN "
                       "KOŞULMAZ (189'un muhur_Hkeskin.json'u — 6/6 TUTTU "
                       "— aynı, DEĞİŞTİRİLMEMİŞ 189c koduyla üretildiği "
                       "için bit-bit geçerli kalır; sha256 denetimiyle "
                       "doğrulanır)."),
            "makine_muhru_kanit": ("sha256(189_configs/189c_zincir.py) ve "
                                   "sha256(189_configs/189d_hukum.py) K0'da "
                                   "kaydedilir; K2/K3 çalışırken yeniden "
                                   "hesaplanıp EŞİT olduğu doğrulanır — "
                                   "değişmemiş makine kanıtı."),
            "sha_189c": sha(QM / "189_configs" / "189c_zincir.py"),
            "sha_189d": sha(QM / "189_configs" / "189d_hukum.py"),
            "uzerine_yazma_yasak": ("scratchpad/189/*_Hkeskin.*, "
                                    "*_Hderin110.*, *_Hderin120.* dosyaları "
                                    "ÜZERİNE YAZILMAZ; yalnız YENİ "
                                    "*_Hderin130.* dosyaları eklenir."),
        },
        "jackknife": "8-blok loo (184-189 AYNEN), se = √(7/8·Σ(θ_i−θ̄)²)",
        "kurallar": KURALLAR,
        "zaman": time.strftime("%a %b %d %H:%M:%S %z %Y"),
        "sha256": None,
    }
    sha_self = hashlib.sha256(
        (QM / "191_configs" / "191a_onkayit.py").read_bytes()).hexdigest()
    ONKAYIT["sha256"] = sha_self
    yol = S191 / "ONKAYIT_191.json"
    json.dump(ONKAYIT, open(yol, "w"), indent=1, ensure_ascii=False,
              default=float)
    print(f"\nONKAYIT_191.json yazıldı -> {yol}")
    print(f"  sha256 = {sha_self}")
    print(f"  damga  = {ONKAYIT['zaman']}")
    print(f"  sha256(189c_zincir.py) = {ONKAYIT['zincir']['sha_189c']}")
    print(f"  sha256(189d_hukum.py)  = {ONKAYIT['zincir']['sha_189d']}")
    print(f"  süre {time.time()-t0:.1f}s")
