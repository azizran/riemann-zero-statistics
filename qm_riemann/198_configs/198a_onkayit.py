# -*- coding: utf-8 -*-
"""
198a — K0: KURAL-ÖNCE ÖN-KAYIT (κ_p'NİN YÜKSEKLİK BAĞIMLILIĞI) — TASLAK, ÇALIŞTIRMA
=================================================================================
KALEM_KAPPA_YUKSEKLIK_25EYL2026.md son hâli. Ölçümden ÖNCE donan: KALEM sha256 + metni,
girdi sha256'ları (üç pencerenin zinciri, 197/198 kodları, 197 HUKUM/harita, kapı dosyaları),
ölçüm parametreleri (geometri eşit-ΔL, havuz maskesi M7'den, pencere başına kapsama modu,
τ'_üst, blok sayıları, L_W, N — kinematik, sonuç DEĞİL), f_γ/f_Λ ve tohum kümeleri, kapı
sonuçları (M1, M6, M7, M8), eşikler, damga. W_alt/W_son'da hiçbir κ/K/P OKUNMAZ.
Yazmayı REDDEDER: M1/M6(H-198a)/M7/M8 geçmemişse, kapı dosyaları güncel 198c sürümüyle
üretilmemişse, ya da W_alt kapsama modu (KALEM "TAM kapsama"; W_alt eşit-ΔL'de B_j min 59/72)
kaptanca belirlenmemişse.
Çıktı: scratchpad/198/ONKAYIT_198.json (sha256 = bu betiğin sha'sı + damga).
"""
import hashlib
import json
import subprocess
import time
from pathlib import Path

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = QM / "scratchpad"
S190, S197, S198 = SCR / "190", SCR / "197", SCR / "198"
KALEM = QM / "KALEM_KAPPA_YUKSEKLIK_25EYL2026.md"
GEOMETRI = "esitL"
# KAPTAN KARARI (KALEM "MAKİNE RAPORU"): W_alt kapsaması KISMİ — 192 kuralı (dilim j'ye yalnız onu
# TAM kapsayan bloklar; çizgi biçimi aynı B_j ile); 13/72 düşük-L blok alt uçta kapsamıyor, B_j ≥ 59/72.
KAPSAMA = {"alt": "kismi", "dusuk": "tam", "son": "tam"}


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def rel(p):
    return str(Path(p).relative_to(QM))


def kalem_bolum(onek):
    s = KALEM.read_text()
    i = s.index(f"\n## {onek}") + 1
    j = s.find("\n## ", i + 3)
    return s[i:j if j > 0 else None].strip()


ZINCIRLER = {"alt": (S198 / "zincir_alt", "gercek_alt", "eta_alt_t0.4_c4000.npz"),
             "dusuk": (S190 / "zincir_dusuk", "gercek_dusuk", "eta_dusuk_t0.4_c4000.npz"),
             "son": (S190 / "muhur_son", "gercek_sonM", "eta_son_t0.4_c4000.npz")}
GIRDI = []
for W, (Z, et, eta) in ZINCIRLER.items():
    GIRDI += [Z / eta, Z / f"K1_{et}.npz", Z / f"OZ_{et}.npz", Z / f"G1_proj_{et}.npz"]
GIRDI += [QM / "198_configs" / f for f in ("198k0_zincir.py", "198b_harita.py", "198c_hukum.py")]
GIRDI += [QM / "197_configs" / f for f in ("197b_harita.py", "197c_hukum.py")]
GIRDI += [QM / "190_configs" / "190b_harita.py", QM / "188_configs" / "188b_harita.py",
          QM / "187_configs" / "187b_katman_defteri.py", QM / "185_configs" / "185b_oz_muhasebe.py",
          S197 / "HUKUM_197.json", S197 / "harita_omega32.npz", SCR / "192" / "ONKAYIT_192.json",
          S198 / "F_198.json", S198 / "M1_198.json", S198 / "M6_on_198.json",
          S198 / "M7_198.json", S198 / "M8_198.json",
          S198 / "dusuk_197" / "harita_omega.npz", S198 / "dusuk_esitL" / "harita_omega.npz"]
GIRDI += [S198 / f"{W}_{GEOMETRI}" / "plan_198.json" for W in ZINCIRLER]


def kapi(ad):
    return json.load(open(S198 / ad))


ONKAYIT = {
    "gorev": "198 — κ_p'nin yükseklik bağımlılığı (W_alt kör, W_düşük referans, W_son yarı-kör)",
    "kalem": KALEM.name + " AYNEN", "kalem_sha256": None, "girdi_sha256": None,
    "kalem_baslik": None, "kalem_metni": None,
    "olcum_parametreleri": None, "pencereler": None, "f": None, "kapilar": None,
    "esikler": {
        "M1": "W_düşük (197 geometrisi) s ve İKİNCİL genlikleri ≤ 1e-10 bağıl",
        "M2": "her pencere: ana-6 P > 3se; 2A(−log2)/3A(−log3) ∈ [0.85, 1.15]",
        "M3": "tarama ∩ HAVUZ' = ∅ (her pencere, assert)",
        "M4": "σ_eff(γ) = f_γ·σ_jk ≤ 0.35; değilse H-198a belirsiz",
        "M5p": "her pencere: 197 M5 (rms ≤ 0.35·std, tepe korr ≥ 0.85; kalibrasyonda −log2, + "
               "yanda +log2) VE R_2 = 2A_2/s, R_3 = 3A_3/s ∈ [0.9, 1.1]; kalan pencere karardan çıkar",
        "M6": "H-198a: doğru ≥ 0.80, yanlış ≤ 0.05 (iki doğru); H-198b gücü < 0.80 ⇒ yalnız KAYIT",
        "M7": "W_düşük havuz yarıları χ²_4 p > 0.01; kalırsa ortak mutlak havuz q ∈ [224, 992]",
        "M8": "H_C-doğru enjeksiyon, pencere geometrisi (W_alt kısmi kapsama); KARAR NİCELİĞİ "
              "|bias_γ| ≤ 0.25σ_eff(γ) (makine raporu sonrası); tek tek κ_p ve Λ yanlılıkları KAYIT",
        "H_198a": "H_S TUTAR ⇔ |γ̂−1| ≤ 2σ_eff VE |γ̂| > 2σ_eff; H_C TUTAR ⇔ |γ̂| ≤ 2σ_eff VE "
                  "|γ̂−1| > 2σ_eff; İKİSİ DE ÖLÜR ⇔ |γ̂| > 3σ_eff VE |γ̂−1| > 3σ_eff; diğer belirsiz",
        "H_198b": "aynı kural Λ_alt için (H_S L_d/L_alt; H_C 1)"},
    "kaptan_kararlari": {"W_alt_kapsama": "KISMİ (192 kuralı)", "M8": "karar niceliği γ̂ yanlılığı",
                         "H_198b": "KAYIT (M6 gücü < 0.80)",
                         "kabul_edilen_yorumlar": ["eşit-ΔL bloklar zincirin 8 jk grubu içinde",
                                                   "ayna kapısı R_2 ve R_3 ayrı ayrı [0.9, 1.1]",
                                                   "+ yan M5 +log2 tepesiyle", "M7 köşegen χ²",
                                                   "M8'de p ≥ 11 için κ = p^{−1/3}"]},
    "kurallar": "TEK DALGA; kurtarma yok; git'e ajan dokunmaz; K0 push edilmeden ölçüm başlamaz; "
                "W_alt/W_son'da ölçüm-öncesi κ/K/P görülmez.",
    "zaman": None, "sha256": None,
}

if __name__ == "__main__":
    if KAPSAMA["alt"] not in ("tam", "kismi"):
        raise SystemExit("W_alt kapsama modu kaptanca belirlenmedi — ön-kayıt yazılmaz")
    kod = sha(QM / "198_configs" / "198c_hukum.py")
    M1, M6, M7, M8 = (kapi(a) for a in ("M1_198.json", "M6_on_198.json", "M7_198.json",
                                         "M8_198.json"))
    for ad, d in (("M1", M1), ("M6", M6), ("M7", M7), ("M8", M8)):
        if d.get("kod_sha256") != kod:
            raise SystemExit(f"{ad} güncel 198c sürümüyle üretilmemiş — ön-kayıt yazılmaz")
    if not (M1["gecti"] and M6["gecti"] and M8["gecti"]):
        raise SystemExit("M1/M6/M8 geçmedi — KALEM: sınav yapılmaz; ön-kayıt yazılmaz")
    F_ = json.load(open(S198 / "F_198.json"))
    if F_["maske"] != M7["maske"]:
        raise SystemExit("f kalibrasyon havuzu ≠ M7 kararı")
    plan = {W: json.load(open(S198 / f"{W}_{GEOMETRI}" / "plan_198.json")) for W in ZINCIRLER}
    ONKAYIT["kalem_sha256"] = sha(KALEM)
    ONKAYIT["kalem_baslik"] = KALEM.read_text().split("\n## ")[0].strip()
    ONKAYIT["kalem_metni"] = {b: kalem_bolum(b) for b in
                              ("Durum", "Kalem cebiri", "Pencereler", "Ölçüm tanımı", "Kapılar",
                               "Hipotezler", "KAYIT", "Ölçüm notları")}
    ONKAYIT["girdi_sha256"] = {rel(p): sha(p) for p in GIRDI}
    ONKAYIT["olcum_parametreleri"] = {
        "geometri": GEOMETRI, "maske": M7["maske"], "kapsama": KAPSAMA,
        "havuz": [0.45, 0.74], "tau_alt": 0.74, "dw": 0.025, "dw_aralik": [-1.40, 2.40],
        "n_blok_kurali": "max(32, ceil(ΔL_W/0.03)); zincir 8 grubu içinde eşit-L",
        "kalibrasyon_penceresi": [-1.30, -0.55], "ikincil_penceresi": [0.15, 2.35],
        "L_d": 10.48392954102207}
    ONKAYIT["pencereler"] = {W: {k: plan[W][k] for k in (
        "aralik", "L_W", "N", "L_min", "L_max", "n_blok", "tau_ust", "havuz_ussu", "B_j_min",
        "tam_kapsama", "kapsamayan_bloklar", "kapsanan_gap_payi_min", "M3")} for W in plan}
    ONKAYIT["f"] = {"f_gamma": F_["f_gamma"], "f_Lambda": F_["f_Lambda"],
                    "K_KAL": F_["tohumlar"], "K_GUC": M6["tohumlar"], "kaynak": "F_198.json"}
    ONKAYIT["kapilar"] = {
        "M1": {k: M1[k] for k in ("s_maks_bagil", "IKINCIL_A_maks_bagil", "gecti")},
        "M6": {"gecti": M6["gecti"], "H198b_yalniz_KAYIT": M6["H198b_yalniz_KAYIT"],
               "dogrular": M6["dogrular"]},
        "M7": {k: M7[k] for k in ("chi2", "p", "gecti", "karar")},
        "M8": {"gecti": M8["gecti"], "gamma_bias": M8["gamma_bias"],
               "gamma_bias_mc_se": M8["gamma_bias_mc_se"],
               "sigma_eff_gamma_medyan": M8["sigma_eff_gamma_medyan"],
               "gamma_bias_bolu_sigma": M8["gamma_bias_bolu_sigma"],
               "KAYIT_kappa_bias": {W: M8["pencereler"][W]["KAYIT_bias"] for W in M8["pencereler"]},
               "KAYIT_Lambda_bias": M8["KAYIT_Lambda_bias"]}}
    try:
        ONKAYIT["zaman"] = subprocess.check_output(["date"], text=True).strip()
    except Exception:
        ONKAYIT["zaman"] = time.strftime("%a %b %d %H:%M:%S %Z %Y")
    ONKAYIT["sha256"] = sha(Path(__file__).resolve())
    yol = S198 / "ONKAYIT_198.json"
    if yol.exists():
        raise SystemExit(f"ONKAYIT_198 ZATEN VAR — üzerine yazılmaz: {yol}")
    json.dump(ONKAYIT, open(yol, "w"), indent=1, ensure_ascii=False)
    print(f"198a / K0 ÖN-KAYIT yazıldı: sha256 {ONKAYIT['sha256']}  damga {ONKAYIT['zaman']}\n"
          f"  -> {yol}")
