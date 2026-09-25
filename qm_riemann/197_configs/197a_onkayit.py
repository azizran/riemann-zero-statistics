# -*- coding: utf-8 -*-
"""
197a — K0: KURAL-ÖNCE ÖN-KAYIT (AYNA YASASI — BK ayna-yanı uydularının kör profil sınavı)
=====================================================================================
TASLAK — KALEM_AYNA_YASASI_25EYL2026.md'nin SON hâline göre; ÖLÇÜMDEN ÖNCE (M1 geçtikten
sonra) ÇALIŞTIRILIR. Donan: KALEM'in sha256'sı ve ilgili bölümleri METİN olarak (AYNEN),
girdi dosyalarının sha256'ları (zincir_dusuk, 190 τ'-dilimleri, referans harita, makine
kodu), S1-S3 parametreleri (32 blok kenarı ve L_b — kinematik, sonuç DEĞİL), hedef listeleri
(kör / kalibrasyon / M2 / KAYIT) ve uyum pencereleri, eşikler (yapılandırılmış + KALEM
metni), M1 kapı sonucu (varsa), zaman damgası. κ/K/P OKUNMAZ.

Çıktı: scratchpad/197/ONKAYIT_197.json (sha256 = bu betiğin sha'sı + damga).
197b 'olcum' modu bu dosyayı, sha'ları ve parametreleri denetlemeden başlamaz.
"""
import hashlib
import json
import subprocess
import time
from fractions import Fraction as F
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = QM / "scratchpad"
S190, S192, S197 = SCR / "190", SCR / "192", SCR / "197"
ZD = S190 / "zincir_dusuk"
KALEM = QM / "KALEM_AYNA_YASASI_25EYL2026.md"
TWO_PI = 2 * np.pi
NJACK = 8
DW = 0.025
EPS = 1e-9


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def rel(p):
    return str(Path(p).relative_to(QM))


def kalem_bolum(onek):
    """KALEM'den '## <onek>…' ile başlayan bölümü AYNEN (metin) al."""
    s = KALEM.read_text()
    i = s.index(f"\n## {onek}") + 1
    j = s.find("\n## ", i + 3)
    return s[i:j if j > 0 else None].strip()


# ---------------- kinematik (sonuç DEĞİL: mid'den L, 32 blok L_b) ----------------
E = np.load(ZD / "eta_dusuk_t0.4_c4000.npz")
MID = np.asarray(E["mid"], float)
L = float(E["L"])
N = len(MID)
KB32 = np.linspace(0, N, 33).astype(int)
KJ8 = np.linspace(0, N, 9).astype(int)
assert np.array_equal(KB32[::4], KJ8), "linspace(0,N,33)[4k] ≠ linspace(0,N,9)[k]"
L_B32 = [float(np.log(MID[KB32[b]:KB32[b + 1]] / TWO_PI).mean()) for b in range(32)]
LW = np.log(MID / TWO_PI)
L_B32_YAYILIM = [float(LW[KB32[b + 1] - 1] - LW[KB32[b]]) for b in range(32)]

# ---------------- S1 HAVUZ' ----------------
K1 = np.load(ZD / "K1_gercek_dusuk.npz")
TAU = np.asarray(K1["tau"], float)
HAVUZ = [0.45, 0.74]
N_HAVUZ = int(((TAU >= HAVUZ[0]) & (TAU < HAVUZ[1])).sum())

# ---------------- S2 τ'-dilimleri ----------------
ONK190 = json.load(open(S190 / "ONKAYIT_190.json"))
KENAR_ESKI = ONK190["dilim_izgara"]
KENAR_YENI = [round(0.74 + 0.005 * k, 3) for k in range(25)]
assert KENAR_YENI[-1] == KENAR_ESKI[0] == 0.86
DW_ARALIK = [-2.2, 2.4]                      # KALEM "Ölçüm notları"
JLO, JHI = int(round(DW_ARALIK[0] / DW)), int(round(DW_ARALIK[1] / DW))


# ---------------- hedefler (KALEM) ----------------
def pencere_j(h, yari=0.03):
    return [j for j in range(-200, 200) if abs(j * DW - h) <= yari + EPS]


def ayna(x, rol):
    """Ayna yanı uydu x (Δω = −log x; r = 1/x). c_BK: ½ℕ → 1/x; b/3, b/6 → 0.25/x."""
    x = F(x)
    c = 1 / x if (2 * x).denominator == 1 else F(1, 4) / x
    return {"x": str(x), "x_ondalik": float(x), "delta_omega": -float(np.log(float(x))),
            "log_r": -float(np.log(float(x))), "c_BK": float(c), "rol": rol}


def liste(xs, rol):
    return {f"x={F(x)}": ayna(F(x), rol) for x in xs}


L_G = 10.48392954102207
_g = lambda om: float(np.exp(om / 2.0) / om)
_KALC = [F(11, 6), F(13, 6), F(7, 3), F(8, 3), F(17, 6), F(19, 6), F(10, 3), F(11, 3), F(5, 3),
         F(23, 6)]
G_CAL = sum(0.25 / float(x) * _g(L_G - np.log(float(x))) for x in _KALC) / \
    sum(0.25 / float(x) for x in _KALC)
KOR = {
    "pencere": [-2.12, -1.58],
    "ana": liste([5, F(11, 2), 6, F(13, 2), 7, F(15, 2), 8], "ana (serbest A_k)"),
    "sikinti_dusen": ["9/2", "17/2", "9"],
    "ceyrek": liste([F(31, 6), F(16, 3), F(17, 3), F(35, 6), F(37, 6), F(19, 3), F(20, 3),
                     F(41, 6), F(43, 6), F(22, 3), F(23, 3), F(47, 6), F(49, 6)],
                    "çeyrek (SABİT: q_f^cal·0.25/x [·g/ḡ_cal M_Rg'de])"),
    "ceyrek_kenar": liste([F(29, 6), F(14, 3), F(25, 3)], "çeyrek kenar (SABİT)"),
    "taban": "a + bΔω serbest",
    "modeller": {"M_BK": "genlik_i = q_f^cal·0.25/x_i",
                 "M_Rg": "genlik_i = q_f^cal·0.25/x_i·g(L + Δω_i)/ḡ_cal"},
    "g_cal": G_CAL,
    "g_cal_tanim": "kalibrasyon çeyrek tepelerinde (listelenen 10) 0.25/x ağırlıklı "
                   "g(L − log x) ortalaması",
    "qf_cal_jk": "replika v'de q_f^cal^{(v)} aynı grupla kalibrasyonda yeniden uydurulur ve "
                 "kör uyuma taşınır",
    "Z_ailesi": ["5", "6", "7", "8"], "H_ailesi": ["11/2", "13/2", "15/2"]}
KALIB = {
    "pencere": [-1.30, -0.55],
    "ana": liste([2, F(5, 2), 3, F(7, 2)], "ana (serbest A_k)"),
    "ceyrek": liste([F(11, 6), F(13, 6), F(7, 3), F(8, 3), F(17, 6), F(19, 6), F(10, 3),
                     F(11, 3)], "çeyrek (ortak q_f)"),
    "ceyrek_kenar": liste([F(5, 3), F(23, 6)], "çeyrek kenar (ortak q_f)"),
    "s_hedefleri": ["2", "3"]}
ANA6 = {"+log2": np.log(2), "+log3": np.log(3), "+log6": np.log(6),
        "+log(3/2)": np.log(1.5), "-log2": -np.log(2), "-log3": -np.log(3)}
M2_ORAN_KALEM = {"+log3/+log2": 0.310, "+log6/+log2": 0.602, "-log2/+log2": 0.240,
                 "-log3/+log3": 0.466}
H192 = json.load(open(S192 / "HUKUM_192.json"))
M2_ORAN_192 = {k: H192["A2"][k]["oran"] for k in M2_ORAN_KALEM}
KAT192 = json.load(open(S192 / "ONKAYIT_192.json"))["B"]["katalog"]
IKINCIL = {ad: {"delta_omega": v["delta_omega"], "liste": v["liste"]}
           for ad, v in KAT192.items()}
RAKIP_RG = {"g": "g(ω') = e^{ω'/2}/ω', ω' = L + Δω, L = 10.48393",
            "R_g": "R^g_x = g(L − log x)/ḡ, ḡ = Σ_a w_a g(L − log x_a)/Σ_a w_a, a ∈ {2,3}, "
                   "w_a = s tahminindeki ağırlıklar",
            "esit_agirlik_KALEM": {"5": 0.754, "6": 0.702, "7": 0.662, "8": 0.629}}

# ---------------- girdi dosyaları ----------------
GIRDI = ([ZD / f for f in sorted(p.name for p in ZD.iterdir() if p.is_file())] +
         [S190 / f"dilim_{i}.npz" for i in range(len(KENAR_ESKI) - 1)] +
         [S190 / "ONKAYIT_190.json", SCR / "188" / "ONKAYIT_188.json",
          S190 / "harita_omega_dusuk.npz", S192 / "ONKAYIT_192.json", S192 / "HUKUM_192.json",
          QM / "185_configs" / "185b_oz_muhasebe.py",
          QM / "187_configs" / "187b_katman_defteri.py",
          QM / "188_configs" / "188b_harita.py", QM / "190_configs" / "190b_harita.py",
          QM / "190_configs" / "190a_onkayit.py", QM / "197_configs" / "197b_harita.py",
          QM / "197_configs" / "197c_hukum.py", S197 / "F_KALIBRASYON.json",
          S197 / "M6_on_deneme.json", S197 / "m1" / "M1_sonuc.json"])


def m1_sonuc():
    out = {}
    for ad in ("M1_sonuc.json", "M1_sonuc_alt.json"):
        y = S197 / "m1" / ad
        if y.exists():
            d = json.load(open(y))
            out[ad] = {k: d.get(k) for k in ("mod", "n_karsilastirma", "maks_bagil_eleman",
                                             "M1a_K_ham", "M1b_K_duz",
                                             "M1c_kesik_katki_maks_bagil", "ncz_esit",
                                             "gecti", "kod_sha256", "zaman")}
    y = S197 / "M6_on_deneme.json"
    if y.exists():                                   # ölçüm-öncesi M6 provası (KAYIT)
        d = json.load(open(y))
        out["M6_on_deneme.json"] = {"gecti": d["M6"]["gecti"], "f": d["M6"].get("f"),
                                    "kod_sha256": d["kod_sha256"], "zaman": d["zaman"]}
    return out or "BEKLİYOR"


ONKAYIT = {
    "gorev": "197 — AYNA YASASI: BK ayna-yanı uydularının kör profil sınavı (düşük pencere)",
    "kalem": KALEM.name + " AYNEN",
    "kalem_sha256": None,
    "girdi_sha256": None,
    "pencere": {"ad": "dusuk", "L": L, "N": N, "L_min": float(LW.min()),
                "L_max": float(LW.max())},
    "kalem_baslik": KALEM.read_text().split("\n## ")[0].strip(),
    "kalem_metni": {b: kalem_bolum(b) for b in
                    ("Kalem cebiri", "Kör bant", "Ölçüm tanımı", "Kapılar", "Hipotezler",
                     "KAYIT kalemleri", "Ölçüm notları")},
    "olcum_parametreleri": {"havuz": HAVUZ, "tau_yeni": [0.74, 0.86], "dtau": 0.005,
                            "n_blok": 32, "dw": DW, "om_chunk": 8000,
                            "dw_aralik": DW_ARALIK},
    "S1": {"havuz": HAVUZ, "n_cizgi": N_HAVUZ,
           "tanim": "HAVUZ' = K1 pencere çizgileri τ = w/L ∈ [0.45, 0.74); mix = "
                    "188b.karisim(P, OZ, aq, loo)[HAVUZ'] AYNEN (çizgi başına; 'dusuk' "
                    "zinciri, tam örneklem); tek maske (HAVUZ'); K_b(j) = Σ C conj(mix)/Σ|mix|²."},
    "S1b": {"birincil": "K_duz",
            "tanim": "ω-dilimi (b, j) için mix_j = mix − m_{S_bj}; m_S = S dilimindeki c^kesik "
                     "evreni çizgilerinin (K1: q ≤ e^{0.86L}) kesik seri terimi 2a·sin(ωg/2)·"
                     "cos(ωm)'nin HAVUZ' çizgilerine TAM-örneklem izdüşümü (çizgi başına, 8-grup "
                     "loo, 188b.c_proj deseni); K_düz = K_matris(C_j, mix_j). τ'-haritasında "
                     "aynı kural τ'-dilimi başına. K_ham (çıkarmasız) KAYIT.",
            "denetim": "M1b: 190 düzeninde K_düz ≡ K bit-bit; M1c/M-kesik: Σ_{tüm K1} katkı "
                       "≡ c^kesik (≤ 1e-10 bağıl, tam + 8 loo)"},
    "S2": {"kenar_yeni": KENAR_YENI, "kenar_eski": KENAR_ESKI,
           "tanim": "τ'-dilim s = (e_s, e_{s+1}] (188 AYNEN); yeni 24 dilim (0.74, 0.86] "
                    "188b.dilimleri_kur AYNEN → scratchpad/197/dilim_<i>.npz; 190'ın 88 "
                    "dilimi (0.86, 1.30] yeniden kullanılır. ω-dilim: j = floor((log q' − "
                    "L_b)/0.025 + 0.5); dilim serisi YALNIZ blok noktalarında 188b.seri_ve_G "
                    "(190b D AYNEN, b190._om_is).",
           "dw_aralik": DW_ARALIK, "j_tam": [JLO, JHI],
           "tarama_kurali": "blok b için [(jlo−½)DW, (jhi+½)DW)'ya değen τ'-dilimlerinin "
                            "TÜM çizgileri taranır (doğrusallık denetimi kesin); aralık "
                            "dışındaki uç ω-dilimleri kısmi olabilir, kapsamaya girmez."},
    "S3": {"n_blok": 32, "kenar": KB32.tolist(), "L_b": L_B32,
           "blok_ici_L_yayilimi": L_B32_YAYILIM,
           "jk_gruplar": [list(range(4 * k, 4 * k + 4)) for k in range(NJACK)],
           "jk_grup_kenar": KJ8.tolist(),
           "tanim": "kenar = linspace(0, N, 33).astype(int); L_b = blok içi mean "
                    "log(m_n/2π); jk: 8 grup × 4 ardışık blok (grup kenarı = 190 8-blok "
                    "kenarı); replika i: grup i blokları B_j'den çıkar + mix loo-i (188b)."},
    "profil": {
        "kappa_b": "κ_b(j) = −Re K_b(j) (192b AYNEN)",
        "kappa_sigma": "κ_Σ(j) = Σ_{b∈B_j} κ_b(j)·n_b / Σ_{b∈B_j} n_b; B_j = j dilimini TAM "
                       "kapsayan bloklar: [c−½DW, c+½DW] ⊂ [0.74L − L_b, 1.30L − L_b] ve "
                       "j ∈ [jlo, jhi] (192 kuralı, tarama alt ucu 0.74); 185 tam dilimin hepsinde B_j = 32 (197b assert)",
        "cizgi_bicimi": "h_b(j; r) = (1/n_b)#{n∈b: floor((L_n − L_b + log r)/0.025 + 0.5) = j}; "
                        "h = Σ_{B_j} n_b h_b/Σ n_b (197b.cizgi_bicimi)",
        "uyum": "doğrusal EKK, pencere = dilim merkezi ∈ [lo, hi]; her jk replikası baştan "
                "(197b.profil_uyumu)",
        "M2_P": "P(h) = Σ_{|c−h| ≤ 0.03 (+1e-9)} κ_Σ(c); se: 8-grup loo jk (197b.parlaklik)"},
    "hedefler": {"kor_uyum": KOR, "kalibrasyon_uyum": KALIB,
                 "M2_ana6": {k: {"delta_omega": float(v), "P_pencere_j": pencere_j(v)}
                             for k, v in ANA6.items()},
                 "ikincil_192_23": {"yan_pencere": [0.15, 2.35], "konumlar": IKINCIL,
                                    "tepe_kurali": "BK c ≥ 0.05 + μ=0 konumları serbest (KALEM)"},
                 "rakip_Rg": RAKIP_RG},
    "esikler": {
        "M1": {"bagil_fark_maks": 1e-10,
               "kapsam": "tam harita ya da Δω ∈ [−1.3, 3.5]'e yayılmış ≥ 12 ω-dilimi, 50 K"},
        "M2": {"ana6_P_bolu_se_min": 3.0, "oran_merkez_KALEM": M2_ORAN_KALEM,
               "oran_merkez_192_tam": M2_ORAN_192, "oran_tolerans_goreli": 0.30,
               "basarisizsa": "SINAV GEÇERSİZ (BK ölümü değil)"},
        "M3": "tarama ∩ HAVUZ' = ∅ (197b assert)",
        "M4": {"sigma_eff_maks": 0.10, "tanim": "σ_eff = f_Z·σ_jk ≤ 0.10, M_BK VE M_Rg",
               "degilse": "tüm kör hükümler 'belirsiz'"},
        "M5": {"artik_rms_bolu_pencere_std_maks": 0.35,
               "mlog2_profil_korelasyon_min": 0.85, "mlog2_dilimleri": "|Δω + log2| ≤ 0.06",
               "basarisizsa": "çizgi-biçimi modeli GEÇERSİZ ⇒ hüküm yalnız KAYIT"},
        "M6": {"dogrular": ["BK", "R-g"],
               "i_gurultusuz_goreli_hata_maks": 0.05,
               "ii_kapsama_min": 0.88,
               "ii_nicelikler": "KARARA GİREN: R̄_Z ve ρ |Â − A| ≤ 2σ_eff (f_Z, f_ρ; doğru-eşleşik "
                                "model); kalibrasyon çapaları A(−log2), A(−log3), q_f^cal ve kör "
                                "ana A_x (x = 5..8, 11/2, 13/2, 15/2) |Â − A| ≤ 2.36σ_jk",
               "iii_dogru_kazanma_min": 0.80, "iii_yanlis_kazanma_maks": 0.05,
               "iii_replika": "K_GUC tohumları (f kalibrasyonundan AYRI), 200'er, σ_eff ile",
               "sinyal": "192 ölçeği S = ters-varyans ağ. ort. {2P(−log2), 3P(−log3)} "
                         "(HUKUM_192 A2_P); R-g doğrusu: genlik S·c·g(L+Δω)/ḡ_s, "
                         "ḡ_s = (g(L−log2)+g(L−log3))/2",
               "gurultu": "e_g ~ N(0, 8·se_j²); ölçüm-öncesi se_j M1 (görülmüş, 8 blok, "
                          "HAVUZ' sütunları; kör pencerede medyan), hükümde ölçülen jk se",
               "ol_oncesi_kalirsa": "sınav YAPILMAZ (tasarım yeniden) — 197a ön-kaydı yazmaz",
               "hukum_aninda_kalirsa": "hüküm yalnız KAYIT"},
        "kapi_onceligi": "M1/M3 > M2 > M6 > M5 > M4 (kaptan kabulü)",
        "H_197a": {"sigma": "σ_eff = f_Z·σ_jk (her iki model)",
                   "BK_oz_tutarli": "|R̄_Z^BK − 1| ≤ 2σ_eff^BK (M_BK uyumu)",
                   "Rg_oz_tutarli": "|R̄_Z^Rg − R̄^g| ≤ 2σ_eff^Rg (M_Rg uyumu; R̄^g M_Rg ağırlıklarıyla)",
                   "BK_TUTAR": "BK öz-tutarlı, R-g DEĞİL (β^BK ve se'si KAYIT)",
                   "RG_KAZANIR": "R-g öz-tutarlı, BK DEĞİL",
                   "IKISI_DE_OLUR": "|R̄_Z^BK − 1| > 3σ_eff^BK VE |R̄_Z^Rg − R̄^g| > 3σ_eff^Rg",
                   "diger": "belirsiz (KAYIT)",
                   "H197b_c": "M_BK uyumundan; M_Rg ile aynı hükümler KAYIT"},
        "H_197b": {"rho": "R̄_H/R̄_Z", "TUTAR": "|ρ − 1| ≤ max(0.25, 2σ_ρ), σ_ρ = f_ρ·σ_jk",
                   "OLUM": "ρ < 0.50 ya da üç hedeften ikisinde A < 0 (ÖLÜM önce; en az iki)",
                   "tolerans_min": 0.25, "olum_rho": 0.50, "olum_negatif_A": 2},
        "H_197c": {"psi": "A_8/Â_8, Â_8 = x=5,6,7'ye uydurulan A ∝ x^{−β} zarfı x=8'de",
                   "TUTAR": [0.60, 1.60], "OLUM_alt": 0.30}},
    "makine_muhurleri": None,
    "f": None,
    "M6_ol_oncesi_ozet": None,
    "kurallar": "TEK DALGA; ölümler kurtarmasız; eşik gevşetme yok; ölçülemeyen "
                "'erişilemedi'; git'e ajan dokunmaz; K0 push edilmeden ölçüm başlamaz; "
                "sonuç ORTAK TEFTİŞE.",
    "zaman": None,
    "sha256": None,
}

if __name__ == "__main__":
    ONKAYIT["kalem_sha256"] = sha(KALEM)
    ONKAYIT["girdi_sha256"] = {rel(p): sha(p) for p in GIRDI}
    ONKAYIT["makine_muhurleri"] = m1_sonuc()
    F = json.load(open(S197 / "F_KALIBRASYON.json"))
    ONKAYIT["f"] = {"f_Z": F["f_Z"], "f_rho": F["f_rho"],
                    "tanim": "f = K_KAL replikalarında z = (tahmin − doğru)/σ_jk standart "
                             "sapması (ddof=1), iki doğrunun büyüğü; BK-doğru M_BK, R-g-doğru "
                             "M_Rg; SABİT (yeniden kalibre edilmez)",
                    "K_KAL": F["tohumlar"],
                    "K_GUC": {"BK": [1972000, 1972199], "Rg": [1973000, 1973199]},
                    "kalibrasyon_dogrular": F["dogrular"]}
    M6d = json.load(open(S197 / "M6_on_deneme.json"))["M6"]
    ONKAYIT["M6_ol_oncesi_ozet"] = {
        "gecti": M6d["gecti"], "nrep": M6d["nrep"], "S": M6d["S"],
        "gurultu_se_medyan": M6d["gurultu_se_medyan"],
        "dogrular": {d: {"i_maks_hata": max(r["i_gurultusuz_goreli_hata"].values()),
                         "ii_min_kapsama": min(r["ii_kapsama_karara_giren"].values()),
                         "ii_kapsama": r["ii_kapsama_karara_giren"],
                         "iii_dogru": r["iii_dogru_kazanma"], "iii_yanlis": r["iii_yanlis_kazanma"],
                         "sigma_eff_BK_medyan": r["sigma_eff_BK_medyan"],
                         "sigma_eff_Rg_medyan": r["sigma_eff_Rg_medyan"], "M4_gecme": r["M4_gecme"]}
                     for d, r in M6d["dogrular"].items()}}
    try:
        ONKAYIT["zaman"] = subprocess.check_output(["date"], text=True).strip()
    except Exception:
        ONKAYIT["zaman"] = time.strftime("%a %b %d %H:%M:%S %Z %Y")
    s = hashlib.sha256(Path(__file__).resolve().read_bytes()).hexdigest()
    ONKAYIT["sha256"] = s
    # K0 kapıları: M1 bu 197b sürümüyle geçmiş; ölçüm-öncesi M6 bu 197c sürümüyle geçmiş
    kb_, kc_ = sha(QM / "197_configs" / "197b_harita.py"), sha(QM / "197_configs" / "197c_hukum.py")
    m1y = S197 / "m1" / "M1_sonuc.json"
    if not (m1y.exists() and json.load(open(m1y)).get("gecti") and
            json.load(open(m1y)).get("kod_sha256") == kb_):
        raise SystemExit("M1 bu 197b sürümüyle GEÇİLMEDİ — ön-kayıt yazılmaz")
    m6y = S197 / "M6_on_deneme.json"
    if json.load(open(m6y))["M6"].get("f") != {"Z": ONKAYIT["f"]["f_Z"], "rho": ONKAYIT["f"]["f_rho"]}:
        raise SystemExit("M6'da kullanılan f ≠ F_KALIBRASYON — ön-kayıt yazılmaz")
    if not (m6y.exists() and json.load(open(m6y))["M6"]["gecti"] and
            json.load(open(m6y)).get("kod_sha256") == kc_):
        raise SystemExit("ölçüm-öncesi M6 bu 197c sürümüyle GEÇİLMEDİ — KALEM: sınav "
                         "yapılmaz (tasarım yeniden); ön-kayıt yazılmaz")
    yol = S197 / "ONKAYIT_197.json"
    if yol.exists():
        raise SystemExit(f"ONKAYIT_197 ZATEN VAR — üzerine yazılmaz: {yol}")
    S197.mkdir(exist_ok=True)
    json.dump(ONKAYIT, open(yol, "w"), indent=1, ensure_ascii=False)
    print("=" * 78)
    print("197a / K0 ÖN-KAYIT yazıldı")
    print("=" * 78)
    print(f"  sha256 = {s}")
    print(f"  KALEM sha256 = {ONKAYIT['kalem_sha256']}")
    print(f"  damga  = {ONKAYIT['zaman']}")
    print(f"  girdi dosyası = {len(ONKAYIT['girdi_sha256'])}")
    print(f"  HAVUZ' çizgi = {N_HAVUZ}; 32 blok L_b ∈ [{min(L_B32):.4f}, {max(L_B32):.4f}]")
    print(f"  M1 = {ONKAYIT['makine_muhurleri']}")
    print(f"  -> {yol}")
