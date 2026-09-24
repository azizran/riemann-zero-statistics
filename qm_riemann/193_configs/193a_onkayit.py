# -*- coding: utf-8 -*-
"""
193a — K0: KURAL-ÖNCE ÖN-KAYIT (KALINTI-SINIFI YASASI — kör, kirliliksiz sınav)
================================================================================
KALEM_KALEM_SINIF_YASASI_24EYL2026 AYNEN. ÖLÇÜMDEN ÖNCE donan: öngörü tablosu
(pay yasası s_r = cos(2π r b/a)/μ(a)), δ = 0.03, 8-blok tanımı (L_b kinematik,
son penceresi mid'inden — sonuç DEĞİL), sınıf tanımları (a, b, sınıflar,
"bölünen" payadan HARİÇ), H-193a/b eşikleri ve ölüm koşulları, KAYIT
tanımları (+log5, +log3, −log3), ve KOMŞU RAPORU: her hedef için a,b ≤ 20
sade rasyonellerin log(a/b) konumlarından |fark| < 0.10 olanlar (yalnız
aritmetik; profillere/verilere BAKILMADAN).

Çıktı: scratchpad/193/ONKAYIT_193.json (sha256 = bu betiğin sha'sı + damga).
"""
import hashlib
import json
import subprocess
import time
from fractions import Fraction
from math import cos, gcd, log, pi
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S155, S188, S193 = SCR / "155", SCR / "188", SCR / "193"
S193.mkdir(exist_ok=True)
TWO_PI = 2 * np.pi
DELTA = 0.03          # uydu seçim penceresi (Δω, blok-yerel)
KOMSU_YARI = 0.10      # komşu raporu arama yarı-genişliği
KIRLENME_ESIK = 2 * DELTA   # 0.06 — kirlenme riski eşiği


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def mu(n):
    n = int(n)
    if n == 1:
        return 1
    k, p, m = 0, 2, n
    while p * p <= m:
        if m % p == 0:
            m //= p
            if m % p == 0:
                return 0
            k += 1
        p += 1
    if m > 1:
        k += 1
    return (-1) ** k


def s_pred(a, b, r):
    return cos(2 * pi * r * b / a) / mu(a)


# ---------------- kinematik: L_son, 8 eşit blok, L_b (sonuç DEĞİL) ----------------
ONK188 = json.load(open(S188 / "ONKAYIT_188.json"))
L_S = float(ONK188["L"])
MID_S = np.asarray(np.load(S155 / "eta_son_t0.4_c4000.npz")["mid"], float)
N_S = len(MID_S)
KJ_S = np.linspace(0, N_S, 9).astype(int)
LW_S = np.log(MID_S / TWO_PI)
L_B = [float(LW_S[KJ_S[b]:KJ_S[b + 1]].mean()) for b in range(8)]
L_B_YAYILIM = [float(LW_S[KJ_S[b + 1] - 1] - LW_S[KJ_S[b]]) for b in range(8)]

# ---------------- sınıf tanımları (KALEM tablosu AYNEN; pay yasasıyla türetildi) --
HEDEFLER = {
    "+log10": {"n": 10, "a": 10, "b": 1, "siniflar": [1, 3, 7, 9],
               "statu": "KÖR, birincil", "nicel_tol": 0.25},
    "+log7": {"n": 7, "a": 7, "b": 1, "siniflar": [1, 2, 3, 4, 5, 6],
              "statu": "KÖR, ikincil", "nicel_tol": 0.30},
    "+log5": {"n": 5, "a": 5, "b": 1, "siniflar": [1, 2, 3, 4],
              "statu": "doğrulayıcı KAYIT", "nicel_tol": None},
    "+log3": {"n": 3, "a": 3, "b": 1, "siniflar": [1, 2],
              "statu": "kontrol (½/½)", "nicel_tol": 0.10},
    "-log3": {"n": "1/3", "a": 3, "b": None, "a_teori": 1,
              "siniflar": [1, 2],
              "statu": "kontrol, sıfır-yapı (½/½, keyfi mod-3)", "nicel_tol": 0.10},
}
for ad, h in HEDEFLER.items():
    n = Fraction(h["n"]) if isinstance(h["n"], str) else Fraction(h["n"])
    log_n = float(log(float(n)))
    h["log_n"] = log_n
    if ad == "-log3":
        # a_teori = 1 (tek koherent sınıf); mod-3 bölme KEYFİ (sıfır-yapı kontrolü,
        # kapalı-form tahmin YOK) — beklenti yalnız eşit pay (Dirichlet).
        h["s_pred"] = {str(r): 0.5 for r in h["siniflar"]}
    elif ad == "+log3":
        h["s_pred"] = {str(r): round(s_pred(h["a"], h["b"], r), 6)
                       for r in h["siniflar"]}
    else:
        h["s_pred"] = {str(r): round(s_pred(h["a"], h["b"], r), 6)
                       for r in h["siniflar"]}
    # gereken τ' üst sınırı (8 blok; pencere-ötesi evreni τ' ≤ 1.30 içinde mi?)
    h["gereken_tau_ust_blok"] = [float((lb + log_n + DELTA) / L_S) for lb in L_B]
    h["evren_ici_mi"] = bool(max(h["gereken_tau_ust_blok"]) <= 1.30)

# KALEM tablosuyla çapraz-doğrulama (üretilen formülün metinle uyumu)
KALEM_TABLO = {
    "+log10": {"1": 0.809, "9": 0.809, "3": -0.309, "7": -0.309},
    "+log7": {"1": -0.623, "6": -0.623, "2": 0.223, "5": 0.223,
              "3": 0.901, "4": 0.901},
    "+log5": {"1": -0.309, "4": -0.309, "2": 0.809, "3": 0.809},
    "+log3": {"1": 0.5, "2": 0.5},
}
for ad, tbl in KALEM_TABLO.items():
    for r, v in tbl.items():
        fark = abs(HEDEFLER[ad]["s_pred"][r] - v)
        assert fark < 5e-4, (ad, r, HEDEFLER[ad]["s_pred"][r], v, fark)

# ---------------- H-193a / H-193b eşikleri (KALEM AYNEN) ----------------
H_193A = {
    "hipotez": "+log10 (birincil): r∈{3,7} negatif, r∈{1,9} pozitif, her biri ≥2σ",
    "nicel": "her s_r öngörünün ±0.25 içinde (±0.20 + vekil payı 0.05)",
    "kosul_yon": {"negatif_2s": [3, 7], "pozitif_2s": [1, 9]},
    "olum": ("r∈{3,7}'den biri s_r ≥ +2·se (pozitife anlamlı) YA DA "
             "r∈{1,9}'dan biri s_r ≤ −2·se (negatife anlamlı)"),
    "kurtarma": "YOK"}
H_193B = {
    "hipotez": "+log7 (ikincil): r1,r6 negatif; r3,r4 pozitif, ≥2σ",
    "nicel": "her s_r öngörünün ±0.30 içinde",
    "kosul_yon": {"negatif_2s": [1, 6], "pozitif_2s": [3, 4]},
    "olum": "r1 ya da r6 s_r ≥ +2·se (pozitife anlamlı)",
    "erisilemedi_kosulu": "κ_top < 4·se_kappa_top (güç yetersiz)",
    "kurtarma": "YOK"}
KAYIT_TANIM = {
    "+log5": "doğrulayıcı; s_r kaydedilir, eşiksiz (yalnız KAYIT)",
    "+log3": "kontrol; her iki sınıf s_r = 0.5 ± 0.10 bandında mı (eşiksiz KAYIT)",
    "-log3": ("kontrol (sıfır-yapı); keyfi mod-3 bölme s_r = 0.5 ± 0.10 bandında "
              "mı (eşiksiz KAYIT) — yapı YOKLUĞUNUN sınavı")}

# ---------------- ölçüm/blok tanımı (KALEM AYNEN) ----------------
BLOK_TANIM = {
    "tanim": ("8 eşit bitişik blok, kenar = linspace(0, N, 9).astype(int) "
              "(184-192 jk blokları AYNEN, N = 299 999 'son'); L_b = blok içi "
              "mean log(m_n/2π) (188f(c)/190a AYNEN)."),
    "kenar": KJ_S.tolist(), "L_b": L_B, "blok_ici_L_yayilimi": L_B_YAYILIM}
OLCUM_TANIM = {
    "pencere": "SON (L_son = 12.029593241726252; 188 makinesi AYNEN)",
    "delta": DELTA,
    "secim": ("Her blok b için, pencere-ötesi çizgi q' (asal kuvvet, "
              "187b.asal_kuvvetler(e^{1.30·L_son}, L_son) AYNEN, τ'≤1.30 evreni), "
              "Δω_b(q') = log q' − L_b; uydu n seçimi |Δω_b(q') − log n| < "
              f"δ = {DELTA}. Seçilen çizgiler q' mod a sınıflarına ayrılır "
              "(gcd(q',a)>1 olanlar ayrı 'bölünen' sınıfı — PAYADAN HARİÇ, "
              "yalnız raporlanır)."),
    "seri": ("188b.seri_ve_G AYNEN, blok b'nin YALNIZ kendi noktalarında "
             "(g_b, mid_b, bid_b≡b) — 190b'nin 'tek-blok izdüşüm deseni' "
             "(izdusum_blok); alt-örneklem n=0,3,6,… (188/192 ile aynı, 100 000); "
             "bir hedef-sınıfta (+log10 r=1) TAM örneklem kontrolü."),
    "kappa": ("κ_r = −Re K_r; K_r = 188 K(b,s) tanımı (K_matris; karışım = "
              "c^kesik − c^öz, 186/187 defterleri; HAVUZ bandı, 188b "
              "bant_maskeleri). Bloklar toplanır (188b.c_proj: Σ_b re/Σ_b im, "
              "blok-loo destekli). 8-blok jackknife (184-192 AYNEN): "
              "se = √(7/8·Σ(θ_i−θ̄)²)."),
    "pay_yasasi": "s_r := κ_r/κ_top = cos(2π r b/a) / μ(a); κ_top = Σ_{r∈sınıflar} κ_r "
                  "(bölünen HARİÇ)."}

# ---------------- KOMŞU RAPORU (yalnız aritmetik; profillere/verilere bakılmadan) --
def komsu_raporu(hedef_log, hedef_ab):
    out = []
    for p in range(1, 21):
        for q in range(1, 21):
            if gcd(p, q) != 1:
                continue
            if hedef_ab is not None and (p, q) == hedef_ab:
                continue
            pos = log(p / q)
            fark = pos - hedef_log
            if abs(fark) < KOMSU_YARI:
                out.append({
                    "delta_omega": round(pos, 6), "a": p, "b": q, "mu_a": mu(p),
                    "fark_hedefe": round(fark, 6),
                    "kirlenme_riskli": bool(abs(fark) < KIRLENME_ESIK - 1e-12)})
    return sorted(out, key=lambda x: abs(x["fark_hedefe"]))


HEDEF_AB = {"+log10": (10, 1), "+log7": (7, 1), "+log5": (5, 1),
            "+log3": (3, 1), "-log3": (1, 3)}
KOMSU = {}
for ad, h in HEDEFLER.items():
    KOMSU[ad] = komsu_raporu(h["log_n"], HEDEF_AB[ad])

# ---------------- ÖN-KAYIT ----------------
ONKAYIT = {
    "gorev": "193 — KALINTI-SINIFI YASASI: kör, kirliliksiz sınav",
    "kalem": "KALEM_SINIF_YASASI_24EYL2026.md AYNEN",
    "girdi_sha256": {
        "KALEM_SINIF_YASASI_24EYL2026.md": sha(QM / "KALEM_SINIF_YASASI_24EYL2026.md"),
        "192c_sinif.py": sha(QM / "192_configs" / "192c_sinif.py"),
        "190b_harita.py": sha(QM / "190_configs" / "190b_harita.py"),
        "190e_kesif.py": sha(QM / "190_configs" / "190e_kesif.py"),
        "188b_harita.py": sha(QM / "188_configs" / "188b_harita.py"),
        "187b_katman_defteri.py": sha(QM / "187_configs" / "187b_katman_defteri.py"),
        "185b_oz_muhasebe.py": sha(QM / "185_configs" / "185b_oz_muhasebe.py"),
        "155/eta_son_t0.4_c4000.npz": sha(S155 / "eta_son_t0.4_c4000.npz"),
        "188/ONKAYIT_188.json": sha(S188 / "ONKAYIT_188.json"),
        "188/dilim_44.npz": sha(S188 / "dilim_44.npz"),
        "188/harita_K_gercek.npz": sha(S188 / "harita_K_gercek.npz")},
    "kinematik": {"L_son": L_S, "N_son": N_S,
                  "not": "L_b yalnız mid'den (190a/192a AYNEN); sonuç değil."},
    "bloklar": BLOK_TANIM,
    "olcum": OLCUM_TANIM,
    "hedefler": HEDEFLER,
    "H_193a": H_193A,
    "H_193b": H_193B,
    "KAYIT": KAYIT_TANIM,
    "komsu_raporu": KOMSU,
    "jackknife": "8-blok loo (184-192 AYNEN), se = √(7/8·Σ(θ_i−θ̄)²)",
    "kurallar": ("TEK DALGA; ölümler kurtarmasız; eşik gevşetme yok; ölçülemeyen "
                 "'erişilemedi'; ön-kayıtlı kapılar (makine mührü dahil) "
                 "ATLANMAZ; git'e dokunulmaz; sonuç ORTAK TEFTİŞE."),
    "zaman": None,
    "sha256": None,
}

if __name__ == "__main__":
    try:
        ONKAYIT["zaman"] = subprocess.check_output(["date"], text=True).strip()
    except Exception:
        ONKAYIT["zaman"] = time.strftime("%a %b %d %H:%M:%S %Z %Y")
    s = hashlib.sha256(Path(__file__).resolve().read_bytes()).hexdigest()
    ONKAYIT["sha256"] = s
    yol = S193 / "ONKAYIT_193.json"
    if yol.exists():
        raise SystemExit(f"ONKAYIT_193 ZATEN VAR — üzerine yazılmaz: {yol}")
    json.dump(ONKAYIT, open(yol, "w"), indent=1, ensure_ascii=False)
    print("=" * 78)
    print("193a / K0 ÖN-KAYIT yazıldı")
    print("=" * 78)
    print(f"  sha256 = {s}")
    print(f"  damga  = {ONKAYIT['zaman']}")
    print(f"  L_son = {L_S:.9f}  N = {N_S}")
    print(f"  L_b = {np.round(L_B, 5).tolist()}")
    for ad, h in HEDEFLER.items():
        print(f"  {ad:>7} (n={h['n']}, a={h['a']}, statü={h['statu']}): "
              f"s_pred={h['s_pred']}  evren_ici={h['evren_ici_mi']}")
        riskli = [k for k in KOMSU[ad] if k["kirlenme_riskli"]]
        print(f"      komşu ({len(KOMSU[ad])}, riskli {len(riskli)}): "
              f"{[(k['a'], k['b'], k['fark_hedefe']) for k in riskli]}")
    print(f"  -> {yol}")
