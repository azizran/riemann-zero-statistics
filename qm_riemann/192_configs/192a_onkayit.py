# -*- coding: utf-8 -*-
"""
192a — K0: KURAL-ÖNCE ÖN-KAYIT (UYDULARIN TEORİSİ: seçim kuralı + kalıntı dönmesi)
================================================================================
KALEM_UYDU_TEORISI_23EYL2026 AYNEN. PROFİLLERE / HARİTALARA BAKILMADAN donan:
  B  — yanmalı / sönmeli katalog (Δω konumları, n = a/b sade, μ(a)),
       "yanar" / "söner" / "belirsiz" ölçütlerinin SAYISAL tanımı (tepe araması
       ±0.03, yerel taban halkası 0.05 < |Δ| ≤ 0.15, yayılım = 1.4826·MAD),
       havuz profili + 8-blok jk tanımı (188b K_hesap AYNEN: C ve karışım loo),
       dilim kapsaması (kinematik), arama/halka dilim listeleri, çakışma notları.
  A1 — beş uydu penceresi (τ', L_son), sınıflar q' mod a, "bölünen" sınıfı,
       K_r (HAVUZ; 188 tanımı), açı/oran eşikleri, yön-tutarlılık KAYDI,
       makine mührü (188 dilim_44 bit-bit), tam-örneklem kontrolü.
  A2 — parlaklık = uydu çevresi Σκ (|c−h| ≤ 0.0625 → 5 dilim; 188 uydu
       penceresi ±0.006 τ'nün düşük penceredeki ω karşılığı ≈ ±0.063),
       oran bantları (×/÷ 2).
  H-192a/b/c ölüm / mühür kuralları (kurtarma YOK).
Yalnız kinematik bilgi okunur (L, L_b, çizgi evreni = asal sayma); κ, K, profil
değeri OKUNMAZ. Girdi dosyalarının sha256'ları kaydedilir.

Çıktı: scratchpad/192/ONKAYIT_192.json (sha256 = bu betiğin sha'sı + damga).
"""
import hashlib
import importlib.util
import json
import subprocess
import sys
import time
from fractions import Fraction
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S = {d: SCR / d for d in ["155", "184", "185", "186", "188", "190", "192"]}
ZD = S["190"] / "zincir_dusuk"
TWO_PI = 2 * np.pi


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def yukle(ad, yol):
    spec = importlib.util.spec_from_file_location(ad, yol)
    m = importlib.util.module_from_spec(spec)
    sys.modules[ad] = m
    spec.loader.exec_module(m)
    return m


# ---------------- kinematik (sonuç DEĞİL) ----------------
ONK188 = json.load(open(S["188"] / "ONKAYIT_188.json"))
ONK190 = json.load(open(S["190"] / "ONKAYIT_190.json"))
L_D = float(ONK190["pencere"]["L"])
L_B_D = [float(x) for x in ONK190["bloklar"]["L_b"]]
L_S = float(ONK188["L"])
MID_S = np.asarray(np.load(S["155"] / "eta_son_t0.4_c4000.npz")["mid"], float)
KJ_S = np.linspace(0, len(MID_S), 9).astype(int)
L_B_S = [float(np.log(MID_S[KJ_S[b]:KJ_S[b + 1]] / TWO_PI).mean()) for b in range(8)]

DW = 0.025
ARAMA = 0.03          # tepe araması yarı-genişliği (merkez |c − h| ≤ 0.03)
HALKA_IC = 0.05       # yerel taban halkası iç sınırı (hariç: |c − h| ≤ 0.05)
HALKA_DIS = 0.15      # yerel taban halkası dış sınırı
YAY_K = 2.0           # taban bandı = medyan ± 2·yayılım
YAY_MAD = 1.4826      # yayılım = 1.4826·MAD (gürbüz σ)
A2_YARI = 0.0625      # A2 parlaklık penceresi yarı-genişliği (5 dilim)
MIN_KAPSAMA = 4       # bir dilim, en az 4 blok onu TAM kapsıyorsa kullanılır
EPS = 1e-9


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


def phi(n):
    n = int(n)
    r, p, m = n, 2, n
    while p * p <= m:
        if m % p == 0:
            while m % p == 0:
                m //= p
            r -= r // p
        p += 1
    if m > 1:
        r -= r // m
    return r


# ---------------- B KATALOĞU (KALEM AYNEN) ----------------
YANMALI = [("+log2", Fraction(2)), ("+log3", Fraction(3)), ("+log5", Fraction(5)),
           ("+log6", Fraction(6)), ("+log7", Fraction(7)), ("+log10", Fraction(10)),
           ("+log(3/2)", Fraction(3, 2)), ("+log(5/3)", Fraction(5, 3)),
           ("+log(5/2)", Fraction(5, 2)), ("+log(5/4)", Fraction(5, 4)),
           ("-log2", Fraction(1, 2)), ("-log3", Fraction(1, 3)),
           ("-log(3/2)", Fraction(2, 3)), ("-log(4/3)", Fraction(3, 4))]
SONMELI = [("+log4", Fraction(4)), ("+log8", Fraction(8)), ("+log9", Fraction(9)),
           ("+log(4/3)", Fraction(4, 3)), ("+log(8/3)", Fraction(8, 3)),
           ("+log(8/5)", Fraction(8, 5)), ("+log(9/2)", Fraction(9, 2)),
           ("+log(9/4)", Fraction(9, 4)), ("+log(9/5)", Fraction(9, 5))]
ANA6 = ["+log2", "-log2", "+log3", "-log3", "+log6", "+log(3/2)"]


def kapsama_dusuk(c):
    """Dilim [c−½DW, c+½DW) blok b'nin Δω aralığı (0.86L−L_b, 1.30L−L_b] içinde mi?"""
    return [b for b in range(8)
            if (c - 0.5 * DW >= 0.86 * L_D - L_B_D[b] - EPS)
            and (c + 0.5 * DW <= 1.30 * L_D - L_B_D[b] + EPS)]


def kapsama_son(c):
    """Son (τ'-dilim) blok profili: dilim merkezleri τ'∈[0.8625,1.2975];
    c ± ½DW, blok b'nin merkez-Δω aralığı içinde mi (np.interp dolgu yok)."""
    return [b for b in range(8)
            if (c - 0.5 * DW >= 0.8625 * L_S - L_B_S[b] - EPS)
            and (c + 0.5 * DW <= 1.2975 * L_S - L_B_S[b] + EPS)]


def jgrid(h, lo, hi):
    """|j·DW − h| ∈ (lo, hi] (lo<0 → [0, hi]) olan merkezler."""
    j0 = int(np.floor((h - hi) / DW)) - 1
    j1 = int(np.ceil((h + hi) / DW)) + 1
    out = []
    for j in range(j0, j1 + 1):
        d = abs(j * DW - h)
        if d <= hi + EPS and (lo < 0 or d > lo + EPS):
            out.append(j)
    return out


KATALOG = {}
for liste, grup in (("yanmali", YANMALI), ("sonmeli", SONMELI)):
    for ad, n in grup:
        a, b = n.numerator, n.denominator
        h = float(np.log(a) - np.log(b))
        ara = jgrid(h, -1, ARAMA)
        hal = jgrid(h, HALKA_IC, HALKA_DIS)
        a2 = jgrid(h, -1, A2_YARI)
        KATALOG[ad] = {
            "liste": liste, "n": f"{a}/{b}", "a": a, "b": b, "mu_a": mu(a),
            "phi_a": phi(a), "koherens_mu_bolu_phi": mu(a) / phi(a),
            "delta_omega": h,
            "ev_dilimi_j": int(np.floor(h / DW + 0.5)),
            "arama_j": ara, "arama_merkez": [round(j * DW, 4) for j in ara],
            "halka_j": hal, "halka_merkez": [round(j * DW, 4) for j in hal],
            "A2_j": a2,
            "kapsama_dusuk": {str(j): kapsama_dusuk(j * DW) for j in sorted(set(ara + hal + a2))},
            "kapsama_son": {str(j): kapsama_son(j * DW) for j in sorted(set(ara + hal + a2))},
            "ana6": ad in ANA6}

# çakışma notları: bir hedefin arama dilimi başka bir katalog hedefinin ev
# dilimi mi / arama kümeleri kesişiyor mu (veri-öncesi, ızgaradan)
CAKISMA = []
adlar = list(KATALOG)
for i, x in enumerate(adlar):
    for y in adlar[i + 1:]:
        ortak = sorted(set(KATALOG[x]["arama_j"]) & set(KATALOG[y]["arama_j"]))
        if ortak:
            CAKISMA.append({"hedefler": [x, y], "ortak_arama_j": ortak,
                            "ortak_merkez": [round(j * DW, 4) for j in ortak],
                            "listeler": [KATALOG[x]["liste"], KATALOG[y]["liste"]]})
EV_ICINDE = []
for x in adlar:
    for y in adlar:
        if x != y and KATALOG[y]["ev_dilimi_j"] in KATALOG[x]["arama_j"]:
            EV_ICINDE.append({"arayan": x, "ev_dilimi_icerilen": y,
                              "j": KATALOG[y]["ev_dilimi_j"]})

# ---------------- A1 (K2) pencereleri + sınıf sayıları (asal sayma) ----------------
A1 = {
    "+log3": {"tau": [1.078, 1.104], "a": 3, "b": 1, "siniflar": [1, 2]},
    "-log3": {"tau": [0.896, 0.922], "a": 3, "b": 3, "a_teori": 1,
              "siniflar": [1, 2]},
    "+log5": {"tau": [1.121, 1.146], "a": 5, "b": 1, "siniflar": [1, 2, 3, 4]},
    "+log6": {"tau": [1.137, 1.161], "a": 6, "b": 1, "siniflar": [1, 5]},
    "+log2": {"tau": [1.045, 1.070], "a": 2, "b": 1, "siniflar": [1],
              "ek_KAYIT_mod4": ("aynı çizgiler q' mod 4 ∈ {1, 3} diye ikiye "
                                "bölünür (ek kontrol; teori a=2 ⇒ mod-4 sınıfları "
                                "AYNI yönde, her biri |K_r|/|K_top| ≈ 0.5); "
                                "hükme girmez.")},
}
b187 = yukle("b187_192a", QM / "187_configs" / "187b_katman_defteri.py")
t0 = time.time()
QA, LAM, TAUA, AQA = b187.asal_kuvvetler(np.exp(1.30 * L_S), L_S)
for ad, w in A1.items():
    m = (TAUA >= w["tau"][0]) & (TAUA <= w["tau"][1])
    q = QA[m].astype(np.int64)
    a = w["a"]
    say = {}
    for r in w["siniflar"]:
        sel = (q % a == r)
        say[str(r)] = {"n": int(sel.sum()), "sum_a": float(AQA[m][sel].sum())}
    bol = (np.gcd(q, a) > 1)
    say["bolunen"] = {"n": int(bol.sum()), "sum_a": float(AQA[m][bol].sum()),
                      "q": q[bol].tolist()}
    kalan = ~bol
    for r in w["siniflar"]:
        kalan &= ~(q % a == r)
    say["siniflanamayan"] = int(kalan.sum())
    w["cizgi_sayisi"] = int(m.sum())
    w["sinif_sayilari"] = say
    # kinematik: blok b uydu konumu τ'_b = (L_b ± log n)/L_son
    ln = float(np.log(w["a"] if ad != "-log3" else 1 / 3))
    w["blok_uydu_tau"] = [float((lb + ln) / L_S) for lb in L_B_S]
    w["cizgi_alt_ust_q"] = [int(q.min()), int(q.max())]
T_SAY = time.time() - t0

# ---------------- ÖN-KAYIT ----------------
ONKAYIT = {
    "gorev": "192 — UYDULARIN TEORİSİ: tarak uydusu kalıntı sınıfını okur (μ(a)/φ(a))",
    "kalem": "KALEM_UYDU_TEORISI_23EYL2026.md AYNEN",
    "girdi_sha256": {
        "KALEM_UYDU_TEORISI_23EYL2026.md": sha(QM / "KALEM_UYDU_TEORISI_23EYL2026.md"),
        "188b_harita.py": sha(QM / "188_configs" / "188b_harita.py"),
        "187b_katman_defteri.py": sha(QM / "187_configs" / "187b_katman_defteri.py"),
        "185b_oz_muhasebe.py": sha(QM / "185_configs" / "185b_oz_muhasebe.py"),
        "190b_harita.py": sha(QM / "190_configs" / "190b_harita.py"),
        "190c_hukum.py": sha(QM / "190_configs" / "190c_hukum.py"),
        "190e_kesif.py": sha(QM / "190_configs" / "190e_kesif.py"),
        "190/harita_omega_dusuk.npz": sha(S["190"] / "harita_omega_dusuk.npz"),
        "190/profiller_190.npz": sha(S["190"] / "profiller_190.npz"),
        "190/zincir_dusuk/K1_gercek_dusuk.npz": sha(ZD / "K1_gercek_dusuk.npz"),
        "190/zincir_dusuk/OZ_gercek_dusuk.npz": sha(ZD / "OZ_gercek_dusuk.npz"),
        "190/zincir_dusuk/G1_proj_gercek_dusuk.npz": sha(ZD / "G1_proj_gercek_dusuk.npz"),
        "188/harita_proj_gercek.npz": sha(S["188"] / "harita_proj_gercek.npz"),
        "188/harita_K_gercek.npz": sha(S["188"] / "harita_K_gercek.npz"),
        "188/dilim_44.npz": sha(S["188"] / "dilim_44.npz"),
        "184/K1_gercek.npz": sha(S["184"] / "K1_gercek.npz"),
        "185/OZ_gercek.npz": sha(S["185"] / "OZ_gercek.npz"),
        "186/G1_proj_gercek.npz": sha(S["186"] / "G1_proj_gercek.npz"),
        "155/eta_son_t0.4_c4000.npz": sha(S["155"] / "eta_son_t0.4_c4000.npz")},
    "kinematik": {"L_dusuk": L_D, "L_b_dusuk": L_B_D, "L_son": L_S, "L_b_son": L_B_S,
                  "not": "L ve L_b yalnız mid'den (190a / 188f(c) AYNEN); sonuç değil."},
    # ======================= B =======================
    "B": {
        "katalog": KATALOG,
        "ana6": ANA6,
        "profil_birincil": (
            "DÜŞÜK pencere ω-dilim (0.025) Δω profili, 190b D çıktısı "
            "scratchpad/190/omega/b<b>_c<k>.npz (blok b, dilim j: re, im = Σ_{n∈b} "
            "dds_{b,j}(n)·(cos, −sin)(ω_q m_n), q ∈ pencere [0.45,0.86)). HAVUZ "
            "profili, dilimi TAM kapsayan bloklar B_j üzerinden: C_j(q) = 2·Σ_{b∈B_j} "
            "(re+i·im)_{b,j}(q) / Σ_{b∈B_j} n_b; K_j = Σ_{q∈HAVUZ} C_j conj(mix_q) / "
            "Σ|mix_q|² (188b K_matris AYNEN, HAVUZ maskesi 188b bant_maskeleri); "
            "κ(j) = −Re K_j. mix = 188b.karisim(P, OZ, aq, disari)[pencere] "
            "('dusuk' zinciri 190/zincir_dusuk)."),
        "jackknife": (
            "8-blok loo (188b K_hesap AYNEN): replika i için C_j'de blok i "
            "çıkar (B_j \\ {i}, pay/payda yeniden), mix de loo (disari=i); "
            "se = √(7/8·Σ(θ_i − θ̄)²). i ∉ B_j ise yalnız mix loo değişir."),
        "kapsama": (
            f"blok b dilimi TAM kapsar ⇔ [c−½·0.025, c+½·0.025) ⊂ (0.86·L − L_b, "
            f"1.30·L − L_b]; bir dilim ancak |B_j| ≥ {MIN_KAPSAMA} ise kullanılır "
            "(arama dilimlerinden biri < 4 ise hedef 'erişilemedi'; halkada < 4 "
            "olan dilim halkadan düşer). Kapsama listeleri katalogda (kinematik)."),
        "profil_ikincil_son": (
            "SON penceresi (188 dosyaları; yeniden koşu YOK): blok başına τ'-dilim "
            "(0.005) κ_b(s) = −Re[C_b(s)·conj(mix)]/Σ|mix|², C_b = 2(re+i·im)_b/n_b "
            "(188/harita_proj_gercek.npz; 190c tau_blok_profili AYNEN), Δω_b(s) = "
            "τ'_s·L_son − L_b,son; yoğunluk κ_b/(0.005·L_son) 0.025 ızgarasına "
            "np.interp (190e son_yog AYNEN), dilimi kapsayan bloklarda n_b-ağırlıklı "
            "ortalama. mix: 184/185/186 son zinciri, 188b.karisim (loo). Aynı "
            "tepe / taban / jk kuralları. HÜKME GİRMEZ (KAYIT)."),
        "tepe": (
            "arama kümesi = merkezi |c_j − h| ≤ 0.03 (+1e-9) olan dilimler "
            "(|B_j| ≥ 4); j* = argmax κ(j); κ_tepe = κ(j*); se = κ(j*)'nin jk "
            "se'si (j* tam veride seçilir, replikalarda sabit)."),
        "yerel_taban": (
            "halka = 0.05 < |c_j − h| ≤ 0.15 (+1e-9) olan dilimler (|B_j| ≥ 4); "
            "taban_medyan = median κ(halka); yayılım = 1.4826·median|κ(halka) − "
            "taban_medyan| (gürbüz σ); üst = taban_medyan + 2·yayılım; "
            "taban bandı = [medyan − 2·yayılım, medyan + 2·yayılım]."),
        "yanar": "κ_tepe > 3·se VE κ_tepe > üst (yerel tabanın üstünde)",
        "soner": ("|κ_tepe| < 2·se YA DA |κ_tepe − taban_medyan| ≤ 2·yayılım "
                  "(yerel taban içinde)"),
        "belirsiz": ("ikisi de değil (ör. üstte ama 2–3 se; ya da tabanın ALTINDA "
                     "'çukur': κ_tepe < medyan − 2·yayılım ve |κ_tepe| ≥ 2·se — "
                     "çukur ayrıca işaretlenir)"),
        "yanar_4se": "κ_tepe > 4·se VE κ_tepe > üst (H-192a ölüm ölçütü)",
        "cakisma_notu": {
            "arama_kumesi_kesisimleri": CAKISMA,
            "ev_dilimi_arama_icinde": EV_ICINDE,
            "okuma": ("Veri-öncesi ızgara notu: listelenen çiftlerde ±0.03 arama "
                      "kümeleri ortak dilim içerir; hüküm LİTERAL kuralla verilir "
                      "(ortak dilim her iki hedefin aramasında kalır). KAYIT: "
                      "çakışmasız arama (başka bir katalog hedefinin ev dilimi "
                      "aramadan çıkarılır) yan yana raporlanır; HÜKMÜ DEĞİŞTİRMEZ.")},
        "KAYIT": ("κ_tepe/se, taban (medyan, yayılım), arama kümesindeki min κ, "
                  "çakışmasız arama sonucu, yalın 'κ_tepe > 4·se' sayımı, son "
                  "penceresi aynı tablo.")},
    # ======================= A1 =======================
    "A1": {
        "pencere": "SON (L_son = 12.029593…; 188 makinesi AYNEN)",
        "uydular": A1,
        "sinif_tanimi": (
            "çizgi q' (asal kuvvet, 187b.asal_kuvvetler(e^{1.30·L_son}, L_son) AYNEN), "
            "τ' = log q'/L_son ∈ [lo, hi] (iki uç dahil); sınıf r = q' mod a "
            "(+log3, −log3: a=3, r∈{1,2}; +log5: a=5, r∈{1,2,3,4}; +log6: a=6, "
            "r∈{1,5}; +log2: a=2, r=1 (tek)); gcd(q', a) > 1 olanlar 'bölünen' "
            "(a'nın asal çarpanlarının kuvvetleri) ayrı sınıf. Toplam = tüm "
            "penceredeki çizgiler (bölünen dahil)."),
        "seri": ("188b.seri_ve_G AYNEN, ALT-ÖRNEKLEM n = 0,3,6,… (100 000; 188 ile "
                 "aynı); bid = orijinal indisin jk bloğu; sınıf başına ayrı seri."),
        "K_r": ("izdüşüm 188b.izdusum AYNEN (alt-örneklem, pencere çizgileri 184 "
                "K1_gercek τ∈[0.45,0.86)); K_r = 188b.K_matris(c_proj(re,im,nb,N,"
                "disari), karisim(P,OZ,aq,disari)[pencere], maskeler)[HAVUZ] — "
                "Σ_{q∈HAVUZ} c^{(r)}_q conj(mix_q)/Σ|mix_q|²; loo 8 replika."),
        "acilar": ("θ_r = arg K_r (derece); göreli θ_r − θ_top (sarılı, (−180,180]); "
                   "sınıf farkı Δ_{r→s} = wrap(arg K_s − arg K_r); jk se: sarılı "
                   "replika sapmaları (188b se_aci deseni)."),
        "oran": "|K_r| / |K_top| (+ jk se)",
        "kosullar": {
            "+log3": "|Δ_{1→2}| ∈ [105°, 135°] VE |K_1|/|K_top|, |K_2|/|K_top| ∈ [0.7, 1.4]",
            "-log3": "|Δ_{1→2}| < 15° VE |K_1|/|K_top|, |K_2|/|K_top| ∈ [0.35, 0.65]",
            "+log5 (KAYIT, H-192b dışı)": ("ardışık |Δ_{r→r+1}| ∈ [57°, 87°] (r=1,2,3) "
                                           "VE dört oran ∈ [0.6, 1.5]"),
            "+log6 (KAYIT, H-192b dışı)": ("|Δ_{1→5}| ∈ [105°, 135°]; yan KAYIT: "
                                           "θ_1 − θ_top, θ_5 − θ_top ≈ ∓60° (±15°)"),
            "+log2 (KAYIT, kontrol)": "tek sınıf: |θ_1 − θ_top| < 15° ve oran ≈ 1 (dönme yok)"},
        "yon_tutarliligi_KAYIT": (
            "faz 2π q' b/a ⇒ ortak bir işaret s = ±1 ile: Δ_{1→2}(+log3) = s·120°, "
            "Δ_{r→r+1}(+log5) = s·72°, Δ_{1→5}(+log6) = −s·120° (mod 360). "
            "KAYIT: işaretlerin bu desene uyup uymadığı."),
        "makine_muhuru": (
            "188/dilim_44.npz (τ' ∈ (1.080, 1.085], +log3 penceresi içinde): aynı "
            "asal-kuvvet listesi ve seri_ve_G ile dds, Gre, Gim, q, a, τ BİT-BİT; "
            "ayrıca tek-dilim K_HAVUZ vs 188/harita_K_gercek.npz K[HAVUZ, 44] "
            "(maks |Δ|; BLAS toplama sırası yüzünden ~1e-16 beklenir)."),
        "dogrusallik": ("dilim 44 çizgilerinin sınıf serileri toplamı = dilim 44 "
                        "serisi (maks|Δ|, kayan-nokta)"),
        "tam_orneklem_kontrolu": (
            "+log3 r=1 ve −log3 r=1 sınıfları TAM örneklemde (N = 299 999) "
            "yeniden: |K_alt − K_tam| / se_alt (KAYIT; 188b C deseni)."),
        "sure_notu": f"asal sayma {T_SAY:.0f}s"},
    # ======================= A2 =======================
    "A2": {
        "pencere": "DÜŞÜK (birincil profil, B ile aynı havuz profili ve jk)",
        "parlaklik": ("P(h) = Σ κ(j), |c_j − h| ≤ 0.0625 (+1e-9) → 5 dilim "
                      "(|B_j| ≥ 4); ham toplam (188 pencere dökümü deseni)."),
        "oranlar": {
            "+log3/+log2": 0.5, "+log5/+log2": 0.25, "+log6/+log2": 0.5,
            "+log7/+log2": 1 / 6, "+log10/+log2": 0.25,
            "-log2/+log2": 0.5, "-log3/+log3": 2 / 3},
        "oran_kalem_yazimi": {"+log7/+log2": 0.17, "-log3/+log3": 0.67},
        "bant": "[öngörü/2, öngörü·2] (×/÷ 2); oran se: jk (replika oranı)",
        "KAYIT": ("taban-düzeltmeli P' = Σ(κ(j) − taban_medyan(h)); tepe oranı "
                  "κ_tepe(h)/κ_tepe(+log2); SON penceresi aynı oranlar "
                  "(yoğunluk × 0.025).")},
    # ======================= HİPOTEZLER =======================
    "H_192a": {
        "hipotez": "SEÇİM KURALI: μ(a)=0 uyduları söner, μ(a)≠0 ana-6 yanar (DÜŞÜK).",
        "mühür": "9 sönmelinin HEPSİ 'söner' VE ana-6'nın HEPSİ 'yanar'",
        "olum": "sönmelilerden herhangi biri 'yanar_4se' (κ_tepe > 4·se VE üstte)",
        "ara": "aksi halde KAYIT (kısmi; hangi kalemin tutmadığı yazılır)",
        "son": "son penceresi aynı kuralla KAYIT; hükme girmez"},
    "H_192b": {
        "hipotez": "KALINTI DÖNMESİ: +log3 sınıfları 120° ayrık, −log3 sınıfları aynı yönde.",
        "mühür": "+log3 koşulu VE −log3 koşulu (A1.kosullar) ikisi birden",
        "olum": "+log3'te |Δ_{1→2}| < 45° YA DA −log3'te |Δ_{1→2}| > 45°",
        "ara": "aksi halde KAYIT",
        "not": "nokta kestirimleriyle karar; se'ler raporlanır"},
    "H_192c": {
        "hipotez": "PARLAKLIK |μ(n)|/φ(n) (artı) ve 1/m (eksi) — A2 bantları (DÜŞÜK).",
        "hukum": "KAYIT (bant sınavı): 7 oranın kaçının bantta olduğu yazılır"},
    "kurallar": ("TEK DALGA; ölümler kurtarmasız; eşik gevşetme yok; ölçülemeyen "
                 "'erişilemedi'; git'e dokunulmaz; derin ikiz inşasına (191) "
                 "dokunulmaz; nice -n 19, tek süreç, OMP/BLAS=1; sonuç ORTAK TEFTİŞE."),
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
    yol = S["192"] / "ONKAYIT_192.json"
    kuru = "--kuru" in sys.argv          # kuru koşu: yazmadan göster
    if yol.exists() and not kuru:
        raise SystemExit(f"ONKAYIT_192 ZATEN VAR — üzerine yazılmaz: {yol}")
    if not kuru:
        json.dump(ONKAYIT, open(yol, "w"), indent=1, ensure_ascii=False)
    print("=" * 78)
    print("192a / K0 ÖN-KAYIT " + ("KURU KOŞU (yazılmadı)" if kuru else "yazıldı"))
    print("=" * 78)
    print(f"  sha256 = {s}")
    print(f"  damga  = {ONKAYIT['zaman']}")
    print(f"  L_dusuk = {L_D:.9f}  L_son = {L_S:.9f}")
    print(f"  L_b_son = {np.round(L_B_S, 4).tolist()}")
    for ad, k in KATALOG.items():
        kd = [len(k['kapsama_dusuk'][str(j)]) for j in k['arama_j']]
        print(f"  {k['liste']:>7} {ad:>10} n={k['n']:>4} μ(a)={k['mu_a']:+d} "
              f"Δω={k['delta_omega']:+.4f} arama={k['arama_merkez']} "
              f"kapsama(düşük)={kd} halka={len(k['halka_j'])}")
    print(f"  çakışmalar: {CAKISMA}")
    print(f"  ev-dilimi arama içinde: {EV_ICINDE}")
    for ad, w in A1.items():
        print(f"  A1 {ad}: τ'{w['tau']} çizgi {w['cizgi_sayisi']} "
              f"sınıflar { {k: v['n'] if isinstance(v, dict) else v for k, v in w['sinif_sayilari'].items()} }"
              f" bölünen q={w['sinif_sayilari']['bolunen']['q']}")
        print(f"      blok uydu τ' = {np.round(w['blok_uydu_tau'], 4).tolist()}")
    print(f"  -> {yol}")
