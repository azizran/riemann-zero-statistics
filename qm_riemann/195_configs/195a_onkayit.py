# -*- coding: utf-8 -*-
"""
195a — K0d ÖN-KAYIT (ADA ÖLÇÜMLERİNDEN ÖNCE)
===========================================
KALEM_UYDU_KARAKTERI_24EYL2026 (commit 50df21e) AYNEN. K0a/K0b/K0c sonuçlarını
okur; biri KALDIYSA yazmayı REDDEDER (DUR). ONKAYIT_195.json zaten varsa
yazmayı REDDEDER. Yalnız kinematik (sıfır konumları → L, bloklar) ve K0
dosyaları okunur; ada A/B ölçümü YOK.

Donanlar: kalem öngörü tablosu + asal-çarpan kuralıyla tam yasak/izinli
matrisi; tanımlar (seri, izdüşüm, öz-terim, K̃ havuzu, jackknife); blok/bant
tanımları ve ada blok listeleri (kinematik); A ızgarası/pencere/halka;
B uydu listesi (−log5 pencere-içi dışlaması, aritmetik); H-195A/B1/B2/B3
karar kuralları; R1/R2/R3 imzaları; İŞARET REFERANSI (K0c); L-eşli zeta
bant referansları (K0c); K0a/b/c özetleri; ajan veri-öncesi notları (hükme
girmez); girdi sha256'ları; betiğin kendi sha256'sı + damga.
"""
import hashlib
import importlib.util
import json
import time
from math import gcd
from pathlib import Path

import numpy as np
from scipy.special import j0 as J0

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("o195", HERE / "195o_ortak.py")
o = importlib.util.module_from_spec(spec)
spec.loader.exec_module(o)

YOL = o.S195 / "ONKAYIT_195.json"


def phi(n):
    return sum(1 for r in range(1, n + 1) if gcd(r, n) == 1)


def log(*a):
    print(*a, flush=True)


if __name__ == "__main__":
    if YOL.exists():
        raise SystemExit(f"ONKAYIT_195.json ZATEN VAR — yazma reddedildi ({YOL})")
    K0ab = json.load(open(o.S195 / "K0ab.json"))
    K0k = json.load(open(o.S195 / "K0a_kusur.json"))
    K0c = json.load(open(o.S195 / "K0c.json"))
    if not K0ab["K0b_GECTI"]:
        raise SystemExit("K0b KALDI — DUR")
    if not K0k["K0a_duzeltilmis_GECTI"]:
        raise SystemExit("K0a (düzeltilmiş) KALDI — DUR")
    for ad in o.ADALAR:
        if any(kk["L"] >= o.L_UST for kk in K0k[ad]["kusurlar"]):
            raise SystemExit(f"K0a: {ad} dosya kusuru ÜST BÖLGEDE — DUR")
    for p in ["zeta_son", "zeta_dusuk"]:
        if not K0ab["K0a"][p]["kapi_0.05"]:
            raise SystemExit(f"K0a {p} KALDI — DUR")
    if not K0c["K0c_GECTI"]:
        raise SystemExit("K0c KALDI — DUR")

    # ---------------- kalem öngörü tablosu ----------------
    kalem = {
        "beta": {"k": 4, "parite": "tek", "sqrtk_phik": 1.000,
                 "sonen": ["log2", "log6", "log10", "log(3/2)"],
                 "yanan": ["±log3", "±log5"]},
        "chi3": {"k": 3, "parite": "tek", "sqrtk_phik": 0.866,
                 "sonen": ["log3", "log6", "log(3/2)"],
                 "yanan": ["±log2", "±log5", "log10"]},
        "chi5e": {"k": 5, "parite": "çift", "sqrtk_phik": 0.559,
                  "sonen": ["log5", "log10"], "yanan": ["±log2", "±log3", "log6"]},
        "chi8e": {"k": 8, "parite": "çift", "sqrtk_phik": 0.707,
                  "sonen": ["log2", "log6", "log10"], "yanan": ["±log3", "±log5"]},
        "chi8o": {"k": 8, "parite": "tek", "sqrtk_phik": 0.707,
                  "sonen": ["log2", "log6", "log10"], "yanan": ["±log3", "±log5"]},
    }
    matris = {}
    for ad in o.ADALAR:
        k = kalem[ad]["k"]
        assert abs(np.sqrt(k) / phi(k) - kalem[ad]["sqrtk_phik"]) < 5e-4
        matris[ad] = {u: ("YASAK" if o.yasak_mi(u, k) else "izinli") for u in o.UYDU_SIRA}
        # kalem listeleriyle tutarlılık (kalem ± yazımı; '+' ön ekli adlar)
        for s in kalem[ad]["sonen"]:
            assert matris[ad]["+" + s] == "YASAK", (ad, s)
        for s in kalem[ad]["yanan"]:
            if s.startswith("±"):
                for sg in "+-":
                    assert matris[ad][sg + s[1:]] == "izinli", (ad, s)
            else:
                assert matris[ad]["+" + s] == "izinli", (ad, s)

    # ---------------- ada blokları (kinematik) ----------------
    bloklar = {}
    for ad in o.ADALAR:
        K = o.ada_kin(ad)
        bl = o.dL_bloklari(K["Lm"], o.L_UST)
        Lb = [float(K["Lm"][idx].mean()) for (_, _, _, idx) in bl]
        bant = [j for (j, lo, hi, idx) in bl if lo >= o.BANT[0] - 1e-9 and lo < o.BANT[1] - 1e-9]
        bloklar[ad] = {"n_blok": len(bl), "j": [int(x[0]) for x in bl],
                       "kenar_lo": [x[1] for x in bl], "kenar_hi": [x[2] for x in bl],
                       "L_b": Lb, "N_b": [int(len(x[3])) for x in bl],
                       "N_ust_toplam": int(sum(len(x[3]) for x in bl)),
                       "L_max": float(K["Lm"].max()),
                       "bant_bloklari_j": [int(x) for x in bant],
                       "bant_U": float(min(o.BANT[1], K["Lm"].max())),
                       "bant_N": int(sum(len(x[3]) for x in bl if x[0] in bant))}
        # −log5 pencere-içi aritmetiği; diğer uyduların pencere-dışı olduğu
        ic = {}
        for u in o.UYDU_SIRA:
            lo_t = [(L + o.log_n(u) - o.DELTA) / L for L in Lb]   # en küçük τ'
            ic[u] = bool(min(lo_t) < o.TAU_HI)
        bloklar[ad]["uydu_pencere_ici_mi"] = ic
        log(f"  {ad:>6}: {len(bl)} blok, L_b {Lb[0]:.3f}…{Lb[-1]:.3f}, N_üst={bloklar[ad]['N_ust_toplam']}, "
            f"bant blokları {bant} (U={bloklar[ad]['bant_U']:.4f}, N={bloklar[ad]['bant_N']}); "
            f"pencere-içi uydular: {[u for u, v in ic.items() if v]}")
    B_DISLANAN = sorted({u for ad in o.ADALAR for u, v in bloklar[ad]["uydu_pencere_ici_mi"].items() if v})
    assert B_DISLANAN == ["-log5"], B_DISLANAN

    # ---------------- A halka/pencere nokta sayıları ----------------
    A_maske = {}
    for u in o.A_LISTE + o.A_KAYIT:
        W, R = o.pencere_halka_maskeleri(u)
        A_maske[u] = {"pencere_nokta": int(W.sum()), "halka_nokta": int(R.sum())}
        assert R.sum() >= 60, (u, "halka çok dar")

    # ---------------- işaret referansı ----------------
    ISR = K0c["ISARET_REFERANSI"]

    # ---------------- ajan veri-öncesi notları (hükme GİRMEZ) ----------------
    tab, _ = o.karakter_tablolari()
    notlar = {}
    for ad in o.ADALAR:
        k = kalem[ad]["k"]
        pj = 1.0
        for p in sorted(o.asal_carpanlar(k)):
            for j in range(1, 40):
                q = float(p) ** j
                a = np.log(p) / (np.pi * np.sqrt(q) * np.log(q))
                pj *= J0(2 * np.pi * a)
        notlar[ad] = {"kalem_sqrtk_phik": float(np.sqrt(k) / phi(k)),
                      "alt1_k_phik": float(k / phi(k)),
                      "Pi_J0_p_bolen_k": float(pj),
                      "alt2_k_phik_bolu_PiJ0": float(k / phi(k) / pj)}

    # ---------------- girdi sha256 ----------------
    girdiler = [o.QM / "KALEM_UYDU_KARAKTERI_24EYL2026.md",
                o.QM / "118a_taperli_L_motoru.py", o.QM / "118b_yedi_ada_taperli_dokum.py",
                o.QM / "105b_yeni_ada_kampanya.py", o.QM / "128_odl_zeros6_2e6_zeros.npz",
                o.S155 / "eta_son_t0.4_c4000.npz", o.ZD / "eta_dusuk_t0.4_c4000.npz",
                o.S188 / "dilim_44.npz", o.S195 / "K0ab.json", o.S195 / "K0a_kusur.json",
                o.S195 / "K0c.json", o.S195 / "K0c_bloklar.npz",
                HERE / "195o_ortak.py", HERE / "195k0ab_kapilar.py",
                HERE / "195k0a2_kusur.py", HERE / "195k0c_zeta.py"] + \
        [o.QM / f"118b_{ad}_zeros.npz" for ad in o.ADALAR]
    girdi_sha = {str(p.name): o.sha(p) for p in girdiler}

    ONK = {
        "gorev": "195 — UYDULARIN KARAKTERİ: L-fonksiyonu adalarında tarak iptali",
        "kalem": "KALEM_UYDU_KARAKTERI_24EYL2026.md (commit 50df21e; türetim teftişinden geçti). "
                 "İŞARET: hipotezler 'K̃_ζ ile AYNI işaret' ile ifade edilir (kalemin B "
                 "paragrafındaki 'K̃<0' yerine hipotezler bölümü geçerli — teftiş (g)).",
        "ongoru_tablosu_kalem": kalem,
        "yasak_izinli_matrisi_asal_carpan_kurali": matris,
        "tanimlar": {
            "L_chi": "log(k t/2π)", "C_chi": "−χ(−1)/8 (çift −1/8, tek +1/8); zeta 7/8",
            "a_q_chi": "χ(q)·Λ(q)/(π√q log q), q asal-kuvvet; χ(q)=0 çizgisi yok",
            "seri": "s(n) = Σ_{q'} 2a_{q'}^χ sin(ω'g_n/2) cos(ω'm_n) (188b.seri_ve_G dds yolu; K0c M1 bit-bit)",
            "izdusum": "c_{q,b} = 2⟨s·e^{−iω_q m}⟩_b (190b.izdusum_blok; K0c M2a bit-bit)",
            "oz_terim": "c_{q,b}^öz = 2⟨2a_q^χ sin(ω_q g/2) cos(ω_q m) e^{−iω_q m}⟩_b (185b.c_oz; K0c M2b 1.5e-15)",
            "pencere_cizgileri": "blok b: τ = log q / L_b ∈ [0.45, 0.86), χ(q) ≠ 0",
            "uydu_cizgileri": "blok b: |log q' − L_b − log n| < 0.03, χ(q') ≠ 0",
            "K_tilde": "K̃ = Σ_b N_b·num_b / Σ_b N_b·den_b; num_b = Σ_q c^{(uydu)}_{q,b} conj(c^öz_{q,b}); "
                       "den_b = Σ_q |c^öz_{q,b}|²; N_b = bloğun TAM nokta sayısı",
            "jackknife": "blok-loo (tek blok dışarıda), se = √((B−1)/B·Σ(θ_j−θ̄)²)",
            "isaretli_istatistik": "Re K̃ (z = Re K̃/se_jk(Re K̃)); Im K̃ ve açı KAYIT (R2 imzası)",
            "L_b": "blok içi mean L_χ(m_n)",
            "orneklem": "adalar ve zeta düşük: TAM; zeta son: seri/izdüşüm ALT=3 (193 AYNEN), öz-terim tam blok",
        },
        "bloklar_tanim": {
            "adalar": f"üst bölge L_χ ≥ {o.L_UST}; eşit ΔL = {o.DL} ızgarası e_j = 8.5 + 0.05 j; "
                      f"üstteki kısmi blok ancak L-genişliği ≥ {o.KISMI_MIN} ise",
            "zeta_son": "193 blokları AYNEN (8 eşit sayım; L_b 193 ön-kaydı)",
            "zeta_dusuk_K0c_oran": "190 blokları AYNEN (8 eşit sayım)",
            "zeta_dusuk_izgara": "adalarla AYNI 8.5-çapalı ΔL=0.05 ızgarası (17 blok) — A kontrolü ve B3 referansı",
            "bant": f"L-eşli bant [{o.BANT[0]}, {o.BANT[1]}): ızgara blokları j=30..37; ada verisi 10.4'ten önce "
                    "biterse üst sınır U_χ = L_max,χ (kısmi blok) ve zeta referansı AYNI kısmi blokla",
        },
        "ada_bloklari": bloklar,
        "A_tanim": {
            "G_hat": "Ĝ_b(Δω) = mean_{n∈b} e^{i(L_b+Δω)m_n} (MUTLAK m_n)",
            "izgara": "Δω ∈ [−2.0, 2.6], adım 0.001",
            "pencere": f"|Δω − log n| < {o.PENCERE_YARI}",
            "halka": f"{o.HALKA_IC} ≤ |Δω − log n| ≤ {o.HALKA_DIS}, a,b ≤ {o.KATALOG_H} sade rasyonellerin "
                     f"(Δω=0 dahil, kendisi hariç) ±{o.DISLAMA} çevresi dışlanır",
            "istatistik": "D_b = mean_pencere|Ĝ_b|² − mean_halka|Ĝ_b|²; D = Σ N_b D_b/Σ N_b; z = D/se_jk(D)",
            "uydular": o.A_LISTE, "kayit_uydular": o.A_KAYIT, "nokta_sayilari": A_maske,
            "zeta_kontrol": "son (193 blokları) ve düşük (ΔL ızgarası) aynı ölçümle — KAYIT",
        },
        "B_tanim": {
            "uydular": [u for u in o.UYDU_SIRA if u not in B_DISLANAN],
            "dislanan": {"-log5": "uydu çizgileri pencere içinde (τ' = log q'/L_b < 0.86 her ada bloğunda) "
                                  "⇒ öz-terim kirliliği; B'de ERİŞİLEMEDİ, A'da ölçülür"},
        },
        "hipotezler": {
            "H_195A": {
                "tanim": "her ada: yasak uydular (n'nin bir asal çarpanı k'yı böler; A listesi) z < 2; "
                         "izinli güçlüler (±log2, ±log3 izinliyse) z > 3",
                "ada_MUHUR": "tüm yasak z < 2 VE tüm izinli-güçlü z > 3",
                "OLUM": "herhangi bir adada herhangi bir yasak uydu z ≥ 4",
                "genel": "ÖLÜM varsa ÖLDÜ; beş ada da MÜHÜR ise MÜHÜR; aksi KAYIT",
                "kayit": "izinli-diğer (+log4, ±log5, +log6, +log10) z; +log(3/2); zeta kontrolleri",
            },
            "H_195B1": {
                "izinli_birincil": [["chi3", "+log2"], ["chi5e", "+log2"], ["chi5e", "+log3"],
                                    ["beta", "+log3"]],
                "yasak_birincil": [["chi3", "+log3"], ["chi5e", "+log5"], ["beta", "+log2"]],
                "referans_isaret": "sign Re K̃_ζ,son (K0c; aşağıda); düşük ile uyumsuzsa o uydu ERİŞİLEMEDİ",
                "MUHUR": "tüm izinli-birincil: işaret = referans VE |z_Re| ≥ 3; tüm yasak-birincil |z_Re| < 2",
                "OLUM": "izinli-birincillerden biri |z_Re| < 2 (R1 lehine) YA DA referansa ZIT işaretli |z_Re| ≥ 3",
                "aksi": "KAYIT",
            },
            "H_195B2": {
                "birincil": {"beta": ["+log3"], "chi3": ["+log2"], "chi8o": ["+log3"],
                             "chi5e": ["+log2", "+log3"], "chi8e": ["+log3"]},
                "tek": ["beta", "chi3", "chi8o"], "cift": ["chi5e", "chi8e"],
                "i_isaret": "tüm birinciller (tek ve çift) referans işaretli VE |z_Re| ≥ 2",
                "ii_buyukluk": "ρ₈ = Re K̃(χ₈ₒ,+log3)/Re K̃(χ₈ₑ,+log3) ∈ [0.60, 1.40]",
                "MUHUR": "(i) VE (ii)",
                "OLUM_R2": "tüm tek-birinciller |z_Re| < 2 VE tüm çift-birinciller referans işaretli |z_Re| ≥ 3",
                "aksi": "KAYIT",
                "kayit": "tek adalarda Im K̃ z-değerleri (R2: 90° faz ⇒ Im ≠ 0, Re ≈ 0)",
            },
            "H_195B3": {
                "ciftler": "B2 birincilleri (6 çift)",
                "bant": "ada: bant blokları havuzu; zeta: K0c 'bant_ref' (aynı bloklar, düşük pencere)",
                "olculebilir": "|z_Re(ada, bant)| ≥ 3 VE |z_Re(ζ, bant)| ≥ 3",
                "oran": "R = Re K̃_χ,bant / Re K̃_ζ,bant; se_R = |R|·√((se_χ/K_χ)² + (se_ζ/K_ζ)²)",
                "ongoru": "√k/φ(k); bant [0.65, 1.35]·öngörü",
                "MUHUR": "≥1 ölçülebilir çift VE tüm ölçülebilirler bant İÇİNDE",
                "OLDU": "bir ölçülebilir çift bant DIŞINDA VE banda uzaklığı ≥ 2·se_R",
                "aksi": "KAYIT (dışarıda ama < 2 se_R)",
                "erisilemedi": "hiç ölçülebilir çift yoksa",
            },
        },
        "rakip_resimler": {
            "R1_iletkensiz": "adalarda iptal YOK: izinli birinciller |z_Re| < 2",
            "R2_paritesiz": "tek adalarda Re ≈ 0 (faz 90°, |Im| anlamlı), çiftlerde var",
            "R3_kaptan": "her iki paritede K̃_ζ ile aynı işaret, √k/φ(k) ölçekli",
        },
        "ISARET_REFERANSI_K0c": ISR,
        "zeta_bant_referansi_K0c": K0c["Y_yeni_Ktilde"]["bant_ref"],
        "zeta_Ktilde_K0c": {k_: K0c["Y_yeni_Ktilde"][k_] for k_ in ["son", "dusuk190", "dusuk_izgara_tum"]},
        "K0a": {"ozet": {ad: {"j0_bagimsiz": K0k[ad]["j0_bagimsiz"],
                              "kusurlar": K0k[ad]["kusurlar"],
                              "duzeltilmis_TUM": K0k[ad]["duzeltilmis_ort_TUM"],
                              "duzeltilmis_UST": K0k[ad]["duzeltilmis_ort_UST"],
                              "literal_duzeltmesiz": K0k[ad]["literal_ort_TUM_duzeltmesiz"],
                              "ham_ortalama_yuvarlama_j0": K0ab["K0a"][ad]["j0_birincil_ortalama_yuvarlama"],
                              "ham_ortalama_yuvarlama_r": K0ab["K0a"][ad]["ort_Nbar_eksi_n_TUM"]}
                         for ad in o.ADALAR},
                "zeta": {p: K0ab["K0a"][p] for p in ["zeta_son", "zeta_dusuk"]},
                "hukum": "GEÇTİ (C_χ kalibrasyonu); χ₃ ve χ₈ₒ dosyalarında L≈4.8'de birer eksik yakın "
                         "sıfır çifti (üst bölge dışında; A/B n kullanmaz)"},
        "K0b": {"hukum": "GEÇTİ", "adalar": {ad: K0ab["K0b"][ad] for ad in o.ADALAR}},
        "K0c": {"M1_seri_bitbit": K0c["M1_seri_bitbit"], "M2a": K0c["M2a_izdusum_bitbit"],
                "M2b": K0c["M2b_coz_goreli_maks_fark"],
                "R_193_yeniden_uretim": K0c["R_193_yeniden_uretim"],
                "oran_kapisi": K0c["G_oran_kapisi"], "hukum": "GEÇTİ"},
        "ajan_veri_oncesi_notlari_HUKME_GIRMEZ": {
            "metin": "(i) Durağan-faz genliği √(2πt*)/Δt_blok: L-eşli kıyasta t*_χ = t*_ζ/k ve blok "
                     "Δt = t·ΔL ⇒ blok-ortalamalı Ĝ çizgi başına √k büyür ⇒ K̃_χ/K̃_ζ ≈ k/φ(k) "
                     "(kalem √k/φ(k) yazıyor). (ii) Jacobi–Anger arka plan çarpımı: χ(q)=0 "
                     "çizgilerinin J_0(2πa_q) çarpanları adada 1 olur ⇒ uydu genliği ek "
                     "1/Π_{p|k,j} J_0(2πa_{p^j}). Bunlar kaptana ORTAK TEFTİŞ için not; "
                     "H-195B3 kalemdeki gibi (√k/φ(k)) sınanır; alternatifler yalnız KAYIT.",
            "sayilar": notlar},
        "girdi_sha256": girdi_sha,
        "kurallar": "ölümler kurtarmasız; eşik gevşetme yok; ölçülemeyene 'erişilemedi'; "
                    "ajan git'e DOKUNMAZ; sonuç ORTAK TEFTİŞE",
    }
    ONK["zaman"] = time.strftime("%a %b %d %H:%M:%S %z %Y")
    ONK["sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    json.dump(ONK, open(YOL, "w"), indent=1, ensure_ascii=False)
    log(f"-> {YOL}\n   sha256 {ONK['sha256']}\n   damga {ONK['zaman']}")
