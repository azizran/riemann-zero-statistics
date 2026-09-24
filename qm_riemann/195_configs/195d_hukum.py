# -*- coding: utf-8 -*-
"""
195d — HÜKÜM (ONKAYIT_195 kuralları AYNEN; eşik gevşetme yok, kurtarma yok)
=========================================================================
Girdi: ONKAYIT_195.json (sha denetimi), A_guc.json, B_cekirdek.json.
Çıktı: scratchpad/195/HUKUM_195.json + ekran tabloları.
"""
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("o195", HERE / "195o_ortak.py")
o = importlib.util.module_from_spec(spec)
spec.loader.exec_module(o)


def log(*a):
    print(*a, flush=True)


if __name__ == "__main__":
    ONK = json.load(open(o.S195 / "ONKAYIT_195.json"))
    s = hashlib.sha256((HERE / "195a_onkayit.py").read_bytes()).hexdigest()
    assert s == ONK["sha256"], "ON-KAYIT SHA UYUMSUZ"
    A = json.load(open(o.S195 / "A_guc.json"))
    B = json.load(open(o.S195 / "B_cekirdek.json"))
    assert A["sha_onkayit"] == ONK["sha256"] and B["sha_onkayit"] == ONK["sha256"]
    ISR = ONK["ISARET_REFERANSI_K0c"]
    ZB = ONK["zeta_bant_referansi_K0c"]
    notlar = ONK["ajan_veri_oncesi_notlari_HUKME_GIRMEZ"]["sayilar"]
    H = {"sha_onkayit": ONK["sha256"]}

    # ======================= H-195A =======================
    log("H-195A — YAPI ÇARPANI GÜCÜ (z = D/se_jk):")
    ada_huk, olum_A = {}, False
    for ad in o.ADALAR:
        k = A["sonuc"][ad]["k"]
        res = A["sonuc"][ad]["uydular"]
        yasak = [u for u in o.A_LISTE if o.yasak_mi(u, k)]
        guclu = [u for u in ["+log2", "-log2", "+log3", "-log3"] if not o.yasak_mi(u, k)]
        y_ok = all(res[u]["z"] < 2 for u in yasak)
        g_ok = all(res[u]["z"] > 3 for u in guclu)
        olum = any(res[u]["z"] >= 4 for u in yasak)
        olum_A |= olum
        h = "ÖLÜM" if olum else ("MÜHÜR" if (y_ok and g_ok) else "KAYIT")
        ada_huk[ad] = {"yasak": {u: res[u]["z"] for u in yasak},
                       "izinli_guclu": {u: res[u]["z"] for u in guclu},
                       "izinli_diger_KAYIT": {u: res[u]["z"] for u in o.A_LISTE
                                              if u not in yasak and u not in guclu},
                       "kayit_3_2": res["+log(3/2)"]["z"], "hukum": h,
                       "yasak_maks_z": max(res[u]["z"] for u in yasak),
                       "guclu_min_z": min(res[u]["z"] for u in guclu)}
        log(f"  {ad:>6}: yasak z maks {ada_huk[ad]['yasak_maks_z']:+.2f} (<2?), "
            f"izinli-güçlü z min {ada_huk[ad]['guclu_min_z']:+.2f} (>3?) → {h}")
    genel_A = "ÖLDÜ" if olum_A else ("MÜHÜR" if all(v["hukum"] == "MÜHÜR"
                                                    for v in ada_huk.values()) else "KAYIT")
    H["H_195A"] = {"adalar": ada_huk, "hukum": genel_A,
                   "zeta_kontrol_KAYIT": {p: {u: A["sonuc"][p]["uydular"][u]["z"]
                                              for u in o.A_LISTE + o.A_KAYIT}
                                          for p in ["zeta_son", "zeta_dusuk"]}}
    log(f"  H-195A: {genel_A}")

    # ======================= H-195B1 =======================
    log("\nH-195B1 — ÇEKİRDEK VARLIK + İŞARET (üst bölge; Re K̃, z_Re):")
    def ust(ad, u):
        return B["adalar"][ad]["ust"][u]
    b1 = ONK["hipotezler"]["H_195B1"]
    izin, yas = [], []
    olum_B1 = False
    erisilemedi_B1 = []
    for ad, u in b1["izinli_birincil"]:
        v = ust(ad, u)
        ref = ISR[u]
        if not ref["uyum"]:
            erisilemedi_B1.append((ad, u))
        ayni = int(np.sign(v["re"])) == ref["isaret_son"]
        ok = ayni and abs(v["z_re"]) >= 3
        ol = abs(v["z_re"]) < 2 or ((not ayni) and abs(v["z_re"]) >= 3)
        olum_B1 |= ol
        izin.append({"ada": ad, "uydu": u, "re": v["re"], "se": v["se_re"], "z": v["z_re"],
                     "ref_isaret": ref["isaret_son"], "ayni_isaret": ayni, "gecti": ok,
                     "olum_kosulu": ol})
        log(f"  izinli {ad:>6} {u:>6}: Re K̃={v['re']:+.5f}±{v['se_re']:.5f} z={v['z_re']:+.2f} "
            f"(ζ işaret {ref['isaret_son']:+d}) → {'✓' if ok else '✗'}")
    for ad, u in b1["yasak_birincil"]:
        v = ust(ad, u)
        ok = abs(v["z_re"]) < 2
        yas.append({"ada": ad, "uydu": u, "re": v["re"], "se": v["se_re"], "z": v["z_re"],
                    "gecti": ok, "isaret_zeta_ile": "ZIT" if np.sign(v["re"]) != ISR[u]["isaret_son"] else "aynı"})
        log(f"  yasak  {ad:>6} {u:>6}: Re K̃={v['re']:+.5f}±{v['se_re']:.5f} z={v['z_re']:+.2f} "
            f"(|z|<2?) → {'✓' if ok else '✗'}")
    if erisilemedi_B1:
        hB1 = "erişilemedi"
    elif olum_B1:
        hB1 = "ÖLDÜ"
    elif all(x["gecti"] for x in izin) and all(x["gecti"] for x in yas):
        hB1 = "MÜHÜR"
    else:
        hB1 = "KAYIT"
    H["H_195B1"] = {"izinli_birincil": izin, "yasak_birincil": yas, "hukum": hB1}
    log(f"  H-195B1: {hB1}")

    # ======================= H-195B2 =======================
    log("\nH-195B2 — PARİTE:")
    b2 = ONK["hipotezler"]["H_195B2"]
    birincil = b2["birincil"]
    satir = []
    for ad, us in birincil.items():
        for u in us:
            v = ust(ad, u)
            satir.append({"ada": ad, "uydu": u, "parite": "tek" if ad in b2["tek"] else "çift",
                          "re": v["re"], "se": v["se_re"], "z": v["z_re"], "im": v["im"],
                          "z_im": v["z_im"], "aci": v["aci"],
                          "ref_ayni": int(np.sign(v["re"])) == ISR[u]["isaret_son"]})
    i_ok = all(x["ref_ayni"] and abs(x["z"]) >= 2 for x in satir)
    r8 = ust("chi8o", "+log3")["re"] / ust("chi8e", "+log3")["re"]
    # jk se (bağımsız bloklar varsayımı; bilgi)
    e8o, e8e = ust("chi8o", "+log3"), ust("chi8e", "+log3")
    se_r8 = abs(r8) * np.sqrt((e8o["se_re"] / e8o["re"]) ** 2 + (e8e["se_re"] / e8e["re"]) ** 2)
    ii_ok = 0.60 <= r8 <= 1.40
    tek = [x for x in satir if x["parite"] == "tek"]
    cift = [x for x in satir if x["parite"] == "çift"]
    olum_R2 = all(abs(x["z"]) < 2 for x in tek) and all(x["ref_ayni"] and abs(x["z"]) >= 3 for x in cift)
    hB2 = "ÖLDÜ" if olum_R2 else ("MÜHÜR" if (i_ok and ii_ok) else "KAYIT")
    for x in satir:
        log(f"  {x['parite']:>4} {x['ada']:>6} {x['uydu']:>6}: Re={x['re']:+.5f}±{x['se']:.5f} "
            f"z={x['z']:+.2f}  Im z={x['z_im']:+.2f}  ∠{x['aci']:+.1f}°")
    log(f"  (i) işaret: {i_ok}   (ii) ρ₈ = {r8:.4f}±{se_r8:.4f} ∈ [0.60,1.40]: {ii_ok}   "
        f"R2 ölüm deseni: {olum_R2} → H-195B2: {hB2}")
    H["H_195B2"] = {"birinciller": satir, "i_isaret": i_ok, "rho8": r8, "se_rho8": se_r8,
                    "ii_buyukluk": ii_ok, "olum_R2": olum_R2, "hukum": hB2}

    # ======================= H-195B3 =======================
    log("\nH-195B3 — L-EŞLİ BANT ÖLÇEĞİ (R = Re K̃_χ,bant / Re K̃_ζ,bant):")
    ciftler = []
    olcul = []
    for ad, us in birincil.items():
        for u in us:
            v = B["adalar"][ad]["bant"][u]
            zr = ZB[ad]["K"][u]
            R = v["re"] / zr["re"]
            seR = abs(R) * np.sqrt((v["se_re"] / v["re"]) ** 2 + (zr["se_re"] / zr["re"]) ** 2)
            pred = notlar[ad]["kalem_sqrtk_phik"]
            lo, hi = 0.65 * pred, 1.35 * pred
            olc = abs(v["z_re"]) >= 3 and abs(zr["z_re"]) >= 3
            icinde = lo <= R <= hi
            uzak = 0.0 if icinde else ((R - hi) / seR if R > hi else (lo - R) / seR)
            c = {"ada": ad, "uydu": u, "K_ada_bant": v["re"], "se_ada": v["se_re"], "z_ada": v["z_re"],
                 "K_zeta_bant": zr["re"], "se_zeta": zr["se_re"], "z_zeta": zr["z_re"],
                 "R": R, "se_R": seR, "ongoru": pred, "bant": [lo, hi], "olculebilir": olc,
                 "icinde": icinde, "banda_uzaklik_se": uzak,
                 "KAYIT_R_bolu_alt1_k_phik": R / notlar[ad]["alt1_k_phik"],
                 "KAYIT_R_bolu_alt2": R / notlar[ad]["alt2_k_phik_bolu_PiJ0"]}
            ciftler.append(c)
            if olc:
                olcul.append(c)
            log(f"  {ad:>6} {u:>6}: K_χ={v['re']:+.5f}±{v['se_re']:.5f}  K_ζ={zr['re']:+.5f}±{zr['se_re']:.5f}"
                f"  R={R:.3f}±{seR:.3f}  öngörü {pred:.3f} bant [{lo:.3f},{hi:.3f}] "
                f"{'İÇİNDE' if icinde else f'DIŞINDA ({uzak:.1f} se)'}  | R/alt1={R/notlar[ad]['alt1_k_phik']:.3f} "
                f"R/alt2={R/notlar[ad]['alt2_k_phik_bolu_PiJ0']:.3f}")
    if not olcul:
        hB3 = "erişilemedi"
    elif all(c["icinde"] for c in olcul):
        hB3 = "MÜHÜR"
    elif any((not c["icinde"]) and c["banda_uzaklik_se"] >= 2 for c in olcul):
        hB3 = "ÖLDÜ"
    else:
        hB3 = "KAYIT"
    H["H_195B3"] = {"ciftler": ciftler, "hukum": hB3}
    log(f"  H-195B3: {hB3}")

    # ======================= rakip resimler =======================
    R1 = all(abs(x["z"]) < 2 for x in satir)
    R2 = olum_R2 or all(abs(x["z_im"]) >= 3 and abs(x["z"]) < 2 for x in tek)
    R3_isaret = all(x["ref_ayni"] and abs(x["z"]) >= 3 for x in satir)
    H["rakip_resimler"] = {
        "R1_iletkensiz_deseni": R1, "R2_paritesiz_deseni": R2,
        "R3_isaret_parite_deseni": R3_isaret, "R3_olcek_B3": hB3,
        "kazanan": ("R3 (işaret/parite/seçim); ölçek √k/φ(k) " + hB3) if (R3_isaret and not R1 and not R2)
        else ("R1" if R1 else ("R2" if R2 else "belirsiz"))}
    log(f"\nRAKİP RESİMLER: R1 {R1}  R2 {R2}  R3(işaret/parite) {R3_isaret} → {H['rakip_resimler']['kazanan']}")

    # ======================= KAYIT: tam tablo + tüm-izinli bant oranları =======================
    tam = {}
    for ad in o.ADALAR:
        k = B["adalar"][ad]["k"]
        tam[ad] = {}
        for u in o.UYDU_SIRA:
            v = B["adalar"][ad]["ust"][u]
            w = B["adalar"][ad]["bant"][u]
            zr = ZB[ad]["K"][u]
            tam[ad][u] = {"yasak": bool(o.yasak_mi(u, k)), "ust_re": v["re"], "ust_se": v["se_re"],
                          "ust_z": v["z_re"], "ust_im_z": v["z_im"], "aci": v["aci"],
                          "bant_re": w["re"], "bant_se": w["se_re"], "bant_z": w["z_re"],
                          "zeta_bant_re": zr["re"], "zeta_bant_se": zr["se_re"],
                          "R_bant": w["re"] / zr["re"],
                          "se_R_bant": abs(w["re"] / zr["re"]) * np.sqrt(
                              (w["se_re"] / w["re"]) ** 2 + (zr["se_re"] / zr["re"]) ** 2),
                          "B_dislanan": u == "-log5"}
    H["KAYIT_tam_tablo"] = tam
    # yasak/boş konum arka planı (ortalama)
    arka = [tam[ad][u]["ust_re"] for ad in o.ADALAR for u in o.UYDU_SIRA
            if tam[ad][u]["yasak"] and u != "-log5"]
    H["KAYIT_yasak_arka_plan"] = {"ortalama": float(np.mean(arka)), "std": float(np.std(arka)),
                                  "n": len(arka)}
    log(f"\nKAYIT yasak konum arka planı: Re K̃ ort {np.mean(arka):+.5f} (std {np.std(arka):.5f}, n={len(arka)})")
    json.dump(H, open(o.S195 / "HUKUM_195.json", "w"), indent=1, ensure_ascii=False,
              default=lambda x: bool(x) if isinstance(x, np.bool_) else float(x))
    log("-> HUKUM_195.json")
