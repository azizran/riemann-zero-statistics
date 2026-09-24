# -*- coding: utf-8 -*-
"""
193c — K2: HÜKÜM (H-193a/b; ONKAYIT_193 AYNEN, kurtarma yok)
=============================================================
K1_sinif.json'daki blok-yerel s_r = κ_r/κ_top (κ_top BÖLÜNEN HARİÇ) ölçümlerini
ön-kayıtlı eşiklerle karşılaştırır:
 H-193a (+log10, birincil, KÖR): r∈{3,7} negatif, r∈{1,9} pozitif, her biri
   ≥2σ; nicel ±0.25. ÖLÜM: r∈{3,7}'den biri ≥2σ pozitif YA DA r∈{1,9}'dan
   biri ≥2σ negatif.
 H-193b (+log7, ikincil, KÖR): r1,r6 negatif; r3,r4 pozitif, ≥2σ; nicel ±0.30.
   ÖLÜM: r1 ya da r6 ≥2σ pozitif. κ_top < 4·se_kappa_top ⇒ 'erişilemedi'.
 KAYIT: +log5 (doğrulayıcı), +log3/−log3 kontrolleri (½/½ ± 0.10).
Çıktı: 193/HUKUM_193.json + ekran.
"""
import hashlib
import json
from pathlib import Path

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S193 = SCR / "193"


def onkayit():
    o = json.load(open(S193 / "ONKAYIT_193.json"))
    s = hashlib.sha256((QM / "193_configs" / "193a_onkayit.py").read_bytes()).hexdigest()
    if s != o["sha256"]:
        raise SystemExit(f"ON-KAYIT SHA UYUMSUZ: {s} != {o['sha256']}")
    return o


if __name__ == "__main__":
    ONK = onkayit()
    K1 = json.load(open(S193 / "K1_sinif.json"))
    assert K1["sha_onkayit"] == ONK["sha256"]
    HED = ONK["hedefler"]
    print("=" * 78)
    print(f"193c / K2 HÜKÜM  [on-kayit {ONK['zaman']} sha {ONK['sha256'][:12]}]")
    print("=" * 78)
    out = {"sha_onkayit": ONK["sha256"], "hedefler": {}}

    # ======================= H-193a (+log10) =======================
    h10 = K1["hedefler"]["+log10"]
    S = h10["siniflar"]
    tol_a = HED["+log10"]["nicel_tol"]
    neg_set, pos_set = (3, 7), (1, 9)
    yon_ok = {}
    nicel_ok = {}
    olum_evren = []
    for r in neg_set:
        c = S[f"r{r}"]
        yon_ok[r] = c["z"] <= -2.0
        nicel_ok[r] = abs(c["s"] - c["s_pred"]) <= tol_a
        if c["z"] >= 2.0:
            olum_evren.append(f"r{r} pozitife anlamlı (z={c['z']:+.2f})")
    for r in pos_set:
        c = S[f"r{r}"]
        yon_ok[r] = c["z"] >= 2.0
        nicel_ok[r] = abs(c["s"] - c["s_pred"]) <= tol_a
        if c["z"] <= -2.0:
            olum_evren.append(f"r{r} negatife anlamlı (z={c['z']:+.2f})")
    hepsi_yon = all(yon_ok.values())
    hepsi_nicel = all(nicel_ok.values())
    if olum_evren:
        hA = "ÖLDÜ"
    elif hepsi_yon and hepsi_nicel:
        hA = "MÜHÜR"
    else:
        eksik = [f"r{r}" for r in (*neg_set, *pos_set)
                 if not (yon_ok[r] and nicel_ok[r])]
        hA = "KAYIT (kısmi: " + ", ".join(eksik) + " eşiği tam sağlamadı)"
    out["hedefler"]["+log10"] = {
        "hukum": hA, "kappa_top": h10["kappa_top"], "se_kappa_top": h10["se_kappa_top"],
        "siniflar": {f"r{r}": {"s": S[f"r{r}"]["s"], "se_s": S[f"r{r}"]["se_s"],
                               "z": S[f"r{r}"]["z"], "s_pred": S[f"r{r}"]["s_pred"],
                               "yon_ok": bool(yon_ok[r]), "nicel_ok": bool(nicel_ok[r])}
                    for r in (*neg_set, *pos_set)},
        "bolunen": h10["bolunen"], "tau_aralik": h10["tau_aralik"],
        "olum_evren": olum_evren}
    print(f"\nH-193a (+log10, birincil): κ_top={h10['kappa_top']:.5f}"
          f"±{h10['se_kappa_top']:.5f}")
    for r in (1, 3, 7, 9):
        c = S[f"r{r}"]
        print(f"   r{r}: s={c['s']:+.4f}±{c['se_s']:.4f}  z={c['z']:+.2f}  "
              f"öngörü={c['s_pred']:+.4f}  yön {'OK' if yon_ok[r] else 'HAYIR'}  "
              f"nicel {'OK' if nicel_ok[r] else 'HAYIR'}")
    print(f"   ölüm evreni: {olum_evren if olum_evren else 'yok'}  -> {hA}")

    # ======================= H-193b (+log7) =======================
    h7 = K1["hedefler"]["+log7"]
    S7 = h7["siniflar"]
    tol_b = HED["+log7"]["nicel_tol"]
    erisilemedi = h7["kappa_top"] < 4 * h7["se_kappa_top"]
    neg7, pos7 = (1, 6), (3, 4)
    yon7 = {}
    nicel7 = {}
    olum7 = []
    for r in neg7:
        c = S7[f"r{r}"]
        yon7[r] = c["z"] <= -2.0
        nicel7[r] = abs(c["s"] - c["s_pred"]) <= tol_b
        if c["z"] >= 2.0:
            olum7.append(f"r{r} pozitife anlamlı (z={c['z']:+.2f})")
    for r in pos7:
        c = S7[f"r{r}"]
        yon7[r] = c["z"] >= 2.0
        nicel7[r] = abs(c["s"] - c["s_pred"]) <= tol_b
    if erisilemedi:
        hB = "erişilemedi"
    elif olum7:
        hB = "ÖLDÜ"
    elif all(yon7.values()) and all(nicel7.values()):
        hB = "MÜHÜR"
    else:
        eksik = [f"r{r}" for r in (*neg7, *pos7) if not (yon7[r] and nicel7[r])]
        hB = "KAYIT (kısmi: " + ", ".join(eksik) + " eşiği tam sağlamadı)"
    yan_kayit = {f"r{r}": {"s": S7[f"r{r}"]["s"], "se_s": S7[f"r{r}"]["se_s"],
                           "z": S7[f"r{r}"]["z"], "s_pred": S7[f"r{r}"]["s_pred"]}
                for r in (2, 5)}
    out["hedefler"]["+log7"] = {
        "hukum": hB, "kappa_top": h7["kappa_top"], "se_kappa_top": h7["se_kappa_top"],
        "erisilemedi_kosulu": bool(erisilemedi),
        "siniflar": {f"r{r}": {"s": S7[f"r{r}"]["s"], "se_s": S7[f"r{r}"]["se_s"],
                               "z": S7[f"r{r}"]["z"], "s_pred": S7[f"r{r}"]["s_pred"],
                               "yon_ok": bool(yon7[r]), "nicel_ok": bool(nicel7[r])}
                    for r in (*neg7, *pos7)},
        "yan_KAYIT_r2_r5": yan_kayit, "bolunen": h7["bolunen"],
        "tau_aralik": h7["tau_aralik"], "olum_evren": olum7}
    print(f"\nH-193b (+log7, ikincil): κ_top={h7['kappa_top']:.5f}±{h7['se_kappa_top']:.5f}"
          f"  (4σ eşik={4*h7['se_kappa_top']:.5f}, erişilemedi={erisilemedi})")
    for r in (1, 6, 3, 4):
        c = S7[f"r{r}"]
        print(f"   r{r}: s={c['s']:+.4f}±{c['se_s']:.4f}  z={c['z']:+.2f}  "
              f"öngörü={c['s_pred']:+.4f}  yön {'OK' if yon7[r] else 'HAYIR'}  "
              f"nicel {'OK' if nicel7[r] else 'HAYIR'}")
    for r in (2, 5):
        c = S7[f"r{r}"]
        print(f"   r{r} (yan KAYIT): s={c['s']:+.4f}±{c['se_s']:.4f}  z={c['z']:+.2f}  "
              f"öngörü={c['s_pred']:+.4f}")
    print(f"   ölüm evreni: {olum7 if olum7 else 'yok'}  -> {hB}")

    # ======================= KAYIT: +log5, +log3, -log3 =======================
    kayit = {}
    h5 = K1["hedefler"]["+log5"]
    S5 = h5["siniflar"]
    kayit["+log5"] = {
        "kappa_top": h5["kappa_top"], "se_kappa_top": h5["se_kappa_top"],
        "siniflar": {f"r{r}": {"s": S5[f"r{r}"]["s"], "se_s": S5[f"r{r}"]["se_s"],
                               "z": S5[f"r{r}"]["z"], "s_pred": S5[f"r{r}"]["s_pred"]}
                    for r in (1, 2, 3, 4)},
        "bolunen": h5["bolunen"], "tau_aralik": h5["tau_aralik"]}
    print(f"\nKAYIT +log5 (doğrulayıcı): κ_top={h5['kappa_top']:.5f}±{h5['se_kappa_top']:.5f}")
    for r in (1, 2, 3, 4):
        c = S5[f"r{r}"]
        print(f"   r{r}: s={c['s']:+.4f}±{c['se_s']:.4f}  z={c['z']:+.2f}  "
              f"öngörü={c['s_pred']:+.4f}")

    for ad in ("+log3", "-log3"):
        h = K1["hedefler"][ad]
        S_ = h["siniflar"]
        bantta = {f"r{r}": bool(abs(S_[f"r{r}"]["s"] - 0.5) <= 0.10) for r in (1, 2)}
        kayit[ad] = {
            "kappa_top": h["kappa_top"], "se_kappa_top": h["se_kappa_top"],
            "siniflar": {f"r{r}": {"s": S_[f"r{r}"]["s"], "se_s": S_[f"r{r}"]["se_s"],
                                   "bantta_0.5pm0.1": bantta[f"r{r}"]}
                        for r in (1, 2)},
            "hepsi_bantta": all(bantta.values()), "bolunen": h["bolunen"],
            "tau_aralik": h["tau_aralik"]}
        print(f"\nKAYIT {ad} (kontrol, ½/½±0.10): κ_top={h['kappa_top']:.5f}"
              f"±{h['se_kappa_top']:.5f}")
        for r in (1, 2):
            c = S_[f"r{r}"]
            print(f"   r{r}: s={c['s']:+.4f}±{c['se_s']:.4f}  "
                  f"bantta {'EVET' if bantta[f'r{r}'] else 'HAYIR'}")

    out["KAYIT"] = kayit
    out["hukum"] = {"H-193a": out["hedefler"]["+log10"]["hukum"],
                    "H-193b": out["hedefler"]["+log7"]["hukum"]}
    out["tam_kontrol_plog10_r1"] = K1["tam_kontrol_plog10_r1"]
    print(f"\nTAM-ÖRNEKLEM KONTROLÜ +log10 r1: {K1['tam_kontrol_plog10_r1']}")

    print("\n" + "=" * 78)
    print(f"HÜKÜMLER: H-193a = {out['hukum']['H-193a']}  |  H-193b = {out['hukum']['H-193b']}")
    print("=" * 78)
    json.dump(out, open(S193 / "HUKUM_193.json", "w"), indent=1, ensure_ascii=False)
    print(f"-> {S193/'HUKUM_193.json'}  BİTTİ")
