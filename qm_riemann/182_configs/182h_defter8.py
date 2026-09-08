# -*- coding: utf-8 -*-
"""
182h — K1(c)/K3: ÇAPASIZ NEDENSEL DEFTERİN n_VF = 8 HÂLİ (O1 AYNEN)
====================================================================
ÖN-MÜHÜR: bütün formüller 180a'da, okunurluk ölçütü O1 (k = 2.0) 181a'da,
n=8 öngörüleri (P2) 182a'da — hepsi ölçümden ÖNCE donduruldu.

MAKİNE İTHAL EDİLİR, KOPYALANMAZ: `181c_defter_n4.ayrisim` (o da
`180f_defter.istat/HAM/ORAN/AD`'ı ithal eder).

ÜREME KAPISI: VS = VS1..VS4, VF = VF1..VF4 ile koşulduğunda
`181/K3_DEFTER_n4.json`'un HER alanını TAM SIFIR farkla üretmeli.
Üretmezse betik ÖLÜR.

Kullanım: 182h_defter8.py
"""
import importlib
import json
import math
import sys
import time
from pathlib import Path

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
for _p in ("181_configs", "180_configs"):
    sys.path.insert(0, str(QM / _p))
C181 = importlib.import_module("181c_defter_n4")
F180 = importlib.import_module("180f_defter")
ayrisim, HAM, ORAN, AD = C181.ayrisim, F180.HAM, F180.ORAN, F180.AD

SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S182, S181, S180 = SCR / "182", SCR / "181", SCR / "180"


def main():
    t0 = time.time()
    OK180 = json.load(open(S180 / "ONKAYIT_K0.json"))
    OK182 = json.load(open(S182 / "ONKAYIT_182.json"))
    D4 = json.load(open(S181 / "K3_DEFTER_n4.json"))
    KS = OK182["capalar"]["k_sigma"]
    capa = OK180["capa_defter"]
    son, Hk = capa["son"], capa["Hkeskin"]
    VS = [json.load(open(S180 / f"G_{g}.json"))
          for g in ("VS1", "VS2", "VS3", "VS4")]
    VF4 = [capa[f"VF{i}"] for i in (1, 2, 3, 4)]
    VF8 = VF4 + [json.load(open(S182 / f"G_VF{i}.json")) for i in (5, 6, 7, 8)]

    print("=" * 78, flush=True)
    print("182h — ÇAPASIZ NEDENSEL DEFTER, n_VS = 4, n_VF = 8 (O1 AYNEN)",
          flush=True)
    print(f"    182 ön-kayıt sha = {OK182['sha256'][:16]}…   O1 k = {KS}",
          flush=True)
    print("=" * 78, flush=True)

    # ------------------------------------------------------- R-KAPISI
    R4 = ayrisim(son, Hk, VS, VF4)
    kotu, ncmp = [], 0
    for konv in ("dogrusal", "log"):
        for k in HAM + ORAN:
            for alan in ("delta", "ZARF", "KILIT", "se_ZARF", "p_zarf",
                         "p_kilit", "VS_ort", "VF_ort", "VS_se", "VF_se",
                         "VS_sd", "VF_sd", "se_p"):
                ncmp += 1
                if R4[konv][k][alan] != D4[konv][k][alan]:
                    kotu.append((konv, k, alan, R4[konv][k][alan],
                                 D4[konv][k][alan]))
    for k in HAM + ORAN:
        for alan in ("logK_VF", "logK_VS", "fark", "KILITlog"):
            ncmp += 1
            if R4["K4"][k][alan] != D4["K4"][k][alan]:
                kotu.append(("K4", k, alan))
    print(f"  R-KAPISI (VF1..VF4 ile 181/K3_DEFTER_n4.json'u üretme): "
          f"{ncmp} alan, {len(kotu)} fark", flush=True)
    if kotu:
        for c in kotu[:12]:
            print("   ✗", c, flush=True)
        raise SystemExit("R-KAPISI TUTMADI — n=8 defteri basılmaz.")
    print("  R-KAPISI ✓ TAM SIFIR FARK", flush=True)

    # -------------------------------------------------------- n_VF = 8
    R8 = ayrisim(son, Hk, VS, VF8)

    print("\n  T2 — DEFTER (HAM esas | oran EK), VF ailesi n = 8")
    print("  gaz        " + "".join(f"{AD[k]:>14s}" for k in HAM)
          + " | " + "".join(f"{AD[k]:>10s}" for k in ORAN))
    for nm, d in ([("son", son), ("Hkeskin", Hk)]
                  + [(f"VF{i+1}", VF8[i]) for i in range(8)]
                  + [(f"VS{i+1}", VS[i]) for i in range(4)]):
        print(f"  {nm:10s} " + "".join(f"{d[k]:14.6f}" for k in HAM)
              + " | " + "".join(f"{d[k]:10.6f}" for k in ORAN), flush=True)

    for etiket, anah in (("DOĞRUSAL", "dogrusal"), ("LOG", "log")):
        print("\n" + "=" * 78)
        print(f"  T3/T4 — AYRIŞIM ({etiket}), n_VS = 4, n_VF = 8")
        print("=" * 78)
        print("  nicelik    Δ          ⟨⟩_VS      ⟨⟩_VF       ZARF±se"
              "            KİLİT±se          p_zarf±se       |Z|/se  hüküm")
        for k in HAM + ORAN:
            r = R8[anah][k]
            t = abs(r["ZARF"]) / r["se_ZARF"]
            r["t"], r["hukum"] = t, ("OKUNUR" if t >= KS else "HÜKÜMSÜZ")
            print(f"  {AD[k]:10s} {r['delta']:+10.6f} {r['VS_ort']:10.6f} "
                  f"{r['VF_ort']:10.6f}  {r['ZARF']:+9.6f}±{r['se_ZARF']:.6f}"
                  f"  {r['KILIT']:+9.6f}±{r['se_ZARF']:.6f} "
                  f"{r['p_zarf']:+8.4f}±{r['se_p']:.4f} {t:8.3f}  "
                  f"{r['hukum']}", flush=True)

    # -------------------- P2 öngörülerinin denetimi -------------------
    print("\n" + "=" * 78)
    print("  T6 — ÖN-KAYITLI ÖNGÖRÜ P2'nin DENETİMİ (182a, veriden önce)")
    print("=" * 78)
    print("  satır      konv.     |Z|/se(4)  |Z|/se(8)ö  |Z|/se(8)  "
          "sd_VF(4)  sd_VF(8)  hüküm(4)  hüküm(8)  ÖNGÖRÜ   ")
    DEG, sapma = {}, 0
    for k in HAM + ORAN:
        for konv in ("dogrusal", "log"):
            a, b = D4[konv][k], R8[konv][k]
            o = OK182["P2_ongoru"][f"{konv}/{k}"]
            t4 = abs(a["ZARF"]) / a["se_ZARF"]
            h4 = "OKUNUR" if t4 >= KS else "HÜKÜMSÜZ"
            tut = b["hukum"] == o["ongorulen_hukum"]
            sapma += 0 if tut else 1
            print(f"  {AD[k]:10s} {konv:9s} {t4:9.3f} {o['t_n8_ongoru']:11.3f} "
                  f"{b['t']:10.3f}  {a['VF_sd']:8.6f} {b['VF_sd']:8.6f}  "
                  f"{h4:9s} {b['hukum']:9s} {o['ongorulen_hukum']}"
                  f"{'' if tut else '  ← ÖNGÖRÜ TUTMADI'}", flush=True)
            DEG[f"{konv}/{k}"] = dict(
                t_n4=t4, t_n8=b["t"], t_n8_ongoru=o["t_n8_ongoru"],
                hukum_n4=h4, hukum_n8=b["hukum"],
                ongorulen_hukum=o["ongorulen_hukum"], ongoru_tuttu=tut,
                ZARF_n4=a["ZARF"], ZARF_n8=b["ZARF"],
                d_merkez=b["p_zarf"] - a["p_zarf"],
                VF_sd_n4=a["VF_sd"], VF_sd_n8=b["VF_sd"],
                VF_se_n4=a["VF_se"], VF_se_n8=b["VF_se"],
                se_ZARF_n4=a["se_ZARF"], se_ZARF_n8=b["se_ZARF"],
                se_p_n4=a["se_p"], se_p_n8=b["se_p"],
                cubuk_orani=b["se_ZARF"] / a["se_ZARF"],
                ulasilabilir_VF_ile=o["ulasilabilir_VF_ile"])
    yeni = [k for k, v in DEG.items()
            if v["hukum_n4"] == "HÜKÜMSÜZ" and v["hukum_n8"] == "OKUNUR"]
    kayip = [k for k, v in DEG.items()
             if v["hukum_n4"] == "OKUNUR" and v["hukum_n8"] == "HÜKÜMSÜZ"]
    print(f"\n  n=4 → n=8 AÇILAN satırlar: {yeni if yeni else 'HİÇBİRİ'}")
    print(f"  n=4 → n=8 KAPANAN satırlar: {kayip if kayip else 'HİÇBİRİ'}")
    print(f"  P2 öngörüsünden sapma: {sapma}/18", flush=True)

    # ----------------------------------------------------------- K4
    print("\n" + "=" * 78)
    print("  T7 — K4: KARIŞTIRMA ÇÖKÜŞÜ (n_VF=8) ile ÇAPRAZ TUTARLILIK")
    print("=" * 78)
    print("  nicelik    K_N(%) n=8   K_N(%) n=4   K'_N(%)     logK−logK'"
          "   KİLİTlog     kalıntı")
    for k in HAM + ORAN:
        r, r4 = R8["K4"][k], D4["K4"][k]
        print(f"  {AD[k]:10s} {r['cokus_VF_yuzde']:+11.2f} "
              f"{r4['cokus_VF_yuzde']:+12.2f} {r['cokus_VS_yuzde']:+11.2f} "
              f"{r['fark']:+13.6f} {r['KILITlog']:+11.6f}  "
              f"{r['kalinti']:+.2e}", flush=True)

    rec = dict(n_VS=4, n_VF=8, k_sigma=KS, dogrusal=R8["dogrusal"],
               log=R8["log"], K4=R8["K4"], degisim=DEG, acilan=yeni,
               kapanan=kayip, P2_sapma=sapma, R_kapisi_fark=len(kotu),
               onkayit_182_sha=OK182["sha256"], zaman=time.ctime(),
               sure_s=time.time() - t0)
    p = S182 / "K3_DEFTER_n8.json"
    p.write_text(json.dumps(rec, indent=1, ensure_ascii=False, default=float))
    print(f"\n-> {p}   ({time.time()-t0:.1f} s)", flush=True)


if __name__ == "__main__":
    main()
