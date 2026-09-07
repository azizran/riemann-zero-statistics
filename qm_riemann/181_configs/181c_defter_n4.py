# -*- coding: utf-8 -*-
"""
181c — K2: ÇAPASIZ NEDENSEL DEFTERİN n=4 HÂLİ
==============================================
ÖN-MÜHÜR: bütün formüller 180a'da (ve okunurluk ölçütü O1, 181a'da)
ölçümden ÖNCE donduruldu. Bu betik hiçbir seçim yapmaz, iki
konvansiyonu da basar.

MAKİNE: kestirici `istat` ve nicelik listeleri (HAM/ORAN/AD)
`180_configs/180f_defter.py`'den İTHAL EDİLİR — kopyalanmaz. 180f'in
`main()`'i VS listesini ("VS1","VS2") gövdesinde sabitlediği için n=4
döngüsü burada kuruludur; buna karşılık **ÜREME KAPISI** koşulur:

  R-KAPISI: aynı döngü VS=(VS1,VS2) ile koşulduğunda 180'in
  `K3_DEFTER.json`'undaki HER alan (delta, ZARF, KİLİT, se, paylar,
  K4) ile **TAM SIFIR** farkla örtüşmeli. Örtüşmezse betik ÖLÜR.
  Yani n=4 tablosu, 180'in makinesinin bit-bit aynısıyla basılmıştır.

Kullanım: 181c_defter_n4.py
"""
import importlib
import json
import math
import sys
import time
from pathlib import Path

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
sys.path.insert(0, str(QM / "180_configs"))
F180 = importlib.import_module("180f_defter")
istat, HAM, ORAN, AD = F180.istat, F180.HAM, F180.ORAN, F180.AD

SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S180, S181 = SCR / "180", SCR / "181"


def ayrisim(son, Hk, VS, VF):
    """180a'nın dondurduğu iki konvansiyon; 180f ile birebir aynı satırlar."""
    OUT = {"dogrusal": {}, "log": {}, "K4": {}}
    for LOG in (False, True):
        anah = "log" if LOG else "dogrusal"
        for k in HAM + ORAN:
            f = (lambda x: math.log(x)) if LOG else (lambda x: x)
            ns, nh = f(son[k]), f(Hk[k])
            vs = [f(d[k]) for d in VS]
            vf = [f(d[k]) for d in VF]
            D = ns - nh
            ms, sds, ses, ns_n = istat(vs)
            mf, sdf, sef, _ = istat(vf)
            Z = ms - mf
            KI = D - Z
            seZ = math.sqrt(ses ** 2 + sef ** 2)
            OUT[anah][k] = dict(
                son=ns, Hk=nh, delta=D, VS_ort=ms, VS_sd=sds, VS_se=ses,
                VF_ort=mf, VF_sd=sdf, VF_se=sef, ZARF=Z, KILIT=KI,
                se_ZARF=seZ, p_zarf=Z / D, p_kilit=KI / D, se_p=seZ / abs(D),
                VS_tek=vs, VF_tek=vf, n_VS=ns_n)
    for k in HAM + ORAN:
        r = OUT["log"][k]
        lK = r["VF_ort"] - r["Hk"]
        lKp = r["VS_ort"] - r["son"]
        OUT["K4"][k] = dict(logK_VF=lK, logK_VS=lKp,
                            cokus_VF_yuzde=100 * (math.exp(lK) - 1),
                            cokus_VS_yuzde=100 * (math.exp(lKp) - 1),
                            fark=lK - lKp, KILITlog=r["KILIT"],
                            kalinti=(lK - lKp) - r["KILIT"])
    return OUT


def main():
    t0 = time.time()
    OK180 = json.load(open(S180 / "ONKAYIT_K0.json"))
    OK181 = json.load(open(S181 / "ONKAYIT_181_K0.json"))
    K3_180 = json.load(open(S180 / "K3_DEFTER.json"))
    KS = OK181["okunurluk"]["k_sigma"]
    capa = OK180["capa_defter"]
    son, Hk = capa["son"], capa["Hkeskin"]
    VF = [capa[g] for g in ("VF1", "VF2", "VF3", "VF4")]

    G = {}
    for g in ("VS1", "VS2", "VS3", "VS4"):
        p = S180 / f"G_{g}.json"
        if not p.exists():
            raise SystemExit(f"{p} yok — ölçüm bitmeden defter basılmaz.")
        G[g] = json.load(open(p))

    print("=" * 78, flush=True)
    print("181c — K2: ÇAPASIZ NEDENSEL DEFTER, n=4", flush=True)
    print(f"    181 ön-kayıt sha = {OK181['sha256'][:16]}…   "
          f"180 ön-kayıt sha = {OK180['sha256'][:16]}…", flush=True)
    print("=" * 78, flush=True)

    # ------------------------------------------------------ R-KAPISI
    R2 = ayrisim(son, Hk, [G["VS1"], G["VS2"]], VF)
    kotu, ncmp = [], 0
    for konv in ("dogrusal", "log"):
        for k in HAM + ORAN:
            for alan in ("delta", "ZARF", "KILIT", "se_ZARF", "p_zarf",
                         "p_kilit", "VS_ort", "VF_ort", "VS_se", "VF_se",
                         "VS_sd", "VF_sd", "se_p"):
                a, b = R2[konv][k][alan], K3_180[konv][k][alan]
                ncmp += 1
                if a != b:
                    kotu.append((konv, k, alan, a, b))
    for k in HAM + ORAN:
        for alan in ("logK_VF", "logK_VS", "fark", "KILITlog"):
            a, b = R2["K4"][k][alan], K3_180["K4"][k][alan]
            ncmp += 1
            if a != b:
                kotu.append(("K4", k, alan, a, b))
    print(f"  R-KAPISI (n=2 ile 180f'in K3_DEFTER.json'unu üretme): "
          f"{ncmp} alan karşılaştırıldı, {len(kotu)} fark", flush=True)
    if kotu:
        for c in kotu[:12]:
            print("   ✗", c, flush=True)
        raise SystemExit("R-KAPISI TUTMADI — n=4 tablosu basılmaz.")
    print("  R-KAPISI ✓ TAM SIFIR FARK (makine 180f'inkiyle birebir)",
          flush=True)

    # -------------------------------------------------------- n=4 defter
    VS4 = [G[g] for g in ("VS1", "VS2", "VS3", "VS4")]
    R4 = ayrisim(son, Hk, VS4, VF)

    print("\n  T2 — DEFTER (HAM esas | oran EK), n=4")
    print("  gaz        " + "".join(f"{AD[k]:>14s}" for k in HAM)
          + " | " + "".join(f"{AD[k]:>10s}" for k in ORAN))
    for nm, d in ([("son", son), ("Hkeskin", Hk)]
                  + [(f"VF{i+1}", VF[i]) for i in range(4)]
                  + [(g, G[g]) for g in ("VS1", "VS2", "VS3", "VS4")]):
        print(f"  {nm:10s} " + "".join(f"{d[k]:14.6f}" for k in HAM)
              + " | " + "".join(f"{d[k]:10.6f}" for k in ORAN), flush=True)

    for etiket, anah in (("DOĞRUSAL", "dogrusal"), ("LOG", "log")):
        print("\n" + "=" * 78)
        print(f"  T3/T4 — AYRIŞIM ({etiket} konvansiyon), n_VS = 4")
        print("=" * 78)
        print("  nicelik        N(son)     N(Hk)        Δ        ⟨⟩_VS      "
              "⟨⟩_VF      ZARF±se           KİLİT±se        p_zarf±se    "
              "p_kilit±se   |Z|/se  hüküm")
        for k in HAM + ORAN:
            r = R4[anah][k]
            t = abs(r["ZARF"]) / r["se_ZARF"]
            h = "OKUNUR" if t >= KS else "HÜKÜMSÜZ"
            r["t"] = t
            r["hukum"] = h
            print(f"  {AD[k]:10s} {r['son']:10.6f} {r['Hk']:10.6f} "
                  f"{r['delta']:+10.6f} {r['VS_ort']:10.6f} "
                  f"{r['VF_ort']:10.6f}  {r['ZARF']:+9.6f}±{r['se_ZARF']:.6f}"
                  f"  {r['KILIT']:+9.6f}±{r['se_ZARF']:.6f} "
                  f"{r['p_zarf']:+8.4f}±{r['se_p']:.4f} "
                  f"{r['p_kilit']:+8.4f}±{r['se_p']:.4f} {t:7.3f}  {h}",
                  flush=True)

    # ------------------------------------------------ n=2 → n=4 kıyas
    print("\n" + "=" * 78)
    print("  T8 — n=2 → n=4: MERKEZ Mİ KAYDI, ÇUBUK MU SIKIŞTI?")
    print("  (ön-kayıtlı öngörü 181a/H1: ZARF sabit + sd sabit varsayımı)")
    print("=" * 78)
    print("  nicelik    konv.     p_zarf(2)   p_zarf(4)  Δmerkez |  se_p(2) "
          "se_p(4)ö  se_p(4)   çubuk×  |  hüküm(2)  hüküm(4)   öngörü")
    DEG = {}
    for k in HAM + ORAN:
        for konv in ("dogrusal", "log"):
            a, b = K3_180[konv][k], R4[konv][k]
            o = OK181["ongoruler"][f"{konv}/{k}"]
            h2 = o["simdiki_hukum"]
            print(f"  {AD[k]:10s} {konv:9s} {a['p_zarf']:+10.4f} "
                  f"{b['p_zarf']:+10.4f} {b['p_zarf']-a['p_zarf']:+8.4f} | "
                  f"{a['se_p']:8.4f} {o['se_p_n4_ongoru']:8.4f} "
                  f"{b['se_p']:8.4f} {b['se_p']/a['se_p']:7.3f}  | "
                  f"{h2:9s} {b['hukum']:9s}  {o['ongorulen_hukum']}"
                  f"{'' if b['hukum']==o['ongorulen_hukum'] else '  ← ÖNGÖRÜ TUTMADI'}",
                  flush=True)
            DEG[f"{konv}/{k}"] = dict(
                p_zarf_n2=a["p_zarf"], p_zarf_n4=b["p_zarf"],
                d_merkez=b["p_zarf"] - a["p_zarf"], se_p_n2=a["se_p"],
                se_p_n4_ongoru=o["se_p_n4_ongoru"], se_p_n4=b["se_p"],
                cubuk_orani=b["se_p"] / a["se_p"], hukum_n2=h2,
                hukum_n4=b["hukum"], ongoru=o["ongorulen_hukum"],
                ongoru_tuttu=bool(b["hukum"] == o["ongorulen_hukum"]),
                t_n2=o["t_n2"], t_n4=b["t"], t_n4_ongoru=o["t_n4_ongoru"],
                sd_VS_n2=o["sd_VS_n2"], sd_VS_n4=b["VS_sd"],
                ulasilabilir=o["ulasilabilir"])

    # ---------------------------------------------------------------- K4
    print("\n" + "=" * 78)
    print("  T7 — K4: KARIŞTIRMA ÇÖKÜŞLERİYLE ÇAPRAZ TUTARLILIK (n=4)")
    print("=" * 78)
    print("  nicelik      K_N (%)      K'_N (%)     logK−logK'   KİLİTlog"
          "     kalıntı     [180 n=2: K'_N %]")
    for k in HAM + ORAN:
        r = R4["K4"][k]
        print(f"  {AD[k]:10s} {r['cokus_VF_yuzde']:+9.2f}  "
              f"{r['cokus_VS_yuzde']:+9.2f}   {r['fark']:+11.6f} "
              f"{r['KILITlog']:+11.6f}  {r['kalinti']:+.2e}   "
              f"{K3_180['K4'][k]['cokus_VS_yuzde']:+9.2f}", flush=True)

    # ----------------------------------------- zarf aktarımı sağlaması
    print("\n" + "=" * 78)
    print("  T9 — ZARF AKTARIMININ SAĞLAMASI (tanı, kapı değil), n=4")
    print("=" * 78)
    tg = OK180["zarf"]["tau_eff"]
    rr = OK180["zarf"]["r"]
    import numpy as _np
    S167 = SCR / "167"
    CVF = [json.load(open(S167 / f"C_{g}.json"))
           for g in ("VF1", "VF2", "VF3", "VF4")]
    bVF = [[b for b in d["bant"] if b.get("olculdu")] for d in CVF]
    BANT = []
    print("  τ_eff   ⟨R⟩_VS(n=4)  ⟨R⟩_VF(n=4)   ölçülen oran   tarif r(τ)"
          "   fark")
    for i, b in enumerate(G["VS1"]["bant"]):
        rs = [G[g]["bant"][i]["R_bant"] for g in ("VS1", "VS2", "VS3", "VS4")]
        rf = [x[i]["R_bant"] for x in bVF]
        ms, mf = sum(rs) / len(rs), sum(rf) / len(rf)
        tar = float(_np.interp(b["tau_eff"], tg, rr))
        BANT.append(dict(tau_eff=b["tau_eff"], R_VS_ort=ms, R_VF_ort=mf,
                         oran=ms / mf, r_tarif=tar,
                         fark_yuzde=100 * (ms / mf / tar - 1),
                         R_VS_tek=rs, R_VF_tek=rf))
        print(f"  {b['tau_eff']:.4f}   {ms:9.4f}   {mf:9.4f}   "
              f"{ms/mf:11.4f}   {tar:10.4f}  {100*(ms/mf/tar-1):+6.1f}%",
              flush=True)

    OUT = dict(zaman=time.ctime(), betik="181c_defter_n4.py",
               onkayit_181_sha=OK181["sha256"],
               onkayit_180_sha=OK180["sha256"],
               R_kapisi=dict(alan=ncmp, fark=0, gecti=True),
               n_VS=4, VS=["VS1", "VS2", "VS3", "VS4"],
               k_sigma=KS, dogrusal=R4["dogrusal"], log=R4["log"],
               K4=R4["K4"], degisim=DEG, bant=BANT,
               sure_s=time.time() - t0)
    p = S181 / "K3_DEFTER_n4.json"
    p.write_text(json.dumps(OUT, indent=1, ensure_ascii=False, default=float))
    print(f"\n-> {p}   ({time.time()-t0:.2f} s)", flush=True)


if __name__ == "__main__":
    main()
