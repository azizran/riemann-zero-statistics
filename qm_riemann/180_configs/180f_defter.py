# -*- coding: utf-8 -*-
"""
180f — K3/K4: NİHAİ ÇAPASIZ NEDENSEL DEFTER
============================================
ÖN-MÜHÜR: bütün formüller `180/ONKAYIT_K0.json`'da, ölçümden önce
dondurulmuştu. Bu betik hiçbir seçim yapmaz; iki konvansiyonu da
hesaplar ve ikisini de yazar.

K3:  ΔN(son↔Hk) = ZARF_N + KİLİT_N   (doğrusal ve log), HAM defter esas
     (M, Q_E, ρ_E, Q_X, ρ_X), oran dili (g_E, g_X, θ, g_cal) yalnız EK;
     tohum-saçılımı hata çubuğu; eski çapa-bağıl %39 ile karşılaştırma.
K4:  ΔQ_E/Δρ_E dilinde 176'nın karıştırma çöküşleriyle çapraz tutarlılık.
     ÖZDEŞLİK (türetim, varsayımsız):
        KİLİTlog_N = log K_N − log K'_N ,
        K_N  := exp(⟨log N⟩_VF) / N(Hk)   [Hk zarfında karıştırma çöküşü]
        K'_N := exp(⟨log N⟩_VS) / N(son)  [son zarfında aynı çöküş]
     yani log dilinde KİLİT PAYI, iki zarfta ölçülen karıştırma
     çöküşlerinin FARKIDIR — ayrı bir kestirici değil, aynı sayıdır.

Kullanım: 180f_defter.py
"""
import json
import math
import time
from pathlib import Path

import numpy as np

SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S180 = SCR / "180"
HAM = ("M", "Q_E", "rho_E", "Q_X", "rho_X")
ORAN = ("g_E", "g_X", "theta", "g_cal")
AD = {"M": "M (=KALİB_u2)", "Q_E": "Q_E", "rho_E": "ρ_E", "Q_X": "Q_X",
      "rho_X": "ρ_X", "g_E": "g_E", "g_X": "g_X", "theta": "θ",
      "g_cal": "g_cal"}


def istat(vals):
    v = np.asarray(vals, float)
    n = len(v)
    sd = float(np.std(v, ddof=1)) if n > 1 else 0.0
    return float(v.mean()), sd, sd / math.sqrt(n), n


def main():
    OK = json.load(open(S180 / "ONKAYIT_K0.json"))
    capa = OK["capa_defter"]
    print("=" * 74, flush=True)
    print("180f — K3/K4: NİHAİ ÇAPASIZ NEDENSEL DEFTER", flush=True)
    print(f"    ön-kayıt {OK['zaman']}  sha={OK['sha256'][:16]}…", flush=True)
    print("=" * 74, flush=True)

    VF = [capa[g] for g in ("VF1", "VF2", "VF3", "VF4")]
    VS = []
    for g in ("VS1", "VS2"):
        p = S180 / f"G_{g}.json"
        if not p.exists():
            raise SystemExit(f"{p} yok — K1 bitmeden K3 koşulmaz.")
        VS.append(json.load(open(p)))
    son, Hk = capa["son"], capa["Hkeskin"]

    print("\n  T2 — DEFTER (HAM esas | oran EK)")
    print("  gaz        " + "".join(f"{AD[k]:>14s}" for k in HAM)
          + " | " + "".join(f"{AD[k]:>10s}" for k in ORAN))
    for nm, d in [("son", son), ("Hkeskin", Hk)] + \
            [(f"VF{i+1}", VF[i]) for i in range(4)] + \
            [(VS[i]["ad"], VS[i]) for i in range(len(VS))]:
        print(f"  {nm:10s} " + "".join(f"{d[k]:14.6f}" for k in HAM)
              + " | " + "".join(f"{d[k]:10.6f}" for k in ORAN), flush=True)

    OUT = {"zaman": time.ctime(), "onkayit_sha": OK["sha256"],
           "dogrusal": {}, "log": {}, "K4": {}}

    for etiket, LOG in (("DOĞRUSAL", False), ("LOG", True)):
        print("\n" + "=" * 74)
        print(f"  T3/T4 — AYRIŞIM ({etiket} konvansiyon)")
        print("=" * 74)
        print("  nicelik        N(son)     N(Hk)        Δ        ⟨⟩_VS      "
              "⟨⟩_VF      ZARF±se           KİLİT±se        p_zarf±se    "
              "p_kilit±se")
        for k in HAM + ORAN:
            f = (lambda x: math.log(x)) if LOG else (lambda x: x)
            try:
                ns, nh = f(son[k]), f(Hk[k])
                vs = [f(d[k]) for d in VS]
                vf = [f(d[k]) for d in VF]
            except ValueError:
                print(f"  {AD[k]:10s}  (log tanımsız — atlandı)")
                continue
            D = ns - nh
            ms, sds, ses, _ = istat(vs)
            mf, sdf, sef, _ = istat(vf)
            Z = ms - mf
            KI = D - Z
            seZ = math.sqrt(ses ** 2 + sef ** 2)
            pz, pk = Z / D, KI / D
            sep = seZ / abs(D)
            print(f"  {AD[k]:10s} {ns:10.6f} {nh:10.6f} {D:+10.6f} "
                  f"{ms:10.6f} {mf:10.6f}  {Z:+9.6f}±{seZ:.6f}  "
                  f"{KI:+9.6f}±{seZ:.6f}  {pz:+8.4f}±{sep:.4f} "
                  f"{pk:+8.4f}±{sep:.4f}", flush=True)
            OUT["log" if LOG else "dogrusal"][k] = dict(
                son=ns, Hk=nh, delta=D, VS_ort=ms, VS_sd=sds, VS_se=ses,
                VF_ort=mf, VF_sd=sdf, VF_se=sef, ZARF=Z, KILIT=KI,
                se_ZARF=seZ, p_zarf=pz, p_kilit=pk, se_p=sep,
                VS_tek=vs, VF_tek=vf)

    # ---------------------------------------------------------------- K4
    print("\n" + "=" * 74)
    print("  T7 — K4: KARIŞTIRMA ÇÖKÜŞLERİYLE ÇAPRAZ TUTARLILIK")
    print("  K_N  = exp⟨log N⟩_VF / N(Hk)   [Hk zarfında karıştırma çöküşü]")
    print("  K'_N = exp⟨log N⟩_VS / N(son)  [son zarfında aynı çöküş]")
    print("  ÖZDEŞLİK: KİLİTlog_N = log K_N − log K'_N   (kalıntı yazılır)")
    print("=" * 74)
    print("  nicelik      K_N (%)      K'_N (%)     logK−logK'   KİLİTlog"
          "     kalıntı")
    for k in HAM + ORAN:
        r = OUT["log"].get(k)
        if r is None:
            continue
        lK = r["VF_ort"] - r["Hk"]
        lKp = r["VS_ort"] - r["son"]
        art = (lK - lKp) - r["KILIT"]
        print(f"  {AD[k]:10s} {100*(math.exp(lK)-1):+9.2f}  "
              f"{100*(math.exp(lKp)-1):+9.2f}   {lK-lKp:+11.6f} "
              f"{r['KILIT']:+11.6f}  {art:+.2e}", flush=True)
        OUT["K4"][k] = dict(logK_VF=lK, logK_VS=lKp,
                            cokus_VF_yuzde=100 * (math.exp(lK) - 1),
                            cokus_VS_yuzde=100 * (math.exp(lKp) - 1),
                            fark=lK - lKp, KILITlog=r["KILIT"], kalinti=art)

    # ------------------------------------------------- eski çapa-bağıl
    E = OK["eski_capa"]
    print("\n" + "=" * 74)
    print("  T6 — ESKİ ÇAPA-BAĞIL SAYILAR ↔ 180'in ÇAPASIZ KARŞILIĞI")
    print("=" * 74)
    pm_lin = OUT["dogrusal"]["M"]
    pm_log = OUT["log"]["M"]
    print(f"  ΔM(son↔Hk) doğrusal = {pm_lin['delta']:+.6f}   "
          f"log = {pm_log['delta']:+.6f}  (175/179: "
          f"{E['DlogM_son_Hk']:+.6f})")
    print(f"  ESKİ (çapa: Hkeskin+HA4): kesim %{100*E['mansetkesim_pay']:.2f}"
          f" / kilit %{100*E['manset_kilit_pay']:.2f}   "
          f"[log konvansiyonu, HA4 çapası]")
    print(f"  180  (ÇAPASIZ):  zarf %{100*pm_log['p_zarf']:.2f} / kilit "
          f"%{100*pm_log['p_kilit']:.2f}  (log)   |   zarf "
          f"%{100*pm_lin['p_zarf']:.2f} / kilit %{100*pm_lin['p_kilit']:.2f}"
          f"  (doğrusal)")
    print(f"  ESKİ E-kanal ara-konumu f: {E['f_gE_dogrusal']:.4f} (doğrusal "
          f"g_E) / {E['f_gE_log']:.4f} (log g_E) / {E['f_mu2E']:.4f} (μ̂²_E)"
          f" / {E['f_b_tau070_080']:.4f},{E['f_b_tau080_095']:.4f} (f_b)")
    for k in ("Q_E", "rho_E"):
        print(f"  180 çapasız {AD[k]:5s}: zarf payı doğrusal "
              f"{OUT['dogrusal'][k]['p_zarf']:+.4f} / log "
              f"{OUT['log'][k]['p_zarf']:+.4f}")
    OUT["eski"] = E

    p = S180 / "K3_DEFTER.json"
    p.write_text(json.dumps(OUT, indent=1, ensure_ascii=False, default=float))
    print(f"\n-> {p}", flush=True)


if __name__ == "__main__":
    main()
