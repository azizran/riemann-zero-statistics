# -*- coding: utf-8 -*-
"""
183f — K1(hüküm) + K3 (ÇATAL DEFTERİ) + K4 (180 defterine etki)
================================================================
ÖN-MÜHÜR: 183/ONKAYIT_183.json. Bu betik HİÇBİR YENİ ALAN ÖLÇMEZ;
yalnız diskteki ölçümleri ön-kayıtlı kurallarla birleştirir.

Girdi: 183/DERINLIK_t1.json, 183/KORR_*.json, 183/SURB_VF*.json,
       176/K1_*.json, 174/K1_*.json, 169/K2b_*.json, 182/ANATOMI_VF*.json
Çıktı: 183/CATAL.json + ekrana defter

Kullanım: 183f_catal.py
"""
import json
import sys
from pathlib import Path

import numpy as np

SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S183, S176, S174, S169, S182 = (SCR / "183", SCR / "176", SCR / "174",
                                SCR / "169", SCR / "182")
VF = [f"VF{i}" for i in range(1, 9)]
VS = [f"VS{i}" for i in range(1, 5)]
VD = ["VD0", "VD1", "VD2", "VD4"]


def jload(p):
    return json.load(open(p)) if Path(p).exists() else None


def se(a):
    a = np.asarray(a, float)
    return float(a.std(ddof=1) / np.sqrt(len(a)))


def main():
    ON = jload(S183 / "ONKAYIT_183.json")
    CAPA = ON["CAPA"]
    out = dict(onkayit_damga=ON["damga"], sha_183a=ON["sha_183a"])

    # ==================================================================
    # K1 — DERİNLİK MERDİVENİ + H-İ1
    # ==================================================================
    print("=" * 78)
    print("K1 — DERİNLİK MERDİVENİ (TEŞHİS gazları; inşa kapıları uygulanmaz)")
    print("=" * 78)
    D = jload(S183 / "DERINLIK_t1.json")
    print(f"  R-KAPISI (d=∞ ↔ 176/z_VF1.npy): maks|Δz| = {D['dz_ref']:.3e}  "
          f"{'✓ BİT-BİT' if D['R_KAPISI'] else '✗ TUTMADI'}")
    print(f"  aktif-küme defteri tutarlı: {D['defter_tutarli']}   "
          f"yineleme: {D['n_yineleme']}   aktif dizi: {D['aktif_dizi']}")
    sat = []
    for d in ["0", "1", "2", "4", "inf"]:
        b = D["basamak"][str(d)]
        ad = b["ad"]
        kr = jload(S183 / f"KORR_{ad}.json")
        k1 = jload(S176 / f"K1_{ad}.json")
        kb = jload(S169 / f"K2b_{ad}.json")
        row = dict(d=d, ad=ad, etiket=b["etiket"],
                   n_ihlal=b["n_ihlal"], min_dz=b["min_dz"],
                   maks_kayma=b.get("maks_kayma"),
                   rms_kayma=b.get("rms_kayma"),
                   sigma_ds=b.get("sigma_ds"), L=b.get("L"),
                   korr=kr["korr_XaXb"] if kr else None,
                   korr_sur=kr["korr_XaXb_sur"] if kr else None,
                   R_eta=k1["eta"]["R"] if k1 else None,
                   kovoran=(k1["eta"]["kov_art_x"] / k1["eta"]["Var"]
                            if k1 else None),
                   kov_art=k1["eta"]["kov_art_x"] if k1 else None,
                   Var_eta=k1["eta"]["Var"] if k1 else None,
                   kov_ciz_art=(k1["eta"]["P"] - k1["eta"]["Var_cizgi"]
                                if k1 else None),
                   m3_cizgi=k1["m3"]["cizgi"] if k1 else None,
                   rho3=kb["ortalama"]["111_hepsi"] if kb else None,
                   rhoE=kb["ortalama"]["100_E"] if kb else None)
        sat.append(row)
    print(f"\n  {'d':>4} {'gaz':>5} {'sıra ihl':>9} {'rms|z−z∞|':>11} "
          f"{'korr(Xa,Xb)':>12} {'R_η':>9} {'Kov-oranı':>11} {'ρ₃':>9} "
          f"{'ρ(E)':>8}")
    for r in sat:
        f = lambda v, w, p: (f"{v:>{w}.{p}f}" if v is not None else " " * w)
        print(f"  {r['d']:>4} {r['ad']:>5} {r['n_ihlal']:>9d} "
              f"{r['rms_kayma']:>11.3e} {f(r['korr'],12,5)} {f(r['R_eta'],9,5)} "
              f"{f(r['kovoran'],11,5)} {f(r['rho3'],9,5)} {f(r['rhoE'],8,4)}")
    out["merdiven"] = sat

    # ---- H-İ1 hükmü (ön-kayıtlı) -------------------------------------
    ref_k = ON["HI1"]["referans"]["korr_XaXb(inf)"]
    ref_v = ON["HI1"]["referans"]["kovoran(inf)"]
    tol_k = ON["HI1"]["tol"]["korr_XaXb"]
    tol_v = ON["HI1"]["tol"]["kovoran"]
    kk = [r["korr"] for r in sat]
    vv = [r["kovoran"] for r in sat]
    hi1 = {}
    if all(x is not None for x in kk) and all(x is not None for x in vv):
        f0k, f0v = kk[0] / ref_k, vv[0] / ref_v
        mon_k = all(kk[i + 1] >= kk[i] - tol_k for i in range(len(kk) - 1))
        mon_v = all(vv[i + 1] >= vv[i] - tol_v for i in range(len(vv) - 1))
        yasa = (f0k <= 0.25 and f0v <= 0.25 and mon_k and mon_v)
        olum = (f0k >= 0.75 or f0v >= 0.75)
        hukum = "ÖLDÜ" if olum else ("YAŞIYOR" if yasa else "HÜKÜMSÜZ")
        hi1 = dict(f0_korr=f0k, f0_kovoran=f0v, tekduze_korr=mon_k,
                   tekduze_kovoran=mon_v, hukum=hukum,
                   ref_korr=ref_k, ref_kovoran=ref_v,
                   tol_korr=tol_k, tol_kovoran=tol_v)
        print(f"\n  H-İ1: f(0)_korr = {f0k:+.4f}   f(0)_Kov-oranı = {f0v:+.4f}")
        print(f"        tekdüze? korr {mon_k}   Kov-oranı {mon_v}   "
              f"(tol {tol_k:.5f} / {tol_v:.5f})")
        print(f"        ÖN-KAYITLI HÜKÜM: **{hukum}**")
    out["H_I1"] = hi1

    # ==================================================================
    # K2 — SURROGATE TABANI
    # ==================================================================
    print("\n" + "=" * 78)
    print("K2 — ALAN-SURROGATE TABANI (yeniden çözüm YOK)")
    print("=" * 78)
    AN = [jload(S182 / f"ANATOMI_{g}.json") for g in VF]
    SB = [jload(S183 / f"SURB_{g}.json") for g in VF[:4]]
    SB = [x for x in SB if x]
    taban = {}
    # SUR-A (diskten, n=8)
    a_korr = [a["korr_XaXb_sur"] for a in AN]
    a_r3 = [a["rho_sur"]["hepsi"] for a in AN]
    a_rE = [a["rho_sur"]["E"] for a in AN]
    taban["korr_XaXb"] = dict(kaynak="SUR-A (182, n=8)", n=8,
                              ort=float(np.mean(a_korr)), se=se(a_korr))
    taban["rho3"] = dict(kaynak="SUR-A (182, n=8)", n=8,
                         ort=float(np.mean(a_r3)), se=se(a_r3))
    taban["rho_E"] = dict(kaynak="SUR-A (182, n=8)", n=8,
                          ort=float(np.mean(a_rE)), se=se(a_rE))
    # ρ_J para birimi: VF ve SUR AYNI kestirimciyle (169'un ρ'su değil)
    taban["rho3_J"] = dict(kaynak="SUR-A ρ_J (182, n=8)", n=8,
                           ort=float(np.mean(a_r3)), se=se(a_r3))
    taban["rhoE_J"] = dict(kaynak="SUR-A ρ_J (182, n=8)", n=8,
                           ort=float(np.mean(a_rE)), se=se(a_rE))
    if SB:
        for ad, key in [("R_eta", "R"), ("kovoran", "kovoran"),
                        ("kov_ciz_art", "kov_ciz_art"),
                        ("kov_art", "kov_art_x"), ("Var_eta", "Var")]:
            v = [s["surrogate"][0]["eta"][key] for s in SB]
            taban[ad] = dict(kaynak=f"SUR-B (183e, n={len(SB)})", n=len(SB),
                             ort=float(np.mean(v)), se=se(v))
        v = [s["surrogate"][0]["m3_cizgi"] for s in SB]
        taban["m3_cizgi"] = dict(kaynak=f"SUR-B (183e, n={len(SB)})",
                                 n=len(SB), ort=float(np.mean(v)), se=se(v))
    # SUR-C (türetim)
    taban["mu2_E"] = dict(kaynak="SUR-C (türetim, DEJENERE)", n=0,
                          ort=-CAPA["rho_var_E_VF"], se=0.0)
    for k, v in taban.items():
        print(f"  {k:14s} SUR = {v['ort']:+.6f} ± {v['se']:.6f}   [{v['kaynak']}]")
    out["taban"] = taban

    # ==================================================================
    # K3 — ÇATAL DEFTERİ
    # ==================================================================
    print("\n" + "=" * 78)
    print("K3 — ÇATAL DEFTERİ  [GK_alt = ikiz−VF ,  GK_üst = ikiz−SUR]")
    print("=" * 78)
    # YOL DÜZELTMESİ (kaptan, 9 Eyl): bütün K1_<gaz>.json dosyaları 176d_K1
    # tarafından S176'ya yazılır (S174 boş çıktı); Hkeskin/son da 176d ile
    # üretildi. İçerik/formül değişikliği YOKTUR.
    K1v = {g: jload(S176 / f"K1_{g}.json")
           for g in VF + ["Hkeskin", "son"]}
    KB = {g: jload(S169 / f"K2b_{g}.json") for g in VF}
    vf = {}
    vf["korr_XaXb"] = [a["korr_XaXb"] for a in AN]
    vf["R_eta"] = [K1v[g]["eta"]["R"] for g in VF]
    vf["kovoran"] = [K1v[g]["eta"]["kov_art_x"] / K1v[g]["eta"]["Var"] for g in VF]
    vf["kov_art"] = [K1v[g]["eta"]["kov_art_x"] for g in VF]
    vf["kov_ciz_art"] = [K1v[g]["eta"]["P"] - K1v[g]["eta"]["Var_cizgi"] for g in VF]
    vf["m3_cizgi"] = [K1v[g]["m3"]["cizgi"] for g in VF]
    vf["rho3"] = [KB[g]["ortalama"]["111_hepsi"] for g in VF]
    vf["rho_E"] = [KB[g]["ortalama"]["100_E"] for g in VF]
    vf["Var_eta"] = [K1v[g]["eta"]["Var"] for g in VF]
    vf["mu2_E"] = [json.load(open(S176 / f"G_{g}.json"))["E"]["mu2"] for g in VF]
    vf["rho3_J"] = [a["rho_J"]["hepsi"] for a in AN]
    vf["rhoE_J"] = [a["rho_J"]["E"] for a in AN]

    kson = jload(S183 / "KORR_son.json")      # gerçek gazın korr(Xa,Xb)'si
    ikiz = dict(
        korr_XaXb=None,           # ikiz (Hkeskin) diskte yok — bkz. BORÇ
        rho3_J=None, rhoE_J=None,  # ρ_J ikizde ölçülmedi
        R_eta=CAPA["R_eta_Hkeskin"], kovoran=CAPA["kovoran_Hkeskin"],
        kov_art=CAPA["kov_art_eta_Hkeskin"],
        kov_ciz_art=CAPA["kov_ciz_art_Hkeskin"],
        m3_cizgi=CAPA["m3_cizgi_Hkeskin"], rho3=CAPA["rho3_Hkeskin"],
        rho_E=0.7947, Var_eta=CAPA["Var_eta_Hkeskin"],
        mu2_E=CAPA["mu2E_Hkeskin"])
    gercek = dict(R_eta=CAPA["R_eta_son"], kovoran=CAPA["kovoran_son"],
                  kov_art=CAPA["kov_art_eta_son"],
                  kov_ciz_art=CAPA["kov_ciz_art_son"],
                  m3_cizgi=CAPA["m3_cizgi_son"], rho3=CAPA["rho3_son"],
                  rho_E=0.8229, Var_eta=CAPA["Var_eta_son"],
                  mu2_E=CAPA["mu2E_son"], rho3_J=None, rhoE_J=None,
                  korr_XaXb=(kson["korr_XaXb"] if kson else None))

    catal = {}
    for tepe, TP in [("ikiz(Hkeskin)", ikiz), ("gerçek(son)", gercek)]:
        print(f"\n  --- TEPE = {tepe} ---")
        print(f"  {'O':13s} {'tepe':>10} {'VF±se':>19} {'SUR±se':>19} "
              f"{'GK_alt':>10} {'IK':>10} {'GK_üst':>10} {'IK/GK_üst':>10}")
        for k in ["korr_XaXb", "rho3", "rho3_J", "rho_E", "rhoE_J", "R_eta",
                  "kovoran", "kov_art", "kov_ciz_art", "m3_cizgi", "mu2_E"]:
            if k not in taban or k not in vf:
                continue
            V = float(np.mean(vf[k]))
            Vse = se(vf[k])
            S = taban[k]["ort"]
            Sse = taban[k]["se"]
            I = TP.get(k)
            row = dict(tepe_ad=tepe, tepe=I, VF=V, VF_se=Vse, SUR=S,
                       SUR_se=Sse, IK=V - S, kaynak_SUR=taban[k]["kaynak"])
            if I is not None:
                row.update(GK_alt=I - V, GK_ust=I - S,
                           pay_insa=(V - S) / (I - S) if (I - S) != 0 else None)
            catal.setdefault(tepe, {})[k] = row
            fi = f"{I:+10.5f}" if I is not None else " " * 10
            fa = f"{row['GK_alt']:+10.5f}" if I is not None else " " * 10
            fu = f"{row['GK_ust']:+10.5f}" if I is not None else " " * 10
            fp = (f"{row['pay_insa']:+10.4f}" if row.get("pay_insa") is not None
                  else " " * 10)
            print(f"  {k:13s} {fi} {V:+11.5f}±{Vse:7.5f} {S:+11.5f}±{Sse:7.5f} "
                  f"{fa} {row['IK']:+10.5f} {fu} {fp}")
    out["catal"] = catal
    catal = catal.get("ikiz(Hkeskin)", {})
    out["gercek_son"] = gercek

    # ---- İŞARET TURNUSOLÜ --------------------------------------------
    print("\n  --- İŞARET TURNUSOLÜ (ön-kayıt K3) ---")
    tur = {}
    if "kovoran" in catal and "kov_art" in catal:
        o, h = catal["kovoran"], catal["kov_art"]
        for nm, (a, b) in [("GK_alt", (h["GK_alt"], o["GK_alt"])),
                           ("IK", (h["IK"], o["IK"])),
                           ("GK_ust", (h["GK_ust"], o["GK_ust"]))]:
            Vimp = a / b if b != 0 else float("nan")
            tur[nm] = dict(ham=a, oran=b, Var_ima=Vimp)
            print(f"    {nm:8s}: HAM {a:+.6f}  ORAN {b:+.6f}  ⇒ ima edilen "
                  f"Var(η) = {Vimp:+.6f}")
        print(f"    ölçülen Var(η): ikiz {CAPA['Var_eta_Hkeskin']:.6f}   "
              f"VF {float(np.mean(vf['Var_eta'])):.6f}   "
              f"SUR {taban.get('Var_eta',{}).get('ort',float('nan')):.6f}")
        tur["Var_olculen"] = dict(ikiz=CAPA["Var_eta_Hkeskin"],
                                  VF=float(np.mean(vf["Var_eta"])),
                                  SUR=taban.get("Var_eta", {}).get("ort"))
    out["turnusol"] = tur

    # ==================================================================
    # K4 — 180 DEFTERİNE ETKİ: inşa-kilidi VS ailesinde de var mı?
    # ==================================================================
    print("\n" + "=" * 78)
    print("K4 — İNŞA-KİLİDİ İKİ AİLEDE ORTAK MI? (ZARF = ⟨N⟩_VS − ⟨N⟩_VF)")
    print("=" * 78)
    k4 = {}
    for nm, fam in [("VF", VF), ("VS", VS)]:
        kr = [jload(S183 / f"KORR_{g}.json") for g in fam]
        kr = [x for x in kr if x]
        if nm == "VF":
            kv = vf["korr_XaXb"]
        else:
            kv = [x["korr_XaXb"] for x in kr]
        k1 = [jload(S176 / f"K1_{g}.json") for g in fam]
        k1 = [x for x in k1 if x]
        ko = [x["eta"]["kov_art_x"] / x["eta"]["Var"] for x in k1]
        k4[nm] = dict(n_korr=len(kv), korr=float(np.mean(kv)),
                      korr_se=se(kv) if len(kv) > 1 else 0.0,
                      n_kov=len(ko), kovoran=float(np.mean(ko)),
                      kovoran_se=se(ko) if len(ko) > 1 else 0.0,
                      korr_hepsi=[float(x) for x in kv],
                      kovoran_hepsi=[float(x) for x in ko])
        print(f"  {nm}: korr(Xa,Xb) = {k4[nm]['korr']:+.5f} ± "
              f"{k4[nm]['korr_se']:.5f} (n={k4[nm]['n_korr']})   "
              f"Kov-oranı = {k4[nm]['kovoran']:+.5f} ± "
              f"{k4[nm]['kovoran_se']:.5f} (n={k4[nm]['n_kov']})")
    if "VS" in k4 and k4["VS"]["n_korr"]:
        for nm in ["korr", "kovoran"]:
            d = k4["VF"][nm] - k4["VS"][nm]
            s = float(np.hypot(k4["VF"][nm + "_se"], k4["VS"][nm + "_se"]))
            ortak = bool(abs(d) <= 2 * s)
            k4[f"fark_{nm}"] = dict(delta=d, se=s, z=d / s if s else None,
                                    ORTAK=ortak)
            print(f"  Δ{nm} (VF−VS) = {d:+.5f} ± {s:.5f}  "
                  f"(|z| = {abs(d)/s:.2f})  ⇒ "
                  f"{'ORTAK (fark ≤ 2se) ⇒ ZARF payında SADELEŞİR' if ortak else 'ORTAK DEĞİL'}")
    out["K4"] = k4

    S183.mkdir(parents=True, exist_ok=True)
    p = S183 / "CATAL.json"
    p.write_text(json.dumps(out, indent=1, ensure_ascii=False, default=float))
    print(f"\n  -> {p}")
    return out


if __name__ == "__main__":
    main()
