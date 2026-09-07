# -*- coding: utf-8 -*-
"""
180e — K2 (HÜKÜM): H-180a ÇİFTE-SAYIM HAKEMİ
=============================================
ÖN-MÜHÜR: karar kuralı `180/ONKAYIT_K0.json` içinde, ölçümden önce
dondurulmuştu. Bu betik SADECE o kuralı diskteki `169/K2b_*.json`
çıktılarına uygular; eşik, gaz kümesi ya da istatistik SEÇİLMEZ.

Kullanım: 180e_hakem.py
"""
import json
import sys
import time
from pathlib import Path

import numpy as np

SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S180, S169 = SCR / "180", SCR / "169"


def main():
    OK = json.load(open(S180 / "ONKAYIT_K0.json"))
    H = OK["H180a"]
    G3 = H["capalar"]["GAUSS3"]
    RHK = H["capalar"]["rho3_Hkeskin"]
    RSON = H["capalar"]["rho3_son"]
    gazlar = H["gazlar"]
    print("=" * 74, flush=True)
    print("180e — H-180a ÇİFTE-SAYIM HAKEMİ (ön-kayıtlı kural)", flush=True)
    print(f"    ön-kayıt {OK['zaman']}  sha={OK['sha256'][:16]}…", flush=True)
    print("=" * 74, flush=True)

    ad = ["100_E", "010_Xa", "001_Xb", "011_XaXb", "110_EXa", "111_hepsi"]
    T = {}
    for g in ("Hkeskin", "son") + tuple(gazlar):
        p = S169 / f"K2b_{g}.json"
        if not p.exists():
            print(f"  ({g}: K2b yok)"); continue
        T[g] = json.load(open(p))["ortalama"]
    print("\n  T5 — kırpma aktarımı tablosu (169_k2b, lo∈[0.52,0.68] ort.)")
    print("  gaz       " + " ".join(f"{a.split('_')[1]:>8s}" for a in ad))
    for g, v in T.items():
        print(f"  {g:9s} " + " ".join(f"{v[a]:8.4f}" for a in ad), flush=True)
    print(f"  GAUSS     " + " ".join(
        f"{(2/np.pi)**(0.5*(1 if a.split('_')[0].count('1')==1 else 0)+0):8.4f}"
        if False else f"{np.sqrt(2/np.pi)**a.split('_')[0].count('1'):8.4f}"
        for a in ad), flush=True)

    r3 = np.array([T[g]["111_hepsi"] for g in gazlar if g in T])
    n = len(r3)
    if n < 2:
        raise SystemExit("VF ölçümleri eksik — hüküm verilmez.")
    bar = float(r3.mean())
    sd = float(np.std(r3, ddof=1))
    se = sd / np.sqrt(n)
    C = (RHK - bar) / (RHK - G3)
    zG = (bar - G3) / se if se > 0 else float("inf")
    sacilim = float(r3.max() - r3.min())
    esik_sacilim = 0.5 * abs(RHK - G3)

    print(f"\n  ⟨ρ₃⟩_VF = {bar:.6f}   sd = {sd:.6f}   se = {se:.6f}   "
          f"(n = {n}; tek tek: {' '.join('%.4f' % x for x in r3)})")
    print(f"  ρ₃(Hkeskin) = {RHK:.4f}   ρ₃(son) = {RSON:.4f}   "
          f"GAUSS3 = {G3:.6f}")
    print(f"  ÇÖKÜŞ KESRİ  C = ({RHK:.4f} − {bar:.6f})/({RHK:.4f} − "
          f"{G3:.6f}) = {C:+.4f}")
    print(f"  z_G = ({bar:.6f} − {G3:.6f})/{se:.6f} = {zG:+.3f}")
    print(f"  tohum saçılımı = {sacilim:.6f}   (HÜKÜMSÜZ eşiği "
          f"{esik_sacilim:.6f})")

    if sacilim > esik_sacilim:
        dal, hukum = "HUKUMSUZ", ("tohum saçılımı çözünürlüğü aşıyor — "
                                  + OK["H180a"]["HUKUMSUZ"]["hukum"])
    elif C >= 0.70 and abs(zG) <= 3:
        dal, hukum = "A", OK["H180a"]["DAL_A"]["hukum"]
    elif C <= 0.30:
        dal, hukum = "B", OK["H180a"]["DAL_B"]["hukum"]
    else:
        dal, hukum = "HUKUMSUZ", OK["H180a"]["HUKUMSUZ"]["hukum"]

    print(f"\n  >>> ÖN-KAYITLI DAL: {dal}")
    print(f"  >>> {hukum}", flush=True)

    # ek: bacak-sayısı ortalamaları (kayıt için, karar dışı)
    NB = {"100_E": 1, "010_Xa": 1, "001_Xb": 1, "011_XaXb": 2,
          "110_EXa": 2, "111_hepsi": 3}
    ek = {}
    for g, v in T.items():
        ek[g] = {k: float(np.mean([v[a] for a in ad if NB[a] == k]))
                 for k in (1, 2, 3)}
    print("\n  bacak-sayısı ortalamaları (kayıt; karar YALNIZ ρ₃'e bağlı)")
    for k in (1, 2, 3):
        gv = np.sqrt(2 / np.pi) ** k
        print(f"   {k} bacak (Gauss {gv:.4f}): " + "  ".join(
            f"{g}={ek[g][k]:.4f}({100*(ek[g][k]/gv-1):+.1f}%)"
            for g in ek), flush=True)

    rec = dict(dal=dal, hukum=hukum, rho3=dict(zip(
        [g for g in gazlar if g in T], r3.tolist())),
        rho3_bar=bar, sd=sd, se=se, C=C, z_G=zG, sacilim=sacilim,
        esik_sacilim=esik_sacilim, rho3_Hk=RHK, rho3_son=RSON, GAUSS3=G3,
        tablo=T, bacak_ort=ek, zaman=time.ctime())
    p = S180 / "K2_HAKEM.json"
    p.write_text(json.dumps(rec, indent=1, ensure_ascii=False, default=float))
    print(f"\n-> {p}", flush=True)


if __name__ == "__main__":
    main()
