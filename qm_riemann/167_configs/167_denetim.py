"""
167 — DENETİM
=============
V1  167'nin ölçüm/öngörü zinciri = 166'nın kayıtlı JSON'u mu?  (BİT DÜZEYİ)
    Karşılaştırılan: Ms0..Ms3, Mu0..Mu3, Ps0..Ps3, Pu0..Pu3, KALİB_s2,
    KALİB_u2, τ_eff, C1, C2, capla + jackknife'lar; üç gazın BÜTÜN
    ölçülen bantlarında.
V2  W çarpanları = 166_bacak'ın bant kayıtları mı? (aynı gp-ağırlığı)
V3  ablasyon gazlarının SADAKAT denetimi (164 ölçütü) — 167_analiz'de
    bant bant basılır; burada gaz özeti.
"""
import importlib
import json
import sys
from pathlib import Path

import numpy as np

SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S166, S167 = SCR / "166", SCR / "167"
ANAHTAR = ([f"M{c}{k}" for c in "su" for k in range(4)]
           + [f"P{c}{k}" for c in "su" for k in range(4)]
           + ["KALIB_s2", "KALIB_u2", "tau_eff", "C1", "C2", "C3", "capla",
              "sMs2", "sPs2", "sMu2", "sPu2", "sKALIB_s2", "sKALIB_u2"])


def main(gazlar=("son", "Hkeskin", "HA4")):
    print("=== V1 — 167 zinciri ↔ 166'nın kayıtlı T1 JSON'u (bit düzeyi) ===")
    genel = 0.0
    for g in gazlar:
        p6 = S166 / f"T1_{g}_olculen_tc0.95.json"
        p7 = S167 / f"C_{g}.json"
        if not (p6.exists() and p7.exists()):
            print(f"  {g}: dosya yok ({p6.exists()}/{p7.exists()})")
            continue
        A = json.load(open(p6))
        B = json.load(open(p7))
        a = {round(b["lo"], 4): b for b in A["bant"] if b.get("olculdu")}
        b_ = {round(b["lo"], 4): b for b in B["bant"] if b.get("olculdu")}
        ort = sorted(set(a) & set(b_))
        mx, nk, arg = 0.0, 0, ""
        for lo in ort:
            for k in ANAHTAR:
                if k in a[lo] and k in b_[lo]:
                    d = abs(float(a[lo][k]) - float(b_[lo][k]))
                    nk += 1
                    if d > mx:
                        mx, arg = d, f"{k}@lo={lo}"
        genel = max(genel, mx)
        print(f"  {g:9s}: {len(ort)} ortak bant, {nk} karşılaştırma, "
              f"maks fark = {mx:.1e}  ({arg})")
    print(f"  GENEL MAKS FARK = {genel:.1e}")

    print("\n=== V2 — W çarpanları ↔ 166_bacak'ın bant kayıtları ===")
    for g in gazlar:
        p6 = S166 / f"BACAK_{g}.json"
        p7 = S167 / f"C_{g}.json"
        if not (p6.exists() and p7.exists()):
            continue
        A = {round(x["lo"], 4): x for x in json.load(open(p6))["bant"]}
        B = json.load(open(p7))
        mx, arg = 0.0, ""
        for b in B["bant"]:
            if not b.get("olculdu"):
                continue
            lo = round(b["lo"], 4)
            if lo not in A:
                continue
            for k in ("W_amp", "W_X", "W_pos"):
                d = abs(A[lo][k] - b[k])
                if d > mx:
                    mx, arg = d, f"{k}@lo={lo}"
        print(f"  {g:9s}: maks fark = {mx:.1e}  ({arg})")

    print("\n=== V3 — SADAKAT (164 ölçütü) gaz özeti ===")
    print("  gaz        bant(lo≥0.52)   R_bant aralığı        SNR aralığı")
    for p in sorted(S167.glob("C_*.json")):
        d = json.load(open(p))
        R = [b["R_bant"] for b in d["bant"]
             if b.get("olculdu") and b["lo"] >= 0.52 - 1e-9]
        S = [b["SNR"] for b in d["bant"]
             if b.get("olculdu") and b["lo"] >= 0.52 - 1e-9]
        if not R:
            continue
        print(f"  {p.stem[2:]:10s} {len(R):3d}          "
              f"{min(R):.3f} – {max(R):.3f}      {min(S):8.1f} – {max(S):8.1f}")


if __name__ == "__main__":
    main(tuple(sys.argv[1:]) or ("son", "Hkeskin", "HA4"))
