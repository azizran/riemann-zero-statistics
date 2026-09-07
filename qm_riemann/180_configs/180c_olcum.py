# -*- coding: utf-8 -*-
"""
180c — K1 (ÖLÇÜM): VS GAZINDA 167 ZİNCİRİ + 172b DEFTERİ
=========================================================
ÖN-MÜHÜR. `176c_olcum.kos` ile AYNI çağrı zinciri:
    167_olcum.kos(ad, 0.40, 0.95, kule=0, duz=0)
sonra 172b'nin `oku`/`kusur`'uyla HAM defter satırı (M, Q_E, ρ_E, Q_X,
ρ_X) ve EK oran satırı (g_E, g_X, θ, g_cal).

NOT (180a'da dondurulmuş): `167_olcum.kos` VS'nin künyesini bulamayınca
λ=1.00, τ_ust=1.00, pen=None varsayar ⇒ `b_nom` HKESKİN'in nominalidir.
`KALIB_u2` (=M) b_nom'a BAĞLI DEĞİLDİR (b_nom yalnız R_bant'a girer), bu
yüzden defter satırı VF ile birebir karşılaştırılabilir. R_bant burada
KAPI değil TANIdır: R_bant(VS)/R_bant(VF) ≈ r(τ) beklenir — zarf
aktarımının bağımsız sağlaması.

Kullanım: 180c_olcum.py VS1 VS2
"""
import importlib
import json
import sys
import time
from pathlib import Path

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
for _p in ("180_configs", "176_configs", "172_configs", "167_configs"):
    sys.path.insert(0, str(QM / _p))
O167 = importlib.import_module("167_olcum")
B172 = importlib.import_module("172b_gE_yasasi")

SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S180 = SCR / "180"
S167 = SCR / "167"


def kos(ad):
    t0 = time.time()
    print("=" * 74, flush=True)
    print(f"180c — VS ÖLÇÜMÜ {ad}  (167_olcum.kos, taban 0.40, τ_c 0.95, "
          f"kule 0, düz 0 — VF ile AYNI çağrı)", flush=True)
    print("=" * 74, flush=True)
    O167.kos(ad, 0.40, 0.95, 0, 0)

    v = B172.oku(ad)
    E = B172.kusur(v["varE_mod"], v["varE_olc"], v["varE_res"], v["korE"],
                   v["gE"])
    X = B172.kusur(v["varX_mod"], v["varX_olc"], v["varX_res"], v["korX"],
                   v["gX"])
    rec = dict(ad=ad, lam=1.0, N=v["N"], T=v["T"], nline=v["nline"],
               M=v["KAL"], theta=v["th"], g_E=v["gE"], g_X=v["gX"],
               g_cal=v["gcal"], W_X=v["WX"], sigX=v["sigX"],
               sigds=v["sigds"], sigC=v["sigC"],
               Q_E=E["Q"], rho_E=E["rho"], pi_E=E["pi"], mu2_E=E["mu2"],
               r_E=E["r"], dK_E=E["dK"],
               Q_X=X["Q"], rho_X=X["rho"], pi_X=X["pi"], mu2_X=X["mu2"],
               dK_X=X["dK"])

    d = json.load(open(S167 / f"C_{ad}.json"))
    rec["bant"] = [dict(lo=b["lo"], tau_eff=b["tau_eff"],
                        R_bant=b["R_bant"], KALIB_u2=b["KALIB_u2"])
                   for b in d["bant"] if b.get("olculdu")]
    rec["c_gaz"] = d.get("c_gaz")
    rec["tau_c_faz"] = d["tau_c_faz"]
    rec["tau_c_gap"] = d["tau_c_gap"]

    S180.mkdir(parents=True, exist_ok=True)
    p = S180 / f"G_{ad}.json"
    p.write_text(json.dumps(rec, indent=1, ensure_ascii=False, default=float))

    print("\n" + "-" * 74, flush=True)
    print(f"  HAM DEFTER ({ad}):", flush=True)
    print(f"    M = {rec['M']:.6f}   Q_E = {rec['Q_E']:.6f}   "
          f"ρ_E = {rec['rho_E']:.6f}   Q_X = {rec['Q_X']:.6f}   "
          f"ρ_X = {rec['rho_X']:.6f}", flush=True)
    print(f"  EK (oran dili): g_E = {rec['g_E']:.6f}  g_X = "
          f"{rec['g_X']:.6f}  θ = {rec['theta']:.6f}  g_cal = "
          f"{rec['g_cal']:.6f}", flush=True)
    print(f"    π_E = {rec['pi_E']:.6f}  μ̂²_E = {rec['mu2_E']:+.6f}  "
          f"r_E = {rec['r_E']:.6f}   |ΔK|/K: E {rec['dK_E']:.2e} "
          f"X {rec['dK_X']:.2e}", flush=True)
    print(f"    R_bant: " + " ".join("%.4f" % b["R_bant"]
                                     for b in rec["bant"]), flush=True)
    print(f"  -> {p}   ({(time.time()-t0)/60:.1f} dk)", flush=True)
    return rec


if __name__ == "__main__":
    for g in (sys.argv[1:] or ["VS1"]):
        kos(g)
