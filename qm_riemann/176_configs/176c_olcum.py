# -*- coding: utf-8 -*-
"""
176c — K1 (ÖLÇÜM): VEKİL GAZDA 167 ZİNCİRİ + ÇARPAN DEFTERİ
============================================================
ÖN-MÜHÜR (koşudan önce yazıldı).

`167_olcum.kos` AYNEN import edilir ve vekil gazda koşulur
(taban 0.40, τ_c 0.95, kule 0, düz 0 — λ/kesim ailesinin gazlarıyla
AYNI çağrı). Sonra 172b'nin `oku`/`kusur` fonksiyonlarıyla — yine
AYNEN import edilerek — gazın G1-satırı hesaplanır:
    g_E, g_X, θ = KALİB_u2/g_cal (lo = 0.60 bandı), M = KALİB_u2,
    π_E = K_E/V_O , ρ_E , Q_E , μ̂²_E , r_E .
Vekilin künyesi 167_ortak.KUNYE'de YOKTUR; `kos` bu durumda
λ = 1.00, τ_ust = 1.00, pen = None varsayar — yani tam olarak
`Hkeskin`'in zarfı. Bu, ön-kaydın "vekilin zarfı gerçeğin zarfıdır"
maddesinin ölçüm tarafındaki karşılığıdır (b_nom aynı formülden çıkar).

KAPI (176a/F0): hüküm bantlarında (lo ∈ [0.52, 0.68]) R_bant ≥ 0.98.

Kullanım: 176c_olcum.py <ad>
"""
import importlib
import json
import sys
import time
from pathlib import Path

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
for _p in ("176_configs", "172_configs", "167_configs"):
    sys.path.insert(0, str(QM / _p))
O167 = importlib.import_module("167_olcum")
B172 = importlib.import_module("172b_gE_yasasi")

SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S176 = SCR / "176"
LO_LO, LO_HI = 0.52, 0.68          # 169 filtresi / hüküm penceresi


def kos(ad):
    t0 = time.time()
    print("=" * 74, flush=True)
    print(f"176c — VEKİL ÖLÇÜMÜ {ad}   (167_olcum.kos, taban 0.40, "
          f"τ_c 0.95, kule 0, düz 0)", flush=True)
    print("=" * 74, flush=True)
    O167.kos(ad, 0.40, 0.95, 0, 0)

    # --- 172b'nin defter satırı ---------------------------------------
    v = B172.oku(ad)
    E = B172.kusur(v["varE_mod"], v["varE_olc"], v["varE_res"], v["korE"],
                   v["gE"])
    X = B172.kusur(v["varX_mod"], v["varX_olc"], v["varX_res"], v["korX"],
                   v["gX"])
    rec = dict(lam=1.0, N=v["N"], T=v["T"], nline=v["nline"], th=v["th"],
               gE=v["gE"], gX=v["gX"], KAL=v["KAL"], WX=v["WX"],
               gcal=v["gcal"], sigX=v["sigX"], sigds=v["sigds"],
               sigC=v["sigC"], E=E, X=X)

    # --- F0 kapısı: R_bant --------------------------------------------
    d = json.load(open(SCR / "167" / f"C_{ad}.json"))
    H = [b for b in d["bant"] if b.get("olculdu")
         and LO_LO - 1e-9 <= b["lo"] <= LO_HI + 1e-9]
    Rb = [b["R_bant"] for b in H]
    rec["R_bant_hukum"] = Rb
    rec["R_bant_min"] = min(Rb) if Rb else float("nan")
    rec["F0_R_bant"] = bool(Rb and min(Rb) >= 0.98)
    rec["c_gaz"] = d.get("c_gaz")
    rec["tau_c_faz"] = d["tau_c_faz"]
    rec["tau_c_gap"] = d["tau_c_gap"]
    rec["kule"] = d["kule"]
    rec["duz"] = d["duz"]

    S176.mkdir(parents=True, exist_ok=True)
    p = S176 / f"G_{ad}.json"
    p.write_text(json.dumps(rec, indent=1, ensure_ascii=False))

    print("\n" + "-" * 74, flush=True)
    print(f"  DEFTER SATIRI ({ad}):", flush=True)
    print(f"    g_E = {rec['gE']:.6f}   g_X = {rec['gX']:.6f}   "
          f"g_cal = {rec['gcal']:.6f}", flush=True)
    print(f"    M (KALİB_u2, lo=0.60) = {rec['KAL']:.6f}   "
          f"θ = M/g_cal = {rec['th']:.6f}", flush=True)
    print(f"    π_E = {E['pi']:.6f}  ρ_E = {E['rho']:.6f}  "
          f"Q_E = {E['Q']:.6f}  μ̂²_E = {E['mu2']:+.6f}  r_E = {E['r']:.6f}",
          flush=True)
    print(f"    π_X = {X['pi']:.6f}  μ̂²_X = {X['mu2']:+.3e}   "
          f"σ_ds = {rec['sigds']:.6f}  σ_Ĉ = {rec['sigC']:.6f}", flush=True)
    print(f"    özdeşlik denetimi |ΔK|/K:  E {E['dK']:.2e}   X {X['dK']:.2e}",
          flush=True)
    print(f"    F0 R_bant (lo∈[0.52,0.68]): "
          f"{' '.join('%.4f' % r for r in Rb)}  ⇒ min {rec['R_bant_min']:.4f}"
          f"  {'✓' if rec['F0_R_bant'] else '✗ KAPI TUTMADI'}", flush=True)
    print(f"  -> {p}   ({(time.time()-t0)/60:.1f} dk)", flush=True)
    return rec


if __name__ == "__main__":
    for g in (sys.argv[1:] or ["VF1"]):
        kos(g)
