"""
165 — HÜKÜM TABLOLARI (K_*.json → rapor satırları)

Kullanım:  165_analiz.py
"""
import importlib
import json
import sys
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
sys.path.insert(0, str(QM / "165_configs"))
K = importlib.import_module("165_cekirdek")
SCR = K.SCR
GAZ = ("Hkeskin", "HA4", "son")


def yukle(v):
    p = SCR / f"K_{v}_t0.4.json"
    return json.loads(p.read_text()) if p.exists() else None


def saglik(v, b):
    """160'ın sağlık kuralı: A²s2 ölçülen anlamlı mı (|Γ| patlaması)."""
    return not (v == "HA4" and b["tau_eff"] > 0.80)


def main():
    D = {v: yukle(v) for v in GAZ}
    print("=" * 96)
    print("T1 — (b)/(a) BANT BANT  [τ∈0.52–0.80, taban 0.40, τ_c=0.95]")
    for kaynak in ("olculen", "nominal"):
        et = f"{kaynak}_tc0.95"
        print(f"\n--- kaynak = {kaynak} ---")
        print(" gaz       τ_eff   A²s2 ölç    HAM      ham/ölç   NRM"
              "      nrm/ölç  ±%25?  Ç1/tot  Ç3/tot   Ç2/tot")
        say = {}
        for v in GAZ:
            if not D[v]:
                continue
            for b in D[v]["kosum"][et]["bant"]:
                if not b.get("olculdu") or b["lo"] < 0.52 - 1e-9:
                    continue
                if not saglik(v, b):
                    print(f" {v:9s} {b['tau_eff']:.4f}  ‡ sağlıksız bant "
                          f"(A²s2={b['A2s2']:+.3f}) — hükme girmiyor")
                    continue
                m, p, pn = b["A2s2"], b["A2s2p"], b["A2s2p_nrm"]
                r = pn / m if m else float("nan")
                ok = "✓" if 0.75 <= r <= 1.25 else " "
                say[v] = say.get(v, [0, 0])
                say[v][1] += 1
                say[v][0] += 1 if ok == "✓" else 0
                tot = p if abs(p) > 1e-300 else float("nan")
                print(f" {v:9s} {b['tau_eff']:.4f} {m:+.6f} {p:+.6f} "
                      f"{p/m if m else float('nan'):+8.3f} {pn:+.6f} "
                      f"{r:+8.3f}   {ok}   {b['A2s2_C1']/tot:+6.3f} "
                      f"{b['A2s2_C3']/tot:+6.3f} {b['A2s2_C2']/tot:+9.1e}")
        tt = [0, 0]
        for v, s in say.items():
            print(f"   {v}: ±%25 içinde {s[0]}/{s[1]}")
            tt[0] += s[0]
            tt[1] += s[1]
        print(f"   TOPLAM: ±%25 içinde {tt[0]}/{tt[1]}")

    print("=" * 96)
    print("KANAL PAYLARI (τ_c=0.95, ölçülen çizgiler) — Ç1 baskın mı?")
    print(" gaz       τ_eff   Ç1        Ç1(T_a)   Ç1(B+C)   Ç2         Ç3"
          "        Ç1/tot  Ç3/tot")
    for v in GAZ:
        if not D[v]:
            continue
        for b in D[v]["kosum"]["olculen_tc0.95"]["bant"]:
            if not b.get("olculdu") or b["lo"] < 0.52 - 1e-9:
                continue
            if not saglik(v, b):
                continue
            tot = b["A2s2p"]
            print(f" {v:9s} {b['tau_eff']:.4f} {b['A2s2_C1']:+.6f} "
                  f"{b['A2s2_C1_Ta']:+.2e} {b['A2s2_C1_BC']:+.6f} "
                  f"{b['A2s2_C2']:+.2e} {b['A2s2_C3']:+.6f} "
                  f"{b['A2s2_C1']/tot:+6.3f} {b['A2s2_C3']/tot:+6.3f}")

    print("=" * 96)
    print("τ_c TARAMASI — çizgi gösteriminin geçerlilik sınırı")
    print(" gaz      kaynak    τ_c   nline  Var(E)mod/ölç Var(X)mod/ölç"
          "  g_E    g_X    g      ⟨nrm/ölç⟩(0.54–0.70)")
    for v in GAZ:
        if not D[v]:
            continue
        for et, r in sorted(D[v]["kosum"].items()):
            A = r["artik"]
            rr = []
            for b in r["bant"]:
                if (b.get("olculdu") and 0.52 <= b["lo"] and
                        b["tau_eff"] < 0.72 and saglik(v, b) and b["A2s2"]):
                    rr.append(b["A2s2p_nrm"] / b["A2s2"])
            print(f" {v:9s} {r['kaynak']:8s} {r['tau_c']:.2f} {A['nline']:6d} "
                  f"{A['varE_mod']/A['varE_olc']:11.3f} "
                  f"{A['varX_mod']/A['varX_olc']:12.3f} {A['gE']:6.3f} "
                  f"{A['gX']:6.3f} {A['gcal']:6.3f}   "
                  f"{np.mean(rr) if rr else float('nan'):+.3f}")

    print("=" * 96)
    print("T4 — GERÇEK ↔ SADAKATLİ: Ç3 payı farkı (HL üst sınırı)")
    print("  τ_eff    Ç3/tot son  Ç3/tot Hkeskin  Ç3/tot HA4   |fark| (son-sent)")
    et = "olculen_tc0.95"
    if all(D[v] for v in GAZ):
        Bs = {v: {round(b["tau"], 3): b for b in D[v]["kosum"][et]["bant"]
                  if b.get("olculdu")} for v in GAZ}
        for t in sorted(Bs["son"]):
            if t < 0.53:
                continue
            try:
                a = Bs["son"][t]["A2s2_C3"] / Bs["son"][t]["A2s2p"]
                b1 = Bs["Hkeskin"][t]["A2s2_C3"] / Bs["Hkeskin"][t]["A2s2p"]
                b2 = Bs["HA4"][t]["A2s2_C3"] / Bs["HA4"][t]["A2s2p"]
            except (KeyError, ZeroDivisionError):
                continue
            if not saglik("HA4", Bs["HA4"][t]):
                b2 = float("nan")
            print(f"  {Bs['son'][t]['tau_eff']:.4f}  {a:+.4f}      "
                  f"{b1:+.4f}         {b2:+.4f}      "
                  f"{abs(a-b1):.4f} / {abs(a-b2):.4f}")


if __name__ == "__main__":
    main()
