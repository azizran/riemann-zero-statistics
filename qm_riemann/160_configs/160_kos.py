"""
160 — KOŞU SÜRÜCÜSÜ.   Kullanım: 160_kos.py <veri> <taban> <izgara>

izgara: t1   → izgara(0.44,0.80,0.04)   (bant bant muhasebe, τ∈(0.42,0.80))
        g158 → izgara(0.28,0.64,0.02)   (158'in ızgarası: denetim + a/b)

veri  : son | orta (gerçek) · keskin | A4 | P1 (sentetik)
Taban altı bantlar elenir (η regresyonu τ ≤ taban çizgilerini çıkarır).
"""
import importlib
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
C = importlib.import_module("160_cekirdek")

IZG = {"t1": C.IZGARA_T1, "g158": C.IZGARA158}

if __name__ == "__main__":
    veri, taban, gad = sys.argv[1], float(sys.argv[2]), sys.argv[3]
    t0 = time.time()
    C.SCR.mkdir(parents=True, exist_ok=True)
    z = C.KOS155.veri_yukle(veri)
    bant = [b for b in IZG[gad] if b[0] >= taban - 1e-9]
    print(f"=== 160 / veri={veri} / taban={taban} / izgara={gad} ===")
    print(f"    n={len(z)}  t∈[{z[0]:.1f},{z[-1]:.1f}]  {len(bant)} bant "
          f"({bant[0][0]}–{bant[-1][1]}); {len(IZG[gad])-len(bant)} elendi",
          flush=True)
    r = C.olc160(z, f"{veri}-t{taban}-{gad}", bant, anahtar=veri, taban=taban,
                 cap=4000)
    r.update(veri=veri, izgara=gad, t_lo=float(z[0]), t_hi=float(z[-1]),
             sure_s=time.time() - t0)
    p = C.SCR / f"S_{veri}_t{taban}_{gad}.json"
    p.write_text(json.dumps(r, indent=1))
    print(f"\n-> {p}   ({(time.time()-t0)/60:.1f} dk)", flush=True)
