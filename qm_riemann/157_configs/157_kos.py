"""
157 — KOŞU SÜRÜCÜSÜ
===================
Kullanım:  157_kos.py <veri> <taban> <bantset>

veri    : son | orta     (gerçek; 155/154'ün pencereleri, birebir)
            son  = zeros6 son 300k              L = 12.030
            orta = zeros6 Z[850000:1150000]     L = 11.464
taban   : η regresyonunun asal-çizgi üst sınırı (0.34 / 0.40 / 0.46 / 0.52)
bantset : ince  = 0.02 ızgara, τ ∈ (0.42, 0.80]   (19 bant; 155'in
                  `ince`+`yuksek` ızgaralarını birebir kapsar)
          kaba  = 0.03 ızgara, τ ∈ (0.43, 0.79]   (12 bant; kenarlar hem
                  0.46'ya hem 0.52'ye düşer → taban değişince ızgara
                  KAYMAZ, yalnız alt bantlar elenir)

TABAN KISITI: zincir3, τ ≤ taban çizgilerini regresyonda ÇIKARIR; bir
bant ancak lo ≥ taban ise ölçülebilir. Sürücü alt bantları otomatik eler
ve kaçının elendiğini yazar.
"""
import importlib
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
K = importlib.import_module("157_cekirdek")

HERE = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
OUT = K.SCR
PENCERE = {"son": (-300000, None), "orta": (850000, 1150000)}


def veri_yukle(ad):
    d = np.load(HERE / "128_odl_zeros6_2e6_zeros.npz")
    Z = np.sort(np.asarray(d["zeros"], dtype=float))
    a, b = PENCERE[ad]
    return Z[len(Z) + a:] if b is None else Z[a:b]


if __name__ == "__main__":
    veri, taban, bset = sys.argv[1], float(sys.argv[2]), sys.argv[3]
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    z = veri_yukle(veri)
    bant = {"ince": K.INCE, "kaba": K.KABA, "tau0": K.TAU0}[bset]
    elenen = [b for b in bant if b[0] < taban - 1e-9]
    bant = [b for b in bant if b[0] >= taban - 1e-9]
    print(f"=== 157 / veri={veri} / taban={taban} / bant={bset} ===")
    print(f"    n={len(z)}  t∈[{z[0]:.1f},{z[-1]:.1f}]  {len(bant)} bant "
          f"({bant[0][0]}–{bant[-1][1]}); {len(elenen)} bant taban altında "
          f"elendi", flush=True)
    r = K.olc157(z, f"{veri}-t{taban}-{bset}", bant, anahtar=veri,
                 taban=taban, cap=4000)
    r["veri"] = veri; r["bantset"] = bset
    r["t_lo"] = float(z[0]); r["t_hi"] = float(z[-1])
    r["sure_s"] = time.time() - t0
    p = OUT / f"R_{veri}_t{taban}_{bset}.json"
    p.write_text(json.dumps(r, indent=1))
    print(f"\n-> {p}   ({(time.time()-t0)/60:.1f} dk)", flush=True)
