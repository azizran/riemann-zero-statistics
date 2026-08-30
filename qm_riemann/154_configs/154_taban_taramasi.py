"""
154 — TABAN TARAMASI: R(τ) τ'nun mu, (τ − taban)'ın mı fonksiyonu?
==================================================================
GEREKÇE: ölçülen φ_Γ düşük τ'da neredeyse doğrusal ve sıfırı τ≈0.514'e,
yani REGRESYON TABANINA (0.52) düşüyor. Aday (c) zaten (τ−0.52) motifini
varsayıyor. O motif fizik mi konvansiyon mu — tek ayırt edici deney:
tabanı kaydır.
  taban 0.46 / 0.52 / 0.58, ÜÇÜNDE DE ortak bant ızgarası (0.03 adım).
  Örtüşme bölgesi τ ∈ (0.585, 0.795].
  R(τ) örtüşmede AYNI kalırsa  → R gerçek bir τ fonksiyonu; (τ−taban)
                                  motifi (aday c) ÖLÜR.
  R(τ) tabanla kayarsa         → R taban-göreli; kapalı form ancak
                                  (τ−taban) değişkeninde yazılabilir.
Not: cap=4000 (standart 720 cap'i taban>0.55'te devreye girip "taban
değişti" iddiasını sahtelerdi).
Kullanım: 154_taban_taramasi.py <taban>   (0.46 | 0.52 | 0.58)
"""
import sys, json, time
import numpy as np
from pathlib import Path
import importlib

sys.path.insert(0, str(Path(__file__).resolve().parent))
CEK = importlib.import_module("154_cekirdek")

HERE = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
OUT = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/154")

ORTAK = [(round(0.465 + 0.03 * i, 3), round(0.495 + 0.03 * i, 3))
         for i in range(11)]          # 0.465 … 0.795

if __name__ == "__main__":
    tb = float(sys.argv[1])
    t0 = time.time()
    d = np.load(HERE / "128_odl_zeros6_2e6_zeros.npz")
    Z = np.sort(np.asarray(d["zeros"], dtype=float))
    zz = Z[len(Z) - 300000:]
    bant = [b for b in ORTAK if b[0] >= tb - 1e-9]
    print(f"=== TABAN {tb} / son-300k / {len(bant)} bant "
          f"({bant[0][0]}–{bant[-1][1]}) ===", flush=True)
    r = CEK.olc(zz, f"taban{tb}", bantlar=bant, taban=tb, cap=4000)
    r["taban"] = tb
    r["sure_s"] = time.time() - t0
    (OUT / f"R_taban_{tb}.json").write_text(json.dumps(r, indent=1))
    print(f"\nbitti — {(time.time()-t0)/60:.1f} dk", flush=True)
