"""
154 — GERÇEK GAZDA R(τ): İKİ PENCERE (L-değişmezlik + hata)
==========================================================
Kullanım: 154_gercek.py <son|orta> [std|ince]
  son  : zeros6'nın son 300k'sı  (150'nin penceresi, L≈12.03)
  orta : Z[850000:1150000]       (bağımsız pencere, L≈11.47)
Bantlar τ=ω/L ile tanımlı ⇒ iki pencerede τ̄ AYNI ama çizgi kümesi ve L
farklı: R(τ) gerçekten τ'nun fonksiyonu mu, yoksa L'ye mi asılı — testi.
"""
import sys, json, time
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import importlib
CEK = importlib.import_module("154_cekirdek")

HERE = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
OUT = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/154")

PENCERE = {"son": (-300000, None), "orta": (850000, 1150000)}

if __name__ == "__main__":
    ad = sys.argv[1]
    bset = sys.argv[2] if len(sys.argv) > 2 else "std"
    t0 = time.time()
    d = np.load(HERE / "128_odl_zeros6_2e6_zeros.npz")
    Z = np.sort(np.asarray(d["zeros"], dtype=float))
    a, b = PENCERE[ad]
    zz = Z[len(Z) + a:] if b is None else Z[a:b]
    print(f"=== GERÇEK / pencere '{ad}' / bant '{bset}' ===")
    print(f"    t ∈ [{zz[0]:.1f}, {zz[-1]:.1f}]  n={len(zz)}", flush=True)
    bant = CEK.BANTLAR_STD if bset == "std" else CEK.BANTLAR_INCE
    r = CEK.olc(zz, f"gercek-{ad}-{bset}", bantlar=bant)
    r["pencere"] = ad; r["bantset"] = bset
    r["t_lo"] = float(zz[0]); r["t_hi"] = float(zz[-1])
    r["sure_s"] = time.time() - t0
    (OUT / f"R_gercek_{ad}_{bset}.json").write_text(json.dumps(r, indent=1))
    print(f"\nbitti — {(time.time()-t0)/60:.1f} dk", flush=True)
