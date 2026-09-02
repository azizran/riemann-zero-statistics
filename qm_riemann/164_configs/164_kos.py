"""
164 — ÖLÇÜM SÜRÜCÜSÜ: 158'in zinciri, sadakatli gazlar için
===========================================================
158_kos.py'nin AYNISI; tek farkı çıktının `scratchpad/164`'e yazılması
(158'in kaydı KİRLETİLMİYOR). Ölçüm parçası KOPYALANMADI:

    155_cekirdek.olc155   ← bant döngüsü (çizgi kaydı + τ_eff + R + n_eff)
    155_kos.veri_yukle    ← gaz yükleyici (scratchpad/155/z_<ad>.npy)
    158_kos.IZGARA158     ← izgara(0.28, 0.64, 0.02)

Kullanım:  164_kos.py <veri> <taban>
"""
import importlib
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, "/Users/ugursezen/Desktop/arin/deney/qm_riemann/155_configs")
sys.path.insert(0, "/Users/ugursezen/Desktop/arin/deney/qm_riemann/158_configs")
K = importlib.import_module("155_cekirdek")
KOS155 = importlib.import_module("155_kos")
K158 = importlib.import_module("158_kos")

SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/164")

if __name__ == "__main__":
    veri, taban = sys.argv[1], float(sys.argv[2])
    t0 = time.time()
    SCR.mkdir(parents=True, exist_ok=True)
    z = KOS155.veri_yukle(veri)
    bant = [b for b in K158.IZGARA158 if b[0] >= taban - 1e-9]
    elenen = len(K158.IZGARA158) - len(bant)
    print(f"=== 164 / veri={veri} / taban={taban} / ızgara 0.02 (158) ===")
    print(f"    n={len(z)}  t∈[{z[0]:.1f},{z[-1]:.1f}]  "
          f"{len(bant)} bant ({bant[0][0]}–{bant[-1][1]}); "
          f"{elenen} bant taban altında elendi", flush=True)
    r = K.olc155(z, f"{veri}-t{taban}-164", bant,
                 anahtar=veri, taban=taban, cap=4000)
    r["veri"] = veri
    r["bantset"] = "158"
    r["t_lo"] = float(z[0])
    r["t_hi"] = float(z[-1])
    r["sure_s"] = time.time() - t0
    p = SCR / f"F_{veri}_t{taban}.json"
    p.write_text(json.dumps(r, indent=1))
    print(f"\n-> {p}   ({(time.time()-t0)/60:.1f} dk)", flush=True)
