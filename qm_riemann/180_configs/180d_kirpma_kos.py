# -*- coding: utf-8 -*-
"""
180d — K2 (KOŞU): 3-BACAK KIRPMA AKTARIMI, KARIŞIK VEKİLLERDE
==============================================================
ÖN-MÜHÜR. Hiçbir yeni makine yazılmaz: `169_k2b.main(gaz, 0.40, 0.95)`
AYNEN çağrılır (169/174e'nin kullandığı çağrının birebir aynısı) ve
çıktısı `169/K2b_<gaz>.json`'a düşer. Bu betik yalnız koşturucudur;
HÜKÜM 180e'de, 180a'da dondurulmuş kuralla verilir.

Gazlar: VF1..VF4 — Hkeskin ZARFLI, fazları U(0,2π) karıştırılmış,
sadakatli zincirle sıfırdan kurulmuş vekiller (176/177'de kurulu).

Kullanım: 180d_kirpma_kos.py VF1 VF2 VF3 VF4
"""
import importlib
import sys
import time
from pathlib import Path

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
sys.path.insert(0, str(QM / "169_configs"))
K2B = importlib.import_module("169_k2b")

if __name__ == "__main__":
    t0 = time.time()
    for g in (sys.argv[1:] or ["VF1"]):
        print("\n" + "#" * 74, flush=True)
        print(f"# 180d — 169_k2b.main({g!r}, 0.40, 0.95)", flush=True)
        print("#" * 74, flush=True)
        K2B.main(g, 0.40, 0.95)
    print(f"\n180d bitti ({(time.time()-t0)/60:.1f} dk)", flush=True)
