"""
167 — PENCERE EKSENİ: aynı gazın alt-pencereleri (T/2, T/4)
===========================================================
Yeni bir gaz İNŞA EDİLMEZ; var olan tekne dizisinin bir dilimi alınır.
Böylece merdiven, λ ve kesim SABİT kalır, yalnız pencere uzunluğu T
(ve dolayısıyla dres = 2π/T, çizgi çözünürlüğü, N) değişir.

Kullanım: 167_pencere.py            (bütün dilimleri yazar)
"""
import importlib
import sys
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
for _p in ("167_configs", "155_configs"):
    sys.path.insert(0, str(QM / _p))
KOS155 = importlib.import_module("155_kos")
ORT = importlib.import_module("167_ortak")

OUT = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/155")


def main():
    kaynaklar = {}
    for ad, (kay, a, b) in ORT.DILIM.items():
        if kay not in kaynaklar:
            kaynaklar[kay] = np.sort(KOS155.veri_yukle(kay))
        z = kaynaklar[kay]
        n = len(z)
        i0, i1 = int(round(a * n)), int(round(b * n))
        zz = z[i0:i1]
        np.save(OUT / f"z_{ad}.npy", zz)
        print(f"{ad:8s} <- {kay}[{i0}:{i1}]  n={len(zz)}  "
              f"t∈[{zz[0]:.1f},{zz[-1]:.1f}]  T={zz[-1]-zz[0]:.1f}",
              flush=True)


if __name__ == "__main__":
    main()
