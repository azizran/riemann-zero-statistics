"""
170 — EK KAPI (K2') ÖLÇÜMÜ: L115 (λ = 1.15) gazında `c`
========================================================
`167_olcum.kos` AYNEN çağrılır (`167_olcum.py L115 0.40 0.95 0 0` ile
birebir); L115 künyesi 167_ortak'ta zaten tanımlı.

ÖN KAYIT bu koşudan ÖNCE yazıldı: scratchpad/170/ONKAYIT_L115.json

Kullanım: 170k_olcum115.py
Çıktı:    scratchpad/167/C_L115.json
"""
import importlib
import sys
from pathlib import Path

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
for _p in ("167_configs", "166_configs", "165_configs", "163_configs",
           "160_configs", "159_configs", "155_configs", "154_configs"):
    sys.path.insert(0, str(QM / _p))
O167 = importlib.import_module("167_olcum")

if __name__ == "__main__":
    O167.kos("L115", 0.40, 0.95, 0, 0)
