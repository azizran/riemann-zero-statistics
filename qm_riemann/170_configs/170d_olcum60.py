"""
170 — K2 ÖLÇÜMÜ: L060 (λ=0.60) gazında `c`
==========================================
ÖLÇÜM ZİNCİRİ KOPYALANMAZ: `167_olcum.kos` AYNEN çağrılır (169'un L070
ölçümüyle birebir aynı yol: `167_olcum.py <gaz> 0.40 0.95 0 0`). Tek
yapılan `167_ortak.KUNYE`'ye L060 künyesini eklemektir.

ÖN KAYIT bu koşudan ÖNCE yazıldı:
    scratchpad/170/ONKAYIT_L060.json   (3 Eylül 2026, 20:48:26)
    P1a=0.8865  P1b=0.4645  P2=0.4035  P3=0.3657  P4=0.3757  P5=0.3292

Kullanım: 170d_olcum60.py
Çıktı:    scratchpad/167/C_L060.json  (167'nin kendi dosya düzeni)
"""
import importlib
import sys
from pathlib import Path

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
for _p in ("167_configs", "166_configs", "165_configs", "163_configs",
           "160_configs", "159_configs", "155_configs", "154_configs"):
    sys.path.insert(0, str(QM / _p))
ORT = importlib.import_module("167_ortak")
ORT.KUNYE["L060"] = dict(gercek=False, lam=0.60, tau_ust=1.00, pen=None,
                         aile="lam", T=1.0)
O167 = importlib.import_module("167_olcum")

if __name__ == "__main__":
    O167.kos("L060", 0.40, 0.95, 0, 0)
