# -*- coding: utf-8 -*-
"""
182c — K1 (ÖLÇÜM): VF5..VF8'de 167 ZİNCİRİ + 172b DEFTERİ
==========================================================
ÖN-MÜHÜR. Yeni makine yazılmaz:
  · ölçüm  = `176c_olcum.kos(ad)` AYNEN (o da `167_olcum.kos(ad, 0.40,
    0.95, kule=0, düz=0)`'i aynen çağırır — VF1..VF4 ile BİREBİR aynı
    çağrı, aynı taban/τ_c/kule/düz);
  · defter satırı = `180a_onkayit.defter(ad)` AYNEN (172b'nin
    `oku`/`kusur`'u) — böylece satır, ON KAYIT'taki VF1..VF4
    satırlarıyla ve G_VS*.json ile ALAN ALAN aynı biçimdedir.

Çıktı: 182/G_<ad>.json  (180c/ONKAYIT biçimi: M, Q_E, rho_E, Q_X, rho_X,
g_E, g_X, theta, g_cal, …) + 176c'nin F0 R_bant kapısı da raporlanır.

Kullanım: 182c_olcum.py VF5 VF6 VF7 VF8
"""
import importlib
import json
import sys
import time
from pathlib import Path

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
for _p in ("182_configs", "180_configs", "176_configs", "172_configs",
           "167_configs", "164_configs"):
    sys.path.insert(0, str(QM / _p))
O176C = importlib.import_module("176c_olcum")
A180 = importlib.import_module("180a_onkayit")

SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S182, S176 = SCR / "182", SCR / "176"


def kos(ad):
    t0 = time.time()
    rec176 = O176C.kos(ad)                 # 167 zinciri + F0 kapısı
    row = A180.defter(ad)                  # HAM + ORAN defter satırı
    row["ad"] = ad
    row["F0_R_bant"] = rec176["F0_R_bant"]
    row["R_bant_min"] = rec176["R_bant_min"]
    row["R_bant_hukum"] = rec176["R_bant_hukum"]
    row["sigds"] = rec176["sigds"]
    row["sigC"] = rec176["sigC"]
    row["sigX"] = rec176["sigX"]
    row["T"] = rec176["T"]
    S182.mkdir(parents=True, exist_ok=True)
    p = S182 / f"G_{ad}.json"
    p.write_text(json.dumps(row, indent=1, ensure_ascii=False, default=float))
    print(f"\n  HAM DEFTER ({ad}): M={row['M']:.6f}  Q_E={row['Q_E']:.6f}  "
          f"ρ_E={row['rho_E']:.6f}  Q_X={row['Q_X']:.6f}  "
          f"ρ_X={row['rho_X']:.6f}")
    print(f"  EK: g_E={row['g_E']:.6f}  g_X={row['g_X']:.6f}  "
          f"θ={row['theta']:.6f}  g_cal={row['g_cal']:.6f}")
    print(f"  -> {p}   ({(time.time()-t0)/60:.1f} dk)", flush=True)
    return row


if __name__ == "__main__":
    for g in (sys.argv[1:] or ["VF5"]):
        kos(g)
