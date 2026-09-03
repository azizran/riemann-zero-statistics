"""
171g — T2c ÖLÇÜMÜ: λ = 0.50 (L050) ve λ = 1.30 (L130) gazlarında `c`
=====================================================================
ÖLÇÜM ZİNCİRİ KOPYALANMAZ: `167_olcum.kos` AYNEN çağrılır — 170'in
L060/L115 ölçümleriyle birebir aynı yol (`kos(<gaz>, 0.40, 0.95, 0, 0)`).
Tek yapılan `167_ortak.KUNYE`'ye künyeyi eklemektir.

ÖN KAYITLAR bu koşulardan ÖNCE yazıldı:
   scratchpad/171/ONKAYIT_L050.json   (2026-09-04 00:00:05)
   scratchpad/171/ONKAYIT_L130.json   (2026-09-04 00:03:56)
   L050:  M0 0.3098 | M2 0.3541 | M3 0.3515 | M8 0.3504 | M1 0.3470 |
          M7 0.3598 | M9 0.3791
   L130:  M0/M2/M3 0.4618 | M8 0.4654 | M1 0.4339 | M7 0.4467 | M9 0.4699

Kullanım: 171g_olcum.py <L050|L130>
Çıktı:    scratchpad/167/C_<gaz>.json  (167'nin kendi dosya düzeni)

SONUÇ (gerçek koşudan, 3.8 dk her biri — 167'nin kendi üyesi c_ampX ile):
 c_ampX(L050) = 0.4369 (log-sd 0.0261, 6 bant, τ-eğimi −0.307)
 c_ampX(L130) = 0.7624 (log-sd 0.2049, 6 bant, τ-eğimi +2.730)
 R_bant: L050 1.50–1.82, L130 1.05–1.16 (aile eğiliminin düzgün devamı;
 Hkeskin 1.22–1.40, L060 1.47–1.74, L115 1.14–1.28).
 HÜKÜM ÜYESİ c_WX = KALİB_u2/W_X ve 5-bant penceresi 171h'de kurulur:
 c(L050) = 0.3800 ± 0.0072, c(L130) = 0.4540 ± 0.0136.
"""
import importlib
import sys
from pathlib import Path

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
for _p in ("167_configs", "166_configs", "165_configs", "163_configs",
           "160_configs", "159_configs", "155_configs", "154_configs"):
    sys.path.insert(0, str(QM / _p))
ORT = importlib.import_module("167_ortak")
ORT.KUNYE["L050"] = dict(gercek=False, lam=0.50, tau_ust=1.00, pen=None,
                         aile="lam", T=1.0)
ORT.KUNYE["L130"] = dict(gercek=False, lam=1.30, tau_ust=1.00, pen=None,
                         aile="lam", T=1.0)
O167 = importlib.import_module("167_olcum")

if __name__ == "__main__":
    ad = sys.argv[1]
    if ad not in ("L050", "L130"):
        raise SystemExit("kullanım: 171g_olcum.py <L050|L130>")
    O167.kos(ad, 0.40, 0.95, 0, 0)
