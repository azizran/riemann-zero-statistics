"""
170 — L060 GAZININ İNŞASI (λ = 0.60): K2'nin ÖRNEKLEM-DIŞI NOKTASI
==================================================================
Hiçbir çözüm/inşa parçası KOPYALANMAZ: `167_insa.main` aynen çağrılır;
tek yaptığımız `167_insa.KONFIG`'e λ=0.60 satırını EKLEMEK.

    S(t) = −Σ_q (λ a_q w_q) sin(ω_q t) ,  a_q = 1/(π m √q) ,  λ = 0.60
    keskin kesim τ ≤ 1.00  (Hkeskin/L085/L070 ile AYNI merdiven, yalnız
    genlik ölçeği farklı)  ;  seviye konvansiyonu c = −½.

ÖN-MÜHÜR (koşudan ÖNCE yazıldı, 3 Eylül 2026 ~20:30).
  λ-serisinden ölçülmüş marjinaller:
      λ      1.00     0.85     0.70
      σ_ds   0.4313   0.3976   0.3559
      σ_X̃    0.24204  0.22270  0.19897
      σ_Ĉ    0.27768  0.25562  0.22886
  σ_ds λ'da DOĞRUSAL DEĞİL (0.4313·0.70 = 0.3019 olurdu, ölçülen 0.3559).
  İki parametreli uyum σ² = A λ² + B (λ=1.00 ve 0.70'ten):
      A = 0.11639, B = 0.06963  ⇒  σ_ds(0.60) = 0.3390
  Kuvvet yasası σ ∝ λ^p (p = 0.539, aynı iki noktadan):
      σ_ds(0.60) = 0.3277
  ÖNGÖRÜ (koşudan önce):  σ_ds(0.60) ∈ [0.327, 0.340],
                          σ_X̃(0.60) ∈ [0.183, 0.190],
                          rms S' = 1.1364, Σa_qω_q = 155.9 (kalemle: λ·λ=1 değeri)
  Ayrıca: sıralılık TAM, maks|F| < 1e−8, ilk-kök hücreleri benzersiz 300000.

SONUÇ bloğu YALNIZ gerçek koşu çıktısındandır (log_insa_L060.txt).

Kullanım: 170_insa60.py [h] [nz] [nwork]
"""
import importlib
import sys
from pathlib import Path

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
sys.path.insert(0, str(QM / "167_configs"))
sys.path.insert(0, str(QM / "164_configs"))
I167 = importlib.import_module("167_insa")
I164 = importlib.import_module("164_insa")

I167.KONFIG["L060"] = dict(lam=0.60)

if __name__ == "__main__":
    I167.main("L060",
              float(sys.argv[1]) if len(sys.argv) > 1 else I164.HIZGARA,
              int(sys.argv[2]) if len(sys.argv) > 2 else I164.NZERO,
              int(sys.argv[3]) if len(sys.argv) > 3 else I164.NWORK)
