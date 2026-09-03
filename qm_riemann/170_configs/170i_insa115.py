"""
170 — L115 GAZININ İNŞASI (λ = 1.15): λ ekseninin DERİN ucu
===========================================================
167'nin `KONFIG`'inde L115 zaten tanımlıdır ama inşası 167'de yarım
kalmıştı (bkz. `log_insa_L115.txt`, 3 Eylül 10:45). Burada AYNEN
`167_insa.main("L115")` çağrılır (kod kopyalanmaz).

NEDEN: K2'nin ölçtüğü `c(λ)` λ ≤ 0.70'te bir TABANA oturuyor
(0.3690 / 0.3689) ve λ = 1.00'de 0.4035'e çıkıyor. İki okuma var:
  (Y1) `c` λ ile ARTIYOR ve derin (λ→büyük) doyum limitine, yani
       4/π² = 0.40528'e ALTTAN yaklaşıyor  ⇒  c(1.15) ∈ (0.4035, 0.4053]
  (Y2) H-D1'in T-yasası: σ_Ĉ λ ile büyür ⇒ Π_b T(σ_b) KÜÇÜLÜR ⇒
       c(1.15) < 0.4035
λ = 1.15 bu ikisini AYIRIR (aynı λ ekseninin 1.00'ın ÖTESİNDEKİ ucu).

ÖN-MÜHÜR (koşudan ÖNCE):
  merdiven (kalemle, λ-doğrusal): rms S' = 2.1781, Σa_qω_q = 298.8
  σ_ds ∝ λ^p, p ≈ 0.50–0.58 ⇒ σ_ds(1.15) ≈ 0.4313·1.15^0.55 = **0.4646**
  σ_X̃(1.15) ≈ 0.5601·σ_ds = **0.2602**;  σ_Ĉ(1.15) ≈ 0.2990
  ızgarada ΔG<0 kesri: 0.1915 (λ=1) → **≈0.24** beklenir
  sıralılık TAM, maks|F| < 1e−8, ilk-kök hücreleri benzersiz 300000.
  (λ=1.15'te rms S' = 2.178 > N̄' = 1.9147: F artık ORTALAMA olarak da
   monoton değil — inşa daha zorlu, ama İLK-KÖK tanımı sıralılığı
   ÖZDEŞ garanti eder.)

SONUÇ bloğu YALNIZ gerçek koşu çıktısındandır (log_insa_L115.txt).
Kullanım: 170i_insa115.py [h] [nz] [nwork]
"""
import importlib
import sys
from pathlib import Path

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
sys.path.insert(0, str(QM / "167_configs"))
sys.path.insert(0, str(QM / "164_configs"))
I167 = importlib.import_module("167_insa")
I164 = importlib.import_module("164_insa")

if __name__ == "__main__":
    I167.main("L115",
              float(sys.argv[1]) if len(sys.argv) > 1 else I164.HIZGARA,
              int(sys.argv[2]) if len(sys.argv) > 2 else I164.NZERO,
              int(sys.argv[3]) if len(sys.argv) > 3 else I164.NWORK)
