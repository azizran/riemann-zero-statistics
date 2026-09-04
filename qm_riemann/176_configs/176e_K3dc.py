# -*- coding: utf-8 -*-
"""
176e — K3 (DC KAÇAĞI): VEKİL GAZDA 174d'NİN MAKİNESİ
======================================================
ÖN-MÜHÜR (koşudan önce yazıldı).

`174d_K3_dc.kos` AYNEN import edilir (özdeşlik
`ort(E) = Σ_q Re[hp_q κ(ω_q)]`, `ξ_q = |hp||κ|cosΔφ`). İki teknik
uyarlama, ikisi de ölçüme dokunmuyor:
  (1) çıktı dizini `176/`;
  (2) `174d` referans μ̂²'yi `172/G1.json`'dan okur; vekil orada
      yoktur, bu yüzden `176c`'nin yazdığı `176/G_<ad>.json` satırı
      modülün G1 sözlüğüne ENJEKTE edilir. Bu, yalnızca Ö2 denetiminin
      (bağımsız yol: Σξ ↔ μ̂²) referansıdır; ölçülen hiçbir sayıyı
      değiştirmez.

SINANAN ÖN-KAYIT MADDESİ: F6 (T-3) — "DC kaçağı faz kilidine KÖRDÜR":
ort(E)(vekil) > 0 ve |ort(E)/0.0621390 − 1| ≤ 0.30; ⟨cosΔφ⟩ 0.50–0.80
bantlarında ≥ 0.90. Tutmazsa T-3 ÖLÜR ve öyle yazılır.

Kullanım: 176e_K3dc.py <ad> [<ad> ...]
"""
import importlib
import json
import sys
from pathlib import Path

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
for _p in ("176_configs", "174_configs"):
    sys.path.insert(0, str(QM / _p))
D174 = importlib.import_module("174d_K3_dc")

SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S176 = SCR / "176"
D174.S174 = S176                # çıktılar 176/ altına


if __name__ == "__main__":
    for g in (sys.argv[1:] or ["VF1"]):
        D174.G1[g] = json.load(open(S176 / f"G_{g}.json"))
        D174.kos(g)
