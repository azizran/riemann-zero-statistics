# -*- coding: utf-8 -*-
"""
191b — DERİN İKİZ τ≤1.30 İNŞASI (189b ve 164_insa DÜZENLENMEDEN).
189b modülü importlib ile yüklenir, AD sözlüğüne "1.30" → "Hderin130"
eklenir, 189b'nin izgara_sinavi / insa_et fonksiyonları AYNEN çağrılır.
Çıktılar 189b'nin yazdığı yerlere gider (scratchpad/189, 164, 155).

Kullanım:
  191b_insa.py izgara 1.30
  191b_insa.py insa 1.30
"""
import importlib.util
import sys
from pathlib import Path

QM = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(QM / "189_configs"))
spec = importlib.util.spec_from_file_location("b189", QM / "189_configs" / "189b_insa.py")
b189 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(b189)
b189.AD["1.30"] = "Hderin130"

if __name__ == "__main__":
    mod, D = sys.argv[1], float(sys.argv[2])
    if f"{D:.2f}" != "1.30":
        raise SystemExit("191b yalnız 1.30 içindir")
    {"izgara": b189.izgara_sinavi, "insa": b189.insa_et}[mod](D)
