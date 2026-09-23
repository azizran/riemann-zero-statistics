# -*- coding: utf-8 -*-
"""
188c — K2: İKİZ (Hkeskin) HARİTASI — H-188d için dilim-profilleri
================================================================
ONKAYIT_188 dondurdu: aynı K/Ĝ/S/W tanımları (188b makinesi AYNEN, import),
Hkeskin kinematiği (b185.kinematik AYNEN), τ' ∈ (0.86, 1.20], dilim 0.01
(34 dilim), alt-örneklem her 3. nokta; karışım = Hkeskin'in kendi
c^kesik − c^öz'ü (187c AYNEN); maskeler gerçek evrenin τ'su (187c AYNEN).

Çıktı: 188/ikiz_dilim_<i>.npz, harita_proj_Hkeskin.npz,
       harita_K_Hkeskin.npz, K1_harita_Hkeskin.json
Kullanım: 188c_ikiz_harita.py [işçi_sayısı]
"""
import importlib.util
import os
import sys
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
           "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
spec = importlib.util.spec_from_file_location(
    "b188", QM / "188_configs" / "188b_harita.py")
b188 = importlib.util.module_from_spec(spec)
sys.modules["b188"] = b188          # fork'lu havuzda işçi fonksiyonu için
spec.loader.exec_module(b188)

if __name__ == "__main__":
    nw = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    ONK = b188.onkayit()
    print("=" * 78)
    print(f"188c / K2 İKİZ HARİTASI (Hkeskin)  [on-kayit {ONK['zaman']} "
          f"sha {ONK['sha256'][:12]}]")
    print("=" * 78, flush=True)
    b188.kos("Hkeskin", "dilim_izgara_ikiz", "ikiz_dilim_", "K1_Hkeskin.npz",
             "OZ_Hkeskin.npz", "G1_proj_Hkeskin.npz", "Hkeskin", nw,
             tam_kontrol=None)
