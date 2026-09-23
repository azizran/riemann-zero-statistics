# -*- coding: utf-8 -*-
"""188-K0 yardımcı: gerçek denizin eta_son_t0.4_c4000.npz önbelleğini 155'in
KENDİ üreticisiyle yeniden kurar (155_kos.veri_yukle('son') +
155_cekirdek.eta_onbellek(z, 'son', 0.4, 4000)). İçerik değişikliği yok."""
import importlib.util
import sys
from pathlib import Path

import numpy as np

QM = Path(__file__).resolve().parents[1] / "qm_riemann"  # taşınabilir
spec = importlib.util.spec_from_file_location("k155", QM / "155_configs" / "155_kos.py")
k155 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(k155)
K = k155.K

z = k155.veri_yukle("son")
print(f"son: n={len(z)} t∈[{z[0]:.3f},{z[-1]:.3f}]", flush=True)
C = K.eta_onbellek(z, "son", 0.4, 4000)
mid, ds, L = C["mid"], C["ds"], C["L"]
g = np.diff(z)
g2 = (ds + 1.0) * 2 * np.pi / np.log(mid / (2 * np.pi))
print(f"L={L:.9f} N={len(mid)} nq={C['nq']} var(ds)={np.var(ds):.5f}")
print(f"tutarlılık maks|g − (ds+1)·2π/log(mid/2π)| = {np.max(np.abs(g - g2)):.2e}")
