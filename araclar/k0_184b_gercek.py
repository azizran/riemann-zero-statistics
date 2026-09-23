# -*- coding: utf-8 -*-
"""188-K0 sürücü: 184b_K1_zarf.py'yi YALNIZ 'gercek'(son) gazı için koşar.
Neden: ONKAYIT_184 'ikiz'=keskin (152'nin z_keskin'i) de listeler; 152 önbelleği
temizlendi ve 188 zincirinde kullanılmıyor (zincirin ikizi 184b2 ile Hkeskin).
184b'nin kos() fonksiyonu AYNEN çağrılır — içerik değişikliği yok."""
import importlib.util
from pathlib import Path

QM = Path(__file__).resolve().parents[1] / "qm_riemann"  # taşınabilir
spec = importlib.util.spec_from_file_location("b184", QM / "184_configs" / "184b_K1_zarf.py")
b184 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(b184)
b184.S184.mkdir(parents=True, exist_ok=True)
b184.ONK = b184.onkayit()
print(f"184b (yalnız gercek) [on-kayit {b184.ONK['zaman']} sha {b184.ONK['sha256'][:12]}]", flush=True)
print(b184.kos("gercek", "son"), flush=True)
