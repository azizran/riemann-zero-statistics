# -*- coding: utf-8 -*-
"""
176d — K1 (GİRİŞİM): VEKİL GAZDA 162'NİN MAKİNESİ
===================================================
ÖN-MÜHÜR (koşudan önce yazıldı).

`174b_K1_girisim.kos` AYNEN import edilir (o da `162_cekirdek.Taban162`'yi
aynen kurar: taban 0.40, cap 4000, τ_çizgi 0.86, NITER 0, BLOK 2000).
Tek değişiklik çıktı dizinidir: 176'nın ürünleri `176/` altına yazılır,
174/175'in defteri kirletilmez.

ÖLÇÜLEN: R_η, R_Ĉ, R_X; bant-bant defteri; üçüncü momentler
(m3, m3_çizgi, skew e1/x1/η_çiz/ds) — rastgele fazlı çizgi alanında
ÖZDEŞ SIFIR olan nesneler (176a/F5).

DENETİM (ön-mühür): L, nline = 3425, m-kimliği ve defter kalıntısı
174b'nin sekiz gazındaki değerlerle aynı mertebede olmalı.

Kullanım: 176d_K1.py <ad> [<ad> ...]
"""
import importlib
import sys
from pathlib import Path

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
for _p in ("176_configs", "174_configs"):
    sys.path.insert(0, str(QM / _p))
B174 = importlib.import_module("174b_K1_girisim")

SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
B174.SCR = SCR / "176"          # çıktılar 176/ altına

if __name__ == "__main__":
    for g in (sys.argv[1:] or ["VF1"]):
        B174.kos(g)
