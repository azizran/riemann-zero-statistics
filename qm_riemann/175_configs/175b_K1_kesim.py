# -*- coding: utf-8 -*-
"""
175b — K1: GİRİŞİM ORANLARI, KESİM AİLESİNDE (162 makinesi)
============================================================
174b'nin makinesi (= 162'nin `Taban162`'si, taban 0.40, cap 4000,
τ_çizgi 0.86, NITER 0) HİÇ DEĞİŞTİRİLMEDEN kesim ailesine uygulanır.
Tek satır bile kopyalanmaz: `174b_K1_girisim` import edilir, `kos()`
çağrılır, çıktı yine `scratchpad/174/K1_<gaz>.json`'a yazılır (aynı
makine, aynı şema; 175 onu okur).

ÖN-MÜHÜR (175b, koşudan önce — 175a'nın ONKAYIT_K2.json'uyla tutarlı):
  Ö-b1  Dört kesim gazında da L = 12.029593242, nline = 3425 ve
        m-kimliği 0.0e+00 çıkmalı (162 makinesi bit düzeyinde aynı).
  Ö-b2  Bant defteri ÖZDEŞ kapanmalı: bağıl kalıntı ≤ 1e−12 (Ö4, 174b).
  Ö-b3  R_η = π_E × r̄ öngörüsü (175a-P1): |R_η/π_E − 1| ≤ %2 dört gazda.
  Ö-b4  R_η(son) = 1.28829 kesim ailesinin [min, max] aralığında OLMALI
        (175a-P2). Değilse H-K1'in η-kanalı ölür.

Kullanım: 175b_K1_kesim.py <gaz> [<gaz> ...]
"""
import importlib
import json
import math
import sys
from pathlib import Path

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
sys.path.insert(0, str(QM / "174_configs"))
B174 = importlib.import_module("174b_K1_girisim")

SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
G1 = json.load(open(SCR / "172" / "G1.json"))
ONK = json.load(open(SCR / "175" / "ONKAYIT_K2.json"))

print("=" * 74)
print("175b — K1 kesim ailesi   [ön-kayıt %s  sha %s]"
      % (ONK["zaman"], ONK["sha256"][:16]))
print("=" * 74)
for g in (sys.argv[1:] or ["HA4"]):
    d = B174.kos(g)
    piE = G1[g]["E"]["pi"]
    R = d["eta"]["R"]
    ong = ONK["kural"]["P1"]["ongoru"].get(g, float("nan"))
    print("  [ÖN-MÜHÜR DENETİMİ %s]" % g)
    print("    Ö-b1  L=%.9f  nline=%d  m-kimlik=%.1e   %s"
          % (d["L"], d["nline"], d["m_kimlik"],
             "✓" if (abs(d["L"] - 12.029593242) < 1e-8 and d["nline"] == 3425)
             else "✗"))
    print("    Ö-b2  defter bağıl kalıntı: η %.1e  Ĉ %.1e  X̃ %.1e   %s"
          % (d["eta"]["defter_bagil"], d["Chat"]["defter_bagil"],
             d["Xtil"]["defter_bagil"],
             "✓" if max(d["eta"]["defter_bagil"], d["Chat"]["defter_bagil"],
                        d["Xtil"]["defter_bagil"]) <= 1e-12 else "✗"))
    print("    Ö-b3  R_η = %.5f   π_E = %.5f   fark %+.2f%%   "
          "(ön-kayıt öngörüsü %.5f, sapma %+.2f%%)   %s"
          % (R, piE, 100 * (R / piE - 1), ong, 100 * (R / ong - 1),
             "✓" if abs(100 * (R / piE - 1)) <= 2.0 else "✗ KÖPRÜ KIRIK"))
    print("=" * 74, flush=True)
