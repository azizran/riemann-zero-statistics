# -*- coding: utf-8 -*-
"""
175c — K3: DC KAÇAĞININ ÇİZGİ-ÇİZGİ DEFTERİ, KESİM AİLESİNDE
==============================================================
174d'nin makinesi (özdeşlik `ort(E) = Σ_q Re[hp_q κ(ω_q)]`) hiç
değiştirilmeden kesim ailesine uygulanır: `174d_K3_dc` import edilir.
Çıktı yine `scratchpad/174/K3_<gaz>.json` (aynı şema).

Not (kod denetimi): 174d'de `A_q = λ·a_q·1{τ≤τ_ust}` yazar ve yorum
"çıplak a_q (λ ve kesim HARİÇ)" der; ama `165_cekirdek.merdiven` erfc
penceresini `M["a"]`ya ZATEN uygular (PENCERE sözlüğü 174d'de
`167_ortak.pencere_dict()` ile HA4+E060 için güncellenir). Dolayısıyla
erfc gazlarında A_q = a_q·w_erfc DOĞRU kurulur; keskin gazlarda maske
τ_ust'u uygular. Bu koşuda ikisi de sınanır (Ö-c1).

ÖN-MÜHÜR (175c, koşudan önce; 175a-P3 ile aynı):
  Ö-c1  Σ_q ξ_q'dan çıkan ort(E), 175a-P3'ün özdeşlik öngörüsünü
        (√(μ̂²_E·K_E/π_E)) ‰5 içinde vermeli; İŞARET POZİTİF.
  Ö-c2  Her gazda μ̂²(ξ) ile 172/G1'in μ̂²_E'si ‰5 içinde (174d'nin Ö2'si).
  Ö-c3  τ = 0.5 işaret dönüşü (⟨cosΔφ⟩ < 0 iken τ<0.5, > 0 iken τ>0.5)
        kesim ailesinin HEPSİNDE görülmeli — orta-nokta tarağı gazdan
        bağımsızdır.

Kullanım: 175c_K3_kesim.py <gaz> [<gaz> ...]
"""
import importlib
import json
import math
import sys
from pathlib import Path

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
sys.path.insert(0, str(QM / "174_configs"))
D174 = importlib.import_module("174d_K3_dc")

SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
ONK = json.load(open(SCR / "175" / "ONKAYIT_K2.json"))
P3 = ONK["kural"]["P3"]["ongoru"]

print("=" * 74)
print("175c — K3 kesim ailesi   [ön-kayıt %s  sha %s]"
      % (ONK["zaman"], ONK["sha256"][:16]))
print("=" * 74)
for g in (sys.argv[1:] or ["HA4"]):
    d = D174.kos(g)
    ong = P3.get(g, float("nan"))
    sap = abs(d["ortE"] / ong - 1) if ong == ong else float("nan")
    print("  [ÖN-MÜHÜR DENETİMİ %s]" % g)
    print("    Ö-c1  ort(E) = %+.6f   ön-kayıt öngörüsü %+.6f   sapma %.2e"
          "   %s" % (d["ortE"], ong, sap,
                     "✓" if (d["ortE"] > 0 and sap <= 5e-3) else "✗"))
    print("    Ö-c2  μ̂²(ξ) / μ̂²(G1) − 1 = %.2e   %s"
          % (abs(d["mu2"] / d["mu2_G1"] - 1),
             "✓" if abs(d["mu2"] / d["mu2_G1"] - 1) <= 5e-3 else "✗"))
    alt = [b for b in d["bant"] if b["hi"] <= 0.50 and b["lo"] >= 0.40]
    ust = [b for b in d["bant"] if b["lo"] >= 0.50]
    okA = all(b["cos"] < 0 for b in alt) if alt else False
    okU = all(b["cos"] > 0 for b in ust) if ust else False
    print("    Ö-c3  ⟨cosΔφ⟩ işareti: τ<0.5 %s | τ>0.5 %s   %s"
          % (" ".join("%+.3f" % b["cos"] for b in alt),
             " ".join("%+.3f" % b["cos"] for b in ust),
             "✓" if (okA and okU) else "✗"))
    print("=" * 74, flush=True)
