"""
164 — ×1.4 FAZ AÇIĞININ YENİDEN YARGILANMASI (152'nin ÖZGÜN ÖLÇÜSÜ)
===================================================================
152'nin "faz oranı" tam olarak şudur:

    oran(bant) = arg Γ_rot(sentetik) / arg Γ_rot(gerçek)
    Γ_rot = Σ(cr_on − cr_off)/Σ(pow_on − pow_off),
    cr = zp·conj(zc)·e^{+i·2πW/L}

yani 155/158'in `phi` alanının ta kendisi. 152'nin konvansiyonu:
    BANTLAR = [(0.525,0.55),(0.55,0.62),(0.62,0.70),(0.70,0.78),(0.78,0.85)]
    η zinciri: taban 0.52, cap 720          (152_gaz.py'nin satırları)
    aday çizgiler: pk(e^{0.86L}), 220 örnek, tohum 21, gap ≥ 2.5·dres

Bu sürücü ÖLÇÜMÜ KOPYALAMAZ: 155_cekirdek.olc155'i o konvansiyonla
çağırır. `keskin` ve `son` koşuları 152'nin kayıtlı sayılarını yeniden
üretmek zorundadır (MEŞRUİYET sınavı — §Dürüstlük).

Kullanım:  164_faz.py <veri> [<veri> ...]
"""
import importlib
import json
import sys
import time
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
sys.path.insert(0, str(QM / "155_configs"))
K = importlib.import_module("155_cekirdek")
KOS155 = importlib.import_module("155_kos")

SCR164 = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
              "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/164")

BANTLAR152 = [(0.525, 0.55), (0.55, 0.62), (0.62, 0.70),
              (0.70, 0.78), (0.78, 0.85)]
# 152_gaz.py'nin GERCEK sözlüğü (Re, Im, R) — TABLODAN-OKUNAN referans,
# yalnız meşruiyet sınavı için; bu koşu `son`'u kendisi ölçüyor.
GERCEK152 = {0.5375: (0.787, 0.206, 0.40), 0.585: (0.550, 0.488, 1.07),
             0.660: (0.118, 0.614, 1.85), 0.740: (-0.148, 0.498, 2.12),
             0.815: (-0.103, 0.655, 1.50)}
# 152'nin kayıtlı keskin/A4 satırları (scratchpad/152/ozet_*.json'dan)
REF152 = {}
for _ad in ("keskin", "A4"):
    _p = SCR164.parent / "152" / f"ozet_{_ad}.json"
    if _p.exists():
        REF152[_ad] = {round(b["tau"], 4): b
                       for b in json.load(open(_p))["bantlar"]}


def kos(veri):
    tb = time.time()
    z = KOS155.veri_yukle(veri)
    print(f"=== 164-FAZ {veri} (152 konvansiyonu: taban 0.52, cap 720) ===",
          flush=True)
    r = K.olc155(z, f"{veri}-152konv", BANTLAR152, anahtar=veri,
                 taban=0.52, cap=720)
    r["veri"] = veri
    r["sure_s"] = time.time() - tb
    SCR164.mkdir(parents=True, exist_ok=True)
    (SCR164 / f"faz152_{veri}.json").write_text(json.dumps(r, indent=1))
    # meşruiyet: 152'nin kayıtlı Γ'sıyla karşılaştır
    ref = REF152.get(veri)
    if ref:
        print("  [MEŞRUİYET] 152'nin kayıtlı Γ'sıyla fark:", flush=True)
        for b in r["bantlar"]:
            if not b.get("olculdu"):
                continue
            q = ref.get(round(b["tau"], 4))
            if q and q.get("re") is not None:
                print(f"    τ={b['tau']:.4f}  ΔRe={b['Gre']-q['re']:+.2e}  "
                      f"ΔIm={b['Gim']-q['im']:+.2e}", flush=True)
    if veri in ("son",):
        print("  [MEŞRUİYET] 152'nin GERCEK satırıyla fark:", flush=True)
        for b in r["bantlar"]:
            g = GERCEK152.get(round(b["tau"], 4))
            if g and b.get("olculdu"):
                print(f"    τ={b['tau']:.4f}  ΔRe={b['Gre']-g[0]:+.2e}  "
                      f"ΔIm={b['Gim']-g[1]:+.2e}  ΔR={b['R']-g[2]:+.3f}",
                      flush=True)
    print(f"-> faz152_{veri}.json  ({(time.time()-tb)/60:.1f} dk)\n",
          flush=True)
    return r


if __name__ == "__main__":
    for v in sys.argv[1:]:
        kos(v)
