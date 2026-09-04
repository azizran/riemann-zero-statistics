# -*- coding: utf-8 -*-
"""
177d — `n = 4` PENCERESİ (ölçülmüş tohumlar verildiğinde)
==========================================================
ÖLÇÜMDEN SONRA yazıldı (176'nın `176g_figur.py`'si gibi). **Kuralı
değiştirmez**; yalnız 177a'da dondurulmuş kuralın geometrisini,
ölçülmüş θ'lar verildiğinde yeniden raporlar. Yeni ölçüm yoktur.

Koşma anı: `VF3` ölçülmüş, `VF4` **henüz kurulmakta** (θ₄ YOKTUR).
Sorulan: dondurulmuş kural, hangi θ₄ değerleri için `n = 4`'te
kesinleşir?

Çıktı: 177/PENCERE_n4.json
"""
import json
import math
from pathlib import Path

import numpy as np

SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S176, S177 = SCR / "176", SCR / "177"

ONK = json.load(open(S177 / "ONKAYIT_177.json"))
K = ONK["kural"]
TH_A, TH_B = K["F3"]["dal_a_yasar"], K["F3"]["dal_b_olur"]
TH_HK, KIL = K["F3"]["capa"]["Hkeskin"], K["H2"]["kilit_th"]
VEK = ["VF1", "VF2", "VF3"]
Y = [json.load(open(S176 / f"G_{v}.json"))["th"] for v in VEK]


def sac(y):
    y = np.asarray(y, float)
    return float(2.0 * np.std(y, ddof=1) / math.sqrt(len(y)))


def kes_a(v):
    y = Y + [v]
    o = float(np.mean(y))
    return bool(o <= TH_A and (TH_A - o) >= sac(y))


def kes_b(v):
    y = Y + [v]
    o = float(np.mean(y))
    return bool(o >= TH_B and (o - TH_B) >= sac(y))


def kes_H2(v):
    y = Y + [v]
    d = [math.log(TH_HK) - math.log(x) for x in y]
    db = float(np.mean(d))
    return bool(db >= KIL and (db - KIL) >= sac(d))


def seg(m, g):
    out, i = [], 0
    while i < len(m):
        if m[i]:
            j = i
            while j + 1 < len(m) and m[j + 1]:
                j += 1
            out.append([float(g[i]), float(g[j])])
            i = j + 1
        else:
            i += 1
    return out


print("=" * 78)
print("177d — n = 4 PENCERESİ   [ön-kayıt %s]" % ONK["zaman"])
print("=" * 78)
print("  ölçülmüş: %s"
      % "  ".join("%s %.7f" % (v, t) for v, t in zip(VEK, Y)))
g = np.linspace(0.40, 0.99, 118001)
A = seg(np.array([kes_a(v) for v in g]), g)
B = seg(np.array([kes_b(v) for v in g]), g)
H = seg(np.array([kes_H2(v) for v in g]), g)
print("  n = 4'te dal (a) KESİN olan θ₄ kümesi : %s" % (A or "BOŞ"))
print("  n = 4'te dal (b) KESİN olan θ₄ kümesi : %s" % (B or "BOŞ"))
print("  n = 4'te H2 (ω_θ ≤ 1) KESİN kümesi    : %s" % (H or "BOŞ"))
bos = not (A or B)
print("\n  ⇒ %s" % ("HİÇBİR θ₄ değeri n = 4'te H1'i kesinleştiremez; "
                    "ön-kayıtlı üst sınır n = 4 olduğu için hüküm "
                    "KALICI HÜKÜMSÜZ olacaktır (eşik gevşetme YOK)."
                    if bos else
                    "n = 4'te kesin bölge vardır; θ₄ karar verir."))
print("  (Bu bir ölçüm değil, dondurulmuş kuralın cebridir; θ₄ HENÜZ "
      "YOKTUR — VF4 kurulmakta.)")

p = S177 / "PENCERE_n4.json"
p.write_text(json.dumps(dict(zaman_notu="177d ölçümden sonra yazıldı; "
                             "kural değişmedi. VF4 kurulurken koştu, "
                             "θ₄ yokken.",
                             olculen={v: t for v, t in zip(VEK, Y)},
                             dal_a=A, dal_b=B, H2=H, bos=bos),
                        indent=1, ensure_ascii=False))
print("  -> %s" % p)
