# -*- coding: utf-8 -*-
"""
177b — DONDURULMUŞ KURALIN GÜÇ GEOMETRİSİ (ölçümden ÖNCE)
==========================================================
ÖN-MÜHÜR: `tohum = 3` vekilinin HİÇBİR niceliği var olmadan koşar
(176b VF3 inşası sürerken; 176c hiç koşmamışken). Yeni ölçüm YOKTUR.

Bu betik yalnızca **177a'nın dondurulmuş kuralının kendi geometrisini**
sorar: 176'nın iki ölçülmüş θ'sı sabitken, üçüncü (ve dördüncü) tohumun
hangi değerleri hükmü KESİN kılar? Kural değişmez; yalnız gücü ölçülür.
Amaç, sonucu görmeden önce "bu merdiven kaç basamak sürer?" sorusunu
kayda geçirmektir.

Çıktı: 177/GUC_GEOMETRISI.json
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
TH_A = K["F3"]["dal_a_yasar"]
TH_B = K["F3"]["dal_b_olur"]
TH_HK = K["F3"]["capa"]["Hkeskin"]
KIL = K["H2"]["kilit_th"]
Y12 = [K["F3"]["olculen176"]["VF1"], K["F3"]["olculen176"]["VF2"]]


def sac(y):
    y = np.asarray(y, float)
    return float(2.0 * np.std(y, ddof=1) / math.sqrt(len(y)))


def kesin_a(ek):
    y = Y12 + list(ek)
    o = float(np.mean(y))
    return bool(o <= TH_A and (TH_A - o) >= sac(y))


def kesin_b(ek):
    y = Y12 + list(ek)
    o = float(np.mean(y))
    return bool(o >= TH_B and (o - TH_B) >= sac(y))


def kesin_H2(ek):
    y = Y12 + list(ek)
    d = [math.log(TH_HK) - math.log(x) for x in y]
    db = float(np.mean(d))
    return bool(db >= KIL and (db - KIL) >= sac(d))


def araliklar(mask, g):
    out, i = [], 0
    while i < len(mask):
        if mask[i]:
            j = i
            while j + 1 < len(mask) and mask[j + 1]:
                j += 1
            out.append([float(g[i]), float(g[j])])
            i = j + 1
        else:
            i += 1
    return out


print("=" * 78)
print("177b — DONDURULMUŞ KURALIN GÜÇ GEOMETRİSİ (ölçümden önce)")
print("=" * 78)
print("  ön-kayıt %s  sha %s" % (ONK["zaman"], ONK["sha256"][:16]))
print("  sabit: θ(VF1) = %.7f , θ(VF2) = %.7f" % tuple(Y12))

g = np.linspace(0.55, 0.98, 86001)
a3 = np.array([kesin_a((v,)) for v in g])
b3 = np.array([kesin_b((v,)) for v in g])
h3 = np.array([kesin_H2((v,)) for v in g])
A3, B3, H3 = araliklar(a3, g), araliklar(b3, g), araliklar(h3, g)
print("\n  n = 3 (θ₃ tek bilinmeyen, 0.55–0.98 taraması):")
print("    dal (a) KESİN olan θ₃ kümesi : %s" % (A3 or "BOŞ"))
print("    dal (b) KESİN olan θ₃ kümesi : %s" % (B3 or "BOŞ"))
print("    H2 (ω_θ ≤ 1) KESİN kümesi    : %s" % (H3 or "BOŞ"))

c = np.linspace(0.55, 0.98, 8601)
a4 = np.array([kesin_a((v, v)) for v in c])
h4 = np.array([kesin_H2((v, v)) for v in c])
A4, H4 = araliklar(a4, c), araliklar(h4, c)
print("\n  n = 4, θ₃ = θ₄ = c kesitinde:")
print("    dal (a) KESİN c aralığı : %s" % (A4 or "BOŞ"))
print("    H2 KESİN c aralığı      : %s" % (H4 or "BOŞ"))

gg = np.linspace(0.70, 0.95, 501)
cnt_a = int(sum(kesin_a((x, y)) for x in gg for y in gg))
cnt_h = int(sum(kesin_H2((x, y)) for x in gg for y in gg))
print("    (0.70–0.95)² ızgarasında KESİN kesir: dal(a) %%%.1f  H2 %%%.1f"
      % (100.0 * cnt_a / gg.size ** 2, 100.0 * cnt_h / gg.size ** 2))

sonuc = ("n = 3 HİÇBİR θ₃ için kesinleşemez ⇒ merdiven zorunlu olarak "
         "n = 4'e gider." if not (A3 or B3) else
         "n = 3 bazı θ₃ değerleri için kesinleşebilir.")
print("\n  ⇒ %s" % sonuc)
print("  (Bu bir ölçüm değil, dondurulmuş kuralın cebridir; θ₃ henüz "
      "YOKTUR.)")

p = S177 / "GUC_GEOMETRISI.json"
p.write_text(json.dumps(dict(
    onkayit=ONK["zaman"], sha=ONK["sha256"], y12=Y12,
    esikler=dict(a=TH_A, b=TH_B, kilit_th=KIL, th_Hk=TH_HK),
    n3=dict(dal_a=A3, dal_b=B3, H2=H3),
    n4_kesit=dict(dal_a=A4, H2=H4),
    n4_izgara=dict(dal_a_kesir=cnt_a / gg.size ** 2,
                   H2_kesir=cnt_h / gg.size ** 2),
    sonuc=sonuc), indent=1, ensure_ascii=False))
print("  -> %s" % p)
