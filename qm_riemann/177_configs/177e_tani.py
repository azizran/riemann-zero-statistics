# -*- coding: utf-8 -*-
"""
177e — TANI: θ'NIN TOHUM GÜRÜLTÜSÜ NEREDEN GELİYOR? + MALİYET
==============================================================
ÖLÇÜMDEN SONRA yazıldı. **Hiçbir hükmü değiştirmez**, hiçbir eşiğe
dokunmaz, yeni ölçüm yapmaz. Yalnız kayda geçmiş dört defter satırının
aritmetiğidir:

  (1) Özdeşlik `θ ≡ M / (g_E · g_X²)` vekillerde SAYISAL olarak sınanır.
  (2) `log θ`'nın tohumlar arası varyansı, `log M`, `log g_E`,
      `2 log g_X` paylarına ayrıştırılır (bağımsızlık VARSAYIMIYLA,
      yalnız büyüklük sıralaması için).
  (3) Dondurulmuş kuralın (2 standart hata) kesinleşmesi için gereken
      tohum sayısı, ölçülmüş `s`'ten hesaplanır.

Çıktı: 177/TANI.json
"""
import json
import math
from pathlib import Path

import numpy as np

SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S176, S177 = SCR / "176", SCR / "177"
VEK = ["VF1", "VF2", "VF3", "VF4"]
G = {v: json.load(open(S176 / f"G_{v}.json")) for v in VEK}
ONK = json.load(open(S177 / "ONKAYIT_177.json"))
K = ONK["kural"]
TH_A, TH_HK = K["F3"]["dal_a_yasar"], K["F3"]["capa"]["Hkeskin"]
D_SON, KIL = K["F3"]["capa"]["dlog_th_son_Hk"], K["H2"]["kilit_th"]

print("=" * 78)
print("177e — TANI (ölçümden sonra; hüküm DEĞİŞMEZ)")
print("=" * 78)

# (1) özdeşlik
res = [G[v]["th"] - G[v]["KAL"] / (G[v]["gE"] * G[v]["gX"] ** 2)
       for v in VEK]
print("\n  (1) θ ≡ M/(g_E g_X²) kalıntısı: %s"
      % "  ".join("%.1e" % r for r in res))

# (2) varyans ayrışımı
lab = dict(th="log θ", KAL="log M", gE="log g_E", gX="log g_X")
d = {k: np.array([math.log(G[v][k]) for v in VEK]) for k in lab}
var = {k: float(np.var(d[k], ddof=1)) for k in lab}
pay = dict(M=var["KAL"], g_E=var["gE"], g_X2=4.0 * var["gX"])
tot = sum(pay.values())
print("\n  (2) tohumlar arası saçılım (n = 4)")
print("      %-8s %-42s %8s %8s" % ("nicelik", "değerler", "CV%", "pay%"))
for k, ad in (("KAL", "M"), ("gE", "g_E"), ("gX", "g_X"), ("th", "θ")):
    y = np.array([G[v][k] for v in VEK])
    cv = 100.0 * y.std(ddof=1) / y.mean()
    p = ("%8.1f" % (100.0 * pay[{"KAL": "M", "gE": "g_E",
                                 "gX": "g_X2"}[k]] / tot)) if k != "th" \
        else "       —"
    print("      %-8s %-42s %8.2f %s"
          % (ad, " ".join("%.6f" % x for x in y), cv, p))
print("      Var(log θ) ölçülen = %.3e ;  Var(log M)+Var(log g_E)"
      "+4Var(log g_X) = %.3e" % (var["th"], tot))
print("      ⇒ θ'nın tohum gürültüsünün **%%%.0f'i g_X'ten** gelir "
      "(g_E'nin payı %%%.1f)."
      % (100.0 * pay["g_X2"] / tot, 100.0 * pay["g_E"] / tot))

# (3) maliyet
th = np.array([G[v]["th"] for v in VEK])
dl = np.array([math.log(TH_HK) - math.log(x) for x in th])
s_th, s_dl = float(th.std(ddof=1)), float(dl.std(ddof=1))
print("\n  (3) 2-standart-hata kesinliği için gereken tohum sayısı")
for ad, ort, s, esik in (("H1 dal (a): θ̄ ≤ %.7f" % TH_A, float(th.mean()),
                          s_th, TH_A),
                         ("H2 üst kenar: Δ̄ ≥ %.7f" % KIL, float(dl.mean()),
                          s_dl, KIL),
                         ("H2 alt kenar: Δ̄ > 0", float(dl.mean()),
                          s_dl, 0.0)):
    marj = abs(ort - esik)
    n_ger = (2.0 * s / marj) ** 2 if marj > 0 else float("inf")
    print("      %-38s ort %+.7f  s %.7f  marj %.7f  ⇒ n ≳ %s"
          % (ad, ort, s, marj,
             ("%.0f" % math.ceil(n_ger)) if np.isfinite(n_ger) else "∞"))
    if ad.startswith("H1") and ort > esik:
        print("        (ortalama eşiğin YANLIŞ tarafında: bu dal hiçbir n "
              "ile YAŞAyamaz)")

print("\n  (4) kilidin θ'da taşıdığı toplam / gerçeğin fazlası = "
      "%%%.1f  (tohum başına %s)"
      % (100.0 * dl.mean() / D_SON,
         " / ".join("%%%.1f" % (100.0 * x / D_SON) for x in dl)))

p = S177 / "TANI.json"
p.write_text(json.dumps(dict(
    vekiller=VEK, ozdeslik_kalinti=res,
    var_log={k: var[k] for k in var}, pay=pay,
    gX_payi=pay["g_X2"] / tot, gE_payi=pay["g_E"] / tot,
    M_payi=pay["M"] / tot,
    theta=list(map(float, th)), delta=list(map(float, dl)),
    s_theta=s_th, s_delta=s_dl,
    yuzde_fazla=100.0 * float(dl.mean()) / D_SON,
    yuzde_fazla_tohum=[100.0 * float(x) / D_SON for x in dl]),
    indent=1, ensure_ascii=False))
print("  -> %s" % p)
