# -*- coding: utf-8 -*-
"""
192f — ÖN-KAYITSIZ KEŞİF (hükümler görüldükten SONRA; hüküm DIŞI, yalnız KAYIT)
==============================================================================
 (a) A1: sınıf katkılarının toplam yönündeki İMZALI payı f_r = Re(K_r·conj K_top)/|K_top|²
     (+ jk se) ve dik payı g_r = Im(…)/|K_top|²; aday "reel-kısım" deseni
     f_r ∝ cos(2π r b/a) (Σ_r cos = μ(a)) ile yan yana — SINANMADI.
 (b) B: sönmeli konumlarda ±0.05 içindeki EN DÜŞÜK κ (çukur derinliği, κ/se) iki
     pencerede; "söner" hükmünün |κ|<2se'den mi yoksa geniş tabandan mı geldiği.
 (c) A2: 2m / m oranları (P(+log6)/P(+log3), P(+log10)/P(+log5)) ve düşük ↔ son
     oran farkları (jk se'ler bağımsız sayılarak).
Çıktı: 192/K_kesif_192.json + ekran.
"""
import json
from pathlib import Path

import numpy as np

SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S192 = SCR / "192"
NJACK = 8
ONK = json.load(open(S192 / "ONKAYIT_192.json"))
KAT = ONK["B"]["katalog"]
A = np.load(S192 / "K2_sinif.npz")
H = json.load(open(S192 / "HUKUM_192.json"))
PF = np.load(S192 / "profiller_192.npz")


def jk(r):
    r = np.asarray(r, float)
    return float(np.sqrt((NJACK - 1) / NJACK * np.sum((r - r.mean()) ** 2)))


out = {}
# (a)
print("(a) A1 — sınıfların toplam yönündeki imzalı payı f_r (Σ f_r = 1) ve dik pay g_r:")
ad_, sn_ = A["ad"], A["sinif"]
Kt, Kr = A["Kt"], A["Kr"]
bab = {"+log3": (3, 1), "-log3": (3, 3), "+log5": (5, 1), "+log6": (6, 1), "+log2": (4, 1)}
a_ = {}
for uy, (a, b) in bab.items():
    ix = [i for i in range(len(ad_)) if ad_[i] == uy]
    T = Kt[ix].sum()
    Tr = Kr[:, ix].sum(1)
    rs = {}
    cz = {}
    for i in ix:
        s = str(sn_[i])
        if s == "bolunen":
            continue
        r = int(s[1]) if uy != "+log2" else int(s[1])
        f = float((Kt[i] * np.conj(T)).real / abs(T) ** 2)
        fr = (Kr[:, i] * np.conj(Tr)).real / np.abs(Tr) ** 2
        gg = float((Kt[i] * np.conj(T)).imag / abs(T) ** 2)
        gr = (Kr[:, i] * np.conj(Tr)).imag / np.abs(Tr) ** 2
        cz[s] = float(np.cos(2 * np.pi * r * b / a))
        rs[s] = {"f": f, "se_f": jk(fr), "g": gg, "se_g": jk(gr)}
    # cos-deseni normalize: Σ cos ile böl (işaret toplamın yönüne göre)
    sc = sum(cz.values())
    for s in rs:
        rs[s]["cos_deseni"] = cz[s] / sc if abs(sc) > 1e-12 else None
    a_[uy] = rs
    print(f"   {uy:>6} (a={a}, b={b}): " + " | ".join(
        f"{s}: f={v['f']:+.3f}±{v['se_f']:.3f} g={v['g']:+.3f}±{v['se_g']:.3f}"
        f" [cos {v['cos_deseni']:+.3f}]" if v["cos_deseni"] is not None else
        f"{s}: f={v['f']:+.3f}±{v['se_f']:.3f}" for s, v in rs.items()))
out["a_A1_imzali_pay"] = a_

# (b)
print("\n(b) sönmeli konumlarda ±0.05 içindeki en düşük κ (çukur):")
b_ = {}
m = PF["merkez"]
for pad, kk, rr in (("dusuk", PF["d_kap"], PF["d_kap_reps"]),
                    ("son", PF["s_kap"], PF["s_kap_reps"])):
    b_[pad] = {}
    for ad, k in KAT.items():
        if k["liste"] != "sonmeli":
            continue
        h = k["delta_omega"]
        s = np.where(np.abs(m - h) <= 0.05 + 1e-9)[0]
        i = s[np.argmin(kk[s])]
        se = jk(rr[:, i])
        b_[pad][ad] = {"merkez": float(m[i]), "kappa_min": float(kk[i]), "se": se,
                       "kappa_bolu_se": float(kk[i] / se)}
    print(f"   [{pad}] " + "; ".join(f"{a}: {v['merkez']:.3f} {v['kappa_bolu_se']:+.1f}σ"
                                     for a, v in b_[pad].items()))
    # taban oranı: Δω ∈ [−1.3, 2.3] düzgün ızgarasında (0.001) rastgele bir h için
    # ±0.05 içindeki en düşük κ/se ≤ −3σ / −5σ olasılığı (profil geneli)
    se_all = np.sqrt((NJACK - 1) / NJACK * np.sum((rr - rr.mean(0)) ** 2, 0))
    z = kk / se_all
    hh = np.arange(-1.3, 2.3 + 1e-9, 0.001)
    zmin = np.array([np.nanmin(z[np.abs(m - h) <= 0.05 + 1e-9]) for h in hh])
    b_[pad]["TABAN_rastgele_h"] = {"P(zmin<=-3)": float(np.mean(zmin <= -3)),
                                   "P(zmin<=-5)": float(np.mean(zmin <= -5))}
    ns3 = sum(v["kappa_bolu_se"] <= -3 for k_, v in b_[pad].items() if k_ in KAT)
    ns5 = sum(v["kappa_bolu_se"] <= -5 for k_, v in b_[pad].items() if k_ in KAT)
    b_[pad]["sonmeli_sayim"] = {"<=-3σ": int(ns3), "<=-5σ": int(ns5), "n": 9}
    print(f"      taban (rastgele h): P(min≤−3σ) = {np.mean(zmin <= -3):.2f}, "
          f"P(min≤−5σ) = {np.mean(zmin <= -5):.2f} | sönmeli: ≤−3σ {ns3}/9, ≤−5σ {ns5}/9")
out["b_sonmeli_cukur"] = b_

# (c)
print("\n(c) A2 — 2m/m oranları ve düşük ↔ son:")
c_ = {}


def P(ad, kk, rr):
    ix = np.array(KAT[ad]["A2_j"]) - int(PF["J"][0])
    return float(kk[ix].sum()), rr[:, ix].sum(1)


for pad, kk, rr in (("dusuk", PF["d_kap"], PF["d_kap_reps"]),
                    ("son", PF["s_kap"], PF["s_kap_reps"])):
    c_[pad] = {}
    for n2, n1 in (("+log6", "+log3"), ("+log10", "+log5"), ("-log2", "+log2"),
                   ("-log3", "+log3")):
        a1, a1r = P(n2, kk, rr)
        b1, b1r = P(n1, kk, rr)
        c_[pad][f"{n2}/{n1}"] = {"oran": a1 / b1, "se": jk(a1r / b1r)}
    print(f"   [{pad}] " + "; ".join(f"{k}: {v['oran']:.3f}±{v['se']:.3f}"
                                     for k, v in c_[pad].items()))
fk = {}
for k_, v in H["A2"].items():
    if isinstance(v.get("KAYIT_son"), dict):
        d = v["oran"] - v["KAYIT_son"]["oran"]
        s = float(np.hypot(v["se_oran"], v["KAYIT_son"]["se"]))
        fk[k_] = {"dusuk_eksi_son": d, "se": s, "sigma": d / s}
print("   düşük − son (A2 oranları): " + "; ".join(
    f"{k}: {v['dusuk_eksi_son']:+.3f}±{v['se']:.3f} ({v['sigma']:+.1f}σ)" for k, v in fk.items()))
c_["dusuk_eksi_son"] = fk
out["c_A2"] = c_
json.dump(out, open(S192 / "K_kesif_192.json", "w"), indent=1, ensure_ascii=False)
print(f"\n-> {S192/'K_kesif_192.json'}  BİTTİ")
