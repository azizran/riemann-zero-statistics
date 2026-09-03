"""
170 — K0: MUHASEBE TABLOSU (168/169'daki her `c` değerinin KOŞUL DEFTERİ)
=========================================================================
Yeni ölçüm YOK. Girdi yalnız `scratchpad/167/C_<gaz>.json` (168 ve 169 ile
AYNI dosyalar) + `scratchpad/167/insa_<gaz>.json`. Hesap parçası
kopyalanmaz: `169_k1.agg_uyeler / band_jk / saglikli / phi_of / ladder`
(o da `166_T1.bant_agg/_jk`) AYNEN import edilir.

Defterin sütunları — bir `c` sayısını tek başına belirleyen HER koşul:
  * GAZ      : gerçek ζ mi, sentetik mi; λ (genlik ölçeği); kesim
               (keskin τ_ust ya da erfc(τ_c,Δ)); φ = BOŞ ÇİZGİ KESRİ
  * PENCERE  : sıfır dilimi (N, T), bant penceresi (lo aralığı, bant sayısı)
  * ÜYE      : c_WX = KALİB_u2/W_X  (168 §A2.8'in üyesi)
               c_ampX = KALİB_u2/(W_amp·W_X)  (JSON'daki `c_u2`, 167'nin üyesi)
  * AYRIŞTIRMA: KALİB_u2 = g_cal · θ ,  g_cal = g_E g_X²

ÖN-MÜHÜR (koşudan ÖNCE): bu betik YENİDEN ÜRETİM betiğidir; 169 §K1.1/K1.2
tablosunun her satırını dört hanede vermelidir:
   c_WX(5 bant): Hkeskin 0.4035, son 0.4122, L085 0.3817, L070 0.3690,
                 K090 0.3699, K070 0.3480, HA4 0.3492, E060 0.3701
   c_WX(9 bant): Hkeskin 0.4051, son 0.4053, L085 0.3751, K090 0.3709,
                 K070 0.3480, HA4 0.3345, E060 0.3701
   θ/θ₀        : K090 0.9080, K070 0.8038, HA4 0.7975, E060 0.7945
Sapma olursa GİZLENMEZ, tabloda "!" ile işaretlenir.

Çıktı: scratchpad/170/K0.json + ekrana tablo
"""
import importlib
import json
import os
import sys
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
for _p in ("169_configs", "167_configs", "166_configs", "165_configs",
           "163_configs", "160_configs", "159_configs", "155_configs",
           "154_configs"):
    sys.path.insert(0, str(QM / _p))
K1 = importlib.import_module("169_k1")
ORT = importlib.import_module("167_ortak")

SCR = ("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
       "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
SCR167, SCR170 = SCR + "/167", SCR + "/170"
C_HIP = 4.0 / np.pi ** 2

# 169 §K1.1/K1.2 ve 168 §A2.6/A2.7'nin YAZILI sayıları (yeniden üretim sınavı)
REF = {   # gaz: (c_WX 5 bant, c_WX 9 bant, g_cal, theta, c_ampX)
    "Hkeskin": (0.4035, 0.4051, 0.2895, 0.8884, 0.5774),
    "son":     (0.4122, 0.4053, 0.2971, 0.9142, 0.5664),
    "L085":    (0.3817, 0.3751, 0.2912, 0.8959, 0.5171),
    "L070":    (0.3690, None,   None,   None,   None),
    "K090":    (0.3699, 0.3709, 0.2877, 0.8067, 0.5393),
    "K070":    (0.3480, 0.3480, 0.3128, 0.7142, 0.5154),
    "HA4":     (0.3492, 0.3345, 0.3194, 0.7086, 0.5061),
    "E060":    (0.3701, 0.3701, 0.3433, 0.7059, 0.5216),
    "HkT2a":   (0.3871, None, 0.2701, 0.9166, 0.5524),
    "HkT2b":   (0.3708, None, 0.2652, 0.8873, 0.5324),
    "HkT4a":   (0.3309, None, 0.2208, 0.9595, 0.4717),
    "HkT4b":   (0.3314, None, 0.2157, 0.9832, 0.4728),
}
EKSEN = {"Hkeskin": "TABAN", "son": "TABAN (GERÇEK ζ)",
         "L085": "λ", "L070": "λ", "L060": "λ", "L115": "λ",
         "K090": "kesim", "K070": "kesim", "HA4": "kesim/taban", "E060": "kesim",
         "HkT2a": "pencere", "HkT2b": "pencere",
         "HkT4a": "pencere", "HkT4b": "pencere"}


def gaz_satiri(ad, d, tau_l):
    """Bir gazın TAM koşul defteri satırı."""
    B5 = K1.saglikli(d, lo_max=0.68)
    B9 = K1.saglikli(d, lo_max=0.80)
    r = dict(gaz=ad, eksen=EKSEN.get(ad, "?"),
             gercek=bool(d["kunye"]["gercek"]),
             lam=d["lam"], tau_ust=d["tau_ust"], pen=d["pen"],
             N=d["N"], T=d["T"], L=d["L"], taban=d["taban"], tau_c=d["tau_c"],
             dilim=d["kunye"].get("dilim"),
             sigds=d["sigds"], sigX=d["sigX"], sigC=d["sigC"],
             gE=d["artik"]["gE"], gX=d["artik"]["gX"],
             gcal=d["artik"]["gE"] * d["artik"]["gX"] ** 2,
             korE=d["artik"]["korE"], korX=d["artik"]["korX"],
             nline=d["artik"]["nline"],
             phi=K1.phi_of(d, tau_l), nb5=len(B5), nb9=len(B9))
    if not B5:
        return r
    rows = [K1.band_jk(b["cizgi"]) for b in B5]
    c = np.array([a["c_WX"] for a in rows])
    s = np.array([a["sc_WX"] for a in rows])
    lg = np.log(c)
    cg = float(np.exp(lg.mean()))
    s_jk = float(np.sqrt(np.sum((s / c) ** 2)) / len(c)) * cg
    s_sc = (float(np.std(lg, ddof=1) / np.sqrt(len(lg))) * cg
            if len(lg) > 1 else float("nan"))
    kal = float(np.exp(np.mean(np.log([a["KALIB_u2"] for a in rows]))))
    r.update(
        c_WX=cg, s_jk=s_jk, s_bant=s_sc, s_tot=float(np.hypot(s_jk, s_sc)),
        c_ampX=float(np.exp(np.mean(np.log([a["c_ampX"] for a in rows])))),
        KALIB=kal, W_X=float(np.mean([a["W_X"] for a in rows])),
        W_amp=float(np.mean([a["W_amp"] for a in rows])),
        theta=kal / (d["artik"]["gE"] * d["artik"]["gX"] ** 2),
        tau_eff=[float(a["tau_eff"]) for a in rows],
        c_bant=[float(x) for x in c],
        c_WX9=(float(np.exp(np.mean(np.log(
            [b["KALIB_u2"] / b["W_X"] for b in B9])))) if B9 else float("nan")))
    return r


def main():
    os.makedirs(SCR170, exist_ok=True)
    D = K1.yukle()
    tau_l = K1.ladder(D["Hkeskin"]["L"], 0.95)
    R = {ad: gaz_satiri(ad, D[ad], tau_l) for ad in sorted(D)}
    th0 = R["Hkeskin"]["theta"]

    print("=" * 132)
    print("K0 — MUHASEBE: HER `c` HANGİ KOŞULDA ÖLÇÜLDÜ  "
          "[ortak bant penceresi lo∈[0.52,0.68], SNR≥3, R_bant≥0.98, "
          "τ_eff<0.85]")
    print("=" * 132)
    print(f"{'gaz':9s} {'eksen':17s} {'λ':>5s} {'kesim':>13s} {'φ':>6s} "
          f"{'N':>7s} {'T/T₀':>5s} {'nb':>3s} | {'σ_ds':>6s} {'σ_X̃':>6s} "
          f"{'σ_Ĉ':>6s} | {'g_cal':>6s} {'θ':>6s} {'KALİB':>6s} {'W_amp':>6s} "
          f"{'W_X':>6s} | {'c_WX':>7s} {'±tot':>6s} {'c_ampX':>6s}")
    print("-" * 132)
    T0 = R["Hkeskin"]["T"]
    for ad in ("L115", "Hkeskin", "son", "L085", "L070", "L060",
               "K090", "K070", "HA4", "E060", "HkT2a", "HkT2b",
               "HkT4a", "HkT4b"):
        if ad not in R or "c_WX" not in R[ad]:
            continue
        r = R[ad]
        kes = ("erfc%.2f/%.3f" % tuple(r["pen"]) if r["pen"]
               else "τ≤%.2f" % (r["tau_ust"] or 1.0))
        if r["gercek"]:
            kes = "GERÇEK ζ"
        print(f"{ad:9s} {r['eksen']:17s} "
              f"{('—' if r['lam'] is None else '%.2f' % r['lam']):>5s} "
              f"{kes:>13s} {r['phi']:6.4f} {r['N']:7d} {r['T']/T0:5.2f} "
              f"{r['nb5']:3d} | {r['sigds']:6.4f} {r['sigX']:6.4f} "
              f"{r['sigC']:6.4f} | {r['gcal']:6.4f} {r['theta']:6.4f} "
              f"{r['KALIB']:6.4f} {r['W_amp']:6.4f} {r['W_X']:6.4f} | "
              f"{r['c_WX']:7.4f} {r['s_tot']:6.4f} {r['c_ampX']:6.4f}")

    # ---- YENİDEN ÜRETİM SINAVI ---------------------------------------
    print("\n" + "=" * 132)
    print("K0(b) — 168/169'UN YAZILI SAYILARIYLA YÜZLEŞME (sapma > %0.15 → '!')")
    print("=" * 132)
    print(f"{'gaz':9s} | {'c_WX5 bu':>9s} {'169':>7s} | {'c_WX9 bu':>9s} "
          f"{'168':>7s} | {'g_cal bu':>9s} {'168':>7s} | {'θ bu':>7s} "
          f"{'168':>7s} | {'c_ampX bu':>10s} {'168':>7s}")
    sapma = []
    for ad, (c5, c9, gc, th, ca) in REF.items():
        if ad not in R or "c_WX" not in R[ad]:
            continue
        r = R[ad]
        line = f"{ad:9s} | {r['c_WX']:9.4f} {c5:7.4f} | "
        for bu, ref in ((r["c_WX9"], c9), (r["gcal"], gc), (r["theta"], th),
                        (r["c_ampX"], ca)):
            if ref is None:
                line += f"{bu:9.4f} {'—':>7s} | "
            else:
                bad = abs(bu / ref - 1) > 0.0015
                if bad:
                    sapma.append((ad, ref, bu))
                line += f"{bu:9.4f} {ref:7.4f}{'!' if bad else ' '}| "
        bad5 = abs(r["c_WX"] / c5 - 1) > 0.0015
        if bad5:
            sapma.append((ad + "/c5", c5, r["c_WX"]))
        print(line + ("  <-- c5 SAPMA" if bad5 else ""))
    print(f"\n  θ/θ₀ (Hkeskin = {th0:.4f}):  " + "  ".join(
        f"{ad}={R[ad]['theta']/th0:.4f}" for ad in
        ("son", "L085", "L070", "K090", "K070", "HA4", "E060")
        if ad in R and "theta" in R[ad]))
    print(f"  168 §A2.7:                  son=1.0290  L085=1.0084  L070=—     "
          f"K090=0.9080  K070=0.8038  HA4=0.7975  E060=0.7945")
    print(f"\n  YENİDEN ÜRETİM: {'TAM (sapma yok)' if not sapma else str(sapma)}")

    # ---- λ EKSENİNİN ÇARPAN AYRIŞTIRMASI ------------------------------
    print("\n" + "=" * 132)
    print("K0(c) — λ EKSENİ: c_WX = (g_E · g_X² · θ)/W_X  ÇARPAN ÇARPAN "
          "(Hkeskin'e oranla)")
    print("=" * 132)
    h = R["Hkeskin"]
    print(f"{'gaz':7s} {'λ':>5s} | {'σ_ds/σ₀':>8s} {'σ_X̃/σ₀':>8s} "
          f"{'σ_Ĉ/σ₀':>8s} | {'g_E/g₀':>7s} {'g_X/g₀':>7s} {'θ/θ₀':>7s} "
          f"{'W_X/W₀':>7s} {'W_amp/W₀':>8s} | {'c/c₀ ÇARPIM':>11s} "
          f"{'c/c₀ ÖLÇÜLEN':>12s}")
    lam_rows = []
    for ad in ("L115", "Hkeskin", "L085", "L070", "L060"):
        if ad not in R or "c_WX" not in R[ad]:
            continue
        r = R[ad]
        carp = ((r["gE"] / h["gE"]) * (r["gX"] / h["gX"]) ** 2
                * (r["theta"] / h["theta"]) / (r["W_X"] / h["W_X"]))
        print(f"{ad:7s} {r['lam']:5.2f} | {r['sigds']/h['sigds']:8.4f} "
              f"{r['sigX']/h['sigX']:8.4f} {r['sigC']/h['sigC']:8.4f} | "
              f"{r['gE']/h['gE']:7.4f} {r['gX']/h['gX']:7.4f} "
              f"{r['theta']/h['theta']:7.4f} {r['W_X']/h['W_X']:7.4f} "
              f"{r['W_amp']/h['W_amp']:8.4f} | {carp:11.4f} "
              f"{r['c_WX']/h['c_WX']:12.4f}")
        lam_rows.append(dict(gaz=ad, lam=r["lam"], carp=carp,
                             olc=r["c_WX"] / h["c_WX"]))

    json.dump(dict(C_HIP=C_HIP, gaz=R, lam=lam_rows, sapma=sapma),
              open(SCR170 + "/K0.json", "w"), indent=1, default=float)
    print(f"\n-> {SCR170}/K0.json")


if __name__ == "__main__":
    main()
