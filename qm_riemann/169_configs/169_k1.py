"""
169 — K1 (HASSASİYET): c'nin jackknife'lı ölçümü vs 4/π² = 0.4052847
====================================================================
Girdi: scratchpad/167/C_<gaz>.json (YENİ GAZ/ÖLÇÜM YOK; 168 ile aynı dosyalar)
Ölçüm parçası kopyalanmaz: bant birleştirmesi `166_T1.bant_agg`, jackknife
`166_T1._jk` AYNEN import edilir (167_olcum.py'nin kullandığı fonksiyonlar).

ÜYE TANIMI (168 §A2.8'in hükmü):   c_WX ≡ KALİB_u2 / W_X
  (`c_u2` alanı JSON'da KALİB_u2/(W_amp·W_X)'tir — W_amp çift sayımdır.)
Jackknife: bant içi 8 grup (160'ın round-robin'i), silme-1; W_amp/W_X de
her silme için gp-ağırlığıyla YENİDEN hesaplanır (oran tahmincisinde
doğrusal hata yayılımı yok) — 167_olcum.py'nin `_agg`'iyle aynı cebir.

Çıktı: scratchpad/169/K1.json  + ekrana tablo
"""
import glob
import importlib
import json
import os
import sys
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
for _p in ("166_configs", "165_configs", "163_configs", "160_configs",
           "159_configs", "155_configs", "154_configs"):
    sys.path.insert(0, str(QM / _p))
T166 = importlib.import_module("166_T1")

SCR = ("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
       "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
SCR167, SCR169 = SCR + "/167", SCR + "/169"
LO_MIN, LO_MAX, SNR_MIN, R_MIN = 0.52, 0.68, 3.0, 0.98
C_HIP = 4.0 / np.pi ** 2                       # = 0.40528473...
NJACK = 8


def yukle():
    D = {}
    for p in sorted(glob.glob(SCR167 + "/C_*.json")):
        D[os.path.basename(p)[2:-5]] = json.load(open(p))
    return D


def agg_uyeler(L, msk=None):
    """166_T1.bant_agg + gp-ağırlıklı W çarpanları + üye tanımları."""
    out = T166.bant_agg(L, msk)
    if msk is None:
        msk = np.ones(len(L), bool)
    gp = np.array([r["gp"] for r in L])[msk]
    gw = gp / gp.sum()
    for k in ("W_amp", "W_X", "W_pos"):
        out[k] = float((np.array([r[k] for r in L])[msk] * gw).sum())
    out["c_WX"] = out["KALIB_u2"] / out["W_X"]
    out["c_ampX"] = out["KALIB_u2"] / (out["W_amp"] * out["W_X"])
    return out


def band_jk(L, njack=NJACK):
    """Bant-içi 8-grup silme-1 jackknife; her üye için ayrı."""
    tam = agg_uyeler(L)
    grp = np.array([r["grup"] for r in L])
    keys = ("KALIB_u2", "c_WX", "c_ampX", "W_X", "W_amp", "tau_eff")
    jk = {k: [] for k in keys}
    for g in range(njack):
        m = grp != g
        if not m.any():
            continue
        try:
            a = agg_uyeler(L, m)
        except (ZeroDivisionError, FloatingPointError):
            continue
        for k in keys:
            jk[k].append(a[k])
    for k in keys:
        tam[f"s{k}"] = T166._jk(jk[k])
    return tam


def saglikli(d, lo_max=LO_MAX):
    out = []
    for b in d["bant"]:
        if not b.get("olculdu"):
            continue
        if b["lo"] < LO_MIN - 1e-9 or b["lo"] > lo_max + 1e-9:
            continue
        if b["tau_eff"] > 0.85 or b["Ms2"] <= 0:
            continue
        if b["R_bant"] < R_MIN or b["SNR"] < SNR_MIN:
            continue
        out.append(b)
    return out


def phi_of(d, tau_ladder):
    """BOŞ ÇİZGİ KESRİ (168 §A2.7 ile birebir)."""
    from math import erfc
    pen, tu = d["pen"], d["tau_ust"] or 1.0
    if pen is not None:
        w = np.array([0.5 * erfc((t - pen[0]) / pen[1]) for t in tau_ladder])
    else:
        w = (tau_ladder <= tu + 1e-12).astype(float)
    return float(1 - w.mean())


def ladder(L=12.0296, tau_c=0.95):
    """168_olcek.ladder ile aynı merdiven (τ listesi yeter)."""
    qmax = int(np.exp(tau_c * L)) + 1
    sieve = np.ones(qmax + 1, bool)
    sieve[:2] = False
    for i in range(2, int(qmax ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = False
    pr = np.nonzero(sieve)[0]
    q = []
    for p in pr:
        v = int(p)
        while v <= qmax:
            q.append(v)
            v *= int(p)
    tau = np.log(np.array(q, float)) / L
    return tau[tau <= tau_c + 1e-12]


def main():
    D = yukle()
    os.makedirs(SCR169, exist_ok=True)
    tau_l = ladder(D["Hkeskin"]["L"], 0.95)
    print("=" * 104)
    print(f"K1 — c (W_X üyesi) vs 4/π² = {C_HIP:.7f}   "
          f"[ortak pencere lo ∈ [{LO_MIN},{LO_MAX}], SNR≥{SNR_MIN}, R≥{R_MIN}]")
    print("=" * 104)
    R = {}
    for ad in sorted(D):
        d = D[ad]
        B = saglikli(d)
        if not B:
            continue
        rows = []
        print(f"\n--- {ad}  (λ={d['lam']} τ_ust={d['tau_ust']} pen={d['pen']} "
              f"N={d['N']} T={d['T']:.0f}) ---")
        print("   τ_eff   KALİB_u2±jk       W_X      W_amp    "
              "**c_WX ± jk**        (c−4/π²)/σ   c_ampX")
        for b in B:
            a = band_jk(b["cizgi"])
            z = (a["c_WX"] - C_HIP) / a["sc_WX"] if a["sc_WX"] > 0 else np.nan
            rows.append(dict(tau_eff=a["tau_eff"], c=a["c_WX"], s=a["sc_WX"],
                             K=a["KALIB_u2"], sK=a["sKALIB_u2"],
                             WX=a["W_X"], Wa=a["W_amp"], z=z,
                             c_ampX=a["c_ampX"], lo=b["lo"]))
            print(f"  {a['tau_eff']:.4f}  {a['KALIB_u2']:.4f}±"
                  f"{a['sKALIB_u2']:.4f}   {a['W_X']:.4f}   {a['W_amp']:.4f}   "
                  f"{a['c_WX']:.4f} ± {a['sc_WX']:.4f}    {z:+7.2f}      "
                  f"{a['c_ampX']:.4f}")
        c = np.array([r["c"] for r in rows])
        s = np.array([r["s"] for r in rows])
        lg = np.log(c)
        cg = float(np.exp(lg.mean()))
        # gaz-düzeyi hata: (i) jackknife'ların ortalamaya taşınması
        #                  (ii) bantlar arası saçılım (sd/√n)
        s_jk = float(np.sqrt(np.sum((s / c) ** 2)) / len(c)) * cg
        s_sc = (float(np.std(lg, ddof=1) / np.sqrt(len(lg))) * cg
                if len(lg) > 1 else float("nan"))
        s_tot = float(np.hypot(s_jk, s_sc))
        te = np.array([r["tau_eff"] for r in rows])
        eg = float(np.polyfit(te, lg, 1)[0]) if len(te) > 2 else float("nan")
        z_jk = (cg - C_HIP) / s_jk
        z_tot = (cg - C_HIP) / s_tot
        R[ad] = dict(c=cg, s_jk=s_jk, s_sc=s_sc, s_tot=s_tot, nb=len(c),
                     z_jk=z_jk, z_tot=z_tot, egim=eg,
                     phi=phi_of(d, tau_l),
                     gcal=d["artik"]["gE"] * d["artik"]["gX"] ** 2,
                     WX=float(np.mean([r["WX"] for r in rows])),
                     T=d["T"], N=d["N"], lam=d["lam"] or 1.0,
                     nline=d["artik"]["nline"])
        print(f"  ⇒ c({ad}) = {cg:.4f}  ± {s_jk:.4f}(jk) ± {s_sc:.4f}(bant) "
              f"= ±{s_tot:.4f}  [{len(c)} bant, τ-eğim {eg:+.3f}]")
        print(f"     4/π² uzaklığı:  {(cg-C_HIP)/C_HIP*100:+.2f}%  = "
              f"{z_jk:+.2f}σ_jk = {z_tot:+.2f}σ_tot")

    # ---------------- 168 TABLOSUNUN YENİDEN ÜRETİMİ (dürüstlük) --------
    print("\n" + "=" * 104)
    print("DÜRÜSTLÜK — 168 §A2.6/§A2.8 tablosunun yeniden üretimi "
          "(aynı pencere, bağımsız kod yolu)")
    print("=" * 104)
    print("  gaz        c_ampX(bu)  168'de   |  c_WX(bu)  168'de  |  "
          "g_cal(bu)  168'de  |  θ(bu)   168'de")
    ref168 = {  # 168 raporunun A2.6 + A2.8 sayıları
        "Hkeskin": (0.5774, 0.4051, 0.2895, 0.8884),
        "son": (0.5664, 0.4053, 0.2971, 0.9142),
        "L085": (0.5171, 0.3751, 0.2912, 0.8959),
        "K090": (0.5393, 0.3709, 0.2877, 0.8067),
        "K070": (0.5154, 0.3480, 0.3128, 0.7142),
        "HA4": (0.5061, 0.3345, 0.3194, 0.7086),
        "E060": (0.5216, 0.3701, 0.3433, 0.7059)}
    for ad, (ca, cw, gc, th) in ref168.items():
        if ad not in R:
            continue
        d, B = D[ad], saglikli(D[ad])
        camp = float(np.exp(np.mean(np.log([b["c_u2"] for b in B]))))
        kal = float(np.exp(np.mean(np.log([b["KALIB_u2"] for b in B]))))
        g = R[ad]["gcal"]
        print(f"  {ad:10s} {camp:.4f}   {ca:.4f}   |  {R[ad]['c']:.4f}  "
              f"{cw:.4f}  |  {g:.4f}   {gc:.4f}  |  {kal/g:.4f}  {th:.4f}")

    # ---------------- KESİM YASASININ φ→0 KESİŞİMİ ----------------------
    print("\n" + "=" * 104)
    print("K1(b) — KESİM YASASININ φ→0 KESİŞİMİ (mutlak normalizasyon)")
    print("=" * 104)
    kes = ["Hkeskin", "K090", "K070", "HA4", "E060"]
    kes = [a for a in kes if a in R]
    phi = np.array([R[a]["phi"] for a in kes])
    cw = np.array([R[a]["c"] for a in kes])
    sw = np.array([R[a]["s_tot"] for a in kes])
    g0 = R["Hkeskin"]["gcal"]
    cw_g = cw * (g0 / np.array([R[a]["gcal"] for a in kes]))  # g_cal düzeltmeli
    print("  gaz        φ       c_WX ± σ        c_WX·(g₀/g)   (=θ-normalize)")
    for i, a in enumerate(kes):
        print(f"  {a:10s} {phi[i]:.4f}  {cw[i]:.4f}±{sw[i]:.4f}    "
              f"{cw_g[i]:.4f}")
    out_fit = {}
    for et, y in (("ham c_WX", cw), ("g_cal-düzeltmeli", cw_g)):
        w = 1.0 / sw ** 2
        A = np.vstack([np.ones_like(phi), phi]).T
        Aw = A * w[:, None]
        cov = np.linalg.inv(A.T @ Aw)
        p = cov @ (Aw.T @ y)
        se = np.sqrt(np.diag(cov))
        c0, sl = p[0], p[1]
        beta = -sl / c0
        z = (c0 - C_HIP) / se[0]
        print(f"\n  [{et}]  c(φ) = c₀ (1 − β φ):")
        print(f"     c₀ = {c0:.4f} ± {se[0]:.4f}   β = {beta:.4f}   "
              f"(4/π² = {C_HIP:.4f} → {(c0-C_HIP)/C_HIP*100:+.2f}% = {z:+.2f}σ)")
        out_fit[et] = dict(c0=float(c0), sc0=float(se[0]), beta=float(beta),
                           z=float(z))

    # ---------------- BANT PENCERESİ DUYARLILIĞI ------------------------
    print("\n" + "=" * 104)
    print("K1(c) — BANT PENCERESİ DUYARLILIĞI (168 §A2.8 DOKUZ bant kullandı; "
          "168 §A2.6/A2.7 ise lo≤0.68)")
    print("=" * 104)
    print("  gaz         c_WX(lo≤0.68)  c_WX(lo≤0.80, 'dokuz bant')  "
          "168 §A2.8'in sayısı")
    for ad in ("Hkeskin", "son", "L085", "K090", "K070", "HA4", "E060"):
        if ad not in R:
            continue
        B9 = saglikli(D[ad], lo_max=0.80)
        c9 = float(np.exp(np.mean(np.log(
            [b["KALIB_u2"] / b["W_X"] for b in B9]))))
        print(f"  {ad:10s}  {R[ad]['c']:.4f}         {c9:.4f}  ({len(B9)} bant)"
              f"                {ref168[ad][1]:.4f}")
        R[ad]["c9"] = c9
        R[ad]["nb9"] = len(B9)

    # ---------------- TAM DÜZELTİLMİŞ ÇÖKME TABLOSU ---------------------
    beta = out_fit["g_cal-düzeltmeli"]["beta"]
    print("\n" + "=" * 104)
    print(f"K1(d) — TAM DÜZELTME: c* = c_WX·(g₀/g_cal)/(1−βφ),  β = {beta:.4f}"
          "   [H-C1: c* = 4/π² her gazda]")
    print("=" * 104)
    print("  gaz         φ      g₀/g_cal   c_WX     **c***    c*/(4/π²)−1   "
          "eksen")
    for ad in sorted(R):
        r = R[ad]
        f = (g0 / r["gcal"]) / (1 - beta * r["phi"])
        cs = r["c"] * f
        eks = ("kesim" if r["phi"] > 0.1 else
               ("pencere" if r["T"] < 1e5 else
                ("λ" if r["lam"] != 1.0 else "taban")))
        print(f"  {ad:10s} {r['phi']:.4f}  {g0/r['gcal']:.4f}   "
              f"{r['c']:.4f}   {cs:.4f}   {100*(cs/C_HIP-1):+7.2f}%    {eks}")
        R[ad]["c_star"] = cs

    json.dump(dict(C_HIP=C_HIP, gaz=R, fit=out_fit, beta=beta),
              open(SCR169 + "/K1.json", "w"), indent=1, default=float)
    print(f"\n-> {SCR169}/K1.json")


if __name__ == "__main__":
    main()
