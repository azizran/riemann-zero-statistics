"""
159 — ANALİZ: T0 (kinematik omurga) · T1 (mekanizma denklemi) ·
T2 (S'nin kanal ayrışımı) · T3 (a_pred, b işareti)

Fit konvansiyonu 158/157'nin `bolum_E`'siyle birebir:
  apsis τ_eff · ağırlık 1/σ_jk · polyfit(x,y,2,w=1/σ) ·
  τ₀ = köklerden x ortalamasına en yakını · a = c₁+2c₀τ₀ · b = c₀
"""
import importlib
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
C = importlib.import_module("159_cekirdek")
TWO_PI = 2 * np.pi
FOUR_PI = 4 * np.pi
SCR = C.SCR

GAZ = ["son", "orta", "keskin", "A4", "P1"]
ETI = {"son": "gerçek-son", "orta": "gerçek-orta", "keskin": "keskin",
       "A4": "A4", "P1": "P1(GUE)"}


def yukle(veri, taban, g):
    p = SCR / f"S_{veri}_t{taban}_{g}.json"
    if not p.exists():
        return None
    return json.load(open(p))


def bantlar(d, lo=None, hi=None):
    B = [b for b in d["bantlar"] if b.get("olculdu")]
    if lo is not None:
        B = [b for b in B if b["tau"] >= lo - 1e-9]
    if hi is not None:
        B = [b for b in B if b["tau"] <= hi + 1e-9]
    return B


def saglik(B):
    """Bant sağlığı: |Γ|≤1 (Γ normalize korelasyon), |φ|≤2.8 (π'ye sarma
    payı), |τ_eff−τ̄| ≤ yarım bant. İLK bozulmadan sonraki bantlar da
    düşer (tahminci bir kez bozulunca üst bantlar bozulmayı devralır)."""
    ok = []
    kir = True
    for b in B:
        yw = 0.5 * (b["hi"] - b["lo"])
        iyi = (abs(b["absG"]) <= 1.0 and abs(b["phi"]) <= 2.8
               and abs(b["tau_eff"] - b["tau"]) <= yw + 1e-9)
        if kir and iyi:
            ok.append(True)
        else:
            kir = False
            ok.append(False)
    return ok


def fit2(x, y, s=None):
    """158'in kuadratik fiti: (τ₀, a=dy/dx|τ₀, b=c₀, χ²/dof)."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    if s is None or not np.all(np.isfinite(s)) or np.any(np.asarray(s) <= 0):
        w = np.ones_like(x)
    else:
        w = 1.0 / np.asarray(s, float)
    c = np.polyfit(x, y, 2, w=w)
    r = np.roots(c)
    r = r[np.abs(r.imag) < 1e-9].real
    t0 = float(r[np.argmin(np.abs(r - x.mean()))]) if len(r) else float("nan")
    a = float(c[1] + 2 * c[0] * t0)
    rez = y - np.polyval(c, x)
    dof = max(len(x) - 3, 1)
    chi = float(np.sum((rez * w) ** 2) / dof) if s is not None else float("nan")
    return t0, a, float(c[0]), chi


def T0():
    print("=" * 78)
    print("T0 — KİNEMATİK OMURGA: φ_Γ = (4π·τ_eff − 2π) + δ")
    print("=" * 78)
    d = yukle("son", 0.4, "g158")
    B = bantlar(d)
    print("  τ̄     τ_eff    φ_meas   4πτe−2π    δ      δ/φ    |Γ|   |Γ_arın|")
    for b in B:
        print(f"  {b['tau']:.2f}  {b['tau_eff']:.4f} {b['phi']:+8.4f} "
              f"{b['kinematik']:+9.4f} {b['delta']:+8.4f} "
              f"{b['delta']/b['phi']:+7.3f} {b['absG']:.3f}  {b['absGc']:.3f}")
    y = np.array([b["phi"] for b in B]); x = np.array([b["tau_eff"] for b in B])
    k = np.array([b["kinematik"] for b in B])
    dl = np.array([b["delta"] for b in B])
    print(f"  → korel(φ, omurga) = {np.corrcoef(y, k)[0,1]:.6f}   "
          f"sd(φ)={y.std():.4f}  sd(omurga)={k.std():.4f}  sd(δ)={dl.std():.4f}"
          f"   |δ|/|φ| medyan = {np.median(np.abs(dl)/np.abs(y)):.3f}")
    print("\n  τ₀ AYRIŞIMI:  omurga tek başına τ₀ = 1/2 verir; ölçülen kayma "
          "−δ(½)/a")
    print("  gaz          τ₀(ölçülen)   δ(τ=½)   ½ − δ(½)/a   fark")
    for v in ("son", "orta", "keskin", "A4", "P1"):
        d2 = yukle(v, 0.4, "g158")
        if d2 is None:
            continue
        Bv = bantlar(d2, 0.43, 0.61)
        xv = [b["tau_eff"] for b in Bv]
        t0, a, _, _ = fit2(xv, [b["phi"] for b in Bv],
                           [b["sPhi_jk"] for b in Bv])
        cd = np.polyfit(xv, [b["delta"] for b in Bv], 2)
        dhalf = float(np.polyval(cd, 0.5))
        tah = 0.5 - dhalf / a
        print(f"  {ETI.get(v,v):11s}  {t0:.4f}      {dhalf:+.4f}   {tah:.4f}"
              f"     {t0-tah:+.4f}")


def T1(veri="son", taban=0.4, g="t1"):
    d = yukle(veri, taban, g)
    if d is None:
        print(f"  [{veri} t{taban} {g}] YOK")
        return
    B = bantlar(d)
    print()
    print("=" * 78)
    print(f"T1 — {ETI.get(veri,veri)} / taban {taban} / ızgara {g} "
          f"({len(B)} bant)")
    print("=" * 78)
    sg = saglik(B)
    print("  τ_eff   |Γ|   φ_meas     A·S      φ/(A·S)  ||  δ_meas    A·S    "
          "δ/(A·S)   M1     δ/M1   ReΓc  δ/(A·S/ReΓc)  sağ")
    sat = []
    for b, s in zip(B, sg):
        AS = b["AS"]["S_tam"]
        ReGc = b["absGc"] * np.cos(b["delta"])
        r1 = b["phi"] / AS if AS else np.nan
        r2 = b["delta"] / AS if AS else np.nan
        r3 = b["delta"] / b["M1"] if b["M1"] else np.nan
        r4 = b["delta"] / (AS / ReGc) if AS else np.nan
        sat.append((r1, r2, r3, r4, s))
        print(f"  {b['tau_eff']:.4f} {b['absG']:.3f} {b['phi']:+8.4f} "
              f"{AS:+8.4f} {r1:+8.2f}  || {b['delta']:+8.4f} {AS:+8.4f} "
              f"{r2:+7.3f} {b['M1']:+8.4f} {r3:+6.3f}  {ReGc:+.3f} "
              f"{r4:+8.3f}   {'✓' if s else '✗'}")
    sat = np.array(sat, float)
    m = sat[:, 4] > 0
    ad = ["φ/(A·S)  [KALEM birebir]", "δ/(A·S)  [kinematikten arınmış]",
          "δ/M1     [ρ, bütün mertebeler]", "δ/(A·S/ReΓ) [sönüm düzeltmeli]"]
    print(f"\n  ORAN ÖZETİ — SAĞLIKLI {int(m.sum())}/{len(B)} bant "
          f"(ideal 1.00, MÜHÜR eşiği ±%25 ⇒ 0.75–1.25)")
    for i, a in enumerate(ad):
        v = sat[m, i]
        v = v[np.isfinite(v)]
        if not len(v):
            continue
        ic = int(np.sum((v > 0.75) & (v < 1.25)))
        print(f"    {a:32s} ort={v.mean():+8.3f} medyan={np.median(v):+8.3f} "
              f"menzil=[{v.min():+.2f},{v.max():+.2f}]  ±%25 içinde "
              f"{ic}/{len(v)}")
    return B


def T1_tanim(veri="son", taban=0.4, g="t1"):
    """P_loc TANIM sınavı + kayan-pencere duyarlılığı + off-line sistematiği."""
    d = yukle(veri, taban, g)
    if d is None:
        return
    B = bantlar(d)
    print()
    print("-" * 78)
    print(f"T1b — P_loc TANIMI, PENCERE ve OFF-LINE SİSTEMATİĞİ "
          f"({ETI.get(veri,veri)}, taban {taban}, {g})")
    print("-" * 78)
    print("  τ_eff    δ     A·S(ρ)   A·S(|ĉ|²,W=32/64/128)        "
          "A·S(⟨η²⟩_W64)  pedestal(W=64)")
    for b in B:
        print(f"  {b['tau_eff']:.4f} {b['delta']:+7.4f} {b['AS']['S_tam']:+8.4f}  "
              f"{b['AS']['S_dem32']:+8.5f} {b['AS']['S_dem64']:+8.5f} "
              f"{b['AS']['S_dem128']:+8.5f}   {b['AS']['S_tot64']:+8.5f}   "
              f"{b['ped']['ped64']:9.1f}")
    print("\n  KAYAN PENCERE (ağırlığı pencerelemek ≡ X'i pencerelemek):")
    print("  τ_eff   A·S(X bond)  A·S(X̄_32)  A·S(X̄_64)  A·S(X̄_128)   "
          "X̄_64/bond")
    for b in B:
        a0 = b["AS"]["S_tam"]
        print(f"  {b['tau_eff']:.4f} {a0:+11.4f} {b['AS']['S_xw32']:+11.5f} "
              f"{b['AS']['S_xw64']:+11.5f} {b['AS']['S_xw128']:+11.5f}   "
              f"{b['AS']['S_xw64']/a0:+8.4f}")
    print("\n  ARA-NOKTA (off-line) REFERANSI: S_on − S_off  vs  yalnız S_on")
    print("  τ_eff   A·S(on−off)   A·S(yalnız on)   fark      %")
    for b in B:
        a0, a1 = b["AS"]["S_tam"], b["AS_on"]["S_tam"]
        print(f"  {b['tau_eff']:.4f} {a0:+12.5f} {a1:+15.5f} {a1-a0:+10.5f} "
              f"{100*(a1-a0)/abs(a0):+8.1f}")
    print("\n  X DEĞİŞKENİ: dsΔ (150/155'in) vs X̃ (tam etkin bond adımı)")
    print("  τ_eff   A·S(dsΔ)     A·S(X̃)      fark %")
    for b in B:
        a0, a1 = b["AS"]["S_tam"], b["AS"]["S_xtil"]
        print(f"  {b['tau_eff']:.4f} {a0:+12.5f} {a1:+12.5f} "
              f"{100*(a1-a0)/abs(a0):+8.2f}")


def T2(veri="son", taban=0.4, g="t1"):
    d = yukle(veri, taban, g)
    if d is None:
        return
    B = bantlar(d)
    print()
    print("=" * 78)
    print(f"T2 — S'nin KANAL AYRIŞIMI ({ETI.get(veri,veri)}, taban {taban}, {g})")
    print(f"     bond payı Cov(X_k,X_tam)/σΔ²: "
          f"lad={d['kov']['lad']/d['sA2']:+.4f} "
          f"eta={d['kov']['eta']/d['sA2']:+.4f} "
          f"dri={d['kov']['dri']/d['sA2']:+.6f}")
    print("=" * 78)
    print("  τ_eff   A·S_tam    A·S_lad    A·S_eta    A·S_dri   kapanış   "
          "lad%    eta%")
    for b in B:
        t, l, e, dr = (b["AS"]["S_tam"], b["AS"]["S_lad"],
                       b["AS"]["S_eta"], b["AS"]["S_dri"])
        print(f"  {b['tau_eff']:.4f} {t:+10.5f} {l:+10.5f} {e:+10.5f} "
              f"{dr:+10.2e} {abs(t-(l+e+dr)):.1e} {100*l/t:+7.1f} "
              f"{100*e/t:+7.1f}")


def T3(gazlar=("son", "orta", "keskin", "A4", "P1"), taban=0.4,
       pen=(0.43, 0.61)):
    print()
    print("=" * 78)
    print(f"T3 — a_pred ve b işareti (ızgara g158, taban {taban}, "
          f"fit penceresi τ̄∈[{pen[0]},{pen[1]}], 158'in W-A'sı)")
    print("=" * 78)
    print(f"  kinematik omurgan eğimi: 4π = {FOUR_PI:.4f} (τ'da tam DOĞRUSAL "
          f"⇒ b'ye katkısı SIFIR)")
    print()
    print("  gaz         n   τ₀(φ)    a(φ)    b(φ)  |  τ₀(δ)   dδ/dτ   b(δ)  "
          "| 4π+dδ/dτ | dA·S/dτ  b(A·S) | 4π+dA·S/dτ | 4π+dM1/dτ")
    out = {}
    for v in gazlar:
        d = yukle(v, taban, "g158")
        if d is None:
            continue
        B = bantlar(d, pen[0], pen[1])
        if len(B) < 4:
            continue
        x = [b["tau_eff"] for b in B]
        sp = [b["sPhi_jk"] for b in B]
        sd = [b["sDelta_jk"] for b in B]
        sa = [b["sAS_jk"] for b in B]
        t0p, ap, bp, _ = fit2(x, [b["phi"] for b in B], sp)
        t0d, ad_, bd, _ = fit2(x, [b["delta"] for b in B], sd)
        t0s, as_, bs, _ = fit2(x, [b["AS"]["S_tam"] for b in B], sa)
        t0m, am, bm, _ = fit2(x, [b["M1"] for b in B], sd)
        # τ₀(φ)'da eğimler (kinematikten arınmış eğriler için)
        def egim(y, s):
            c = np.polyfit(x, y, 2, w=1.0 / np.asarray(s))
            return float(c[1] + 2 * c[0] * t0p)
        dd = egim([b["delta"] for b in B], sd)
        ds = egim([b["AS"]["S_tam"] for b in B], sa)
        dm = egim([b["M1"] for b in B], sd)
        out[v] = dict(t0p=t0p, ap=ap, bp=bp, t0d=t0d, bd=bd, dd=dd,
                      t0s=t0s, as_=as_, bs=bs, ds=ds, dm=dm, bm=bm, n=len(B),
                      a_AS=FOUR_PI + ds, a_M1=FOUR_PI + dm, a_dl=FOUR_PI + dd)
        print(f"  {ETI.get(v,v):11s} {len(B):2d} {t0p:.4f} {ap:7.3f} {bp:+7.2f} "
              f"| {t0d:.4f} {dd:+7.3f} {bd:+7.2f} | {FOUR_PI+dd:9.3f} "
              f"| {ds:+8.3f} {bs:+7.2f} | {FOUR_PI+ds:10.3f} "
              f"| {FOUR_PI+dm:9.3f}")
    print()
    print("  KALEM'in birebir T3'ü: a_pred = d(A·S)/dτ | (A·S)'nin SIFIR "
          "GEÇİŞİNDE")
    print("  gaz          τ₀(A·S)  d(A·S)/dτ|τ₀(A·S)   a_meas   a_pred/a_meas")
    for v, r in out.items():
        print(f"  {ETI.get(v,v):11s}  {r['t0s']:.4f} {r['as_']:+16.3f} "
              f"{r['ap']:9.3f} {r['as_']/r['ap']:14.3f}")
    return out


def T3_ref(out):
    print()
    print("  REFERANSLAR: gerçek a = 10.713 ± 0.059 (158, 32 fit), "
          "b = −7.22 ± 0.79 (16/16 negatif); σ_konv(gerçek) = 0.075")
    print()
    print("  a ÖNGÖRÜSÜ — hata bütçesi (σ_konv birimlerinde, gerçek için)")
    print("  gaz         a_meas   4π+dδ/dτ   4π+dM1/dτ  hata%   4π+dA·S/dτ "
          " hata%   4π (çıplak)  hata%")
    for v, r in out.items():
        e1 = 100 * (r["a_M1"] - r["ap"]) / r["ap"]
        e2 = 100 * (r["a_AS"] - r["ap"]) / r["ap"]
        e3 = 100 * (FOUR_PI - r["ap"]) / r["ap"]
        print(f"  {ETI.get(v,v):11s} {r['ap']:7.3f} {r['a_dl']:10.3f} "
              f"{r['a_M1']:11.3f} {e1:+6.1f} {r['a_AS']:11.3f} {e2:+7.1f} "
              f"{FOUR_PI:11.3f} {e3:+7.1f}")
    ks = list(out)
    am = np.array([out[v]["ap"] for v in ks])
    for nm, key in (("4π+dδ/dτ", "a_dl"), ("4π+dM1/dτ", "a_M1"),
                    ("4π+dA·S/dτ", "a_AS")):
        pr = np.array([out[v][key] for v in ks])
        sr = (np.argsort(np.argsort(am)) == np.argsort(np.argsort(pr))).all()
        print(f"    korel(a_meas, {nm}) = {np.corrcoef(am, pr)[0,1]:+.4f}  "
              f"sıralama {'AYNI' if sr else 'FARKLI'}  "
              f"({' < '.join(ETI.get(k,k) for k in np.array(ks)[np.argsort(pr)])})")
    print()
    print("  b İŞARETİ (kinematik omurga doğrusal ⇒ b(φ) ≡ b(δ) olmalı):")
    print("  gaz         b(φ)     b(δ)   |Δ|     b(M1)    b(A·S)  "
          "işaret(M1) işaret(A·S)")
    for v, r in out.items():
        u1 = "✓" if np.sign(r["bm"]) == np.sign(r["bp"]) else "✗"
        u2 = "✓" if np.sign(r["bs"]) == np.sign(r["bp"]) else "✗"
        print(f"  {ETI.get(v,v):11s} {r['bp']:+7.2f} {r['bd']:+7.2f} "
              f"{abs(r['bp']-r['bd']):.3f} {r['bm']:+8.2f} {r['bs']:+8.2f}"
              f"      {u1}          {u2}")
    print("  H-b (gerçek b<0 ↔ GUE-boyalı P1 b>0): "
          f"ölçülen {out['son']['bp']:+.2f} / {out['P1']['bp']:+.2f}; "
          f"A·S öngörüsü {out['son']['bs']:+.2f} / {out['P1']['bs']:+.2f}; "
          f"M1 öngörüsü {out['son']['bm']:+.2f} / {out['P1']['bm']:+.2f}")


def T4_taban():
    print()
    print("=" * 78)
    print("T4 — TABAN DAYANIKLILIĞI (0.40 vs 0.52), örtüşen bantlar")
    print("=" * 78)
    for v, g in (("son", "g158"), ("son", "t1"), ("A4", "g158")):
        a = yukle(v, 0.4, g); b = yukle(v, 0.52, g)
        if a is None or b is None:
            continue
        A = {x["tau"]: x for x in bantlar(a)}
        Bd = {x["tau"]: x for x in bantlar(b)}
        ort = sorted(set(A) & set(Bd))
        if not ort:
            continue
        print(f"\n  {ETI.get(v,v)} / {g} — {len(ort)} örtüşen bant")
        print("  τ̄     δ(0.40)  δ(0.52)   A·S(0.40) A·S(0.52)  "
              "δ/A·S(0.40) δ/A·S(0.52)")
        for t in ort:
            x, y = A[t], Bd[t]
            print(f"  {t:.2f} {x['delta']:+8.4f} {y['delta']:+8.4f} "
                  f"{x['AS']['S_tam']:+10.4f} {y['AS']['S_tam']:+9.4f} "
                  f"{x['delta']/x['AS']['S_tam']:+11.3f} "
                  f"{y['delta']/y['AS']['S_tam']:+11.3f}")


if __name__ == "__main__":
    T0()
    for v in ("son", "A4", "keskin", "P1"):
        T1(v, 0.4, "t1")
    T1("son", 0.52, "t1")
    T1_tanim("son", 0.4, "t1")
    T1_tanim("A4", 0.4, "t1")
    T2("son", 0.4, "t1")
    T2("A4", 0.4, "t1")
    T2("keskin", 0.4, "t1")
    o = T3()
    T3_ref(o)
    T4_taban()
