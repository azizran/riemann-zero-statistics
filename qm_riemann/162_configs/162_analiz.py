"""
162 — ANALİZ
B1  DENETİM: `gercek` koşusu 160'ın kayıtlı JSON'unu bit düzeyinde veriyor mu
B2  T2 ANA TABLO: bant bant δ_gerçek / δ_vekil (+ K1, K2, A·S, K2/δ)
B3  ÇEŞİT AYRIŞIMI: V0 / VL (çizgi fazı) / VC (süreklilik) — hangisi taşıyor
B4  ⟨σX̃²⟩ SINAVI: A²s2/2 ve K2 — Gauss vekilde sıfıra düşüyor mu
B5  a ve b (g158 ızgarası, 158'in 4 penceresi) — b'nin işareti vekilde
B6  İKİNCİ-MOMENT BÜTÇESİ: vekil gerçekten aynı ikinci momentlere mi sahip
B7  A·σ_X EŞLEŞTİRİLMİŞ karşılaştırma (vekilin σ_X'i şiştiği için)
B8  ÖZET
"""
import importlib
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
C = importlib.import_module("162_cekirdek")
sys.path.insert(0, str(C.QM / "160_configs"))
A160 = importlib.import_module("160_analiz")

SCR = C.SCR
SCR160 = A160.SCR
TWO_PI = 2 * np.pi
FOUR_PI = 4 * np.pi
GAZ = ("son", "keskin", "A4")
ETI = {"son": "gerçek-son", "keskin": "keskin", "A4": "A4"}
VG = ("VG101", "VG102", "VG103", "VG104", "VG105")


def yukle(veri, cesit, g="t1", taban=0.4):
    p = SCR / f"S_{veri}_t{taban}_{g}_{cesit}.json"
    return json.load(open(p)) if p.exists() else None


def vg_listesi(veri, g="t1"):
    return [c for c in VG if yukle(veri, c, g) is not None]


def bant_dizisi(d, ok_sadece=True):
    B = A160.bantlar(d)
    if ok_sadece:
        s = A160.saglik(B)
        B = [b for b, o in zip(B, s) if o]
    return B


def eslesen(veri, cesitler, g="t1"):
    """Aynı τ'daki bantları eşleştirir; GERÇEĞİN sağlık kuralını uygular."""
    dd = {c: yukle(veri, c, g) for c in cesitler}
    if any(v is None for v in dd.values()):
        return None, None
    ref = bant_dizisi(dd["gercek"])
    tt = [b["tau"] for b in ref]
    out = {}
    for c, d in dd.items():
        ix = {b["tau"]: b for b in A160.bantlar(d)}
        out[c] = [ix.get(t) for t in tt]
    return tt, out


def _ort(vals):
    v = np.array([x for x in vals if x is not None and np.isfinite(x)], float)
    if len(v) == 0:
        return float("nan"), float("nan")
    return float(v.mean()), float(v.std(ddof=1) if len(v) > 1 else 0.0)


# ------------------------------------------------------------------ B1
def B1():
    print("\n" + "=" * 78)
    print("B1. DENETİM — 162'nin `gercek` koşusu 160'ın kaydını yeniden "
          "üretiyor mu?")
    print("=" * 78)
    kol = ["phi", "delta", "tau_eff", "absG", "K1", "K2", "Kfull", "M1", "M0",
           "sPhi_jk", "sDelta_jk", "sK1_jk", "sKf_jk", "rho", "absE1"]
    print("  gaz      ızgara  ortak bant   maks fark (bütün sütunlar)")
    for v in GAZ:
        for g in ("t1", "g158"):
            d = yukle(v, "gercek", g)
            p160 = SCR160 / f"S_{v}_t0.4_{g}.json"
            if d is None or not p160.exists():
                continue
            r = json.load(open(p160))
            ix = {b["tau"]: b for b in r["bantlar"] if b.get("olculdu")}
            mx, n = 0.0, 0
            for b in A160.bantlar(d):
                o = ix.get(b["tau"])
                if o is None:
                    continue
                n += 1
                for k in kol:
                    if k in b and k in o and np.isfinite(b[k]) and np.isfinite(o[k]):
                        mx = max(mx, abs(b[k] - o[k]))
                for k in ("S_tam", "S_xtil"):
                    mx = max(mx, abs(b["AS"][k] - o["AS"][k]))
            print(f"  {ETI[v]:<10} {g:<6} {n:>6}        {mx:.3e}")


# ------------------------------------------------------------------ B2
def B2(g="t1", lo=0.50, hi=0.82):
    print("\n" + "=" * 78)
    print(f"B2. T2 — GAUSS-TAYFSAL KAPANIŞ SINAVI ({g}, τ∈({lo},{hi}))")
    print("=" * 78)
    ozet = {}
    for v in GAZ:
        vgl = vg_listesi(v, g)
        tt, B = eslesen(v, ["gercek", "V0"] + vgl, g)
        if tt is None:
            continue
        print(f"\n  --- {ETI[v]} --- (vekil tohumları: {', '.join(vgl)})")
        print("   τ_eff  A·σ_X  |Γ|    δ_gerçek  δ_V0      δ_vekil(ort±sd)   "
              "  δ_g/δ_v   K1_g/K1_v  K2/δ g   K2/δ v   |Γ|_v")
        satir = []
        for i, t in enumerate(tt):
            bg = B["gercek"][i]
            if bg is None or not (lo <= bg["tau_eff"] <= hi):
                continue
            b0 = B["V0"][i]
            dv = [B[c][i]["delta"] for c in vgl if B[c][i] is not None]
            k1v = [B[c][i]["K1"] for c in vgl if B[c][i] is not None]
            k2v = [B[c][i]["K2"] for c in vgl if B[c][i] is not None]
            gv = [B[c][i]["absG"] for c in vgl if B[c][i] is not None]
            axv = [TWO_PI * B[c][i]["tau_eff"] * np.sqrt(B[c][i]["sA2"])
                   for c in vgl if B[c][i] is not None]
            m, s = _ort(dv)
            m1, _ = _ort(k1v)
            m2, _ = _ort(k2v)
            mg, _ = _ort(gv)
            mx, _ = _ort(axv)
            aX = TWO_PI * bg["tau_eff"] * np.sqrt(bg["sA2"])
            print(f"   {bg['tau_eff']:.4f} {aX:5.3f} {bg['absG']:.3f} "
                  f"{bg['delta']:+9.4f} {b0['delta'] if b0 else np.nan:+9.4f} "
                  f"{m:+8.4f}±{s:.4f}  {bg['delta']/m if m else np.nan:+8.2f}  "
                  f"{bg['K1']/m1 if m1 else np.nan:+9.2f}  "
                  f"{bg['K2']/bg['delta']:+7.3f}  {m2/m if m else np.nan:+7.3f}  "
                  f"{mg:.3f}")
            satir.append(dict(tau=bg["tau_eff"], aX=aX, aXv=mx,
                              dg=bg["delta"], d0=b0["delta"] if b0 else None,
                              dv=m, dvs=s, k1g=bg["K1"], k1v=m1,
                              k2g=bg["K2"], k2v=m2, Gg=bg["absG"], Gv=mg,
                              ASg=bg["AS"]["S_xtil"],
                              ASv=_ort([B[c][i]["AS"]["S_xtil"]
                                        for c in vgl if B[c][i]])[0]))
        if satir:
            o = np.array([r["dg"] / r["dv"] for r in satir if r["dv"]])
            o0 = np.array([r["dg"] / r["d0"] for r in satir
                           if r["d0"]])
            print(f"   ORAN δ_gerçek/δ_vekil: medyan {np.median(o):+.2f}  "
                  f"menzil [{o.min():+.2f}, {o.max():+.2f}]   "
                  f"|oran−1| medyan {np.median(np.abs(o-1))*100:.0f} %")
            print(f"   KONTROL δ_gerçek/δ_V0 : medyan {np.median(o0):+.5f}  "
                  f"menzil [{o0.min():+.5f}, {o0.max():+.5f}]  "
                  f"(1.0 olmalı — boru hattı kimliği)")
            ozet[v] = satir
    return ozet


# ------------------------------------------------------------------ B3
def B3(g="t1", lo=0.50, hi=0.82):
    print("\n" + "=" * 78)
    print("B3. ÇEŞİT AYRIŞIMI — δ'yı hangi faz yapısı taşıyor?")
    print("    VL = yalnız ÇİZGİ fazları karıştırıldı (süreklilik gerçek)")
    print("    VC = yalnız SÜREKLİLİK karıştırıldı (çizgi fazları gerçek)")
    print("=" * 78)
    for v in GAZ:
        cl = [c for c in ("V0", "VL101", "VC101", "VG101")
              if yukle(v, c, g) is not None]
        tt, B = eslesen(v, ["gercek"] + cl, g)
        if tt is None:
            continue
        print(f"\n  --- {ETI[v]} ---")
        print("   τ_eff   δ_gerçek   " + "  ".join(f"δ_{c:<7}" for c in cl)
              + "   || oran gerçek/·")
        for i, t in enumerate(tt):
            bg = B["gercek"][i]
            if bg is None or not (lo <= bg["tau_eff"] <= hi):
                continue
            vs = [B[c][i]["delta"] if B[c][i] else np.nan for c in cl]
            print(f"   {bg['tau_eff']:.4f} {bg['delta']:+10.4f}   " +
                  "  ".join(f"{x:+9.4f}" for x in vs) + "   || " +
                  "  ".join(f"{bg['delta']/x if x else np.nan:+6.2f}"
                            for x in vs))


# ------------------------------------------------------------------ B4
def B4(g="t1", lo=0.50, hi=0.82):
    print("\n" + "=" * 78)
    print("B4. ⟨σX̃²⟩ SINAVI — K2'yi taşıyan ÜÇÜNCÜ-MOMENT nesnesi")
    print("    s2 = ⟨σX̃0²⟩/⟨ρ⟩ ; Re[Es]'in ilk X-bağımlı terimi −A²s2/2")
    print("=" * 78)
    for v in GAZ:
        vgl = vg_listesi(v, g)
        tt, B = eslesen(v, ["gercek", "V0"] + vgl, g)
        if tt is None:
            continue
        print(f"\n  --- {ETI[v]} ---")
        print("   τ_eff   −A²s2/2 gerçek   V0        vekil(ort±sd)      "
              "oran g/v   K2 gerçek   K2 vekil    işaret")
        for i, t in enumerate(tt):
            bg = B["gercek"][i]
            if bg is None or not (lo <= bg["tau_eff"] <= hi):
                continue
            sg = -0.5 * bg["KUM"]["A2s2"]
            s0 = -0.5 * B["V0"][i]["KUM"]["A2s2"] if B["V0"][i] else np.nan
            sv = [-0.5 * B[c][i]["KUM"]["A2s2"] for c in vgl if B[c][i]]
            m, s = _ort(sv)
            k2v, _ = _ort([B[c][i]["K2"] for c in vgl if B[c][i]])
            isr = "AYNI" if np.sign(sg) == np.sign(m) else "TERS"
            print(f"   {bg['tau_eff']:.4f}  {sg:+.6f}     {s0:+.6f}  "
                  f"{m:+.6f}±{s:.6f}  {sg/m if m else np.nan:+8.2f}   "
                  f"{bg['K2']:+.5f}   {k2v:+.5f}   {isr}")


# ------------------------------------------------------------------ B5
def B5():
    print("\n" + "=" * 78)
    print("B5. a ve b — 158/160'ın fit konvansiyonu (g158, 4 pencere)")
    print("    a = dφ/dτ|τ₀ , b = ½d²φ/dτ²|τ₀ ; b(δ)=b(K1+K2), b(K1), b(A·S)")
    print("=" * 78)
    for v in GAZ:
        vgl = vg_listesi(v, "g158")
        print(f"\n  --- {ETI[v]} --- (vekiller: {', '.join(vgl) or 'yok'})")
        print("   çeşit    τ₀      a(φ)      b(φ)      b(δ)      b(K1)     "
              "b(A·S)    4π+dδ/dτ  n")
        for c in ["gercek", "V0"] + vgl:
            d = yukle(v, c, "g158")
            if d is None:
                continue
            F = [A160.fit_seti(d, p) for p in ("W-A", "W-D", "W-B", "W-C")]
            F = [f for f in F if f]
            if not F:
                print(f"   {c:<8} fit yok")
                continue
            def ms(k):
                a = np.array([f[k] for f in F], float)
                return a.mean(), (a.std(ddof=1) if len(a) > 1 else 0.0)
            t0 = np.mean([f["t0"] for f in F])
            am, asd = ms("a"); bm, bsd = ms("b")
            bd, _ = ms("bw_K1+K2"); bk, _ = ms("bw_K1"); ba, _ = ms("bw_A·S")
            dd, _ = ms("d_K1+K2")
            print(f"   {c:<8} {t0:.4f} {am:8.3f}  {bm:+8.3f}  {bd:+8.3f}  "
                  f"{bk:+8.3f}  {ba:+8.3f}  {FOUR_PI+dd:8.3f}  {len(F)}")


# ------------------------------------------------------------------ B6
def B6(g="t1"):
    print("\n" + "=" * 78)
    print("B6. İKİNCİ-MOMENT BÜTÇESİ — vekil gerçekten aynı mı?")
    print("=" * 78)
    print("   gaz     çeşit    σΔ²=Var(dsΔ)  Var(η)     ⟨X̃⟩       "
          "çizgi-payı η / Ĉ   varyans kapanışı η / Ĉ")
    for v in GAZ:
        for c in ["gercek", "V0"] + vg_listesi(v, g) + ["VL101", "VC101"]:
            d = yukle(v, c, g)
            if d is None:
                continue
            tb = d.get("taban162", {})
            kp = tb.get("kapanis", [np.nan, np.nan])
            print(f"   {v:<7} {c:<8} {d['sA2']:.6f}      {d['s_eta']:.6f}  "
                  f"{d['Xort']:+.3e}   {tb.get('pay_eta', np.nan):+.4f} / "
                  f"{tb.get('pay_C', np.nan):+.4f}    "
                  f"{kp[0]:.4f} / {kp[1]:.4f}")


# ------------------------------------------------------------------ B7
def B7(g="t1", lo=0.50, hi=0.82):
    print("\n" + "=" * 78)
    print("B7. A·σ_X EŞLEŞTİRİLMİŞ KARŞILAŞTIRMA")
    print("    Vekilin σΔ²'si şişiyor; 160 δ'nın A·σ_X'e güçlü bağlı olduğunu")
    print("    ölçmüştü. Bu tablo, GERÇEĞİN δ(A·σ_X) eğrisini vekilin KENDİ")
    print("    A·σ_X'inde değerlendirir — şişme tek başına açıklıyor mu?")
    print("=" * 78)
    for v in GAZ:
        vgl = vg_listesi(v, g)
        tt, B = eslesen(v, ["gercek"] + vgl, g)
        if tt is None:
            continue
        R = []
        for i in range(len(tt)):
            bg = B["gercek"][i]
            if bg is None:
                continue
            R.append((TWO_PI * bg["tau_eff"] * np.sqrt(bg["sA2"]), bg["delta"]))
        R.sort()
        xr = np.array([r[0] for r in R]); yr = np.array([r[1] for r in R])
        print(f"\n  --- {ETI[v]} --- (gerçek eğri A·σ_X ∈ "
              f"[{xr.min():.3f},{xr.max():.3f}])")
        print("   τ_eff   A·σ_X_v  δ_vekil    δ_gerçek(A·σ_X_v)  oran   "
              "δ_gerçek(kendi τ)  ham oran")
        for i in range(len(tt)):
            bg = B["gercek"][i]
            if bg is None or not (lo <= bg["tau_eff"] <= hi):
                continue
            av = _ort([TWO_PI * B[c][i]["tau_eff"] * np.sqrt(B[c][i]["sA2"])
                       for c in vgl if B[c][i]])[0]
            dv = _ort([B[c][i]["delta"] for c in vgl if B[c][i]])[0]
            dg_i = float(np.interp(av, xr, yr))
            print(f"   {bg['tau_eff']:.4f} {av:7.3f}  {dv:+9.4f}  "
                  f"{dg_i:+13.4f}      {dg_i/dv if dv else np.nan:+6.2f}   "
                  f"{bg['delta']:+13.4f}      "
                  f"{bg['delta']/dv if dv else np.nan:+6.2f}")


# ------------------------------------------------------------------ B8
def B8(g="t1", lo=0.50, hi=0.82):
    print("\n" + "=" * 78)
    print("B8. ÖZET — T2/T3 hükmü")
    print("=" * 78)
    print("   gaz       bant  δ_g/δ_V0 (kontrol)  δ_g/δ_VG medyan [menzil]  "
          " K1_g/K1_VG  s2_g/s2_VG  δ_g/δ_VL  δ_g/δ_VC")
    for v in GAZ:
        vgl = vg_listesi(v, g)
        cl = ["gercek", "V0"] + vgl + [c for c in ("VL101", "VC101")
                                       if yukle(v, c, g)]
        tt, B = eslesen(v, cl, g)
        if tt is None:
            continue
        o, o0, ok1, os2, ol, oc = [], [], [], [], [], []
        for i in range(len(tt)):
            bg = B["gercek"][i]
            if bg is None or not (lo <= bg["tau_eff"] <= hi):
                continue
            m, _ = _ort([B[c][i]["delta"] for c in vgl if B[c][i]])
            m1, _ = _ort([B[c][i]["K1"] for c in vgl if B[c][i]])
            m2, _ = _ort([B[c][i]["KUM"]["A2s2"] for c in vgl if B[c][i]])
            if m:
                o.append(bg["delta"] / m)
            if B["V0"][i] and B["V0"][i]["delta"]:
                o0.append(bg["delta"] / B["V0"][i]["delta"])
            if m1:
                ok1.append(bg["K1"] / m1)
            if m2:
                os2.append(bg["KUM"]["A2s2"] / m2)
            for ad, lst in (("VL101", ol), ("VC101", oc)):
                if ad in B and B[ad][i] and B[ad][i]["delta"]:
                    lst.append(bg["delta"] / B[ad][i]["delta"])
        def q(a):
            a = np.array(a, float)
            return f"{np.median(a):+6.2f}" if len(a) else "   —  "
        o = np.array(o, float)
        print(f"   {ETI[v]:<9} {len(o):>4}   "
              f"{np.median(o0) if len(o0) else np.nan:+.5f}          "
              f"{np.median(o):+6.2f} [{o.min():+.2f},{o.max():+.2f}]   "
              f"{q(ok1)}      {q(os2)}      {q(ol)}    {q(oc)}")


if __name__ == "__main__":
    B1()
    B2()
    B3()
    B4()
    B5()
    B6()
    B7()
    B8()
