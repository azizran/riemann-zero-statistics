# -*- coding: utf-8 -*-
"""
175d — K2: DONDURULMUŞ KURALIN UYGULANMASI ve AYRIŞIM DEFTERİ
==============================================================
Girdi: 175/ONKAYIT_K2.json (175a, zaman damgalı, sha256'lı),
       174/K1_<gaz>.json (162 makinesi; 174b + 175b),
       174/K3_<gaz>.json (DC kaçağı; 174d + 175c),
       172/G1.json + G2.json + G3.json.
YENİ ÖLÇÜM YOK — yalnız dondurulmuş kural uygulanır.
"""
import json
import math
import sys
from pathlib import Path

import numpy as np

SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S175, S174, S172 = SCR / "175", SCR / "174", SCR / "172"
G1 = json.load(open(S172 / "G1.json"))
G2 = json.load(open(S172 / "G2.json"))
G3 = json.load(open(S172 / "G3.json"))
ONK = json.load(open(S175 / "ONKAYIT_K2.json"))
KES5 = ["E060", "HA4", "K070", "K090", "Hkeskin"]     # τ̄_A artan sırada
LAM7 = ["L050", "L060", "L070", "L085", "Hkeskin", "L115", "L130"]
LV = np.array([.50, .60, .70, .85, 1.00, 1.15, 1.30])
BANT = [(0.40, 0.50), (0.50, 0.60), (0.60, 0.70), (0.70, 0.80), (0.80, 0.95)]

K1 = {p.stem[3:]: json.load(open(p)) for p in sorted(S174.glob("K1_*.json"))}
K3 = {p.stem[3:]: json.load(open(p)) for p in sorted(S174.glob("K3_*.json"))}
print("=" * 74)
print("175d — K2 YÜZLEŞME   [ön-kayıt %s  sha %s]"
      % (ONK["zaman"], ONK["sha256"][:16]))
print("=" * 74)
print("  K1 (162 makinesi) ölçülü: %s" % ", ".join(sorted(K1)))
print("  K3 (DC kaçağı)    ölçülü: %s" % ", ".join(sorted(K3)))
OUT = dict(onkayit=ONK["zaman"], sha=ONK["sha256"])

# ═════════════════════════════════════════════════════════════════════
# K2-1  P1: KÖPRÜ KESİM EKSENİNDE DE TUTUYOR MU?
# ═════════════════════════════════════════════════════════════════════
print("\n" + "=" * 74)
print("K2-1 (P1) — KÖPRÜ KİMLİĞİ kesim ekseninde:  R_η(162) ↔ π_E(172)")
print("=" * 74)
print("%-9s %10s %10s %9s | %10s %9s | %9s %9s"
      % ("gaz", "R_η", "π_E", "fark%", "ön-kayıt", "sapma%", "R_Ĉ", "R_X"))
P1 = {}
for g in KES5 + ["son"]:
    if g not in K1:
        continue
    R, piE = K1[g]["eta"]["R"], G1[g]["E"]["pi"]
    ong = ONK["kural"]["P1"]["ongoru"].get(g, float("nan"))
    P1[g] = dict(R=R, piE=piE, fark=100 * (R / piE - 1), ongoru=ong,
                 sapma=100 * (R / ong - 1) if ong == ong else float("nan"),
                 R_C=K1[g]["Chat"]["R"], R_X=K1[g]["Xtil"]["R"])
    print("%-9s %10.5f %10.5f %+8.2f%% | %10.5f %+8.2f%% | %9.5f %9.5f"
          % (g, R, piE, P1[g]["fark"], ong, P1[g]["sapma"],
             P1[g]["R_C"], P1[g]["R_X"]))
mx = max(abs(P1[g]["fark"]) for g in KES5 if g in P1)
print("  ön-kayıt P1: |R_η/π_E - 1| <= 2%% (kesim gazlari) => olculen maks "
      "%.2f%%  → **%s**" % (mx, "✓ GEÇTİ" if mx <= 2.0 else "✗ KÖPRÜ KIRIK"))
OUT["P1"] = dict(tablo=P1, maks_fark=mx, gecti=bool(mx <= 2.0))

# ═════════════════════════════════════════════════════════════════════
# K2-2  P2: GERÇEK GAZ KESİM AİLESİNİN İÇİNDE Mİ? (aile-dışı fazlanın adresi)
# ═════════════════════════════════════════════════════════════════════
print("\n" + "=" * 74)
print("K2-2 (P2) — AİLE-DIŞI FAZLANIN ADRESİ:  λ ailesi ↔ KESİM ailesi")
print("=" * 74)
TB = {g: ONK["kesim_ekseni"][g]["tau_bar"] for g in KES5}


def kesim_es(y):
    """τ̄_A koordinatında kesim ailesinin y'yi verdiği yer (tekdüze kolda)."""
    yy = np.array([y[g] for g in KES5])
    tt = np.array([TB[g] for g in KES5])
    d = np.diff(yy)
    v = y["son"]
    if np.all(d > 0) and yy[0] <= v <= yy[-1]:
        return float(np.interp(v, yy, tt)), True
    if np.all(d < 0) and yy[-1] <= v <= yy[0]:
        return float(np.interp(v, yy[::-1], tt[::-1])), True
    ic = bool(yy.min() <= v <= yy.max())
    return float("nan"), ic


NIC = {
    "R_η":          lambda g: K1[g]["eta"]["R"],
    "R_Ĉ":          lambda g: K1[g]["Chat"]["R"],
    "R_X":          lambda g: K1[g]["Xtil"]["R"],
    "m3":           lambda g: K1[g]["m3"]["olc"],
    "m3_çizgi":     lambda g: K1[g]["m3"]["cizgi"],
    "m3/m3_çizgi":  lambda g: K1[g]["m3"]["olc"] / K1[g]["m3"]["cizgi"],
    "skew(e1)":     lambda g: K1[g]["m3"]["skew_e1"],
    "skew(x1)":     lambda g: K1[g]["m3"]["skew_x1"],
    "skew(ds)":     lambda g: K1[g]["m3"]["skew_ds"],
    "skew(η_çiz)":  lambda g: K1[g]["eta"]["skew_cizgi"],
    "P_η":          lambda g: K1[g]["eta"]["P"],
    "Var(η)":       lambda g: K1[g]["eta"]["Var"],
    "çizgi payı η": lambda g: K1[g]["eta"]["pay"],
    # ort(E): ÖZDEŞLİKTEN (|ort E| = √(μ̂²_E·K_E/π_E)); 175c bunu altı gazda
    # 4.1e−13'te doğruladı, dolayısıyla λ merdiveni için de kullanılabilir
    # (K3 makinesi λ merdiveninde yalnız L085+Hkeskin'de koşmuştu).
    "ort(E)":       lambda g: math.sqrt(G1[g]["E"]["mu2"] * G1[g]["E"]["K"]
                                        / G1[g]["E"]["pi"]),
    "μ̂²_E":         lambda g: G1[g]["E"]["mu2"],
    "g_E":          lambda g: G1[g]["gE"],
    "g_X":          lambda g: G1[g]["gX"],
    "θ":            lambda g: G2["theta"][g],
    "M = g_E g_X²θ": lambda g: G1[g]["KAL"],
    "σ_ds":         lambda g: G1[g]["sigds"],
}
print("%-14s %10s | %-46s | %8s %8s %s"
      % ("nicelik", "son", "kesim ailesi (E060→Hkeskin)", "λ-ailesi",
         "τ̄_A^eş", "kesim ailesi"))
ATL = {}
for ad, f in NIC.items():
    try:
        y = {g: float(f(g)) for g in KES5 + ["son"]}
    except (KeyError, TypeError):
        continue
    if any(v != v for v in y.values()):
        continue
    tes, ic = kesim_es(y)
    # λ ailesi durumu (174'ün atlası ile aynı ölçüt)
    try:
        yl = np.array([float(f(g)) for g in LAM7])
        icl = bool(yl.min() <= y["son"] <= yl.max())
        lam_dur = "içeride" if icl else "**DIŞINDA**"
    except (KeyError, TypeError):
        yl, icl, lam_dur = None, None, "—"
    ATL[ad] = dict(son=y["son"], kesim={g: y[g] for g in KES5},
                   tau_es=tes, ic=ic, lam_ic=icl,
                   lam=[float(v) for v in yl] if yl is not None else None)
    print("%-14s %10.5f | %-46s | %8s %8s %s"
          % (ad, y["son"], " ".join("%8.4f" % y[g] for g in KES5), lam_dur,
             ("%.4f" % tes) if tes == tes else "  —  ",
             "**İÇERİDE**" if ic else "dışında"))
OUT["ATLAS"] = ATL
kurtarilan = [a for a, v in ATL.items() if v["lam_ic"] is False and v["ic"]]
kacan = [a for a, v in ATL.items() if v["lam_ic"] is False and not v["ic"]]
print("\n  λ ailesinin DIŞINDA olup KESİM ailesinin İÇİNDE olanlar: %s"
      % (", ".join(kurtarilan) or "—"))
print("  iki ailenin de dışında kalanlar: %s" % (", ".join(kacan) or "—"))
OUT["P2"] = dict(kurtarilan=kurtarilan, kacan=kacan)

# ═════════════════════════════════════════════════════════════════════
# K2-3  P3: DC KAÇAĞI ÖZDEŞLİK DENETİMİ
# ═════════════════════════════════════════════════════════════════════
print("\n" + "=" * 74)
print("K2-3 (P3) — DC KAÇAĞI: özdeşlik öngörüsünün denetimi")
print("=" * 74)
print("%-9s %12s %12s %10s %12s" % ("gaz", "ort(E) ölç", "ön-kayıt",
                                    "sapma", "μ̂²(ξ)/μ̂²(G1)−1"))
P3 = {}
for g in KES5 + ["son"]:
    if g not in K3:
        continue
    ong = ONK["kural"]["P3"]["ongoru"][g]
    P3[g] = dict(ortE=K3[g]["ortE"], ongoru=ong,
                 sapma=abs(K3[g]["ortE"] / ong - 1),
                 mu2=abs(K3[g]["mu2"] / K3[g]["mu2_G1"] - 1))
    print("%-9s %+12.6f %+12.6f %10.2e %12.2e"
          % (g, K3[g]["ortE"], ong, P3[g]["sapma"], P3[g]["mu2"]))
mxp = max(v["sapma"] for v in P3.values())
print("  ön-kayıt P3 (‰5, işaret +): maks sapma %.2e → **%s**"
      % (mxp, "✓ GEÇTİ" if mxp <= 5e-3 else "✗ ÖLDÜ"))
OUT["P3"] = dict(tablo=P3, maks=mxp, gecti=bool(mxp <= 5e-3))

# ═════════════════════════════════════════════════════════════════════
# K2-4  P4: MÜHÜR KURALI (KALEM, literal)
# ═════════════════════════════════════════════════════════════════════
def bant_xi(g, lo, hi):
    for b in K3[g]["bant"]:
        if abs(b["lo"] - lo) < 1e-9 and abs(b["hi"] - hi) < 1e-3:
            return b["xi"], b["n"], b["mod"], b["cos"]
    return float("nan"), 0, float("nan"), float("nan")


print("\n" + "=" * 74)
print("K2-4 (P4) — MÜHÜR: τ > 0.70 açığı erfc-ikizde ≥%50 küçülüyor mu?")
print("=" * 74)
hi70 = [(0.70, 0.80), (0.80, 0.95)]
a_son = sum(bant_xi("son", *b)[0] for b in hi70)
a_hk = sum(bant_xi("Hkeskin", *b)[0] for b in hi70)
a_er = sum(bant_xi("HA4", *b)[0] for b in hi70)
g_hk, g_er = a_son - a_hk, a_son - a_er
esik = ONK["kural"]["P4"]["esik"]
print("  Σξ(τ>0.70):  son %+.6e | keskin-ikiz %+.6e | erfc-ikiz %+.6e"
      % (a_son, a_hk, a_er))
print("  açık(gerçek−keskin) = %+.6e   açık(gerçek−erfc) = %+.6e"
      % (g_hk, g_er))
print("  |açık_erfc| / |açık_keskin| = %.3f   (mühür eşiği ≤ 0.50)"
      % (abs(g_er) / abs(g_hk)))
muhur = abs(g_er) <= esik
print("  ⇒ ön-kayıtlı MÜHÜR KURALI: **%s**"
      % ("H-K1 YAŞAR" if muhur else "H-K1'in literal biçimi ÖLDÜ "
         "(ön-kayıt bunu önceden yazmıştı)"))
print("  işaret: açık_keskin %s , açık_erfc %s  ⇒ %s"
      % ("−" if g_hk < 0 else "+", "−" if g_er < 0 else "+",
         "**SANDVİÇ: gerçek gaz iki ikizin ARASINDA**"
         if g_hk * g_er < 0 else "aynı yönde (sandviç YOK)"))
OUT["P4"] = dict(xi_son=a_son, xi_keskin=a_hk, xi_erfc=a_er,
                 acik_keskin=g_hk, acik_erfc=g_er,
                 oran=abs(g_er) / abs(g_hk), muhur=bool(muhur),
                 sandvic=bool(g_hk * g_er < 0))

# ═════════════════════════════════════════════════════════════════════
# K2-5  P5: KESİM KESRİ f — BANT BANT (kesim payı ↔ kilit payı)
# ═════════════════════════════════════════════════════════════════════
print("\n" + "=" * 74)
print("K2-5 (P5) — KESİM KESRİ  f = [y(son)−y(keskin)] / [y(erfc)−y(keskin)]")
print("=" * 74)
print("  DC kaçağı bant bant:")
print("   %-11s %6s %12s %12s %12s %9s" % ("τ bandı", "n", "Σξ(son)",
                                           "Σξ(keskin)", "Σξ(erfc)", "f_b"))
FB = {}
for lo, hi in BANT:
    xs, n, _, _ = bant_xi("son", lo, hi)
    xk, _, _, _ = bant_xi("Hkeskin", lo, hi)
    xe, _, _, _ = bant_xi("HA4", lo, hi)
    fb = (xs - xk) / (xe - xk) if xe != xk else float("nan")
    FB["%.2f-%.2f" % (lo, hi)] = dict(n=n, son=xs, keskin=xk, erfc=xe, f=fb)
    print("   %.2f–%.2f   %6d %+12.4e %+12.4e %+12.4e %+9.4f"
          % (lo, hi, n, xs, xk, xe, fb))
fv = np.array([v["f"] for v in FB.values()])
fw = np.array([abs(v["erfc"] - v["keskin"]) for v in FB.values()])
f_ag = float((fv * fw).sum() / fw.sum())
yay = float(fv.max() - fv.min())
print("   TOPLAM ort(E):  %+.6f  %+.6f  %+.6f   f = %+.4f"
      % (K3["son"]["ortE"], K3["Hkeskin"]["ortE"], K3["HA4"]["ortE"],
         (K3["son"]["ortE"] - K3["Hkeskin"]["ortE"]) /
         (K3["HA4"]["ortE"] - K3["Hkeskin"]["ortE"])))
print("   ağırlıklı f = %+.4f   bantlar arası YAYILIM = %.4f  (ön-kayıt "
      "eşiği ±0.15 ⇒ %s)" % (f_ag, yay,
                             "SAF KESİM" if yay <= 0.30 else
                             "**KESİM TEK BAŞINA TAŞIMIYOR**"))
print("\n  gözlenebilir bazında f (162 makinesi + 172 defteri):")
FY = {}
for ad, f in NIC.items():
    try:
        ys, yk, ye = float(f("son")), float(f("Hkeskin")), float(f("HA4"))
    except (KeyError, TypeError):
        continue
    if ys != ys or yk != yk or ye != ye or ye == yk:
        continue
    FY[ad] = (ys - yk) / (ye - yk)
    print("     f(%-13s) = %+8.4f   [son %.5f | keskin %.5f | erfc %.5f]"
          % (ad, FY[ad], ys, yk, ye))
OUT["P5"] = dict(bant=FB, f_agirlikli=f_ag, yayilim=yay, f_nicelik=FY)

# ═════════════════════════════════════════════════════════════════════
# K2-6  P6: ΔM AYRIŞIM DEFTERİ — kesim + kilit + θ
# ═════════════════════════════════════════════════════════════════════
print("\n" + "=" * 74)
print("K2-6 (P6) — ΔM AYRIŞIMI:  ΔM = kesim payı + kilit payı + θ payı")
print("=" * 74)
fyz = FY["R_η"]                       # f* — ÖLÇÜLEN kesim kesri (η kanalı)


def dlog(f, a, b):
    return math.log(f(a) / f(b))


gE = lambda g: G1[g]["gE"]
gX = lambda g: G1[g]["gX"]
th = lambda g: G1[g]["th"]
M = lambda g: G1[g]["KAL"]
dM = dlog(M, "son", "Hkeskin")
dE, dX2, dT = (dlog(gE, "son", "Hkeskin"), 2 * dlog(gX, "son", "Hkeskin"),
               dlog(th, "son", "Hkeskin"))
eM = dlog(M, "HA4", "Hkeskin")
eE, eX2, eT = (dlog(gE, "HA4", "Hkeskin"), 2 * dlog(gX, "HA4", "Hkeskin"),
               dlog(th, "HA4", "Hkeskin"))
kesim = fyz * eM
kilit = (dE + dX2) - fyz * (eE + eX2)
tpay = dT - fyz * eT
print("  f* = f(R_η) = %+.5f   (162 makinesinin ÖLÇTÜĞÜ kesim kesri)" % fyz)
print("  Δlog M(gerçek ← keskin ikiz) = %+.6f   (+%.2f%%)"
      % (dM, 100 * (math.exp(dM) - 1)))
print("    = Δlog g_E %+.6f + 2Δlog g_X %+.6f + Δlog θ %+.6f"
      % (dE, dX2, dT))
print("  Δlog M(erfc ikiz ← keskin ikiz) = %+.6f  = %+.6f + %+.6f + %+.6f"
      % (eM, eE, eX2, eT))
print("\n  %-24s %12s %10s" % ("pay", "Δlog", "% (ΔM'nin)"))
for ad, v in (("KESİM payı  (f*·erfc yolu)", kesim),
              ("KİLİT payı  (E,X artığı)", kilit),
              ("θ payı      (θ artığı)", tpay)):
    print("  %-24s %+12.6f %+10.1f%%" % (ad, v, 100 * v / dM))
print("  %-24s %+12.6f %+10.1f%%   (kalıntı %.1e)"
      % ("TOPLAM", kesim + kilit + tpay, 100 * (kesim + kilit + tpay) / dM,
         abs(kesim + kilit + tpay - dM)))
# ayrıca: kanal kanal kesim payları
print("\n  kanal kanal (aynı f* ile):")
for ad, dv, ev in (("g_E", dE, eE), ("2·g_X", dX2, eX2), ("θ", dT, eT)):
    print("     %-6s Δlog(son)=%+.6f  kesim=%+.6f  artık=%+.6f  "
          "(kesimin kapattığı %+.0f%%)"
          % (ad, dv, fyz * ev, dv - fyz * ev,
             100 * fyz * ev / dv if dv else float("nan")))
OUT["P6"] = dict(f_yildiz=fyz, dM=dM, dE=dE, dX2=dX2, dT=dT,
                 eM=eM, eE=eE, eX2=eX2, eT=eT,
                 kesim=kesim, kilit=kilit, theta=tpay,
                 kalinti=abs(kesim + kilit + tpay - dM))

# ═════════════════════════════════════════════════════════════════════
# K2-7  P8: θ DEFTERİ (K4)
# ═════════════════════════════════════════════════════════════════════
print("\n" + "=" * 74)
print("K2-7 (P8 = K4) — θ DEFTERİ: keskin ↔ erfc ↔ gerçek")
print("=" * 74)
print("  %-9s %9s %9s %9s %9s %9s %9s"
      % ("gaz", "θ(G1)", "θ(G2)", "g_cal", "g_E", "g_X", "M"))
for g in KES5 + ["son"]:
    print("  %-9s %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f"
          % (g, G1[g]["th"], G2["theta"][g], G1[g]["gcal"], G1[g]["gE"],
             G1[g]["gX"], G1[g]["KAL"]))
fth1 = FY.get("θ", float("nan"))
print("  f(θ) = %+.4f  → **%s**"
      % (fth1, "H-K4 ÖLDÜ: kesim ekseni θ'yı TERS yöne taşıyor"
         if fth1 < 0 else "kesim θ'yı doğru yöne taşıyor"))
OUT["P8"] = dict(f_theta=fth1,
                 tablo={g: dict(th1=G1[g]["th"], th2=G2["theta"][g],
                                gcal=G1[g]["gcal"], gE=G1[g]["gE"],
                                gX=G1[g]["gX"], M=G1[g]["KAL"])
                        for g in KES5 + ["son"]})

json.dump(OUT, open(S175 / "K2.json", "w"), indent=1, ensure_ascii=False)
print("\n-> %s" % (S175 / "K2.json"))
