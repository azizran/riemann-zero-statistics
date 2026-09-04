# -*- coding: utf-8 -*-
"""
174f — EK TANILAR (ÖN-MÜHÜRSÜZ, İKİNCİ TUR) + FİGÜR
=====================================================
Bu betikteki her sayı ölçümden SONRA sorulmuş sorulardır; hiçbiri
174a'nın dondurulmuş kuralında yoktu ve hiçbiri bir ölümü kurtarmaz.
Amaç: K2'nin ölen maddelerinin ARDINDAN neyin ayakta kaldığını görmek.

E1  λ_eş ATLASI — K1'in ürettiği bütün nicelikler için 172e'nin ters
    çevirmesiyle λ_eş (7-nokta merdiveni).
E2  θ'nın DOĞRU vekili: θ tanımı gereği bir ORANDIR (KALİB = ölçülen /
    model). K2-E'nin ön-kayıtlı vekili yalnız ÖLÇÜLEN üçüncü momentti
    (m3) ve öldü. Oran vekili m3/m3_çizgi denenir.
E3  FİGÜR (174_gercek_imza.png).
"""
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S174, S172 = SCR / "174", SCR / "172"
QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
LAM7 = ["L050", "L060", "L070", "L085", "Hkeskin", "L115", "L130"]
LV = np.array([.50, .60, .70, .85, 1.00, 1.15, 1.30])
G1 = json.load(open(S172 / "G1.json"))
G2 = json.load(open(S172 / "G2.json"))
G3 = json.load(open(S172 / "G3.json"))
K1 = {json.load(open(p))["veri"]: json.load(open(p))
      for p in S174.glob("K1_*.json")}
K3 = {json.load(open(p))["veri"]: json.load(open(p))
      for p in S174.glob("K3_*.json")}


def lam_es(y, v):
    for kol in (slice(0, 6), slice(0, 5), slice(None)):
        yy, ll = y[kol], LV[kol]
        d = np.diff(yy)
        if np.all(d > 0):
            if yy[0] <= v <= yy[-1]:
                return float(np.interp(v, yy, ll)), True
        elif np.all(d < 0):
            if yy[-1] <= v <= yy[0]:
                return float(np.interp(v, yy[::-1], ll[::-1])), True
    return float("nan"), False


NIC = {
 "R_η (girişim, η)": lambda d: d["eta"]["R"],
 "R_Ĉ (girişim, Ĉ)": lambda d: d["Chat"]["R"],
 "R_X (girişim, X̃)": lambda d: d["Xtil"]["R"],
 "P_η = Σ|c|²/2": lambda d: d["eta"]["P"],
 "Var(η)": lambda d: d["eta"]["Var"],
 "Var(η_çizgi)": lambda d: d["eta"]["Var_cizgi"],
 "çizgi payı η": lambda d: d["eta"]["pay"],
 "çizgi payı X̃": lambda d: d["Xtil"]["pay"],
 "bantlar-arası Kov(η)": lambda d: d["eta"]["capraz"],
 "iç-bant ΣV_b(η)": lambda d: d["eta"]["ic_bant"],
 "Kov(η_çiz,η_artık)": lambda d: d["eta"]["P"] - d["eta"]["Var_cizgi"],
 "m3 = ⟨e1x1²⟩/σσ²": lambda d: d["m3"]["olc"],
 "m3_çizgi (model)": lambda d: d["m3"]["cizgi"],
 "**m3/m3_çizgi**": lambda d: d["m3"]["olc"] / d["m3"]["cizgi"],
 "skew(e1)": lambda d: d["m3"]["skew_e1"],
 "skew(η_çizgi)": lambda d: d["eta"]["skew_cizgi"],
 "skew(x1)": lambda d: d["m3"]["skew_x1"],
 "skew(ds)": lambda d: d["m3"]["skew_ds"],
 "σ_ds": lambda d: d["sigma"]["ds"],
 "σ_Ĉ": lambda d: d["sigma"]["Chat"],
}

print("=" * 78)
print("174f / E1 — λ_eş ATLASI (ÖN-MÜHÜRSÜZ ikinci tur)")
print("=" * 78)
print("  hedefler (172e): λ_eş(E)=0.7768  λ_eş(X)=0.9089  λ_eş(marj)=0.9354"
      "  λ_eş(θ)=0.6916  çapa=0.9363")
print("%-22s %10s | %-52s | %9s %s"
      % ("nicelik", "son", "λ merdiveni", "λ_eş", "durum"))
ATL = {}
for ad, f in NIC.items():
    y = np.array([f(K1[g]) for g in LAM7])
    v = f(K1["son"])
    le, ok = lam_es(y, v)
    d = np.diff(y)
    tek = "↑" if np.all(d[:5] > 0) else ("↓" if np.all(d[:5] < 0) else "∩∪")
    ic = y.min() <= v <= y.max()
    ATL[ad] = dict(son=v, lam=le, merdiven=[float(x) for x in y],
                   tekduze=tek, icerde=bool(ic))
    print("%-22s %10.5f | %-52s | %9s %s%s"
          % (ad, v, " ".join("%8.4f" % x for x in y),
             ("%.4f" % le) if le == le else "  YOK  ", tek,
             "" if ic else "  **ARALIK DIŞI (%s)**"
             % ("üstünde" if v > y.max() else "altında")))

print("\n" + "=" * 78)
print("174f / E2 — θ'nın ORAN vekili (K2-E'nin ölümünden sonra)")
print("=" * 78)
print("  θ tanımı: θ = KALİB_u2/g_cal , KALİB = ⟨ölçülen 3. moment⟩ /"
      " ⟨model 3. moment⟩")
print("  ön-kayıtlı vekil m3 (yalnız ölçülen)  → λ_eş = %s  [ÖLDÜ]"
      % (("%.4f" % ATL["m3 = ⟨e1x1²⟩/σσ²"]["lam"])
         if ATL["m3 = ⟨e1x1²⟩/σσ²"]["lam"] ==
         ATL["m3 = ⟨e1x1²⟩/σσ²"]["lam"] else "YOK (aralık dışı, üstünde)"))
r = ATL["**m3/m3_çizgi**"]
print("  oran vekili m3/m3_çizgi              → λ_eş = %s   (θ hedefi "
      "0.6916, fark %+.4f)"
      % (("%.4f" % r["lam"]) if r["lam"] == r["lam"] else "YOK",
         r["lam"] - 0.6916 if r["lam"] == r["lam"] else float("nan")))
th = np.array([G2["theta"][g] for g in LAM7])
print("  θ merdiveni: " + " ".join("%.4f" % x for x in th)
      + "   θ(son) = %.4f" % G2["theta"]["son"])
print("  [not] oran vekili θ ile AYNI YÖNDE (λ < çapa) ama nicel ıska;"
      " ön-mühürsüz olduğu için hiçbir ölümü kurtarmaz.")

# ---------------------------------------------------------------- figür
fig, ax = plt.subplots(2, 2, figsize=(13.5, 9.6))
c_son, c_ik, c_l = "#c1121f", "#1d3557", "#2a9d8f"

# (a) R_η(λ) tümseği
a = ax[0, 0]
y = np.array([K1[g]["eta"]["R"] for g in LAM7])
a.plot(LV, y, "o-", color=c_ik, lw=2, ms=6, label=r"$R_\eta(\lambda)$ ailesi")
a.axhline(K1["son"]["eta"]["R"], color=c_son, lw=2, ls="--",
          label=r"gerçek $\zeta$: $R_\eta=1.2883$")
a.axhline(max(y), color="0.6", lw=1, ls=":")
a.text(0.03, 0.62, "aile tavanı %.4f (λ=0.70)\ngerçek gaz %.2f%% üstünde"
       % (max(y), 100 * (K1["son"]["eta"]["R"] / max(y) - 1)),
       transform=a.transAxes, fontsize=9, color="0.35", va="top")
a.plot([1.0], [K1["Hkeskin"]["eta"]["R"]], "s", ms=10, mfc="none",
       mec=c_ik, mew=2)
a.annotate("sadakatli ikiz", (1.0, K1["Hkeskin"]["eta"]["R"]),
           textcoords="offset points", xytext=(-12, 12), fontsize=9,
           color=c_ik)
a.set_xlabel(r"$\lambda$ (merdiven genlik ölçeği)")
a.set_ylabel(r"$R_\eta=\Sigma|c_q(\eta)|^2/2 \;\div\; \mathrm{Var}(\eta)$")
a.set_title("(a) η-kanalı girişim oranı: gerçek gaz AİLENİN DIŞINDA\n"
            r"$\lambda_{eş}$ YOK — K2-A öldü", fontsize=11)
a.legend(fontsize=9, loc="lower left")
a.grid(alpha=.3)

# (b) R_X ve R_Ĉ tekdüze
b = ax[0, 1]
for ad, key, col, hed in (("$R_X$ (X̃=ΔĈ)", "Xtil", c_l, 0.9089),
                          ("$R_{\\hat C}$ (bellek)", "Chat", "#e76f51", None)):
    yy = np.array([K1[g][key]["R"] for g in LAM7])
    b.plot(LV, yy, "o-", color=col, lw=2, ms=5, label=ad)
    v = K1["son"][key]["R"]
    le, _ = lam_es(yy, v)
    b.axhline(v, color=col, ls="--", lw=1.2, alpha=.8)
    b.plot([le], [v], "*", ms=18, color=c_son, zorder=5)
    b.annotate(r"$\lambda_{eş}=%.4f$" % le, (le, v), textcoords="offset points",
               xytext=(8, 8), fontsize=10, color=c_son)
b.axvline(0.9089, color="0.5", ls=":", lw=1.5)
b.annotate("172e:\nλ_eş(X)=0.9089", (0.9089, 0.955), fontsize=9,
           color="0.35", ha="left", va="bottom")
b.set_xlabel(r"$\lambda$")
b.set_ylabel("girişim oranı")
b.set_title("(b) X̃ kanalı TEKDÜZE ve gerçek gaz İÇERİDE:\n"
            r"$\lambda_{eş}(R_X)=0.9026$ ↔ defterin 0.9089'u — K2-G tam isabet",
            fontsize=11)
b.legend(fontsize=9)
b.grid(alpha=.3)

# (c) DC kaçağının bant defteri
c = ax[1, 0]
if "son" in K3 and "Hkeskin" in K3:
    B = [0.40, 0.50, 0.60, 0.70, 0.80, 0.95]
    mid = [0.5 * (B[i] + B[i + 1]) for i in range(len(B) - 1)]
    ss, tt = [], []
    for i in range(len(B) - 1):
        for gz, acc in (("son", ss), ("Hkeskin", tt)):
            tau = np.array(K3[gz]["tau"])
            xi = np.array(K3[gz]["xi"])
            m = (tau > B[i]) & (tau <= B[i + 1])
            acc.append(float(xi[m].sum()))
    w = 0.035
    c.bar(np.array(mid) - w / 2, ss, w, color=c_son, label="gerçek ζ")
    c.bar(np.array(mid) + w / 2, tt, w, color=c_ik, label="sadakatli ikiz")
    c.axhline(0, color="k", lw=.8)
    c.set_xlabel(r"$\tau = \log q / L$")
    c.set_ylabel(r"$\sum_{q\in bant}\xi_q$,   $\xi_q=\mathrm{Re}[h_q\kappa(\omega_q)]$")
    c.set_title("(c) DC kaçağının bant defteri: işaret τ=0.5'te dönüyor,\n"
                "gerçek−ikiz açığının %83'ü τ>0.70'te", fontsize=11)
    c.legend(fontsize=9)
    c.grid(alpha=.3, axis="y")

# (d) κ yasası
d_ = ax[1, 1]
for gz, col, mk in (("son", c_son, "o"), ("Hkeskin", c_ik, "s"),
                    ("L085", c_l, "^")):
    if gz not in K3:
        continue
    tau = np.array(K3[gz]["tau"])
    ka = np.array(K3[gz]["kap_abs"])
    dd = json.load(open(S174 / f"K3_{gz}.json"))
    kc = np.abs(np.array(dd["kap_cip"]))            # πA τ  (çıplak)
    e = np.arange(0.40, 0.96, 0.05)
    xs, ys = [], []
    for lo, hi in zip(e[:-1], e[1:]):
        m = (tau > lo) & (tau <= hi)
        if m.sum() < 3:
            continue
        xs.append(0.5 * (lo + hi))
        ys.append(float(np.median(ka[m] / kc[m])))
    d_.plot(xs, ys, mk + "-", color=col, lw=1.6, ms=5, label=gz)
d_.axvline(0.5, color="0.5", ls=":", lw=1.5)
d_.annotate("τ=0.5: κ işaret değiştirir\n(orta-nokta tarağı, cos πτ)",
            (0.55, 0.05), fontsize=9, color="0.35")
d_.set_xlabel(r"$\tau$")
d_.set_ylabel(r"$|\kappa|_{ölç}\;/\;\pi A_q\tau_q$  (medyan)")
d_.set_title("(d) tarak rezonansı: çıplak doğrusal tepkiye oran\n"
             "τ<0.5'te cos πτ çukuru, τ>0.5'te çarpımsal kanalın fazlası",
             fontsize=11)
d_.legend(fontsize=9)
d_.grid(alpha=.3)

fig.suptitle("174 — GERÇEĞİN İMZASI: girişim oranının λ-ailesi ve DC "
             "kaçağının anatomisi", fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.965])
p = QM / "174_gercek_imza.png"
fig.savefig(p, dpi=140)
print("\n-> %s" % p)
json.dump(ATL, open(S174 / "ATLAS.json", "w"), indent=1, ensure_ascii=False)
print("-> %s" % (S174 / "ATLAS.json"))
