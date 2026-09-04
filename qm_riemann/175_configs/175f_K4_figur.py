# -*- coding: utf-8 -*-
"""
175f — K4 θ DEFTERİ (tam tablo) + FİGÜR
=========================================
YENİ ÖLÇÜM YOK: 175d'nin K2.json'u, 175e'nin K3M.json'u, 174/K1_*, K3_*
ve 172/G1-G2 okunur. Çıktı: 175/K4.json + 175_kesim_kilit.png
"""
import json
import math
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
S175, S174, S172 = SCR / "175", SCR / "174", SCR / "172"
G1 = json.load(open(S172 / "G1.json"))
G2 = json.load(open(S172 / "G2.json"))
K2 = json.load(open(S175 / "K2.json"))
ONK = json.load(open(S175 / "ONKAYIT_K2.json"))
K3M = json.load(open(S175 / "K3M.json")) if (S175 / "K3M.json").exists() else {}
K1 = {p.stem[3:]: json.load(open(p)) for p in sorted(S174.glob("K1_*.json"))}
K3 = {p.stem[3:]: json.load(open(p)) for p in sorted(S174.glob("K3_*.json"))}

KES5 = ["E060", "HA4", "K070", "K090", "Hkeskin"]
LAM7 = ["L050", "L060", "L070", "L085", "Hkeskin", "L115", "L130"]
LV = [.50, .60, .70, .85, 1.00, 1.15, 1.30]
TB = [ONK["kesim_ekseni"][g]["tau_bar"] for g in KES5]
BANT = [(0.40, 0.50), (0.50, 0.60), (0.60, 0.70), (0.70, 0.80), (0.80, 0.95)]

print("=" * 74)
print("175f — K4 θ DEFTERİ + FİGÜR")
print("=" * 74)

# ── K4: θ tam defteri ────────────────────────────────────────────────
print("\n  θ ve çarpanlar — λ ailesi ↔ KESİM ailesi ↔ gerçek")
print("  %-9s %7s %9s %9s %9s %9s %9s %9s"
      % ("gaz", "λ/τ̄_A", "θ(G1)", "θ(G2)", "g_E", "g_X", "g_cal", "M"))
T4 = {}
for g in LAM7:
    T4[g] = dict(aile="lam", x=LV[LAM7.index(g)], th1=G1[g]["th"],
                 th2=G2["theta"][g], gE=G1[g]["gE"], gX=G1[g]["gX"],
                 gcal=G1[g]["gcal"], M=G1[g]["KAL"])
    print("  %-9s %7.2f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f"
          % (g, T4[g]["x"], T4[g]["th1"], T4[g]["th2"], T4[g]["gE"],
             T4[g]["gX"], T4[g]["gcal"], T4[g]["M"]))
print("  " + "-" * 70)
for g, tb in zip(KES5, TB):
    T4[g] = dict(aile="kesim", x=tb, th1=G1[g]["th"], th2=G2["theta"][g],
                 gE=G1[g]["gE"], gX=G1[g]["gX"], gcal=G1[g]["gcal"],
                 M=G1[g]["KAL"])
    print("  %-9s %7.4f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f"
          % (g, tb, T4[g]["th1"], T4[g]["th2"], T4[g]["gE"], T4[g]["gX"],
             T4[g]["gcal"], T4[g]["M"]))
print("  " + "-" * 70)
T4["son"] = dict(aile="gerçek", x=float("nan"), th1=G1["son"]["th"],
                 th2=G2["theta"]["son"], gE=G1["son"]["gE"],
                 gX=G1["son"]["gX"], gcal=G1["son"]["gcal"], M=G1["son"]["KAL"])
print("  %-9s %7s %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f"
      % ("son", "—", T4["son"]["th1"], T4["son"]["th2"], T4["son"]["gE"],
         T4["son"]["gX"], T4["son"]["gcal"], T4["son"]["M"]))
for ad, k in (("θ(G1)", "th1"), ("θ(G2)", "th2"), ("g_cal", "gcal"),
              ("M", "M")):
    kv = [T4[g][k] for g in KES5]
    lv = [T4[g][k] for g in LAM7]
    print("  %-6s: kesim ailesi [%.5f, %.5f] %s ; λ ailesi [%.5f, %.5f] %s"
          % (ad, min(kv), max(kv),
             "İÇERİR" if min(kv) <= T4["son"][k] <= max(kv) else "**DIŞINDA**",
             min(lv), max(lv),
             "içerir" if min(lv) <= T4["son"][k] <= max(lv) else "**DIŞINDA**"))
json.dump(dict(T4=T4, f_theta=K2["P8"]["f_theta"]),
          open(S175 / "K4.json", "w"), indent=1, ensure_ascii=False)

# ── FİGÜR ────────────────────────────────────────────────────────────
fig, ax = plt.subplots(2, 2, figsize=(14.5, 10.2))
plt.rcParams.update({"font.size": 9})

# (a) R_η: λ ailesi (tümsek, gerçek ÜSTÜNDE) ↔ kesim ailesi (gerçek İÇİNDE)
a = ax[0, 0]
Rl = [K1[g]["eta"]["R"] for g in LAM7]
Rk = [K1[g]["eta"]["R"] for g in KES5]
Rs = K1["son"]["eta"]["R"]
a.plot(LV, Rl, "o-", color="#888", label="λ ailesi (174)")
a.axhline(Rs, color="crimson", lw=1.6, ls="--",
          label=f"gerçek gaz  R_η = {Rs:.4f}")
a.axhline(max(Rl), color="#888", lw=0.8, ls=":")
a.set_xlabel("λ (genlik ölçeği)")
a.set_ylabel(r"$R_\eta$  (girişim oranı, 162 makinesi)")
a2 = a.twiny()
a2.plot(TB, Rk, "s-", color="#1f77b4", label="KESİM ailesi (175)")
for x, y, g in zip(TB, Rk, KES5):
    a2.annotate(g, (x, y), textcoords="offset points", xytext=(3, -11),
                fontsize=7, color="#1f77b4")
a2.set_xlabel(r"$\bar\tau_A = \sum A^2\tau/\sum A^2$  (kesim ekseni)",
              color="#1f77b4")
a2.tick_params(axis="x", colors="#1f77b4")
h1, l1 = a.get_legend_handles_labels()
h2, l2 = a2.get_legend_handles_labels()
a.legend(h1 + h2, l1 + l2, fontsize=8, loc="lower center")
a.set_title("(a) aile-dışı fazlanın adresi: λ ailesi TAŞIMIYOR, "
            "KESİM ailesi TAŞIYOR", fontsize=9.5)

# (b) DC kaçağı bant defteri — SANDVİÇ
b = ax[0, 1]
xs = np.arange(len(BANT))
wd = 0.2
for i, (g, c) in enumerate((("son", "crimson"), ("Hkeskin", "#333"),
                            ("HA4", "#1f77b4"), ("K090", "#7fbf7f"))):
    if g not in K3:
        continue
    v = []
    for lo, hi in BANT:
        v.append(next((x["xi"] for x in K3[g]["bant"]
                       if abs(x["lo"] - lo) < 1e-9 and abs(x["hi"] - hi) < 1e-3),
                      np.nan))
    b.bar(xs + (i - 1.5) * wd, v, wd, color=c,
          label={"son": "gerçek", "Hkeskin": "keskin ikiz",
                 "HA4": "erfc ikiz (0.68/0.125)", "K090": "keskin τ≤0.90"}[g])
b.set_xticks(xs)
b.set_xticklabels(["%.2f–%.2f" % t for t in BANT])
b.set_xlabel(r"$\tau$ bandı")
b.set_ylabel(r"$\sum_q \xi_q$   (DC kaçağı)")
b.axhline(0, color="k", lw=0.6)
b.legend(fontsize=8)
b.set_title("(b) DC kaçağı: gerçek gaz iki ikizin ARASINDA (sandviç)",
            fontsize=9.5)

# (c) kesim kesri f
c = ax[1, 0]
FY = K2["P5"]["f_nicelik"]
sec = [k for k in ["R_η", "π_E", "ort(E)", "P_η", "Var(η)", "g_E", "g_X",
                   "m3", "skew(ds)", "μ̂²_E", "θ", "M = g_E g_X²θ", "σ_ds"]
       if k in FY]
vals = [FY[k] for k in sec]
LO, HI = -1.0, 1.6
cols = ["#1f77b4" if 0 <= v <= 1.2 else "crimson" for v in vals]
c.barh(range(len(sec)), [min(max(v, LO), HI) for v in vals], color=cols)
for i, v in enumerate(vals):
    if not (LO < v < HI):
        c.annotate(" %+.2f →" % v if v > 0 else "← %+.2f " % v,
                   (HI if v > 0 else LO, i), fontsize=7.5, va="center",
                   ha="right" if v > 0 else "left", color="white",
                   fontweight="bold")
c.set_yticks(range(len(sec)))
c.set_yticklabels(sec, fontsize=8)
c.axvline(0, color="k", lw=0.8)
c.axvline(1, color="k", lw=0.8, ls=":")
c.axvline(K2["P5"]["f_agirlikli"], color="green", lw=1.4,
          label="f (DC kaçağı, ağırlıklı) = %.3f" % K2["P5"]["f_agirlikli"])
c.set_xlabel(r"kesim kesri  $f=[y_{gerçek}-y_{keskin}]/[y_{erfc}-y_{keskin}]$"
             "   (ölçek kırpıldı)")
c.set_xlim(LO, HI)
c.legend(fontsize=8, loc="lower right")
c.set_title("(c) kesim payı ↔ kilit payı: f'ler AYNI DEĞİL "
            "(θ ve M'de TERS işaret)", fontsize=9.5)

# (d) κ oranı: 174e yasası ↔ W_pos'suz yasa ↔ çarpımsal düzeltme
d = ax[1, 1]
K3G = json.load(open(S175 / "K3G.json")) if (S175 / "K3G.json").exists() else {}
if "son" in K3M:
    bs = K3M["son"]["bant"]
    xc = [0.5 * (r["lo"] + r["hi"]) for r in bs]
    d.plot(xc, [r["med"] for r in bs], "o-", color="#888", ms=4,
           label=r"174e yasası ($\hat\kappa \propto W_{pos}$) — tüm çizgiler")
if "son" in K3G:
    bs = K3G["son"]["bant"]
    xc = [0.5 * (r["lo"] + r["hi"]) for r in bs]
    d.plot(xc, [r["a0"] for r in bs], "o-", color="crimson", ms=4,
           label=r"ASAL: $|\kappa|/|\hat\kappa_0|$  ($W_{pos}$ YOK)")
    d.plot(xc, [r["k0"] for r in bs], "^--", color="darkorange", ms=5,
           label=r"KULE: $|\kappa|/|\hat\kappa_0|$  (×2 fazla)")
    d.plot(xc, [r["k2"] for r in bs], "*-", color="#1f77b4", ms=10,
           label=r"KULE: $|\kappa|/|\hat\kappa_2|$  (+ çarpımsal $q_1q_2{=}Q$)")
d.axhline(1, color="k", lw=0.8)
d.axhspan(0.8, 1.25, color="green", alpha=0.12)
d.axvline(0.5, color="k", lw=0.6, ls=":")
d.set_yscale("log")
d.set_ylim(0.3, 12)
d.set_xlabel(r"$\tau$")
d.set_ylabel(r"medyan  $|\kappa| / |\hat\kappa|$   (log)")
d.legend(fontsize=7, ncol=1, loc="upper left")
d.set_title("(d) κ: fazla W_pos'un aşırı sönümüydü; KULE fazlası "
            "ÇARPIMSAL kanaldır", fontsize=9.5)

fig.suptitle("175 — FAZLANIN AYRIŞIMI: kesim payı + kilit payı  "
             "(ön-kayıt %s, sha %s)" % (ONK["zaman"][:16], ONK["sha256"][:12]),
             fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.97])
p = QM / "175_kesim_kilit.png"
fig.savefig(p, dpi=130)
print("\n-> %s" % p)
print("-> %s" % (S175 / "K4.json"))
