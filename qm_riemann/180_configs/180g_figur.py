# -*- coding: utf-8 -*-
"""
180g — FİGÜR: NEDENSEL DEFTER (4 panel)
========================================
(a) H-180a: 3-bacak kırpma aktarımı ρ₃ — son, Hkeskin, VF1..4, Gauss;
    çöküş kesri C.
(b) ÇAPASIZ AYRIŞIM (doğrusal): her HAM nicelik için ΔN = ZARF + KİLİT.
(c) Aynı ayrışım LOG konvansiyonunda + eski çapa-bağıl %8.6/%91.4 çubuğu.
(d) Zarf sağlaması: r(τ) tarifi ↔ ölçülen R_bant(VS)/R_bant(VF).

Kullanım: 180g_figur.py
"""
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S180 = SCR / "180"
HAM = ("M", "Q_E", "rho_E", "Q_X", "rho_X")
ET = {"M": "M", "Q_E": "$Q_E$", "rho_E": r"$\rho_E$", "Q_X": "$Q_X$",
      "rho_X": r"$\rho_X$"}

OK = json.load(open(S180 / "ONKAYIT_K0.json"))
HK = json.load(open(S180 / "K2_HAKEM.json"))
DF = json.load(open(S180 / "K3_DEFTER.json"))

fig, ax = plt.subplots(2, 2, figsize=(15.5, 11))

# ---------------------------------------------------------------- (a)
a = ax[0, 0]
G3 = HK["GAUSS3"]
names = ["son", "Hkeskin"] + list(HK["rho3"].keys())
vals = [HK["rho3_son"], HK["rho3_Hk"]] + list(HK["rho3"].values())
cols = ["#c0392b", "#2c3e50"] + ["#2980b9"] * len(HK["rho3"])
a.bar(range(len(vals)), vals, color=cols)
a.axhline(G3, color="k", ls="--", lw=1.6,
          label=f"Gauss $(2/\\pi)^{{3/2}}$ = {G3:.4f}")
a.axhline(HK["rho3_bar"], color="#2980b9", ls=":", lw=1.4,
          label=f"$\\langle\\rho_3\\rangle_{{VF}}$ = {HK['rho3_bar']:.4f}")
a.set_xticks(range(len(vals)))
a.set_xticklabels(names, rotation=30)
a.set_ylim(min(min(vals), G3) - 0.006, max(vals) + 0.006)
a.set_ylabel(r"$\rho_3$  (3-bacak kırpma aktarımı)")
a.set_title(f"(a) H-180a — çöküş kesri C = {HK['C']:+.3f},  "
            f"$z_G$ = {HK['z_G']:+.2f}  ⇒ DAL {HK['dal']}")
a.legend(fontsize=8)
a.grid(alpha=0.3, axis="y")

# ---------------------------------------------------------------- (b,c)
for j, (key, ttl) in enumerate((("dogrusal", "(b) ÇAPASIZ AYRIŞIM — "
                                 "DOĞRUSAL konvansiyon"),
                                ("log", "(c) ÇAPASIZ AYRIŞIM — "
                                 "LOG konvansiyon"))):
    b = ax[0, 1] if j == 0 else ax[1, 0]
    D = DF[key]
    x = np.arange(len(HAM))
    pz = [D[k]["p_zarf"] for k in HAM]
    pk = [D[k]["p_kilit"] for k in HAM]
    se = [D[k]["se_p"] for k in HAM]
    b.bar(x - 0.2, pz, 0.38, yerr=se, capsize=4, color="#16a085",
          label="ZARF payı")
    b.bar(x + 0.2, pk, 0.38, yerr=se, capsize=4, color="#8e44ad",
          label="KİLİT payı")
    b.axhline(0, color="k", lw=0.9)
    b.axhline(1, color="k", lw=0.7, ls=":")
    b.set_xticks(x)
    b.set_xticklabels([ET[k] for k in HAM])
    b.set_ylabel("pay  (ΔN'ye oranla)")
    b.set_title(ttl + f"\nΔ{'log ' if key=='log' else ''}M = "
                f"{D['M']['delta']:+.6f}")
    for xi, (p, q) in enumerate(zip(pz, pk)):
        b.text(xi - 0.2, p, f"{100*p:+.0f}%", ha="center",
               va="bottom" if p >= 0 else "top", fontsize=8)
        b.text(xi + 0.2, q, f"{100*q:+.0f}%", ha="center",
               va="bottom" if q >= 0 else "top", fontsize=8)
    if key == "log":
        b.axhline(OK["eski_capa"]["manset_kilit_pay"], color="#c0392b",
                  ls="--", lw=1.3,
                  label="eski çapa-bağıl kilit payı %91.4")
    b.set_ylim(-2.6, 3.6)          # Q_X taşıyor: değerler etiketlerde
    b.legend(fontsize=8, loc="lower left")
    b.grid(alpha=0.3, axis="y")

# ---------------------------------------------------------------- (d)
d = ax[1, 1]
tg = np.array(OK["zarf"]["tau_eff"])
rr = np.array(OK["zarf"]["r"])
d.plot(tg, rr, "o-", color="#2c3e50", label=r"tarif  $r(\tau)$ "
       "= $R_{bant}$(son)/$R_{bant}$(Hkeskin)")
VS = [json.load(open(S180 / f"G_{g}.json")) for g in ("VS1", "VS2")]
VF = [json.load(open(SCR / "176" / f"G_{g}.json"))
      for g in ("VF1", "VF2", "VF3", "VF4")]
tv = np.array([b["tau_eff"] for b in VS[0]["bant"]])
rvs = np.mean([[b["R_bant"] for b in v["bant"]] for v in VS], axis=0)
rvf = np.mean([[r for r in v["R_bant_hukum"]] for v in VF], axis=0) \
    if "R_bant_hukum" in VF[0] else None
# VF'nin tam bant dizisi 167/C_VF*.json'dan
rvf_full = []
for g in ("VF1", "VF2", "VF3", "VF4"):
    bb = json.load(open(SCR / "167" / f"C_{g}.json"))["bant"]
    rvf_full.append([b["R_bant"] for b in bb if b.get("olculdu")])
rvf_full = np.mean(rvf_full, axis=0)
n = min(len(tv), len(rvf_full), len(rvs))
d.plot(tv[:n], rvs[:n] / rvf_full[:n], "s--", color="#e67e22",
       label=r"ölçülen  $\langle R_{bant}\rangle_{VS}/\langle R_{bant}"
             r"\rangle_{VF}$")
d.axhline(1.0, color="k", lw=0.8, ls=":")
d.set_xlabel(r"$\tau_{eff}$")
d.set_ylabel("oran")
d.set_title("(d) ZARF AKTARIMININ SAĞLAMASI (kapı değil, tanı)")
d.legend(fontsize=8)
d.grid(alpha=0.3)

fig.suptitle("180 — NEDENSEL DEFTER: çapasız zarf/kilit ayrışımı ve "
             "çifte-sayım hakemi", fontsize=13)
fig.tight_layout(rect=(0, 0, 1, 0.97))
p = QM / "180_nedensel_defter.png"
fig.savefig(p, dpi=125)
print(f"-> {p}")
