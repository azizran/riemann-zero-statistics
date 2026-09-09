# -*- coding: utf-8 -*-
"""
183j — FİGÜR: derinlik merdiveni + üç-nokta çatalı
===================================================
Yalnız 183/CATAL.json'u çizer; hiçbir yeni sayı üretmez.
Kullanım: 183j_figur.py
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
C = json.load(open(SCR / "183" / "CATAL.json"))
ON = json.load(open(SCR / "183" / "ONKAYIT_183.json"))
CAPA = ON["CAPA"]
M = C["merdiven"]
X = np.arange(len(M))
ETK = [m["d"] if m["d"] != "inf" else "∞" for m in M]

fig, ax = plt.subplots(2, 3, figsize=(16.5, 9.2))

# --- 1 korr(Xa,Xb) vs derinlik -------------------------------------
a = ax[0, 0]
y = [m["korr"] for m in M]
a.plot(X, y, "o-", color="C0", lw=2, ms=8)
a.axhline(CAPA["korr_XaXb_SUR"], color="C3", ls="--",
          label=f"SUR-A tabanı {CAPA['korr_XaXb_SUR']:+.4f}")
a.axhspan(CAPA["korr_XaXb_VF"] - CAPA["korr_XaXb_VF_sd"],
          CAPA["korr_XaXb_VF"] + CAPA["korr_XaXb_VF_sd"], color="C0",
          alpha=0.13, label="VF ailesi ±sd (n=8)")
a.axhline(0, color="k", lw=0.6)
a.set_xticks(X); a.set_xticklabels(ETK)
a.set_xlabel("çözücü derinliği d (Newton güncellemesi)")
a.set_ylabel("korr(Xa, Xb)")
a.set_title("(1) İnşa-kilidi imzası derinlikle değişmiyor")
a.legend(fontsize=8); a.grid(alpha=0.3)

# --- 2 Kov-oranı vs derinlik ---------------------------------------
a = ax[0, 1]
y = [m["kovoran"] for m in M]
a.plot(X, y, "s-", color="C2", lw=2, ms=8)
sur = C["taban"].get("kovoran", {}).get("ort")
if sur is not None:
    a.axhline(sur, color="C3", ls="--", label=f"SUR-B tabanı {sur:+.4f}")
a.axhline(CAPA["kovoran_Hkeskin"], color="C1", ls=":",
          label=f"ikiz (Hkeskin) {CAPA['kovoran_Hkeskin']:+.4f}")
a.axhline(0, color="k", lw=0.6)
a.set_xticks(X); a.set_xticklabels(ETK)
a.set_xlabel("çözücü derinliği d"); a.set_ylabel("Kov(η_artık, η)/Var(η)")
a.set_title("(2) İşaret turnusolü: VF pozitif, gerçek negatif")
a.legend(fontsize=8); a.grid(alpha=0.3)

# --- 3 ρ₃ vs derinlik ----------------------------------------------
a = ax[0, 2]
y = [m["rho3"] for m in M]
a.plot(X, y, "^-", color="C4", lw=2, ms=8)
a.axhline(CAPA["GAUSS3"], color="k", ls="--", label="GAUSS3 = (2/π)^{3/2}")
a.axhline(CAPA["rho3_sur_VF"], color="C3", ls="--",
          label=f"SUR-A tabanı {CAPA['rho3_sur_VF']:.4f}")
a.axhline(CAPA["rho3_Hkeskin"], color="C1", ls=":", label="ikiz 0.5284")
a.set_xticks(X); a.set_xticklabels(ETK)
a.set_xlabel("çözücü derinliği d"); a.set_ylabel("ρ₃ (3-bacak kırpma)")
a.set_title("(3) Gauss-altı derinliğe bağlı değil")
a.legend(fontsize=8); a.grid(alpha=0.3)

# --- 4 derinlik düğmesinin dinamik aralığı --------------------------
a = ax[1, 0]
r = [m["rms_kayma"] for m in M]
a.semilogy(X, [max(v, 1e-16) for v in r], "o-", color="C5", lw=2, ms=8)
a.set_xticks(X); a.set_xticklabels(ETK)
a.set_xlabel("çözücü derinliği d")
a.set_ylabel("rms |z(d) − z(∞)|")
a.set_title("(4) Düğmenin dinamik aralığı: braket zaten çözmüş")
a.grid(alpha=0.3, which="both")

# --- 5 üç-nokta çatalı ---------------------------------------------
a = ax[1, 1]
CT = C["catal"]["ikiz(Hkeskin)"]
ks = [k for k in ["kovoran", "kov_art", "kov_ciz_art", "m3_cizgi", "rho3",
                  "rho3_J"]
      if k in CT and CT[k].get("GK_ust") is not None]
for i, k in enumerate(ks):
    d = CT[k]
    n = max(abs(d["GK_alt"]), abs(d["GK_ust"])) or 1.0
    a.plot([d["GK_alt"] / n, d["GK_ust"] / n], [i, i], "-",
           color="C7", lw=6, alpha=0.5, solid_capstyle="butt")
    a.plot(d["GK_alt"] / n, i, "o", color="C0", ms=10,
           label="GK_alt = ikiz−VF" if i == 0 else None)
    a.plot(d["GK_ust"] / n, i, "D", color="C3", ms=9,
           label="GK_üst = ikiz−SUR" if i == 0 else None)
a.axvline(0, color="k", lw=0.8)
a.set_yticks(range(len(ks)))
a.set_yticklabels([f"{k}\n(IK payı {CT[k]['pay_insa']:+.3f})"
                   for k in ks], fontsize=8)
a.set_xlabel("çatal, en büyük kola göre normalize")
a.set_title("(5) ÇATAL: [alt sınır, üst kol]")
a.legend(fontsize=8, loc="lower right"); a.grid(alpha=0.3, axis="x")

# --- 6 K4: iki aile -------------------------------------------------
a = ax[1, 2]
K4 = C.get("K4", {})
if "VS" in K4 and K4["VS"].get("n_korr"):
    lab = ["korr(Xa,Xb)", "Kov-oranı"]
    vf = [K4["VF"]["korr"], K4["VF"]["kovoran"]]
    vfe = [K4["VF"]["korr_se"], K4["VF"]["kovoran_se"]]
    vs = [K4["VS"]["korr"], K4["VS"]["kovoran"]]
    vse = [K4["VS"]["korr_se"], K4["VS"]["kovoran_se"]]
    w = 0.34
    xx = np.arange(2)
    a.bar(xx - w / 2, vf, w, yerr=vfe, capsize=5, color="C0",
          label=f"VF (Hkeskin zarfı, n={K4['VF']['n_korr']})")
    a.bar(xx + w / 2, vs, w, yerr=vse, capsize=5, color="C1",
          label=f"VS (son zarfı, n={K4['VS']['n_korr']})")
    a.set_xticks(xx); a.set_xticklabels(lab)
    a.axhline(0, color="k", lw=0.6)
    a.set_title("(6) K4: inşa-kilidi iki ailede de var mı?")
    a.legend(fontsize=8); a.grid(alpha=0.3, axis="y")
else:
    a.text(0.5, 0.5, "K4: VS ölçümü yok", ha="center", va="center")
    a.axis("off")

fig.suptitle("183 — İNŞA-KİLİDİ: derinlik merdiveni ve üç-nokta çatalı "
             f"(ön-kayıt {ON['damga']})", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.96])
p = QM / "183_insa_kilidi.png"
fig.savefig(p, dpi=125)
print(f"-> {p}")
