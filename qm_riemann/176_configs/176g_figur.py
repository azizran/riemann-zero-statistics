# -*- coding: utf-8 -*-
"""
176g — FİGÜR (ölçümden SONRA yazıldı; hiçbir hüküm buradan çıkmaz)
==================================================================
Dört panel:
 (a) R_η ekseni: λ ailesi, KESİM ailesi, gerçek/ikizler ve VEKİLLER;
     ön-kayıtlı merkez (0.6949), bant [0.60,1.05] ve ölüm çizgisi 1.15.
 (b) Δlog(·←Hkeskin) çubukları: gerçek, erfc-ikiz ve VEKİL sütunu
     (= gazın TOPLAM kilit içeriği) — g_E, 2g_X, θ, M.
 (c) DC kaçağı: ⟨cosΔφ⟩ bant bant (F6/T-3) + Σξ toplamları.
 (d) H-F2 kaba kalem (F8): τ>0.70'te Δ_amp / Δ_res / Δ_faz ve
     kilit^HF2, 0.61 ± 0.10 hedef bandıyla.
"""
import json
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                        # noqa: E402

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S176, S174, S172 = SCR / "176", SCR / "174", SCR / "172"
VEK = sys.argv[1:] or ["VF1", "VF2"]

K2 = json.load(open(S176 / "K2.json"))
G1 = json.load(open(S172 / "G1.json"))


def k1(g, d=S174):
    return json.load(open(d / f"K1_{g}.json"))


def k3(g, d=S174):
    return json.load(open(d / f"K3_{g}.json"))


LAM7 = ["L050", "L060", "L070", "L085", "Hkeskin", "L115", "L130"]
KES = ["E060", "HA4", "K070", "K090", "Hkeskin"]

fig, ax = plt.subplots(2, 2, figsize=(15.5, 10.2))
fig.suptitle("176 — KİLİDİN FATURASI: faz-karıştırma sınavı "
             "(zarf birebir korunur, yalnız çizgi fazları karışır)",
             fontsize=13, fontweight="bold")

# ────────────────────────────────────────────────────────── (a) R_η
a = ax[0, 0]
Rl = [k1(g)["eta"]["R"] for g in LAM7]
Rk = [k1(g)["eta"]["R"] for g in KES]
Rs = k1("son")["eta"]["R"]
Rh = k1("Hkeskin")["eta"]["R"]
Rv = [k1(v, S176)["eta"]["R"] for v in VEK]
a.axhspan(0.60, 1.05, color="#cfe8ff", alpha=.6, zorder=0,
          label="ön-kayıtlı bant [0.60, 1.05]")
a.axhline(K2["F1"]["merkez"], color="#1f6fb4", ls="--", lw=1.6,
          label="ön-kayıtlı merkez %.4f" % K2["F1"]["merkez"])
a.axhline(1.15, color="#b40426", ls=":", lw=1.6, label="ölüm eşiği 1.15")
a.plot(np.full(len(Rl), 0.0) + np.linspace(-.16, .16, len(Rl)), Rl, "o",
       color="#888", ms=6, label="λ ailesi (7 gaz)")
a.plot(np.full(len(Rk), 1.0) + np.linspace(-.16, .16, len(Rk)), Rk, "s",
       color="#e08b00", ms=6, label="KESİM ailesi (5 gaz)")
a.plot([2.0], [Rh], "D", color="#2b7a2b", ms=10, label="Hkeskin (ikiz)")
a.plot([2.0], [Rs], "*", color="k", ms=18, label="son (gerçek)")
for i, (v, r) in enumerate(zip(VEK, Rv)):
    a.plot([3.0 + .18 * (i - .5)], [r], "v", color="#b40426", ms=13,
           label="VEKİL (faz karışık)" if i == 0 else None)
    a.annotate("%s %.4f" % (v, r), (3.0 + .18 * (i - .5), r),
               textcoords="offset points", xytext=(-16, 13 - 26 * i),
               fontsize=8, color="#b40426", ha="center")
a.annotate("", xy=(3.0, np.mean(Rv)), xytext=(2.0, Rh),
           arrowprops=dict(arrowstyle="->", color="#b40426", lw=2.0,
                           alpha=.8))
a.text(2.5, .5 * (Rh + np.mean(Rv)), "Λ_R = %+.4f" % K2["F1"]["Lambda_R"],
       ha="center", va="bottom", fontsize=10, color="#b40426",
       fontweight="bold")
a.set_xticks([0, 1, 2, 3])
a.set_xticklabels(["λ ailesi", "kesim ailesi", "ikiz / gerçek", "VEKİL"])
a.set_ylabel(r"$R_\eta = \Sigma|c_q(\eta)|^2/2\ \div\ \mathrm{Var}(\eta)$")
a.set_title("(a) girişim oranı: faz kilidi sökülünce", fontsize=11)
a.legend(fontsize=7.5, loc="center left", ncol=1)
a.grid(alpha=.25)

# ─────────────────────────────────────────── (b) Δlog(·←Hkeskin)
b = ax[0, 1]
T = K2["F9"]["tablo"]
ad = [t["ad"] for t in T]
x = np.arange(len(T))
w = .26
b.bar(x - w, [t["d_son"] for t in T], w, color="k", label="gerçek ← ikiz")
b.bar(x, [t["d_erfc"] for t in T], w, color="#e08b00",
      label="erfc-ikiz ← ikiz (KESİM ekseni)")
b.bar(x + w, [-t["d_vekil"] for t in T], w, color="#b40426",
      label="VEKİL ← ikiz (kilit sökülmüş)")
for i, t in enumerate(T):
    b.annotate("%+.4f" % t["d_son"], (i - w, t["d_son"]), ha="center",
               va="bottom" if t["d_son"] > 0 else "top", fontsize=7.5)
    b.annotate("%+.4f" % (-t["d_vekil"]), (i + w, -t["d_vekil"]),
               ha="center", va="top" if t["d_vekil"] > 0 else "bottom",
               fontsize=7.5, color="#b40426")
b.axhline(0, color="k", lw=.8)
b.set_xticks(x)
b.set_xticklabels(ad, fontsize=9)
b.set_ylabel(r"$\Delta\log$")
b.set_title("(b) çarpan defteri: gerçeğin fazlası, kesim ekseni, "
            "ve kilidin TOPLAM içeriği", fontsize=11)
b.legend(fontsize=8)
b.grid(alpha=.25, axis="y")

# ────────────────────────────────────────────── (c) DC kaçağı
c = ax[1, 0]
GAZ = [("son", "k", "*"), ("Hkeskin", "#2b7a2b", "D"),
       ("HA4", "#e08b00", "s")] + [(v, "#b40426", "v") for v in VEK]
for g, col, mk in GAZ:
    d = k3(g, S176 if g in VEK else S174)
    bb = [x for x in d["bant"] if x["lo"] >= 0.40 - 1e-9]
    xs = [.5 * (x["lo"] + x["hi"]) for x in bb]
    ys = [x["cos"] for x in bb]
    c.plot(xs, ys, mk + "-", color=col, ms=7, lw=1.4,
           label="%s   ort(E)=%+.4f" % (g, d["ortE"]))
c.axhline(0.90, color="#1f6fb4", ls="--", lw=1.4,
          label="F6 eşiği ⟨cosΔφ⟩ ≥ 0.90")
c.axhline(0, color="k", lw=.8)
c.axvspan(0.50, 0.80, color="#cfe8ff", alpha=.45, zorder=0)
c.set_xlabel(r"$\tau$")
c.set_ylabel(r"$\langle\cos\Delta\varphi\rangle$")
c.set_title("(c) DC kaçağı faz kilidine kör mü? (F6 / T-3)", fontsize=11)
c.legend(fontsize=7.5, loc="lower left")
c.grid(alpha=.25)

# ────────────────────────────────────────────── (d) H-F2 (F8)
d_ = ax[1, 1]
F8 = K2["F8"]
lab = [r"$\Delta_{\rm amp}$" + "\n(genlik)", r"$\Delta_{\rm res}$" + "\n(κ)",
       r"$\Delta_{\rm faz}$" + "\n(faz uyumu)"]
val = [F8["d_amp"], F8["d_res"], F8["d_faz"]]
col = ["#888", "#1f6fb4", "#b40426"]
pay = [v / F8["toplam"] for v in val]
d_.barh(np.arange(3), pay, .5, color=col)
d_.axvspan(F8["hedef"] - F8["tolerans"], F8["hedef"] + F8["tolerans"],
           color="#cfe8ff", alpha=.65, zorder=0,
           label="H-F2 hedef bandı 0.61 ± 0.10")
d_.axvline(F8["hedef"], color="#1f6fb4", ls="--", lw=1.6)
for i, (v, p) in enumerate(zip(val, pay)):
    d_.annotate("  %.1f%%   (Δlog %+.4f)" % (100 * p, v), (p, i),
                va="center", ha="left", fontsize=9.5)
d_.axvline(0, color="k", lw=.8)
d_.set_yticks(np.arange(3))
d_.set_yticklabels([r"$\Delta_{\rm amp}$  genlik", r"$\Delta_{\rm res}$  κ",
                    r"$\Delta_{\rm faz}$  faz uyumu"], fontsize=10)
d_.set_xlim(0, 0.92)
d_.invert_yaxis()
d_.set_xlabel(r"$\tau>0.70$ açığının payı  "
              r"($\Delta\log$ / toplam $-0.4091$)")
d_.set_title("(d) H-F2 kaba kalem: kilit$^{HF2}$ = %.4f  ⇒ %s\n"
             "(açığın %%%.0f'i GENLİK, %%%.0f'i faz uyumu — ön-kayıt "
             "tersini bekliyordu)"
             % (F8["kilit_HF2"], "YAŞADI" if F8["gecti"] else "ÖLDÜ",
                100 * pay[0], 100 * pay[2]), fontsize=10.5)
d_.legend(fontsize=8.5, loc="lower right")
d_.grid(alpha=.25, axis="x")

fig.tight_layout(rect=(0, 0, 1, 0.965))
p = QM / "176_kilit_faturasi.png"
fig.savefig(p, dpi=135)
print("-> %s" % p)
