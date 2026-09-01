"""159 — figür: kinematik omurga, mekanizma fazı, T1 merdiveni, T2, T3."""
import importlib
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
C = importlib.import_module("159_cekirdek")
AN = importlib.import_module("159_analiz")
FOUR_PI = 4 * np.pi
OUT = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann/159_phi_AS_sinavi.png")

REN = {"son": "#1a1a1a", "orta": "#666666", "keskin": "#1f77b4",
       "A4": "#d62728", "P1": "#2ca02c"}
ETI = AN.ETI


def B(v, t, g):
    d = AN.yukle(v, t, g)
    return (AN.bantlar(d), AN.saglik(AN.bantlar(d))) if d else ([], [])


fig, ax = plt.subplots(2, 3, figsize=(19, 10.5))

# --- P1: φ ham tayfı + kinematik omurga -------------------------------
a = ax[0, 0]
for v in ("son", "keskin", "A4", "P1"):
    bb, sg = B(v, 0.4, "t1")
    x = np.array([b["tau_eff"] for b in bb]); y = np.array([b["phi"] for b in bb])
    s = np.array(sg)
    a.plot(x[s], y[s], "o-", color=REN[v], ms=4, lw=1.3, label=ETI[v])
    if (~s).any():
        a.plot(x[~s], y[~s], "x", color=REN[v], ms=6, alpha=.45)
xx = np.linspace(0.44, 0.80, 50)
a.plot(xx, FOUR_PI * xx - 2 * np.pi, "k--", lw=2,
       label="kinematik omurga 4πτ−2π")
a.axhline(0, color="0.8", lw=0.6); a.axvline(0.5, color="0.8", lw=0.6)
a.set_xlabel("τ_eff"); a.set_ylabel("φ_Γ = arg Γ")
a.set_title("P1  ölçülen φ_Γ ve kinematik omurga\n(158'in 'ilkel gözlenebilir'i)")
a.legend(fontsize=8); a.grid(alpha=.25)

# --- P2: δ = φ − omurga ------------------------------------------------
a = ax[0, 1]
for v in ("son", "keskin", "A4", "P1"):
    bb, sg = B(v, 0.4, "t1")
    x = np.array([b["tau_eff"] for b in bb]); y = np.array([b["delta"] for b in bb])
    s = np.array(sg)
    a.plot(x[s], y[s], "o-", color=REN[v], ms=4, lw=1.3, label=ETI[v])
    if (~s).any():
        a.plot(x[~s], y[~s], "x", color=REN[v], ms=6, alpha=.5)
a.axhline(0, color="0.8", lw=0.6); a.axvline(0.5, color="0.8", lw=0.6)
a.set_xlabel("τ_eff"); a.set_ylabel("δ = arg[Σ cr·e^{−2iA}] (mekanizma fazı)")
a.set_title("P2  omurga çıkarıldıktan sonra kalan faz\n(× = bozuk bant)")
a.legend(fontsize=8); a.grid(alpha=.25)

# --- P3: mekanizma merdiveni (gerçek gaz) ------------------------------
a = ax[0, 2]
bb, sg = B("son", 0.4, "t1")
x = [b["tau_eff"] for b in bb]
a.plot(x, [b["delta"] for b in bb], "ko-", ms=5, lw=1.6, label="δ (ölçülen)")
a.plot(x, [b["M1"] for b in bb], "s-", color="#ff7f0e", ms=4, lw=1.3,
       label="M1 = arg⟨ρe^{iAX}⟩/⟨ρ⟩  (ρ, tüm mertebeler)")
a.plot(x, [b["AS"]["S_tam"] for b in bb], "^-", color="#9467bd", ms=4, lw=1.3,
       label="A·S = A·Cov(ρ,dsΔ)/⟨ρ⟩  (KALEM, 1. mertebe)")
a.plot(x, [b["AS"]["S_dem64"] for b in bb], "v-", color="#8c564b", ms=4, lw=1.1,
       label="A·S, P_loc=|ĉ|² (W=64 kayan pencere)")
a.axhline(0, color="0.8", lw=0.6)
a.set_xlabel("τ_eff"); a.set_ylabel("faz [rad]")
a.set_title("P3  mekanizma merdiveni — gerçek gaz (taban 0.40)")
a.legend(fontsize=7.5); a.grid(alpha=.25)

# --- P4: T1 oranları ---------------------------------------------------
a = ax[1, 0]
for v in ("son", "keskin", "A4", "P1"):
    bb, sg = B(v, 0.4, "t1")
    s = np.array(sg)
    x = np.array([b["tau_eff"] for b in bb])
    r2 = np.array([b["delta"] / b["AS"]["S_tam"] for b in bb])
    r3 = np.array([b["delta"] / b["M1"] if b["M1"] else np.nan for b in bb])
    a.plot(x[s], r2[s], "o-", color=REN[v], ms=4, lw=1.3, label=f"{ETI[v]}  δ/(A·S)")
    a.plot(x[s], r3[s], "s--", color=REN[v], ms=3, lw=1.0, alpha=.6,
           label=f"{ETI[v]}  δ/M1")
a.axhspan(0.75, 1.25, color="#2ca02c", alpha=.13)
a.axhline(1.0, color="#2ca02c", lw=1.2)
a.set_ylim(-1.5, 5.5)
a.set_xlabel("τ_eff"); a.set_ylabel("oran (ideal 1.00)")
a.set_title("P4  T1 HÜKMÜ — bant bant oran\nyeşil şerit = MÜHÜR eşiği ±%25")
a.legend(fontsize=6.5, ncol=2); a.grid(alpha=.25)

# --- P5: T2 kanal ayrışımı --------------------------------------------
a = ax[1, 1]
for v, mrk in (("son", "o"), ("A4", "s"), ("keskin", "^")):
    bb, sg = B(v, 0.4, "t1")
    s_ = np.array(sg)
    x = np.array([b["tau_eff"] for b in bb])
    yl = np.array([b["AS"]["S_lad"] for b in bb])
    ye = np.array([b["AS"]["S_eta"] for b in bb])
    a.plot(x[s_], yl[s_], mrk + "-", color=REN[v], ms=4,
           lw=1.4, label=f"{ETI[v]}  merdiven payı")
    a.plot(x[s_], ye[s_], mrk + "--", color=REN[v], ms=3,
           lw=1.0, alpha=.6, label=f"{ETI[v]}  η payı")
a.axhline(0, color="0.8", lw=0.6)
a.set_ylim(-0.8, 0.35)
a.set_xlabel("τ_eff"); a.set_ylabel("A·S kanal payı")
a.set_title("P5  T2 — S'nin TAM kanal ayrışımı\nA·S_tam = A·S_lad + A·S_η + A·S_drift")
a.legend(fontsize=6.5, ncol=1); a.grid(alpha=.25)

# --- P6: T3 a öngörüsü -------------------------------------------------
a = ax[1, 2]
o = {}
for v in ("son", "orta", "keskin", "A4", "P1"):
    d = AN.yukle(v, 0.4, "g158")
    if d is None:
        continue
    bb = AN.bantlar(d, 0.43, 0.61)
    x = [b["tau_eff"] for b in bb]
    t0p, ap, bp, _ = AN.fit2(x, [b["phi"] for b in bb],
                             [b["sPhi_jk"] for b in bb])

    def eg(y, s):
        c = np.polyfit(x, y, 2, w=1.0 / np.asarray(s))
        return float(c[1] + 2 * c[0] * t0p)
    sd = [b["sDelta_jk"] for b in bb]; sa = [b["sAS_jk"] for b in bb]
    o[v] = (ap, FOUR_PI + eg([b["delta"] for b in bb], sd),
            FOUR_PI + eg([b["M1"] for b in bb], sd),
            FOUR_PI + eg([b["AS"]["S_tam"] for b in bb], sa))
am = np.array([o[v][0] for v in o])
for i, (nm, mk, col) in enumerate(
        ((("4π + dδ/dτ  (ölçülen mekanizma)"), "o", "#1a1a1a"),
         ("4π + dM1/dτ  (ρ, tüm mertebeler)", "s", "#ff7f0e"),
         ("4π + d(A·S)/dτ  (KALEM 1. mertebe)", "^", "#9467bd"))):
    a.plot(am, [o[v][i + 1] for v in o], mk, color=col, ms=9, label=nm)
a.plot([10, 13.5], [10, 13.5], "k--", lw=1, label="1:1")
a.axhline(FOUR_PI, color="#2ca02c", lw=1.2, ls=":", label="çıplak 4π = 12.566")
for v, dy in zip(o, (-13, +7, -13, -13, -13)):
    a.annotate(ETI[v], (o[v][0], o[v][1]), fontsize=7,
               xytext=(5, dy), textcoords="offset points")
a.set_xlabel("a_meas = dφ_Γ/dτ|τ₀  (158'in sayısı)")
a.set_ylabel("a_pred")
a.set_title("P6  T3 — a öngörüsü\n(gerçek referans a = 10.713 ± 0.059)")
a.legend(fontsize=7); a.grid(alpha=.25)

fig.suptitle("159 — mekanizma denklemi φ_Γ(τ) = A·S(τ) sınavı:  "
             "φ_Γ = (4πτ_eff − 2π) + δ,  δ ≈ A·S", fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.965])
fig.savefig(OUT, dpi=125)
print("->", OUT)
