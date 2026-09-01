"""160 — figür: δ'nın tam muhasebesi, iki çarpan, σ'nın kimliği, a/b ensemble."""
import importlib
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
C = importlib.import_module("160_cekirdek")
AN = importlib.import_module("160_analiz")
FOUR_PI = 4 * np.pi
OUT = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann/160_kuadratur.png")

REN = {"son": "#1a1a1a", "orta": "#666666", "keskin": "#1f77b4",
       "A4": "#d62728", "P1": "#2ca02c"}
ETI = AN.ETI


def B(v, t, g):
    d = AN.yukle(v, t, g)
    if not d:
        return [], []
    bb = AN.bantlar(d)
    return bb, AN.saglik(bb)


fig, ax = plt.subplots(2, 3, figsize=(19, 10.5))

# --- P1: muhasebe merdiveni (gerçek gaz) ------------------------------
a = ax[0, 0]
bb, sg = B("son", 0.4, "t1")
x = np.array([b["tau_eff"] for b in bb])
for ad, key, st in (("δ (ölçülen)", "delta", "o-"),
                    ("K1+K2 (TAM muhasebe)", "Kfull", "s--"),
                    ("K1 (reel ağırlık, tam kar.f.)", "K1", "^-"),
                    ("A·S (KALEM, 1. mertebe)", None, "v-")):
    y = np.array([b["AS"]["S_tam"] if key is None else b[key] for b in bb])
    a.plot(x, y, st, ms=6, lw=1.6, label=ad, alpha=.9 if key else .8)
a.axhline(0, color="k", lw=.6)
a.set_xlabel("τ_eff"); a.set_ylabel("faz [rad]")
a.set_title("P1  δ'nın muhasebesi — gerçek gaz (son, taban 0.40)\n"
            "K1+K2 δ'nın ÜSTÜNE oturuyor (bağıl ≤1.2e−05)", fontsize=10)
a.legend(fontsize=8); a.grid(alpha=.3)

# --- P2: bağıl hata, üç seviye, beş gaz -------------------------------
a = ax[0, 1]
for v in ("son", "keskin", "A4", "P1"):
    bb, sg = B(v, 0.4, "t1")
    s = np.array(sg)
    if not s.any():
        continue
    x = np.array([b["tau_eff"] for b in bb])[s]
    dl = np.array([b["delta"] for b in bb])[s]
    AS = np.array([b["AS"]["S_tam"] for b in bb])[s]
    K1 = np.array([b["K1"] for b in bb])[s]
    KF = np.array([b["Kfull"] for b in bb])[s]
    a.semilogy(x, np.abs((AS - dl) / dl) * 100, "v:", color=REN[v], ms=6,
               alpha=.75)
    a.semilogy(x, np.abs((K1 - dl) / dl) * 100, "^--", color=REN[v], ms=6,
               alpha=.9)
    a.semilogy(x, np.maximum(np.abs((KF - dl) / dl) * 100, 1e-6), "o-",
               color=REN[v], ms=5, label=ETI[v])
a.axhline(10, color="r", ls="--", lw=1.2)
a.text(0.445, 12, "±%10 hedefi", color="r", fontsize=8)
a.set_xlabel("τ_eff"); a.set_ylabel("|tahmin − δ| / |δ|  [%]")
a.set_title("P2  kapanış seviyesi (▽ A·S · △ K1 · ○ K1+K2)\n"
            "yalnız K1+K2 hedefin altında — ve 5 mertebe altında", fontsize=10)
a.legend(fontsize=8); a.grid(alpha=.3, which="both")

# --- P3: iki çarpan ---------------------------------------------------
a = ax[0, 2]
for v, t, ls in (("son", 0.4, "-"), ("son", 0.52, "--"),
                 ("keskin", 0.4, "-"), ("A4", 0.4, "-")):
    bb, sg = B(v, t, "t1")
    s = np.array(sg)
    if not s.any():
        continue
    x = np.array([b["tau_eff"] for b in bb])[s]
    K1 = np.array([b["K1"] for b in bb])[s]
    dl = np.array([b["delta"] for b in bb])[s]
    AS = np.array([b["AS"]["S_tam"] for b in bb])[s]
    m = (np.abs(dl) > 0.2 * np.abs(dl).max()) & (np.abs(AS) > 0.02)
    a.plot(x[m], (K1 / AS)[m], "^" + ls, color=REN[v], ms=6, alpha=.8,
           label=f"{ETI[v]} t{t} — kesme K1/(A·S)")
    a.plot(x[m], (dl / K1)[m], "o" + ls, color=REN[v], ms=6, mfc="none",
           label=f"{ETI[v]} t{t} — kuadratür δ/K1")
a.axhline(1, color="k", lw=.8)
a.set_ylim(0.3, 3.2)
a.set_xlabel("τ_eff"); a.set_ylabel("çarpan")
a.set_title("P3  iki çarpan: sonlu-kesme (△) × kuadratür (○)\n"
            "kuadratür çarpanı tabana duyarsız (%1), kesme çarpanı değil",
            fontsize=10)
a.legend(fontsize=7); a.grid(alpha=.3)

# --- P4: a ensemble (32 fit) ------------------------------------------
a = ax[1, 0]
acc = {k: [] for k in ("a", "δ", "K1", "A·S")}
for v in ("son", "orta"):
    for t in (0.28, 0.34, 0.4, 0.46):
        d = AN.yukle(v, t, "g158")
        if d is None:
            continue
        for pen in ("W-A", "W-D", "W-B", "W-C"):
            f = AN.fit_seti(d, pen)
            if f is None:
                continue
            acc["a"].append(f["a"])
            acc["δ"].append(FOUR_PI + f["d_δ"])
            acc["K1"].append(FOUR_PI + f["d_K1"])
            acc["A·S"].append(FOUR_PI + f["d_A·S"])
pos = [0, 1, 2, 3]
lab = ["a(φ)\nÖLÇÜLEN", "4π+dδ/dτ\n(K1+K2)", "4π+dK1/dτ\n(σ yok)",
       "4π+d(A·S)/dτ\n(1. mertebe)"]
for i, k in enumerate(("a", "δ", "K1", "A·S")):
    v = np.array(acc[k])
    a.errorbar(pos[i], v.mean(), yerr=v.std(ddof=1), fmt="o", ms=9, capsize=6,
               color=["#1a1a1a", "#2ca02c", "#ff7f0e", "#d62728"][i])
    a.scatter(np.full(len(v), pos[i]) + np.random.default_rng(3).normal(0, .05, len(v)),
              v, s=9, alpha=.35,
              color=["#1a1a1a", "#2ca02c", "#ff7f0e", "#d62728"][i])
a.axhspan(10.713 - 0.059, 10.713 + 0.059, color="#2ca02c", alpha=.14)
a.axhline(10.713, color="#2ca02c", lw=1.2, ls="--")
a.axhline(FOUR_PI, color="gray", lw=1, ls=":")
a.text(3.35, FOUR_PI, "4π (çıplak)", fontsize=8, color="gray", va="center")
a.text(1.05, 10.713 + 0.10, "158 referansı: 10.713 ± 0.059", fontsize=8,
       color="#2ca02c")
a.set_xticks(pos); a.set_xticklabels(lab, fontsize=8)
a.set_ylabel("a"); a.set_xlim(-0.5, 3.9)
a.set_title("P4  a — gerçek gazın 32 fiti (2 veri × 4 taban × 4 pencere)\n"
            "σ kanalı öngörüyü 6.2σ'dan 0.6σ'ya indiriyor", fontsize=10)
a.grid(alpha=.3, axis="y")

# --- P5: b kapanışı, gaz gaz ------------------------------------------
a = ax[1, 1]
gz = ["son", "orta", "keskin", "A4", "P1"]
w = 0.2
for j, (key, ad, c) in enumerate((("b", "b(φ) ölçülen", "#1a1a1a"),
                                  ("b_δ", "b(δ) = b(K1+K2)", "#2ca02c"),
                                  ("b_K1", "b(K1) — σ yok", "#ff7f0e"),
                                  ("b_A·S", "b(A·S)", "#d62728"))):
    vals = []
    for v in gz:
        d = AN.yukle(v, 0.4, "g158")
        f = AN.fit_seti(d, "W-A") if d else None
        vals.append(f[key] if f else np.nan)
    a.bar(np.arange(len(gz)) + (j - 1.5) * w, vals, w, color=c, label=ad,
          alpha=.9)
a.axhline(0, color="k", lw=.8)
a.annotate("A4: işareti\nSADECE σ kanalı\ndoğru veriyor", xy=(3.0, 3.2),
           xytext=(2.05, 9.6), fontsize=8, ha="center",
           arrowprops=dict(arrowstyle="->", lw=1.1))
a.set_ylim(-9.5, 12.5)
a.set_xticks(range(len(gz))); a.set_xticklabels([ETI[g] for g in gz], fontsize=8)
a.set_ylabel("b = eğrilik"); a.legend(fontsize=8)
a.set_title("P5  b — kuadratür kanalı beş gazın BEŞİNDE işareti veriyor\n"
            "(159'un tek karşı-örneği A4 çözüldü)", fontsize=10)
a.grid(alpha=.3, axis="y")

# --- P6: σ'nın kimliği ------------------------------------------------
a = ax[1, 2]
a2 = a.twinx()
for v in ("son", "keskin", "A4", "P1"):
    bb, sg = B(v, 0.4, "t1")
    s = np.array(sg)
    if not s.any() or "TANI" not in bb[0]:
        continue
    x = np.array([b["tau_eff"] for b in bb])[s]
    a.plot(x, np.array([b["TANI"]["c_sdr"] for b in bb])[s], "o-",
           color=REN[v], ms=5, label=ETI[v])
    a.plot(x, np.array([b["TANI"]["c_sr"] for b in bb])[s], "s:",
           color=REN[v], ms=4, alpha=.55)
    a2.plot(x, np.array([b["TANI"]["c_sx"] for b in bb])[s], "^--",
            color=REN[v], ms=4, alpha=.55)
a.axhline(0, color="k", lw=.6)
a.set_xlabel("τ_eff"); a.set_ylabel("korel(σ,Δρ) [○]  ·  korel(σ,ρ) [□]")
a2.set_ylabel("korel(σ, X̃)  [△, sağ eksen]", fontsize=9)
a.set_title("P6  σ'nın kimliği: ○ korel(σ,Δρ) = 0.02…0.38 · □ korel(σ,ρ) ≈ 0\n"
            "TAM kimlik (ρ+iσ)·e^{iA(n+1)} = η_{n+1}·e^{−iA·C_n}·K  "
            "(180 çizgide ≤1.2e−08)", fontsize=9)
a.legend(fontsize=7, loc="upper left"); a.grid(alpha=.3)

fig.suptitle("160 — δ'nın TAM muhasebesi: sonlu-kesme (K1) × kuadratür (K2). "
             "δ = arg⟨(ρ+iσ)e^{iAX̃}⟩ bağıl ≤7e−04 ile kapanıyor.",
             fontsize=13, y=0.995)
fig.tight_layout(rect=[0, 0, 1, 0.975])
fig.savefig(OUT, dpi=115)
print(f"-> {OUT}")
