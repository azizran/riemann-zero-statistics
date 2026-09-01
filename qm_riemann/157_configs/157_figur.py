"""157 — 6 panelli figür (157_Rlad_yarisi.png)."""
import importlib
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
A = importlib.import_module("157_analiz")

SCR = A.SCR
HERE = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
REN = {0.28: "#7b3294", 0.34: "#c2a5cf", 0.40: "#2c7fb8",
       0.46: "#1a9641", 0.52: "#d95f02", 0.58: "#a6611a"}


def main():
    D = A.yukle()
    kayit = {}
    for v in A.VERI:
        for t in A.TABAN:
            if (v, t, "ince") in D:
                bs = A.bant(D[(v, t, "ince")])
            elif (v, t, "tau0") in D:
                bs = A.bant(D[(v, t, "tau0")])
            else:
                continue
            ust = 0.57 if len([b for b in bs if b["tau"] <= 0.57]) >= 5 \
                else 0.63
            kayit[(v, t)] = dict(
                seta=D[(v, t, "ince")]["s_eta"] if (v, t, "ince") in D
                else D[(v, t, "tau0")]["s_eta"],
                phi=A.tau0_bir(bs, "phi", ust)[0],
                tam=A.tau0_bir(bs, "tam", ust)[0],
                lad=A.tau0_bir(bs, "lad", ust)[0])

    fig, ax = plt.subplots(2, 3, figsize=(17.5, 9.6))

    # --- P1: tayf, taban 0.46, iki pencere
    a = ax[0, 0]
    for v, mk in (("son", "o"), ("orta", "s")):
        if (v, 0.46, "kaba") not in D:
            continue
        bs = A.bant(D[(v, 0.46, "kaba")])
        x = [b["tau_eff"] for b in bs]
        for kan, c in (("tam", "#666666"), ("lad", "#d7191c")):
            y = [b["kanal"][kan]["Rham"] for b in bs]
            e = [b["kanal"][kan]["sRham_jk"] for b in bs]
            a.errorbar(x, y, yerr=e, marker=mk, ms=4, lw=1.1, color=c,
                       ls="-" if v == "son" else "--",
                       label=f"R_{kan} ({v})")
    a.axhline(0, color="k", lw=0.6); a.axvline(0.5, color="b", lw=0.6, ls=":")
    a.set_xlabel("τ_eff"); a.set_ylabel("R (ham)")
    a.set_title("P1 — R_tam ve R_lad tayfı (taban 0.46, iki pencere)")
    a.legend(fontsize=7); a.grid(alpha=.3)

    # --- P2: taban yelpazesi (R_lad)
    a = ax[0, 1]
    for t in A.TABAN:
        if ("son", t, "kaba") not in D:
            continue
        bs = A.bant(D[("son", t, "kaba")])
        a.plot([b["tau_eff"] for b in bs],
               [b["kanal"]["lad"]["Rham"] for b in bs],
               "o-", ms=3.5, lw=1.1, color=REN[t], label=f"taban {t:.2f}")
    a.axhline(0, color="k", lw=0.6)
    a.set_xlabel("τ_eff"); a.set_ylabel("R_lad (ham)")
    a.set_title("P2 — R_lad TABAN YELPAZESİ (son)\n"
                "156 ±%3'ü yalnız 0.46–0.58'de ölçmüştü")
    a.legend(fontsize=7); a.grid(alpha=.3)

    # --- P3: R(taban) sabit τ'da
    a = ax[0, 2]
    for tau, ls in ((0.655, "-"), (0.715, "--")):
        for kan, c in (("tam", "#666666"), ("lad", "#d7191c")):
            xs, ys = [], []
            for t in A.TABAN:
                if ("son", t, "kaba") not in D:
                    continue
                bb = [b for b in A.bant(D[("son", t, "kaba")])
                      if b["tau"] == tau]
                if bb and bb[0]["kanal"][kan]["artik"] < 0.02:
                    xs.append(t); ys.append(bb[0]["kanal"][kan]["Rham"])
            if xs:
                a.plot(xs, ys, "o" + ls, color=c, ms=5,
                       label=f"R_{kan}, τ̄={tau}")
    a.axvspan(0.46, 0.58, color="#ffe08a", alpha=.45, zorder=0)
    a.text(0.52, a.get_ylim()[0], "156'nın penceresi", ha="center",
           fontsize=7, va="bottom")
    a.set_xlabel("regresyon tabanı"); a.set_ylabel("R (ham)")
    a.set_title("P3 — R'nin TABAN ekseni: R_lad'ın düzlüğü\n"
                "durağan bir noktanın yerel düzlüğü")
    a.legend(fontsize=7); a.grid(alpha=.3)

    # --- P4: τ₀ vs σ_η²  (üç nicelik üst üste = KİMLİK)
    a = ax[1, 0]
    for ad, c, mk in (("phi", "#2c7fb8", "o"), ("tam", "#666666", "s"),
                      ("lad", "#d7191c", "^")):
        for v, fs in (("son", "full"), ("orta", "none")):
            xs = [kayit[(v, t)]["seta"] for t in A.TABAN if (v, t) in kayit
                  and np.isfinite(kayit[(v, t)][ad])]
            ys = [kayit[(v, t)][ad] for t in A.TABAN if (v, t) in kayit
                  and np.isfinite(kayit[(v, t)][ad])]
            a.plot(xs, ys, mk, color=c, ms=7, fillstyle=fs,
                   label=f"τ₀({ad}) {v}")
    xg = np.linspace(0, 0.11, 20)
    a.plot(xg, 0.519 - 0.19 * xg, "k:", lw=1.2,
           label="155: 0.519 − 0.19σ_η²")
    a.axhline(0.5, color="b", lw=0.7, ls="--")
    a.set_xlabel("σ_η²  (regresyon tabanının vekili)")
    a.set_ylabel("τ₀")
    a.set_title("P4 — τ₀'ın TABAN kayması: üç tanım da aynı yasayı izliyor\n"
                "(kanal değiştirmek 0.19·σ_η² kaymasını kaldırmıyor)")
    a.legend(fontsize=6.5, ncol=2); a.grid(alpha=.3)

    # --- P5: yarış
    a = ax[1, 1]
    P = A.veri_kur(D, "lad", "kaba", 0.46)
    sg = A.sigma(P, "II")
    t0m = float(np.nanmean([kayit[(v, 0.46)]["lad"] for v in A.VERI
                            if (v, 0.46) in kayit]))
    AD = A.adaylar(P, sg, t0m)
    res = []
    for ad, (tanim, k, chi2, th0, mfn) in AD.items():
        c2, th = A._opt(chi2, th0)
        res.append((c2 + 2 * k, ad, c2, k, th, mfn))
    res.sort()
    x = np.array([p["x"] for p in P]); y = np.array([p["R"] for p in P])
    a.errorbar(x, y, yerr=sg, fmt="ko", ms=4, lw=1, label="ölçüm (R_lad)")
    for (aic, ad, c2, k, th, mfn), c in zip(res[:3],
                                            ("#d7191c", "#2c7fb8", "#1a9641")):
        o = np.argsort(x)
        a.plot(x[o], np.array(mfn(th))[o], "-", color=c, lw=1.4,
               label=f"{ad}: AIC={aic:.1f} (k={k})")
    a.axhline(0, color="k", lw=0.6)
    a.set_xlabel("τ_eff"); a.set_ylabel("R_lad (ham)")
    a.set_title("P5 — kapalı-form yarışı (taban 0.46, iki pencere ortak)\nAIC: hata modeli II = σ_jk ⊕ σ_pencere")
    a.legend(fontsize=7); a.grid(alpha=.3)

    # --- P6: φ_Γ'nın taban yelpazesi (ilkel nesnenin konvansiyonu)
    a = ax[1, 2]
    for t in A.TABAN:
        k = ("son", t, "ince") if ("son", t, "ince") in D else \
            ("son", t, "tau0")
        if k not in D:
            continue
        bs = [b for b in A.bant(D[k]) if b["tau"] <= 0.62]
        if len(bs) < 3:
            continue
        a.errorbar([b["tau_eff"] for b in bs], [b["phi"] for b in bs],
                   yerr=[b["sPhi_jk"] for b in bs], marker="o", ms=3.5,
                   lw=1.1, color=REN[t], label=f"taban {t:.2f}")
    a.axhline(0, color="k", lw=0.6); a.axvline(0.5, color="b", lw=0.7, ls=":")
    a.set_xlabel("τ_eff"); a.set_ylabel("φ_Γ")
    a.set_title("P6 — İLKEL nesne φ_Γ ve sıfırının taban kayması")
    a.legend(fontsize=7); a.grid(alpha=.3)

    fig.suptitle("157 — R_lad(τ): taban-bağımsız nesne sınavı ve "
                 "kapalı-form yarışı", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    p = HERE / "157_Rlad_yarisi.png"
    fig.savefig(p, dpi=135)
    print("->", p)


if __name__ == "__main__":
    main()
