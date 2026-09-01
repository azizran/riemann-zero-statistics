"""
161 — FİGÜR: δ haritası (6 panel)
=================================
P1  δ(τ_eff) tayfları — dokuz gaz, taban 0.40   (haritanın kendisi)
P2  δ tayfı, GERÇEK gaz × beş taban              (konvansiyon ötelemesi)
P3  ANA PANEL: (δ(½), dδ/dτ|½) düzlemi — konvansiyon bulutları + elipsler
P4  üçüncü eksen: (δ(½), b) düzlemi
P5  tamamlayıcılık: hangi eksen hangi gazı yakalıyor (σ_konv birimi)
P6  merdiven: A4 → N5 → P1 → gerçek, dört nicelikte normalize
Çıktı: 161_delta_harita.png
"""
import importlib
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse

sys.path.insert(0, str(Path(__file__).resolve().parent))
K = importlib.import_module("161_cekirdek")
A = importlib.import_module("161_analiz")

REN = {"keskin": "#7b3294", "A4": "#c2a5cf", "J14": "#2166ac",
       "J26": "#67a9cf", "N5": "#1a9850", "N5z": "#a6d96a",
       "P1": "#e08214", "son": "#b2182b", "orta": "#d6604d"}
ETI = {"keskin": "keskin (saf merdiven)", "A4": "A4 (erfc-0.68)",
       "J14": "J14 (yapısız titreşim)", "J26": "J26 (titreşim ×2)",
       "N5": "N5 (kısa-menzil itme)", "N5z": "N5z (itme ×8)",
       "P1": "P1 (GUE aralık)", "son": "GERÇEK son (L=12.03)",
       "orta": "GERÇEK orta (L=11.46)"}
TABAN_IYI = A.TABAN_IYI


def elips(ax, X, renk, n=2.0, **kw):
    """Bulutun n-σ kovaryans elipsi."""
    if len(X) < 3:
        return
    m, C = X.mean(0), np.cov(X.T, ddof=1)
    w, V = np.linalg.eigh(C)
    w = np.maximum(w, 1e-18)
    ang = np.degrees(np.arctan2(V[1, -1], V[0, -1]))
    ax.add_patch(Ellipse(m, 2 * n * np.sqrt(w[-1]), 2 * n * np.sqrt(w[0]),
                         angle=ang, fc=renk, ec=renk, alpha=0.13, lw=1.2, **kw))


def main():
    D = K.yukle()
    P = {k: K.delta_bant(K.bant(d)) for k, d in D.items() if k[0] != "P0"}
    R = A.topla(P)

    fig, ax = plt.subplots(2, 3, figsize=(21, 12.5))
    fig.suptitle("161 — rotor düzeltmesi sonrası δ haritası:  "
                 r"$\delta(\tau)=\varphi_\Gamma-(4\pi\tau_{\rm eff}-2\pi)$"
                 "   ·   parmak-izi üçlüsü (δ(½), dδ/dτ|½, b)",
                 fontsize=15, weight="bold")

    # ---------------------------------------------------------------- P1
    a0 = ax[0, 0]
    for v in K.MERDIVEN:
        p = P[(v, 0.40)]
        a0.errorbar(p["x"], p["d"], yerr=p["s_d"], marker="o", ms=3.5, lw=1.3,
                    capsize=2, color=REN[v], label=ETI[v])
    a0.axhline(0, color="k", lw=0.7)
    a0.axvline(0.5, color="k", lw=0.7, ls=":")
    a0.set_xlabel(r"$\tau_{\rm eff}$"); a0.set_ylabel(r"$\delta$")
    a0.set_title("P1 · δ tayfları (taban 0.40)\n"
                 "gerçek TEK BAŞINA aşağı dalıyor; sentetikler yukarıda kalıyor",
                 fontsize=10.5)
    a0.legend(fontsize=7.2, ncol=2, loc="upper left", framealpha=0.92)
    a0.grid(alpha=0.25)

    # ---------------------------------------------------------------- P2
    a1 = ax[0, 1]
    cm = plt.cm.viridis(np.linspace(0.05, 0.85, len(K.TABAN)))
    for i, t in enumerate(K.TABAN):
        p = P[("son", t)]
        a1.plot(p["x"], p["d"], "o-", ms=3.5, lw=1.3, color=cm[i],
                label=f"gerçek, taban {t:.2f}")
        p = P[("A4", t)]
        m = p["tau"] <= 0.615
        a1.plot(p["x"][m], p["d"][m], "s--", ms=3, lw=1.0,
                color=str(0.72 - 0.13 * i), alpha=0.95)
    a1.axhline(0, color="k", lw=0.7); a1.axvline(0.5, color="k", lw=0.7, ls=":")
    a1.set_xlabel(r"$\tau_{\rm eff}$"); a1.set_ylabel(r"$\delta$")
    a1.set_title("P2 · taban konvansiyonu δ'yı ÖTELİYOR\n"
                 "düz = gerçek (öteleme 0.13), kesik = A4 (öteleme ~0.01)",
                 fontsize=10.5)
    a1.legend(fontsize=7.6, loc="lower left"); a1.grid(alpha=0.25)

    # ---------------------------------------------------------------- P3
    a2 = ax[0, 2]
    for v in K.MERDIVEN:
        pts = A.sec(R, v, TABAN_IYI, A.PENCERELER, A.AGIRLIKLAR)
        X = np.array([[f["d_half"], f["dd"]] for f in pts], float)
        if not len(X):
            continue
        a2.scatter(X[:, 0], X[:, 1], s=13, color=REN[v], alpha=0.55,
                   edgecolors="none")
        elips(a2, X, REN[v])
        a2.scatter([X[:, 0].mean()], [X[:, 1].mean()], s=130, marker="*",
                   color=REN[v], edgecolors="k", linewidths=0.6, zorder=5,
                   label=ETI[v])
    a2.axhline(0, color="k", lw=0.7); a2.axvline(0, color="k", lw=0.7)
    a2.set_xlabel(r"$\delta(\frac{1}{2})$")
    a2.set_ylabel(r"$d\delta/d\tau|_{1/2}$   $(=a_{1/2}-4\pi)$")
    a2.set_title("P3 · ANA PANEL: (δ(½), dδ/dτ) düzlemi\n"
                 "nokta = bir konvansiyon seçimi (4 taban × 4 pencere × 2 ağırlık);"
                 " elips = 2σ_konv", fontsize=10.5)
    a2.annotate("N5z: dδ/dτ ekseninde gerçekle ÖRTÜŞÜYOR (0.9σ),\n"
                "δ(½) ekseninde AYRIŞIYOR (2.7σ)  →  2B: 2.9σ",
                xy=(0.069, -1.57), xytext=(-0.055, -1.20), fontsize=8.2,
                color="#4d7a1e",
                arrowprops=dict(arrowstyle="->", color="#4d7a1e", lw=1.1))
    a2.legend(fontsize=7.2, loc="upper left", framealpha=0.92); a2.grid(alpha=0.25)

    # ---------------------------------------------------------------- P4
    a3 = ax[1, 0]
    for v in K.MERDIVEN:
        pts = A.sec(R, v, TABAN_IYI, A.PENCERELER, A.AGIRLIKLAR)
        X = np.array([[f["d_half"], f["b"]] for f in pts], float)
        if not len(X):
            continue
        a3.scatter(X[:, 0], X[:, 1], s=11, color=REN[v], alpha=0.45,
                   edgecolors="none")
        elips(a3, X, REN[v])
        a3.scatter([X[:, 0].mean()], [X[:, 1].mean()], s=130, marker="*",
                   color=REN[v], edgecolors="k", linewidths=0.6, zorder=5,
                   label=ETI[v])
    a3.axhline(0, color="k", lw=0.9); a3.axvline(0, color="k", lw=0.7)
    a3.set_xlabel(r"$\delta(\frac{1}{2})$"); a3.set_ylabel(r"$b$  (eğrilik)")
    a3.set_title("P4 · üçüncü eksen b — b'nin İŞARETİ P1'i ayırıyor\n"
                 "(b, rotor düzeltmesinden ETKİLENMEZ: b(δ) ≡ b(φ))", fontsize=10.5)
    a3.legend(fontsize=7.2, loc="upper left"); a3.grid(alpha=0.25)

    # ---------------------------------------------------------------- P5
    a4 = ax[1, 1]
    eks = [("d_half", r"$\delta(\frac{1}{2})$"), ("dd", r"$d\delta/d\tau|_{1/2}$"),
           ("b", r"$b$")]
    Xg = {alan: np.array([f[alan] for f in A.sec(R, ("son", "orta"), TABAN_IYI,
                                                 A.PENCERELER, A.AGIRLIKLAR)])
          for alan, _ in eks}
    G3 = np.array([[f["d_half"], f["dd"], f["b"]]
                   for f in A.sec(R, ("son", "orta"), TABAN_IYI,
                                  A.PENCERELER, A.AGIRLIKLAR)], float)
    mg, Cg = G3.mean(0), np.cov(G3.T, ddof=1)
    w = 0.2
    xs = np.arange(len(K.SENTETIK))
    for j, (alan, ad) in enumerate(eks):
        vals = []
        for v in K.SENTETIK:
            Xv = np.array([f[alan] for f in A.sec(R, v, TABAN_IYI,
                                                  A.PENCERELER, A.AGIRLIKLAR)])
            vals.append(abs(Xv.mean() - Xg[alan].mean())
                        / np.sqrt(Xv.var(ddof=1) + Xg[alan].var(ddof=1)))
        a4.bar(xs + (j - 1.5) * w, vals, w, label=ad,
               color=["#4575b4", "#d73027", "#fdae61"][j], alpha=0.85)
    v3 = []
    for v in K.SENTETIK:
        Xv = np.array([[f["d_half"], f["dd"], f["b"]]
                       for f in A.sec(R, v, TABAN_IYI, A.PENCERELER,
                                      A.AGIRLIKLAR)], float)
        dmu = Xv.mean(0) - mg
        v3.append(float(np.sqrt(dmu @ np.linalg.solve(np.cov(Xv.T, ddof=1) + Cg, dmu))))
    a4.bar(xs + 1.5 * w, v3, w, label="3B birleşim", color="#1a9850", alpha=0.95)
    a4.axhline(3, color="k", ls="--", lw=1.0)
    a4.text(-0.45, 3.15, "3σ_konv", fontsize=8)
    a4.axhline(1.1, color="gray", ls=":", lw=1.0)
    a4.text(3.3, 1.25, "158'in a-uzayındaki en zayıf ayrımı: N5z = 1.1σ",
            fontsize=7.8, color="dimgray")
    a4.set_xticks(xs); a4.set_xticklabels(K.SENTETIK, fontsize=9)
    a4.set_ylabel("gerçekten ayrım  (σ_konv)")
    a4.set_title("P5 · TAMAMLAYICILIK: üç eksen AYRI gazları yakalıyor\n"
                 "dδ/dτ N5z'yi kaçırır · δ(½) yakalar · b P1'i yakalar",
                 fontsize=10.5)
    a4.legend(fontsize=8.5); a4.grid(alpha=0.25, axis="y")

    # ---------------------------------------------------------------- P6
    a5 = ax[1, 2]
    zin = ["A4", "N5", "P1", "son"]
    t158 = {v: float(np.mean([K.REF158_WA[(v, t)][0] for t in TABAN_IYI]))
            for v in K.MERDIVEN}
    a158 = {v: float(np.mean([K.REF158_WA[(v, t)][1] for t in TABAN_IYI]))
            for v in K.MERDIVEN}
    oz = {v: {alan: A.ort(R, v, alan)[0] for alan in ("d_half", "dd", "b")}
          for v in K.MERDIVEN}
    seri = {r"$\delta(\frac{1}{2})$": [oz[v]["d_half"] for v in zin],
            r"$d\delta/d\tau|_{1/2}$": [oz[v]["dd"] for v in zin],
            r"$\tau_0$ (158)": [t158[v] for v in zin],
            r"$a$ (158)": [a158[v] for v in zin]}
    mk = ["o", "s", "^", "D"]
    for i, (ad, vv) in enumerate(seri.items()):
        vv = np.array(vv, float)
        nrm = (vv - vv[0]) / (vv[-1] - vv[0]) * 100
        a5.plot(range(4), nrm, mk[i] + "-", ms=8, lw=1.8, label=ad, alpha=0.85)
    a5.set_xticks(range(4))
    a5.set_xticklabels(["A4\n(başlangıç)", "N5\n(itme)", "P1\n(GUE)",
                        "gerçek\n(hedef)"], fontsize=9)
    a5.set_ylabel("A4 → gerçek yolunun yüzdesi")
    a5.axhline(0, color="k", lw=0.7); a5.axhline(100, color="k", lw=0.7)
    a5.set_title("P6 · MERDİVEN: δ(½) ve dδ/dτ, τ₀ ile a'nın\n"
                 "merdiven kesirlerini yeniden üretiyor", fontsize=10.5)
    a5.legend(fontsize=9); a5.grid(alpha=0.25)

    fig.tight_layout(rect=(0, 0, 1, 0.955))
    out = Path(__file__).resolve().parent.parent / "161_delta_harita.png"
    fig.savefig(out, dpi=125)
    print("[yazıldı]", out)


if __name__ == "__main__":
    main()
