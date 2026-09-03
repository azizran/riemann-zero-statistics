"""
171k — FİGÜR + ν MERDİVENİ + c'NİN MİNİMUMU
============================================
Yeni ölçüm YOK; 171c/d/e/h'nin JSON'larından çizim ve iki türev hesabı:
  (i)  ν_λ = Δlog KALİB / Δlog W_X  merdiveni (yedi λ noktası)
  (ii) c(λ)'nın minimumunun yeri (üç düşük-λ noktasından parabol)

ÖN-MÜHÜR (koşudan ÖNCE, 4 Eylül 2026 00:25):
 W1 c(λ) minimumu λ* ≈ 0.649 (0.62–0.68), c_min ≈ 0.368 — ve λ* tam olarak
    ν_λ = 1 noktasıdır (özdeşlik: dc/dλ = 0 ⟺ ν_λ = 1).
 W2 ν merdiveni: 0.50←0.60 için ≈ 1.50; 1.30←1.15 için ≈ +0.35 (yani
    λ ≥ 1.15'te ν SIFIRDAN GERİ TIRMANIYOR — 170'in "ν tekdüze tırmanır"
    resmi λ ∈ [0.6, 1.15] penceresinin yasasıymış).

Çıktı: 171_nu.png + scratchpad/171/NU_MERDIVEN.json

SONUÇ (gerçek koşudan):
 W1 ✓✓ c(λ)'nın minimumu λ* = **0.6487**, c_min = **0.3675** (öngörü 0.649).
       Özdeşlik doğrulandı: ν'nün 1'i kestiği yer ile aynı.
 W2 ✓✓ ν merdiveni (5-bant ort.): 0.60→0.50 **+1.489 ± 0.063** (ν−1 = +7.7σ),
       0.70→0.60 +0.993, 0.85→0.70 +0.563, 1.00→0.85 +0.205,
       1.15→1.00 −0.034, **1.30→1.15 +0.271 ± 0.073** (öngörü ≈ +0.35).
       ⇒ ν(λ) TEKDÜZE DEĞİL: λ ≈ 1.08'de minimumu var; 170 §K0.3'ün
       "tekdüze tırmanış" resmi o minimumun sol kanadıymış.
 Doğrudan c farkı c(0.50) − c(0.60) = +0.0111 ± 0.0102 = +1.1σ; hükmü
 taşıyan istatistik ν'dür (+7.7σ), çünkü bantlar arası saçılım oranda
 sadeleşir.
"""
import importlib
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
for _p in ("169_configs", "167_configs", "166_configs", "165_configs",
           "163_configs", "160_configs", "159_configs", "155_configs",
           "154_configs"):
    sys.path.insert(0, str(QM / _p))
K1 = importlib.import_module("169_k1")

SCR = ("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
       "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S171 = SCR + "/171"
C_HIP = 4.0 / np.pi ** 2
SIRA = ["L130", "L115", "Hkeskin", "L085", "L070", "L060", "L050"]
LAMS = [1.30, 1.15, 1.00, 0.85, 0.70, 0.60, 0.50]


def main():
    D = K1.yukle()
    AO = json.load(open(S171 + "/A_ONKAYIT.json"))
    AY = json.load(open(S171 + "/A_YUZLESME.json"))
    MO = json.load(open(S171 + "/M_ONKAYIT.json"))
    YZ = json.load(open(S171 + "/T2C_YUZLESME.json"))

    # ---- bant defterleri ----------------------------------------------
    B = {}
    for g in SIRA + ["son"]:
        rows = [K1.band_jk(b["cizgi"]) for b in K1.saglikli(D[g], 0.68)]
        B[g] = dict(K=np.array([r["KALIB_u2"] for r in rows]),
                    sK=np.array([r["sKALIB_u2"] for r in rows]),
                    W=np.array([r["W_X"] for r in rows]),
                    t=np.array([r["tau_eff"] for r in rows]))
    defter = {d["gaz"]: d for d in YZ["defter"]}

    # ---- (i) ν merdiveni ------------------------------------------------
    print("=" * 88)
    print("171k — ν MERDİVENİ  (ν_λ = Δlog KALİB / Δlog W_X, orta bant ve "
          "5-bant ortalaması)")
    print("=" * 88)
    nu = []
    for a, b in zip(SIRA[:-1], SIRA[1:]):
        dk = np.log(B[b]["K"] / B[a]["K"])
        dw = np.log(B[b]["W"] / B[a]["W"])
        nu_b = float(dk[2] / dw[2])
        nu_o = float(dk.mean() / dw.mean())
        # jackknife hatası: KALİB'lerin bant-içi jk'leri (W_X marjinal, ~tam)
        sdk = np.sqrt((B[a]["sK"] / B[a]["K"]) ** 2
                      + (B[b]["sK"] / B[b]["K"]) ** 2)
        s_b = float(sdk[2] / abs(dw[2]))
        s_o = float(np.sqrt((sdk ** 2).sum()) / len(sdk) / abs(dw.mean()))
        nu.append(dict(dan=a, ye=b, lam_a=LAMS[SIRA.index(a)],
                       lam_b=LAMS[SIRA.index(b)], nu_orta=nu_b, nu_ort=nu_o,
                       s_orta=s_b, s_ort=s_o,
                       z1=float((nu_o - 1.0) / s_o)))
        print("  %-8s → %-8s (λ %.2f → %.2f):  ν = %+.3f ± %.3f (orta bant)"
              "   %+.3f ± %.3f (5 bant)   [ν−1 = %+.1fσ]"
              % (a, b, LAMS[SIRA.index(a)], LAMS[SIRA.index(b)],
                 nu_b, s_b, nu_o, s_o, (nu_o - 1.0) / s_o))
    print("  [170 §K0.3 aynı merdiveni λ ∈ [0.60, 1.15]'te −0.02 → 1.02 "
          "vermişti]")

    # ---- (ii) c'nin minimumu -------------------------------------------
    lam = np.array(LAMS)
    c = np.array([defter[g]["c"] for g in SIRA])
    s = np.array([defter[g]["s"] for g in SIRA])
    i = int(np.argmin(c))
    x = lam[i - 1:i + 2][::-1]
    y = c[i - 1:i + 2][::-1]
    p = np.polyfit(x, y, 2)
    lstar = float(-p[1] / (2 * p[0]))
    cmin = float(np.polyval(p, lstar))
    print("\n  c(λ) MİNİMUMU (üç nokta λ = %s):  λ* = %.4f   c_min = %.4f"
          % (list(x), lstar, cmin))
    print("  c defteri: " + "  ".join("%.2f:%.4f" % (l, v)
                                      for l, v in zip(lam, c)))
    print("  ⇒ 170 §K2'nin 'λ ≤ 0.70'te TABAN' okuması ÖRNEKLEM-DIŞI ÖLDÜ: "
          "c(0.50) = %.4f > c(0.60) = %.4f (%+.1fσ)"
          % (c[6], c[5], (c[6] - c[5]) / np.hypot(s[6], s[5])))

    # ---- FİGÜR ---------------------------------------------------------
    fig, ax = plt.subplots(2, 2, figsize=(13.5, 9.6))
    fig.suptitle("171 — ν üyesi: KALİB_u2(λ,τ) = A(τ)·M(λ) — A'nın kimliği "
                 "ve M'nin iki ucu", fontsize=13)

    # (a) A(τ)
    a0 = ax[0, 0]
    tA = np.array(AO["tau"])
    A5 = np.array(AO["A5"])
    a0.errorbar(tA, A5, yerr=A5 * np.array(AO["sA"]), fmt="ko", ms=6,
                capsize=3, label="A(τ) ölçülen (5 gaz, uyum penceresi)")
    a0.errorbar([AY["t6"], AY["t7"]], [AY["A6"], AY["A7"]],
                yerr=[AY["A6"] * AY["sA6"], AY["A7"] * AY["sA7"]],
                fmt="rs", ms=8, capsize=4,
                label="ÖRNEKLEM-DIŞI (9 bant)")
    tt = np.linspace(0.52, 0.80, 200)
    sX_hk = D["Hkeskin"]["sigX"]
    t0 = tA[2]
    A2 = np.exp(-2 * np.pi ** 2 * (tt ** 2 - t0 ** 2) * sX_hk ** 2)
    a0.plot(tt, A2, "b-", lw=2.2,
            label="A2Hk: $W_X(τ;σ_{\\tilde X}(λ{=}1))$  **0 parametre**")
    for nm, sty in (("A1  DW-tam Hk", "g--"), ("A3s  165-F ailesi (c₂/c₁)",
                                               "m-."),
                    ("A4s  doymuş tarak-sayım 1/(1+κe^{τL})", "c:")):
        r = [q for q in AO["adaylar"] if q["ad"] == nm][0]
        if r["param"] is None:
            if nm.startswith("A1"):
                sd = D["Hkeskin"]["sigds"]
                v = np.exp(-0.5 * np.pi ** 2 * (tt ** 2 - t0 ** 2)
                           * (sd ** 2 + 4 * sX_hk ** 2))
        elif nm.startswith("A3s"):
            rr = r["param"]
            f = (np.sin(2 * np.pi * tt) ** 2 + rr * np.sin(np.pi * tt) ** 4)
            f0 = (np.sin(2 * np.pi * t0) ** 2 + rr * np.sin(np.pi * t0) ** 4)
            v = f / f0
        else:
            k = r["param"]
            v = (1 / (1 + k * np.exp(tt * 12.0296))) / \
                (1 / (1 + k * np.exp(t0 * 12.0296)))
        a0.plot(tt, v, sty, lw=1.4, label=nm.split("[")[0].strip())
    a0.axvspan(0.52, 0.71, color="0.92", zorder=0)
    a0.text(0.615, 0.63, "uyum penceresi", ha="center", fontsize=8,
            color="0.35")
    a0.set_xlabel("τ (bant)")
    a0.set_ylabel("A(τ)  (orta bantta 1)")
    a0.set_title("(a) A(τ)'nun kimliği — kazanan SIFIR parametreli", fontsize=11)
    a0.legend(fontsize=7.5, loc="lower left")
    a0.set_ylim(0.55, 1.20)
    a0.grid(alpha=0.3)

    # (b) M(λ)
    a1 = ax[0, 1]
    Mv = np.array([defter[g]["M"] for g in SIRA])
    sMv = np.array([MO["sM"].get(g, YZ.get(g, {}).get("sM", 0.007))
                    for g in SIRA])
    a1.errorbar(lam, Mv, yerr=sMv, fmt="ko", ms=7, capsize=3,
                label="M(λ) ölçülen")
    a1.errorbar([1.00], [defter["son"]["M"]], yerr=[MO["sM"]["son"]],
                fmt="r*", ms=15, label="GERÇEK ζ (σ_X̃ ⇒ λ_eş = 0.93)")
    ll = np.linspace(0.45, 1.35, 200)
    LC = MO["LAM_C"]
    a1.plot(ll, np.where(ll >= LC, 1.0, (LC / ll) ** 0.1898), "b--", lw=1.8,
            label="M2 kırık kuvvet (λ_c = 1.011, türetimli)")
    p9 = np.polyfit(np.log(lam[1:6]), np.log(Mv[1:6]), 2)
    a1.plot(ll, np.exp(np.polyval(p9, np.log(ll))), "g-", lw=1.6,
            label="M9 log-λ kuadratik (EMPİRİK ÇIPA)")
    for g, xx in (("L050", 0.50), ("L130", 1.30)):
        for nm, col in (("M2 kırık kuvvet (λ_c=1.0109)", "b"),
                        ("M9 log-λ kuadratik [EMPİRİK]", "g")):
            a1.plot([xx], [YZ[g]["hukum"][nm]["M_on"]], col + "v", ms=9,
                    mfc="none", mew=1.6)
    a1.axvspan(0.45, 0.55, color="#ffe9e9", zorder=0)
    a1.axvspan(1.22, 1.35, color="#ffe9e9", zorder=0)
    a1.text(0.50, 1.26, "ÖN-KAYITLI", ha="center", fontsize=8, color="#a33")
    a1.text(1.285, 1.26, "ÖN-KAYITLI", ha="center", fontsize=8, color="#a33")
    a1.set_xlabel("λ")
    a1.set_ylabel("M(λ) = geo.ort_b KALİB_b/KALİB_b(Hk)")
    a1.set_title("(b) M(λ): türetimli aileler düşük uçta ÖLDÜ", fontsize=11)
    a1.legend(fontsize=7.5, loc="upper right")
    a1.grid(alpha=0.3)

    # (c) c(λ)
    a2 = ax[1, 0]
    a2.errorbar(lam, c, yerr=s, fmt="ko-", ms=7, capsize=3, lw=1.2,
                label="c = KALİB_u2/W_X")
    a2.errorbar([1.00], [defter["son"]["c"]], yerr=[defter["son"]["s"]],
                fmt="r*", ms=15, label="GERÇEK ζ")
    a2.axhline(C_HIP, color="0.4", ls=":", lw=1.5, label="4/π² = 0.40528")
    xx = np.linspace(0.48, 0.74, 100)
    a2.plot(xx, np.polyval(p, xx), "b--", lw=1.4,
            label="parabol ⇒ λ* = %.3f, c_min = %.4f" % (lstar, cmin))
    a2.axvline(lstar, color="b", ls=":", lw=1)
    for g, xx2, col in (("L050", 0.50, "#c33"), ("L130", 1.30, "#c33")):
        a2.plot([xx2], [YZ[g]["hukum"]["M2 kırık kuvvet (λ_c=1.0109)"]["c_on"]],
                "v", color=col, ms=9, mfc="none", mew=1.6)
    a2.annotate("ön-kayıtlı M2 öngörüsü", xy=(0.50, 0.3541),
                xytext=(0.58, 0.335), fontsize=8, color="#c33",
                arrowprops=dict(arrowstyle="->", color="#c33", lw=1))
    a2.set_xlabel("λ")
    a2.set_ylabel("c")
    a2.set_title("(c) c(λ)'nın TABANI YOK: minimum λ* ≈ %.2f" % lstar,
                 fontsize=11)
    a2.legend(fontsize=7.5, loc="upper left")
    a2.grid(alpha=0.3)

    # (d) şekil genişliği
    a3 = ax[1, 1]
    gen = MO["genislik"]
    grup = [("λ ailesi", ["L115", "Hkeskin", "L085", "L070", "L060"], "o",
             "tab:blue"),
            ("yeni uçlar", ["L050", "L130"], "D", "tab:red"),
            ("gerçek ζ", ["son"], "*", "k"),
            ("kesim", ["K090", "K070", "HA4", "E060"], "s", "tab:orange"),
            ("pencere", ["HkT2a", "HkT2b", "HkT4a", "HkT4b"], "^",
             "tab:green")]
    xk = 0
    for et, gs, mk, col in grup:
        for g in gs:
            v = (YZ[g]["sX_eff"] if g in ("L050", "L130")
                 else gen[g]["sX_eff"])
            a3.plot([xk], [v], mk, color=col, ms=11 if mk == "*" else 7)
            a3.text(xk, v + 0.004, g, rotation=90, fontsize=6.5,
                    ha="center", va="bottom")
            xk += 1
        xk += 0.7
    ge5 = [gen[g]["sX_eff"] for g in ("L115", "Hkeskin", "L085", "L070",
                                      "L060")]
    a3.axhspan(min(ge5), max(ge5), color="#dce8f7", zorder=0,
               label="λ-ailesi bandı")
    a3.axhline(sX_hk, color="b", ls="--", lw=1.5,
               label="σ_X̃(λ = 1.00) = %.5f" % sX_hk)
    a3.set_xticks([])
    a3.set_ylabel("gaz-başına şekil genişliği σ*/2")
    a3.set_title("(d) A(τ) λ-değişmezliğinin GEÇERLİLİK PENCERESİ", fontsize=11)
    a3.legend(fontsize=8, loc="upper left")
    a3.grid(alpha=0.3, axis="y")
    a3.set_ylim(0.17, 0.36)

    fig.tight_layout(rect=[0, 0, 1, 0.965])
    p_ = QM / "171_nu.png"
    fig.savefig(p_, dpi=140)
    print("\n-> %s" % p_)
    json.dump(dict(nu=nu, lam_star=lstar, c_min=cmin,
                   lam=[float(x) for x in lam], c=[float(x) for x in c],
                   s=[float(x) for x in s], M=[float(x) for x in Mv]),
              open(S171 + "/NU_MERDIVEN.json", "w"), indent=1)
    print("-> %s/NU_MERDIVEN.json" % S171)


if __name__ == "__main__":
    main()
