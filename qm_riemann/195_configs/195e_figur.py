# -*- coding: utf-8 -*-
"""
195e — FİGÜR: 195_uydu_karakteri.png
Üst sıra: her ada için yapı çarpanı gücü N_b·|Ĝ_b(Δω)|² (gürültü tabanı birimi,
blok-ağırlıklı ortalama; log ölçek); zeta düşük profili açık gri arka plan;
yasak uydu konumları turuncu kesikli, izinliler mavi noktalı çizgi. [195b]
Alt sıra: her ada için L-eşli bantta K̃_χ/K̃_ζ çubukları (± jk se), izinli
mavi / yasak turuncu; kalem öngörüsü √k/φ(k) siyah düz çizgi (bant ±%35 gri
şerit); ajan veri-öncesi notu (k/φ(k))/ΠJ₀ gri kesikli (HÜKME GİRMEZ). [195c/195d]
"""
import importlib.util
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("o195", HERE / "195o_ortak.py")
o = importlib.util.module_from_spec(spec)
spec.loader.exec_module(o)

MAVI, TURUNCU = "#2a78d6", "#eb6834"
INK, INK2, GRI = "#0b0b0b", "#52514e", "#b9b8b2"

if __name__ == "__main__":
    P = np.load(o.S195 / "A_profiller.npz")
    H = json.load(open(o.S195 / "HUKUM_195.json"))
    ONK = json.load(open(o.S195 / "ONKAYIT_195.json"))
    A = json.load(open(o.S195 / "A_guc.json"))
    notlar = ONK["ajan_veri_oncesi_notlari_HUKME_GIRMEZ"]["sayilar"]
    dw = P["dw"]
    uy_bar = ["+log2", "-log2", "+log3", "-log3", "+log5", "+log6", "+log10", "+log(3/2)"]
    et_bar = ["+2", "−2", "+3", "−3", "+5", "+6", "+10", "+3/2"]

    plt.rcParams.update({"font.size": 8.5, "axes.edgecolor": INK2, "axes.labelcolor": INK,
                         "xtick.color": INK2, "ytick.color": INK2, "axes.linewidth": 0.6})
    fig, ax = plt.subplots(2, 5, figsize=(17, 7.2), gridspec_kw={"height_ratios": [1, 1.05]})
    zref = P["zeta_dusuk_Pn"]
    YMAX = 1.25 * max(H["KAYIT_tam_tablo"][ad][u]["R_bant"] + H["KAYIT_tam_tablo"][ad][u]["se_R_bant"]
                      for ad in o.ADALAR for u in uy_bar)
    for i, ad in enumerate(o.ADALAR):
        k = ONK["K0b"]["adalar"][ad]["k"]
        a0 = ax[0, i]
        a0.plot(dw, zref, color=GRI, lw=0.6, label="ζ düşük (ΔL ızgarası)")
        a0.plot(dw, P[f"{ad}_Pn"], color=INK, lw=0.7, label=o.ADA_ETIKET[ad])
        for u in o.A_LISTE + o.A_KAYIT:
            x = o.log_n(u)
            if o.yasak_mi(u, k):
                a0.axvline(x, color=TURUNCU, lw=1.0, ls="--", alpha=0.9, zorder=0)
            else:
                a0.axvline(x, color=MAVI, lw=0.8, ls=":", alpha=0.9, zorder=0)
        a0.set_yscale("log")
        a0.set_ylim(0.05, 200)
        a0.set_xlim(-2.0, 2.6)
        par = "tek" if ONK["K0b"]["adalar"][ad]["parite"] == 1 else "çift"
        zmax_y = max(A["sonuc"][ad]["uydular"][u]["z"] for u in o.A_LISTE if o.yasak_mi(u, k))
        zmin_i = min(A["sonuc"][ad]["uydular"][u]["z"] for u in ["+log2", "-log2", "+log3", "-log3"]
                     if not o.yasak_mi(u, k))
        a0.set_title(f"{o.ADA_ETIKET[ad]}  (k={k}, {par})\nyasak z maks {zmax_y:+.1f} · "
                     f"izinli-güçlü z min {zmin_i:+.1f}", fontsize=8.5, color=INK)
        a0.set_xlabel("Δω = ω − L_b")
        if i == 0:
            a0.set_ylabel("N_b·|Ĝ_b|²  (gürültü tabanı birimi)")
            a0.legend(loc="upper left", fontsize=7, frameon=False)
        a0.grid(axis="y", color="#e6e5e0", lw=0.5)
        a0.spines[["top", "right"]].set_visible(False)

        a1 = ax[1, i]
        T = H["KAYIT_tam_tablo"][ad]
        R = np.array([T[u]["R_bant"] for u in uy_bar])
        se = np.array([T[u]["se_R_bant"] for u in uy_bar])
        renk = [TURUNCU if T[u]["yasak"] else MAVI for u in uy_bar]
        xs = np.arange(len(uy_bar))
        a1.bar(xs, R, width=0.72, color=renk, edgecolor="white", linewidth=1.0, zorder=2)
        a1.errorbar(xs, R, yerr=se, fmt="none", ecolor=INK, elinewidth=0.8, capsize=2, zorder=3)
        pred = notlar[ad]["kalem_sqrtk_phik"]
        alt2 = notlar[ad]["alt2_k_phik_bolu_PiJ0"]
        a1.axhspan(0.65 * pred, 1.35 * pred, color="#dcdbd5", alpha=0.6, zorder=0)
        a1.axhline(pred, color=INK, lw=1.2, zorder=1)
        a1.axhline(alt2, color=INK2, lw=1.0, ls="--", zorder=1)
        a1.axhline(0, color=INK2, lw=0.6)
        from matplotlib.lines import Line2D
        a1.legend(handles=[Line2D([], [], color=INK, lw=1.2, label=f"kalem √k/φ(k) = {pred:.3f} (±%35 şerit)"),
                           Line2D([], [], color=INK2, lw=1.0, ls="--",
                                  label=f"ajan notu (k/φ(k))/ΠJ₀ = {alt2:.2f} (hükme girmez)")],
                  loc="upper right", fontsize=6.5, frameon=True, facecolor="white", edgecolor="none",
                  framealpha=0.9)
        a1.set_xticks(xs)
        a1.set_xticklabels([f"log {e}" for e in et_bar], rotation=45, fontsize=7)
        a1.set_ylim(-0.7, YMAX)
        U = ONK["ada_bloklari"][ad]["bant_U"]
        a1.set_title(f"L-eşli bant [10.0, {U:.3f}):  K̃_χ / K̃_ζ,düşük", fontsize=8.5, color=INK)
        if i == 0:
            a1.set_ylabel("K̃_χ / K̃_ζ  (Re; ± jk se)")
        a1.spines[["top", "right"]].set_visible(False)
        a1.grid(axis="y", color="#e6e5e0", lw=0.5, zorder=0)
    from matplotlib.patches import Patch
    from matplotlib.legend import Legend
    ax[1, 0].add_artist(Legend(ax[1, 0], [Patch(color=MAVI), Patch(color=TURUNCU)],
                               ["izinli", "yasak (p | k)"], loc="upper left", fontsize=7,
                               frameon=False, bbox_to_anchor=(0.0, 0.80)))
    fig.suptitle("195 — Uyduların karakteri: L-fonksiyonu adalarında tarak iptali  "
                 f"(H-195A {H['H_195A']['hukum']} · B1 {H['H_195B1']['hukum']} · "
                 f"B2 {H['H_195B2']['hukum']} · B3 {H['H_195B3']['hukum']})", fontsize=10.5, color=INK)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(o.QM / "195_uydu_karakteri.png", dpi=150, facecolor="white")
    print("-> 195_uydu_karakteri.png")
