"""
158 — FİGÜR: a'nın kimlik sınavı (4 panel)
==========================================
P1  φ_Γ(τ_eff) tayfı, taban 0.40 — bütün gazlar, negatif dal dahil
P2  ÇÖKÜŞ SINAVI: φ_Γ(τ_eff − τ₀) — eğriler üst üste biniyor mu?
    (biniyorsa a evrensel; binmiyorsa şekil gaza bağlı)
P3  a vs τ₀ — korelasyon merdiveni; hata çubukları taban yayılımı
P4  a vs taban — taban-bağımsızlık her gazda ayrı ayrı
Çıktı: 158_a_kimlik.png
"""
import importlib
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
A = importlib.import_module("158_analiz")

REN = {"keskin": "#7b3294", "A4": "#c2a5cf", "J14": "#2166ac",
       "J26": "#67a9cf", "N5": "#1a9850", "N5z": "#a6d96a",
       "P1": "#e08214", "son": "#b2182b", "orta": "#d6604d"}
ETI = {"keskin": "keskin (saf merdiven)", "A4": "A4 (erfc-0.68)",
       "J14": "J14 (yapısız titreşim)", "J26": "J26 (titreşim ×2)",
       "N5": "N5 (kısa-menzil itme)", "N5z": "N5z (itme ×8)",
       "P1": "P1 (GUE aralık)", "son": "GERÇEK son (L=12.03)",
       "orta": "GERÇEK orta (L=11.46)"}


def main():
    D = A.yukle()
    RA = {}
    for v in A.MERDIVEN:
        for t in A.TABAN:
            if (v, t) in D:
                f = A.olc(A.bant(D[(v, t)]), "A", 2)
                if f:
                    RA[(v, t)] = f

    fig, ax = plt.subplots(2, 2, figsize=(15.5, 11))
    fig.suptitle("158 — konvansiyonsuz sayı a = dφ_Γ/dτ|$_{τ_0}$: "
                 "makine-evrenseli mi, korelasyon-duyarlı mı?",
                 fontsize=14, weight="bold")

    # --- P1: ham tayf
    a0 = ax[0, 0]
    for v in A.MERDIVEN:
        if (v, 0.40) not in D:
            continue
        bs = A.bant(D[(v, 0.40)])
        x = [b["tau_eff"] for b in bs]
        y = [b["phi"] for b in bs]
        e = [b["sPhi_jk"] for b in bs]
        a0.errorbar(x, y, yerr=e, marker="o", ms=3.5, lw=1.3, capsize=2,
                    color=REN[v], label=ETI[v])
    a0.axhline(0, color="k", lw=0.7)
    a0.axvline(0.5, color="k", lw=0.7, ls=":")
    a0.set_xlabel("τ_eff  (güç-ağırlıklı bant frekansı)")
    a0.set_ylabel("φ_Γ = arg Γ_rot   [rad]")
    a0.set_title("P1 — ilkel gözlenebilir, taban 0.40 (negatif dal dahil)")
    a0.legend(fontsize=7.5, loc="upper left")
    a0.grid(alpha=0.25)

    # --- P2: çöküş sınavı
    a1 = ax[0, 1]
    for v in A.MERDIVEN:
        if (v, 0.40) not in RA or (v, 0.40) not in D:
            continue
        t0 = RA[(v, 0.40)]["t0"]
        bs = [b for b in A.bant(D[(v, 0.40)])
              if 0.43 - 1e-9 <= b["tau"] <= 0.61 + 1e-9]
        x = np.array([b["tau_eff"] for b in bs]) - t0
        y = np.array([b["phi"] for b in bs])
        a1.plot(x, y, marker="o", ms=3.5, lw=1.3, color=REN[v],
                label=f"{v}: a={RA[(v,0.40)]['a']:.2f}")
    u = np.linspace(-0.09, 0.11, 50)
    a1.plot(u, A.REF_A * u + A.REF_B * u**2, "k--", lw=2,
            label=f"157 gerçek: a={A.REF_A:.2f}, b={A.REF_B:.2f}")
    a1.axhline(0, color="k", lw=0.7); a1.axvline(0, color="k", lw=0.7)
    a1.set_xlabel("τ_eff − τ₀(gazın kendi sıfırı)")
    a1.set_ylabel("φ_Γ   [rad]")
    a1.set_title("P2 — ÇÖKÜŞ SINAVI: sıfıra kaydırılınca eğriler bindi mi?")
    a1.legend(fontsize=7.5, loc="upper left")
    a1.grid(alpha=0.25)

    # --- P3: a vs τ₀
    a2 = ax[1, 0]
    for v in A.MERDIVEN:
        ks = [(v, t) for t in A.TABAN if (v, t) in RA
              and RA[(v, t)]["n_neg"] >= 1]
        if not ks:
            continue
        t0 = np.array([RA[k]["t0"] for k in ks])
        aa = np.array([RA[k]["a"] for k in ks])
        a2.errorbar(t0.mean(), aa.mean(),
                    xerr=(t0.max() - t0.min()) / 2 or None,
                    yerr=(aa.max() - aa.min()) / 2 or None,
                    marker="s", ms=9, capsize=4, color=REN[v], label=ETI[v])
    a2.axhspan(A.REF_A - A.REF_EA, A.REF_A + A.REF_EA, color="#b2182b",
               alpha=0.15)
    a2.axhline(A.REF_A, color="#b2182b", ls="--", lw=1.4)
    a2.text(0.4885, A.REF_A + 0.06, "157: a = 10.759 ± 0.114 (gerçek)",
            color="#b2182b", fontsize=8)
    a2.axvline(0.5, color="k", lw=0.7, ls=":")
    a2.set_xlabel("τ₀ (gazın kendi sıfır-geçişi)")
    a2.set_ylabel("a = dφ_Γ/dτ|$_{τ_0}$")
    a2.set_title("P3 — a, τ₀'ın korelasyon merdivenini izliyor mu?")
    a2.legend(fontsize=7.5)
    a2.grid(alpha=0.25)

    # --- P4: a vs taban
    a3 = ax[1, 1]
    for v in A.MERDIVEN:
        ks = [t for t in A.TABAN if (v, t) in RA]
        if not ks:
            continue
        ic = [t for t in ks if RA[(v, t)]["n_neg"] >= 1]
        dis = [t for t in ks if RA[(v, t)]["n_neg"] < 1]
        a3.errorbar(ic, [RA[(v, t)]["a"] for t in ic],
                    yerr=[RA[(v, t)]["e_a"] for t in ic],
                    marker="o", ms=5, lw=1.4, capsize=3, color=REN[v],
                    label=ETI[v])
        if dis:
            a3.plot(dis, [RA[(v, t)]["a"] for t in dis], "x", ms=7,
                    color=REN[v])
    ic_hepsi = [RA[k]["a"] for k in RA if RA[k]["n_neg"] >= 1]
    a3.set_ylim(min(ic_hepsi) - 0.35, max(ic_hepsi) + 0.35)
    a3.axhline(A.REF_A, color="#b2182b", ls="--", lw=1.2)
    a3.set_xlabel("regresyon tabanı  (× = τ₀ pencere dışı, ekstrapolasyon)")
    a3.set_ylabel("a = dφ_Γ/dτ|$_{τ_0}$")
    a3.set_title("P4 — a'nın taban konvansiyonuna bağımlılığı (her gaz ayrı)")
    a3.legend(fontsize=7.5, ncol=2)
    a3.grid(alpha=0.25)

    fig.tight_layout(rect=[0, 0, 1, 0.965])
    p = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann/158_a_kimlik.png")
    fig.savefig(p, dpi=135)
    print(f"-> {p}")


if __name__ == "__main__":
    main()
