"""
155 — FİGÜR: τ₀'ın kimliği (6 panel)
"""
import importlib
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
A = importlib.import_module("155_analiz")
RP = importlib.import_module("155_rapor")
HERE = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")

REN = {"son": "#1f4e9c", "orta": "#3f7fd0", "dusuk": "#7fb2f0",
       "keskin": "#c0392b", "A4": "#e07b39", "N5": "#7e57c2",
       "N5z": "#b39ddb", "J14": "#2e7d32", "J26": "#81c784",
       "P1": "#8d6e63", "P0": "#bdbdbd"}


def main():
    D = RP.tum_dosyalar()
    ciftler = sorted({(v, t) for (v, t, _) in D})
    HAV = {c: RP.havuz(D, *c) for c in ciftler}
    fig, ax = plt.subplots(2, 3, figsize=(16.5, 9.2))

    # --- 1: φ(τ) sıfır civarı, taban 0.40
    a = ax[0, 0]
    for v in ("son", "orta", "dusuk", "keskin", "A4", "N5", "J14", "P1"):
        j = RP.sec(D, v, 0.40)
        if not j:
            continue
        B = [b for b in j["bantlar"] if b.get("olculdu") and
             0.41 < b["tau"] < 0.60]
        if not B:
            continue
        a.errorbar([b["tau_eff"] for b in B], [b["phi"] for b in B],
                   yerr=[b["sPhi_jk"] for b in B], marker="o", ms=3.5, lw=1.2,
                   color=REN.get(v, "k"), label=v, capsize=2)
    a.axhline(0, color="k", lw=0.8)
    a.axvline(0.5, color="k", ls=":", lw=1.0)
    a.set_xlim(0.41, 0.60); a.set_ylim(-1.1, 1.0)
    a.set_xlabel("τ"); a.set_ylabel("φ_Γ  [rad]")
    a.set_title("1) φ_Γ(τ) — taban 0.40; nokta çizgi τ=½")
    a.legend(fontsize=7, ncol=2); a.grid(alpha=0.25)

    # --- 2: sıfırın yakın planı
    a = ax[0, 1]
    for v in ("son", "orta", "dusuk", "keskin", "A4"):
        j = RP.sec(D, v, 0.40)
        if not j:
            continue
        B = [b for b in j["bantlar"] if b.get("olculdu") and
             0.465 < b["tau"] < 0.535]
        a.errorbar([b["tau_eff"] for b in B], [b["phi"] for b in B],
                   yerr=[b["sPhi_jk"] for b in B], marker="o", ms=5, lw=1.4,
                   color=REN.get(v, "k"), label=v, capsize=3)
    a.axhline(0, color="k", lw=0.8); a.axvline(0.5, color="k", ls=":", lw=1.0)
    a.set_xlabel("τ"); a.set_ylabel("φ_Γ"); a.grid(alpha=0.25)
    a.set_title("2) sıfır geçişi — gerçek (mavi) ≠ sentetik (kırmızı)")
    a.legend(fontsize=8)

    # --- 3: taban taraması
    a = ax[0, 2]
    for v in ("son", "orta", "dusuk", "keskin", "A4", "N5"):
        pts = []
        for (vv, t) in ciftler:
            if vv != v or t > 0.47:
                continue
            d = RP.t0(HAV[(vv, t)], 0.46, 0.56, 0.02, 2)
            if d:
                pts.append((t, d["tau0"], d["s_jk"]))
        if len(pts) < 2:
            continue
        pts.sort()
        a.errorbar([p[0] for p in pts], [p[1] for p in pts],
                   yerr=[p[2] for p in pts], marker="s", ms=5, lw=1.5,
                   color=REN.get(v, "k"), label=v, capsize=3)
    a.axhline(0.5, color="k", ls=":", lw=1.0)
    a.set_xlabel("regresyon tabanı τ_taban"); a.set_ylabel("τ₀")
    a.set_title("3) EN BÜYÜK SİSTEMATİK: τ₀ ↔ taban\n(gerçekte kayıyor, "
                "sentetikte kaymıyor)")
    a.legend(fontsize=8); a.grid(alpha=0.25)

    # --- 4: H-L sınavı
    a = ax[1, 0]
    for tb, mk in ((0.40, "o"), (0.46, "s")):
        pts = []
        for v in ("son", "orta", "dusuk"):
            j = RP.sec(D, v, tb)
            if not j:
                continue
            d = RP.t0(HAV[(v, tb)], 0.46, 0.56, 0.02, 2)
            if d:
                pts.append((1 / j["L"], d["tau0"], max(d["s_jk"], 2e-4)))
        if not pts:
            continue
        pts.sort()
        a.errorbar([p[0] for p in pts], [p[1] for p in pts],
                   yerr=[10 * p[2] for p in pts], marker=mk, ms=6, lw=1.5,
                   capsize=3, label=f"ölçüm, taban {tb} (hata ×10)")
    x = np.linspace(0.080, 0.098, 50)
    a.plot(x, 0.5 + 0.18 * x, "r--", lw=1.6, label="H-L: ½ + 0.18/L")
    a.axhline(0.5, color="k", ls=":", lw=1.0)
    a.set_xlabel("1/L"); a.set_ylabel("τ₀")
    a.set_title("4) H-L SINAVI — τ₀, L'ye asılı DEĞİL")
    a.legend(fontsize=8); a.grid(alpha=0.25)

    # --- 5: H-N — δ vs σ_ds²; N5 ↔ J14 eşleşmiş çifti vurgulu
    a = ax[1, 1]
    P = {}
    for v in [x for x in RP.SIRA if "S0." not in x]:
        j = RP.sec(D, v, 0.40)
        if not j or (v, 0.40) not in HAV:
            continue
        d = RP.t0(HAV[(v, 0.40)], 0.42, 0.56, 0.02, 2)
        if not d:
            continue
        P[v] = (j["s_ds"], d["tau0"] - 0.5, max(d["s_jk"], 2e-4))
        a.errorbar(P[v][0], P[v][1], yerr=P[v][2], marker="o", ms=8,
                   color=REN.get(v, "k"), capsize=3)
        dy = {"son": 6, "orta": -14, "dusuk": -24}.get(v, 4)
        a.annotate(v, (P[v][0], P[v][1]), fontsize=8, xytext=(6, dy),
                   textcoords="offset points")
    if "N5" in P and "J14" in P:
        a.annotate("", xy=P["N5"][:2], xytext=P["J14"][:2],
                   arrowprops=dict(arrowstyle="<->", color="k", lw=1.6))
        a.text(0.5 * (P["N5"][0] + P["J14"][0]),
               0.5 * (P["N5"][1] + P["J14"][1]) + 0.0012,
               "AYNI σ_ds²\nfarklı τ₀", fontsize=8, ha="center")
    a.axhline(0, color="k", ls=":", lw=1.0)
    a.axvline(0.1674, color="#1f4e9c", ls="--", lw=1.0, alpha=0.6)
    a.set_xlabel("σ_ds²  (adım varyansı; kesikli mavi = gerçek)")
    a.set_ylabel("δ = τ₀ − ½")
    a.set_title("5) H-N NİCEL — δ momente asılı DEĞİL\n(N5 ile J14 aynı "
                "momentte, τ₀ farklı)")
    a.grid(alpha=0.25)

    # --- 6: τ₀ vs σ_η² (taban etkisinin ölçeklenmesi)
    a = ax[1, 2]
    for v in ("son", "keskin", "A4"):
        pts = []
        for (vv, t) in ciftler:
            if vv != v or t > 0.47:
                continue
            d = RP.t0(HAV[(vv, t)], 0.46, 0.56, 0.02, 2)
            j = RP.sec(D, vv, t)
            if d:
                pts.append((j["s_eta"], d["tau0"], t))
        if len(pts) < 2:
            continue
        pts.sort()
        a.plot([p[0] for p in pts], [p[1] for p in pts], "o-", ms=6,
               color=REN.get(v, "k"), label=v)
        for p in pts:
            a.annotate(f"{p[2]:.2f}", (p[0], p[1]), fontsize=7,
                       xytext=(3, -10), textcoords="offset points")
        if len(pts) >= 3:
            b_, a_ = np.polyfit([p[0] for p in pts], [p[1] for p in pts], 1)
            xx = np.linspace(0, max(p[0] for p in pts), 20)
            a.plot(xx, a_ + b_ * xx, "--", lw=1.0, color=REN.get(v, "k"),
                   alpha=0.6)
            a.plot([0], [a_], "*", ms=14, color=REN.get(v, "k"))
    a.axhline(0.5, color="k", ls=":", lw=1.0)
    a.set_xlabel("σ_η²  (artık arka plan; etiketler = taban)")
    a.set_ylabel("τ₀")
    a.set_title("6) τ₀ ↔ arka plan seviyesi\n(★ = σ_η²→0 ekstrapolasyonu)")
    a.legend(fontsize=8); a.grid(alpha=0.25)

    fig.suptitle("155 — τ₀ (bond-transfer fazının sıfır-geçişi) kimlik "
                 "kampanyası", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    p = HERE / "155_tau0_kampanya.png"
    fig.savefig(p, dpi=125)
    print(f"-> {p}")


if __name__ == "__main__":
    main()
