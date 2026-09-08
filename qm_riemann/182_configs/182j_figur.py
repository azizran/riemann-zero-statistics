# -*- coding: utf-8 -*-
"""
182j — FİGÜR: GAUSS-ALTININ ANATOMİSİ (6 panel)
===============================================
(a) 3-bacak kırpma aktarımı ρ₃, n=8 + çözünürlük saçılımı ve eşik.
(b) Bacak bacak ρ: E kör olmayan tek bacak (Hkeskin / son / VF n=8).
(c) ÇİZGİ / BANT ayrışımı: log açığın κ_E (çizgi) ve β_E (bant) payları.
(d) ρ(E) vs τ_eff — bant profili ve eğim (E ↔ Xa/Xb).
(e) SURROGATE merdiveni: ρ(E) → ρ*(E) (yalnız E) → ρ°(E) (üçü de).
(f) Anti-girişimin adresi: R_η, Kov(artık,η)/Var, Var(çizgi)/P.

Kullanım: 182j_figur.py
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
S182, S180, S169 = SCR / "182", SCR / "180", SCR / "169"

HK8 = json.load(open(S182 / "K2_HAKEM_n8.json"))
OK180 = json.load(open(S180 / "ONKAYIT_K0.json"))["H180a"]["capalar"]
HYP = json.load(open(S182 / "K2_HIPOTEZ.json"))
gz = [f"VF{i}" for i in range(1, 9)]
A = {g: json.load(open(S182 / f"ANATOMI_{g}.json")) for g in gz}
G1, G3 = HYP["GAUSS1"], HYP["GAUSS3"]
R8 = HK8["R8"]

fig, ax = plt.subplots(2, 3, figsize=(19.5, 11))

# ------------------------------------------------------------- (a)
a = ax[0, 0]
nm = ["son", "Hkeskin"] + gz
vl = [OK180["rho3_son"], OK180["rho3_Hkeskin"]] + \
     [R8["rho3"][g] for g in gz]   # son/Hkeskin: 180a'da dondurulmuş çapa
cl = ["#c0392b", "#2c3e50"] + ["#2980b9"] * 8
a.bar(range(len(vl)), vl, color=cl)
a.axhline(G3, color="k", ls="--", lw=1.8,
          label=f"Gauss $(2/\\pi)^{{3/2}}$ = {G3:.4f}")
a.axhline(R8["rho3_bar"], color="#2980b9", ls=":", lw=1.5,
          label=f"$\\langle\\rho_3\\rangle_{{VF}}$ = {R8['rho3_bar']:.4f}"
                f"  ($z_G$ = {R8['z_G']:+.1f})")
a.set_xticks(range(len(nm)))
a.set_xticklabels(nm, rotation=45, ha="right", fontsize=8)
a.set_ylim(0.40, 0.58)
a.set_ylabel(r"$\rho_3$ (3 bacak kırpma aktarımı)")
a.set_title(f"(a) n=8: Gauss'un %{100*(1-R8['rho3_bar']/G3):.1f} ALTI\n"
            f"saçılım {R8['sacilim']:.4f} > eşik "
            f"{R8['esik_sacilim']:.4f} ⇒ hakem HÜKÜMSÜZ (aynen)",
            fontsize=10)
a.legend(fontsize=8, loc="upper right")
a.grid(alpha=.3, axis="y")

# ------------------------------------------------------------- (b)
a = ax[0, 1]
LEG = ["100_E", "010_Xa", "001_Xb", "011_XaXb", "110_EXa", "111_hepsi"]
ETI = ["E", "Xa", "Xb", "XaXb", "EXa", "hepsi"]
x = np.arange(len(LEG))
for g, c, mk in (("Hkeskin", "#2c3e50", "s"), ("son", "#c0392b", "o")):
    if g in R8["tablo"]:
        a.plot(x, [R8["tablo"][g][k] for k in LEG], mk + "-", color=c,
               label=g, ms=7)
mu = np.array([HK8["bacak_n8"][k]["ort"] for k in LEG])
sd = np.array([HK8["bacak_n8"][k]["sd"] for k in LEG])
a.errorbar(x, mu, yerr=sd, fmt="^-", color="#2980b9", ms=8, capsize=4,
           label="VF ailesi (n=8, ±sd)")
a.plot(x, [G1 ** k.split("_")[0].count("1") for k in LEG], "k--", lw=1.8,
       label=r"Gauss $(\sqrt{2/\pi})^k$")
a.set_xticks(x); a.set_xticklabels(ETI)
a.set_ylabel(r"$\rho$(varyant)")
a.set_title("(b) Bütün hareket E bacağında; Xa/Xb karıştırmaya kör\n(Hkeskin/son bacak satırları 180-T5'te; 182 yalnız VF ailesini ölçtü)",
            fontsize=10)
a.legend(fontsize=8); a.grid(alpha=.3)

# ------------------------------------------------------------- (c)
a = ax[0, 2]
kE = np.array(HYP["kappa_E"]); rE = np.array(HYP["rho_E"])
lc = np.log(kE) - np.log(G1)          # ÇİZGİ payı (log)
lt = np.log(rE) - np.log(G1)          # toplam açık (log)
lb = np.log(np.array(HYP["beta_E"]))  # BANT payı (log)
w = np.arange(8)
a.bar(w - .2, lc, .4, color="#e67e22", label=r"ÇİZGİ: $\log\kappa_E-\log\sqrt{2/\pi}$")
a.bar(w + .2, lb, .4, color="#8e44ad", label=r"BANT: $\log\beta_E$")
a.plot(w, lt, "ko-", ms=6, label=r"TOPLAM: $\log\rho(E)-\log\sqrt{2/\pi}$")
a.axhline(0, color="k", lw=.8)
a.set_xticks(w); a.set_xticklabels(gz, rotation=45, fontsize=8)
a.set_ylabel("log açık")
a.set_title(f"(c) ÇİZGİ/BANT ayrışımı — PAY_ÇİZGİ = "
            f"{np.mean(HYP['pay_cizgi']):+.3f}\naçığın tamamı BANT'ta "
            f"(çizgi ters yönde çalışıyor)", fontsize=10)
a.legend(fontsize=8); a.grid(alpha=.3, axis="y")

# ------------------------------------------------------------- (d)
a = ax[1, 0]
for g in gz:
    d = A[g]
    a.plot(d["tau"], d["bant_rhoE"], "o-", ms=4, alpha=.7, lw=1.2,
           color="#2980b9")
for g, c in (("Hkeskin", "#2c3e50"), ("son", "#c0392b")):
    p = S169 / f"K2b_{g}.json"
    if p.exists():
        dd = json.load(open(p))
        a.plot([b["tau_eff"] for b in dd["bant"]],
               [b["oran"]["100_E"] for b in dd["bant"]], "s-", color=c,
               ms=6, lw=2, label=g)
a.axhline(G1, color="k", ls="--", lw=1.8, label=r"Gauss $\sqrt{2/\pi}$")
a.plot([], [], "o-", color="#2980b9", label="VF1..VF8")
a.set_xlabel(r"$\tau_{\rm eff}$"); a.set_ylabel(r"$\rho(E)$")
a.set_title(f"(d) $\\rho(E)$ bant profili — eğim "
            f"$s_E$ = {np.mean(HYP['egim_E']):+.3f}\n"
            f"($s_{{Xa}}$ = {np.mean(HYP['egim_Xa']):+.3f}, "
            f"$s_{{Xb}}$ = {np.mean(HYP['egim_Xb']):+.3f})", fontsize=10)
a.legend(fontsize=8); a.grid(alpha=.3)

# ------------------------------------------------------------- (e)
a = ax[1, 1]
seri = [("ρ(E)\n(169 kestirimcisi)", HYP["rho_E"], "#2980b9"),
        ("ρ_J(E)\n(bağımsız kest.)", HYP["rho_J_E"], "#16a085"),
        ("ρ*(E)\nyalnız E surrogate", HYP["rho_Esur"], "#8e44ad"),
        ("ρ°(E)\nüçü de surrogate", HYP["rho_sur_E"], "#e67e22")]
for i, (et, v, c) in enumerate(seri):
    v = np.array(v)
    a.errorbar([i], [v.mean()], yerr=[v.std(ddof=1)], fmt="o", ms=11,
               color=c, capsize=6)
    a.scatter([i] * 8, v, s=16, color=c, alpha=.45)
a.axhline(G1, color="k", ls="--", lw=1.8, label=r"Gauss $\sqrt{2/\pi}$")
a.set_xticks(range(4)); a.set_xticklabels([s[0] for s in seri], fontsize=8)
a.set_ylabel(r"$\rho(E)$")
a.set_title(f"(e) Surrogate merdiveni — geri kazanım "
            f"Rec = {np.mean(HYP['Rec']):+.2f}\n"
            f"E'nin fazı silinince açık kapanıyor ⇒ H-G2", fontsize=10)
a.legend(fontsize=8); a.grid(alpha=.3, axis="y")

# ------------------------------------------------------------- (f)
a = ax[1, 2]
Re = np.array(HYP["R_eta"]); kv = np.array(HYP["kov_art_over_var"])
n = len(Re)
a.bar(np.arange(n) - .2, Re, .4, color="#2980b9", label=r"$R_\eta$")
a.bar(np.arange(n) + .2, kv, .4, color="#c0392b",
      label=r"Kov$(\eta_{\rm artık},\eta)/$Var$(\eta)\ \equiv 1-R_\eta$")
a.axhline(1.0, color="k", ls="--", lw=1.5, label=r"$R_\eta=1$ (girişimsiz)")
a.axhline(0.0, color="k", lw=.8)
a.axhline(1.26504, color="#2c3e50", ls=":", lw=1.6,
          label=r"$R_\eta$(Hkeskin) = 1.265")
a.set_xticks(range(n)); a.set_xticklabels(gz[:n], rotation=45, fontsize=8)
a.set_title(f"(f) ANTİ-GİRİŞİMİN ADRESİ: "
            f"$\\langle R_\\eta\\rangle$ = {Re.mean():.3f} < 1\n"
            f"artık η ile POZİTİF ilişkili; "
            f"Var(çizgi)/P = "
            f"{np.mean([A[g]['R_eta']['Var_cizgi']/A[g]['R_eta']['P'] for g in gz if A[g]['R_eta']]):.3f} > 1",
            fontsize=10)
a.legend(fontsize=8); a.grid(alpha=.3, axis="y")

fig.suptitle("182 — GAUSS-ALTININ ANATOMİSİ: karışık (kilitsiz) vekil gaz, "
             "n = 8 tohum  —  H-G1 %s / H-G2 %s / H-G3 %s"
             % (HYP["hukum"]["H_G1"], HYP["hukum"]["H_G2"],
                HYP["hukum"]["H_G3"]), fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.965])
p = QM / "182_gauss_alti.png"
fig.savefig(p, dpi=125)
print("->", p)
