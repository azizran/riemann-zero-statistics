"""
170 — FİGÜR: T(σ) yasası, arcsine ara-değer noktaları, c(λ) çukuru,
             çarpan ayrıştırması
====================================================================
Yeni ölçüm YOK: scratchpad/170/{K1_T.json, K2_yuzlesme.json, K0.json,
ONKAYIT_L060.json} ve scratchpad/169/K2b_Hkeskin.json okunur.
Çıktı: 170_sonlu_doyum.png
"""
import importlib
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
sys.path.insert(0, str(QM / "170_configs"))
TB = importlib.import_module("170b_T_yasasi")

SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
T_INF = np.sqrt(2 / np.pi)
C_HIP = 4 / np.pi ** 2

K1T = json.load(open(SCR / "170/K1_T.json"))
YZ = json.load(open(SCR / "170/K2_yuzlesme.json"))
K0 = json.load(open(SCR / "170/K0.json"))["gaz"]
OK = json.load(open(SCR / "170/ONKAYIT_L060.json"))

fig, ax = plt.subplots(2, 2, figsize=(13.2, 9.2))

# --- (a) T(σ) --------------------------------------------------------
a = ax[0, 0]
s = np.linspace(0.01, 5.0, 600)
a.plot(s, TB.T_sigma(s), "k-", lw=2.2, label=r"$T(\sigma)=\sqrt{\rho/\arcsin\rho}$,"
       "\n" r"$\rho=1-e^{-\sigma^2}$   (T2)+(T4)")
a.plot(s, T_INF * (1 + (np.sqrt(2) / np.pi) * np.exp(-s ** 2 / 2)), "--",
       color="tab:orange", lw=1.6,
       label=r"derin açılım (T5): $\sqrt{2/\pi}\,[1+\frac{\sqrt{2}}{\pi}e^{-\sigma^2/2}]$")
a.axhline(T_INF, color="tab:red", ls=":", lw=1.5,
          label=r"$T_\infty=\sqrt{2/\pi}=0.7979$")
ren = {"Hkeskin": "tab:blue", "L085": "tab:green", "L070": "tab:purple",
       "son": "tab:brown"}
for g, r in K1T["gaz"].items():
    for b in r["bacak"][:3]:
        a.plot(b["sigma_b"], b["T"], "o", color=ren.get(g, "grey"), ms=6,
               mec="k", mew=0.5)
    a.plot(min(5.0, r["sigma_toplam"]), T_INF, "*", color=ren.get(g, "grey"),
           ms=15, mec="k", mew=0.5)
a.annotate("ölçülen tek BACAKLAR\n(σ_b = 1.0–1.26 rad,\nmod 2π DÜZGÜN DEĞİL)",
           xy=(1.15, 0.935), xytext=(1.55, 0.975), fontsize=8.5,
           arrowprops=dict(arrowstyle="->", lw=1))
a.annotate("Ç4 TOPLAMI (σ_Σ = 4.0–4.8 rad)\nmod 2π DÜZGÜN  ⇒  $T=\\sqrt{2/\\pi}$",
           xy=(4.5, T_INF), xytext=(2.4, 0.845), fontsize=8.5,
           arrowprops=dict(arrowstyle="->", lw=1))
a.set_xlabel(r"$\sigma_b = 2\pi\tau_b\sigma_{\hat C}$  [rad]")
a.set_ylabel(r"bacak aktarımı  $T$")
a.set_title("(a) K1 — türetilen aktarım yasası ve ölçülen $\\sigma_b$'ler",
            fontsize=11)
a.set_ylim(0.77, 1.02)
a.legend(fontsize=8, loc="upper right")
a.grid(alpha=0.3)

# --- (b) arcsine ara-değer eğrisi ------------------------------------
b_ = ax[0, 1]
r = np.linspace(0.02, 1.0, 400)
b_.plot(r, TB.T2_r(r), "k-", lw=2.2,
        label=r"$T_2(r)=\frac{2}{\pi}\arcsin(r)/r$   (T3)")
b_.axhline(2 / np.pi, color="tab:red", ls=":", lw=1.5,
           label=r"$2/\pi=T_\infty^2=0.6366$")
xs = [x["r"] for x in K1T["arcsine"]]
ys = [x["T2_olc"] for x in K1T["arcsine"]]
b_.plot(xs, ys, "o", color="tab:blue", ms=8, mec="k",
        label="169 §K2.3 ÖLÇÜLEN (5 bant)")
for i, (x, y, d) in enumerate(zip(xs, ys, K1T["arcsine"])):
    b_.annotate(f"{d['fark']:+.1f}%", (x, y), textcoords="offset points",
                xytext=(-14, 10 if i % 2 else -18), fontsize=8)
b_.set_xlabel(r"$r$ = korr(bant-sınırlı model alanı, $\eta$)")
b_.set_ylabel(r"iki-bacak aktarımı  $\langle sgn\,sgn\rangle/r$")
b_.set_title("(b) K1 — T eğrisi 169'un ara-değer noktalarını vuruyor",
             fontsize=11)
b_.set_xlim(0, 1.02)
b_.set_ylim(0.60, 1.02)
b_.legend(fontsize=8.5, loc="upper left")
b_.grid(alpha=0.3)

# --- (c) c(λ) --------------------------------------------------------
c_ = ax[1, 0]
lam = [r["lam"] for r in YZ["lam"]]
cc = [r["c"] for r in YZ["lam"]]
ss = [r["s"] for r in YZ["lam"]]
hd = [r["hd1"] for r in YZ["lam"]]
c_.errorbar(lam, cc, yerr=ss, fmt="o-", color="tab:blue", ms=8, capsize=4,
            lw=2, label="ÖLÇÜLEN $c_{W_X}$ (±σ_tot)")
c_.plot(lam, hd, "s--", color="tab:red", ms=7,
        label="H-D1 (oranlı): $c(Hk)\\,\\Pi T(\\lambda)/\\Pi T(1)$")
c_.axhline(C_HIP, color="k", ls=":", lw=1.4, label=r"$4/\pi^2=0.40528$")
for k, m, col in (("P3", "^", "tab:green"), ("P4", "v", "tab:olive"),
                  ("P5", "x", "grey")):
    c_.plot(0.60, YZ["onkayit"][k], m, color=col, ms=9,
            label=f"ön kayıt {k} = {YZ['onkayit'][k]:.4f}")
p115 = SCR / "170/ONKAYIT_L115.json"
if p115.exists() and "L115" in K0 and "c_WX" in K0["L115"]:
    o115 = json.load(open(p115))["ongoru"]
    c_.errorbar([1.15], [K0["L115"]["c_WX"]], yerr=[K0["L115"]["s_tot"]],
                fmt="o", color="tab:blue", ms=8, capsize=4)
    c_.plot([1.00, 1.15], [cc[0], K0["L115"]["c_WX"]], "-", color="tab:blue",
            lw=2)
    c_.plot(1.15, o115["Q1"], "s", color="tab:red", ms=7)
    c_.plot([1.00, 1.15], [hd[0], o115["Q1"]], "--", color="tab:red", lw=1.5)
    for k, m, col in (("Q4", "^", "tab:green"), ("Q5", "v", "tab:olive"),
                      ("Q6", "x", "grey")):
        c_.plot(1.15, o115[k], m, color=col, ms=9)
    c_.annotate("ÖN KAYITLI\nλ=1.15", xy=(1.15, K0["L115"]["c_WX"]),
                xytext=(1.03, 0.372), fontsize=9,
                arrowprops=dict(arrowstyle="->", lw=1.2))
c_.annotate("ÖN KAYITLI\nÖRNEKLEM-DIŞI", xy=(0.60, cc[-1]),
            xytext=(0.585, 0.335), fontsize=9,
            arrowprops=dict(arrowstyle="->", lw=1.2))
c_.annotate("iniş DURUYOR:\nc(0.60) = c(0.70)", xy=(0.66, 0.3690),
            xytext=(0.70, 0.347), fontsize=9,
            arrowprops=dict(arrowstyle="->", lw=1))
c_.set_xlabel(r"$\lambda$  (merdiven genlik ölçeği)")
c_.set_ylabel(r"$c_{W_X}=\mathrm{KALİB}_{u2}/W_X$")
c_.set_title("(c) K2/K2′ — ön kayıtlı sınav: H-D1 −13.2σ; taban λ≈0.6, "
             "4/π² λ=1.00'de KESİLİYOR", fontsize=10.5)
c_.set_xlim(0.555, 1.20)
c_.set_ylim(0.315, 0.50)
c_.legend(fontsize=8, loc="upper right")
c_.grid(alpha=0.3)

# --- (d) çarpan ayrıştırması -----------------------------------------
d_ = ax[1, 1]
gz = ["L115", "Hkeskin", "L085", "L070", "L060"]
lm = [1.15, 1.00, 0.85, 0.70, 0.60]
h = K0["Hkeskin"]
for key, lab, col, mk in (("KALIB", r"$\mathrm{KALİB}_{u2}$ (PAY)", "tab:blue", "o"),
                          ("W_X", r"$W_X$ (PAYDA)", "tab:red", "s"),
                          ("c_WX", r"$c_{W_X}$ = pay/payda", "k", "D"),
                          ("gE", r"$g_E$ (Gram şişmesi)", "tab:green", "^"),
                          ("theta", r"$\theta$", "tab:orange", "v")):
    y = [K0[g][key] / h[key] for g in gz]
    d_.plot(lm, y, mk + "-", color=col, ms=7, lw=1.8, label=lab)
d_.plot(lm, [K0[g]["sigC"] / h["sigC"] for g in gz], "--", color="grey",
        lw=1.5, label=r"$\sigma_{\hat C}/\sigma_0$  ($\approx\lambda^{0.55}$)")
d_.axhline(1.0, color="k", lw=0.8)
d_.set_xlabel(r"$\lambda$")
d_.set_ylabel("Hkeskin'e oran")
d_.set_title("(d) K0 — λ ekseninin ADRESİ: pay artıyor, payda daha hızlı",
             fontsize=11)
d_.legend(fontsize=8, loc="center left")
d_.grid(alpha=0.3)

fig.suptitle("170 — EJDERHA: sonlu-doyum düzeltme yasası c(λ)  "
             "[H-D1 ÖLDÜ; λ-ekseni bir PAYDA (W_X) olgusudur]",
             fontsize=12.5)
fig.tight_layout(rect=(0, 0, 1, 0.965))
out = QM / "170_sonlu_doyum.png"
fig.savefig(out, dpi=135)
print(f"-> {out}")
