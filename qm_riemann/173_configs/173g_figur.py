"""
173g — FİGÜR: dokuz/on noktalı λ merdiveni ve hakem
====================================================
Yeni ölçüm YOK; `173/EGRI.json`, `173/HAKEM.json`, `173/EK.json`,
`172/G4.json` okunur.

(a) c(λ): U eğrisi, vadi λ* = 0.6486, 4/π² çizgisi, G4'ün iki c mührü
(b) M(λ): sekiz sağlıklı nokta + log-λ parabolü; G4'ün iki M mührü ve
    171'in M ailesinin λ = 0.40 öngörüleri
(c) α(λ) ve Q_E(λ): iki tümsek/plato, λ* ve λ_c
(d) ν(λ) merdiveni (ν = 1 çizgisi) ve R_bant(λ) sadakat sınırı

Kullanım: 173g_figur.py     Çıktı: 173_hakem.png
"""
import json

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt            # noqa: E402

SCR = ("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
       "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/")
QM = "/Users/ugursezen/Desktop/arin/deney/qm_riemann/"
E = json.load(open(SCR + "173/EGRI.json"))
HK = json.load(open(SCR + "173/HAKEM.json"))
EK = json.load(open(SCR + "173/EK.json"))
G4 = json.load(open(SCR + "172/G4.json"))
D9 = E["defter"]
SIRA = [g for g in ["L040", "L050", "L060", "L070", "L085", "Hkeskin",
                    "L115", "L130", "L140", "L145"] if g in D9]
lam = np.array([D9[g]["lam"] for g in SIRA])
sag = np.array([D9[g]["saglikli"] for g in SIRA])
c = np.array([D9[g]["c"] for g in SIRA])
sc = np.array([D9[g]["s_tot"] for g in SIRA])
M = np.array([D9[g]["M"] for g in SIRA])
sM = np.array([D9[g]["sM"] for g in SIRA])
al = np.array([D9[g]["alfa"] for g in SIRA])
QE = np.array([D9[g]["QE"] for g in SIRA])
Rlo = np.array([min(D9[g]["R_bant"]) for g in SIRA])
Rhi = np.array([max(D9[g]["R_bant"]) for g in SIRA])
C_HIP = 4 / np.pi ** 2
LSTAR = float(E["P2"]["kesisim"][0]) if E["P2"]["kesisim"] else 0.6486
LAMC = 1.9147 / 1.894
S, U = sag, ~sag

fig, ax = plt.subplots(2, 2, figsize=(13.5, 9.5))

# ---------------------------------------------------------------- (a) c
a = ax[0, 0]
a.errorbar(lam[S], c[S], yerr=sc[S], fmt="o-", color="#1b6ca8", ms=6,
           capsize=3, label=u"ölçülen c (sağlıklı)")
if U.any():
    a.errorbar(lam[U], c[U], yerr=sc[U], fmt="x--", color="#95a5a6", ms=8,
               capsize=3, label=u"SAĞLIKSIZ (R_bant < 0.98)")
a.axhline(C_HIP, ls=":", c="k", lw=1)
a.text(1.24, C_HIP + .004, r"$4/\pi^2$", fontsize=9)
a.axvline(LSTAR, ls="--", c="#c0392b", lw=1)
a.text(LSTAR + .02, 0.352, u"vadi λ* = %.4f\n(171: 0.6487)" % LSTAR,
       fontsize=8, color="#c0392b")
a.plot([0.40], [G4["ongoru"]["c"][0]], "v", color="#e67e22", ms=11,
       label=u"172/G4 ÇARPAN mührü")
a.plot([0.40], [G4["ongoru"]["c_dogrudan"][0]], "^", color="#16a085", ms=11,
       label=u"172/G4 DOĞRUDAN mührü")
a.plot([1.45], [G4["ongoru"]["c"][1]], "v", color="#e67e22", ms=11)
a.plot([1.45], [G4["ongoru"]["c_dogrudan"][1]], "^", color="#16a085", ms=11)
a.set_xlabel(r"$\lambda$"); a.set_ylabel("c")
a.set_title(u"(a) c(λ): U eğrisi — 4/π² bir limit değil, İKİ KEZ kesilen "
            u"bir seviye", fontsize=10)
a.legend(fontsize=7.5, loc="upper center"); a.grid(alpha=.3)

# ---------------------------------------------------------------- (b) M
a = ax[0, 1]
a.errorbar(lam[S], M[S], yerr=sM[S], fmt="o", color="#1b6ca8", ms=6,
           capsize=3, label=u"ölçülen M (sağlıklı)")
if U.any():
    a.errorbar(lam[U], M[U], yerr=sM[U], fmt="x", color="#95a5a6", ms=9,
               capsize=3, label=u"SAĞLIKSIZ")
xx = np.linspace(0.38, 1.34, 200)
p = np.array(EK["T1"]["kat"])
a.plot(xx, np.exp(np.polyval(p, np.log(xx))), "-", color="#1b6ca8", lw=1.2,
       label=u"8 sağlıklı nokta: log-λ parabolü (rms %.2f%%)" % EK["T1"]["rms"])
a.plot([0.40], [G4["ongoru"]["M"][0]], "v", color="#e67e22", ms=12,
       label=u"172/G4 ÇARPAN mührü  (+4.9σ ⇒ ÖLDÜ)")
a.plot([0.40], [G4["ongoru"]["M_dogrudan"][0]], "^", color="#16a085", ms=12,
       label=u"172/G4 DOĞRUDAN mührü (+0.9σ ⇒ AYAKTA)")
a.plot([1.45], [G4["ongoru"]["M"][1]], "v", color="#e67e22", ms=12)
a.plot([1.45], [G4["ongoru"]["M_dogrudan"][1]], "^", color="#16a085", ms=12)
for nm, mk, col in ((u"M2 kırık kuvvet (+9.3σ)", 1.1924, "#c0392b"),
                    (u"M3 varyans açığı (+10.4σ)", 1.1674, "#8e44ad"),
                    (u"M0 düz (+17.6σ)", 1.0000, "#7f8c8d")):
    a.plot([0.40], [mk], "s", color=col, ms=7, alpha=.9)
    a.annotate(nm, (0.40, mk), (0.55, mk + 0.015), fontsize=7, color=col,
               arrowprops=dict(arrowstyle="-", color=col, lw=.6))
if 1.40 in list(lam):
    j = list(lam).index(1.40)
    a.annotate(u"λ = 1.40 SAĞLIKLI ama\nparabolün %.1f%% altında (−4.6σ)"
               % abs(100 * (M[j] / np.exp(np.polyval(p, np.log(1.40))) - 1)),
               (1.40, M[j]), (1.02, 0.855), fontsize=7.5, color="#b03a2e",
               arrowprops=dict(arrowstyle="->", color="#b03a2e", lw=.8))
a.set_xlabel(r"$\lambda$"); a.set_ylabel("M")
a.set_title(u"(b) HAKEM: λ = 0.40'ta iki mühür %6.7 ayrışıyordu", fontsize=10)
a.legend(fontsize=7.5, loc="upper right"); a.grid(alpha=.3)

# ------------------------------------------------------------ (c) α, Q_E
a = ax[1, 0]
a.plot(lam[S], al[S], "o-", color="#16a085", label=r"$\alpha$ (A-genişliği)")
if U.any():
    a.plot(lam[U], al[U], "x--", color="#16a085", alpha=.45, ms=9)
a.axvline(LSTAR, ls="--", c="#c0392b", lw=1)
a.axvline(LAMC, ls="--", c="k", lw=1)
a.text(LAMC + .01, 1.18, r"$\lambda_c = 1.0109$" "\n(inşa doyumu)",
       fontsize=8)
a.axvspan(0.50, 1.15, color="#16a085", alpha=.07)
a.text(0.72, 0.28, u"171'in A-çarpanlaşma\npenceresi [0.50, 1.15]\n"
       u"— iki yandan doğrulandı", fontsize=8, color="#0e6655", ha="center")
a2 = a.twinx()
a2.plot(lam[S], QE[S], "s-", color="#1b6ca8", label=r"$Q_E$")
if U.any():
    a2.plot(lam[U], QE[U], "x--", color="#1b6ca8", alpha=.45, ms=9)
a2.set_ylabel(r"$Q_E$ (izdüşüm kusuru)", color="#1b6ca8")
a.set_xlabel(r"$\lambda$")
a.set_ylabel(r"$\alpha$", color="#16a085")
a.set_title(u"(c) α: iki yandan uçurumlu plato;  Q_E tümseği λ_c'de",
            fontsize=10)
a.grid(alpha=.3)

# ------------------------------------------------------------- (d) ν, R
a = ax[1, 1]
nu = E["nu"]
lmid = np.array([r["lam_orta"] for r in nu])
nuv = np.array([r["nu_ort"] for r in nu])
nus = np.array([r["s_ort"] for r in nu])
saglik_adim = np.array([D9[r["dan"]]["saglikli"] and D9[r["ye"]]["saglikli"]
                        for r in nu])
a.errorbar(lmid[saglik_adim], nuv[saglik_adim], yerr=nus[saglik_adim],
           fmt="o-", color="#8e44ad", capsize=3, label=r"$\nu$ (5 bant)")
if (~saglik_adim).any():
    a.errorbar(lmid[~saglik_adim], nuv[~saglik_adim], yerr=nus[~saglik_adim],
               fmt="x--", color="#95a5a6", capsize=3, ms=9,
               label=u"ν (sağlıksız adım)")
a.axhline(1.0, ls=":", c="k", lw=1)
a.text(0.98, 1.10, r"$\nu = 1 \Leftrightarrow dc/d\lambda = 0$", fontsize=8)
a.axvline(LSTAR, ls="--", c="#c0392b", lw=1)
a.set_xlabel(r"$\lambda$"); a.set_ylabel(r"$\nu$", color="#8e44ad")
a3 = a.twinx()
a3.plot(lam, Rlo, "v-", color="#d35400", ms=5, lw=1, label=u"R_bant (min)")
a3.plot(lam, Rhi, "^-", color="#d35400", ms=5, lw=1, alpha=.6)
a3.axhline(0.98, ls=":", c="#d35400", lw=1)
a3.set_ylabel(u"R_bant (sadakat)", color="#d35400")
kk = EK["T3"]["kesis1"]
if kk:
    a3.text(1.02, 1.05, u"R_bant = 1 → λ ≈ %.2f–%.2f\n"
            u"(sadakat sınırı)" % (kk[0][1], kk[-1][1]), fontsize=8,
            color="#d35400")
a.set_title(u"(d) ν(λ) ve SADAKAT SINIRI R_bant(λ)", fontsize=10)
a.legend(fontsize=7.5, loc="upper left"); a.grid(alpha=.3)

fig.suptitle(u"173 — HAKEM GAZLARI: λ = 0.40 doğrudan yolu seçti; "
             u"λ ≈ 1.4'te inşa sadakati bitiyor", fontsize=13)
fig.tight_layout(rect=(0, 0, 1, 0.96))
fig.savefig(QM + "173_hakem.png", dpi=130)
print("-> %s173_hakem.png" % QM)
