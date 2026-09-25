# -*- coding: utf-8 -*-
"""
199fig_4 — NOT 7 FİGÜRÜ 4: fig_n7_mirror_en.{png,pdf}
=======================================================
Mühürlü 197 sonucunu (H-197a BK TUTAR, H-197b TUTAR, H-197c TUTAR) gösterir; YENİ ÖLÇÜM YOK.
197c_hukum.py importlib ile yüklenir (DÜZENLENMEZ); yalnız saf fonksiyonları çağrılır:
kinematik32, Bicim, f_oku, hesap(tam=False) — yani mühürlü harita_omega32.npz profili
üzerinde ön-kayıtlı doğrusal profil uyumu (deterministik, M6 simülasyonu YOK, dosya YAZMAZ).
Yeniden hesaplanan R_x, R̄_Z, ρ, ψ HUKUM_197.json ile assert edilir.
(a) kör bant [−2.12, −1.58]: κ_Σ ± jk se (32 blok, K_düz), M_BK uyumu, BK öngörüsü
    (A_x = s/x; taban ve sabit çeyrek aile M_BK uyumundan — 197c.ciz ile AYNI kurgu).
(b) R_x = A_x/(s/x): tam sayı (x = 5…8) ve buçuklu (11/2, 13/2, 15/2) aileler, M_BK uyumu
    (dolu) ve aynı veri M_Rg ile (gri, KAYIT); BK = 1 ve rakip R^g(x) = g(L − log x)/ḡ.
Girdi: scratchpad/197/{harita_omega32.npz, HUKUM_197.json, F_KALIBRASYON.json},
       scratchpad/190/zincir_dusuk/eta_dusuk_t0.4_c4000.npz (kinematik)
"""
import importlib.util
import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True           # içe aktarılan betikler için __pycache__ YAZMA
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
S197 = QM / "scratchpad" / "197"
CIKTI = QM / "fig_n7_mirror_en"

spec = importlib.util.spec_from_file_location("h197c", QM / "197_configs" / "197c_hukum.py")
c = importlib.util.module_from_spec(spec)
sys.modules["h197c"] = c
spec.loader.exec_module(c)

MAVI, KIZIL, YESIL, MOR = "#0072B2", "#D55E00", "#009E73", "#AA3377"
INK, GRI = "#1a1a1a", "#8c8c8c"
plt.rcParams.update({"font.size": 9, "axes.labelsize": 10, "xtick.labelsize": 9,
                     "ytick.labelsize": 9, "legend.fontsize": 9, "axes.linewidth": 0.8,
                     "mathtext.fontset": "dejavusans", "savefig.dpi": 300,
                     "pdf.fonttype": 42, "ps.fonttype": 42})


def harf(ax, s, x=0.015, y=0.975):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=10, fontweight="bold",
            va="top", ha="left", bbox=dict(fc="white", ec="none", pad=1.2), zorder=6)


# ---- mühürlü profil + ön-kayıtlı uyum (197c.mod_hukum ile AYNI girdi hazırlığı) ----
HK = json.load(open(S197 / "HUKUM_197.json"))
Hm = np.load(S197 / "harita_omega32.npz")
assert str(Hm["birincil"]) == "duz"
tam = Hm["tam"]
J = Hm["J"][tam]
prof = Hm["kappa_sigma_duz"][:, tam]
mid, kb, Lb = c.kinematik32()
assert np.array_equal(kb, Hm["kb"]) and np.array_equal(Lb, Hm["L_b"])
H = c.Bicim(mid, kb, Lb, Hm["kapsama"][:, tam], J)
f = c.f_oku()
o = c.hesap(prof, H, J, tam=False, f=f)

# ---- denetim: HUKUM_197 ile birebir ----
nk = HK["nicelikler"]
assert abs(o["_s"][0] - nk["s"]["deger"]) < 1e-12
for mo in ("BK", "Rg"):
    for ad, v in nk["modeller"][mo]["R"].items():
        assert abs(o["modeller"][mo]["R"][ad]["deger"] - v["deger"]) < 1e-10, (mo, ad)
        assert abs(o["modeller"][mo]["R"][ad]["se"] - v["se"]) < 1e-10, (mo, ad)
    for q in ("R_Z", ):
        assert abs(o["modeller"][mo][q]["deger"] - nk["modeller"][mo][q]["deger"]) < 1e-10
mB, mR = o["modeller"]["BK"], o["modeller"]["Rg"]
fZ, frho = f["Z"], f["rho"]
RZ, sZ = mB["R_Z"]["deger"], fZ * mB["R_Z"]["sigma_Z"]
rho, srho = mB["rho"]["deger"], frho * mB["rho"]["se"]
psi, spsi = mB["psi"]["deger"], mB["psi"]["se"]
print(f"s = {o['_s'][0]:.5f} ± {c.jk(o['_s'][1:]):.5f}")
print("R_x (M_BK):", {a: f"{v['deger']:.3f}±{v['se']:.3f}" for a, v in mB["R"].items()})
print(f"R̄_Z^BK = {RZ:.3f} ± {sZ:.3f} (σ_eff; jk {mB['R_Z']['sigma_Z']:.4f}); "
      f"R̄_Z^Rg = {mR['R_Z']['deger']:.3f} vs R̄^g = {mR['R_g_Z']:.3f}; "
      f"(R̄^g − R̄_Z^Rg)/σ_eff^Rg = {(mR['R_Z']['deger'] - mR['R_g_Z']) / (fZ * mR['R_Z']['sigma_Z']):.1f}")
print(f"ρ = {rho:.3f} ± {srho:.3f} (σ_eff); ψ = {psi:.3f} ± {spsi:.3f}; "
      f"β^BK = {mB['beta']['deger']:.2f} ± {mB['beta']['se']:.2f}")
print("hüküm (HUKUM_197):", HK["hukum"])

fk = o["_fit"]
s0 = o["_s"][0]
cc = J * c.DW
se = np.array([c.jk(prof[1:, j]) for j in range(len(J))])
m = fk["BK"]["m"]
thB = fk["BK"]["th"][0]
ong = thB[-2] + thB[-1] * cc[m] + fk["BK"]["sabit"][0]
for x in c.KOR["ana"]:
    ong = ong + s0 / float(x) * H(-np.log(float(x)))[0, m]
fitB = fk["BK"]["model_toplam"][0]
rms_ong = float(np.sqrt(np.mean((prof[0, m] - ong) ** 2)))
rms_fit = float(np.sqrt(np.mean((prof[0, m] - fitB) ** 2)))
print(f"kör bant: {int(m.sum())} dilim; rms(veri − M_BK) = {rms_fit:.2e}, "
      f"rms(veri − BK öngörüsü) = {rms_ong:.2e}, medyan jk se = {np.median(se[m]):.2e}")

# ---- figür ----
fig, (ax, bx) = plt.subplots(2, 1, figsize=(6.5, 6.3), layout="constrained",
                             gridspec_kw={"height_ratios": [1.1, 1.0]})
K3 = 1e3
ax.axhline(0, color=GRI, lw=0.6)
for x in c.KOR["ana"]:
    ax.axvline(-np.log(float(x)), color=INK, lw=0.6, alpha=0.35,
               ls="-" if x.denominator == 1 else (0, (1, 2)), zorder=0)
ax.errorbar(cc[m], K3 * prof[0, m], yerr=K3 * se[m], fmt="o", ms=3.6, color=INK,
            elinewidth=0.9, capsize=0, zorder=4, label=r"measured $\pm$ jk s.e.")
ax.plot(cc[m], K3 * fitB, color=MAVI, lw=1.6, zorder=3, label=r"fit $M_{\rm BK}$ (free $A_x$)")
ax.plot(cc[m], K3 * ong, color=KIZIL, lw=1.5, ls=(0, (4, 1.8)), zorder=3,
        label=r"BK prediction $A_x=s/x$")
ax.set_xlim(-2.135, -1.565)
ax.set_ylim(-2.75, 5.0)
ax.set_xlabel(r"$\Delta\omega=\omega'-L_{\rm loc}$")
ax.set_ylabel(r"$\kappa_\Sigma\;(\times10^{-3})$")
ha_, la_ = ax.get_legend_handles_labels()
ax.legend([ha_[i] for i in (2, 0, 1)], [la_[i] for i in (2, 0, 1)], loc="lower center", ncol=3,
          frameon=True, framealpha=0.95, edgecolor="none", borderaxespad=0.3,
          handlelength=2.2, columnspacing=1.2)
top = ax.twiny()
top.set_xlim(ax.get_xlim())
top.set_xticks([-np.log(float(x)) for x in c.KOR["ana"]])
top.set_xticklabels([str(x) for x in c.KOR["ana"]])
top.set_xlabel(r"$x$  ($\Delta\omega={-}\log x$)", labelpad=3)
harf(ax, "(a)")

xs = [float(x) for x in c.KOR["ana"]]
Zi = [i for i, x in enumerate(c.KOR["ana"]) if x.denominator == 1]
Hi = [i for i, x in enumerate(c.KOR["ana"]) if x.denominator == 2]
RB = np.array([mB["R"][c.ad_x(x)]["deger"] for x in c.KOR["ana"]])
eB = np.array([mB["R"][c.ad_x(x)]["se"] for x in c.KOR["ana"]])
RR = np.array([mR["R"][c.ad_x(x)]["deger"] for x in c.KOR["ana"]])
eR = np.array([mR["R"][c.ad_x(x)]["se"] for x in c.KOR["ana"]])
xx = np.linspace(4.7, 8.3, 200)
Rg = c.g_rg(c.L_G - np.log(xx)) / mR["g_bar"]
bx.axhline(1.0, color=KIZIL, lw=1.5, ls=(0, (4, 1.8)), zorder=1, label=r"BK: $R_x=1$")
bx.plot(xx, Rg, color=YESIL, lw=1.6, ls=(0, (1, 1.4)), zorder=1,
        label=r"line-density rival $R^g(x)$")
X = np.array(xs)
bx.errorbar(X[Zi] + 0.07, RR[Zi], yerr=eR[Zi], fmt="s", ms=3.6, color=GRI, mfc="white",
            mew=0.9, elinewidth=0.8, capsize=0, zorder=2, label=r"same data, fit $M_{\rm Rg}$")
bx.errorbar(X[Hi] + 0.07, RR[Hi], yerr=eR[Hi], fmt="s", ms=3.6, color=GRI, mfc="white",
            mew=0.9, elinewidth=0.8, capsize=0, zorder=2)
bx.errorbar(X[Zi], RB[Zi], yerr=eB[Zi], fmt="o", ms=5.4, color=MAVI, mec="white", mew=0.6,
            elinewidth=1.2, capsize=2.5, zorder=4, label=r"integer $x$ ($M_{\rm BK}$)")
bx.errorbar(X[Hi], RB[Hi], yerr=eB[Hi], fmt="D", ms=5.0, color=MOR, mfc="white", mew=1.2,
            elinewidth=1.2, capsize=2.5, zorder=4, label=r"half-integer $x$ ($M_{\rm BK}$)")
bx.set_xlim(4.7, 8.3)
bx.set_ylim(0.52, 1.24)
bx.set_xticks(xs)
bx.set_xticklabels([str(x) for x in c.KOR["ana"]])
bx.set_xlabel(r"$x$  (satellite at $\Delta\omega={-}\log x$)")
bx.set_ylabel(r"$R_x=A_x\,/\,(s/x)$")
bx.text(0.985, 0.955, rf"$\bar R_{{\mathbb{{Z}}}}={RZ:.3f}\pm{sZ:.3f}$" + "\n" +
        rf"$\rho={rho:.3f}\pm{srho:.3f}$", transform=bx.transAxes, ha="right", va="top",
        fontsize=9, bbox=dict(fc="white", ec="none", pad=1.5))
hh, ll = bx.get_legend_handles_labels()
order = [3, 4, 2, 0, 1]
bx.legend([hh[i] for i in order], [ll[i] for i in order], loc="upper center",
          bbox_to_anchor=(0.5, -0.24), ncol=3, frameon=False, borderaxespad=0.0,
          handlelength=2.0, columnspacing=1.1)
harf(bx, "(b)")

fig.savefig(CIKTI.with_suffix(".png"), dpi=300)
fig.savefig(CIKTI.with_suffix(".pdf"))
print(f"-> {CIKTI}.png / .pdf")
