"""154 — figür: R(τ) iki pencere, adaylar, ρ, taban taraması, sentetik tanık."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.special import erfcx

OUT = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/154")
HERE = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
TWO_PI = 2 * np.pi


def yuk(p):
    j = json.loads((OUT / p).read_text())
    return j, [x for x in j["bantlar"] if x.get("olculdu")]


def iyi(b):
    # geçerli ölçüm: faz eşleşti, |Γ| patlamadı VE faz dal kesiğinden
    # SARMADI. Sarma imzası: φ monoton artarken aniden büyük NEGATİF
    # (sentetik gazın üst bantlarında −2.9 / −2.7 / −2.5).
    return [x for x in b if x["artik"] < 0.02 and x["absG"] < 1.5
            and x["phi"] > -2.0]


def R_erfc(t, tc=0.68, dl=0.125, p=2.0):
    return p / (np.pi**1.5 * dl) / erfcx((np.asarray(t, float) - tc) / dl)


Y = json.loads((OUT / "yarisma.json").read_text())
ft = np.array(Y["ince_tau"])
tg = np.linspace(0.52, 0.80, 300)
fig, ax = plt.subplots(2, 3, figsize=(18, 9.5))

# ---- (1) R(τ) + sıfır-parametreli adaylar ----
a = ax[0, 0]
a.errorbar(ft, Y["ince_R_son"], yerr=Y["ince_sig_son"], fmt="o", ms=5,
           color="C0", label="R ölçülen — son 300k")
a.errorbar(ft + .002, Y["ince_R_orta"], yerr=Y["ince_sig_orta"], fmt="s",
           ms=5, mfc="none", color="C1", label="R ölçülen — orta 300k")
a.plot(Y["tau"], Y["Rg_son"], "kv", ms=7, mfc="none",
       label=r"$R_{gauss}=\varphi/(A^2\sigma_\Delta^2)$")
a.plot(tg, R_erfc(tg), "--", color="C3", lw=1.6,
       label=r"(d) $-d\ln\rho_{erfc}/dA$  [0 par]")
a.axhline(6.20 / TWO_PI, ls=":", color="C2",
          label=r"(a) sabit $6.20/2\pi=0.99$")
c = np.array(Y["lnrho_kubik_son"])
a.plot(tg, -np.polyval(np.polyder(c), tg) / TWO_PI, "-.", color="C4",
       label=r"(f1) $-d\ln\rho_{ölç}/dA$  [0 par]")
a.plot(tg, -np.polyval(np.polyder(c), tg) / np.pi, "-", color="C5", lw=1.3,
       label=r"(f) $2\times$ aynısı [0 par]")
a.set_ylim(-0.2, 5); a.set_xlabel(r"$\tau$"); a.set_ylabel("R")
a.set_title("(1) sıfır-parametreli adaylar — hepsi ıskalıyor")
a.legend(fontsize=7.5, loc="upper left"); a.grid(alpha=.3)

# ---- (2) hayatta kalan 2-parametreli formlar + artıklar ----
a = ax[0, 1]
a.errorbar(ft, Y["ince_R_son"], yerr=Y["ince_sig_son"], fmt="o", ms=5,
           color="C0", label="R ölçülen (son)")
a.errorbar(ft + .002, Y["ince_R_orta"], yerr=Y["ince_sig_orta"], fmt="s",
           ms=5, mfc="none", color="C1", label="R ölçülen (orta)")
S = {s["ad"][:2].strip(): s for s in (Y.get("skor_syst") or Y["skor18"])}
if "c" in S:
    R8, w8 = S["c"]["theta"]
    a.plot(tg, R8 * (1 - np.exp(-(tg - 0.52) / w8)), "-", color="C3", lw=1.6,
           label=rf"(c) $R_\infty(1-e^{{-(\tau-0.52)/w}})$, "
                 rf"$R_\infty$={R8:.2f}, w={w8:.3f}")
if "g" in S:
    cg, t0 = S["g"]["theta"]
    a.plot(tg, cg * (tg - t0) / tg**2, "--", color="C2", lw=1.6,
           label=rf"(g) $c(\tau-\tau_0)/\tau^2$, c={cg:.2f}, "
                 rf"$\tau_0$={t0:.4f}")
a.axvline(0.5153, color="k", ls=":", lw=1)
a.text(0.518, 0.15, r"doğrudan ölçülen $R{=}0$: $\tau_0=0.5153\pm0.0008$",
       fontsize=7.5, rotation=90, va="bottom")
a.set_xlabel(r"$\tau$"); a.set_ylabel("R")
a.set_title("(2) hayatta kalan 2-parametreli formlar")
a.legend(fontsize=8, loc="upper left"); a.grid(alpha=.3)

# ---- (3) φ_Γ(τ) ve işaret değişimi ----
a = ax[0, 2]
try:
    for i, tb in enumerate(["0.46", "0.52", "0.58"]):
        j, b = yuk(f"R_taban_{tb}.json")
        a.plot([x["tau"] for x in b], [x["phi"] for x in b], "o-", ms=4,
               color=f"C{i}", label=f"taban {tb}")
except FileNotFoundError:
    pass
a.axhline(0, color="k", lw=.8); a.axvline(0.5134, color="k", ls=":", lw=1)
a.set_xlabel(r"$\tau$"); a.set_ylabel(r"$\varphi_\Gamma$")
a.set_title(r"(3) ilkel gözlenebilir $\varphi_\Gamma$: sıfırı $\tau_0$"
            "\nTABANLA KAYMIYOR")
a.legend(fontsize=8); a.grid(alpha=.3)

# ---- (4) taban taraması: R ----
a = ax[1, 0]
try:
    for i, tb in enumerate(["0.46", "0.52", "0.58"]):
        j, b = yuk(f"R_taban_{tb}.json")
        g = iyi(b)
        a.errorbar([x["tau"] for x in g], [x["R"] for x in g],
                   yerr=[max(x["sR_jk"], .02) for x in g], fmt="o-", ms=4,
                   color=f"C{i}", label=f"taban {tb}")
    a.annotate("(τ−taban) motifi\ndoğru olsaydı bu üç eğri\nyatayda 0.06'şar "
               "kayardı", xy=(0.60, 1.4), xytext=(0.63, 0.5), fontsize=8,
               arrowprops=dict(arrowstyle="->", lw=.8))
except FileNotFoundError:
    pass
a.axhline(0, color="k", lw=.8)
a.set_xlabel(r"$\tau$"); a.set_ylabel("R")
a.set_title("(4) TABAN TARAMASI — R, (τ−taban)'ın fonksiyonu DEĞİL")
a.legend(fontsize=8); a.grid(alpha=.3)

# ---- (5) sentetik tanık: R ----
a = ax[1, 1]
a.errorbar(Y["tau"], Y["R_son"], yerr=Y["sig_son"], fmt="ko-", lw=2,
           label="GERÇEK (son 300k)")
for i, (ad, rn) in enumerate([("keskin — soğurma YOK", "keskin"),
                              ("A4 — erfc 0.68/0.125", "A4")]):
    p = OUT / f"R_sentetik_{rn}_std.json"
    if not p.exists():
        continue
    j, b = yuk(p.name); g = iyi(b)
    a.plot([x["tau"] for x in g], [x["R"] for x in g], "s--", ms=6,
           color=f"C{i+2}", label=f"sentetik {ad}")
    bd = [x for x in b if x not in g]
    if bd:
        a.plot([x["tau"] for x in bd], [x["R"] for x in bd], "x", ms=10,
               color=f"C{i+2}", label="  ↑ elenen (faz/patlak)")
a.set_ylim(-0.5, 5.5); a.set_xlabel(r"$\tau$"); a.set_ylabel("R")
a.set_title("(5) sentetik tanık — R şekli 'makine'den")
a.legend(fontsize=8); a.grid(alpha=.3)

# ---- (6) ρ(τ): gerçek düşüyor, sentetik DÜZ ----
a = ax[1, 2]
a.semilogy(ft, Y["ince_rho_son"], "o-", color="C0", label="ρ GERÇEK (son)")
a.semilogy(ft, Y["ince_rho_orta"], "s--", color="C1", mfc="none",
           label="ρ GERÇEK (orta)")
a.semilogy(tg, np.exp(2.62 - 6.20 * tg), ":", color="k",
           label="144: exp(2.62−6.20τ)")
for i, rn in enumerate(["keskin", "A4"]):
    p = OUT / f"R_sentetik_{rn}_std.json"
    if not p.exists():
        continue
    j, b = yuk(p.name)
    if b and b[0].get("rho") is not None:
        a.semilogy([x["tau"] for x in b], [x["rho"] for x in b], "^-",
                   color=f"C{i+2}", label=f"ρ sentetik {rn}")
a.set_xlabel(r"$\tau$"); a.set_ylabel(r"$\rho$")
a.set_title("(6) keskin gazda ρ DÜZ ama R yine tırmanıyor\n"
            r"⇒ R'nin şekli soğurmadan gelmiyor")
a.legend(fontsize=8); a.grid(alpha=.3, which="both")

plt.tight_layout()
plt.savefig(HERE / "154_R_kapali_form.png", dpi=130)
print("-> 154_R_kapali_form.png")
