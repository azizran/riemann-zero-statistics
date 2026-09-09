# -*- coding: utf-8 -*-
"""184e — ZARFIN ANATOMİSİ figürü (K1 profil + K2 yarışı + K3 köprüsü)."""
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.special import erfc

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S184 = SCR / "184"
S180 = SCR / "180"


def load(ad):
    p = S184 / f"K1_{ad}.npz"
    return np.load(p) if p.exists() else None


def band_prof(K2, tw):
    sat = [s for s in K2["bant"][tw] if not s["etiket"].startswith("KUYRUK")]
    return (np.array([s["tau"] for s in sat]),
            np.array([s["w_g"] for s in sat]),
            np.array([s["w_tw"] for s in sat]),
            np.array([s["r"] for s in sat]),
            np.array([s["se_r"] for s in sat]))


def main():
    K2 = json.load(open(S184 / "K2_sonuc.json"))
    tau, wg, wHk, r, se = band_prof(K2, "Hkeskin")
    yar = K2["fit"]["yaris"]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.2))

    # --- panel 1: w(τ) profilleri ---
    ax1.axhline(1.0, color="0.6", lw=1, ls=":", label="analitik nominal (w=1)")
    ax1.plot(tau, wg, "o-", color="tab:red", ms=6, label="gerçek (son)")
    ax1.plot(tau, wHk, "s-", color="tab:blue", ms=5,
             label="sadakatli-keskin ikiz (Hkeskin)")
    for ad, c, mk, lb in (("HA4", "tab:green", "^", "erfc-ikiz (HA4, 0.68)"),
                          ("ikiz", "0.5", "x", "keskin (152 kaba)")):
        d = load(ad)
        if d is not None:
            t = d["tau"]; wq = d["wq"]; aeff = d["aq_eff"]
            edg = [0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.86]
            xs, ys = [], []
            for i in range(len(edg) - 1):
                m = (t > edg[i]) & (t <= edg[i + 1])
                if m.sum():
                    xs.append(np.average(t[m], weights=aeff[m]))
                    ys.append(d["ahat"][m].sum() / aeff[m].sum())
            ax1.plot(xs, ys, mk + "-", color=c, ms=5, label=lb, alpha=0.8)
    ax1.set_xlabel(r"$\tau=\log q/L$")
    ax1.set_ylabel(r"$w(\tau)=\hat a_q/(2a_q\sin\pi\tau)$")
    ax1.set_title("K1 — çizgi-genlik profilleri (etkin genlik / nominal)")
    ax1.legend(fontsize=8, loc="upper left")
    ax1.grid(alpha=0.3)

    # --- panel 2: saf zarf r(τ) + hipotezler ---
    ax2.errorbar(tau, r, yerr=se, fmt="o", color="k", ms=6, capsize=3,
                 label=r"saf zarf $r=\hat a_{gerçek}/\hat a_{Hkeskin}$", zorder=5)
    tt = np.linspace(0.44, 0.87, 300)
    p2 = yar["H-Z2_guc"]["par"]
    ax2.plot(tt, np.clip(1 - p2[0] * tt ** p2[1], 0, None), "-",
             color="tab:green", lw=2,
             label=rf"H-Z2 güç $1-{p2[0]:.3f}\tau^{{{p2[1]:.2f}}}$ "
                   rf"(χ²/dof={yar['H-Z2_guc']['chi2_dof']:.2f}) KAZANIR")
    p1 = yar["H-Z1_erfc"]["par"]
    ax2.plot(tt, 0.5 * erfc((tt - p1[0]) / p1[1]), "--", color="tab:orange",
             lw=1.6, label=rf"H-Z1 erfc $\tau_c={p1[0]:.2f},\Delta={p1[1]:.2f}$ "
                          rf"(χ²/dof={yar['H-Z1_erfc']['chi2_dof']:.2f})")
    p3 = yar["H-Z3_gaussDW"]["par"]
    ax2.plot(tt, np.exp(-(2 * np.pi * tt) ** 2 * p3[0] ** 2 / 2), ":",
             color="tab:purple", lw=1.6,
             label=rf"H-Z3 Gauss-DW $\sigma={p3[0]:.3f}$ "
                   rf"(χ²/dof={yar['H-Z3_gaussDW']['chi2_dof']:.1f}) ÖLDÜ")
    if (S180 / "r_tau_profil.npy").exists():
        tg, rr = np.load(S180 / "r_tau_profil.npy")
        ax2.plot(tg, rr, "d", color="tab:blue", ms=5, mfc="none",
                 label="180 R_bant zarfı (özdeş)")
    ax2.axvline(0.70, color="0.7", lw=0.8, ls=":")
    ax2.set_xlabel(r"$\tau=\log q/L$")
    ax2.set_ylabel(r"$r(\tau)=\hat a_{gerçek}/\hat a_{Hkeskin}$")
    ax2.set_title("K2/K3 — saf zarf: yumuşak güç-yasası (kesim DEĞİL)")
    ax2.legend(fontsize=8, loc="lower left")
    ax2.grid(alpha=0.3)

    fig.suptitle("184 — ZARFIN ANATOMİSİ: gerçeği ikizinden ayıran ~%80'in "
                 "kimliği = doğrudan-ölçülen çizgi-genlik zarfı", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out = QM / "184_zarf.png"
    fig.savefig(out, dpi=130)
    print(f"-> {out}")


if __name__ == "__main__":
    main()
