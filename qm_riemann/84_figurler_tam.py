"""
84 — REVİZYON FİGÜRLERİ + ODLYZKO BATARYASI (TAM TABAN) (20 Ağustos 2026)
==========================================================================
83'ün tam-taban kampanyasından Not 2-3 revizyon figürleri:
  fig_law_en      — korunum/kısıt yasası (yeni doğru 1.017 − 0.884v)
  fig_crossing_en — 132-nokta tam-taban eğrisi + termal model (kesiş 0.443)
  fig_tau_en      — Not 2: w(τ) çökmesi, mavi=12 birincil pencere,
                    kırmızı=Odlyzko (L=24.475), kesikli=yeni kılavuz
Ayrıca: Odlyzko tam-taban bataryası (p=2..13 w/v) yazdırılır ve yasa
çiftleri 84_yasa_ciftleri.npz'e kaydedilir. fig_tent_en DEĞİŞMEZ
(noktasal yoğun-ızgara ölçümü, örnekleme sistematiğinden muaf).
"""

import numpy as np
import matplotlib.pyplot as plt

exec(open("83_buyuk_yeniden_olcum.py").read().split("# ---- pencereler")[0])

LAW = []
d41 = np.load(HERE / "41_bigT_windows.npz")
K41 = sorted({x.split("_")[1] for x in d41.files}, key=lambda s: int(s[:-1]))
for k in K41[:6]:
    LAW.append((d41[f"gaps_{k}"], d41[f"amps_{k}"], d41[f"tmid_{k}"], False))
for f in ["55_win_1e+08.npz", "55_win_1e+09.npz", "55_win_1e+10.npz", "55_win_1e+11.npz"]:
    d = np.load(HERE / f)
    LAW.append((d["gaps"], d["amps"], d["tmid"], False))
d53 = np.load(HERE / "53_odlyzko_amps.npz")
LAW.append((d53["gaps"], d53["max_amps"], d53["t_mid"], True))

tr, oos, odl = [], [], []
for wi, (gaps, amps, tmid, anch) in enumerate(LAW):
    out, guc, L = kanal_tam(gaps, amps, tmid, P4 + P2, anch, powers_too=False)
    for p in P4 + P2:
        tau, w_, s_, wsin, v_ = out[p]
        if anch:
            odl.append((p, tau, w_, s_, v_))
        if w_ > 0:
            (tr if p in P4 else oos).append((v_, np.sqrt(w_), tau, w_))
tr = np.array(tr); oos = np.array(oos)
np.savez(HERE / "84_yasa_ciftleri.npz", tr=tr, oos=oos, odl=np.array(odl))

print("ODLYZKO TAM-TABAN BATARYASI (L=24.475):")
d83 = np.load(HERE / "83_tam_taban_egri.npz")
tauE, wE, swE = d83["tau"], d83["w"], d83["sw"]
mg = tauE < 0.30
cg = np.polyfit(tauE[mg], wE[mg], 1, w=1/swE[mg])
print(f"  yeni kılavuz (12 pencere, τ<0.30): w = {cg[1]:.3f} {cg[0]:+.3f}·τ")
for p, tau, w_, s_, v_ in odl:
    print(f"  p={int(p):>2} τ={tau:.4f}: w = {w_:.4f} ± {s_:.4f} "
          f"(kılavuz {cg[1]+cg[0]*tau:.4f}, fark {w_-(cg[1]+cg[0]*tau):+.4f})  v = {v_:.4f}")

# ---- fig_law_en
fig, ax = plt.subplots(figsize=(8.2, 5.4))
ax.plot(tr[:, 0], tr[:, 1], "o", ms=4, c="steelblue", alpha=0.6,
        label=r"$p \in \{2,3,5,7\}$ (44 pairs, complete-basis)")
ax.plot(oos[:, 0], oos[:, 1], "s", ms=6, c="firebrick", zorder=5,
        label=r"$p \in \{11,13\}$ (out of sample)")
xx = np.linspace(0, max(tr[:, 0].max(), oos[:, 0].max()) * 1.05, 50)
ax.plot(xx, 1.0171 - 0.8835 * xx, "k-", lw=1.2,
        label=r"$\sqrt{w} = 1.017 - 0.884\,v$")
ax.plot(xx, 1 - 1.0683 * xx, "k--", lw=0.9, alpha=0.45,
        label="previous convention: $1 - 1.068\,v$")
ax.set_xlabel(r"$v$ — zero-displacement channel")
ax.set_ylabel(r"$\sqrt{w}$ — amplitude transmission")
ax.legend(fontsize=9); ax.grid(alpha=0.3)
plt.tight_layout(); plt.savefig(HERE / "fig_law_en.png", dpi=110); plt.close()
print("fig_law_en ok")

# ---- fig_crossing_en
c_jit, B_new, A_new = 1.17, -0.0924, 1.020
fig, ax = plt.subplots(figsize=(8.6, 5.4))
ax.errorbar(tauE, wE, yerr=swE, fmt="o", ms=3, c="firebrick", alpha=0.55,
            label=r"$w(\tau)$ (132 measurements, complete-basis)")
tt = np.linspace(0.001, 0.62, 500)
ax.plot(tt, A_new*(1-2*tt)*np.exp(-c_jit*tt**2)*(tt < 0.5) + B_new, "k-", lw=1.4,
        label=rf"thermal Bragg: $A(1-2\tau)e^{{-c\tau^2}} + B$ ($A$={A_new:.2f})")
ax.axvline(0.5, color="gray", ls="--", lw=1, label=r"ideal mirror $\tau = 1/2$")
ax.axvline(0.447, color="steelblue", ls=":", lw=1.2,
           label=r"measured crossing $0.447 \pm 0.005$")
ax.axhline(0, color="gray", lw=0.6)
ax.set_xlabel(r"$\tau = \log p / L$"); ax.set_ylabel(r"$w$")
ax.legend(fontsize=9); ax.grid(alpha=0.3)
plt.tight_layout(); plt.savefig(HERE / "fig_crossing_en.png", dpi=110); plt.close()
print("fig_crossing_en ok")

# ---- fig_tau_en (Not 2)
fig, ax = plt.subplots(figsize=(8.2, 5.4))
ax.errorbar(tauE, wE, yerr=swE, fmt="o", ms=3.5, c="steelblue", alpha=0.6,
            label=r"$t \leq 1.6\times10^6$ (complete-basis)")
oarr = np.array(odl)
ax.errorbar(oarr[:, 1], oarr[:, 2], yerr=oarr[:, 3], fmt="s", ms=6,
            c="firebrick", zorder=5, label=r"$t = 2.677\times10^{11}$ (out of sample)")
tt = np.linspace(0, 0.35, 50)
ax.plot(tt, cg[1] + cg[0]*tt, "k--", lw=1.1,
        label=rf"guide fitted to blue only: ${cg[1]:.3f} {cg[0]:+.3f}\,\tau$")
ax.set_xlabel(r"$\tau = \log p / L$")
ax.set_ylabel(r"$w$ — amplitude transmission")
ax.set_xlim(0, 0.65)
ax.legend(fontsize=9); ax.grid(alpha=0.3)
plt.tight_layout(); plt.savefig(HERE / "fig_tau_en.png", dpi=110); plt.close()
print("fig_tau_en ok")
