"""
103 — NOT 4 FİGÜRÜ v2: KIRINIM (taperli + sıfır ızgarası) (21 Ağustos 2026)
==========================================================================
fig_diffraction_en'i 102a'nın taperli yeniden ölçümüyle günceller.

SOL PANEL  (94'ten değişmedi): aritmetik benekler, mutlak açık-formül
           öngörüsü vs ölçüm.
SAĞ PANEL  (YENİ): çizgi-dışı karanlık alan, 102a_karanlik.npz'den
           |Ĝ| = √(I/n) biriminde:
             • orta noktalar, dikdörtgen pencere  (86'nın ölçümü)
             • orta noktalar, Hann                 (gerçek taban)
             • karıştırılmış-gap vekili, Hann
             • SIFIRLAR z_n, Hann                  (8-10 kademe altta)
             • kusursuz örgü, Hann                 (alet tabanı)
             • 1/√n atım gürültüsü
           İÇ KUTU: örnekleme-fazı eğrisi I(c), t(c) = z + ḡ/2 +
           c(g−ḡ)/2 — c=0 ve c=2 sıfır örgüsü (karanlık), c=1 tepe.

102a_karanlik.npz ızgara indeksleri (GRIDS sırası):
  0 gerçek orta | 1 SIFIRLAR z_n | 2 karıştırılmış-gap |
  3 RvM pürüzsüz | 4 kusursuz örgü | 5 Poisson
  son ek _0 = tapersiz, _1 = Hann.
"""

import numpy as np
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
P11 = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]


def Ghat(t, omegas, chunk=30000):
    out = np.zeros(len(omegas), dtype=complex)
    for s0 in range(0, len(t), chunk):
        tt = t[s0:s0 + chunk]
        out += np.exp(1j * np.outer(omegas, tt)).sum(axis=1)
    return out / len(t)


def rvm_N(t):
    x = t / TWO_PI
    return x * np.log(x / np.e) + 7 / 8


def sigma_t(gaps, tmid):
    t0 = tmid[0] - gaps[0] / 2
    kk = np.arange(len(tmid)) + 0.5
    ts = t0 + kk * TWO_PI / np.log(t0 / TWO_PI)
    for _ in range(6):
        ts = ts - (rvm_N(ts) - rvm_N(t0) - kk) / (np.log(ts / TWO_PI) / TWO_PI)
    return float((tmid - ts).std()), ts


def hann(t):
    return 0.5 * (1 - np.cos(TWO_PI * (t - t[0]) / (t[-1] - t[0])))


def I_omega(t, omegas, taper=True, chunk=256):
    """102a ile birebir aynı kestirimci: I = |Σ w e^{iω(t−c)}|²/Σw²."""
    omegas = np.asarray(omegas, float)
    w = hann(t) if taper else np.ones_like(t)
    nrm = float((w ** 2).sum())
    tc = t - 0.5 * (t[0] + t[-1])
    out = np.empty(len(omegas))
    for s0 in range(0, len(omegas), chunk):
        ob = omegas[s0:s0 + chunk]
        S = (np.exp(1j * np.outer(ob, tc)) * w[None, :]).sum(axis=1)
        out[s0:s0 + chunk] = np.abs(S) ** 2 / nrm
    return out


d41 = np.load(HERE / "41_bigT_windows.npz")
K41 = sorted({x.split("_")[1] for x in d41.files}, key=lambda s: int(s[:-1]))
WNDS = [(d41[f"gaps_{k}"], d41[f"amps_{k}"], d41[f"tmid_{k}"]) for k in K41[:6]]

# ---------------- SOL PANEL: benekler (94 ile birebir) ----------------
meas, pred, taus = [], [], []
for p in P11:
    ms, prs, ts_ = [], [], []
    for gaps, amps, tmid in WNDS:
        L = float(np.log(tmid / TWO_PI).mean())
        tau = np.log(p) / L
        sig, _ = sigma_t(gaps, tmid)
        dw = np.exp(-np.log(p) ** 2 * sig ** 2 / 2)
        ms.append(abs(Ghat(tmid, np.array([np.log(p)]))[0]))
        prs.append(tau * p ** -0.5 * np.cos(np.pi * tau) * dw)
        ts_.append(tau)
    meas.append(np.mean(ms)); pred.append(np.mean(prs)); taus.append(np.mean(ts_))
print("sol panel ok", flush=True)

# ---------------- SAĞ PANEL: 102a'nın taperli taraması ----------------
d = np.load(HERE / "102a_karanlik.npz")
om, n = d["om"], int(d["n"])
amp = lambda key: np.sqrt(d[key] / n)          # I → |Ĝ|

# ---------------- İÇ KUTU: örnekleme-fazı eğrisi ----------------
gaps, amps, tmid = WNDS[1]                     # 102a'nın penceresi
z = np.empty(len(gaps) + 1)
z[0] = tmid[0] - gaps[0] / 2
z[1:] = z[0] + np.cumsum(gaps)
gbar = gaps.mean()
cs = np.array([0.0, 0.25, 0.5, 1.0, 1.5, 2.0])
Ic = []
for c in cs:
    t_c = z[:-1] + gbar / 2 + c * (gaps - gbar) / 2
    Ic.append(np.median(I_omega(np.sort(t_c), om, taper=True)))
    print(f"  c={c:.2f}  I={Ic[-1]:.3e}", flush=True)
Ic = np.array(Ic)

# ---------------- figür ----------------
fig, axes = plt.subplots(1, 3, figsize=(17.4, 5.0),
                         gridspec_kw={"width_ratios": [1.05, 1.35, 0.72]})

ax = axes[0]
ax.plot(taus, pred, "s-", c="steelblue", ms=6,
        label=r"explicit formula: $\tau p^{-1/2}\cos(\pi\tau)\cdot\mathrm{DW}$")
ax.plot(taus, meas, "o", c="firebrick", ms=6,
        label=r"measured spot $|\hat G(\log p)|$")
ax.set_xlabel(r"$\tau = \log p / L$"); ax.set_ylabel(r"$|\hat G|$")
ax.set_title("Arithmetic spots: parameter-free prediction")
ax.legend(fontsize=9); ax.grid(alpha=0.3)

ax = axes[1]
ax.semilogy(om, amp("I_0_0"), ".", ms=2, c="darkorange", alpha=0.45,
            label="midpoints, rectangular window")
ax.semilogy(om, amp("I_2_1"), ".", ms=2, c="steelblue", alpha=0.45,
            label="shuffled-gap surrogate (Hann)")
ax.semilogy(om, amp("I_0_1"), ".", ms=2, c="firebrick", alpha=0.55,
            label="midpoints, Hann")
ax.semilogy(om, amp("I_1_1"), ".", ms=2, c="seagreen", alpha=0.55,
            label=r"zeros $z_n$, Hann")
ax.semilogy(om, amp("I_4_1"), ".", ms=2, c="0.6", alpha=0.45,
            label="perfect lattice, Hann (instrument)")
ax.axhline(1 / np.sqrt(n), color="k", lw=0.9, ls="--", label=r"shot noise $1/\sqrt{n}$")
ax.set_xlabel(r"$\omega$ (off-line)"); ax.set_ylabel(r"$|\hat G(\omega)|$")
ax.set_ylim(1e-13, 3e-2)
ax.set_title("Dark field: the midpoint grid vs the zeros")
ax.legend(fontsize=7.3, loc="lower center", ncol=3, markerscale=3.5,
          framealpha=0.92)
ax.grid(alpha=0.3)

ax = axes[2]
ax.semilogy(cs, Ic, "o-", c="firebrick", ms=6, lw=1.4)
ax.set_xticks([0, 0.5, 1, 1.5, 2])
ax.set_xlabel(r"sampling phase $c$")
ax.set_ylabel(r"median $I(\omega)$ off the lines")
ax.set_title("The dark field is a sampling phase", fontsize=11)
ax.grid(alpha=0.3)
ax.annotate("zero lattice", xy=(0.0, Ic[0]), xytext=(0.12, 3e-12),
            fontsize=8, color="seagreen")
ax.annotate("zero lattice", xy=(2.0, Ic[-1]), xytext=(1.05, 3e-12),
            fontsize=8, color="seagreen", ha="left")
ax.annotate("midpoints", xy=(1.0, Ic[3]), xytext=(1.0, 4e-3),
            fontsize=8, color="firebrick", ha="center")
ax.set_ylim(3e-15, 6e-2)

fig.tight_layout()
fig.savefig(HERE / "fig_diffraction_en.png", dpi=110)
print("fig_diffraction_en ok (v2)")
