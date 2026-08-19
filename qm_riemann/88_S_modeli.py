"""
88 — TAM S(ω) MODELİ: CETVELİN TAYFI (21 Ağustos 2026)
==========================================================================
86-87'nin bileşenlerini tek yapı-çarpanı modelinde birleştirir:
  S(ω) = |Ĝ(ω)|²,  Ĝ(ω) = ⟨e^{iωt_n}⟩ (gap ortaları)

  M1  BENEKLER + GEOMETRİ: gap regresyonu dalganın FARKINI ölçer
      (2 sin(πτ) çarpanı), ortalar ORTALAMASINI taşır (cos(πτ)):
        |Ĝ(log p)| = (v_p p^{-1/2}/2) · f(τ),  f(τ) = πτ·cot(πτ).
      86'nın benek-sönümü bilmecesinin (1.07→0.70) aday çözümü —
      parametresiz; f(1/2)=0 → beneklerin katta sönmesi = Nyquist.
  M2  TARAK = KARAKTERİSTİK FONKSİYON: |Ĝ_x(2πk)| tam olarak orta-nokta
      faz dağılımının kar. fonksiyonudur; Gauss yaklaşımı e^{-c·k²}
      0.885 çarpanıyla ıskalıyordu → Edgeworth (ölçülen kurtosis) testi.
  M3  KARANLIK ALAN vs MONTGOMERY RAMPASI: F(α) = n|Ĝ|², α = ω/L.
      Berry resmi: GUE rampası F=α, asal çizgilerin δ-tepelerinden
      İNŞA edilir → bin-ortalama (çizgiler DAHİL) ≈ rampa; çizgi-arası
      taban ≪ rampa. İnce-ızgara taramayla doğrudan test.
  M4  SENTEZ FİGÜRÜ: ölçülen tayf + model bileşenleri tek resimde.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from sympy import primerange

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

d41 = np.load(HERE / "41_bigT_windows.npz")
K41 = sorted({x.split("_")[1] for x in d41.files}, key=lambda s: int(s[:-1]))
WNDS = [(d41[f"gaps_{k}"], d41[f"amps_{k}"], d41[f"tmid_{k}"]) for k in K41[:6]]
d83 = np.load(HERE / "83_tam_taban_egri.npz")

# ============ M1: BENEK GEOMETRİSİ ============
print("M1 — BENEK GEOMETRİSİ: oran = |Ĝ|ölç / (v·p^{-1/2}/2) vs f(τ)=πτ·cot(πτ)")
print(f"{'p':>3} {'⟨τ⟩':>6} {'oran':>6} {'f(τ)':>6} {'oran/f':>7}")
M1 = []
for p in P11:
    rats, fs, taus = [], [], []
    for gaps, amps, tmid in WNDS:
        L = float(np.log(tmid / TWO_PI).mean())
        m = (np.abs(d83["L"] - L) < 0.01) & (np.abs(d83["tau"] - np.log(p)/L) < 1e-9)
        if m.sum() != 1:
            continue
        v = float(d83["v"][m][0])
        tau = np.log(p) / L
        Gm = abs(Ghat(tmid, np.array([np.log(p)]))[0])
        rats.append(Gm / (v * p**-0.5 / 2))
        fs.append(np.pi * tau / np.tan(np.pi * tau))
        taus.append(tau)
    M1.append((np.mean(taus), np.mean(rats), np.mean(fs)))
    print(f"{p:>3} {np.mean(taus):>6.3f} {np.mean(rats):>6.3f} {np.mean(fs):>6.3f} "
          f"{np.mean(rats)/np.mean(fs):>7.3f}")
M1 = np.array(M1)
print(f"  oran/f ortalaması = {np.mean(M1[:,1]/M1[:,2]):.3f} ± "
      f"{np.std(M1[:,1]/M1[:,2]):.3f}  (1'e yakın ve DÜZ ise geometri çözüyor)")

# ============ M2: TARAK = KARAKTERİSTİK FONKSİYON ============
print("\nM2 — TARAK (k=1,2): Gauss vs Edgeworth(ölçülen kurtosis) vs ölçüm")
print(f"{'L':>6} {'k':>2} {'ölçüm':>8} {'Gauss':>8} {'Edgeworth':>9}")
for gaps, amps, tmid in WNDS[:3]:
    L = float(np.log(tmid / TWO_PI).mean())
    x = rvm_N(tmid)
    P1 = np.exp(2j * np.pi * x)
    delta = np.angle(P1 * np.conj(P1.mean() / abs(P1.mean()))) / TWO_PI
    sig = delta.std()
    k4 = float(((delta - delta.mean())**4).mean() / sig**4 - 3)
    for k in [1, 2]:
        meas = abs(np.exp(2j * np.pi * k * x).mean())
        kk = TWO_PI * k
        gauss = np.exp(-kk**2 * sig**2 / 2)
        edge = gauss * (1 + k4 * sig**4 * kk**4 / 24)
        print(f"{L:>6.2f} {k:>2} {meas:>8.4f} {gauss:>8.4f} {edge:>9.4f}"
              f"   (σ_δ={sig:.3f}, kurt={k4:+.2f})")

# ============ M3: MONTGOMERY RAMPASI (Berry ayrışımı) ============
print("\nM3 — F(α) = n|Ĝ|² vs GUE rampası F=α (L=10.37; ince ızgara):")
gaps, amps, tmid = WNDS[1]
L = float(np.log(tmid / TWO_PI).mean())
n = len(tmid)
span = tmid[-1] - tmid[0]
step = TWO_PI / span / 2.2
om_fine = np.arange(0.10, 4.4, step)
G_fine = np.abs(Ghat(tmid, om_fine))**2 * n
lines = [np.log(q) for q in
         [2,3,4,5,7,8,9,11,13,16,17,19,23,25,27,29,31,32,37,41,43,47,49,53,59,61,64,67,71,73,79]]
online = np.zeros(len(om_fine), dtype=bool)
for l in lines:
    online |= np.abs(om_fine - l) < 3 * step
print(f"  ızgara adımı {step:.5f}, {len(om_fine)} nokta; çizgi bandında "
      f"{online.sum()} nokta")
print(f"  {'α-bin':>13} {'F̄ (hepsi)':>10} {'F̄ (çizgisiz)':>12} {'rampa α':>8}")
bins = np.arange(0.10, 4.4, 0.55)
M3 = []
for lo in bins:
    m = (om_fine >= lo) & (om_fine < lo + 0.55)
    if m.sum() < 50:
        continue
    a_mid = (lo + 0.275) / L
    F_all = G_fine[m].mean()
    F_off = G_fine[m & ~online].mean()
    M3.append(((lo + 0.275), F_all, F_off, a_mid))
    print(f"  [{lo:.2f},{lo+0.55:.2f}) {F_all:>10.3f} {F_off:>12.3f} {a_mid:>8.3f}")

# ============ M4: SENTEZ FİGÜRÜ ============
fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.0))
ax = axes[0]
tt = np.linspace(0.02, 0.49, 200)
ax.plot(tt, np.pi * tt / np.tan(np.pi * tt), "k-", lw=1.3,
        label=r"$f(\tau) = \pi\tau\cot(\pi\tau)$ (geometri)")
ax.plot(M1[:, 0], M1[:, 1], "o", c="firebrick", ms=6, label="ölçülen oran")
ax.axvline(0.5, color="gray", ls="--", lw=1)
ax.set_xlabel(r"$\tau$"); ax.set_ylabel("benek / naif öngörü")
ax.set_title("M1: benek-sönümü = ayrık geometri")
ax.legend(fontsize=9); ax.grid(alpha=0.3)
ax = axes[1]
M3a = np.array(M3)
ax.plot(M3a[:, 3], M3a[:, 1], "o-", c="firebrick", label="F̄ bin (çizgiler dahil)")
ax.plot(M3a[:, 3], M3a[:, 2], "s-", c="steelblue", label="F̄ bin (çizgi-arası)")
aa = np.linspace(0, M3a[:, 3].max() * 1.05, 50)
ax.plot(aa, aa, "k--", lw=1.1, label="GUE rampası F = α")
ax.set_xlabel(r"$\alpha = \omega/L$"); ax.set_ylabel(r"$F = n\,|\hat G|^2$")
ax.set_title("M3: rampa, çizgilerden inşa (Berry)")
ax.legend(fontsize=9); ax.grid(alpha=0.3)
ax = axes[2]
ax.semilogy(om_fine, np.sqrt(G_fine / n), ",", c="firebrick", alpha=0.5)
for l in lines:
    if l < 4.4:
        ax.axvline(l, color="steelblue", lw=0.4, alpha=0.4)
ax.axhline(1/np.sqrt(n), color="k", lw=0.8, ls="--", label=r"$1/\sqrt{n}$")
ax.set_xlabel(r"$\omega$"); ax.set_ylabel(r"$|\hat G(\omega)|$")
ax.set_title("M4: cetvelin tayfı (mavi: p^k çizgileri)")
ax.legend(fontsize=9); ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(HERE / "88_S_modeli.png", dpi=110)
print("\nFigür: 88_S_modeli.png")
