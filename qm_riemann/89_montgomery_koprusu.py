"""
89 — MONTGOMERY KÖPRÜSÜ: v-TOPLAM KURALI TÜRETİMİ VE TESTİ (21 Ağustos)
==========================================================================
TÜRETİM (ilk kez kapalı form):
  S(t) çizgisi: Λ(q)/(π√q log q)  ⟹  u-çizgisi U_q = 2Λ(q)/(L√q log q)
  Ortalar ortalama alır (cos πτ), gap'ler fark alır (sin πτ):
    BENEK:  |Ĝ(log q)| = Λ(q)/(L√q) · cos(πτ) · DW      [v'siz, mutlak]
    GAP:    v-katsayısı = U_q·(L/π)·sin(πτ)
  Asallar (Λ=log p):   ★ v(τ) = (2/π)·sin(πτ) ★
    → onset v ≈ 2τ (ölçülen 2.014 ✓ türedi), doyum 2/π = 0.6366
    → kuvvetler: v(p^k) = (2/π)sin(πτ)/k  (YENİ öngörü: 1/k bastırma)
  Pythagoras: benek ∝ cos²(πτ), gap ∝ sin²(πτ) — Montgomery rampası
  (F=α teoremi) iki kanala sin²+cos²=1 ile bölüşülür.

TESTLER:
  T1  v-yasası: 132-eğrinin v'leri vs (2/π)sin(πτ); serbest genlik A_v.
  T2  MUTLAK benek: |Ĝ(log p)| vs τ·p^{-1/2}·cos(πτ)·e^{-ω²σ_t²/2}
      (hiçbir kanal ölçümü girmiyor — explicit formülden doğrudan).
  T3  ENJEKSİYON KALİBRASYONU: gerçek sıfır dizisine bilinen sentetik
      dalga ekle (aritmetik-dışı ω*), v-regresyonu ve beneği aynı boru
      hattıyla ölç → a_v(τ), a_G(τ) zayıflama çarpanları → k₀'ın ve
      v_ölç/v_gerçek ~0.86-0.92 açığının ölçüm-geometrisi payı.
  T4  RAMPA BÖLÜŞÜMÜ: bin-F̄ (çizgiler) vs α·cos²(πα)·DW·k₀².
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

def sigma_t(gaps, tmid):
    t0 = tmid[0] - gaps[0] / 2
    kk = np.arange(len(tmid)) + 0.5
    ts = t0 + kk * TWO_PI / np.log(t0 / TWO_PI)
    for _ in range(6):
        fdel = rvm_N(ts) - rvm_N(t0) - kk
        ts = ts - fdel / (np.log(ts / TWO_PI) / TWO_PI)
    return float((tmid - ts).std())

# ============ T1: v-YASASI ============
print("T1 — v-TOPLAM KURALI: v(τ) vs (2/π)·sin(πτ)  [tüm 132-eğri v'leri]")
tauv, vv = d83["tau"], d83["v"]
pred = (2 / np.pi) * np.sin(np.pi * tauv)
rat = vv / pred
A_v = np.sum(vv * np.sin(np.pi * tauv)) / np.sum(np.sin(np.pi * tauv)**2)
print(f"  n = {len(vv)} nokta; ⟨v/öngörü⟩ = {rat.mean():.3f} ± {rat.std():.3f}")
print(f"  serbest genlik fiti A_v·sin(πτ): A_v = {A_v:.4f}  "
      f"(2/π = {2/np.pi:.4f}; oran {A_v/(2/np.pi):.3f})")
for lo, hi in [(0, 0.1), (0.1, 0.2), (0.2, 0.3), (0.3, 0.45), (0.45, 0.62)]:
    m = (tauv >= lo) & (tauv < hi)
    if m.sum():
        print(f"    τ∈[{lo:.2f},{hi:.2f}): ⟨v/öngörü⟩ = {rat[m].mean():.3f} (n={m.sum()})")

# ============ T2: MUTLAK BENEK ============
print("\nT2 — MUTLAK BENEK: |Ĝ| vs τ·p^{-1/2}·cos(πτ)·DW  [kanal ölçümü YOK]")
print(f"{'p':>3} {'ölçüm':>8} {'öngörü':>8} {'oran':>6}")
sigs = {id(w[2]): sigma_t(w[0], w[2]) for w in WNDS}
T2 = []
for p in P11:
    ms, prs = [], []
    for gaps, amps, tmid in WNDS:
        L = float(np.log(tmid / TWO_PI).mean())
        tau = np.log(p) / L
        om = np.log(p)
        dw = np.exp(-om**2 * sigs[id(tmid)]**2 / 2)
        ms.append(abs(Ghat(tmid, np.array([om]))[0]))
        prs.append(tau * p**-0.5 * np.cos(np.pi * tau) * dw)
    T2.append((np.log(p) / 10.4, np.mean(ms), np.mean(prs)))
    print(f"{p:>3} {np.mean(ms):>8.4f} {np.mean(prs):>8.4f} {np.mean(ms)/np.mean(prs):>6.3f}")
T2 = np.array(T2)
print(f"  oran ortalaması = {(T2[:,1]/T2[:,2]).mean():.3f} ± {(T2[:,1]/T2[:,2]).std():.3f}")

# ============ T3: ENJEKSİYON KALİBRASYONU ============
print("\nT3 — ENJEKSİYON (gerçek diziye bilinen dalga; L=10.37):")
gaps, amps, tmid = WNDS[1]
L = float(np.log(tmid / TWO_PI).mean())
z = np.empty(len(gaps) + 1)
z[0] = tmid[0] - gaps[0] / 2
z[1:] = z[0] + np.cumsum(gaps)
gmean = gaps.mean()
print(f"  {'τ*':>6} {'a_v = v_ölç/v_naif':>18} {'sin-öngörü':>11} "
      f"{'a_G = Ĝ_ölç/Ĝ_naif':>18} {'cos-öngörü':>11}")
for taustar in [0.10, 0.25, 0.40]:
    om = taustar * L
    U = 0.02
    zp = z + U * np.cos(om * z)
    gp = np.diff(zp)
    mp_ = 0.5 * (zp[:-1] + zp[1:])
    Lw = np.log(mp_ / TWO_PI)
    yg = np.log(gp * Lw / TWO_PI)
    X = np.vstack([np.ones_like(mp_), np.cos(om * mp_), np.sin(om * mp_)]).T
    b, *_ = np.linalg.lstsq(X, yg, rcond=None)
    Lw0 = np.log(tmid / TWO_PI)
    yg0 = np.log(gaps * Lw0 / TWO_PI)
    X0 = np.vstack([np.ones_like(tmid), np.cos(om * tmid), np.sin(om * tmid)]).T
    b0, *_ = np.linalg.lstsq(X0, yg0, rcond=None)
    coef = np.hypot(b[1] - b0[1], b[2] - b0[2])
    a_v = coef / (U * om)
    sin_pred = 2 * np.sin(om * gmean / 2) / (om * gmean)
    dG = abs(Ghat(mp_, np.array([om]))[0] - Ghat(tmid, np.array([om]))[0])
    a_G = dG / (om * U / 2)
    cos_pred = np.cos(om * gmean / 2)
    print(f"  {taustar:>6.2f} {a_v:>18.3f} {sin_pred:>11.3f} "
          f"{a_G:>18.3f} {cos_pred:>11.3f}")

# ============ T4: RAMPA BÖLÜŞÜMÜ ============
print("\nT4 — RAMPA: bin-F̄(çizgiler) vs α·cos²(πα)·DW  (L=10.37)")
n = len(tmid)
span = tmid[-1] - tmid[0]
step = TWO_PI / span / 2.2
om_fine = np.arange(0.10, 4.4, step)
G_fine = np.abs(Ghat(tmid, om_fine))**2 * n
sig = sigs[id(tmid)]
print(f"  {'α':>6} {'F̄ ölçüm':>9} {'model':>7}")
for lo in np.arange(0.65, 4.4, 0.55):
    m = (om_fine >= lo) & (om_fine < lo + 0.55)
    a_mid = (lo + 0.275) / L
    om_mid = lo + 0.275
    model = a_mid * np.cos(np.pi * a_mid)**2 * np.exp(-om_mid**2 * sig**2)
    print(f"  {a_mid:>6.3f} {G_fine[m].mean():>9.3f} {model:>7.3f}")

# ============ FİGÜR ============
fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.2))
ax = axes[0]
ax.plot(tauv, vv, "o", ms=4, c="steelblue", alpha=0.6, label="ölçülen v (132 nokta)")
tt = np.linspace(0, 0.62, 200)
ax.plot(tt, (2/np.pi) * np.sin(np.pi * tt), "k-", lw=1.4,
        label=r"türetim: $v = \frac{2}{\pi}\sin(\pi\tau)$")
ax.plot(tt, A_v * np.sin(np.pi * tt), "k--", lw=1,
        label=f"serbest genlik: {A_v:.3f}·sin(πτ)")
ax.plot(tt, 2 * tt, ":", c="gray", lw=1, label="onset 2τ (türedi)")
ax.axhline(2/np.pi, color="firebrick", lw=0.8, ls=":", label="doyum 2/π")
ax.set_xlabel(r"$\tau$"); ax.set_ylabel("v")
ax.set_title("Montgomery köprüsü: v-toplam kuralı")
ax.legend(fontsize=9); ax.grid(alpha=0.3)
ax = axes[1]
ax.plot(T2[:, 0], T2[:, 2], "s-", c="steelblue", ms=6,
        label=r"explicit formül: $\tau p^{-1/2}\cos(\pi\tau)\cdot$DW")
ax.plot(T2[:, 0], T2[:, 1], "o", c="firebrick", ms=6, label="ölçülen benek")
ax.set_xlabel(r"$\tau$ (≈, ort. pencere)"); ax.set_ylabel(r"$|\hat G(\log p)|$")
ax.set_title("Mutlak benek: kanal ölçümsüz öngörü")
ax.legend(fontsize=9); ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(HERE / "89_kopru.png", dpi=110)
print("\nFigür: 89_kopru.png")
