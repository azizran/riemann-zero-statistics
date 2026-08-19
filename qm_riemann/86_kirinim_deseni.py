"""
86 — SIFIR ÖRGÜSÜNÜN KIRINIM DESENİ: CETVEL TEORİSİ İLK ÖLÇÜM (20 Ağustos)
==========================================================================
Fikir: "cetvel bükülmesi"nin tamamı tek nesnede yaşar — örneklem
ızgarasının kendi Fourier'i (yapı çarpanı):
    Ĝ(ω) = (1/n) Σ_n exp(i ω t_n),   t_n = gap ortaları.
İki kolonun Gram-çiftlenimi Ĝ'nin fark/toplam frekansındaki değeridir;
81-83'ün bütün "bükülme" fenomenolojisi bu desenin gölgesiydi.

BİRİNCİ-İLKE TEORİSİ (u = yerdeğiştirme alanı, t_n = t^s_n + u_n):
  Ĝ(ω) ≈ [pürüzsüz ≈ 0] + iω·⟨e^{iωt^s} u⟩ →
  P1  ASAL BENEKLERİ: |Ĝ(log p)| = v_p · p^{-1/2} / 2
      (v_p 83'ün tam-taban ölçümü; U_p = v_p p^{-1/2}/log p, çarpan iω)
      — PARAMETRESİZ öngörü.
  P2  FAZ KİLİDİ: v kanalı cos-kilitli (faz 0) ⟹ u sin-kilitli ⟹
      iω çarpanıyla Ĝ(log p) REEL-kilitli olmalı (Im ≈ 0), işaret tutarlı.
  P3  TARAK (katlanmamış çerçeve x = N_RvM(t)): birinci mertebe benek
      |Ĝ_x(2π)| ≈ e^{-c_jit(L)}  (Debye-Waller; (2πσ_u)²/2 = c_jit)
      — sıcaklığın ÜÇÜNCÜ bağımsız ölçümü.
  T4  GENİŞBANT: çizgi-dışı |Ĝ(ω)|² tabanı = u-alanının güç tayfı
      (S-Alanı Denizi'nin ilk haritası); kontrol: karıştırılmış-gap
      vekil ızgara (aritmetik silinir, taban kalır) + RvM pürüzsüz
      (her şey söner).
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from sympy import primerange

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
P11 = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]

def Ghat(t, omegas, chunk=40000):
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

# ---- P1 + P2: asal benekleri, pencere pencere
print("P1/P2 — ASAL BENEKLERİ (öngörü |Ĝ| = v_p·p^{-1/2}/2; faz reel-kilit):")
print(f"{'p':>3} {'⟨|Ĝ|ölç⟩':>9} {'⟨öngörü⟩':>9} {'oran':>6} {'⟨Re⟩':>8} {'⟨Im⟩':>8} {'taban':>7}")
oranlar = []
for p in P11:
    Gm, Gpred, Re, Im, flo = [], [], [], [], []
    for gaps, amps, tmid in WNDS:
        L = float(np.log(tmid / TWO_PI).mean())
        m = (np.abs(d83["L"] - L) < 0.01) & (np.abs(d83["tau"] - np.log(p)/L) < 1e-9)
        if m.sum() != 1:
            continue
        v = float(d83["v"][m][0])
        om = np.log(p)
        G = Ghat(tmid, np.array([om, om * 1.0173]))   # ikincisi çizgi-dışı taban
        Gm.append(abs(G[0])); Gpred.append(v * p**-0.5 / 2)
        Re.append(G[0].real); Im.append(G[0].imag); flo.append(abs(G[1]))
    Gm, Gpred = np.mean(Gm), np.mean(Gpred)
    oranlar.append(Gm / Gpred)
    print(f"{p:>3} {Gm:>9.4f} {Gpred:>9.4f} {Gm/Gpred:>6.2f} "
          f"{np.mean(Re):>+8.4f} {np.mean(Im):>+8.4f} {np.mean(flo):>7.4f}")
print(f"  oran ortalaması = {np.mean(oranlar):.3f} ± {np.std(oranlar):.3f}")

# ---- P3: tarak (katlanmamış çerçeve), pencere pencere sıcaklık
print("\nP3 — TARAK BENEĞİ (öngörü |Ĝ_x(2π)| ≈ e^{-c_jit}):")
print(f"{'L':>6} {'|Ĝ_x(2π)|':>10} {'e^-c_jit':>9} {'oran':>6}")
for gaps, amps, tmid in WNDS:
    L = float(np.log(tmid / TWO_PI).mean())
    x = rvm_N(tmid)
    t0 = tmid[0] - gaps[0] / 2
    kk = np.arange(len(tmid)) + 0.5
    ts = t0 + kk * TWO_PI / np.log(t0 / TWO_PI)
    for _ in range(6):
        fdel = rvm_N(ts) - rvm_N(t0) - kk
        ts = ts - fdel / (np.log(ts / TWO_PI) / TWO_PI)
    c_j = (L * (tmid - ts).std()) ** 2 / 2
    om_grid = TWO_PI + np.linspace(-0.02, 0.02, 41)
    Gx = Ghat(x, om_grid)
    pk = np.abs(Gx).max()
    print(f"{L:>6.2f} {pk:>10.4f} {np.exp(-c_j):>9.4f} {pk/np.exp(-c_j):>6.2f}")

# ---- T4: genişbant taban — üç ızgara
print("\nT4 — GENİŞBANT (L=10.37; çizgi-dışı ω taraması):")
gaps, amps, tmid = WNDS[1]
L = float(np.log(tmid / TWO_PI).mean())
rng = np.random.default_rng(86)
g_sh = rng.permutation(gaps)
t_sh = tmid[0] + np.cumsum(g_sh) - g_sh / 2
t0 = tmid[0] - gaps[0] / 2
kk = np.arange(len(tmid)) + 0.5
ts = t0 + kk * TWO_PI / np.log(t0 / TWO_PI)
for _ in range(6):
    fdel = rvm_N(ts) - rvm_N(t0) - kk
    ts = ts - fdel / (np.log(ts / TWO_PI) / TWO_PI)
lines = set()
for q in range(2, 60):
    fs = [f for f in [q] if all(q % pp or q == pp for pp in range(2, q))]
lines = [np.log(q) for q in [2,3,4,5,7,8,9,11,13,16,17,19,23,25,27,29,31,32,37,41,43,47,49,53,59]]
om_scan = np.linspace(0.15, 4.2, 1200)
om_scan = np.array([o for o in om_scan if min(abs(o - l) for l in lines) > 0.015])
G_real = np.abs(Ghat(tmid, om_scan))
G_shuf = np.abs(Ghat(t_sh, om_scan))
G_rvm = np.abs(Ghat(ts, om_scan))
n = len(tmid)
print(f"  n = {n}; beyaz-gürültü tabanı 1/√n = {1/np.sqrt(n):.4f}")
for isim, G in [("gerçek ızgara", G_real), ("karıştırılmış-gap", G_shuf),
                ("RvM pürüzsüz", G_rvm)]:
    print(f"  {isim:>18}: medyan |Ĝ| = {np.median(G):.4f}  "
          f"ω<1'de {np.median(G[om_scan < 1]):.4f}, ω>3'te {np.median(G[om_scan > 3]):.4f}")

# ---- figür
fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.2))
ax = axes[0]
ps = np.array(P11)
meas = []
pred = []
for p in P11:
    Gm, Gp = [], []
    for gaps_, amps_, tmid_ in WNDS:
        Lw = float(np.log(tmid_ / TWO_PI).mean())
        m = (np.abs(d83["L"] - Lw) < 0.01) & (np.abs(d83["tau"] - np.log(p)/Lw) < 1e-9)
        if m.sum() != 1:
            continue
        Gm.append(abs(Ghat(tmid_, np.array([np.log(p)]))[0]))
        Gp.append(float(d83["v"][m][0]) * p**-0.5 / 2)
    meas.append(np.mean(Gm)); pred.append(np.mean(Gp))
ax.plot(np.log(ps), pred, "s-", c="steelblue", ms=6, label="öngörü: v·p^{-1/2}/2 (83'ün v'si)")
ax.plot(np.log(ps), meas, "o", c="firebrick", ms=6, label="ölçülen |Ĝ(log p)|")
ax.set_xlabel("ω = log p"); ax.set_ylabel("|Ĝ|")
ax.set_title("Asal benekleri: parametresiz öngörü vs ölçüm")
ax.legend(fontsize=9); ax.grid(alpha=0.3)
ax = axes[1]
ax.semilogy(om_scan, G_real, ".", ms=2, c="firebrick", alpha=0.6, label="gerçek ızgara")
ax.semilogy(om_scan, G_shuf, ".", ms=2, c="steelblue", alpha=0.5, label="karıştırılmış-gap")
ax.semilogy(om_scan, G_rvm, ".", ms=2, c="gray", alpha=0.5, label="RvM pürüzsüz")
ax.axhline(1/np.sqrt(n), color="k", lw=0.8, ls="--", label="1/√n")
ax.set_xlabel("ω (çizgi-dışı)"); ax.set_ylabel("|Ĝ(ω)|")
ax.set_title("Genişbant ışıltı: üç ızgara")
ax.legend(fontsize=9); ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(HERE / "86_kirinim.png", dpi=110)
print("\nFigür: 86_kirinim.png")
