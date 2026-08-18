"""
79 — B'NİN EĞİMİ: ARKA ODA ÇADIRIN YAMACINI İZLİYOR MU? (19 Ağustos 2026)
==========================================================================
Fısıltı (Not 3 çizim seansı): B = −0.16 platosu ayna-kanalı sızıntısıysa
düz olmamalı — çadırı izlemeli: B(τ) ∝ −2(1−τ) (kat ötesi).

GECENİN SONUÇLARI (bu script tam kaydıdır):
  T1  Havuzlanmış kat-ötesi eğim NEGATİF çıkar — ama bu SIMPSON TUZAĞI:
      yüksek-τ noktaları yalnız küçük-L pencerelerinden gelir.
  T1b L-kontrollü eğim hafif POZİTİF (+0.18±0.11): çadır yönünde ama
      düzle de uyumlu — biçim çözülmedi. Termal-çadır (DW×tent) fiti
      anti-sönüm ister (s<0) → basit termal-çadır ÖLÜ.
  T2  Tam-eğri yarışı (132 nokta): sabit-B ve çadır-B eşdeğer (kolinerite);
      ikisi birlikte fit edilirse dejenere patlar (M3).
  T3  PLASEBO DERSİ: p+0.5 sahteleri asal çizgisine 3-16 çözünürlük
      mesafesinde — SIZINTI (+0.07 yalancı taban). Geometrik-orta sahtelerle
      taban ±0.02-0.05'e iner; sinyal bunun 7-17 katı.
  T4  AYNI-BAND TARAMASI [0.505,0.55]: kısıtlı tabanla B̂ = −0.31…−0.36,
      L=5.6→16.6 boyunca SABİT (36+41+55 veri setleri). Daha önceki
      "eriyen B(L)" pencere-ortalamasının τ-karışım artefaktıydı.
  T5  TABAN BUKALEMUNU (2×2): aynı bandda küçük asallar kovaryat olarak
      eklenince B̂ = −0.12'ye iner — o da L-SABİT. Çıplak −0.35 vs
      giydirilmiş −0.12 (×2.8); Not 3'ün −0.16'sı kendi konvansiyonunda
      tutarlı, arada. Giydirme mekanizması (71'in etkileşim kanalı?) AÇIK.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.optimize import least_squares
from sympy import primerange

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
A_CG = (np.e**2 - 5) / 2
B0, B1 = 2.7580, -0.0543
P11 = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
P_EXT = [37, 41, 43, 47, 53, 59]
FAKES = [2.31, 6.7, 23.7, 34.7, 44.3, 50.9, 57.3]
TAU_CAP = 0.78

def channels(gaps, amps, tmid, primes, fakes):
    Lw = np.log(tmid / TWO_PI)
    g_u = gaps * Lw / TWO_PI
    a_u = amps / np.sqrt(A_CG * Lw + B0 + B1 / Lw)
    a_u /= np.sqrt((a_u**2).mean())
    ya = np.log(a_u)
    cw = [np.ones_like(ya), g_u, g_u**2]
    qs = list(primes) + list(fakes)
    for q in qs:
        arg = tmid * np.log(q)
        cw += [np.cos(arg), np.sin(arg)]
    Xw = np.vstack(cw).T
    bw, *_ = np.linalg.lstsq(Xw, ya, rcond=None)
    sew = np.sqrt((ya - Xw @ bw).var() * np.diag(np.linalg.inv(Xw.T @ Xw)))
    L = float(Lw.mean())
    out = []
    for i, q in enumerate(qs):
        out.append((np.log(q) / L, bw[3 + 2*i] / q**-0.5,
                    sew[3 + 2*i] / q**-0.5, q not in primes, q, L))
    return out

def wfit(X, y, s):
    Xw_, yw = X / s[:, None], y / s
    c, *_ = np.linalg.lstsq(Xw_, yw, rcond=None)
    cov = np.linalg.inv(Xw_.T @ Xw_)
    return c, np.sqrt(np.diag(cov)), float(((yw - Xw_ @ c) ** 2).sum())

def bhat(rows, keep):
    pw = np.array([(w, s) for t, w, s, f, q, _ in rows if q in keep])
    return (np.sum(pw[:, 0] / pw[:, 1]**2) / np.sum(1 / pw[:, 1]**2),
            1 / np.sqrt(np.sum(1 / pw[:, 1]**2)))

WNDS = []
d36 = np.load(HERE / "36_T100k.npz")
edges = np.geomspace(d36["t_mid"][0], d36["t_mid"][-1] * 1.0001, 13)
for i in range(12):
    m = (d36["t_mid"] >= edges[i]) & (d36["t_mid"] < edges[i + 1])
    if m.sum() >= 500:
        WNDS.append((d36["intervals"][m], d36["max_amps"][m], d36["t_mid"][m]))
d41 = np.load(HERE / "41_bigT_windows.npz")
for k in sorted({x.split("_")[1] for x in d41.files}, key=lambda s: int(s[:-1])):
    WNDS.append((d41[f"gaps_{k}"], d41[f"amps_{k}"], d41[f"tmid_{k}"]))

# ---- taban A (Not 3: 11 asal) ve taban B (uzatılmış + uzak-olmayan sahteler)
ptsA, ptsB = [], []
for gaps, amps, tmid in WNDS:
    L = float(np.log(tmid / TWO_PI).mean())
    ptsA += channels(gaps, amps, tmid, P11, [])
    ext = [p for p in P_EXT if np.log(p) / L <= TAU_CAP]
    fks = [q for q in FAKES if np.log(q) / L <= TAU_CAP]
    ptsB += channels(gaps, amps, tmid, P11 + ext, fks)
A_ = np.array([(t, w, s) for t, w, s, f, q, L in ptsA])
BpL = np.array([(t, w, s, L) for t, w, s, f, q, L in ptsB if not f])
Bf = np.array([(t, w, s) for t, w, s, f, q, L in ptsB if f])
print(f"Taban A: {len(A_)} nokta (τ≤{A_[:,0].max():.2f}); "
      f"taban B: {len(BpL)} asal + {len(Bf)} sahte (τ≤{BpL[:,0].max():.2f})")

# ---- T1: havuzlanmış kat-ötesi eğim (SIMPSON TUZAĞI — uyarıyla)
print("\nT1 — KAT-ÖTESİ HAVUZLANMIŞ (τ>0.5; DİKKAT: yüksek τ yalnız küçük L'den):")
for isim, P in [("Taban A", A_[:, :3]), ("Taban B", BpL[:, :3])]:
    m = P[:, 0] > 0.5
    t, w, s = P[m].T
    n = int(m.sum())
    cF, eF, x2F = wfit(np.ones((n, 1)), w, s)
    cL, eL, x2L = wfit(np.vstack([np.ones(n), t]).T, w, s)
    print(f"  {isim} (n={n}): düz B={cF[0]:+.4f}±{eF[0]:.4f} (χ²/dof {x2F/(n-1):.2f}) | "
          f"eğim {cL[1]:+.3f}±{eL[1]:.3f} (çadır öngörüsü ≈ +{abs(cF[0])/(1-t.mean()):.2f})")

# ---- T1b: L-kontrollü ve pencere-içi eğimler
m = BpL[:, 0] > 0.5
t, w, s, L = BpL[m].T
n = int(m.sum())
cM, eM, _ = wfit(np.vstack([np.ones(n), t, L - L.mean()]).T, w, s)
print(f"\nT1b — L-KONTROLLÜ eğim = {cM[1]:+.3f} ± {eM[1]:.3f} "
      f"(dL katsayısı {cM[2]:+.4f}±{eM[2]:.4f})")
for Lv in np.unique(L):
    mm = L == Lv
    if mm.sum() >= 4:
        c_, e_, _ = wfit(np.vstack([np.ones(int(mm.sum())), t[mm]]).T, w[mm], s[mm])
        print(f"  pencere-içi L={Lv:5.2f}: eğim = {c_[1]:+.3f} ± {e_[1]:.3f}")
r3 = least_squares(lambda p: (-p[0]*2*(1-t)*np.exp(-p[1]*t**2) - w) / s, [0.4, 2.0])
print(f"  termal-çadır fiti: k={r3.x[0]:.3f}, sönüm s={r3.x[1]:+.3f} "
      f"(s<0 = anti-sönüm → fiziksel değil, model ölü)")

# ---- T3: plasebo dersi
mF = Bf[:, 0] > 0.45
print(f"\nT3 — PLASEBO: uzak-sahte |w| ort = {np.abs(Bf[mF,1]).mean():.4f} "
      f"(sinyal −0.35'in ~%5'i). p+0.5 sahteleri SIZDIRIR (çözünürlük ~3e-4).")

# ---- T4: aynı-band taraması (kısıtlı taban) + T5: 2×2 taban etkisi
print("\nT4/T5 — AYNI BAND [0.505,0.55], iki taban:")
print(f"{'L':>6} {'kısıtlı-B̂':>12} {'giydirilmiş-B̂':>14}")
scanL, scanB1, scanE1, scanB2, scanE2 = [], [], [], [], []
for gaps, amps, tmid in WNDS:
    Lv = float(np.log(tmid / TWO_PI).mean())
    ps = [p for p in primerange(int(np.exp(0.505 * Lv)), int(np.exp(0.55 * Lv)))][:10]
    if len(ps) < 2:
        continue
    small = [p for p in P11 if np.log(p) / Lv < 0.45]
    b1, e1 = bhat(channels(gaps, amps, tmid, ps, []), set(ps))
    b2, e2 = bhat(channels(gaps, amps, tmid, small + ps, []), set(ps))
    print(f"{Lv:>6.2f} {b1:>+9.4f}±{e1:.3f} {b2:>+10.4f}±{e2:.3f}")
    scanL.append(Lv); scanB1.append(b1); scanE1.append(e1)
    scanB2.append(b2); scanE2.append(e2)
for f in ["55_win_1e+08.npz"]:
    d = np.load(HERE / f)
    gaps, amps, tmid = d["gaps"], d["amps"], d["tmid"]
    Lv = float(np.log(tmid / TWO_PI).mean())
    ps = [p for p in primerange(int(np.exp(0.505 * Lv)), int(np.exp(0.55 * Lv)))][:10]
    b1, e1 = bhat(channels(gaps, amps, tmid, ps, []), set(ps))
    print(f"{Lv:>6.2f} {b1:>+9.4f}±{e1:.3f}        (55, aynı-motor kontrol)")
    scanL.append(Lv); scanB1.append(b1); scanE1.append(e1)
scanL = np.array(scanL); scanB1 = np.array(scanB1); scanE1 = np.array(scanE1)
scanB2 = np.array(scanB2); scanE2 = np.array(scanE2)
c1 = np.sum(scanB1/scanE1**2)/np.sum(1/scanE1**2)
c2 = np.sum(scanB2/scanE2**2)/np.sum(1/scanE2**2)
print(f"\nkombine: çıplak B = {c1:+.4f}, giydirilmiş B = {c2:+.4f} "
      f"(oran {c1/c2:.2f}) — İKİSİ DE L-SABİT")

# ---- grafik
fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.2))
ax = axes[0]
mB = BpL[:, 0] > 0.42
ax.errorbar(BpL[mB, 0], BpL[mB, 1], yerr=BpL[mB, 2], fmt="o", ms=4,
            c="firebrick", alpha=0.7, label="asal (taban B)")
ax.plot(Bf[mF, 0], Bf[mF, 1], "x", ms=6, c="gray", label="uzak-sahte taban")
tt = np.linspace(0.42, 0.80, 100)
ax.plot(tt, np.full_like(tt, -0.181), "k-", lw=1.2, label="düz B (havuz)")
ax.plot(tt, -0.201 * 2 * (1 - tt), "b--", lw=1.2, label="çadır −k·2(1−τ)")
ax.axvline(0.5, color="gray", ls="--", lw=1)
ax.axhline(0, color="gray", lw=0.6)
ax.set_xlabel("τ"); ax.set_ylabel("w")
ax.set_title("Kat-ötesi: biçim çözülmedi (Simpson uyarısı metinde)")
ax.legend(fontsize=9); ax.grid(alpha=0.3)
ax = axes[1]
ax.errorbar(scanL[:len(scanB2)], scanB2, yerr=scanE2, fmt="s", ms=5, c="steelblue",
            label=f"giydirilmiş (küçük asallar kovaryat): {c2:+.3f}")
ax.errorbar(scanL, scanB1, yerr=scanE1, fmt="o", ms=5, c="firebrick",
            label=f"çıplak (yalnız band asalları): {c1:+.3f}")
ax.axhline(c1, color="firebrick", lw=0.8, ls=":")
ax.axhline(c2, color="steelblue", lw=0.8, ls=":")
ax.set_xlabel("L"); ax.set_ylabel("B̂ (band [0.505, 0.55])")
ax.set_title("Arka oda: iki okuma, ikisi de L-sabit")
ax.legend(fontsize=9); ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(HERE / "79_B_egimi.png", dpi=110)
print("Grafik: 79_B_egimi.png")
