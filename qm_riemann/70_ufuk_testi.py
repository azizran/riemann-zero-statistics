"""
70 — RIEMANN-SIEGEL UFKU TESTİ: PERDE NEDEN TAM KAPANMIYOR? (18 Ağustos)
==========================================================================

Hipotez: genlik (RS ana toplamı) yalnız p ≤ √(t/2π) asallarını içerir →
τ = log p / L = 1/2'de SERT UFUK. Sıfır konumları (explicit formula) için
ufuk yok. Öngörüler:
  P1: w(τ) ufka yaklaşırken çöker (~0), yasanın v-kesiminden ÖNCE
  P2: v(τ) τ=1/2'yi pürüzsüz geçer
  P3: yasa (√w + 1.07v = 1) ufuk bölgesinde bükülür — kapanmama çözülür

Ölçüm: 36'nın 6 + 41'in 6 penceresi × p ∈ {2..31} (11 asal), tmid fazları.
Uyarı payı: büyük τ'da tepe-konumu yayılımı w'yi ~e^{-0.8τ²} kadar
söndürebilir (τ=0.55'te ~0.79 çarpanı) — çöküş bununla KARIŞMAZ (o %20,
bu ~%100).
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
A = (np.e**2 - 5) / 2
B0, B1 = 2.7580, -0.0543
PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
FAKES = [2.31, 6.7]

def channels(gaps, amps, tmid):
    Lw = np.log(tmid / TWO_PI)
    g_u = gaps * Lw / TWO_PI
    a_u = amps / np.sqrt(A * Lw + B0 + B1 / Lw)
    a_u /= np.sqrt((a_u**2).mean())
    ya, yg = np.log(a_u), np.log(g_u)
    cw = [np.ones_like(ya), g_u, g_u**2]
    cv = [np.ones_like(yg)]
    for q in PRIMES + FAKES:
        arg = tmid * np.log(q)
        cw += [np.cos(arg), np.sin(arg)]
        cv += [np.cos(arg), np.sin(arg)]
    Xw = np.vstack(cw).T
    bw, *_ = np.linalg.lstsq(Xw, ya, rcond=None)
    sew = np.sqrt((ya - Xw @ bw).var() * np.diag(np.linalg.inv(Xw.T @ Xw)))
    Xv = np.vstack(cv).T
    bv, *_ = np.linalg.lstsq(Xv, yg, rcond=None)
    sev = np.sqrt((yg - Xv @ bv).var() * np.diag(np.linalg.inv(Xv.T @ Xv)))
    L = float(Lw.mean())
    out = []
    for i, p in enumerate(PRIMES):
        w_p = bw[3 + 2 * i] / p**-0.5          # işaretli (cos) bileşen
        sw_p = sew[3 + 2 * i] / p**-0.5
        v_p = np.hypot(bv[1 + 2 * i], bv[2 + 2 * i]) / p**-0.5
        sv_p = sev[1 + 2 * i] / p**-0.5
        out.append((np.log(p) / L, p, L, w_p, sw_p, v_p, sv_p))
    return out

rows = []
d36 = np.load(HERE / "36_T100k.npz")
edges = np.geomspace(d36["t_mid"][0], d36["t_mid"][-1] * 1.0001, 13)
for i in range(12):
    m = (d36["t_mid"] >= edges[i]) & (d36["t_mid"] < edges[i + 1])
    if m.sum() < 500:
        continue
    rows += channels(d36["intervals"][m], d36["max_amps"][m], d36["t_mid"][m])
d41 = np.load(HERE / "41_bigT_windows.npz")
for k in sorted({x.split("_")[1] for x in d41.files}, key=lambda s: int(s[:-1])):
    rows += channels(d41[f"gaps_{k}"], d41[f"amps_{k}"], d41[f"tmid_{k}"])

rows.sort()
print(f"{len(rows)} (τ, w, v) üçlüsü — τ aralığı: {rows[0][0]:.3f} … {rows[-1][0]:.3f}\n")
print("UFUK BÖLGESİ (τ ≥ 0.40):")
print(f"{'τ':>6} {'p':>3} {'L':>6} {'w':>8} {'±':>6} {'v':>7} {'yasa √w-öngörüsü²':>17}")
for tau, p, L, w_p, sw_p, v_p, sv_p in rows:
    if tau >= 0.40:
        law = max(1 - 1.0683 * v_p, 0) ** 2
        print(f"{tau:>6.3f} {p:>3} {L:>6.2f} {w_p:>8.4f} {sw_p:>6.4f} "
              f"{v_p:>7.4f} {law:>17.4f}")

# bin özetleri
arr = np.array([(t, w_, sw_, v_, sv_) for t, p, L, w_, sw_, v_, sv_ in rows])
print("\nBİN ÖZETLERİ:")
print(f"{'τ-bin':>13} {'n':>3} {'w_ort':>8} {'v_ort':>7}")
for lo, hi in [(0.30,0.40),(0.40,0.45),(0.45,0.50),(0.50,0.55),(0.55,0.62)]:
    m = (arr[:, 0] >= lo) & (arr[:, 0] < hi)
    if m.sum():
        print(f"[{lo:.2f},{hi:.2f}) {m.sum():>3} {arr[m,1].mean():>8.4f} "
              f"{arr[m,3].mean():>7.4f}")

# grafik
fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.4))
ax = axes[0]
ax.errorbar(arr[:, 0], arr[:, 1], yerr=arr[:, 2], fmt="o", ms=3.5,
            c="firebrick", alpha=0.7, label="w (işaretli cos bileşeni)")
ax.errorbar(arr[:, 0], arr[:, 3], yerr=arr[:, 4], fmt="s", ms=3.5,
            c="teal", alpha=0.7, label="v")
ax.axvline(0.5, color="k", ls="--", lw=1.2, label="RS ufku τ = 1/2")
ax.axvline(0.40, color="gray", ls=":", lw=1, label="τ* ≈ 0.40")
ax.axhline(0, color="gray", lw=0.6)
ax.set_xlabel("τ = log p / L"); ax.set_ylabel("kanal katsayısı")
ax.set_title("İki kanal ufka doğru"); ax.legend(fontsize=9); ax.grid(alpha=0.3)
ax = axes[1]
m_ok = arr[:, 1] > 0
sq = np.sqrt(arr[m_ok, 1])
col = np.where(arr[m_ok, 0] < 0.40, "steelblue", "firebrick")
ax.scatter(arr[m_ok, 3], sq, s=16, c=col)
xx = np.linspace(0, 0.7, 50)
ax.plot(xx, np.clip(1 - 1.0683 * xx, 0, None), "k--", lw=1.1,
        label="yasa: √w = 1 − 1.07·v")
ax.scatter([], [], c="steelblue", label="τ < 0.40")
ax.scatter([], [], c="firebrick", label="τ ≥ 0.40 (ufuk bölgesi)")
ax.set_xlabel("v"); ax.set_ylabel("√w")
ax.set_title("Yasa ufuk bölgesinde bükülüyor mu?")
ax.legend(fontsize=9); ax.grid(alpha=0.3)
plt.tight_layout()
out = HERE / "70_ufuk_testi.png"
plt.savefig(out, dpi=110)
print(f"\nGrafik: {out.name}")
