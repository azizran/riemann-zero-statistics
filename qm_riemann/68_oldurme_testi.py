"""
68 — ADAY YASANIN ÖLDÜRME TESTİ: p=11,13 OUT-OF-SAMPLE (18 Ağustos 2026)
==========================================================================

Aday yasa (67): √w + b·v = 1, b = 1.0683 (a≡1 sabitli fit; yalnız
p∈{2,3,5,7}, τ≤0.197 çiftlerinden).

Test: p=11 ve 13 için w ŞİMDİ ölçülüyor (yasa türetiminde yok) —
τ=0.20-0.26 bölgesi de menzil dışı. Öngörü: √w = 1 − 1.0683·v.

Ölçüm: 11 pencere (41'in 6'sı + 55'in 4'ü + Odlyzko), 6-asal regresyon
(g̃, g̃² kontrollü, tmid fazları — 56 konvansiyonu).
"""

import numpy as np
import mpmath as mp
import matplotlib.pyplot as plt
from pathlib import Path

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
A = (np.e**2 - 5) / 2
B0, B1 = 2.7580, -0.0543
T0_ODL = 267653395647
PRIMES6 = [2, 3, 5, 7, 11, 13]
B_LAW = 1.0683

mp.mp.dps = 30

def measure_wv(gaps, amps, tmid, anchored=False):
    Lw = np.log(((T0_ODL + tmid) if anchored else tmid) / TWO_PI)
    g_u = gaps * Lw / TWO_PI
    a_u = amps / np.sqrt(A * Lw + B0 + B1 / Lw)
    a_u /= np.sqrt((a_u**2).mean())
    ya, yg = np.log(a_u), np.log(g_u)
    cols_w = [np.ones_like(ya), g_u, g_u**2]
    cols_v = [np.ones_like(yg)]
    for p in PRIMES6:
        ph0 = float(mp.fmod(T0_ODL * mp.log(p), 2 * mp.pi)) if anchored else 0.0
        arg = ph0 + tmid * np.log(p)
        cols_w += [np.cos(arg), np.sin(arg)]
        cols_v += [np.cos(arg), np.sin(arg)]
    Xw = np.vstack(cols_w).T
    bw, *_ = np.linalg.lstsq(Xw, ya, rcond=None)
    sew = np.sqrt((ya - Xw @ bw).var() * np.diag(np.linalg.inv(Xw.T @ Xw)))
    Xv = np.vstack(cols_v).T
    bv, *_ = np.linalg.lstsq(Xv, yg, rcond=None)
    out = []
    L = float(Lw.mean())
    for i, p in enumerate(PRIMES6):
        w_p = bw[3 + 2 * i] / p**-0.5
        se_p = sew[3 + 2 * i] / p**-0.5
        v_p = np.hypot(bv[1 + 2 * i], bv[2 + 2 * i]) / p**-0.5
        out.append((p, np.log(p) / L, w_p, se_p, v_p))
    return L, out

WINDOWS = []
d41 = np.load(HERE / "41_bigT_windows.npz")
for k in sorted({x.split("_")[1] for x in d41.files}, key=lambda s: int(s[:-1])):
    WINDOWS.append((d41[f"gaps_{k}"], d41[f"amps_{k}"], d41[f"tmid_{k}"], False))
for f in ["55_win_1e+08.npz", "55_win_1e+09.npz", "55_win_1e+10.npz", "55_win_1e+11.npz"]:
    d = np.load(HERE / f)
    WINDOWS.append((d["gaps"], d["amps"], d["tmid"], False))
d53 = np.load(HERE / "53_odlyzko_amps.npz")
WINDOWS.append((d53["gaps"], d53["max_amps"], d53["t_mid"], True))

new_pairs, old_pairs = [], []
for gaps, amps, tmid, anch in WINDOWS:
    L, rows = measure_wv(gaps, amps, tmid, anchored=anch)
    for p, tau, w_p, se_p, v_p in rows:
        rec = (tau, w_p, se_p, v_p, p, L)
        (new_pairs if p in (11, 13) else old_pairs).append(rec)

print("OUT-OF-SAMPLE: p=11,13 çiftleri (yasa bunları hiç görmedi)")
print(f"{'p':>3} {'L':>6} {'τ':>6} {'v':>7} {'√w ölç':>8} {'√w öngörü':>10} {'fark':>7} {'σ':>5}")
devs = []
for tau, w_p, se_p, v_p, p, L in sorted(new_pairs):
    if w_p <= 0:
        print(f"{p:>3} {L:>6.2f} {tau:>6.3f}  w<0 ölçüldü ({w_p:.3f}) — atlandı")
        continue
    sq_m = np.sqrt(w_p)
    sq_pred = 1 - B_LAW * v_p
    se_sq = se_p / (2 * sq_m)
    dev = sq_m - sq_pred
    devs.append((dev, se_sq))
    print(f"{p:>3} {L:>6.2f} {tau:>6.3f} {v_p:>7.3f} {sq_m:>8.4f} {sq_pred:>10.4f} "
          f"{dev:>+7.4f} {abs(dev)/se_sq:>5.1f}")

dv = np.array([d for d, s in devs])
print(f"\nYENİ çiftler: ortalama fark {dv.mean():+.4f}, RMS {np.sqrt((dv**2).mean()):.4f}")
print(f"(Eski 44 çiftin taban RMS'i 0.0091 idi — kıyas ölçütü)")

# grafik: tüm çiftler tek doğru üstünde mi?
fig, ax = plt.subplots(figsize=(8.5, 5.6))
op = np.array([(v_p, np.sqrt(w_p)) for tau, w_p, se_p, v_p, p, L in old_pairs if w_p > 0])
npr = np.array([(v_p, np.sqrt(w_p)) for tau, w_p, se_p, v_p, p, L in new_pairs if w_p > 0])
ax.plot(op[:, 0], op[:, 1], "o", ms=4, c="steelblue", alpha=0.6,
        label="p ∈ {2,3,5,7} (yasanın türetildiği 44)")
ax.plot(npr[:, 0], npr[:, 1], "s", ms=7, c="firebrick", zorder=5,
        label="p ∈ {11,13} (OUT-OF-SAMPLE)")
xx = np.linspace(0, 0.65, 50)
ax.plot(xx, 1 - B_LAW * xx, "k--", lw=1.2, label=f"√w = 1 − {B_LAW}·v")
ax.set_xlabel("v — konum kanalı"); ax.set_ylabel("√w — genlik-iletim genliği")
ax.set_title("Genlik korunumu adayı: öldürme testi")
ax.legend(fontsize=9); ax.grid(alpha=0.3)
plt.tight_layout()
out = HERE / "68_oldurme_testi.png"
plt.savefig(out, dpi=110)
print(f"\nGrafik: {out.name}")
