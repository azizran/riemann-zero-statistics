"""
74 — KAT KAPISI: (1−2τ) İKİ-TOPLAM GİRİŞİMİNDEN TÜRETİLDİ + TEST
==================================================================
(18 Ağustos 2026)

TÜRETİM (3 satır): Z = 2Σ_{n≤N} n^{-1/2} cos(θ − t log n), N=√(t/2π).
|Z|²'de cos(t log p) frekansını yalnız (m, mp) çiftleri üretir (mp ≤ N):
  katsayı = 4 p^{-1/2} Σ_{m≤N/p} 1/m ;  DC = 2 Σ_{n≤N} 1/n
  → bağıl genlik R_pred(p) = 2 p^{-1/2} · S₁(⌊N/p⌋)/S₁(N)
  → sürekli limitte p^{-1/2}·2(1−2τ);  p>N'de TAM SIFIR (keskin ufuk).

Bu bir ÖZDEŞLİK — jitter/DW dokunamaz. Test: noktasal |Z|² modülasyonu
(1.27M grid noktası, L=9.86 penceresi) parametresiz tahminle kıyas.
Beklenti: (i) fold altında üst üste; (ii) τ=1/2'de keskin kesiliş;
(iii) fazlar saf-cos. Ardından max-tabanlı w ile fark = kristal fiziği.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from sympy import primerange

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi

def theta(t):
    return t / 2 * np.log(t / TWO_PI) - t / 2 - np.pi / 8 + 1 / (48 * t)

def Z_rs(t, chunk=200000):
    t = np.asarray(t, dtype=np.float64)
    out = np.empty_like(t)
    for s in range(0, len(t), chunk):
        tt = t[s:s + chunk]
        a = np.sqrt(tt / TWO_PI)
        N = a.astype(np.int64)
        th = theta(tt)
        z = np.zeros_like(tt)
        for Nv in np.unique(N):
            m = N == Nv
            n = np.arange(1, Nv + 1)
            ph = th[m, None] - tt[m, None] * np.log(n)[None, :]
            z[m] = 2 * (np.cos(ph) @ (n**-0.5))
        p = a - N
        cp = np.cos(TWO_PI * p)
        cp = np.where(np.abs(cp) < 1e-8, 1e-8, cp)
        psi = np.cos(TWO_PI * (p**2 - p - 1 / 16)) / cp
        z += (-1) ** (N - 1) * (tt / TWO_PI) ** -0.25 * psi
        out[s:s + chunk] = z
    return out

# ---- pencere ve grid ----
T_LO, T_HI, STEP = 107252.0, 132748.0, 0.02
grid = np.arange(T_LO, T_HI, STEP)
print(f"Z gridi: {len(grid)} nokta hesaplanıyor...")
Z = Z_rs(grid)
Y = Z**2
Y = Y / Y.mean()
L = float(np.log(0.5 * (T_LO + T_HI) / TWO_PI))
N_RS = int(np.sqrt(0.5 * (T_LO + T_HI) / TWO_PI))
print(f"L = {L:.3f}, N_RS = {N_RS} (fold: p = {N_RS})")

# ---- test asalları (fold'un iki yakası) ----
PS = [2, 3, 5, 7, 11, 13, 17, 23, 31, 43, 59, 79, 101, 113, 127, 131, 137,
      139, 149, 163, 181, 199, 251, 307, 397, 499]

cols = [np.ones_like(Y)]
for p in PS:
    arg = grid * np.log(p)
    cols += [np.cos(arg), np.sin(arg)]
X = np.vstack(cols).T
b, *_ = np.linalg.lstsq(X, Y, rcond=None)
se = np.sqrt((Y - X @ b).var() * np.diag(np.linalg.inv(X.T @ X)))

S1 = lambda M: np.sum(1.0 / np.arange(1, max(int(M), 1) + 1))
S1N = S1(N_RS)

print(f"\n{'p':>4} {'τ':>6} | {'ölçülen R':>10} {'±':>7} | {'TAHMİN':>8} | {'sin':>7}")
rows = []
for i, p in enumerate(PS):
    tau = np.log(p) / L
    meas = b[1 + 2 * i] / p**-0.5
    err = se[1 + 2 * i] / p**-0.5
    sinp = b[2 + 2 * i] / p**-0.5
    pred = 2 * S1(N_RS / p) / S1N if p <= N_RS else 0.0
    rows.append((tau, meas, err, pred))
    tag = " ◄fold" if 120 < p < 160 else ""
    print(f"{p:>4} {tau:>6.3f} | {meas:>10.4f} {err:>7.4f} | {pred:>8.4f} | "
          f"{sinp:>7.4f}{tag}")

rows = np.array(rows)
m_lo = rows[:, 0] < 0.47
resid = rows[m_lo, 1] - rows[m_lo, 3]
print(f"\nFold-altı uyum: artık RMS = {np.sqrt(np.mean(resid**2)):.4f} "
      f"(tahmin ölçeği ~0.3-1.9, PARAMETRESİZ)")
m_hi = rows[:, 0] > 0.52
print(f"Fold-ötesi: ölçülen |R| ort = {np.abs(rows[m_hi, 1]).mean():.4f} "
      f"(tahmin: TAM 0)")

# grafik
fig, ax = plt.subplots(figsize=(9, 5.6))
tt = np.linspace(0.03, 0.68, 400)
pp = np.exp(tt * L)
pr = np.array([2 * S1(N_RS / p) / S1N if p <= N_RS else 0.0 for p in pp])
ax.plot(tt, pr, "k-", lw=1.4, label="TÜRETİM: 2·S₁(N/p)/S₁(N) → 2(1−2τ)")
ax.errorbar(rows[:, 0], rows[:, 1], yerr=rows[:, 2], fmt="o", ms=5,
            c="firebrick", label="ölçülen |Z|² modülasyonu")
ax.axvline(0.5, color="gray", ls="--", lw=1.1, label="kat çizgisi τ=1/2")
ax.axhline(0, color="gray", lw=0.6)
ax.set_xlabel("τ = log p / L"); ax.set_ylabel("bağıl modülasyon / p^(−1/2)")
ax.set_title("Kat kapısı: parametresiz aritmetik tahmin vs ölçüm (L=9.86)")
ax.legend(fontsize=9); ax.grid(alpha=0.3)
plt.tight_layout()
out = HERE / "74_kat_turetimi.png"
plt.savefig(out, dpi=110)
print(f"\nGrafik: {out.name}")
