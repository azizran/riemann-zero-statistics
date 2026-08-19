"""
90 — KUVVETLERİN v'Sİ: 1/k TESTİ (21 Ağustos 2026)
==========================================================================
89'un türetimi: u-çizgisi U_q = 2Λ(q)/(L√q log q) ⟹ kuvvetler için
  v(p^k) = (2/π)sin(πτ)·D-perde/k  ⟹  k·v(p^k) = asal v-eğrisi(τ_q).
Karşı hipotez (log q ağırlığı, bastırma yok): v(p^k) = asal eğrisi → R≈k.

T1  KIRINIM BENEĞİ (regresyonsuz, mutlak): |Ĝ(log q)| vs
    (τ/k)·q^{-1/2}·cos(πτ)·DW  — 1/k'li ve 1/k'siz öngörü yarışı.
T2  GAP KANALI: tam-taban gap regresyonundan v_q; oran testi
    R = k·v_q / v_asal(τ_q)  [aynı pencere, aynı regresyon, interp] —
    κ ve D biçiminden bağımsız. Plasebo: tamsayı-olmayan frekanslar.
"""

import numpy as np
from pathlib import Path

exec(open("83_buyuk_yeniden_olcum.py").read().split("# ---- pencereler")[0])

QPOW = [(4, 2, 2), (8, 2, 3), (9, 3, 2), (16, 2, 4),
        (25, 5, 2), (27, 3, 3), (32, 2, 5), (49, 7, 2)]
FAKES = [4.53, 8.61, 26.4]

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
        fdel = rvm_N(ts) - rvm_N(t0) - kk
        ts = ts - fdel / (np.log(ts / TWO_PI) / TWO_PI)
    return float((tmid - ts).std())

d41 = np.load(HERE / "41_bigT_windows.npz")
K41 = sorted({x.split("_")[1] for x in d41.files}, key=lambda s: int(s[:-1]))
W41 = [(d41[f"gaps_{k}"], d41[f"amps_{k}"], d41[f"tmid_{k}"]) for k in K41[:6]]

d36 = np.load(HERE / "36_T100k.npz")
edges = np.geomspace(d36["t_mid"][0], d36["t_mid"][-1] * 1.0001, 13)
EXT = []
for i in range(12):
    m = (d36["t_mid"] >= edges[i]) & (d36["t_mid"] < edges[i + 1])
    if m.sum() >= 500:
        EXT.append((d36["intervals"][m], d36["max_amps"][m], d36["t_mid"][m]))
EXT += W41

# ============ T1: KIRINIM BENEKLERİ ============
print("T1 — KUVVET BENEKLERİ (6 pencere; öngörü (τ/k)·q^{-1/2}·cos·DW):")
print(f"{'q':>3} {'k':>2} {'ölçüm':>8} {'1/k öngörü':>10} {'oran':>6} {'k-siz oran':>10}")
for q, pp, kk_ in QPOW:
    ms, prs = [], []
    for gaps, amps, tmid in W41:
        L = float(np.log(tmid / TWO_PI).mean())
        tau = np.log(q) / L
        if tau > 0.55:
            continue
        om = np.log(q)
        dw = np.exp(-om**2 * sigma_t(gaps, tmid)**2 / 2)
        ms.append(abs(Ghat(tmid, np.array([om]))[0]))
        prs.append((tau / kk_) * q**-0.5 * np.cos(np.pi * tau) * dw)
    if ms:
        r = np.mean(ms) / np.mean(prs)
        print(f"{q:>3} {kk_:>2} {np.mean(ms):>8.4f} {np.mean(prs):>10.4f} "
              f"{r:>6.3f} {r/kk_:>10.3f}")

# ============ T2: GAP KANALI ============
print("\nT2 — GAP KANALI (12+6 pencere; R = k·v_q / v_asal(τ_q)):")
acc = {q: [] for q, _, _ in QPOW}
fake_acc = []
for gaps, amps, tmid in EXT:
    ya, yg, g_u, L = unfold(gaps, amps, tmid)
    targets = [q for q, _, _ in QPOW if np.log(q) / L <= 0.55]
    qs = sorted(set(P11) | set(targets) | set(tam_taban(L, 0.55)))
    qs_all = qs + FAKES
    bv, sev, _ = chunked_reg(yg, tmid, qs_all, False, g_u=None)
    prim = []
    for p in P11:
        if np.log(p) / L > 0.60:
            continue
        i = qs_all.index(p)
        prim.append((np.log(p) / L,
                     np.hypot(bv[1 + 2*i], bv[2 + 2*i]) / p**-0.5))
    prim = np.array(sorted(prim))
    for q, pp, kk_ in QPOW:
        tau = np.log(q) / L
        if q not in targets or tau < prim[0, 0] or tau > prim[-1, 0]:
            continue
        i = qs_all.index(q)
        vq = np.hypot(bv[1 + 2*i], bv[2 + 2*i]) / q**-0.5
        se = sev[1 + 2*i] / q**-0.5
        vpr = np.interp(tau, prim[:, 0], prim[:, 1])
        acc[q].append((kk_ * vq / vpr, kk_ * se / vpr, vq, se))
    for f in FAKES:
        i = qs_all.index(f)
        fake_acc.append(np.hypot(bv[1 + 2*i], bv[2 + 2*i]) / f**-0.5)

print(f"{'q':>3} {'k':>2} {'n':>3} {'R = k·v/v_asal':>14} {'k-siz R/k':>10} {'⟨v_q⟩':>7}")
Rs = []
for q, pp, kk_ in QPOW:
    a = np.array(acc[q])
    if len(a) == 0:
        continue
    Rm = np.sum(a[:, 0] / a[:, 1]**2) / np.sum(1 / a[:, 1]**2)
    Re_ = 1 / np.sqrt(np.sum(1 / a[:, 1]**2))
    Rs.append((q, kk_, Rm, Re_))
    print(f"{q:>3} {kk_:>2} {len(a):>3} {Rm:>8.3f}±{Re_:.3f} {Rm/kk_:>10.3f} "
          f"{a[:,2].mean():>7.3f}")
print(f"\nplasebo ⟨v⟩ = {np.mean(fake_acc):.4f} (sinyal ölçeği ~0.1-0.4)")
k1 = [r for q, k_, r, e in Rs if k_ >= 2]
print(f"KARAR: k≥2 kuvvetlerinde ⟨R⟩ = {np.mean(k1):.3f} — "
      f"1/k hipotezi R≈1 der, k'sız hipotez R≈k (2-5) der.")
