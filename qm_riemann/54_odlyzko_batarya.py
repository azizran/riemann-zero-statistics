"""
54 — ODLYZKO BATARYASI: L=24.48'DE TÜM ÖLÇÜMLER (17 Ağustos 2026)
===================================================================

Kaldıraç kolu 12.45 → 24.48 (t: 1.6×10⁶ → 2.7×10¹¹). Testler:

  T1  N_eff − L sabitliği: +0.93 (P) / +1.86 (S) burada da tutuyor mu?
  T2  w(2) @ τ=0.0283: 45b kılavuzu ~0.856 der; tam saydamlık 1 derdi.
  T3  Sum rule u(q) ≈ 1 bu yükseklikte de mi?
  T4  τ* öngörüsü: L=24.5'te Q_tepe ≈ e^{0.42L} ~ 29000 → Q≤300 tamamen
      tepe-ÖNCESİ → r*(Q) DÜŞMEDEN tırmanmalı (12.45'te düşüyordu!)
  T5  V_res(L=24.5): lineer uzantı ~0.244+0.031×12 ≈ 0.62 öngörür — bakalım.

Not: 10⁴ aralıkla Q=300 soyması 158 sütun → plasebo tabanı ~%1.6 artık
görünür olur; her şey eş-plaseboyla raporlanır.
"""

import numpy as np
import mpmath as mp
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import pearsonr, spearmanr
from sympy import primerange

rng = np.random.default_rng(54)
HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
A = (np.e**2 - 5) / 2
B0, B1 = 2.7580, -0.0543
T0 = 267653395647

d = np.load(HERE / "53_odlyzko_amps.npz")
gaps, amps, tmid, L = d["gaps"], d["max_amps"], d["t_mid"], float(d["L"])
g_u = gaps * L / TWO_PI
a_u = amps / np.sqrt(A * L + B0 + B1 / L)
a_u /= np.sqrt((a_u**2).mean())
y_a, y_g = np.log(a_u), np.log(g_u)
n = len(g_u)
rP, _ = pearsonr(g_u, a_u); rS, _ = spearmanr(g_u, a_u)
print(f"L = {L:.4f}, {n} aralık")
print(f"HAM: r_P = {rP:.4f} ± {(1-rP**2)/np.sqrt(n):.4f}   r_S = {rS:.4f}")
print(f"unfold kontrol: mean(g̃)={g_u.mean():.4f}, mean(ã²)={(a_u**2).mean():.4f}\n")

# ---------------------------------------------------------------
# T1: CUE N=24..27 → N_eff
# ---------------------------------------------------------------
def haar_unitary_batch(M, N):
    G = (rng.standard_normal((M, N, N)) + 1j * rng.standard_normal((M, N, N))) / np.sqrt(2)
    Q, R = np.linalg.qr(G)
    diag = np.einsum("mii->mi", R)
    return Q * (diag / np.abs(diag))[:, None, :]

def cue_r(N, n_gaps=150000, grid=24, chunk=1000):
    M = int(np.ceil(n_gaps / N))
    gs, ms = [], []
    for s in range(0, M, chunk):
        m = min(chunk, M - s)
        U = haar_unitary_batch(m, N)
        ph = np.sort(np.angle(np.linalg.eigvals(U)), axis=1)
        gp = np.diff(np.concatenate([ph, ph[:, :1] + TWO_PI], axis=1), axis=1)
        uu = np.arange(1, grid + 1) / (grid + 1)
        th = ph[:, :, None] + gp[:, :, None] * uu[None, None, :]
        df = th[:, :, :, None] - ph[:, None, None, :]
        amp = np.prod(2 * np.abs(np.sin(df / 2)), axis=-1)
        gs.append(gp.ravel()); ms.append(amp.max(axis=-1).ravel())
    g = np.concatenate(gs); a = np.concatenate(ms)
    return pearsonr(g, a)[0], spearmanr(g, a)[0]

print("T1 — CUE referansları:")
Ns = np.arange(23, 29)
cueP, cueS = [], []
for N in Ns:
    p_, s_ = cue_r(N)
    cueP.append(p_); cueS.append(s_)
    print(f"  N={N}: P {p_:.4f}  S {s_:.4f}")
NP = np.interp(-rP, -np.array(cueP), Ns.astype(float))
NS = np.interp(-rS, -np.array(cueS), Ns.astype(float))
print(f"  → N_eff−L: Pearson {NP-L:+.2f} (küçük-T: +0.93±0.15)")
print(f"             Spearman {NS-L:+.2f} (küçük-T: +1.86±0.15)\n")

# ---------------------------------------------------------------
# T2/T3: kanallar ve sum rule (çapalı fazlar)
# ---------------------------------------------------------------
mp.mp.dps = 30
def anchored_phase(q):
    return float(mp.fmod(T0 * mp.log(q), 2 * mp.pi))

def wave_cols(qlist):
    cols = []
    for q in qlist:
        ph0 = anchored_phase(q)
        arg = ph0 + tmid * np.log(q)
        cols += [np.cos(arg), np.sin(arg)]
    return cols

def ef_weight(q):
    for p in [2, 3, 5, 7]:
        k = round(np.log(q) / np.log(p))
        if k >= 2 and abs(p**k - q) < 0.5:
            return p**(-k / 2) / k
    return q**-0.5

PRIMES = [2, 3, 5, 7, 11, 13]
SUMRULE_Q = PRIMES + [4, 9, 25, 8, 27, 49]

# w kanalı (g̃ kontrollü) ve toplam u (kontrolsüz)
pc = wave_cols(PRIMES)
Xw = np.vstack([np.ones(n), g_u, g_u**2] + pc).T
bw = np.linalg.lstsq(Xw, y_a, rcond=None)[0]
resw = y_a - Xw @ bw
sew = np.sqrt(resw.var() * np.diag(np.linalg.inv(Xw.T @ Xw)))
Xg = np.vstack([np.ones(n)] + pc).T
bg = np.linalg.lstsq(Xg, y_g, rcond=None)[0]

print("T2 — kanallar (L=24.48):")
tau_list, w_list = [], []
for i, p in enumerate(PRIMES):
    w_p = bw[3 + 2 * i] / p**-0.5
    se_p = sew[3 + 2 * i] / p**-0.5
    v_p = np.hypot(bg[1 + 2 * i], bg[2 + 2 * i]) / p**-0.5
    tau = np.log(p) / L
    guide = 0.940 - 2.967 * tau
    tau_list.append(tau); w_list.append(w_p)
    print(f"  p={p:>2} (τ={tau:.4f}): w = {w_p:+.3f} ± {se_p:.3f}"
          f"   [45b kılavuzu: {guide:.3f}]   |v| = {v_p:.3f}")

pcs = wave_cols(SUMRULE_Q)
Xu = np.vstack([np.ones(n)] + pcs).T
bu = np.linalg.lstsq(Xu, y_a, rcond=None)[0]
print("\nT3 — sum rule:")
for i, q in enumerate(SUMRULE_Q):
    u = np.hypot(bu[1 + 2 * i], bu[2 + 2 * i]) / ef_weight(q)
    print(f"  q={q:>2}: |u| = {u:.3f}")

# ---------------------------------------------------------------
# T4/T5: soyma dizisi + plasebo
# ---------------------------------------------------------------
def prime_powers_upto(Q):
    out = []
    for p in primerange(2, Q + 1):
        pk = p
        while pk <= Q:
            out.append(float(pk)); pk *= p
    return sorted(out)

def strip_r(qlist_or_freqs, is_freq=False):
    if is_freq:
        cols = []
        for om in qlist_or_freqs:
            arg = tmid * om
            cols += [np.cos(arg), np.sin(arg)]
    else:
        cols = wave_cols(qlist_or_freqs)
    C = np.vstack([np.ones(n)] + cols).T
    ba_ = np.linalg.lstsq(C, y_a, rcond=None)[0]
    bg_ = np.linalg.lstsq(C, y_g, rcond=None)[0]
    ya = y_a - C[:, 1:] @ ba_[1:]
    yg = y_g - C[:, 1:] @ bg_[1:]
    ar = np.exp(ya); ar /= np.sqrt((ar**2).mean())
    gr = np.exp(yg); gr *= g_u.mean() / gr.mean()
    va = float((C[:, 1:] @ ba_[1:]).var())
    return pearsonr(gr, ar)[0], spearmanr(gr, ar)[0], va

print("\nT4 — soyma dizisi (τ* testi: hep TIRMANMALI, düşüş olmamalı):")
print(f"  {'Q':>5} {'τ_Q':>6} {'r*_P':>8} {'r*_S':>8} {'V_soyulan':>10} {'plasebo r_P':>12}")
V_tot = float(y_a.var())
for Q in [13, 50, 150, 300]:
    qs = prime_powers_upto(Q)
    rp_, rs_, va = strip_r(qs)
    # eş sayıda plasebo frekansı
    used = np.log(np.array(qs))
    plc = []
    while len(plc) < len(qs):
        c = rng.uniform(np.log(2) * 0.9, np.log(320))
        if np.abs(used - c).min() > 0.01 and all(abs(x - c) > 0.01 for x in plc):
            plc.append(c)
    rpp, _, vap = strip_r(np.array(plc), is_freq=True)
    print(f"  {Q:>5} {np.log(Q)/L:>6.3f} {rp_:>8.4f} {rs_:>8.4f} {va:>10.4f} {rpp:>12.4f}")

_, _, va300 = strip_r(prime_powers_upto(300))
V_res = V_tot - va300
print(f"\nT5 — V_res(L=24.48) = {V_tot:.4f} − {va300:.4f} = {V_res:.4f}")
print(f"  küçük-T lineer uzantısı öngörüsü: 0.244 + 0.031×(24.48−12.45) ≈ "
      f"{0.244 + 0.031*(24.48-12.45):.3f}")
print(f"  (plasebo tabanı ~{vap:.4f} — V_res buna göre yorumlanmalı)")
