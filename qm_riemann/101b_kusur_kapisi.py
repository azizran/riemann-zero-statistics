"""
101b — KUSUR KAPISI: 101'İN ADA "ÇÖZÜLMEMİŞ" SONUCU GERÇEK Mİ? (20 Ağu)
==========================================================================
101 dört adada τ=0.04'te D_spont ≈ 0.96-0.98 ölçtü (ζ: 0.041 — donmuş).
ŞÜPHE: motorun ~%0.3 kaçık yakın-çifti örgü KUSURU üretir; kusur hem
ρ'ya hem ds'e korele sıçrama sokar → sahte "çözülme".
KONTROLLER (ön-mühürlü beklentiler):
  C1  ζ 120k penceresi BENİM estimatörümle → 93 ile aynıysa ~0.04-0.05
      çıkmalı (kod-eşdeğerlik kapısı).
  C2  Aynı pencereye %0.3 rastgele sıfır SİLME enjeksiyonu → kusur
      hipotezi doğruysa τ=0.04'te D 0.04 → O(1)'e fırlamalı.
  C3  Poisson süreci (yapısız) → kinematik kilit referansı (~1/cos(κ/2)).
  C4  χ₃ DAR pencere (üst %20 log) → genişlik artefaktı mı kontrolü.
"""

import numpy as np
from pathlib import Path

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
rng = np.random.default_rng(101)

PKS = [2,3,4,5,7,8,9,11,13,16,17,19,23,25,27,29,31,32,37,41,43,47,49,53,
       59,61,64,67,71,73,79,81,83,89,97,101,103,107,109,113,121,125,127,128]
LINES_ALL = [np.log(q) for q in PKS]

def D_est(z, dens_log_arg, taus):
    """z: sıfırlar; dens_log_arg: q_eff öyle ki yoğunluk log(q_eff·t/2π)/2π."""
    gaps = np.diff(z)
    mids = 0.5 * (z[:-1] + z[1:])
    L = float(np.log(dens_log_arg * mids / TWO_PI).mean())
    ds = gaps * np.log(dens_log_arg * mids / TWO_PI) / TWO_PI - 1
    tt = (mids - mids.mean()) / (mids[-1] - mids[0])
    P = np.vstack([np.ones_like(tt), tt, tt**2, tt**3]).T
    ds = ds - P @ np.linalg.lstsq(P, ds, rcond=None)[0]
    out = []
    for tau0 in taus:
        oms = tau0 * L + np.linspace(-0.018 * L, 0.018 * L, 140)
        oms = np.array([o for o in oms
                        if min(abs(o - l) for l in LINES_ALL) > 0.01])
        num = 0.0 + 0j; den = 0.0
        for s0 in range(0, len(oms), 35):
            ob = oms[s0:s0 + 35]
            rr = np.exp(1j * np.outer(ob, z)).sum(axis=1)
            GG = (np.exp(1j * np.outer(ob, mids)) * ds[None, :]).sum(axis=1)
            num += (GG * np.conj(rr)).sum()
            den += (np.abs(rr)**2).sum()
        kap = 2 * np.pi * tau0
        out.append(abs(num) / ((2 * np.sin(kap / 2) / kap) * den))
    return L, out

TAUS = [0.04, 0.08, 0.15, 0.30]
print(f"{'kontrol':<28} {'L':>5} " + " ".join(f"τ={t:.2f}" for t in TAUS))

# C1 — ζ 120k penceresi (93'ün verisi, benim kodum)
d41 = np.load(HERE / "41_bigT_windows.npz")
gz, tm = d41["gaps_120k"], d41["tmid_120k"]
zz = np.empty(len(gz) + 1)
zz[0] = tm[0] - gz[0] / 2
zz[1:] = zz[0] + np.cumsum(gz)
L, D = D_est(zz, 1.0, TAUS)
print(f"{'C1 ζ-120k (temiz)':<28} {L:>5.2f} " + " ".join(f"{d:>6.3f}" for d in D))

# C2 — aynı pencere, %0.3 rastgele silme (kusur enjeksiyonu)
for frac in [0.001, 0.003]:
    kill = rng.random(len(zz)) < frac
    zk = zz[~kill]
    L, D = D_est(zk, 1.0, TAUS)
    print(f"{f'C2 ζ-120k sil %{100*frac:.1f}':<28} {L:>5.2f} "
          + " ".join(f"{d:>6.3f}" for d in D))

# C3 — Poisson (yapısız, kinematik referans); yoğunluğu ζ-120k'ya benzet
n = 30000
lam = np.log(tm.mean() / TWO_PI) / TWO_PI
zp = tm[0] + np.cumsum(rng.exponential(1 / lam, n))
L, D = D_est(zp, 1.0, TAUS)
kapv = [1 / np.cos(np.pi * t) for t in TAUS]
print(f"{'C3 Poisson':<28} {L:>5.2f} " + " ".join(f"{d:>6.3f}" for d in D))
print(f"{'    (kinematik 1/cos(κ/2))':<28} {'':>5} "
      + " ".join(f"{k:>6.3f}" for k in kapv))

# C4 — χ₃ dar pencere (üst %20 log)
zc = np.load(HERE / "99_chi3_zeros.npz")["zeros"]
lo = np.exp(np.log(zc[0] + 1) + 0.8 * (np.log(zc[-1]) - np.log(zc[0] + 1)))
zn = zc[zc >= lo]
L, D = D_est(zn, 3.0, TAUS)
print(f"{'C4 χ₃ dar (n=%d)' % len(zn):<28} {L:>5.2f} "
      + " ".join(f"{d:>6.3f}" for d in D))
