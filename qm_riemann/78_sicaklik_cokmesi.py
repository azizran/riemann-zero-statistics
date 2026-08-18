"""
78 — SICAKLIK ÇÖKMESİ (72'nin V2 analizi, repoya kayıt) (19 Ağustos)
======================================================================
Denetim S5: pencere-bazlı c(L)·s modeli scriptte yoktu. Bu script Not 3'ün
"window temperature absorbs the fixed-τ scatter" sayısını üretir:
model w = A·(1−2τ)·exp(−s·c_jitter(L)·τ²)·[τ<½] + B, A ve ORTAK s serbest
(c_jitter pencere-bazlı ölçüm; s fit edilir — metinde de öyle yazar).
Çıktı: RMS ve soğuk/sıcak pencere artık-yarılması, sabit-c kıyasıyla.
"""

import numpy as np
from pathlib import Path
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
exec(open(HERE / "72_sicak_kristal.py").read().split("# ---- T2")[0]
     .replace('print(f"  L={L:5.2f}', '_ = (f"  L={L:5.2f}'))

cs_arr, tau_l, w_l, sw_l = [], [], [], []
for (gaps, amps, tmid), c_w in zip(WNDS, cs):
    for t_, w_, s_ in w_channel(gaps, amps, tmid):
        cs_arr.append(c_w); tau_l.append(t_); w_l.append(w_); sw_l.append(s_)
cs_arr = np.array(cs_arr); tau = np.array(tau_l)
w = np.array(w_l); sw = np.array(sw_l)

def mod(A_, cvec):
    return A_ * (1 - 2 * tau) * np.exp(-cvec * tau**2) * (tau < 0.5) + B_PLATO

def split(resid):
    lo = resid[cs_arr < np.median(cs_arr)].mean()
    hi = resid[cs_arr >= np.median(cs_arr)].mean()
    return lo, hi

f0 = least_squares(lambda p: (mod(p[0], np.full_like(tau, c_jit := np.mean(cs_arr))) - w) / sw, [1.05])
r0 = mod(f0.x[0], np.full_like(tau, np.mean(cs_arr))) - w
lo0, hi0 = split(r0)
print(f"SABİT-c  : A={f0.x[0]:.3f}  RMS={np.sqrt(np.mean(r0**2)):.4f}  "
      f"soğuk/sıcak artık: {lo0:+.4f}/{hi0:+.4f}")

f2 = least_squares(lambda p: (mod(p[0], p[1] * cs_arr) - w) / sw, [1.05, 3.0],
                   bounds=([0.3, 0.2], [3, 8]))
r2 = mod(f2.x[0], f2.x[1] * cs_arr) - w
lo2, hi2 = split(r2)
print(f"PENCERE-c: A={f2.x[0]:.3f}, ortak s={f2.x[1]:.3f}  "
      f"RMS={np.sqrt(np.mean(r2**2)):.4f}  soğuk/sıcak artık: {lo2:+.4f}/{hi2:+.4f}")
print("(Not 3 iddiası: pencere-sıcaklığı sabit-τ yarılmasını ±0.003'e indirir)")
