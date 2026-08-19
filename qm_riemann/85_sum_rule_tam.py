"""
85 — SUM RULE, KOŞULSUZ TAM-TABAN (revizyon denetimi K1) (20 Ağustos 2026)
==========================================================================
Denetçi bulgusu: u(4)=0.70 (83) ile sum rule u≈1 (51) arasındaki fark
taban değil KOŞULLAMA farkı — 83'ün w-regresyonu g_u, g_u² kovaryatlı
("sabit boşlukta" içerik), 51'in sum rule'u koşulsuz (toplam içerik).
Bu script sum rule'un doğru nesnesini üretir: KOŞULSUZ + TAM-TABAN
u_q = √(a²+b²) / (p^{-k/2}/k), altı 41 penceresi, 23 frekans.
Beklenti (denetçi ön-ölçümü): seviye ≈0.95–1.00, driftsiz — eski
"+1..+12% drift"in çoğu taban transferiymiş.
"""

import numpy as np
from sympy import primerange, factorint

exec(open("83_buyuk_yeniden_olcum.py").read().split("# ---- pencereler")[0])

FREQS = [2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 17, 19, 23, 25, 27, 29, 31, 32,
         37, 41, 43, 47, 49]

d41 = np.load(HERE / "41_bigT_windows.npz")
K41 = sorted({x.split("_")[1] for x in d41.files}, key=lambda s: int(s[:-1]))
WNDS = [(d41[f"gaps_{k}"], d41[f"amps_{k}"], d41[f"tmid_{k}"]) for k in K41[:6]]

acc = {q: [] for q in FREQS}
for gaps, amps, tmid in WNDS:
    ya, yg, g_u, L = unfold(gaps, amps, tmid)
    qs = sorted(set(FREQS) | set(tam_taban(L, 0.45)))
    bw, sew, _ = chunked_reg(ya, tmid, qs, False, g_u=None)   # KOŞULSUZ
    for q in FREQS:
        i = qs.index(q)
        fac = factorint(q)
        (pp, kk), = fac.items()
        wgt = pp ** (-kk / 2) / kk
        amp = float(np.hypot(bw[1 + 2*i], bw[2 + 2*i]))
        se = float(sew[1 + 2*i])
        acc[q].append((amp / wgt, se / wgt))

print("KOŞULSUZ TAM-TABAN SUM RULE (altı 41 penceresi, ağırlıklı ort):")
print(f"{'q':>3} {'u_q':>7} {'±':>6}")
us = []
for q in FREQS:
    a = np.array(acc[q])
    um = np.sum(a[:, 0] / a[:, 1]**2) / np.sum(1 / a[:, 1]**2)
    ue = 1 / np.sqrt(np.sum(1 / a[:, 1]**2))
    us.append((q, um, ue))
    print(f"{q:>3} {um:>7.3f} {ue:>6.3f}")
arr = np.array([(u, e) for q, u, e in us])
asal = np.array([(u, e) for q, u, e in us if q in
                 (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47)])
print(f"\nasal aralığı: {asal[:,0].min():.3f}–{asal[:,0].max():.3f}  "
      f"(23 frekans: {arr[:,0].min():.3f}–{arr[:,0].max():.3f})")
lo = [u for q, u, e in us if q <= 8]
hi = [u for q, u, e in us if q >= 29]
print(f"drift kontrolü: küçük-q ort {np.mean(lo):.3f} vs büyük-q ort {np.mean(hi):.3f}")
