"""
61 — ÇERÇEVE KONTROLÜ: SUM RULE tmid vs PÜRÜZSÜZ FAZLAR (K4 doğrulaması)
==========================================================================

Denetçi bulgusu (K4): tmid fazları sıfır yer-değiştirmelerinin jitter'ını
taşıyor; fazlar pürüzsüz konumlarda değerlendirilince u seviyesi ~0.96'ya
iniyor ve +%1→%12 drift kayboluyor. Kendi doğrulamamız:

Pürüzsüz çerçeve: Riemann-von Mangoldt akışı — pencerenin ilk sıfırından
başlayıp N(t_k) − N(t_0) = k denklemini Newton ile çöz (S(t) jitter'ı yok,
t_0'daki ofset sabit). Fazlar t_smooth(k)'da; u'lar iki çerçevede kıyas.

Pencere: L=12.45 (denetçinin kullandığı) + L=9.86 (tekrar kontrolü).
"""

import numpy as np
from pathlib import Path

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
A = (np.e**2 - 5) / 2
B0, B1 = 2.7580, -0.0543

def ef_weight(q):
    for p in [2, 3, 5, 7]:
        k = round(np.log(q) / np.log(p))
        if k >= 2 and abs(p**k - q) < 0.5:
            return p**(-k / 2) / k
    return q**-0.5

QLIST = [2, 3, 5, 7, 11, 13, 4, 8, 9, 16, 25, 27, 32, 49, 17, 19, 23, 29, 31, 37, 41, 43, 47]

def rvm_N(t):
    x = t / TWO_PI
    return x * np.log(x / np.e) + 7 / 8

def smooth_positions(g_lo_first, n_gaps_total, k_arr):
    """N(t_k) - N(t_0) = k'yi Newton'la çöz (k: her aralığın orta-indeksi)."""
    t0 = g_lo_first
    N0 = rvm_N(t0)
    t = t0 + k_arr * TWO_PI / np.log(t0 / TWO_PI)  # başlangıç tahmini
    for _ in range(6):
        f = rvm_N(t) - N0 - k_arr
        fp = np.log(t / TWO_PI) / TWO_PI
        t = t - f / fp
    return t

d41 = np.load(HERE / "41_bigT_windows.npz")
for wkey in ["1600k", "120k"]:
    gaps = d41[f"gaps_{wkey}"]; tmid = d41[f"tmid_{wkey}"]; amps = d41[f"amps_{wkey}"]
    Lw = np.log(tmid / TWO_PI); L = float(Lw.mean())
    a_u = amps / np.sqrt(A * Lw + B0 + B1 / Lw)
    a_u /= np.sqrt((a_u**2).mean())
    y = np.log(a_u)
    n = len(y)

    # pürüzsüz konumlar: aralık ortası indeksleri k = 0.5, 1.5, ...
    g_lo_first = tmid[0] - gaps[0] / 2
    k_arr = np.arange(n) + 0.5
    t_sm = smooth_positions(g_lo_first, n, k_arr)
    jitter = tmid - t_sm
    print(f"\n=== L={L:.2f} ({wkey}): jitter rms = {jitter.std():.4f} "
          f"(ort. aralık {gaps.mean():.4f}) ===")
    print(f"{'q':>4} {'u(tmid)':>9} {'u(smooth)':>10}")

    for frame, tv in [("tmid", tmid), ("smooth", t_sm)]:
        cols = [np.ones_like(y)]
        for q in QLIST:
            arg = tv * np.log(q)
            cols += [np.cos(arg), np.sin(arg)]
        X = np.vstack(cols).T
        b, *_ = np.linalg.lstsq(X, y, rcond=None)
        if frame == "tmid":
            u_tmid = [np.hypot(b[1 + 2*i], b[2 + 2*i]) / ef_weight(q)
                      for i, q in enumerate(QLIST)]
        else:
            u_sm = [np.hypot(b[1 + 2*i], b[2 + 2*i]) / ef_weight(q)
                    for i, q in enumerate(QLIST)]
    for i, q in enumerate(QLIST):
        print(f"{q:>4} {u_tmid[i]:>9.3f} {u_sm[i]:>10.3f}")
    ut, us = np.array(u_tmid), np.array(u_sm)
    print(f"  ortalama: tmid {ut.mean():.3f} (aralık {ut.min():.2f}-{ut.max():.2f})"
          f" | smooth {us.mean():.3f} (aralık {us.min():.2f}-{us.max():.2f})")
    # drift: log q'ya karşı eğim
    lq = np.log(np.array(QLIST, dtype=float))
    for name, uu in [("tmid", ut), ("smooth", us)]:
        sl = np.polyfit(lq, uu, 1)[0]
        print(f"  {name} drift eğimi (u vs log q): {sl:+.4f}")
