"""
101e — KALDIRILMIŞ DİP KURTARMASI: HİBRİT MOTOR+MPMATH (20 Ağustos)
==========================================================================
101c bulgusu: ince ızgara (+54-72) yetmedi — sayım sertifikası ada
başına ~285-620 kayıp sıfır gösteriyor ve kurtarma turları düzeltmiyor.
TEŞHİS: kayıplar yakın-çift değil KALDIRILMIŞ DİP — motorun RS-düzeltmesiz
O((qt)^{-1/4}) ≈ 0.05-0.1 hatası dar çift arasındaki |Z| çukurunu sıfır
üstüne kaldırıyor; işaret-taraması hiçbir ızgarada göremez.
ÇÖZÜM (hibrit): kusur bölgelerinde motor-|Z| dipleri bul (işaretsiz yerel
minimum) → her dipte 3 GERÇEK-Z (mpmath Hurwitz) değerlendirmesi →
parabol min < 0 ise çift, kökleri sıfır olarak ekle.
SERTİFİKA HEDEFİ: kalan sürüklenme ±2 içinde (saf S-salınımı).
"""

import numpy as np
import mpmath as mp
import time
from pathlib import Path

exec(open("98_L_motoru.py").read().split('CHI4 = ')[0])

HERE = Path(".").resolve()
mp.mp.dps = 12
CHI4 = {0: 0, 1: 1, 2: 0, 3: -1}
CHI3 = {0: 0, 1: 1, 2: -1}
CHI5 = {0: 0, 1: 1, 2: 1j, 3: -1j, 4: -1}
z6 = np.exp(1j * np.pi / 3)
CHI7 = {0: 0, 1: 1, 3: z6, 2: z6**2, 6: z6**3, 4: z6**4, 5: z6**5}
ISLANDS = [("chi3", 3, CHI3), ("beta", 4, CHI4),
           ("chi5", 5, CHI5), ("chi7", 7, CHI7)]

def gercek_Z(M, tab, t):
    """mpmath Hurwitz ile gerçek Z (motor theta konvansiyonuyla)."""
    s = mp.mpc(0.5, t)
    Lv = mp.mpc(0)
    for a in range(1, M.q):
        if tab[a % M.q] != 0:
            Lv += tab[a % M.q] * mp.zeta(s, mp.mpf(a) / M.q)
    Lv *= mp.power(M.q, -s)
    th = float(M.theta(np.array([t]))[0])
    return float((mp.e**(1j * th) * Lv).real)

def basamak_bul(M, zz, win=80, esik=0.6):
    th = M.theta(zz)
    d = np.arange(len(zz)) - (th - th[0]) / np.pi
    from numpy.lib.stride_tricks import sliding_window_view
    sw = sliding_window_view(d, win)
    med = np.median(sw, axis=1)
    adim = med[win + 1:] - med[:-(win + 1)]
    kotu = np.where(np.abs(adim) > esik)[0] + win // 2
    bolgeler = []
    for i in kotu:
        t0 = zz[max(0, i - win)]
        t1 = zz[min(len(zz) - 1, i + 2 * win)]
        if bolgeler and t0 <= bolgeler[-1][1]:
            bolgeler[-1] = (bolgeler[-1][0], max(bolgeler[-1][1], t1))
        else:
            bolgeler.append((t0, t1))
    return bolgeler, d

for etiket, q, tab in ISLANDS:
    out = HERE / f"101e_{etiket}_zeros.npz"
    if out.exists():
        print(f"[{etiket}] önbellek var, geçiliyor", flush=True); continue
    M = Lmotor(q, tab, 1)
    zz = np.load(HERE / f"101c_{etiket}_zeros.npz")["zeros"]
    t0c = time.time()
    n_eval = 0
    for tur in range(4):
        bolgeler, d = basamak_bul(M, zz)
        dd = d - np.median(d[:200])
        if not bolgeler:
            break
        yeni = []
        for (a, b) in bolgeler:
            gbar = TWO_PI / np.log(q * 0.5 * (a + b) / TWO_PI)
            ts = np.arange(a, b, 0.002 * gbar)
            v = M.Z(ts)
            av = np.abs(v)
            # işaretsiz yerel minimumlar (dip adayları)
            imin = np.where((av[1:-1] < av[:-2]) & (av[1:-1] <= av[2:]))[0] + 1
            def uzak(t):
                j = np.searchsorted(zz, t)
                dl = t - zz[j - 1] if j > 0 else 1e9
                dr = zz[j] - t if j < len(zz) else 1e9
                return min(dl, dr)
            imin = [i for i in imin if av[i] < 0.8
                    and np.sign(v[max(0, i - 60)]) == np.sign(v[min(len(v)-1, i + 60)])
                    and uzak(ts[i]) > 0.05 * gbar]
            imin = sorted(imin, key=lambda i: av[i])[:8]
            for i in imin:
                tstar = ts[i]; h = 0.10 * gbar
                z3 = [gercek_Z(M, tab, tstar - h),
                      gercek_Z(M, tab, tstar),
                      gercek_Z(M, tab, tstar + h)]
                n_eval += 3
                c2 = (z3[0] - 2 * z3[1] + z3[2]) / (2 * h * h)
                c1 = (z3[2] - z3[0]) / (2 * h)
                c0 = z3[1]
                disc = c1 * c1 - 4 * c2 * c0
                # çift kabulü: parabol ekseni kesiyor (disc>0) ve kökler
                # dip civarında; gerçek-min>0 (Lehmer-dip) → disc<0 → ret
                if abs(c2) > 1e-12 and disc > 0:
                    r1 = tstar + (-c1 - np.sqrt(disc)) / (2 * c2)
                    r2 = tstar + (-c1 + np.sqrt(disc)) / (2 * c2)
                    if (abs(r1 - tstar) < 2 * h and abs(r2 - tstar) < 2 * h
                            and abs(r2 - r1) < 0.8 * gbar):
                        yeni += [min(r1, r2), max(r1, r2)]
        if not yeni:
            print(f"  tur {tur+1}: kurtarılacak dip yok, kalan bölge "
                  f"{len(bolgeler)}", flush=True)
            break
        zz = np.unique(np.concatenate([zz, np.array(yeni)]))
        keep = np.concatenate([[True], np.diff(zz) > 5e-4])
        zz = zz[keep]
        print(f"  tur {tur+1}: {len(bolgeler)} bölge, +{len(yeni)} sıfır "
              f"(mpmath {n_eval} eval, {time.time()-t0c:.0f} sn)", flush=True)
    bolgeler, d = basamak_bul(M, zz)
    dd = d - np.median(d[:200])
    print(f"[{etiket}] SERTİFİKA: n={len(zz)}; kalan bölge={len(bolgeler)}; "
          f"sürüklenme [{dd.min():.2f}, {dd.max():.2f}]", flush=True)
    np.savez(out, zeros=zz)
print("BİTTİ", flush=True)
