"""
101f — SERTİFİKA TEMİZLİĞİ: KOPYA BUDAMA + DOĞRULAMA + 2. KURTARMA (20 Ağu)
==========================================================================
101e yaraları: (1) inatçı diplerde tekrar-ekleme (chi3 +5, chi5 +6
pozitif sürüklenme = sahte kopyalar), (2) chi7'de ~16 kayıp (mutlak 0.8
eşiği adanın |Z| ölçeğine dar).
İŞLEM: pozitif basamaklı bölgelerde eklenen kümeleri TEK çifte indir +
mpmath orta-nokta işaret testiyle doğrula (sahteyse at); negatif
basamaklı bölgelerde uyarlanır eşikli (0.45×medyan|Z|) ikinci kurtarma.
ÇIKTI: 101f_{ada}_zeros.npz (zeros + kalan şüpheli bölgeler 'bolgeler'
dizisi — 101d bu bölgeleri KESER, analize sokmaz).
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

def yerel_adim(zz, d, a, b):
    ia, ib = np.searchsorted(zz, a), np.searchsorted(zz, b)
    sol = np.median(d[max(0, ia - 120):max(1, ia - 10)])
    sag = np.median(d[min(len(d) - 1, ib + 10):min(len(d), ib + 120)])
    return sag - sol

for etiket, q, tab in ISLANDS:
    out = HERE / f"101f_{etiket}_zeros.npz"
    if out.exists():
        print(f"[{etiket}] önbellek var", flush=True); continue
    M = Lmotor(q, tab, 1)
    zc = np.load(HERE / f"101c_{etiket}_zeros.npz")["zeros"]
    zz = np.load(HERE / f"101e_{etiket}_zeros.npz")["zeros"]
    t0c = time.time(); n_eval = 0
    for tur in range(4):
        bolgeler, d = basamak_bul(M, zz)
        if not bolgeler:
            break
        degisti = False
        sil, ekle = [], []
        for (a, b) in bolgeler:
            adim = yerel_adim(zz, d, a, b)
            gbar = TWO_PI / np.log(q * 0.5 * (a + b) / TWO_PI)
            if adim > 0.5:
                # eklenen kümeleri bul (101c'de olmayanlar)
                m = (zz >= a) & (zz <= b)
                zw = zz[m]
                j = np.searchsorted(zc, zw)
                dist = np.minimum(np.abs(zw - zc[np.clip(j, 0, len(zc)-1)]),
                                  np.abs(zw - zc[np.clip(j - 1, 0, len(zc)-1)]))
                ins = zw[dist > 1e-6]
                if len(ins) == 0:
                    continue
                kum = np.split(ins, np.where(np.diff(ins) > 0.5 * gbar)[0] + 1)
                for c in kum:
                    if len(c) < 2:
                        sil += list(c); degisti = True; continue
                    orta = len(c) // 2
                    cift = c[max(0, orta - 1):max(0, orta - 1) + 2]
                    # doğrulama: motor omuz işareti vs mpmath orta işaret
                    tl, th_ = cift[0], cift[-1]
                    sh = np.sign(M.Z(np.array([tl - 0.3 * gbar,
                                               th_ + 0.3 * gbar])))
                    zmid = gercek_Z(M, tab, 0.5 * (tl + th_)); n_eval += 1
                    fazla = [x for x in c if x not in cift]
                    sil += fazla
                    if sh[0] == sh[1] and np.sign(zmid) == sh[0]:
                        sil += list(cift)          # sahte çift
                    degisti = True
            elif adim < -0.5:
                # uyarlanır eşikli kurtarma
                ts = np.arange(a, b, 0.002 * gbar)
                v = M.Z(ts); av = np.abs(v)
                esik = 0.45 * np.median(av)
                imin = np.where((av[1:-1] < av[:-2]) &
                                (av[1:-1] <= av[2:]))[0] + 1
                def uzak(t):
                    j = np.searchsorted(zz, t)
                    dl = t - zz[j - 1] if j > 0 else 1e9
                    dr = zz[j] - t if j < len(zz) else 1e9
                    return min(dl, dr)
                imin = [i for i in imin if av[i] < esik
                        and np.sign(v[max(0, i - 60)]) ==
                            np.sign(v[min(len(v) - 1, i + 60)])
                        and uzak(ts[i]) > 0.05 * gbar]
                imin = sorted(imin, key=lambda i: av[i])[:12]
                for i in imin:
                    tstar = ts[i]; h = 0.10 * gbar
                    z5t = [tstar - h, tstar - h / 2, tstar,
                           tstar + h / 2, tstar + h]
                    z5 = [gercek_Z(M, tab, x) for x in z5t]; n_eval += 5
                    A = np.vstack([np.ones(5), np.array(z5t) - tstar,
                                   (np.array(z5t) - tstar)**2]).T
                    c0, c1, c2 = np.linalg.lstsq(A, np.array(z5),
                                                 rcond=None)[0]
                    disc = c1 * c1 - 4 * c2 * c0
                    if abs(c2) > 1e-12 and disc > 0:
                        r1 = tstar + (-c1 - np.sqrt(disc)) / (2 * c2)
                        r2 = tstar + (-c1 + np.sqrt(disc)) / (2 * c2)
                        if (abs(r1 - tstar) < 2 * h and abs(r2 - tstar) < 2 * h
                                and 1e-4 < abs(r2 - r1) < 0.8 * gbar):
                            ekle += [min(r1, r2), max(r1, r2)]
                            degisti = True
        if not degisti:
            break
        if sil:
            zz = zz[~np.isin(zz, np.array(sil))]
        if ekle:
            zz = np.unique(np.concatenate([zz, np.array(ekle)]))
            keep = np.concatenate([[True], np.diff(zz) > 5e-4])
            zz = zz[keep]
        print(f"  tur {tur+1}: −{len(sil)} +{len(ekle)} "
              f"(eval {n_eval}, {time.time()-t0c:.0f} sn)", flush=True)
    bolgeler, d = basamak_bul(M, zz)
    dd = d - np.median(d[:200])
    print(f"[{etiket}] SON SERTİFİKA: n={len(zz)}; kalan bölge="
          f"{len(bolgeler)}; sürüklenme [{dd.min():.2f}, {dd.max():.2f}]",
          flush=True)
    np.savez(out, zeros=zz,
             bolgeler=np.array(bolgeler if bolgeler else []).reshape(-1, 2))
print("BİTTİ", flush=True)
