"""
105b — YENİ-ADA KAMPANYASI + SERTİFİKA BORU HATTI (25 Ağustos 2026)
==========================================================================
Üç yeni L-fonksiyon adası türetilir; 101c/101e/101f boru hattının
TAMAMI uygulanır (ince ızgara + sayım-sertifikası + hibrit mpmath
dip-kurtarma + temizlik + kalan-bölge kaydı).

ADALAR (105a kapısından geçti):
  chi5e  q=5, a=0 (ÇİFT), reel Legendre: χ(1,2,3,4)=(+1,−1,−1,+1)
         ölü: 5, 25 | ilk sağ kalan çizgi log 2
  chi8e  q=8, a=0 (ÇİFT), reel: +1 için n≡±1(8), −1 için n≡±3(8)
         ölü: 2,4,8,... | ilk sağ kalan çizgi log 3
  chi8o  q=8, a=1 (TEK),  reel: +1 için n≡1,3, −1 için n≡5,7
         ölü: 2,4,8,... | ilk sağ kalan çizgi log 3
T1 seçimi ~60-75k sıfır hedefine göre (chi5e 48000, chi8* 42000).

BORU HATTI (101'in dersleri, birebir aynı makine):
  A) ince ızgara grid_frac=0.03, T0=200 → 105b_{ada}_ham.npz
  B) SAYIM SERTİFİKASI d_i = i − (θ(z_i)−θ(z_0))/π; basamak bölgelerinde
     grid_frac=0.005 kurtarma (maks 3 tur)
  C) KALDIRILMIŞ DİP kurtarması (101e): işaretsiz |Z| minimumu +
     3 gerçek-Z (mpmath Hurwitz) + parabol kökleri
  D) TEMİZLİK (101f): pozitif basamakta sahte-kopya budama (mpmath
     orta-nokta işaret testiyle), negatif basamakta uyarlanır eşikli
     (0.45×medyan|Z|) ikinci kurtarma
  ÇIKTI: 105b_{ada}_zeros.npz  (zeros + bolgeler)  — 105c/105d bunu okur;
  105d kalan bölgeleri KESER (düzlük segmentasyonu).

NOT (dürüstlük): motorun RS-düzeltmesiz konum kayması (105a G5:
|Δγ|/⟨g⟩ ~1e-2 @t≈1e3, ~1e-3 @t≈3e4) yayımlanmış dört adada da
mevcuttur; benek/faz ölçümleri duyarsız (98-S2), donma ölçümü
sertifika+segmentasyonla korunur.

==========================================================================
SONUÇ (25 Ağustos — üç ada da sertifikalandı; toplam ~13 dk)
==========================================================================
ada    T1     ham(0.03)  B(+ince)  C(+dip)  SON n   kalan bölge  sürüklenme
chi5e  48000    72451     +2         +376   72829       2        [-4.53,+0.83]
chi8e  42000    65816     +2         +130   65948       2        [-0.63,+4.96]
chi8o  42000    65503     +8         +432   65943       0        [-1.01,+2.41]
mpmath eval: 703 / 217 / 741. Süre: ~10 / ~3 / ~9 dk.
101'in dersi BİREBİR tekrarlandı: ince ızgara yalnız +2/+2/+8 buldu,
sayım sertifikası ise 379/128/435 kayıp gösterdi; hepsi KALDIRILMIŞ DİP
çıktı ve hibrit mpmath kurtarmasıyla geri geldi.
DÜRÜST NOT (105e'de yakalandı): sertifikanın basamak dedektörü (80'lik
medyan) KISA çukur kusurlarını göremiyor — chi8e'de t≈19352'de 2 sıfır
kaybı + hemen telafi, net basamak bırakmadan geçti ve 105d'nin ilk
koşusunda eşik-altı D'yi 0.02→0.48 şişirdi. Kapı 105d/105e'ye eklendi;
sıfır listeleri DEĞİŞMEDİ (kusurlu bölge analizden kesiliyor).
"""

import numpy as np
import mpmath as mp
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
exec(open(HERE / "98_L_motoru.py").read().split('CHI4 = ')[0])
mp.mp.dps = 12

CHI5E = {0: 0, 1: 1, 2: -1, 3: -1, 4: 1}
CHI8E = {0: 0, 1: 1, 2: 0, 3: -1, 4: 0, 5: -1, 6: 0, 7: 1}
CHI8O = {0: 0, 1: 1, 2: 0, 3: 1, 4: 0, 5: -1, 6: 0, 7: -1}

ISLANDS = [("chi5e", 5, CHI5E, 0, 48000.0),
           ("chi8e", 8, CHI8E, 0, 42000.0),
           ("chi8o", 8, CHI8O, 1, 42000.0)]


def gercek_Z(M, tab, t):
    s = mp.mpc(0.5, t)
    Lv = mp.mpc(0)
    for r in range(1, M.q):
        if tab[r % M.q] != 0:
            Lv += tab[r % M.q] * mp.zeta(s, mp.mpf(r) / M.q)
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


for etiket, q, tab, a, T1 in ISLANDS:
    son = HERE / f"105b_{etiket}_zeros.npz"
    if son.exists():
        print(f"[{etiket}] nihai önbellek var, geçiliyor", flush=True)
        continue
    M = Lmotor(q, tab, a)
    print(f"\n===== {etiket} (q={q}, a={a}, T1={T1:.0f}) =====", flush=True)

    # ---------- A) ince ızgara ----------
    ham = HERE / f"105b_{etiket}_ham.npz"
    if ham.exists():
        zz = np.load(ham)["zeros"]
        print(f"  A) ham önbellek: {len(zz)} sıfır", flush=True)
    else:
        t0 = time.time()
        zz = M.sifir_bul(200.0, T1, grid_frac=0.03)
        np.savez(ham, zeros=zz)
        print(f"  A) ince ızgara (0.03): {len(zz)} sıfır, "
              f"{time.time()-t0:.0f} sn", flush=True)
    n_ham = len(zz)

    # ---------- B) sayım sertifikası + ince kurtarma ----------
    t0 = time.time()
    for tur in range(3):
        bolgeler, d = basamak_bul(M, zz)
        if not bolgeler:
            break
        yeni = [zz]
        for (x, y) in bolgeler:
            yeni.append(M.sifir_bul(x - 0.5, y + 0.5, grid_frac=0.005))
        zz = np.unique(np.concatenate(yeni))
        keep = np.concatenate([[True], np.diff(zz) > 1e-3])
        zz = zz[keep]
        print(f"     B tur {tur+1}: {len(bolgeler)} şüpheli bölge → n={len(zz)}",
              flush=True)
    bolgeler, d = basamak_bul(M, zz)
    dd = d - np.median(d[:200])
    zc = zz.copy()          # 101c muadili (kopya-budamada referans)
    print(f"  B) SERTİFİKA-1: n={len(zz)} ({len(zz)-n_ham:+d}); "
          f"kalan bölge={len(bolgeler)}; sürüklenme "
          f"[{dd.min():.2f}, {dd.max():.2f}]  ({time.time()-t0:.0f} sn)",
          flush=True)

    # ---------- C) kaldırılmış dip kurtarması (101e) ----------
    t0c = time.time(); n_eval = 0; n_c0 = len(zz)
    for tur in range(4):
        bolgeler, d = basamak_bul(M, zz)
        if not bolgeler:
            break
        yeni = []
        for (x, y) in bolgeler:
            gbar = TWO_PI / np.log(q * 0.5 * (x + y) / TWO_PI)
            ts = np.arange(x, y, 0.002 * gbar)
            v = M.Z(ts); av = np.abs(v)
            imin = np.where((av[1:-1] < av[:-2]) & (av[1:-1] <= av[2:]))[0] + 1

            def uzak(t, _zz=zz):
                j = np.searchsorted(_zz, t)
                dl = t - _zz[j - 1] if j > 0 else 1e9
                dr = _zz[j] - t if j < len(_zz) else 1e9
                return min(dl, dr)
            imin = [i for i in imin if av[i] < 0.8
                    and np.sign(v[max(0, i - 60)]) ==
                        np.sign(v[min(len(v) - 1, i + 60)])
                    and uzak(ts[i]) > 0.05 * gbar]
            imin = sorted(imin, key=lambda i: av[i])[:8]
            for i in imin:
                tstar = ts[i]; h = 0.10 * gbar
                z3 = [gercek_Z(M, tab, tstar - h), gercek_Z(M, tab, tstar),
                      gercek_Z(M, tab, tstar + h)]
                n_eval += 3
                c2 = (z3[0] - 2 * z3[1] + z3[2]) / (2 * h * h)
                c1 = (z3[2] - z3[0]) / (2 * h)
                c0 = z3[1]
                disc = c1 * c1 - 4 * c2 * c0
                if abs(c2) > 1e-12 and disc > 0:
                    r1 = tstar + (-c1 - np.sqrt(disc)) / (2 * c2)
                    r2 = tstar + (-c1 + np.sqrt(disc)) / (2 * c2)
                    if (abs(r1 - tstar) < 2 * h and abs(r2 - tstar) < 2 * h
                            and abs(r2 - r1) < 0.8 * gbar):
                        yeni += [min(r1, r2), max(r1, r2)]
        if not yeni:
            print(f"     C tur {tur+1}: kurtarılacak dip yok, kalan bölge "
                  f"{len(bolgeler)}", flush=True)
            break
        zz = np.unique(np.concatenate([zz, np.array(yeni)]))
        keep = np.concatenate([[True], np.diff(zz) > 5e-4])
        zz = zz[keep]
        print(f"     C tur {tur+1}: {len(bolgeler)} bölge, +{len(yeni)} sıfır "
              f"(mpmath {n_eval} eval, {time.time()-t0c:.0f} sn)", flush=True)
    bolgeler, d = basamak_bul(M, zz)
    dd = d - np.median(d[:200])
    print(f"  C) SERTİFİKA-2: n={len(zz)} ({len(zz)-n_c0:+d} dip-kurtarma); "
          f"kalan bölge={len(bolgeler)}; sürüklenme "
          f"[{dd.min():.2f}, {dd.max():.2f}]", flush=True)

    # ---------- D) temizlik (101f) ----------
    t0d = time.time(); n_d0 = len(zz)
    for tur in range(4):
        bolgeler, d = basamak_bul(M, zz)
        if not bolgeler:
            break
        degisti = False
        sil, ekle = [], []
        for (x, y) in bolgeler:
            adim = yerel_adim(zz, d, x, y)
            gbar = TWO_PI / np.log(q * 0.5 * (x + y) / TWO_PI)
            if adim > 0.5:
                m = (zz >= x) & (zz <= y)
                zw = zz[m]
                j = np.searchsorted(zc, zw)
                dist = np.minimum(
                    np.abs(zw - zc[np.clip(j, 0, len(zc) - 1)]),
                    np.abs(zw - zc[np.clip(j - 1, 0, len(zc) - 1)]))
                ins = zw[dist > 1e-6]
                if len(ins) == 0:
                    continue
                kum = np.split(ins, np.where(np.diff(ins) > 0.5 * gbar)[0] + 1)
                for c in kum:
                    if len(c) < 2:
                        sil += list(c); degisti = True; continue
                    orta = len(c) // 2
                    cift = c[max(0, orta - 1):max(0, orta - 1) + 2]
                    tl, th_ = cift[0], cift[-1]
                    sh = np.sign(M.Z(np.array([tl - 0.3 * gbar,
                                               th_ + 0.3 * gbar])))
                    zmid = gercek_Z(M, tab, 0.5 * (tl + th_)); n_eval += 1
                    fazla = [w for w in c if w not in cift]
                    sil += fazla
                    if sh[0] == sh[1] and np.sign(zmid) == sh[0]:
                        sil += list(cift)
                    degisti = True
            elif adim < -0.5:
                ts = np.arange(x, y, 0.002 * gbar)
                v = M.Z(ts); av = np.abs(v)
                esik = 0.45 * np.median(av)
                imin = np.where((av[1:-1] < av[:-2]) &
                                (av[1:-1] <= av[2:]))[0] + 1

                def uzak2(t, _zz=zz):
                    j = np.searchsorted(_zz, t)
                    dl = t - _zz[j - 1] if j > 0 else 1e9
                    dr = _zz[j] - t if j < len(_zz) else 1e9
                    return min(dl, dr)
                imin = [i for i in imin if av[i] < esik
                        and np.sign(v[max(0, i - 60)]) ==
                            np.sign(v[min(len(v) - 1, i + 60)])
                        and uzak2(ts[i]) > 0.05 * gbar]
                imin = sorted(imin, key=lambda i: av[i])[:12]
                for i in imin:
                    tstar = ts[i]; h = 0.10 * gbar
                    z5t = [tstar - h, tstar - h / 2, tstar,
                           tstar + h / 2, tstar + h]
                    z5 = [gercek_Z(M, tab, w) for w in z5t]; n_eval += 5
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
        print(f"     D tur {tur+1}: −{len(sil)} +{len(ekle)} "
              f"(eval {n_eval}, {time.time()-t0d:.0f} sn)", flush=True)

    bolgeler, d = basamak_bul(M, zz)
    dd = d - np.median(d[:200])
    Leff = float(np.log(q * zz.mean() / TWO_PI))
    print(f"  D) SON SERTİFİKA: n={len(zz)} ({len(zz)-n_d0:+d} temizlik); "
          f"KALAN BÖLGE={len(bolgeler)}; sürüklenme "
          f"[{dd.min():.2f}, {dd.max():.2f}]; L_eff≈{Leff:.2f}; "
          f"toplam mpmath eval={n_eval}", flush=True)
    np.savez(son, zeros=zz,
             bolgeler=np.array(bolgeler if bolgeler else []).reshape(-1, 2),
             q=q, a=a, T1=T1)
print("\nBİTTİ", flush=True)
