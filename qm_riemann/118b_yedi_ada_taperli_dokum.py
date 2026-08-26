"""
118b — YEDİ GL(1) ADASININ TAPER'LI YENİDEN DÖKÜMÜ + SERTİFİKA KIYASI
       (26 Ağustos 2026)
==========================================================================
ÖN-MÜHÜR (ölçümden ÖNCE yazıldı)
==========================================================================
118a taper'lı GL(1) motorunu kapılardan geçirdi. Bu script yedi adayı
(chi3, beta, chi5, chi7 — 101c/101e/101f'nin dörtlüsü; chi5e, chi8e,
chi8o — 105b'nin üçlüsü) TAM SERTİFİKA BORU HATTINDAN geçirir. Boru hattı
105b'nin A/B/C/D adımlarıyla BİREBİR aynıdır; TEK DEĞİŞEN MOTORDUR.

  A) ince ızgara grid_frac = 0.03, T0 = 200          [118a'da yapıldı]
  B) sayım sertifikası d_i = i − (θ(z_i)−θ(z_0))/π; basamak bölgelerinde
     grid_frac = 0.005 kurtarma taraması (maks 3 tur)
  C) KALDIRILMIŞ DİP kurtarması (101e): işaretsiz |Z| minimumu +
     3 gerçek-Z (mpmath Hurwitz) + parabol kökleri
  D) TEMİZLİK (101f): pozitif basamakta sahte-kopya budama (mpmath
     orta-nokta işaret testi), negatif basamakta uyarlanır eşikli
     (0.45×medyan|Z|) ikinci kurtarma
  ÇIKTI: 118b_{ada}_zeros.npz (zeros + bolgeler + q + a + T1)

Ardından 107c'nin SERTİFİKA DEFTERİ birebir yeniden koşulur (K1 düzlük/
basamak, K2 kısa-çukur, K3 sıçrama + kopya; lo-kesimi, log-üçleme,
n≥3000 penceresi) ve 107c_sertifika.npz ile YAN YANA basılır.

ÖN-MÜHÜR (beklenti):
  B1  DİP-KURTARMA İHTİYACI ~0'A İNER. 101/105'te ada başına 128-435
      kayıp dip vardı ve HEPSİ mpmath'le kurtarıldı; 117'nin Δ'sında
      taper'la bu sayı 0 çıktı. GL(1)'de de C adımının katkısı
      onlarca-mertebeden birler-mertebesine düşmeli.
  B2  KALAN ŞÜPHELİ BÖLGE AZALIR (107c: ζ 0, chi8o 0, β 1, chi3 2,
      chi5e 2, chi8e 2, chi5 3, chi7 8 — toplam 18).
  B3  MASKELENEN % DÜŞER (107c: %3.7-10.0; χ₇ %10.0 ile en kötü).
  B4  SÜRÜKLENME ARALIĞI DARALIR.
  RİSK / DÜRÜST KAYIT: taper motorunun sıfır KONUMLARI keskininkinden
  ~1e-2⟨g⟩ mertebesinde farklıdır (117-S3'te medyan 5.8e-3). Yani yeni
  liste eskisinin "düzeltilmiş kopyası" değil, BAĞIMSIZ BİR DÖKÜMdür;
  118c'nin donma sayıları ondalıklarına kadar aynı çıkmak ZORUNDA
  DEĞİLDİR. Beklenen: hükümler aynı, sayılar bandın içinde oynar.
==========================================================================
SONUÇ (26 Ağustos, koşu 59 s — B1 ✓✓✓ B2 ✓✓✓ B3 ✓ B4 ✓✓)
==========================================================================
DÖKÜM (KESKİN → TAPER):
  ada   keskin ham  taper ham   Δham   B    C    D   taper n  keskin n
  chi3      79923      80199   +276   +4   +2   +0    80205     80211
  beta      74110      74420   +310   +0   +8   +0    74428     74430
  chi5      72514      72819   +305   +4   +8   +0    72831     72836
  chi7      64436      65048   +612   +0   +8   +0    65056     65055
  chi5e     72451      72823   +372   +0   +8   +0    72831     72829
  chi8e     65816      65938   +122   +4   +4   +0    65946     65948
  chi8o     65503      65939   +436   +2   +2   +0    65943     65943

B1 ✓✓✓ DİP-KURTARMA İHTİYACI ÇÖKTÜ. C adımı eskiden ada başına
  128-435 sıfır kurtarıyordu; şimdi 2-8. mpmath eval sayısı
  217-741 → 3-12 (60 kat azalma). Toplam sıfır sayıları eski
  dökümle 0-6 sıfır farkla örtüşüyor (chi8o birebir 65943) — iki
  BAĞIMSIZ yol aynı sayıya varıyor, sertifikanın en güçlü delili.
  117'nin Δ'da gördüğü "taper'la dip-kurtarma 0" olgusu GL(1)'de de
  (tam sıfır değil ama 60 kat azalmayla) doğrulandı.

B2 ✓✓✓ KALAN ŞÜPHELİ BÖLGE 18 → 0. Yedi adanın YEDİSİNDE DE nihai
  sertifikada kalan bölge SIFIR (eski: χ₇ 8, χ₅ 3, χ₃ 2, chi5e 2,
  chi8e 2, β 1, chi8o 0).
  ÜÇ KUSUR KAPISI DA SUSTU: K1 ham ihlal 85/254/789/169/84/85/0 → 0
  (yedisinde de), K2 ham 24/23/121/20/0/12/0 → 0, K3 yalnız χ₅'te 2.
  maks |med₂₀−med₂₀₀|: 1.85/2.01/1.97/1.90/0.37/1.40/0.20 →
  0.24/0.29/0.27/0.24/0.21/0.20/0.20 (eşik 0.7'nin çok altında,
  hepsi tertemiz tabanda).

B3 ✓ MASKELENEN ORAN düştü ama az: χ₇ %10.0 → %4.0 (en büyük kazanç),
  χ₃ %5.2 → %3.4, χ₈ᵉ %5.4 → %4.0, β %4.1 → %3.6; χ₅ %5.5 → %5.4,
  chi5e/chi8o değişmedi (%3.7 / %4.0). Kalan maskeleme artık
  KUSURDAN DEĞİL, n < 3000 düşen en alt log-penceresinden geliyor.
  Pencere sayısı da düştü (χ₇ 6 → 2, χ₅ 4 → 3, ötekiler 2): kusur
  kesikleri kalkınca segmentler birleşti.

B4 ✓✓ SÜRÜKLENME ARALIĞI: eski dökümlerde [-4.53,+0.83] / [-0.63,+4.96]
  gibi genişlikler vardı; yenide yedi adada da [-1.17, +2.37] içinde
  (saf S(t) salınımı; 65-80 bin sıfırda 4σ ≈ 2 normaldir).

ADLI KUSUR DENETİMİ — SÜRPRİZ DEĞİL AMA ÖLÇÜLDÜ: 104d/105e/107a'nın
  ELLE bulduğu YEDİ kusurun YEDİSİ DE artık analiz penceresinin
  İÇİNDE. Bu bir kapı kaçağı DEĞİL: aynı noktaların ±400 sıfırlık
  komşuluğunda ölçülen sürüklenme istatistikleri
     χ₇ 3081.67:  maks|Δmed₈₀| 0.49 → 0.05 ; maks|m₂₀−m₂₀₀| 1.97 → 0.13
     χ₃ 3872.54:  0.18 → 0.10 ; 1.85 → 0.15
     χ₅ 5820.00:  2.00 → 0.10 ; 2.01 → 0.19
     χ₈ᵉ 19352 :  0.13 → 0.09 ; 1.40 → 0.16
  yani KUSURUN KENDİSİ YOK OLDU. 104/105/107'nin üç seferde elle
  kovaladığı bütün kusurlar TEK BİR ALET DEĞİŞİKLİĞİNİN eseriymiş.
  χ₇'nin altı "kopya sıfırı" (u<0.005) ve chi8e'nin biri de gitti.

TEK YENİ KUSUR (dürüst kayıt): χ₅'te t = 34633.3216/34633.3246
  kopya çifti (u < 0.005). Eski dökümde orada 34633.3127/34633.3355
  vardı ve TAPER TARAMASI bu çifti bulamamıştı (118a'nın bulunamayan
  listesinde). Yani B adımının ince kurtarması (grid 0.005) çifti
  yakalarken çok yakın iki kök üretmiş. K3 kapısı onu yakaladı ve
  χ₅'i 2 yerine 3 pencereye böldü — analiz dışında.
==========================================================================
"""
import numpy as np
import mpmath as mp
import time
from pathlib import Path
from numpy.lib.stride_tricks import sliding_window_view

T00 = time.time()
HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
exec(open(HERE / "118a_taperli_L_motoru.py").read().split(
    'def kesin_Z')[0].split('if __name__')[0])
HERE = Path(__file__).resolve().parent
mp.mp.dps = 12


def gercek_Z(M, tab, t):
    s = mp.mpc(0.5, t)
    Lv = mp.mpc(0)
    for r in range(1, M.q):
        if tab[r % M.q] != 0:
            Lv += tab[r % M.q] * mp.zeta(s, mp.mpf(r) / M.q)
    Lv *= mp.power(M.q, -s)
    th = float(M.theta(np.array([t]))[0])
    return float((mp.e ** (1j * th) * Lv).real)


def basamak_bul(M, zz, win=80, esik=0.6):
    th = M.theta(zz)
    d = np.arange(len(zz)) - (th - th[0]) / np.pi
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


# ======================================================================
#  BORU HATTI  (105b'nin A/B/C/D adımları, TEK DEĞİŞEN MOTOR)
# ======================================================================
print("=" * 78, flush=True)
print("118b — YEDİ ADA, TAPER'LI MOTORLA YENİDEN DÖKÜM", flush=True)
print("=" * 78, flush=True)

ISTAT = {}
for etiket, q, tab, a, T1, _, _ in ADALAR:
    son = HERE / f"118b_{etiket}_zeros.npz"
    if son.exists():
        d = np.load(son)
        ISTAT[etiket] = dict(nham=int(d["nham"]), nB=int(d["nB"]),
                             nC=int(d["nC"]), nD=int(d["nD"]),
                             bolge=int(np.asarray(d["bolgeler"]
                                                  ).reshape(-1, 2).shape[0]),
                             sur=(float(d["surlo"]), float(d["surhi"])),
                             eval=int(d["neval"]), sure=float(d["sure"]))
        print(f"[{etiket}] nihai önbellek var, geçiliyor", flush=True)
        continue
    tA = time.time()
    M = LmotorT(q, tab, a, taper_c=0.5)
    print(f"\n===== {etiket} (q={q}, a={a}, T1={T1:.0f}) =====", flush=True)

    # ---------- A) ince ızgara (118a'da yapıldı, önbellekten) ----------
    ham = HERE / f"118b_{etiket}_ham.npz"
    if ham.exists():
        zz = np.load(ham)["zeros"]
        print(f"  A) taper'lı ince ızgara (118a önbelleği): {len(zz)} sıfır",
              flush=True)
    else:
        t0 = time.time()
        zz = M.sifir_bul(200.0, T1, grid_frac=0.03)
        np.savez(ham, zeros=zz, q=q, a=a, T1=T1, taper_c=0.5)
        print(f"  A) taper'lı ince ızgara (0.03): {len(zz)} sıfır, "
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
        print(f"     B tur {tur+1}: {len(bolgeler)} şüpheli bölge → "
              f"n={len(zz)}", flush=True)
    bolgeler, d = basamak_bul(M, zz)
    dd = d - np.median(d[:200])
    zc = zz.copy()                 # kopya-budamada referans
    n_B = len(zz)
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
    n_C = len(zz)
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
                    A_ = np.vstack([np.ones(5), np.array(z5t) - tstar,
                                    (np.array(z5t) - tstar) ** 2]).T
                    c0, c1, c2 = np.linalg.lstsq(A_, np.array(z5),
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
    sure = time.time() - tA
    print(f"  D) SON SERTİFİKA: n={len(zz)} ({len(zz)-n_d0:+d} temizlik); "
          f"KALAN BÖLGE={len(bolgeler)}; sürüklenme "
          f"[{dd.min():.2f}, {dd.max():.2f}]; L_eff≈{Leff:.2f}; "
          f"mpmath eval={n_eval}; {sure:.0f} sn", flush=True)
    np.savez(son, zeros=zz,
             bolgeler=np.array(bolgeler if bolgeler else []).reshape(-1, 2),
             q=q, a=a, T1=T1, taper_c=0.5, nham=n_ham, nB=n_B, nC=n_C,
             nD=len(zz), neval=n_eval, surlo=dd.min(), surhi=dd.max(),
             sure=sure)
    ISTAT[etiket] = dict(nham=n_ham, nB=n_B, nC=n_C, nD=len(zz),
                         bolge=len(bolgeler), sur=(dd.min(), dd.max()),
                         eval=n_eval, sure=sure)

# ======================================================================
#  DÖKÜM KARŞILAŞTIRMASI (105b/101 tablosu vs 118b)
# ======================================================================
ESKI_DOKUM = {   # (ham0.03, B kazancı, C kazancı, son n, kalan bölge)
    "chi3":  (79923, None, None, 80211, 2),
    "beta":  (None,  None, None, 74430, 1),
    "chi5":  (None,  None, None, 72836, 3),
    "chi7":  (None,  None, None, 65055, 8),
    "chi5e": (72451, 2, 376, 72829, 2),
    "chi8e": (65816, 2, 130, 65948, 2),
    "chi8o": (65503, 8, 432, 65943, 0),
}
print("\n" + "=" * 92, flush=True)
print("DÖKÜM KARŞILAŞTIRMASI — KESKİN (101c/101e/101f, 105b) vs TAPER (118b)",
      flush=True)
print("=" * 92, flush=True)
print(f"{'ada':>7} {'keskin ham':>10} {'taper ham':>10} {'Δham':>7} "
      f"{'B(+ince)':>9} {'C(+dip)':>8} {'D(temizlik)':>12} "
      f"{'taper n':>8} {'keskin n':>9} {'kalan böl.':>10} "
      f"{'sürüklenme':>18} {'eval':>6}", flush=True)
for etiket, q, tab, a, T1, ham_f, son_f in ADALAR:
    r = ISTAT[etiket]
    e = ESKI_DOKUM[etiket]
    hk = len(np.load(HERE / ham_f)["zeros"])
    print(f"{etiket:>7} {hk:>10} {r['nham']:>10} {r['nham']-hk:>+7} "
          f"{r['nB']-r['nham']:>+9} {r['nC']-r['nB']:>+8} "
          f"{r['nD']-r['nC']:>+12} {r['nD']:>8} {e[3]:>9} "
          f"{r['bolge']:>4} (eski {e[4]}) "
          f"[{r['sur'][0]:>+6.2f},{r['sur'][1]:>+6.2f}] {r['eval']:>6}",
          flush=True)


# ======================================================================
#  SERTİFİKA DEFTERİ (107c birebir) — YENİ SIFIRLARLA
# ======================================================================
def _roll(x, w, f):
    if len(x) <= w:
        return np.full(len(x), float(f(x)))
    m = f(sliding_window_view(x, w), axis=1)
    out = np.empty(len(x))
    out[w // 2:w // 2 + len(m)] = m
    out[:w // 2] = m[0]
    out[w // 2 + len(m):] = m[-1]
    return out


def segmentle_ix(zz, bad, min_n=3000):
    seg, kes, s0 = [], 0, 0
    for i in range(1, len(zz) + 1):
        if i == len(zz) or bad[i] != bad[i - 1]:
            if not bad[s0] and i - s0 >= min_n:
                seg.append((s0, i))
            if i < len(zz) and bad[i]:
                kes += 1
            s0 = i
    return seg, kes


def kapilari_kos(zc, sayim, q):
    """107c'nin kapı zinciri: K1 (basamak) + K2 (kısa çukur) → log-üçleme
    → K3 (sıçrama + kopya), n ≥ 3000 pencereleri."""
    n = len(zc)
    d = np.arange(n) - (sayim(zc) - sayim(zc[0]))
    med = np.median(sliding_window_view(d, 80), axis=1)
    adim = med[81:] - med[:-81]
    ham1 = np.zeros(n, bool)
    ham1[np.clip(np.where(np.abs(adim) > 0.5)[0] + 40, 0, n - 1)] = True
    g1 = np.zeros(n, bool)
    for i in np.where(ham1)[0]:
        g1[max(0, i - 160):i + 320] = True
    sap = np.abs(_roll(d, 20, np.median) - _roll(d, 200, np.median))
    ham2 = sap > 0.7
    g2 = np.zeros(n, bool)
    for i in np.where(ham2)[0]:
        g2[max(0, i - 160):i + 160] = True
    seg, kes12 = segmentle_ix(zc, g1 | g2)
    kenar = np.exp(np.linspace(np.log(zc[0]), np.log(zc[-1] * 1.0001), 4))
    parcalar = []
    for (x, y) in seg:
        s = zc[x:y]
        for i in range(3):
            m = (s >= kenar[i]) & (s < kenar[i + 1])
            if m.sum() >= 3000:
                j = np.where(m)[0]
                parcalar.append((x + j[0], x + j[-1] + 1))
    son, ham3, dup3, rmax = [], 0, [], []
    for (x, y) in parcalar:
        p = zc[x:y]
        dd = np.arange(len(p)) - (sayim(p) - sayim(p[0]))
        r = _roll(dd, 21, np.mean) - _roll(dd, 801, np.median)
        gp = np.diff(p)
        mid = 0.5 * (p[:-1] + p[1:])
        u = gp * np.log(q * mid / TWO_PI) / TWO_PI
        bad = np.abs(r) > 0.30
        for i in np.where(u < 0.005)[0]:
            bad[i] = True
            bad[min(i + 1, len(bad) - 1)] = True
            dup3.append(float(p[i]))
        ham3 += int(bad.sum())
        rmax.append(float(np.abs(r).max()))
        kes = np.zeros(len(p), bool)
        for i in np.where(bad)[0]:
            kes[max(0, i - 600):min(len(p), i + 601)] = True
        s0 = 0
        for i in range(1, len(p) + 1):
            if i == len(p) or kes[i] != kes[i - 1]:
                if not kes[s0] and i - s0 >= 3000:
                    son.append((x + s0, x + i))
                s0 = i
    gp = np.diff(zc)
    mid = 0.5 * (zc[:-1] + zc[1:])
    u_all = gp * np.log(q * mid / TWO_PI) / TWO_PI
    kop_all = [float(zc[i]) for i in np.where(u_all < 0.005)[0]]
    return dict(n=n, ham1=int(ham1.sum()), ham2=int(ham2.sum()), ham3=ham3,
                sapmax=float(sap.max()), rmax=rmax, dup_pencere=dup3,
                dup_ada=kop_all, son=son,
                nson=sum(b - a for (a, b) in son), npen=len(son))


ADLI = [("chi7", 3081.6688, "104d kopya çifti #1"),
        ("chi7", 3081.8885, "104d kopya çifti #2"),
        ("chi7", 3095.0, "104d |d|>1 bloğunun ortası"),
        ("chi3", 3872.54, "104d χ₃ w1 kusuru (baş)"),
        ("chi3", 3890.24, "104d χ₃ w1 kusuru (son)"),
        ("chi5", 5820.0, "107a χ₅ kısa çukuru"),
        ("chi8e", 19352.0, "105e kısa çukuru")]

print("\n" + "=" * 118, flush=True)
print("SERTİFİKA DEFTERİ — 107c (KESKİN) vs 118b (TAPER)", flush=True)
print("=" * 118, flush=True)
E = np.load(HERE / "107c_sertifika.npz", allow_pickle=True)
DEF, PENCERE_T = {}, {}
for etiket, q, tab, a, T1, _, _ in ADALAR:
    dd = np.load(HERE / f"118b_{etiket}_zeros.npz")
    zraw = dd["zeros"]
    bolge = int(np.asarray(dd["bolgeler"]).reshape(-1, 2).shape[0])
    lo = np.exp(np.log(zraw[0] + 1) + 0.25 *
                (np.log(zraw[-1]) - np.log(zraw[0] + 1)))
    zc = zraw[zraw >= lo]
    M = LmotorT(q, tab, a, 0.5)
    say = (lambda t, M=M: M.theta(t) / np.pi)
    r = kapilari_kos(zc, say, float(q))
    r["nham"] = len(zraw)
    r["bolge"] = bolge
    r["qilk"] = 3 if etiket in ("beta", "chi8e", "chi8o") else 2
    Ls, Ns = [], []
    PENCERE_T[etiket] = []
    for (x, y) in r["son"]:
        p = zc[x:y]
        PENCERE_T[etiket].append((float(p[0]), float(p[-1])))
        mids = 0.5 * (p[:-1] + p[1:])
        Ls.append(float(np.log(q * mids / TWO_PI).mean()))
        Ns.append(len(p))
    r["Lort"] = float(np.average(Ls, weights=np.array(Ns, float)))
    DEF[etiket] = r

SIRA = ["chi3", "chi5", "chi7", "chi5e", "beta", "chi8e", "chi8o"]
print(f"{'küme':>6} {'':>4} {'n_ham':>7} {'n_analiz':>8} {'kalan':>6} "
      f"{'K1':>5} {'K2':>5} {'K3':>5} {'kopya':>6} {'maskelenen':>10} "
      f"{'%':>5} {'n_pen':>5} {'⟨L⟩':>6} {'τ_ilk':>7}", flush=True)
for a in SIRA:
    r = DEF[a]
    mask = r["n"] - r["nson"]
    tilk = np.log(r["qilk"]) / r["Lort"]
    ek = float(E[f"Lort_{a}"])
    print(f"{a:>6} {'ESKİ':>4} {int(E[f'nham_{a}']):>7} "
          f"{int(E[f'n_{a}']):>8} {int(E[f'bolge_{a}']):>6} "
          f"{int(E[f'ham1_{a}']):>5} {int(E[f'ham2_{a}']):>5} "
          f"{int(E[f'ham3_{a}']):>5} {len(E[f'dupada_{a}']):>6} "
          f"{int(E[f'n_{a}'])-int(E[f'nson_{a}']):>10} "
          f"{100*(int(E[f'n_{a}'])-int(E[f'nson_{a}']))/int(E[f'n_{a}']):>5.1f} "
          f"{int(E[f'npen_{a}']):>5} {ek:>6.2f} "
          f"{np.log(r['qilk'])/ek:>7.4f}", flush=True)
    print(f"{'':>6} {'YENİ':>4} {r['nham']:>7} {r['n']:>8} {r['bolge']:>6} "
          f"{r['ham1']:>5} {r['ham2']:>5} {r['ham3']:>5} "
          f"{len(r['dup_ada']):>6} {mask:>10} {100*mask/r['n']:>5.1f} "
          f"{r['npen']:>5} {r['Lort']:>6.2f} {tilk:>7.4f}", flush=True)

print(f"\n{'küme':>6} {'maks|med₂₀−med₂₀₀| eski→yeni':>30} "
      f"{'K3 maks|r| pencere başına (yeni)':>34}", flush=True)
for a in SIRA:
    r = DEF[a]
    print(f"{a:>6} {float(E[f'sapmax_{a}']):>13.2f} → {r['sapmax']:>6.2f}"
          f"          " + " ".join(f"{x:.2f}" for x in r["rmax"]), flush=True)

print("\n--- ADLI KUSUR DENETİMİ: ESKİ DÖKÜMÜN KUSURLARI YENİ VERİDE VAR MI? ---",
      flush=True)
print("  Kusur ARTIK YOKSA kapı işaret vermez ve nokta analize GİRER. Bu bir",
      flush=True)
print("  kapı kaçağı DEĞİL, kusurun ortadan kalkmasıdır — ama iddia edilmez,",
      flush=True)
print("  ÖLÇÜLÜR: her noktanın ±400 sıfırlık komşuluğunda sayım-sürüklenmesi",
      flush=True)
print("  istatistikleri eski ve yeni dökümde yan yana verilir.", flush=True)
print(f"  {'küme':>6} {'t':>11} {'ESKİ maks|Δmed₈₀|':>17} "
      f"{'ESKİ maks|m₂₀−m₂₀₀|':>19} {'YENİ maks|Δmed₈₀|':>17} "
      f"{'YENİ maks|m₂₀−m₂₀₀|':>19}  hüküm", flush=True)
ESKI_DOSYA = {e[0]: e[6] for e in ADALAR}


def yerel_kusur(zc, sayim, t, pad=400):
    i = int(np.searchsorted(zc, t))
    a, b = max(0, i - pad), min(len(zc), i + pad)
    p = zc[a:b]
    if len(p) < 250:
        return np.nan, np.nan
    d = np.arange(len(p)) - (sayim(p) - sayim(p[0]))
    med = np.median(sliding_window_view(d, 80), axis=1)
    adim = np.abs(med[81:] - med[:-81]).max() if len(med) > 82 else np.nan
    sap = np.abs(_roll(d, 20, np.median) - _roll(d, 200, np.median)).max()
    return float(adim), float(sap)


for (et, t, ad) in ADLI:
    ic = any(x <= t <= y for (x, y) in PENCERE_T[et])
    q = [e[1] for e in ADALAR if e[0] == et][0]
    tab = [e[2] for e in ADALAR if e[0] == et][0]
    aa = [e[3] for e in ADALAR if e[0] == et][0]
    Mn = LmotorT(q, tab, aa, 0.5)
    Me = Lmotor(q, tab, aa)
    sy_n = (lambda x, M=Mn: M.theta(x) / np.pi)
    sy_e = (lambda x, M=Me: M.theta(x) / np.pi)
    z_es = np.load(HERE / ESKI_DOSYA[et])["zeros"]
    z_ye = np.load(HERE / f"118b_{et}_zeros.npz")["zeros"]
    ae, se = yerel_kusur(z_es, sy_e, t)
    an, sn = yerel_kusur(z_ye, sy_n, t)
    hk = ("KUSUR YOK (analize girdi) ✓" if (an < 0.5 and sn < 0.7)
          else ("kusur SÜRÜYOR ama kapı kaçırdı ✗" if ic
                else "kusur sürüyor, kapı yakaladı"))
    print(f"  {et:>6} {t:>11.4f} {ae:>17.2f} {se:>19.2f} "
          f"{an:>17.2f} {sn:>19.2f}  {hk}", flush=True)

print("\n--- NİHAİ ANALİZ PENCERELERİ (yeni) ---", flush=True)
for a in SIRA:
    print(f"  [{a}] {DEF[a]['npen']} pencere: " + "  ".join(
        f"[{x:.0f},{y:.0f}]" for (x, y) in PENCERE_T[a]), flush=True)

print("\n--- KOPYA SIFIR ENVANTERİ (u < 0.005), yeni döküm ---", flush=True)
for a in SIRA:
    ka = DEF[a]["dup_ada"]
    print(f"  [{a}] " + (f"{len(ka)} kopya: "
                         + " ".join(f"{t:.4f}" for t in ka)
                         if ka else "kopya YOK")
          + f"   (eski: {len(E[f'dupada_{a}'])})", flush=True)

np.savez(HERE / "118b_sertifika.npz", adlar=np.array(SIRA),
         **{f"{k}_{a}": DEF[a][k] for a in SIRA
            for k in ["nham", "n", "bolge", "ham1", "ham2", "ham3",
                      "nson", "npen", "sapmax", "Lort", "qilk"]},
         **{f"dupada_{a}": np.array(DEF[a]["dup_ada"]) for a in SIRA},
         **{f"pen_{a}": np.array(PENCERE_T[a]).reshape(-1, 2) for a in SIRA})
print(f"\n118b_sertifika.npz yazıldı. Toplam süre {time.time()-T00:.0f} s.",
      flush=True)
