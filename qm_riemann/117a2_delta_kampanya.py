"""
117a2 — Δ KAMPANYASI + SERTİFİKA BORU HATTI (26 Ağustos 2026)
==========================================================================
ÖN-MÜHÜR (ölçümden ÖNCE)
==========================================================================
117a'nın kapılarından geçen GL(2) motoruyla Ramanujan Δ'nın sıfır
örgüsü dökülür. 101c/101e/101f + 104e + 105e boru hattının TAMAMI
uygulanır — ama BİR FARKLA, ve bu fark dürüstçe kayda geçer:

  GL(1) kampanyalarında dip-kurtarmanın hakemi mpmath'in Hurwitz-ζ
  tabanlı GERÇEK Z'siydi (motordan BAĞIMSIZ bir oracle). Δ için böyle
  bir oracle YOK: 117a-G3'ün kesin Hecke integral formülü e^{πt/2}
  iptali yüzünden t ≲ 250'de kalıyor, kampanya ise t ~ 2e4'e gidiyor.
  Bu yüzden dip-kurtarma hakemi UZLAŞI olur: aynı ana toplamın FARKLI
  kesim profilleriyle (taper c = 0.35 / 0.50 / 0.80 ve keskin) bağımsız
  değerlendirmeleri. Bunlar bağımsız FORMÜLLER değildir — yalnız hata
  terimleri farklıdır (117a-G3'te keskin %4.4, c=0.5 %0.31 hata verdi).
  ⟹ Δ'nın dip-kurtarması GL(1)'inkinden ZAYIF bir sertifikadır;
  kurtarılan çiftler ayrıca sayılır ve analiz onlarsız da tekrarlanır.

ADIMLAR
  A) ince ızgara grid_frac=0.03, T0=200 → 117a_delta_ham.npz
  B) SAYIM SERTİFİKASI d_i = i − (θ(z_i)−θ(z_0))/π; basamak bölgelerinde
     grid_frac=0.004 ile kurtarma (maks 3 tur)
  C) KALDIRILMIŞ DİP kurtarması — UZLAŞI hakemiyle (yukarıdaki not)
  D) TEMİZLİK: pozitif basamakta sahte-kopya budama, negatif basamakta
     uyarlanır eşikli ikinci kurtarma
  ÇIKTI: 117a_delta_zeros.npz (zeros, bolgeler, kurtarilan)

ÖNGÖRÜ (mühür): taper motorunun hatası GL(1) motorundan ~14 kat küçük
olduğu için KALDIRILMIŞ DİP sayısı GL(1) adalarındakinden (ada başına
130-430) BELİRGİN AZ olmalı. Çıkmazsa: derece-2'nin daha uzun toplamı
hatayı geri getiriyordur — o da bir bulgudur.

==========================================================================
SONUÇ (26 Ağustos — T1 = 36000, toplam 3.7 dk)
==========================================================================
adım                       n        değişim   kalan bölge   sürüklenme
A) ince ızgara 0.03     87541          —            —            —
B) sayım sertifikası    87543         +2            1       [−3.15, +1.39]
C) dip kurtarma         87543         +0            1       [−3.15, +1.39]
D) temizlik             87543         +0            1       [−3.15, +1.39]
L_eff ≈ 16.05;  A adımı 193 sn.

ÖN-MÜHÜRLÜ ÖNGÖRÜ İSABET ETTİ — VE BEKLENENDEN ÇOK DAHA GÜÇLÜ:
GL(1) kampanyalarında sayım sertifikası ada başına 128-435 KAYIP sıfır
(kaldırılmış dip) gösteriyordu; burada 87.5k sıfırda toplam sürüklenme
yayılımı 4.5 basamak ve UZLAŞAN TEK BİR DİP YOK. Sebep 117a-G3'te
ölçüldü: taper (c=0.5) kesimi motorun hatasını keskin kesime göre
11-14 kat düşürüyor ⟹ dar çiftlerin |Z| çukuru artık sıfır üstüne
kaldırılmıyor. Yani "kaldırılmış dip" GL(1)'de de ALETSELDİ, derecenin
değil KESİM PROFİLİNİN sonucuymuş. (GL(1) adalarının taper'lı yeniden
dökümü — 101/105 verilerinin bağımsız denetimi — AÇIK İŞ olarak kayda
geçiyor; oradaki mpmath kurtarmaları muhtemelen gerçek sıfırları
buluyordu ama gerek olmayabilirdi.)

BAĞIMSIZ SAĞLAMA (kuram): N(T) = (T/π)·log(T/2πe) + O(log T) ⟹
N(36000) ≈ 87 700; ölçülen 87 543 (%0.2 içinde). Bu, hem derece-2
yoğunluğunu hem de L_eff = 2·log(t/2π) konvansiyonunu doğruluyor
(kaptanın mühründeki (2/π)·log yoğunluğu 2 kat fazlaydı).

DÜRÜST SINIR: kalan 1 kusurlu bölge analizden KESİLİR (düzlük
segmentasyonu, 117b). Dip-kurtarma hakemi bu seferde hiç devreye
girmediği için "uzlaşı hakeminin zayıflığı" pratikte test EDİLMEDİ —
Δ'da gerekmedi; başka bir GL(2) formunda gerekirse ilk sınanacak yer.
"""

import numpy as np
import time
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
exec(open(HERE / "117a_delta_motoru.py").read().split('if __name__')[0])

T0 = 200.0
T1 = float(sys.argv[1]) if len(sys.argv) > 1 else 20000.0

print("=" * 74, flush=True)
print(f"117a2 — Δ KAMPANYASI  T ∈ [{T0:.0f}, {T1:.0f}]", flush=True)
print("=" * 74, flush=True)

t00 = time.time()
TAU = tau_tablosu(NMAX_TAU)
M = DeltaMotor(TAU, taper_c=0.5)                      # ana motor
VAR = [DeltaMotor(TAU, taper_c=c) for c in (0.35, 0.8)] + \
      [DeltaMotor(TAU, taper_c=0.0)]                  # uzlaşı hakemleri
print(f"τ tablosu hazır ({time.time()-t00:.1f} sn); "
      f"X(T1) = {T1/TWO_PI:.0f} ≤ {NMAX_TAU} ✓", flush=True)


def basamak_bul(M, zz, win=80, esik=0.6):
    d = M.sayim_sapmasi(zz)
    from numpy.lib.stride_tricks import sliding_window_view
    med = np.median(sliding_window_view(d, win), axis=1)
    adim = med[win + 1:] - med[:-(win + 1)]
    kotu = np.where(np.abs(adim) > esik)[0] + win // 2
    bolgeler = []
    for i in kotu:
        a = zz[max(0, i - win)]
        b = zz[min(len(zz) - 1, i + 2 * win)]
        if bolgeler and a <= bolgeler[-1][1]:
            bolgeler[-1] = (bolgeler[-1][0], max(bolgeler[-1][1], b))
        else:
            bolgeler.append((a, b))
    return bolgeler, d


def yerel_adim(zz, d, a, b):
    ia, ib = np.searchsorted(zz, a), np.searchsorted(zz, b)
    sol = np.median(d[max(0, ia - 120):max(1, ia - 10)])
    sag = np.median(d[min(len(d) - 1, ib + 10):min(len(d), ib + 120)])
    return sag - sol


def parabol_kok(mot, tstar, h):
    """3 nokta parabolü → çift kök (varsa)."""
    z3 = mot.Z(np.array([tstar - h, tstar, tstar + h]))
    c2 = (z3[0] - 2 * z3[1] + z3[2]) / (2 * h * h)
    c1 = (z3[2] - z3[0]) / (2 * h)
    c0 = z3[1]
    disc = c1 * c1 - 4 * c2 * c0
    if abs(c2) < 1e-14 or disc <= 0:
        return None
    r1 = tstar + (-c1 - np.sqrt(disc)) / (2 * c2)
    r2 = tstar + (-c1 + np.sqrt(disc)) / (2 * c2)
    return (min(r1, r2), max(r1, r2))


def uzlasi_cift(tstar, h, gbar):
    """UZLAŞI HAKEMİ: ana motor + 3 varyant. Hepsi çift kök vermeli ve
    kökler 0.25⟨g⟩ içinde uyuşmalı. Ana motorun kökleri döner."""
    ana = parabol_kok(M, tstar, h)
    if ana is None or abs(ana[1] - ana[0]) > 0.8 * gbar:
        return None
    uy = 0
    for mv in VAR:
        r = parabol_kok(mv, tstar, h)
        if r is not None and abs(r[0] - ana[0]) < 0.25 * gbar \
                and abs(r[1] - ana[1]) < 0.25 * gbar:
            uy += 1
    return ana if uy >= 2 else None


# ---------------- A) ince ızgara ----------------
ham = HERE / "117a_delta_ham.npz"
if ham.exists():
    zz = np.load(ham)["zeros"]
    print(f"A) ham önbellek: {len(zz)} sıfır", flush=True)
else:
    t0 = time.time()
    zz = M.sifir_bul(T0, T1, grid_frac=0.03)
    np.savez(ham, zeros=zz, T0=T0, T1=T1)
    print(f"A) ince ızgara (0.03): {len(zz)} sıfır, {time.time()-t0:.0f} sn",
          flush=True)
n_ham = len(zz)

# ---------------- B) sayım sertifikası ----------------
t0 = time.time()
for tur in range(3):
    bolgeler, d = basamak_bul(M, zz)
    if not bolgeler:
        break
    yeni = [zz]
    for (x, y) in bolgeler:
        yeni.append(M.sifir_bul(x - 0.5, y + 0.5, grid_frac=0.004))
    zz = np.unique(np.concatenate(yeni))
    zz = zz[np.concatenate([[True], np.diff(zz) > 1e-3])]
    print(f"   B tur {tur+1}: {len(bolgeler)} şüpheli bölge → n={len(zz)}",
          flush=True)
bolgeler, d = basamak_bul(M, zz)
dd = d - np.median(d[:200])
zc = zz.copy()
print(f"B) SERTİFİKA-1: n={len(zz)} ({len(zz)-n_ham:+d}); "
      f"kalan bölge={len(bolgeler)}; sürüklenme [{dd.min():.2f},{dd.max():.2f}]"
      f"  ({time.time()-t0:.0f} sn)", flush=True)

# ---------------- C) kaldırılmış dip kurtarması (uzlaşı hakemi) --------
t0 = time.time()
n_c0 = len(zz)
kurtarilan = []
for tur in range(4):
    bolgeler, d = basamak_bul(M, zz)
    if not bolgeler:
        break
    yeni = []
    for (x, y) in bolgeler:
        gbar = np.pi / np.log(0.5 * (x + y) / TWO_PI)
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
            r = uzlasi_cift(ts[i], 0.10 * gbar, gbar)
            if r is not None:
                yeni += [r[0], r[1]]
    if not yeni:
        print(f"   C tur {tur+1}: uzlaşan dip yok, kalan bölge {len(bolgeler)}",
              flush=True)
        break
    kurtarilan += yeni
    zz = np.unique(np.concatenate([zz, np.array(yeni)]))
    zz = zz[np.concatenate([[True], np.diff(zz) > 5e-4])]
    print(f"   C tur {tur+1}: {len(bolgeler)} bölge, +{len(yeni)} sıfır "
          f"({time.time()-t0:.0f} sn)", flush=True)
bolgeler, d = basamak_bul(M, zz)
dd = d - np.median(d[:200])
print(f"C) SERTİFİKA-2: n={len(zz)} ({len(zz)-n_c0:+d} dip-kurtarma); "
      f"kalan bölge={len(bolgeler)}; sürüklenme [{dd.min():.2f},{dd.max():.2f}]",
      flush=True)

# ---------------- D) temizlik ----------------
t0 = time.time()
n_d0 = len(zz)
for tur in range(4):
    bolgeler, d = basamak_bul(M, zz)
    if not bolgeler:
        break
    degisti = False
    sil, ekle = [], []
    for (x, y) in bolgeler:
        adim = yerel_adim(zz, d, x, y)
        gbar = np.pi / np.log(0.5 * (x + y) / TWO_PI)
        if adim > 0.5:                      # fazla sıfır → sahte kopya
            zw = zz[(zz >= x) & (zz <= y)]
            j = np.searchsorted(zc, zw)
            dist = np.minimum(np.abs(zw - zc[np.clip(j, 0, len(zc) - 1)]),
                              np.abs(zw - zc[np.clip(j - 1, 0, len(zc) - 1)]))
            ins = zw[dist > 1e-6]
            if len(ins) == 0:
                continue
            for c in np.split(ins, np.where(np.diff(ins) > 0.5 * gbar)[0] + 1):
                if len(c) < 2:
                    sil += list(c); degisti = True; continue
                o = len(c) // 2
                cift = c[max(0, o - 1):max(0, o - 1) + 2]
                tl, th_ = cift[0], cift[-1]
                sh = np.sign(M.Z(np.array([tl - 0.3 * gbar, th_ + 0.3 * gbar])))
                zmid = M.Z(np.array([0.5 * (tl + th_)]))[0]
                sil += [w for w in c if w not in cift]
                if sh[0] == sh[1] and np.sign(zmid) == sh[0]:
                    sil += list(cift)
                degisti = True
        elif adim < -0.5:                   # eksik sıfır → ikinci kurtarma
            ts = np.arange(x, y, 0.002 * gbar)
            v = M.Z(ts); av = np.abs(v)
            esik = 0.45 * np.median(av)
            imin = np.where((av[1:-1] < av[:-2]) & (av[1:-1] <= av[2:]))[0] + 1

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
                r = uzlasi_cift(ts[i], 0.10 * gbar, gbar)
                if r is not None and abs(r[1] - r[0]) > 1e-4:
                    ekle += [r[0], r[1]]; degisti = True
    if not degisti:
        break
    if sil:
        zz = zz[~np.isin(zz, np.array(sil))]
    if ekle:
        kurtarilan += ekle
        zz = np.unique(np.concatenate([zz, np.array(ekle)]))
        zz = zz[np.concatenate([[True], np.diff(zz) > 5e-4])]
    print(f"   D tur {tur+1}: −{len(sil)} +{len(ekle)} "
          f"({time.time()-t0:.0f} sn)", flush=True)

bolgeler, d = basamak_bul(M, zz)
dd = d - np.median(d[:200])
Leff = float(2 * np.log(zz.mean() / TWO_PI))
print(f"D) SON SERTİFİKA: n={len(zz)} ({len(zz)-n_d0:+d} temizlik); "
      f"KALAN BÖLGE={len(bolgeler)}; sürüklenme "
      f"[{dd.min():.2f},{dd.max():.2f}]; L_eff≈{Leff:.2f}", flush=True)
print(f"   kurtarılan (uzlaşı hakemli, zayıf sertifika) toplam: "
      f"{len(kurtarilan)} sıfır", flush=True)
print(f"   TOPLAM SÜRE: {(time.time()-t00)/60:.1f} dk", flush=True)

np.savez(HERE / "117a_delta_zeros.npz", zeros=zz,
         bolgeler=np.array(bolgeler if bolgeler else []).reshape(-1, 2),
         kurtarilan=np.array(sorted(kurtarilan)), T0=T0, T1=T1)
print("\n117a_delta_zeros.npz yazıldı.", flush=True)
