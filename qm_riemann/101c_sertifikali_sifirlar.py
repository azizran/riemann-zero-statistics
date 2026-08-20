"""
101c — SERTİFİKALI SIFIRLAR: SAYIM-DENETİMLİ YENİDEN TÜRETİM (20 Ağustos)
==========================================================================
101b kapısı: %0.1 kaçık bile donma ölçümünü 0.09→0.55'e bozuyor.
Motorun grid_frac=0.12 taraması ~%0.3 kaçırıyordu → dört adayı
grid_frac=0.03 ile yeniden türet + SAYIM SERTİFİKASI:
   d_i = i − (θ(z_i) − θ(z_0))/π  sürüklenmesi izlenir (S(t) ±O(1)
   salınır; kaçık sıfır kalıcı −1 basamağı bırakır). İleri/geri medyan
   (pencere 80) farkı > 0.6 olan bölgelere grid_frac=0.005 kurtarma
   taraması; birleştir, tekilleştir, yeniden denetle (maks 3 tur).
Çıktı: 101c_{ada}_zeros.npz + sertifika raporu (eklenen sıfır sayısı,
kalan basamak sayısı — SIFIR olmalı).
"""

import numpy as np
import time
from pathlib import Path

exec(open("98_L_motoru.py").read().split('CHI4 = ')[0])

HERE = Path(".").resolve()
CHI4 = {0: 0, 1: 1, 2: 0, 3: -1}
CHI3 = {0: 0, 1: 1, 2: -1}
CHI5 = {0: 0, 1: 1, 2: 1j, 3: -1j, 4: -1}
z6 = np.exp(1j * np.pi / 3)
CHI7 = {0: 0, 1: 1, 3: z6, 2: z6**2, 6: z6**3, 4: z6**4, 5: z6**5}

ISLANDS = [("chi3", 3, CHI3, 55000.0), ("beta", 4, CHI4, 50000.0),
           ("chi5", 5, CHI5, 48000.0), ("chi7", 7, CHI7, 42000.0)]

def basamak_bul(M, zz, win=80, esik=0.6):
    """Sayım sürüklenmesinde basamakları bul → şüpheli t-bölgeleri."""
    th = M.theta(zz)
    d = np.arange(len(zz)) - (th - th[0]) / np.pi
    # ileri/geri medyan farkı
    from numpy.lib.stride_tricks import sliding_window_view
    if len(d) < 2 * win + 1:
        return [], d
    sw = sliding_window_view(d, win)
    med = np.median(sw, axis=1)          # med[i] = medyan d[i:i+win]
    # basamak i civarında: med[i+1] − med[i−win] farkı
    adim = med[win + 1:] - med[:-(win + 1)]
    kotu = np.where(np.abs(adim) > esik)[0] + win // 2
    bolgeler = []
    for i in kotu:
        t0 = zz[max(0, i - win // 2)]
        t1 = zz[min(len(zz) - 1, i + 3 * win // 2)]
        if bolgeler and t0 <= bolgeler[-1][1]:
            bolgeler[-1] = (bolgeler[-1][0], max(bolgeler[-1][1], t1))
        else:
            bolgeler.append((t0, t1))
    return bolgeler, d

for etiket, q, tab, T1 in ISLANDS:
    out = HERE / f"101c_{etiket}_zeros.npz"
    if out.exists():
        print(f"[{etiket}] önbellek var, geçiliyor"); continue
    M = Lmotor(q, tab, 1)
    t0 = time.time()
    zz = M.sifir_bul(200.0, T1, grid_frac=0.03)
    n0 = len(zz)
    eski = np.load(HERE / f"99_{etiket}_zeros.npz")["zeros"]
    print(f"[{etiket}] ince ızgara: {n0} sıfır ({n0 - len(eski):+d} vs 99; "
          f"{time.time()-t0:.0f} sn)")
    for tur in range(3):
        bolgeler, d = basamak_bul(M, zz)
        if not bolgeler:
            break
        print(f"  tur {tur+1}: {len(bolgeler)} şüpheli bölge, kurtarma...")
        yeni = [zz]
        for (a, b) in bolgeler:
            yeni.append(M.sifir_bul(a - 0.5, b + 0.5, grid_frac=0.005))
        zz = np.unique(np.concatenate(yeni))
        # tekilleştir: 1e-3'ten yakınları birleştir
        keep = np.concatenate([[True], np.diff(zz) > 1e-3])
        zz = zz[keep]
    bolgeler, d = basamak_bul(M, zz)
    dd = d - np.median(d[:200])
    print(f"  SERTİFİKA: n={len(zz)} ({len(zz)-n0:+d} kurtarılan); "
          f"kalan basamak bölgesi = {len(bolgeler)}; "
          f"sürüklenme aralığı [{dd.min():.2f}, {dd.max():.2f}] "
          f"(S-salınımı beklenir, basamaksız)")
    np.savez(out, zeros=zz)
print("BİTTİ")
