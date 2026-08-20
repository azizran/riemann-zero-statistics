"""
101d — DONMA SINIRI ADALARDA, SERTİFİKALI SIFIRLARLA (20 Ağustos)
==========================================================================
101'in ön-mühürlü öngörüleri (P1-P3) burada test edilir; 101 ham hâliyle
ARTEFAKTTI (101b: kusur kapısı), veri 101c sertifikalı sıfırlar.

ÖN-MÜHÜR (101'den devir + keskinleştirme):
  P1  χ₃, χ₅, χ₇ çözülme eğrisi ζ'nınkini izler: τ=0.04-0.05'te donuk
      (<0.15), ilk çizgi log2/L ≈ 0.068-0.073'te yükselişte, τ≥0.3 CUE.
  P2  AYRIŞTIRICI — β (2-çizgisi yok, ilk çizgi log3/L ≈ 0.113):
      çözülme GECİKİR; τ=0.08-0.10'da ötekiler ~0.45-0.65 iken β
      belirgin düşük (≲0.3); yarı-çözülme ~×1.6 kaymalı.
  P3  τ≥0.3'te dört ada + ζ tek CUE eğrisinde.

Makine: 93 estimatörü — yerel katlama + kübik detrend, ω-bandı
τ0·Lw ± 0.02Lw (160 frek, çizgi-dışı >0.01), segmentler-arası pay
havuzu, 40-permütasyonlu vekil taban. ζ referansı aynı estimatörle.

VERİ DİSİPLİNİ — DÜZLÜK SEGMENTASYONU (101b kusur kapısının gereği):
donma estimatörü kayıp/sahte sıfıra aşırı duyarlı (101b: %0.1 → 0.55).
Küresel tamlık yerine SERTİFİKALI-DÜZ segmentler: sayım sürüklenmesi
d_i = i − Δθ/π'nin 80'lik medyanı bir segment içinde |Δ|<0.5 düz
kalmalı; her basamakta kes, kesim çevresini (±160 sıfır) at. Kusurlar
segment sınırlarında dışarıda kalır — analize hiçbiri giremez.
ζ'ya da aynı disiplin uygulanır (pürüzsüz sayımla).

SONUÇ (20 Ağustos — P1, P2, P3 ÜÇÜ DE DOĞRULANDI; figür 101g):
  τ=0.068 (log2/L): ζ 0.341, χ₃ 0.394, χ₅ 0.350, χ₇ 0.465 — dördü
  çözülmede. β: 0.040 (vekil taban 0.028!) — HÂLÂ DONUK. β 0.085-0.113
  arasında dik çözülür: 0.116 → 0.460 → 0.705 tam log3/L=0.113'te;
  0.14'te ortak eğriye biner. τ=0.3'te beş aile 1.007-1.055 (tek CUE).
  Yarı-çözülme kayması ~×1.4-1.5 (öngörü log3/log2 = 1.585).
  HÜKÜM: donma sınırı adanın İLK SAĞ KALAN ÇİZGİSİDİR — ölü çizgi
  (χ₄(2)=0) çözülmeyi log3/L'e erteler. Kırınımdaki ölü benek
  dinamikte de ölüdür: benek yoksa çözülme yok. "Gazı çözen
  aritmetiktir" iddiasının kontrollü (karakter-programlı) deneyi.
"""

import numpy as np
from pathlib import Path

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
rng = np.random.default_rng(1010)

PKS = [2,3,4,5,7,8,9,11,13,16,17,19,23,25,27,29,31,32,37,41,43,47,49,53,
       59,61,64,67,71,73,79,81,83,89,97,101,103,107,109,113,121,125,127,128]
LINES_ALL = [np.log(qq) for qq in PKS]

def pencere_hazirla(z, qeff):
    gaps = np.diff(z)
    mids = 0.5 * (z[:-1] + z[1:])
    Lw = float(np.log(qeff * mids / TWO_PI).mean())
    ds = gaps * np.log(qeff * mids / TWO_PI) / TWO_PI - 1
    tt = (mids - mids.mean()) / (mids[-1] - mids[0])
    P = np.vstack([np.ones_like(tt), tt, tt**2, tt**3]).T
    ds = ds - P @ np.linalg.lstsq(P, ds, rcond=None)[0]
    return (z, mids, ds, Lw)

def egri(WIN, taus):
    """93-usulü: pencereler-arası havuzlanmış D(τ) + vekil taban."""
    Dv, Vv = [], []
    for tau0 in taus:
        num = 0.0 + 0j; den = 0.0
        Gs, rs = [], []
        for (zz, tm, ds, Lw) in WIN:
            oms = tau0 * Lw + np.linspace(-0.02 * Lw, 0.02 * Lw, 160)
            oms = np.array([o for o in oms
                            if min(abs(o - l) for l in LINES_ALL) > 0.01])
            for s0 in range(0, len(oms), 40):
                ob = oms[s0:s0 + 40]
                rr = np.exp(1j * np.outer(ob, zz)).sum(axis=1)
                GG = (np.exp(1j * np.outer(ob, tm)) * ds[None, :]).sum(axis=1)
                num += (GG * np.conj(rr)).sum()
                den += (np.abs(rr)**2).sum()
                Gs.append(GG); rs.append(rr)
        kap = 2 * np.pi * tau0
        c = 2 * np.sin(kap / 2) / kap
        G_all, r_all = np.concatenate(Gs), np.concatenate(rs)
        fl = []
        for _ in range(40):
            perm = rng.permutation(len(r_all))
            fl.append(abs((G_all * np.conj(r_all[perm])).sum()) / (c * den))
        Dv.append(abs(num) / (c * den)); Vv.append(np.mean(fl))
    return Dv, Vv

TAUS = [0.04, 0.05, 0.06, 0.068, 0.075, 0.085, 0.095, 0.105, 0.113,
        0.125, 0.14, 0.16, 0.20, 0.30]

def duz_segmentler(zz, sayim, min_n=3000, win=80, esik=0.5, pad=160):
    """Sayım-sürüklenmesi düz parçalara böl; basamak çevresini at."""
    d = np.arange(len(zz)) - (sayim(zz) - sayim(zz[0]))
    from numpy.lib.stride_tricks import sliding_window_view
    med = np.median(sliding_window_view(d, win), axis=1)
    adim = med[win + 1:] - med[:-(win + 1)]
    bad = np.zeros(len(zz), bool)
    for i in np.where(np.abs(adim) > esik)[0] + win // 2:
        bad[max(0, i - pad):i + 2 * pad] = True
    seg, kes = [], 0
    s0 = 0
    for i in range(1, len(zz) + 1):
        if i == len(zz) or bad[i] != bad[i - 1]:
            if not bad[s0] and i - s0 >= min_n:
                seg.append(zz[s0:i])
            if i < len(zz) and bad[i]:
                kes += 1
            s0 = i
    return seg, kes

# ζ referansı — 41'in 6 penceresi (aynı düzlük disiplini, pürüzsüz sayım)
d41 = np.load(HERE / "41_bigT_windows.npz")
K41 = sorted({x.split("_")[1] for x in d41.files}, key=lambda s: int(s[:-1]))
def rvm_sayim(t):
    x = t / TWO_PI
    return x * np.log(x / np.e)
WINZ = []
nk_z = 0
for k in K41[:6]:
    gz, tm = d41[f"gaps_{k}"], d41[f"tmid_{k}"]
    zz = np.empty(len(gz) + 1)
    zz[0] = tm[0] - gz[0] / 2
    zz[1:] = zz[0] + np.cumsum(gz)
    seg, kes = duz_segmentler(zz, rvm_sayim)
    nk_z += kes
    WINZ += [pencere_hazirla(s, 1.0) for s in seg]
print(f"[zeta] {len(WINZ)} segment ({nk_z} kesim)", flush=True)
SONUC = {}
SONUC["zeta"] = egri(WINZ, TAUS)

exec(open(HERE / "98_L_motoru.py").read().split('CHI4 = ')[0])
CHI4 = {0: 0, 1: 1, 2: 0, 3: -1}
CHI3 = {0: 0, 1: 1, 2: -1}
CHI5 = {0: 0, 1: 1, 2: 1j, 3: -1j, 4: -1}
z6c = np.exp(1j * np.pi / 3)
CHI7 = {0: 0, 1: 1, 3: z6c, 2: z6c**2, 6: z6c**3, 4: z6c**4, 5: z6c**5}

for etiket, q, tab in [("chi3", 3, CHI3), ("beta", 4, CHI4),
                       ("chi5", 5, CHI5), ("chi7", 7, CHI7)]:
    zc = np.load(HERE / f"101f_{etiket}_zeros.npz")["zeros"]
    lo = np.exp(np.log(zc[0] + 1) + 0.25 * (np.log(zc[-1]) - np.log(zc[0] + 1)))
    zc = zc[zc >= lo]
    M = Lmotor(q, tab, 1)
    seg, kes = duz_segmentler(zc, lambda t: M.theta(t) / np.pi)
    # segmentleri 3 log-pencere kenarından da böl (τ çözünürlüğü)
    kenar = np.exp(np.linspace(np.log(zc[0]), np.log(zc[-1] * 1.0001), 4))
    WIN = []
    for s in seg:
        for i in range(3):
            p = s[(s >= kenar[i]) & (s < kenar[i + 1])]
            if len(p) >= 3000:
                WIN.append(pencere_hazirla(p, float(q)))
    SONUC[etiket] = egri(WIN, TAUS)
    print(f"[{etiket}] n={len(zc)}, {len(WIN)} segment ({kes} kesim), "
          f"L'ler: " + " ".join(f"{w[3]:.2f}" for w in WIN), flush=True)

adlar = ["zeta", "chi3", "beta", "chi5", "chi7"]
print(f"\n{'τ':>6} " + " ".join(f"{a:>6}" for a in adlar) + "   vekil(maks)")
for i, t in enumerate(TAUS):
    vek = max(SONUC[a][1][i] for a in adlar)
    print(f"{t:>6.3f} " + " ".join(f"{SONUC[a][0][i]:>6.3f}" for a in adlar)
          + f"   {vek:>6.3f}")
print("\nEşikler: log2/L ≈ 0.068-0.073 (ζ,χ₃,χ₅,χ₇) | log3/L ≈ 0.113 (β)")
