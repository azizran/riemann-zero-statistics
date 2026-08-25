"""
108b — ADA PERDESİ, TAM-TABAN ALETİYLE (24 Ağustos, kaptan kuyusu)
==========================================================================
108'in tek-çizgi okuması tuzağa düştü (D≈1.05-1.15, ζ'da bile — 81-83'ün
örnekleme-ızgarası şişmesi; ζ kontrol sütunu tuzağı ele verdi → ders:
perde ancak TAM-TABAN ortak regresyonuyla ölçülür). Burada 83'ün
chunked normal-denklem makinesi adalara uyarlanır.

KAVRAM DENETİMİ (docstring'e mühür): 92-T1b D'yi v/(2τ·a_v) ile,
89 ise rigid=(2/π)sin(πτ) ile yazmış — iki norm τ büyüdükçe ayrışır
(0.55'te 1.75×). Burada İKİSİ de basılır; fizik tartışması sinüs-normda
(89'un türetilmiş katı biçimi) yürütülür.

ÖN-MÜHÜR (P2 revize): sinüs-normda gerçek perde daha yatıksa ayrışma
küçülür; korunan öngörü SIRALAMA — aynı τ'da ölü-2 üçlüsü (soğuk)
canlı-2'den YÜKSEK D; ln-oran ≈ σ² oranı 0.80±0.2. P1: her ada D<1.
P3: kervan-içi çökme.

SONUÇ (24 Ağustos):
  ALET DOĞRULANDI: ζ sütunu 2τ-normda bilinen perdeyi üretti
  (0.94→0.54 ≈ kampanyanın 0.91→0.48'i; 108'in tek-çizgi şişmesi
  tam-tabanla kayboldu). KAVRAM: sinüs-normda perde 0.955→0.823 —
  "0.48'e düşüş"ün çoğu 2τ/sinüs geometrisiymiş (89'un türetilmiş
  katı biçimi).
  P1 ✓ (her ailede D<1, τ ile iner). P3 ✓✓✓ — SEKİZ SÜTUN TEK EĞRİ
  (sabit p'de saçılma ±0.01-0.02): perde EVRENSEL. P2 RET — ölü-2
  üçlüsü canlı-2 ile binde-birkaç içinde özdeş (ln-oran ~1.0, 0.80
  değil): perde SICAKLIKTAN BAĞIMSIZ → sıcaklık-DW hipotezi öldü;
  perde yalnız yerel/evrensel istatistiğe bağlı olabilir. a_v ve
  fiziksellik: 108c (a_v≈1.00; D_fiz 0.958→0.797; seçicilik bonusu).
"""

import numpy as np
from pathlib import Path
from sympy import primerange

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
Q_CAP = 720

def pk_list(lim):
    out = []
    for p in primerange(2, int(lim) + 1):
        q = p
        while q <= lim:
            out.append(q); q *= p
    return sorted(set(out))

def chunked_v(yg, tmid, qs, chunk=40000):
    """83'ün v-kanalı: sabit + 2 sürüklenme + cos/sin(q) ortak regresyon."""
    tt = (tmid - tmid.mean()) / (tmid[-1] - tmid[0])
    C = 3 + 2 * len(qs)
    XtX = np.zeros((C, C)); Xty = np.zeros(C)
    n = len(yg)
    def cols(sl):
        c = [np.ones(sl.stop - sl.start), tt[sl], tt[sl]**2]
        for q in qs:
            arg = tmid[sl] * np.log(q)
            c += [np.cos(arg), np.sin(arg)]
        return np.vstack(c).T
    for s0 in range(0, n, chunk):
        sl = slice(s0, min(s0 + chunk, n))
        Xc = cols(sl)
        XtX += Xc.T @ Xc; Xty += Xc.T @ yg[sl]
    b = np.linalg.solve(XtX, Xty)
    return b

def ada_perde(zlist, qeff):
    """Sertifikalı sıfır listeleri → çizgi başına (p, τ, v/p^{-1/2})."""
    acc = {}
    for zz in zlist:
        g = np.diff(zz)
        m = 0.5 * (zz[:-1] + zz[1:])
        Lw = np.log(qeff * m / TWO_PI)
        L = float(Lw.mean())
        yg = np.log(g * Lw / TWO_PI)
        qs = pk_list(min(np.exp(0.52 * L), Q_CAP))
        b = chunked_v(yg, m, qs)
        for i, q in enumerate(qs):
            iv = 3 + 2 * i
            v = np.hypot(b[iv], b[iv + 1]) / q**-0.5
            acc.setdefault(q, []).append((np.log(q) / L, v))
    return {q: (np.mean([t for t, v in vs]), np.mean([v for t, v in vs]))
            for q, vs in acc.items()}

SETS = [
    ("zeta",  1, None, [], "41"),
    ("chi3",  3, "101f_chi3_zeros.npz",  {3}, None),
    ("beta",  4, "101f_beta_zeros.npz",  {2}, None),
    ("chi5",  5, "101f_chi5_zeros.npz",  {5}, None),
    ("chi7",  7, "101f_chi7_zeros.npz",  {7}, None),
    ("chi5e", 5, "105b_chi5e_zeros.npz", {5}, None),
    ("chi8e", 8, "105b_chi8e_zeros.npz", {2}, None),
    ("chi8o", 8, "105b_chi8o_zeros.npz", {2}, None),
]
OLU2 = {"beta", "chi8e", "chi8o"}
ASAL = set(primerange(2, 200))

SONUC = {}
for ad, qeff, dosya, olu, src in SETS:
    if src == "41":
        d41 = np.load(HERE / "41_bigT_windows.npz")
        zlist = []
        for k in ["120k", "200k"]:
            gz, tm = d41[f"gaps_{k}"], d41[f"tmid_{k}"]
            zz = np.empty(len(gz) + 1)
            zz[0] = tm[0] - gz[0] / 2
            zz[1:] = zz[0] + np.cumsum(gz)
            zlist.append(zz)
    else:
        zc = np.load(HERE / dosya)["zeros"]
        lo = np.exp(np.log(zc[0] + 1) + 0.25 * (np.log(zc[-1]) - np.log(zc[0] + 1)))
        zlist = [zc[zc >= lo]]
    SONUC[ad] = ada_perde(zlist, float(qeff))
    print(f"[{ad}] {len(SONUC[ad])} q", flush=True)

def D_sin(t, v):
    return v / ((2 / np.pi) * np.sin(np.pi * t))
def D_2t(t, v):
    return v / (2 * t)

print(f"\nAsal çizgiler, D_sin = v/((2/π)sin(πτ)) [parantez: D_2τ = v/2τ]:")
print(f"{'p':>4}", end="")
for ad, *_ in SETS:
    print(f" {ad:>7}", end="")
print()
plist = sorted({q for ad, *_ in SETS for q in SONUC[ad]
                if q in ASAL and 0.06 < SONUC[ad][q][0] < 0.52})
for p in plist:
    print(f"{p:>4}", end="")
    for ad, qeff, dosya, olu, src in SETS:
        if p in olu or p not in SONUC[ad]:
            print(f" {'—':>7}", end="")
        else:
            t, v = SONUC[ad][p]
            print(f" {D_sin(t, v):>7.3f}", end="")
    print()

print("\nτ-binli kervan ortalamaları (yalnız asallar; her iki norm):")
print(f"{'τ-bin':>12} {'canlı D_sin':>11} {'ölü D_sin':>10} {'ln-oran':>8} "
      f"{'canlı D_2τ':>10} {'ölü D_2τ':>9}")
for lo_, hi_ in [(0.06, 0.15), (0.15, 0.25), (0.25, 0.35), (0.35, 0.45),
                 (0.45, 0.52)]:
    cs, os_, c2, o2 = [], [], [], []
    for ad, qeff, dosya, olu, src in SETS:
        for q, (t, v) in SONUC[ad].items():
            if q in ASAL and q not in olu and lo_ <= t < hi_:
                (os_ if ad in OLU2 else cs).append(D_sin(t, v))
                (o2 if ad in OLU2 else c2).append(D_2t(t, v))
    if cs and os_:
        mc, mo = np.mean(cs), np.mean(os_)
        lr = np.log(mo) / np.log(mc) if 0 < mc < 1 else np.nan
        print(f"[{lo_:.2f},{hi_:.2f}) {mc:>11.3f} {mo:>10.3f} {lr:>8.3f} "
              f"{np.mean(c2):>10.3f} {np.mean(o2):>9.3f}  (n={len(cs)}/{len(os_)})")
np.savez(HERE / "108b_perde.npz",
         **{ad: np.array([(q, t, v) for q, (t, v) in sorted(SONUC[ad].items())])
            for ad in SONUC})
print("\nKayıt: 108b_perde.npz")
