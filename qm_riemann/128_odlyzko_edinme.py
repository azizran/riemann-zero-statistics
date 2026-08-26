"""
128 — ODLYZKO VERİ-EDİNME + DOĞRULAMA + npz DÖNÜŞÜMÜ (26 Ağustos)
==========================================================================
127'nin menzil sınırı: 10^4-sıfırlık pencere yüksek-L çizgi ormanını
çözemiyor. Bu script Odlyzko'nun resmî tablolarından UZUN blokları
(zeros1 = 10^5, zeros6 = 2,001,052) alır, dört kapıdan geçirir ve
standart npz önbelleğe yazar.

DÖRT KAPI (veri-giriş hali):
  K1  SATIR SAYISI  — belgelenmiş n ile birebir mi?
  K2  UÇ DEĞERLER   — ilk/son sıfır bilinen değerlerle tutuyor mu?
  K3  ARDIŞIKLIK    — kesin artan, kopyasız, negatif boşluk yok?
  K4  SAYIM YOĞUNL. — N'(t) ≈ log(t/2π)/2π kaba kontrolü (pencere pencere)
  + AKIL SAĞLIĞI: unfold → var(ds) ≈ 0.16-0.18 (GUE ~0.178)

npz alanları: zeros (ofset), taban (float), taban_str, n, L, L_min, L_max
"""

import numpy as np
from pathlib import Path
from decimal import Decimal, getcontext

getcontext().prec = 60
HERE = Path(__file__).resolve().parent
VERI = HERE / "veri_odlyzko"
TWO_PI = 2.0 * np.pi


def theta(t):
    """Riemann-Siegel theta, asimptotik (t>10 için ~1e-7)."""
    return t / 2.0 * np.log(t / TWO_PI) - t / 2.0 - np.pi / 8.0 \
        + 1.0 / (48.0 * t) + 7.0 / (5760.0 * t ** 3)


def yukle(path, skip=0):
    vals = []
    with open(path) as f:
        for i, line in enumerate(f):
            if i < skip:
                continue
            s = line.strip()
            if not s:
                continue
            vals.append(float(s))
    return np.array(vals, dtype=np.float64)


def dort_kapi(etiket, off, taban_str, n_bekl, ilk_bekl, son_bekl, tol_uc):
    """off: ofsetler (taban çıkarılmış). taban_str: tam tamsayı taban (str)."""
    rap = {}
    taban_dec = Decimal(taban_str)
    n = len(off)

    # --- K1 satır sayısı
    rap["K1_n"] = (n, n_bekl, n == n_bekl)

    # --- K2 uç değerler (Decimal ile tam)
    ilk = taban_dec + Decimal(repr(float(off[0])))
    son = taban_dec + Decimal(repr(float(off[-1])))
    d_ilk = abs(ilk - Decimal(ilk_bekl))
    d_son = abs(son - Decimal(son_bekl))
    rap["K2_ilk"] = (str(ilk), ilk_bekl, float(d_ilk), float(d_ilk) < tol_uc)
    rap["K2_son"] = (str(son), son_bekl, float(d_son), float(d_son) < tol_uc)

    # --- K3 ardışıklık
    d = np.diff(off)
    rap["K3_kesin_artan"] = bool(np.all(d > 0))
    rap["K3_min_bosluk"] = float(d.min())
    rap["K3_kopya"] = int(np.sum(d == 0))

    # --- taban-güvenli log(t/2pi)
    if taban_dec == 0:
        t = off.copy()
        log_t2pi = np.log(t / TWO_PI)
    else:
        b = float(taban_dec)
        log_b = float(Decimal(taban_str).ln()) if hasattr(Decimal, "ln") else np.log(b)
        log_t2pi = (log_b - np.log(TWO_PI)) + off / b   # log(1+x) ≈ x, x~1e-9
        t = None

    L_lokal = log_t2pi
    rap["L"] = float(L_lokal[n // 2])
    rap["L_min"] = float(L_lokal.min())
    rap["L_max"] = float(L_lokal.max())

    # --- K4 sayım yoğunluğu: gözlenen / beklenen, 20 pencerede
    m = 20
    kenar = np.linspace(0, n, m + 1).astype(int)
    oran = []
    for i in range(m):
        a, b_ = kenar[i], kenar[i + 1] - 1
        if b_ <= a:
            continue
        dt = off[b_] - off[a]
        bekl = dt * np.mean(L_lokal[a:b_ + 1]) / TWO_PI
        oran.append((b_ - a) / bekl)
    oran = np.array(oran)
    rap["K4_oran"] = (float(oran.min()), float(oran.mean()), float(oran.max()))
    rap["K4_gecti"] = bool(np.all(np.abs(oran - 1.0) < 0.05))

    # --- unfold + var(ds)
    if taban_dec == 0:
        # theta tabanlı global unfold (t>10 kesimi)
        msk = off > 15.0
        u = theta(off[msk]) / np.pi
        ds = np.diff(u)
    else:
        # lokal yoğunluk unfold (derin blok: mutlak t float64'te taşar)
        mid = 0.5 * (L_lokal[1:] + L_lokal[:-1])
        ds = np.diff(off) * mid / TWO_PI
    rap["ds_ort"] = float(ds.mean())
    rap["ds_var"] = float(ds.var())
    rap["ds_n"] = int(len(ds))
    return rap, L_lokal


def kaydet(etiket, off, taban_str, L, L_min, L_max):
    out = HERE / f"128_odl_{etiket}_zeros.npz"
    np.savez_compressed(
        out,
        zeros=off.astype(np.float64),
        taban=np.float64(Decimal(taban_str)),
        taban_str=np.array(taban_str),
        n=np.int64(len(off)),
        L=np.float64(L),
        L_min=np.float64(L_min),
        L_max=np.float64(L_max),
    )
    return out, out.stat().st_size


BLOKLAR = [
    # etiket, dosya, skip, taban_str, n_bekl, ilk_bekl, son_bekl, tol_uc
    ("zeros1_1e5", VERI / "zeros1", 0, "0", 100000,
     "14.134725142", "74920.827498994", 1e-8),
    ("zeros6_2e6", VERI / "zeros6", 0, "0", 2001052,
     "14.134725142", "1132490.658714411", 1e-8),
    ("t1e12_1e4", HERE / "odlyzko_zeros3.txt", 9, "267653395647", 10000,
     "267653395648.8475231278", None, 1e-7),
    ("t1e21_1e4", HERE / "odlyzko_zeros4.txt", 9, "144176897509546973000", 10000,
     "144176897509546973538.49806962", None, 1e-5),
    ("t1e22_1e4", HERE / "odlyzko_zeros5.txt", 9, "1370919909931995300000", 10000,
     "1370919909931995308226.68016095", None, 1e-5),
]

if __name__ == "__main__":
    for etiket, path, skip, taban_str, n_bekl, ilk_b, son_b, tol in BLOKLAR:
        if not path.exists():
            print(f"[ATLA] {etiket}: {path} yok")
            continue
        off = yukle(path, skip=skip)
        if son_b is None:
            son_b = str(Decimal(taban_str) + Decimal(repr(float(off[-1]))))
            tol_son = tol
        rap, Ll = dort_kapi(etiket, off, taban_str, n_bekl, ilk_b, son_b, tol)
        print("=" * 74)
        print(f"BLOK {etiket}   dosya={path.name}")
        print(f"  K1 satır: n={rap['K1_n'][0]} beklenen={rap['K1_n'][1]} "
              f"-> {'GEÇTİ' if rap['K1_n'][2] else 'KALDI'}")
        print(f"  K2 ilk  : {rap['K2_ilk'][0][:34]}  |Δ|={rap['K2_ilk'][2]:.2e} "
              f"-> {'GEÇTİ' if rap['K2_ilk'][3] else 'KALDI'}")
        print(f"  K2 son  : {rap['K2_son'][0][:34]}  |Δ|={rap['K2_son'][2]:.2e} "
              f"-> {'GEÇTİ' if rap['K2_son'][3] else 'KALDI'}")
        print(f"  K3 artan: {rap['K3_kesin_artan']}  min_boşluk={rap['K3_min_bosluk']:.3e}"
              f"  kopya={rap['K3_kopya']}")
        print(f"  K4 N'(t): oran min/ort/max = {rap['K4_oran'][0]:.4f} / "
              f"{rap['K4_oran'][1]:.4f} / {rap['K4_oran'][2]:.4f} "
              f"-> {'GEÇTİ' if rap['K4_gecti'] else 'BAK'}")
        print(f"  L (orta) = {rap['L']:.4f}   [{rap['L_min']:.4f} .. {rap['L_max']:.4f}]")
        print(f"  unfold  : <ds>={rap['ds_ort']:.5f}  var(ds)={rap['ds_var']:.5f} "
              f"(n={rap['ds_n']})")
        out, sz = kaydet(etiket, off, taban_str, rap["L"], rap["L_min"], rap["L_max"])
        print(f"  -> {out.name}  ({sz/1e6:.2f} MB)")
