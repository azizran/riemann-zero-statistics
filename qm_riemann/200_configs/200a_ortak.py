"""
200a_ortak.py — KALEM 200-A: ortak yardımcı fonksiyonlar (k-istatistikleri, blok
jackknife, khi-kare seçim hattı). 200a_olcum.py, 200a_analiz.py ve 200a_m6_sentetik.py
tarafından ortak kullanılır. Hiçbir gerçek veri okumaz / yazmaz — yalnız sayısal
yardımcı fonksiyonlar.
"""
import json
import numpy as np
from scipy import special
from pathlib import Path

HERE = Path(__file__).resolve().parent
WINDOWS = ("W1", "W2", "W3")
RUNGS = (0, 1)
MODELS = ("CUE", "a_k", "hyb2", "hyb3", "hyb5", "hyb7", "hyb11", "Gauss")
ARITHMETIC_MODELS = ("a_k", "hyb2", "hyb3", "hyb5", "hyb7", "hyb11")
SIGMA_TEORI = {2: 0.0188, 3: 0.0334}   # sigma_teori,r  (önceden ilan, KALEM)

# 12-bileşenli artık vektörünün sırası: (pencere, basamak, kümülant)
STAT_ORDER = [(w, b, r) for w in WINDOWS for b in RUNGS for r in (2, 3)]


def load_tahmin(path=None):
    path = path or (HERE / "200a_tahmin.json")
    with open(path) as f:
        return json.load(f)


# ---------------------------------------------------------------------------
# Vektörize kap(r,N,b) — eğik CUE log-kümülantları, tam-sayı olmayan N için analitik
# devam (200t_merdiven.py'nin mpmath sürümüyle 1e-14 mertebesinde örtüşür; bkz.
# scratchpad/k200a_makine doğrulaması). scipy.special (digamma/polygamma) kullanır,
# tek tek mpmath çağrısından ~1e5x hızlı — M7 ve sürekli-X uyumu (200a_analiz.py)
# için gerekli.
# ---------------------------------------------------------------------------
def Sr_vec(r, n, a):
    x = a + 1
    if r == 1:
        return (x + n - 1) * special.digamma(x + n) - (x - 1) * special.digamma(x) - n
    if r == 2:
        return (special.digamma(x + n) + (x + n - 1) * special.polygamma(1, x + n)
                - special.digamma(x) - (x - 1) * special.polygamma(1, x))
    if r == 3:
        return (2 * special.polygamma(1, x + n) + (x + n - 1) * special.polygamma(2, x + n)
                - 2 * special.polygamma(1, x) - (x - 1) * special.polygamma(2, x))
    raise ValueError(r)


def kap_vec(r, N, b):
    n = N - b
    return Sr_vec(r, n, 2 * b) - 2.0 ** (1 - r) * Sr_vec(r, n, b)


# Asal başına aritmetik/hibrit kümülant katkıları (KALEM: kappa2_p = 1/2 Li2(1/p),
# kappa3_p = 3/2 S_{1,2}(1/p)), mpmath ile bir kez, X<=13 için önbelleğe alınır.
_PRIME_CUM_CACHE = {}


def prime_cumulant(p):
    if p in _PRIME_CUM_CACHE:
        return _PRIME_CUM_CACHE[p]
    import mpmath as mp
    mp.mp.dps = 30
    x = mp.mpf(1) / p
    k2 = float(mp.polylog(2, x) / 2)
    s = mp.mpf(0)
    H = mp.mpf(0)
    m = 1
    t = x
    while True:
        m += 1
        H += mp.mpf(1) / (m - 1)
        t *= x
        term = H * t / m ** 2
        s += term
        if term < mp.mpf(10) ** -25:
            break
    k3 = float(1.5 * s)
    _PRIME_CUM_CACHE[p] = (k2, k3)
    return k2, k3


def primes_leq(X):
    from sympy import primerange
    return list(primerange(2, int(np.floor(X)) + 1))


# ---------------------------------------------------------------------------
# k-istatistikleri (yansız), ham güç toplamlarından (MathWorld "k-statistic")
# n, S1=sum x, S2=sum x^2, S3=sum x^3, S4=sum x^4
# ---------------------------------------------------------------------------
def kstat_from_sums(n, s1, s2, s3):
    """Yansız k2, k3 — ham güç toplamlarından (n>=3 gerekir)."""
    n = np.asarray(n, dtype=np.float64)
    s1 = np.asarray(s1, dtype=np.float64)
    s2 = np.asarray(s2, dtype=np.float64)
    s3 = np.asarray(s3, dtype=np.float64)
    k2 = (n * s2 - s1 ** 2) / (n * (n - 1))
    k3 = (2 * s1 ** 3 - 3 * n * s1 * s2 + n ** 2 * s3) / (n * (n - 1) * (n - 2))
    return k2, k3


def power_sums_of(x):
    """Bir x örnekleminden (n, S1, S2, S3, S4) — yalnız TOPLAMLAR döner, x atılabilir."""
    x = np.asarray(x, dtype=np.float64)
    n = x.size
    s1 = float(np.sum(x))
    s2 = float(np.sum(x ** 2))
    s3 = float(np.sum(x ** 3))
    s4 = float(np.sum(x ** 4))
    return np.array([n, s1, s2, s3, s4], dtype=np.float64)


def merge_sums(list_of_sums):
    """Birkaç (n,S1..S4) satırını topla (blok birleştirme / jackknife 'hariç tut')."""
    arr = np.asarray(list_of_sums, dtype=np.float64)
    return arr.sum(axis=0)


def merge_adjacent_blocks(block_sums, factor=2):
    """(...,B,5) dizisindeki bitişik `factor` bloğu birleştirerek (...,B/factor,5) üretir."""
    block_sums = np.asarray(block_sums, dtype=np.float64)
    B = block_sums.shape[-2]
    assert B % factor == 0, f"blok sayısı {B} {factor}'e bölünmeli"
    newshape = block_sums.shape[:-2] + (B // factor, factor, block_sums.shape[-1])
    return block_sums.reshape(newshape).sum(axis=-2)


# ---------------------------------------------------------------------------
# Bir-blok-dışarıda (leave-one-block-out) jackknife: pencere başına 4x4 kovaryans
# blocks_win: shape (B, 2, 5)  [blok, basamak(0/1), (n,S1,S2,S3,S4)]
# döndürür: point estimate (k2_b0,k3_b0,k2_b1,k3_b1), 4x4 kovaryans, jk replikaları
# ---------------------------------------------------------------------------
def jackknife_4x4(blocks_win):
    blocks_win = np.asarray(blocks_win, dtype=np.float64)
    B = blocks_win.shape[0]
    assert blocks_win.shape[1] == 2 and blocks_win.shape[2] == 5

    total = blocks_win.sum(axis=0)  # (2,5)
    point = np.empty(4)
    for bi, b in enumerate(RUNGS):
        n, s1, s2, s3, s4 = total[bi]
        k2, k3 = kstat_from_sums(n, s1, s2, s3)
        point[2 * bi] = k2
        point[2 * bi + 1] = k3

    reps = np.empty((B, 4))
    for i in range(B):
        loo = total - blocks_win[i]  # (2,5): i'inci blok hariç toplam
        for bi in range(2):
            n, s1, s2, s3, s4 = loo[bi]
            k2, k3 = kstat_from_sums(n, s1, s2, s3)
            reps[i, 2 * bi] = k2
            reps[i, 2 * bi + 1] = k3

    rep_bar = reps.mean(axis=0)
    d = reps - rep_bar[None, :]
    cov = (B - 1) / B * (d.T @ d)
    return point, cov, reps


def jackknife_var_iid(blocks_win):
    """Basit iid varsayımlı (bağımsız-özdeş dağılım) varyans tahmini — Var_jk/Var_iid
    tanısı için. Delta yöntemiyle k2,k3'ün örnekleme varyansı (yaklaşık, normal-teori):
      Var(k2) ~ 2 k2^2/n ;  Var(k3) ~ 6 k2^3/n   (standart 4./6. kümülant terimleri ihmal)
    Yalnızca TANI amaçlı (M8/rapor), karar hesaplarında KULLANILMAZ (jackknife kullanılır).
    """
    blocks_win = np.asarray(blocks_win, dtype=np.float64)
    total = blocks_win.sum(axis=0)
    out = np.empty(4)
    for bi in range(2):
        n, s1, s2, s3, s4 = total[bi]
        k2, k3 = kstat_from_sums(n, s1, s2, s3)
        out[2 * bi] = 2 * k2 ** 2 / n
        out[2 * bi + 1] = 6 * max(k2, 1e-12) ** 3 / n
    return out


# ---------------------------------------------------------------------------
# Khi-kare seçim hattı
# ---------------------------------------------------------------------------
def build_full_cov(cov_per_window, f=None):
    """3 tane 4x4 kovaryanstan (pencere başına) blok-köşegen 12x12 kurar.
    f: None ya da uzunluk-12 kalibrasyon vektörü (C_ij -> f_i f_j C_ij)."""
    C = np.zeros((12, 12))
    for wi in range(3):
        C[4 * wi:4 * wi + 4, 4 * wi:4 * wi + 4] = cov_per_window[wi]
    if f is not None:
        f = np.asarray(f, dtype=np.float64)
        C = C * np.outer(f, f)
    return C


def sigma_teori_vector():
    return np.array([SIGMA_TEORI[r] for (_, _, r) in STAT_ORDER])


def model_prediction_vector(tahmin, model):
    """tahmin.json'dan model'in 12-bileşenli tahmin vektörü (STAT_ORDER sırasıyla)."""
    out = np.empty(12)
    for i, (w, b, r) in enumerate(STAT_ORDER):
        pred = tahmin[w]["pred"][f"b{b}"][model]
        out[i] = pred[0] if r == 2 else pred[1]
    return out


def chi2(resid, C_full):
    Cinv = np.linalg.inv(C_full)
    return float(resid @ Cinv @ resid)


def chi2_per_window(resid, cov_per_window_full):
    """Diagnostik: her pencerenin kendi 4x4 bloğundan khi-kare katkısı (toplamı = tam khi-kare
    yalnızca C blok-köşegen olduğunda tam eşittir; f-kalibrasyonu blok-köşegenliği bozmaz)."""
    out = []
    for wi in range(3):
        r = resid[4 * wi:4 * wi + 4]
        C = cov_per_window_full[4 * wi:4 * wi + 4, 4 * wi:4 * wi + 4]
        out.append(float(r @ np.linalg.inv(C) @ r))
    return out
