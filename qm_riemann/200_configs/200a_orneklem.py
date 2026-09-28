"""
200a_orneklem.py — KALEM 200-A: gerçekleştirilebilir (realizable) sentetik örnekleyiciler
================================================================================================

Hiçbir gerçek Z/Z' verisi YOK — yalnız model tanımlarından (KALEM Model ailesi F) çekilen
sentetik örnekler. Kullanılır: 200a_m6_sentetik.py (M6: sentetik güç + SE kalibrasyonu,
seçim gücü) ve 200a_analiz.py'nin sentetik-blok uçtan-uca testinde (fixture üretimi).

Ürün yasası (BHNY Prop 2.2, |Λ|^{2b}-eğik genellemesi, TEOREM_KUCUK_ARALIK_ISKELET §3, §6b):
  X_n = Π_{j=1}^n |1-ξ_j|,  n = N - b (N: CUE boyutu, b: eğim üsteli / 2)
  - ξ_1 (çember): yoğunluk ∝ |1-e^{iα}|^{2b}  (b=0: düzgün)
  - ξ_j (j>=2, disk): yoğunluk ∝ (1-|z|²)^{j-2} |1-z|^{2b}  (b=0: |z|²~Beta(1,j-1), düzgün faz)
  Genel b için ret-örnekleme (rejection sampling): b=0 tabanından öner, |1-z|^{2b}/4^b ile kabul et
  (4^b = sup_{|z|<=1}|1-z|^{2b}, z→-1'de erişilir).

  N tam sayı değilse: floor(N-b) ve ceil(N-b) çarpan sayılarının KARIŞIMI, ağırlık = kesir kısmı
  (KALEM M6(a): "non-integer N: mixture of floor/ceil sizes with weight = fractional part").

hyb(X): yukarıdaki tilted-CUE(N_X) örneklemesi + bağımsız rastgele Euler çarpanları
  (asal p<=X, faz düzgün [0,2π)): −log|1 − p^{−1/2} e^{iα_p}| toplamı.

Gauss: b=0 log|N(0,1)|  (= log chi_1);  b=1 log chi_2 (Rayleigh, Kac–Rice eğilmiş).
"""
import numpy as np
from sympy import primerange

EULER_GAMMA = np.euler_gamma


def _oversample_mult(b):
    return {0: 1.15, 1: 6.0, 2: 24.0}.get(b, 8.0 * 4.0 ** max(b - 1, 0))


def sample_circle_factor(b, n, rng):
    """log|1-xi_1|, xi_1 ~ |1-e^{i alpha}|^{2b}-tilted uniform on circle."""
    if n == 0:
        return np.empty(0)
    if b == 0:
        alpha = rng.uniform(0, 2 * np.pi, n)
        return np.log(np.abs(1 - np.exp(1j * alpha)))
    M = 4.0 ** b
    mult = _oversample_mult(b)
    out = np.empty(n)
    got = 0
    while got < n:
        m = int((n - got) * mult) + 16
        alpha = rng.uniform(0, 2 * np.pi, m)
        w = np.abs(1 - np.exp(1j * alpha)) ** (2 * b)
        u = rng.uniform(0, 1, m) * M
        keep = alpha[u < w]
        take = min(len(keep), n - got)
        if take:
            out[got:got + take] = np.log(np.abs(1 - np.exp(1j * keep[:take])))
            got += take
    return out


def sample_disk_factor(j, b, n, rng):
    """log|1-xi_j| (j>=2), xi_j ~ (1-|z|^2)^{j-2}|1-z|^{2b}-tilted on unit disk."""
    if n == 0:
        return np.empty(0)
    if b == 0:
        U = rng.uniform(0, 1, n)
        r2 = 1 - U ** (1.0 / (j - 1))
        r = np.sqrt(r2)
        theta = rng.uniform(0, 2 * np.pi, n)
        z = r * np.exp(1j * theta)
        return np.log(np.abs(1 - z))
    M = 4.0 ** b
    mult = _oversample_mult(b)
    out = np.empty(n)
    got = 0
    while got < n:
        m = int((n - got) * mult) + 16
        U = rng.uniform(0, 1, m)
        r2 = 1 - U ** (1.0 / (j - 1))
        r = np.sqrt(r2)
        theta = rng.uniform(0, 2 * np.pi, m)
        z = r * np.exp(1j * theta)
        w = np.abs(1 - z) ** (2 * b)
        u = rng.uniform(0, 1, m) * M
        keep_mask = u < w
        zk = z[keep_mask]
        take = min(len(zk), n - got)
        if take:
            out[got:got + take] = np.log(np.abs(1 - zk[:take]))
            got += take
    return out


def sample_log_tilted_product(N, b, n, rng):
    """log X_{N-b} under the 2b-tilted product law; N may be non-integer (floor/ceil mixture,
    weight = fractional part -> matches E[#factors]=N-b)."""
    nf = N - b
    if nf < 0:
        nf = 0.0
    n_lo = int(np.floor(nf))
    frac = nf - n_lo
    n_hi = n_lo + 1

    use_hi = rng.uniform(0, 1, n) < frac
    idx_lo = np.where(~use_hi)[0]
    idx_hi = np.where(use_hi)[0]
    out = np.empty(n)
    for idx_set, nfac in ((idx_lo, n_lo), (idx_hi, n_hi)):
        m = len(idx_set)
        if m == 0:
            continue
        acc = np.zeros(m)
        if nfac >= 1:
            acc += sample_circle_factor(b, m, rng)
        for j in range(2, nfac + 1):
            acc += sample_disk_factor(j, b, m, rng)
        out[idx_set] = acc
    return out


def sample_euler_factor(primes, n, rng):
    out = np.zeros(n)
    for p in primes:
        alpha = rng.uniform(0, 2 * np.pi, n)
        out += -np.log(np.abs(1 - p ** -0.5 * np.exp(1j * alpha)))
    return out


def primes_upto(X):
    return list(primerange(2, int(np.floor(X)) + 1))


def sample_hybrid(X, L, b, n, rng):
    NX = L / (np.e ** EULER_GAMMA * np.log(X))
    primes = primes_upto(X)
    cue = sample_log_tilted_product(NX, b, n, rng)
    eul = sample_euler_factor(primes, n, rng) if primes else np.zeros(n)
    return cue + eul, NX


def sample_gauss(b, n, rng):
    if b == 0:
        return np.log(np.abs(rng.standard_normal(n)))
    if b == 1:
        r = np.sqrt(rng.standard_normal(n) ** 2 + rng.standard_normal(n) ** 2)
        return np.log(r)
    raise NotImplementedError(f"b={b} Gauss örneklemesi 200-A kapsamında değil")


def sample_model(model, L, b, n, rng):
    """Genel arayüz: model in {'CUE','hyb2','hyb3','hyb5','hyb7','hyb11','Gauss'} (a_k HARİÇ:
    a_k'nın 200t_tam_yasa.py'ye göre sonlu-L olasılık yasası yok, bkz KALEM Bulgu-1)."""
    if model == "CUE":
        return sample_log_tilted_product(L, b, n, rng), L
    if model == "Gauss":
        return sample_gauss(b, n, rng), None
    if model.startswith("hyb"):
        X = int(model[3:])
        return sample_hybrid(X, L, b, n, rng)
    raise ValueError(f"'{model}' gerçekleştirilebilir bir örnekleyici değil (a_k analitik-yalnız)")
