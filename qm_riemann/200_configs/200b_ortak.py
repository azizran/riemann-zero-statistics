"""
200b_ortak.py — KALEM 200-B: ortak yardımcı fonksiyonlar (pencere tanımları, tepe
(hump) ölçümü — ızgara + parabolik rafinasyon, tek-parça (b=2, tek "rung") k-istatistiği
jackknife'ı). 200b_olcum.py, 200b_analiz.py, 200b_m1_m2_kapilar.py, 200b_m6_sentetik.py,
200b_m7_karisim.py tarafından ortak kullanılır.

Hiçbir gerçek Z/M/log(M/δ̃²) DAĞILIM ÖZETİ burada hesaplanmaz — yalnız sayısal
yardımcı fonksiyonlar (200a_ortak.py'nin b=2 sürümü). 200a_motor.py / 200a_ortak.py
DEĞİŞTİRİLMEDİ, yalnız İÇE AKTARILIR (görev kuralı: "Reusable code (do not edit)").
"""
import importlib.util
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi

# Pencereler — KALEM ile AYNI (200-A ile aynı sıfır-indeks aralıkları; görev metninde
# birebir verildi): W1 [20868,198239), W2 [198239,598742), W3 [598742,2001052).
WINDOW_IDX = {
    "W1": (20868, 198239),
    "W2": (198239, 598742),
    "W3": (598742, 2001052),
}
WINDOWS = ("W1", "W2", "W3")
EPSILONS = (0.1, 0.2, 0.3)
EPS_PRIMARY = 0.2
B_RUNG = 2  # yakın çift basamağı (bu KALEM'in konusu)

# Aritmetik sabitler (KALEM / 200t_aritmetik_kesin.py, yakınsamış)
A2_ARITH = -0.088124
A3_ARITH = 0.233653


def _load(modname):
    spec = importlib.util.spec_from_file_location(modname, HERE / f"{modname}.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


# ---------------------------------------------------------------------------
# Gerçek pencere sınırları (yalnız KONUM — zeros dosyasından t_a,t_b,idx; M/Z YOK)
# ---------------------------------------------------------------------------
def load_real_windows(zeros_path=None):
    zeros_path = zeros_path or (HERE.parent / "128_odl_zeros6_2e6_zeros.npz")
    zeros = np.load(zeros_path)["zeros"]
    out = {}
    for w, (i0, i1) in WINDOW_IDX.items():
        out[w] = {"idx": (i0, i1), "t_a": float(zeros[i0]), "t_b": float(zeros[i1 - 1]),
                  "n_zeros_window": int(i1 - i0)}
    return out


def load_sayimlar(path=None):
    path = path or (HERE / "200b_sayimlar.json")
    with open(path) as f:
        return json.load(f)


# ---------------------------------------------------------------------------
# Olay (event) çıkarımı — YALNIZ KONUM (γ_n, γ_{n+1}, m_n, L_n, δ̃_n). Hiçbir Z/M yok.
# ---------------------------------------------------------------------------
def events_from_zeros(zero_t, idx0, eps_max=0.3):
    """zero_t: bu pencerenin TÜM sıfırları (1B, artan). idx0: bu dilimin global başlangıç
    indeksi (n'nin mutlak sıfır-indeksi olarak raporlanması için).
    n+1 < idx1 koşulu ZATEN sağlanmış olmalı: zero_t, tam pencere aralığı olmalı ve
    çağıran, son elemanı ÇİFT OLUŞTURMAK için kullanmamalı (bu fonksiyon zero_t'nin
    ardışık TÜM çiftlerini üretir: len(zero_t)-1 tane çift; pencere sınırını AŞAN çift
    yoktur çünkü zero_t zaten [idx0,idx1) dilimidir ve n+1<idx1 otomatik sağlanır).
    Döner: dict(n=global sıfır indeksi (n), m=orta nokta, L=log(m/2pi), delta_tilde, delta)
    """
    zero_t = np.asarray(zero_t, dtype=np.float64)
    ga = zero_t[:-1]
    gb = zero_t[1:]
    n_idx = idx0 + np.arange(len(zero_t) - 1, dtype=np.int64)
    m = 0.5 * (ga + gb)
    L = np.log(m / TWO_PI)
    delta = gb - ga
    delta_tilde = delta * L / TWO_PI
    mask = delta_tilde < eps_max
    return {"n": n_idx[mask], "ga": ga[mask], "gb": gb[mask], "m": m[mask], "L": L[mask],
            "delta": delta[mask], "delta_tilde": delta_tilde[mask]}


# ---------------------------------------------------------------------------
# Tepe (hump) ölçümü: [ga,gb] aralığında |Z| maksimumu, n_grid iç nokta ızgarası +
# 3-nokta parabolik rafinasyon (|Z| üzerinde, log değil — KALEM: "for |Z|").
# Z_func: t (1B dizi) -> Z (1B dizi) (motor, ya da sahte motor --fake modunda).
# Vektörize: TÜM (ga,gb) çiftleri için AYNI ANDA.
# ---------------------------------------------------------------------------
def peak_grid_parabola(Z_func, ga, gb, n_grid=33, chunk=None):
    ga = np.asarray(ga, dtype=np.float64)
    gb = np.asarray(gb, dtype=np.float64)
    n_ev = len(ga)
    if n_ev == 0:
        return np.empty(0)
    # iç ızgara: k/(n_grid+1), k=1..n_grid (uç noktalar HARİÇ — sıfırlarda Z=0)
    frac = np.arange(1, n_grid + 1, dtype=np.float64) / (n_grid + 1)
    t_grid = ga[:, None] + (gb - ga)[:, None] * frac[None, :]   # (n_ev, n_grid)
    flat_t = t_grid.reshape(-1)

    if chunk is None:
        Z_flat = Z_func(flat_t)
    else:
        Z_flat = np.empty_like(flat_t)
        for s in range(0, len(flat_t), chunk):
            Z_flat[s:s + chunk] = Z_func(flat_t[s:s + chunk])
    absZ = np.abs(Z_flat).reshape(n_ev, n_grid)

    k = np.clip(np.argmax(absZ, axis=1), 1, n_grid - 2)
    r = np.arange(n_ev)
    y0, y1, y2 = absZ[r, k - 1], absZ[r, k], absZ[r, k + 1]
    den = y0 - 2.0 * y1 + y2
    safe_den = np.where(den == 0, -1e-300, den)
    M = y1 - 0.125 * (y2 - y0) ** 2 / safe_den
    # taban güvenliği: parabolik rafinasyon 3-nokta ızgara maksimumunun altına
    # düşmemeli (M9/tutarlılık: den>0 durumunda -parabol tepe yerine çukur- yakalarsa)
    M = np.maximum(M, y1)
    return M


# ---------------------------------------------------------------------------
# Genel (tek-rung, 2 istatistik: k2,k3) blok jackknife — 200a_ortak.jackknife_4x4'ün
# tek-rung sürümü. blocks: (B,5) [n,S1,S2,S3,S4].
# ---------------------------------------------------------------------------
def kstat_from_sums_1(n, s1, s2, s3):
    n = float(n)
    k2 = (n * s2 - s1 ** 2) / (n * (n - 1))
    k3 = (2 * s1 ** 3 - 3 * n * s1 * s2 + n ** 2 * s3) / (n * (n - 1) * (n - 2))
    return k2, k3


def jackknife_2(blocks):
    """blocks: (B,5) float64 [n,S1,S2,S3,S4] pencere-içi t-blokları. Döner:
    point (2,) [k2,k3], cov (2,2), reps (B,2)."""
    blocks = np.asarray(blocks, dtype=np.float64)
    B = blocks.shape[0]
    total = blocks.sum(axis=0)
    n, s1, s2, s3, s4 = total
    k2, k3 = kstat_from_sums_1(n, s1, s2, s3)
    point = np.array([k2, k3])

    reps = np.empty((B, 2))
    for i in range(B):
        loo = total - blocks[i]
        nn, ss1, ss2, ss3, ss4 = loo
        k2i, k3i = kstat_from_sums_1(nn, ss1, ss2, ss3)
        reps[i] = [k2i, k3i]

    rep_bar = reps.mean(axis=0)
    d = reps - rep_bar[None, :]
    cov = (B - 1) / B * (d.T @ d)
    return point, cov, reps


def block_index(t, t_a, t_b, nblocks):
    idx = ((t - t_a) / (t_b - t_a) * nblocks).astype(np.int64)
    return np.clip(idx, 0, nblocks - 1)


def power_sums_blocked(x, blk_idx, nblocks):
    """x,blk_idx (aynı uzunluk) -> (nblocks,5) [n,S1,S2,S3,S4]."""
    x = np.asarray(x, dtype=np.float64)
    counts = np.bincount(blk_idx, minlength=nblocks).astype(np.float64)
    s1 = np.bincount(blk_idx, weights=x, minlength=nblocks)
    s2 = np.bincount(blk_idx, weights=x ** 2, minlength=nblocks)
    s3 = np.bincount(blk_idx, weights=x ** 3, minlength=nblocks)
    s4 = np.bincount(blk_idx, weights=x ** 4, minlength=nblocks)
    return np.stack([counts, s1, s2, s3, s4], axis=1)


def merge_adjacent_blocks_1(block_sums, factor):
    block_sums = np.asarray(block_sums, dtype=np.float64)
    B = block_sums.shape[0]
    assert B % factor == 0
    return block_sums.reshape(B // factor, factor, 5).sum(axis=1)


# ---------------------------------------------------------------------------
# Doğrusal ara değer (N ızgarasında, M6b(i) sonlu-ε tablosu için)
# ---------------------------------------------------------------------------
def lin_interp_table(N_grid, val_grid, N_query):
    """N_grid artan; N_query'de doğrusal ara değer (uçlarda kelepçelenir/extrapolate
    edilmez -> sınırda tutulur, KALEM aralığı [9,12] içinde L̄'lerin hepsi zaten var)."""
    return float(np.interp(N_query, N_grid, val_grid))
