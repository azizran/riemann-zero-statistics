"""
200c_ortak.py — KALEM 200-C: ortak yardımcılar (pencere/dosya tanımları, sıfırların
yüklenmesi, olay çıkarımı, çapa tablosu önbelleği). 200c_m*.py, 200c_olcum.py,
200c_analiz.py tarafından ortak kullanılır. Hiçbir gerçek Z/Z'/M dağılım özeti
burada hesaplanmaz — yalnız pencere/dosya META VERİSİ ve sayısal yardımcı
fonksiyonlar (200a_ortak, 200b_ortak, 200c_motor, 200c_veri İÇE AKTARILIR,
DEĞİŞTİRİLMEZ).
"""
import hashlib
import importlib.util
import sys
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
TWO_PI = 2 * np.pi


def _load(modname, where=HERE):
    if modname in sys.modules:
        return sys.modules[modname]
    spec = importlib.util.spec_from_file_location(modname, where / f"{modname}.py")
    m = importlib.util.module_from_spec(spec)
    sys.modules[modname] = m
    spec.loader.exec_module(m)
    return m


ortak_a = _load("200a_ortak")
ortak_b = _load("200b_ortak")
motor = _load("200c_motor")
veri = _load("200c_veri")
kap = ortak_a.kap_vec  # kap(r, N, b) — vektörize, b=0,1,2 hepsi burada kullanılacak

WINDOWS = ("C1", "C2", "C3", "C4")
FILES = {
    "C1": "zeros_8846000.dat",
    "C2": "zeros_99146000.dat",
    "C3": "zeros_997946000.dat",
    "C4": "zeros_30599546000.dat",
}
L_NOMINAL = {"C1": 14.194, "C2": 16.577, "C3": 18.884, "C4": 22.306}
N_MAX = 1_500_000
NBLOCKS = 128
EPS_B2_SAVE = 0.3     # olcum.py'de saklanan en genis esik (analiz.py sonra kirpar)
A2_ARITH = ortak_b.A2_ARITH   # -0.088124
A3_ARITH = ortak_b.A3_ARITH   # +0.233653


def data_path(w):
    return ROOT / "veri_lmfdb" / FILES[w]


def sha256_of(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ---------------------------------------------------------------------------
# Pencere yükleme — YALNIZ KONUM (T0 tam sayı taban, off float64 ofsetler).
# ---------------------------------------------------------------------------
class Window:
    __slots__ = ("name", "path", "T0", "off", "n_read", "bloklar", "t_a_off", "t_b_off")


def load_window(w, n_max=N_MAX):
    path = data_path(w)
    taban, off, bloklar = veri.oku(str(path), n_max=n_max)
    win = Window()
    win.name = w
    win.path = path
    win.T0 = int(taban)
    win.off = off
    win.n_read = len(off)
    win.bloklar = bloklar
    win.t_a_off = float(off[0])
    win.t_b_off = float(off[-1])
    return win


def L_of_offset(win, off_arr):
    """L(t) = log(t/2pi), t = T0+off (float64 — mutlak buyukluk icin guvenli,
    bkz 200c_motor modul notu: mutlak t'nin float64 hassasiyeti ~1e-6, log'un
    hassasiyetine etkisi ihmal edilebilir)."""
    t_abs = win.T0 + np.asarray(off_arr, dtype=np.float64)
    return np.log(t_abs / TWO_PI)


# ---------------------------------------------------------------------------
# Çapa tablosu (motor) — pencere başına ÖNBELLEKLENİR (tekrar tekrar inşa
# edilmesin diye; her M/olcum betiği kendi sureci icinde bir kez kurar).
# ---------------------------------------------------------------------------
_ANCHOR_CACHE = {}


def get_anchor_table(win, n_workers=8, margin=8):
    key = (win.name, n_workers)
    if key in _ANCHOR_CACHE:
        return _ANCHOR_CACHE[key]
    tab = motor.build_anchor_table(win.T0, win.t_a_off, win.t_b_off, n_workers=n_workers, margin=margin)
    _ANCHOR_CACHE[key] = tab
    return tab


def Z_func_for(win, tab, chunk=None):
    return lambda t_off: motor.Z_vec(t_off, tab, chunk=chunk)


# ---------------------------------------------------------------------------
# Olay (yakın çift) çıkarımı — YALNIZ KONUM (200b_ortak.events_from_zeros'un
# T0+ofset temsiline uyarlanmış sürümü; L, absolute t üzerinden hesaplanır).
# ---------------------------------------------------------------------------
def events_from_window(win, eps_max=EPS_B2_SAVE):
    off = win.off
    ga = off[:-1]
    gb = off[1:]
    n_idx = np.arange(len(off) - 1, dtype=np.int64)  # dosya-ici 0-tabanli indeks
    m = 0.5 * (ga + gb)
    L = L_of_offset(win, m)
    delta_off = gb - ga  # ofsetlerin farki = mutlak farkla AYNI (T0 sabit)
    delta_tilde = delta_off * L / TWO_PI
    mask = delta_tilde < eps_max
    return {"n": n_idx[mask], "ga": ga[mask], "gb": gb[mask], "m": m[mask], "L": L[mask],
            "delta": delta_off[mask], "delta_tilde": delta_tilde[mask]}


# ---------------------------------------------------------------------------
# Ortak t-blokları (128 esit-t-uzunluklu, pencere t_a_off..t_b_off araliginda)
# ---------------------------------------------------------------------------
def block_index(t_off, win, nblocks=NBLOCKS):
    idx = ((t_off - win.t_a_off) / (win.t_b_off - win.t_a_off) * nblocks).astype(np.int64)
    return np.clip(idx, 0, nblocks - 1)


def power_sums_blocked(x, blk_idx, nblocks=NBLOCKS):
    return ortak_b.power_sums_blocked(x, blk_idx, nblocks)


def load_tahmin_200c(path=None):
    path = path or (HERE / "200c_tahmin.json")
    with open(path) as f:
        return json.load(f)
