"""
200c_motor.py — KALEM 200-C: yüksek-t ÇAPALI Riemann-Siegel motoru, Z(t) ve Z'(t)
=====================================================================================

Görev metninden ("Task 1"):

  Zamanlar (tamsayı taban T0 + float64 ofset) olarak temsil edilir: t = T0 + t_ofset,
  T0 dosyanın (LMFDB/Platt) tabanı (tam sayı), t_ofset küçük (pencere genişliği
  mertebesinde, ~1e6'nın altında) — böylece t büyük olsa bile (t~3.1e10) faz
  hesaplarında kullanılan farklar (dt) her zaman küçük ve float64'te KESİN kalır.

  Faz (n'inci terim, t = t_a + dt, t_a bir ÇAPA):
    [theta(t_a) - t_a*log(n) mod 2pi]        (mpmath, YÜKSEK HASSASİYET, çapa başına BİR KEZ)
    + (theta'(t_a) - log n) * dt              (float64 — dt küçük, KESİN)
    + theta''(t_a) * dt^2 / 2                 (float64)

  Atılan kübik terim: theta'''(t_a)*dt^3/6 ~ -dt^3/(12*t_a^2) (theta'''(t)~-1/(2t^2)
  baskın terim) — HER değerlendirme noktasında |bu terim| <= 1e-10 rad olacak
  şekilde çapa aralığı (spacing) pencere başına HESAPLANIR (bkz anchor_spacing()),
  en küçük t (pencerenin alt ucu) kullanılarak MUHAFAZAKAR/worst-case.

  Kalan terim (C0+C1): 200a_motor.py'nin AYNI Chebyshev katsayılarıyla (DEĞİŞTİRİLMEDİ,
  yalnız İÇE AKTARILIR) — p = frac(sqrt(t/2pi)) çapa-göreli hesaplanır (p_base çapada
  mpmath ile, sonra dt ile float64 doğrusal kayma a_slope*dt).

  Z'(t): ana-toplamın ANALİTİK türevi (yukarıdaki faz açılımının dt'ye göre türevi
  zaten theta'(t_a)+theta''(t_a)*dt - log n verir — bu TUTARLI bir yaklaşım, çünkü
  fazın kendisi de aynı Taylor açılımını kullanıyor) + kalan terimin MERKEZİ FARK
  türevi (200a_motor ile aynı yöntem, h~1e-3, R(t) zaten C0+C1 toplamı).

Bu dosya SADECE motor tanımlarını içerir; hiçbir dosyaya yazma / hiçbir istatistik
hesabı YAPMAZ. Motor kördür (girdiye göre davranış değişmez); körlük M0c-M9c ve
ölçüm betiklerinde uygulanır.

Performans notu: çapa başına taban-faz dizisi (n=1..N_max_pencere) HESAPLANMASI
mpmath gerektirir (t_a büyük olduğunda t_a*log(n)'nin float64'te KESİN
hesaplanamaması — göreli hassasiyet t_a ile çarpılınca mutlak hassasiyeti
t_a*eps~t_a*2.2e-16 mertebesine düşürür, t_a~3e10'da ~7e-6 rad >> 1e-10 hedefi).
Kıyaslama (bu makinede, dps=30, saf-Python mpmath — gmpy2 YOK): ~4.4-5 us/çağrı
(mp.fmod(t*log(n),2pi) mertebesinde) — toplam ~1.8e7 (anchor,n) çifti tüm 4
pencerede ⇒ seri tahmini ~90-100 s; burada multiprocessing (çapa başına
BAĞIMSIZ) ile ek güvenlik payı için paralelleştirilir.
"""
import importlib.util
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import mpmath as mp

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
# DİKKAT (bulunan hata, düzeltildi): "MP_TWO_PI = 2*mp.pi" modül YÜKLENİRKEN bir
# KEZ hesaplanırsa, mp.pi o ANKİ (modül import sırasındaki, mpmath VARSAYILANI
# dps=15) hassasiyette DONAR — sonraki mp.mp.dps=30 ayarları bu donmuş sabiti
# ETKİLEMEZ. fmod(theta_a, 2*pi) gibi büyük-sayı indirgemelerinde bu, ~1e-5
# rad mertebesinde SİSTEMATİK faz hatasına yol açtı (teşhis: bkz
# scratchpad/k200c_makine/debug_c4_phase*.py). Çözüm: 2*pi HER ZAMAN
# mp.mp.dps ayarlandıktan SONRA, fonksiyon içinde taze hesaplanır (bkz
# _anchor_worker: `two_pi = 2 * mp.pi` dps set edildikten hemen sonra) —
# hiçbir modül-seviyesi dondurulmuş mpmath sabiti YOK.
DEFAULT_DPS = 30
DEFAULT_CUBIC_TOL = 1e-10


def _load(modname):
    """sys.modules'e KAYDEDER (multiprocessing 'spawn' ile pickle-by-reference'ın
    çapa işçi fonksiyonlarını doğru çözmesi için gerekli — aksi halde her
    dinamik-yükleme farklı bir modül nesnesi üretir ve Pool.map pickling'de
    'not the same object' hatası verir)."""
    if modname in sys.modules:
        return sys.modules[modname]
    spec = importlib.util.spec_from_file_location(modname, HERE / f"{modname}.py")
    m = importlib.util.module_from_spec(spec)
    sys.modules[modname] = m
    spec.loader.exec_module(m)
    return m


# 200a_motor.py DEĞİŞTİRİLMEDİ — yalnız C0/C1 Chebyshev değerlendirme fonksiyonları
# (_C0, _C1) İÇE AKTARILIYOR (KALEM Görev 1: "Remainder: C0 + C1 exactly as in
# 200a_motor, reuse its Chebyshev coefficients").
_motor_a = _load("200a_motor")
_C0 = _motor_a._C0
_C1 = _motor_a._C1


# ---------------------------------------------------------------------------
# theta(t), theta'(t), theta''(t) — float64, mutlak t (büyüklük/eğim terimleri
# için; bunlar YAVAŞ değişen düzgün fonksiyonlar, float64 göreli hassasiyeti
# (~1e-16) burada yeterli — yalnız FAZIN KENDİSİ (t*log n mertebesinde büyük
# sayılar) mpmath gerektiriyor, bkz modül docstring).
# ---------------------------------------------------------------------------
def theta_f64(t):
    t = np.asarray(t, dtype=np.float64)
    return t / 2 * np.log(t / TWO_PI) - t / 2 - np.pi / 8 + 1 / (48 * t) + 7 / (5760 * t ** 3)


def theta_prime_f64(t):
    t = np.asarray(t, dtype=np.float64)
    return 0.5 * np.log(t / TWO_PI) - 1 / (48 * t ** 2) - 7 / (1920 * t ** 4)


def theta_pprime_f64(t):
    """theta''(t) = d/dt[theta'(t)] = 1/(2t) + 1/(24 t^3) + 7/(480 t^5)
    (theta'(t) = 0.5 log(t/2pi) - 1/(48t^2) - 7/(1920 t^4)'ün türevi; elle
    doğrulandı, ayrıca __main__ duman testinde merkezi farka karşı sınandı)."""
    t = np.asarray(t, dtype=np.float64)
    return 1 / (2 * t) + 1 / (24 * t ** 3) + 7 / (480 * t ** 5)


def theta_tprime_f64(t):
    """theta'''(t) ~ baskın terim -1/(2t^2) (atılan kübik terimin büyüklüğünü
    raporlamak için; tam ifade -1/(2t^2) - 1/(8t^4) - 7/(96 t^6))."""
    t = np.asarray(t, dtype=np.float64)
    return -1 / (2 * t ** 2) - 1 / (8 * t ** 4) - 7 / (96 * t ** 6)


def anchor_spacing(t_min, cubic_tol=DEFAULT_CUBIC_TOL):
    """h: |theta'''(t_min)|*h^3/6 <= cubic_tol  =>  h ~ (12*tol*t_min^2)^(1/3)
    (theta'''(t)~-1/(2t^2) baskın terimiyle, worst-case t_min kullanılarak —
    pencerenin ALT UCUNDA en sıkı, t arttıkça daha gevşer). spacing = 2h,
    tamsayıya YUVARLANARAK AŞAĞI (güvenlik payı için) döner."""
    t_min = float(t_min)
    h = (12.0 * cubic_tol * t_min ** 2) ** (1.0 / 3.0)
    spacing = max(int(np.floor(2 * h)), 1)
    return spacing, h


# ---------------------------------------------------------------------------
# Çapa tablosu inşası (mpmath, çapa başına BİR KEZ)
# ---------------------------------------------------------------------------
def _anchor_worker(args):
    T0_int, offsets_chunk, N_max_window, dps = args
    mp.mp.dps = dps
    two_pi = 2 * mp.pi  # dps set edildikten SONRA taze hesaplanır (bkz modül üstü not)
    logn_mp = [None] + [mp.log(n) for n in range(1, N_max_window + 1)]  # 1-indeksli
    n_a = len(offsets_chunk)
    phase = np.empty((n_a, N_max_window), dtype=np.float64)
    thp = np.empty(n_a)
    thpp = np.empty(n_a)
    pbase = np.empty(n_a)
    aslope = np.empty(n_a)
    tln = np.empty(N_max_window)
    for i, off in enumerate(offsets_chunk):
        t_a = mp.mpf(T0_int) + mp.mpf(int(off))
        theta_a = t_a / 2 * mp.log(t_a / two_pi) - t_a / 2 - mp.pi / 8 + 1 / (48 * t_a) + 7 / (5760 * t_a ** 3)
        theta_a_mod = float(mp.fmod(theta_a, two_pi))
        for n in range(1, N_max_window + 1):
            tln[n - 1] = float(mp.fmod(t_a * logn_mp[n], two_pi))
        phase[i] = (theta_a_mod - tln) % TWO_PI

        a0m = mp.sqrt(t_a / two_pi)
        N_a = int(a0m)
        pbase[i] = float(a0m - N_a)
        aslope[i] = float(1 / (4 * mp.pi * a0m))
        t_a_f = float(t_a)
        thp[i] = 0.5 * np.log(t_a_f / TWO_PI) - 1 / (48 * t_a_f ** 2) - 7 / (1920 * t_a_f ** 4)
        thpp[i] = 1 / (2 * t_a_f) + 1 / (24 * t_a_f ** 3) + 7 / (480 * t_a_f ** 5)
    return phase, thp, thpp, pbase, aslope


class AnchorTable:
    __slots__ = ("T0_int", "offsets", "spacing", "phase", "thp", "thpp", "pbase",
                 "aslope", "logn", "N_max_window", "dps", "cubic_tol", "h",
                 "offset_min", "offset_max", "build_time_s")


def build_anchor_table(T0_int, offset_min, offset_max, cubic_tol=DEFAULT_CUBIC_TOL,
                        dps=DEFAULT_DPS, n_workers=8, margin=8):
    """Bir pencere için tam çapa tablosu. offset_min/offset_max: pencerenin
    sorgulanacak t_ofset aralığı (>=0). Çapalar [offset_min - spacing, offset_max
    + spacing] aralığını (kenar payıyla) kapsayacak şekilde eşit aralıklı
    tamsayı ofsetlerde kurulur."""
    t0 = time.time()
    t_min = T0_int + offset_min
    t_max = T0_int + offset_max
    spacing, h = anchor_spacing(t_min, cubic_tol)

    a0 = int(np.floor(offset_min / spacing)) * spacing - spacing
    a1 = int(np.ceil(offset_max / spacing)) * spacing + spacing
    offsets = np.arange(a0, a1 + spacing, spacing, dtype=np.float64)

    N_max_window = int(np.sqrt(t_max / TWO_PI)) + margin

    n_workers = max(1, min(n_workers, len(offsets)))
    chunks = np.array_split(offsets, n_workers)
    jobs = [(T0_int, ch, N_max_window, dps) for ch in chunks if len(ch) > 0]
    if n_workers > 1 and len(jobs) > 1:
        with Pool(len(jobs)) as pool:
            results = pool.map(_anchor_worker, jobs)
    else:
        results = [_anchor_worker(j) for j in jobs]

    phase = np.concatenate([r[0] for r in results], axis=0)
    thp = np.concatenate([r[1] for r in results])
    thpp = np.concatenate([r[2] for r in results])
    pbase = np.concatenate([r[3] for r in results])
    aslope = np.concatenate([r[4] for r in results])

    mp.mp.dps = dps
    logn = np.array([float(mp.log(n)) for n in range(1, N_max_window + 1)], dtype=np.float64)

    tab = AnchorTable()
    tab.T0_int = int(T0_int)
    tab.offsets = offsets
    tab.spacing = spacing
    tab.phase = phase
    tab.thp = thp
    tab.thpp = thpp
    tab.pbase = pbase
    tab.aslope = aslope
    tab.logn = logn
    tab.N_max_window = N_max_window
    tab.dps = dps
    tab.cubic_tol = cubic_tol
    tab.h = h
    tab.offset_min = offset_min
    tab.offset_max = offset_max
    tab.build_time_s = time.time() - t0
    return tab


# ---------------------------------------------------------------------------
# Değerlendirme: Z(t), Z'(t) — çapa tablosundan, vektörize, N'ye göre gruplanmış
# ---------------------------------------------------------------------------
def _default_chunk(N_max_window):
    return max(200, int(6.0e7 / max(N_max_window, 1)))


def Z_and_Zprime(t_offset, table, chunk=None, h_fd=1e-3):
    """t_offset: pencerenin T0'ından float64 ofset (1B dizi). table: build_anchor_table
    çıktısı. Döner: (Z, Zprime), t_offset ile aynı şekilde float64."""
    t_offset = np.asarray(t_offset, dtype=np.float64)
    if t_offset.ndim != 1:
        raise ValueError("t_offset bir boyutlu olmalı")
    if chunk is None:
        chunk = _default_chunk(table.N_max_window)

    out_Z = np.empty_like(t_offset)
    out_Zp = np.empty_like(t_offset)

    off0 = table.offsets[0]
    spacing = table.spacing
    n_anchors = len(table.offsets)

    for s in range(0, len(t_offset), chunk):
        to = t_offset[s:s + chunk]
        aidx = np.clip(np.round((to - off0) / spacing).astype(np.int64), 0, n_anchors - 1)
        dt = to - table.offsets[aidx]
        t_abs = table.T0_int + to
        b = t_abs / TWO_PI
        a = np.sqrt(b)
        N = a.astype(np.int64)

        thp_a = table.thp[aidx]
        thpp_a = table.thpp[aidx]
        pbase_a = table.pbase[aidx]
        aslope_a = table.aslope[aidx]

        Zm = np.zeros_like(to)
        Zpm = np.zeros_like(to)
        for Nv in np.unique(N):
            if Nv <= 0:
                continue
            m = N == Nv
            logn = table.logn[:Nv]
            wts = np.arange(1, Nv + 1, dtype=np.float64) ** -0.5
            base_m = table.phase[aidx[m], :Nv]
            dt_m = dt[m][:, None]
            ph = base_m + (thp_a[m][:, None] - logn[None, :]) * dt_m + thpp_a[m][:, None] * dt_m ** 2 / 2.0
            cosph = np.cos(ph)
            sinph = np.sin(ph)
            Zm[m] = 2.0 * (cosph @ wts)
            coef = thp_a[m][:, None] + thpp_a[m][:, None] * dt_m - logn[None, :]
            Zpm[m] = -2.0 * ((coef * sinph) @ wts)

        # kalan terim (C0+C1), merkezi fark türeviyle (200a_motor ile aynı yöntem)
        sign_c0 = (-1.0) ** (N - 1)

        def remainder(dt_eval):
            p = (pbase_a + aslope_a * dt_eval) % 1.0
            c0 = _C0(p)
            c1 = _C1(p)
            return sign_c0 * b ** -0.25 * (c0 + c1 * b ** -0.5)

        R = remainder(dt)
        Rp = (remainder(dt + h_fd) - remainder(dt - h_fd)) / (2 * h_fd)

        out_Z[s:s + chunk] = Zm + R
        out_Zp[s:s + chunk] = Zpm + Rp

    return out_Z, out_Zp


def Z_vec(t_offset, table, chunk=None):
    Z, _ = Z_and_Zprime(t_offset, table, chunk=chunk)
    return Z


def Zprime_vec(t_offset, table, chunk=None, h_fd=1e-3):
    _, Zp = Z_and_Zprime(t_offset, table, chunk=chunk, h_fd=h_fd)
    return Zp


if __name__ == "__main__":
    # Küçük duman testi: C1 penceresinin başındaki dar bir dilim, mpmath'e karşı.
    T0 = 8846000
    tab = build_anchor_table(T0, offset_min=0.0, offset_max=3000.0, n_workers=4)
    print(f"duman testi: spacing={tab.spacing} h={tab.h:.2f} n_anchors={len(tab.offsets)} "
          f"N_max_window={tab.N_max_window} build_time={tab.build_time_s:.2f}s")

    rng = np.random.default_rng(1)
    to = rng.uniform(0.0, 3000.0, 6)
    Z, Zp = Z_and_Zprime(to, tab)

    mp.mp.dps = 25
    for tv, zv, zpv in zip(to, Z, Zp):
        tm = mp.mpf(T0) + mp.mpf(tv)
        zmp = float(mp.siegelz(tm))
        zpmp = float(mp.siegelz(tm, derivative=1))
        print(f"t_ofs={tv:12.4f}  Z={zv:+.8f} (mp {zmp:+.8f}, d={zv-zmp:+.2e})  "
              f"Z'={zpv:+.6f} (mp {zpmp:+.6f}, d={zpv-zpmp:+.2e})")

    # theta'' formülü duman testi (merkezi farka karşı)
    tt = np.array([1e7, 1e9, 3e10])
    hh = 1.0
    num = (theta_prime_f64(tt + hh) - theta_prime_f64(tt - hh)) / (2 * hh)
    ana = theta_pprime_f64(tt)
    print("theta'' num vs analitik fark:", np.abs(num - ana))
