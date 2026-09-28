"""
200a_motor.py — KALEM 200-A: Motor B (vektörize Riemann–Siegel, float64) + Z'(t)
====================================================================================

Z(t) = 2 * Sum_{n<=N} n^{-1/2} cos(theta(t) - t log n) + R(t),   N = floor(sqrt(t/2pi))
R(t) = (-1)^(N-1) (t/2pi)^(-1/4) [C0(p) + C1(p)*(t/2pi)^(-1/2)],  p = frac(sqrt(t/2pi))
       C0(p) = cos(2*pi*(p^2-p-1/16)) / cos(2*pi*p)
       C1(p) = -C0'''(p) / (96*pi^2)
       (C1 EKLENDİ 26 Eyl 2026: W1 penceresinin alt ucunda M1 kapısı C0-yalnız motorla
       KALDI — |ΔZ| max 6.2e-5 > 2e-5 eşiği, bkz 200a_MAKINE_RAPORU.md §"pre-seal
       deviation". Koordinatörün talimatıyla doğrulanmış C1 terimi eklendi.)

C0(p), C1(p) p=1/4, 3/4'te KALDIRILABİLİR (removable) tekilliğe sahip (num/den
ikisi de sıfıra gider, oran sonlu bir limite yakınsar — bkz aşağıdaki YÜKSEK
HASSASİYET NOTU). Float64'te (ve hatta mpmath'ta TAM p=1/4 ya da 3/4'te, ör.
dps=40'ta bile) doğrudan cos(..)/cos(..) formülünü hesaplamak KATASTROFİK İPTAL
(catastrophic cancellation) verir — inşa sırasında doğrulandı: p=0.75 TAM'da
mpmath dps=40 ile num~2e-43, den~2e-41, oran=0.0093 (YANLIŞ; doğrusu 0.5 limit).
Bu yüzden C0, C1 burada doğrudan HESAPLANMAZ: p in [0,1] üzerinde, mpmath dps=40
ile (tekilliklere TAM denk düşmeyecek şekilde küçük kaydırmalarla) üretilmiş,
doğrulanmış (yoğun ızgarada max|hata| < 2e-15, hedef 1e-12'nin çok altında)
derece-28 Chebyshev polinom katsayılarıyla (CHEB_C0, CHEB_C1 aşağıda) SABİT
olarak saklanır ve numpy.polynomial.chebyshev.chebval ile değerlendirilir —
üretim kodu hiçbir zaman 0/0 formunu görmez (bkz
scratchpad/k200a_makine/build_c0_c1_chebyshev.py, bir-kerelik üretici).

Z'(t) = -2 * Sum_{n<=N} n^{-1/2} (theta'(t) - log n) sin(theta(t) - t log n) + R'(t)
theta(t)  = (t/2) log(t/2pi) - t/2 - pi/8 + 1/(48t) + 7/(5760 t^3)
theta'(t) = (1/2) log(t/2pi) - 1/(48 t^2) - 7/(1920 t^4)
R'(t): merkezi sonlu fark, yalnız kalan terimin kendisi üzerinde (h ~ 1e-3):
       R'(t) ~ [R(t+h) - R(t-h)] / (2h)   (C1'in katkısı otomatik dahil olur,
       çünkü R(t) zaten C0+C1 toplamı).

Bu dosya SADECE motor tanımlarını içerir; hiçbir dosyaya yazma / hiçbir istatistik
hesabı YAPMAZ. t argümanı olarak gerçek sıfır konumları da, sentetik t de verilebilir —
motor kördür (girdiye göre davranış değişmez); körlük M1/M2/olcum betiklerinde
uygulanır (bu dosyada hiçbir gerçek veri kullanılmaz / okunmaz).
"""
import numpy as np

TWO_PI = 2 * np.pi

# Derece-28 Chebyshev katsayıları, p in [0,1] (x = 2p-1 in [-1,1]) üzerinde C0(p), C1(p).
# Üretim: scratchpad/k200a_makine/build_c0_c1_chebyshev.py (mpmath dps=40, Lobatto
# düğümleri; p=1/4,3/4'e tam denk düşen düğümler 1e-7 kaydırıldı). Doğrulama: 2001+8
# noktalık yoğun ızgarada (0.25/0.75'e 1e-6 ve 1e-9 yakın noktalar dahil) mpmath'a
# karşı max|hata C0|=1.887e-15, max|hata C1|=1.631e-16 (hedef < 1e-12, rahatça geçti).
CHEB_C0 = np.array([
    0.6426672862397692, -2.58235202249818e-17, 0.27197299999785546, -7.250531671243924e-17,
    0.01073860581934049, 9.147388699029754e-18, -0.001374381529633682, 6.808888670445681e-18,
    -0.00012468221880335763, 3.239169947603857e-17, -5.764599703287179e-07, 2.016929676023804e-16,
    2.728067430441586e-07, -1.380168966982116e-16, 8.077952904788358e-09, 5.4703632024314076e-17,
    -2.088461313955312e-10, 1.099785003465083e-16, -1.3115522807486417e-11, 6.265479655871297e-17,
    -1.4103940039395554e-14, 1.0513393443275229e-16, 1.0170682753359706e-14, -3.7641387286761194e-17,
    2.1216060462161778e-16, 7.920553719517148e-18, -1.1529845140177923e-17, 1.352186862200935e-16,
    -4.7326498504079754e-17,
])
CHEB_C1 = np.array([
    7.876222252623142e-19, 0.010697913921003017, 1.6140147479705943e-17, 0.01717065124337789,
    1.6761028241798142e-17, 0.002793211149788478, 2.314129791151925e-18, -3.6375653719265965e-05,
    -6.047103174638284e-19, -2.7108955231160883e-05, 2.1675319090777567e-17, -1.0483749866510241e-06,
    -1.7408563499797703e-18, 5.886467164373442e-08, -1.7915795287294143e-17, 4.3229672572099e-09,
    -1.9825267591548533e-17, -1.1369627126060306e-11, -3.0565935306960177e-18, -6.699839200700353e-12,
    -7.120777572929532e-18, -1.0079587468841147e-13, -6.370485752314593e-18, 5.14758557014337e-15,
    2.4399847195196155e-18, 1.4926756051141325e-16, 4.23233910260304e-18, -4.328965908666393e-18,
    -9.608540811497415e-19,
])


def _C0(p):
    x = 2.0 * p - 1.0
    return np.polynomial.chebyshev.chebval(x, CHEB_C0)


def _C1(p):
    x = 2.0 * p - 1.0
    return np.polynomial.chebyshev.chebval(x, CHEB_C1)


def theta(t):
    t = np.asarray(t, dtype=np.float64)
    return t / 2 * np.log(t / TWO_PI) - t / 2 - np.pi / 8 + 1 / (48 * t) + 7 / (5760 * t ** 3)


def theta_prime(t):
    t = np.asarray(t, dtype=np.float64)
    return 0.5 * np.log(t / TWO_PI) - 1 / (48 * t ** 2) - 7 / (1920 * t ** 4)


def _remainder(t):
    """C0+C1 düzeltme terimi R(t) (kendi başına, ana toplamdan bağımsız).
    R(t) = (-1)^(N-1) (t/2pi)^(-1/4) [C0(p) + C1(p)*(t/2pi)^(-1/2)]
    C0, C1 Chebyshev-polinom değerlendirmesiyle (bkz üstteki not) — hiçbir 0/0
    formu YOK, tekillik koruması gerekmiyor (polinom her yerde düzgün)."""
    t = np.asarray(t, dtype=np.float64)
    b = t / TWO_PI
    a = np.sqrt(b)
    N = a.astype(np.int64)
    p = a - N
    c0 = _C0(p)
    c1 = _C1(p)
    return ((-1.0) ** (N - 1)) * b ** -0.25 * (c0 + c1 * b ** -0.5)


def _default_chunk(t):
    Nmax = int(np.sqrt(float(np.max(t)) / TWO_PI)) + 1
    return max(200, int(8.0e7 / max(Nmax, 1)))


def Z_and_Zprime(t, chunk=None, h=1e-3):
    """Vektörize Z(t) ve Z'(t), bellek sınırlamak için parçalar (chunk) hâlinde.

    Girdi: t (1B dizi, pozitif, [1.8e4, 1.2e6] aralığında güvenli).
    Çıktı: (Z, Zprime) — her ikisi de t ile aynı şekilde float64 dizi.
    """
    t = np.asarray(t, dtype=np.float64)
    if t.ndim != 1:
        raise ValueError("t bir boyutlu dizi olmalı")
    if chunk is None:
        chunk = _default_chunk(t)
    out_Z = np.empty_like(t)
    out_Zp = np.empty_like(t)
    for s in range(0, len(t), chunk):
        tt = t[s:s + chunk]
        a = np.sqrt(tt / TWO_PI)
        N = a.astype(np.int64)
        th = theta(tt)
        thp = theta_prime(tt)
        Zm = np.zeros_like(tt)
        Zpm = np.zeros_like(tt)
        for Nv in np.unique(N):
            m = N == Nv
            if Nv <= 0:
                continue
            n = np.arange(1, Nv + 1)
            logn = np.log(n)
            wts = n ** -0.5
            ph = th[m, None] - tt[m, None] * logn[None, :]
            cosph = np.cos(ph)
            sinph = np.sin(ph)
            Zm[m] = 2.0 * (cosph @ wts)
            coef = thp[m, None] - logn[None, :]
            Zpm[m] = -2.0 * ((coef * sinph) @ wts)
        R = _remainder(tt)
        Rp = (_remainder(tt + h) - _remainder(tt - h)) / (2 * h)
        out_Z[s:s + chunk] = Zm + R
        out_Zp[s:s + chunk] = Zpm + Rp
    return out_Z, out_Zp


def Z_vec(t, chunk=None):
    """Yalnız Z(t) (rung-0 için biraz daha ucuz yol: yine de Z' hesaplanır ve atılır,
    kod sadeliği için; maliyeti kabul edilebilir, bkz. MAKINE_RAPORU)."""
    Z, _ = Z_and_Zprime(t, chunk=chunk)
    return Z


def Zprime_vec(t, chunk=None, h=1e-3):
    _, Zp = Z_and_Zprime(t, chunk=chunk, h=h)
    return Zp


if __name__ == "__main__":
    # Küçük, kendi kendine yeten duman testi (sentetik t; gerçek veri YOK).
    rng = np.random.default_rng(12345)
    ts = rng.uniform(2.0e4, 1.1e6, 5)
    Z, Zp = Z_and_Zprime(ts)
    for tv, zv, zpv in zip(ts, Z, Zp):
        print(f"t={tv:14.3f}  Z={zv:+.6f}  Z'={zpv:+.6f}")
    print("motor duman testi tamam (mpmath karşılaştırması 200a_kapilar_m1_m2.py'de).")
