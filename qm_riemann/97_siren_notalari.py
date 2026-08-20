"""
97 — SİREN NOTALARI: χ₅ (KOMPLEKS) VE χ₃ KIRINIMLARI (21 Ağustos, gece)
==========================================================================
Not 5 seferi: kompleks karakterde FAZ-KADRANI testi + üçüncü sıcaklık.

ÖNCEDEN MÜHÜRLENEN ÖNGÖRÜLER (ölçümden önce):
  χ₅ (mod 5, χ(2)=i; tek karakter, a=1, iletken 5):
    P1  5 ve 25 benekleri ÖLÜ (χ=0).
    P2  DÖRT-KONUM KADRANI (çapa: χ=+1 → 180°): kalıntı sınıfına göre
        q≡1 → 180° (11) | q≡4 → 0° (4, 9, 19) |
        q≡2 → 270° (2, 7, 17) | q≡3 → 90° (3, 8, 13, 23)
        (küresel yönelim ±arg χ konvansiyonuna bağlı: keskin iddia,
        dört sınıfın DÖRT AYRI konuma derece-düzeyinde kümelenmesi).
    P3  Genlikler mutlak yasada (L = log(5t/2π), tarak-DW);
        kuvvetler (4, 8, 9) kuvvet-açığı bandında (~0.75-0.9).
  χ₃ (mod 3, χ(2)=−1; tek, gerçek, iletken 3):
    P4  3, 9, 27 ÖLÜ; 2,5,8,11 → 0°; 4,7,13 → 180°.
  P5  SICAKLIK EVRENSELLİĞİ: her iki adanın σ_u'su (tarak-termometre)
      zeta'nın log-ısınma eğrisiyle tek çizgide (ζ pencereleri + β +
      χ₃ + χ₅ → eğrinin ilk çok-L hali).

Makine: 96'nın tamamlanmış-fonksiyon yöntemi; kompleks χ için kök-sayısı
fazı Gauss toplamından + ampirik faz-kalibrasyonu (çifte sağlama).
"""

import numpy as np
import mpmath as mp
from pathlib import Path

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
mp.mp.dps = 15

def make_L(q_cond, chi_table):
    """chi_table: {a mod q: kompleks değer}; L(s,χ) Hurwitz toplamı."""
    ays = [(a, chi_table[a]) for a in range(1, q_cond) if chi_table[a] != 0]
    def L(s):
        tot = mp.mpc(0)
        for a, c in ays:
            tot += c * mp.zeta(s, mp.mpf(a) / q_cond)
        return mp.power(q_cond, -s) * tot
    return L

def kok_fazi(q_cond, chi_table):
    tau = mp.mpc(0)
    for a in range(1, q_cond):
        tau += chi_table[a] * mp.e**(2j * mp.pi * a / q_cond)
    eps = tau / (1j * mp.sqrt(q_cond))
    return float(mp.arg(eps)), abs(complex(eps))

def sifir_uret(q_cond, chi_table, T0, T1, etiket):
    cache = HERE / f"97_{etiket}_zeros.npz"
    if cache.exists():
        return np.load(cache)["zeros"]
    Lf = make_L(q_cond, chi_table)
    aeps, meps = kok_fazi(q_cond, chi_table)
    print(f"[{etiket}] kök-sayısı: |ε| = {meps:.6f} (1 olmalı), arg ε = {aeps:.4f}")
    def Z(t):
        s = mp.mpc(0.5, t)
        psi = (t / 2) * mp.log(q_cond / mp.pi) + mp.im(mp.loggamma((s + 1) / 2))
        z = mp.e**(1j * (psi - aeps / 2)) * Lf(s)
        return float(mp.re(z)), float(mp.im(z))
    # ampirik faz sağlaması
    ims = [abs(Z(t)[1]) for t in np.linspace(T0 + 3, T0 + 60, 12)]
    print(f"[{etiket}] faz sağlaması: maks|Im Z| örneklem = {max(ims):.2e}")
    ts = []
    t = T0
    while t < T1:
        ts.append(t)
        t += 0.25 * TWO_PI / np.log(q_cond * t / TWO_PI)
    print(f"[{etiket}] ızgara {len(ts)} nokta — değerlendiriliyor...")
    vals = np.empty(len(ts))
    for i, tt in enumerate(ts):
        vals[i], _ = Z(float(tt))
        if i % 3000 == 0:
            print(f"  [{etiket}] {i}/{len(ts)}")
    sc = np.where(np.sign(vals[:-1]) * np.sign(vals[1:]) < 0)[0]
    print(f"[{etiket}] {len(sc)} işaret değişimi — inceltiliyor...")
    zeros = []
    for i in sc:
        a, b = float(ts[i]), float(ts[i + 1])
        fa = vals[i]
        for _ in range(22):
            m = 0.5 * (a + b)
            fm, _ = Z(m)
            if fa * fm <= 0:
                b = m
            else:
                a, fa = m, fm
        zeros.append(0.5 * (a + b))
    zz = np.array(zeros)
    np.savez(cache, zeros=zz)
    print(f"[{etiket}] {len(zz)} sıfır önbelleğe alındı")
    return zz

def Ghat(t, omegas, chunk=30000):
    out = np.zeros(len(omegas), dtype=complex)
    for s0 in range(0, len(t), chunk):
        tt = t[s0:s0 + chunk]
        out += np.exp(1j * np.outer(omegas, tt)).sum(axis=1)
    return out / len(t)

def analiz(zz, q_cond, chi_table, LAMk, etiket):
    gaps = np.diff(zz)
    mids = 0.5 * (zz[:-1] + zz[1:])
    n = len(mids)
    Leff = float(np.log(q_cond * mids / TWO_PI).mean())
    def Nsm(t):
        x = t / TWO_PI
        return x * np.log(q_cond * x / np.e)
    x_unf = Nsm(mids) - Nsm(zz[0])
    comb1 = abs(np.exp(2j * np.pi * x_unf).mean())
    c_jit = -np.log(comb1)
    sig_u = np.sqrt(2 * c_jit) / TWO_PI
    sig_t = sig_u * TWO_PI / np.log(q_cond * mids.mean() / TWO_PI)
    print(f"\n===== {etiket}: n={n}, L_eff={Leff:.3f} =====")
    print(f"tarak: |Ĝ_x(2π)| = {comb1:.4f} → σ_u = {sig_u:.3f}  ← SICAKLIK NOKTASI")
    print(f"{'q':>3} {'χ-sınıf':>8} {'ölçüm|Ĝ|':>9} {'öngörü':>8} {'oran':>6} "
          f"{'faz':>8} {'öngörü-faz':>10}")
    for q, lam, kk in LAMk:
        om = np.log(q)
        tau = om / Leff
        G = Ghat(mids, np.array([om]))[0]
        chi = chi_table[q % q_cond]
        ph = (np.degrees(np.angle(G)) + 360) % 360
        if chi == 0:
            print(f"{q:>3} {'0 (ölü)':>8} {abs(G):>9.4f} {'YOK':>8} {'—':>6} "
                  f"{ph:>7.0f}° {'(taban)':>10}")
        else:
            dw = np.exp(-om**2 * sig_t**2 / 2)
            pred = lam / (Leff * np.sqrt(q)) * np.cos(np.pi * tau) * dw
            pred_ph = (180 + np.degrees(np.angle(complex(chi)))) % 360
            print(f"{q:>3} {str(np.round(complex(chi),2)):>8} {abs(G):>9.4f} "
                  f"{pred:>8.4f} {abs(G)/pred:>6.3f} {ph:>7.1f}° {pred_ph:>9.0f}°")
    fakes = np.array([np.log(v) for v in [2.31, 6.7, 10.4]])
    print(f"plasebo ort = {np.abs(Ghat(mids, fakes)).mean():.4f};  "
          f"1/√n = {1/np.sqrt(n):.4f}")
    return Leff, sig_u

# ---- χ₅: χ(2)=i (2 üreteç: 2→i, 4→−1, 3→−i, 1→1)
CHI5 = {0: 0, 1: 1, 2: 1j, 3: -1j, 4: -1}
LAMK5 = [(5, np.log(5), 1), (25, np.log(5), 2),
         (2, np.log(2), 1), (7, np.log(7), 1), (17, np.log(17), 1),
         (3, np.log(3), 1), (8, np.log(2), 3), (13, np.log(13), 1),
         (23, np.log(23), 1), (4, np.log(2), 2), (9, np.log(3), 2),
         (19, np.log(19), 1), (11, np.log(11), 1)]
zz5 = sifir_uret(5, CHI5, 200.0, 4200.0, "chi5")

# ---- χ₃: χ(2)=−1
CHI3 = {0: 0, 1: 1, 2: -1}
LAMK3 = [(3, np.log(3), 1), (9, np.log(3), 2),
         (2, np.log(2), 1), (5, np.log(5), 1), (8, np.log(2), 3),
         (11, np.log(11), 1), (4, np.log(2), 2), (7, np.log(7), 1),
         (13, np.log(13), 1)]
zz3 = sifir_uret(3, CHI3, 200.0, 4600.0, "chi3")

L5, s5 = analiz(zz5, 5, CHI5, LAMK5, "χ₅ ADASI (kompleks — kadran testi)")
L3, s3 = analiz(zz3, 3, CHI3, LAMK3, "χ₃ ADASI")

print("\n===== SICAKLIK EVRENSELLİĞİ (P5) =====")
print("L_eff   σ_u    kaynak")
print(f"{L3:5.2f}  {s3:.3f}  χ₃")
print(f"7.20   0.161  β (96)")
print(f"{L5:5.2f}  {s5:.3f}  χ₅")
print("9.86   0.198  ζ (76-ölçümü)   12.45  0.272  ζ")
print("→ dört L-fonksiyonu tek log-ısınma eğrisinde mi? (rapor)")
