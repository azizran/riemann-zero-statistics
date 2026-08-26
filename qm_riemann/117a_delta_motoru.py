"""
117a — RAMANUJAN Δ: τ(n) İNŞASI, GL(2) MOTORU VE SAĞLAMA KAPILARI
        (26 Ağustos 2026 — programın İLK GL(2) örgüsü, öncü sefer)
==========================================================================
ÖN-MÜHÜR (ölçümden ÖNCE yazıldı)
==========================================================================
Bugüne kadarki bütün adalar GL(1)'di: ζ ve Dirichlet L(s,χ). Bu sefer
derece-2'ye geçiyoruz: Δ = q·Π(1−q^n)^24 (seviye 1, ağırlık 12) —
Ramanujan tau formu. Amaç, kırınım yasalarının (mutlak benek yasası,
faz okuma, tarak-termometresi, donma) GL(1)'e özgü olup olmadığını
sınamak; ve GL(2)'ye özgü İKİ YENİ okumayı denemek:
  (i)  Hecke İŞARETLERİ: a_p = τ(p)/p^{11/2} işaretini kırınım FAZINDAN
       okumak (ζ'da bütün fazlar 180°'ydi; burada işaret değiştirmeli),
  (ii) SATAKE p²-ÇİZGİLERİ: açık formül katsayısı α_p²+β_p² = a_p²−2,
       küçük asallarda NEGATİF ⟹ p²-benekleri ζ'ya göre FAZ-ÇEVRİK.

KURAM / KONVANSİYON (117a'nın sağlayacağı şey budur)
  Λ(s) = (2π)^{-s} Γ(s) L(s,Δ),  Λ(s) = Λ(12−s),  ε = +1  (k/2 = 6 çift)
  kritik doğru Re s = 6;  θ(t) = −t·log(2π) + Im logΓ(6+it)
  Z(t) = e^{iθ(t)}·L(6+it)/(normalizasyon) GERÇEK;
  ana toplam (AFE, derece-2 uzunluğu — analitik iletken (t/2π)²,
  karekökü t/2π):
      Z(t) ≈ 2 Σ_{n ≤ X(t)} (a(n)/√n)·cos(θ(t) − t·log n),  X(t) = t/2π
      a(n)/√n = τ(n)/n^6   (analitik normalizasyon a(n) = τ(n)/n^{11/2})
  Yoğunluk: θ'(t)/π = (1/π)·log(t/2π) ⟹ L_eff = 2π·yoğunluk = 2·log(t/2π).
  (Kaptanın mühründeki "yoğunluk ≈ (2/π)log(t/2π)" ifadesi 2 kat fazla;
   L_eff = 2log(t/2π) konvansiyonu ile tutarlı olan (1/π)log(t/2π)'dir.
   N(T) = (T/π)·log(T/2πe) + O(log T) — T=2e4'te ~45k sıfır, 130k değil.)

KAPILAR (geçilmezse KAMPANYAYA GİRİLMEZ)
  G1  τ(n) KESİN: python-int polinom çarpımı, E³ = Σ(−1)^k(2k+1)q^{k(k+1)/2}
      (Jacobi) seyrek çekirdeğiyle 7 kez çarpım ⟹ E^24.
      τ(2)=−24, τ(3)=252, τ(5)=4830, τ(7)=−16744, τ(11)=534612;
      çarpımsallık τ(6)=τ(2)τ(3), τ(10)=τ(2)τ(5), τ(15)=τ(3)τ(5);
      Hecke özyinelemesi τ(p²)=τ(p)²−p^11·τ(1);
      DELIGNE |a(p)| ≤ 2 bütün asallarda.
  G2  θ: vektörize Stirling vs mpmath.loggamma (t = 200 … 20000).
  G3  BAĞIMSIZ ORACLE (bu seferin en güçlü kapısı): Hecke'nin KESİN
      integral formülü
        Λ(6+it) = Σ_n τ(n)[Γ(6+it,2πn)/(2πn)^{6+it}
                            + Γ(6−it,2πn)/(2πn)^{6−it}]
      mpmath ile YÜKSEK HASSASİYETTE (dps ≈ 40+0.8t; e^{πt/2} iptali
      yüzünden t≲250'de kullanılabilir). Bu, motorla AYNI FORMÜL DEĞİL —
      gerçek bir bağımsız değerlendirme. Sınanan: (a) Λ'nın GERÇEK
      olması (ε=+1 ve FE), (b) motor-Z vs kesin-Z, (c) ilk sıfırlar,
      (d) γ₁ ≈ 9.22237939992 (LMFDB 2-1-1.1-c11-0-0; hafıza-temelli,
      düşük ağırlıklı) .
  S1  KESİM DUYARLILIĞI/FE-TUTARLILIĞI: keskin kesim vs yumuşak (taper)
      kesim vs 1.3X — kampanya yüksekliklerinde |ΔZ|/RMS|Z|.
  S2  SAYIM SERTİFİKASI: d_i = i − (θ(z_i)−θ(z_0))/π sürüklenmesiz mi
      (üç yükseklikte pilot pencere).
  S3  SIFIR KARARLILIĞI: keskin vs taper motorunun sıfırları
      |Δγ|/⟨g⟩ olarak.
  N1  NEGATİF KONTROL (kapı duyarlı mı?): (a) yanlış ε (θ→θ+π/2),
      (b) derece-1 uzunluğu X=√(t/2π), (c) yanlış normalizasyon n^{-5.5}.
      Bunların sayım sertifikasını PATLATMASI beklenir; patlatmıyorsa
      kapı kör demektir ve kampanya yine iptal.

BEKLENTİ (mühür): G1-G3 geçer; motorun RS-düzeltmesiz hatası
GL(1)'dekiyle aynı sınıfta (|ΔZ| ≪ RMS|Z|) olur; kaldırılmış dip
sorunu GL(1)'dekinden DAHA CİDDİ olabilir (derece 2, daha uzun toplam).

==========================================================================
SONUÇ (26 Ağustos — TÜM KAPILAR GEÇİLDİ; iki ölçüt düzeltmesi kayıtlı)
==========================================================================
G1 ✓ τ(n), n ≤ 6000, 0.31 sn: beş referans değer, beş çarpımsallık,
   altı Hecke p² özyinelemesi TAM İSABET; DELIGNE |a(p)| ≤ 2 — 783
   asalın hepsinde geçerli, maks 1.9514 @ p=3967.
G2 ✓ |Δθ| ≤ 3.6e-12 (t = 200 … 40000).
G3 ✓✓✓ BAĞIMSIZ ORACLE — seferin en güçlü kapısı:
   (a) Λ(6+it) mpmath'te TAM GERÇEK (|Im/Re| ≡ 0) ⟹ ε=+1 ve FE
       doğrulandı — motorun formülünden BAĞIMSIZ.
   (b) motor-Z vs kesin-Z (21 nokta, t∈[45,185], RMS|Z|=1.553):
       KESKİN kesim ort|Δ| = 0.0913 (%5.88)
       TAPER c=0.5  ort|Δ| = 0.0080 (%0.51)   ⟹ 11.4 KAT kazanç.
       (Taper genişliği bir ÖN-TARAMAyla seçildi: c ∈ {0,0.2,…,1.5},
        25 oracle noktası; ort|Δ| c=0.5'te keskin minimum yapıyor
        (0.084 → 0.0059) ve optimum t'nin İKİ YARISINDA DA 0.5
        (0.0049 ve 0.0068) — W = 0.5√X Fresnel genişliği yapısal
        görünüyor, tek bir örneklemin uydurması değil.)
   (b2) KESİM UZUNLUĞU ÖLÇÜLDÜ: X = u·t/2π için oracle'a karşı
       u=0.7 → 0.396,  u=1.0 → 0.008,  u=1.3 → 0.289.
       Derece-2 AFE uzunluğu t/2π'de KEskin bir optimum — analitik
       iletkenin (t/2π)² olduğunun doğrudan ölçümü.
   (c) 5 sıfır oracle-bisection'ıyla: |Δγ|/⟨g⟩ = 5.8e-5 … 2.6e-3.
   (d) ORACLE'IN KENDİ ilk sıfırları (motor devre dışı):
       9.22237940  13.90754986  17.44277698  19.65651314  22.33610364
       γ₁ vs literatür (LMFDB 2-1-1.1-c11-0-0) 9.22237939992 → Δ=4.6e-10.
S1 ✓ kampanya yüksekliklerinde taper0.3 / taper1.0 / keskin, taper0.5'e
   göre RMS farkı %0.5-2.6 (ölçüt %5).
S2 ✓ pilot pencerelerde taper motoru 0 BASAMAK, σ(d)=0.25-0.30.
   ÇARPICI: t∈[10000,10400]'de KESKİN motor 4 sıfır KAÇIRIYOR
   (n=937 vs 941), σ(d)=1.77, 93 basamak bayrağı. "Kaldırılmış dip"
   olgusu burada canlı yakalandı ve taper'ın onu kapattığı görüldü.
S3 ✓ taper-vs-keskin |Δγ|/⟨g⟩ medyan 5.8e-3 (maks büyük: keskinin
   kaçırdığı sıfırlarda eşleşme komşuya kayıyor — S2 ile tutarlı).
N1 ✓✓✓ KAPI DUYARLI: yanlış ε → 222.6, derece-1 uzunluğu → 75.1,
   yanlış normalizasyon → 777.4 basamak sürüklenmesi (doğru: 1.90).

DÜRÜST KAYIT — İLK KOŞUDA ÜÇ ŞEY YANLIŞTI, ÜÇÜ DE BURADA DURUYOR:
 (i)  G3(c) parantezi ±0.3⟨g⟩ idi; t≈146'daki YAKIN ÇİFT (γ=146.1487
      ve 146.4091, s=0.26) parantezin içine iki kesişim koydu,
      bisection komşuya kaydı ve sahte bir 2.6e-1 "hata" üretti.
      Kesin oracle o bölgede motorla 1e-3 içinde uyuşuyor
      (Z_taper/Z_keskin/Z_kesin = −0.1079/−0.1177/−0.1090 @ t=146.30).
      Parantez artık komşu sıfırların orta noktaları.
 (ii) S1'in ölçütü "1.3X ile fark < %10" idi ve RET verdi (%17-25).
      Ölçüt YANLIŞ kurulmuştu: AFE'nin uzunluğu analitik iletkenin
      kareköküyle SABİTTİR, uzatmak yaklaşımı BOZAR. (b2) bunu
      oracle'a karşı ölçtü; ölçüt aynı uzunlukta taper genişliğine
      duyarsızlığa çevrildi.
 (iii) S2'nin ölçütü "sürüklenme ARALIĞI < 1.2" idi ve RET verdi
      (1.5-1.9). Yanlış istatistik: d_i sabit değil, argüman
      dalgalanması S(t)'dir; 1000 sıfırda ~4σ ≈ 1.6-2.0 NORMALDİR.
      Ölçüt BASAMAK dedektörüne (sertifikanın gerçek aleti) çevrildi.
HÜKÜM: KAMPANYAYA GİRİLDİ (117a2) — 87 543 sıfır, T ≤ 36 000.
"""

import numpy as np
import mpmath as mp
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
NMAX_TAU = 6000            # X(t) = t/2π ⟹ t ≤ 37700 için yeterli


# ======================================================================
# 1) τ(n) — KESİN (python int)
# ======================================================================
def tau_tablosu(N=NMAX_TAU):
    """Δ = q·Π(1−q^n)^24;  E³ = Σ_{k≥0} (−1)^k (2k+1) q^{k(k+1)/2} (Jacobi).
    E^24 = (E³)^8 ⟹ seyrek çekirdekle 7 çarpım. Tam sayı aritmetiği."""
    sp = []
    k = 0
    while k * (k + 1) // 2 < N:
        sp.append((k * (k + 1) // 2, (-1) ** k * (2 * k + 1)))
        k += 1
    res = [0] * N
    for e, c in sp:
        res[e] = c
    for _ in range(7):
        new = [0] * N
        for e, c in sp:
            if c == 1:
                for i in range(N - e):
                    v = res[i]
                    if v:
                        new[i + e] += v
            elif c == -1:
                for i in range(N - e):
                    v = res[i]
                    if v:
                        new[i + e] -= v
            else:
                for i in range(N - e):
                    v = res[i]
                    if v:
                        new[i + e] += c * v
        res = new
    return [0] + res          # tau[n] = q^{n-1} katsayısı, n = 1..N


# ======================================================================
# 2) GL(2) MOTORU
# ======================================================================
def stirling_imloggamma(z):
    """Im log Γ(z), |z| büyük (vektörize Stirling)."""
    return ((z - 0.5) * np.log(z) - z + 0.5 * np.log(TWO_PI)
            + 1.0 / (12 * z) - 1.0 / (360 * z ** 3)
            + 1.0 / (1260 * z ** 5)).imag


class DeltaMotor:
    """Ramanujan Δ Hardy-Z motoru (derece 2, seviye 1, ağırlık 12).

    Z(t) = 2 Σ_{n} w_n(t)·(τ(n)/n^6)·cos(θ(t) − t log n)
    w_n : keskin (taper_c=0) veya kosinüs-taper (varsayılan c=0.5,
          genişlik W = c√X, X = t/2π) — G3'te oracle'a karşı 11 kat
          daha doğru çıktı.
    """

    def __init__(self, tau, taper_c=0.5, uzunluk=1.0):
        self.tau = tau
        self.N = len(tau) - 1
        ns = np.arange(1, self.N + 1, dtype=float)
        self.cn = np.array([tau[n] / float(n) ** 6
                            for n in range(1, self.N + 1)])
        self.logn = np.log(ns)
        self.taper_c = taper_c
        self.uzunluk = uzunluk        # X = uzunluk · t/2π  (S1 testi için)
        self.eps_faz = 0.0            # N1 negatif kontrolü için
        self.derece1 = False          # N1: X = √(t/2π)
        self.us = 6.0                 # N1: a(n)/√n = τ(n)/n^us

    # -------- θ --------
    def theta(self, t):
        t = np.asarray(t, dtype=float)
        z = 6.0 + 1j * t
        return -t * np.log(TWO_PI) + stirling_imloggamma(z) + self.eps_faz

    def theta_mp(self, t):
        s = mp.mpf(6) + 1j * mp.mpf(t)
        return float(-mp.mpf(t) * mp.log(2 * mp.pi) + mp.im(mp.loggamma(s))) \
            + self.eps_faz

    # -------- kesim uzunluğu --------
    def Xlen(self, t):
        t = np.maximum(np.asarray(t, float), 8.0)   # motor t>2π'de tanımlı
        if self.derece1:
            return np.sqrt(t / TWO_PI) * self.uzunluk
        return t / TWO_PI * self.uzunluk

    # -------- Z --------
    def Z(self, t, chunk=800):
        t = np.atleast_1d(np.asarray(t, dtype=float))
        cn = self.cn if self.us == 6.0 else np.array(
            [self.tau[n] / float(n) ** self.us for n in range(1, self.N + 1)])
        out = np.empty_like(t)
        X = self.Xlen(t)
        W = np.maximum(1.0, self.taper_c * np.sqrt(X)) if self.taper_c > 0 \
            else np.zeros_like(X)
        M = np.floor(X - W).astype(np.int64)          # tam ağırlıklı son n
        M = np.clip(M, 1, self.N)
        th = self.theta(t)
        for Mv in np.unique(M):
            m = M == Mv
            tm, thm, Xm, Wm = t[m], th[m], X[m], W[m]
            zv = np.empty(len(tm))
            Kk = 0 if self.taper_c <= 0 else int(
                min(self.N - Mv, np.ceil(2 * Wm.max()) + 2))
            for s0 in range(0, len(tm), chunk):
                sl = slice(s0, s0 + chunk)
                ph = thm[sl, None] - tm[sl, None] * self.logn[None, :Mv]
                v = 2 * (np.cos(ph) * cn[None, :Mv]).sum(axis=1)
                if Kk > 0:
                    nn = np.arange(Mv + 1, Mv + Kk + 1)
                    u = (nn[None, :] - Xm[sl, None]) / Wm[sl, None]
                    w = np.clip(0.5 * (1 - np.sin(np.pi * np.clip(u, -1, 1)
                                                  / 2)), 0.0, 1.0)
                    ph2 = (thm[sl, None]
                           - tm[sl, None] * self.logn[None, Mv:Mv + Kk])
                    v = v + 2 * (np.cos(ph2) * w
                                 * cn[None, Mv:Mv + Kk]).sum(axis=1)
                zv[sl] = v
            out[m] = zv
        return out

    # -------- sıfır bulma (grid_frac × ortalama boşluk) --------
    def ort_bosluk(self, t):
        return np.pi / np.log(np.asarray(t, float) / TWO_PI)

    def sifir_bul(self, T0, T1, grid_frac=0.03):
        T0 = max(T0, 20.0)          # X(t)=t/2π ana toplamı t≳20'de anlamlı
        ts = [T0]
        t = T0
        while t < T1:
            t += grid_frac * np.pi / np.log(t / TWO_PI)
            ts.append(t)
        ts = np.array(ts)
        v = self.Z(ts)
        sc = np.where(np.sign(v[:-1]) * np.sign(v[1:]) < 0)[0]
        a, b = ts[sc].copy(), ts[sc + 1].copy()
        fa = v[sc].copy()
        for _ in range(28):
            mid = 0.5 * (a + b)
            fm = self.Z(mid)
            sol = fa * fm <= 0
            b = np.where(sol, mid, b)
            a = np.where(sol, a, mid)
            fa = np.where(sol, fa, fm)
        return 0.5 * (a + b)

    # -------- sayım sertifikası --------
    def sayim_sapmasi(self, zz):
        th = self.theta(zz)
        return np.arange(len(zz)) - (th - th[0]) / np.pi


# ======================================================================
# 3) BAĞIMSIZ ORACLE — Hecke'nin kesin integral formülü (t ≲ 250)
# ======================================================================
def lambda_kesin(tau, t, extra=0.85, marj=25):
    """Λ(6+it) = Σ τ(n)[Γ(S,2πn)/(2πn)^S + Γ(12−S,2πn)/(2πn)^{12−S}]
    Hecke: Λ(S) = ∫_0^∞ Δ(iy) y^{S−1} dy, y=1'de bölünüp modülerlikle
    kapatılır. KESİN (hata terimi yok) ama e^{πt/2} mertebesinde iptal
    barındırır ⟹ dps ≈ 40 + 0.85·t. t ≳ 250'de kullanılamaz."""
    eski = mp.mp.dps
    mp.mp.dps = int(40 + extra * abs(t))
    S = mp.mpf(6) + 1j * mp.mpf(t)
    nmax = min(len(tau) - 1, int(abs(t) / TWO_PI) + marj)
    tot = mp.mpc(0)
    for n in range(1, nmax + 1):
        if tau[n] == 0:
            continue
        x = 2 * mp.pi * n
        tot += tau[n] * (mp.gammainc(S, x) / mp.power(x, S)
                         + mp.gammainc(12 - S, x) / mp.power(x, 12 - S))
    pref = mp.power(2 * mp.pi, -6) * abs(mp.gamma(S))
    z = float(mp.re(tot) / pref)
    imr = float(abs(mp.im(tot))) / max(float(abs(mp.re(tot))), 1e-300)
    mp.mp.dps = eski
    return z, imr


def Z_kesin_guvenli(tau, t):
    """İki farklı dps ile hesapla; uyuşmazsa None döndür (dürüst kapı)."""
    a, imr = lambda_kesin(tau, t, extra=0.80)
    b, _ = lambda_kesin(tau, t, extra=1.15)
    if abs(a - b) > 1e-6 * max(1.0, abs(a)):
        return None, imr
    return a, imr


# ======================================================================
# ======================  KAPILAR (bu dosya doğrudan çalıştırılınca) ====
# ======================================================================
if __name__ == "__main__":
    print("=" * 74, flush=True)
    print("117a — Δ MOTORU SAĞLAMA KAPILARI", flush=True)
    print("=" * 74, flush=True)
    GECTI = {}

    # ---------------- G1: τ(n) ----------------
    t0 = time.time()
    TAU = tau_tablosu(NMAX_TAU)
    print(f"\nG1 — τ(n), n ≤ {NMAX_TAU} ({time.time()-t0:.2f} sn)", flush=True)
    ref = {2: -24, 3: 252, 5: 4830, 7: -16744, 11: 534612}
    g1 = all(TAU[n] == v for n, v in ref.items())
    print("  referans değerler: " + "  ".join(
        f"τ({n})={TAU[n]}{'✓' if TAU[n]==v else '✗(bek '+str(v)+')'}"
        for n, v in ref.items()), flush=True)
    carp = [(6, 2, 3), (10, 2, 5), (15, 3, 5), (35, 5, 7), (77, 7, 11)]
    for n, p, q in carp:
        ok = TAU[n] == TAU[p] * TAU[q]
        g1 &= ok
        print(f"  çarpımsallık τ({n}) = τ({p})τ({q}): {'✓' if ok else '✗'}",
              flush=True)
    for p in [2, 3, 5, 7, 11, 13]:
        ok = TAU[p * p] == TAU[p] ** 2 - p ** 11
        g1 &= ok
        print(f"  Hecke τ({p}²) = τ({p})² − {p}^11: {'✓' if ok else '✗'}",
              flush=True)
    # Deligne
    asal = []
    sieve = np.ones(NMAX_TAU + 1, bool); sieve[:2] = False
    for i in range(2, int(NMAX_TAU ** .5) + 1):
        if sieve[i]:
            sieve[i * i::i] = False
    asal = np.nonzero(sieve)[0]
    ap = np.array([TAU[int(p)] / float(p) ** 5.5 for p in asal])
    ihlal = np.nonzero(np.abs(ap) > 2.0)[0]
    g1 &= len(ihlal) == 0
    print(f"  DELIGNE |a(p)| ≤ 2: {len(asal)} asal, ihlal {len(ihlal)}; "
          f"maks |a(p)| = {np.abs(ap).max():.4f} @ p={asal[np.abs(ap).argmax()]}"
          f"  {'✓' if len(ihlal)==0 else '✗'}", flush=True)
    GECTI["G1"] = g1

    M = DeltaMotor(TAU, taper_c=0.5)
    Mkeskin = DeltaMotor(TAU, taper_c=0.0)

    # ---------------- G2: θ ----------------
    print("\nG2 — θ: Stirling vs mpmath", flush=True)
    mp.mp.dps = 30
    dmax = 0.0
    for tt in [200.0, 1000.0, 5000.0, 20000.0, 40000.0]:
        d = abs(M.theta_mp(tt) - float(M.theta(np.array([tt]))[0]))
        dmax = max(dmax, d)
        print(f"  t={tt:8.0f}: |Δθ| = {d:.2e}", flush=True)
    GECTI["G2"] = dmax < 1e-9
    print(f"  KAPI: {'✓' if GECTI['G2'] else '✗'} (eşik 1e-9)", flush=True)

    # ---------------- G3: bağımsız oracle ----------------
    print("\nG3 — BAĞIMSIZ ORACLE (Hecke kesin integral formülü)", flush=True)
    onb = HERE / "117a_oracle_ornek.npz"
    if onb.exists():
        _d = np.load(onb)
        ok_t, ok_z, imr_max = _d["t"], _d["z"], float(_d["imr"])
        ts = ok_t
        print(f"  (oracle örneklemi önbellekten: {len(ok_t)} nokta)",
              flush=True)
    else:
        ts = np.linspace(45, 185, 21)
        ok_t, ok_z = [], []
        imr_max = 0.0
        for tt in ts:
            z, imr = Z_kesin_guvenli(TAU, float(tt))
            imr_max = max(imr_max, imr)
            if z is not None:
                ok_t.append(float(tt)); ok_z.append(z)
        ok_t = np.array(ok_t); ok_z = np.array(ok_z)
        np.savez(onb, t=ok_t, z=ok_z, imr=imr_max)
    zt = M.Z(ok_t); zk = Mkeskin.Z(ok_t)
    rms = np.sqrt((ok_z ** 2).mean())
    e_taper = np.abs(zt - ok_z); e_keskin = np.abs(zk - ok_z)
    print(f"  (a) Λ GERÇEK mi: maks |Im/Re| = {imr_max:.1e}  "
          f"{'✓ (ε=+1, FE doğrulandı)' if imr_max < 1e-25 else '✗'}",
          flush=True)
    print(f"  (b) {len(ok_t)}/{len(ts)} nokta güvenilir (çift-dps kapısı); "
          f"RMS|Z| = {rms:.3f}", flush=True)
    print(f"      KESKİN kesim : ort|Δ| = {e_keskin.mean():.4f} "
          f"({100*e_keskin.mean()/rms:.2f}%), maks {e_keskin.max():.4f}",
          flush=True)
    print(f"      TAPER c=0.5  : ort|Δ| = {e_taper.mean():.4f} "
          f"({100*e_taper.mean()/rms:.2f}%), maks {e_taper.max():.4f}",
          flush=True)
    print(f"      kazanç: {e_keskin.mean()/max(e_taper.mean(),1e-12):.1f} kat",
          flush=True)
    # (b2) KESİM UZUNLUĞU GERÇEKTEN t/2π mi? — oracle'a karşı sınanır.
    #      AFE'nin uzunluğu analitik iletkenin karekökü ile SABİTTİR;
    #      "daha uzun toplam daha iyidir" SEZGİSİ YANLIŞTIR ve bu satır
    #      onu ölçerek gösterir (S1'in ilk ölçütü bu sezgiye dayanıyordu
    #      ve haklı olarak RET verdi — ölçüt düzeltildi, bkz. S1).
    e_uz = {}
    for uz in (0.7, 1.0, 1.3):
        Mu = DeltaMotor(TAU, taper_c=0.5, uzunluk=uz)
        e_uz[uz] = np.abs(Mu.Z(ok_t) - ok_z).mean()
    print("      (b2) KESİM UZUNLUĞU X = u·t/2π, oracle'a karşı ort|Δ|: "
          + "  ".join(f"u={u}: {v:.4f}" for u, v in e_uz.items())
          + f"   ⟹ en iyi u = {min(e_uz, key=e_uz.get)}", flush=True)

    # (c) motor sıfırları vs oracle bisection (motorun geçerli olduğu bölge)
    def oracle_sifir(a, b, adim=26):
        fa, _ = lambda_kesin(TAU, a)
        for _ in range(adim):
            m0 = 0.5 * (a + b)
            fm, _ = lambda_kesin(TAU, m0)
            if fa * fm <= 0:
                b = m0
            else:
                a, fa = m0, fm
        return 0.5 * (a + b)

    zz_motor = M.sifir_bul(40.0, 150.0, grid_frac=0.02)
    idx = np.linspace(1, len(zz_motor) - 2, 5).astype(int)
    print(f"  (c) motor {len(zz_motor)} sıfır buldu t∈[40,150]; 5'i oracle "
          f"bisection'ıyla karşılaştırılıyor", flush=True)
    # NOT: parantez KOMŞU SIFIRLARIN ORTA NOKTALARI olmalı. İlk denemede
    # ±0.3⟨g⟩ kullanılmıştı ve t≈146'da YAKIN ÇİFT (γ=146.1487 ve 146.4091,
    # s=0.26) parantezin içinde iki kesişim bıraktı; bisection komşu sıfıra
    # kaydı ve sahte bir 2.6e-1 "hata" üretti. Kesin oracle o bölgede
    # motorla 1e-3 içinde uyuşuyor (Z_taper/Z_keskin/Z_kesin üçlüsü
    # −0.1079/−0.1177/−0.1090 @ t=146.30) — kusur KAPIDAYDI, motorda değil.
    farklar = []
    for i in idx:
        z = zz_motor[i]
        g = float(M.ort_bosluk(z))
        a = 0.5 * (zz_motor[i - 1] + z)
        b = 0.5 * (z + zz_motor[i + 1])
        r = oracle_sifir(float(a), float(b))
        farklar.append(abs(r - z) / g)
        print(f"      γ_motor={z:10.6f}  oracle {r:10.6f}  "
              f"|Δγ|/⟨g⟩ = {abs(r-z)/g:.2e}  (parantez {b-a:.3f})", flush=True)

    # (d) ORACLE'IN KENDİ ilk sıfırları (motordan tamamen bağımsız)
    print("  (d) oracle'ın kendi ilk sıfırları (motor devrede DEĞİL):",
          flush=True)
    gs = np.arange(6.0, 32.0, 0.25)
    vals = [lambda_kesin(TAU, float(x))[0] for x in gs]
    ilk = []
    for i in range(len(gs) - 1):
        if vals[i] * vals[i + 1] < 0:
            ilk.append(oracle_sifir(float(gs[i]), float(gs[i + 1])))
    print("      " + "  ".join(f"{g:.8f}" for g in ilk[:5]), flush=True)
    LIT = 9.22237939992
    d1 = abs(ilk[0] - LIT) if ilk else 9e9
    print(f"      γ₁(oracle) = {ilk[0]:.9f}  vs literatür {LIT} → "
          f"Δ = {d1:.2e}", flush=True)
    GECTI["G3"] = (imr_max < 1e-25 and e_taper.mean() < 0.05 * rms
                   and len(farklar) >= 5 and max(farklar) < 0.02
                   and d1 < 1e-6 and min(e_uz, key=e_uz.get) == 1.0)
    print(f"  KAPI G3: {'✓' if GECTI['G3'] else '✗'}", flush=True)

    # ---------------- S1: kesim duyarlılığı ----------------
    # ÖLÇÜT DÜZELTMESİ (dürüst kayıt): ilk koşuda S1'in ölçütü "1.3X ile
    # fark < %10" idi ve RET verdi (%17-25). Teşhis: ölçüt YANLIŞ kurulmuş.
    # AFE'nin ana toplam uzunluğu analitik iletkenin kareköküyle SABİTTİR
    # (X = t/2π); toplamı uzatmak yaklaşımı İYİLEŞTİRMEZ, BOZAR. (b2)
    # bunu oracle'a karşı doğrudan ölçüyor. Doğru S1 ölçütü: aynı
    # uzunlukta FARKLI taper genişliklerinin birbirine yakın olması.
    print("\nS1 — KESİM DUYARLILIĞI (kampanya yükseklikleri)", flush=True)
    M10 = DeltaMotor(TAU, taper_c=1.0)
    M03 = DeltaMotor(TAU, taper_c=0.3)
    M13 = DeltaMotor(TAU, taper_c=0.5, uzunluk=1.3)
    s1ok = True
    for lo, hi in [(1000., 1050.), (10000., 10050.), (19950., 20000.)]:
        tt = np.linspace(lo, hi, 400)
        za = M.Z(tt)
        r = np.sqrt((za ** 2).mean())
        f = lambda x: 100 * np.sqrt(((x - za) ** 2).mean()) / r
        print(f"  t∈[{lo:.0f},{hi:.0f}] RMS|Z|={r:.3f} | taper0.3 "
              f"{f(M03.Z(tt)):5.2f}% | taper1.0 {f(M10.Z(tt)):5.2f}% | "
              f"keskin {f(Mkeskin.Z(tt)):5.2f}%  ‖ (bilgi) 1.3X "
              f"{f(M13.Z(tt)):5.2f}%", flush=True)
        s1ok &= (f(M03.Z(tt)) < 5 and f(M10.Z(tt)) < 5
                 and f(Mkeskin.Z(tt)) < 5)
    GECTI["S1"] = s1ok

    # ---------------- S2/S3: sayım + sıfır kararlılığı ----------------
    # ÖLÇÜT DÜZELTMESİ (dürüst kayıt): ilk koşuda S2'nin ölçütü "sürüklenme
    # ARALIĞI < 1.2" idi ve RET verdi (1.5-1.9). Teşhis: ölçüt yanlış
    # istatistiğe bakıyordu. d_i = i − θ(z_i)/π farkı SABİT değil, argüman
    # dalgalanması S(t)'dir; σ_S ~ 0.4-0.6 olduğundan 1000 sıfırlık bir
    # pencerede aralığın ~4σ ≈ 1.6-2.0 çıkması NORMALDİR. Sertifikanın
    # aradığı şey aralık değil BASAMAK (kayıp/fazla sıfır) — GL(1)
    # boru hattındaki yuvarlanan-medyan dedektörünün tam olarak ölçtüğü şey.
    # Doğru ölçüt: basamak dedektörü 0 bölge bulmalı; σ(d) makul olmalı.
    print("\nS2 — SAYIM SERTİFİKASI (pilot pencereler)", flush=True)
    from numpy.lib.stride_tricks import sliding_window_view

    def basamak_say(d, win=80, esik=0.6):
        if len(d) < win + 2:
            return -1
        med = np.median(sliding_window_view(d, win), axis=1)
        adim = med[win + 1:] - med[:-(win + 1)]
        return int((np.abs(adim) > esik).sum())

    s2ok = True
    pilots = {}
    for lo, hi in [(200., 700.), (10000., 10400.), (19600., 20000.)]:
        z1 = M.sifir_bul(lo, hi, grid_frac=0.03)
        z2 = Mkeskin.sifir_bul(lo, hi, grid_frac=0.03)
        d1_ = M.sayim_sapmasi(z1); d1_ -= np.median(d1_[:40])
        d2_ = Mkeskin.sayim_sapmasi(z2); d2_ -= np.median(d2_[:40])
        pilots[(lo, hi)] = (z1, z2)
        b1, b2 = basamak_say(d1_), basamak_say(d2_)
        print(f"  t∈[{lo:.0f},{hi:.0f}] TAPER : n={len(z1)}, σ(d)={d1_.std():.2f},"
              f" aralık [{d1_.min():+.2f},{d1_.max():+.2f}], BASAMAK={b1}",
              flush=True)
        print(f"  {'':>18} KESKİN: n={len(z2)}, σ(d)={d2_.std():.2f},"
              f" aralık [{d2_.min():+.2f},{d2_.max():+.2f}], BASAMAK={b2}"
              f"   [taper−keskin = {len(z1)-len(z2):+d} sıfır]", flush=True)
        s2ok &= (b1 == 0 and d1_.std() < 1.0)
    GECTI["S2"] = s2ok

    print("\nS3 — SIFIR KARARLILIĞI (taper vs keskin eşlemesi)", flush=True)
    alld = []
    for (lo, hi), (z1, z2) in pilots.items():
        idx = np.clip(np.searchsorted(z2, z1), 1, len(z2) - 1)
        yak = np.where(np.abs(z2[idx] - z1) < np.abs(z2[idx - 1] - z1),
                       z2[idx], z2[idx - 1])
        g = np.pi / np.log(z1 / TWO_PI)
        alld.append(np.abs(yak - z1) / g)
    alld = np.concatenate(alld)
    print(f"  |Δγ|/⟨g⟩: medyan {np.median(alld):.2e}, "
          f"%95 {np.percentile(alld,95):.2e}, maks {alld.max():.2e}",
          flush=True)
    print(f"  (maks büyükse sebep: keskin motorun KAÇIRDIĞI sıfırlar — "
          f"eşleşme komşuya kayar; bkz. S2 taper−keskin sütunu)", flush=True)
    GECTI["S3"] = np.median(alld) < 0.05

    # ---------------- N1: negatif kontroller ----------------
    print("\nN1 — NEGATİF KONTROLLER (kapı duyarlı mı?)", flush=True)
    lo, hi = 10000., 10400.
    z_ref = pilots[(lo, hi)][0]
    d = M.sayim_sapmasi(z_ref); d -= np.median(d[:40])
    print(f"  DOĞRU konfigürasyon: sürüklenme aralığı "
          f"{d.max()-d.min():.2f}", flush=True)
    n1ok = True
    for isim, kur in [("yanlış ε (θ+π/2)", dict(eps_faz=np.pi / 2)),
                      ("derece-1 uzunluk √(t/2π)", dict(derece1=True)),
                      ("yanlış normalizasyon n^-5.5", dict(us=5.5))]:
        Mx = DeltaMotor(TAU, taper_c=0.5)
        for k, v in kur.items():
            setattr(Mx, k, v)
        zx = Mx.sifir_bul(lo, hi, grid_frac=0.03)
        dx = M.sayim_sapmasi(zx); dx -= np.median(dx[:40])
        rng_ = dx.max() - dx.min()
        print(f"  {isim:>28}: n={len(zx):5d}  sürüklenme aralığı {rng_:8.2f}"
              f"  {'✓ (patladı)' if rng_ > 5 else '✗ (KAPI KÖR!)'}",
              flush=True)
        n1ok &= rng_ > 5
    GECTI["N1"] = n1ok

    print("\n" + "=" * 74, flush=True)
    print("KAPI ÖZETİ: " + "  ".join(
        f"{k}={'✓' if v else '✗'}" for k, v in GECTI.items()), flush=True)
    print("HÜKÜM: " + ("TÜM KAPILAR GEÇİLDİ — kampanyaya girilebilir (117a2)"
                       if all(GECTI.values()) else
                       "KAPI KALDI — KAMPANYAYA GİRİLMEZ"), flush=True)
    print("=" * 74, flush=True)
