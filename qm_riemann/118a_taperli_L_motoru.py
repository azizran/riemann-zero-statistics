"""
118a — TAPER'LI GL(1) MOTORU + YER-GERÇEĞİ SINAVI (26 Ağustos 2026)
==========================================================================
ÖN-MÜHÜR (ölçümden ÖNCE yazıldı)
==========================================================================
117 (Ramanujan Δ, GL(2)) seferinin sürprizi: "KALDIRILMIŞ DİP" olgusunun
asıl nedeni RS-düzeltmesinin eksikliği DEĞİL, ana toplamın KESKİN
KESİLMESİYMİŞ. Kosinüs-taper'lı kesim (W = 0.5√X) oracle'a karşı hatayı
11-14 KAT düşürdü ve t∈[10000,10400]'de keskin motorun kaçırdığı 4 sıfırı
geri getirdi. 117'nin kapanış cümlesi: "AÇIK İŞ: 101/105'in yedi GL(1)
adasının taper'lı yeniden dökümü."

Bu script o işin BİRİNCİ adımıdır: taper'ı GL(1) motoruna (98_L_motoru)
taşır ve motoru YER GERÇEĞİNE karşı sınar.

TAPER (117a kalıbı birebir; tek değişen X'in tanımı):
  GL(1) AFE ana toplamı  X(t) = √(qt/2π)   (GL(2)'de X = t/2π idi)
  W(t) = max(1, c·√X),  c = 0.5
  n ≤ X−W          : tam ağırlık
  X−W < n < X+W    : w(n) = ½[1 − sin(π·u/2)],  u = (n−X)/W
  n ≥ X+W          : 0
  (keskin motor: w = 1 for n ≤ ⌊X⌋, 0 sonra — 98'in aynısı)

KAPILAR
  G1  θ VE ε DEĞİŞMEMELİ. Taper yalnız toplamın ağırlıklarına dokunur;
      θ(t), arg ε, |ε| BİT-BİT aynı çıkmalı (98-S1 mpmath karşılaştırması
      da tekrarlanır). Ayrıca c=0 kurulumu 98'in Z'siyle BİREBİR aynı
      olmalı (kalıtımın doğru kurulduğunun kanıtı).
  G2  BAĞIMSIZ ORACLE — GL(1)'DE 117'DEN DAHA GÜÇLÜ. Dirichlet L için
      mpmath Hurwitz-ζ kesin Z verir (117'de oracle t≲250'de ölüyordu;
      burada TÜM yükseklikte çalışır). Yedi adada da:
        (a) Z gerçek mi (|Im/Re|) — FE + ε konvansiyonu denetimi
        (b) keskin vs taper ort|Δ| / RMS|Z|
        (c) c ∈ {0, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5} taraması: 117'nin
            c = 0.5'i GL(1)'de de iyi mi?
  G3  SIFIR KONUMU DOĞRULUĞU — MPMATH YER GERÇEĞİNE KARŞI. Sertifikalı
      önbellek sıfırlarından örneklem alınır; her biri KESİN Z ile
      (3 mpmath eval + parabol) GERÇEK köke rafine edilir; sonra
      keskin ve taper motorlarının kendi sıfırlarının o gerçek köke
      uzaklığı |Δγ|/⟨g⟩ olarak kıyaslanır.
      [ÖLÇÜT DÜZELTMESİ — dürüst kayıt: ilk koşuda G3 "önbellek
       sıfırlarında medyan |Z|" idi. O ölçüt DAİRESELDİR: önbellek
       sıfırları zaten KESKİN motorun kökleridir, dolayısıyla keskin
       motor orada tanım gereği 1e-11 verir (ölçüldü: 1.7-2.1e-11),
       taper 1.5e-2 — ve "keskin daha iyi" gibi görünür. Ölçülen şey
       motorun doğruluğu değil, listenin hangi motordan geldiğiydi.
       Ölçüt BAĞIMSIZ yer gerçeğine (mpmath kökü) çevrildi.]
  G4  ★ ASIL SINAV — YER GERÇEĞİ ★
      101e/101f (ve 105b-C/D) mpmath dip-kurtarmasıyla, işaret-taraması
      GÖREMEDİĞİ için elle kurtarılmış sıfırlar var:
        eski dörtlü : R = 101f_zeros \\ 101c_zeros
        yeni üçlü   : R = 105b_zeros \\ 105b_ham
      Bu R kümesi, KESKİN motorun işaret-taramasıyla bulamadığı
      sıfırların TAM LİSTESİDİR (tanım gereği keskin taramada %0).
      SORU: taper'lı motor SALT İŞARET-TARAMASIYLA (grid_frac 0.03,
      mpmath YOK) bunların yüzde kaçını buluyor?
      Eşleme: karşılıklı-tekil (injektif) en-yakın atama, tolerans
      0.25⟨g⟩ (kurtarılanlar dar ÇİFTLER hâlinde geldiği için tekillik
      şart: yoksa çiftin iki üyesi tek sıfıra eşleşir ve sayı şişer).

BEKLENTİ (mühür): G1-G3 geçer. G4'te BÜYÜK ÇOĞUNLUK (>%80) beklenir —
117'nin t∈[10000,10400] pilotu keskinin kaçırdığı 4 sıfırın 4'ünü de
taper'la bulmuştu. Bulunamayanlar DÜRÜSTÇE listelenir; bulunamayan bir
sıfır iki şeyden biridir: (i) taper'ın da kapatamadığı gerçek dip,
(ii) 101f'nin budayamadığı SAHTE kurtarma (o zaman taper'ın bulmaması
bir KUSUR değil bir KAZANÇtır) — ikisi ayırt edilemez, ikisi de yazılır.

YAN ÜRÜN: taper'lı tam-menzil taramalar 118b_{ada}_ham.npz olarak
önbelleklenir; 118b boru hattı onları A adımı olarak kullanır.
==========================================================================
SONUÇ (26 Ağustos — DÖRT KAPI DA GEÇİLDİ; bir ölçüt düzeltmesi kayıtlı)
==========================================================================
G1 ✓ YEDİ ADADA DA |Δθ| = 0, |Δε| = 0, |ε| = 1.000000000, ve c = 0
   kurulumu 98'in Z'siyle BİT-BİT aynı (|ΔZ| = 0). mpmath θ farkı 0.
   Taper yalnız ağırlıklara dokunuyor; konvansiyon değişmedi.

G2 ✓✓ BAĞIMSIZ ORACLE (mpmath Hurwitz, kampanya yüksekliklerinde):
   Λ gerçek: maks |Im/Re| = 2.0e-11 … 5.4e-11 (FE + ε ✓, yedi adada).
     ada    RMS|Z|   keskin   c=0.3   c=0.5   c=0.7   c=1.0   kazanç(0.5)
     chi3    1.798    3.38%   0.50%   0.46%   0.29%   0.34%     7.4×
     beta    1.958    3.26%   0.79%   0.29%   0.17%   0.30%    11.1×
     chi5    4.972    1.34%   0.24%   0.14%   0.13%   0.10%     9.7×
     chi7    3.061    2.99%   0.64%   0.34%   0.26%   0.25%     8.7×
     chi5e   3.210    1.92%   0.52%   0.23%   0.15%   0.17%     8.4×
     chi8e   1.835    2.22%   0.56%   0.45%   0.27%   0.24%     4.9×
     chi8o   2.302    4.34%   0.70%   0.35%   0.21%   0.21%    12.5×
   117'nin GL(2)'de ölçtüğü 11-14 kat, GL(1)'de 4.9-12.5 kat olarak
   TEKRARLANDI. NOT: GL(1)'de optimum c 0.7-1.0'a kayıyor (117'de 0.5
   idi); c = 0.5 yine de bir mertebe kazandırıyor ve 117 KALIBINA
   SADIK kalmak için (yedi adanın tek bir konvansiyonla dökülmesi)
   korundu. c = 0.7'ye geçmek ~1.5× daha kazandırırdı — AÇIK İŞ.

G3 ✓✓ KONUM DOĞRULUĞU (mpmath köküne karşı, ada başına 30 örneklem):
   medyan |Δγ|/⟨g⟩  keskin 3.5e-3 … 1.6e-2  →  taper 2.1e-4 … 6.6e-4
   kazanç 11-50 kat (chi5e 11.1, chi5 14.5, chi8e 14.8, beta 22.3,
   chi3 26.6, chi7 32.4, chi8o 50.0). 105a-G5'in "~1e-2 @t≈1e3"
   rakamı taper'la ~1e-4'e iniyor.

G4 ✓✓✓ YER GERÇEĞİ — SEFERİN ASIL SINAVI:
     ada    kurtarılmış  taper buldu    %     keskin    %
     chi3        288         284      98.6      0     0.0
     beta        320         312      97.5      0     0.0
     chi5        322         316      98.1      0     0.0
     chi7        619         610      98.5      0     0.0
     chi5e       378         374      98.9      0     0.0
     chi8e       132         131      99.2      1     0.8
     chi8o       440         440     100.0      0     0.0
     TOPLAM     2499        2467      98.7      1     0.04
   101/105'te mpmath dip-kurtarmasıyla ELLE çıkarılmış 2499 sıfırın
   %98.7'si taper'lı motorda SALT İŞARET-TARAMASIYLA (mpmath YOK)
   geri geliyor. Keskin motor aynı kümede 1/2499 buluyor.

BULUNAMAYAN 32 SIFIRIN TEŞHİSİ (dürüst kayıt; hepsi KESİN Z ile VE
taper-Z'nin kendi 0.001⟨g⟩ örneklemesiyle tek tek incelendi — ve ilk
hipotezim ÇÜRÜDÜ, aşağıda):
  Hepsi dar ÇİFTLER hâlinde, 18 küme; ayrım 0.014-0.118⟨g⟩. Kesin Z
  her kümede 2-3 işaret değişimi gösteriyor: sıfırlar GERÇEK.
  ÖN HİPOTEZ (ÇÜRÜDÜ): "taper dibi indirdi, sorun sadece 0.03⟨g⟩'lik
  IZGARA ADIMI; grid_frac küçültülünce kapanır."
  ÖLÇÜM: taper-Z'nin KENDİSİ 18 kümenin yalnız 5'inde işaret
  değiştiriyor; 13'ünde DİP HÂLÂ KALDIRILMIŞ (min|taper-Z| =
  2e-4 … 1.5e-2, sıfırın üstünde). Ve ızgara sıklaştırma testi
  (grid_frac 0.03 → 0.01, üç adada) yalnız χ₃'te +2 sıfır getirdi,
  β ve χ₇'de HİÇ (312/320 ve 610/619 aynı kaldı).
  ⟹ DOĞRU HÜKÜM: TAPER, KALDIRILMIŞ DİP MEKANİZMASINI YOK ETMİYOR;
  NÜFUSUNU ~60 KAT AZALTIYOR (2499 → ~40). En dar çiftlerde taper'ın
  kalan ~3e-4⟨g⟩ hatası dibi hâlâ sıfırın üstünde tutuyor. Bu, işi
  bitiren şeyin BORU HATTI olduğunu gösteriyor: 118b'nin C adımı
  (mpmath dip-kurtarma) ada başına 2-8 sıfırla bu kalıntıyı kapatıyor
  ve kalan şüpheli bölge 0'a iniyor. Taper'ın kazancı, dip-kurtarmayı
  GEREKSİZ kılmak değil, onu SAYIM SERTİFİKASININ YAKALAYABİLECEĞİ
  kadar seyrekleştirmektir (117'nin Δ'da 0 çıkması, GL(2)'nin daha
  büyük boşluklarının şansıymış).
  AYRICA 5'i (χ₇: 935.4419, 1531.0633, 1913.7452, 3081.6688,
  3081.8891) 107c'nin KOPYA ENVANTERİNDEKİ (u < 0.005) çiftler —
  eski dökümün SAHTE kopyaları. Taper'ın onları üretmemesi bir
  kaçırma değil, bir DÜZELTMEdir.

DÜRÜST KAYIT — ÖLÇÜT DÜZELTMESİ (G3): ilk koşuda G3 "sertifikalı
önbellek sıfırlarında medyan |Z|" idi ve RET verdi (taper 1.5e-2 vs
keskin 1.7e-11). Teşhis: ölçüt DAİRESEL. O sıfırlar zaten KESKİN
motorun kökleri; keskin motor orada tanım gereği makine sıfırı verir.
Ölçülen şey motorun doğruluğu değil, listenin hangi motordan geldiğiydi.
Ölçüt bağımsız yer gerçeğine (mpmath kökü) çevrildi ve 11-50 kat
kazanç ortaya çıktı. (117'nin S1/S2 ölçüt düzeltmeleriyle aynı sınıf
hata: kapı yanlış istatistiğe bakıyordu.)

HÜKÜM: TAPER'LI GL(1) MOTORU KAMPANYAYA HAZIR (118b).
"""
import numpy as np
import mpmath as mp
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
exec(open(HERE / "98_L_motoru.py").read().split('CHI4 = ')[0])
HERE = Path(__file__).resolve().parent
mp.mp.dps = 15

CHI3 = {0: 0, 1: 1, 2: -1}
CHI4 = {0: 0, 1: 1, 2: 0, 3: -1}
CHI5 = {0: 0, 1: 1, 2: 1j, 3: -1j, 4: -1}
z6c = np.exp(1j * np.pi / 3)
CHI7 = {0: 0, 1: 1, 3: z6c, 2: z6c**2, 6: z6c**3, 4: z6c**4, 5: z6c**5}
CHI5E = {0: 0, 1: 1, 2: -1, 3: -1, 4: 1}
CHI8E = {0: 0, 1: 1, 2: 0, 3: -1, 4: 0, 5: -1, 6: 0, 7: 1}
CHI8O = {0: 0, 1: 1, 2: 0, 3: 1, 4: 0, 5: -1, 6: 0, 7: -1}

# (etiket, q, tablo, a, T1, eski_ham_dosya, eski_son_dosya)
ADALAR = [
    ("chi3",  3, CHI3,  1, 55000.0, "101c_chi3_zeros.npz",  "101f_chi3_zeros.npz"),
    ("beta",  4, CHI4,  1, 50000.0, "101c_beta_zeros.npz",  "101f_beta_zeros.npz"),
    ("chi5",  5, CHI5,  1, 48000.0, "101c_chi5_zeros.npz",  "101f_chi5_zeros.npz"),
    ("chi7",  7, CHI7,  1, 42000.0, "101c_chi7_zeros.npz",  "101f_chi7_zeros.npz"),
    ("chi5e", 5, CHI5E, 0, 48000.0, "105b_chi5e_ham.npz",   "105b_chi5e_zeros.npz"),
    ("chi8e", 8, CHI8E, 0, 42000.0, "105b_chi8e_ham.npz",   "105b_chi8e_zeros.npz"),
    ("chi8o", 8, CHI8O, 1, 42000.0, "105b_chi8o_ham.npz",   "105b_chi8o_zeros.npz"),
]


# ======================================================================
# TAPER'LI GL(1) MOTORU
# ======================================================================
class LmotorT(Lmotor):
    """98'in Lmotor'u + kosinüs-taper'lı AFE kesimi (117a kalıbı).

    taper_c = 0  →  98'in KESKİN motoruyla BİT-BİT aynı (G1'de sınanır).
    taper_c > 0  →  W = taper_c·√X genişliğinde kosinüs penceresi.
    θ, ε, karakter fazları DEĞİŞMEZ (kalıtım).
    """

    def __init__(self, q, chi_table, a_par, taper_c=0.5):
        Lmotor.__init__(self, q, chi_table, a_par)
        self.N = len(self.logn)
        self.taper_c = taper_c
        self.amp = self.abschi * self.wn

    def Xlen(self, t):
        """GL(1) AFE uzunluğu: analitik iletken qt/2π, karekökü."""
        return np.sqrt(self.q * np.asarray(t, float) / TWO_PI)

    def Z(self, t, chunk=4000):
        t = np.atleast_1d(np.asarray(t, dtype=float))
        out = np.empty_like(t)
        X = self.Xlen(t)
        if self.taper_c <= 0:
            W = np.zeros_like(X)
            M = X.astype(np.int64)            # 98: N = ⌊√(qt/2π)⌋
        else:
            W = np.maximum(1.0, self.taper_c * np.sqrt(X))
            M = np.floor(X - W).astype(np.int64)
        M = np.clip(M, 1, self.N)
        th = self.theta(t)
        amp = self.amp
        for Mv in np.unique(M):
            m = M == Mv
            tm, thm, Xm, Wm = t[m], th[m], X[m], W[m]
            zv = np.empty(len(tm))
            Kk = 0 if self.taper_c <= 0 else int(
                min(self.N - Mv, np.ceil(2 * Wm.max()) + 2))
            for s0 in range(0, len(tm), chunk):
                sl = slice(s0, s0 + chunk)
                ph = (thm[sl, None] - tm[sl, None] * self.logn[None, :Mv]
                      + self.phichi[None, :Mv])
                v = 2 * (np.cos(ph) * amp[None, :Mv]).sum(axis=1)
                if Kk > 0:
                    nn = np.arange(Mv + 1, Mv + Kk + 1, dtype=float)
                    u = (nn[None, :] - Xm[sl, None]) / Wm[sl, None]
                    w = np.clip(0.5 * (1 - np.sin(np.pi * np.clip(u, -1, 1)
                                                  / 2)), 0.0, 1.0)
                    ph2 = (thm[sl, None]
                           - tm[sl, None] * self.logn[None, Mv:Mv + Kk]
                           + self.phichi[None, Mv:Mv + Kk])
                    v = v + 2 * (np.cos(ph2) * w
                                 * amp[None, Mv:Mv + Kk]).sum(axis=1)
                zv[sl] = v
            out[m] = zv
        return out

    def ort_bosluk(self, t):
        return TWO_PI / np.log(self.q * np.asarray(t, float) / TWO_PI)

    def sayim_sapmasi(self, zz):
        th = self.theta(zz)
        return np.arange(len(zz)) - (th - th[0]) / np.pi


def kesin_Z(M, tab, t):
    """BAĞIMSIZ ORACLE: L(s,χ) = q^{-s} Σ_a χ(a) ζ(s, a/q) (Hurwitz).
    Motorun AFE formülünden bağımsız; Z = e^{iθ}L GERÇEK olmalı."""
    s = mp.mpc(0.5, t)
    Lv = mp.mpc(0)
    for a in range(1, M.q):
        if tab[a % M.q] != 0:
            Lv += tab[a % M.q] * mp.zeta(s, mp.mpf(a) / M.q)
    Lv *= mp.power(M.q, -s)
    th = float(M.theta(np.array([t]))[0])
    z = mp.e ** (1j * th) * Lv
    re, im = float(z.real), float(z.imag)
    return re, abs(im) / max(abs(re), 1e-300)


def injektif_esle(R, Z, tol):
    """R'nin her elemanına Z'den TEKİL en yakın eşi ata (greedy, artan
    mesafe). DÖNÜŞ: bulundu maskesi + mesafeler.

    Tekillik ŞART: kurtarılan sıfırlar dar ÇİFTLER hâlinde gelir; tekil
    olmayan eşlemede çiftin iki üyesi de tek bir motor sıfırına eşleşip
    başarı oranını YAPAY OLARAK ŞİŞİRİR."""
    if len(Z) == 0:
        return np.zeros(len(R), bool), np.full(len(R), np.inf)
    aday = []
    idx = np.searchsorted(Z, R)
    for i, r in enumerate(R):
        for j in range(max(0, idx[i] - 2), min(len(Z), idx[i] + 3)):
            aday.append((abs(Z[j] - r), i, j))
    aday.sort()
    rk = np.zeros(len(R), bool)
    zk = set()
    dist = np.full(len(R), np.inf)
    for d, i, j in aday:
        if rk[i] or j in zk or d > tol[i]:
            continue
        rk[i] = True
        zk.add(j)
        dist[i] = d
    return rk, dist


if __name__ == "__main__":
    T00 = time.time()
    print("=" * 78, flush=True)
    print("118a — TAPER'LI GL(1) MOTORU: SAĞLAMA KAPILARI", flush=True)
    print("=" * 78, flush=True)
    GECTI = {}
    rng = np.random.default_rng(118)

    # ---------------- G1: θ / ε / keskin-özdeşlik ----------------
    print("\nG1 — θ, ε DEĞİŞMEDİ Mİ? (+ c=0 ⇔ 98 özdeşliği)", flush=True)
    g1 = True
    for et, q, tab, a, T1, _, _ in ADALAR:
        M0 = Lmotor(q, tab, a)
        Mt = LmotorT(q, tab, a, taper_c=0.5)
        Mk = LmotorT(q, tab, a, taper_c=0.0)
        ts = np.array([500.3, 2000.7, 8000.11, 0.8 * T1 + 0.37])
        dth = np.abs(Mt.theta(ts) - M0.theta(ts)).max()
        deps = abs(Mt.arg_eps - M0.arg_eps) + abs(Mt.abs_eps - M0.abs_eps)
        dZ = np.abs(Mk.Z(ts) - M0.Z(ts)).max()
        # mpmath θ (98-S1)
        t0 = float(ts[1])
        s = mp.mpc(0.5, t0)
        ps_mp = float((t0 / 2) * mp.log(q / mp.pi)
                      + mp.im(mp.loggamma((s + a) / 2))) - M0.arg_eps / 2
        dmp = abs(ps_mp - float(Mt.theta(np.array([t0]))[0]))
        ok = (dth == 0.0 and deps == 0.0 and dZ == 0.0 and dmp < 1e-9
              and abs(M0.abs_eps - 1.0) < 1e-9)
        g1 &= ok
        print(f"  {et:>6} q={q} a={a}: |Δθ|={dth:.1e}  |Δε|={deps:.1e}  "
              f"|ΔZ(c=0 vs 98)|={dZ:.1e}  |ε|={M0.abs_eps:.9f}  "
              f"|Δθ_mpmath|={dmp:.1e}  {'✓' if ok else '✗'}", flush=True)
    GECTI["G1"] = g1

    # ---------------- G2: bağımsız oracle ----------------
    print("\nG2 — BAĞIMSIZ ORACLE (mpmath Hurwitz; her yükseklikte geçerli)",
          flush=True)
    CS = [0.0, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5]
    g2 = True
    ORACLE = {}
    for et, q, tab, a, T1, _, _ in ADALAR:
        Mt = LmotorT(q, tab, a, 0.5)
        tt = np.exp(np.linspace(np.log(500.0), np.log(T1), 24))
        tt = tt + rng.uniform(0.0, 1.0, len(tt))     # sıfırlara denk gelmesin
        t0 = time.time()
        ex, imr = [], 0.0
        for x in tt:
            z, r = kesin_Z(Mt, tab, float(x))
            ex.append(z)
            imr = max(imr, r)
        ex = np.array(ex)
        rms = float(np.sqrt((ex ** 2).mean()))
        err = {}
        for c in CS:
            Mx = LmotorT(q, tab, a, c)
            err[c] = float(np.abs(Mx.Z(tt) - ex).mean())
        en_iyi = min(err, key=err.get)
        kaz = err[0.0] / max(err[0.5], 1e-300)
        ORACLE[et] = (rms, err, imr, en_iyi, kaz)
        ok = (imr < 1e-9 and err[0.5] < 0.02 * rms and kaz > 3.0)
        g2 &= ok
        print(f"  {et:>6}: |Im/Re|maks={imr:.1e}  RMS|Z|={rms:.3f}  "
              f"({time.time()-t0:.0f} s)", flush=True)
        print("          ort|Δ|/RMS  " + "  ".join(
            f"c={c}:{100*err[c]/rms:5.2f}%" for c in CS)
            + f"   en iyi c={en_iyi}", flush=True)
        print(f"          KAZANÇ keskin/taper0.5 = {kaz:.1f} kat  "
              f"{'✓' if ok else '✗'}", flush=True)
    GECTI["G2"] = g2

    # ---------------- G3/G4: önbellek + yer gerçeği ----------------
    print("\nG3/G4 — KONUM DOĞRULUĞU (mpmath yer gerçeği) ve KURTARILMIŞ "
          "SIFIRLARIN YENİDEN BULUNMASI", flush=True)
    RAPOR = {}
    for et, q, tab, a, T1, ham_f, son_f in ADALAR:
        t0 = time.time()
        Mt = LmotorT(q, tab, a, 0.5)
        Mk = LmotorT(q, tab, a, 0.0)
        z_ham = np.load(HERE / ham_f)["zeros"]      # keskin işaret-taraması
        z_son = np.load(HERE / son_f)["zeros"]      # mpmath ile kurtarılmış
        # --- G4: yer gerçeği kümesi R = z_son \ z_ham ---
        j = np.searchsorted(z_ham, z_son)
        d1 = np.abs(z_son - z_ham[np.clip(j, 0, len(z_ham) - 1)])
        d2 = np.abs(z_son - z_ham[np.clip(j - 1, 0, len(z_ham) - 1)])
        R = z_son[np.minimum(d1, d2) > 1e-6]
        # --- taper'lı TAM MENZİL işaret taraması (mpmath YOK) ---
        onb = HERE / f"118b_{et}_ham.npz"
        if onb.exists():
            z_new = np.load(onb)["zeros"]
            tar = -1.0
        else:
            tt0 = time.time()
            z_new = Mt.sifir_bul(200.0, T1, grid_frac=0.03)
            tar = time.time() - tt0
            np.savez(onb, zeros=z_new, q=q, a=a, T1=T1, taper_c=0.5)
        gbar = Mt.ort_bosluk(R) if len(R) else np.array([])
        bulundu, dist = injektif_esle(R, z_new, 0.25 * gbar)
        # keskin taramanın aynı kümedeki başarısı (tanım gereği ~0)
        bul_k, _ = injektif_esle(R, z_ham, 0.25 * gbar)
        print(f"\n  [{et}] keskin-ham {len(z_ham)} | sertifikalı {len(z_son)} "
              f"| TAPER ham {len(z_new)} ({len(z_new)-len(z_ham):+d} vs keskin"
              + (f", {tar:.0f} s tarama)" if tar >= 0 else ", önbellek)"),
              flush=True)

        # --- G3: konum doğruluğu, mpmath yer gerçeğine karşı ---
        # Örneklem KESKİN taramadan alınır (o listenin her üyesi keskin
        # motorun köküdür); gerçek kök 3 kesin-Z + parabol ile bulunur.
        NS = 30
        orn = z_ham[np.linspace(0.02 * len(z_ham), 0.98 * len(z_ham),
                                NS).astype(int)]
        ek, tk = [], []
        for z0 in orn:
            g = float(Mt.ort_bosluk(z0))
            h = 0.06 * g
            y = [kesin_Z(Mt, tab, float(z0 + dx))[0] for dx in (-h, 0.0, h)]
            c2 = (y[0] - 2 * y[1] + y[2]) / (2 * h * h)
            c1 = (y[2] - y[0]) / (2 * h)
            c0 = y[1]
            if abs(c2) < 1e-14:
                continue
            disc = c1 * c1 - 4 * c2 * c0
            if disc <= 0:
                continue
            rr = [z0 + (-c1 - np.sqrt(disc)) / (2 * c2),
                  z0 + (-c1 + np.sqrt(disc)) / (2 * c2)]
            tstar = min(rr, key=lambda x: abs(x - z0))
            if abs(tstar - z0) > 0.3 * g:
                continue                      # rafine güvenilmez, atla
            j = np.searchsorted(z_new, tstar)
            en = min(abs(z_new[k] - tstar)
                     for k in range(max(0, j - 1), min(len(z_new), j + 2)))
            ek.append(abs(z0 - tstar) / g)
            tk.append(en / g)
        ek, tk = np.array(ek), np.array(tk)
        oran = dict(keskin=(float(np.median(ek)), float(np.percentile(ek, 90))),
                    taper=(float(np.median(tk)), float(np.percentile(tk, 90))),
                    n=len(ek))
        RAPOR[et] = dict(nham=len(z_ham), nson=len(z_son), nnew=len(z_new),
                         R=R, bulundu=bulundu, dist=dist,
                         nk=int(bul_k.sum()), oran=oran, tar=tar)
        print(f"        G3 |Δγ|/⟨g⟩ (mpmath köküne, {len(ek)} örneklem): "
              f"KESKİN medyan {oran['keskin'][0]:.2e} (%90 "
              f"{oran['keskin'][1]:.2e})  |  TAPER medyan "
              f"{oran['taper'][0]:.2e} (%90 {oran['taper'][1]:.2e})  "
              f"⟹ {oran['keskin'][0]/max(oran['taper'][0],1e-99):.1f} kat",
              flush=True)
        print(f"        G4 YER GERÇEĞİ: kurtarılmış {len(R)} sıfır → "
              f"taper işaret-taraması {int(bulundu.sum())} buldu "
              f"(%{100*bulundu.mean() if len(R) else 0:.1f});  "
              f"keskin tarama {int(bul_k.sum())} (%"
              f"{100*bul_k.mean() if len(R) else 0:.1f})", flush=True)
        kayip = R[~bulundu]
        if len(kayip):
            print(f"        BULUNAMAYAN {len(kayip)} sıfır (dürüst liste"
                  + (", ilk 30):" if len(kayip) > 30 else "):"), flush=True)
            print("          " + "  ".join(f"{x:.4f}" for x in kayip[:30]),
                  flush=True)
        print(f"        ({time.time()-t0:.0f} s)", flush=True)

    g3 = all(RAPOR[e]["oran"]["taper"][0] < RAPOR[e]["oran"]["keskin"][0]
             for e in RAPOR)
    tot_R = sum(len(RAPOR[e]["R"]) for e in RAPOR)
    tot_B = sum(int(RAPOR[e]["bulundu"].sum()) for e in RAPOR)
    GECTI["G3"] = g3
    GECTI["G4"] = tot_B / max(tot_R, 1) > 0.80

    print("\n" + "=" * 78, flush=True)
    print("YER GERÇEĞİ ÖZET TABLOSU", flush=True)
    print("=" * 78, flush=True)
    print(f"{'ada':>7} {'kurtarılmış':>11} {'taper buldu':>11} {'%':>6} "
          f"{'keskin':>7} {'%':>6} {'Δham':>7} {'oracle kaz.':>11} "
          f"{'|Δγ|/⟨g⟩ keskin→taper':>24}", flush=True)
    for et, q, tab, a, T1, _, _ in ADALAR:
        r = RAPOR[et]
        nR = len(r["R"])
        nB = int(r["bulundu"].sum())
        o = r["oran"]
        print(f"{et:>7} {nR:>11} {nB:>11} {100*nB/max(nR,1):>6.1f} "
              f"{r['nk']:>7} {100*r['nk']/max(nR,1):>6.1f} "
              f"{r['nnew']-r['nham']:>+7d} {ORACLE[et][4]:>10.1f}x "
              f"{o['keskin'][0]:>10.1e} → {o['taper'][0]:>8.1e}", flush=True)
    print(f"{'TOPLAM':>7} {tot_R:>11} {tot_B:>11} "
          f"{100*tot_B/max(tot_R,1):>6.1f}", flush=True)

    np.savez(HERE / "118a_yer_gercegi.npz",
             adlar=np.array([e[0] for e in ADALAR]),
             cs=np.array(CS),
             **{f"err_{e}": np.array([ORACLE[e][1][c] for c in CS])
                for e in ORACLE},
             **{f"rms_{e}": ORACLE[e][0] for e in ORACLE},
             **{f"kaz_{e}": ORACLE[e][4] for e in ORACLE},
             **{f"R_{e}": RAPOR[e]["R"] for e in RAPOR},
             **{f"bul_{e}": RAPOR[e]["bulundu"] for e in RAPOR},
             **{f"nham_{e}": RAPOR[e]["nham"] for e in RAPOR},
             **{f"nnew_{e}": RAPOR[e]["nnew"] for e in RAPOR},
             **{f"nson_{e}": RAPOR[e]["nson"] for e in RAPOR},
             **{f"dgk_{e}": RAPOR[e]["oran"]["keskin"][0] for e in RAPOR},
             **{f"dgt_{e}": RAPOR[e]["oran"]["taper"][0] for e in RAPOR})

    print("\nKAPI ÖZETİ: " + "  ".join(f"{k}={'✓' if v else '✗'}"
                                      for k, v in GECTI.items()), flush=True)
    print(f"Süre {time.time()-T00:.0f} s. 118a_yer_gercegi.npz + "
          f"118b_*_ham.npz yazıldı.", flush=True)
