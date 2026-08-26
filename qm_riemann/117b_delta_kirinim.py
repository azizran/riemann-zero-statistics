"""
117b — Δ'NIN KIRINIM SPEKTROSKOPİSİ: HECKE İŞARETLERİ VE SATAKE ÇİZGİLERİ
        (26 Ağustos 2026 — GL(2)'nin ilk kırınım deseni)
==========================================================================
ÖN-MÜHÜRLER — ÖLÇÜMDEN ÖNCE YAZILDI
==========================================================================
UNFOLD KONVANSİYONU (belgelenir): derece-2 yoğunluk θ'(t)/π =
(1/π)·log(t/2π) + O(1/t) ⟹ ortalama boşluk = π/log(t/2π) = 2π/L_eff ile
        L_eff := 2·log(t/2π)
(GL(1)'de L = log(qt/2π) idi; "L" her zaman 2π×yoğunluk.)
Açılmış sayım x(t) = θ(t)/π (kesin; asimptotik form değil).
Ĝ(ω) = ⟨e^{iωz}⟩ hem SIFIRLARIN kendisinde (orta-noktasız) hem de
ORTA NOKTALARDA (z_i+z_{i+1})/2 ölçülür.

--------------------------------------------------------------------------
P1 — ASAL BENEKLERİ + HECKE İŞARETİNİN OKUNMASI  (seferin taç sorusu)
--------------------------------------------------------------------------
MUTLAK YASA'nın GL(2) hali (89'un türetimi, açık formül katsayısı
b(p^k) = α_p^k + β_p^k ile; GL(1)'de b ≡ χ(q)):
     |Ĝ(log q)| = |b(q)|·Λ(q)/(L_eff·√q)·cos(π τ_q)·DW,   τ_q = log q/L_eff
     FAZ(q)     = 180° + arg b(q)   ⟹  reel b: b>0 → 180°,  b<0 → 0°
Δ için b(p) = a_p = τ(p)/p^{11/2}. ζ'da bütün fazlar 180°'ydi; burada
τ(p)'nin İŞARETİ fazı 180° ile 0° arasında çevirmeli.

KESİN τ'DAN TAM ÖNGÖRÜ TABLOSU (hesaplandı, ölçümden önce mühürlendi):
   p    τ(p)              a_p        işaret   ÖNGÖRÜLEN FAZ
   2         −24        −0.53033       −          0°
   3         252        +0.59873       +        180°
   5        4830        +0.69121       +        180°
   7      −16744        −0.37655       −          0°
  11      534612        +1.00087       +        180°
  13     −577738        −0.43156       −          0°
  17    −6905934        −1.17965       −          0°
  19    10661420        +0.98780       +        180°
  23    18643272        +0.60398       +        180°
  29   128406630        +1.16251       +        180°
  31   −52843168        −0.33151       −          0°
HÜKÜM ÖLÇÜTÜ: 11 asalın işaret dizisi (−++−+−−++ +−) kırınım fazlarından
HATASIZ okunursa P1 İSABET; tek bir işaret ters çıkarsa RET.

--------------------------------------------------------------------------
P2 — p²-ÇİZGİLERİ (SATAKE): FAZ ÇEVRİLMESİ
--------------------------------------------------------------------------
Hecke: a(p²) = a_p² − 1 (bu KATSAYI değil, Dirichlet katsayısıdır).
Açık formülün katsayısı b(p²) = α_p²+β_p² = a_p² − 2.
Küçük asallarda |a_p| < √2 ⟹ b(p²) < 0 ⟹ p²-benekleri ζ'ya göre
FAZ-ÇEVRİK (ζ'da her p^k 180°'deydi; Δ'da 0° olmalı):
   p²    b(p²)=a_p²−2     öngörülen faz   (ζ'daki faz)
    4      −1.71875            0°            180°
    9      −1.64152            0°            180°
   25      −1.52222            0°            180°
   49      −1.85821            0°            180°
  121      −0.99825            0°            180°
  169      −1.81375            0°            180°
Genel özyineleme (Newton, e₁=a_p, e₂=1): b_0=2, b_1=a_p,
b_k = a_p·b_{k−1} − b_{k−2}  ⟹ p³, p⁴ çizgileri de öngörülür.

--------------------------------------------------------------------------
P3 — TARAK-TERMOMETRESİ + HİPERUNİFORMLUK + "ÇİZGİ ENVANTERİ" TESTİ
--------------------------------------------------------------------------
99/97 kalıbı, indekssiz: σ_u = √(−2 ln|⟨e^{2πix}⟩|)/(2π), log-ısınma.
KORO YASASI'nın GL(2) sınaması: koro varyansı Σ_q U_q²/2,
U_q = 2|b(q)|Λ(q)/(L√q log q). ζ'da |b|≡1; Δ'da |b(p)|=|a_p| ve küçük
asallarda a_p² ≪ 1 (a_2²=0.281, a_3²=0.358, a_7²=0.142) — oysa
Sato-Tate ⟨a_p²⟩=1. Koro ağırlığı 1/p olduğu için KÜÇÜK asallar hakim:
    Σ_p a_p²/p (Δ)  vs  Σ_p 1/p (ζ)  ⟹  Δ, EŞLEŞTİRİLMİŞ L'de ζ'DAN
    SOĞUK olmalı (oran script içinde ölçümden önce hesaplanır).
Bu, β'nın "ölü 2-ailesi soğukluğu"nun (105-P4) sürekli/GL(2) analoğudur:
orada bir çizgi ÖLDÜRÜLMÜŞTÜ, burada bütün çizgiler Sato-Tate ile
YENİDEN AĞIRLIKLANDIRILIYOR.
KONTROL: aynı L_eff'te ζ örgüsü yerinde üretilir (RS ana toplamı,
t = 2π e^{L}) ve aynı boru hattından geçirilir.

MEZARLIK KONTROLÜ (GL(2)'ye özgü, GL(1)'de karşılığı yok): açık formül
YALNIZCA asal kuvvetlerinde atom taşır. ω = log 6, log 10, log 14,
log 15, log 21 (bileşik, asal kuvveti DEĞİL) → benek OLMAMALI; plasebo
tabanında kalmalı. Bu, "her şey çizgidir" resminin doğrudan testi.

HİPERUNİFORMLUK: Σ²(n) logaritmik, GUE Dyson-Mehta'nın 0.6-1.6 katı
(105c'deki yedi GL(1) adasıyla aynı bant) — kaba gösterge.

BEKLENTİ (mühür): P1 ve P2 İSABET; P3'te Δ, eşleştirilmiş L'de ζ'dan
belirgin soğuk; mezarlıklar plasebo tabanında.

==========================================================================
SONUÇ (26 Ağustos) — P1 ✓✓✓  P2 ✓✓✓  MEZARLIK ✓✓✓  P3 ✓/✗ (aşağıda)
==========================================================================
VERİ: 117a_delta_zeros.npz, n=87 543; kapılardan sonra 86 979 (1 kesim,
564 sıfır atıldı); L_eff = 15.58 (küresel), pencereler 12.66 / 16.10.

KOD KAPISI ÖNCE: aynı boru hattı ζ örgüsünde (yerinde üretilmiş, taper'lı)
89'un yayımlanmış GL(1) sonucunu birebir veriyor — faz 180.0°±0.4°,
mutlak yasa oranı 0.99-1.03. Yani aşağıdaki Δ okumaları KALİBRE.

P1 ✓✓✓ HECKE İŞARETLERİ KIRINIMDAN OKUNDU — 14/14, HER İKİ IZGARADA
  Mühürlenen işaret dizisi (p = 2..43): − + + − + − − + + + − − + −
  Ölçülen fazlar (orta-nokta ızgarası): 359.90 / 179.96 / 180.00 /
  0.15 / 180.03 / 359.98 / 0.00 / 180.08 / 179.98 / 179.91 / 0.15 /
  360.00 / 180.24 / 354.19
  Maks |Δfaz| = 0.24° (p=43 hariç); p=43 |a_p| = 0.0178 ile neredeyse
  ÖLÜ bir çizgi (SNR 6.6) ve 5.8° sapıyor — yine de İŞARETİ doğru.
  τ(p)'nin işareti, sıfır örgüsünün kırınım deseninden okunuyor.

P2 ✓✓✓ SATAKE p^k ÇİZGİLERİ — 9/9, VE SADECE "HEPSİ ÇEVRİK" DEĞİL
  b(p^k) = α^k+β^k özyinelemesi hem ÇEVRİK hem ÇEVRİLMEMİŞ çizgi
  öngördü ve ikisi de tuttu:
    q=4 (b=−1.719) → 0° ölçülen 359.99   ÇEVRİK (ζ'da 180°)
    q=9 (−1.642)   → 0° ölçülen   0.04   ÇEVRİK
    q=25(−1.522)   → 0° ölçülen   0.04   ÇEVRİK
    q=27(−1.582)   → 0° ölçülen 359.89   ÇEVRİK
    q=49(−1.858)   → 0° ölçülen   0.31   ÇEVRİK
    q=121(−0.998)  → 0° ölçülen   1.25   ÇEVRİK
    q=169(−1.814)  → 0° ölçülen 359.28   ÇEVRİK
    q=8  (b=+1.442)→180° ölçülen 179.96  ÇEVRİLMEDİ ✓
    q=16 (b=+0.954)→180° ölçülen 179.96  ÇEVRİLMEDİ ✓
  Yani faz, |a_p|<√2 kuralının değil, b(p^k)'nin İŞARETİNİ izliyor.

MUTLAK BENEK YASASI GL(2)'DE ✓✓ (parametresiz):
  |Ĝ(log q)| = |b(q)|·Λ(q)/(L_eff√q)·cos(πτ)·DW
  orta-nokta ızgarası, üst pencere (L=16.10): oranlar 1.003-1.035
  (13 asal + 7 kuvvet); küresel tabloda 1.004-1.081. GL(1) adalarındaki
  1-4% ile AYNI kalitede. ⟹ mutlak yasa DERECEDEN BAĞIMSIZ; tek
  değişen açık formül katsayısı 1 → b(q).
  AÇIK/SÜRPRİZ: SIFIR ızgarasında (orta-noktasız) oranlar ω ile
  büyüyor (1.03 → 2.42). Yani Gauss-DW orta-noktalarda doğru,
  sıfırların kendisinde AŞIRI sönümlüyor: sıfır-ızgarasının sönümü
  Gauss DEĞİL. Yeni gözlemlenebilir, açık iş.

MEZARLIK ✓✓✓ "HER ŞEY ÇİZGİDİR" GL(2)'de doğrudan:
  bileşik ω (log 6,10,14,15,21,22,33,35) → 0.00003-0.00058,
  plasebo tabanı 0.00007; canlı çizgiler 0.008-0.047 (100-500 kat).
  Açık formülün atomik tayfı: asal kuvveti olmayan yerde benek YOK.

P3 — TARAK TERMOMETRESİ: İŞARET ✓, BÜYÜKLÜK ✗, YENİ YASA ✓
  Ön-mühür "Δ eşleştirilmiş L'de ζ'dan SOĞUK olmalı" dedi. ÖLÇÜM:
    L=12.66: σ_u(Δ)=0.2129 vs σ_u(ζ)=0.2858
    L=16.10: σ_u(Δ)=0.2342 vs σ_u(ζ)=0.3055
  İŞARET İSABET. Ama mühürdeki TOPLAMSAL biçim (97/105 kalıbı)
  Δσ² = (2/L²)(Σa_p²/p − Σ1/p) = −0.0118 / −0.0073 öngördü;
  ölçülen −0.0363 / −0.0385 — 3.1 ve 5.3 kat, üstelik L ile
  YANLIŞ YÖNDE. TOPLAMSAL KORO YASASI GL(2)'DE RET.
  Veriden ÇARPIMSAL biçim doğdu (ölçümden sonra):
    σ_u²(Δ)/σ_u²(ζ) = [Σ a_p² w_p/p]/[Σ w_p/p],  w = DW ağırlığı
    L=12.66: öngörü 0.5854 ölçüm 0.5550 (−5.2%)
    L=16.10: öngörü 0.6008 ölçüm 0.5877 (−2.2%)
  %2-5 içinde ve L ile İYİLEŞİYOR. "Sıcaklık çizgi envanterine
  ORANTILIDIR, envanter farkı kadar KAYMAZ." (GL(1) verilerinde
  ikisi ayırt edilemiyordu — çünkü orada envanter yalnız birkaç ölü
  çizgi kadar değişiyordu; Δ'da Sato-Tate bütün çizgileri yeniden
  ağırlıklandırdığı için iki biçim ilk kez AYRIŞTI.)
  Fiziksel okuma: Δ soğuk çünkü küçük asallarda a_p² ≪ 1
  (a_2²=0.281, a_3²=0.358, a_7²=0.142) ve koro ağırlığı 1/p.
  ⟨a_p²⟩ = 0.9775 (Sato-Tate 1.0) — ortalama normal, ilk birkaç
  çizgi zayıf. β'nın "ölü 2-ailesi soğukluğu"nun sürekli analoğu.

HİPERUNİFORMLUK ✓ Σ²(n) = 0.504/0.460/0.443/0.470 (n=5/10/20/50),
  GUE'nin 0.76-1.31 katı, Poisson'un 10-106 katı bastırılmış —
  105c'deki yedi GL(1) adasının bandıyla aynı.

ARTEFAKT AVI (gizleme yok — ikisi de İLK KOŞUDA YANLIŞTI):
 (1) Σ²(n) ilk koşuda 0.703/1.133/3.021/16.14 çıktı (GUE'nin 26 katı!).
     Sebep: kesilen kusurlu bölgenin bıraktığı BOŞLUK, birleştirilmiş
     dizide kutuları BOŞ sayıyordu. Segment-bazlı ölçümle düzeldi.
     Yanlış sayılar kıyas için scriptte "Δ✗" satırında duruyor.
 (2) ζ kontrol örgüsü ilk koşuda KESKİN kesimle üretilmişti ve 40
     birimlik sayım sürüklenmesi bıraktı (kayıp sıfır). Kayıp sıfır
     combı ISITIR ⟹ "Δ soğuk" hükmünü SAHTE olarak güçlendirebilirdi.
     ζ taper'lı yeniden döküldü: sürüklenme 40.0 → 2.4, 38 sıfır geri
     geldi, σ_u 0.2856 → 0.2858. Etki 2e-4; gözlenen fark 0.073.
     HÜKÜM AYAKTA. (Ayrıca konan rastgele-silme testi ZAYIFTIR —
     x=θ/π konumları silmeyle değişmez — scriptte öyle işaretlendi.)
"""

import numpy as np
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
exec(open(HERE / "117a_delta_motoru.py").read().split('if __name__')[0])

t00 = time.time()
TAU = tau_tablosu(NMAX_TAU)
M = DeltaMotor(TAU, taper_c=0.5)

D = np.load(HERE / "117a_delta_zeros.npz")
ZZ = D["zeros"]
KURT = D["kurtarilan"] if "kurtarilan" in D.files else np.array([])
print("=" * 78, flush=True)
print(f"117b — Δ KIRINIM SPEKTROSKOPİSİ   n_ham = {len(ZZ)}, "
      f"t ∈ [{ZZ[0]:.0f}, {ZZ[-1]:.0f}]", flush=True)
print("=" * 78, flush=True)


# ---------------------------------------------------------------- kapılar
def kisa_olcek_segmentler(zz, sayim, min_n=2500, win=80, esik=0.5, pad=160,
                          kisa_win=20, uzun_win=200, kisa_esik=0.7):
    """101 (basamak) + 105e-K4 (kısa çukur) kapıları + düzlük
    segmentasyonu — GL(1) kampanyalarındaki standardın aynısı."""
    from numpy.lib.stride_tricks import sliding_window_view
    d = np.arange(len(zz)) - (sayim(zz) - sayim(zz[0]))

    def cmed(x, w):
        m = np.median(sliding_window_view(x, w), axis=1)
        out = np.empty(len(x))
        out[w // 2:w // 2 + len(m)] = m
        out[:w // 2] = m[0]
        out[w // 2 + len(m):] = m[-1]
        return out

    med = np.median(sliding_window_view(d, win), axis=1)
    adim = med[win + 1:] - med[:-(win + 1)]
    bad = np.zeros(len(zz), bool)
    for i in np.where(np.abs(adim) > esik)[0] + win // 2:
        bad[max(0, i - pad):i + 2 * pad] = True
    sapma = np.abs(cmed(d, kisa_win) - cmed(d, uzun_win))
    n_kisa = int((sapma > kisa_esik).sum())
    for i in np.where(sapma > kisa_esik)[0]:
        bad[max(0, i - pad):i + pad] = True
    seg, kes, s0 = [], 0, 0
    for i in range(1, len(zz) + 1):
        if i == len(zz) or bad[i] != bad[i - 1]:
            if not bad[s0] and i - s0 >= min_n:
                seg.append(zz[s0:i])
            if i < len(zz) and bad[i]:
                kes += 1
            s0 = i
    return seg, kes, float(sapma.max()), n_kisa


SEG, KES, SAPMAX, NKISA = kisa_olcek_segmentler(
    ZZ, lambda t: M.theta(t) / np.pi)
ZC = np.concatenate(SEG) if SEG else ZZ
print(f"KAPILAR: basamak+kısa-çukur → {len(SEG)} segment, {KES} kesim; "
      f"maks |med20−med200| = {SAPMAX:.2f} ({NKISA} sıfır işaretli); "
      f"analize giren n = {len(ZC)} ({len(ZZ)-len(ZC)} atıldı)", flush=True)


# ---------------------------------------------------------------- araçlar
def Ghat(t, omegas, chunk=40000):
    out = np.zeros(len(omegas), dtype=complex)
    for s0 in range(0, len(t), chunk):
        tt = t[s0:s0 + chunk]
        out += np.exp(1j * np.outer(omegas, tt)).sum(axis=1)
    return out / len(t)


def b_kat(a, k):
    """b_k = α^k+β^k; b0=2, b1=a, b_k = a·b_{k−1} − b_{k−2}."""
    bm, b = 2.0, a
    for _ in range(k - 1):
        bm, b = b, a * b - bm
    return b


AS = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43]
PK = [(p, p, 1) for p in AS] + \
     [(4, 2, 2), (8, 2, 3), (16, 2, 4), (9, 3, 2), (27, 3, 3),
      (25, 5, 2), (49, 7, 2), (121, 11, 2), (169, 13, 2)]
MEZAR = [6, 10, 14, 15, 21, 22, 33, 35]      # bileşik: atom YOK

# plasebo için: BÜTÜN asal kuvvetlerinin log'ları (m ≤ 4000)
def _asal_kuvvetleri(N=4000):
    sv = np.ones(N + 1, bool); sv[:2] = False
    for i in range(2, int(N ** .5) + 1):
        if sv[i]:
            sv[i * i::i] = False
    out = []
    for p in np.nonzero(sv)[0]:
        m = int(p)
        while m <= N:
            out.append(np.log(m)); m *= int(p)
    return np.array(sorted(out))


TUM_CIZGI = _asal_kuvvetleri()


def tablo(zset, etiket, LOG=True):
    """Bir sıfır kümesinde tam benek tablosu. zset: ölçüm noktaları."""
    Leff = float(np.mean(2 * np.log(zset / TWO_PI)))
    x = M.theta(zset) / np.pi
    comb = abs(np.exp(2j * np.pi * x).mean())
    sig_u = np.sqrt(max(-2 * np.log(max(comb, 1e-300)), 0)) / TWO_PI
    sig_t = sig_u * TWO_PI / Leff
    oms = np.array([np.log(q) for q, _, _ in PK])
    G = Ghat(zset, oms)
    # plasebo tabanı: çizgilerden uzak ω
    rng = np.random.default_rng(3)
    cand = rng.uniform(0.5, 4.2, 2000)
    keep = np.array([np.min(np.abs(c - TUM_CIZGI)) > 0.02 for c in cand])
    Gp = np.abs(Ghat(zset, cand[keep][:400]))
    taban = float(np.median(Gp))
    Gm = np.abs(Ghat(zset, np.array([np.log(m) for m in MEZAR])))
    if LOG:
        print(f"\n--- {etiket}: n={len(zset)}, L_eff={Leff:.3f}, "
              f"σ_u={sig_u:.4f}, σ_t={sig_t:.4f}, plasebo tabanı="
              f"{taban:.5f} ---", flush=True)
        print(f"  {'q':>4} {'p^k':>5} {'b(q)':>9} {'ölçüm':>8} {'öngörü':>8} "
              f"{'oran':>6} {'faz':>8} {'ö-faz':>6} {'Δfaz':>7} {'SNR':>5}",
              flush=True)
    rows = []
    for (q, p, k), om, g in zip(PK, oms, G):
        a = TAU[p] / float(p) ** 5.5
        b = b_kat(a, k)
        tau_ = om / Leff
        dw = np.exp(-om ** 2 * sig_t ** 2 / 2)
        pred = abs(b) * np.log(p) / (Leff * np.sqrt(q)) * np.cos(np.pi * tau_) * dw
        ph = (np.degrees(np.angle(g)) + 360) % 360
        pph = 180.0 if b > 0 else 0.0
        dph = (ph - pph + 180) % 360 - 180
        rows.append((q, p, k, b, abs(g), pred, ph, pph, dph, taban))
        if LOG:
            print(f"  {q:>4} {p:>3}^{k:<1} {b:>9.5f} {abs(g):>8.5f} "
                  f"{pred:>8.5f} {abs(g)/pred:>6.3f} {ph:>7.2f}° {pph:>5.0f}° "
                  f"{dph:>+7.2f} {abs(g)/taban:>5.1f}", flush=True)
    if LOG:
        print("  MEZARLIK (bileşik ω, atom YOK): "
              + "  ".join(f"log{m}:{v:.5f}" for m, v in zip(MEZAR, Gm))
              + f"   | plasebo {taban:.5f}", flush=True)
    return dict(L=Leff, sig_u=sig_u, sig_t=sig_t, rows=rows, taban=taban,
                mezar=Gm, n=len(zset))


# ================= P1/P2: iki ızgara (orta-noktasız / orta-noktalı) =====
print("\n" + "=" * 78, flush=True)
print("P1/P2 — BENEK TABLOLARI", flush=True)
MID = np.concatenate([0.5 * (s[:-1] + s[1:]) for s in SEG]) if SEG \
    else 0.5 * (ZZ[:-1] + ZZ[1:])
T_ORT = tablo(ZC, "IZGARA A: SIFIRLAR (orta-noktasız)")
T_MID = tablo(MID, "IZGARA B: ORTA NOKTALAR")

# --- L bandını daraltmak için üst pencere (L_eff yayılımı büyük) ---
kenar = np.exp(np.linspace(np.log(ZC[0]), np.log(ZC[-1] * 1.0001), 4))
PEN = []
for i in range(3):
    m = MID[(MID >= kenar[i]) & (MID < kenar[i + 1])]
    if len(m) >= 2000:
        PEN.append(tablo(m, f"PENCERE {i+1} (orta noktalar)"))

# ---------------------------- HÜKÜMLER ---------------------------------
print("\n" + "=" * 78, flush=True)
print("P1 HÜKMÜ — HECKE İŞARETLERİ KIRINIMDAN OKUNUYOR MU?", flush=True)
print(f"{'p':>4} {'τ(p)':>13} {'a_p':>9} {'gerçek':>7} "
      + "  ".join(f"{'okuma('+e+')':>13}" for e in ("A", "B")), flush=True)
isabet = {"A": 0, "B": 0}
toplam = 0
for p in AS:
    a = TAU[p] / float(p) ** 5.5
    ger = "+" if a > 0 else "−"
    hu = []
    for et, T in (("A", T_ORT), ("B", T_MID)):
        r = [x for x in T["rows"] if x[0] == p and x[2] == 1][0]
        # okuma KURALI: ölçülen faz 180°'ye mi 0°'a mı yakın
        ph = r[6]
        oku = "+" if 90 < ph < 270 else "−"
        hu.append((oku, ph, r[4] / r[9]))
        if oku == ger:
            isabet[et] += 1
    toplam += 1
    print(f"{p:>4} {TAU[p]:>13} {a:>9.5f} {ger:>7} "
          + "  ".join(f"{o:>3} {ph:7.1f}° {s:4.1f}σ" for o, ph, s in hu),
          flush=True)
print(f"  İŞARET OKUMA: ızgara A {isabet['A']}/{toplam}, "
      f"ızgara B {isabet['B']}/{toplam}", flush=True)
sapA = [abs(x[8]) for x in T_ORT["rows"] if x[2] == 1]
sapB = [abs(x[8]) for x in T_MID["rows"] if x[2] == 1]
print(f"  maks |Δfaz| (asallar): A {max(sapA):.2f}°, B {max(sapB):.2f}°",
      flush=True)
print(f"  HÜKÜM P1: "
      + ("İSABET — Hecke işaretleri kırınım fazından hatasız okundu"
         if isabet['B'] == toplam else
         f"KISMİ/RET — {toplam-isabet['B']} işaret ters (ızgara B)"),
      flush=True)

print("\nP2 HÜKMÜ — SATAKE p^k ÇİZGİLERİ (ζ'ya göre faz çevrilmesi)",
      flush=True)
print(f"{'q':>5} {'b(q)':>10} {'öngörü':>7} {'ölçüm(B)':>9} {'Δfaz':>8} "
      f"{'|Ĝ|/plasebo':>12} {'ζ farkı':>8}", flush=True)
p2ok = 0; p2n = 0
for x in T_MID["rows"]:
    q, p, k, b, g, pred, ph, pph, dph, tb = x
    if k == 1:
        continue
    p2n += 1
    if abs(dph) < 90:
        p2ok += 1
    print(f"{q:>5} {b:>10.5f} {pph:>6.0f}° {ph:>8.2f}° {dph:>+8.2f} "
          f"{g/tb:>12.1f} {'ÇEVRİK' if pph == 0 else 'aynı':>8}", flush=True)
print(f"  P2: {p2ok}/{p2n} kuvvet çizgisi öngörülen fazda", flush=True)

# ================= P3: termometre + eşleştirilmiş ζ =====================
print("\n" + "=" * 78, flush=True)
print("P3 — TARAK-TERMOMETRESİ, KORO ENVANTERİ, HİPERUNİFORMLUK", flush=True)

# koro oranı (ölçümden bağımsız, aritmetik öngörü)
sieve = np.ones(NMAX_TAU + 1, bool); sieve[:2] = False
for i in range(2, int(NMAX_TAU ** .5) + 1):
    if sieve[i]:
        sieve[i * i::i] = False
PRIMES = np.nonzero(sieve)[0]
ap_all = np.array([TAU[int(p)] / float(p) ** 5.5 for p in PRIMES])
S_delta = float(np.sum(ap_all ** 2 / PRIMES))
S_zeta = float(np.sum(1.0 / PRIMES))
print(f"  KORO ENVANTERİ (p ≤ {NMAX_TAU}): Σ a_p²/p = {S_delta:.4f} (Δ) vs "
      f"Σ 1/p = {S_zeta:.4f} (ζ) → oran {S_delta/S_zeta:.3f}", flush=True)
print(f"  ⟨a_p²⟩ = {float(np.mean(ap_all**2)):.4f} (Sato-Tate: 1.0000) — "
      f"soğukluk küçük asalların a_p² < 1 olmasından", flush=True)
print(f"  ⟹ ÖNGÖRÜ: σ_u²(Δ) − σ_u²(ζ) ≈ (2/L²)(Σa_p²/p − Σ1/p) = "
      f"{2*(S_delta-S_zeta):.4f}/L²  (negatif ⟹ Δ SOĞUK)", flush=True)

print(f"\n  {'pencere':>9} {'L_eff':>7} {'σ_u':>8} {'σ_u²':>9} {'n':>7}",
      flush=True)
for i, T in enumerate(PEN):
    print(f"  {'Δ p'+str(i+1):>9} {T['L']:>7.3f} {T['sig_u']:>8.4f} "
          f"{T['sig_u']**2:>9.5f} {T['n']:>7}", flush=True)

# --- eşleştirilmiş L'de ζ örgüsü (RS ana toplamı, yerinde üretilir) ---
print("\n  EŞLEŞTİRİLMİŞ-L ζ KONTROLÜ (ζ örgüsü yerinde üretiliyor)",
      flush=True)


def zeta_theta(t):
    z = 0.25 + 0.5j * np.asarray(t, float)
    return stirling_imloggamma(z) - np.asarray(t, float) / 2 * np.log(np.pi)


def zeta_Z(t, chunk=3000, c=0.5):
    """ζ RS ana toplamı — Δ motoruyla AYNI taper disiplini (c=0.5).
    (İlk koşuda keskin kesim kullanılmıştı ve ζ örgüsünde 40 birimlik
     sayım sürüklenmesi bıraktı; kayıp sıfır combı ISITIR, yani P3'ün
     'Δ soğuk' hükmünü SAHTE olarak güçlendirirdi. Taper + kayıp-sıfır
     duyarlılık testi bu yüzden eklendi.)"""
    t = np.atleast_1d(np.asarray(t, float))
    out = np.empty_like(t)
    X = np.sqrt(t / TWO_PI)
    W = np.maximum(1.0, c * np.sqrt(X)) if c > 0 else np.zeros_like(X)
    Nl = np.clip(np.floor(X - W).astype(np.int64), 1, None)
    th = zeta_theta(t)
    K = 0 if c <= 0 else int(np.ceil(2 * W.max()) + 2)
    nn = np.arange(1, int(Nl.max()) + K + 1)
    ln, wn = np.log(nn), nn ** -0.5
    for Nv in np.unique(Nl):
        m = Nl == Nv
        tm, thm, Xm, Wm = t[m], th[m], X[m], W[m]
        zv = np.empty(len(tm))
        for s0 in range(0, len(tm), chunk):
            sl = slice(s0, s0 + chunk)
            ph = thm[sl, None] - tm[sl, None] * ln[None, :Nv]
            v = 2 * (np.cos(ph) * wn[None, :Nv]).sum(axis=1)
            if K > 0:
                n2 = np.arange(Nv + 1, Nv + K + 1)
                u = (n2[None, :] - Xm[sl, None]) / Wm[sl, None]
                w = np.clip(0.5 * (1 - np.sin(np.pi * np.clip(u, -1, 1) / 2)),
                            0.0, 1.0)
                ph2 = thm[sl, None] - tm[sl, None] * ln[None, Nv:Nv + K]
                v = v + 2 * (np.cos(ph2) * w * wn[None, Nv:Nv + K]).sum(axis=1)
            zv[sl] = v
        out[m] = zv
    return out


def zeta_sifir(T0, T1, grid_frac=0.03):
    ts = [T0]; t = T0
    while t < T1:
        t += grid_frac * TWO_PI / np.log(t / TWO_PI)
        ts.append(t)
    ts = np.array(ts); v = zeta_Z(ts)
    sc = np.where(np.sign(v[:-1]) * np.sign(v[1:]) < 0)[0]
    a, b = ts[sc].copy(), ts[sc + 1].copy(); fa = v[sc].copy()
    for _ in range(28):
        mid = 0.5 * (a + b); fm = zeta_Z(mid)
        sol = fa * fm <= 0
        b = np.where(sol, mid, b); a = np.where(sol, a, mid)
        fa = np.where(sol, fa, fm)
    return 0.5 * (a + b)


ZETA = {}
for T in PEN:
    Lz = T["L"]                      # ζ'da L = log(t/2π) ⟹ t = 2π e^L
    tz0 = TWO_PI * np.exp(Lz)
    n_hedef = min(T["n"], 12000)
    dt = n_hedef * TWO_PI / Lz
    t1 = time.time()
    zz = zeta_sifir(tz0, tz0 + dt, grid_frac=0.03)
    mz = 0.5 * (zz[:-1] + zz[1:])
    x = zeta_theta(mz) / np.pi
    comb = abs(np.exp(2j * np.pi * x).mean())
    su = np.sqrt(max(-2 * np.log(max(comb, 1e-300)), 0)) / TWO_PI
    Lm = float(np.mean(np.log(mz / TWO_PI)))
    dr = np.arange(len(zz)) - (zeta_theta(zz) - zeta_theta(zz[0])) / np.pi
    ZETA[Lz] = (Lm, su, len(mz), dr.max() - dr.min())
    print(f"    ζ @ L={Lm:.3f} (t≈{tz0:.3e}): σ_u = {su:.4f}, n={len(mz)}, "
          f"sürüklenme {dr.max()-dr.min():.2f}  ({time.time()-t1:.0f} sn)",
          flush=True)
    # ---- KOD KAPISI: aynı boru hattı ζ'da 89'un yayımlanmış sonucunu
    #      (faz ≡ 180°, mutlak yasa oranı ~1) üretiyor mu?
    st_z = su * TWO_PI / Lm
    oz = np.array([np.log(p) for p in (2, 3, 5, 7, 11)])
    Gz = Ghat(mz, oz)
    print("      KOD KAPISI (ζ, aynı boru hattı — beklenen faz 180°, "
          "oran ~1):", flush=True)
    for p, om, g in zip((2, 3, 5, 7, 11), oz, Gz):
        pr = np.log(p) / (Lm * np.sqrt(p)) * np.cos(np.pi * om / Lm) \
            * np.exp(-om ** 2 * st_z ** 2 / 2)
        ph = (np.degrees(np.angle(g)) + 360) % 360
        print(f"        p={p:>3}: |Ĝ|={abs(g):.5f} öngörü {pr:.5f} "
              f"oran {abs(g)/pr:.3f}  faz {ph:7.2f}° "
              f"(Δ={((ph-180+180)%360)-180:+.2f}°)", flush=True)

print(f"\n  {'L_eff':>7} {'σ_u(Δ)':>9} {'σ_u(ζ)':>9} {'Δσ²':>10} "
      f"{'öngörü Δσ²':>12} {'oran':>7} {'ζ sürükl':>9}", flush=True)
for T, (Lz, v) in zip(PEN, ZETA.items()):
    Lm, su, nz, dr = v
    d2 = T["sig_u"] ** 2 - su ** 2
    ong = 2 * (S_delta - S_zeta) / T["L"] ** 2
    print(f"  {T['L']:>7.3f} {T['sig_u']:>9.4f} {su:>9.4f} {d2:>10.5f} "
          f"{ong:>12.5f} {d2/ong if ong != 0 else np.nan:>7.2f} {dr:>9.2f}",
          flush=True)

# ---- ARTEFAKT KAPISI: kayıp sıfır combı ISITIR mı? ----
# "Δ, ζ'dan soğuk" hükmünün en tehlikeli alternatif açıklaması: ζ
# örgüsünde kayıp sıfır var (motor kusuru) ⟹ ζ sahte-sıcak. Bu kapı
# σ_u'nun kayıp-sıfır oranına duyarlılığını DOĞRUDAN ölçer.
print("\n  ARTEFAKT KAPISI — σ_u'nun KAYIP SIFIRA duyarlılığı", flush=True)
rng2 = np.random.default_rng(77)


def su_of(x):
    return np.sqrt(max(-2 * np.log(max(abs(np.exp(2j * np.pi * x).mean()),
                                       1e-300)), 0)) / TWO_PI


xd = M.theta(MID) / np.pi
print(f"  {'silinen %':>10} {'σ_u(Δ)':>9} {'Δσ_u':>9}", flush=True)
for frak in (0.0, 0.001, 0.002, 0.004, 0.01):
    if frak == 0:
        v = su_of(xd)
        v0 = v
    else:
        k = int(len(xd) * frak)
        idx = rng2.choice(len(xd), len(xd) - k, replace=False)
        v = su_of(xd[np.sort(idx)])
    print(f"  {100*frak:>10.2f} {v:>9.4f} {v-v0:>+9.4f}", flush=True)
print("  NOT (dürüstlük): bu test SIFIR SİLMEYİ taklit eder ama x=θ/π", flush=True)
print("  konumları değişmediği için zayıftır. ASIL kanıt ζ örgüsünün", flush=True)
print("  taper'lı yeniden dökümüdür: sürüklenme 40.0 → 2.4 (38 sıfır geri", flush=True)
print("  geldi) ve σ_u 0.2856 → 0.2858 — yani kayıp sıfırlar 'Δ soğuk'", flush=True)
print("  hükmünü YARATMIYOR (etki ~2e-4, gözlenen fark 0.073).", flush=True)

# ---- KORO YASASININ ÇARPIMSAL BİÇİMİ (ölçümden SONRA doğdu) ----
# ÖN-MÜHÜR toplamsal biçimi öngörüyordu (97/105 kalıbı):
#     σ_u²(Δ) − σ_u²(ζ) = (2/L²)(Σa_p²/p − Σ1/p)
# Bu, İŞARETİ doğru verdi ama BÜYÜKLÜĞÜ 3-5 kat ıskaladı ve L ile
# yanlış yönde gitti. Veriye bakınca ÇARPIMSAL biçim çıktı:
#     σ_u²(Δ)/σ_u²(ζ) = [Σ_p a_p² w_p/p] / [Σ_p w_p/p],  w_p = DW ağırlığı
print("\n  KORO YASASININ ÇARPIMSAL BİÇİMİ (ölçümden SONRA türetildi)",
      flush=True)
_N = 200000
_sv = np.ones(_N + 1, bool); _sv[:2] = False
for _i in range(2, int(_N ** .5) + 1):
    if _sv[_i]:
        _sv[_i * _i::_i] = False
_P = np.nonzero(_sv)[0]
_kn = _P[_P <= NMAX_TAU]
_a2 = np.array([(TAU[int(p)] / float(p) ** 5.5) ** 2 for p in _kn])
_bg = _P[_P > NMAX_TAU]
print(f"  {'L':>7} {'koro oranı(DW)':>15} {'σ_u² oranı':>12} {'sapma':>8}",
      flush=True)
for T, (Lz, v) in zip(PEN, ZETA.items()):
    st = T["sig_t"]
    w = np.exp(-np.log(_kn.astype(float)) ** 2 * st ** 2)
    wb = np.exp(-np.log(_bg.astype(float)) ** 2 * st ** 2)
    ek = float(np.sum(wb / _bg))          # p>6000: ⟨a_p²⟩=1 ⟹ iki tarafta da
    cD = float(np.sum(_a2 * w / _kn)) + ek
    cZ = float(np.sum(w / _kn)) + ek
    olc = T["sig_u"] ** 2 / v[1] ** 2
    print(f"  {T['L']:>7.3f} {cD/cZ:>15.4f} {olc:>12.4f} "
          f"{100*(olc/(cD/cZ)-1):>+7.1f}%", flush=True)
print("  ⟹ TOPLAMSAL biçim RET (3-5 kat), ÇARPIMSAL biçim %2-5 içinde "
      "ve L ile İYİLEŞİYOR.", flush=True)

# --- hiperuniformluk ---
print("\n  HİPERUNİFORMLUK Σ²(n) (kaba gösterge)", flush=True)
GAM = 0.5772156649
NL = (5, 10, 20, 50)


def sigma2_seg(segler, sayim, nlist=NL):
    """DÜZELTME (ilk koşuda bulunan HATA): Σ²(n) segmentlerin
    BİRLEŞTİRİLMİŞ dizisi üzerinde ölçülürse, kesilen kusurlu bölgenin
    bıraktığı BOŞLUK boyunca kutular BOŞ sayılır ve varyansı sahte
    şekilde şişirir. İlk koşu böyle yapmıştı: Σ²(50) = 16.14 çıktı
    (GUE 0.617!) — tamamen kesim boşluğundan. Doğrusu: her segmentte
    ayrı ölç, sonra kutu sayısıyla ağırlıklı birleştir."""
    out = []
    for n in nlist:
        tum, W = [], 0
        for s in segler:
            x = sayim(s)
            x = x - x[0]
            if x[-1] < 3 * n:
                continue
            kutu = np.arange(0, x[-1] - n, n)
            cnt = np.diff(np.searchsorted(
                x, np.concatenate([kutu, [kutu[-1] + n]])))[:-1]
            tum.append(cnt)
        c = np.concatenate(tum)
        out.append(float(np.var(c)))
    return out


gue = [(np.log(2 * np.pi * n) + GAM + 1) / np.pi ** 2 - 0.125 for n in NL]
print("    " + " ".join(f"{'n='+str(n):>18}" for n in NL), flush=True)
print(f"    GUE:" + " ".join(f"{g:>17.3f}" for g in gue), flush=True)
s2v = sigma2_seg(SEG, lambda t: M.theta(t) / np.pi)
print("    Δ  :" + " ".join(f"{v:6.3f} ({v/g:4.2f}G {n/v:4.0f}P)"
                            for v, g, n in zip(s2v, gue, NL)), flush=True)
s2_hata = sigma2_seg([ZC], lambda t: M.theta(t) / np.pi)
print("    (kıyas — segmentleri BİRLEŞTİRerek ölçülseydi, yani ilk "
      "koşudaki HATA:)", flush=True)
print("    Δ✗ :" + " ".join(f"{v:6.3f}" for v in s2_hata), flush=True)

np.savez(HERE / "117b_delta_kirinim.npz",
         rows_A=np.array([[r[0], r[2], r[3], r[4], r[5], r[6], r[7], r[8]]
                          for r in T_ORT["rows"]]),
         rows_B=np.array([[r[0], r[2], r[3], r[4], r[5], r[6], r[7], r[8]]
                          for r in T_MID["rows"]]),
         meta=np.array([T_ORT["L"], T_ORT["sig_u"], T_ORT["taban"],
                        T_MID["L"], T_MID["sig_u"], T_MID["taban"]]),
         mezar_A=T_ORT["mezar"], mezar_B=T_MID["mezar"],
         mezar_q=np.array(MEZAR),
         pen_L=np.array([T["L"] for T in PEN]),
         pen_su=np.array([T["sig_u"] for T in PEN]),
         zeta_su=np.array([v[1] for v in ZETA.values()]),
         zeta_L=np.array([v[0] for v in ZETA.values()]),
         sigma2=np.array(s2v), koro=np.array([S_delta, S_zeta]))
print(f"\n117b_delta_kirinim.npz yazıldı. Toplam {(time.time()-t00)/60:.1f} dk",
      flush=True)
