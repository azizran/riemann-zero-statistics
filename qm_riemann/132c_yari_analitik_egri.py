"""
132c — YARI-ANALİTİK EĞRİ: R(τ) SÜREKLİ Mİ, ÇİZGİYE Mİ ÖZGÜ? (27 Ağu)
==========================================================================
132a: Taylor serisi ıraksak (rms S'/ρ̄ = 1.17), kesilemez.
132b: yasa "kuyruğun, p-dalgasınca oynatılmış noktalarda örneklenişi"
      kanalından doğuyor (C1 yaşıyor, C2 ölü, C4p tek çizgiyle %66-92).
Bu lab kalemin YAZACAĞI TOPLAMI arıyor — ve bir KESİN ÖZDEŞLİKle
başlıyor (pertürbasyon YOK):

  ρ̄ ḡ_t = 1 (ızgara N̄'yi çevirerek kurulur) ⟹
  ds_n = ρ̄(z_{n+1}−z_n) − 1 = S(z_n) − S(z_{n+1})
       = Σ_Q 2 a_Q sin(ω_Q g_n/2) · cos(ω_Q m_n)          (★)
  a_Q = Λ(Q)/(π√Q logQ), ω_Q = logQ, g_n = z_{n+1}−z_n,
  m_n = (z_n+z_{n+1})/2,  ve ω_Q ḡ/2 = π·τ_Q  (τ_Q = ω_Q/L̄).

(★) TEK TOPLAMDIR — çift toplam değil. Bütün τ-trigonometrisi buradan:
  • ÇİZGİ GENLİĞİ:  A₁(q) = 2 a_q sin(π τ_q)            [sin(πτ)]
  • İKİ-NOKTA ÇEKİRDEĞİ: η_n ve η_{n+1} orta noktaları mm ∓ ḡ/2'de;
    p-dalgası her iki boşluğu da modüle eder ⟹ cos(ω_p ḡ/2) = cos(πτ)
    ortak çarpanı çıkar:
        R_nn(τ) = κ · cos(π τ),   κ = ḡ (∂C/∂g)/C(ḡ)
        C(g) = Σ_{Q∈kuyruk} 2 a_Q² sin²(ω_Q g/2) cos(ω_Q g)
        V(g) = Σ_{Q∈kuyruk} 2 a_Q² sin²(ω_Q g/2)   [= σ_η², R_p=ḡV'/2V]
  Yani "−2cos(πτ)"nin cos(πτ)'si KİNEMATİK (iki-nokta geometrisi),
  "−2" ise kuyruk çekirdeğinin ḡ'deki logaritmik türevi.

SINAVLAR:
  Y1  ÖZDEŞLİK: (★) yakınsamış lab'da nokta-nokta doğrulanır mı?
  Y2  GENLİK YASASI: ölçülen A₁(q) ?= 2a_q⟨sin(ω_q g_n/2)⟩ (19 çizgi).
  Y3  τ-TARAMASI: R_p(τ), R_nn(τ) tüm fit çizgilerinde (τ=0.10..0.52);
      R_nn(τ)/cos(πτ) SABİT mi (çarpanlaşma)? δ(τ) = R_nn/(−2cosπτ)−1.
  Y4  SÜREKLİLİK: ÇİZGİ-DIŞI boyalı sondaj frekanslarında (τ_x, hiçbir
      logq/L̄'ye denk değil) aynı yasa okunuyor mu? 126'nın "seçicilik"i
      vuruş-defterinden geliyordu: o defter R(τ_x)≈0 der. (★) defteri
      ise SÜREKLİ der — çarpanlaşma kinematik. AYIRT EDİCİ SINAV.
  Y5  ÇEKİRDEK: κ_kuram (sıfır parametre) vs κ_ölçüm.

SONUÇ (27 Ağustos) — Y1✓ Y2✓ Y3✓ Y4 ÇARPICI Y5~ Y6 AYIRICI:
lab: L̄=8.1113, ḡ=0.7771 (2π/L̄=0.7746 ✓), σ_η=0.1634, c₁/σ²=−0.572.

Y1 ÖZDEŞLİK ✓: korelasyon(ds, S(z_n)−S(z_{n+1})) = 0.9753; çarpım
  biçimi Σ2a·sin(ωg/2)cos(ωm) aynı korelasyon, rms oranı 1.024.
  (%22'lik artık: ρ'nun boşluk boyunca değişimi + n_iter=5'in tam sabit
  nokta olmaması.) (★) PERTÜRBASYON DEĞİL, KİMLİK — kalemin zemini.

Y2 GENLİK YASASI ✓✓ SIFIR PARAMETRE, 28 ÇİZGİ:
  A₁ölçüm / [2a_q⟨sin(ω_q g/2)⟩] = %84-108 (asallarda τ ile düzgün
  %97→%86 düşüş = öz-tutarlı bastırma). 2a·sin(πτ) ile 2a·⟨sin(ωg/2)⟩
  zaten %1-5 içinde. ⟹ sin(πτ) yarısı ÖLÇÜLDÜ.

Y3 τ-TARAMASI ✓✓ — "−2" SABİT OLARAK AYRIŞTIRILDI:
  κ_ölçüm = R_nn/cos(πτ):
    ASAL & τ≤0.25 : −1.996 ± 0.086  (%4.3 saçılma, n=4)  ← −2 !
    ASAL & τ≤0.35 : −2.093 ± 0.158  (%7.5, n=7)
    τ>0.35'te çarpanlaşma KOPUYOR: R_nn −2cosπτ'yi geçiyor ve
    τ→0.5'te sıfırlanmıyor (R_nn ≈ −0.5 platosu; hedef 0'a gider).
  δ(τ) = R_nn/(−2cosπτ) − 1 (asallar): −0.061(.085) −0.022(.135)
    +0.021(.198) +0.053(.240) +0.031(.296) +0.103(.316) +0.199(.349)
    +0.279(.363) +0.399(.387) … ⟹ LAB'DA δ, τ≳0.26'da POZİTİF.
    129'un GERÇEK VERİDEKİ δ(0.26-0.30) ≈ −0.10'u ile İŞARET ZIT —
    açık uyuşmazlık (lab kuyruğu τ≤1.4'te budalı; gerçek kuyruk sonsuz).
    ARTEFAKT ŞÜPHESİ KAYITLI, kalem bunu δ'nın kaynağı olarak kullanmaz.
  Asal KUVVETLERİ (4,8,9,16,25,27,32,49,64) sapıyor (A₁ küçük → gürültü).

Y5 ÇEKİRDEK (kuyruk TEK toplamı, boşluk-modülasyonu kanalı, 0 parametre):
  V(ḡ)=Σ2a²sin²(ωḡ/2) → σ_η_kuram 0.272 [ölçüm 0.163 — 1.7× fazla]
  c₁/σ²_kuram = −0.4427 [ölçüm −0.5719, gerçek veri −0.52..−0.59] ✓ %23
  R_p_kuram = ḡV'/(2V) = −0.813 [ölçüm −0.86..−0.94] ✓ %6-14 — ÇARPICI
  κ_kuram = ḡ(∂C/∂g_n)/C = −3.22 [ölçüm −2.00] — 1.6× fazla.
  ⟹ Tek-toplam çekirdeği R_p'yi ve c₁/σ²'yi neredeyse doğru veriyor;
  κ'nın büyüklüğünü fazla veriyor (öz-tutarlı bastırma hesaba katılmalı).

Y4 ÇİZGİ-DIŞI BOYALI SONDAJLAR — 126'NIN SEÇİCİLİĞİ TAM LAB'DA SAĞLAM:
  9 çizgi-dışı τ_x'te (asal-benzeri ağırlıkla boyanmış, tam öz-tutarlı
  lab): R_nn = −0.115/−0.055/+0.281/+0.150/+0.234/+0.221/+0.473/+0.129/
  +0.231 — YASA YOK (hedefin %3-6'sı ya da İŞARET TERS); R_p ≈ +0.5..+0.7.
  Üstelik A₁ölçüm/A₁kuram = %48-50 (dokuzunda da!) — boyalı dalga
  ÖZDEŞLİĞİN genlik yasasına da uymuyor; gerçek çizgide %84-108.
  ÇİZGİ-İÇİ kontrol (var olan çizgi ikizlenir, ağırlık ×2): R_nn gerçek
  çizginin %65/%59/%57'sine düşüyor — tam da "boyalı pay yalnız A₁'i
  şişirir, R üretmez" öngörüsünün verdiği oran (A₁gerçek/A₁toplam=0.65;
  öngörü −1.157/−0.808/−0.639, ölçüm −1.231/−0.728/−0.558; %6-13) ✓✓
  ⟹ ÇİZGİ-KİLİDİ NİCEL: R yalnızca ARİTMETİK DESTEKLİ genlikten doğar.

Y6 KANAL ÇAPRAZI — SEÇİCİLİK KANALDA DEĞİL, İPTALDE:
  132b'nin C4 kanalı (kuyruk YALNIZ tek bir dalganın oynattığı noktalarda
  örneklenir) BOYALI dalgayla da çalışıyor:
    çizgi (q=3/11/23):  R_nn = −1.177/−1.346/−1.267
    boyalı (τ=.17/.30/.38): R_nn = −1.124/−1.147/−0.962
  Neredeyse AYNI. ⟹ C4 kanalı KİNEMATİK ve süreklidir; tek başına
  seçici DEĞİLDİR. Tam lab'da çizgi-dışı frekansın R'si sıfıra
  gidiyorsa, bunu yapan şey KARŞIT İŞARETLİ İKİNCİ BİR KANALIN
  İPTALİDİR — ve bu iptal, aritmetik desteği olan frekanslarda
  TAMAMLANMAZ. Faz çaprazıyla da tutarlı: 132b'de FZ-C3 (kinematik
  kanal) faz karıştırmayı %55-65 sağ atlatıyor, ama FZ-C0 (tam lab)
  ölüyor — yani faz-duyarlı olan iptali bozan aritmetik artıktır.
  130'un kule-reti ile birlikte: iptali bozan şey TAM rezonans
  (Q=pQ' — yalnız aynı-taban kuleleri, ihmal) OLAMAZ; geriye
  YAKIN-rezonans aritmetik artık kalıyor — kalemin bir sonraki hedefi.

ÇİZİM: 132c_egri.png (çizgiler, asal kuvvetleri, çizgi-dışı sondajlar,
κ_ölçüm·cosπτ ve −2cosπτ üst üste).
"""

import numpy as np
import time
from sympy import primerange, factorint

TWO_PI = 2 * np.pi
L0 = 7.0
NZ = 40000
N_ITER = 5


def pk_list(lim):
    out = []
    for p in primerange(2, int(lim) + 1):
        q = p
        while q <= lim:
            out.append(q); q *= p
    return sorted(set(out))


def rvm_N(t):
    x = t / TWO_PI
    return x * np.log(x / np.e) + 7 / 8


t0 = TWO_PI * np.exp(L0)
idx = np.arange(NZ + 1, dtype=float)
tg = t0 + idx * TWO_PI / L0
for _ in range(8):
    tg = tg - (rvm_N(tg) - rvm_N(t0) - idx) / (np.log(tg / TWO_PI) / TWO_PI)
rho = np.log(tg / TWO_PI) / TWO_PI

QCAP = int(np.exp(1.4 * L0))
QS = pk_list(QCAP)
QARR = np.array(QS, float)
LOGQ = np.log(QARR)
LAMv = np.array([float(np.log(list(factorint(q).items())[0][0])) for q in QS])
W0 = LAMv / (np.pi * np.sqrt(QARR) * LOGQ)
KUY0 = LOGQ / L0 > 0.52


def S_of(t, lq, w, chunk=8000):
    s = np.zeros_like(t)
    if w.size == 0:
        return s
    for a in range(0, len(t), chunk):
        tt = t[a:a + chunk]
        s[a:a + chunk] = -(np.sin(np.outer(tt, lq)) * w[None, :]).sum(axis=1)
    return s


def lab(lq, w, n_iter=N_ITER):
    u = np.zeros_like(tg)
    for _ in range(n_iter):
        u = 0.5 * u + 0.5 * (-S_of(tg + u, lq, w) / rho)
    return np.sort(tg + u)


def olc_ext(z, ekstra=()):
    """131 zinciri + isteğe bağlı çizgi-dışı sondaj frekansları."""
    g = np.diff(z)
    m = 0.5 * (z[:-1] + z[1:])
    Lw = np.log(m / TWO_PI)
    Lb = float(Lw.mean())
    ds = g * Lw / TWO_PI - 1
    qs_fit = pk_list(min(int(np.exp(0.52 * Lb)), 720))
    # ekstra: (etiket, ω) — AÇISAL FREKANS DOĞRUDAN verilir; τ'dan
    # yeniden türetilmez (lab'ın L̄'si ile fit'in L̄'si arasındaki 10⁻³'lük
    # fark, pencere çözünürlüğü 2π/T ≈ 1.7e−4 mertebesinde olduğu için
    # sondajı ıskalatır — 27 Ağu'da yakalanan tuzak).
    etiket = [("q", q, np.log(q)) for q in qs_fit] + \
             [("x", ad, wx) for ad, wx in ekstra]
    freqs = [e[2] for e in etiket]
    tt = (m - m.mean()) / (m[-1] - m[0])
    cols = [np.ones_like(m), tt, tt**2]
    for f in freqs:
        cols += [np.cos(f * m), np.sin(f * m)]
    X = np.vstack(cols).T
    b1, *_ = np.linalg.lstsq(X, ds, rcond=None)
    eta = ds - X @ b1
    s2 = float((eta**2).mean())
    b2, *_ = np.linalg.lstsq(X, eta**2 - s2, rcond=None)
    ee = eta[:-1] * eta[1:]
    c1 = float(ee.mean())
    mm = 0.5 * (m[:-1] + m[1:])
    tt2 = (mm - mm.mean()) / (mm[-1] - mm[0])
    cols2 = [np.ones_like(mm), tt2, tt2**2]
    for f in freqs:
        cols2 += [np.cos(f * mm), np.sin(f * mm)]
    X2 = np.vstack(cols2).T
    b3, *_ = np.linalg.lstsq(X2, ee - c1, rcond=None)
    out = []
    for i, (kind, lab_, f) in enumerate(etiket):
        cA, sA = b1[3 + 2 * i], b1[3 + 2 * i + 1]
        A1 = np.hypot(cA, sA); ph = np.arctan2(sA, cA)
        Rp = (b2[3 + 2 * i] * np.cos(ph) +
              b2[3 + 2 * i + 1] * np.sin(ph)) / (2 * s2 * A1)
        Rn = (b3[3 + 2 * i] * np.cos(ph) +
              b3[3 + 2 * i + 1] * np.sin(ph)) / (2 * c1 * A1)
        out.append(dict(kind=kind, lab=lab_, f=f, tau=f / Lb,
                        A1=A1, Rp=Rp, Rn=Rn))
    return out, np.sqrt(s2), c1 / s2, Lb, g, m, ds, eta


# ================= YAKINSAMIŞ LAB ===================================
ts = time.time()
z = lab(LOGQ, W0)
out, se, c1r, Lb, g, m, ds, eta = olc_ext(z)
gbar = float(g.mean())
print(f"lab: {len(QS)} çizgi, n_iter={N_ITER}, {time.time()-ts:.0f} sn", flush=True)
print(f"L̄ = {Lb:.4f}   ḡ = {gbar:.4f}   2π/L̄ = {TWO_PI/Lb:.4f}   "
      f"σ_η = {se:.4f}   c₁/σ² = {c1r:+.3f}", flush=True)

# ---- Y1: KESİN ÖZDEŞLİK (★) ----------------------------------------
print("\n########## Y1: KESİN ÖZDEŞLİK ds = S(z_n) − S(z_{n+1}) ##########",
      flush=True)
Sz = S_of(z, LOGQ, W0)
ds_ozd = Sz[:-1] - Sz[1:]


def korel(a, b):
    a = a - a.mean(); b = b - b.mean()
    return float((a * b).mean() / np.sqrt((a**2).mean() * (b**2).mean()))


print(f"  rms(ds ölçüm) = {np.sqrt((ds**2).mean()):.5f}   "
      f"rms(ds özdeşlik) = {np.sqrt((ds_ozd**2).mean()):.5f}")
print(f"  korelasyon = {korel(ds, ds_ozd):.6f}   "
      f"rms(fark)/rms(ds) = {np.sqrt(((ds-ds_ozd)**2).mean())/np.sqrt((ds**2).mean()):.4f}",
      flush=True)
# çarpım biçimi: Σ 2a sin(ωg/2) cos(ωm)
acc = np.zeros_like(m)
for a0 in range(0, len(m), 4000):
    mm_ = m[a0:a0 + 4000]; gg = g[a0:a0 + 4000]
    acc[a0:a0 + 4000] = (2 * W0[None, :]
                         * np.sin(np.outer(gg, LOGQ) / 2)
                         * np.cos(np.outer(mm_, LOGQ))).sum(axis=1)
print(f"  ÇARPIM BİÇİMİ Σ2a·sin(ωg/2)cos(ωm): korelasyon(ds) = {korel(ds, acc):.6f}"
      f"   rms oranı = {np.sqrt((acc**2).mean())/np.sqrt((ds**2).mean()):.4f}",
      flush=True)

# ---- Y2: GENLİK YASASI A₁ = 2a_q ⟨sin(ω_q g/2)⟩ ---------------------
print("\n########## Y2: GENLİK YASASI  A₁(q) = 2a_q·⟨sin(ω_q g/2)⟩ ##########",
      flush=True)
print(f"{'q':>4} {'τ':>6} {'A₁ ölçüm':>10} {'2a·sin(πτ)':>11} "
      f"{'2a·⟨sin(ωg/2)⟩':>15} {'ölç/kuram':>10}", flush=True)
for r in out:
    q = r["lab"]; i = QS.index(q); w = W0[i]; f = r["f"]
    pred_s = 2 * w * np.sin(np.pi * r["tau"])
    pred_g = 2 * w * float(np.sin(f * g / 2).mean())
    print(f"{q:>4} {r['tau']:>6.3f} {r['A1']:>10.4f} {pred_s:>11.4f} "
          f"{pred_g:>15.4f} {r['A1']/pred_g:>9.1%}", flush=True)

# ---- Y5: ÇEKİRDEK (sıfır parametre) --------------------------------
KUYF = LOGQ / Lb > 0.52          # fit tabanının dışı (η'nın taşıyıcısı)
wk = LOGQ[KUYF]; ak = W0[KUYF]
h = 1e-5


def Vf(gg):
    return float(np.sum(2 * ak**2 * np.sin(wk * gg / 2)**2))


def Cf(gn, gp, D):
    return float(np.sum(2 * ak**2 * np.sin(wk * gn / 2) *
                        np.sin(wk * gp / 2) * np.cos(wk * D)))


Vg = Vf(gbar)
Vp = (Vf(gbar + h) - Vf(gbar - h)) / (2 * h)
Cg = Cf(gbar, gbar, gbar)
dCn = (Cf(gbar + h, gbar, gbar + h / 2) -
       Cf(gbar - h, gbar, gbar - h / 2)) / (2 * h)
kap_kuram = gbar * dCn / Cg
print("\n########## Y5: ÇEKİRDEK (sıfır parametre, kuyruk tek toplamı) ##########")
print(f"  V(ḡ) = σ_η²_kuram = {Vg:.5f}  (σ = {np.sqrt(Vg):.4f})  "
      f"[ÖLÇÜM {se**2:.5f} / {se:.4f}]")
print(f"  C(ḡ) = c₁_kuram   = {Cg:.5f}   c₁/σ²_kuram = {Cg/Vg:+.4f}  "
      f"[ÖLÇÜM {c1r:+.4f}]")
print(f"  R_p_kuram = ḡV'/(2V) = {gbar*Vp/(2*Vg):+.4f}")
print(f"  κ_kuram   = ḡ(∂C/∂g_n)/C = {kap_kuram:+.4f}   "
      f"(hedef yasa κ = −2)", flush=True)

# ---- Y3: τ-TARAMASI (çizgiler) -------------------------------------
print("\n########## Y3: τ-TARAMASI — ÇİZGİLER ##########", flush=True)
print(f"{'q':>4} {'τ':>6} {'R_p':>8} {'R_nn':>8} {'−2cosπτ':>9} "
      f"{'R_nn/hedef':>11} {'δ(τ)':>8} {'R_nn/cosπτ':>11}", flush=True)
ASAL = set(primerange(2, 4000))
tl, rl, kl, is_p = [], [], [], []
for r in out:
    tau = r["tau"]; hed = -2 * np.cos(np.pi * tau)
    kk = r["Rn"] / np.cos(np.pi * tau)
    tl.append(tau); rl.append(r["Rn"]); kl.append(kk)
    is_p.append(r["lab"] in ASAL)
    print(f"{r['lab']:>4}{'' if r['lab'] in ASAL else '*'} {tau:>6.3f} "
          f"{r['Rp']:>+8.3f} {r['Rn']:>+8.3f} "
          f"{hed:>+9.3f} {r['Rn']/hed:>10.1%} {r['Rn']/hed-1:>+8.3f} "
          f"{kk:>+11.3f}", flush=True)
kl = np.array(kl); tl = np.array(tl); rl = np.array(rl)
is_p = np.array(is_p)
print("  (* = asal kuvveti; A₁ küçük olduğu için R gürültülü)")
for ad, msk in [("TÜM çizgiler", np.ones_like(is_p)),
                ("yalnız ASALLAR", is_p),
                ("ASAL & τ≤0.35", is_p & (tl <= 0.35)),
                ("ASAL & τ≤0.25", is_p & (tl <= 0.25))]:
    v = kl[msk]
    print(f"  κ_ölçüm [{ad:>14}] n={msk.sum():>2}  ort {v.mean():+.3f}  "
          f"std {v.std():.3f}  ({v.std()/abs(v.mean()):.1%} saçılma)", flush=True)

# ---- Y4: ÇİZGİ-DIŞI BOYALI SONDAJLAR -------------------------------
print("\n########## Y4: ÇİZGİ-DIŞI BOYALI SONDAJLAR (süreklilik sınavı) ##########",
      flush=True)
# Boyalı sondaj: S alanına ω_x frekanslı, ASAL-BENZERİ ağırlıklı bir çizgi
# eklenir (a_x = 1/(π e^{ω_x/2}) — Q=e^{ω_x} asalmış gibi). İki tür:
#   DIŞI : ω_x hiçbir log q'ya denk değil  → 126'nın "seçicilik" sınavı
#   İÇİ  : ω_x = log q (var olan çizgi ikizlenir) → makine kontrolü
TX_DIS = [0.12, 0.17, 0.22, 0.26, 0.30, 0.34, 0.38, 0.44, 0.50]
Q_ICI = [3, 11, 23]
cizgi_w = LOGQ.copy()
sondalar = [("dış", tx, tx * Lb) for tx in TX_DIS] + \
           [("iç", np.log(q) / Lb, float(np.log(q))) for q in Q_ICI]
print(f"{'tür':>4} {'τ_x':>6} {'Δω(en yakın çizgi)':>19} {'A₁ ölçüm':>9} "
      f"{'A₁ kuram':>9} {'oran':>7} {'R_p':>8} {'R_nn':>8} {'−2cosπτ':>9} "
      f"{'R/hedef':>8} {'σ_η':>7}", flush=True)
sonda = []
for tur, tx, wx in sondalar:
    ax = 1.0 / (np.pi * np.exp(wx / 2))
    lq2 = np.append(LOGQ, wx); w2 = np.append(W0, ax)
    zx = lab(lq2, w2)
    ox, sex, _, Lbx, gx, *_ = olc_ext(zx, ekstra=[("x", wx)])
    r = ox[-1]
    tau_x = wx / Lbx
    hed = -2 * np.cos(np.pi * tau_x)
    dmin = float(np.abs(cizgi_w - wx).min())
    A1k = 2 * ax * float(np.sin(wx * gx / 2).mean())
    if tur == "dış":
        sonda.append((tau_x, r["Rn"]))
    print(f"{tur:>4} {tau_x:>6.3f} {dmin:>19.4f} {r['A1']:>9.4f} {A1k:>9.4f} "
          f"{r['A1']/A1k:>6.1%} {r['Rp']:>+8.3f} {r['Rn']:>+8.3f} "
          f"{hed:>+9.3f} {r['Rn']/hed:>7.1%} {sex:>7.4f}", flush=True)

# ---- Y6: KANAL ÇAPRAZI — 132b'nin C4 hücresi BOYALI dalgayla ---------
# Y4 "yasa çizgiye özgü" dedi; 132b-C4p ise "tek çizginin dalgası kuyruğu
# oynatınca yasanın %66-92'si doğuyor" demişti. C4 kanalı KİNEMATİK
# görünüyor (herhangi bir dalga kuyruğu oynatabilir) — o hâlde boyalı
# dalga da aynı R'yi vermeli. VERMEZSE seçicilik kanalın İÇİNDE.
print("\n########## Y6: C4 KANALI — GERÇEK ÇİZGİ vs BOYALI DALGA ##########",
      flush=True)
KUY = KUY0
lq_k, w_k = LOGQ[KUY], W0[KUY]
lq_t, w_t = LOGQ[~KUY], W0[~KUY]
S_tab = S_of(tg, lq_t, w_t)
print(f"{'tür':>5} {'τ':>6} {'A₁':>8} {'R_p':>8} {'R_nn':>8} {'−2cosπτ':>9} "
      f"{'R/hedef':>8} {'kesişme':>8}", flush=True)
for tur, etk, wx, ax in (
        [("çizgi", q, float(np.log(q)), W0[QS.index(q)]) for q in [3, 11, 23]] +
        [("boyalı", tx, tx * Lb, 1.0 / (np.pi * np.exp(tx * Lb / 2)))
         for tx in [0.170, 0.300, 0.380]]):
    lqx = np.array([wx]); wxa = np.array([ax])
    u_x = -S_of(tg, lqx, wxa) / rho
    if tur == "çizgi":
        zc = tg + (-(S_of(tg + u_x, lq_k, w_k) + S_tab) / rho)
    else:
        zc = tg + (-(S_of(tg + u_x, lq_k, w_k) + S_tab
                     + S_of(tg, lqx, wxa)) / rho)
    kes = int((np.diff(zc) <= 0).sum())
    oc, sec, _, Lbc, *_ = olc_ext(np.sort(zc), ekstra=[("x", wx)])
    r = oc[-1] if tur == "boyalı" else \
        [d for d in oc if d["kind"] == "q" and d["lab"] == etk][0]
    hed = -2 * np.cos(np.pi * r["tau"])
    print(f"{tur:>5} {r['tau']:>6.3f} {r['A1']:>8.4f} {r['Rp']:>+8.3f} "
          f"{r['Rn']:>+8.3f} {hed:>+9.3f} {r['Rn']/hed:>7.1%} {kes:>8d}",
          flush=True)

# ---- ÖZET EĞRİ + PNG ------------------------------------------------
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    kfit = kl[is_p & (tl <= 0.35)].mean()
    tt = np.linspace(0.05, 0.55, 400)
    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.plot(tt, -2 * np.cos(np.pi * tt), "k-", lw=2, label=r"$-2\cos\pi\tau$")
    ax.plot(tt, kfit * np.cos(np.pi * tt), "--", color="tab:red", lw=1.6,
            label=rf"$\kappa_{{ölçüm}}\cos\pi\tau$  ($\kappa={kfit:+.2f}$, asal τ≤0.35)")
    ax.plot(tl[is_p], rl[is_p], "o", color="tab:blue", ms=6, label="lab: asal çizgiler")
    ax.plot(tl[~is_p], rl[~is_p], "^", color="tab:cyan", ms=5, label="lab: asal kuvvetleri")
    st = np.array([s[0] for s in sonda]); sr = np.array([s[1] for s in sonda])
    ax.plot(st, sr, "s", color="tab:orange", ms=7, label="lab: çizgi-DIŞI boyalı sondaj")
    ax.set_xlabel(r"$\tau = \log q / \bar L$"); ax.set_ylabel(r"$R_{nn}$")
    ax.axhline(0, color="0.7", lw=0.8)
    ax.set_title("132c — R_nn(τ): çizgiler, çizgi-dışı sondajlar ve çekirdek")
    ax.legend(fontsize=9); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig("132c_egri.png", dpi=130)
    print("\n132c_egri.png yazıldı", flush=True)
except Exception as e:
    print(f"\n(çizim atlandı: {e})", flush=True)
