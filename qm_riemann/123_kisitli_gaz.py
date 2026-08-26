"""
123 — KISITLI GAZ: SABİTLENMİŞ-MOD LOG-GAZININ NEFESİ (26 Ağustos)
==========================================================================
Kısıt-hipotezi kaleminin belirleyici deneyi. İki mini-no-go (Gaussian
koşullama faz-kilitli nefes üretemez; lineer-yanıt rejimi R_eq=+1'i
kopyalar [116]) → −2 ancak AŞIRI-TERMAL sabitleme rejiminde doğabilir.
DENEY: dairesel β=2 log-gazı (CUE ağırlığı) + m-modunun İKİ kuadratürü
sert cezayla mıhlanır (Re ρ_m = C₀, Im ρ_m = 0); Metropolis örneklemesi;
standart zincirle R_p ölçümü; sabitleme gücü taranır.

ÖN-MÜHÜR:
  K0  C₀=0 (kısıtsız): R ≈ +1 civarı (116'nın öz-dalga dengesiyle uyum).
  KARAR EĞRİSİ: C₀/σ_termal = ~2 / ~5 / ~10 arttıkça R_p +1'den −'ye
  KAYARSA → −2 = kısıtlı log-gazın evrensel eyer-dalgalanma fiziği
  (aritmetiğin rolü: kısıtı sağlamak — büyük birleşme).
  KAYMAZSA (+1 kalırsa) → aritmetik, kısıtlı gazdan da derin.
N=256, m=26 (τ≈0.10); ceza λ[(C−C₀)²+S²], λ=2 (C-dalgalanması ~0.5 ≪
σ_termal≈√(m/2)≈3.6 → etkin mıh). Süpürme: 1500 ısınma + 3000 örnekleme
(20'de bir örnek). Gerçekleşen dalga genliği A1 raporlanır.

SONUÇ (26 Ağustos) — H-KISIT RET: SABİTLENMİŞ GAZ İNATLA ADYABATİK:
  C₀/σ = 1.9 / 5 / 10 → A1 = 0.053 / 0.138 / 0.289; R_p = +1.07 /
  +0.96 / +0.95. On kat termal sabitlemede ve %29'luk dev dalgada bile
  R = +1'de çakılı — −2'ye kayış YOK. (C₀=0 satırı A1≈0 → R tanımsız,
  beklenen.) HÜKÜM: −2, dengede-gaz + mod-mıhlama fiziğinden ÇIKMAZ —
  büyük-sapma eyeri bile adyabatik. DÖRDÜNCÜ DARALMA: anomali
  (1) evrensel, (2) denge-dışı, (3) sektör-seçici, (4) sabitlenmiş-
  gazla üretilemez. KALAN ADAYLAR: (a) ÇOKLU-MOD sabitleme — gerçek
  gaz ~60 çizgiyi p^k-kule yapısıyla BİRDEN taşır; modlar-arası yapı
  tek-mod deneyinde yok (sıradaki basamak adayı); (b) determinizm/
  spektral yapı — ζ örgüsü bir Gibbs/koşullu ölçü değil, deterministik
  bir tayfın ergodik izi; −2 iz-formülünün 2. mertebe yapısı olabilir.
"""

import numpy as np
import time

TWO_PI = 2 * np.pi
rng = np.random.default_rng(123)
N = 256
m = 26
gbar = TWO_PI / N
LAM = 2.0
STEP = 0.35 * gbar

def enerji_fark(th, i, yeni):
    d_es = th - th[i]
    d_ye = th - yeni
    d_es[i] = 1.0; d_ye[i] = 1.0
    le = np.log(np.abs(np.sin(0.5 * d_es)))
    ly = np.log(np.abs(np.sin(0.5 * d_ye)))
    le[i] = 0.0; ly[i] = 0.0
    return -2.0 * (ly.sum() - le.sum())

def kos(C0, n_burn=1500, n_samp=3000, her=20):
    th = np.sort(rng.uniform(0, TWO_PI, N))
    C = np.cos(m * th).sum(); S = np.sin(m * th).sum()
    kabul = 0; toplam = 0
    ornekler = []
    for sw in range(n_burn + n_samp):
        for _ in range(N):
            i = rng.integers(N)
            yeni = th[i] + rng.normal(0, STEP)
            dE = enerji_fark(th, i, yeni)
            dC = np.cos(m * yeni) - np.cos(m * th[i])
            dS = np.sin(m * yeni) - np.sin(m * th[i])
            dV = LAM * (((C + dC - C0)**2 - (C - C0)**2) +
                        ((S + dS)**2 - S**2))
            toplam += 1
            if dE + dV < 0 or rng.random() < np.exp(-(dE + dV)):
                th[i] = yeni % TWO_PI
                C += dC; S += dS
                kabul += 1
        if sw >= n_burn and (sw - n_burn) % her == 0:
            ornekler.append(np.sort(th.copy()))
    return ornekler, kabul / toplam

def R_olc(ornekler):
    XtX = np.zeros((5, 5)); Xty1 = np.zeros(5)
    DATA = []
    for th in ornekler:
        dth = np.diff(np.concatenate([th, [th[0] + TWO_PI]]))
        mid = th + dth / 2
        ds = dth * N / TWO_PI - 1
        X = np.vstack([np.ones(N), np.cos(m * mid), np.sin(m * mid),
                       np.cos(2 * m * mid), np.sin(2 * m * mid)]).T
        XtX += X.T @ X; Xty1 += X.T @ ds
        DATA.append((X, ds))
    b1 = np.linalg.solve(XtX, Xty1)
    A1 = np.hypot(b1[1], b1[2]); ph = np.arctan2(b1[2], b1[1])
    s2acc = 0.0; n_t = 0
    Xty2 = np.zeros(5)
    ETA = []
    for X, ds in DATA:
        eta = ds - X @ b1
        ETA.append((X, eta))
        s2acc += (eta**2).sum(); n_t += len(eta)
    s2 = s2acc / n_t
    for X, eta in ETA:
        Xty2 += X.T @ (eta**2 - s2)
    b2 = np.linalg.solve(XtX, Xty2)
    Rp = (b2[1] * np.cos(ph) + b2[2] * np.sin(ph)) / (2 * s2 * A1)
    return A1, Rp, np.sqrt(s2)

sig_term = np.sqrt(m / 2.0)
print(f"N={N}, m={m} (τ={m/N:.3f}); σ_termal(kuadratür) ≈ {sig_term:.2f}")
print(f"{'C₀':>5} {'C₀/σ':>6} {'A1':>7} {'R_p':>7} {'σ_η':>6} {'kabul':>6}")
for C0 in [0.0, 7.0, 18.0, 36.0]:
    t0 = time.time()
    orn, acc = kos(C0)
    A1, Rp, se = R_olc(orn)
    print(f"{C0:>5.0f} {C0/sig_term:>6.1f} {A1:>7.4f} {Rp:>+7.3f} "
          f"{se:>6.3f} {acc:>6.2f}  [{time.time()-t0:.0f} sn]", flush=True)
print("\nKARAR: C₀/σ arttıkça R_p +1'den −'ye kayıyor mu? "
      "(gerçek gaz: −2; denge: +1)")
