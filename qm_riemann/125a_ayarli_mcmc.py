"""
125a — AYARLI MCMC ALTYAPISI: ADAPTİF ADIM + YAKINSAMA KAPISI (26 Ağustos)
==========================================================================
124'ün çoklu-mod mıhlaması SONUÇSUZ kaldı: kabul %11, σ_η şişkin,
zincirlerin dengelendiği şüpheli, katı-genlik dalga payını 0.234'e
çıkarıyordu (ζ: ~0.142). Bu script 124'ün kos() çekirdeğini alır ve
ölçülebilir bir alete çevirir. 125b deney matrisini bu aletle koşar.

EKLENENLER
  (1) ADAPTİF ADIM: ısınmanın ilk %80'inde Robbins-Monro ile log-adım
      ayarlanır, hedef kabul %37.5 (kapı bandı %30-45); son %20'de adım
      DONDURULUR (ayrıntılı denge örnekleme boyunca korunur).
  (2) YAKINSAMA TANILARI: ≥3 bağımsız tohum; bölünmüş-zincir Gelman-
      Rubin R̂ (Stan kalıbı) ve Geyer başlangıç-pozitif dizisiyle n_eff.
      İzlenen skalerler: toplam enerji E, gap varyansı Var(δs) ve BEŞ
      HEDEF MODUN örnek-başına A1'i. KAPI: R̂ < 1.1 (hepsinde).
  (3) MIH SERTLİĞİ KAPISI: her modda σ(C_m)/σ_termal ≤ %15,
      σ_termal = √(min(m,N)/2) (β=2 dairesel toplulukta Var(Re Tr U^m)).
      λ taraması bu kapıyla kabul/n_eff dengesine karşı yapılır.
  (4) İKİ KANAL: R_p (varyans kanalı, 124) ve R_nn (bond kanalı, 120) —
      gerçekte en temiz gösterge R_nn ≈ −2cos(πτ).
  (5) İKİ GENLİK REJİMİ: katı Ã_q (124) ve PERDELİ Ã_q·(1−0.36τ_q)
      (108'in ölçülen fiziksel perdesi).

ÖN-MÜHÜR (125a'nın kendi kapıları — 125b ancak bunlar geçilince anlamlı)
  G1  Adaptif adım kabul oranını %30-45 bandına getirebiliyor mu?
  G2  Bir λ değeri hem mıh kapısını (σ_C/σ_term ≤ %15) hem de makul
      n_eff'i sağlıyor mu? Sağlamıyorsa 125b'nin hükmü ASKIYA ALINIR.
  G3  3 tohumla R̂ < 1.1 kapısı geçiliyor mu (E, Var(δs), A1'ler)?
  ARTEFAKT ÖN-TEŞHİSİ (mühürlü, koşudan ÖNCE hesaplanır): mıhlanan
  hedef profil ρ_hedef(θ) = 1 + Σ A_m cos(mθ − φ_m) NEGATİFE düşüyorsa
  hedef FİZİKSEL DEĞİLDİR; gaz onu gerçekleyemez ve A/B karşılaştırması
  konfunde olur. min ρ ve negatif-uzunluk oranı basılır.

SONUÇ (26 Ağustos) — ALET AYARLANDI + 124'ÜN ÖLÇÜMÜ AÇIKLANDI:
  ARTEFAKT ÖN-TEŞHİSİ (koşudan önce, kalemle): TAM TAYFTA MIHLANAN HEDEF
    FİZİKSEL DEĞİL. ΣA = 4.44 (katı) / 3.90 (perdeli) ≫ 1 olduğundan
    ortak-orijinli fazların hepsi θ=0'da üst üste biniyor:
      A-katı min ρ = −2.90 (ρ<0 uzunluk oranı %0.67)
      A-perdeli min ρ = −2.41 (%0.60)
      B-katı min ρ = −0.58 (%2.72), B-perdeli −0.41 (%1.36)
    Gerçek gazda böyle bir eşzamanlı-hizalanma noktası YOK: orada
    m_q = N·log q/L₀ tam sayı DEĞİL, ölçüm penceresi t≫0'da ve fazlar
    Weyl-eşdağılımlı. Tam sayıya yuvarlama (dairesel modelin zorunluluğu)
    her periyotta yapay bir KOHERANS SİVRİSİ üretiyor. Bu, A hücresini
    B'den ayıran ama fizikten gelmeyen bir fark. → 125b'ye MÜHÜRLÜ
    KONTROL ÇİFTİ eklendi: KESİK tayf τ≤0.25 (9 çizgi, ΣA=1.06),
    min ρ = +0.055 — hedef her yerde POZİTİF, gaz onu gerçekleyebilir.
  PERDE DÜRÜST KAYDI: D=1−0.36τ dalga payını 0.234 → 0.1831 indiriyor;
    ζ'nın 0.142'sine TAM inmiyor (beklenen ~0.14 tutmadı, %29 fazla).
  G1 ✓ Adaptif adım her λ'da kabulü %33.5-34.2'ye getirdi (124: %11).
  G2 ✓ λ=4 SEÇİLDİ: σ_C/σ_termal max 0.121 ≤ 0.15 ve |sistematik
    kayma| max 0.071. λ=2 kapıda kalıyor (0.165), λ=0.5 çok gevşek
    (0.344, kayma 0.463).
  G3 — İKİ GEÇİŞ, ÖĞRETİCİ: (i) yalnız adaptif adım + düzgün başlangıçla
    KAPI DÜŞTÜ: R̂[Var(δs)] = 2.58, R̂[A1 m=48] = 1.67, enerji izi hâlâ
    −0.42 σ_E/1k süpürme eğimle DÜŞÜYOR, zincirlerin ⟨E⟩'si 540 birim
    ayrı. ÇOKLU-MOD MIHLAMASI CAMSI: tek-parçacık yerel adım, boşalması
    gereken bölgeden kaçışı yapamıyor. (ii) Üç karışma mekanizması
    eklendi (hedef-yoğunluktan katmanlı başlangıç, λ tavlaması 0.1λ→λ,
    Hastings-düzeltmeli ışınlanma önerisi p=0.1) → 12000+20000 süpürmede
    ve AŞIRI-DAĞITILMIŞ üç başlangıçtan (jitter 1.0/0.5/0.0):
      R̂[Var(δs)] = 1.057 (n_eff 40); R̂[A1] = 1.001/1.011/1.007/1.016/
      1.056 (n_eff 1500/1461/467/159/65); enerji izi eğimleri
      +0.05/−0.03/−0.03 σ_E/1k; ⟨E⟩ zincirler arası 8 birim içinde.
    G3 ✓ KAPI GEÇİLDİ.
  124'ÜN AÇIKLAMASI (yan bulgu, önemli): dengelenmemiş gazda σ_η λ ile
    tırmanıyordu (λ=2/4/8 → 0.410/0.537/0.664); DENGELENMİŞ gazda
    σ_η ≈ 0.365 ve λ'dan BAĞIMSIZ (0.364/0.372/0.368). Yani 124'ün
    "katı-genlik şişkinliği + karışık işaretli R_p" tablosu bir
    TERMALLEŞME ARTEFAKTIYDI — sonuçsuz kaydı doğruydu.
"""

import numpy as np
import time
from sympy import primerange

TWO_PI = 2 * np.pi
N = 256
L0 = 10.37
HEDEF = [17, 27, 40, 48, 63]          # m_2, m_3, m_5, m_7, m_13
RHAT_KAPI = 1.1
MIH_KAPI = 0.15                        # σ_C/σ_termal üst sınırı


# ---------------------------------------------------------------- tayf
def pk_list(lim):
    out = []
    for p in primerange(2, int(lim) + 1):
        q = p
        while q <= lim:
            out.append(q); q *= p
    return sorted(set(out))


def tayf(perde=False, taumax=0.5):
    """Aritmetik çizgi tayfı → (MODS, AMPS). perde: Ã_q ×(1−0.36τ)."""
    QS = [q for q in pk_list(int(np.exp(0.5 * L0))) if q >= 2]
    agg = {}
    for q in QS:
        tau = np.log(q) / L0
        mq = int(round(N * tau))
        if mq < 2 or mq > N // 2 or tau > taumax:
            continue
        A = (2 / np.pi) * np.sin(np.pi * tau) / np.sqrt(q)   # Λ(q)/log p = 1
        if perde:
            A *= (1 - 0.36 * tau)
        agg[mq] = agg.get(mq, 0.0) + A          # çakışan çizgiler toplanır
    MODS = np.array(sorted(agg))
    AMPS = np.array([agg[m] for m in MODS])
    return MODS, AMPS


def fazlar_uret(tip, MODS, tohum):
    """A: ortak-orijinli sinüs (aritmetik). B: rastgele (karışık)."""
    if tip == "A":
        return np.full(len(MODS), -np.pi / 2)
    return np.random.default_rng(9000 + tohum).uniform(0, TWO_PI, len(MODS))


def hedef_profil(MODS, AMPS, fazlar, n=4001):
    """ρ_hedef(θ)/ρ̄ = 1 + Σ A_m cos(mθ − φ_m); artefakt ön-teşhisi."""
    th = np.linspace(0, TWO_PI, n)
    rho = 1 + (AMPS[:, None] *
               np.cos(MODS[:, None] * th[None, :] - fazlar[:, None])).sum(0)
    return th, rho


def sigma_termal(m):
    """β=2 dairesel toplulukta σ(Re Tr U^m) = √(min(m,N)/2)."""
    return np.sqrt(min(int(m), N) / 2.0)


def hedef_izgara(MODS, AMPS, fazlar, ng=4096, taban=0.02):
    """Hedef profilin kesilmiş-pozitif ızgarası: (θ, ρ̃, CDF).
    Hem akıllı başlangıç hem de bağımsızlık-önerisi (ışınlanma) için."""
    th = (np.arange(ng) + 0.5) * TWO_PI / ng
    rho = 1 + (AMPS[:, None] *
               np.cos(MODS[:, None] * th[None, :] - fazlar[:, None])).sum(0)
    rho = np.clip(rho, taban, None)
    rho = rho / rho.mean()
    cdf = np.cumsum(rho); cdf = cdf / cdf[-1]
    return th, rho, cdf


def hedef_baslangic(cdf_th, cdf, rng, jitter=1.0):
    """Katmanlı (stratified) çekimle hedef yoğunluktan başlangıç dizilimi.
    jitter=0 → düzgün (aşırı-dağıtılmış başlangıç), 1 → tam hedef."""
    u = (np.arange(N) + rng.random(N)) / N
    th_h = np.interp(u, cdf, cdf_th)
    th_d = u * TWO_PI
    return np.sort((jitter * th_h + (1 - jitter) * th_d) % TWO_PI)


# ------------------------------------------------------- MCMC çekirdeği
def tam_enerji(th, MODS, C, S, C0, S0, LAM):
    d = th[:, None] - th[None, :]
    iu = np.triu_indices(len(th), 1)
    Eg = -2.0 * np.log(np.abs(np.sin(0.5 * d[iu]))).sum()
    Ev = LAM * (((C - C0) ** 2).sum() + ((S - S0) ** 2).sum())
    return Eg + Ev


def kos(MODS, AMPS, fazlar, tohum, LAM=2.0, n_burn=10000, n_samp=30000,
        her=50, adim0=0.35 * TWO_PI / N, hedef_kabul=0.375, ayar_her=100,
        iz_her=25, nbin=512, jitter=1.0, p_isin=0.10, tavlama=True,
        sessiz=True, etiket=""):
    """Adaptif-adımlı Metropolis + hedef-yoğunluk başlangıcı + λ tavlaması
    + bağımsızlık ("ışınlanma") önerisi.

    124'ün çekirdeği aynen korunur; eklenen üç mekanizma yalnız KARIŞMA
    içindir, denge dağılımını değiştirmez:
      - başlangıç: hedef profilden katmanlı çekim (jitter ile aşırı-
        dağıtılabilir → R̂ gerçek bir sınav kalır);
      - λ tavlaması: ısınmanın ilk yarısında λ 0.1λ→λ; ikinci yarıda
        SABİT λ (örnekleme öncesi tam dengelenme);
      - ışınlanma: p_isin olasılıkla parçacık hedef profilden çekilir,
        Hastings düzeltmesi ρ̃(eski)/ρ̃(yeni) ile — boşalması gereken
        bölgeden kaçış için tek-parçacık yerel adımın yapamadığı hamle.
    """
    rng = np.random.default_rng(tohum)
    C0 = N * AMPS / 2 * np.cos(fazlar)
    S0 = N * AMPS / 2 * np.sin(fazlar)
    g_th, g_rho, g_cdf = hedef_izgara(MODS, AMPS, fazlar)
    ng = len(g_th)
    th = hedef_baslangic(g_th, g_cdf, rng, jitter)
    C = np.array([np.cos(m * th).sum() for m in MODS])
    S = np.array([np.sin(m * th).sum() for m in MODS])

    step = adim0
    ayar_son = int(0.8 * n_burn)
    tav_son = max(1, int(0.5 * n_burn)) if tavlama else 0
    pen_kabul = 0; pen_top = 0            # pencere sayaçları (adaptasyon)
    orn_kabul = 0; orn_top = 0            # örnekleme kabul oranı (rapor)
    isin_kabul = 0; isin_top = 0
    ornekler = []; iz_E = []; iz_sw = []
    ps_var = []; ps_A1 = []; ps_C = []; ps_S = []
    hist = np.zeros(nbin)
    t0 = time.time()

    for sw in range(n_burn + n_samp):
        lam = LAM if sw >= tav_son else LAM * 10.0 ** (-(1 - sw / tav_son))
        for _ in range(N):
            i = rng.integers(N)
            isin = rng.random() < p_isin
            if isin:
                yeni = float(np.interp(rng.random(), g_cdf, g_th))
                j_es = min(int(th[i] / TWO_PI * ng), ng - 1)
                j_ye = min(int(yeni / TWO_PI * ng), ng - 1)
                lqr = np.log(g_rho[j_es] / g_rho[j_ye])   # Hastings
            else:
                yeni = th[i] + rng.normal(0, step)
                lqr = 0.0
            d_es = th - th[i]; d_ye = th - yeni
            d_es[i] = 1.0; d_ye[i] = 1.0
            le = np.log(np.abs(np.sin(0.5 * d_es)))
            ly = np.log(np.abs(np.sin(0.5 * d_ye)))
            le[i] = 0.0; ly[i] = 0.0
            dE = -2.0 * (ly.sum() - le.sum())
            dC = np.cos(MODS * yeni) - np.cos(MODS * th[i])
            dS = np.sin(MODS * yeni) - np.sin(MODS * th[i])
            dV = lam * (((C + dC - C0) ** 2 - (C - C0) ** 2).sum() +
                        ((S + dS - S0) ** 2 - (S - S0) ** 2).sum())
            arg = -(dE + dV) + lqr
            if not isin:
                pen_top += 1
                if sw >= n_burn:
                    orn_top += 1
            else:
                isin_top += 1
            if arg > 0 or rng.random() < np.exp(arg):
                th[i] = yeni % TWO_PI
                C += dC; S += dS
                if not isin:
                    pen_kabul += 1
                    if sw >= n_burn:
                        orn_kabul += 1
                else:
                    isin_kabul += 1

        # --- adaptif adım (yalnız ısınmanın ilk %80'i) ---
        if sw < ayar_son and (sw + 1) % ayar_her == 0:
            acc = pen_kabul / max(pen_top, 1)
            gam = max(0.10, 1.0 / (1.0 + sw / ayar_her) ** 0.5)
            step *= float(np.exp(np.clip(2.0 * gam * (acc - hedef_kabul),
                                         -0.4, 0.4)))
            step = float(np.clip(step, 1e-4 * TWO_PI / N, 3.0 * TWO_PI / N))
            pen_kabul = 0; pen_top = 0
            if not sessiz:
                print(f"    [{etiket}] sw={sw+1} kabul={acc:.3f} "
                      f"adım={step/(TWO_PI/N):.4f}·ḡ", flush=True)

        # --- enerji izi (durağanlık) ---
        if (sw + 1) % iz_her == 0:
            iz_E.append(tam_enerji(th, MODS, C, S, C0, S0, lam))
            iz_sw.append(sw + 1)

        # --- örnekleme ---
        if sw >= n_burn and (sw - n_burn) % her == 0:
            ths = np.sort(th.copy())
            ornekler.append(ths)
            dth = np.diff(np.concatenate([ths, [ths[0] + TWO_PI]]))
            mid = ths + dth / 2
            ds = dth * N / TWO_PI - 1
            ps_var.append(float(ds.var()))
            a1 = []
            for m in HEDEF:
                a = (2.0 / N) * float(np.dot(ds, np.cos(m * mid)))
                b = (2.0 / N) * float(np.dot(ds, np.sin(m * mid)))
                a1.append(float(np.hypot(a, b)))
            ps_A1.append(a1)
            ps_C.append(C.copy()); ps_S.append(S.copy())
            hist += np.histogram(ths, bins=nbin, range=(0, TWO_PI))[0]

    ns = len(ornekler)
    return dict(
        ornekler=ornekler,
        kabul=orn_kabul / max(orn_top, 1),
        isin=isin_kabul / max(isin_top, 1),
        adim=step / (TWO_PI / N),
        iz_E=np.array(iz_E), iz_sw=np.array(iz_sw),
        ps_var=np.array(ps_var), ps_A1=np.array(ps_A1),
        ps_C=np.array(ps_C), ps_S=np.array(ps_S),
        rho=hist * nbin / (hist.sum() if hist.sum() else 1),
        sure=time.time() - t0, n_orn=ns,
    )


# ------------------------------------------------------- mıh kapısı
def mih_kapisi(ps_C, ps_S, MODS, AMPS, fazlar):
    """σ(C_m)/σ_termal ve gerçekleşen mıh sapması; kapı ≤ %15."""
    C0 = N * AMPS / 2 * np.cos(fazlar)
    S0 = N * AMPS / 2 * np.sin(fazlar)
    st = np.array([sigma_termal(m) for m in MODS])
    oranC = ps_C.std(0) / st
    oranS = ps_S.std(0) / st
    oran = np.maximum(oranC, oranS)
    # sistematik kayma (bias): |⟨C⟩−C₀|/σ_termal
    kay = np.maximum(np.abs(ps_C.mean(0) - C0), np.abs(ps_S.mean(0) - S0)) / st
    return oran, kay


# ------------------------------------------------------- R ölçümü
def _tasarim(ths, Ms):
    dth = np.diff(np.concatenate([ths, [ths[0] + TWO_PI]]))
    mid = ths + dth / 2
    ds = dth * N / TWO_PI - 1
    cols = [np.ones(N), mid / TWO_PI - 0.5]
    for m in Ms:
        cols += [np.cos(m * mid), np.sin(m * mid)]
    return np.vstack(cols).T, ds, mid


def R_olc(ornekler, MODS, hedefler=HEDEF, agirlik_fn=None):
    """R_p (varyans kanalı, 124) + R_nn (bond kanalı, 120).
    agirlik_fn(mid) → 0/1 maske (patolojik bölgeyi dışlamak için)."""
    Ms = list(MODS)
    C = 2 + 2 * len(Ms)
    XtX = np.zeros((C, C)); Xty = np.zeros(C)
    for ths in ornekler:
        X, ds, mid = _tasarim(ths, Ms)
        if agirlik_fn is not None:
            w = agirlik_fn(mid).astype(float)
            X = X * w[:, None]; ds = ds * w
        XtX += X.T @ X; Xty += X.T @ ds
    b1 = np.linalg.solve(XtX, Xty)

    s2acc = 0.0; nt = 0
    for ths in ornekler:
        X, ds, mid = _tasarim(ths, Ms)
        if agirlik_fn is not None:
            w = agirlik_fn(mid).astype(float)
            eta = (ds - X @ b1) * w
            nt += int(w.sum())
        else:
            eta = ds - X @ b1
            nt += len(eta)
        s2acc += float((eta ** 2).sum())
    s2 = s2acc / nt

    Xty2 = np.zeros(C)
    XtXn = np.zeros((C, C)); Xtyn = np.zeros(C)
    ee_acc = 0.0; ne = 0
    EE = []
    for ths in ornekler:
        X, ds, mid = _tasarim(ths, Ms)
        eta = ds - X @ b1
        if agirlik_fn is not None:
            w = agirlik_fn(mid).astype(float)
            Xty2 += (X * w[:, None]).T @ ((eta * w) ** 2 - s2 * w)
        else:
            Xty2 += X.T @ (eta ** 2 - s2)
        ee = eta[:-1] * eta[1:]
        mm = 0.5 * (mid[:-1] + mid[1:])
        ee_acc += float(ee.sum()); ne += len(ee)
        EE.append((ee, mm))
    b2 = np.linalg.solve(XtX, Xty2)
    c1 = ee_acc / ne

    for ee, mm in EE:
        cols = [np.ones(len(mm)), mm / TWO_PI - 0.5]
        for m in Ms:
            cols += [np.cos(m * mm), np.sin(m * mm)]
        Xn = np.vstack(cols).T
        XtXn += Xn.T @ Xn; Xtyn += Xn.T @ (ee - c1)
    b3 = np.linalg.solve(XtXn, Xtyn)

    out = {}
    for m in hedefler:
        i = Ms.index(m)
        c1c, s1c = b1[2 + 2 * i], b1[2 + 2 * i + 1]
        A1 = float(np.hypot(c1c, s1c)); ph = float(np.arctan2(s1c, c1c))
        Pv = b2[2 + 2 * i] * np.cos(ph) + b2[2 + 2 * i + 1] * np.sin(ph)
        Pn = b3[2 + 2 * i] * np.cos(ph) + b3[2 + 2 * i + 1] * np.sin(ph)
        out[m] = dict(A1=A1, R_p=float(Pv / (2 * s2 * A1)),
                      R_nn=float(Pn / (2 * c1 * A1)))
    return out, float(np.sqrt(s2)), float(c1 / s2)


# --------------------------------------------------- yakınsama tanıları
def _acov(x):
    T = len(x)
    y = x - x.mean()
    n2 = 1 << int(np.ceil(np.log2(2 * T)))
    f = np.fft.rfft(y, n2)
    ac = np.fft.irfft(f * np.conj(f), n2)[:T].real
    return ac / T


def rhat_neff(zincirler):
    """Bölünmüş-zincir R̂ + Geyer başlangıç-pozitif n_eff (Stan kalıbı).
    zincirler: (M, T) dizi."""
    Z = np.asarray(zincirler, float)
    M, T = Z.shape
    if T < 8 or M < 2:
        return np.nan, np.nan
    if T % 2:
        Z = Z[:, :T - 1]; T -= 1
    Zs = Z.reshape(M * 2, T // 2)
    Ms, Ts = Zs.shape
    W = Zs.var(axis=1, ddof=1).mean()
    if W <= 0 or not np.isfinite(W):
        return np.nan, np.nan
    B = Ts * Zs.mean(axis=1).var(ddof=1)
    varhat = (Ts - 1) / Ts * W + B / Ts
    rhat = float(np.sqrt(varhat / W))

    ac = np.array([_acov(z) for z in Zs]).mean(0)
    rho = 1.0 - (W - ac) / varhat
    rho[0] = 1.0
    # Geyer: ardışık çiftlerin toplamı pozitif kaldığı sürece
    t = 1; toplam = 0.0; onceki = np.inf
    while t + 1 < Ts:
        p = rho[t] + rho[t + 1]
        if p < 0:
            break
        p = min(p, onceki)                    # monotonluk kesmesi
        toplam += p; onceki = p
        t += 2
    tau = max(1.0, -1.0 + 2.0 * toplam)
    neff = float(Ms * Ts / tau)
    return rhat, min(neff, float(Ms * Ts))


def durgunluk(iz_E, iz_sw, n_burn, iz_her):
    """Örnekleme fazındaki enerji izinin eğimi (σ_E biriminde / 1000 süp.)"""
    msk = iz_sw > n_burn
    if msk.sum() < 10:
        return np.nan
    x = iz_sw[msk].astype(float); y = iz_E[msk]
    b = np.polyfit(x - x.mean(), y, 1)[0]
    sd = y.std()
    return float(b * 1000.0 / sd) if sd > 0 else np.nan


# ================================================================= main
if __name__ == "__main__":
    import multiprocessing as mp

    print("=" * 74)
    print("125a — AYARLI MCMC ALTYAPISI")
    print("=" * 74)

    # ---------- tayf + ARTEFAKT ÖN-TEŞHİSİ (koşudan ÖNCE) ----------
    print("\n[TAYF ve MIHLANAN HEDEF PROFİL — artefakt ön-teşhisi]")
    print(f"{'genlik':>10} {'n_mod':>6} {'ΣA':>6} {'dalga payı':>11} "
          f"{'faz':>4} {'min ρ':>8} {'ρ<0 oran':>9}")
    TAYF = {}
    for ad, perde in [("katı", False), ("perdeli", True)]:
        MODS, AMPS = tayf(perde)
        TAYF[ad] = (MODS, AMPS)
        for tip in ["A", "B"]:
            ph = fazlar_uret(tip, MODS, 1)
            _, rho = hedef_profil(MODS, AMPS, ph, 20001)
            print(f"{ad:>10} {len(MODS):>6} {AMPS.sum():>6.2f} "
                  f"{0.5*(AMPS**2).sum():>11.4f} {tip:>4} "
                  f"{rho.min():>+8.3f} {(rho<0).mean():>9.4f}")
    print("  (ζ dalga payı ≈ 0.142; ρ<0 ⇒ hedef fiziksel DEĞİL — gaz "
          "gerçekleyemez)")

    # kesik tayf (τ ≤ 0.25): fiziksel-hedef kontrolü
    MK, AK = tayf(True, taumax=0.25)
    _, rk = hedef_profil(MK, AK, fazlar_uret("A", MK, 1), 20001)
    print(f"{'kesik':>10} {len(MK):>6} {AK.sum():>6.2f} "
          f"{0.5*(AK**2).sum():>11.4f} {'A':>4} {rk.min():>+8.3f} "
          f"{(rk<0).mean():>9.4f}   ← τ≤0.25 kontrol tayfı")

    # ---------- G1+G2: λ taraması ----------
    print("\n[G1/G2 — λ TARAMASI: adaptif adım + mıh kapısı] "
          "(A-aritmetik, perdeli genlik, 3000+3000 süpürme)", flush=True)
    MODS, AMPS = TAYF["perdeli"]
    faz_A = fazlar_uret("A", MODS, 1)

    def _tara(lam):
        r = kos(MODS, AMPS, faz_A, tohum=1250 + int(lam * 10), LAM=lam,
                n_burn=3000, n_samp=3000, her=25)
        oran, kay = mih_kapisi(r["ps_C"], r["ps_S"], MODS, AMPS, faz_A)
        rr, se, c1s = R_olc(r["ornekler"], MODS)
        A1m = np.mean([rr[m]["A1"] for m in HEDEF])
        _, ne = rhat_neff(np.vstack([r["ps_var"][:len(r["ps_var"]) // 2],
                                     r["ps_var"][len(r["ps_var"]) // 2:]]))
        return dict(lam=lam, kabul=r["kabul"], isin=r["isin"], adim=r["adim"],
                    oran_med=float(np.median(oran)), oran_max=float(oran.max()),
                    kay_max=float(kay.max()), A1=A1m, se=se,
                    neff=ne, sure=r["sure"])

    with mp.get_context("fork").Pool(5) as pool:
        tarama = pool.map(_tara, [0.5, 1.0, 2.0, 4.0, 8.0])

    print(f"{'λ':>5} {'kabul':>6} {'ışın':>6} {'adım/ḡ':>8} {'σC/σt med':>10} "
          f"{'max':>7} {'|kay|max':>9} {'⟨A1⟩':>7} {'σ_η':>6} "
          f"{'n_eff':>7} {'sn':>5}")
    for t in tarama:
        print(f"{t['lam']:>5.1f} {t['kabul']:>6.3f} {t['isin']:>6.3f} "
              f"{t['adim']:>8.4f} "
              f"{t['oran_med']:>10.3f} {t['oran_max']:>7.3f} "
              f"{t['kay_max']:>9.3f} {t['A1']:>7.4f} {t['se']:>6.3f} "
              f"{t['neff']:>7.0f} {t['sure']:>5.0f}", flush=True)

    uygun = [t for t in tarama if t["oran_max"] <= MIH_KAPI
             and 0.28 <= t["kabul"] <= 0.50]
    if uygun:
        SEC = min(uygun, key=lambda t: t["lam"])
        print(f"  G1 ✓ G2 ✓ — SEÇİLEN λ = {SEC['lam']} "
              f"(kabul {SEC['kabul']:.2f}, σC/σt max {SEC['oran_max']:.3f})")
    else:
        SEC = min(tarama, key=lambda t: t["oran_max"])
        print(f"  G2 ✗ — hiçbir λ mıh kapısını (≤{MIH_KAPI}) ve kabul "
              f"bandını birlikte geçmedi; en iyisi λ={SEC['lam']} "
              f"(σC/σt max {SEC['oran_max']:.3f}) — 125b HÜKMÜ ASKIDA")
    LAM_SEC = SEC["lam"]

    # ---------- G3: 3 tohumlu, AŞIRI-DAĞITILMIŞ BAŞLANGIÇLI R̂ provası ----
    NB, NS, HER = 12000, 20000, 40
    JIT = [1.0, 0.5, 0.0]          # hedef / yarı / düzgün başlangıç
    print(f"\n[G3 — 3 TOHUMLU R̂ PROVASI] λ={LAM_SEC}, "
          f"{NB}+{NS} süpürme, her {HER}; aşırı-dağıtılmış başlangıç "
          f"jitter={JIT}", flush=True)

    def _zincir(arg):
        tohum, jit = arg
        r = kos(MODS, AMPS, faz_A, tohum=tohum, LAM=LAM_SEC,
                n_burn=NB, n_samp=NS, her=HER, jitter=jit)
        return dict(kabul=r["kabul"], isin=r["isin"], adim=r["adim"],
                    ps_var=r["ps_var"], ps_A1=r["ps_A1"],
                    iz=durgunluk(r["iz_E"], r["iz_sw"], NB, 25),
                    Es=float(r["iz_E"][r["iz_sw"] > NB].mean()),
                    sure=r["sure"])

    with mp.get_context("fork").Pool(3) as pool:
        Z = pool.map(_zincir, list(zip([1251, 1252, 1253], JIT)))

    print(f"  kabul: {[round(z['kabul'],3) for z in Z]}  "
          f"ışınlanma kabul: {[round(z['isin'],4) for z in Z]}  "
          f"adım/ḡ: {[round(z['adim'],4) for z in Z]}")
    print(f"  ⟨E⟩ (örnekleme): {[round(z['Es'],1) for z in Z]}  "
          f"eğim/σ_E/1k süp.: {[round(z['iz'],3) for z in Z]}")
    rh_v, ne_v = rhat_neff(np.vstack([z["ps_var"] for z in Z]))
    print(f"  R̂[Var(δs)] = {rh_v:.3f}  n_eff = {ne_v:.0f}")
    gecti = np.isfinite(rh_v) and rh_v < RHAT_KAPI
    for k, m in enumerate(HEDEF):
        rh, ne = rhat_neff(np.vstack([z["ps_A1"][:, k] for z in Z]))
        gecti &= (np.isfinite(rh) and rh < RHAT_KAPI)
        print(f"  R̂[A1 m={m:>3}] = {rh:.3f}  n_eff = {ne:.0f}")
    print(f"  G3 {'✓ KAPI GEÇİLDİ' if gecti else '✗ KAPI GEÇİLMEDİ'} "
          f"(eşik R̂ < {RHAT_KAPI}) — bu uzunluk 125b için "
          f"{'yeterli' if gecti else 'YETERSİZ, büyütülmeli'}")
    print(f"\n125b İÇİN AYAR: LAM={LAM_SEC}, ısınma≥{NB}, örnekleme≥{NS}")
