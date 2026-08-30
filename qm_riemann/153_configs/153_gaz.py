"""
153 — İTMELİ GAZ  (152_gaz.py'nin parametrik kopyası)
=====================================================
152'nin ölçüm zinciri BİREBİR korunmuştur (satır satır aynı kod yolu:
`g = np.diff(z)`den itibaren tek karakter değişmedi). Eklenen tek şey,
o zincire GİREN z'nin nasıl üretildiği — ve zincirin sonunda raporlanan
küçük-s göstergesi (ölçümü etkilemez, yalnız yazdırılır).

SORU: saf-merdiven gazının anormal dispersiyon fazı gerçeğin ~1.4 katı
ve 152'nin sekiz konfigürasyonunda KATI kaldı. Eksik aday malzeme:
GUE SEVİYE İTMESİ — birinci-mertebe merdiven gazında yok. İki yoldan
sınanıyor:

  İ1 (tip "itme")  — YEREL İTME GEVŞETMESİ
      erfc-0.68 gazının Newton çözümünden SONRA, kısa-menzilli itme
      süpürmeleri:   z_n ← z_n + ε·ḡ_n·[ f(g_{n−1}) − f(g_n) ],
      f(s) = (ḡ/s)²,  g_n = z_{n+1} − z_n,  yalnız komşu çift.
      Görev notundaki f(s)=c·ḡ²/s² ile aynı; c dimensiyonsuz olduğu için
      ε'nun uzunluk taşıması gerekir — burada ε ḡ birimindedir (ε·ḡ
      çarpanı) ve c ≡ 1 ε'ya emilmiştir. TEK KADRAN: ε (+ süpürme sayısı).
      Yer değiştirme ±λ·min(g_{n−1},g_n) ile kelepçelenir (λ=0.4): bu
      kelepçe SIRALAMANIN BOZULMAMASINI GARANTİ EDER (komşu iki nokta en
      fazla 0.4g kadar birbirine yaklaşabildiğinden yeni aralık ≥ 0.2g).
      Kelepçenin ne sıklıkta bağladığı ve pencere yoğunluğunun ne kadar
      kaydığı her süpürmede raporlanır — z YENİDEN ÖLÇEKLENMEZ/ÖTELENMEZ
      (hareket zaten teleskopik olduğu için yoğunluk korunumu kendiliğinden
      gelir; kelepçe onu bozarsa sayıyla görülsün diye).
      ε'nun İKİ İŞARETİ de koşulur. Sebep ölçümdür, tercih değil: taban
      gazında P(s<0.3)=0.00013 iken gerçekte 0.02420 — sentetik gaz kısa
      menzilde gerçeğin ~190 katı KATI. Görevin kalibrasyon hedefi (küçük-s
      kuyruğunu gerçeğe yaklaştırmak) bu yüzden ε<0 (kısa-menzilli
      yumuşatma) ister; ε>0 (yazıldığı gibi itme) hedeften uzaklaştırır.
      İkisi de tabloda, ikisi de etiketli.

  İ2 (tip "boya")  — İTMELİ-TABAN + BOYALI MERDİVEN
      1) bağımsız gap'ler s_k (üç taban seçeneği, aşağıda) → açılmış
         konumlar u_k = n0 + Σ s_j
      2) pürüzsüz sayımın tersi:  rvm_N(x_k) = u_k  (merdiven YOK)
      3) merdiven yer değiştirmesi TEK GEÇİŞ boyanır (öz-tutarlılık YOK):
         z = x − S(x)/N'(x) = x + ḡ(x)·Σ a_q sin(ω_q x)
         (S(z) = −Σ a_q sin(ω_q z) olduğundan iki yazım özdeştir; kodda
          152'nin S fonksiyonu kullanılır ki merdiven bit-bit aynı olsun.)
      taban seçenekleri:
         gue      : P(s) = (32/π²)s²e^{−4s²/π} (GUE Wigner surmise, ort=1)
                    ters-CDF örneklemesi; CDF(x)=erf(2x/√π)−(4/π)xe^{−4x²/π}
         poisson  : s = −ln U  (itme YOK) — İ2 yolunda itmenin payını
                    ayırmak için kontrol
         det      : s ≡ 1  (ne itme ne düzensizlik) — "tek geçiş mi,
                    öz-tutarlılık mı?" sorusunu ayıran kontrol

Merdiven her koşuda erfc(τ_c=0.68, Δ=0.125) kesimlidir — 152'nin S3'ü
gerçeğe oturtan TABAN konfigürasyonu.

Kullanım:  python 153_gaz.py <konfig>
"""

import sys
import time
import math
import json
import numpy as np
from pathlib import Path
import multiprocessing as mp

HERE = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCRATCH = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
               "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/153")
TWO_PI = 2 * np.pi
NEWTON_IT = 20          # 152 ile aynı
NWORK = 6               # S(z) havuzu
NPT = 25000             # süreç başına nokta bloğu
TAU_C, DELTA = 0.68, 0.125      # TABAN kesim (152'nin A4'ü)
LAMBDA_KLP = 0.4                # itme kelepçesi: ±λ·min(sol,sağ) aralık

# konfig: (tip, parametreler)
KONFIG = {
    # --- kontrol: 152'nin A4'ü; z diske yazılır, İ1 bunu okur ---
    "taban":  ("erfc", {}),
    # --- İ1a: ε>0, görevde YAZILDIĞI GİBİ itme (komşuları uzaklaştırır) ---
    "R1":     ("itme", {"eps": 0.003, "nsup": 5}),
    "R2":     ("itme", {"eps": 0.010, "nsup": 5}),
    "R3":     ("itme", {"eps": 0.030, "nsup": 5}),
    # --- İ1b: ε<0, KALİBRASYON HEDEFİNİN istediği yön ---
    # Ön ölçüm (kalib.py): taban gazında P(s<0.3)=0.00013, gerçek 0.02420 —
    # yani sentetik gaz kısa menzilde zaten GEREĞİNDEN İTİCİ. ε>0 hedefi
    # daha da uzaklaştırıyor (P(s<0.3) → 0). Kuyruk oranını gerçeğe
    # YAKLAŞTIRAN tek yön ε<0'dır (kısa-menzilli yumuşatma). ε aşağıdaki
    # üç değerde P(s<0.3)≈0.024'e kalibre edildi; üçü farklı süpürme
    # sayısında aynı hedefi tutturur (sonuç süpürme sayısına bağlı mı?).
    "N3":     ("itme", {"eps": -0.00700, "nsup": 3}),
    "N5":     ("itme", {"eps": -0.00370, "nsup": 5}),
    "N10":    ("itme", {"eps": -0.00180, "nsup": 10}),
    # kasıtlı aşırı doz merdiveni: faz doza duyarlı mı, yoksa hiç mi
    # kıpırdamıyor? (N5x'te faz oranı belirgin düşünce doz artırıldı.)
    "N5x":    ("itme", {"eps": -0.01000, "nsup": 5}),
    "N5y":    ("itme", {"eps": -0.02000, "nsup": 5}),
    "N5z":    ("itme", {"eps": -0.03000, "nsup": 5}),
    # --- J: AYIRT EDİCİ KONTROL ---
    # N-serisi fazı düşürüyor ama σ_ds²'yi de büyütüyor. Faz hangisine
    # bağlı: kısa-menzil YAPISINA mı, yoksa yalnız σ_ds²'ye mi? J'ler
    # yapısız (inkoherent, bağımsız) titreşimle AYNI σ_ds²'yi üretir:
    #   σ_ds² ≈ 0.1128 + 2σ_j²  → σ_j=0.141 ≈ N5,  σ_j=0.260 ≈ N5x
    # Faz oranları N5/N5x ile aynı çıkarsa faz yalnız σ_ds²'yi görüyordur.
    "J14":    ("titresim", {"sigma_j": 0.141}),
    "J26":    ("titresim", {"sigma_j": 0.260}),
    # --- İ2: itmeli taban + boyalı merdiven ---
    "P1":     ("boya", {"taban": "gue"}),
    "P0":     ("boya", {"taban": "poisson"}),
    "Pd":     ("boya", {"taban": "det"}),
}

GERCEK = {0.5375: (0.787, 0.206, 0.40), 0.585: (0.550, 0.488, 1.07),
          0.660: (0.118, 0.614, 1.85), 0.740: (-0.148, 0.498, 2.12),
          0.815: (-0.103, 0.655, 1.50)}
BANTLAR = [(0.525, 0.55), (0.55, 0.62), (0.62, 0.70),
           (0.70, 0.78), (0.78, 0.85)]
# gerçek zeros6 son-300k'dan ÖLÇÜLDÜ (aynı açılma: s = g·log(mid/2π)/2π)
GERCEK_PS = {0.1: 0.000960, 0.2: 0.007327, 0.3: 0.024203,
             0.5: 0.103847, 0.7: 0.247217}


def rvm_N(t):
    x = t / TWO_PI
    return x * np.log(x / np.e) + 7 / 8


def rvm_d(t):
    return np.log(t / TWO_PI) / TWO_PI


def pk(lim):
    from sympy import primerange
    out = []
    for p in primerange(2, lim + 1):
        q = p
        while q <= lim:
            out.append(q); q *= p
    return sorted(set(out))


# ---------- S(z) ve S'(z): 151b/152 ile aynı formül, bloklu ----------
def S_ve_Sp(z, om, a, blok=800, npt=40000):
    S = np.zeros_like(z); Sp = np.zeros_like(z)
    for b0 in range(0, len(om), blok):
        w = om[b0:b0 + blok]; aa = a[b0:b0 + blok]
        for s0 in range(0, len(z), npt):
            sl = slice(s0, min(s0 + npt, len(z)))
            arg = np.outer(z[sl], w)
            S[sl] += -np.sin(arg) @ aa
            Sp[sl] += -np.cos(arg) @ (aa * w)
            del arg
    return S, Sp


_G = {}


def _init(om, a):
    _G["om"] = om; _G["a"] = a


def _work(zc):
    return S_ve_Sp(zc, _G["om"], _G["a"], npt=NPT)


# ---------- İ1: kısa-menzilli itme süpürmesi ----------
def itme_supur(z, eps, nsup, log=print):
    """z_n ← z_n + ε·ḡ_n·[f(g_{n−1}) − f(g_n)],  f(s)=(ḡ/s)².

    Uçlarda eksik aralık ḡ ile doldurulur (f=1, nötr). Yer değiştirme
    ±λ·min(sol,sağ) ile kelepçelenir → hiçbir aralık sıfıra inemez.
    z yeniden ölçeklenmez/ötelenmez.
    """
    gb = 1.0 / rvm_d(z)                 # yerel ortalama aralık, nokta başına
    W = np.array_split(np.arange(len(z)), 30)

    def pencere_g(zz):
        return np.array([(zz[w[-1]] - zz[w[0]]) / (len(w) - 1) for w in W])

    pen0 = pencere_g(z)                 # süpürmelerden ÖNCEKİ yoğunluk
    for it in range(nsup):
        gp = np.concatenate(([gb[0]], np.diff(z), [gb[-1]]))  # padli aralıklar
        f = (gb / gp[:-1]) ** 2 - (gb / gp[1:]) ** 2   # f(g_{n−1}) − f(g_n)
        ham = eps * gb * f
        lim = LAMBDA_KLP * np.minimum(gp[:-1], gp[1:])
        d = np.clip(ham, -lim, lim)
        nklp = int(np.sum(np.abs(ham) > lim + 1e-18))
        z = z + d
        capraz = int(np.sum(np.diff(z) <= 0))
        kay = float(np.max(np.abs(pencere_g(z) - pen0)) / gb.mean())
        log(f"  itme süpürme {it+1}/{nsup}: |Δz|ort={np.mean(np.abs(d)):.5f} "
            f"maks={np.max(np.abs(d)):.5f} (ḡ={gb.mean():.4f}), "
            f"kelepçe {nklp} ({100.0*nklp/len(z):.3f}%), "
            f"sıra bozan {capraz}, pencere-ḡ kayması {kay:.2e}·ḡ", flush=True)
        if capraz:
            raise SystemExit("HATA: itme sıralamayı bozdu — kelepçe garantisi "
                             "çalışmadı, koşu düşürüldü.")
    return z


# ---------- İ2: bağımsız gap örnekleyicileri ----------
def gap_ornekle(taban, n, rng):
    if taban == "det":
        return np.ones(n)
    if taban == "poisson":
        return -np.log(rng.random(n))
    if taban == "gue":
        # ters-CDF: CDF(x) = erf(2x/√π) − (4/π)·x·e^{−4x²/π}
        from scipy.special import erf
        xs = np.linspace(0.0, 6.0, 600001)
        cdf = erf(2 * xs / np.sqrt(np.pi)) - (4 / np.pi) * xs * np.exp(
            -4 * xs ** 2 / np.pi)
        cdf[0] = 0.0
        cdf = np.maximum.accumulate(cdf)
        u = rng.random(n) * cdf[-1]
        return np.interp(u, cdf, xs)
    raise ValueError(taban)


def main(ad):
    tip, par = KONFIG[ad]
    t_bas = time.time()
    print(f"=== KONFIG {ad}  ({tip} {par}) ===", flush=True)

    d = np.load(HERE / "128_odl_zeros6_2e6_zeros.npz")
    Z = np.sort(np.asarray(d["zeros"], dtype=float))
    zr = Z[len(Z) - 300000:]
    t0, t1 = float(zr[0]), float(zr[-1])
    L_hedef = float(np.log(0.5 * (t0 + t1) / TWO_PI))

    # ---- merdiven (τ ≤ 1.00) ----
    from sympy import primerange
    lim = int(np.exp(1.00 * L_hedef))
    lad = []
    for p in primerange(2, lim + 1):
        q, m = p, 1
        while q <= lim:
            lad.append((np.log(q), 1.0 / (np.pi * m * np.sqrt(q))))
            q *= p; m += 1
    om_l = np.array([w for w, _ in lad])
    a_l = np.array([a for _, a in lad])
    print(f"merdiven: {len(lad)} çizgi (τ≤1.00, L={L_hedef:.3f})", flush=True)

    # ---- TABAN kesim: erfc(0.68, 0.125) — her konfigürasyonda ----
    tau = om_l / L_hedef
    wgt = np.array([0.5 * math.erfc((x - TAU_C) / DELTA) for x in tau])
    a_ham = a_l.copy()
    a_l = a_l * wgt
    print(f"  erfc-kesim τ_c={TAU_C} Δ={DELTA}: Σa·w/Σa="
          f"{a_l.sum()/a_ham.sum():.4f}  (w>0.01 olan çizgi: "
          f"{int((wgt>0.01).sum())})", flush=True)

    z_onbellek = SCRATCH / "z_taban.npy"

    if tip == "itme":
        # İ1: taban çözümünü oku, Newton'u tekrarlama (merdiven özdeş)
        if not z_onbellek.exists():
            raise SystemExit("HATA: önce 'taban' konfigürasyonu koşulmalı "
                             f"({z_onbellek} yok).")
        z = np.load(z_onbellek)
        print(f"  taban (erfc-0.68) Newton çözümü okundu ({len(z)} nokta).",
              flush=True)
        z = itme_supur(z, par["eps"], par["nsup"])

    elif tip == "titresim":
        # AYIRT EDİCİ KONTROL: yapısız (inkoherent) titreşim — 152'nin H-B'si,
        # ama TABAN erfc-0.68 çözümü üzerinde ve σ_ds² N-serisiyle eşleşen
        # dozda. Aynı σ_ds², kısa-menzil yapısı YOK.
        if not z_onbellek.exists():
            raise SystemExit("HATA: önce 'taban' konfigürasyonu koşulmalı "
                             f"({z_onbellek} yok).")
        z = np.load(z_onbellek)
        sj = par["sigma_j"]
        gbar_t = 1.0 / rvm_d(z.mean())
        rj = np.random.default_rng(7)
        z = z + rj.normal(0.0, sj * gbar_t, size=len(z))
        capraz = int(np.sum(np.diff(z) < 0))
        z = np.sort(z)
        print(f"  titreşim σ_j={sj} (σ={sj*gbar_t:.4f}, ḡ_t={gbar_t:.4f}), "
              f"rng(7); sıra bozan çift: {capraz} "
              f"({100.0*capraz/len(z):.2f}%)", flush=True)

    elif tip == "boya":
        # İ2: bağımsız gap'lerle taban + tek geçiş merdiven boyası
        n0 = int(np.ceil(rvm_N(t0)))
        rng_g = np.random.default_rng(11)
        s = gap_ornekle(par["taban"], 299999, rng_g)
        s = s / s.mean()                     # ortalama tam 1'e normalize
        u = np.concatenate(([0.0], np.cumsum(s))) + n0     # 300000 hedef
        print(f"  taban gap'leri: {par['taban']}, n={len(s)}, "
              f"ort={s.mean():.6f}, var={s.var():.4f}, "
              f"P(s<0.3)={np.mean(s<0.3):.5f}", flush=True)
        # pürüzsüz sayımın tersi (merdiven YOK) — 152'nin 1. adımıyla aynı
        x = np.full_like(u, 0.5 * (t0 + t1))
        for _ in range(60):
            F = rvm_N(x) - u
            x = np.clip(x - F / rvm_d(x), 100.0, None)
            if np.max(np.abs(F)) < 1e-9:
                break
        print(f"  pürüzsüz taban: maks|F|={np.max(np.abs(rvm_N(x)-u)):.2e}, "
              f"t∈[{x[0]:.1f},{x[-1]:.1f}]", flush=True)
        # tek geçiş merdiven boyası: z = x − S(x)/N'(x)
        ctx = mp.get_context("spawn")
        pool = ctx.Pool(NWORK, initializer=_init, initargs=(om_l, a_l))
        try:
            parts = np.array_split(x, NWORK)
            res = pool.map(_work, parts)
            S = np.concatenate([r[0] for r in res])
        finally:
            pool.close(); pool.join()
        dz = -S / rvm_d(x)
        z = x + dz
        capraz = int(np.sum(np.diff(z) <= 0))
        print(f"  merdiven boyası (tek geçiş): |Δz|ort={np.mean(np.abs(dz)):.5f}"
              f" maks={np.max(np.abs(dz)):.5f} (ḡ={1/rvm_d(x.mean()):.4f}); "
              f"sıra bozan çift {capraz} ({100.0*capraz/len(z):.2f}%)",
              flush=True)

    else:
        # kontrol/taban: 152'nin A4'ü — öz-tutarlı Newton
        n0 = int(np.ceil(rvm_N(t0))); n1 = n0 + 300000
        ns = np.arange(n0, n1, dtype=float)
        z = np.full_like(ns, 0.5 * (t0 + t1))
        for _ in range(30):
            F = rvm_N(z) - ns
            z = np.clip(z - F / rvm_d(z), 100.0, None)
            if np.max(np.abs(F)) < 1e-9:
                break
        print(f"  pürüzsüz çözüm: maks|F|="
              f"{np.max(np.abs(rvm_N(z)-ns)):.2e}", flush=True)
        gbar_t = 1.0 / rvm_d(z.mean())
        ctx = mp.get_context("spawn")
        pool = ctx.Pool(NWORK, initializer=_init, initargs=(om_l, a_l))
        try:
            for it in range(NEWTON_IT):
                parts = np.array_split(z, NWORK)
                res = pool.map(_work, parts)
                S = np.concatenate([r[0] for r in res])
                Sp = np.concatenate([r[1] for r in res])
                F = rvm_N(z) + S - ns
                payda = np.maximum(rvm_d(z) + Sp, 0.3 * rvm_d(z))
                adim = np.clip(0.8 * F / payda, -1.0 * gbar_t, 1.0 * gbar_t)
                z = np.clip(z - adim, 100.0, None)
                aF = np.abs(F)
                mf = float(aF.max())
                if it % 2 == 0 or it == NEWTON_IT - 1 or mf < 1e-3:
                    print(f"  Newton {it:2d}: maks|F|={mf:.3e}  "
                          f"medyan|F|={np.median(aF):.3e}  "
                          f"%99|F|={np.percentile(aF, 99):.3e}  "
                          f"({time.time()-t_bas:.0f}s)", flush=True)
                if mf < 1e-3:
                    break
        finally:
            pool.close(); pool.join()

    z = np.sort(z)
    if tip == "erfc":
        np.save(z_onbellek, z)
        print(f"  z kaydedildi -> {z_onbellek}", flush=True)

    # ================= 152/151b ZİNCİRİNİN AYNISI =================
    g = np.diff(z)
    mid = 0.5 * (z[:-1] + z[1:])
    Lw = np.log(mid / TWO_PI)
    L = float(Lw.mean())
    ds = g * Lw / TWO_PI - 1
    Nn = len(ds)
    tt = (mid - mid.mean()) / (mid[-1] - mid[0])
    qs = pk(min(int(np.exp(0.52 * L)), 720))
    fr = np.array([np.log(q) for q in qs])
    C = 3 + 2 * len(fr)
    XtX = np.zeros((C, C)); Xty = np.zeros(C)
    for s0 in range(0, Nn, 40000):
        sl = slice(s0, min(s0 + 40000, Nn))
        arg = np.outer(mid[sl], fr)
        Xc = np.empty((sl.stop - sl.start, C))
        Xc[:, 0] = 1; Xc[:, 1] = tt[sl]; Xc[:, 2] = tt[sl]**2
        Xc[:, 3::2] = np.cos(arg); Xc[:, 4::2] = np.sin(arg)
        XtX += Xc.T @ Xc; Xty += Xc.T @ ds[sl]
        del Xc, arg
    b = np.linalg.solve(XtX, Xty)
    eta = np.empty(Nn)
    for s0 in range(0, Nn, 40000):
        sl = slice(s0, min(s0 + 40000, Nn))
        arg = np.outer(mid[sl], fr)
        eta[sl] = ds[sl] - (b[0] + b[1]*tt[sl] + b[2]*tt[sl]**2 +
                            np.cos(arg) @ b[3::2] + np.sin(arg) @ b[4::2])
        del arg
    c1 = float(np.mean(eta[:-1] * eta[1:]))
    s_ds = float(np.var(ds)); s_eta = float(np.var(eta))
    print(f"\nS3 — {ad}: σ_ds²={s_ds:.4f}  σ_η²={s_eta:.4f}  c₁={c1:+.5f}"
          f"   [gerçek: 0.1674 / 0.0227 / −0.01158]", flush=True)

    # --- küçük-s göstergesi (ölçüm zincirine DOKUNMAZ, yalnız rapor) ---
    ss = ds + 1.0            # zincirin kendi açılmış aralığı
    ps = {k: float(np.mean(ss < k)) for k in (0.1, 0.2, 0.3, 0.5, 0.7)}
    print("P(s) küçük-s: " + "  ".join(
        f"P(s<{k})={ps[k]:.5f}[ger {GERCEK_PS[k]:.5f}]" for k in
        (0.1, 0.2, 0.3, 0.5)), flush=True)

    dsA = 0.5 * (ds[:-1] + ds[1:]); dsA -= dsA.mean()
    allq = pk(int(np.exp(0.86 * L)))
    allw = np.array([np.log(q) for q in allq])
    T = mid[-1] - mid[0]; dres = TWO_PI / T
    e0, e1 = eta[:-1], eta[1:]
    m0 = mid[:-1]
    rng = np.random.default_rng(21)
    print("\nS1/S2 — bant   sentetik Γrot(Re,Im)  R_sent  | gerçek Γrot  R_ger",
          flush=True)
    satirlar = []
    for lo, hi in BANTLAR:
        tumu = [(q, w) for q, w in zip(allq, allw) if lo < w / L <= hi]
        cand = tumu
        if len(cand) > 220:
            idx = rng.choice(len(cand), 220, replace=False)
            cand = [cand[i] for i in idx]
        crN = 0j; crO = 0j; on0 = off0 = 0.0
        nkul = 0
        for q, w in cand:
            j = np.searchsorted(allw, w)
            koms = [allw[k] for k in (j - 1, j + 1)
                    if 0 <= k < len(allw) and abs(allw[k] - w) > 1e-12]
            gap = min(abs(w - k) for k in koms)
            if gap < 2.5 * dres:
                continue
            nkul += 1
            for W, hedef in ((w, True), (w + gap / 2, False)):
                cw, sw = np.cos(W * m0), np.sin(W * m0)
                zc = complex(2 * np.mean(e0 * cw), -2 * np.mean(e0 * sw))
                zp = complex(2 * np.mean(e1 * cw), -2 * np.mean(e1 * sw))
                cr = zp * np.conj(zc) * np.exp(1j * TWO_PI * W / L)
                if hedef:
                    crN += cr; on0 += abs(zc)**2
                else:
                    crO += cr; off0 += abs(zc)**2
        tb = 0.5 * (lo + hi)
        gr = GERCEK[round(tb, 4)]
        if nkul == 0 or abs(on0 - off0) < 1e-300:
            print(f"  {tb:.4f}  (—,—)        —      "
                  f"| ({gr[0]:+.3f},{gr[1]:+.3f})  {gr[2]:.2f}    "
                  f"[ölçülemedi: kullanılan aday {nkul}, payda "
                  f"{on0-off0:.2e}]", flush=True)
            satirlar.append({"tau": round(tb, 4), "re": None, "im": None,
                             "R": None, "R_artik": None, "n_aday": nkul})
            continue
        G = (crN - crO) / (on0 - off0)
        A = TWO_PI * tb
        phm = float(np.angle(G))
        best = None
        for R in np.linspace(0.0, 3.0, 121):
            wg = np.exp(-A * R * dsA)
            M = np.mean(wg * np.exp(-1j * A * dsA)) / np.mean(wg)
            if best is None or abs(np.angle(M) - phm) < best[0]:
                best = (abs(np.angle(M) - phm), R)
        artik = best[0]
        print(f"  {tb:.4f}  ({G.real:+.3f},{G.imag:+.3f})   {best[1]:4.2f}   "
              f"| ({gr[0]:+.3f},{gr[1]:+.3f})  {gr[2]:.2f}    "
              f"[R-artık {artik:.3f}, aday {nkul}]", flush=True)
        satirlar.append({"tau": round(tb, 4), "re": G.real, "im": G.imag,
                         "R": best[1], "R_artik": artik, "n_aday": nkul})

    ozet = {"ad": ad, "tip": tip, "par": par, "newton_it": NEWTON_IT,
            "sigma_ds2": s_ds, "sigma_eta2": s_eta, "c1": c1,
            "ps": ps, "bantlar": satirlar, "sure_s": time.time() - t_bas}
    (SCRATCH / f"ozet_{ad}.json").write_text(json.dumps(ozet, indent=1))
    print(f"\n[{ad}] bitti — {(time.time()-t_bas)/60:.1f} dk", flush=True)


if __name__ == "__main__":
    if len(sys.argv) != 2 or sys.argv[1] not in KONFIG:
        raise SystemExit(f"kullanım: 153_gaz.py <{'|'.join(KONFIG)}>")
    main(sys.argv[1])
