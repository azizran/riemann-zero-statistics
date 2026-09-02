"""
164 — GENLİK-SADAKATLİ SENTETİK GAZ İNŞASI
==========================================
152'nin sönümlü/kelepçeli Newton'u  N̄(z) + S(z) = n  denklemini EKSİK
gerçekliyor: 163 §1c'nin ölçtüğü gerçekleşen/nominal merdiven genliği
keskin'de 0.54, A4'te 0.85 (gerçek veri 1.00–1.04). Bu modül denklemi
SADAKATLE çözer.

────────────────────────────────────────────────────────────────────────
NEDEN 152 EKSİK KURUYOR — TEŞHİS (bu koşuda sayıyla doğrulanıyor)
────────────────────────────────────────────────────────────────────────
F(z) ≡ N̄(z) + S(z) − n ,  S(t) = −Σ_q a_q sin(ω_q t) ,  a_q = 1/(πm√q)

    N̄'(t) = log(t/2π)/2π = 1.9147          (bu pencerede)
    rms S' = sqrt(½ Σ (a_q ω_q)²) = 1.8940   ← N̄' ile AYNI mertebe!
    max |S'| ≤ Σ a_q ω_q = 259.9

Yani F MONOTON DEĞİLDİR: F' = N̄' + S' işaret değiştirir. 152'nin
Newton'u tam bu yüzden kelepçelidir:

    payda = max(N̄' + S', 0.3·N̄')      ← negatif türevi 0.3·N̄'e kırpar
    adım  = clip(0.8·F/payda, ±ḡ)      ← sönüm 0.8 + ±ḡ kelepçe
    20 iterasyon, |F| < 1e-3'te dur

Kelepçe, F'nin dik/negatif-türevli bölgelerinde adımı SİSTEMATİK olarak
küçültür; çözüm merdivenin tepe-dip genliğini tam kuramadan durur. Sonuç:
ds'nin çizgi genlikleri nominalin altında kalır (keskin 0.54).

────────────────────────────────────────────────────────────────────────
SADAKATLİ ÇÖZÜM — SIRALI İLK-KÖK
────────────────────────────────────────────────────────────────────────
F monoton olmadığı için "z_n = F(z)=n'in kökü" tanımı tek başına
sıralılığı garanti etmez. Doğru tanım:

    z_n = (verilen bir alt sınırdan itibaren) F(z) = n'in İLK kökü

Bu tanım sıralılığı ÖZDEŞ olarak garanti eder: seviye n+1 > n olduğundan
n+1'in ilk kökü n'inkinden küçük olamaz.

ALGORİTMA
  1. İnce ızgara. Adım h ≤ (en kısa dalga boyu)/10 = (2π/ω_max)/10
     = 0.0522; ızgara-yakınsama sınavı h = 0.015 dedi (~10.4M nokta).
     S ızgarada BİR KEZ, TAM (yaklaşıksız) hesaplanır — 7 süreçli havuz.
  2. G_k = N̄(z_k) + S(z_k). Koşan maksimum M_k = max(G_0..G_k) MONOTONDUR,
     ve M_k ≥ n'i sağlayan İLK k'da G_k = M_k'dir. Dolayısıyla
     k_n = searchsorted(M, n) her n için İLK ızgara-geçişini verir;
     braket [z_{k−1}, z_k] içinde G_{k−1} < n ≤ G_k.
  3. Braket içinde KORUMALI NEWTON: doğrusal ara değerden başlanır, her
     adımda TAM S ve S' ile Newton denenir; adım braketten çıkarsa ya da
     F' ≤ 0 ise ikiye bölme yapılır. Braket her adımda daraltılır, yani
     yakınsama GARANTİLİDİR (en kötü hâlde ikiye bölme).
     Yalnız yakınsamamış noktalar üzerinde çalışılır (aktif küme küçülür).
  4. Yakınsama ölçütü: |F| ≤ 1e-8 HER teknede (raporlanır, gizlenmez).

Kübik interpolasyon YERİNE doğrusal ara değer + TAM Newton kullanıldı:
kübik spline hatası h⁴·rms(S⁗)/384 = 2.7e−05 (rms S⁗ = 1665) yani
zaten TAM değerlendirme gerekiyordu; doğrusal ara değer + 2–3 TAM Newton
adımı hem daha ucuz hem de |F| ölçütünü İNTERPOLANTTA değil GERÇEK
F'de sağlıyor.

────────────────────────────────────────────────────────────────────────
Kullanım
  164_insa.py Skeskin            # sadakatli keskin (τ≤1.00, pencere yok)
  164_insa.py SA4                # sadakatli A4 (erfc 0.68/0.125)
  164_insa.py SA1                # sadakatli A1 (erfc 0.75/0.10) — yedek
  164_insa.py eski_keskin        # 152'nin Newton'unu AYNEN yeniden üret
  164_insa.py eski_A4            # 152'nin Newton'u + erfc 0.68/0.125
  164_insa.py NKkeskin           # kök-seçimi kontrolü: EN YAKIN kök
"""
import json
import math
import multiprocessing as mp
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCRR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
            "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
SCR164 = SCRR / "164"
SCR155 = SCRR / "155"
SCR152 = SCRR / "152"
TWO_PI = 2 * np.pi

NWORK = 7
NPT = 10000          # süreç içi nokta bloğu
BLOK = 800           # süreç içi çizgi bloğu
NZERO = 300000       # üretilen tekne sayısı (152 ile aynı)
HIZGARA = 0.015      # ızgara adımı (en kısa dalga boyunun 1/35'i)
# h SEÇİMİ ÖLÇÜLDÜ (V-ızgara, §7): 5000 teknelik smoke koşusunda
#   h = 0.0400 → h = 0.00375'e göre 4/5000 tekne farklı (alt-ızgara çift
#                geçişi kaçırılıyor; maks |Δz| = 0.28)
#   h = 0.0150 → h = 0.00375 ile 0/5000 fark (maks |Δz| = 0.0)
#   h = 0.0075 → h = 0.00375 ile 0/5000 fark
# yani h = 0.015 ızgara-yakınsamıştır. 1/10-dalga-boyu ölçütü (0.0522)
# TEK BAŞINA yetmiyor; ilk-kök seçimi daha ince ızgara istiyor.
TOL = 1e-8           # |F| ölçütü

KONFIG = {
    "Skeskin":     ("sadakatli", {}),
    "SA4":         ("sadakatli", {"tau_c": 0.68, "delta": 0.125}),
    "SA1":         ("sadakatli", {"tau_c": 0.75, "delta": 0.10}),
    "eski_keskin": ("eski", {}),
    "eski_A4":     ("eski", {"tau_c": 0.68, "delta": 0.125}),
    "NKkeskin":    ("enyakin", {}),
    # SEVİYE KONVANSİYONU c = −½ (§3d/§6'nın bulgusu): gerçek sayma
    # fonksiyonu sıfırda +1 atlar, asal-toplam onun DÜZGÜN sürümüdür ⇒
    # N̄(γ_n) + S_düz(γ_n) = n − ½. 152 (ve ilk koşu) c = 0 kullandı.
    "Hkeskin":     ("sadakatli", {"c": -0.5}),
    "HA4":         ("sadakatli", {"tau_c": 0.68, "delta": 0.125, "c": -0.5}),
}


def rvm_N(t):
    x = t / TWO_PI
    return x * np.log(x / np.e) + 7 / 8


def rvm_d(t):
    return np.log(t / TWO_PI) / TWO_PI


# --------------------------------------------------------------- merdiven
def merdiven(L_hedef, tau_c=None, delta=None, tau_ust=1.00):
    """152/153/155 ile BİREBİR aynı merdiven (+ isteğe bağlı erfc penceresi).

    Döner: om, a (pencereli), a_ham (penceresiz), w (pencere)."""
    from sympy import primerange
    lim = int(np.exp(tau_ust * L_hedef))
    lad = []
    for p in primerange(2, lim + 1):
        q, m = p, 1
        while q <= lim:
            lad.append((np.log(q), 1.0 / (np.pi * m * np.sqrt(q))))
            q *= p
            m += 1
    om = np.array([x for x, _ in lad])
    a_ham = np.array([x for _, x in lad])
    if tau_c is None:
        return om, a_ham, a_ham, np.ones_like(a_ham)
    w = np.array([0.5 * math.erfc((x - tau_c) / delta) for x in om / L_hedef])
    return om, a_ham * w, a_ham, w


# ------------------------------------------------------- S(z) ve S'(z)
def S_ve_Sp(z, om, a, deriv=True, blok=BLOK, npt=NPT):
    S = np.zeros_like(z)
    Sp = np.zeros_like(z) if deriv else None
    for b0 in range(0, len(om), blok):
        w = om[b0:b0 + blok]
        aa = a[b0:b0 + blok]
        aw = aa * w if deriv else None
        for s0 in range(0, len(z), npt):
            sl = slice(s0, min(s0 + npt, len(z)))
            arg = np.outer(z[sl], w)
            S[sl] += -np.sin(arg) @ aa
            if deriv:
                Sp[sl] += -np.cos(arg) @ aw
            del arg
    return (S, Sp) if deriv else (S, None)


_G = {}


def _init(om, a):
    _G["om"] = om
    _G["a"] = a


def _work_S(zc):
    return S_ve_Sp(zc, _G["om"], _G["a"], deriv=False)[0]


def _work_SP(zc):
    return S_ve_Sp(zc, _G["om"], _G["a"], deriv=True)


def _par(pool, fn, z, nw=NWORK):
    if len(z) == 0:
        return None
    parts = np.array_split(z, min(nw, max(1, len(z) // 200 + 1)))
    return pool.map(fn, parts)


def S_par(pool, z):
    return np.concatenate(_par(pool, _work_S, z))


def SSp_par(pool, z):
    r = _par(pool, _work_SP, z)
    return (np.concatenate([x[0] for x in r]),
            np.concatenate([x[1] for x in r]))


# --------------------------------------------------------- SADAKATLİ ÇÖZÜM
def _log(*a):
    print(*a, flush=True)


def coz_sadakatli(om, a, ns, t0, pool, h=HIZGARA, log=_log):
    """N̄(z)+S(z)=n'in SIRALI İLK-KÖK çözümü (ızgara braketi + korumalı Newton)."""
    tb = time.time()
    n0, n1 = float(ns[0]), float(ns[-1])
    # ızgara aralığı: pürüzsüz köklerin ±pay'ı (|S| ≤ Σa_q)
    Sam = float(np.abs(a).sum())
    pay = (Sam + 5.0) / rvm_d(t0) + 5.0
    zlo = t0
    for _ in range(60):                     # pürüzsüz N̄(z)=n0
        zlo -= (rvm_N(zlo) - n0) / rvm_d(zlo)
    zhi = t0
    for _ in range(60):
        zhi -= (rvm_N(zhi) - n1) / rvm_d(zhi)
    zlo -= pay
    zhi += pay
    ng = int(np.ceil((zhi - zlo) / h)) + 1
    log(f"  ızgara: [{zlo:.1f}, {zhi:.1f}] adım h={h} → {ng} nokta "
        f"(pay={pay:.1f}, Σa={Sam:.3f})")

    zg = zlo + h * np.arange(ng)
    S = S_par(pool, zg)
    G = rvm_N(zg) + S
    log(f"  S ızgarada hesaplandı ({time.time()-tb:.0f}s); "
        f"G∈[{G.min():.1f},{G.max():.1f}]")

    # ızgara tanısı: F'nin monotonluk ihlali
    dG = np.diff(G)
    fr_neg = float(np.mean(dG < 0))
    log(f"  ızgara tanısı: ΔG<0 kesri = {fr_neg:.4f} "
        f"(F' işaret değiştiriyor); min ΔG = {dG.min():+.4f}, "
        f"maks ΔG = {dG.max():+.4f}")

    M = np.maximum.accumulate(G)
    k = np.searchsorted(M, ns.astype(float), side="left")
    if k.max() >= ng or k.min() < 1:
        raise SystemExit("ızgara aralığı yetersiz (k sınırda)")
    lo = zg[k - 1]
    hi = zg[k]
    Glo = G[k - 1]
    Ghi = G[k]
    # doğrusal ara değer (braket içinde)
    z = lo + (ns - Glo) / np.maximum(Ghi - Glo, 1e-300) * h
    z = np.clip(z, lo, hi)
    log(f"  braketler kuruldu; ilk-kök hücreleri: benzersiz {len(np.unique(k))}"
        f" / {len(k)}")

    # ---- korumalı Newton (aktif küme) --------------------------------
    # KAYAN-NOKTA TABANI: F = (N̄(z) − n) + S(z); N̄ ≈ 1.7e6 olduğundan
    # çıkarma iptali |F|'yi ulp(n) ≈ 2.3e−10'un altına indiremez. İç eşik
    # bu tabanın hemen üstüne (2e−9) konur; RAPOR ÖLÇÜTÜ TOL = 1e−8 ve
    # taban ondan 4 kat küçüktür.
    TOLIN = max(2e-9, 8.0 * np.spacing(float(ns[-1])))
    log(f"  kayan-nokta tabanı: ulp(n)={np.spacing(float(ns[-1])):.2e}; "
        f"iç eşik TOLIN={TOLIN:.2e}, rapor ölçütü TOL={TOL:.0e}")
    akt = np.arange(len(ns))
    nbis = 0
    for it in range(60):
        za = z[akt]
        Sa, Spa = SSp_par(pool, za)
        F = rvm_N(za) + Sa - ns[akt]
        Fp = rvm_d(za) + Spa
        pos = F > 0                                  # braket güncelle
        hi[akt] = np.where(pos, za, hi[akt])
        lo[akt] = np.where(pos, lo[akt], za)
        kal = np.abs(F) > TOLIN
        log(f"  Newton {it:2d}: aktif={len(akt):7d}  "
            f"maks|F|={np.abs(F).max():.3e}  "
            f"medyan|F|={np.median(np.abs(F)):.3e}  "
            f"kalan={int(kal.sum()):7d}  ({time.time()-tb:.0f}s)")
        if not kal.any():
            break
        akt, za, F, Fp = akt[kal], za[kal], F[kal], Fp[kal]
        # Newton adımı; braket dışına çıkarsa ya da Fp ≤ 0 ise ikiye böl
        with np.errstate(divide="ignore", invalid="ignore"):
            zn = za - F / Fp
        kotu = ((Fp <= 0) | ~np.isfinite(zn)
                | (zn <= lo[akt]) | (zn >= hi[akt]))
        nbis += int(kotu.sum())
        z[akt] = np.where(kotu, 0.5 * (lo[akt] + hi[akt]), zn)
    # son bir TAM değerlendirme (interpolant DEĞİL, gerçek S)
    Sf = S_par(pool, z)
    Ffin = rvm_N(z) + Sf - ns
    nas = int((np.abs(Ffin) > TOL).sum())
    log(f"  BİTTİ: maks|F| = {np.abs(Ffin).max():.3e}  "
        f"(|F|>{TOL:g} olan tekne: {nas})  "
        f"ikiye-bölme adımı: {nbis}  ({time.time()-tb:.0f}s)")
    tani = dict(fr_dG_neg=fr_neg, ng=ng, h=h,
                maxF=float(np.abs(Ffin).max()),
                medF=float(np.median(np.abs(Ffin))),
                ulp=float(np.spacing(float(ns[-1]))), tolin=float(TOLIN),
                nF_asan=nas, nbis=nbis,
                hucre_benzersiz=int(len(np.unique(k))))
    return z, Ffin, tani


# -------------------------------- KÖK-SEÇİMİ KONTROLÜ: "EN YAKIN KÖK"
def coz_enyakin(om, a, ns, t0, t1, pool, tohum=None, log=_log):
    """SİSTEMATİK KONTROLÜ — İLK-KÖK yerine 152'nin çözümüne EN YAKIN kök.

    F monoton olmadığında 'z_n = F=n'in kökü' çok değerlidir. İlk-kök
    kuralı sıralılığı garanti eder ama seçimin kendisi bir MODEL
    kararıdır. Bu kontrol aynı denklemi 152'nin (sönümlü) çözümünden
    başlatılan SÖNÜMSÜZ Newton ile çözer: yakınsadığı kök 'ilk' değil
    'başlangıca en yakın' köktür. Böylece ölçülen farkların ne kadarı
    ÇÖZÜM SADAKATİNDEN, ne kadarı KÖK SEÇİMİNDEN geliyor ayrılır.
    (Sıralılık burada GARANTİ DEĞİL — bozulursa raporlanır.)"""
    tb = time.time()
    z0 = SCR152 / "z_keskin.npy"
    z = np.sort(np.load(z0)) if z0.exists() else None
    if z is None or len(z) != len(ns):
        raise SystemExit("en-yakın kontrolü 152'nin z_keskin.npy'sini ister")
    log(f"  başlangıç: 152'nin sönümlü çözümü ({len(z)} nokta)")
    for it in range(30):
        S, Sp = SSp_par(pool, z)
        F = rvm_N(z) + S - ns
        Fp = rvm_d(z) + Sp
        mf = float(np.abs(F).max())
        log(f"  Newton(enyakın) {it:2d}: maks|F|={mf:.3e}  "
            f"medyan|F|={np.median(np.abs(F)):.3e}  ({time.time()-tb:.0f}s)")
        if mf < 1e-9:
            break
        with np.errstate(divide="ignore", invalid="ignore"):
            adim = F / np.where(np.abs(Fp) > 1e-9, Fp, np.nan)
        adim = np.where(np.isfinite(adim), adim, 0.0)
        z = z - np.clip(adim, -0.5, 0.5)
    S = S_par(pool, z)
    F = rvm_N(z) + S - ns
    sira = int((np.diff(z) <= 0).sum())
    log(f"  BİTTİ(enyakın): maks|F|={np.abs(F).max():.3e}  "
        f"sıra bozan çift: {sira}")
    return np.sort(z), F, dict(maxF=float(np.abs(F).max()),
                               nF_asan=int((np.abs(F) > TOL).sum()),
                               sira_bozan=sira)


# ------------------------------------------------------- 152'NİN NEWTON'U
def coz_eski(om, a, ns, t0, t1, pool, log=_log):
    """152_gaz.py'nin sönümlü/kelepçeli Newton'u — SATIR SATIR AYNI."""
    tb = time.time()
    NEWTON_IT = 20
    z = np.full_like(ns, 0.5 * (t0 + t1))
    for _ in range(30):
        F = rvm_N(z) - ns
        z = np.clip(z - F / rvm_d(z), 100.0, None)
        if np.max(np.abs(F)) < 1e-9:
            break
    log(f"  pürüzsüz çözüm: maks|F|={np.max(np.abs(rvm_N(z)-ns)):.2e}")
    gbar_t = 1.0 / rvm_d(z.mean())
    for it in range(NEWTON_IT):
        S, Sp = SSp_par(pool, z, )
        F = rvm_N(z) + S - ns
        payda = np.maximum(rvm_d(z) + Sp, 0.3 * rvm_d(z))
        adim = np.clip(0.8 * F / payda, -1.0 * gbar_t, 1.0 * gbar_t)
        z = np.clip(z - adim, 100.0, None)
        aF = np.abs(F)
        mf = float(aF.max())
        if it % 4 == 0 or it == NEWTON_IT - 1 or mf < 1e-3:
            log(f"  Newton(152) {it:2d}: maks|F|={mf:.3e}  "
                f"medyan|F|={np.median(aF):.3e}  ({time.time()-tb:.0f}s)")
        if mf < 1e-3:
            break
    S = S_par(pool, z)
    F = rvm_N(z) + S - ns
    log(f"  152-Newton bitti: maks|F|={np.abs(F).max():.3e}  "
        f"medyan|F|={np.median(np.abs(F)):.3e}")
    return np.sort(z), F, dict(maxF=float(np.abs(F).max()),
                               medF=float(np.median(np.abs(F))))


# ----------------------------------------------------------------- ana
def main(ad, h=HIZGARA, nz=NZERO):
    tip, par = KONFIG[ad]
    tbas = time.time()
    SCR164.mkdir(parents=True, exist_ok=True)
    SCR155.mkdir(parents=True, exist_ok=True)
    print(f"=== 164 İNŞA {ad} ({tip} {par}) h={h} nz={nz} ===", flush=True)

    d = np.load(HERE / "128_odl_zeros6_2e6_zeros.npz")
    Z = np.sort(np.asarray(d["zeros"], dtype=float))
    zr = Z[len(Z) - 300000:]
    t0, t1 = float(zr[0]), float(zr[-1])
    L_hedef = float(np.log(0.5 * (t0 + t1) / TWO_PI))

    om, a, a_ham, wgt = merdiven(L_hedef, par.get("tau_c"), par.get("delta"))
    print(f"merdiven: {len(om)} çizgi (τ≤1.00, L={L_hedef:.4f})", flush=True)
    if par.get("tau_c") is not None:
        print(f"  erfc-kesim τ_c={par['tau_c']} Δ={par['delta']}: "
              f"Σa·w/Σa={a.sum()/a_ham.sum():.4f}  "
              f"(w>0.01 çizgi: {int((wgt>0.01).sum())})", flush=True)
    print(f"  N̄'={rvm_d(0.5*(t0+t1)):.4f}  rms S'="
          f"{np.sqrt(0.5*np.sum((a*om)**2)):.4f}  Σa_qω_q="
          f"{np.sum(a*om):.1f}  rms S={np.sqrt(0.5*np.sum(a**2)):.4f}",
          flush=True)

    n0 = int(np.ceil(rvm_N(t0)))
    ns = np.arange(n0, n0 + nz, dtype=float) + par.get("c", 0.0)
    if par.get("c"):
        print(f"  SEVİYE KONVANSİYONU: N̄+S = n {par['c']:+g}", flush=True)

    ctx = mp.get_context("spawn")
    pool = ctx.Pool(NWORK, initializer=_init, initargs=(om, a))
    try:
        if tip == "sadakatli":
            z, F, tani = coz_sadakatli(om, a, ns, t0, pool, h=h)
        elif tip == "enyakin":
            z, F, tani = coz_enyakin(om, a, ns, t0, t1, pool)
        else:
            z, F, tani = coz_eski(om, a, ns, t0, t1, pool)
    finally:
        pool.close()
        pool.join()

    dz = np.diff(z)
    sirali = bool(np.all(dz > 0))
    print(f"  sıralılık: {'TAM' if sirali else 'BOZUK'}  "
          f"min Δz={dz.min():.6f}  (sıra bozan çift: {int((dz<=0).sum())})",
          flush=True)

    g = dz
    mid = 0.5 * (z[:-1] + z[1:])
    Lw = np.log(mid / TWO_PI)
    L = float(Lw.mean())
    ds = g * Lw / TWO_PI - 1
    print(f"  σ_ds² = {np.var(ds):.4f}   L={L:.4f}   "
          f"P(s<0.3)={np.mean((g/(TWO_PI/Lw))<0.3):.5f}   "
          f"[gerçek: σ_ds²=0.1674]", flush=True)

    p164 = SCR164 / f"z_{ad}.npy"
    np.save(p164, z)
    np.save(SCR155 / f"z_{ad}.npy", z)     # 155_kos.veri_yukle buradan okur
    tani.update(ad=ad, tip=tip, par=par, L=L, L_hedef=L_hedef,
                sirali=sirali, min_dz=float(dz.min()),
                sigma_ds2=float(np.var(ds)), n=int(len(z)),
                nline=int(len(om)), sure_s=time.time() - tbas)
    (SCR164 / f"insa_{ad}.json").write_text(json.dumps(tani, indent=1))
    print(f"-> {p164}  ({(time.time()-tbas)/60:.1f} dk)", flush=True)

    # eski inşanın MEŞRUİYET sınavı: 152'nin kayıtlı z'siyle karşılaştır
    ref = {"eski_keskin": SCR152 / "z_keskin.npy",
           "eski_A4": SCRR / "154" / "z_A4.npy"}.get(ad)
    if ref is not None and ref.exists():
        zr2 = np.sort(np.load(ref))
        dd = float(np.max(np.abs(zr2 - z))) if len(zr2) == len(z) else float("nan")
        print(f"  [MEŞRUİYET] {ref.name} ile maks|Δz| = {dd:.3e}", flush=True)
        tani["ref_maks_dz"] = dd
        (SCR164 / f"insa_{ad}.json").write_text(json.dumps(tani, indent=1))


if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] not in KONFIG:
        raise SystemExit(f"kullanım: 164_insa.py <{'|'.join(KONFIG)}> [h] [nz]")
    hh = float(sys.argv[2]) if len(sys.argv) > 2 else HIZGARA
    nn = int(sys.argv[3]) if len(sys.argv) > 3 else NZERO
    main(sys.argv[1], hh, nn)
