"""
156 — ÜÇ-KANALLI AĞIRLIK AYRIŞIMI ÇEKİRDEĞİ
==========================================
Soru: Γ_rot'un fazını üreten "yerel sağkalım-ağırlık bağlaşımı"
   Γ_rot(τ) = ⟨w·e^{−iA·dsΔ}⟩ / ⟨w⟩,   w = e^{−A·R·dsΔ}
hangi YEREL büyüklükle korelasyonlu? Tam adımla mı (dsΔ), yoksa adımın
yalnız bir PAYIYLA mı — merdiven dalgasıyla (koherent, çizgi içeriğiyle
paylaşılan) ya da η artığıyla?

AYRIŞIM (0.52-taban regresyonundan, birebir 150/154 zinciri):
   ds = drift + lad + eta,   drift = b0 + b1·tt + b2·tt²
                              lad   = Σ (b_c cos ωt + b_s sin ωt)   (τ≤0.52)
                              eta   = artık
Bond adımı (n ↔ n+1 bağı) her kanal için ayrı:
   X_tam = (ds_n+ds_{n+1})/2,  X_lad, X_eta, X_dri aynı biçimde.
Hepsi ortalaması çıkarılmış (ağırlıkta sabit oran içinde sadeleşir;
X_tam'ın ortalama çıkarması 150 ile birebir aynıdır).

MODEL (tek kanallı, k ∈ {tam, lad, eta}):
   w_k = e^{−A·R_k·X_k},   M_k(R) = ⟨w_k·e^{−iA·X_tam}⟩ / ⟨w_k⟩
FAZ FAKTÖRÜ HER ZAMAN TAM ADIMDIR (fiziksel faz ilerlemesi tam adımdan
gelir); yalnız AĞIRLIK kanalı değişir. R_k, ölçülen fazı eşleyecek
şekilde çözülür; sonra AŞIRI-BELİRLEME skoru |Re M_k − Re Γ| bakılır —
tek serbest parametre fazı yediği için Re bağımsız bir sınavdır.

ÖLÇEK NORMALİZASYONU (sayısal, fiziği değiştirmez): kanalların
varyansları çok farklı olduğundan ham R_k'ler karşılaştırılamaz.
Gauss limitinde faz = A²·R_k·Cov(X_k, X_tam) olduğundan
   Y_k = X_k · σΔ² / Cov(X_k, X_tam)
tanımlanır; Y_k üzerinden çözülen R̃_k, Gauss limitinde kanaldan bağımsız
olarak R_tam'a eşittir — yani R̃ doğrudan karşılaştırılabilir ve arama
aralığı üç kanalda da aynı olur. Ham karşılık: R_k = R̃_k · s_k,
s_k = σΔ²/Cov(X_k, X_tam). (Gauss limiti burada YAKLAŞIK değil, sadece
ölçek seçimi; kök ampirik dağılımla çözülür.)

İKİ KANALLI SIRT: w = e^{−A(R̃_l·Y_lad + R̃_e·Y_eta)}. R̃_l ızgarası
üzerinde her sütunda fazı eşleyen R̃_e kökü bulunur (sırt eğrisi), sonra
sırt üzerinde |Re M − Re Γ| en küçük olan çift seçilir.

KOPYA-KAYMASI DENETİMİ: Γ_rot ölçüm döngüsü 154_cekirdek.olc'den
kopyadır; R kök-bulucu ve M hesabı KOPYA DEĞİL, 154_cekirdek'in kendi
fonksiyonları import edilerek çağrılır. Ayrıca her koşuda 154'ün
kaydettiği Γ/R değerleriyle karşılaştırma yazdırılır.
"""

import importlib
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "154_configs"))
CEK = importlib.import_module("154_cekirdek")

TWO_PI = 2 * np.pi

BANTLAR_STD = CEK.BANTLAR_STD          # [(0.525,0.55),(0.55,0.62),
                                       #  (0.62,0.70),(0.70,0.78),(0.78,0.85)]
# görevin istediği üç bant = std'nin 2., 3., 4. bandı
GOREV_BANTLARI = [(0.55, 0.62), (0.62, 0.70), (0.70, 0.78)]

# R̃ arama aralığı: 154'ün (−3, 8) aralığından geniş, çünkü normalize
# edilmiş kanallarda Gauss ölçeği ampirik kökten uzak düşebiliyor.
RMIN, RMAX, RADIM = -8.0, 20.0, 0.05
# sırt taramasında R̃_e aralığı (ham karşılığı gerçekte ×21, A4'te ×6.6)
SIRT_EMIN, SIRT_EMAX, SIRT_EADIM = -6.0, 10.0, 0.05
# faz kapasitesi taraması (kanal tek başına hangi fazları üretebiliyor)
KAP_GRID = np.arange(-4.0, 8.0001, 0.05)


def zincir3(z, taban=0.52, cap=720):
    """150/154'ün ds → 0.52-taban regresyonu; drift/lad/eta AYRI döner."""
    g = np.diff(z)
    mid = 0.5 * (z[:-1] + z[1:])
    Lw = np.log(mid / TWO_PI)
    L = float(Lw.mean())
    ds = g * Lw / TWO_PI - 1
    Nn = len(ds)
    tt = (mid - mid.mean()) / (mid[-1] - mid[0])

    qs = CEK.pk(min(int(np.exp(taban * L)), cap))
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

    drift = b[0] + b[1] * tt + b[2] * tt**2
    lad = np.empty(Nn)
    for s0 in range(0, Nn, 40000):
        sl = slice(s0, min(s0 + 40000, Nn))
        arg = np.outer(mid[sl], fr)
        lad[sl] = np.cos(arg) @ b[3::2] + np.sin(arg) @ b[4::2]
        del arg
    eta = ds - drift - lad

    kapanis = float(np.max(np.abs(ds - (drift + lad + eta))))
    return dict(mid=mid, L=L, ds=ds, drift=drift, lad=lad, eta=eta,
                nq=len(qs), kapanis=kapanis,
                s_ds=float(np.var(ds)), s_eta=float(np.var(eta)),
                s_lad=float(np.var(lad)), s_dri=float(np.var(drift)),
                c1=float(np.mean(eta[:-1] * eta[1:])))


def _bond(x):
    """Bond adımı: komşu ortalaması, ortalaması çıkarılmış."""
    y = 0.5 * (x[:-1] + x[1:])
    return y - y.mean()


def _kapasite(A, Yk, ex):
    """Kanal TEK BAŞINA hangi faz aralığını üretebiliyor? (arg M(R̃) menzili)

    R̃ ızgarası üzerinde arg M taranır; en büyük/en küçük değer ve nerede
    olduğu döner. Ölçülen faz bu menzilin DIŞINDAYSA o kanal fazı hiçbir
    R ile üretemez — tek-kanallı model o bantta ELENİR.
    """
    a = np.empty(len(KAP_GRID)); ne = np.empty(len(KAP_GRID))
    for i, R in enumerate(KAP_GRID):
        M, n = _Mw(R * Yk, A, ex)
        a[i] = np.angle(M); ne[i] = n
    i_mx, i_mn = int(np.argmax(a)), int(np.argmin(a))
    return dict(arg_max=float(a[i_mx]), R_at_max=float(KAP_GRID[i_mx]),
                neff_at_max=float(ne[i_mx]),
                arg_min=float(a[i_mn]), R_at_min=float(KAP_GRID[i_mn]))


def _Mw(x_lin, A, ex):
    """M ve n_eff; x_lin = ağırlık üstelinin R'siz kısmı (w = e^{−A·x_lin})."""
    x = -A * x_lin
    w = np.exp(x - x.max())
    sw = w.sum()
    return (w @ ex) / sw, float(sw * sw / (w @ w))


def _sirt(A, Yl, Ye, ex, phm, Rl_grid, rmin=SIRT_EMIN, rmax=SIRT_EMAX,
          adim=SIRT_EADIM):
    """Faz-eşlenmiş (R̃_l, R̃_e) sırtı: her R̃_l için fazı çözen R̃_e."""
    def M(Rl, Re):
        x = -A * (Rl * Yl + Re * Ye)
        w = np.exp(x - x.max())
        sw = w.sum()
        return (w @ ex) / sw, float(sw * sw / (w @ w))

    sirt = []
    for Rl in Rl_grid:
        def f(Re):
            return np.angle(M(Rl, Re)[0]) - phm
        # R̃_e = 0'dan dışa doğru ilk kök (154'ün sahte-kök kuralı)
        kok = None
        for yon, sinir in ((+1, rmax), (-1, rmin)):
            if yon * sinir <= 0:
                continue
            Rg = np.arange(0.0, sinir + yon * 1e-9, yon * adim)
            d = np.array([f(R) for R in Rg])
            aday = None
            for i in range(len(Rg) - 1):
                if d[i] == 0.0:
                    aday = float(Rg[i]); break
                if d[i] * d[i + 1] < 0:
                    lo, hi = sorted((float(Rg[i]), float(Rg[i + 1])))
                    flo = f(lo)
                    for _ in range(50):
                        mdl = 0.5 * (lo + hi)
                        if flo * f(mdl) <= 0:
                            hi = mdl
                        else:
                            lo = mdl; flo = f(lo)
                    aday = 0.5 * (lo + hi); break
            if aday is None:
                continue
            if kok is None or abs(aday) < abs(kok):
                kok = aday
        if kok is None:
            continue
        Mv, neff = M(Rl, kok)
        sirt.append(dict(Rl=float(Rl), Re=float(kok),
                         artik=float(abs(np.angle(Mv) - phm)),
                         ReM=float(Mv.real), ImM=float(Mv.imag),
                         absM=float(abs(Mv)), neff=neff))
    return sirt


def olc3(z, etiket, bantlar=None, taban=0.52, cap=720, tohum=21,
         Rl_grid=None, ayrinti=True):
    """Üç-kanallı ağırlık ayrışımı ölçümü. z: sıralı nokta dizisi."""
    if bantlar is None:
        bantlar = BANTLAR_STD
    if Rl_grid is None:
        Rl_grid = np.round(np.arange(-1.0, 4.001, 0.125), 4)

    C = zincir3(z, taban, cap)
    mid, L, ds = C["mid"], C["L"], C["ds"]
    X = {"tam": _bond(ds), "lad": _bond(C["lad"]),
         "eta": _bond(C["eta"]), "dri": _bond(C["drift"])}
    Xt = X["tam"]
    sA2 = float(np.var(Xt))
    kov = {k: float(np.mean(v * Xt)) for k, v in X.items()}
    var = {k: float(np.var(v)) for k, v in X.items()}
    # toplama denetimi: X_tam = X_lad + X_eta + X_dri (ortalamalar çıkık)
    top_hata = float(np.max(np.abs(Xt - (X["lad"] + X["eta"] + X["dri"]))))
    s = {k: sA2 / kov[k] for k in ("tam", "lad", "eta")}
    Y = {k: X[k] * s[k] for k in ("tam", "lad", "eta")}

    print(f"[{etiket}] N={len(z)}  L={L:.4f}  regresör={C['nq']}  "
          f"ds-kapanış={C['kapanis']:.2e}  bond-kapanış={top_hata:.2e}",
          flush=True)
    print(f"  nokta varyansları : σ_ds²={C['s_ds']:.5f}  σ_lad²={C['s_lad']:.5f}"
          f"  σ_η²={C['s_eta']:.5f}  σ_drift²={C['s_dri']:.2e}  "
          f"c₁(η)={C['c1']:+.5f}", flush=True)
    print(f"  bond varyansları  : σΔ²={sA2:.5f}  Var(X_lad)={var['lad']:.5f}"
          f"  Var(X_eta)={var['eta']:.5f}  Var(X_dri)={var['dri']:.2e}",
          flush=True)
    print(f"  PAY YAPISI Cov(X_k,X_tam)/σΔ² : lad={kov['lad']/sA2:+.4f}  "
          f"eta={kov['eta']/sA2:+.4f}  dri={kov['dri']/sA2:+.4f}   "
          f"(toplam {(kov['lad']+kov['eta']+kov['dri'])/sA2:+.4f})", flush=True)
    print(f"  Cov(X_lad,X_eta)={float(np.mean(X['lad']*X['eta'])):+.6f}  "
          f"ölçekler s_lad={s['lad']:.3f} s_eta={s['eta']:.3f}", flush=True)

    qm = CEK.pk_m(int(np.exp(0.86 * L)))
    allq = sorted(qm)
    allw = np.array([np.log(q) for q in allq])
    T = mid[-1] - mid[0]; dres = TWO_PI / T
    eta = C["eta"]
    e0, e1 = eta[:-1], eta[1:]
    m0 = mid[:-1]
    rng = np.random.default_rng(tohum)

    satir = []
    for lo, hi in bantlar:
        tumu = [(q, w) for q, w in zip(allq, allw) if lo < w / L <= hi]
        cand = tumu
        if len(cand) > 220:
            idx = rng.choice(len(cand), 220, replace=False)
            cand = [cand[i] for i in idx]
        crN = 0j; crO = 0j; on0 = off0 = 0.0
        kul = 0
        for q, w in cand:
            j = np.searchsorted(allw, w)
            koms = [allw[k] for k in (j - 1, j + 1)
                    if 0 <= k < len(allw) and abs(allw[k] - w) > 1e-12]
            gap = min(abs(w - k) for k in koms)
            if gap < 2.5 * dres:
                continue
            kul += 1
            for W, hedef in ((w, True), (w + gap / 2, False)):
                cw, sw_ = np.cos(W * m0), np.sin(W * m0)
                zc = complex(2 * np.mean(e0 * cw), -2 * np.mean(e0 * sw_))
                zp = complex(2 * np.mean(e1 * cw), -2 * np.mean(e1 * sw_))
                cr = zp * np.conj(zc) * np.exp(1j * TWO_PI * W / L)
                if hedef:
                    crN += cr; on0 += abs(zc)**2
                else:
                    crO += cr; off0 += abs(zc)**2

        tb = 0.5 * (lo + hi)
        A = TWO_PI * tb
        ex = np.exp(-1j * A * Xt)
        Me = complex(np.mean(ex))
        payda = on0 - off0
        if kul == 0 or abs(payda) < 1e-300:
            satir.append(dict(tau=round(tb, 4), lo=lo, hi=hi, kul=kul,
                              olculdu=False))
            print(f"  τ={tb:.4f}  ÖLÇÜLEMEDİ (aday {kul}, payda {payda:.1e})",
                  flush=True)
            continue
        G = (crN - crO) / payda
        phm = float(np.angle(G))
        patlak = bool(abs(G) > 1.5)

        rec = dict(tau=round(tb, 4), lo=lo, hi=hi, kul=kul, olculdu=True,
                   L=L, sA2=sA2, kov={k: kov[k] for k in kov},
                   var={k: var[k] for k in var}, olcek=s,
                   Gre=float(G.real), Gim=float(G.imag), absG=float(abs(G)),
                   phi=phm, absMe=float(abs(Me)), patlak=patlak, kanal={})
        if ayrinti:
            print(f"  τ={tb:.4f}  aday={kul}  |Γ|={abs(G):.4f}  "
                  f"φ={phm:+.4f}  ReΓ={G.real:+.4f}  |M_emp|={abs(Me):.4f}"
                  f"{'   << |Γ|>1.5 PATLAK' if patlak else ''}", flush=True)
        # ---- üç tek-kanallı model ----
        for k in ("tam", "lad", "eta"):
            kap = _kapasite(A, Y[k], ex)
            yeter = bool(kap["arg_min"] - 1e-9 <= phm <= kap["arg_max"] + 1e-9)
            Rt, art, M, neff = CEK._R_ile_faz_es(phm, A, Y[k], ex,
                                                 Rmin=RMIN, Rmax=RMAX,
                                                 adim=RADIM)
            dRe = abs(M.real - G.real)
            rec["kanal"][k] = dict(Rn=Rt, Rham=Rt * s[k], artik=art,
                                   ReM=float(M.real), ImM=float(M.imag),
                                   absM=float(abs(M)), dRe=float(dRe),
                                   neff=neff, kapasite=kap, faz_yeter=yeter,
                                   korel=float(kov[k] / np.sqrt(var[k] * sA2)))
            if ayrinti:
                bay = ("" if art < 0.02 else
                       f" << FAZ ÜRETİLEMİYOR (menzil {kap['arg_min']:+.3f}"
                       f"..{kap['arg_max']:+.3f})")
                print(f"      {k:3s}: R̃={Rt:+7.3f}  R_ham={Rt*s[k]:+9.3f}  "
                      f"artık={art:.4f}  ReM={M.real:+.4f}  "
                      f"|ΔRe|={dRe:.4f}  n_eff={neff:9.0f}  ρ_kanal="
                      f"{kov[k]/np.sqrt(var[k]*sA2):+.3f}"
                      f"  {'✓' if (art < 0.02 and dRe < 0.06) else '✗'}{bay}",
                      flush=True)
        # ---- iki kanallı sırt ----
        sirt = _sirt(A, Y["lad"], Y["eta"], ex, phm, Rl_grid)
        for p in sirt:
            p["Rl_ham"] = p["Rl"] * s["lad"]
            p["Re_ham"] = p["Re"] * s["eta"]
            p["dRe"] = abs(p["ReM"] - G.real)
        sirt_ok = [p for p in sirt if p["artik"] < 0.02]
        if sirt_ok:
            en = min(sirt_ok, key=lambda p: abs(p["ReM"] - G.real))
            rec["sirt"] = sirt
            rec["en_iyi2"] = dict(en, dRe=abs(en["ReM"] - G.real),
                                  Rl_ham=en["Rl"] * s["lad"],
                                  Re_ham=en["Re"] * s["eta"])
            if ayrinti:
                print(f"      2K : R̃_l={en['Rl']:+.3f} R̃_e={en['Re']:+.3f}"
                      f"  (ham {en['Rl']*s['lad']:+.3f} / "
                      f"{en['Re']*s['eta']:+.3f})  ReM={en['ReM']:+.4f}  "
                      f"|ΔRe|={abs(en['ReM']-G.real):.4f}  "
                      f"n_eff={en['neff']:.0f}  sırt-nokta={len(sirt_ok)}",
                      flush=True)
        else:
            rec["sirt"] = sirt
            rec["en_iyi2"] = None
            if ayrinti:
                print("      2K : sırt boş (hiçbir R̃_l'de faz kökü yok)",
                      flush=True)
        satir.append(rec)

    return dict(etiket=etiket, L=L, sA2=sA2, s_ds=C["s_ds"], s_lad=C["s_lad"],
                s_eta=C["s_eta"], s_dri=C["s_dri"], c1=C["c1"],
                kapanis=C["kapanis"], bond_kapanis=top_hata,
                kov=kov, var=var, olcek=s, bantlar=satir)
