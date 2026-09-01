"""
157 — R_lad(τ) ÖLÇÜM ÇEKİRDEĞİ (taban-bağımsız nesne adayı)
==========================================================
156, Γ_rot'un fazını üreten sağkalım-ağırlık bağlaşımının MERDİVEN
payında yaşadığını buldu ve R_lad'ın regresyon TABANINDAN neredeyse
bağımsız (±%3) olduğunu ölçtü. 154'ün kapalı-form yarışı ise taban ile
%6–27 kayan R_tam üzerinde koşulmuştu. Bu çekirdek yarışı R_lad üstünde
tekrarlamak için gerekli tayfı üretir.

ÜÇ RAPORUN ZİNCİRİ TEK YERDE (kopya değil, IMPORT):
   154_cekirdek : pk, pk_m, _M, _R_150grid          (ölçüm/kök çekirdeği)
   156_cekirdek : zincir3 (drift/lad/eta ayrışımı), _bond
Bu dosyanın kendi yazdığı tek şey BANT DÖNGÜSÜdür; o da 155'in
çizgi-bazında kayıt + τ_eff apsisi eklentisiyle (F1/F2) genişletilmiştir.

NE ÖLÇÜLÜR
----------
Her bantta, 150/154/156'nın modeli:
    Γ_rot(τ) = ⟨w_k·e^{−iA·X_tam}⟩ / ⟨w_k⟩ ,  w_k = e^{−A·R_k·X_k}
FAZ FAKTÖRÜ HER ZAMAN TAM ADIMDIR; yalnız AĞIRLIK kanalı k değişir
(k ∈ {tam, lad, eta}). R_k ölçülen fazı eşleyecek biçimde çözülür; Re
kısmı bağımsız bir AŞIRI-BELİRLEME sınavı olarak kalır (|ΔRe| < 0.06).

156'nın ölçek normalizasyonu birebir korunur:
    s_k = σΔ² / Cov(X_k, X_tam),   Y_k = X_k·s_k,   R_k(ham) = R̃_k·s_k
R̃ üç kanalda aynı aralıkta aranır; ham R ise doğrudan X_k'nin
katsayısıdır (w = e^{−A·R_ham·X_k}) ve gazlar/tabanlar arası
karşılaştırmalar HAM R üzerinden yapılır.

157'NİN EKLERİ
--------------
 (G1) τ_eff apsisi (155/F2): τ_eff = Σ τ_i·(pow_on−pow_off)_i /
      Σ(pow_on−pow_off)_i. 155, bant-ortası apsisin τ₀'a SAHTE bir
      L-bağımlılığı imal ettiğini ölçtü; bu koşuda apsis her yerde τ_eff.
 (G2) arg M_k(R) IZGARASI JSON'a yazılır. Böylece "φ-ilkel" adayı (h)
      — φ(τ) modellenip R'nin ondan TÜRETİLMESİ — ölçümü yeniden
      koşmadan, aynı ampirik ağırlık dağılımıyla sınanabilir.
 (G3) Jackknife φ, R_tam, R_lad, τ_eff, ρ için EŞZAMANLI (aynı 8 grup,
      bütün bantlarda aynı k indisi) — τ₀ fiti replika-replika
      tekrarlanabilsin diye (155'in hata mimarisi).

KİMLİK UYARISI (raporun anahtar noktası, koda not olarak da düşülüyor):
R_k = 0 ⟺ M_k(0) = ⟨e^{−iA·X_tam}⟩ = M_emp, ve M_emp KANALDAN
BAĞIMSIZDIR. Dolayısıyla R_lad'ın sıfır-geçişi ile R_tam'ınki AYNI τ'da
olmak zorundadır (φ_Γ = arg M_emp koşulu). Bu çekirdek her iki kanalın
sıfırını da ölçer ki bu kimlik SAYIYLA doğrulansın (ya da çürütülsün).
"""

import importlib
import sys
from pathlib import Path

import numpy as np

_H = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_H / "154_configs"))
sys.path.insert(0, str(_H / "156_configs"))
C154 = importlib.import_module("154_cekirdek")
C156 = importlib.import_module("156_cekirdek")

TWO_PI = 2 * np.pi
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/157")

KANALLAR = ("tam", "lad", "eta")


def izgara(lo, hi, w):
    n = int(round((hi - lo) / w))
    return [(round(lo + w * i, 4), round(lo + w * (i + 1), 4))
            for i in range(n)]


# 0.02 ızgara: 155'in `ince`(0.42–0.56) ve `yuksek`(0.52–0.68) ızgaralarını
# BİREBİR kapsar → kopya-kayması denetimi doğrudan yapılabilir.
INCE = izgara(0.42, 0.80, 0.02)            # 19 bant
# 0.03 ızgara, 0.43'e demirli: kenarları hem 0.46'ya hem 0.52'ye DÜŞER,
# yani her iki taban da ızgarayı kaydırmadan koşulabilir. (Görevin
# istediği 8–10 bantlık tayf; taban 0.40'ta 12, 0.46'da 11, 0.52'de 9.)
KABA = izgara(0.43, 0.79, 0.03)            # 12 bant
# yalnız sıfır-geçişi bölgesi (τ₀ fiti için; INCE'nin alt parçasıyla BİREBİR
# aynı bantlar — pahalı yüksek-τ bantlarını koşmadan)
TAU0 = izgara(0.42, 0.62, 0.02)            # 10 bant

# R̃ arama ızgarası. 156: (−8, 20); 155: (−6, 8). Burada negatif dal
# (τ<τ₀) DA istendiği için alt uç açık; ızgara 0'ı tam üstünde taşır.
RG_MIN, RG_MAX, RG_ADIM = -10.0, 14.0, 0.05
# JSON'a yazılan (h)-adayı ızgarası (daha dar, daha yoğun)
KAY_MIN, KAY_MAX, KAY_ADIM = -8.0, 10.0, 0.05


# ---------------------------------------------------------------- önbellek
def kanal_onbellek(z, anahtar, taban, cap=4000):
    """156.zincir3'ün drift/lad/eta ayrışımı, diske önbelleklenmiş.

    Regresyon (300k × ~250 sütun) tek pahalı adım ve aynı (pencere,taban)
    için birkaç bant ızgarası koşulacak.
    """
    SCR.mkdir(parents=True, exist_ok=True)
    yol = SCR / f"kanal_{anahtar}_t{taban}_c{cap}.npz"
    if yol.exists():
        d = np.load(yol)
        print(f"  [kanal önbellek] {yol.name}", flush=True)
        return {k: (d[k] if d[k].ndim else float(d[k])) for k in d.files}
    C = C156.zincir3(z, taban, cap)
    np.savez_compressed(
        yol, mid=C["mid"], L=C["L"], ds=C["ds"], drift=C["drift"],
        lad=C["lad"], eta=C["eta"], nq=C["nq"], kapanis=C["kapanis"],
        s_ds=C["s_ds"], s_eta=C["s_eta"], s_lad=C["s_lad"],
        s_dri=C["s_dri"], c1=C["c1"])
    return C


# --------------------------------------------------------------- R çözücü
class Cozucu:
    """arg M_k(R̃) ızgarasını BİR KEZ tarar; birçok φ için kök verir.

    Kök mantığı birebir 154._R_ile_faz_es: R=0'DAN DIŞA doğru aranır ve
    ilk (|R|'si en küçük) kök alınır — arg M büyük |R|'de monoton
    olmadığından sahte kökler ancak böyle elenir.
    """

    def __init__(self, A, Y, ex, rmin=RG_MIN, rmax=RG_MAX, adim=RG_ADIM):
        self.A, self.Y, self.ex = A, Y, ex
        n0 = int(round(-rmin / adim))
        n1 = int(round(rmax / adim))
        self.Rg = np.round(adim * np.arange(-n0, n1 + 1), 8)
        self.i0 = n0                                # Rg[i0] == 0.0
        ang = np.empty(len(self.Rg)); re = np.empty(len(self.Rg))
        ne = np.empty(len(self.Rg))
        for i, R in enumerate(self.Rg):
            M, n = C154._M(R, A, Y, ex)
            ang[i] = np.angle(M); re[i] = M.real; ne[i] = n
        self.ang, self.re, self.ne = ang, re, ne

    def _f(self, R, phm):
        return np.angle(C154._M(R, self.A, self.Y, self.ex)[0]) - phm

    def _bisek(self, lo, hi, phm):
        flo = self._f(lo, phm)
        for _ in range(50):
            mdl = 0.5 * (lo + hi)
            if flo * self._f(mdl, phm) <= 0:
                hi = mdl
            else:
                lo = mdl; flo = self._f(lo, phm)
        return 0.5 * (lo + hi)

    def coz(self, phm):
        kok = None
        for yon in (+1, -1):
            idx = (np.arange(self.i0, len(self.Rg)) if yon > 0
                   else np.arange(self.i0, -1, -1))
            d = self.ang[idx] - phm
            aday = None
            for i in range(len(idx) - 1):
                if d[i] == 0.0:
                    aday = float(self.Rg[idx[i]]); break
                if d[i] * d[i + 1] < 0:
                    lo, hi = sorted((float(self.Rg[idx[i]]),
                                     float(self.Rg[idx[i + 1]])))
                    aday = self._bisek(lo, hi, phm); break
            if aday is None:
                continue
            if kok is None or abs(aday) < abs(kok):
                kok = aday
        tavan = kok is None
        if tavan:                       # kök yok → en yakın ızgara noktası
            kok = float(self.Rg[int(np.argmin(np.abs(self.ang - phm)))])
        M, neff = C154._M(kok, self.A, self.Y, self.ex)
        return (float(kok), float(abs(np.angle(M) - phm)), M, neff, tavan)

    def kayit(self):
        """(h)-adayı için: arg M / Re M / n_eff ızgarası (daha dar aralık)."""
        m = (self.Rg >= KAY_MIN - 1e-9) & (self.Rg <= KAY_MAX + 1e-9)
        return dict(R=[float(v) for v in self.Rg[m]],
                    arg=[float(v) for v in self.ang[m]],
                    re=[float(v) for v in self.re[m]],
                    neff=[float(v) for v in self.ne[m]])


# ----------------------------------------------------------------- ölçüm
def olc157(z, etiket, bantlar, anahtar, taban=0.52, cap=4000, njack=8,
           tohum=21, ayrinti=True, kanallar=KANALLAR, kayit_izgara=True):
    """R_tam / R_lad / R_eta tayfı, τ_eff apsisi ve çizgi-bazında kayıt."""
    C = kanal_onbellek(z, anahtar, taban, cap)
    mid, L, ds = C["mid"], float(C["L"]), C["ds"]
    eta = C["eta"]
    X = {"tam": C156._bond(ds), "lad": C156._bond(C["lad"]),
         "eta": C156._bond(C["eta"]), "dri": C156._bond(C["drift"])}
    Xt = X["tam"]
    sA2 = float(np.var(Xt))
    kov = {k: float(np.mean(v * Xt)) for k, v in X.items()}
    var = {k: float(np.var(v)) for k, v in X.items()}
    top_hata = float(np.max(np.abs(Xt - (X["lad"] + X["eta"] + X["dri"]))))
    s = {k: sA2 / kov[k] for k in KANALLAR}
    Y = {k: X[k] * s[k] for k in KANALLAR}
    korel = {k: float(kov[k] / np.sqrt(var[k] * sA2)) for k in KANALLAR}

    m2 = sA2
    m3 = float(np.mean(Xt**3)); m4 = float(np.mean(Xt**4))
    mom = dict(sA2=sA2, k3=m3, k4=m4 - 3 * m2 * m2, carp=m3 / m2**1.5,
               sds2=float(np.var(ds)))

    print(f"[{etiket}] N={len(z)} L={L:.4f} regresör={int(C['nq'])} "
          f"σΔ²={sA2:.5f} σ_ds²={float(C['s_ds']):.5f} "
          f"σ_lad²={float(C['s_lad']):.5f} σ_η²={float(C['s_eta']):.5f} "
          f"c₁={float(C['c1']):+.5f}", flush=True)
    print(f"  ds-kapanış={float(C['kapanis']):.2e} bond-kapanış={top_hata:.2e}"
          f"  PAY lad={kov['lad']/sA2:+.4f} eta={kov['eta']/sA2:+.4f} "
          f"dri={kov['dri']/sA2:+.4f}  ölçek s_lad={s['lad']:.4f} "
          f"s_eta={s['eta']:.3f}  korel_lad={korel['lad']:+.4f}", flush=True)

    qm = C154.pk_m(int(np.exp(0.86 * L)))
    allq = sorted(qm)
    allw = np.array([np.log(q) for q in allq])
    T = mid[-1] - mid[0]
    dres = TWO_PI / T
    e0, e1 = eta[:-1], eta[1:]
    m0 = mid[:-1]
    rng = np.random.default_rng(tohum)

    satir = []
    if ayrinti:
        print("   τ̄     τ_eff   N  kul   |Γ|    φ_Γ     argMe   "
              "R_tam   R_lad(ham)  ±jk   artık  |ΔRe|   n_eff", flush=True)
    for lo, hi in bantlar:
        tumu = [(q, w) for q, w in zip(allq, allw) if lo < w / L <= hi]
        cand = tumu
        if len(cand) > 220:
            idx = rng.choice(len(cand), 220, replace=False)
            cand = [cand[i] for i in idx]
        Ls = []
        kul = 0
        for q, w in cand:
            j = np.searchsorted(allw, w)
            koms = [allw[k] for k in (j - 1, j + 1)
                    if 0 <= k < len(allw) and abs(allw[k] - w) > 1e-12]
            gap = min(abs(w - k) for k in koms)
            if gap < 2.5 * dres:
                continue
            gi = kul % njack
            kul += 1
            aq = 1.0 / (np.pi * qm[q] * np.sqrt(q))
            rec = dict(q=int(q), w=float(w), tau=float(w / L), grup=gi,
                       gp=float((2 * aq * np.sin(np.pi * w / L))**2))
            for W, hedef in ((w, True), (w + gap / 2, False)):
                cw, sw = np.cos(W * m0), np.sin(W * m0)
                zc = complex(2 * np.mean(e0 * cw), -2 * np.mean(e0 * sw))
                zp = complex(2 * np.mean(e1 * cw), -2 * np.mean(e1 * sw))
                cr = zp * np.conj(zc) * np.exp(1j * TWO_PI * W / L)
                et = "on" if hedef else "off"
                rec[f"cr_{et}_re"] = float(cr.real)
                rec[f"cr_{et}_im"] = float(cr.imag)
                rec[f"pow_{et}"] = float(abs(zc)**2)
            Ls.append(rec)

        tb = 0.5 * (lo + hi)
        A = TWO_PI * tb
        ex = np.exp(-1j * A * Xt)
        Me = complex(np.mean(ex))
        if kul == 0:
            satir.append(dict(tau=round(tb, 4), lo=lo, hi=hi, N=len(tumu),
                              kul=0, olculdu=False))
            continue

        crN = np.array([complex(r["cr_on_re"], r["cr_on_im"]) for r in Ls])
        crO = np.array([complex(r["cr_off_re"], r["cr_off_im"]) for r in Ls])
        pN = np.array([r["pow_on"] for r in Ls])
        pO = np.array([r["pow_off"] for r in Ls])
        gpv = np.array([r["gp"] for r in Ls])
        tv = np.array([r["tau"] for r in Ls])
        grp = np.array([r["grup"] for r in Ls])
        u = pN - pO
        payda = float(u.sum())
        if abs(payda) < 1e-300:
            satir.append(dict(tau=round(tb, 4), lo=lo, hi=hi, N=len(tumu),
                              kul=kul, olculdu=False))
            continue
        G = complex((crN - crO).sum() / payda)
        phm = float(np.angle(G))
        rho = payda / float(gpv.sum())
        tau_eff = float((tv * u).sum() / payda)
        patlak = bool(abs(G) > 1.5)

        rec = dict(tau=round(tb, 4), lo=lo, hi=hi, N=len(tumu), kul=kul,
                   olculdu=True, L=L, sA2=sA2, taban=taban,
                   Gre=float(G.real), Gim=float(G.imag), absG=float(abs(G)),
                   phi=phm, absMe=float(abs(Me)), argMe=float(np.angle(Me)),
                   rho=rho, tau_eff=tau_eff, tau_ari=float(tv.mean()),
                   patlak=patlak, kov=kov, var=var, olcek=s, korel=korel,
                   kanal={}, cizgi=Ls)

        nanv = float("nan")
        phj = [nanv] * njack; tej = [nanv] * njack; rhoj = [nanv] * njack
        phi_k = {}
        for k in range(njack):
            m = grp != k
            if not m.any():
                continue
            pk_ = float(u[m].sum())
            if abs(pk_) < 1e-300:
                continue
            Gk = complex((crN[m] - crO[m]).sum() / pk_)
            if abs(Gk) > 1.5:
                continue
            phj[k] = float(np.angle(Gk))
            rhoj[k] = pk_ / float(gpv[m].sum())
            tej[k] = float((tv[m] * u[m]).sum() / pk_)
            phi_k[k] = phj[k]

        def jkerr(v):
            v = np.asarray(v, float); v = v[np.isfinite(v)]
            if len(v) < 4:
                return float("nan")
            return float(np.sqrt((len(v) - 1) / len(v)
                                 * np.sum((v - v.mean())**2)))

        for kan in kanallar:
            cz = Cozucu(A, Y[kan], ex)
            Rt, art, M, neff, tavan = cz.coz(phm)
            Rj = [nanv] * njack
            for k, pk_ang in phi_k.items():
                rk, ak, _, _, tv_ = cz.coz(pk_ang)
                if ak < 0.02 and not tv_:
                    Rj[k] = rk
            d = dict(Rn=Rt, Rham=Rt * s[kan], artik=art,
                     ReM=float(M.real), ImM=float(M.imag),
                     absM=float(abs(M)), dRe=float(abs(M.real - G.real)),
                     neff=neff, tavan=bool(tavan), korel=korel[kan],
                     olcek=s[kan], R_jk=Rj,
                     sR_jk=jkerr(Rj), sRham_jk=jkerr(Rj) * abs(s[kan]),
                     njk=int(np.isfinite(np.array(Rj, float)).sum()))
            if kayit_izgara and kan in ("tam", "lad"):
                d["izgara"] = cz.kayit()
            rec["kanal"][kan] = d

        rec.update(sPhi_jk=jkerr(phj), sTeff_jk=jkerr(tej),
                   sRho_jk=jkerr(rhoj), phi_jk=phj, teff_jk=tej,
                   rho_jk=rhoj)
        satir.append(rec)
        if ayrinti:
            kt, kl = rec["kanal"]["tam"], rec["kanal"]["lad"]
            bay = "" if kl["artik"] < 0.02 else " << LAD FAZ TUTMADI"
            pat = " << |Γ|>1.5" if patlak else ""
            print(f"  {tb:.4f} {tau_eff:.4f} {len(tumu):4d} {kul:3d} "
                  f"{abs(G):.3f} {phm:+.4f} {np.angle(Me):+.4f} "
                  f"{kt['Rham']:+7.3f} {kl['Rham']:+9.3f} "
                  f"{kl['sRham_jk']:6.3f} {kl['artik']:6.4f} "
                  f"{kl['dRe']:6.4f} {kl['neff']:8.0f}{bay}{pat}", flush=True)

    return dict(etiket=etiket, L=L, sA2=sA2, taban=taban, cap=cap,
                s_ds=float(C["s_ds"]), s_eta=float(C["s_eta"]),
                s_lad=float(C["s_lad"]), s_dri=float(C["s_dri"]),
                c1=float(C["c1"]), nq=int(C["nq"]),
                kapanis=float(C["kapanis"]), bond_kapanis=top_hata,
                kov=kov, var=var, olcek=s, korel=korel, mom=mom,
                bantlar=satir)
