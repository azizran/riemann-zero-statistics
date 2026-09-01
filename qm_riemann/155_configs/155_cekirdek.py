"""
155 — τ₀ KAMPANYASI: ÖLÇÜM ÇEKİRDEĞİ
====================================
154_cekirdek.py'nin ölçüm zincirini AYNEN kullanır — kopya DEĞİL, import:

    eta_zinciri, pk_m, _M, _R_ile_faz_es, _R_150grid   ← 154'ten import

Yalnız BANT DÖNGÜSÜ yeniden yazıldı; sebebi üç ek çıktı:

  (F1) ÇİZGİ-BAZINDA kayıt. Her kullanılan çizgi için (τ_i, cr_on_i,
       cr_off_i, pow_on_i, pow_off_i, gp_i) saklanır. Böylece φ, ρ,
       τ_eff ve her türlü jackknife/yeniden-gruplama ÖLÇÜM SONRASI,
       tekrar koşmadan hesaplanabilir. (154 yalnız 8 gruba toplanmış
       hâli saklıyordu.)
  (F2) τ_eff — bandın GÜÇ-AĞIRLIKLI etkin frekansı:
           τ_eff = Σ_i τ_i·(pow_on_i − pow_off_i) / Σ_i (pow_on_i − pow_off_i)
       Γ = Σ_i cr_i / Σ_i pow_i olduğundan, φ küçükken
       arg Γ ≈ Σ pow_i·φ(τ_i)/Σ pow_i = φ(τ_eff) (φ doğrusalsa TAM).
       Bir SIFIR-GEÇİŞİ ölçümünde apsisin bant ORTASI mı τ_eff mi
       olduğu doğrudan τ₀'a yazılır: bu, 154'te ölçülmemiş bir
       sistematiktir.
  (F3) Adım momentleri: dsΔ'nın σ², κ₃, κ₄ (H-N'in alias-kayması
       modelleri için).

Ayrıca hız: bir bantta A, dsΔ, ex sabit olduğundan arg M(R) ızgarası
BİR KEZ hesaplanıp merkez + 8 jackknife kökü için yeniden kullanılır
(154 her kök için ızgarayı baştan tarıyordu). Kök bulma mantığı ve
bisection birebir 154'ünkidir; sonuçlar bit düzeyinde aynı çıkmalı —
155_dogrulama.py bunu 154'ün kayıtlı JSON'una karşı sınar.
"""

import importlib
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "154_configs"))
C154 = importlib.import_module("154_cekirdek")

TWO_PI = 2 * np.pi
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/155")

# --- bant ızgaraları -------------------------------------------------
def izgara(lo, hi, w):
    n = int(round((hi - lo) / w))
    return [(round(lo + w * i, 4), round(lo + w * (i + 1), 4))
            for i in range(n)]

# negatif dal ince haritası (görev 1): τ∈(0.42,0.56), 7 bant × 0.02
INCE_NEG = izgara(0.42, 0.56, 0.02)
# taban 0.46 ile ölçülebilir alt küme (taban altındaki çizgiler
# regresyonda ÇIKARILDIĞI için lo ≥ taban olmalı)
INCE_NEG_046 = [b for b in INCE_NEG if b[0] >= 0.46 - 1e-9]
# düşük-L penceresi seyrek çizgili: 0.03'lük ızgara da koşulur
KABA_NEG = izgara(0.42, 0.60, 0.03)
KABA_NEG_046 = [b for b in KABA_NEG if b[0] >= 0.46 - 1e-9]
# taban taraması için derin harita: φ'nin sıfırın ÇOK altındaki şekli
GENIS = izgara(0.30, 0.60, 0.02)
# 154'ün ekstrapolasyon bölgesi (0.525-0.615) tam kapsansın diye
YUKSEK = izgara(0.52, 0.68, 0.02)


# --- SAHTE TABAN: aynı sayıda regresör, ASAL OLMAYAN frekanslarda ----
def _eta_fr(z, fr):
    """154.eta_zinciri ile AYNI algoritma, ama regresör frekansları dışarıdan.

    154'ün fonksiyonu frekansları kendi üretir (fr = log q, q ≤ e^{taban·L});
    burada fr parametre — böylece 'aynı sayıda regresör ama ASAL OLMAYAN
    frekans' kontrolü kurulabilir. `dogrula_fr` bu kopyanın gerçek fr ile
    154'ü bit düzeyinde yeniden ürettiğini sınar.
    """
    g = np.diff(z)
    mid = 0.5 * (z[:-1] + z[1:])
    Lw = np.log(mid / TWO_PI)
    L = float(Lw.mean())
    ds = g * Lw / TWO_PI - 1
    Nn = len(ds)
    tt = (mid - mid.mean()) / (mid[-1] - mid[0])
    fr = np.asarray(fr, float)
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
    return dict(mid=mid, L=L, ds=ds, eta=eta, nq=len(fr),
                s_ds=float(np.var(ds)), s_eta=float(np.var(eta)),
                c1=float(np.mean(eta[:-1] * eta[1:])))


def sahte_fr(L, taban_gercek, taban_ust, cap=4000):
    """[0,taban_gercek] GERÇEK asal-kuvvet frekansları +
       (taban_gercek,taban_ust] aralığında AYNI SAYIDA ama asal olmayan
       (ardışık çizgi aralıklarının ORTA noktaları) frekanslar."""
    qs = C154.pk(min(int(np.exp(taban_ust * L)), cap))
    w = np.array([np.log(q) for q in qs])
    gercek = w[w <= taban_gercek * L]
    ust = np.where(w > taban_gercek * L)[0]
    sahte = []
    for i in ust:
        nxt = w[i + 1] if i + 1 < len(w) else w[i] + (w[i] - w[i - 1])
        sahte.append(0.5 * (w[i] + nxt))
    return np.concatenate([gercek, np.array(sahte)]), len(gercek), len(sahte)


def dogrula_fr(z, taban, cap=4000):
    """_eta_fr, gerçek frekans kümesiyle 154.eta_zinciri'ni yeniden üretir mi?"""
    a = C154.eta_zinciri(z, taban, cap)
    qs = C154.pk(min(int(np.exp(taban * a["L"])), cap))
    b = _eta_fr(z, [np.log(q) for q in qs])
    return dict(d_eta=float(np.max(np.abs(a["eta"] - b["eta"]))),
                d_seta=abs(a["s_eta"] - b["s_eta"]),
                d_c1=abs(a["c1"] - b["c1"]), nq=a["nq"])


# --- η zinciri önbelleği ---------------------------------------------
def eta_onbellek(z, anahtar, taban, cap, sahte=None):
    """154.eta_zinciri'nin AYNISI, diske önbelleklenmiş.

    Regresyon (300k × ~500 sütun) tek pahalı adım; aynı (pencere,taban)
    için birçok bant ızgarası koşacağımızdan bir kez hesaplanır.
    sahte=(taban_gercek, taban_ust) verilirse SAHTE-TABAN kontrolü.
    """
    SCR.mkdir(parents=True, exist_ok=True)
    ek = "" if sahte is None else f"_sahte{sahte[0]}-{sahte[1]}"
    yol = SCR / f"eta_{anahtar}_t{taban}_c{cap}{ek}.npz"
    if yol.exists():
        d = np.load(yol)
        print(f"  [η önbellek] {yol.name}", flush=True)
        return dict(mid=d["mid"], L=float(d["L"]), ds=d["ds"], eta=d["eta"],
                    nq=int(d["nq"]), s_ds=float(d["s_ds"]),
                    s_eta=float(d["s_eta"]), c1=float(d["c1"]))
    if sahte is not None:
        Lk = float(np.log(0.5 * (z[:-1] + z[1:]) / TWO_PI).mean())
        fr, ng, ns = sahte_fr(Lk, sahte[0], sahte[1], cap)
        print(f"  [SAHTE TABAN] {ng} gerçek (τ≤{sahte[0]}) + {ns} sahte "
              f"(τ∈({sahte[0]},{sahte[1]}], ardışık çizgi ORTA noktaları) "
              f"= {len(fr)} regresör", flush=True)
        Ç = _eta_fr(z, fr)
    else:
        Ç = C154.eta_zinciri(z, taban, cap)
    np.savez_compressed(yol, mid=Ç["mid"], L=Ç["L"], ds=Ç["ds"], eta=Ç["eta"],
                        nq=Ç["nq"], s_ds=Ç["s_ds"], s_eta=Ç["s_eta"],
                        c1=Ç["c1"])
    return Ç


# --- hızlandırılmış R kökü (mantık birebir 154._R_ile_faz_es) --------
class RCozucu:
    """arg M(R) ızgarasını bir kez tarar; birçok φ için kök verir."""

    def __init__(self, A, dsA, ex, Rmin=-6.0, Rmax=8.0, adim=0.02):
        self.A, self.dsA, self.ex = A, dsA, ex
        self.Rmin, self.Rmax, self.adim = Rmin, Rmax, adim
        self.yon = {}
        for yon, sinir in ((+1, Rmax), (-1, Rmin)):
            if yon * sinir <= 0:
                continue
            Rg = np.arange(0.0, sinir + yon * 1e-9, yon * adim)
            ang = np.array([np.angle(C154._M(R, A, dsA, ex)[0]) for R in Rg])
            self.yon[yon] = (Rg, ang)
        Rg = np.arange(Rmin, Rmax + 1e-9, adim)
        self.tam = (Rg, np.array([np.angle(C154._M(R, A, dsA, ex)[0])
                                  for R in Rg]))

    def _f(self, R, phm):
        return np.angle(C154._M(R, self.A, self.dsA, self.ex)[0]) - phm

    def _bisek(self, lo, hi, phm):
        flo = self._f(lo, phm)
        for _ in range(60):
            mdl = 0.5 * (lo + hi)
            if flo * self._f(mdl, phm) <= 0:
                hi = mdl
            else:
                lo = mdl
                flo = self._f(lo, phm)
        return 0.5 * (lo + hi)

    def coz(self, phm):
        kok = None
        for yon in (+1, -1):
            if yon not in self.yon:
                continue
            Rg, ang = self.yon[yon]
            d = ang - phm
            aday = None
            for i in range(len(Rg) - 1):
                if d[i] == 0.0:
                    aday = float(Rg[i]); break
                if d[i] * d[i + 1] < 0:
                    aday = self._bisek(*sorted((Rg[i], Rg[i + 1])), phm); break
            if aday is None:
                continue
            if kok is None or abs(aday) < abs(kok):
                kok = aday
        if kok is None:
            Rg, ang = self.tam
            kok = float(Rg[int(np.argmin(np.abs(ang - phm)))])
        M, neff = C154._M(kok, self.A, self.dsA, self.ex)
        return float(kok), float(abs(np.angle(M) - phm)), M, neff


# --- ana ölçüm --------------------------------------------------------
def olc155(z, etiket, bantlar, anahtar, taban=0.52, cap=4000, njack=8,
           tohum=21, ayrinti=True, Rmin=-6.0, sahte=None):
    """154.olc ile AYNI ölçüm; ek olarak çizgi-bazında kayıt + τ_eff.

    Rmin=-6.0 (154'te -3.0): kök R=0'dan DIŞA doğru arandığı ve İLK
    (en küçük |R|) kök alındığı için aralığı genişletmek zaten bulunan
    hiçbir kökü DEĞİŞTİREMEZ; yalnız |R|>3 gereken (çok negatif φ)
    bantlarda kök bulunmasını sağlar.
    """
    Ç = eta_onbellek(z, anahtar, taban, cap, sahte)
    mid, L, ds, eta = Ç["mid"], Ç["L"], Ç["ds"], Ç["eta"]
    dsA = 0.5 * (ds[:-1] + ds[1:])
    dsA = dsA - dsA.mean()
    sA2 = float(np.var(dsA))
    m2 = sA2
    m3 = float(np.mean(dsA**3))
    m4 = float(np.mean(dsA**4))
    momds = dict(sA2=sA2, k3=m3, k4=m4 - 3 * m2 * m2,
                 carp=m3 / m2**1.5, bas=m4 / m2**2 - 3.0,
                 sds2=float(np.var(ds)),
                 ds_k3=float(np.mean((ds - ds.mean())**3)),
                 ds_carp=float(np.mean((ds - ds.mean())**3) / np.var(ds)**1.5))
    print(f"[{etiket}] N={len(z)}  L={L:.4f}  σΔ²={sA2:.4f}  κ₃Δ={m3:+.5f} "
          f"(çarp {momds['carp']:+.3f})  σ_ds²={Ç['s_ds']:.4f}  "
          f"σ_η²={Ç['s_eta']:.4f}  c₁={Ç['c1']:+.5f}  (regresör {Ç['nq']})",
          flush=True)

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
        print("  τ̄      τ_eff    N  kul   |Γ|    φ_Γ      R      ±jk   "
              "artık   n_eff     ρ", flush=True)
    for lo, hi in bantlar:
        tumu = [(q, w) for q, w in zip(allq, allw) if lo < w / L <= hi]
        cand = tumu
        if len(cand) > 220:
            idx = rng.choice(len(cand), 220, replace=False)
            cand = [cand[i] for i in idx]
        Ls = []          # çizgi-bazında kayıt (F1)
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
            gpi = (2 * aq * np.sin(np.pi * w / L))**2
            rec = dict(q=int(q), w=float(w), tau=float(w / L), grup=gi,
                       gp=float(gpi))
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
        ex = np.exp(-1j * A * dsA)
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
        u = pN - pO                       # Γ'nın paydasının çizgi payları
        payda = float(u.sum())
        if abs(payda) < 1e-300:
            satir.append(dict(tau=round(tb, 4), lo=lo, hi=hi, N=len(tumu),
                              kul=kul, olculdu=False))
            continue
        G = complex((crN - crO).sum() / payda)
        phm = float(np.angle(G))
        rho = payda / float(gpv.sum())
        tau_eff = float((tv * u).sum() / payda)        # (F2)
        tau_on = float((tv * pN).sum() / pN.sum())
        tau_ari = float(tv.mean())

        coz = RCozucu(A, dsA, ex, Rmin=Rmin)
        Rin, artin, Min, neff = coz.coz(phm)
        R150, art150, M150 = C154._R_150grid(phm, A, dsA, ex)

        # --- jackknife (grup sil) : φ, R, ρ, τ_eff hepsi AYNI gruplardan.
        # Diziler njack UZUNLUĞUNDA; atlanan replika nan'dır — böylece
        # k indisi BANTLAR ARASINDA aynı grubu gösterir ve τ₀ fiti
        # replika-replika (bütün bantlar birlikte) tekrarlanabilir.
        nanv = float("nan")
        phj = [nanv] * njack
        Rj = [nanv] * njack
        rhoj = [nanv] * njack
        tej = [nanv] * njack
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
            pk_ang = float(np.angle(Gk))
            phj[k] = pk_ang
            rhoj[k] = pk_ / float(gpv[m].sum())
            tej[k] = float((tv[m] * u[m]).sum() / pk_)
            rk, ak, _, _ = coz.coz(pk_ang)
            if ak < 0.02:
                Rj[k] = rk

        def jkerr(v):
            v = np.asarray(v, float)
            v = v[np.isfinite(v)]
            if len(v) < 4:
                return float("nan")
            return float(np.sqrt((len(v) - 1) / len(v)
                                 * np.sum((v - v.mean())**2)))

        rec = dict(tau=round(tb, 4), lo=lo, hi=hi, N=len(tumu), kul=kul,
                   olculdu=True, L=L, sA2=sA2,
                   Gre=float(G.real), Gim=float(G.imag), absG=float(abs(G)),
                   phi=phm, sPhi_jk=jkerr(phj), absMe=float(abs(Me)),
                   argMe=float(np.angle(Me)),
                   R=Rin, artik=artin, Re_M=float(Min.real), neff=neff,
                   sR_jk=jkerr(Rj),
                   njk=int(np.isfinite(np.array(Rj, float)).sum()),
                   R150=R150, art150=art150,
                   rho=rho, sRho_jk=jkerr(rhoj),
                   tau_eff=tau_eff, sTeff_jk=jkerr(tej),
                   tau_on=tau_on, tau_ari=tau_ari,
                   phi_jk=phj, R_jk=Rj, teff_jk=tej,
                   cizgi=Ls)
        satir.append(rec)
        if ayrinti:
            bay = "" if artin < 0.02 else "  << FAZ TUTMADI"
            pat = "  << |Γ|>1.5" if abs(G) > 1.5 else ""
            print(f"  {tb:.4f} {tau_eff:.4f} {len(tumu):4d} {kul:3d} "
                  f"{abs(G):.3f} {phm:+.4f} {Rin:+7.3f} {rec['sR_jk']:6.3f} "
                  f"{artin:7.4f} {neff:8.0f} {rho:7.4f}"
                  f"  ±φ={rec['sPhi_jk']:.4f}{bay}{pat}", flush=True)
    return dict(etiket=etiket, L=L, sA2=sA2, s_ds=Ç["s_ds"], s_eta=Ç["s_eta"],
                c1=Ç["c1"], taban=taban, cap=cap, mom=momds, bantlar=satir)
