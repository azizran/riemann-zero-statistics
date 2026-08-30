"""
154 — R(τ) ÖLÇÜM ÇEKİRDEĞİ (tek kod yolu)
=========================================
150_gamma_insasi.py'nin zinciri BİREBİR buraya taşındı ve tek bir
fonksiyona (`olc`) kapatıldı; gerçek pencereler DE sentetik gazlar DA
aynı kodu çağırır. Amaç: 152 raporunun uyardığı "kopya kayması"
riskini sıfırlamak — gerçek/sentetik farkı yalnız GİRDİ z dizisinden
gelsin.

150'ye göre KASITLI ve tek tek işaretli eklemeler (ölçümü değiştirmez,
yalnız hata/denetim üretir):
  (E1) R taraması 150'nin linspace(0,3,121) ızgarasıyla AYNEN de
       yapılır (R_150) — tekrar-üretim kontrolü. Ayrıca tavanı 8'e
       çekilmiş ince tarama + bisection (R_ince) eklenir; sentetik
       gazda 150'nin 3.0 tavanına yapışan bantlar için gerekli.
  (E2) Faz artığı |arg M(R) − arg Γ| kaydedilir. > 0.02 rad ise R
       eşleşmesi TUTMAMIŞTIR (faz sarması / tavan) — o R anlamsız.
  (E3) Jackknife: aday çizgiler K gruba bölünür, birer grup dışarıda
       bırakılarak R yeniden ölçülür → bant-içi istatistiksel hata.
  (E4) n_eff = (Σw)²/Σw² kaydedilir: ağırlık birkaç uç noktaya
       çökmüşse ortalama anlamsızdır.
  (E5) exp taşmasına karşı w = exp(x − max x) (oran alındığı için
       sabiti sadeleşir; sayısal olarak özdeş).
  (E6) AYNI aday çizgilerden ρ(τ) da ölçülür (144'ün tanımı birebir:
       ρ = (Σon − Σoff)/Σ(2a·sin(πω/L))², a = 1/(πm√q)). Γ döngüsü
       zaten on/off biriktiriyor; tek ek çıplak-güç paydası. Böylece
       R(τ) ve ρ(τ) AYNI bantlardan, AYNI koşudan çıkar — "R =
       −dlnρ/dA" adayı ölçülmüş ρ ile SIFIR parametreyle sınanabilir.
"""

import numpy as np
from sympy import primerange

TWO_PI = 2 * np.pi

BANTLAR_STD = [(0.525, 0.55), (0.55, 0.62), (0.62, 0.70),
               (0.70, 0.78), (0.78, 0.85)]
# şekil taraması için daha ince (ve daha gürültülü) ızgara
BANTLAR_INCE = [(0.525, 0.555), (0.555, 0.585), (0.585, 0.615),
                (0.615, 0.645), (0.645, 0.675), (0.675, 0.705),
                (0.705, 0.735), (0.735, 0.765), (0.765, 0.795)]


def pk(lim):
    out = []
    for p in primerange(2, lim + 1):
        q = p
        while q <= lim:
            out.append(q); q *= p
    return sorted(set(out))


def pk_m(lim):
    """pk ile AYNI küme, ek olarak kule basamağı m (q = p^m)."""
    out = {}
    for p in primerange(2, lim + 1):
        q, m = p, 1
        while q <= lim:
            out[q] = m
            q *= p; m += 1
    return out


def eta_zinciri(z, taban=0.52, cap=720):
    """150'nin ds → 0.52-taban regresyonu → η zinciri (birebir).

    taban/cap YALNIZ taban-taraması için parametredir; varsayılan
    (0.52, 720) 150/149/144'ün standart konvansiyonudur ve o değerlerde
    kod birebir 150'dir. NOT: taban>0.55'te exp(taban·L) 720'yi aşar,
    yani standart cap sessizce devreye girer — tarama için cap
    yükseltilmelidir (aksi halde "taban" değişmiş sayılmaz).
    """
    g = np.diff(z)
    mid = 0.5 * (z[:-1] + z[1:])
    Lw = np.log(mid / TWO_PI)
    L = float(Lw.mean())
    ds = g * Lw / TWO_PI - 1
    Nn = len(ds)
    tt = (mid - mid.mean()) / (mid[-1] - mid[0])

    qs = pk(min(int(np.exp(taban * L)), cap))
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
    return dict(mid=mid, L=L, ds=ds, eta=eta, nq=len(qs),
                s_ds=float(np.var(ds)), s_eta=float(np.var(eta)),
                c1=float(np.mean(eta[:-1] * eta[1:])))


def _M(R, A, dsA, ex):
    """⟨w·e^{−iA·dsΔ}⟩/⟨w⟩ (E5 kaydırmalı), n_eff ile."""
    x = -A * R * dsA
    w = np.exp(x - x.max())
    sw = w.sum()
    return (w @ ex) / sw, float(sw * sw / (w @ w))


def _R_ile_faz_es(phm, A, dsA, ex, Rmin=-3.0, Rmax=8.0, adim=0.02):
    """arg M(R) = phm'i çözen R (kaba tarama + bisection). Artığı da döner.

    Taban 0.46 taramasında τ<0.52 bantlarında φ_Γ NEGATİF çıkıyor; bunun
    kökü NEGATİF R'dedir (ağırlık ters korelasyonlu).

    DİKKAT: arg M(R) BÜYÜK |R|'de monoton DEĞİL — ağırlık dağılımın bir
    kuyruğuna çökünce faz sarıyor ve SAHTE kökler doğuyor (ilk denemede
    τ=0.69 bandı −2.87 gibi anlamsız bir köke yapıştı, artık 1.56).
    Bu yüzden kök R=0'DAN DIŞA DOĞRU aranır ve ilk rastlanan (|R|'si en
    küçük) kök alınır. Rmin=0 ile davranış 150'nin özgün taramasıyla
    özdeştir; negatif dal yalnız φ<arg M_emp olan bantlarda devreye girer.
    """
    def f(R):
        return np.angle(_M(R, A, dsA, ex)[0]) - phm

    def bisek(lo, hi):
        flo = f(lo)
        for _ in range(60):
            mdl = 0.5 * (lo + hi)
            if flo * f(mdl) <= 0:
                hi = mdl
            else:
                lo = mdl; flo = f(lo)
        return 0.5 * (lo + hi)

    kok = None
    for yon, sinir in ((+1, Rmax), (-1, Rmin)):
        if yon * sinir <= 0:
            continue
        Rg = np.arange(0.0, sinir + yon * 1e-9, yon * adim)
        d = np.array([f(R) for R in Rg])
        for i in range(len(Rg) - 1):
            if d[i] == 0.0:
                aday = float(Rg[i]); break
            if d[i] * d[i + 1] < 0:
                aday = bisek(*sorted((Rg[i], Rg[i + 1]))); break
        else:
            continue
        if kok is None or abs(aday) < abs(kok):
            kok = aday
    if kok is None:                       # kök yok → en yakın ızgara noktası
        Rg = np.arange(Rmin, Rmax + 1e-9, adim)
        d = np.array([f(R) for R in Rg])
        kok = float(Rg[int(np.argmin(np.abs(d)))])
    M, neff = _M(kok, A, dsA, ex)
    return float(kok), float(abs(np.angle(M) - phm)), M, neff


def _R_150grid(phm, A, dsA, ex):
    """150'nin AYNEN ızgarası: linspace(0,3,121), en yakın faz."""
    best = None
    for R in np.linspace(0.0, 3.0, 121):
        M, _ = _M(R, A, dsA, ex)
        d = abs(np.angle(M) - phm)
        if best is None or d < best[0]:
            best = (d, float(R), M)
    return best[1], best[0], best[2]


def olc(z, etiket, bantlar=None, njack=8, tohum=21, ayrinti=True,
        taban=0.52, cap=720):
    """150'nin tam ölçümü + E1..E6. z: sıralı sıfır/nokta dizisi."""
    if bantlar is None:
        bantlar = BANTLAR_STD
    Ç = eta_zinciri(z, taban, cap)
    mid, L, ds, eta = Ç["mid"], Ç["L"], Ç["ds"], Ç["eta"]
    dsA = 0.5 * (ds[:-1] + ds[1:]); dsA = dsA - dsA.mean()
    sA2 = float(np.var(dsA))
    print(f"[{etiket}] N={len(z)}  L={L:.4f}  σΔ²={sA2:.4f}  "
          f"σ_ds²={Ç['s_ds']:.4f}  σ_η²={Ç['s_eta']:.4f}  "
          f"c₁={Ç['c1']:+.5f}  (regresör {Ç['nq']})", flush=True)

    qm = pk_m(int(np.exp(0.86 * L)))
    allq = sorted(qm)
    allw = np.array([np.log(q) for q in allq])
    T = mid[-1] - mid[0]; dres = TWO_PI / T
    e0, e1 = eta[:-1], eta[1:]
    m0 = mid[:-1]
    rng = np.random.default_rng(tohum)

    satir = []
    if ayrinti:
        print("  bant τ̄    N   kul   |Γ|   |M_emp|   φ_Γ    R_150  R_ince "
              " ±jk    artık  n_eff   Re_Γ  Re_M", flush=True)
    for lo, hi in bantlar:
        tumu = [(q, w) for q, w in zip(allq, allw) if lo < w / L <= hi]
        cand = tumu
        if len(cand) > 220:
            idx = rng.choice(len(cand), 220, replace=False)
            cand = [cand[i] for i in idx]
        # jackknife grupları (round-robin, kullanılan adaylar üzerinden)
        gN = np.zeros(njack, dtype=complex); gO = np.zeros(njack, dtype=complex)
        gn = np.zeros(njack); go = np.zeros(njack); gp = np.zeros(njack)
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
            # (E6) çıplak (soğurmasız) çizgi gücü — 144'ün paydası
            aq = 1.0 / (np.pi * qm[q] * np.sqrt(q))
            gp[gi] += (2 * aq * np.sin(np.pi * w / L))**2
            for W, hedef in ((w, True), (w + gap / 2, False)):
                cw, sw = np.cos(W * m0), np.sin(W * m0)
                zc = complex(2 * np.mean(e0 * cw), -2 * np.mean(e0 * sw))
                zp = complex(2 * np.mean(e1 * cw), -2 * np.mean(e1 * sw))
                cr = zp * np.conj(zc) * np.exp(1j * TWO_PI * W / L)
                if hedef:
                    gN[gi] += cr; gn[gi] += abs(zc)**2
                else:
                    gO[gi] += cr; go[gi] += abs(zc)**2

        tb = 0.5 * (lo + hi)
        A = TWO_PI * tb
        ex = np.exp(-1j * A * dsA)
        Me = complex(np.mean(ex))
        payda = (gn.sum() - go.sum())
        if kul == 0 or abs(payda) < 1e-300:
            satir.append(dict(tau=round(tb, 4), lo=lo, hi=hi, N=len(tumu),
                              kul=kul, olculdu=False))
            if ayrinti:
                print(f"  {tb:.4f}  {len(tumu):4d} {kul:3d}   — (ölçülemedi)",
                      flush=True)
            continue
        G = (gN.sum() - gO.sum()) / payda
        phm = float(np.angle(G))
        rho = float(payda / gp.sum()) if gp.sum() > 0 else float("nan")
        rj = [((gn.sum() - gn[k]) - (go.sum() - go[k])) / (gp.sum() - gp[k])
              for k in range(njack) if gp.sum() - gp[k] > 0]
        srho = (float(np.sqrt((len(rj) - 1) / len(rj) *
                              np.sum((np.array(rj) - np.mean(rj))**2)))
                if len(rj) >= 4 else float("nan"))
        R150, art150, M150 = _R_150grid(phm, A, dsA, ex)
        Rin, artin, Min, neff = _R_ile_faz_es(phm, A, dsA, ex)

        # jackknife
        Rj = []
        if kul >= njack:
            for k in range(njack):
                pk_ = (gn.sum() - gn[k]) - (go.sum() - go[k])
                if abs(pk_) < 1e-300:
                    continue
                Gk = ((gN.sum() - gN[k]) - (gO.sum() - gO[k])) / pk_
                if abs(Gk) > 1.5:
                    continue
                rk, ak, _, _ = _R_ile_faz_es(float(np.angle(Gk)), A, dsA, ex)
                if ak < 0.02:
                    Rj.append(rk)
        if len(Rj) >= 4:
            Rj = np.array(Rj)
            sR = float(np.sqrt((len(Rj) - 1) / len(Rj) *
                               np.sum((Rj - Rj.mean())**2)))
        else:
            sR = float("nan")

        rec = dict(tau=round(tb, 4), lo=lo, hi=hi, N=len(tumu), kul=kul,
                   olculdu=True, L=L, sA2=sA2,
                   Gre=float(G.real), Gim=float(G.imag), absG=float(abs(G)),
                   phi=phm, absMe=float(abs(Me)),
                   R150=R150, art150=art150, Re_M150=float(M150.real),
                   R=Rin, artik=artin, Re_M=float(Min.real), neff=neff,
                   sR_jk=sR, njk=len(Rj), rho=rho, sRho_jk=srho)
        satir.append(rec)
        if ayrinti:
            bay = "" if artin < 0.02 else "  << FAZ TUTMADI"
            pat = "  << |Γ|>1.5 PATLAK" if abs(G) > 1.5 else ""
            print(f"  {tb:.4f}  {len(tumu):4d} {kul:3d}  {abs(G):.3f}  "
                  f"{abs(Me):.3f}  {phm:+.3f}  {R150:5.2f}  {Rin:6.3f} "
                  f"{sR:6.3f} {artin:7.4f} {neff:8.0f}  {G.real:+.3f} "
                  f"{Min.real:+.3f}  ρ={rho:.4f}±{srho:.4f}{bay}{pat}",
                  flush=True)
    return dict(etiket=etiket, L=L, sA2=sA2, s_ds=Ç["s_ds"],
                s_eta=Ç["s_eta"], c1=Ç["c1"], bantlar=satir)
