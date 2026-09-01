"""
159 — MEKANİZMA DENKLEMİ φ_Γ(τ) = A·S(τ) SINAVI: ÖLÇÜM ÇEKİRDEĞİ
================================================================
KALEM (01 Eylül):
    Γ_rot = ⟨w·e^{−iA·dsΔ}⟩/⟨w⟩ , w = yerel bant gücü
    ⇒ φ_Γ(τ) = A·S(τ),  S = Cov(P_yerel, dsΔ)/⟨P_yerel⟩,  A = 2πτ

Bu modül denklemin iki tarafını da bağımsız ölçer ve ARADAKİ HER ADIMI
ayrı ayrı kaydeder. Ölçüm zinciri kopya değil: çizgi seçimi, tohum,
boşluk filtresi, on/off ara-nokta referansı ve jackknife grupları
155_cekirdek.olc155 ile birebir aynıdır (159_analiz V1 bunu 158'in
kayıtlı φ tablosuna karşı bit düzeyinde sınar).

────────────────────────────────────────────────────────────────────
TAHMİNCİNİN CEBRİ (bu koşuda 1e−15 ile doğrulandı; 159_analiz V2)
────────────────────────────────────────────────────────────────────
Tanımlar (n = 0 … Nn−2):
    c_n = η_n·e^{−iW·m_n}
    zc  = 2⟨c_n⟩                     (kod: pow_on = |zc|²)
    zp  = 2⟨η_{n+1}·e^{−iW·m_n}⟩
    X̃_n = (m_{n+1}−m_n)·L/2π − 1     (tam etkin bond adımı)
    dsΔ_n = (ds_n+ds_{n+1})/2 − ort   (150/154/155'in adım değişkeni)

m_{n+1} − m_n = (2π/L)(1+X̃_n) olduğundan ÖZDEŞ olarak

    zp·conj(zc)·e^{−iA}/4 = ⟨(ρ_n + i·σ_n)·e^{+iA·X̃_n}⟩
    ρ_n + i·σ_n ≡ c_{n+1}·conj(⟨c⟩)        ⟨ρ⟩ = |⟨c⟩|² = pow/4

Yani HAM korelatör zaten e^{+iA} taşır; "bilinen faz ilerlemesini
çıkarmak" için e^{−iA} ile çarpmak gerekir. 148'den beri bütün zincir
(148/149/150/151/152/153/154/155/156/157) e^{+i·2πW/L} ile ÇARPIYOR:

    φ_Γ = arg Γ = (2A − 2π) + δ ,   δ = arg[zp·conj(zc)·e^{−iA}]

2A − 2π = 4π(τ−½): KİNEMATİK OMURGA (gazdan bağımsız).
δ: mekanizma fazı — KALEM'in φ_Γ'sının karşılığı olan nicelik.

────────────────────────────────────────────────────────────────────
MEKANİZMA MERDİVENİ (her basamak ayrı ölçülür)
────────────────────────────────────────────────────────────────────
    δ        = arg⟨(ρ+iσ)e^{iAX}⟩                       (TAM, ölçülen)
    M1       = arg[⟨ρ·e^{iAX}⟩/⟨ρ⟩]      σ kanalı atılır (mekanizma,
                                          bütün mertebeler)
    A·S      = A·Cov(ρ,X)/⟨ρ⟩            + birinci mertebe (KALEM)

YEREL GÜÇ TANIMI — gerekçe. Görev "u=η cos(ωm), v=η sin(ωm),
P_loc = kayan-pencere ortalaması (u²+v²)" diyor. u²+v² ≡ η² olduğundan
bu okuma frekanstan bağımsızdır (bandın hiçbir çizgisini ayırmaz);
demodüle okuma P_dem = |ĉ|² ise ⟨P_dem⟩ = pow/4 + σ_η²/W, yani
W-bond penceresine sızan geniş-bant GÜRÜLTÜ PEDESTALI taşır (bu
koşuda pedestal/koherent = 4.6 … 500). Tahmincinin cebrinden çıkan —
ve mekanizma denklemini ÖZDEŞ yapan — ağırlık

    ρ_n = Re[c_{n+1}·conj⟨c⟩]        (çizginin yerel KOHERENT gücü;
                                      ⟨ρ⟩ = pow/4 tam olarak)

Üç okuma da ölçülür ve tabloya girer:
    S_koh   (ρ, birincil)   S_dem^W (|ĉ|²)   S_tot^W (⟨η²⟩_W)
Kayan pencere duyarlılığı ρ tarafında X'i pencerelemekle ölçülür:
Cov(K*ρ, X) = Cov(ρ, K*X) (simetrik çekirdek) — yani "ağırlığı
pencerele" ile "X'i pencerele" AYNI sınavdır.

KANALLAR (T2): ds = drift + lad + η (156_cekirdek.zincir3). Kovaryans
doğrusal olduğundan S_tam = S_lad + S_eta + S_dri ÖZDEŞ; T2 bir model
seçimi değil TAM bir ayrışımdır.
"""

import importlib
import sys
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
for _p in ("155_configs", "156_configs", "154_configs"):
    sys.path.insert(0, str(QM / _p))
K155 = importlib.import_module("155_cekirdek")
KOS155 = importlib.import_module("155_kos")
K156 = importlib.import_module("156_cekirdek")
C154 = importlib.import_module("154_cekirdek")

TWO_PI = 2 * np.pi
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/159")

IZGARA158 = K155.izgara(0.28, 0.64, 0.02)      # 158'in ızgarası (denetim, a/b)
IZGARA_T1 = K155.izgara(0.44, 0.80, 0.04)      # T1'in geniş taraması (9 bant)

PENCERELER = (32, 64, 128)
KANALLAR = ("tam", "lad", "eta", "dri")


def kayan(x, W):
    """W-bond kayan ORTALAMA, geçerli bölge. ort[j] ↔ merkez j+(W−1)//2."""
    cs = np.concatenate((x[:1] * 0, np.cumsum(x)))
    return (cs[W:] - cs[:-W]) / W, (W - 1) // 2


def _S(w, x):
    """Cov(w,x)/⟨w⟩ — S'nin tanımı."""
    wm = w.mean()
    return float((np.dot(w, x) / len(w) - wm * x.mean()) / wm)


class Yerel:
    """Bir (veri, taban) için kanal dizileri, pencere dilimleri, çizgi motoru."""

    def __init__(self, z, anahtar, taban, cap=4000):
        C = K155.eta_onbellek(z, anahtar, taban, cap)          # 155/158 ile aynı
        self.mid, self.L, self.ds, self.eta = C["mid"], C["L"], C["ds"], C["eta"]
        self.s_ds, self.s_eta, self.c1, self.nq = (C["s_ds"], C["s_eta"],
                                                   C["c1"], C["nq"])
        D = K156.zincir3(z, taban, cap)                        # 156 ile aynı
        self.eta_kapanis = float(np.max(np.abs(D["eta"] - self.eta)))
        self.s_lad, self.s_dri = D["s_lad"], D["s_dri"]
        X = {"tam": K156._bond(D["ds"]), "lad": K156._bond(D["lad"]),
             "eta": K156._bond(D["eta"]), "dri": K156._bond(D["drift"])}
        self.bond_kapanis = float(np.max(np.abs(
            X["tam"] - (X["lad"] + X["eta"] + X["dri"]))))
        self.X = X
        self.sA2 = float(np.var(X["tam"]))
        self.kov = {k: float(np.mean(v * X["tam"])) for k, v in X.items()}
        self.e0, self.e1 = self.eta[:-1], self.eta[1:]
        self.m0 = self.mid[:-1]
        self.N = len(self.m0)
        # tam etkin bond adımı (kinematik kimliği ÖZDEŞ yapan değişken)
        self.Xtil = (self.mid[1:] - self.mid[:-1]) * self.L / TWO_PI - 1.0
        self.Xtil0 = self.Xtil - self.Xtil.mean()
        self.dXtil = float(np.max(np.abs(self.Xtil0 - X["tam"])))
        # pencere dilimleri
        self.dil = {}
        for W in PENCERELER:
            o = (W - 1) // 2
            n = self.N - W + 1
            xw, _ = kayan(X["tam"], W)
            self.dil[W] = dict(o=o, n=n, tam=np.ascontiguousarray(X["tam"][o:o + n]),
                               xw=xw, Ptot=kayan(self.e0 ** 2, W)[0])

    def cizgi(self, Wf, tam_merdiven=True):
        """Bir frekans için tam kayıt: zc/zp/pow + mekanizma merdiveni."""
        cw = np.cos(Wf * self.m0)
        sw = np.sin(Wf * self.m0)
        zc = complex(2 * np.mean(self.e0 * cw), -2 * np.mean(self.e0 * sw))
        zp = complex(2 * np.mean(self.e1 * cw), -2 * np.mean(self.e1 * sw))
        A = TWO_PI * Wf / self.L
        out = dict(pow=float(abs(zc) ** 2), zc=zc, zp=zp, A=A,
                   argraw=float(np.angle(zp * np.conj(zc))))
        # c_{n+1} = η_{n+1}·e^{−iW·m_{n+1}}
        cw1 = np.cos(Wf * self.mid[1:])
        sw1 = np.sin(Wf * self.mid[1:])
        c1r, c1i = self.e1 * cw1, -self.e1 * sw1
        zb = zc / 2.0                       # ⟨c⟩
        rho = c1r * zb.real + c1i * zb.imag           # Re[c₁ conj⟨c⟩]
        sig = c1i * zb.real - c1r * zb.imag           # Im[c₁ conj⟨c⟩]
        rm = float(rho.mean())
        out["rho_ort"] = rm
        out["rho_ort_norm"] = rm / (out["pow"] / 4.0) if out["pow"] else np.nan
        Xt = self.X["tam"]
        for k in KANALLAR:                            # T2: TAM ayrışım
            out[f"S_{k}"] = _S(rho, self.X[k])
        out["S_xtil"] = _S(rho, self.Xtil0)           # X̃ değişkeniyle
        for W in PENCERELER:                          # pencere duyarlılığı
            d = self.dil[W]
            o, n = d["o"], d["n"]
            out[f"S_xw{W}"] = _S(rho[o:o + n], d["xw"])   # X pencerelenmiş
            cb, _ = kayan(self.e0 * cw - 1j * (self.e0 * sw), W)
            Pd = cb.real ** 2 + cb.imag ** 2
            out[f"S_dem{W}"] = _S(Pd, d["tam"])
            out[f"ped{W}"] = float(Pd.mean() / (out["pow"] / 4.0)) \
                if out["pow"] else np.nan
            out[f"S_tot{W}"] = _S(d["Ptot"], d["tam"])
        if tam_merdiven:
            ph = np.exp(1j * A * Xt)
            out["M1_re"] = float(np.dot(rho, ph.real) / len(Xt) / rm)
            out["M1_im"] = float(np.dot(rho, ph.imag) / len(Xt) / rm)
            out["M0_re"] = float(np.dot(sig, ph.real) / len(Xt) / rm)
            out["M0_im"] = float(np.dot(sig, ph.imag) / len(Xt) / rm)
        return out


def olcS(z, etiket, bantlar, anahtar, taban=0.40, cap=4000, njack=8,
         tohum=21, ayrinti=True, tam_merdiven=True):
    """158'in bant/çizgi döngüsü + her çizgide mekanizma merdiveni."""
    Y = Yerel(z, anahtar, taban, cap)
    L = Y.L
    print(f"[{etiket}] N={len(z)}  L={L:.4f}  σΔ²={Y.sA2:.5f}  "
          f"σ_ds²={Y.s_ds:.5f}  σ_η²={Y.s_eta:.5f}  "
          f"η-kapanış={Y.eta_kapanis:.1e}  bond-kapanış={Y.bond_kapanis:.1e}  "
          f"maks|X̃−dsΔ|={Y.dXtil:.2e}", flush=True)
    print(f"  pay yapısı Cov(X_k,X_tam)/σΔ²: lad={Y.kov['lad']/Y.sA2:+.4f} "
          f"eta={Y.kov['eta']/Y.sA2:+.4f} dri={Y.kov['dri']/Y.sA2:+.4f}",
          flush=True)

    qm = C154.pk_m(int(np.exp(0.86 * L)))
    allq = sorted(qm)
    allw = np.array([np.log(q) for q in allq])
    dres = TWO_PI / (Y.mid[-1] - Y.mid[0])
    rng = np.random.default_rng(tohum)

    satir = []
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
                       gp=float((2 * aq * np.sin(np.pi * w / L)) ** 2))
            for Wf, et in ((w, "on"), (w + gap / 2, "off")):
                r = Y.cizgi(Wf, tam_merdiven)
                cr = r["zp"] * np.conj(r["zc"]) * np.exp(1j * TWO_PI * Wf / L)
                rec[f"A_{et}"] = r["A"]
                rec[f"cr_{et}_re"] = float(cr.real)
                rec[f"cr_{et}_im"] = float(cr.imag)
                rec[f"pow_{et}"] = r["pow"]
                rec[f"argraw_{et}"] = r["argraw"]
                for kk, vv in r.items():
                    if kk[0] in "SMp" and kk not in ("pow",):
                        rec[f"{kk}_{et}"] = vv
                rec[f"rho_norm_{et}"] = r["rho_ort_norm"]
            Ls.append(rec)

        tb = 0.5 * (lo + hi)
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
        Aon = np.array([r["A_on"] for r in Ls])
        Aof = np.array([r["A_off"] for r in Ls])
        u = pN - pO
        payda = float(u.sum())
        if abs(payda) < 1e-300:
            satir.append(dict(tau=round(tb, 4), lo=lo, hi=hi, N=len(tumu),
                              kul=kul, olculdu=False))
            continue
        G = complex((crN - crO).sum() / payda)
        phm = float(np.angle(G))
        tau_eff = float((tv * u).sum() / payda)
        Gc = complex((crN * np.exp(-2j * Aon) - crO * np.exp(-2j * Aof)).sum()
                     / payda)
        dlt = float(np.angle(Gc))

        def col(key, et):
            return np.array([r[f"{key}_{et}"] for r in Ls])

        def AS(key, msk=None):
            """Σ(pN·A_on·S_on − pO·A_off·S_off)/Σ(pN−pO)."""
            m = slice(None) if msk is None else msk
            d = payda if msk is None else float(u[msk].sum())
            return float(((pN[m] * Aon[m] * col(key, "on")[m]
                           - pO[m] * Aof[m] * col(key, "off")[m]).sum()) / d)

        def AS_on(key):
            return float((pN * Aon * col(key, "on")).sum() / float(pN.sum()))

        def Sb(key):
            return float((pN * col(key, "on") - pO * col(key, "off")).sum()
                         / payda)

        ASt, Sbt, ASo = {}, {}, {}
        anahtarlar = ([f"S_{k}" for k in KANALLAR] + ["S_xtil"] +
                      [f"S_xw{W}" for W in PENCERELER] +
                      [f"S_dem{W}" for W in PENCERELER] +
                      [f"S_tot{W}" for W in PENCERELER])
        for key in anahtarlar:
            ASt[key] = AS(key)
            Sbt[key] = Sb(key)
        for key in ("S_tam", "S_dem64"):
            ASo[key] = AS_on(key)
        ped = {f"ped{W}": float(np.mean(col(f"ped{W}", "on")))
               for W in PENCERELER}

        # tüm-mertebeli mekanizma (σ kanalı ayrı)
        M1 = M0 = float("nan")
        if tam_merdiven:
            z1 = ((pN * (col("M1_re", "on") + 1j * col("M1_im", "on")))
                  - (pO * (col("M1_re", "off") + 1j * col("M1_im", "off")))
                  ).sum() / payda
            zs = ((pN * (col("M0_re", "on") + 1j * col("M0_im", "on")))
                  - (pO * (col("M0_re", "off") + 1j * col("M0_im", "off")))
                  ).sum() / payda
            M1 = float(np.angle(z1))
            M0 = float(np.angle(z1 + 1j * zs))

        nanv = float("nan")
        phj = [nanv] * njack; dlj = [nanv] * njack
        tej = [nanv] * njack; asj = [nanv] * njack
        for kk in range(njack):
            m = grp != kk
            if not m.any():
                continue
            pk_ = float(u[m].sum())
            if abs(pk_) < 1e-300:
                continue
            Gk = complex((crN[m] - crO[m]).sum() / pk_)
            if abs(Gk) > 1.5:
                continue
            phj[kk] = float(np.angle(Gk))
            tej[kk] = float((tv[m] * u[m]).sum() / pk_)
            dlj[kk] = float(np.angle(
                (crN[m] * np.exp(-2j * Aon[m])
                 - crO[m] * np.exp(-2j * Aof[m])).sum() / pk_))
            asj[kk] = AS("S_tam", m)

        def jkerr(v):
            v = np.asarray(v, float)
            v = v[np.isfinite(v)]
            if len(v) < 4:
                return float("nan")
            return float(np.sqrt((len(v) - 1) / len(v)
                                 * np.sum((v - v.mean()) ** 2)))

        kin = 2 * (TWO_PI * tau_eff) - TWO_PI
        rec = dict(tau=round(tb, 4), lo=lo, hi=hi, N=len(tumu), kul=kul,
                   olculdu=True, L=L, sA2=Y.sA2,
                   absG=float(abs(G)), Gre=float(G.real), Gim=float(G.imag),
                   phi=phm, sPhi_jk=jkerr(phj),
                   delta=dlt, sDelta_jk=jkerr(dlj), absGc=float(abs(Gc)),
                   delta_b=float(phm - kin), kinematik=float(kin),
                   tau_eff=tau_eff, sTeff_jk=jkerr(tej), tau_ari=float(tv.mean()),
                   A_eff=float(TWO_PI * tau_eff),
                   M1=M1, M0=M0, AS=ASt, S=Sbt, AS_on=ASo, ped=ped,
                   sAS_jk=jkerr(asj), rho=payda / float(gpv.sum()),
                   rho_norm=float(np.mean(col("rho_norm", "on"))),
                   phi_jk=phj, delta_jk=dlj, teff_jk=tej)
        satir.append(rec)
        if ayrinti:
            r1 = dlt / ASt["S_tam"] if ASt["S_tam"] else float("nan")
            print(f"  {tb:.4f} τe={tau_eff:.4f} kul={kul:3d} |Γ|={abs(G):.3f} "
                  f"φ={phm:+.4f} kin={kin:+.4f} δ={dlt:+.4f} M1={M1:+.4f} "
                  f"A·S={ASt['S_tam']:+.4f} δ/A·S={r1:+.3f} "
                  f"ped64={ped['ped64']:.1f}", flush=True)
    return dict(etiket=etiket, L=L, sA2=Y.sA2, s_ds=Y.s_ds, s_eta=Y.s_eta,
                s_lad=Y.s_lad, s_dri=Y.s_dri, c1=Y.c1, taban=taban, cap=cap,
                kov={k: Y.kov[k] for k in Y.kov}, dXtil=Y.dXtil,
                eta_kapanis=Y.eta_kapanis, bond_kapanis=Y.bond_kapanis,
                pencereler=list(PENCERELER), bantlar=satir)
