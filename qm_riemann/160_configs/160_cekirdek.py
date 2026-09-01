"""
160 — δ'NIN TAM MUHASEBESİ: kuadratür (σ) kanalı + sonlu-kesme. ÖLÇÜM ÇEKİRDEĞİ
==============================================================================
159 iki açık çarpan bıraktı: (i) birinci-mertebe kesme (A·σ_X = 0.68…1.14,
küçük DEĞİL) ve (ii) kuadratür (σ) kanalı. Bu modül ikisini de AÇILIM
YAPMADAN, TAM ifadeyle ölçer.

────────────────────────────────────────────────────────────────────────────
MUHASEBE (üç ölçülebilir bileşen; hiçbiri fit değil)
────────────────────────────────────────────────────────────────────────────
159'un V2 ile 7.8e−09'da doğruladığı özdeşlik:

    zp·conj(zc)·e^{−iA}/4 = ⟨(ρ_n + i·σ_n)·e^{+iA·X̃_n}⟩
    ρ+iσ ≡ c_{n+1}·conj⟨c⟩ ,  c_n = η_n e^{−iW m_n} ,  A = 2πW/L
    X̃_n = (m_{n+1}−m_n)·L/2π − 1        (TAM etkin bond adımı)

Buradan δ = arg[zp conj(zc) e^{−iA}] için üç basamaklı TAM muhasebe:

    (K1)  reel-ağırlık, TAM karakteristik fonksiyon
          E1 = ⟨ρ e^{iAX̃}⟩/⟨ρ⟩ ,      K1 = arg E1
          — küçük-faz açılımı YOK; 159'un "sonlu-kesme" çarpanı burada
            kendiliğinden içerilir (A·S, E1'in BİRİNCİ mertebesidir).
    (K2)  kuadratür düzeltmesi
          Es = ⟨σ e^{iAX̃}⟩/⟨ρ⟩ ,      K2 = arg(E1 + i·Es) − arg(E1)
    (K3)  artık:  δ_meas − (K1+K2)

X̃ ile K1+K2 ÖZDEŞ olarak δ'dır (kapanış ~1e−15; V2). Bu bir sınav değil,
bir kimliktir — ve muhasebenin KAPALI olduğunun kanıtıdır. SINAV, KALEM'in
kendi değişkeni dsΔ ile koşulan sürümdür:

    M1 = arg⟨ρ e^{iA·dsΔ}⟩/⟨ρ⟩ ,  M0 = arg[⟨(ρ+iσ)e^{iA·dsΔ}⟩/⟨ρ⟩]
    K3(dsΔ) = δ − M0    ← değişken ikamesi sistematiği, ÖLÇÜLÜR

────────────────────────────────────────────────────────────────────────────
KESMENİN KAPALI FORMU (kümülant açılımı — 159'un 3. sıradaki adımı)
────────────────────────────────────────────────────────────────────────────
X̃0 = X̃ − ⟨X̃⟩ ; ρ-ölçüsünün merkezi momentleri μ_k = ⟨ρ X̃0^k⟩/⟨ρ⟩ ;
kümülantlar κ1 = μ1, κ2 = μ2−μ1², κ3 = μ3−3μ1μ2+2μ1³,
κ4 = μ4−4μ1μ3−3μ2²+12μ1²μ2−6μ1⁴.

    K1 = Im log E1 = A·κ1 − A³·κ3/6 + O(A⁵)      (+ A⟨X̃⟩, ~5e−05)
    log|E1| = −A²κ2/2 + A⁴κ4/24 + O(A⁶)          (SÖNÜM)

A·S ≡ A·Cov(ρ,X̃)/⟨ρ⟩ = A·κ1 (kovaryans merkezlemeye duyarsız), yani
KALEM'in ifadesi kümülant serisinin BİRİNCİ terimidir. −A³κ3/6 terimi
"sonlu-kesme" çarpanının kapalı formudur ve burada ayrı ölçülür.

σ kanalının serisi:  Es = s0 + iA·s1 − A²s2/2 − iA³s3/6 + …,
s_k = ⟨σ X̃0^k⟩/⟨ρ⟩. Re[Es]'in ilk X-bağımlı terimi −A²s2/2 — yani
σ kanalı A'da İKİNCİ mertebede girer; reel-ağırlık + küçük-faz ansatzının
iki ayrı nedenle göremediği bir kanaldır.

────────────────────────────────────────────────────────────────────────────
σ'NIN KİMLİĞİ (kalem): ρ'nun TAŞIYICI FREKANSTAKİ KUADRATÜRÜ
────────────────────────────────────────────────────────────────────────────
u_n ≡ ρ_n + iσ_n = c_{n+1} conj⟨c⟩ için TAM özyineleme:

    u_{n+1} = (η_{n+2}/η_{n+1}) · u_n · e^{−iA(1+X̃_{n+1})}

Yani (ρ,σ) bir FAZÖRDÜR: bond başına −A dönüyor, genliği |η| ile
modüle. Genlik yavaş değişirken
    Δρ_n ≡ ρ_{n+1}−ρ_n ≈ (cos A − 1)·ρ_n + sin A · σ_n
    ⇒  σ_n ≈ tan(A/2)·ρ_n + Δρ_n / sin A
(σ'yı ρ ve gradyanına regres ederek α,β,R² ölçülür; öngörü
α_pred = tan(A/2), β_pred = 1/sin A.)

Demodüle edilmiş fazör TAM olarak
    v_n ≡ u_n·e^{+iA(n+1)} = η_{n+1}·e^{−iA·C_n}·(e^{−iW m_0}conj⟨c⟩),
    C_n = Σ_{k≤n} X̃_k   (BİRİKMİŞ adım sapması)
⇒ **σ ≠ 0 ⟺ çizgi küresel taşıyıcıya göre faz BİRİKTİRMİŞTİR.** σ, yerel
frekans kaymasının (faz difüzyonunun) kuadratür izidir; |Γ|'nın sönümü de
aynı C_n'in dağılımından gelir (|E1| ≈ e^{−A²κ2/2}).

────────────────────────────────────────────────────────────────────────────
HIZ NOTU: e^{iAX̃} YENİ TRİGONOMETRİ İSTEMEZ
────────────────────────────────────────────────────────────────────────────
A·X̃_n = W(m_{n+1}−m_n) − A olduğundan
    e^{iAX̃} = e^{−iA}·(cos Wm_{n+1} + i sin Wm_{n+1})(cos Wm_n − i sin Wm_n)
ve bu dört dizi zaten zc/zp/ρ/σ için hesaplanıyor. 159'un pencere
makinesi (PENCERELER: |ĉ|², ⟨η²⟩_W) BU KOŞUDA KOŞULMAZ — 159 §4a/4b onu
kapattı; yerine kümülantlar ve σ tanıları geliyor.
"""

import importlib
import sys
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
sys.path.insert(0, str(QM / "159_configs"))
C159 = importlib.import_module("159_cekirdek")
K155, KOS155, C154 = C159.K155, C159.KOS155, C159.C154

TWO_PI = 2 * np.pi
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/160")

IZGARA158 = C159.IZGARA158          # izgara(0.28,0.64,0.02)
IZGARA_T1 = C159.IZGARA_T1          # izgara(0.44,0.80,0.04)
KANALLAR = C159.KANALLAR            # ("tam","lad","eta","dri")


def _S(w, x):
    """Cov(w,x)/⟨w⟩ — 159'un tanımı (bit düzeyinde aynı fonksiyon)."""
    return C159._S(w, x)


def _Sr(w, x, rm):
    """Cov(w,x)/⟨ρ⟩ — σ gibi ORTALAMASI SIFIRA YAKIN ağırlıklar için."""
    wm = w.mean()
    return float((np.dot(w, x) / len(w) - wm * x.mean()) / rm)


def _kor(a, b):
    """Pearson korelasyonu (tam dizi)."""
    a = a - a.mean()
    b = b - b.mean()
    da = float(np.dot(a, a)) ** 0.5
    db = float(np.dot(b, b)) ** 0.5
    if da == 0.0 or db == 0.0:
        return float("nan")
    return float(np.dot(a, b) / (da * db))


class Yerel160(C159.Yerel):
    """159'un Yerel'i + merkezi kuvvetler, ρ/X gradyanları (σ tanıları)."""

    def __init__(self, z, anahtar, taban, cap=4000):
        super().__init__(z, anahtar, taban, cap)
        X0 = self.Xtil0                        # X̃ − ⟨X̃⟩
        self.Xort = float(self.Xtil.mean())
        self.P1_ = X0
        self.P2_ = X0 * X0
        self.P3_ = self.P2_ * X0
        self.P4_ = self.P2_ * self.P2_
        # dX̃/dn (ileri fark; son eleman 0 — 300k'da 1 elemanlık sınır etkisi)
        self.dX = np.concatenate((np.diff(X0), np.zeros(1)))

    # ---------------------------------------------------------------
    def cizgi(self, Wf, tanilar=False):
        """Bir frekans için TAM kayıt (159'un sayıları + 160'ın muhasebesi)."""
        # --- 159 ile BİT DÜZEYİNDE aynı blok ------------------------
        cw = np.cos(Wf * self.m0)
        sw = np.sin(Wf * self.m0)
        zc = complex(2 * np.mean(self.e0 * cw), -2 * np.mean(self.e0 * sw))
        zp = complex(2 * np.mean(self.e1 * cw), -2 * np.mean(self.e1 * sw))
        A = TWO_PI * Wf / self.L
        out = dict(pow=float(abs(zc) ** 2), zc=zc, zp=zp, A=A,
                   argraw=float(np.angle(zp * np.conj(zc))))
        cw1 = np.cos(Wf * self.mid[1:])
        sw1 = np.sin(Wf * self.mid[1:])
        c1r, c1i = self.e1 * cw1, -self.e1 * sw1
        zb = zc / 2.0                                   # ⟨c⟩
        rho = c1r * zb.real + c1i * zb.imag             # Re[c₁ conj⟨c⟩]
        sig = c1i * zb.real - c1r * zb.imag             # Im[c₁ conj⟨c⟩]
        rm = float(rho.mean())
        out["rho_ort"] = rm
        out["rho_ort_norm"] = rm / (out["pow"] / 4.0) if out["pow"] else np.nan
        Xt = self.X["tam"]
        for k in KANALLAR:                              # T2 (159 ile aynı)
            out[f"S_{k}"] = _S(rho, self.X[k])
        out["S_xtil"] = _S(rho, self.Xtil0)
        ph = np.exp(1j * A * Xt)
        out["M1_re"] = float(np.dot(rho, ph.real) / len(Xt) / rm)
        out["M1_im"] = float(np.dot(rho, ph.imag) / len(Xt) / rm)
        out["M0_re"] = float(np.dot(sig, ph.real) / len(Xt) / rm)
        out["M0_im"] = float(np.dot(sig, ph.imag) / len(Xt) / rm)
        # --- 160: TAM değişken X̃ ile karakteristik fonksiyon --------
        # e^{iAX̃} = e^{−iA}·e^{iW(m_{n+1}−m_n)}  (yeni trigonometri YOK)
        dr = cw1 * cw + sw1 * sw                        # cos(W Δm)
        di = sw1 * cw - cw1 * sw                        # sin(W Δm)
        ca, sa = float(np.cos(A)), float(np.sin(A))
        Cx = dr * ca + di * sa
        Sx = di * ca - dr * sa
        N = len(rho)
        sm = float(sig.mean())
        out["E1_re"] = float(np.dot(rho, Cx) / N / rm)
        out["E1_im"] = float(np.dot(rho, Sx) / N / rm)
        out["Es_re"] = float(np.dot(sig, Cx) / N / rm)
        out["Es_im"] = float(np.dot(sig, Sx) / N / rm)
        # kapanış denetimi: E1 + i·Es = [zp conj(zc) e^{−iA}/4]/⟨ρ⟩
        lhs = complex(out["E1_re"] - out["Es_im"], out["E1_im"] + out["Es_re"])
        rhs = zp * np.conj(zc) * complex(ca, -sa) / 4.0 / rm
        out["kapanis"] = float(abs(lhs - rhs) / max(abs(rhs), 1e-300))
        # --- kümülantlar (ρ ve σ ölçüleri, X̃0 merkezli) -------------
        u1 = float(np.dot(rho, self.P1_) / N / rm)
        u2 = float(np.dot(rho, self.P2_) / N / rm)
        u3 = float(np.dot(rho, self.P3_) / N / rm)
        u4 = float(np.dot(rho, self.P4_) / N / rm)
        out["k1"] = u1
        out["k2"] = u2 - u1 * u1
        out["k3"] = u3 - 3 * u1 * u2 + 2 * u1 ** 3
        out["k4"] = u4 - 4 * u1 * u3 - 3 * u2 * u2 + 12 * u1 * u1 * u2 \
            - 6 * u1 ** 4
        out["s0"] = sm / rm
        out["s1"] = float(np.dot(sig, self.P1_) / N / rm)
        out["s2"] = float(np.dot(sig, self.P2_) / N / rm)
        out["s3"] = float(np.dot(sig, self.P3_) / N / rm)
        out["Ss_tam"] = _Sr(sig, Xt, rm)          # Cov(σ,dsΔ)/⟨ρ⟩
        out["Ss_xtil"] = _Sr(sig, self.Xtil0, rm)  # Cov(σ,X̃)/⟨ρ⟩
        # --- σ TANILARI (isteğe bağlı; ek ~8 nokta çarpımı) ---------
        if tanilar:
            dro = np.concatenate((np.diff(rho), np.zeros(1)))
            out["c_sr"] = _kor(sig, rho)
            out["c_sdr"] = _kor(sig, dro)
            out["c_sx"] = _kor(sig, self.P1_)
            out["c_sdx"] = _kor(sig, self.dX)
            out["c_rx"] = _kor(rho, self.P1_)
            out["c_rdx"] = _kor(rho, self.dX)
            # σ ≈ α·ρ + β·Δρ  (kuadratür kimliği); öngörü α=tan(A/2), β=1/sinA
            a11 = float(np.dot(rho, rho)); a12 = float(np.dot(rho, dro))
            a22 = float(np.dot(dro, dro))
            b1 = float(np.dot(rho, sig)); b2 = float(np.dot(dro, sig))
            det = a11 * a22 - a12 * a12
            if abs(det) > 0:
                al = (a22 * b1 - a12 * b2) / det
                be = (a11 * b2 - a12 * b1) / det
                ss = float(np.dot(sig, sig))
                out["q_al"] = float(al)
                out["q_be"] = float(be)
                out["q_R2"] = float((al * b1 + be * b2) / ss) if ss > 0 else np.nan
            else:
                out["q_al"] = out["q_be"] = out["q_R2"] = float("nan")
            out["q_al_p"] = float(np.tan(A / 2.0))
            out["q_be_p"] = float(1.0 / sa) if sa != 0 else float("nan")
        return out


# ---------------------------------------------------------------------
# çizgi kayıtlarında saklanan skaler anahtarlar
SKA = ([f"S_{k}" for k in KANALLAR] + ["S_xtil", "Ss_tam", "Ss_xtil",
                                       "k1", "k2", "k3", "k4",
                                       "s0", "s1", "s2", "s3", "kapanis"])
KOMPLEKS = ("M1", "M0", "E1", "Es")
TANI = ("c_sr", "c_sdr", "c_sx", "c_sdx", "c_rx", "c_rdx",
        "q_al", "q_be", "q_R2", "q_al_p", "q_be_p")


def olc160(z, etiket, bantlar, anahtar, taban=0.40, cap=4000, njack=8,
           tohum=21, ayrinti=True, tanilar=True):
    """159'un bant/çizgi döngüsü (aynı tohum, aynı filtre, aynı jackknife)
    + δ'nın TAM muhasebesi (K1, K2, K3) + kümülantlar + σ tanıları."""
    Y = Yerel160(z, anahtar, taban, cap)
    L = Y.L
    print(f"[{etiket}] N={len(z)}  L={L:.4f}  σΔ²={Y.sA2:.5f}  "
          f"σ_ds²={Y.s_ds:.5f}  σ_η²={Y.s_eta:.5f}  "
          f"η-kapanış={Y.eta_kapanis:.1e}  bond-kapanış={Y.bond_kapanis:.1e}  "
          f"maks|X̃−dsΔ|={Y.dXtil:.2e}  ⟨X̃⟩={Y.Xort:.3e}", flush=True)

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
                r = Y.cizgi(Wf, tanilar=(tanilar and et == "on"))
                cr = r["zp"] * np.conj(r["zc"]) * np.exp(1j * TWO_PI * Wf / L)
                rec[f"A_{et}"] = r["A"]
                rec[f"cr_{et}_re"] = float(cr.real)
                rec[f"cr_{et}_im"] = float(cr.imag)
                rec[f"pow_{et}"] = r["pow"]
                rec[f"argraw_{et}"] = r["argraw"]
                for kk in SKA:
                    rec[f"{kk}_{et}"] = r[kk]
                for kk in KOMPLEKS:
                    rec[f"{kk}_re_{et}"] = r[f"{kk}_re"]
                    rec[f"{kk}_im_{et}"] = r[f"{kk}_im"]
                rec[f"rho_norm_{et}"] = r["rho_ort_norm"]
                if tanilar and et == "on":
                    for kk in TANI:
                        rec[kk] = r[kk]
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

        def AS(key, msk=None, p=1):
            """Σ(pN·A_on^p·v_on − pO·A_off^p·v_off)/Σ(pN−pO)."""
            m = slice(None) if msk is None else msk
            d = payda if msk is None else float(u[msk].sum())
            return float(((pN[m] * Aon[m] ** p * col(key, "on")[m]
                           - pO[m] * Aof[m] ** p * col(key, "off")[m]).sum()) / d)

        def CZ(key, msk=None):
            """Kompleks nesnenin bant birleştirmesi (159'un z1/zs'i ile aynı)."""
            m = slice(None) if msk is None else msk
            d = payda if msk is None else float(u[msk].sum())
            vN = col(f"{key}_re", "on") + 1j * col(f"{key}_im", "on")
            vO = col(f"{key}_re", "off") + 1j * col(f"{key}_im", "off")
            return complex(((pN[m] * vN[m] - pO[m] * vO[m]).sum()) / d)

        ASt = {k: AS(k) for k in
               [f"S_{c}" for c in KANALLAR] + ["S_xtil", "Ss_tam", "Ss_xtil"]}
        Sbt = {k: float((pN * col(k, "on") - pO * col(k, "off")).sum() / payda)
               for k in [f"S_{c}" for c in KANALLAR] + ["S_xtil"]}
        ASo = {"S_tam": float((pN * Aon * col("S_tam", "on")).sum()
                              / float(pN.sum()))}
        # kümülant terimleri: A^p·κ_p bant birleştirmesi
        KUM = {f"A{p}k{p}": AS(f"k{p}", p=p) for p in (1, 2, 3, 4)}
        KUM.update({f"A{p}s{p}": AS(f"s{p}", p=p) for p in (1, 2, 3)})
        KUM["s0"] = AS("s0", p=0)

        z1 = CZ("M1"); zs = CZ("M0")            # dsΔ değişkeni (159)
        zE = CZ("E1"); zS = CZ("Es")            # X̃ değişkeni (TAM)
        M1 = float(np.angle(z1))
        M0 = float(np.angle(z1 + 1j * zs))
        K1 = float(np.angle(zE))
        Kf = float(np.angle(zE + 1j * zS))
        K2 = float(np.angle((zE + 1j * zS) / zE))

        nanv = float("nan")
        phj = [nanv] * njack; dlj = [nanv] * njack
        tej = [nanv] * njack; asj = [nanv] * njack
        k1j = [nanv] * njack; kfj = [nanv] * njack; m1j = [nanv] * njack
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
            k1j[kk] = float(np.angle(CZ("E1", m)))
            kfj[kk] = float(np.angle(CZ("E1", m) + 1j * CZ("Es", m)))
            m1j[kk] = float(np.angle(CZ("M1", m)))

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
                   M1=M1, M0=M0, K1=K1, K2=K2, Kfull=Kf,
                   absE1=float(abs(zE)), absEf=float(abs(zE + 1j * zS)),
                   E1_re=float(zE.real), E1_im=float(zE.imag),
                   Es_re=float(zS.real), Es_im=float(zS.imag),
                   AS=ASt, S=Sbt, AS_on=ASo, KUM=KUM,
                   sAS_jk=jkerr(asj), sK1_jk=jkerr(k1j), sKf_jk=jkerr(kfj),
                   sM1_jk=jkerr(m1j),
                   rho=payda / float(gpv.sum()),
                   rho_norm=float(np.mean(col("rho_norm", "on"))),
                   kapanis=float(np.max(np.abs(col("kapanis", "on")))),
                   phi_jk=phj, delta_jk=dlj, teff_jk=tej)
        if tanilar:
            w_ = pN / float(pN.sum())
            rec["TANI"] = {k: float(np.dot(w_, np.array([r[k] for r in Ls])))
                           for k in TANI}
        satir.append(rec)
        if ayrinti:
            AS0 = ASt["S_tam"]
            print(f"  {tb:.4f} τe={tau_eff:.4f} kul={kul:3d} |Γ|={abs(G):.3f} "
                  f"φ={phm:+.4f} δ={dlt:+.4f} K1={K1:+.4f} K2={K2:+.4f} "
                  f"K1+K2={Kf:+.4f} K3={dlt-Kf:+.2e} M0={M0:+.4f} "
                  f"A·S={AS0:+.4f} δ/A·S={dlt/AS0 if AS0 else np.nan:+.3f} "
                  f"kap={rec['kapanis']:.1e}", flush=True)
    return dict(etiket=etiket, L=L, sA2=Y.sA2, s_ds=Y.s_ds, s_eta=Y.s_eta,
                s_lad=Y.s_lad, s_dri=Y.s_dri, c1=Y.c1, taban=taban, cap=cap,
                kov={k: Y.kov[k] for k in Y.kov}, dXtil=Y.dXtil,
                Xort=Y.Xort, eta_kapanis=Y.eta_kapanis,
                bond_kapanis=Y.bond_kapanis, bantlar=satir)
