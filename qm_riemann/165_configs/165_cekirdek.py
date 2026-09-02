"""
165 — Ψ⋆κ: FAZ-KİLİDİ KATMANININ KANAL-ETİKETLİ, PENCERE-FAZLI HESABI
=====================================================================
Ölçüm zinciri KOPYALANMAZ: `163_cekirdek` (o da 160/159/156/155/154'ü)
import edilir. `bant_adaylari`, `olc_cizgi`, `gaz` aynen kullanılır;
165 yalnız ÖNGÖRÜ tarafını değiştirir.

══════════════════════════════════════════════════════════════════════
0. 163'ÜN İHMALİ ve 165'İN DÜZELTMESİ (KALEM 2 Eylül)
══════════════════════════════════════════════════════════════════════
163 üçlü-toplamı TAM REZONANS (ν ≡ 0) koşuluna kesti; ν ≠ 0 üçlüleri
ya attı (T1) ya da KUTU pencereyle birim ağırlıkla topladı (T2 → W ile
doğrusal patlama). Doğru ağırlık pencere dönüşümüdür:

    ölçülen ortalama  ⟨·⟩ = (1/N) Σ_n (·)      ⇒
    e^{iν s_n} teriminin ağırlığı  κ(ν) ≡ (1/N) Σ_n e^{iν s_n}

`κ` üç şeyi birden taşır: (i) sonlu pencere → |κ| ~ |sinc(νT/2)|,
(ii) pencere MERKEZİ → arg κ ≈ ν·s̄  (s̄ ~ 1.05e6 ⇒ ν ~ 1e-6 bile tam tur),
(iii) TARAĞIN kendi yapısı → ν = ±ω_q'da κ'nın rezonans tepeleri
(143'ün 𝒢'si). ν ≡ 0'da κ = 1 ÖZDEŞ.

══════════════════════════════════════════════════════════════════════
1. TEK SİTEYE İNDİRGEME (165'in kilit basitleştirmesi)
══════════════════════════════════════════════════════════════════════
Ölçülen nesne (160 §2, 163 §2b ile ÖZDEŞ):

    c_{n+1} = η_{n+1} e^{−iω_Q m_{n+1}} ,  ⟨c⟩ = h_Q/2 ,
    ρ_n + iσ_n = c_{n+1} conj⟨c⟩ ,  ⟨ρ⟩ = rm (ÖLÇÜLEN),
    s_k = ⟨σ X̃0^k⟩/⟨ρ⟩ = Im[ h̄_Q J_k ] / (2 rm) ,  J_k = ⟨c_{n+1} X̃0_n^k⟩

163 J_k'yı m_n sitesinde yazdı ve m_{n+1} = m_n + ḡ(1+X̃_n) farkını
`Γ(θ,k)` (X̃'nin marjinali) ile taşımak zorunda kaldı — birinci mertebede
kesilen yer orası. 165 BÜTÜN alanları TEK sitede, **s_n ≡ m_{n+1}**'de
yazar; o zaman hiçbir kaydırma çekirdeği kalmaz:

    η_{n+1} = Σ_q Re[ h'_q e^{iω_q s_n} ] + artık ,  h'_q ≡ 2⟨η_{n+1}e^{−iω_q s_n}⟩
    X̃0_n   = Σ_q Re[ y_q  e^{iω_q s_n} ] + artık ,  y_q  ≡ 2⟨X̃0_n e^{−iω_q s_n}⟩

ÇIPLAK (merdiven) DEĞERLERİ — 163 §2a'dan tek satırda:
    ds_n = Σ_q b_q cos(ω_q m_n),  b_q = 2a_q sin(πτ_q)
    ⇒ η_{n+1} ⊃ b_q cos(ω_q s_n)                    ⇒  **h'_q = b_q**  (reel)
    X̃0_n = ½(ds_n+ds_{n+1}) − ort = Σ_q B_q cos(ω_q m_n + πτ_q),
           B_q = b_q cos(πτ_q) = a_q sin(2πτ_q)
    m_n = s_n − ḡ, ω_q ḡ = 2πτ_q ⇒ **y_q = B_q e^{−iπτ_q}**
(probe1: gerçekleşen gazda |y_q|/B_q = 0.99 ve arg y_q + πτ_q = −0.02
τ ≤ 0.2'de — çıplak değerler ölçümle doğrulanmıştır.)

══════════════════════════════════════════════════════════════════════
2. BELİRLENİMLİ DÖRTLÜ TOPLAM — TAM, KANAL-ETİKETLİ
══════════════════════════════════════════════════════════════════════
    J_2 = (1/8) Σ_{(1,ε₁)} Σ_{(2,ε₂)} Σ_{(3,ε₃)}
              h'^{(ε₁)}_1 y^{(ε₂)}_2 y^{(ε₃)}_3 · κ(Δ) ,
    Δ = ε₁ω₁ − ω_Q + ε₂ω₂ + ε₃ω₃          (x^{(+)}=x, x^{(−)}=x̄)

Dört merdiven frekansı (taşıyıcı ω_Q dahil) ⇒ "dördüncü kat".
κ TAM alınır: pencere uzunluğu, pencere merkezi fazı (arg κ ≈ ν s̄) ve
tarağın kendi yapısı içindedir.

**ÖZDEŞLİK (165'in hesap hilesi).** κ(Δ) = (1/N)Σ_n e^{iΔ s_n} olduğundan
toplam ve n-ortalaması yer değiştirir:

    J_2 = ⟨ [Σ_q Re(h'_q e^{iω_q s})] · e^{−iω_Q s} · [Σ_q Re(y_q e^{iω_q s})]² ⟩

yani BÜTÜN kombinasyonların κ-ağırlıklı toplamı = model alanlarının
gerçek tarak üzerindeki ortalaması. Budama YOK, kesme YOK: (2·8981)³ ≈
5.8e12 terim O(N·n_çizgi)'de TAM toplanır. (Alt-merdivende doğrudan
sayım ile 1e−12'de doğrulanır — V4.)

KANAL ETİKETLERİ (görevin tanımı):
  • **Ç1 — çift-iptal, ν ≡ 0 ARİTMETİKSİZ.** {ε₁ω₁, −ω_Q, ε₂ω₂, ε₃ω₃}
    ikişerli birebir sadeleşir. Üç aile:
      (A) q₁=Q,ε₁=+ ; q₂=q₃=q, ε₂=−ε₃  → T_a: (h'_Q/4)Σ_q|y_q|²
      (B) q₂=Q,ε₂=+ ; q₃=q₁, ε₃=−ε₁    → (y_Q/8)·2Σ_{q₁}Re[h'_1 ȳ_1]
      (C) q₃=Q,ε₃=+ ; q₂=q₁, ε₂=−ε₁    → (B) ile aynı
    çakışmalar (üç terim, mertebe a_Q³) tam çıkarılır.
    ÇIPLAK LİMİT: s2^{Ç1} = −½ sin(2πτ_Q) Σ_q a_q² sin²(2πτ_q)
    — **163'ün kapalı formunun ta kendisi** (§2e), burada TÜRETİLEREK
    çıkıyor: Re[h'_1ȳ_1] = b_1B_1cos(πτ_1) = a_1² sin²(2πτ_1) ve
    taşıyıcı Im[h̄_Q y_Q]/(4rm) = −½sin(2πτ_Q)/... (bkz. `c1_capla`).
  • **Ç2 — ν ≡ 0 ARİTMETİKLİ.** q₁^{ε₁}q₂^{ε₂}q₃^{ε₃} = Q'nun
    çift-iptal DIŞI çözümleri. Asal kuvvetlerde bunlar YALNIZ aynı
    asalın kuleleridir (163 §2c teklik lemması); `kule_cozumleri` sayar.
  • **Ç3 — ν ≠ 0, pencere-fazlı.** Geri kalan HER ŞEY. Kapalı biçimde
    Ç3 = TOPLAM − Ç1 − Ç2 (özdeş kapanış), |Δ|-çözünürlüklü dağılımı
    alt-merdiven sayımıyla ölçülür (`165_kanal.py`).

OFF (ara-nokta) frekansında W ladder çizgisi DEĞİLDİR ⇒ Δ = 0 çözümü
yoktur ⇒ **Ç1 = Ç2 = 0 ÖZDEŞ**, off katkısının tamamı Ç3'tür.
"""

import importlib
import math
import sys
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
for _p in ("163_configs", "160_configs", "159_configs"):
    sys.path.insert(0, str(QM / _p))
C163 = importlib.import_module("163_cekirdek")
C160 = importlib.import_module("160_cekirdek")

TWO_PI = 2 * np.pi
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/165")
IZGARA_T1 = C160.IZGARA_T1          # izgara(0.44,0.80,0.04) — 160'ın bantları

# gazların merdiven penceresi (164_insa.PENCERE ile aynı)
PENCERE = {"A4": (0.68, 0.125), "SA4": (0.68, 0.125), "eski_A4": (0.68, 0.125),
           "HA4": (0.68, 0.125), "SA1": (0.75, 0.10)}


# ---------------------------------------------------------------------
def gaz(veri, taban=0.40, cap=4000):
    return C163.gaz(veri, taban, cap)


def merdiven(L, tau_max=0.95, veri=None):
    """163_cekirdek.merdiven + gazın erfc penceresi (varsa)."""
    M = dict(C163.merdiven(L, tau_max))
    pen = PENCERE.get(veri)
    M["pen"] = pen
    if pen is not None:
        w = np.array([0.5 * math.erfc((float(t) - pen[0]) / pen[1])
                      for t in M["tau"]])
        M["w_pen"] = w
        M["a"] = M["a"] * w
        M["b"] = M["b"] * w
        M["B"] = M["B"] * w
    else:
        M["w_pen"] = np.ones_like(M["a"])
    return M


def tayf_s(sites, fields, W, blok=25000, fblok=256):
    """[2⟨f_n e^{−iW s_n}⟩] her f ∈ fields, her W (bellek: fblok×blok)."""
    N = len(sites)
    W = np.asarray(W, float)
    acc = [np.zeros(len(W), dtype=complex) for _ in fields]
    for f0 in range(0, len(W), fblok):
        fs = slice(f0, min(f0 + fblok, len(W)))
        Wc = W[fs]
        for a in range(0, N, blok):
            sl = slice(a, min(a + blok, N))
            arg = np.outer(sites[sl], Wc)
            Cc = np.cos(arg)
            Ss = np.sin(arg)
            del arg
            for i, f in enumerate(fields):
                acc[i][fs] += f[sl] @ Cc - 1j * (f[sl] @ Ss)
            del Cc, Ss
    return [2 * a / N for a in acc]


def sentez(sites, W, amp, blok=25000, fblok=256):
    """Σ_q Re[amp_q e^{iW_q s_n}] — model alanı."""
    N = len(sites)
    out = np.zeros(N)
    ar = np.ascontiguousarray(amp.real)
    ai = np.ascontiguousarray(amp.imag)
    for f0 in range(0, len(W), fblok):
        fs = slice(f0, min(f0 + fblok, len(W)))
        Wc = W[fs]
        for a in range(0, N, blok):
            sl = slice(a, min(a + blok, N))
            arg = np.outer(sites[sl], Wc)
            out[sl] += np.cos(arg) @ ar[fs] - np.sin(arg) @ ai[fs]
            del arg
    return out


def kappa(sites, nu, blok=50000):
    """κ(ν) = ⟨e^{iν s_n}⟩ — TAM pencere/tarak çekirdeği."""
    nu = np.atleast_1d(np.asarray(nu, float))
    N = len(sites)
    acc = np.zeros(len(nu), dtype=complex)
    for a in range(0, N, blok):
        sl = slice(a, min(a + blok, N))
        arg = np.outer(sites[sl], nu)
        acc += np.cos(arg).sum(0) + 1j * np.sin(arg).sum(0)
        del arg
    return acc / N


# ---------------------------------------------------------------------
class Model165:
    """Bir gazın s-sitesindeki çizgi tayfı + model alanları + öngörü.

    ÇİZGİ KESİMİ τ_c (165'in ölçtüğü SINIR). Çizgi (sabit yarım-gap
    genlikli) gösterim yalnız KOHERENT çizgiler için geçerlidir; iki
    bağımsız ölçüt aynı yeri veriyor:
      (i) FAZ SEĞİRMESİ — s_n = s̄ + ḡ(n−n̄) + ḡ Ĉ_n, Ĉ = birikmiş adım
          sapması. e^{iω ḡ Ĉ} açılımı ancak ω ḡ σ_Ĉ ≲ 1 iken bir çizgidir:
          **τ_c = 1/(2π σ_Ĉ)**  (ḡL = 2π).
      (ii) YARIM-GAP MODÜLASYONU — gerçek genlik 2a_q sin(ω_q g_n/2)
          = 2a_q sin(πτ_q(1+ds_n)); sabit b_q = 2a_q sin(πτ_q) yaklaşımı
          πτ_q σ_ds ≲ 1, yani τ_c ≈ 1/(π σ_ds) iken geçerlidir.
    İkisi de bu gazlarda τ_c ≈ 0.55–0.60 veriyor ve BAĞIMSIZ bir sınav
    doğruluyor: model alanının gücü Var(E_mod) ancak τ_c ≈ 0.57'de
    ölçülen Var(η)'yı tutuyor; τ ≤ 0.95'te ÜÇ KAT aşıyor (§tanı).
    """

    def __init__(self, veri, taban=0.40, cap=4000, tau_max=0.95,
                 kaynak="olculen", Y=None, onbellek=True, tau_c=None):
        self.veri, self.taban, self.tau_max = veri, taban, tau_max
        self.kaynak = kaynak
        self.Y = Y if Y is not None else gaz(veri, taban, cap)
        Y = self.Y
        self.s = np.ascontiguousarray(Y.mid[1:])
        self.T = float(self.s[-1] - self.s[0])
        self.dres = TWO_PI / self.T
        self.sbar = float(self.s.mean())
        self.M = merdiven(Y.L, tau_max, veri)
        w = self.M["w"]
        yol = SCR / f"tayf_{veri}_t{taban}_tm{tau_max}.npz"
        if onbellek and yol.exists():
            d = np.load(yol)
            hp, y = d["hp"], d["y"]
        else:
            hp, y = tayf_s(self.s, [Y.e1, Y.Xtil0], w)
            if onbellek:
                SCR.mkdir(parents=True, exist_ok=True)
                np.savez(yol, hp=hp, y=y, w=w)
        self.hp_olc, self.y_olc = hp, y
        # NOMİNAL (sıfır parametre): η regresyonu q ≤ min(e^{taban L},cap)
        # çizgilerini SİLER (154.eta_zinciri: pk(min(exp(taban·L),cap)))
        qlim = min(int(np.exp(taban * Y.L)), cap)
        mask = (self.M["q"] > qlim).astype(float)
        self.hp_nom = self.M["b"] * mask
        self.y_nom = self.M["B"] * np.exp(-1j * np.pi * self.M["tau"])
        self.qlim = qlim
        self.koherans()
        self.tau_c = tau_c if tau_c is not None else tau_max
        self.sec(kaynak)

    # -----------------------------------------------------------------
    def koherans(self):
        """σ_Ĉ (TRENDSİZ) ve iki bağımsız τ_c ölçütü.

        Ĉ_n = Σ_{k<n} X̃0_k = (m_n−m_0)L/2π − n; sabit L kullanıldığı için
        pencere boyunca YAVAŞ bir kayma taşır (log(t/2π) 11.95→12.10).
        Faz-seğirmesi ölçütü yalnız SALINIMLI parçayı ister; bu yüzden
        Ĉ'den n'de kübik trend çıkarılır (163 §6d'nin σ_Ĉ ≈ 0.27'si)."""
        Y = self.Y
        C = np.cumsum(Y.Xtil0)
        n = np.arange(len(C), dtype=float)
        n = (n - n.mean()) / (n[-1] if len(n) > 1 else 1.0)
        V = np.vstack([np.ones_like(n), n, n * n, n ** 3]).T
        c, *_ = np.linalg.lstsq(V, C, rcond=None)
        Cd = C - V @ c
        self.sigC = float(np.std(Cd))
        self.sigds = float(np.std(Y.ds))
        self.tau_c_faz = 1.0 / (TWO_PI * self.sigC) if self.sigC else np.inf
        self.tau_c_gap = 1.0 / (np.pi * self.sigds) if self.sigds else np.inf
        return self

    def sec(self, kaynak=None, tau_c=None):
        if kaynak is not None:
            self.kaynak = kaynak
        if tau_c is not None:
            self.tau_c = tau_c
        m = self.M["tau"] <= self.tau_c + 1e-12
        self.msk = m
        hp = self.hp_olc if self.kaynak == "olculen" else self.hp_nom
        y = self.y_olc if self.kaynak == "olculen" else self.y_nom
        self.hp = np.where(m, hp, 0.0)
        self.y = np.where(m, y, 0.0)
        return self

    def alanlar(self, kmax=3):
        """G_k = E·X^k (E = η_{n+1} modeli, X = X̃0 modeli)."""
        w = self.M["w"][self.msk]
        E = sentez(self.s, w, self.hp[self.msk])
        X = sentez(self.s, w, self.y[self.msk])
        X = X - X.mean()
        self.E, self.X = E, X
        e1, x1 = self.Y.e1 - self.Y.e1.mean(), self.Y.Xtil0
        # TEK-SKALER DEKONVOLÜSYON (yalnız İKİNCİ momentler; üçüncü
        # momentten hiçbir girdi yok): model alanının ölçülen alana en
        # iyi ölçek çarpanı. Çizgi gösterimi çözülmemiş/koherent-olmayan
        # yüksek-τ çizgilerinde alanı ŞİŞİRİYOR; g bunu bir sayıyla
        # düzeltir. BÜTÜN kanallar aynı g_E g_X² ile ölçeklendiği için
        # KANAL PAYLARI g'den BAĞIMSIZDIR.
        gE = float(np.dot(E, e1) / np.dot(E, E)) if np.dot(E, E) else 1.0
        gX = float(np.dot(X, x1) / np.dot(X, X)) if np.dot(X, X) else 1.0
        self.gE, self.gX = gE, gX
        self.gcal = gE * gX * gX
        self.artik = dict(
            varE_mod=float(np.var(E)), varE_olc=float(np.var(e1)),
            varE_res=float(np.var(e1 - E)),
            korE=float(np.corrcoef(E, e1)[0, 1]),
            varX_mod=float(np.var(X)), varX_olc=float(np.var(x1)),
            varX_res=float(np.var(x1 - X)),
            korX=float(np.corrcoef(X, x1)[0, 1]),
            gE=gE, gX=gX, gcal=self.gcal,
            sigC=self.sigC, sigds=self.sigds,
            tau_c_faz=self.tau_c_faz, tau_c_gap=self.tau_c_gap,
            nline=int(self.msk.sum()), tau_c=float(self.tau_c))
        G = [E]
        for _ in range(kmax):
            G.append(G[-1] * X)
        self.G = G
        return self

    def ongor(self, Wler, kmax=3, fblok=256):
        """Her W için J_k = ⟨G_k e^{−iW s}⟩ (TAM κ-ağırlıklı dörtlü toplam)."""
        J = tayf_s(self.s, self.G[:kmax + 1], np.asarray(Wler, float),
                   fblok=fblok)
        return [j / 2.0 for j in J]     # tayf_s 2⟨·⟩ döner; J_k = ⟨·⟩

    # -----------------------------------------------------------------
    def kanal1(self, iQ):
        """Ç1 (çift-iptal, ν≡0 aritmetiksiz) — TAM kapalı toplam.

        iQ: taşıyıcının merdiven indeksi (on-çizgi). J_2^{Ç1} döner.
        """
        hp, y = self.hp, self.y
        SA = 0.25 * hp[iQ] * float(np.sum(np.abs(y) ** 2))      # (A) = T_a
        SB = 0.25 * y[iQ] * float(np.sum((hp * np.conj(y)).real))  # (B)
        cak = (0.125 * hp[iQ] * abs(y[iQ]) ** 2 * 2            # A∩B, A∩C
               + 0.125 * np.conj(hp[iQ]) * y[iQ] ** 2)         # B∩C
        return SA + 2 * SB - cak, SA, 2 * SB, cak

    def kule_cozumleri(self, iQ):
        """Ç2: q₁^{ε₁}q₂^{ε₂}q₃^{ε₃} = Q'nun çift-iptal DIŞI çözümleri.

        Asal kuvvetlerde yalnız aynı asalın kuleleri; ayrıca q₁ η'da
        var olmalı (τ₁ > taban). J_2^{Ç2} döner."""
        q = self.M["q"]
        mult = self.M["mult"]
        Q = int(q[iQ])
        m = int(mult[iQ])
        p = int(round(Q ** (1.0 / m)))
        # aynı asalın merdivendeki kuleleri
        idx = [i for i in range(len(q)) if int(q[i]) == p ** int(mult[i])]
        if not idx:
            return 0j, 0
        hp, y = self.hp, self.y
        tot = 0j
        n = 0
        for i1 in idx:
            for e1 in (+1, -1):
                for i2 in idx:
                    for e2 in (+1, -1):
                        for i3 in idx:
                            for e3 in (+1, -1):
                                k = (e1 * mult[i1] + e2 * mult[i2]
                                     + e3 * mult[i3])
                                if k != m:
                                    continue
                                # çift-iptal mi? (Ç1'de sayıldı mı)
                                c1 = (i1 == iQ and e1 == 1 and i2 == i3
                                      and e2 == -e3)
                                c2 = (i2 == iQ and e2 == 1 and i3 == i1
                                      and e3 == -e1)
                                c3 = (i3 == iQ and e3 == 1 and i2 == i1
                                      and e2 == -e1)
                                if c1 or c2 or c3:
                                    continue
                                v = (hp[i1] if e1 > 0 else np.conj(hp[i1]))
                                v = v * (y[i2] if e2 > 0 else np.conj(y[i2]))
                                v = v * (y[i3] if e3 > 0 else np.conj(y[i3]))
                                tot += v / 8.0
                                n += 1
        return tot, n


# ---------------------------------------------------------------------
def s_den_J(hQ, J, rm):
    """s_k = Im[h̄_Q J_k]/(2⟨ρ⟩) ; u_k = Re[h̄_Q J_k]/(2⟨ρ⟩)."""
    v = np.conj(hQ) * J
    return float(v.imag) / (2 * rm), float(v.real) / (2 * rm)


def c1_capla(tauQ, tau, a, taban=0.0, tau_c=1e9):
    """ÇIPLAK Ç1 (163 §2e): s2 = −½ sin(2πτ_Q)·Σ_{taban<τ≤τ_c} a²sin²(2πτ).

    Toplam η'nın çizgi kümesi üzerindedir: regresyon τ ≤ taban'ı siler."""
    m = (tau > taban + 1e-12) & (tau <= tau_c + 1e-12)
    return -0.5 * np.sin(TWO_PI * tauQ) * float(
        np.sum(a[m] ** 2 * np.sin(TWO_PI * tau[m]) ** 2))
