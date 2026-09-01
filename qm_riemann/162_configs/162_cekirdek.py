"""
162 — GAZIN BELLEĞİ C_n: KİMLİK + GAUSS-TAYFSAL KAPANIŞ. ÖLÇÜM ÇEKİRDEĞİ
========================================================================
KALEM (01 Eylül, akşam):  C_n = Σ_{k<n} ds_k ,  ds_k = S(z_k) − S(z_{k+1})
⇒ C_n = −S(z_n) + sabit.   160'ın D1 kimliği:
    (ρ_n + iσ_n)·e^{+iA(n+1)} = η_{n+1}·e^{−iA·C_n}·K ,  C_n = Σ_{k≤n} X̃_k
Buradan δ'nın TAM indirgemesi (hiçbir yaklaşım yok):
    zc/2 = ⟨η_n     e^{−iW m_n}⟩ ,  zp/2 = ⟨η_{n+1} e^{−iW m_n}⟩
    δ    = arg[zp·conj(zc)·e^{−iA}] ,  m_n = m_0 + (2π/L)(n + C_{n−1})
Yani δ, (η_n , C_n) ORTAK SÜRECİNİN bir fonksiyonelidir; başka girdi yok.

SORU (T2): δ bu sürecin YALNIZ İKİNCİ MOMENTLERİNDEN mi geliyor?
SINAV: aynı oto- ve çapraz-kovaryansa sahip ORTAK-GAUSS vekil üret,
160 zincirini vekile uygula, δ_gerçek/δ_vekil'e bak.

────────────────────────────────────────────────────────────────────────
C'NİN İKİ PARÇASI (tam ayrışım, fit yok)
────────────────────────────────────────────────────────────────────────
    dsΔ_n = (ds_n+ds_{n+1})/2 − ort          (156._bond = 160'ın X_tam'ı)
    r_n   = X̃_n − dsΔ_n                       SÜRÜKLENME (≈ L/Lw_n − 1;
              korel 0.974; gücünün %94'ü en düşük 500 FFT kutusunda)
    u_n   = n + Σ_{k<n} r_k                   sürüklenmiş "kafes" ekseni
    Ĉ_n   = Σ_{k<n} dsΔ_k                     DURAĞAN BELLEK (Var 0.0726,
                                              Ĉ_0 = Ĉ_{son} = 0)
    m_n   = m_0 + (2π/L)(u_n + Ĉ_n)           ← ÖZDEŞ (bu koşuda 0.0e+00)
Var(Σ X̃) = 1.9e+04 iken Var(Ĉ) = 0.073: C'nin görünen devasa varyansı
TÜMÜYLE yoğunluk sürüklenmesidir; ARİTMETİK bellek Ĉ'dir ve odur ki
−S(z_n)'ye eşittir (T1).

────────────────────────────────────────────────────────────────────────
VEKİL REÇETESİ — ve reddedilen iki denemenin kaydı
────────────────────────────────────────────────────────────────────────
(RED-1) Düz n-uzayı FFT faz karıştırması. Ölçüm konum uzayındadır;
  n-uzayında bir asal çizgi SABİT frekansta DEĞİLDİR: sürüklenme
  C̄ = Σr salınımı ≈ 470 birim olduğundan τ=0.45'te çizgi n-tayfında
  ≈ 210 kutuya yayılır. Fazlar karıştırılınca uyumlu toplanma biter.
(RED-2) u ekseninde DÜZGÜN ızgaraya taşı → faz karıştır → geri taşı.
  V0 (faz karıştırmasız) kontrolü ÇALIŞIR (pow/gp 0.628 → 0.619), ama
  vekil τ>0.5 çizgilerini ÖLDÜRÜR (0.628 → 0.008): düzgün u ızgarası
  Nyquist'i τ=0.5'e koyar, oysa ÖLÇÜM τ ≤ 0.80'i düzensiz örneklemenin
  (sürüklenmenin) sayesinde AYIRT EDEBİLİYOR. Bant-sınırlı ara temsil
  bu üstün-Nyquist bilgiyi yok eder. (İkisi de ölçüldü; V-K'da tabloda.)

(KABUL) KONUM-UZAYI FAZ KARIŞTIRMASI. İkinci-moment yapısı konum
uzayında şudur: her asal frekans w_q = log q'daki çizgi genliği |c_q|,
çizgiler-arası çapraz-tayf, ve artık sürekliliğin tayfı. Vekil:

  1) ÇIKARIM: ölçümün gördüğü BÜTÜN çizgiler için
        c_q(x) = 2⟨x_n e^{−i w_q m_n}⟩ ,  q ∈ pk_m(e^{0.86L}) (160 ile aynı küme)
     x ∈ {η, Ĉ}.  c_q, ÖLÇÜMÜN o çizgide gördüğü genliğin ta kendisidir
     (zc = c_q); bu yüzden ham izdüşüm kullanılır (NITER = 0) — vekil
     böylece pow_on'u BİREBİR yeniden üretir. VARYANS BÜTÇESİ ölçülür ve
     rapora girer: Ĉ için Σ|c_q|²/2 = 0.0741 ↔ Var(Ĉ) = 0.0726 (%+2,
     kapanıyor), η için 0.0699 ↔ 0.0543 (%+29, KAPANMIYOR). Yani η'nın
     çizgileri gerçekte YIKICI GİRİŞİM yapıyor; rastgele fazlı vekilde
     bu girişim yok ve vekil gaz "daha gürültülü" oluyor. Bu bir seçim
     değil, ölçülen bir olgudur ve bütün oranlarda açıkça taşınır
     (162_dogrulama V5 bunun bir Gram artefaktı OLMADIĞINI sınar).
  2) ÇİZGİ/SÜREKLİLİK AYRIMI: x_çizgi = Re Σ_q c_q e^{i w_q m_n},
     x_sür = x − x_çizgi.
  3) VEKİL:  c'_q = c_q·e^{iψ_q}  (ψ_q ~ U(0,2π), η ve Ĉ için AYNI ψ_q
     ⇒ çizgi çapraz-tayfı TAM korunur, çizgi güçleri TAM korunur)
     ve x'_sür = n-uzayı ORTAK faz karıştırması (Prichard–Theiler:
     |F_k| korunur, aynı rastgele faz η ve Ĉ artıklarına; k=0 ve
     Nyquist fazı korunur ⇒ ortalama korunur).
     x' = x'_sür + Re Σ_q c'_q e^{i w_q m_n}
  4) ÖZ-TUTARLI YENİDEN KURULUM. Ĉ bir ALANDIR (Ĉ_n = −S(z_n), S konumun
     fonksiyonu) ve konumlar Ĉ'den doğar; vekil gaz AYNI sabit noktayı
     çözer (sürüklenme r ve u ekseni GERÇEĞİN kendisi, dokunulmaz):
        Ĉ'_n = Ĉ'_sür,n + Re Σ_q c'_q e^{i w_q m'_n}
        m'_n = m_0 + (2π/L)(u_n + Ĉ'_n − Ĉ'_0)
        dsΔ' = diff(Ĉ') ,  X̃' = r + dsΔ' ,  η'_n = η'_sür,n + ReΣ c'^η_q e^{iw_q m'_n}
     Tutarlı: X̃'_n ≡ (m'_{n+1}−m'_n)L/2π − 1 = r_n + dsΔ'_n. V0'da m' = m
     TAM sabit noktadır (ölçüldü: 0.0e+00), yani boru hattı kimliği geçer.

(RED-3) Öz-tutarlılık OLMADAN (çizgiler gerçeğin m'sinde sentezlenip
  vekilin m'siyle okunursa) çizgi gücü pow/gp 0.73 → 0.20, 0.63 → 0.09,
  0.47 → 0.04'e düşer ve τ≥0.6'da pow_on < pow_off olur: ölçüm kırılır.
  Sebebi fiziksel — alan, gazın KENDİ konumlarında okunmalıdır.

Çok sayıda (≈3400) rastgele fazlı çizgi + faz karıştırılmış süreklilik
⇒ süreç Gauss'tur (Rice); ikinci momentler korunur; korunmayan tek şey
FAZ İLİŞKİLERİDİR — yani "aritmetik bellek".

ÇEŞİTLER (aynı boru hattı, ne karıştırıldığı değişir):
    V0  ψ≡0, süreklilik dokunulmaz  → gerçeğin ta kendisi (kimlik kontrolü)
    VG  çizgi + süreklilik           → ANA Gauss vekil
    VL  yalnız çizgi fazları         → "çizgiler-arası kilit" tanısı
    VC  yalnız süreklilik            → süreklilik tanısı

ÖLÇÜM ZİNCİRİ KOPYALANMAZ: `olc162`, 160_cekirdek.olc160'ı OLDUĞU GİBİ
çağırır; yalnız `Yerel160` sınıfı geçici olarak hazır nesneyi döndüren
bir fabrikayla değiştirilir ⇒ aday listesi, 220-örnekleme tohumu (21),
gap filtresi, on/off ara-nokta referansı ve 8-grup jackknife BİT
DÜZEYİNDE 160'ınkidir (162_dogrulama V1 bunu sınar).
"""

import copy
import importlib
import sys
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
for _p in ("160_configs", "159_configs", "156_configs"):
    sys.path.insert(0, str(QM / _p))
C160 = importlib.import_module("160_cekirdek")
C159 = importlib.import_module("159_cekirdek")
K156 = importlib.import_module("156_cekirdek")
K155, KOS155, C154 = C160.K155, C160.KOS155, C160.C154

TWO_PI = 2 * np.pi
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/162")

IZGARA158 = C160.IZGARA158
IZGARA_T1 = C160.IZGARA_T1
KANALLAR = C160.KANALLAR
TAU_CIZGI = 0.86          # çizgi evreni: 160'ın aday listesiyle AYNI
BLOK = 2000               # trig blok boyu
NITER = 0                 # Jacobi arındırma adımı sayısı


# --- konum-uzayı çizgi cebri -----------------------------------------
def cizgi_cikar(m, diziler, w, blok=BLOK):
    """c_q(x) = 2⟨x_n e^{−i w_q m_n}⟩ , her x için (tek trig geçişi)."""
    N = len(m)
    re = [np.zeros(len(w)) for _ in diziler]
    im = [np.zeros(len(w)) for _ in diziler]
    for s0 in range(0, N, blok):
        sl = slice(s0, min(s0 + blok, N))
        A = np.outer(m[sl], w)
        Cc = np.cos(A)
        Ss = np.sin(A)
        for k, x in enumerate(diziler):
            re[k] += Cc.T @ x[sl]
            im[k] -= Ss.T @ x[sl]
        del A, Cc, Ss
    return [2.0 * (r + 1j * i) / N for r, i in zip(re, im)]


def cizgi_kur(m, katsayilar, w, blok=BLOK):
    """x_çizgi = Re Σ_q c_q e^{i w_q m_n} = Σ_q |c_q| cos(w_q m_n + arg c_q).
    (Genlik/faz biçimi tek transandantal çağrı kullanır — iki kat hızlı.)"""
    N = len(m)
    gen = [np.abs(c) for c in katsayilar]
    faz = [np.angle(c) for c in katsayilar]
    out = [np.empty(N) for _ in katsayilar]
    for s0 in range(0, N, blok):
        sl = slice(s0, min(s0 + blok, N))
        A0 = np.outer(m[sl], w)
        for k in range(len(katsayilar)):
            out[k][sl] = np.cos(A0 + faz[k]) @ gen[k]
        del A0
    return out


def cizgi_kur_d(m, c, w, blok=BLOK):
    """x_çizgi ve d/dm türevi (Newton için; tek trig geçişinde)."""
    N = len(m)
    g = np.abs(c)
    f = np.angle(c)
    gw = g * w
    y = np.empty(N)
    yp = np.empty(N)
    for s0 in range(0, N, blok):
        sl = slice(s0, min(s0 + blok, N))
        A = np.outer(m[sl], w) + f
        y[sl] = np.cos(A) @ g
        yp[sl] = -(np.sin(A) @ gw)
        del A
    return y, yp


def faz_karistir(diziler, tohum):
    """n-uzayı ORTAK faz karıştırma. |F_k| ve ortalama korunur.
    tohum None ⇒ ψ ≡ 0 (dizi değişmez)."""
    if tohum is None:
        return [x.copy() for x in diziler]
    M = len(diziler[0])
    rng = np.random.default_rng(tohum)
    ph = rng.uniform(0.0, TWO_PI, M // 2 + 1)
    ph[0] = 0.0
    if M % 2 == 0:
        ph[-1] = 0.0
    ps = np.exp(1j * ph)
    return [np.fft.irfft(np.fft.rfft(x) * ps, n=M) for x in diziler]


# ---------------------------------------------------------------------
class Taban162:
    """Bir (veri, taban) için GERÇEK Yerel160 + konum-uzayı vekil donanımı."""

    def __init__(self, veri, taban=0.40, cap=4000, tau_cizgi=TAU_CIZGI,
                 ayrinti=True):
        self.veri, self.taban, self.cap = veri, taban, cap
        self.z = KOS155.veri_yukle(veri)
        self.Y = C160.Yerel160(self.z, veri, taban, cap)
        Y = self.Y
        self.Nn = len(Y.eta)
        dsB = Y.X["tam"]                       # = _bond(ds), uzunluk Nn−1
        self.dsB = dsB
        self.r = Y.Xtil - dsB                  # SÜRÜKLENME (tam ayrışım)
        self.Chat = np.concatenate(([0.0], np.cumsum(dsB)))   # DURAĞAN bellek
        self.u = np.concatenate(([0.0],
                                 np.arange(1, self.Nn) + np.cumsum(self.r)))
        self.m0_ilk = float(Y.mid[0])
        self.dv = TWO_PI / Y.L
        self.m_kimlik = float(np.max(np.abs(
            Y.mid - (self.m0_ilk + self.dv * (self.u + self.Chat)))))
        # çizgi evreni (160'ın aday kümesiyle aynı üretici)
        qm = C154.pk_m(int(np.exp(tau_cizgi * Y.L)))
        self.q = np.array(sorted(qm))
        self.w = np.log(self.q.astype(float))
        self.mq = np.array([qm[int(x)] for x in self.q], dtype=float)
        # çizgi katsayıları (Jacobi arındırmalı) ve çizgi/süreklilik ayrımı
        hedef = [Y.eta, self.Chat]
        c = cizgi_cikar(Y.mid, hedef, self.w)
        self.yakinsama = []
        for it in range(NITER + 1):
            L_ = cizgi_kur(Y.mid, c, self.w)
            art = [h - l for h, l in zip(hedef, L_)]
            gc = [float(np.sum(np.abs(ci) ** 2) / 2) for ci in c]
            va = [float(np.var(a)) for a in art]
            self.yakinsama.append(
                [(g + v) / float(np.var(h)) for g, v, h in zip(gc, va, hedef)])
            if it == NITER:
                break
            d = cizgi_cikar(Y.mid, art, self.w)
            c = [ci + di for ci, di in zip(c, d)]
        self.c_eta, self.c_C = c
        self.eta_sur, self.C_sur = art
        self.pay_eta = float(1.0 - np.var(self.eta_sur) / np.var(Y.eta))
        self.pay_C = float(1.0 - np.var(self.C_sur) / np.var(self.Chat))
        self.kapanis = self.yakinsama[-1]
        if ayrinti:
            print(f"  [162 taban] {veri} t{taban}: Nn={self.Nn} L={Y.L:.4f} "
                  f"m-kimlik={self.m_kimlik:.1e} Var(Ĉ)={self.Chat.var():.5f} "
                  f"Var(ΣX̃)={np.cumsum(Y.Xtil).var():.4g} çizgi={len(self.q)}",
                  flush=True)
            print(f"      çizgi-payı η={self.pay_eta:+.4f} Ĉ={self.pay_C:+.4f}"
                  f"   varyans kapanışı (Σ|c|²/2+Var_sür)/Var: "
                  f"η={self.kapanis[0]:.4f} Ĉ={self.kapanis[1]:.4f}"
                  f"   [it0: {self.yakinsama[0][0]:.3f}/"
                  f"{self.yakinsama[0][1]:.3f}]", flush=True)

    # -----------------------------------------------------------------
    def vekil_yerel(self, tohum, cesit="VG", nsabit=6, esik=1e-9,
                    ayrinti=True):
        """cesit: V0 | VG | VL | VC.  tohum=None ⇒ V0.

        ÖZ-TUTARLILIK: Ĉ, konumların FONKSİYONU olan bir alandır
        (Ĉ_n = −S(z_n)) ve konumlar Ĉ'den doğar. Vekil gaz da aynı
        sabit-nokta denklemini çözer:
            Ĉ'_n = Ĉ'_sür,n + Re Σ_q c'_q e^{i w_q m'_n}
            m'_n = m_0 + (2π/L)(u_n + Ĉ'_n − Ĉ'_0)
        Büzülme çarpanı ≈ |dsΔ| (rms 0.23) ⇒ birkaç adımda yakınsar.
        V0'da m' = m TAM olarak sabit noktadır (m-kimlik = 0).
        """
        if cesit == "V0" or tohum is None:
            cesit = "V0"
        rng = np.random.default_rng(tohum) if tohum is not None else None
        if cesit in ("VG", "VL"):
            rot = np.exp(1j * rng.uniform(0.0, TWO_PI, len(self.w)))
        else:
            rot = np.ones(len(self.w))
        st = (tohum + 500000) if (cesit in ("VG", "VC")) else None
        es, cs = faz_karistir([self.eta_sur, self.C_sur], st)
        ce, cc = self.c_eta * rot, self.c_C * rot
        # NEWTON: h(m) = m − m_0 − (2π/L)(u + Ĉ_sür + F(m)) = 0
        # h'(m) = 1 − (2π/L)F'(m). Adım sınırı 0.3·(2π/L) — 152/154'ün
        # Newton çözücüsünün yaptığı gibi yanlış köke atlamayı engeller.
        sab = self.m0_ilk + self.dv * (self.u + cs - cs[0])
        m = self.Y.mid.copy()
        adim = 0.3 * self.dv
        iz = []
        for it in range(nsabit):
            F, Fp = cizgi_kur_d(m, cc, self.w)
            h = m - (sab + self.dv * (F - F[0]))
            iz.append((float(np.max(np.abs(h))), float(np.std(h))))
            if iz[-1][1] < esik:
                break
            hp = 1.0 - self.dv * Fp
            hp = np.where(np.abs(hp) < 0.25, np.sign(hp) * 0.25 + (hp == 0),
                          hp)
            m = m + np.clip(-h / hp, -adim, adim)
        Ch = cs + cizgi_kur(m, [cc], self.w)[0]
        eta = es + cizgi_kur(m, [ce], self.w)[0]
        if ayrinti:
            print(f"    [{cesit} t={tohum}] Newton {len(iz)} adım, "
                  f"|h| maks {iz[0][0]:.2e}→{iz[-1][0]:.2e}, "
                  f"rms {iz[0][1]:.2e}→{iz[-1][1]:.2e}", flush=True)
        Y = self._kur(eta, Ch, tohum, cesit)
        Y._162["fix_iz"] = iz
        return Y

    def _kur(self, eta, Ch, tohum, cesit):
        Y = copy.copy(self.Y)
        Y.eta = eta
        dsB = np.diff(Ch)
        Y.ds = np.concatenate((dsB, dsB[-1:]))      # yalnız kayıt için
        Xtil = self.r + dsB
        mid = self.m0_ilk + self.dv * (self.u + (Ch - Ch[0]))
        Y.mid = mid
        Y.Chat = Ch
        sif = np.zeros_like(dsB)
        Y.X = {"tam": dsB, "lad": sif, "eta": sif, "dri": sif}
        Y.sA2 = float(np.var(dsB))
        Y.kov = {k: float(np.mean(v * dsB)) for k, v in Y.X.items()}
        Y.Xtil = Xtil
        Y.Xtil0 = Xtil - Xtil.mean()
        Y.Xort = float(Xtil.mean())
        Y.dXtil = float(np.max(np.abs(Y.Xtil0 - dsB)))
        Y.e0, Y.e1 = eta[:-1], eta[1:]
        Y.m0 = mid[:-1]
        Y.N = len(Y.m0)
        X0 = Y.Xtil0
        Y.P1_ = X0
        Y.P2_ = X0 * X0
        Y.P3_ = Y.P2_ * X0
        Y.P4_ = Y.P2_ * Y.P2_
        Y.dX = np.concatenate((np.diff(X0), np.zeros(1)))
        Y.s_ds = float(np.var(Y.ds))
        Y.s_eta = float(np.var(eta))
        Y.c1 = float(np.mean(eta[:-1] * eta[1:]))
        Y.s_lad = Y.s_dri = 0.0
        Y.eta_kapanis = float("nan")
        Y.bond_kapanis = float("nan")
        Y.dil = {}
        for W in C159.PENCERELER:
            o = (W - 1) // 2
            n = Y.N - W + 1
            xw, _ = C159.kayan(dsB, W)
            Y.dil[W] = dict(o=o, n=n,
                            tam=np.ascontiguousarray(dsB[o:o + n]),
                            xw=xw, Ptot=C159.kayan(Y.e0 ** 2, W)[0])
        Y._162 = dict(tohum=tohum, cesit=cesit)
        return Y


class _Fabrika:
    def __init__(self, Y):
        self.Y = Y

    def __call__(self, *a, **k):
        return self.Y


def olc162(Y, etiket, bantlar, **kw):
    """160'ın olc160'ını AYNEN koşar; Yerel160 yerine hazır Y konur."""
    eski = C160.Yerel160
    C160.Yerel160 = _Fabrika(Y)
    try:
        return C160.olc160(np.zeros(len(Y.mid) + 1), etiket, bantlar,
                           anahtar="-", **kw)
    finally:
        C160.Yerel160 = eski


# --- merdiven (T1) ----------------------------------------------------
def merdiven(t, L, tau_max=0.9, blok=4000):
    """S_merdiven(t) = −Σ_{q=p^m, log q ≤ τ_max·L} a_q sin(t log q),
       a_q = 1/(π m p^{m/2}) = 1/(π m √q)   (155/160'ın aq'suyla AYNI)."""
    qm = C154.pk_m(int(np.exp(tau_max * L)))
    qs = np.array(sorted(qm), dtype=float)
    ms = np.array([qm[int(q)] for q in qs], dtype=float)
    w = np.log(qs)
    aq = 1.0 / (np.pi * ms * np.sqrt(qs))
    S = np.zeros(len(t))
    for s0 in range(0, len(t), blok):
        sl = slice(s0, min(s0 + blok, len(t)))
        S[sl] = np.sin(np.outer(t[sl], w)) @ aq
    return -S, dict(nq=len(qs), w=w, aq=aq, tau=w / L, q=qs, m=ms)
