"""
163 — KİLİTLİ-FAZ İSTATİSTİĞİNİN İLK-İLKE SINAVI: ÖLÇÜM + ÜÇLÜ-TOPLAM ÇEKİRDEĞİ
==============================================================================
Hiçbir ölçüm parçası kopyalanmaz: `160_cekirdek` (o da 159/155/156/154)
import edilir; `Yerel160` aynen kullanılır. Bant/çizgi SEÇİMİ (aday listesi,
220-örnekleme tohumu 21, `gap < 2.5·dres` filtresi, on/off ara-nokta
referansı `W' = w + gap/2`, 8-grup round-robin) `olc160`'ın döngüsünden
alınmıştır ve `163_dogrulama` V1 bunu 160'ın kayıtlı JSON'una karşı
BİT DÜZEYİNDE sınar (aynı bantlarda maks fark 0.0e+00 beklenir).

══════════════════════════════════════════════════════════════════════════
1. BELİRLENİMLİ MERDİVEN → ÇİZGİ GENLİKLERİ (yarım-gap çarpanları)
══════════════════════════════════════════════════════════════════════════
Sentetik gaz  N_sm(z) + S(z) = n,  S(t) = −Σ_q a_q sin(ω_q t),
ω_q = log q,  a_q = 1/(π m √q)·(pencere).  N_sm' = 1/ḡ, ḡ = 2π/L ⇒

    ds_n = −[S(z_{n+1})−S(z_n)] = Σ_q 2a_q sin(ω_q g_n/2)·cos(ω_q m_n)

  • η çizgi genliği      b_q = 2a_q sin(πτ_q)          (kodun `gp` = b_q²)
  • X̃ = dsΔ çizgi gen.   B_q = b_q cos(πτ_q) = a_q sin(2πτ_q),  faz +πτ_q
  • Ĉ = Σ dsΔ            a_q cos(πτ_q)  ⇒ 162 §2b'nin cos²(πτ) sütunu

Bu üç YARIM-GAP çarpanı (sin πτ, cos πτ) tek kaynaktan gelir:
sin/cos(ω ḡ/2) = sin/cos(πτ).  143'ün G formülündeki cos(πτ) ile aynı kalıp.

══════════════════════════════════════════════════════════════════════════
2. ⟨σX̃^k⟩'nin ÜÇLÜ-TOPLAMI — TAM TÜRETİM (yaklaşım yalnız 1. mertebe)
══════════════════════════════════════════════════════════════════════════
Tanım (160 §2, hiçbir yaklaşım yok):
    c_n = η_n e^{−iω_Q m_n},  ⟨c⟩ = h_Q/2  (h_Q = zc: ÖLÇÜLEN çizgi genliği)
    ρ_n + iσ_n = c_{n+1}·conj⟨c⟩ ,  ⟨ρ⟩ = |h_Q|²/4
    s_k = ⟨σX̃0^k⟩/⟨ρ⟩ = 2·Im[h̄_Q J_k]/|h_Q|² ,  J_k ≡ ⟨c_{n+1} X̃0_n^k⟩
    u_k = ⟨ρX̃0^k⟩/⟨ρ⟩ = 2·Re[h̄_Q J_k]/|h_Q|²      (u1 = S = 159/160'ın A·S'i)

η'yı KENDİ çizgisi + gerisi diye ayır (site m_{n+1}):
    η_{n+1} = Re[h_Q e^{iω_Q m_{n+1}}] + η^r_{n+1}
    c_{n+1} = h_Q/2 + (h̄_Q/2)e^{−2iω_Q m_{n+1}} + η^r_{n+1}e^{−iω_Q m_{n+1}}

(T_a) KENDİ ÇİZGİSİ:  (h_Q/2)·⟨X̃0^k⟩ — REEL katsayı ⇒ **s_k'ya KATKISI SIFIR**,
      u_k'ya ⟨X̃0^k⟩ verir (u0 = 1 ÖZDEŞ). Gauss/köşegen pay burada ölür:
      kuadratür kanalı kendi çizgisinden BESLENMEZ.
(T_b) 2ω_Q harmoniği: yalnız (q₂,q₃) = (Q,Q) çifti besler ⇒ x_Q² mertebesi.
(T_c) GERİSİ: q₁ ≠ Q çizgileri. m_{n+1} = m_n + ḡ(1+X̃_n) ÖZDEŞ olduğundan
      her q₁ terimi ⟨e^{iν m_n}·e^{iA_ν(1+X̃_n)}X̃0_n^k⟩ biçimindedir,
      ν = ω₁−ω_Q (fark) ya da −(ω₁+ω_Q) (toplam), A_ν = 2π ν/L.
      Tarak-ortalaması X̃0^k e^{iA_ν X̃}'nin −ν bileşenini ister; X̃'nin
      çizgileri x_q olduğundan bu bileşen ancak

          −ν = ε₂ω₂ + ε₃ω₃  ⟺  **q₁ = Q·q₂^{ε₂}q₃^{ε₃}**   (ÇARPIMSAL ÜÇLÜ)

      ile yaşar. Q = p^a, q₁ = r^b (r ≠ p) iken asal-kuvvet çözümü TEKTİR:
      **(q₂,q₃) = (Q, q₁)** — EVRENSEL ÜÇLÜ. (Aynı asalın kuleleri ek çözüm
      verir: §3.) Bileşenin katsayısı X̃'nin İKİNCİ türevidir:
          ⟨e^{iνm}·G(X̃)⟩ = ¼·⟨G''⟩·x_Q·x̄₁   (fark) ya da ¼⟨G''⟩x_Q x₁ (toplam)
      — köşegen (q₂=q₃) "dressing" ⟨G''⟩'nün İÇİNDEDİR (tam X̃ ile ölçülür).

SONUÇ (ÜÇLÜ-TOPLAM):

    Γ(θ,k) ≡ ⟨ d²/dX²[ e^{iθ(1+X)} X0^k ] ⟩
           = e^{iθ}[ k(k−1)Φ_{k−2}(θ) + 2ikθ Φ_{k−1}(θ) − θ² Φ_k(θ) ]
    Φ_j(θ) ≡ ⟨ X0^j e^{iθ X̃} ⟩          (X̃'nin MARJİNALİ; ölçülür)

    J_k = (h_Q/2)Φ_k(0)                                         … T_a
        + (h̄_Q x_Q²/8)·Γ(−4πτ_Q, k)                             … T_b
        + (x_Q/8)·Σ_{q₁≠Q}[ h₁x̄₁·Γ(2π(τ₁−τ_Q),k)
                          + h̄₁x₁·Γ(−2π(τ₁+τ_Q),k) ]            … T_c

ÇIPLAK LİMİT (Φ_0→1, Φ_1→0, Φ_2→σ², θ²σ² ihmal; h=b, x=Be^{iπτ}):

    Γ(θ,2) → 2e^{iθ}  ⇒  J_2 = (b_Q/2)σ² + (x_Q/2)e^{−2πiτ_Q}Σ b₁B₁cos(πτ₁)
    **s2 → −cos(πτ_Q)·sin(πτ_Q)·Σ_{q₁} b₁B₁cos(πτ₁)
          = −½ sin(2πτ_Q) · Σ_{q₁} a₁² sin²(2πτ₁)**

— YARIM-GAP çarpanlarının işaret muhasebesi: −sin(πτ_Q)·cos(πτ_Q) taşıyıcıdan
(T4'ün cevabı), Σ b₁B₁cos(πτ₁) = Σ4a²sin²(πτ)cos²(πτ) çıplakta DAİMA ≥ 0.
Köşegen (q₁=Q) terimi sin(π(τ₁−τ_Q))|_{q₁=Q} = 0 ile TAM ölür ⇒ ikinci
momentlerden gelen pay YOKTUR (162'nin Gauss-ötesi hükmünün kalem karşılığı).

══════════════════════════════════════════════════════════════════════════
3. KULE (aynı asal) ve TARAK-ARACILI (G) katmanları
══════════════════════════════════════════════════════════════════════════
Q = p^a, q₁ = p^c iken q₁ = Q q₂^{ε₂}q₃^{ε₃} denkleminin ek çözümleri
(q₂,q₃) = (p^i, p^j) kuleleridir; `kule_ciftleri` bunları sayar.
Tarak-aracılı katman: ⟨e^{iνm}⟩ = 𝒢(ν) = −πτ_ν a_ν cos(πτ_ν)·1{ν = ±ω_q}
(143'ün G_cc'sinin birinci terimi; burada γ_ν = Ĉ'nin çizgi genliğiyle
𝒢(ν) = iπτ_ν conj(γ_ν)). Beş asal-kuvvetli kısıt Q q₁^{±}q₂^{±}q₃^{±}q₄^{±}=1
genel Q için yine yalnız KULE çözümleri verir; büyüklüğü §T2'de ölçülür.
"""

import importlib
import sys
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
for _p in ("160_configs", "159_configs"):
    sys.path.insert(0, str(QM / _p))
C160 = importlib.import_module("160_cekirdek")

TWO_PI = 2 * np.pi
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/163")
IZGARA_T1 = C160.IZGARA_T1
IZGARA158 = C160.IZGARA158


def gaz(veri, taban=0.40, cap=4000):
    z = C160.KOS155.veri_yukle(veri)
    return C160.Yerel160(z, veri, taban, cap)


def merdiven(L, tau_max=1.00):
    qm = C160.C154.pk_m(int(np.exp(tau_max * L)))
    q = np.array(sorted(qm), dtype=np.int64)
    mult = np.array([qm[int(x)] for x in q], dtype=float)
    w = np.log(q.astype(float))
    tau = w / L
    a = 1.0 / (np.pi * mult * np.sqrt(q.astype(float)))
    return dict(q=q, mult=mult, w=w, tau=tau, a=a,
                b=2 * a * np.sin(np.pi * tau),
                B=a * np.sin(2 * np.pi * tau), qm=qm)


def tayf(Y, W, blok=20000, fblok=400):
    """h_q = 2⟨η_n e^{−iωm_n}⟩ (= zc) , x_q = 2⟨X̃0_n e^{−iωm_n}⟩ ,
       hp_q = 2⟨η_{n+1} e^{−iωm_n}⟩ (tanı için).  Bellek: fblok×blok."""
    m, e0, e1, X0 = Y.m0, Y.e0, Y.e1, Y.Xtil0
    N = len(m)
    W = np.asarray(W, float)
    acc = np.zeros((3, len(W)), dtype=complex)
    for f0 in range(0, len(W), fblok):
        fs = slice(f0, min(f0 + fblok, len(W)))
        Wc = W[fs]
        for s in range(0, N, blok):
            sl = slice(s, min(s + blok, N))
            arg = np.outer(m[sl], Wc)
            C = np.cos(arg)
            S = np.sin(arg)
            del arg
            acc[0][fs] += e0[sl] @ C - 1j * (e0[sl] @ S)
            acc[1][fs] += X0[sl] @ C - 1j * (X0[sl] @ S)
            acc[2][fs] += e1[sl] @ C - 1j * (e1[sl] @ S)
            del C, S
    a = 2 * acc / N
    return a[0], a[1], a[2]


# ---------------------------------------------------------------------
class PhiTablo:
    """Φ_j(θ) = ⟨X̃0^j e^{iθX̃}⟩ — X̃'nin MARJİNALİNDEN, histogram + ızgara.

    Binleme hatası (Δθ)²/24 ≤ 5e−6 (Δ = menzil/nbin, |θ| ≤ 12).
    """

    def __init__(self, Y, jmax=4, nbin=6000, thmin=-13.0, thmax=6.0,
                 dth=0.001):
        X = np.asarray(Y.Xtil, float)
        X0 = X - X.mean()
        lo, hi = X.min(), X.max()
        pad = 1e-9 * max(1.0, hi - lo)
        idx = np.clip(((X - lo) / (hi - lo + pad) * nbin).astype(int), 0, nbin - 1)
        cx = lo + (np.arange(nbin) + 0.5) * (hi - lo + pad) / nbin
        cx0 = cx - X.mean()
        self.W = np.zeros((jmax + 1, nbin))
        for j in range(jmax + 1):
            self.W[j] = np.bincount(idx, weights=X0 ** j, minlength=nbin) / len(X)
        self.cx = cx
        self.th = np.arange(thmin, thmax + dth / 2, dth)
        E = np.exp(1j * np.outer(self.th, cx))
        self.tab = (self.W @ E.T)            # (jmax+1, nth)
        self.dth, self.thmin, self.nth = dth, thmin, len(self.th)
        # kesin kontrol: Φ_j(0) = ⟨X0^j⟩
        self.mom = np.array([float(np.mean(X0 ** j)) for j in range(jmax + 1)])
        i0 = int(round((0.0 - thmin) / dth))
        self.hata0 = float(np.max(np.abs(self.tab[:, i0].real - self.mom)))

    def Phi(self, j, th):
        th = np.asarray(th, float)
        f = (th - self.thmin) / self.dth
        i = np.clip(f.astype(int), 0, self.nth - 2)
        w = f - i
        return self.tab[j][i] * (1 - w) + self.tab[j][i + 1] * w

    def Gam(self, th, k):
        """Γ(θ,k) = e^{iθ}[k(k−1)Φ_{k−2} + 2ikθΦ_{k−1} − θ²Φ_k]"""
        th = np.asarray(th, float)
        t = -(th ** 2) * self.Phi(k, th)
        if k >= 1:
            t = t + 2j * k * th * self.Phi(k - 1, th)
        if k >= 2:
            t = t + k * (k - 1) * self.Phi(k - 2, th)
        return np.exp(1j * th) * t


# ---------------------------------------------------------------------
def kule_ciftleri(Q, q1, qm, L):
    """q₁ = Q·q₂^{ε₂}q₃^{ε₃} denkleminin EVRENSEL DIŞI (kule) çözümleri.

    Döner: [(q2, q3, e2, e3), ...] — evrensel (Q,q₁,+,−)/(q₁,Q,−,+) ve
    toplam kanadı (Q,q₁,+,+)/(q₁,Q,+,+) AYRI ele alınır; burada yalnız
    aynı asalın kuleleri (Q ve q₁ aynı p'nin kuvvetiyse) sayılır.
    """
    from sympy import factorint
    fq, f1 = factorint(int(Q)), factorint(int(q1))
    if len(fq) != 1 or len(f1) != 1:
        return []
    (p, a), = fq.items()
    (r, c), = f1.items()
    if p != r:
        return []
    out = []
    kmax = int(np.log(max(qm)) / np.log(p)) + 1
    for i in range(1, kmax + 1):
        for j in range(1, kmax + 1):
            for e2 in (+1, -1):
                for e3 in (+1, -1):
                    if c == a + e2 * i + e3 * j:
                        q2, q3 = p ** i, p ** j
                        if q2 in qm and q3 in qm:
                            if not (q2 == Q and q3 == q1 and (e2, e3) == (1, -1)) \
                               and not (q2 == q1 and q3 == Q and (e2, e3) == (-1, 1)):
                                out.append((q2, q3, e2, e3))
    return out


# ---------------------------------------------------------------------
def bant_adaylari(Y, bantlar, tohum=21, njack=8):
    """160'ın aday seçimi — BİREBİR (V1 bit düzeyinde sınar)."""
    L = Y.L
    qm = C160.C154.pk_m(int(np.exp(0.86 * L)))
    allq = sorted(qm)
    allw = np.array([np.log(q) for q in allq])
    dres = TWO_PI / (Y.mid[-1] - Y.mid[0])
    rng = np.random.default_rng(tohum)
    cikti = []
    for lo, hi in bantlar:
        tumu = [(q, w) for q, w in zip(allq, allw) if lo < w / L <= hi]
        cand = tumu
        if len(cand) > 220:
            idx = rng.choice(len(cand), 220, replace=False)
            cand = [cand[i] for i in idx]
        Ls, kul = [], 0
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
            Ls.append(dict(q=int(q), w=float(w), tau=float(w / L), grup=gi,
                           gap=float(gap),
                           gp=float((2 * aq * np.sin(np.pi * w / L)) ** 2)))
        cikti.append(dict(lo=lo, hi=hi, tau=round(0.5 * (lo + hi), 4),
                          N=len(tumu), kul=kul, cizgi=Ls))
    return cikti


# ---------------------------------------------------------------------
def olc_cizgi(Y, Wf, kmax=3):
    """160'ın çizgi cebri (ρ, σ, u_k, s_k) + h_Q, x_Q."""
    m0 = Y.m0
    cw, sw = np.cos(Wf * m0), np.sin(Wf * m0)
    zc = complex(2 * np.mean(Y.e0 * cw), -2 * np.mean(Y.e0 * sw))
    A = TWO_PI * Wf / Y.L
    cw1, sw1 = np.cos(Wf * Y.mid[1:]), np.sin(Wf * Y.mid[1:])
    c1r, c1i = Y.e1 * cw1, -Y.e1 * sw1
    zb = zc / 2.0
    rho = c1r * zb.real + c1i * zb.imag
    sig = c1i * zb.real - c1r * zb.imag
    rm = float(rho.mean())
    N = len(rho)
    o = dict(A=float(A), pow=float(abs(zc) ** 2), h=zc, rho_ort=rm,
             x=complex(2 * np.mean(Y.Xtil0 * cw), -2 * np.mean(Y.Xtil0 * sw)),
             hp=complex(2 * np.mean(Y.e1 * cw), -2 * np.mean(Y.e1 * sw)))
    P = np.ones(N)
    for k in range(0, kmax + 1):
        if k:
            P = P * Y.Xtil0
        o[f"u{k}"] = float(np.dot(rho, P) / N / rm)
        o[f"s{k}"] = float(np.dot(sig, P) / N / rm)
    return o


def ongor_cizgi(hQ, xQ, tauQ, tau1, h1, x1, PT, kmax=3, kule=None):
    """ÜÇLÜ-TOPLAM öngörüsü: J_k → (u_k, s_k). tau1/h1/x1 = q₁ ≠ Q kümesi."""
    out = {}
    th_f = TWO_PI * (tau1 - tauQ)          # fark kanadı
    th_s = -TWO_PI * (tau1 + tauQ)         # toplam kanadı
    n2 = abs(hQ) ** 2
    hb = np.conj(hQ)
    for k in range(0, kmax + 1):
        Gf = PT.Gam(th_f, k)
        Gs = PT.Gam(th_s, k)
        Ta = hQ / 2 * PT.mom[k]
        Tb = np.conj(hQ) * xQ ** 2 / 8 * complex(PT.Gam(np.array([-2 * TWO_PI * tauQ]), k)[0])
        Sf = complex(np.sum(h1 * np.conj(x1) * Gf))
        Ss = complex(np.sum(np.conj(h1) * x1 * Gs))
        Tc = xQ / 8 * (Sf + Ss)
        Tk = 0j
        if kule:
            Tk = xQ / 8 * complex(np.sum(kule[k])) if k in kule else 0j
        J = Ta + Tb + Tc + Tk
        out[f"u{k}_pred"] = 2 * float((hb * J).real) / n2
        out[f"s{k}_pred"] = 2 * float((hb * J).imag) / n2
        out[f"s{k}_Ta"] = 2 * float((hb * Ta).imag) / n2
        out[f"s{k}_Tb"] = 2 * float((hb * Tb).imag) / n2
        out[f"s{k}_fark"] = 2 * float((hb * xQ / 8 * Sf).imag) / n2
        out[f"s{k}_topl"] = 2 * float((hb * xQ / 8 * Ss).imag) / n2
        out[f"u{k}_Ta"] = 2 * float((hb * Ta).real) / n2
        out[f"u{k}_fark"] = 2 * float((hb * xQ / 8 * Sf).real) / n2
    return out
