# -*- coding: utf-8 -*-
"""
195o — ORTAK MODÜL: genelleştirilmiş (χ-ağırlıklı) tarak-uydusu makinesi
=======================================================================
KALEM_UYDU_KARAKTERI_24EYL2026 için. 188b/190b/193b makinesinin χ'ye
genelleştirilmiş hâli. χ ≡ 1, k = 1 iken zeta makinesine indirgenir (K0c'de
188b ile bit-bit sınanır).

Tanımlar (kalem + görev metni):
  L_χ(t) = log(k t / 2π);  2πN̄_χ(t) = t·L_χ(t) − t + 2πC_χ,
  C_χ = −χ(−1)/8  (çift −1/8, tek +1/8);  zeta: C = 7/8.
  a_q^χ = χ(q)·Λ(q)/(π √q log q)   (q asal-kuvvet; χ(q)=0 çizgisi YOK).
  seri (bir çizgi kümesi Q):  s(n) = Σ_{q'∈Q} 2 a_{q'}^χ sin(ω' g_n/2) cos(ω' m_n)
     (187b.katman_serisi / 188b.seri_ve_G dds yolu AYNEN — bit-bit).
  izdüşüm:  c_q = 2⟨s·e^{−iω_q m}⟩_blok   (190b.izdusum_blok deseni AYNEN).
  öz-terim: c_q^öz = 2⟨2a_q^χ sin(ω_q g/2) cos(ω_q m) e^{−iω_q m}⟩_blok
     (185b.c_oz ile aynı nicelik; K0c'de sınanır).
  blok b: pencere çizgileri τ = log q / L_b ∈ [0.45, 0.86), χ(q) ≠ 0;
          uydu n çizgileri |log q' − L_b − log n| < δ = 0.03, χ(q') ≠ 0.
  K̃_uydu = Σ_b N_b·num_b / Σ_b N_b·den_b,
     num_b = Σ_q c_{q,b}^{(uydu)} conj(c_{q,b}^öz),  den_b = Σ_q |c_{q,b}^öz|²
     (N_b = bloğun TAM nokta sayısı; blok-loo jackknife).
  Ĝ_b(Δω) = mean_{n∈b} e^{i(L_b+Δω) m_n}  (MUTLAK m_n).
"""
import ast
import hashlib
import os
import re
import sys
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
           "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S155, S188, S190, S193, S195 = (SCR / d for d in ["155", "188", "190", "193", "195"])
ZD = S190 / "zincir_dusuk"
TWO_PI = 2 * np.pi

ADALAR = ["beta", "chi3", "chi5e", "chi8e", "chi8o"]
ADA_ETIKET = {"beta": "β (χ₄)", "chi3": "χ₃", "chi5e": "χ₅ₑ", "chi8e": "χ₈ₑ",
              "chi8o": "χ₈ₒ"}
# 118a'daki tablo adları (üretici kod)
ADA_TABLO = {"beta": "CHI4", "chi3": "CHI3", "chi5e": "CHI5E", "chi8e": "CHI8E",
             "chi8o": "CHI8O"}

DELTA = 0.03          # uydu seçim yarı-genişliği (193 AYNEN)
DL = 0.05             # eşit-ΔL blok genişliği (adalar, zeta düşük ızgara)
L_UST = 8.5           # adalarda üst bölge alt sınırı
KISMI_MIN = 0.0125    # üstteki kısmi blok en az bu kadar L-genişliği taşımalı
TAU_LO, TAU_HI = 0.45, 0.86
BANT = (10.0, 10.4)   # L-eşli bant (H-195B3)

# uydular: ad -> (a, b)  (n = a/b)
UYDULAR = {"+log2": (2, 1), "-log2": (1, 2), "+log3": (3, 1), "-log3": (1, 3),
           "+log4": (4, 1), "+log5": (5, 1), "-log5": (1, 5), "+log6": (6, 1),
           "+log10": (10, 1), "+log(3/2)": (3, 2)}
UYDU_SIRA = list(UYDULAR.keys())
A_LISTE = ["+log2", "-log2", "+log3", "-log3", "+log4", "+log5", "-log5",
           "+log6", "+log10"]          # görev: n ∈ {±2, ±3, 4, ±5, 6, 10}
A_KAYIT = ["+log(3/2)"]


def log_n(ad):
    a, b = UYDULAR[ad]
    return float(np.log(a / b))


def asal_carpanlar(n):
    out, p = set(), 2
    while p * p <= n:
        while n % p == 0:
            out.add(p)
            n //= p
        p += 1
    if n > 1:
        out.add(n)
    return out


def yasak_mi(ad, k):
    """n = a/b'nin bir asal çarpanı k'yı bölüyorsa YASAK (kalem madde 2)."""
    a, b = UYDULAR[ad]
    return any(k % p == 0 for p in asal_carpanlar(a) | asal_carpanlar(b))


def jk_se(reps):
    reps = np.asarray(reps)
    B = reps.shape[0]
    return np.sqrt((B - 1) / B * np.sum((reps - reps.mean(0)) ** 2, 0))


def sha(yol):
    return hashlib.sha256(Path(yol).read_bytes()).hexdigest()


# ---------------------------------------------------------------- karakterler
def karakter_tablolari():
    """118a (Not-5 / 118b üreticisi) metninden χ tabloları + ADALAR satırları."""
    txt = (QM / "118a_taperli_L_motoru.py").read_text()
    tab = {}
    for ad, isim in ADA_TABLO.items():
        m = re.search(rf"^{isim} = (\{{.*?\}})\s*$", txt, re.M)
        tab[ad] = {int(k_): int(v) for k_, v in ast.literal_eval(m.group(1)).items()}
    satir = {}
    for m in re.finditer(r'\("(\w+)",\s*(\d+),\s*(\w+),\s*(\d+),', txt):
        satir[m.group(1)] = (int(m.group(2)), m.group(3), int(m.group(4)))
    return tab, satir


def chi_dizi(tablo, k, q):
    """χ(q) (tamamen çarpımsal; q asal-kuvvet de olsa χ(q mod k))."""
    lut = np.array([tablo[r] for r in range(k)], float)
    return lut[(np.asarray(q, np.int64) % k)]


# ---------------------------------------------------------------- kinematik
def Nbar(t, k, C):
    t = np.asarray(t, float)
    return (t * np.log(k * t / TWO_PI) - t) / TWO_PI + C


def C_chi(parite):
    return -1.0 / 8 if parite == 0 else 1.0 / 8


def ada_kin(ad):
    d = np.load(QM / f"118b_{ad}_zeros.npz")
    z = np.sort(np.asarray(d["zeros"], float))
    k, a = int(d["q"]), int(d["a"])
    g = np.diff(z)
    mid = 0.5 * (z[:-1] + z[1:])
    return dict(z=z, g=g, mid=mid, k=k, a=a, Lm=np.log(k * mid / TWO_PI))


def zeta_kin(pencere):
    """son: 155 eta_son önbelleği (185b 'gercek' kolu AYNEN); dusuk: 190 zinciri."""
    yol = S155 / "eta_son_t0.4_c4000.npz" if pencere == "son" else \
        ZD / "eta_dusuk_t0.4_c4000.npz"
    d = np.load(yol)
    mid = np.asarray(d["mid"], float)
    ds = np.asarray(d["ds"], float)
    g = (ds + 1.0) * TWO_PI / np.log(mid / TWO_PI)
    n0 = {"son": 1701053, "dusuk": 200001}[pencere]   # K0a'da sıfır dosyasıyla doğrulanır
    return dict(g=g, mid=mid, k=1, a=None, Lm=np.log(mid / TWO_PI), n0=n0,
                L=float(d["L"]))


# ---------------------------------------------------------------- bloklar
def esit_sayim_bloklari(N, nblok=8):
    """184-193 jk blokları AYNEN: kenar = linspace(0, N, 9)."""
    kenar = np.linspace(0, N, nblok + 1).astype(int)
    return [np.arange(kenar[b], kenar[b + 1]) for b in range(nblok)]


def dL_bloklari(Lm, alt_sinir, ust_sinir=None, dl=DL, kismi_min=KISMI_MIN):
    """Eşit ΔL ızgarası e_j = alt_sinir + dl·j; blok j = {e_j ≤ L < e_{j+1}}.
    Üstte kısmi blok ancak L-genişliği ≥ kismi_min ise alınır.
    Döner: [(j, e_lo, e_hi_etkin, idx)]"""
    Lmin = float(Lm.min())
    U = float(Lm.max()) if ust_sinir is None else min(float(Lm.max()), ust_sinir)
    out = []
    j = 0
    while True:
        lo = round(alt_sinir + dl * j, 10)
        hi = round(alt_sinir + dl * (j + 1), 10)
        if lo >= U:
            break
        idx = np.where((Lm >= lo) & (Lm < hi) & (Lm <= U))[0]
        span = min(hi, U) - max(lo, Lmin)
        if len(idx) > 0 and span >= kismi_min:
            out.append((j, lo, min(hi, U), idx))
        j += 1
    return out


# ---------------------------------------------------------------- çizgiler
def asal_kuvvetler(qmax):
    """187b.asal_kuvvetler AYNEN (sympy eleği); (q, Λ, a_q). a_q ifadesi 187b ile
    aynı işlem sırası: lams/(π·√q·log q)."""
    from sympy import primerange
    qs, lams = [], []
    for p in primerange(2, int(qmax) + 1):
        qs.append(p)
        lams.append(np.log(p))
        pk = p * p
        while pk <= qmax:
            qs.append(pk)
            lams.append(np.log(p))
            pk *= p
    qs = np.array(qs, float)
    lams = np.array(lams)
    srt = np.argsort(qs)
    qs, lams = qs[srt], lams[srt]
    lq = np.log(qs)
    return qs, lams, lams / (np.pi * np.sqrt(qs) * lq)


# ---------------------------------------------------------------- seri / izdüşüm
PBLOK = 2048      # 187b/188b AYNEN
LBLOK = 4096      # 187b/188b AYNEN
PROJ_BLOK = 1500  # 186b/187b/190b AYNEN


def seri(g, mid, w_l, a_l):
    """dds(n) = Σ 2a_l sin(w g/2) cos(w m) — 188b.seri_ve_G'nin dds yolu AYNEN
    (aynı bloklama ve işlem sırası). a_l = χ·a (χ ≡ 1 ⇒ a_l ≡ a, bit-bit)."""
    N = len(mid)
    dds = np.zeros(N)
    kats = 2.0 * a_l
    nl = len(w_l)
    for s0 in range(0, N, PBLOK):
        sl = slice(s0, min(s0 + PBLOK, N))
        acc = np.zeros(sl.stop - sl.start)
        for l0 in range(0, nl, LBLOK):
            ll = slice(l0, min(l0 + LBLOK, nl))
            P = np.outer(mid[sl], w_l[ll])
            np.cos(P, out=P)
            K = np.sin(0.5 * np.outer(g[sl], w_l[ll]))
            K *= P
            acc += K @ kats[ll]
            del P, K
        dds[sl] = acc
    return dds


def izdusum_blok(D, mid_b, w_win):
    """190b.izdusum_blok AYNEN: re = Σ D cos(wm), im = −Σ D sin(wm)."""
    ns = D.shape[1]
    re_ = np.zeros((ns, len(w_win)))
    im_ = np.zeros((ns, len(w_win)))
    n = len(mid_b)
    for s0 in range(0, n, PROJ_BLOK):
        sl = slice(s0, min(s0 + PROJ_BLOK, n))
        P = np.outer(mid_b[sl], w_win)
        cP = np.cos(P)
        sP = np.sin(P)
        del P
        re_ += D[sl].T @ cP
        im_ -= D[sl].T @ sP
        del cP, sP
    return re_, im_


def c_oz_blok(g_b, mid_b, w, a_chi):
    """c^öz = (4a/N) Σ sin(wg/2) cos(wm) e^{−iwm} (185b.c_oz ile aynı nicelik)."""
    N = len(mid_b)
    Sre = np.zeros(len(w))
    Sim = np.zeros(len(w))
    for s0 in range(0, N, PROJ_BLOK):
        sl = slice(s0, min(s0 + PROJ_BLOK, N))
        P = np.outer(mid_b[sl], w)
        c = np.cos(P)
        s = np.sin(P)
        del P
        X = np.sin(0.5 * np.outer(g_b[sl], w))
        X *= c
        Sre += (X * c).sum(0)
        Sim -= (X * s).sum(0)
        del c, s, X
    return 4.0 * a_chi * (Sre + 1j * Sim) / N


# ---------------------------------------------------------------- çekirdek (blok)
def cekirdek_blok(g_full, mid_full, idx_full, idx_an, Lb, Q, A, CHI, uydular):
    """Tek blok: pencere (τ=log q/L_b∈[0.45,0.86), χ≠0) öz-terim + her uydu için
    blok-yerel seçim, seri (idx_an noktalarında), izdüşüm.
    Döner dict: num (nU complex), den, N_tam, N_an, ncz_pencere, ncz_uydu (nU),
    c_oz (nwin complex), c_uydu (nU, nwin complex), q_pencere."""
    lq = np.log(Q)
    tau = lq / Lb
    pw = (tau >= TAU_LO) & (tau < TAU_HI) & (CHI != 0)
    w_win = lq[pw]
    a_win = CHI[pw] * A[pw]
    g_b, m_b = g_full[idx_full], mid_full[idx_full]
    coz = c_oz_blok(g_b, m_b, w_win, a_win)
    g_a, m_a = g_full[idx_an], mid_full[idx_an]
    cols, ncz = [], []
    for ad in uydular:
        sel = (np.abs(lq - Lb - log_n(ad)) < DELTA) & (CHI != 0)
        ncz.append(int(sel.sum()))
        if sel.any():
            cols.append(seri(g_a, m_a, lq[sel], CHI[sel] * A[sel]))
        else:
            cols.append(np.zeros(len(m_a)))
    D = np.stack(cols, axis=1)
    re_, im_ = izdusum_blok(D, m_a, w_win)
    cu = 2.0 * (re_ + 1j * im_) / len(m_a)
    num = cu @ np.conj(coz)
    den = float(np.sum(np.abs(coz) ** 2))
    return dict(num=num, den=den, N_tam=len(idx_full), N_an=len(idx_an),
                ncz_pencere=int(pw.sum()), ncz_uydu=np.array(ncz), c_oz=coz,
                c_uydu=cu, q_pencere=Q[pw], Lb=Lb)


def havuz(num, den, Nw):
    """num (B, nU), den (B,), Nw (B,) → K̃ (nU), loo replikaları (B, nU)."""
    num = np.asarray(num)
    den = np.asarray(den)
    Nw = np.asarray(Nw, float)
    P = (Nw[:, None] * num).sum(0)
    Q = (Nw * den).sum()
    K = P / Q
    reps = np.array([(P - Nw[j] * num[j]) / (Q - Nw[j] * den[j])
                     for j in range(len(den))])
    return K, reps


def ozet(K, reps):
    """Re/Im K̃, jk se, z; açı."""
    se_re = jk_se(reps.real)
    se_im = jk_se(reps.imag)
    return dict(re=K.real, im=K.imag, se_re=se_re, se_im=se_im,
                z_re=K.real / se_re, z_im=K.imag / se_im,
                aci=np.degrees(np.angle(K)))


# ---------------------------------------------------------------- Ĝ gücü
DW_IZGARA = np.round(np.arange(-2.0, 2.6 + 1e-9, 0.001), 6)
G_PBLOK = 1024


def G_guc_blok(mid_b, Lb, dgrid=DW_IZGARA):
    """|Ĝ_b(Δω)|², Ĝ_b = mean e^{i(L_b+Δω)m} (MUTLAK m)."""
    w = Lb + dgrid
    re_ = np.zeros(len(w))
    im_ = np.zeros(len(w))
    for s0 in range(0, len(mid_b), G_PBLOK):
        P = np.outer(mid_b[s0:s0 + G_PBLOK], w)
        re_ += np.cos(P).sum(0)
        im_ += np.sin(P).sum(0)
        del P
    n = len(mid_b)
    return (re_ / n) ** 2 + (im_ / n) ** 2


def rasyonel_katalog(hmax=6):
    """1 ≤ a, b ≤ hmax sade rasyonellerin log'ları (Δω=0 dahil)."""
    from math import gcd
    out = set()
    for a in range(1, hmax + 1):
        for b in range(1, hmax + 1):
            if gcd(a, b) == 1:
                out.add((a, b))
    return sorted(out, key=lambda x: np.log(x[0] / x[1]))


PENCERE_YARI = 0.03     # uydu penceresi |Δω − log n| < 0.03 (görev)
HALKA_IC, HALKA_DIS = 0.06, 0.25
DISLAMA = 0.035         # halkadan dışlanan komşu rasyonel yarı-genişliği
KATALOG_H = 6


def pencere_halka_maskeleri(ad, dgrid=DW_IZGARA):
    ln = log_n(ad)
    W = np.abs(dgrid - ln) < PENCERE_YARI
    d = np.abs(dgrid - ln)
    R = (d >= HALKA_IC) & (d <= HALKA_DIS)
    a0, b0 = UYDULAR[ad]
    from math import gcd
    gg = gcd(a0, b0)
    kendi = (a0 // gg, b0 // gg)
    for (a, b) in rasyonel_katalog(KATALOG_H):
        if (a, b) == kendi:
            continue
        R &= np.abs(dgrid - np.log(a / b)) >= DISLAMA
    return W, R
