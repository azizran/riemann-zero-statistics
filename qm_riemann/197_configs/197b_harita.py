# -*- coding: utf-8 -*-
"""
197b — AYNA YASASI: 32-BLOK Δω HARİTASI (190b/188b makinesi AYNEN + KALEM'in beyanlı sapmaları)
===========================================================================================
KALEM_AYNA_YASASI_25EYL2026.md "Ölçüm tanımı" (S1, S1b, S2, S3) ve "Kapılar" (M1, M3).
190b_harita.py importlib ile yüklenir (o da 188b'yi yükler); DOSYALAR DÜZENLENMEZ.
AYNEN çağrılanlar: b190.onkayit, b190.kinematik_eta, b190.izdusum_blok,
b190._om_init/_om_is (ω-görev işçisi), b188.jk_kenar, b188.dilimleri_kur, b188.karisim,
b188.bant_maskeleri, b188.K_matris, b188.c_proj, b187.asal_kuvvetler, b187.c_kesik.
Çalışma-anı yamaları: b188.S188 → dilim dizini (190b emsali); b190.S190 → scratchpad/190.

S1b K_düz (BİRİNCİL): ω-dilimi (b, j) için mix_j = mix − m_{S_bj}; m_S = S dilimindeki
çizgilerin (yalnız c^kesik evrenindekiler, K1: q ≤ e^{0.86L}) kesik seriye katkısının
HAVUZ çizgilerine TAM-örneklem izdüşümü (186b merdiven terimi 2a sin(ωg/2)cos(ωm), çizgi
başına; 188b.c_proj deseniyle 8-grup loo). K_düz = K_matris(C_j, mix_j). K_ham (çıkarmasız)
KAYIT. Katkı makinesi Σ_{tüm K1} m ≡ c^kesik ile denetlenir (M1c / M-kesik).

YAPILANDIRMA ANAHTARI (ilk argüman):
  m1         Sapmalar KAPALI: tam HAVUZ τ∈[0.45,0.86), 8 blok, yalnız τ'∈(0.86,1.30];
             50 K satırı. M1a: K_ham ≡ 190 K; M1b: K_düz (aynı dilim-başına çıkarma yolu)
             ≡ 190 K; M1c (yalnız tam): Σ_{K1} katkı ≡ c^kesik (tam + 8 loo).
             KÖRLÜK: alt ucu Δω < −1.3'e değen dilim (j < −51) HESAPLANMAZ/OKUNMAZ.
             --alt: temsilî 12 ω-dilimi (merkez −1.275 … 3.5), tüm 50 K satırı (M1c yok).
  olcum      S1 HAVUZ' τ∈[0.45,0.74); S1b K_düz; S2 yeni τ'-dilimleri (0.74,0.86] (0.005,
             188b.dilimleri_kur AYNEN) + 190'ın (0.86,1.30]; S3 32 blok, jk 8 grup × 4 blok;
             Δω ∈ [−2.2, +2.4] (tam dilimler; 32 bloğun hepsi kapsar — assert).
             ONKAYIT_197 + bu kod sürümüyle M1 geçişi ŞART. EKRANA/LOG'A κ/K/P YAZILMAZ.
  zamanlama  yalnız dilim serisi üretimi + katkı makinesi zamanlaması; κ HESAPLANMAZ.
Kullanım: 197b_harita.py <m1|olcum|zamanlama> [--alt] [--isci N]
Kütüphane (197c): jk, profil_sigma, parlaklik, cizgi_bicimi.
"""
import argparse
import hashlib
import importlib.util
import json
import os
import re
import sys
import time
from pathlib import Path
from types import SimpleNamespace

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
           "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = QM / "scratchpad"
S188, S190, S197 = SCR / "188", SCR / "190", SCR / "197"
ZD = S190 / "zincir_dusuk"
KALEM = QM / "KALEM_AYNA_YASASI_25EYL2026.md"
TWO_PI = 2 * np.pi
NJACK = 8
DW = 0.025            # ω-dilim genişliği (190 AYNEN)
EPS = 1e-9            # 192 kapsama/pencere toleransı AYNEN
OM_CHUNK = 8000       # 190b AYNEN
KK_BLOK = 1500        # katkı makinesi nokta bloğu

# ---- KÖRLÜK: Δω < −1.3 ön-kayıtlı kör sınav bandı (m1/zamanlama'da dokunulmaz) ----
J_GOR = -51           # en düşük görülmüş dilim: alt ucu (−51.5)·0.025 = −1.2875 ≥ −1.3
ALT_J = np.round(np.linspace(J_GOR, 140, 12)).astype(int)   # M1 temsilî: −1.275 … 3.5

HAVUZ_190 = (0.45, 0.86)          # m1 (190 AYNEN)
HAVUZ_197 = (0.45, 0.74)          # S1: HAVUZ' = τ ∈ [0.45, 0.74)
TAU_YENI = (0.74, 0.86)           # S2: yeni τ'-dilimleri (e, e+0.005], 24 dilim
DTAU = 0.005
NBLOK = 32                        # S3
DW_ARALIK = (-2.2, 2.4)           # olcum: tam dilim merkezleri (KALEM ölçüm notu)
T190_OMEGA_YEDEK = 212.0          # log_190b ω aşaması (6 işçi) — süre kalibrasyonu


def yukle(ad, yol):
    spec = importlib.util.spec_from_file_location(ad, yol)
    m = importlib.util.module_from_spec(spec)
    sys.modules[ad] = m
    spec.loader.exec_module(m)
    return m


b190 = yukle("b190", QM / "190_configs" / "190b_harita.py")
b190.S190 = S190
b188 = b190.b188
b187 = b188.b187


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


# ============================ ortak kurulum ============================
def ortak():
    ONK = b190.onkayit()                       # 190a sha denetimi (190b AYNEN)
    ONK188 = json.load(open(S188 / "ONKAYIT_188.json"))
    kenar_eski = np.array(ONK["dilim_izgara"])
    assert np.array_equal(kenar_eski, np.array(ONK188["dilim_izgara_gercek"]))
    assert ONK["ince_bant_kenar"] == ONK188["ince_bant_kenar"]
    assert ONK["bant8_kenar"] == ONK188["bant8_kenar"]
    L = float(ONK["pencere"]["L"])
    G = np.load(ZD / "K1_gercek_dusuk.npz")
    OZ = np.load(ZD / "OZ_gercek_dusuk.npz")
    P = np.load(ZD / "G1_proj_gercek_dusuk.npz")
    g, mid, Lg = b190.kinematik_eta(ZD / "eta_dusuk_t0.4_c4000.npz")
    assert Lg == L, (Lg, L)
    N = len(mid)
    kj8 = b188.jk_kenar(N)
    if not (np.array_equal(np.diff(kj8), G["nb"]) and
            np.array_equal(np.diff(kj8), P["nb"])):
        raise SystemExit("JK BLOK KENARLARI 184 İLE UYUMSUZ")
    bid8 = np.searchsorted(kj8, np.arange(N), side="right") - 1
    t0 = time.time()
    qa, _, taua, aqa = b187.asal_kuvvetler(np.exp(kenar_eski[-1] * L), L)
    print(f"[dusuk] N={N} L={L:.9f}; asal-kuvvetler q≤e^(1.30L): {len(qa)} "
          f"({time.time()-t0:.0f}s)", flush=True)
    q_all = np.asarray(G["q"], float)
    # tarama çizgisi (qa indisi) → K1 satırı (c^kesik evreni), yoksa −1
    k1_satir = np.full(len(qa), -1, int)
    ix = np.searchsorted(qa, q_all)
    assert np.array_equal(qa[ix], q_all), "K1 çizgileri asal-kuvvet listesinde değil"
    k1_satir[ix] = np.arange(len(q_all))
    return SimpleNamespace(
        ONK=ONK, ONK188=ONK188, kenar_eski=kenar_eski, L=L, OZ=OZ, P=P,
        w_all=np.asarray(G["w"], float), tau_all=np.asarray(G["tau"], float),
        aq_all=np.asarray(G["aq"], float), q_all=q_all, k1_satir=k1_satir,
        g=g, mid=mid, N=N, kj8=kj8, bid8=bid8, qa=qa, taua=taua, aqa=aqa)


def bloklar(N, mid, nblok):
    """S3: kenar = linspace(0, N, nblok+1).astype(int); L_b = blok içi mean log(m/2π)."""
    kb = np.linspace(0, N, nblok + 1).astype(int)
    Lb = np.array([float(np.log(mid[kb[b]:kb[b + 1]] / TWO_PI).mean())
                   for b in range(nblok)])
    return kb, Lb


def havuz(c, aralik):
    win = (c.tau_all >= aralik[0]) & (c.tau_all < aralik[1])
    return win, c.w_all[win], c.tau_all[win], c.q_all[win]


def kenar_yeni_kur():
    n = int(round((TAU_YENI[1] - TAU_YENI[0]) / DTAU))
    k = np.round(TAU_YENI[0] + DTAU * np.arange(n + 1), 3)
    assert k[0] == TAU_YENI[0] and k[-1] == TAU_YENI[1]
    return k


def m3_oz_terim(q_h, q_t, etiket):
    """M3 (KALEM): tarama çizgileri ∩ HAVUZ = ∅ (tam sayı q karşılaştırması)."""
    a = np.round(q_h).astype(np.int64)
    b = np.unique(np.round(q_t).astype(np.int64))
    kes = np.intersect1d(a, b)
    assert kes.size == 0, (f"M3 İHLAL [{etiket}]: HAVUZ ∩ tarama ≠ ∅ "
                           f"({kes.size} çizgi; ilk {kes[:5].tolist()})")
    print(f"  M3 [{etiket}] HAVUZ ({a.size} çizgi, q∈[{a.min()},{a.max()}]) ∩ tarama "
          f"({b.size} çizgi, q∈[{b.min()},{b.max()}]) = ∅  ✓", flush=True)
    return {"n_havuz": int(a.size), "n_tarama": int(b.size), "kesisim": 0,
            "q_havuz_maks": int(a.max()), "q_tarama_min": int(b.min())}


# ============================ S1b katkı makinesi ============================
_KK = {}


def _kk_init(g, mid, w_l, a_l, w_h):
    _KK.update(g=g, mid=mid, w_l=w_l, a_l=a_l, w_h=w_h)


def _kk_is(arg):
    gi, lo, hi = arg
    t0 = time.time()
    g, mid, w_l, a_l, w_h = (_KK[k] for k in ("g", "mid", "w_l", "a_l", "w_h"))
    kats = 2.0 * a_l
    re_ = np.zeros((len(w_l), len(w_h)))
    im_ = np.zeros((len(w_l), len(w_h)))
    for s0 in range(lo, hi, KK_BLOK):
        sl = slice(s0, min(s0 + KK_BLOK, hi))
        D = np.sin(0.5 * np.outer(g[sl], w_l))
        D *= np.cos(np.outer(mid[sl], w_l))
        D *= kats                                   # çizgi başına kesik-seri terimi
        Ph = np.outer(mid[sl], w_h)
        re_ += D.T @ np.cos(Ph)
        im_ -= D.T @ np.sin(Ph)
    return gi, re_, im_, time.time() - t0


def kesik_katki(c, w_h, nw):
    """K1 evrenindeki HER çizgi q'' için kesik seri terimi 2a''·sin(ω''g/2)·cos(ω''m)'nin
    HAVUZ çizgilerine izdüşümü, 8 jk grubu toplamı (186b işaret deseni) → re, im (8, n_K1, nh)."""
    import multiprocessing as mp
    t0 = time.time()
    nl, nh = len(c.w_all), len(w_h)
    re_ = np.zeros((NJACK, nl, nh))
    im_ = np.zeros((NJACK, nl, nh))
    isler = [(gi, int(c.kj8[gi]), int(c.kj8[gi + 1])) for gi in range(NJACK)]
    with mp.get_context("fork").Pool(min(nw, NJACK), initializer=_kk_init,
                                     initargs=(c.g, c.mid, c.w_all, c.aq_all, w_h)) as pool:
        for gi, r_, i_, dt in pool.imap_unordered(_kk_is, isler):
            re_[gi], im_[gi] = r_, i_
    print(f"  [katkı] {nl} K1 çizgisi → {nh} HAVUZ çizgisi, 8 grup ({time.time()-t0:.0f}s)",
          flush=True)
    return re_, im_


def kesik_denetim(c, re_, im_, win):
    """Σ_{tüm K1} katkı ≡ c^kesik (187b.c_kesik, tam + 8 loo) — HAVUZ çizgilerinde."""
    nb8 = np.diff(c.kj8)
    out = []
    for v in range(NJACK + 1):
        top = b188.c_proj(re_.sum(1), im_.sum(1), nb8, c.N, v - 1)
        ref = b187.c_kesik(c.P, v - 1)[win]
        out.append(float(np.max(np.abs(top - ref)) / np.max(np.abs(ref))))
    return out


def katki_kur(k1_satir, Mv):
    """(v, qa indisleri) → dilimin c^kesik evrenindeki çizgilerinin toplam katkısı; yoksa None."""
    def katki(v, satir_qa):
        r = k1_satir[satir_qa]
        r = r[r >= 0]
        if r.size == 0:
            return None
        return Mv[v][r].sum(0)
    return katki


# ============================ ω-dilim makinesi (190b D AYNEN) ============================
def plan_kur(qa, aqa, sec, Lb, dizin, jfiltre=None):
    """190b D planı AYNEN: blok b, j = floor((log q' − L_b)/DW + 0.5), OM_CHUNK gruplama."""
    isler, plan, say = [], {}, {}
    for b in range(len(Lb)):
        io = sec[b]
        j = np.floor((np.log(qa[io]) - Lb[b]) / DW + 0.5).astype(int)
        uj = np.unique(j)
        if jfiltre is not None:
            uj = jfiltre(b, uj)
        plan[b] = uj
        say[b] = {int(jj): int((j == jj).sum()) for jj in uj}
        grup, sayac, k = [], 0, 0
        for jj in uj:
            s = io[j == jj]
            grup.append((int(jj), np.log(qa[s]), aqa[s]))
            sayac += len(s)
            if sayac >= OM_CHUNK or jj == uj[-1]:
                yol = dizin / f"b{b}_c{k}.npz"
                if not yol.exists():
                    isler.append((b, k, yol, [x[0] for x in grup],
                                  [x[1] for x in grup], [x[2] for x in grup]))
                grup, sayac = [], 0
                k += 1
    return isler, plan, say


def omega_kos(isler, c, w_win, kb, nw, etiket):
    """190b D işçi havuzu AYNEN (b190._om_init/_om_is); kenar = blok kenarları."""
    if not isler:
        return 0.0
    import multiprocessing as mp
    isler.sort(key=lambda x: -sum(len(w_) for w_ in x[4]))
    t0 = time.time()
    with mp.get_context("fork").Pool(nw, initializer=b190._om_init,
                                     initargs=(c.g, c.mid, c.bid8, w_win, kb)) as pool:
        for kk, (b, k, nl, dt) in enumerate(pool.imap_unordered(b190._om_is, isler), 1):
            print(f"    [{etiket} ω] blok {b} parça {k} bitti ({nl} çizgi, {dt:.0f}s)  "
                  f"{kk}/{len(isler)}  toplam {time.time()-t0:.0f}s", flush=True)
    return time.time() - t0


def omega_birlestir(dizin, plan, say, nblok, mixes, maskeler, sec, Lb, c, katki):
    """190b birleştirme AYNEN (K_ham) + S1b K_düz: C = 2(re+i·im)/n_b;
    K_ham = K_matris(C, mix); K_düz(j) = K_matris(C_j, mix − m_{S_bj}) — katkısı olmayan
    dilimler aynı çağrıda (mix − 0) toplu hesaplanır (190 düzeninde K_düz ≡ K bit-bit).
    → J, K_ham, K_düz (nv, nblok, nmaske, nJ), ncz, n_kesik, nb."""
    dolu = [b for b in range(nblok) if len(plan[b])]
    jmin = min(int(plan[b].min()) for b in dolu)
    jmax = max(int(plan[b].max()) for b in dolu)
    J = np.arange(jmin, jmax + 1)
    nJ = len(J)
    nv, nm = len(mixes), len(maskeler)
    Kh = np.zeros((nv, nblok, nm, nJ), complex)
    Kd = np.zeros((nv, nblok, nm, nJ), complex)
    n = np.zeros((nblok, nJ), int)
    nk = np.zeros((nblok, nJ), int)
    nb = np.zeros(nblok, int)
    sifir = np.zeros(len(mixes[0]), complex)
    for b in dolu:
        io = sec[b]
        jl = np.floor((np.log(c.qa[io]) - Lb[b]) / DW + 0.5).astype(int)
        parca = sorted(dizin.glob(f"b{b}_c*.npz"), key=lambda p: int(p.stem.split("_c")[1]))
        gordu = []
        for p in parca:
            d = np.load(p)
            C = 2.0 * (d["re"] + 1j * d["im"]) / int(d["nb"])
            jlist = d["j"]
            ix = jlist - jmin
            satir = [io[jl == jj] for jj in jlist]
            nk[b, ix] = [int((c.k1_satir[s] >= 0).sum()) for s in satir]
            for v, mx in enumerate(mixes):
                Kh[v, b][:, ix], _ = b188.K_matris(C, mx, maskeler)
                kat = [katki(v, s) for s in satir]
                bos = np.array([k for k, a in enumerate(kat) if a is None], int)
                if bos.size:
                    Kd[v, b][:, ix[bos]], _ = b188.K_matris(C[bos], mx - sifir, maskeler)
                for k, a in enumerate(kat):
                    if a is not None:
                        Kd[v, b][:, ix[k]] = b188.K_matris(C[k:k + 1], mx - a, maskeler)[0][:, 0]
            n[b, ix] = d["ncz"]
            nb[b] = int(d["nb"])
            gordu += [int(x) for x in jlist]
        assert sorted(gordu) == [int(x) for x in plan[b]], f"blok {b}: parça j'leri ≠ plan"
        assert all(n[b, jj - jmin] == say[b][int(jj)] for jj in plan[b]), f"blok {b}: ncz ≠ plan"
    return J, Kh, Kd, n, nk, nb


# ============================ profil / çizgi biçimi (KALEM; 197c için) ============================
def jk(r):
    r = np.asarray(r, float)
    return float(np.sqrt((NJACK - 1) / NJACK * np.sum((r - r.mean()) ** 2)))


def profil_sigma(K, nb, kap):
    """KALEM: κ_Σ(j) = Σ_{b∈B_j} κ_b(j)·n_b / Σ_{b∈B_j} n_b, κ_b = −Re K_b (192b AYNEN).
    K: (1+NJACK, nblok, nJ) — v=0 tam, v=1+i mix loo-i; replika i'de grup i'nin
    (nblok/8 ardışık blok) blokları B_j'den çıkar. → (1+NJACK, nJ)."""
    nblok = K.shape[1]
    grup = np.arange(nblok) // (nblok // NJACK)
    kb_ = -K.real
    w = nb[:, None] * kap
    out = np.full((NJACK + 1, K.shape[2]), np.nan)
    with np.errstate(invalid="ignore", divide="ignore"):
        out[0] = (kb_[0] * w).sum(0) / w.sum(0)
        for i in range(NJACK):
            wi = w * (grup != i)[:, None]
            out[1 + i] = (kb_[1 + i] * wi).sum(0) / wi.sum(0)
    return out


def parlaklik(prof, J, h, yari=0.03):
    """KALEM M2: P(h) = Σ_{|c−h| ≤ 0.03} κ_Σ(c) (+EPS, 192 AYNEN); se = jk(8 replika)."""
    m = np.abs(J * DW - h) <= yari + EPS
    p = float(prof[0, m].sum())
    pr = prof[1:, m].sum(1)
    return p, pr, jk(pr), J[m]


def cizgi_bicimi(mid, kb, Lb, kap, J, log_r):
    """KALEM ÇİZGİ BİÇİMİ (parametresiz, kinematik):
    h_b(j; r) = (1/n_b)·#{n ∈ b : floor((L_n − L_b + log r)/0.025 + 0.5) = j}, L_n = log(m_n/2π);
    h(j; r) = Σ_{b∈B_j} n_b·h_b / Σ_{b∈B_j} n_b; replika i: grup i blokları B_j'den çıkar
    (profil_sigma ile tutarlı). → (1+NJACK, nJ)."""
    nblok = len(Lb)
    nb = np.diff(kb)
    grup = np.arange(nblok) // (nblok // NJACK)
    Ln = np.log(mid / TWO_PI)
    hb = np.zeros((nblok, len(J)))
    for b in range(nblok):
        j = np.floor((Ln[kb[b]:kb[b + 1]] - Lb[b] + log_r) / DW + 0.5).astype(int)
        ok = (j >= J[0]) & (j <= J[-1])
        hb[b] = np.bincount(j[ok] - J[0], minlength=len(J)) / nb[b]
    w = nb[:, None] * kap
    out = np.full((NJACK + 1, len(J)), np.nan)
    with np.errstate(invalid="ignore", divide="ignore"):
        out[0] = (hb * w).sum(0) / w.sum(0)
        for i in range(NJACK):
            wi = w * (grup != i)[:, None]
            out[1 + i] = (hb * wi).sum(0) / wi.sum(0)
    return out


# ============================ m1 ============================
def mod_m1(c, alt, nw):
    t00 = time.time()
    kb, Lb = bloklar(c.N, c.mid, NJACK)
    assert np.array_equal(kb, c.kj8)
    assert np.array_equal(Lb, np.array(c.ONK["bloklar"]["L_b"])), "L_b ÖN-KAYITLA UYUMSUZ"
    win, w_win, tau_win, q_win = havuz(c, HAVUZ_190)
    kenar = c.kenar_eski
    io = np.where((c.taua > kenar[0]) & (c.taua <= kenar[-1]))[0]
    m3 = m3_oz_terim(q_win, c.qa[io], "m1")
    maskeler, et = b188.bant_maskeleri(tau_win, c.ONK188)
    mix_tam = b188.karisim(c.P, c.OZ, c.aq_all)[win]
    m1c, Mv = None, None
    if not alt:                                  # M1c: katkı makinesi ≡ c^kesik
        re_k, im_k = kesik_katki(c, w_win, nw)
        m1c = kesik_denetim(c, re_k, im_k, win)
        print(f"  M1c Σ_K1 katkı ≡ c^kesik (1044 HAVUZ çizgisi): maks bağıl (tam+8 loo) = "
              f"{max(m1c):.2e}", flush=True)
        Mv = [b188.c_proj(re_k, im_k, np.diff(c.kj8), c.N, -1)]
        del re_k, im_k
    katki = katki_kur(c.k1_satir, Mv)
    if alt:
        jf, ad = (lambda b, uj: uj[np.isin(uj, ALT_J)]), "omega_alt"
    else:
        jf, ad = (lambda b, uj: uj[uj >= J_GOR]), "omega"
    dizin = S197 / "m1" / ad
    dizin.mkdir(parents=True, exist_ok=True)
    sec = {b: io for b in range(NJACK)}
    isler, plan, say = plan_kur(c.qa, c.aqa, sec, Lb, dizin, jf)
    for b in range(NJACK):
        assert len(plan[b]) == 0 or plan[b].min() >= J_GOR     # KÖRLÜK
    print(f"  [m1{' alt' if alt else ''}] 8 blok; hesaplanacak dilim sayıları "
          f"{[len(plan[b]) for b in range(NJACK)]} (hepsi j ≥ {J_GOR}); eksik görev "
          f"{len(isler)} (işçi {nw})", flush=True)
    t_om = omega_kos(isler, c, w_win, kb, nw, "m1")
    J, Kh, Kd, n, nk, nb = omega_birlestir(dizin, plan, say, NJACK, [mix_tam], maskeler,
                                           sec, Lb, c, katki)
    Kh, Kd = Kh[0], Kd[0]                           # (8, 50, nJ)
    assert nk.sum() == 0, "m1'de tarama çizgisi c^kesik evreninde olmamalı"

    # ---- referans: YALNIZ görülmüş dilimler (j ≥ J_GOR); kör hücreler okunmaz ----
    H = np.load(S190 / "harita_omega_dusuk.npz")
    Jr = H["J"]
    gor = Jr >= J_GOR
    Jg = Jr[gor]
    Kref = H["K"][:, :, gor]
    nref = H["ncz"][:, gor]
    et_ref = [str(x) for x in H["etiket"]]
    nb_ref, Lb_ref = H["nb"], H["L_b"]
    del H
    assert et_ref == list(et), "etiketler uyumsuz"
    nb_esit = bool(np.array_equal(nb_ref, nb))
    Lb_esit = bool(np.array_equal(Lb_ref, Lb))
    rap = {}
    kume_tam, n_esit, nhucre = True, True, 0
    for ad_k, KK in (("ham", Kh), ("duz", Kd)):
        el_maks, abs_maks, ref_maks, bitbit, sifir_esit = 0.0, 0.0, 0.0, True, True
        for b in range(NJACK):
            if ad_k == "ham" and not alt:
                kume_tam &= set(int(x) for x in plan[b]) == set(int(x) for x in Jg[nref[b] > 0])
            for jj in plan[b]:
                a = KK[b, :, jj - J[0]]
                r = Kref[b, :, jj - Jg[0]]
                if ad_k == "ham":
                    n_esit &= bool(n[b, jj - J[0]] == nref[b, jj - Jg[0]])
                    nhucre += 1
                bitbit &= bool(np.array_equal(a, r))
                d = np.abs(a - r)
                m = np.abs(r) > 0
                sifir_esit &= bool(np.all(a[~m] == 0))
                if m.any():
                    el_maks = max(el_maks, float(np.max(d[m] / np.abs(r[m]))))
                abs_maks = max(abs_maks, float(d.max()))
                ref_maks = max(ref_maks, float(np.abs(r).max()))
        rap[ad_k] = {"maks_bagil_eleman": el_maks,
                     "maks_mutlak_bolu_maks_ref": abs_maks / ref_maks if ref_maks else None,
                     "bit_bit": bitbit, "sifir_hucreler_esit": sifir_esit}
    gecti = bool(all(rap[k]["maks_bagil_eleman"] <= 1e-10 and rap[k]["sifir_hucreler_esit"]
                     for k in rap) and n_esit and nb_esit and Lb_esit and kume_tam
                 and nhucre > 0 and (m1c is None or max(m1c) <= 1e-10))
    sonuc = {
        "kapi": "M1 (KALEM): S1-S3 kapalı (tam HAVUZ [0.45,0.86), 8 blok, yalnız "
                "τ'∈(0.86,1.30]) → 190 harita_omega_dusuk.npz K (50 satır) ≤ 1e-10 bağıl; "
                "M1a K_ham, M1b K_düz (S1b dilim-başına çıkarma yolu; 190'da ≡ K), "
                "M1c Σ_K1 katkı ≡ c^kesik",
        "mod": "alt (12 temsilî ω-dilimi)" if alt else "tam (görülmüş bölge: tüm j ≥ −51)",
        "korluk": f"j < {J_GOR} (dilim alt ucu < −1.2875) HESAPLANMADI/OKUNMADI",
        "alt_j": ALT_J.tolist() if alt else None,
        "n_hucre_blok_x_dilim": nhucre, "n_K_satir": len(et),
        "n_karsilastirma": nhucre * len(et),
        "M1a_K_ham": rap["ham"], "M1b_K_duz": rap["duz"],
        "M1c_kesik_katki_maks_bagil": m1c,
        "maks_bagil_eleman": max(rap["ham"]["maks_bagil_eleman"],
                                 rap["duz"]["maks_bagil_eleman"]),
        "ncz_esit": bool(n_esit), "nb_esit": nb_esit, "L_b_esit": Lb_esit,
        "hucre_kumesi_tam": bool(kume_tam) if not alt else "uygulanmaz (alt)",
        "M3": m3, "esik": 1e-10, "gecti": gecti,
        "kod_sha256": sha(Path(__file__).resolve()),
        "sure_omega_s": t_om, "sure_toplam_s": time.time() - t00, "isci": nw,
        "zaman": time.strftime("%Y-%m-%d %H:%M:%S %z")}
    yol = S197 / "m1" / ("M1_sonuc_alt.json" if alt else "M1_sonuc.json")
    json.dump(sonuc, open(yol, "w"), indent=1, ensure_ascii=False)
    for k in ("ham", "duz"):
        print(f"  M1{'a' if k == 'ham' else 'b'} K_{k}: maks bağıl (eleman) = "
              f"{rap[k]['maks_bagil_eleman']:.2e}; bit-bit {rap[k]['bit_bit']}", flush=True)
    print(f"  M1 {'alt' if alt else 'tam'}: {nhucre} hücre × {len(et)} satır; ncz {n_esit}; "
          f"nb {nb_esit}; L_b {Lb_esit}; hücre kümesi {kume_tam}  → "
          f"{'GEÇTİ' if gecti else 'KALDI'}", flush=True)
    print(f"-> {yol}  ({time.time()-t00:.0f}s)", flush=True)
    return sonuc


# ============================ olcum ============================
def onkayit197():
    yol = S197 / "ONKAYIT_197.json"
    if not yol.exists():
        raise SystemExit("ONKAYIT_197.json YOK — K0 mühürlenmeden ölçüm başlamaz")
    o = json.load(open(yol))
    s = sha(QM / "197_configs" / "197a_onkayit.py")
    if s != o["sha256"]:
        raise SystemExit(f"ON-KAYIT SHA UYUMSUZ: {s} != {o['sha256']}")
    if sha(KALEM) != o["kalem_sha256"]:
        raise SystemExit("KALEM ön-kayıttan sonra DEĞİŞMİŞ")
    for rel, h in o["girdi_sha256"].items():
        if sha(QM / rel) != h:
            raise SystemExit(f"GİRDİ SHA UYUMSUZ: {rel}")
    p = o["olcum_parametreleri"]
    assert p["havuz"] == list(HAVUZ_197), p["havuz"]
    assert p["tau_yeni"] == list(TAU_YENI) and p["dtau"] == DTAU
    assert p["n_blok"] == NBLOK and p["dw"] == DW and p["om_chunk"] == OM_CHUNK
    assert p["dw_aralik"] == list(DW_ARALIK), p["dw_aralik"]
    kod = sha(Path(__file__).resolve())
    m1 = [S197 / "m1" / "M1_sonuc.json", S197 / "m1" / "M1_sonuc_alt.json"]
    if not any((d := json.load(open(y))).get("gecti") and d.get("kod_sha256") == kod
               for y in m1 if y.exists()):
        raise SystemExit("M1 kapısı bu kod sürümüyle GEÇİLMEDİ — ölçüm başlamaz")
    return o


def kapsanan_tau_dilimleri(kenar_tum, L, Lb_b, jlo, jhi):
    """Blok b için tam dilimlerin [(jlo−½)DW, (jhi+½)DW) aralığına değen τ'-dilimleri
    (bitişik). Bu dilimlerin TÜM çizgileri taranır → doğrusallık denetimi kesin olur."""
    lo_s = kenar_tum[:-1] * L - Lb_b
    hi_s = kenar_tum[1:] * L - Lb_b
    s_ix = np.where((hi_s >= (jlo - 0.5) * DW) & (lo_s < (jhi + 0.5) * DW))[0]
    assert s_ix.size and np.all(np.diff(s_ix) == 1)
    return s_ix


def olcum_plani(c):
    """S1-S3 kurulumu (κ YOK): bloklar, HAVUZ', τ'-dilim ızgarası, taranan çizgiler, kapsama."""
    L = c.L
    kb, Lb = bloklar(c.N, c.mid, NBLOK)
    assert np.array_equal(kb[::NBLOK // NJACK], c.kj8), \
        "S3: 32-blok grup kenarları 190'ın 8-blok kenarlarıyla BİREBİR değil"
    nb = np.diff(kb)
    win1, w1, tau1, q1 = havuz(c, HAVUZ_197)
    kenar_yeni = kenar_yeni_kur()
    assert kenar_yeni[-1] == c.kenar_eski[0]
    kenar_tum = np.concatenate([kenar_yeni[:-1], c.kenar_eski])
    tarama = np.where((c.taua > kenar_tum[0]) & (c.taua <= kenar_tum[-1]))[0]
    m3 = m3_oz_terim(q1, c.qa[tarama], "olcum")
    jlo = int(round(DW_ARALIK[0] / DW))
    jhi = int(round(DW_ARALIK[1] / DW))
    # kapsama (192 kuralı, tarama (0.74, 1.30]): tam aralıktaki her dilimde B_j = 32
    mt = np.arange(jlo, jhi + 1) * DW
    kapT = np.array([(mt - 0.5 * DW >= kenar_tum[0] * L - Lb[b] - EPS) &
                     (mt + 0.5 * DW <= kenar_tum[-1] * L - Lb[b] + EPS) for b in range(NBLOK)])
    Bj = kapT.sum(0)
    assert Bj.min() == NBLOK, f"kapsama eksik: min B_j = {Bj.min()}"
    pay_alt = float(((jlo - 0.5) * DW - (kenar_tum[0] * L - Lb)).min())
    pay_ust = float(((kenar_tum[-1] * L - Lb) - (jhi + 0.5) * DW).min())
    sec, Sb = {}, {}
    for b in range(NBLOK):
        s_ix = kapsanan_tau_dilimleri(kenar_tum, L, Lb[b], jlo, jhi)
        Sb[b] = s_ix
        sec[b] = np.where((c.taua > kenar_tum[s_ix[0]]) &
                          (c.taua <= kenar_tum[s_ix[-1] + 1]))[0]
        jt = np.floor((np.log(c.qa[tarama]) - Lb[b]) / DW + 0.5).astype(int)
        gerek = tarama[(jt >= jlo) & (jt <= jhi)]
        assert np.isin(gerek, sec[b]).all(), f"blok {b}: tam aralık çizgileri eksik"
    print(f"  S1 HAVUZ' {int(win1.sum())} çizgi; S2 yeni {len(kenar_yeni)-1} τ'-dilimi; "
          f"kapsama: {jhi-jlo+1} tam dilimin HEPSİNDE B_j = {int(Bj.min())} "
          f"(pay: alt {pay_alt:.3f}, üst {pay_ust:.3f})", flush=True)
    return SimpleNamespace(kb=kb, Lb=Lb, nb=nb, win1=win1, w1=w1, tau1=tau1, q1=q1,
                           kenar_yeni=kenar_yeni, kenar_tum=kenar_tum, tarama=tarama,
                           m3=m3, jlo=jlo, jhi=jhi, sec=sec, Sb=Sb,
                           kapsama={"n_tam_dilim": jhi - jlo + 1, "B_j_min": int(Bj.min()),
                                    "B_j_maks": int(Bj.max()), "pay_alt": pay_alt,
                                    "pay_ust": pay_ust})


def mod_olcum(c, nw):
    t00 = time.time()
    ONK197 = onkayit197()
    print(f"  [on-kayit {ONK197['zaman']} sha {ONK197['sha256'][:12]}]", flush=True)
    L = c.L
    pl = olcum_plani(c)
    assert pl.Lb.tolist() == ONK197["S3"]["L_b"], "L_b32 ÖN-KAYITLA UYUMSUZ"
    assert int(pl.win1.sum()) == ONK197["S1"]["n_cizgi"]
    grup = np.arange(NBLOK) // (NBLOK // NJACK)
    hv = [np.ones(len(pl.w1), bool)]                        # tek maske: HAVUZ'
    mixes = ([b188.karisim(c.P, c.OZ, c.aq_all)[pl.win1]] +
             [b188.karisim(c.P, c.OZ, c.aq_all, i)[pl.win1] for i in range(NJACK)])
    sonuc = {"sha_onkayit": ONK197["sha256"], "M3": pl.m3, "N": c.N, "L": L,
             "n_havuz_ussu": int(pl.win1.sum()), "jlo_jhi": [pl.jlo, pl.jhi],
             "kapsama": pl.kapsama, "birincil": "K_duz (S1b)"}

    # ---------------- A: S2 yeni τ'-dilimleri (188b.dilimleri_kur AYNEN) ----------------
    ns_y = len(pl.kenar_yeni) - 1
    b188.S188 = S197
    t0 = time.time()
    b188.dilimleri_kur("yeni", c.g, c.mid, c.bid8, pl.kenar_yeni, L, c.qa, c.taua,
                       c.aqa, "dilim_", nw)
    sonuc["sure_A_s"] = time.time() - t0
    dilim_yol = ([S197 / f"dilim_{i}.npz" for i in range(ns_y)] +
                 [S190 / f"dilim_{i}.npz" for i in range(len(c.kenar_eski) - 1)])
    ns = len(dilim_yol)
    assert ns == len(pl.kenar_tum) - 1
    dilim_ix = []
    for s, yol in enumerate(dilim_yol):          # dosya ↔ çizgi seçimi (190 dahil)
        idx = np.where((c.taua > pl.kenar_tum[s]) & (c.taua <= pl.kenar_tum[s + 1]))[0]
        d = np.load(yol)
        assert np.array_equal(d["q"], c.qa[idx]) and np.array_equal(d["a"], c.aqa[idx]), yol
        dilim_ix.append(idx)
    print(f"  [A] τ'-dilimleri: {ns_y} yeni + {ns - ns_y} (190) = {ns}", flush=True)

    # ---------------- S1b katkı makinesi (HAVUZ') + M-kesik denetimi ----------------
    yolM = S197 / "kesik_katki32.npz"
    if not yolM.exists():
        re_k, im_k = kesik_katki(c, pl.w1, nw)
        np.savez_compressed(yolM, re=re_k, im=im_k, w_K1=c.w_all, w_havuz=pl.w1,
                            nb8=np.diff(c.kj8))
    KM = np.load(yolM)
    re_k, im_k = KM["re"], KM["im"]
    mkes = kesik_denetim(c, re_k, im_k, pl.win1)
    print(f"  M-kesik Σ_K1 katkı ≡ c^kesik (HAVUZ'): maks bağıl = {max(mkes):.2e}", flush=True)
    assert max(mkes) <= 1e-10, "katkı makinesi c^kesik'i üretmiyor"
    sonuc["M_kesik"] = mkes
    Mv = [b188.c_proj(re_k, im_k, np.diff(c.kj8), c.N, v - 1) for v in range(NJACK + 1)]
    del re_k, im_k
    katki = katki_kur(c.k1_satir, Mv)

    # ---------------- B: τ'-izdüşüm (32 blok, HAVUZ'; 190b.izdusum_blok AYNEN) ----------------
    yolB = S197 / "harita_proj32.npz"
    if not yolB.exists():
        t0 = time.time()
        D = np.stack([np.load(y)["dds"] for y in dilim_yol], axis=1)
        re_ = np.zeros((NBLOK, ns, len(pl.w1)))
        im_ = np.zeros((NBLOK, ns, len(pl.w1)))
        for b in range(NBLOK):
            lo, hi = pl.kb[b], pl.kb[b + 1]
            re_[b], im_[b] = b190.izdusum_blok(D[lo:hi], c.mid[lo:hi], pl.w1)
        del D
        np.savez_compressed(yolB, re=re_, im=im_, nb=pl.nb, kb=pl.kb,
                            kenar=pl.kenar_tum, w_havuz=pl.w1)
        print(f"  [B] τ'-izdüşüm -> {yolB.name} ({time.time()-t0:.0f}s)", flush=True)
    PR = np.load(yolB)
    PRre, PRim, PRnb = PR["re"], PR["im"], PR["nb"]
    assert np.array_equal(PRnb, pl.nb) and PRre.shape == (NBLOK, ns, len(pl.w1))
    Kt_ham = np.zeros((NJACK + 1, NBLOK, ns), complex)
    Kt_duz = np.zeros_like(Kt_ham)
    kat_s = [[katki(v, dilim_ix[s]) for s in range(ns)] for v in range(NJACK + 1)]
    for b in range(NBLOK):
        Cb = 2.0 * (PRre[b] + 1j * PRim[b]) / int(PRnb[b])
        for v, mx in enumerate(mixes):
            Kt_ham[v, b] = b188.K_matris(Cb, mx, hv)[0][0]
            for s in range(ns):
                a = kat_s[v][s]
                Kt_duz[v, b, s] = (Kt_ham[v, b, s] if a is None else
                                   b188.K_matris(Cb[s:s + 1], mx - a, hv)[0][0, 0])
    np.savez_compressed(S197 / "harita_K32.npz", K_duz=Kt_duz, K_ham=Kt_ham,
                        kenar=pl.kenar_tum, nb=pl.nb, L_b=pl.Lb, kb=pl.kb,
                        kapsanan_s=np.array([[pl.Sb[b][0], pl.Sb[b][-1]]
                                             for b in range(NBLOK)]),
                        etiket_v=np.array(["tam"] + [f"loo{i}" for i in range(NJACK)]))

    # ---------------- C: 32-blok ω-dilimleri (190b D AYNEN; Δω ∈ DW_ARALIK) ----------------
    dizin = S197 / "omega32"
    dizin.mkdir(exist_ok=True)
    isler, plan, say = plan_kur(c.qa, c.aqa, pl.sec, pl.Lb, dizin)
    print(f"  [C] {NBLOK} blok; çizgi/blok {min(len(pl.sec[b]) for b in range(NBLOK))}-"
          f"{max(len(pl.sec[b]) for b in range(NBLOK))}; dilim/blok "
          f"{min(len(plan[b]) for b in range(NBLOK))}-{max(len(plan[b]) for b in range(NBLOK))};"
          f" eksik görev {len(isler)} (işçi {nw})", flush=True)
    sonuc["sure_C_s"] = omega_kos(isler, c, pl.w1, pl.kb, nw, "olcum")
    J, Kh, Kd, n, nkes, nb_om = omega_birlestir(dizin, plan, say, NBLOK, mixes, hv,
                                                pl.sec, pl.Lb, c, katki)
    K_ham, K_duz = Kh[:, :, 0, :], Kd[:, :, 0, :]                # (9, 32, nJ)
    del Kh, Kd
    assert np.array_equal(nb_om, pl.nb)

    # doğrusallık (K_ham; tek mix): Σ_j K_ω(b,j) = Σ_{s∈S_b} K_τ'(b,s)
    dog = [float(abs(K_ham[0, b].sum() - Kt_ham[0, b, pl.Sb[b]].sum()) /
                 abs(Kt_ham[0, b, pl.Sb[b]].sum())) for b in range(NBLOK)]
    print(f"  doğrusallık (K_ham) |ΣK_ω − ΣK_τ'|/|ΣK_τ'| (32 blok) maks = {max(dog):.1e}",
          flush=True)
    sonuc["dogrusallik_K_ham"] = dog

    merkez = J * DW
    tam = (J >= pl.jlo) & (J <= pl.jhi)
    kap = np.zeros((NBLOK, len(J)), bool)
    for b in range(NBLOK):                        # 192 kapsama AYNEN, tarama (0.74, 1.30]
        kap[b] = ((merkez - 0.5 * DW >= pl.kenar_tum[0] * L - pl.Lb[b] - EPS) &
                  (merkez + 0.5 * DW <= pl.kenar_tum[-1] * L - pl.Lb[b] + EPS) & tam)
    assert kap[:, tam].all()
    np.savez_compressed(
        S197 / "harita_omega32.npz", birincil=np.array("duz"),
        K_duz=K_duz, kappa_sigma_duz=profil_sigma(K_duz, pl.nb, kap),
        K_ham=K_ham, kappa_sigma_ham=profil_sigma(K_ham, pl.nb, kap),
        ncz=n, n_kesik=nkes, nb=pl.nb, kb=pl.kb, L_b=pl.Lb, grup=grup, kapsama=kap,
        tam=tam, J=J, merkez=merkez, dw=DW, dw_aralik=np.array(DW_ARALIK),
        havuz=np.array(HAVUZ_197), w_havuz=pl.w1,
        etiket_v=np.array(["tam"] + [f"loo{i}" for i in range(NJACK)]))
    sonuc.update({"nJ": int(len(J)), "n_tam_dilim": int(tam.sum()),
                  "n_kesik_hucre": int((nkes > 0).sum()),
                  "sure_toplam_s": time.time() - t00})
    json.dump(sonuc, open(S197 / "K1_harita_197.json", "w"), indent=1, ensure_ascii=False)
    print(f"-> harita_omega32.npz, harita_K32.npz, harita_proj32.npz, kesik_katki32.npz, "
          f"K1_harita_197.json BİTTİ ({time.time()-t00:.0f}s) — κ/P EKRANA YAZILMADI.",
          flush=True)


# ============================ zamanlama ============================
def log190_omega():
    t, n = T190_OMEGA_YEDEK, None
    try:
        s = (S190 / "log_190b.txt").read_text()
        ts = re.findall(r"\[ω\] blok \d+ parça \d+ bitti .*? toplam (\d+)s", s)
        if ts:
            t = float(ts[-1])
        m = re.search(r"pencere-ötesi \(τ'∈\(0\.86,1\.3\]\): (\d+)", s)
        if m:
            n = int(m.group(1))
    except OSError:
        pass
    return t, n


def mod_zamanlama(c, nw, n_dilim=2):
    t00 = time.time()
    pl = olcum_plani(c)
    ncz_y = np.array([int(((c.taua > pl.kenar_yeni[i]) & (c.taua <= pl.kenar_yeni[i + 1])).sum())
                      for i in range(len(pl.kenar_yeni) - 1)])
    nl_b = np.array([len(pl.sec[b]) for b in range(NBLOK)])
    nsl_b, nparca = [], 0
    for b in range(NBLOK):
        j = np.floor((np.log(c.qa[pl.sec[b]]) - pl.Lb[b]) / DW + 0.5).astype(int)
        uj, cnt = np.unique(j, return_counts=True)
        nsl_b.append(len(uj))
        nparca += int(np.ceil(cnt.sum() / OM_CHUNK))
    pt_line = float(np.sum(pl.nb * nl_b))
    n_kes = int((c.k1_satir[pl.tarama] >= 0).sum())
    print(f"  S2 yeni τ'-dilimleri: {len(ncz_y)} dilim, {int(ncz_y.sum())} çizgi "
          f"(dilim başına {ncz_y.min()}-{ncz_y.max()}); tarama ∩ c^kesik evreni: {n_kes}",
          flush=True)
    print(f"  S3 ω: 32 blok, blok başına çizgi {nl_b.min()}-{nl_b.max()}, dilim "
          f"{min(nsl_b)}-{max(nsl_b)}; ~{nparca} görev; Σ n_b·çizgi = {pt_line:.3e}", flush=True)
    # S2 zamanlaması (tam N, 188b.dilimleri_kur AYNEN, ayrı dizin, tek işçi)
    zd = S197 / "zamanlama" / time.strftime("%Y%m%d_%H%M%S")
    zd.mkdir(parents=True, exist_ok=False)
    b188.S188 = zd
    kz = pl.kenar_yeni[-(n_dilim + 1):]
    nl_test = int(sum(((c.taua > kz[i]) & (c.taua <= kz[i + 1])).sum() for i in range(n_dilim)))
    t0 = time.time()
    b188.dilimleri_kur("zaman", c.g, c.mid, c.bid8, kz, c.L, c.qa, c.taua, c.aqa, "zdilim_", 1)
    t_s2 = time.time() - t0
    cpu_s2 = t_s2 / nl_test * ncz_y.sum()
    # katkı makinesi (HAVUZ'; karışım düzeyi — κ değil)
    t0 = time.time()
    kesik_katki(c, pl.w1, nw)
    t_kk = time.time() - t0
    t190, n190 = log190_omega()
    if n190 is None:
        n190 = int(((c.taua > c.kenar_eski[0]) & (c.taua <= c.kenar_eski[-1])).sum())
    duvar_om = pt_line / (c.N * n190 / t190)
    toplam = cpu_s2 / nw + t_kk + duvar_om + 60
    print(f"  [S2] {n_dilim} dilim ({nl_test} çizgi) {t_s2:.1f}s → yeni dilimler CPU ≈ "
          f"{cpu_s2:.0f}s; katkı makinesi {t_kk:.0f}s; ω (190 kalibre) ≈ {duvar_om:.0f}s; "
          f"TOPLAM ≈ {toplam/60:.1f} dk", flush=True)
    out = {"S2_dilim": int(len(ncz_y)), "S2_cizgi": int(ncz_y.sum()), "S2_test_s": t_s2,
           "S2_cpu_s": cpu_s2, "katki_s": t_kk, "omega_duvar_s": duvar_om,
           "nokta_cizgi": pt_line, "cizgi_blok": nl_b.tolist(), "dilim_blok": nsl_b,
           "n_havuz_ussu": int(pl.win1.sum()), "kapsama": pl.kapsama, "M3": pl.m3,
           "tarama_kesik_evreninde": n_kes, "toplam_tahmin_s": toplam,
           "sure_test_s": time.time() - t00}
    json.dump(out, open(zd / "zamanlama.json", "w"), indent=1, ensure_ascii=False)
    print(f"  -> {zd/'zamanlama.json'}", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("mod", choices=["m1", "olcum", "zamanlama"])
    ap.add_argument("--alt", action="store_true", help="m1: 12 temsilî ω-dilimi")
    ap.add_argument("--isci", type=int, default=6)
    a = ap.parse_args()
    S197.mkdir(exist_ok=True)
    print("=" * 78)
    print(f"197b / {a.mod.upper()}{' (alt)' if a.alt else ''}  işçi={a.isci}  "
          f"{time.strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 78, flush=True)
    C = ortak()
    if a.mod == "m1":
        mod_m1(C, a.alt, a.isci)
    elif a.mod == "olcum":
        mod_olcum(C, a.isci)
    else:
        mod_zamanlama(C, a.isci)
