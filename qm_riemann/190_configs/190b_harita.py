# -*- coding: utf-8 -*-
"""
190b — K1: İPTAL ÇEKİRDEĞİ HARİTASI, DÜŞÜK PENCERE (188b makinesi AYNEN)
======================================================================
ONKAYIT_190 dondurdu. 188b_harita.py importlib ile yüklenir (DOSYA DÜZENLENMEZ);
fonksiyonları (seri_ve_G, izdusum, c_proj, K_matris, G_loo, bant_maskeleri,
karisim, dilimleri_kur, jk_kenar) AYNEN çağrılır; yalnız çalışma-anı yaması
b188.S188 → scratchpad/190 (dilim kontrol noktaları). Girdi zinciri 190k0
'dusuk' (scratchpad/190/zincir_dusuk). TAM örneklem (ALT=1).

Aşamalar:
 M1  seri_ve_G ≡ 187b.katman_serisi (düşük, dilim 0, tam örneklem) bit-bit
 M2  188'in dilim_0.npz'i (son, alt-örneklem) bu kod yoluyla bit-bit
 A   88 τ'-dilim serisi (kontrol noktalı; çok süreçli)            [H-190b]
 B   tek geçişte izdüşüm → harita_proj_dusuk.npz
 C   K (41 ince + 8 bant + HAVUZ) tam + loo, S, S_Re, W → harita_K_dusuk.npz
 D   blok-blok Δω ω-dilimleri (0.025; blok içi seri + izdüşüm) [H-190a]
 M3  doğrusallık: Σ_j K_HAVUZ(b,j) = Σ_s K_HAVUZ(b,s) her blokta
Çıktı: 190/dilim_<i>.npz, harita_proj_dusuk.npz, harita_K_dusuk.npz,
       omega/b<b>_c<k>.npz, harita_omega_dusuk.npz, K1_harita_dusuk.json
Kullanım: 190b_harita.py [işçi_sayısı]
"""
import hashlib
import importlib.util
import json
import os
import sys
import time
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
           "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S155, S188, S190 = SCR / "155", SCR / "188", SCR / "190"
ZD = S190 / "zincir_dusuk"
OM = S190 / "omega"
TWO_PI = 2 * np.pi
NJACK = 8
OM_CHUNK = 8000      # ω-görevi başına yaklaşık çizgi sayısı (dilimler bölünmez)


def yukle(ad, yol):
    spec = importlib.util.spec_from_file_location(ad, yol)
    m = importlib.util.module_from_spec(spec)
    sys.modules[ad] = m
    spec.loader.exec_module(m)
    return m


b188 = yukle("b188", QM / "188_configs" / "188b_harita.py")
b185 = b188.b185
b187 = b188.b187          # 187b_katman_defteri (188b'nin kullandığı örnek)


def onkayit():
    o = json.load(open(S190 / "ONKAYIT_190.json"))
    sha = hashlib.sha256(
        (QM / "190_configs" / "190a_onkayit.py").read_bytes()).hexdigest()
    if sha != o["sha256"]:
        raise SystemExit(f"ON-KAYIT SHA UYUMSUZ: {sha} != {o['sha256']}")
    return o


def kinematik_eta(yol):
    """185b 'gercek' kolu AYNEN: g = (ds+1)·2π/log(mid/2π)."""
    d = np.load(yol)
    mid = np.asarray(d["mid"], float)
    ds = np.asarray(d["ds"], float)
    L = float(d["L"])
    g = (ds + 1.0) * TWO_PI / np.log(mid / TWO_PI)
    return g, mid, L


def izdusum_blok(D, mid_b, w_win):
    """188b.izdusum'un tek-blok iç döngüsü AYNEN (PROJ_BLOK nokta bloğu)."""
    ns = D.shape[1]
    re = np.zeros((ns, len(w_win)))
    im = np.zeros((ns, len(w_win)))
    n = len(mid_b)
    for s0 in range(0, n, b188.PROJ_BLOK):
        sl = slice(s0, min(s0 + b188.PROJ_BLOK, n))
        P = np.outer(mid_b[sl], w_win)
        cP = np.cos(P)
        sP = np.sin(P)
        del P
        re += D[sl].T @ cP
        im -= D[sl].T @ sP
        del cP, sP
    return re, im


# ---------------- ω-görevi işçisi ----------------
_O = {}


def _om_init(g, mid, bid, w_win, kenar_jk):
    _O.update(g=g, mid=mid, bid=bid, w_win=w_win, kj=kenar_jk)


def _om_is(arg):
    b, k, yol, jler, w_list, a_list = arg
    t0 = time.time()
    lo, hi = _O["kj"][b], _O["kj"][b + 1]
    g_b, mid_b, bid_b = _O["g"][lo:hi], _O["mid"][lo:hi], _O["bid"][lo:hi]
    D = np.zeros((hi - lo, len(jler)))
    ncz = np.zeros(len(jler), int)
    for c, (w_l, a_l) in enumerate(zip(w_list, a_list)):
        D[:, c], _, _ = b188.seri_ve_G(g_b, mid_b, bid_b, w_l, a_l)
        ncz[c] = len(w_l)
    re, im = izdusum_blok(D, mid_b, _O["w_win"])
    tmp = str(yol) + ".tmp.npz"
    np.savez(tmp, re=re, im=im, j=np.array(jler), ncz=ncz, nb=hi - lo, blok=b)
    os.replace(tmp, yol)
    return b, k, int(ncz.sum()), time.time() - t0


def main(nw):
    ONK = onkayit()
    ONK188 = json.load(open(S188 / "ONKAYIT_188.json"))
    kenar = np.array(ONK["dilim_izgara"])
    assert np.array_equal(kenar, np.array(ONK188["dilim_izgara_gercek"]))
    assert ONK["ince_bant_kenar"] == ONK188["ince_bant_kenar"]
    assert ONK["bant8_kenar"] == ONK188["bant8_kenar"]
    L = float(ONK["pencere"]["L"])
    ns = len(kenar) - 1
    print("=" * 78)
    print(f"190b / K1 HARİTA (düşük pencere)  [on-kayit {ONK['zaman']} sha "
          f"{ONK['sha256'][:12]}]")
    print("=" * 78, flush=True)

    G = np.load(ZD / "K1_gercek_dusuk.npz")
    OZ = np.load(ZD / "OZ_gercek_dusuk.npz")
    P = np.load(ZD / "G1_proj_gercek_dusuk.npz")
    w_all = np.asarray(G["w"], float)
    tau_all = np.asarray(G["tau"], float)
    aq_all = np.asarray(G["aq"], float)
    win = (tau_all >= 0.45) & (tau_all < 0.86)
    w_win, tau_win = w_all[win], tau_all[win]

    g, mid, Lg = kinematik_eta(ZD / "eta_dusuk_t0.4_c4000.npz")
    assert Lg == L, (Lg, L)
    N = len(mid)
    kenar_jk = b188.jk_kenar(N)
    if not (np.array_equal(np.diff(kenar_jk), G["nb"]) and
            np.array_equal(np.diff(kenar_jk), P["nb"])):
        raise SystemExit("JK BLOK KENARLARI 184 İLE UYUMSUZ")
    bid = np.searchsorted(kenar_jk, np.arange(N), side="right") - 1
    L_b = np.array([float(np.log(mid[kenar_jk[b]:kenar_jk[b + 1]] / TWO_PI).mean())
                    for b in range(NJACK)])
    assert np.array_equal(L_b, np.array(ONK["bloklar"]["L_b"])), "L_b ÖN-KAYITLA UYUMSUZ"
    print(f"[dusuk] N={N} L={L:.9f} pencere çizgi={int(win.sum())} dilim={ns} "
          f"(TAM örneklem)", flush=True)

    qmax = np.exp(kenar[-1] * L)
    t0 = time.time()
    qa, lam, taua, aqa = b187.asal_kuvvetler(qmax, L)
    otesi = (taua > kenar[0]) & (taua <= kenar[-1])
    print(f"  asal-kuvvetler q≤{qmax:.0f}: {len(qa)} ({time.time()-t0:.0f}s); "
          f"pencere-ötesi (τ'∈({kenar[0]},{kenar[-1]}]): {int(otesi.sum())}",
          flush=True)
    sonuc = {"sha_onkayit": ONK["sha256"], "L": L, "N": N,
             "n_asal_kuvvet": int(len(qa)), "n_otesi": int(otesi.sum()),
             "n_pencere": int(win.sum()), "L_b": L_b.tolist()}

    # ---------------- M1 ----------------
    idx0 = np.where((taua > kenar[0]) & (taua <= kenar[1]))[0]
    d_ref = b187.katman_serisi(g, mid, np.log(qa[idx0]), aqa[idx0], "M1")
    d_yeni, _, _ = b188.seri_ve_G(g, mid, bid, np.log(qa[idx0]), aqa[idx0])
    M1 = float(np.max(np.abs(d_ref - d_yeni)))
    print(f"  M1 maks|seri_ve_G − 187b.katman_serisi| (dilim 0, {len(idx0)} çizgi) "
          f"= {M1:.1e}", flush=True)
    sonuc["M1"] = M1
    del d_ref, d_yeni

    # ---------------- M2 (son, 188 dilim_0 yeniden üretimi) ----------------
    ONK188L = float(ONK188["L"])
    gs, ms, Ls = kinematik_eta(S155 / "eta_son_t0.4_c4000.npz")
    Ns = len(ms)
    kjs = b188.jk_kenar(Ns)
    bids = np.searchsorted(kjs, np.arange(Ns), side="right") - 1
    alt = np.arange(0, Ns, b188.ALT)
    qs_, _, ts_, as_ = b187.asal_kuvvetler(np.exp(1.30 * ONK188L), ONK188L)
    i0 = np.where((ts_ > 0.86) & (ts_ <= 0.865))[0]
    dd, Gre, Gim = b188.seri_ve_G(gs[alt], ms[alt], bids[alt], np.log(qs_[i0]),
                                  as_[i0])
    R = np.load(S188 / "dilim_0.npz")
    M2 = {"dds": bool(np.array_equal(dd, R["dds"])),
          "Gre": bool(np.array_equal(Gre, R["Gre"])),
          "Gim": bool(np.array_equal(Gim, R["Gim"])),
          "q": bool(np.array_equal(qs_[i0], R["q"])),
          "a": bool(np.array_equal(as_[i0], R["a"])),
          "tau": bool(np.array_equal(ts_[i0], R["tau"]))}
    M2["HEPSI"] = all(M2.values())
    print(f"  M2 188 dilim_0 (son, alt-örneklem, {len(i0)} çizgi) bit-bit: {M2}",
          flush=True)
    sonuc["M2"] = M2
    del gs, ms, bids, qs_, ts_, as_, dd, Gre, Gim, R

    # ---------------- A: τ'-dilim serileri ----------------
    b188.S188 = S190
    b188.dilimleri_kur("dusuk", g, mid, bid, kenar, L, qa, taua, aqa, "dilim_", nw)

    # ---------------- B: izdüşüm ----------------
    yolB = S190 / "harita_proj_dusuk.npz"
    if not yolB.exists():
        t0 = time.time()
        D = np.stack([np.load(S190 / f"dilim_{i}.npz")["dds"] for i in range(ns)],
                     axis=1)
        re, im, nb_a, N_a = b188.izdusum(D, mid, bid, w_win)
        del D
        np.savez_compressed(yolB, re=re, im=im, nb=nb_a, N=N_a, kenar=kenar)
        print(f"  [izdüşüm] -> {yolB.name} ({time.time()-t0:.0f}s)", flush=True)
    PR = np.load(yolB)
    re, im, nb_a, N_a = PR["re"], PR["im"], PR["nb"], int(PR["N"])

    # ---------------- C: K, S, S_Re, W (188b.kos gövdesi AYNEN) ----------------
    Gre_s, Gim_s, a_s, tau_s, ncz = [], [], [], [], []
    for i in range(ns):
        d = np.load(S190 / f"dilim_{i}.npz")
        Gre_s.append(d["Gre"])
        Gim_s.append(d["Gim"])
        a_s.append(d["a"])
        tau_s.append(d["tau"])
        ncz.append(len(d["a"]))

    def SW(disari):
        S = np.zeros(ns)
        SRe = np.zeros(ns)
        for i in range(ns):
            Gh = b188.G_loo(Gre_s[i], Gim_s[i], nb_a, N_a, disari)
            S[i] = np.sum(a_s[i] * np.abs(Gh))
            SRe[i] = np.sum(a_s[i] * np.pi * tau_s[i] *
                            np.cos(np.pi * tau_s[i]) * Gh.real)
        return S, SRe

    W = np.array([a.sum() for a in a_s])
    S_tam, SRe_tam = SW(-1)
    S_reps = np.array([SW(j)[0] for j in range(NJACK)])
    SRe_reps = np.array([SW(j)[1] for j in range(NJACK)])

    maskeler, et = b188.bant_maskeleri(tau_win, ONK188)

    def K_hesap(disari):
        C = b188.c_proj(re, im, nb_a, N_a, disari)
        mix = b188.karisim(P, OZ, aq_all, disari)[win]
        return b188.K_matris(C, mix, maskeler)

    K_tam, pay_tam = K_hesap(-1)
    rep = [K_hesap(j) for j in range(NJACK)]
    K_reps = np.array([r[0] for r in rep])
    pay_reps = np.array([r[1] for r in rep])
    orta = 0.5 * (kenar[:-1] + kenar[1:])
    iH = len(et) - 1
    for ad, m in [("le120", kenar[1:] <= 1.20 + 1e-9), ("le130", np.ones(ns, bool))]:
        z = K_tam[iH][m].sum()
        zr = K_reps[:, iH][:, m].sum(-1)
        se_mod = float(b185.jk_se(np.abs(zr)))
        dfark = (np.degrees(np.angle(zr)) - np.degrees(np.angle(z)) + 180) % 360 - 180
        se_aci = float(np.sqrt((NJACK - 1) / NJACK * np.sum((dfark - dfark.mean()) ** 2)))
        sonuc[f"HAVUZ_toplam_{ad}"] = {"mod": float(abs(z)), "se_mod": se_mod,
                                        "aci": float(np.degrees(np.angle(z))),
                                        "se_aci": se_aci}
        print(f"  HAVUZ Σ_s K ({ad}): {abs(z):.4f}±{se_mod:.4f} "
              f"∠{np.degrees(np.angle(z)):+.2f}°±{se_aci:.2f}", flush=True)
    np.savez_compressed(
        S190 / "harita_K_dusuk.npz", K=K_tam, K_reps=K_reps, pay=pay_tam,
        pay_reps=pay_reps, kenar=kenar, orta=orta, etiket=np.array(et),
        S=S_tam, S_reps=S_reps, SRe=SRe_tam, SRe_reps=SRe_reps, W=W,
        ncizgi=np.array(ncz), ince_kenar=np.array(ONK188["ince_bant_kenar"]))
    sonuc["dilim_cizgi"] = ncz
    print("  -> harita_K_dusuk.npz", flush=True)

    # ---------------- D: blok-blok Δω ω-dilimleri ----------------
    OM.mkdir(exist_ok=True)
    DW = float(ONK["omega_dilim"]["genislik"])
    io = np.where(otesi)[0]
    lq = np.log(qa[io])
    isler, plan = [], {}
    for b in range(NJACK):
        j = np.floor((lq - L_b[b]) / DW + 0.5).astype(int)
        uj = np.unique(j)
        plan[b] = uj
        grup, sayac = [], 0
        k = 0
        for jj in uj:
            sel = io[j == jj]
            grup.append((int(jj), np.log(qa[sel]), aqa[sel]))
            sayac += len(sel)
            if sayac >= OM_CHUNK or jj == uj[-1]:
                yol = OM / f"b{b}_c{k}.npz"
                if not yol.exists():
                    isler.append((b, k, yol, [x[0] for x in grup],
                                  [x[1] for x in grup], [x[2] for x in grup]))
                grup, sayac = [], 0
                k += 1
    print(f"  [ω] {NJACK} blok; dilim sayıları {[len(plan[b]) for b in range(NJACK)]};"
          f" eksik görev {len(isler)} (işçi {nw})", flush=True)
    if isler:
        import multiprocessing as mp
        isler.sort(key=lambda x: -sum(len(w_) for w_ in x[4]))
        t0 = time.time()
        ctx = mp.get_context("fork")
        with ctx.Pool(nw, initializer=_om_init,
                      initargs=(g, mid, bid, w_win, kenar_jk)) as pool:
            kk = 0
            for b, k, nl, dt in pool.imap_unordered(_om_is, isler):
                kk += 1
                print(f"    [ω] blok {b} parça {k} bitti ({nl} çizgi, {dt:.0f}s)  "
                      f"{kk}/{len(isler)}  toplam {time.time()-t0:.0f}s", flush=True)

    # birleştir: ortak j ızgarası
    jmin = min(int(plan[b].min()) for b in range(NJACK))
    jmax = max(int(plan[b].max()) for b in range(NJACK))
    J = np.arange(jmin, jmax + 1)
    nJ = len(J)
    mix_tam = b188.karisim(P, OZ, aq_all)[win]
    K_om = np.zeros((NJACK, len(maskeler), nJ), complex)
    n_om = np.zeros((NJACK, nJ), int)
    nb_om = np.zeros(NJACK, int)
    for b in range(NJACK):
        parca = sorted(OM.glob(f"b{b}_c*.npz"),
                       key=lambda p: int(p.stem.split("_c")[1]))
        for p in parca:
            d = np.load(p)
            C = 2.0 * (d["re"] + 1j * d["im"]) / int(d["nb"])
            Kb, _ = b188.K_matris(C, mix_tam, maskeler)
            ix = d["j"] - jmin
            K_om[b][:, ix] = Kb
            n_om[b, ix] = d["ncz"]
            nb_om[b] = int(d["nb"])
        assert n_om[b].sum() == len(io), (b, n_om[b].sum(), len(io))
    # M3: doğrusallık (τ'-dilim blok izdüşümü ile)
    M3 = []
    for b in range(NJACK):
        Cb = 2.0 * (re[b] + 1j * im[b]) / nb_a[b]
        Ks, _ = b188.K_matris(Cb, mix_tam, maskeler)
        a_ = K_om[b, iH].sum()
        b_ = Ks[iH].sum()
        M3.append(float(abs(a_ - b_) / abs(b_)))
    print(f"  M3 doğrusallık |ΣK_ω − ΣK_τ'|/|ΣK_τ'| (8 blok, HAVUZ) maks = "
          f"{max(M3):.1e}", flush=True)
    sonuc["M3"] = M3
    np.savez_compressed(S190 / "harita_omega_dusuk.npz", K=K_om, ncz=n_om, J=J,
                        dw=DW, merkez=J * DW, L_b=L_b, nb=nb_om,
                        etiket=np.array(et))
    json.dump(sonuc, open(S190 / "K1_harita_dusuk.json", "w"), indent=1,
              ensure_ascii=False)
    print("-> harita_K_dusuk.npz, harita_omega_dusuk.npz, K1_harita_dusuk.json "
          "BİTTİ", flush=True)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 6)
