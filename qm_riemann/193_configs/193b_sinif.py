# -*- coding: utf-8 -*-
"""
193b — K1: KALINTI-SINIFI ÖLÇÜMÜ (blok-yerel Δω seçimi; 192c/188b makinesi AYNEN)
=================================================================================
ONKAYIT_193 dondurdu. 188b_harita.py (seri_ve_G, izdusum, c_proj, karisim,
K_matris, bant_maskeleri, jk_kenar) ve 190b_harita.py (izdusum_blok — tek-blok
izdüşüm deseni) importlib ile yüklenir (DOSYALAR DÜZENLENMEZ), AYNEN çağrılır.

YALNIZ SEÇİM değişti (192c'ye göre): 192c'nin TÜM pencerede ortak (blok-bağımsız)
τ'-penceresi yerine, HER BLOK kendi yerel L_b'siyle Δω_b(q') = log q' − L_b
hesaplar; uydu n için |Δω_b(q') − log n| < δ olan çizgiler O BLOĞUN KENDİ
noktalarıyla (g_b, mid_b — yalnız blok b) sınıflandırılıp seri/izdüşüm
hesaplanır (190b'nin ω-dilim / tek-blok izdüşüm deseni). Sonra 8 bloğun re/im
katkıları 188b.c_proj/K_matris'e AYNEN (blok-loo destekli, jackknife) verilir —
yani blok-yerel D matrisi, biçimce 188b.izdusum'un ürettiği (NJACK, nsınıf,
nwin) diziyle ÖZDEŞ; geri kalan makine (c_proj, K_matris, κ=−Re K, jackknife)
AYNEN 192c/188b.

 M   MAKİNE MÜHRÜ: 188/dilim_44 bit-bit (192c AYNEN) — seri_ve_G'nin 188 ile
     bit-bit aynılığı + tek-dilim K karşılaştırması + mod-3 doğrusallık.
 A   5 hedef × (sınıf + 'bölünen') × 8 blok: blok-yerel seçim + blok-yerel
     seri + blok-yerel izdüşüm → RE/IM (8, nsınıf, nwin) [kontrol noktalı].
 B   K_r = K_matris(c_proj(RE,IM,nb_a,N_a,disari), karışım, maskeler)[HAVUZ];
     κ_r = −Re K_r; κ_top = Σ_{r∈sınıflar} κ_r (BÖLÜNEN HARİÇ, ayrı raporlanır);
     s_r = κ_r/κ_top; 8-blok jackknife (loo).
 C   +log10 r=1 sınıfında TAM-örneklem kontrolü (blok-yerel, ALT=1, N=299999).
Çıktı: 193/proj_<ad>.npz, 193/muhur_193.json, 193/tam_plog10_r1.npz,
       193/K1_sinif.json, 193/K1_sinif.npz
Koşu: nohup (asal-kuvvetler + ~184 blok-yerel seri çağrısı; birkaç dakika
sürebilir — 5 dk'lık ön-plan komutu YASAK).
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
    os.environ[_v] = "1"

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S155, S184, S185, S186, S188, S193 = (SCR / d for d in
                                      ["155", "184", "185", "186", "188", "193"])
NJACK = 8
DELTA = 0.03
TWO_PI = 2 * np.pi


def yukle(ad, yol):
    spec = importlib.util.spec_from_file_location(ad, yol)
    m = importlib.util.module_from_spec(spec)
    sys.modules[ad] = m
    spec.loader.exec_module(m)
    return m


b188 = yukle("b188", QM / "188_configs" / "188b_harita.py")
b185 = b188.b185
b187 = b188.b187
b190 = yukle("b190", QM / "190_configs" / "190b_harita.py")   # izdusum_blok


def onkayit():
    o = json.load(open(S193 / "ONKAYIT_193.json"))
    s = hashlib.sha256((QM / "193_configs" / "193a_onkayit.py").read_bytes()).hexdigest()
    if s != o["sha256"]:
        raise SystemExit(f"ON-KAYIT SHA UYUMSUZ: {s} != {o['sha256']}")
    return o


def jk(r):
    r = np.asarray(r, float)
    return float(np.sqrt((NJACK - 1) / NJACK * np.sum((r - r.mean()) ** 2)))


def log(*a):
    print(*a, flush=True)


def blok_yerel_D(ad, h, L_B, qi, qa, taua, aqa, kenar_jk, g, mid, bid, tamörneklem=False):
    """Her blok kendi Δω_b = log q' − L_b'sinde |Δω_b − log n| < δ seçer, O
    bloğun kendi noktalarıyla (190b tek-blok deseni) sınıf serisi + izdüşüm
    hesaplar. Döner: RE, IM (NJACK, nsınıf, nwin), ncizgi (nsınıf,), sınıf adları,
    taua_aralık (min,max — pencere-ötesi/evren kontrolü)."""
    a = h["a"]
    R_list = h["siniflar"]
    log_n = h["log_n"]
    siniflar_ad = [f"r{r}" for r in R_list] + ["bolunen"]
    nS = len(siniflar_ad)
    RE = None
    IM = None
    ncizgi = np.zeros(nS, int)
    tau_min, tau_max = np.inf, -np.inf
    for b in range(NJACK):
        Lb = L_B[b]
        dW = np.log(qi) - Lb
        sel = np.abs(dW - log_n) < DELTA
        if sel.any():
            tau_min = min(tau_min, float(taua[sel].min()))
            tau_max = max(tau_max, float(taua[sel].max()))
        qb = qi[sel]
        wl_all = np.log(qa[sel])
        al_all = aqa[sel]
        r_masks = [(qb % a == r) for r in R_list]
        bol_mask = np.gcd(qb, a) > 1
        kap = np.zeros(len(qb), bool)
        for m_ in r_masks:
            assert not (kap & m_).any(), (ad, b, "sınıflar örtüşüyor")
            kap |= m_
        assert not (kap & bol_mask).any(), (ad, b, "bölünen sınıflarla örtüşüyor")
        kap |= bol_mask
        assert kap.all(), (ad, b, "sınıflandırılamayan çizgi var")

        lo, hi = kenar_jk[b], kenar_jk[b + 1]
        g_b, mid_b, bid_b = g[lo:hi], mid[lo:hi], bid[lo:hi]
        assert np.all(bid_b == b)

        cols = []
        for k_, m_ in enumerate(r_masks + [bol_mask]):
            ncizgi[k_] += int(m_.sum())
            if m_.any():
                d_, _, _ = b188.seri_ve_G(g_b, mid_b, bid_b, wl_all[m_], al_all[m_])
            else:
                d_ = np.zeros(len(mid_b))
            cols.append(d_)
        D_b = np.stack(cols, axis=1)
        re_b, im_b = b190.izdusum_blok(D_b, mid_b, w_win_G[0])
        if RE is None:
            RE = np.zeros((NJACK, nS, re_b.shape[1]))
            IM = np.zeros((NJACK, nS, re_b.shape[1]))
        RE[b] = re_b
        IM[b] = im_b
    return RE, IM, ncizgi, siniflar_ad, (tau_min, tau_max)


w_win_G = [None]   # basit taşıyıcı (blok_yerel_D içinde kapatma için)


if __name__ == "__main__":
    T0 = time.time()
    ONK = onkayit()
    ONK188 = json.load(open(S188 / "ONKAYIT_188.json"))
    L = float(ONK188["L"])
    HEDEFLER = ONK["hedefler"]
    L_B = np.array(ONK["bloklar"]["L_b"])
    log("=" * 78)
    log(f"193b / K1 SINIF ÖLÇÜMÜ (blok-yerel Δω)  [on-kayit {ONK['zaman']} sha "
        f"{ONK['sha256'][:12]}]")
    log("=" * 78)

    # ---- kinematik + pencere (188b/192c AYNEN) ----
    G = np.load(S184 / "K1_gercek.npz")
    OZ = np.load(S185 / "OZ_gercek.npz")
    P = np.load(S186 / "G1_proj_gercek.npz")
    w_all = np.asarray(G["w"], float)
    tau_all = np.asarray(G["tau"], float)
    aq_all = np.asarray(G["aq"], float)
    win = (tau_all >= 0.45) & (tau_all < 0.86)
    w_win, tau_win = w_all[win], tau_all[win]
    w_win_G[0] = w_win

    g, mid, Lg = b185.kinematik("gercek")
    N = len(mid)
    kenar_jk = b188.jk_kenar(N)
    assert np.array_equal(np.diff(kenar_jk), G["nb"])
    bid_tam = np.searchsorted(kenar_jk, np.arange(N), side="right") - 1
    alt = np.arange(0, N, b188.ALT)
    g_a, mid_a, bid_a = g[alt], mid[alt], bid_tam[alt]
    # ALT dizisindeki GERÇEK blok sınırları: bid_a sıralı (monoton artan) —
    # sınırları doğrudan bid_a'dan çıkar (kendi linspace'imizi İCAT ETMEYİZ).
    assert np.all(np.diff(bid_a) >= 0), "bid_a sıralı değil"
    kenar_jk_a = np.searchsorted(bid_a, np.arange(NJACK + 1))
    assert kenar_jk_a[0] == 0 and kenar_jk_a[-1] == len(alt)
    assert np.array_equal(bid_a, np.searchsorted(kenar_jk_a, np.arange(len(alt)),
                                                  side="right") - 1)
    maskeler, et = b188.bant_maskeleri(tau_win, ONK188)
    iH = len(et) - 1
    assert et[iH] == "HAVUZ"
    nb_a = np.diff(kenar_jk_a)
    N_a = len(alt)
    log(f"  N={N} alt={len(alt)} nb_a={nb_a.tolist()} L={L:.9f}")

    # L_b tutarlılık (K0'daki kinematikle bit-bit aynı olmalı)
    L_B_check = np.array([float(np.log(mid[kenar_jk[b]:kenar_jk[b + 1]] / TWO_PI).mean())
                          for b in range(NJACK)])
    assert np.allclose(L_B_check, L_B, atol=1e-12), "L_B ON-KAYITLA UYUMSUZ"

    t0 = time.time()
    qa, lam, taua, aqa = b187.asal_kuvvetler(np.exp(1.30 * L), L)
    qi = qa.astype(np.int64)
    log(f"  asal-kuvvetler: {len(qa)} ({time.time()-t0:.0f}s)")

    mix_reps = [b188.karisim(P, OZ, aq_all, j)[win] for j in range(NJACK)]
    mix_tam = b188.karisim(P, OZ, aq_all)[win]

    def K_hav(re, im, nb, Nn, disari):
        C = b188.c_proj(re, im, nb, Nn, disari)
        mix = mix_tam if disari < 0 else mix_reps[disari]
        return b188.K_matris(C, mix, maskeler)[0]

    muhur = {}
    # ================= M: makine mührü (dilim 44; 192c AYNEN) =================
    kenar = np.array(ONK188["dilim_izgara_gercek"])
    i44 = 44
    idx = np.where((taua > kenar[i44]) & (taua <= kenar[i44 + 1]))[0]
    t0 = time.time()
    dd, Gre, Gim = b188.seri_ve_G(g_a, mid_a, bid_a, np.log(qa[idx]), aqa[idx])
    R = np.load(S188 / f"dilim_{i44}.npz")
    bit = {"dds": bool(np.array_equal(dd, R["dds"])),
           "Gre": bool(np.array_equal(Gre, R["Gre"])),
           "Gim": bool(np.array_equal(Gim, R["Gim"])),
           "q": bool(np.array_equal(qa[idx], R["q"])),
           "a": bool(np.array_equal(aqa[idx], R["a"])),
           "tau": bool(np.array_equal(taua[idx], R["tau"]))}
    bit["HEPSI"] = all(bit.values())
    log(f"  MÜHÜR dilim_44 ({len(idx)} çizgi) bit-bit: {bit}  "
        f"({time.time()-t0:.0f}s)")
    muhur["dilim44_bitbit"] = bit
    re1, im1, nb_a1, N_a1 = b188.izdusum(dd[:, None], mid_a, bid_a, w_win)
    K1t = K_hav(re1, im1, nb_a1, N_a1, -1)[:, 0]
    K1r = np.array([K_hav(re1, im1, nb_a1, N_a1, j)[:, 0] for j in range(NJACK)])
    HK = np.load(S188 / "harita_K_gercek.npz")
    dK = float(np.max(np.abs(K1t - HK["K"][:, i44])))
    dKr = float(np.max(np.abs(K1r - HK["K_reps"][:, :, i44])))
    muhur["dilim44_K_maks_fark_tam"] = dK
    muhur["dilim44_K_maks_fark_reps"] = dKr
    log(f"  MÜHÜR tek-dilim K vs 188 harita_K: maks|Δ| tam {dK:.1e}, loo {dKr:.1e}")
    q44 = qi[idx]
    top = np.zeros_like(dd)
    for sel in ((q44 % 3 == 1), (q44 % 3 == 2), (np.gcd(q44, 3) > 1)):
        if sel.any():
            d_, _, _ = b188.seri_ve_G(g_a, mid_a, bid_a, np.log(qa[idx][sel]),
                                      aqa[idx][sel])
            top += d_
    dlin = float(np.max(np.abs(top - dd)))
    muhur["dilim44_dogrusallik_maks"] = dlin
    log(f"  DOĞRUSALLIK Σ_r seri_r − seri (dilim 44, mod 3) maks|Δ| = {dlin:.1e} "
        f"(maks|seri| = {np.max(np.abs(dd)):.3f})")
    json.dump(muhur, open(S193 / "muhur_193.json", "w"), indent=1, ensure_ascii=False)
    del dd, Gre, Gim, R, re1, im1, top, HK

    # ================= A/B: blok-yerel sınıf ölçümü =================
    out = {"sha_onkayit": ONK["sha256"], "muhur": muhur, "hedefler": {}}
    for ad, h in HEDEFLER.items():
        t0 = time.time()
        yol = S193 / f"proj_{ad.replace('+', 'p').replace('-', 'm')}.npz"
        if yol.exists():
            D = np.load(yol)
            RE, IM = D["re"], D["im"]
            ncizgi = D["ncizgi"]
            siniflar_ad = [str(x) for x in D["siniflar"]]
            tau_aralik = (float(D["tau_min"]), float(D["tau_max"]))
            log(f"  [{ad}] kontrol noktası var")
        else:
            RE, IM, ncizgi, siniflar_ad, tau_aralik = blok_yerel_D(
                ad, h, L_B, qi, qa, taua, aqa, kenar_jk_a, g_a, mid_a, bid_a)
            np.savez_compressed(yol, re=RE, im=IM, ncizgi=ncizgi,
                                siniflar=np.array(siniflar_ad),
                                tau_min=tau_aralik[0], tau_max=tau_aralik[1])
            log(f"  [{ad}] blok-yerel seri+izdüşüm bitti ({time.time()-t0:.0f}s), "
                f"çizgi/sınıf={dict(zip(siniflar_ad, ncizgi.tolist()))}")

        R_list = h["siniflar"]
        Kt = K_hav(RE, IM, nb_a, N_a, -1)[iH]
        Kr = np.array([K_hav(RE, IM, nb_a, N_a, j)[iH] for j in range(NJACK)])
        kappa_t = -Kt.real
        kappa_r = -Kr.real
        nR = len(R_list)
        kappa_top_t = float(kappa_t[:nR].sum())
        kappa_top_r = kappa_r[:, :nR].sum(1)
        se_kappa_top = jk(kappa_top_r)

        siniflar_out = {}
        for k_, r in enumerate(R_list):
            s_t = float(kappa_t[k_] / kappa_top_t)
            s_r_reps = kappa_r[:, k_] / kappa_top_r
            se_s = jk(s_r_reps)
            z = float(s_t / se_s) if se_s > 0 else float("nan")
            siniflar_out[f"r{r}"] = {
                "n_cizgi": int(ncizgi[k_]), "kappa": float(kappa_t[k_]),
                "se_kappa": jk(kappa_r[:, k_]), "s": s_t, "se_s": float(se_s),
                "z": z, "s_pred": h["s_pred"][str(r)]}
        kb = nR
        bol_top = kappa_top_t + kappa_t[kb]
        bol_reps_top = kappa_top_r + kappa_r[:, kb]
        oran_bol = float(kappa_t[kb] / bol_top) if bol_top else float("nan")
        oran_bol_r = kappa_r[:, kb] / bol_reps_top
        bolunen_out = {"n_cizgi": int(ncizgi[kb]), "kappa": float(kappa_t[kb]),
                       "se_kappa": jk(kappa_r[:, kb]),
                       "oran_dahil_toplam": oran_bol,
                       "se_oran_dahil_toplam": jk(oran_bol_r)}

        out["hedefler"][ad] = {
            "a": h["a"], "n": str(h["n"]), "tau_aralik": tau_aralik,
            "kappa_top": kappa_top_t, "se_kappa_top": se_kappa_top,
            "siniflar": siniflar_out, "bolunen": bolunen_out}
        log(f"\n  [{ad}] τ'∈{tau_aralik}: κ_top={kappa_top_t:.5f}±{se_kappa_top:.5f} "
            f"(4σ eşik={4*se_kappa_top:.5f})")
        for et_, c in siniflar_out.items():
            log(f"     sınıf {et_:>3} ({c['n_cizgi']:6d} çizgi): κ={c['kappa']:+.5f}"
                f"±{c['se_kappa']:.5f}  s={c['s']:+.4f}±{c['se_s']:.4f}  "
                f"z={c['z']:+.2f}  öngörü={c['s_pred']:+.4f}")
        log(f"     bölünen ({bolunen_out['n_cizgi']:6d} çizgi): κ={bolunen_out['kappa']:+.5f}"
            f"±{bolunen_out['se_kappa']:.5f}  (payda-dışı; dahil-toplam oranı "
            f"{bolunen_out['oran_dahil_toplam']:.4f}±{bolunen_out['se_oran_dahil_toplam']:.4f})")

    # ================= C: +log10 r=1 TAM-örneklem kontrolü (blok-yerel) =====
    yolT = S193 / "tam_plog10_r1.npz"
    h10 = HEDEFLER["+log10"]
    a10 = h10["a"]
    log_n10 = h10["log_n"]
    if not yolT.exists():
        t0 = time.time()
        RE_T = np.zeros((NJACK, 1, len(w_win)))
        IM_T = np.zeros((NJACK, 1, len(w_win)))
        ncizgi_T = 0
        for b in range(NJACK):
            Lb = L_B[b]
            dW = np.log(qi) - Lb
            sel = np.abs(dW - log_n10) < DELTA
            qb = qi[sel]
            wl_all = np.log(qa[sel])
            al_all = aqa[sel]
            r1_mask = (qb % a10 == 1)
            ncizgi_T += int(r1_mask.sum())
            lo, hi = kenar_jk[b], kenar_jk[b + 1]
            g_b, mid_b, bid_b = g[lo:hi], mid[lo:hi], bid_tam[lo:hi]
            assert np.all(bid_b == b)
            if r1_mask.any():
                d_, _, _ = b188.seri_ve_G(g_b, mid_b, bid_b, wl_all[r1_mask],
                                          al_all[r1_mask])
            else:
                d_ = np.zeros(len(mid_b))
            re_b, im_b = b190.izdusum_blok(d_[:, None], mid_b, w_win)
            RE_T[b] = re_b
            IM_T[b] = im_b
        nb_T = np.diff(kenar_jk)
        np.savez_compressed(yolT, re=RE_T, im=IM_T, nb=nb_T, N=N, ncizgi=ncizgi_T)
        log(f"  [TAM +log10 r1] ({ncizgi_T} çizgi, {time.time()-t0:.0f}s)")
    TT = np.load(yolT)
    KT = K_hav(TT["re"], TT["im"], TT["nb"], int(TT["N"]), -1)[iH, 0]
    KTr = np.array([K_hav(TT["re"], TT["im"], TT["nb"], int(TT["N"]), j)[iH, 0]
                    for j in range(NJACK)])
    kappa_T = -KT.real
    kappa_Tr = -KTr.real
    se_kappa_T = jk(kappa_Tr)
    k_alt = out["hedefler"]["+log10"]["siniflar"]["r1"]["kappa"]
    se_alt = out["hedefler"]["+log10"]["siniflar"]["r1"]["se_kappa"]
    fark = abs(kappa_T - k_alt)
    out["tam_kontrol_plog10_r1"] = {
        "n_cizgi_tam": int(TT["ncizgi"]), "kappa_alt": k_alt, "se_alt": se_alt,
        "kappa_tam": float(kappa_T), "se_tam": float(se_kappa_T),
        "fark": float(fark), "fark_bolu_se_alt": float(fark / se_alt) if se_alt else None}
    log(f"\n  TAM-ÖRNEKLEM +log10 r1: κ_alt={k_alt:.5f}±{se_alt:.5f} "
        f"κ_tam={kappa_T:.5f}±{se_kappa_T:.5f} |Δ|={fark:.2e} "
        f"(={fark/se_alt:.2f} se_alt)")

    json.dump(out, open(S193 / "K1_sinif.json", "w"), indent=1, ensure_ascii=False)
    log(f"\n-> K1_sinif.json  BİTTİ ({time.time()-T0:.0f}s)")
