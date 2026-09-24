# -*- coding: utf-8 -*-
"""
194b — K1: BOŞ PENCERE + μ=0 PENCERE ÖLÇÜMÜ (193b_sinif.py makinesi AYNEN)
=============================================================================
193b_sinif.py'nin blok-yerel Δω seçimi (δ=0.03), sınıf ayrıştırması, HAVUZ κ
ve 8-blok jackknife makinesi (blok_yerel_D, jk, K_hav düzeni — 188b/190b/187b
AYNEN çağrı zinciri) importlib ile YÜKLENİR (dosya DÜZENLENMEZ) ve AYNEN
çağrılır. YALNIZ pencere merkezleri (ve bazı pencerelerde sınıf modülü)
değişti — ONKAYIT_194'ten (194a) okunur.

 M   MAKİNE MÜHÜRÜ: +log10 hedefini (193'ün KENDİ HEDEFLER tanımıyla, a193'ten
     AYNEN) bu betiğin kendi çağrı zinciriyle yeniden ölçüp κ_top'u 193'ün
     ARŞİVLENMİŞ K1_sinif.json değeriyle bit-bit karşılaştırır.
 A   Her BOŞ pencere için: birincil κ (a=1, tek sınıf — HAVUZ, sınıfsız) +
     mod 3, mod 5, mod 10 sınıf payları (H-194c için).
 B   Her μ=0 penceresi (Δω aralığı İÇİNDEKİ 8 tanesi) için: birincil κ (a=1).
 Her ölçüm 8-blok jackknife (loo) replikalarıyla birlikte saklanır (194c'nin
 havuzlanmış/pooled jackknife ile κ̄_boş, κ̄_μ0, eğim vb. hesaplaması için).

Çıktı: scratchpad/194/proj_<pencere>_<mod>.npz (kontrol noktası, resumable),
       scratchpad/194/muhur_194.json, scratchpad/194/K1_bos_194.json
Koşu: nohup (asal-kuvvetler + ~25 blok-yerel seri çağrısı; birkaç dakika
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
S184, S185, S186, S188, S193, S194 = (SCR / d for d in
                                      ["184", "185", "186", "188", "193", "194"])
NJACK = 8
DELTA = 0.03
TWO_PI = 2 * np.pi


def yukle(ad, yol):
    spec = importlib.util.spec_from_file_location(ad, yol)
    m = importlib.util.module_from_spec(spec)
    sys.modules[ad] = m
    spec.loader.exec_module(m)
    return m


def log(*a):
    print(*a, flush=True)


def onkayit():
    o = json.load(open(S194 / "ONKAYIT_194.json"))
    s = hashlib.sha256((QM / "194_configs" / "194a_onkayit.py").read_bytes()).hexdigest()
    if s != o["sha256"]:
        raise SystemExit(f"ON-KAYIT SHA UYUMSUZ: {s} != {o['sha256']}")
    return o


if __name__ == "__main__":
    T0 = time.time()
    ONK = onkayit()
    log("=" * 78)
    log(f"194b / K1 BOŞ+MU0 PENCERE ÖLÇÜMÜ  [on-kayit {ONK['zaman']} sha "
        f"{ONK['sha256'][:12]}]")
    log("=" * 78)

    # ---- 193b_sinif.py'yi AYNEN yükle (blok_yerel_D, jk, b188/b190/b187) ----
    b193 = yukle("b193", QM / "193_configs" / "193b_sinif.py")
    a193 = yukle("a193", QM / "193_configs" / "193a_onkayit.py")  # +log10 hedef tanımı (mühür)
    b188, b190, b187 = b193.b188, b193.b190, b193.b187
    jk = b193.jk
    blok_yerel_D = b193.blok_yerel_D
    assert b193.DELTA == DELTA == 0.03, "DELTA AYNEN olmalı"

    ONK188 = json.load(open(S188 / "ONKAYIT_188.json"))
    L = float(ONK188["L"])

    # ---- kinematik + pencere (193b AYNEN) ----
    G = np.load(S184 / "K1_gercek.npz")
    OZ = np.load(S185 / "OZ_gercek.npz")
    P = np.load(S186 / "G1_proj_gercek.npz")
    w_all = np.asarray(G["w"], float)
    tau_all = np.asarray(G["tau"], float)
    aq_all = np.asarray(G["aq"], float)
    win = (tau_all >= 0.45) & (tau_all < 0.86)
    w_win, tau_win = w_all[win], tau_all[win]
    b193.w_win_G[0] = w_win

    g, mid, Lg = b188.b185.kinematik("gercek")
    N = len(mid)
    kenar_jk = b188.jk_kenar(N)
    assert np.array_equal(np.diff(kenar_jk), G["nb"])
    bid_tam = np.searchsorted(kenar_jk, np.arange(N), side="right") - 1
    alt = np.arange(0, N, b188.ALT)
    g_a, mid_a, bid_a = g[alt], mid[alt], bid_tam[alt]
    assert np.all(np.diff(bid_a) >= 0), "bid_a sıralı değil"
    kenar_jk_a = np.searchsorted(bid_a, np.arange(NJACK + 1))
    assert kenar_jk_a[0] == 0 and kenar_jk_a[-1] == len(alt)
    maskeler, et = b188.bant_maskeleri(tau_win, ONK188)
    iH = len(et) - 1
    assert et[iH] == "HAVUZ"
    nb_a = np.diff(kenar_jk_a)
    N_a = len(alt)
    log(f"  N={N} alt={len(alt)} nb_a={nb_a.tolist()} L={L:.9f}")

    L_B = np.array([float(np.log(mid[kenar_jk[b]:kenar_jk[b + 1]] / TWO_PI).mean())
                    for b in range(NJACK)])
    log(f"  L_B={np.round(L_B,5).tolist()}")

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

    # ================= M: MAKİNE MÜHÜRÜ (+log10, 193 ile bit-bit) =================
    h10 = a193.HEDEFLER["+log10"]
    t0 = time.time()
    RE10, IM10, ncizgi10, siniflar10, tau10 = blok_yerel_D(
        "+log10(muhur)", h10, L_B, qi, qa, taua, aqa, kenar_jk_a, g_a, mid_a, bid_a)
    Kt10 = K_hav(RE10, IM10, nb_a, N_a, -1)[iH]
    kappa_t10 = -Kt10.real
    nR10 = len(h10["siniflar"])
    kappa_top_194 = float(kappa_t10[:nR10].sum())
    K193 = json.load(open(S193 / "K1_sinif.json"))
    kappa_top_193 = float(K193["hedefler"]["+log10"]["kappa_top"])
    fark_muhur = abs(kappa_top_194 - kappa_top_193)
    muhur = {"kappa_top_194": kappa_top_194, "kappa_top_193": kappa_top_193,
             "fark_mutlak": fark_muhur, "sure_s": time.time() - t0,
             "ncizgi": {k: int(v) for k, v in zip(siniflar10, ncizgi10.tolist())}}
    log(f"\n  MÜHÜR +log10 κ_top: 194={kappa_top_194:.10f}  193={kappa_top_193:.10f}  "
        f"|Δ|={fark_muhur:.3e}  ({muhur['sure_s']:.0f}s)")
    json.dump(muhur, open(S194 / "muhur_194.json", "w"), indent=1, ensure_ascii=False)
    del RE10, IM10

    # ================= yardımcı: bir pencereyi ölç =================
    def olc(ad_dosya, log_n, a, siniflar):
        yol = S194 / f"proj_{ad_dosya}.npz"
        if yol.exists():
            D = np.load(yol)
            RE, IM = D["re"], D["im"]
            ncizgi = D["ncizgi"]
            siniflar_ad = [str(x) for x in D["siniflar"]]
            tau_aralik = (float(D["tau_min"]), float(D["tau_max"]))
            log(f"    [{ad_dosya}] kontrol noktası var")
        else:
            t0 = time.time()
            h = {"a": a, "siniflar": siniflar, "log_n": log_n}
            RE, IM, ncizgi, siniflar_ad, tau_aralik = blok_yerel_D(
                ad_dosya, h, L_B, qi, qa, taua, aqa, kenar_jk_a, g_a, mid_a, bid_a)
            np.savez_compressed(yol, re=RE, im=IM, ncizgi=ncizgi,
                                siniflar=np.array(siniflar_ad),
                                tau_min=tau_aralik[0], tau_max=tau_aralik[1])
            log(f"    [{ad_dosya}] bitti ({time.time()-t0:.0f}s) "
                f"çizgi/sınıf={dict(zip(siniflar_ad, ncizgi.tolist()))}")
        Kt = K_hav(RE, IM, nb_a, N_a, -1)[iH]
        Kr = np.array([K_hav(RE, IM, nb_a, N_a, j)[iH] for j in range(NJACK)])
        kappa_t = -Kt.real
        kappa_r = -Kr.real
        nR = len(siniflar)
        kappa_top_t = float(kappa_t[:nR].sum())
        kappa_top_r = kappa_r[:, :nR].sum(1)
        out = {"tau_aralik": tau_aralik, "kappa_top": kappa_top_t,
               "kappa_top_se": jk(kappa_top_r),
               "kappa_top_reps": kappa_top_r.tolist(),
               "siniflar": {}}
        for k_, r in enumerate(siniflar):
            s_t = float(kappa_t[k_] / kappa_top_t) if kappa_top_t else float("nan")
            s_reps = kappa_r[:, k_] / kappa_top_r
            out["siniflar"][f"r{r}"] = {
                "n_cizgi": int(ncizgi[k_]), "kappa": float(kappa_t[k_]),
                "se_kappa": jk(kappa_r[:, k_]), "s": s_t, "se_s": jk(s_reps),
                "kappa_reps": kappa_r[:, k_].tolist(), "s_reps": s_reps.tolist()}
        kb = nR
        out["bolunen"] = {"n_cizgi": int(ncizgi[kb]), "kappa": float(kappa_t[kb]),
                          "se_kappa": jk(kappa_r[:, kb])}
        del RE, IM
        return out

    OUT = {"sha_onkayit": ONK["sha256"], "muhur": muhur,
           "bos_pencereler": {}, "mu0_pencereler": {}}

    # ================= A: BOŞ PENCERELER (birincil a=1 + mod 3,5,10) =================
    for bp in ONK["bos_pencereler"]:
        ad = bp["ad"]
        merkez = bp["merkez"]
        log(f"\n  === BOŞ {ad} (merkez={merkez:.5f}) ===")
        primer = olc(f"bos_{ad}_a1", merkez, 1, [0])
        mod3 = olc(f"bos_{ad}_mod3", merkez, 3, [1, 2])
        mod5 = olc(f"bos_{ad}_mod5", merkez, 5, [1, 2, 3, 4])
        mod10 = olc(f"bos_{ad}_mod10", merkez, 10, [1, 3, 7, 9])
        OUT["bos_pencereler"][ad] = {
            "merkez": merkez, "komsu_sol": bp["komsu_sol"], "komsu_sag": bp["komsu_sag"],
            "kappa": primer["kappa_top"], "se_kappa": primer["kappa_top_se"],
            "kappa_reps": primer["kappa_top_reps"], "n_cizgi": primer["siniflar"]["r0"]["n_cizgi"],
            "mod3": mod3, "mod5": mod5, "mod10": mod10}
        log(f"    κ={primer['kappa_top']:+.5f}±{primer['kappa_top_se']:.5f}  "
            f"n_cizgi={primer['siniflar']['r0']['n_cizgi']}")

    # ================= B: MU=0 PENCERELERİ (yalnız aralık-içi; birincil a=1) =====
    for mp in ONK["mu0_pencereler"]:
        if mp["durum"] != "aralik_ici":
            continue
        ad = mp["ad"]
        merkez = mp["merkez"]
        ad_dosya = "mu0_" + ad.replace("+", "p").replace("(", "").replace(")", "").replace("/", "_")
        log(f"\n  === MU0 {ad} (merkez={merkez:.5f}) ===")
        primer = olc(f"{ad_dosya}_a1", merkez, 1, [0])
        OUT["mu0_pencereler"][ad] = {
            "merkez": merkez, "a": mp["a"], "b": mp["b"],
            "kappa": primer["kappa_top"], "se_kappa": primer["kappa_top_se"],
            "kappa_reps": primer["kappa_top_reps"],
            "n_cizgi": primer["siniflar"]["r0"]["n_cizgi"]}
        log(f"    κ={primer['kappa_top']:+.5f}±{primer['kappa_top_se']:.5f}  "
            f"n_cizgi={primer['siniflar']['r0']['n_cizgi']}")

    json.dump(OUT, open(S194 / "K1_bos_194.json", "w"), indent=1, ensure_ascii=False)
    log(f"\n-> K1_bos_194.json  BİTTİ ({time.time()-T0:.0f}s)")
