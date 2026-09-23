# -*- coding: utf-8 -*-
"""
188b — K1: İPTAL ÇEKİRDEĞİ HARİTASI K(τ_b, τ') (gerçek kinematik; ANA KOŞU)
===========================================================================
ONKAYIT_188 dondurdu:
  dilim s = τ' ∈ (e_s, e_{s+1}], 88 dilim × 0.005, (0.86, 1.30]; TÜM
  asal-kuvvetler, a_q = Λ(q)/(π√q·log q) (187b.asal_kuvvetler AYNEN).
  dds_s(n) = Σ_{q'∈s} 2a' sin(ω' g_n/2) cos(ω' m_n)  (187b katman_serisi
  AYNEN) — ALT-ÖRNEKLEM n = 0,3,6,… (100k).
  c^{(s)}_q = 2⟨dds_s e^{−iω_q m}⟩_alt (jk-blok toplamları, loo).
  K(b,s) = Σ_b c^{(s)} conj(karışım)/Σ_b |karışım|², karışım = c^kesik − c^öz
  (187c AYNEN, tam örneklem).
  Ĝ_mid(ω') = mean_n e^{+iω' m_n} (MUTLAK m_n) aynı alt-örneklemde, jk-blok
  toplamlarıyla; S(s) = Σ a'|Ĝ|, W(s) = Σ a'.
Aşamalar: A dilim serileri (kontrol noktalı dilim_<i>.npz; çok süreçli);
B tek geçişte izdüşüm; C iki dilimde TAM-örneklem kontrolü; D K + S/W +
K1 tutarlılık kontrolü.

Çıktı: 188/dilim_<i>.npz, harita_proj_gercek.npz, kontrol_tam.npz,
       harita_K_gercek.npz, K1_harita.json
Kullanım: 188b_harita.py [işçi_sayısı]
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
S184, S185, S186, S187, S188 = (SCR / d for d in
                                ["184", "185", "186", "187", "188"])
NJACK = 8
PBLOK = 2048      # 187b AYNEN
LBLOK = 4096      # 187b AYNEN
PROJ_BLOK = 1500  # 186b/187b AYNEN
ALT = 3


def _yukle(ad, yol):
    spec = importlib.util.spec_from_file_location(ad, yol)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


b185 = _yukle("b185", QM / "185_configs" / "185b_oz_muhasebe.py")
b187 = _yukle("b187", QM / "187_configs" / "187b_katman_defteri.py")


def onkayit():
    o = json.load(open(S188 / "ONKAYIT_188.json"))
    sha = hashlib.sha256(
        (QM / "188_configs" / "188a_onkayit.py").read_bytes()).hexdigest()
    if sha != o["sha256"]:
        raise SystemExit(f"ON-KAYIT SHA UYUMSUZ: {sha} != {o['sha256']}")
    return o


def jk_kenar(N):
    return np.linspace(0, N, NJACK + 1).astype(int)


def seri_ve_G(g, mid, bid, w_l, a_l):
    """187b.katman_serisi AYNEN (aynı bloklama/işlem sırası → bit-bit aynı dds)
    + Ĝ_mid jk-blok toplamları (Σ cos(ω'm), Σ sin(ω'm))."""
    N = len(mid)
    dds = np.zeros(N)
    kats = 2.0 * a_l
    nl = len(w_l)
    Gre = np.zeros((NJACK, nl))
    Gim = np.zeros((NJACK, nl))
    for s0 in range(0, N, PBLOK):
        sl = slice(s0, min(s0 + PBLOK, N))
        acc = np.zeros(sl.stop - sl.start)
        bb = bid[sl]
        parca = [(b, np.where(bb == b)[0]) for b in np.unique(bb)]
        for l0 in range(0, nl, LBLOK):
            ll = slice(l0, min(l0 + LBLOK, nl))
            P = np.outer(mid[sl], w_l[ll])
            Sp = np.sin(P)
            for b, r in parca:
                Gim[b, ll] += Sp[r[0]:r[-1] + 1].sum(0)
            del Sp
            np.cos(P, out=P)
            for b, r in parca:
                Gre[b, ll] += P[r[0]:r[-1] + 1].sum(0)
            K = np.sin(0.5 * np.outer(g[sl], w_l[ll]))
            K *= P
            acc += K @ kats[ll]
            del P, K
        dds[sl] = acc
    return dds, Gre, Gim


# ---------------- işçi (çok süreçli) ----------------
_W = {}


def _init(g, mid, bid):
    _W["g"], _W["mid"], _W["bid"] = g, mid, bid


def _is(arg):
    i, yol, q, w_l, a_l, tau_l = arg
    t0 = time.time()
    dds, Gre, Gim = seri_ve_G(_W["g"], _W["mid"], _W["bid"], w_l, a_l)
    tmp = str(yol) + ".tmp.npz"
    np.savez(tmp, dds=dds, Gre=Gre, Gim=Gim, q=q, a=a_l, tau=tau_l,
             n=len(_W["mid"]))
    os.replace(tmp, yol)
    return i, len(w_l), time.time() - t0


def dilimleri_kur(etiket, g, mid, bid, kenar, L, qa, taua, aqa, onek, nw):
    """Dilim serilerini (kontrol noktalı) çok süreçli kur."""
    import multiprocessing as mp
    ns = len(kenar) - 1
    isler = []
    for i in range(ns):
        yol = S188 / f"{onek}{i}.npz"
        if yol.exists():
            continue
        idx = np.where((taua > kenar[i]) & (taua <= kenar[i + 1]))[0]
        isler.append((i, yol, qa[idx], np.log(qa[idx]), aqa[idx], taua[idx]))
    print(f"  [{etiket}] {ns} dilim; eksik {len(isler)} (işçi {nw})",
          flush=True)
    if not isler:
        return
    isler.sort(key=lambda x: -len(x[2]))
    t0 = time.time()
    ctx = mp.get_context("fork")
    with ctx.Pool(nw, initializer=_init, initargs=(g, mid, bid)) as pool:
        k = 0
        for i, nl, dt in pool.imap_unordered(_is, isler):
            k += 1
            print(f"    [{etiket}] dilim {i:2d} bitti ({nl} çizgi, {dt:.0f}s)"
                  f"  {k}/{len(isler)}  toplam {time.time()-t0:.0f}s",
                  flush=True)


def izdusum(D, mid, bid, w_win):
    """D: (N, ns) dilim serileri. jk-blok toplamları re/im (8, ns, nwin)."""
    ns = D.shape[1]
    re = np.zeros((NJACK, ns, len(w_win)))
    im = np.zeros((NJACK, ns, len(w_win)))
    N = len(mid)
    for b in range(NJACK):
        rr = np.where(bid == b)[0]
        lo, hi = rr[0], rr[-1] + 1
        for s0 in range(lo, hi, PROJ_BLOK):
            sl = slice(s0, min(s0 + PROJ_BLOK, hi))
            P = np.outer(mid[sl], w_win)
            cP = np.cos(P)
            sP = np.sin(P)
            del P
            re[b] += D[sl].T @ cP
            im[b] -= D[sl].T @ sP
            del cP, sP
    nb = np.array([(bid == b).sum() for b in range(NJACK)])
    return re, im, nb, N


def c_proj(re, im, nb, N, disari=-1):
    """(ns, nwin) karmaşık izdüşüm, loo destekli."""
    if disari < 0:
        return 2.0 * (re.sum(0) + 1j * im.sum(0)) / N
    return 2.0 * ((re.sum(0) - re[disari]) + 1j * (im.sum(0) - im[disari])) \
        / (N - nb[disari])


def karisim(P, OZ, aq, disari=-1):
    """187c AYNEN: karışım = c^kesik − c^öz (tam örneklem, loo)."""
    return b187.c_kesik(P, disari) - b185.c_oz(OZ, aq, disari)


def bant_maskeleri(tau_win, ONK):
    ince = ONK["ince_bant_kenar"]
    b8 = ONK["bant8_kenar"]
    m = [(tau_win >= ince[b]) & (tau_win < ince[b + 1])
         for b in range(len(ince) - 1)]
    m += [(tau_win >= b8[b]) & (tau_win < b8[b + 1]) for b in range(8)]
    m.append((tau_win >= b8[0]) & (tau_win < b8[-1]))
    et = [f"{ince[b]:.2f}-{ince[b+1]:.2f}" for b in range(len(ince) - 1)]
    et += [f"B{b8[b]:.2f}-{b8[b+1]:.2f}" for b in range(8)] + ["HAVUZ"]
    return m, et


def K_matris(C, mix, maskeler):
    """K(b,s) = Σ_b C[s,q] conj(mix_q) / Σ_b |mix_q|²  → (nbant, ns)."""
    out = np.zeros((len(maskeler), C.shape[0]), complex)
    pay = np.zeros(len(maskeler))
    for k, m in enumerate(maskeler):
        pay[k] = np.sum(np.abs(mix[m]) ** 2)
        out[k] = (C[:, m] @ np.conj(mix[m])) / pay[k]
    return out, pay


def G_loo(Gre, Gim, nb, N, disari=-1):
    if disari < 0:
        return (Gre.sum(0) + 1j * Gim.sum(0)) / N
    return ((Gre.sum(0) - Gre[disari]) + 1j * (Gim.sum(0) - Gim[disari])) \
        / (N - nb[disari])


def kos(gaz, kenar_anahtar, onek, Kad, OZad, Pad, cikti_on, nw,
        tam_kontrol=None):
    ONK = onkayit()
    L = ONK["L"]
    kenar = np.array(ONK[kenar_anahtar])
    ns = len(kenar) - 1
    G = np.load(S184 / "K1_gercek.npz")   # çizgi evreni + maskeler (187c AYNEN)
    Gk = np.load(S184 / Kad)              # yalnız jk-blok kontrolü
    OZ = np.load(S185 / OZad)
    P = np.load(S186 / Pad)
    w_all = np.asarray(G["w"], float)
    tau_all = np.asarray(G["tau"], float)
    aq_all = np.asarray(G["aq"], float)
    win = (tau_all >= 0.45) & (tau_all < 0.86)
    w_win, tau_win = w_all[win], tau_all[win]

    g, mid, Lg = b185.kinematik(gaz)
    N = len(mid)
    kenar_jk = jk_kenar(N)
    if not (np.array_equal(np.diff(kenar_jk), G["nb"]) and
            np.array_equal(np.diff(kenar_jk), Gk["nb"])):
        raise SystemExit("JK BLOK KENARLARI 184 İLE UYUMSUZ")
    bid_tam = np.searchsorted(kenar_jk, np.arange(N), side="right") - 1
    alt = np.arange(0, N, ALT)
    g_a, mid_a, bid_a = g[alt], mid[alt], bid_tam[alt]
    print(f"[{gaz}] N={N} alt={len(alt)} L_kin={Lg:.9f} L_onk={L:.9f} "
          f"pencere çizgi={win.sum()} dilim={ns}", flush=True)

    qmax = np.exp(kenar[-1] * L)
    t0 = time.time()
    qa, lam, taua, aqa = b187.asal_kuvvetler(qmax, L)
    print(f"  asal-kuvvetler q≤{qmax:.0f}: {len(qa)} ({time.time()-t0:.0f}s);"
          f" pencere-ötesi (τ'∈({kenar[0]},{kenar[-1]}]): "
          f"{int(((taua > kenar[0]) & (taua <= kenar[-1])).sum())}", flush=True)

    # --- MAKİNE MÜHRÜ: seri_ve_G dds ≡ 187b.katman_serisi (bit-bit) ---
    idx0 = np.where((taua > kenar[0]) & (taua <= kenar[1]))[0]
    d_ref = b187.katman_serisi(g_a, mid_a, np.log(qa[idx0]), aqa[idx0],
                               "mühür")
    d_yeni, _, _ = seri_ve_G(g_a, mid_a, bid_a, np.log(qa[idx0]), aqa[idx0])
    muhur_bit = float(np.max(np.abs(d_ref - d_yeni)))
    print(f"  MAKİNE MÜHRÜ: maks|dds_188 − dds_187b| (dilim 0) = "
          f"{muhur_bit:.1e}", flush=True)

    # --- A: dilim serileri ---
    dilimleri_kur(gaz, g_a, mid_a, bid_a, kenar, L, qa, taua, aqa, onek, nw)

    # --- B: izdüşüm (tek geçiş) ---
    yolB = S188 / f"harita_proj_{cikti_on}.npz"
    if not yolB.exists():
        t0 = time.time()
        D = np.stack([np.load(S188 / f"{onek}{i}.npz")["dds"]
                      for i in range(ns)], axis=1)
        re, im, nb_a, N_a = izdusum(D, mid_a, bid_a, w_win)
        np.savez_compressed(yolB, re=re, im=im, nb=nb_a, N=N_a, kenar=kenar)
        print(f"  [izdüşüm] -> {yolB.name} ({time.time()-t0:.0f}s)",
              flush=True)
    PR = np.load(yolB)
    re, im, nb_a, N_a = PR["re"], PR["im"], PR["nb"], int(PR["N"])

    # dilim çizgi özellikleri (Ĝ, S, W)
    Gre_s, Gim_s, a_s, tau_s, ncz = [], [], [], [], []
    for i in range(ns):
        d = np.load(S188 / f"{onek}{i}.npz")
        Gre_s.append(d["Gre"])
        Gim_s.append(d["Gim"])
        a_s.append(d["a"])
        tau_s.append(d["tau"])
        ncz.append(len(d["a"]))

    def SW(disari):
        S = np.zeros(ns)
        SRe = np.zeros(ns)
        for i in range(ns):
            Gh = G_loo(Gre_s[i], Gim_s[i], nb_a, N_a, disari)
            S[i] = np.sum(a_s[i] * np.abs(Gh))
            SRe[i] = np.sum(a_s[i] * np.pi * tau_s[i] *
                            np.cos(np.pi * tau_s[i]) * Gh.real)
        return S, SRe

    W = np.array([a.sum() for a in a_s])
    S_tam, SRe_tam = SW(-1)
    S_reps = np.array([SW(j)[0] for j in range(NJACK)])
    SRe_reps = np.array([SW(j)[1] for j in range(NJACK)])

    # --- D: K matrisi (tam + loo) ---
    maskeler, et = bant_maskeleri(tau_win, ONK)
    aq_win = aq_all[win]

    def K_hesap(disari):
        C = c_proj(re, im, nb_a, N_a, disari)
        mix = karisim(P, OZ, aq_all, disari)[win]
        return K_matris(C, mix, maskeler)

    K_tam, pay_tam = K_hesap(-1)
    rep = [K_hesap(j) for j in range(NJACK)]
    K_reps = np.array([r[0] for r in rep])
    pay_reps = np.array([r[1] for r in rep])

    sonuc = {"sha_onkayit": ONK["sha256"], "gaz": gaz, "N_alt": N_a,
             "nb_alt": nb_a.tolist(), "dilim_cizgi": ncz,
             "muhur_bit_dds": muhur_bit, "etiket": et}

    # --- C: iki dilimde TAM-örneklem kontrolü ---
    if tam_kontrol:
        yolC = S188 / f"kontrol_tam_{cikti_on}.npz"
        sec = [int(np.argmin(np.abs(kenar[:-1] - lo))) for lo in tam_kontrol]
        if not yolC.exists():
            t0 = time.time()
            Dt = []
            for i in sec:
                idx = np.where((taua > kenar[i]) & (taua <= kenar[i + 1]))[0]
                d, _, _ = seri_ve_G(g, mid, bid_tam, np.log(qa[idx]),
                                    aqa[idx])
                Dt.append(d)
            Dt = np.stack(Dt, axis=1)
            ret, imt, nbt, Nt = izdusum(Dt, mid, bid_tam, w_win)
            np.savez_compressed(yolC, re=ret, im=imt, nb=nbt, N=Nt,
                                sec=np.array(sec))
            print(f"  [tam kontrol] -> {yolC.name} ({time.time()-t0:.0f}s)",
                  flush=True)
        TC = np.load(yolC)
        Kt, _ = K_matris(c_proj(TC["re"], TC["im"], TC["nb"], int(TC["N"])),
                         karisim(P, OZ, aq_all)[win], maskeler)
        Kt_reps = np.array([K_matris(
            c_proj(TC["re"], TC["im"], TC["nb"], int(TC["N"]), j),
            karisim(P, OZ, aq_all, j)[win], maskeler)[0]
            for j in range(NJACK)])
        tk = {}
        print("\n  TAM-ÖRNEKLEM KONTROLÜ (K_alt vs K_tam; karmaşık):")
        for jj, i in enumerate(sec):
            ka = K_tam[:, i]
            kt = Kt[:, jj]
            fark = np.abs(ka - kt)
            se_c = np.sqrt(b185.jk_se(K_reps[:, :, i].real) ** 2 +
                           b185.jk_se(K_reps[:, :, i].imag) ** 2)
            se_t = np.sqrt(b185.jk_se(Kt_reps[:, :, jj].real) ** 2 +
                           b185.jk_se(Kt_reps[:, :, jj].imag) ** 2)
            ince = slice(0, 41)
            kor = np.real(np.vdot(kt[ince], ka[ince])) / (
                np.linalg.norm(kt[ince]) * np.linalg.norm(ka[ince]))
            tk[f"({kenar[i]:.3f},{kenar[i+1]:.3f}]"] = {
                "K_alt_HAVUZ": [float(ka[-1].real), float(ka[-1].imag)],
                "K_tam_HAVUZ": [float(kt[-1].real), float(kt[-1].imag)],
                "se_alt_HAVUZ": float(se_c[-1]), "se_tam_HAVUZ": float(se_t[-1]),
                "fark_HAVUZ": float(fark[-1]),
                "fark_bolu_se_alt_HAVUZ": float(fark[-1] / se_c[-1]),
                "B8_fark": fark[41:49].tolist(),
                "B8_fark_bolu_se_alt": (fark[41:49] / se_c[41:49]).tolist(),
                "ince41_karmasik_korr": float(kor),
                "ince41_medyan_fark_bolu_se_alt":
                    float(np.median(fark[:41] / se_c[:41]))}
            print(f"   dilim ({kenar[i]:.3f},{kenar[i+1]:.3f}]: HAVUZ "
                  f"K_alt={ka[-1]:.5f} K_tam={kt[-1]:.5f} |Δ|={fark[-1]:.2e}"
                  f" (={fark[-1]/se_c[-1]:.2f} se_alt; se_tam={se_t[-1]:.1e});"
                  f" 41-bant karmaşık korr={kor:.3f}", flush=True)
        sonuc["tam_kontrol"] = tk

    # --- K1 tutarlılık: HAVUZ Σ_{τ'≤1.20} K ---
    orta = 0.5 * (kenar[:-1] + kenar[1:])
    m120 = kenar[1:] <= 1.20 + 1e-9
    iH = len(et) - 1

    def top(Kx, m):
        return Kx[..., m].sum(-1)

    for ad, m in [("le120", m120), ("le130", np.ones(ns, bool))]:
        z = top(K_tam[iH], m)
        zr = top(K_reps[:, iH], m)
        se_mod = float(b185.jk_se(np.abs(zr)))
        dfark = (np.degrees(np.angle(zr)) - np.degrees(np.angle(z)) + 180) \
            % 360 - 180
        se_aci = float(np.sqrt((NJACK - 1) / NJACK *
                               np.sum((dfark - dfark.mean()) ** 2)))
        sonuc[f"HAVUZ_toplam_{ad}"] = {
            "mod": float(abs(z)), "se_mod": se_mod,
            "aci": float(np.degrees(np.angle(z))), "se_aci": se_aci}
        print(f"  HAVUZ Σ_s K ({ad}): {abs(z):.4f}±{se_mod:.4f} "
              f"∠{np.degrees(np.angle(z)):+.2f}°±{se_aci:.2f}", flush=True)

    np.savez_compressed(
        S188 / f"harita_K_{cikti_on}.npz", K=K_tam, K_reps=K_reps,
        pay=pay_tam, pay_reps=pay_reps, kenar=kenar, orta=orta,
        etiket=np.array(et), S=S_tam, S_reps=S_reps, SRe=SRe_tam,
        SRe_reps=SRe_reps, W=W, ncizgi=np.array(ncz),
        ince_kenar=np.array(ONK["ince_bant_kenar"]))
    json.dump(sonuc, open(S188 / f"K1_harita_{cikti_on}.json", "w"),
              indent=1, ensure_ascii=False)
    print(f"-> harita_K_{cikti_on}.npz, K1_harita_{cikti_on}.json  BİTTİ",
          flush=True)


if __name__ == "__main__":
    nw = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    ONK = onkayit()
    print("=" * 78)
    print(f"188b / K1 İPTAL ÇEKİRDEĞİ HARİTASI (gerçek)  [on-kayit "
          f"{ONK['zaman']} sha {ONK['sha256'][:12]}]")
    print("=" * 78, flush=True)
    kos("gercek", "dilim_izgara_gercek", "dilim_", "K1_gercek.npz",
        "OZ_gercek.npz", "G1_proj_gercek.npz", "gercek", nw,
        tam_kontrol=[1.050, 0.900])
