# -*- coding: utf-8 -*-
"""
187b — K1: KATMAN DEFTERİ (gerçek kinematik; ANA KOŞU)
======================================================
ONKAYIT_187 dondurdu:
  katman i = τ' ∈ (ızgara[i-1], ızgara[i]] asal-kuvvetleri; ızgara GERÇEK
  [0.86, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20]; a_q = Λ(q)/(π√q·log q).
  dds_i(n) = Σ_{q∈i} 2a_q sin(ω_q g_n/2) cos(ω_q m_n)   (ρ≡1, nominal)
  c^{Δi}_q = 2⟨dds_i e^{−iω_q m}⟩ ;  c^ya(τ_c) = c^kesik + Σ_{i≤τ_c} c^{Δi}
  m^ya_b(τ_c) = [Σ_b |c^ya| − Σ_b â^öz]/Σ_b ae ;  hedef m_ölç,b
  c^kesik = 186 G1_proj_gercek re_bir/im_bir AYNEN (yeniden hesap YOK).
Aşamalar: A katman serileri (kontrol noktalı, katman_<i>_gercek.npy);
B tek geçişte tüm serilerin bant-izdüşümü + taban momidleri (makine mührü);
C defter + jackknife + H-F1 çatal hükmü.

Çıktı: 187/katman_<i>_gercek.npy, katmanproj_gercek.npz, K1_katman_defteri.json
"""
import hashlib
import importlib.util
import json
import time
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S155, S184, S185, S186, S187 = (SCR / d for d in
                                ["155", "184", "185", "186", "187"])
TWO_PI = 2 * np.pi
NJACK = 8
PBLOK = 2048     # nokta bloğu (katman inşası)
LBLOK = 4096     # çizgi bloğu (katman inşası)
PROJ_BLOK = 1500  # izdüşüm nokta bloğu (186b AYNEN)

spec = importlib.util.spec_from_file_location(
    "b185", QM / "185_configs" / "185b_oz_muhasebe.py")
b185 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(b185)


def onkayit():
    o = json.load(open(S187 / "ONKAYIT_187.json"))
    sha = hashlib.sha256(
        (QM / "187_configs" / "187a_onkayit.py").read_bytes()).hexdigest()
    if sha != o["sha256"]:
        raise SystemExit(f"ON-KAYIT SHA UYUMSUZ: {sha} != {o['sha256']}")
    return o


def asal_kuvvetler(qmax, L):
    """Tüm asal-kuvvetler q ≤ qmax: (q, Λ, τ'=log q/L, a_q) — sympy sieve."""
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
    return qs, lams, lq / L, lams / (np.pi * np.sqrt(qs) * lq)


def katman_serisi(g, mid, w_l, a_l, etiket):
    """dds_i(n) = Σ 2a sin(ωg/2)cos(ωm) — çift-bloklu vektörleme."""
    t0 = time.time()
    N = len(mid)
    dds = np.zeros(N)
    kats = 2.0 * a_l
    for s0 in range(0, N, PBLOK):
        sl = slice(s0, min(s0 + PBLOK, N))
        acc = np.zeros(sl.stop - sl.start)
        for l0 in range(0, len(w_l), LBLOK):
            ll = slice(l0, min(l0 + LBLOK, len(w_l)))
            P = np.outer(mid[sl], w_l[ll])
            np.cos(P, out=P)
            K = np.sin(0.5 * np.outer(g[sl], w_l[ll]))
            K *= P
            acc += K @ kats[ll]
            del P, K
        dds[sl] = acc
        if s0 % (PBLOK * 32) == 0:
            print(f"    [{etiket}] nokta {s0}/{N}  ({time.time()-t0:.0f}s)",
                  flush=True)
    print(f"  [{etiket}] seri bitti: {len(w_l)} çizgi, "
          f"({time.time()-t0:.0f}s)", flush=True)
    return dds


def c_kesik(P, disari=-1):
    """186 G1_proj re_bir/im_bir bloklarından kompleks kesik izdüşüm (loo)."""
    N = int(P["N"])
    nb = P["nb"]
    if disari < 0:
        re, im, n = P["re_bir"].sum(0), P["im_bir"].sum(0), N
    else:
        re = P["re_bir"].sum(0) - P["re_bir"][disari]
        im = P["im_bir"].sum(0) - P["im_bir"][disari]
        n = N - nb[disari]
    return 2.0 * (re + 1j * im) / n


def c_katman(KP, i, disari=-1):
    N = int(KP["N"])
    nb = KP["nb"]
    if disari < 0:
        re, im, n = KP["re"][i].sum(0), KP["im"][i].sum(0), N
    else:
        re = KP["re"][i].sum(0) - KP["re"][i][disari]
        im = KP["im"][i].sum(0) - KP["im"][i][disari]
        n = N - nb[disari]
    return 2.0 * (re + 1j * im) / n


if __name__ == "__main__":
    ONK = onkayit()
    print("=" * 78)
    print(f"187b / K1 KATMAN DEFTERİ (gerçek)  [on-kayit {ONK['zaman']} "
          f"sha {ONK['sha256'][:12]}]")
    print("=" * 78, flush=True)
    L = ONK["L"]
    izg = ONK["katman_izgara_gercek"]          # [0.86 ... 1.20]
    NKAT = len(izg) - 1                        # 7 katman

    G = np.load(S184 / "K1_gercek.npz")
    OZg = np.load(S185 / "OZ_gercek.npz")
    Pg = np.load(S186 / "G1_proj_gercek.npz")
    w_win = np.asarray(G["w"], float)
    tau = np.asarray(G["tau"], float)
    aq = np.asarray(G["aq"], float)
    ae = np.asarray(G["aq_eff"], float)
    g, mid, Lg = b185.kinematik("gercek")
    ds_ref = np.asarray(np.load(S155 / "eta_son_t0.4_c4000.npz")["ds"], float)
    N = len(mid)
    kenar_jk = np.linspace(0, N, NJACK + 1).astype(int)
    if not np.array_equal(np.diff(kenar_jk), G["nb"]):
        raise SystemExit("JK BLOK KENARLARI 184 İLE UYUMSUZ")

    # --- katman çizgi listeleri ---
    qmax = np.exp(izg[-1] * L)
    qa, lam, taua, aqa = asal_kuvvetler(qmax, L)
    katman_idx = [np.where((taua > izg[i]) & (taua <= izg[i + 1]))[0]
                  for i in range(NKAT)]
    print(f"\nKatman çizgi sayıları (q≤{qmax:.0f}, toplam yeni "
          f"{sum(len(k) for k in katman_idx)}):")
    for i in range(NKAT):
        print(f"  katman {i+1} ({izg[i]:.2f},{izg[i+1]:.2f}]: "
              f"{len(katman_idx[i])} çizgi", flush=True)

    # --- AŞAMA A: katman serileri (kontrol noktalı) ---
    for i in range(NKAT):
        yol = S187 / f"katman_{i+1}_gercek.npy"
        if yol.exists():
            print(f"  [katman {i+1}] mevcut, atlanıyor", flush=True)
            continue
        idx = katman_idx[i]
        dds = katman_serisi(g, mid, np.log(qa[idx]), aqa[idx],
                            f"katman {i+1}")
        np.save(yol, dds)
        print(f"  -> {yol}", flush=True)

    # --- AŞAMA B: tek geçişte taban serisi + tüm katman izdüşümleri ---
    yolB = S187 / "katmanproj_gercek.npz"
    if not yolB.exists():
        t0 = time.time()
        katmanlar = [np.load(S187 / f"katman_{i+1}_gercek.npy")
                     for i in range(NKAT)]
        taban = np.zeros(N)
        re = np.zeros((NKAT, NJACK, len(w_win)))
        im = np.zeros((NKAT, NJACK, len(w_win)))
        dyy = dxy = 0.0
        kats_w = 2.0 * aq
        for b in range(NJACK):
            lo, hi = kenar_jk[b], kenar_jk[b + 1]
            for s0 in range(lo, hi, PROJ_BLOK):
                sl = slice(s0, min(s0 + PROJ_BLOK, hi))
                P = np.outer(mid[sl], w_win)
                cP = np.cos(P)
                sP = np.sin(P)
                del P
                Kw = np.sin(0.5 * np.outer(g[sl], w_win))
                Kw *= cP
                tb = Kw @ kats_w
                taban[sl] = tb
                dyy += (tb ** 2).sum()
                dxy += (tb * ds_ref[sl]).sum()
                for i in range(NKAT):
                    v = katmanlar[i][sl]
                    re[i, b] += v @ cP
                    im[i, b] -= v @ sP
                del cP, sP, Kw
            print(f"  [izdüşüm] jk-blok {b+1}/8  ({time.time()-t0:.0f}s)",
                  flush=True)
        # MAKİNE MÜHRÜ 1: taban momentleri 186 G1_proj(bir) ile örtüşmeli
        dyy186 = float(Pg["dyy_bir"].sum())
        dxy186 = float(Pg["dxy_bir"].sum())
        f1 = abs(dyy - dyy186) / dyy186
        f2 = abs(dxy - dxy186) / dxy186
        print(f"  MÜHÜR-1 taban momentleri vs 186: |Δdyy|/dyy={f1:.2e}  "
              f"|Δdxy|/dxy={f2:.2e}  (eşik 1e-6)", flush=True)
        np.save(S187 / "katman_0_gercek.npy", taban)
        np.savez_compressed(yolB, re=re, im=im, N=N, nb=np.diff(kenar_jk),
                            izgara=np.array(izg), dyy_taban=dyy,
                            dxy_taban=dxy, muhur1=np.array([f1, f2]))
        print(f"  -> {yolB}  ({time.time()-t0:.0f}s)", flush=True)
    else:
        print("  [izdüşüm] mevcut, atlanıyor", flush=True)

    KP = np.load(yolB)
    muhur1 = KP["muhur1"]

    # --- AŞAMA C: defter + jackknife + H-F1 ---
    kenar_b = json.load(open(S186 / "ONKAYIT_186.json"))["kenar"]
    maskeler = [(tau >= kenar_b[b]) & (tau < kenar_b[b + 1]) for b in range(8)]
    maskeler.append((tau >= kenar_b[0]) & (tau < kenar_b[-1]))   # HAVUZ
    et = [f"{kenar_b[b]:.2f}-{kenar_b[b+1]:.2f}" for b in range(8)] + ["HAVUZ"]

    def defter(disari):
        c_olc = b185.c_olculu(G, disari)
        a_oz = np.abs(b185.c_oz(OZg, aq, disari))
        ck = c_kesik(Pg, disari)
        ckat = [c_katman(KP, i, disari) for i in range(NKAT)]
        out = np.zeros((9, 2 + NKAT))    # m_olc, m_ya(0.86), m_ya(katman1..7)
        for k, m in enumerate(maskeler):
            sae = ae[m].sum()
            soz = a_oz[m].sum()
            out[k, 0] = (np.abs(c_olc[m]).sum() - soz) / sae
            cya = ck[m].copy()
            out[k, 1] = (np.abs(cya).sum() - soz) / sae
            for i in range(NKAT):
                cya = cya + ckat[i][m]
                out[k, 2 + i] = (np.abs(cya).sum() - soz) / sae
        return out

    tam = defter(-1)
    reps = np.array([defter(i) for i in range(NJACK)])
    se = b185.jk_se(reps)

    m_olc = tam[:, 0]
    mya = tam[:, 1:]                    # (9, 8): τ_c = 0.86..1.20
    dm = np.diff(tam[:, 1:], axis=1)    # (9, 7) katman artımları
    dm_reps = np.diff(reps[:, :, 1:], axis=2)
    se_dm = b185.jk_se(dm_reps)
    acik = mya - m_olc[:, None]         # (9, 8)
    kapanan_reps = 1.0 - (reps[:, :, -1] - reps[:, :, 0]) / \
        (reps[:, :, 1] - reps[:, :, 0])
    kapanan = 1.0 - acik[:, -1] / acik[:, 0]
    se_kapanan = b185.jk_se(kapanan_reps)

    # MAKİNE MÜHRÜ 2: m^ya(0.86) 186'nın g_bir m_pred'ini hane hane vermeli
    g186 = json.load(open(S186 / "G1_sonuc.json"))["g_bir"]["m_pred"]
    muhur2 = float(np.abs(mya[:8, 0] - np.array(g186)).max())
    print(f"\nMÜHÜR-2: m^ya(0.86) vs 186 g_bir m_pred maks|Δ| = {muhur2:.2e}")

    # korelasyon yan sütunu (tam-örneklem)
    taban = np.load(S187 / "katman_0_gercek.npy")
    korr = []
    dsy = taban.copy()
    ss = float((ds_ref ** 2).sum())
    for i in range(NKAT + 1):
        if i > 0:
            dsy += np.load(S187 / f"katman_{i}_gercek.npy")
        korr.append(float((dsy @ ds_ref) /
                          np.sqrt(ss * (dsy ** 2).sum())))

    tau_bar = np.array([np.sum(ae[m] * tau[m]) / ae[m].sum()
                        for m in maskeler])
    print("\nK1 YAKINSAMA DEFTERİ  m^ya(τ_c; bant)  [hedef m_ölç; ±jk se]:")
    bas = "  ".join(f"{t:.2f}" for t in izg)
    print(f"{'bant':>10} {'m_ölç':>8} | τ_c: {bas}")
    for k in range(9):
        s = " ".join(f"{mya[k, j]:.4f}" for j in range(NKAT + 1))
        print(f"{et[k]:>10} {m_olc[k]:8.4f} | {s}")
    print("\nKATMAN ARTIMLARI Δm (HAVUZ; ±jk se):")
    for i in range(NKAT):
        print(f"  katman {i+1} ({izg[i]:.2f},{izg[i+1]:.2f}]: "
              f"Δm = {dm[8, i]:+.4f} ± {se_dm[8, i]:.4f}   "
              f"|Δm|/2se = {abs(dm[8, i])/(2*se_dm[8, i]):.2f}")
    print("\nKORELASYON YAN SÜTUNU korr(ds^ya(τ_c), ds):")
    for j in range(NKAT + 1):
        print(f"  τ_c={izg[j]:.2f}: {korr[j]:.4f}")

    # H-F1 çatalı (HAVUZ)
    Hp = 8
    mono = bool(np.all(np.diff(acik[Hp]) < 0) and np.all(acik[Hp] > 0))
    cat_a_ciftler = [i for i in range(NKAT - 1)
                     if abs(dm[Hp, i]) < 2 * se_dm[Hp, i]
                     and abs(dm[Hp, i + 1]) < 2 * se_dm[Hp, i + 1]]
    cat_a = len(cat_a_ciftler) > 0
    kalan = float(acik[Hp, -1] / acik[Hp, 0])
    cat_b = bool(kalan > 0.5)
    kap_50 = bool(kapanan[Hp] >= 0.5)
    print("\nH-F1 ÇATALI (HAVUZ):")
    print(f"  monoton yaklaşma: {'EVET' if mono else 'HAYIR'}")
    print(f"  açık(0.86) = {acik[Hp, 0]:.4f}  açık(1.20) = {acik[Hp, -1]:.4f}")
    print(f"  kapanan pay = {kapanan[Hp]:.3f} ± {se_kapanan[Hp]:.3f}  "
          f"(≥%50: {'EVET' if kap_50 else 'HAYIR'})")
    print(f"  çatal (a) SIĞ İPTAL (art arda iki |Δm|<2se): "
          f"{'EVET — çiftler ' + str([(i+1, i+2) for i in cat_a_ciftler]) if cat_a else 'HAYIR'}")
    print(f"  çatal (b) DERİN KUYRUK (kalan açık >%50): "
          f"{'EVET' if cat_b else 'HAYIR'}  [kalan = {kalan:.3f}]")

    sonuc = {
        "sha_onkayit": ONK["sha256"], "izgara": izg,
        "katman_cizgi": [int(len(k)) for k in katman_idx],
        "etiket": et, "tau_bar": tau_bar.tolist(),
        "m_olc": m_olc.tolist(), "se_m_olc": se[:, 0].tolist(),
        "m_ya": mya.tolist(), "se_m_ya": se[:, 1:].tolist(),
        "dm": dm.tolist(), "se_dm": se_dm.tolist(),
        "acik": acik.tolist(),
        "kapanan": kapanan.tolist(), "se_kapanan": se_kapanan.tolist(),
        "korr": korr,
        "muhur1_taban_momentleri": muhur1.tolist(),
        "muhur2_maksfark_186": muhur2,
        "HF1": {"monoton": mono, "cat_a": cat_a,
                "cat_a_ciftler": [(i + 1, i + 2) for i in cat_a_ciftler],
                "cat_b": cat_b, "kalan_acik_pay": kalan,
                "kapanan_pay": float(kapanan[Hp]),
                "se_kapanan": float(se_kapanan[Hp])},
    }
    json.dump(sonuc, open(S187 / "K1_katman_defteri.json", "w"), indent=1,
              ensure_ascii=False)
    print(f"\n-> {S187/'K1_katman_defteri.json'}  BİTTİ", flush=True)
