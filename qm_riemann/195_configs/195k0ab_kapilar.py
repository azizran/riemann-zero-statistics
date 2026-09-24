# -*- coding: utf-8 -*-
"""
195k0ab — K0a (C_χ MUTLAK KALİBRASYONU) + K0b (KARAKTER KİMLİĞİ)
================================================================
KALEM_UYDU_KARAKTERI_24EYL2026 kapıları; ölçümden ÖNCE; biri tutmazsa DUR.

K0a: her ada için m_n = (z_n+z_{n+1})/2; dosya t≈200'den başladığı için ilk sıfır
     GERÇEK ilk sıfır DEĞİL ⇒ indeks kaydırması j0 (z_1'den önceki sıfır sayısı):
       (i)  birincil: j0 = round(mean_i[N̄_χ(m_i) − i])  (görev: N̄_χ'den yuvarlama)
       (ii) tek-nokta: round(N̄_χ(z_1) − 1/2)            (KAYIT)
       (iii) BAĞIMSIZ SAYIM (KAYIT): mpmath ile Hardy Z_χ(t) = e^{iθ_χ(t)} L(½+it,χ)
             işaret değişimleri (0, m_1) aralığında (adım 0.02 + yakın-kaçış
             incelemesi 0.001) ⇒ N(m_1) = j0 + 1 olmalı.
     n = j0 + i;  kapı: |mean[N̄_χ(m_n) − n]| < 0.05 (C_χ = −χ(−1)/8).
     Zeta (son, düşük): n = küresel sıfır indeksi (sıfır dosyasıyla eşleşme
     doğrulanır), C = 7/8; aynı kapı.
K0b: 118a (118b'nin exec ettiği üretici) metninden χ tabloları; değerler,
     parite, gerçellik, tam-çarpımsallık, İLKELLİK, Gauss toplamı τ(χ) assert;
     105b (χ₅ₑ/χ₈ₑ/χ₈ₒ'nin ilk üreticisi) tablolarıyla ve npz 'q','a' ile eşleşme.
Çıktı: scratchpad/195/K0ab.json
Kullanım: 195k0ab_kapilar.py   (nohup; mpmath sayımı ~birkaç dk)
"""
import importlib.util
import json
import re
import sys
import time
import ast
from math import gcd
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
spec = importlib.util.spec_from_file_location(
    "o195", Path(__file__).resolve().parent / "195o_ortak.py")
o = importlib.util.module_from_spec(spec)
spec.loader.exec_module(o)


def log(*a):
    print(*a, flush=True)


# ------------------------------------------------------------ mpmath sayımı
def _Z_say(arg):
    ad, k, tab, par, t_ust, z1, adim = arg
    import mpmath as mp
    mp.mp.dps = 20
    chi = [tab[r] for r in range(k)]

    def Z(t):
        s = mp.mpc(0.5, t)
        th = t / 2 * mp.log(mp.mpf(k) / mp.pi) + mp.im(mp.loggamma((s + par) / 2))
        v = mp.exp(1j * th) * mp.dirichlet(s, chi)
        return float(mp.re(v)), float(mp.im(v))

    t0 = time.time()
    ts = np.arange(adim, t_ust, adim)
    vals = np.zeros(len(ts))
    ims = np.zeros(len(ts))
    for i, t in enumerate(ts):
        vals[i], ims[i] = Z(float(t))
    sg = np.sign(vals)
    deg = np.where(sg[:-1] * sg[1:] < 0)[0]
    # yakın-kaçış: |Z| yerel minimumu, komşular aynı işaretli, küçük değer
    av = np.abs(vals)
    inc = 0
    ekstra = 0
    for i in range(1, len(av) - 1):
        if av[i] < av[i - 1] and av[i] < av[i + 1] and sg[i - 1] == sg[i] == sg[i + 1]:
            yerel = max(av[max(0, i - 25):i + 25].max(), 1e-300)
            if av[i] < 0.1 * yerel:
                inc += 1
                tt = np.arange(ts[i - 1], ts[i + 1], 0.001)
                vv = np.array([Z(float(x))[0] for x in tt])
                ekstra += int(np.sum(np.sign(vv[:-1]) * np.sign(vv[1:]) < 0))
    say = int(len(deg) + ekstra)
    rel_im = float(np.max(np.abs(ims)) / max(np.max(np.abs(vals)), 1e-300))
    # z_1 civarında işaret değişimi var mı?
    z1_yakin = bool(np.any(np.abs(ts[deg] + adim / 2 - z1) < adim))
    return ad, dict(sayim_0_m1=say, grid_isaret_degisimi=int(len(deg)),
                    yakin_kacis_incelenen=inc, yakin_kacis_ekstra=ekstra,
                    maks_goreli_Im=rel_im, z1_isaret_degisimi_goruldu=z1_yakin,
                    adim=adim, t_ust=float(t_ust), sure_s=time.time() - t0)


def K0b():
    tab, satir = o.karakter_tablolari()
    out = {}
    beklenen = {
        "beta": (4, 1, {1: 1, 3: -1}),
        "chi3": (3, 1, {1: 1, 2: -1}),
        "chi5e": (5, 0, {r: (1 if pow(r, 2, 5) == 1 or r in (1, 4) else -1)
                         for r in range(1, 5)}),   # Legendre (·/5)
        "chi8e": (8, 0, {1: 1, 3: -1, 5: -1, 7: 1}),
        "chi8o": (8, 1, {1: 1, 3: 1, 5: -1, 7: -1}),
    }
    # Legendre (·/5) Euler ölçütüyle bağımsız
    leg5 = {r: (1 if pow(r, 2, 5) == 1 else -1) for r in range(1, 5)}
    assert leg5 == {1: 1, 2: -1, 3: -1, 4: 1}
    beklenen["chi5e"] = (5, 0, leg5)
    # 105b tabloları (ilk üretici)
    t105 = (o.QM / "105b_yeni_ada_kampanya.py").read_text()
    tab105 = {}
    for ad, isim in [("chi5e", "CHI5E"), ("chi8e", "CHI8E"), ("chi8o", "CHI8O")]:
        m = re.search(rf"^{isim} = (\{{.*?\}})\s*$", t105, re.M)
        tab105[ad] = {int(a): int(b) for a, b in ast.literal_eval(m.group(1)).items()}
    # 118b gerçekten 118a'yı exec ediyor mu
    t118b = (o.QM / "118b_yedi_ada_taperli_dokum.py").read_text()
    exec_118a = '118a_taperli_L_motoru.py").read().split(' in t118b
    for ad in o.ADALAR:
        k, par, deg = beklenen[ad]
        T = tab[ad]
        kk, isim, apar = satir[ad]
        d = np.load(o.QM / f"118b_{ad}_zeros.npz")
        assert kk == k and apar == par, (ad, "118a ADALAR satırı", kk, apar)
        assert int(d["q"]) == k and int(d["a"]) == par, (ad, "npz q/a")
        assert set(T.keys()) == set(range(k)), (ad, "tablo tam değil")
        for r in range(k):
            if gcd(r, k) > 1:
                assert T[r] == 0, (ad, r, "gcd>1 ama χ≠0")
            else:
                assert T[r] == deg[r], (ad, r, T[r], deg[r], "DEĞER UYUŞMAZ")
                assert T[r] in (1, -1), (ad, "gerçel değil")
        # tam çarpımsallık
        for x in range(k):
            for y in range(k):
                assert T[(x * y) % k] == T[x] * T[y], (ad, x, y, "çarpımsal değil")
        # parite
        assert T[k - 1] == (1 if par == 0 else -1), (ad, "parite uyumsuz")
        # ilkellik: her öz bölen d|k (d<k) için {m≡1 (mod d), gcd(m,k)=1} üzerinde χ≢1
        ilkel = True
        for dv in [dd for dd in range(1, k) if k % dd == 0]:
            if all(T[m] == 1 for m in range(1, k) if gcd(m, k) == 1 and m % dv == 1 % dv):
                ilkel = False
        assert ilkel, (ad, "İLKEL DEĞİL")
        # Gauss toplamı
        tau = sum(T[r] * np.exp(2j * np.pi * r / k) for r in range(k))
        beklenen_tau = np.sqrt(k) if par == 0 else 1j * np.sqrt(k)
        assert abs(tau - beklenen_tau) < 1e-12, (ad, tau, beklenen_tau)
        if ad in tab105:
            assert tab105[ad] == T, (ad, "105b tablosu ile uyuşmaz")
        out[ad] = {"tablo_118a": isim, "k": k, "parite": par,
                   "chi(-1)": T[k - 1], "degerler": {str(r): T[r] for r in range(k)},
                   "ilkel": ilkel, "gauss_tau": [float(tau.real), float(tau.imag)],
                   "105b_eslesme": (ad in tab105),
                   "C_chi": -T[k - 1] / 8.0, "sqrt_k_bolu_phi_k":
                       float(np.sqrt(k) / sum(1 for r in range(1, k) if gcd(r, k) == 1))}
        log(f"  K0b {ad:>6}: k={k} parite={par} χ(−1)={T[k-1]:+d} değerler="
            f"{[T[r] for r in range(k)]} ilkel={ilkel} τ(χ)={tau:.6f} "
            f"{'(105b ✓)' if ad in tab105 else ''}")
    out["118b_exec_118a"] = exec_118a
    assert exec_118a
    return out, tab


if __name__ == "__main__":
    T0 = time.time()
    log("=" * 78)
    log("195k0ab — K0a C_χ MUTLAK KALİBRASYONU + K0b KARAKTER KİMLİĞİ")
    log("=" * 78)
    sonuc = {"zaman_baslangic": time.ctime()}

    # ---------------- K0b ----------------
    log("\nK0b — KARAKTER KİMLİĞİ (118a üretici metni; 105b; npz):")
    k0b, tab = K0b()
    sonuc["K0b"] = k0b
    sonuc["K0b_GECTI"] = True
    log("  K0b: TÜM assert'ler GEÇTİ")

    # ---------------- K0a (adalar) ----------------
    log("\nK0a — C_χ MUTLAK KALİBRASYONU (adalar):")
    import multiprocessing as mp
    isler = []
    kin = {}
    for ad in o.ADALAR:
        K = o.ada_kin(ad)
        kin[ad] = K
        isler.append((ad, K["k"], tab[ad], K["a"], float(K["mid"][0]),
                      float(K["z"][0]), 0.02))
    ctx = mp.get_context("fork")
    with ctx.Pool(len(isler)) as pool:
        mp_sonuc = dict(pool.map(_Z_say, isler))

    k0a = {}
    gecti = True
    for ad in o.ADALAR:
        K = kin[ad]
        k, par = K["k"], K["a"]
        C = o.C_chi(par)
        mid, z = K["mid"], K["z"]
        i = np.arange(1, len(mid) + 1)
        off = o.Nbar(mid, k, C) - i
        j0 = int(np.round(off.mean()))
        j0_tek = int(np.round(o.Nbar(z[0], k, C) - 0.5))
        n = j0 + i
        r_all = float(np.mean(o.Nbar(mid, k, C) - n))
        ust = K["Lm"] >= o.L_UST
        r_ust = float(np.mean(o.Nbar(mid[ust], k, C) - n[ust]))
        # yanlış-parite kontrolü (KAYIT): C yerine −C
        r_yanlis = float(np.mean(o.Nbar(mid, k, -C) - n))
        m = mp_sonuc[ad]
        bagimsiz_j0 = m["sayim_0_m1"] - 1
        ok = abs(r_all) < 0.05
        gecti &= ok
        k0a[ad] = {"k": k, "parite": par, "C_chi": C, "z_1": float(z[0]),
                   "ilk_sifir_gercek_ilk_mi": False,
                   "j0_birincil_ortalama_yuvarlama": j0,
                   "ham_ortalama_kayma": float(off.mean()),
                   "j0_tek_nokta_Nbar_z1": j0_tek,
                   "j0_bagimsiz_mpmath": bagimsiz_j0,
                   "mpmath": m,
                   "j0_uyum": bool(j0 == bagimsiz_j0 == j0_tek),
                   "ort_Nbar_eksi_n_TUM": r_all, "ort_Nbar_eksi_n_UST": r_ust,
                   "N_orta_nokta": int(len(mid)), "N_ust": int(ust.sum()),
                   "yanlis_parite_kontrol_KAYIT": r_yanlis,
                   "kapi_0.05": ok}
        log(f"  {ad:>6}: z_1={z[0]:.4f} (ilk sıfır DEĞİL); j0 ort-yuvarlama={j0} "
            f"(ham {off.mean():+.4f}), tek-nokta={j0_tek}, mpmath bağımsız={bagimsiz_j0} "
            f"[{m['grid_isaret_degisimi']} işaret değ. + {m['yakin_kacis_ekstra']} "
            f"yakın-kaçış; maks|Im|/|Z|={m['maks_goreli_Im']:.1e}; {m['sure_s']:.0f}s]")
        log(f"          mean[N̄−n] TÜM={r_all:+.5f}  ÜST(L≥8.5)={r_ust:+.5f}  "
            f"(yanlış parite: {r_yanlis:+.4f})  → {'GEÇTİ' if ok else 'KALDI'}")

    # ---------------- K0a (zeta) ----------------
    log("\nK0a — zeta (C = 7/8):")
    Z = np.load(o.QM / "128_odl_zeros6_2e6_zeros.npz")["zeros"]
    assert abs(Z[0] - 14.134725141734693) < 1e-6, "zeros6 ilk sıfır değil"
    for pen in ["son", "dusuk"]:
        K = o.zeta_kin(pen)
        mid = K["mid"]
        i0 = int(np.searchsorted(Z, mid[0]))
        assert i0 == K["n0"], (pen, i0, K["n0"])
        zz = Z[i0 - 1:i0 - 1 + len(mid) + 1]
        assert np.array_equal(0.5 * (zz[:-1] + zz[1:]), mid), "orta nokta uyuşmaz"
        g_maks = float(np.max(np.abs(np.diff(zz) - K["g"])))
        n = K["n0"] + np.arange(len(mid))
        r = float(np.mean(o.Nbar(mid, 1, 7.0 / 8) - n))
        ok = abs(r) < 0.05
        gecti &= ok
        k0a[f"zeta_{pen}"] = {"n0": K["n0"], "ort_Nbar_eksi_n": r, "kapi_0.05": ok,
                              "maks_fark_g_eta_vs_sifir": g_maks,
                              "N_orta_nokta": int(len(mid))}
        log(f"  zeta {pen:>5}: n0={K['n0']} (sıfır dosyası indeksi ✓), "
            f"mean[N̄−n]={r:+.5f} → {'GEÇTİ' if ok else 'KALDI'}  "
            f"(g: eta vs sıfır maks|Δ|={g_maks:.1e})")

    sonuc["K0a"] = k0a
    sonuc["K0a_GECTI"] = bool(gecti)
    sonuc["zaman_bitis"] = time.ctime()
    json.dump(sonuc, open(o.S195 / "K0ab.json", "w"), indent=1, ensure_ascii=False)
    log(f"\nK0a: {'GEÇTİ' if gecti else 'KALDI — DUR'};  K0b: GEÇTİ")
    log(f"-> {o.S195/'K0ab.json'}  ({time.time()-T0:.0f}s)")
