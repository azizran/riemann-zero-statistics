# -*- coding: utf-8 -*-
"""
198b — κ_p YÜKSEKLİK: 197 MAKİNESİ, PENCERE + BLOK GEOMETRİSİ PARAMETRELİ
=========================================================================
KALEM_KAPPA_YUKSEKLIK_25EYL2026.md "Ölçüm tanımı". 197b_harita.py importlib ile yüklenir
(DÜZENLENMEZ); yapı taşları AYNEN: m3_oz_terim, kesik_katki, kesik_denetim, katki_kur, plan_kur,
omega_kos (b190._om_is), omega_birlestir (K_ham + S1b K_düz); 188b.karisim/jk_kenar/c_proj,
187b.asal_kuvvetler, 190b.kinematik_eta. Yeni olan: pencere kurulumu, blok geometrisi ve
grup-genelleştirilmiş profil_sigma_g / cizgi_bicimi_g (eşit grup geometrisinde 197b'ninkiyle
BİT-BİT aynı — M1'de sınanır).
 GEOMETRİ:
  '197'   32 eşit gap-sayılı blok (linspace(0,N,33)); jk grup = blok//4 (M1 kod yolu).
  'esitL' (KALEM teftiş sonrası) EŞİT-ΔL: n_B = max(32, ⌈ΔL_W/0.03⌉), hedef genişlik
          w = ΔL_W/n_B; jk 8 bitişik grup = zincirin 8 eşit gap-sayılı bloğu (karışım loo'su
          yalnız bu kenarlarla tutarlı); her grup kendi L aralığında ⌈ΔL_g/w⌉ eşit-L bloğa
          bölünür (grup içi eşit genişlik ≤ w; toplam blok ≥ n_B olabilir).
 Her pencerede: HAVUZ' = K1 çizgileri τ = log q/L_W ∈ [0.45, 0.74) (+ M7 yarıları [0.45,0.60),
 [0.60,0.74) ve ortak mutlak havuz q ∈ [224, 992] AYNI koşuda ayrı maskeler); K_düz BİRİNCİL;
 tarama τ' ∈ (0.74, τ'_üst], τ'_üst = max(1.30, +2.40 üst kapsaması için gereken, 0.005↑);
 Δω tam dilimleri j ∈ [−56, 96]; blok başına taranan çizgiler = tam dilimlerin TÜM çizgileri.
 Kapsama 192 kuralı (B_j = dilimi TAM kapsayan bloklar). KALEM TAM kapsama ister: varsayılan
 'tam' (eksikse DURUR); '--kapsama kismi' yalnız kaptan kararıyla (ön-kayda yazılır).
Pencereler: dusuk (190 zinciri; REFERANS), son (190 muhur_son; tam örneklem), alt (198k0).
Modlar:
  plan <W>     kinematik + geometri + kapsama + çizgi sayıları (κ YOK)
  harita <W>   ω-haritası; W ∈ {alt, son}: ONKAYIT_198 + M1 + M6 ŞART; W = dusuk serbest.
               EKRANA κ/K/P YAZILMAZ.
Kullanım: 198b_harita.py <plan|harita> <alt|dusuk|son> [--geo 197|esitL] [--kapsama tam|kismi]
          [--isci N]
"""
import argparse
import hashlib
import importlib.util
import json
import math
import sys
import time
from pathlib import Path
from types import SimpleNamespace

sys.dont_write_bytecode = True
import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = QM / "scratchpad"
S190, S198 = SCR / "190", SCR / "198"
KALEM = QM / "KALEM_KAPPA_YUKSEKLIK_25EYL2026.md"
TWO_PI = 2 * np.pi
NJACK, DW, EPS = 8, 0.025, 1e-9
HAVUZ = (0.45, 0.74)
M7_YARILAR = ((0.45, 0.60), (0.60, 0.74))
MUTLAK_Q = (224, 992)
MASKE_ADLARI = ("tam", "yA", "yB", "mutlak")
TAU_ALT = 0.74
TAU_UST_MIN = 1.30
DW_ARALIK = (-1.40, 2.40)
NB_MIN, DL_HEDEF = 32, 0.03
PENCERELER = {
    "dusuk": {"zincir": S190 / "zincir_dusuk", "etiket": "gercek_dusuk",
              "eta": "eta_dusuk_t0.4_c4000.npz", "aralik": "zeros6 Z[200000:500000]"},
    "son": {"zincir": S190 / "muhur_son", "etiket": "gercek_sonM",
            "eta": "eta_son_t0.4_c4000.npz", "aralik": "zeros6 son 300k"},
    "alt": {"zincir": S198 / "zincir_alt", "etiket": "gercek_alt",
            "eta": "eta_alt_t0.4_c4000.npz", "aralik": "zeros6 Z[20000:200000]"},
}


def yukle(ad, yol):
    spec = importlib.util.spec_from_file_location(ad, yol)
    m = importlib.util.module_from_spec(spec)
    sys.modules[ad] = m
    spec.loader.exec_module(m)
    return m


h197 = yukle("h197b", QM / "197_configs" / "197b_harita.py")
b190, b188, b187 = h197.b190, h197.b188, h197.b187


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def jaralik():
    return int(round(DW_ARALIK[0] / DW)), int(round(DW_ARALIK[1] / DW))


def dizin(W, geo):
    return S198 / f"{W}_{geo}"


# ============================ geometri ============================
def geometri(mid, kj8, geo):
    """→ kb (blok kenarları), Lb (blok içi mean log(m/2π)), grup (blok → jk grubu)."""
    N = len(mid)
    Ln = np.log(mid / TWO_PI)
    if geo == "197":
        kb = np.linspace(0, N, 33).astype(int)
        grup = np.arange(32) // 4
        assert np.array_equal(kb[::4], kj8)
    else:
        dLW = float(Ln[-1] - Ln[0])
        nB = max(NB_MIN, math.ceil(dLW / DL_HEDEF - 1e-12))
        w = dLW / nB
        kenar, grup = [0], []
        for g in range(NJACK):
            lo, hi = int(kj8[g]), int(kj8[g + 1])
            L0, L1 = Ln[lo], Ln[hi - 1]
            ng = max(1, math.ceil((L1 - L0) / w - 1e-9))
            for k in range(1, ng):
                kenar.append(lo + int(np.searchsorted(Ln[lo:hi], L0 + k * (L1 - L0) / ng)))
            kenar.append(hi)
            grup += [g] * ng
        kb = np.array(kenar)
        grup = np.array(grup)
        assert np.all(np.diff(kb) > 0)
    Lb = np.array([float(np.log(mid[kb[b]:kb[b + 1]] / TWO_PI).mean())
                   for b in range(len(kb) - 1)])
    return kb, Lb, grup


def profil_sigma_g(K, nb, kap, grup):
    """197b.profil_sigma AYNEN, jk grubu açık (grup[b]); replika i'de grup i blokları çıkar."""
    kb_ = -K.real
    w = nb[:, None] * kap
    out = np.full((NJACK + 1, K.shape[2]), np.nan)
    with np.errstate(invalid="ignore", divide="ignore"):
        out[0] = (kb_[0] * w).sum(0) / w.sum(0)
        for i in range(NJACK):
            wi = w * (grup != i)[:, None]
            out[1 + i] = (kb_[1 + i] * wi).sum(0) / wi.sum(0)
    return out


def cizgi_bicimi_g(mid, kb, Lb, kap, J, log_r, grup):
    """197b.cizgi_bicimi AYNEN, jk grubu açık."""
    nb = np.diff(kb)
    Ln = np.log(mid / TWO_PI)
    hb = np.zeros((len(Lb), len(J)))
    for b in range(len(Lb)):
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


def kapsama(J, L, Lb, tau_ust, tam=None):
    mt = J * DW
    k = np.array([(mt - 0.5 * DW >= TAU_ALT * L - Lb[b] - EPS) &
                  (mt + 0.5 * DW <= tau_ust * L - Lb[b] + EPS) for b in range(len(Lb))])
    return k if tam is None else k & tam


# ============================ pencere kurulumu ============================
def ortak_W(W, geo):
    """197b.ortak() muadili (aynı öznitelik adları) + geometri."""
    pw = PENCERELER[W]
    Z = pw["zincir"]
    G = np.load(Z / f"K1_{pw['etiket']}.npz")
    OZ = np.load(Z / f"OZ_{pw['etiket']}.npz")
    P = np.load(Z / f"G1_proj_{pw['etiket']}.npz")
    g, mid, L = b190.kinematik_eta(Z / pw["eta"])           # 185b 'gercek' kolu AYNEN
    assert float(G["L"]) == L and float(P["L"]) == L
    N = len(mid)
    kj8 = b188.jk_kenar(N)
    if not (np.array_equal(np.diff(kj8), G["nb"]) and np.array_equal(np.diff(kj8), P["nb"])
            and np.array_equal(np.diff(kj8), OZ["nb"])):
        raise SystemExit("JK BLOK KENARLARI ZİNCİRLE UYUMSUZ")
    bid8 = np.searchsorted(kj8, np.arange(N), side="right") - 1
    kb, Lb, grup = geometri(mid, kj8, geo)
    for gg in range(NJACK):                              # grup kenarı = zincir 8-blok kenarı
        bb = np.where(grup == gg)[0]
        assert kb[bb[0]] == kj8[gg] and kb[bb[-1] + 1] == kj8[gg + 1]
    jlo, jhi = jaralik()
    gerek = float(np.max(Lb) + (jhi + 0.5) * DW) / L
    tau_ust = round(max(TAU_UST_MIN, math.ceil(gerek / 0.005 - 1e-9) * 0.005), 3)
    qa, _, taua, aqa = b187.asal_kuvvetler(np.exp(tau_ust * L), L)
    q_all = np.asarray(G["q"], float)
    k1_satir = np.full(len(qa), -1, int)
    ix = np.searchsorted(qa, q_all)
    assert np.array_equal(qa[ix], q_all), "K1 çizgileri asal-kuvvet listesinde değil"
    k1_satir[ix] = np.arange(len(q_all))
    return SimpleNamespace(
        W=W, geo=geo, L=L, N=N, OZ=OZ, P=P, g=g, mid=mid, kj8=kj8, bid8=bid8, kb=kb, Lb=Lb,
        grup=grup, w_all=np.asarray(G["w"], float), tau_all=np.asarray(G["tau"], float),
        aq_all=np.asarray(G["aq"], float), q_all=q_all, k1_satir=k1_satir,
        qa=qa, taua=taua, aqa=aqa, tau_ust=tau_ust, tau_ust_gerek=gerek)


def maskeler_W(c, win1):
    tw, qw = c.tau_all[win1], np.round(c.q_all[win1])
    m = [np.ones(int(win1.sum()), bool)]
    m += [(tw >= a) & (tw < b) for a, b in M7_YARILAR]
    m.append((qw >= MUTLAK_Q[0]) & (qw <= MUTLAK_Q[1]))
    return m


def plan_W(c):
    win1 = (c.tau_all >= HAVUZ[0]) & (c.tau_all < HAVUZ[1])
    tarama = np.where((c.taua > TAU_ALT) & (c.taua <= c.tau_ust))[0]
    m3 = h197.m3_oz_terim(c.q_all[win1], c.qa[tarama], f"{c.W}/{c.geo}")
    jlo, jhi = jaralik()
    J = np.arange(jlo, jhi + 1)
    kap = kapsama(J, c.L, c.Lb, c.tau_ust)
    Bj = kap.sum(0)
    nb = np.diff(c.kb)
    eksik = {int(b): [round(float((J[~kap[b]] * DW).min()), 3),
                      round(float((J[~kap[b]] * DW).max()), 3)]
             for b in range(len(c.Lb)) if not kap[b].all()}
    wpay = (nb[:, None] * kap).sum(0) / nb.sum()
    sec = {}
    lq = np.log(c.qa[tarama])
    for b in range(len(c.Lb)):
        jt = np.floor((lq - c.Lb[b]) / DW + 0.5).astype(int)
        sec[b] = tarama[(jt >= jlo) & (jt <= jhi)]
    Ln = np.log(c.mid / TWO_PI)
    yay = np.array([Ln[c.kb[b + 1] - 1] - Ln[c.kb[b]] for b in range(len(c.Lb))])
    mk = maskeler_W(c, win1)
    kal = (J * DW >= -1.30 - EPS) & (J * DW <= -0.55 + EPS)
    rap = {"pencere": c.W, "geometri": c.geo, "aralik": PENCERELER[c.W]["aralik"],
           "L_W": c.L, "N": c.N, "L_min": float(Ln.min()), "L_max": float(Ln.max()),
           "n_blok": int(len(c.Lb)), "blok_gap_sayisi": [int(x) for x in nb],
           "grup": c.grup.tolist(), "L_b": c.Lb.tolist(), "blok_ici_L_yayilimi": yay.tolist(),
           "havuz_ussu": {"aralik": list(HAVUZ), "n_cizgi": int(win1.sum()),
                          "maske_cizgi": dict(zip(MASKE_ADLARI, [int(x.sum()) for x in mk]))},
           "K1_evreni": int(len(c.w_all)), "tau_ust": c.tau_ust,
           "tau_ust_gerek_2.4": c.tau_ust_gerek, "n_tarama": int(len(tarama)),
           "tarama_kesik_evreninde": int((c.k1_satir[tarama] >= 0).sum()),
           "M3": m3, "j_tam": [jlo, jhi], "B_j_min": int(Bj.min()), "n_blok_toplam": len(c.Lb),
           "tam_kapsama": bool(Bj.min() == len(c.Lb)), "kapsamayan_bloklar": eksik,
           "kapsanan_gap_payi_min": float(wpay.min()),
           "kapsanan_gap_payi_kalibrasyon_min": float(wpay[kal].min()),
           "cizgi_blok": [int(len(sec[b])) for b in range(len(c.Lb))],
           "nokta_cizgi": float(np.sum(nb * np.array([len(sec[b]) for b in range(len(c.Lb))])))}
    return win1, tarama, sec, J, rap


def mod_plan(W, geo):
    c = ortak_W(W, geo)
    rap = plan_W(c)[-1]
    d = dizin(W, geo)
    d.mkdir(parents=True, exist_ok=True)
    json.dump(rap, open(d / "plan_198.json", "w"), indent=1, ensure_ascii=False)
    print(f"  [{W}/{geo}] L_W={rap['L_W']:.9f} N={rap['N']} L∈[{rap['L_min']:.4f},"
          f"{rap['L_max']:.4f}]; n_blok={rap['n_blok']} (gap/blok {min(rap['blok_gap_sayisi'])}-"
          f"{max(rap['blok_gap_sayisi'])}); blok-içi L yayılımı "
          f"{min(rap['blok_ici_L_yayilimi']):.4f}-{max(rap['blok_ici_L_yayilimi']):.4f}; "
          f"HAVUZ' {rap['havuz_ussu']['maske_cizgi']}; τ'_üst {rap['tau_ust']}", flush=True)
    print(f"      kapsama: B_j min {rap['B_j_min']}/{rap['n_blok']} → "
          f"{'TAM' if rap['tam_kapsama'] else 'EKSİK'}; kapsanan gap payı min "
          f"{rap['kapsanan_gap_payi_min']:.4f} (kalibrasyon penceresinde "
          f"{rap['kapsanan_gap_payi_kalibrasyon_min']:.4f}); kapsamayan bloklar "
          f"{rap['kapsamayan_bloklar']}", flush=True)
    print(f"      Σ n_b·çizgi = {rap['nokta_cizgi']:.3e}  -> {d/'plan_198.json'}", flush=True)


# ============================ harita (197b olcum C + S1b) ============================
def onkayit198():
    yol = S198 / "ONKAYIT_198.json"
    if not yol.exists():
        raise SystemExit("ONKAYIT_198.json YOK — K0 mühürlenmeden kör pencere ölçülmez")
    o = json.load(open(yol))
    if sha(QM / "198_configs" / "198a_onkayit.py") != o["sha256"]:
        raise SystemExit("ON-KAYIT SHA UYUMSUZ")
    if sha(KALEM) != o["kalem_sha256"]:
        raise SystemExit("KALEM ön-kayıttan sonra DEĞİŞMİŞ")
    for rel, h in o["girdi_sha256"].items():
        if sha(QM / rel) != h:
            raise SystemExit(f"GİRDİ SHA UYUMSUZ: {rel}")
    for ad in ("M1_198.json", "M6_on_198.json"):
        if not json.load(open(S198 / ad)).get("gecti"):
            raise SystemExit(f"{ad} GEÇMEDİ — ölçüm başlamaz")
    return o


def mod_harita(W, geo, kapsama_mod, nw):
    t00 = time.time()
    if W != "dusuk":
        o = onkayit198()
        assert o["olcum_parametreleri"]["geometri"] == geo
        assert o["olcum_parametreleri"]["kapsama"][W] == kapsama_mod
    c = ortak_W(W, geo)
    win1, tarama, sec, Jp, rap = plan_W(c)
    if not rap["tam_kapsama"] and kapsama_mod != "kismi":
        raise SystemExit(f"[{W}/{geo}] Δω {list(DW_ARALIK)} TAM kapsanmıyor (B_j min "
                         f"{rap['B_j_min']}/{rap['n_blok']}): {rap['kapsamayan_bloklar']} — "
                         f"KALEM gereği durur ('--kapsama kismi' yalnız kaptan kararıyla)")
    d = dizin(W, geo)
    d.mkdir(parents=True, exist_ok=True)
    nblok = len(c.Lb)
    w1 = c.w_all[win1]
    mk = maskeler_W(c, win1)
    mixes = ([b188.karisim(c.P, c.OZ, c.aq_all)[win1]] +
             [b188.karisim(c.P, c.OZ, c.aq_all, i)[win1] for i in range(NJACK)])
    yolM = d / "kesik_katki.npz"
    if not yolM.exists():
        re_k, im_k = h197.kesik_katki(c, w1, nw)
        np.savez_compressed(yolM, re=re_k, im=im_k, w_K1=c.w_all, w_havuz=w1)
    KM = np.load(yolM)
    re_k, im_k = KM["re"], KM["im"]
    mkes = h197.kesik_denetim(c, re_k, im_k, win1)
    print(f"  M-kesik Σ_K1 katkı ≡ c^kesik (HAVUZ'): maks bağıl = {max(mkes):.2e}", flush=True)
    assert max(mkes) <= 1e-10
    Mv = [b188.c_proj(re_k, im_k, np.diff(c.kj8), c.N, v - 1) for v in range(NJACK + 1)]
    del re_k, im_k
    katki = h197.katki_kur(c.k1_satir, Mv)
    om = d / "omega"
    om.mkdir(exist_ok=True)
    isler, plan, say = h197.plan_kur(c.qa, c.aqa, sec, c.Lb, om)
    print(f"  [{W}/{geo}] {nblok} blok; çizgi/blok {min(len(sec[b]) for b in range(nblok))}-"
          f"{max(len(sec[b]) for b in range(nblok))}; eksik görev {len(isler)} (işçi {nw})",
          flush=True)
    t_om = h197.omega_kos(isler, c, w1, c.kb, nw, f"{W}/{geo}")
    J, Kh, Kd, n, nkes, nb_om = h197.omega_birlestir(om, plan, say, nblok, mixes, mk,
                                                     sec, c.Lb, c, katki)
    nb = np.diff(c.kb)
    assert np.array_equal(nb_om, nb)
    jlo, jhi = jaralik()
    assert J[0] == jlo and J[-1] == jhi, (J[0], J[-1])
    tam = np.ones(len(J), bool)
    kap = kapsama(J, c.L, c.Lb, c.tau_ust, tam)
    kay = {}
    for k, ad in enumerate(MASKE_ADLARI):
        kay[f"kappa_sigma_duz_{ad}"] = profil_sigma_g(Kd[:, :, k, :], nb, kap, c.grup)
        kay[f"kappa_sigma_ham_{ad}"] = profil_sigma_g(Kh[:, :, k, :], nb, kap, c.grup)
    np.savez_compressed(
        d / "harita_omega.npz", birincil=np.array("duz"), pencere=np.array(W),
        geometri=np.array(geo), maske_adlari=np.array(MASKE_ADLARI),
        K_duz=Kd, K_ham=Kh, **kay, ncz=n, n_kesik=nkes, nb=nb, kb=c.kb, L_b=c.Lb,
        grup=c.grup, L_W=c.L, N=c.N, kapsama=kap, tam=tam, J=J, merkez=J * DW, dw=DW,
        dw_aralik=np.array(DW_ARALIK), havuz=np.array(HAVUZ), tau_ust=c.tau_ust, w_havuz=w1,
        kapsama_mod=np.array(kapsama_mod),
        etiket_v=np.array(["tam"] + [f"loo{i}" for i in range(NJACK)]))
    sonuc = {"pencere": W, "geometri": geo, "kapsama_mod": kapsama_mod, "L_W": c.L, "N": c.N,
             "n_blok": nblok, "M3": rap["M3"], "M_kesik": mkes, "tau_ust": c.tau_ust,
             "havuz_ussu": rap["havuz_ussu"], "B_j_min": rap["B_j_min"],
             "sure_omega_s": t_om, "sure_toplam_s": time.time() - t00,
             "kod_sha256": sha(Path(__file__).resolve()),
             "zaman": time.strftime("%Y-%m-%d %H:%M:%S %z")}
    json.dump(sonuc, open(d / "K1_harita_198.json", "w"), indent=1, ensure_ascii=False)
    print(f"-> {d/'harita_omega.npz'} BİTTİ ({time.time()-t00:.0f}s) — κ/P EKRANA YAZILMADI",
          flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("mod", choices=["plan", "harita"])
    ap.add_argument("pencere", choices=list(PENCERELER))
    ap.add_argument("--geo", choices=["197", "esitL"], default="esitL")
    ap.add_argument("--kapsama", choices=["tam", "kismi"], default="tam")
    ap.add_argument("--isci", type=int, default=6)
    a = ap.parse_args()
    print("=" * 78)
    print(f"198b / {a.mod.upper()} [{a.pencere}/{a.geo}]  {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 78, flush=True)
    if a.mod == "plan":
        mod_plan(a.pencere, a.geo)
    else:
        mod_harita(a.pencere, a.geo, a.kapsama, a.isci)
