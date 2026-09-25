# -*- coding: utf-8 -*-
"""
198c — HÜKÜM: κ_p'NİN YÜKSEKLİK BAĞIMLILIĞI (KALEM_KAPPA_YUKSEKLIK_25EYL2026.md)
================================================================================
197c fonksiyonları AYNEN (importlib; DÜZENLEME YOK): uyum_kal, ikincil, par, ozet, jk, c_bk,
bilesenler_yan, kat192, parlaklik, jk_gurultu, IKINCIL_PENCERE, EK_PENCERE; s_W kodu 197c.hesap
AYNEN. Çizgi biçimi 198b.cizgi_bicimi_g (eşit-ΔL geometri, jk grubu açık).
 Pencere başına: s_W (kalibrasyon [−1.30, −0.55]), İKİNCİL + yan [0.15, 2.35] (34 bileşen)
 ⇒ κ_p(W) = A_W(+log p)/(s_W·c_BK(p)), p ∈ {2,3,5,7}, 8-grup jk.
 γ̂: 1 − κ_p(W) = B_p (L_W/L_d)^{−γ}, serbest B_p + ortak γ, ağırlık 1/se² (B_p profillenir,
 γ ızgara + sınırlı 1-B arama). σ_γ = √Σ_W var_jk,W (pencere başına loo, diğerleri tam) × f_γ.
 Λ_alt = Σ_p w_p r_p/Σ w_p, r_p = (1 − κ_p^alt)/(1 − κ_p^d), w_p = 1/var_jk(r_p) (alt + d
 loo'ları); σ_Λ × f_Λ. KAYIT: Λ_son, Λ_son/Λ_alt (çapasız: (1−κ^son)/(1−κ^alt)).
 H-198a (γ), H-198b (Λ_alt; M6 gücü < %80 ise KAYIT). Kapılar: M1, M2, M3, M4, M5', M6, M7, M8.
YORUMLAR (rapor edildi): (a) M5' ayna R̄_W: s_W, {2A_2, 3A_3}'ün ters-varyans ortalaması olduğundan
 ağırlıklı R̄_W ≡ 1; kapı HER BİRİ R_2 = 2A_2/s, R_3 = 3A_3/s ∈ [0.9, 1.1] olarak uygulanır.
 (b) M5' + yan penceresinde "−log2 korelasyonu" → +log2 tepesi (|Δω − log2| ≤ 0.06).
 (c) M7 χ²: d_p = κ_p(yA) − κ_p(yB), χ² = Σ d_p²/var_jk(d_p) (4 sd; köşegen jk), p = χ²_4 kuyruğu.
 (d) M6 gürültüsü: 197 düşük-pencere jk kovaryansı (A_2, A_3, A_5, A_7, s; göreli) × N_d/N_W;
 κ = κ_doğru·(1 + e_A)/(1 + e_s) (s ortak ⇒ pencere içi κ korelasyonu korunur).
 (e) M8 doğrusu: tüm c_BK ≥ 0.02 uyduları [−1.45, 2.45], genlik s·c(r)·Π_{p|pay(r)} κ_p
 (p ≤ 7: 197 κ_p; p ≥ 11: p^{−1/3} — kaptan kabulü); gürültü W_düşük ölçülen dilim jk se ×
 √(N_d/N_W); W_alt kapsaması KISMİ (kaptan kararı; 192 kuralı, B_j ≥ 59/72).
 M8 KARAR NİCELİĞİ (KALEM son hâli): |bias_γ| ≤ 0.25σ_γ; κ_p ve Λ yanlılıkları KAYIT; hüküm
 json'unda pencere başına M8 yanlılığı ve M8-düzeltilmiş κ_p KAYIT olarak yazılır.
Modlar: m1 | m7 | f_kal | m6_on | m8 | hukum   (sıra: m1 → m7 → f_kal → m6_on → m8 → [198a]
        → [198b harita alt/son] → hukum)
"""
import argparse
import hashlib
import importlib.util
import json
import sys
import time
from fractions import Fraction as F
from pathlib import Path

sys.dont_write_bytecode = True
import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = QM / "scratchpad"
S197, S198 = SCR / "197", SCR / "198"
KALEM = QM / "KALEM_KAPPA_YUKSEKLIK_25EYL2026.md"
NJACK, DW, EPS = 8, 0.025, 1e-9
PRIMES = (2, 3, 5, 7)
L_D = 10.48392954102207
PENC = ("alt", "dusuk", "son")
BILESIK = {"+log6": F(6), "+log10": F(10), "+log(3/2)": F(3, 2), "+log(5/2)": F(5, 2),
           "+log(5/4)": F(5, 4), "+log(5/3)": F(5, 3), "-log(3/2)": F(2, 3), "-log(4/3)": F(3, 4)}
K_KAL = {"HS": list(range(1980000, 1980200)), "HC": list(range(1981000, 1981200))}
K_GUC = {"HS": list(range(1982000, 1982200)), "HC": list(range(1983000, 1983200))}
K_M8 = {W: list(range(1984000 + 1000 * i, 1984000 + 1000 * i + 200)) for i, W in enumerate(PENC)}
_t = [set(v) for d in (K_KAL, K_GUC, K_M8) for v in d.values()]
assert all(not (a & b) for i, a in enumerate(_t) for b in _t[i + 1:])


def yukle(ad, yol):
    spec = importlib.util.spec_from_file_location(ad, yol)
    m = importlib.util.module_from_spec(spec)
    sys.modules[ad] = m
    spec.loader.exec_module(m)
    return m


c197 = yukle("c197", QM / "197_configs" / "197c_hukum.py")
b198 = yukle("b198b", QM / "198_configs" / "198b_harita.py")
jk, par, ozet, c_bk = c197.jk, c197.par, c197.ozet, c197.c_bk


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


# ============================ geometri / çizgi biçimi ============================
class BicimG:
    def __init__(self, mid, kb, Lb, kap, J, grup):
        self.a = (mid, kb, Lb, kap, J)
        self.grup = grup
        self.on = {}

    def __call__(self, log_r):
        k = round(float(log_r), 12)
        if k not in self.on:
            self.on[k] = b198.cizgi_bicimi_g(*self.a, float(log_r), self.grup)
        return self.on[k]


def geo_W(W, geo="esitL"):
    """Pencere geometrisi (çekirdek YOK): mid, bloklar, grup, τ'_üst, kapsama, L_W, N."""
    pw = b198.PENCERELER[W]
    _, mid, L = b198.b190.kinematik_eta(pw["zincir"] / pw["eta"])
    kj8 = b198.b188.jk_kenar(len(mid))
    kb, Lb, grup = b198.geometri(mid, kj8, geo)
    jlo, jhi = b198.jaralik()
    J = np.arange(jlo, jhi + 1)
    gerek = float(np.max(Lb) + (jhi + 0.5) * DW) / L
    tau_ust = round(max(b198.TAU_UST_MIN, np.ceil(gerek / 0.005 - 1e-9) * 0.005), 3)
    kap = b198.kapsama(J, L, Lb, tau_ust)
    return {"mid": mid, "kb": kb, "Lb": Lb, "grup": grup, "J": J, "kap": kap, "L": L,
            "N": len(mid), "H": BicimG(mid, kb, Lb, kap, J, grup)}


def harita_oku(W, geo, maske):
    d = np.load(S198 / f"{W}_{geo}" / "harita_omega.npz")
    g = geo_W(W, geo)
    assert np.array_equal(g["kb"], d["kb"]) and np.array_equal(g["kap"], d["kapsama"])
    return {**g, "prof": d[f"kappa_sigma_duz_{maske}"], "prof_ham": d[f"kappa_sigma_ham_{maske}"],
            "L": float(d["L_W"])}


# ============================ pencere analizi (197c AYNEN) ============================
def pencere_analiz(prof, H, J, tam=True):
    kat = c197.kat192()
    kal = c197.uyum_kal(prof, H, J)
    # ---- s_W: 197c.hesap AYNEN ----
    oran = {a: float(a) * par(kal, c197.ad_x(F(a))) for a in (2, 3)}
    se_o = {a: jk(oran[a][1:]) for a in oran}
    w_a = {a: 1.0 / se_o[a] ** 2 for a in oran}
    s = sum(w_a[a] * oran[a] for a in oran) / sum(w_a.values())
    mu0 = [k["r"] for k in kat.values() if k["liste"] == "sonmeli"]
    ik, fik = c197.ikincil(prof, H, J, c197.IKINCIL_PENCERE, mu0, kat, +1)
    kap = {p: par(fik, str(F(p))) / (s * c_bk(F(p))) for p in PRIMES}
    out = {"s": s, "w_a": w_a, "oran": oran, "kal": kal, "fik": fik, "ik": ik, "kappa": kap}
    if not tam:
        return out
    ek, fek = c197.ikincil(prof, H, J, c197.EK_PENCERE, [], kat, -1)
    bil = {}
    for ad, r in BILESIK.items():
        f = fik if r > 1 else fek
        A = par(f, str(r))
        if A is None:
            bil[ad] = "tanımsız"
            continue
        carp = np.ones(NJACK + 1)
        for p in PRIMES:
            if r.numerator % p == 0:
                carp = carp * kap[p]
        kr = A / (s * c_bk(r))
        bil[ad] = {"kappa_r": ozet(kr), "kappa_carpimi": ozet(carp), "oran": ozet(kr / carp)}
    Rx = {x: par(kal, c197.ad_x(F(x))) * float(x) / s for x in (2, 3, F(5, 2), F(7, 2))}
    mu = {str(r): ozet(par(fik, str(r)) / s) for r in mu0 if par(fik, str(r)) is not None}
    # ---- M2 ----
    P = {ad: c197.parlaklik(prof, J, h) for ad, h in c197.ANA6.items()}
    m2y = {ad: bool(P[ad][0] > 3 * P[ad][2]) for ad in P}
    ayna = oran[2] / oran[3]
    M2 = bool(all(m2y.values()) and 0.85 <= ayna[0] <= 1.15)
    # ---- M5' ----
    def m5(fit, pik):
        mk = fit["m"]
        rms = float(np.sqrt(np.mean(fit["artik"] ** 2)))
        std = float(np.std(prof[0, mk]))
        sel = np.abs(J[mk] * DW - pik) <= 0.06 + EPS
        korr = float(np.corrcoef(prof[0, mk][sel], fit["model"][0][sel])[0, 1])
        return {"artik_rms": rms, "std": std, "oran": rms / std, "korr": korr,
                "gecti": bool(rms <= 0.35 * std and korr >= 0.85)}
    m5k, m5i = m5(kal, -np.log(2)), m5(fik, np.log(2))
    R2, R3 = float(Rx[2][0]), float(Rx[3][0])
    M5p = bool(m5k["gecti"] and m5i["gecti"] and 0.9 <= R2 <= 1.1 and 0.9 <= R3 <= 1.1)
    out.update({
        "ek": ek, "bilesik": bil, "R_ayna": {str(x): ozet(v) for x, v in Rx.items()},
        "q_f_bolu_s": ozet(par(kal, "q_f") / s), "mu0_bolu_s": mu,
        "M2": {"gecti": M2, "yanma": {ad: {"P_bolu_se": P[ad][0] / P[ad][2], "yanar": m2y[ad]}
                                      for ad in P},
               "ayna_2A2_bolu_3A3": ozet(ayna)},
        "M5p": {"gecti": M5p, "kalibrasyon": m5k, "arti_yan": m5i, "R_2": R2, "R_3": R3}})
    return out


# ============================ γ̂, Λ, hükümler ============================
_GRID = np.linspace(-4.0, 6.0, 2001)


def gamma_hat(y, w, Lr):
    """1 − κ = B_p (L/L_d)^{−γ}: B_p profillenir (ağırlıklı doğrusal), γ ızgara + sınırlı arama."""
    from scipy.optimize import minimize_scalar

    def chi(g):
        g = np.atleast_1d(g)
        U = Lr[None, :] ** (-g[:, None])                       # (ng, nW)
        num = np.einsum("wp,wp,gw->gp", w, y, U)
        den = np.einsum("wp,gw->gp", w, U ** 2)
        B = num / den
        r = y[None] - U[:, :, None] * B[:, None, :]
        return (w[None] * r ** 2).sum((1, 2)), B
    c, _ = chi(_GRID)
    g0 = _GRID[int(np.argmin(c))]
    res = minimize_scalar(lambda g: chi(g)[0][0], bounds=(g0 - 0.01, g0 + 0.01),
                          method="bounded", options={"xatol": 1e-10})
    return float(res.x), chi(res.x)[1][0]


def gamma_jk(K, Ls):
    """K: {W: (9, 4) κ}, Ls: {W: L_W}. → γ̂, σ_jk (pencere başına loo birleşimi), se_κ."""
    Ws = list(K)
    se = np.array([[jk(K[W][1:, i]) for i in range(4)] for W in Ws])
    w = 1.0 / se ** 2
    Lr = np.array([Ls[W] / L_D for W in Ws])
    y0 = np.array([1 - K[W][0] for W in Ws])
    g0, B0 = gamma_hat(y0, w, Lr)
    var = 0.0
    for a, W in enumerate(Ws):
        gi = []
        for i in range(NJACK):
            y = y0.copy()
            y[a] = 1 - K[W][1 + i]
            gi.append(gamma_hat(y, w, Lr)[0])
        var += (NJACK - 1) / NJACK * np.sum((np.array(gi) - np.mean(gi)) ** 2)
    return g0, float(np.sqrt(var)), B0, se


def lam(Ka, Kb):
    """Λ = Σ w r/Σ w, r_p = (1 − κ_a)/(1 − κ_b); var_jk(r) = a-loo + b-loo; σ_jk birleşik."""
    r0 = (1 - Ka[0]) / (1 - Kb[0])
    ra = (1 - Ka[1:]) / (1 - Kb[0])
    rb = (1 - Ka[0]) / (1 - Kb[1:])
    v = (NJACK - 1) / NJACK * (((ra - ra.mean(0)) ** 2).sum(0) + ((rb - rb.mean(0)) ** 2).sum(0))
    w = 1.0 / v
    L0 = float((w * r0).sum() / w.sum())
    La, Lb_ = (ra * w).sum(1) / w.sum(), (rb * w).sum(1) / w.sum()
    s = float(np.sqrt((NJACK - 1) / NJACK * (((La - La.mean()) ** 2).sum() +
                                             ((Lb_ - Lb_.mean()) ** 2).sum())))
    return L0, s, r0, np.sqrt(v)


def kural(x, s, hS, hC):
    """H-198a/b ortak kuralı (KALEM birebir)."""
    dS, dC = abs(x - hS), abs(x - hC)
    if dS <= 2 * s and dC > 2 * s:
        return "H_S TUTAR"
    if dC <= 2 * s and dS > 2 * s:
        return "H_C TUTAR"
    if dC > 3 * s and dS > 3 * s:
        return "İKİSİ DE ÖLÜR"
    return "belirsiz"


# ============================ M6 (κ_p düzeyi, sentetik) ============================
def m6_gurultu_modeli(maske):
    """197 düşük-pencere jk kovaryansı (göreli; A_2, A_3, A_5, A_7, s)."""
    geo = "197" if maske == "tam" else "esitL"
    h = harita_oku("dusuk", geo, maske)
    a = pencere_analiz(h["prof"], h["H"], h["J"], tam=False)
    X = np.array([par(a["fik"], str(F(p))) for p in PRIMES] + [a["s"]])   # (5, 9)
    d = X[:, 1:] / X[:, :1] - 1
    dm = d - d.mean(1, keepdims=True)
    S = (NJACK - 1) / NJACK * dm @ dm.T
    kd = np.array([a["kappa"][p][0] for p in PRIMES])
    sek = np.array([jk(a["kappa"][p][1:]) for p in PRIMES])
    return S, kd, sek, geo


def pencere_bilgi():
    out = {}
    for W in PENC:
        pw = b198.PENCERELER[W]
        e = np.load(pw["zincir"] / pw["eta"])
        out[W] = {"L": float(e["L"]), "N": int(len(e["mid"]))}
    return out


def sentetik_kappa(hip, kd, S, pb, rng):
    """Bir replika: {W: (9, 4) κ} — jk yapılı gürültü, pencere ölçeği N_d/N_W."""
    g = 1.0 if hip == "HS" else 0.0
    B = 1 - kd
    K = {}
    for W in PENC:
        kt = 1 - B * (pb[W]["L"] / L_D) ** (-g)
        e = rng.multivariate_normal(np.zeros(5), NJACK * S * pb["dusuk"]["N"] / pb[W]["N"],
                                    size=NJACK)
        v = np.vstack([e.mean(0)] + [(e.sum(0) - e[i]) / (NJACK - 1) for i in range(NJACK)])
        K[W] = kt[None, :] * (1 + v[:, :4]) / (1 + v[:, 4:5])
    return K


def m6_replika(hip, kd, S, pb, tohum, f=None):
    K = sentetik_kappa(hip, kd, S, pb, np.random.default_rng(tohum))
    Ls = {W: pb[W]["L"] for W in PENC}
    g, sg, _, _ = gamma_jk(K, Ls)
    La, sa, _, _ = lam(K["alt"], K["dusuk"])
    return g, sg, La, sa


def mod_f_kal(maske):
    yol = S198 / "F_198.json"
    if yol.exists():
        raise SystemExit("f SABİT — yeniden kalibre EDİLMEZ (dosya var)")
    S, kd, sek, geo = m6_gurultu_modeli(maske)
    pb = pencere_bilgi()
    LS = pb["dusuk"]["L"] / pb["alt"]["L"]
    out = {}
    for hip in ("HS", "HC"):
        gt, lt = (1.0, LS) if hip == "HS" else (0.0, 1.0)
        z = np.array([[(g - gt) / sg, (La - lt) / sa]
                      for g, sg, La, sa in (m6_replika(hip, kd, S, pb, t) for t in K_KAL[hip])])
        out[hip] = {"z_gamma_ort": float(z[:, 0].mean()), "z_gamma_sd": float(z[:, 0].std(ddof=1)),
                    "z_Lambda_ort": float(z[:, 1].mean()), "z_Lambda_sd": float(z[:, 1].std(ddof=1)),
                    "n": len(z)}
        print(f"  [{hip}-doğru] z_γ ort {out[hip]['z_gamma_ort']:+.3f} sd "
              f"{out[hip]['z_gamma_sd']:.3f}; z_Λ ort {out[hip]['z_Lambda_ort']:+.3f} sd "
              f"{out[hip]['z_Lambda_sd']:.3f}", flush=True)
    F_ = {"f_gamma": max(out[h]["z_gamma_sd"] for h in out),
          "f_Lambda": max(out[h]["z_Lambda_sd"] for h in out), "dogrular": out,
          "tohumlar": {h: [K_KAL[h][0], K_KAL[h][-1]] for h in K_KAL}, "maske": maske,
          "gurultu_geometri": geo, "kappa_d": kd.tolist(), "se_kappa_d": sek.tolist(),
          "pencere": pb, "kod_sha256": sha(Path(__file__).resolve()),
          "zaman": time.strftime("%Y-%m-%d %H:%M:%S %z")}
    json.dump(F_, open(yol, "w"), indent=1, ensure_ascii=False)
    print(f"  f_γ = {F_['f_gamma']:.4f}   f_Λ = {F_['f_Lambda']:.4f}  -> {yol}", flush=True)


def f_oku():
    d = json.load(open(S198 / "F_198.json"))
    return d["f_gamma"], d["f_Lambda"], d["maske"]


def mod_m6_on():
    fg, fl, maske = f_oku()
    S, kd, sek, geo = m6_gurultu_modeli(maske)
    pb = pencere_bilgi()
    LS = pb["dusuk"]["L"] / pb["alt"]["L"]
    sonuc, gecti = {}, True
    for hip in ("HS", "HC"):
        say_a, say_b, sgs, sls = {}, {}, [], []
        for t in K_GUC[hip]:
            g, sg, La, sa = m6_replika(hip, kd, S, pb, t)
            ka = kural(g, fg * sg, 1.0, 0.0)
            kb = kural(La, fl * sa, LS, 1.0)
            say_a[ka] = say_a.get(ka, 0) + 1
            say_b[kb] = say_b.get(kb, 0) + 1
            sgs.append(fg * sg)
            sls.append(fl * sa)
        n = len(K_GUC[hip])
        oa = {k: v / n for k, v in say_a.items()}
        ob = {k: v / n for k, v in say_b.items()}
        dk, yk = ("H_S TUTAR", "H_C TUTAR") if hip == "HS" else ("H_C TUTAR", "H_S TUTAR")
        ga = oa.get(dk, 0) >= 0.80 and oa.get(yk, 0) <= 0.05
        gecti &= ga
        sonuc[hip] = {"H198a": oa, "H198a_dogru": oa.get(dk, 0), "H198a_yanlis": oa.get(yk, 0),
                      "H198a_gecti": bool(ga), "H198b": ob, "H198b_dogru": ob.get(dk, 0),
                      "H198b_yanlis": ob.get(yk, 0),
                      "sigma_eff_gamma_medyan": float(np.median(sgs)),
                      "sigma_eff_Lambda_medyan": float(np.median(sls)),
                      "M4_gecme": float(np.mean(np.array(sgs) <= 0.35))}
        r = sonuc[hip]
        print(f"  [{hip}-doğru] H-198a doğru {r['H198a_dogru']:.3f} yanlış {r['H198a_yanlis']:.3f} "
              f"{oa}; H-198b doğru {r['H198b_dogru']:.3f} yanlış {r['H198b_yanlis']:.3f}; "
              f"σ_eff(γ) med {r['sigma_eff_gamma_medyan']:.3f}, σ_eff(Λ) med "
              f"{r['sigma_eff_Lambda_medyan']:.4f}; M4 {r['M4_gecme']:.2f}", flush=True)
    guc_b = min(sonuc[h]["H198b_dogru"] for h in sonuc)
    out = {"gecti": bool(gecti), "H198b_guc_min": guc_b, "H198b_yalniz_KAYIT": bool(guc_b < 0.80),
           "f_gamma": fg, "f_Lambda": fl, "maske": maske, "dogrular": sonuc,
           "tohumlar": {h: [K_GUC[h][0], K_GUC[h][-1]] for h in K_GUC},
           "F_198_sha256": sha(S198 / "F_198.json"), "kod_sha256": sha(Path(__file__).resolve()),
           "zaman": time.strftime("%Y-%m-%d %H:%M:%S %z")}
    json.dump(out, open(S198 / "M6_on_198.json", "w"), indent=1, ensure_ascii=False)
    print(f"  M6 (H-198a) → {'GEÇTİ' if gecti else 'KALDI'}; H-198b gücü min {guc_b:.3f} → "
          f"{'yalnız KAYIT' if guc_b < 0.80 else 'hükme girer'}", flush=True)


# ============================ M1 ============================
def mod_m1():
    H197 = json.load(open(S197 / "HUKUM_197.json"))
    h = harita_oku("dusuk", "197", "tam")
    # kod yolu: grup-genelleştirilmiş profil ≡ 197 profili (ortak dilimler)
    P7 = np.load(S197 / "harita_omega32.npz")
    J7 = P7["J"]
    ix = np.searchsorted(J7, h["J"])
    assert np.array_equal(J7[ix], h["J"])
    prof_fark = float(np.max(np.abs(h["prof"] - P7["kappa_sigma_duz"][:, ix])) /
                      np.max(np.abs(P7["kappa_sigma_duz"][:, ix])))
    hb7 = b198.h197.cizgi_bicimi(h["mid"], h["kb"], h["Lb"], h["kap"], h["J"], np.log(2.0))
    cb_fark = float(np.max(np.abs(h["H"](np.log(2.0)) - hb7)))
    a = pencere_analiz(h["prof"], h["H"], h["J"])
    s_ref = H197["nicelikler"]["s"]
    fs = max(abs(ozet(a["s"])["deger"] / s_ref["deger"] - 1), abs(ozet(a["s"])["se"] / s_ref["se"] - 1))
    ref = H197["KAYIT"]["IKINCIL_arti"]["konumlar_192"]
    fA = 0.0
    for ad, v in ref.items():
        if "A" in v:
            m = a["ik"]["konumlar_192"][ad]["A"]
            fA = max(fA, abs(m["deger"] / v["A"]["deger"] - 1), abs(m["se"] / v["A"]["se"] - 1))
    gecti = bool(max(fs, fA, prof_fark) <= 1e-10 and cb_fark == 0.0)
    out = {"kapi": "M1: W_düşük (197 geometrisi) → HUKUM_197 İKİNCİL genlikleri ve s ≤ 1e-10",
           "s_maks_bagil": fs, "IKINCIL_A_maks_bagil": fA, "n_konum": sum("A" in v for v in ref.values()),
           "profil_maks_bagil_vs_197": prof_fark, "cizgi_bicimi_g_vs_197b_maks_fark": cb_fark,
           "kappa_d": {p: ozet(a["kappa"][p]) for p in PRIMES}, "gecti": gecti,
           "kod_sha256": sha(Path(__file__).resolve()),
           "harita_sha256": sha(S198 / "dusuk_197" / "harita_omega.npz"),
           "zaman": time.strftime("%Y-%m-%d %H:%M:%S %z")}
    json.dump(out, open(S198 / "M1_198.json", "w"), indent=1, ensure_ascii=False)
    print(f"  M1: s maks bağıl {fs:.2e}; İKİNCİL A (deger+se, {out['n_konum']} konum) maks bağıl "
          f"{fA:.2e}; profil vs 197 {prof_fark:.2e}; çizgi biçimi_g vs 197b {cb_fark:.1e} → "
          f"{'GEÇTİ' if gecti else 'KALDI'}", flush=True)
    print(f"  κ^d: { {p: (round(v['deger'], 4), round(v['se'], 4)) for p, v in out['kappa_d'].items()} }",
          flush=True)


# ============================ M7 ============================
def mod_m7():
    from scipy.stats import chi2
    out = {}
    K = {}
    for yar in ("yA", "yB"):
        h = harita_oku("dusuk", "esitL", yar)
        K[yar] = np.array([pencere_analiz(h["prof"], h["H"], h["J"], tam=False)["kappa"][p]
                           for p in PRIMES])            # (4, 9)
    d = K["yA"] - K["yB"]
    var = np.array([jk(d[i, 1:]) ** 2 for i in range(4)])
    x2 = float(np.sum(d[:, 0] ** 2 / var))
    p = float(chi2.sf(x2, 4))
    out = {"kapi": "M7: W_düşük (eşit-ΔL) havuz yarıları τ∈[0.45,0.60) / [0.60,0.74); χ²_4 p > 0.01",
           "kappa_yA": {q: ozet(K["yA"][i]) for i, q in enumerate(PRIMES)},
           "kappa_yB": {q: ozet(K["yB"][i]) for i, q in enumerate(PRIMES)},
           "fark": {q: ozet(d[i]) for i, q in enumerate(PRIMES)}, "chi2": x2, "p": p,
           "gecti": bool(p > 0.01),
           "karar": "HAVUZ' τ∈[0.45,0.74)" if p > 0.01 else "ORTAK MUTLAK HAVUZ q∈[224,992]",
           "maske": "tam" if p > 0.01 else "mutlak", "kod_sha256": sha(Path(__file__).resolve()),
           "zaman": time.strftime("%Y-%m-%d %H:%M:%S %z")}
    json.dump(out, open(S198 / "M7_198.json", "w"), indent=1, ensure_ascii=False)
    for q in PRIMES:
        a, b = out["kappa_yA"][q], out["kappa_yB"][q]
        print(f"  κ_{q}: yA {a['deger']:.4f}±{a['se']:.4f}  yB {b['deger']:.4f}±{b['se']:.4f}", flush=True)
    print(f"  M7 χ²_4 = {x2:.2f}, p = {p:.4f} → {'GEÇTİ' if p > 0.01 else 'KALDI'} ({out['karar']})",
          flush=True)


# ============================ M8 (geometri yanlılığı) ============================
def kappa_dogru(p, kd):
    return kd[PRIMES.index(p)] if p in PRIMES else float(p) ** (-1.0 / 3.0)


def m8_dogru_profil(g, S, kd, a0, b0):
    bil = {**c197.bilesenler_yan((-1.45, -0.02), 0.02), **c197.bilesenler_yan((0.02, 2.45), 0.02)}
    t = np.zeros((NJACK + 1, len(g["J"]))) + a0 + b0 * g["J"] * DW
    for r, cc in bil.items():
        k = 1.0
        for p in c197.carpanlar(r.numerator):
            k *= kappa_dogru(p, kd)
        t += S * cc * k * g["H"](np.log(float(r)))
    return t, len(bil)


def mod_m8():
    fg, _, maske = f_oku()
    h = harita_oku("dusuk", "esitL", maske)               # görülmüş: ölçülen dilim se'si
    se_d = np.array([jk(h["prof"][1:, j]) for j in range(len(h["J"]))])
    a_d = pencere_analiz(h["prof"], h["H"], h["J"], tam=False)
    kd = np.array([a_d["kappa"][p][0] for p in PRIMES])
    S = float(a_d["s"][0])
    th = a_d["kal"]["th"][0]
    a0, b0 = float(th[-2]), float(th[-1])
    pb = pencere_bilgi()
    sonuc, Kr = {}, {W: [] for W in PENC}
    for W in PENC:
        g = geo_W(W, "esitL")
        t, nbil = m8_dogru_profil(g, S, kd, a0, b0)
        se_j = se_d * np.sqrt(pb["dusuk"]["N"] / pb[W]["N"])
        a0r = pencere_analiz(t, g["H"], g["J"], tam=False)
        gurultusuz = {p: float(a0r["kappa"][p][0] - kd[i]) for i, p in enumerate(PRIMES)}
        bias = []
        for tohum in K_M8[W]:
            p_ = c197.jk_gurultu(t, se_j, np.random.default_rng(tohum))
            a = pencere_analiz(p_, g["H"], g["J"], tam=False)
            kv = np.array([a["kappa"][p] for p in PRIMES]).T              # (9, 4)
            Kr[W].append(kv)
            bias.append(kv[0] - kd)
        bias = np.array(bias)
        mb, sb = bias.mean(0), bias.std(0, ddof=1) / np.sqrt(len(bias))
        sonuc[W] = {"n_blok": int(len(g["Lb"])), "B_j_min": int(g["kap"].sum(0).min()),
                    "kapsama": "kısmi (192 kuralı)" if g["kap"].sum(0).min() < len(g["Lb"])
                    else "tam", "n_dogru_bileseni": nbil,
                    "KAYIT_gurultusuz_bias": gurultusuz,
                    "KAYIT_bias": {p: float(mb[i]) for i, p in enumerate(PRIMES)},
                    "KAYIT_bias_mc_se": {p: float(sb[i]) for i, p in enumerate(PRIMES)},
                    "KAYIT_bias_0.005_icinde": bool(np.all(np.abs(mb) <= 0.005))}
        print(f"  [{W}] bias(κ_p) = { {p: round(float(mb[i]), 4) for i, p in enumerate(PRIMES)} } "
              f"(MC se ≤ {sb.max():.4f}); gürültüsüz { {p: round(v, 4) for p, v in gurultusuz.items()} }"
              f"; B_j min {sonuc[W]['B_j_min']}/{sonuc[W]['n_blok']}", flush=True)
    Ls = {W: pb[W]["L"] for W in PENC}
    gs, ss = [], []
    for k in range(len(K_M8["alt"])):
        g_, s_, _, _ = gamma_jk({W: Kr[W][k] for W in PENC}, Ls)
        gs.append(g_)
        ss.append(fg * s_)
    gb, sg = float(np.mean(gs)), float(np.median(ss))
    # KARAR NİCELİĞİ (KALEM makine raporu sonrası): |bias_γ| ≤ 0.25σ_γ; κ_p ve Λ yanlılıkları KAYIT
    gecti = bool(abs(gb) <= 0.25 * sg)
    lb = {}
    for ad, (a_, b_) in {"Lambda_alt": ("alt", "dusuk"), "Lambda_son": ("son", "dusuk"),
                         "Lambda_son_bolu_alt": ("son", "alt")}.items():
        v = np.array([lam(Kr[a_][k], Kr[b_][k])[0] - 1.0 for k in range(len(K_M8["alt"]))])
        lb[ad] = {"bias": float(v.mean()), "mc_se": float(v.std(ddof=1) / np.sqrt(len(v)))}
    out = {"kapi": "M8 (KALEM son hâli): H_C-doğru sentetik, pencere geometrisi; KARAR |bias_γ| ≤ "
                   "0.25σ_γ; κ_p ve Λ yanlılıkları KAYIT",
           "pencereler": sonuc, "KAYIT_Lambda_bias": lb,
           "gamma_bias": gb, "gamma_bias_mc_se": float(np.std(gs, ddof=1) / np.sqrt(len(gs))),
           "sigma_eff_gamma_medyan": sg, "gamma_bias_bolu_sigma": gb / sg, "gecti": gecti,
           "maske": maske, "f_gamma": fg, "S": S, "kappa_d": kd.tolist(),
           "kod_sha256": sha(Path(__file__).resolve()), "zaman": time.strftime("%Y-%m-%d %H:%M:%S %z")}
    json.dump(out, open(S198 / "M8_198.json", "w"), indent=1, ensure_ascii=False)
    print(f"  KAYIT Λ yanlılıkları: { {k: round(v['bias'], 4) for k, v in lb.items()} }", flush=True)
    print(f"  γ̂ bias = {gb:+.4f} (MC se {out['gamma_bias_mc_se']:.4f}); σ_eff(γ) medyan {sg:.3f}; "
          f"oran {gb/sg:+.3f} → M8 {'GEÇTİ' if gecti else 'KALDI'}", flush=True)


# ============================ HÜKÜM ============================
def ciz(yol, K, Ls, g0, sg, B0, La, sa, baslik):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    INK, INK2, MUTED = "#0b0b0b", "#52514e", "#a3a29c"
    RENK = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK2,
                         "xtick.color": INK2, "ytick.color": INK2})
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.8), gridspec_kw={"width_ratios": [1.6, 1]})
    Ws = list(K)
    Lx = np.array([Ls[W] for W in Ws])
    xx = np.linspace(Lx.min() - 0.2, Lx.max() + 0.2, 100)
    for i, p in enumerate(PRIMES):
        k = np.array([K[W][0, i] for W in Ws])
        e = np.array([jk(K[W][1:, i]) for W in Ws])
        ax[0].errorbar(Lx, k, yerr=e, fmt="o", color=RENK[i], ms=5, capsize=3, label=f"κ_{p}")
        kd = K["dusuk"][0, i] if "dusuk" in K else k[0]
        ax[0].plot(xx, 1 - (1 - kd) * (xx / L_D) ** (-1.0), color=RENK[i], lw=1, ls="--")
        ax[0].plot(xx, np.full_like(xx, kd), color=RENK[i], lw=1, ls=":")
        ax[0].plot(xx, 1 - B0[i] * (xx / L_D) ** (-g0), color=RENK[i], lw=1.8)
    ax[0].set_xlabel("L_W")
    ax[0].set_ylabel("κ_p = A(+log p)/(s·c_BK)")
    ax[0].set_title(f"γ̂ = {g0:.3f} ± {sg:.3f} (σ_eff)   — düz: uyum, kesikli: H_S, noktalı: H_C")
    ax[0].legend(fontsize=8, frameon=False, ncol=4)
    ax[1].errorbar([0], [La], yerr=[sa], fmt="o", color=INK, capsize=4)
    ax[1].axhline(L_D / Ls["alt"], color=RENK[1], ls="--", label="H_S: L_d/L_alt")
    ax[1].axhline(1.0, color=RENK[2], ls=":", label="H_C: 1")
    ax[1].set_xticks([0])
    ax[1].set_xticklabels(["Λ_alt (kör bacak)"])
    ax[1].legend(fontsize=8, frameon=False)
    ax[1].set_title(f"Λ_alt = {La:.3f} ± {sa:.3f} (σ_eff)")
    fig.suptitle(baslik, fontsize=11, color=INK)
    fig.tight_layout()
    fig.savefig(yol, dpi=150)
    plt.close(fig)


def mod_hukum():
    o = json.load(open(S198 / "ONKAYIT_198.json"))
    if sha(QM / "198_configs" / "198a_onkayit.py") != o["sha256"]:
        raise SystemExit("ON-KAYIT SHA UYUMSUZ")
    if sha(KALEM) != o["kalem_sha256"]:
        raise SystemExit("KALEM ön-kayıttan sonra DEĞİŞMİŞ")
    for rel, hh in o["girdi_sha256"].items():
        if sha(QM / rel) != hh:
            raise SystemExit(f"GİRDİ SHA UYUMSUZ: {rel}")
    fg, fl, maske = o["f"]["f_gamma"], o["f"]["f_Lambda"], o["olcum_parametreleri"]["maske"]
    m6 = json.load(open(S198 / "M6_on_198.json"))
    m8 = json.load(open(S198 / "M8_198.json"))
    Ks, An, Ls = {}, {}, {}
    for W in PENC:
        h = harita_oku(W, o["olcum_parametreleri"]["geometri"], maske)
        a = pencere_analiz(h["prof"], h["H"], h["J"])
        Ks[W] = np.array([a["kappa"][p] for p in PRIMES]).T
        An[W], Ls[W] = a, h["L"]
    gecerli = [W for W in PENC if An[W]["M5p"]["gecti"] and An[W]["M2"]["gecti"]]
    Kd = {W: Ks[W] for W in gecerli}
    g0, sjk, B0, se = gamma_jk(Kd, Ls) if len(gecerli) >= 2 and "dusuk" in gecerli else \
        (np.nan, np.nan, np.full(4, np.nan), None)
    sg = fg * sjk
    La, sla, ra, _ = lam(Ks["alt"], Ks["dusuk"])
    Ls_, sls, _, _ = lam(Ks["son"], Ks["dusuk"])
    Rsa, srsa, _, _ = lam(Ks["son"], Ks["alt"])
    LS = L_D / Ls["alt"]
    M4 = bool(sg <= 0.35)
    ha = kural(g0, sg, 1.0, 0.0) if M4 else "belirsiz (M4)"
    hb = kural(La, fl * sla, LS, 1.0)
    if m6["H198b_yalniz_KAYIT"]:
        hb = f"KAYIT (M6 gücü < 0.80) — hesaplanan: {hb}"
    if "alt" not in gecerli:
        hb = f"karardan çıktı (W_alt M2/M5') — hesaplanan: {hb}"
    son = {"hukum": {"H-198a": ha, "H-198b": hb}, "karara_giren_pencereler": gecerli,
           "gamma": {"deger": g0, "sigma_jk": sjk, "sigma_eff": sg, "f_gamma": fg,
                     "B_p": dict(zip(map(str, PRIMES), np.asarray(B0).tolist()))},
           "Lambda_alt": {"deger": La, "sigma_eff": fl * sla, "H_S": LS, "H_C": 1.0},
           "KAYIT": {"Lambda_son": {"deger": Ls_, "sigma_eff": fl * sls, "H_S": L_D / Ls["son"]},
                     "Lambda_son_bolu_alt": {"deger": Rsa, "sigma_eff": fl * srsa,
                                             "H_S": Ls["alt"] / Ls["son"], "H_C": 1.0},
                     "pencereler": {W: {
                         "L_W": Ls[W], "kappa": {p: ozet(An[W]["kappa"][p]) for p in PRIMES},
                         "p_eksi_1_3": {p: float(p) ** (-1 / 3) for p in PRIMES},
                         "bilesik": An[W]["bilesik"], "R_ayna": An[W]["R_ayna"],
                         "q_f_bolu_s": An[W]["q_f_bolu_s"], "mu0_bolu_s": An[W]["mu0_bolu_s"],
                         "s": ozet(An[W]["s"]),
                         "M8_kappa_bias": m8["pencereler"][W]["KAYIT_bias"],
                         "kappa_M8_duzeltilmis": {p: float(An[W]["kappa"][p][0] -
                                                           m8["pencereler"][W]["KAYIT_bias"][str(p)])
                                                  for p in PRIMES}} for W in PENC},
                     "M8_Lambda_bias": m8["KAYIT_Lambda_bias"],
                     "H198b_etiketi": "KAYIT (M6 gücü < 0.80; KALEM kuralı)"},
           "kapilar": {"M2": {W: An[W]["M2"] for W in PENC}, "M5p": {W: An[W]["M5p"] for W in PENC},
                       "M4": {"gecti": M4, "sigma_eff_gamma": sg, "esik": 0.35}},
           "sha_onkayit": o["sha256"], "kod_sha256": sha(Path(__file__).resolve()),
           "zaman": time.strftime("%Y-%m-%d %H:%M:%S %z")}
    json.dump(son, open(S198 / "HUKUM_198.json", "w"), indent=1, ensure_ascii=False, default=str)
    ciz(S198 / "198_kappa_yukseklik.png", Kd if gecerli else Ks, Ls, g0, sg, B0, La, fl * sla,
        "198 — κ_p'nin yükseklik bağımlılığı (BİRİNCİL K_düz)")
    print(json.dumps(son["hukum"], ensure_ascii=False, indent=1))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("mod", choices=["m1", "m7", "f_kal", "m6_on", "m8", "hukum"])
    ap.add_argument("--maske", choices=["tam", "mutlak"], default=None)
    a = ap.parse_args()
    print("=" * 78)
    print(f"198c / {a.mod.upper()}  {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 78, flush=True)
    if a.mod == "m1":
        mod_m1()
    elif a.mod == "m7":
        mod_m7()
    elif a.mod == "f_kal":
        m = a.maske or json.load(open(S198 / "M7_198.json"))["maske"]
        mod_f_kal(m)
    elif a.mod == "m6_on":
        mod_m6_on()
    elif a.mod == "m8":
        mod_m8()
    else:
        mod_hukum()
