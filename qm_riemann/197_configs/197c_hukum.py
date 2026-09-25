# -*- coding: utf-8 -*-
"""
197c — HÜKÜM: AYNA YASASI (KALEM_AYNA_YASASI_25EYL2026.md AYNEN; "İKİNCİ MAKİNE RAPORU" sürümü)
=============================================================================================
Girdi: scratchpad/197/harita_omega32.npz (197b olcum): κ_Σ^{(v)}(j), v = tam + 8 grup loo;
BİRİNCİL = K_düz (S1b), K_ham KAYIT. Çizgi biçimi h(j; r) veriden (197b.cizgi_bicimi,
parametresiz); her jk replikası v kendi κ_Σ^{(v)} ve h^{(v)}'siyle BAŞTAN uydurulur.
 KALİBRASYON [−1.30, −0.55]: serbest A_k (x = 2, 5/2, 3, 7/2), çeyrek (11/6 … 11/3 + 5/3, 23/6)
   ortak serbest q_f^cal, taban a + bΔω.  s = ters-varyans ağ. ort. {2A_2, 3A_3}.
 KÖR [−2.12, −1.58]: serbest A_k (x = 5, 11/2, 6, 13/2, 7, 15/2, 8), serbest a + bΔω; çeyrek
   (31/6 … 49/6 + 29/6, 14/3, 25/3) SABİT genlikle, iki model:
     M_BK: q_f^cal·0.25/x_i;   M_Rg: q_f^cal·0.25/x_i·g(L + Δω_i)/ḡ_cal,
   ḡ_cal = kalibrasyon çeyrek tepelerinde 0.25/x ağırlıklı g ortalaması (listelenen 10 tepe).
   q_f^cal'ın jk belirsizliği: replika v'de q_f^cal^{(v)} aynı grupla yeniden uydurulur ve kör
   uyuma o replikada taşınır.
 R_x = A_x/(s/x); R̄_Z (ağ. 1/se_x², x = 5..8), σ_Z jk; R̄^g = aynı ağırlıklarla R^g_x
   (R^g_x = g(L − log x)/ḡ, ḡ = s ağırlıklarıyla); β ağırlıklı NLS; ρ; ψ — her model için.
 H-197a öz-tutarlılık (M_BK / M_Rg); H-197b, H-197c M_BK'dan (M_Rg KAYIT).
 Kapılar M2, M4, M5, M6 (+ M1/M3 kayıtları). M6: iki-doğrulu güç sınavı (BK-doğru, R-g-doğru).
 KABUL EDİLEN YORUMLAR (kaptan, 25 Eyl): (2) kenardan taşan ana tepeler modelden düşer
 (KALEM'e işlendi); (3) kapı önceliği M1/M3 > M2 > M6 > M5 > M4; (4) H-197b'de ÖLÜM önce,
 "üç hedeften ikisinde" = en az iki; (5) "q_f/R̄_Z" = (q_f^cal/s)/R̄_Z^BK (BK: 1);
 (6) − yan ikincil pencere [−0.55, −0.15] KALEM DIŞI EK, KAYIT.
 Bu betiğin yorumları: M4 = σ^BK ≤ 0.10 VE σ^Rg ≤ 0.10; R̄^g her modelin kendi 1/se_x²
 ağırlıklarıyla (R-g öz-tutarlılığında M_Rg ağırlıkları).
 ÜÇÜNCÜ MAKİNE RAPORU (kaptan): karar σ'ları σ_eff = f·σ_jk — H-197a (öz-tutarlılık 2σ/3σ) ve
 M4 f_Z ile, H-197b σ_ρ f_ρ ile; β karardan çıktı (KAYIT). f_Z, f_ρ: f_kal modunda K_KAL
 tohumlarıyla (BK-doğru M_BK, R-g-doğru M_Rg; iki doğrunun z sd'lerinin BÜYÜĞÜ) → sabit olarak
 F_KALIBRASYON.json → ONKAYIT_197.json; hüküm ön-kayıttaki f'leri kullanır. M6 (iii) güç
 replikaları K_GUC tohumlarıyla (K_KAL ile kesişimsiz). Sentetik provada M2 ATLANIR (loglanır).
Modlar (sıra: f_kal → m6_on → [197a] → [197b olcum] → hukum):
  f_kal     ölçüm-ÖNCESİ f kalibrasyonu (M1 gürültüsü, 192 sinyali) → F_KALIBRASYON.json.
  m6_on     ölçüm-ÖNCESİ M6 (σ_eff ile; güç tohumları ayrı) → M6_on_deneme.json.
  sentetik  uçtan uca prova (BK-doğru, tüm c ≥ 0.02 uyduları) → scratchpad/197/sentetik_test/.
  hukum     GERÇEK (ölçümden sonra): HUKUM_197.json + 197_ayna_yasasi.png (scratchpad/197).
Kullanım: 197c_hukum.py <f_kal|m6_on|sentetik|hukum>
"""
import argparse
import hashlib
import importlib.util
import json
import sys
import time
import warnings
from fractions import Fraction as F
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = QM / "scratchpad"
S190, S192, S197 = SCR / "190", SCR / "192", SCR / "197"
ZD = S190 / "zincir_dusuk"
KALEM = QM / "KALEM_AYNA_YASASI_25EYL2026.md"
L_G = 10.48392954102207          # KALEM: g(ω') için L
DW = 0.025
EPS = 1e-9
NJACK = 8
T7 = 2.36                        # M6 (ii): t₇, %95
AGIRLIK_ESIK = 1e-3              # pencerede Σh < bu ise bileşen tanımsız (sütun düşer)
MODELLER = ("BK", "Rg")


def yukle(ad, yol):
    spec = importlib.util.spec_from_file_location(ad, yol)
    m = importlib.util.module_from_spec(spec)
    sys.modules[ad] = m
    spec.loader.exec_module(m)
    return m


h197 = yukle("h197b", QM / "197_configs" / "197b_harita.py")
jk, cizgi_bicimi, parlaklik = h197.jk, h197.cizgi_bicimi, h197.parlaklik


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


# ============================ BK cebiri (KALEM 3-4) ============================
def carpanlar(n):
    out, p = {}, 2
    while p * p <= n:
        while n % p == 0:
            out[p] = out.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def c_bk(r):
    """c(r) = Π_p f_p(k_p); f_p(k<0) = p^k, f_p(1) = p/(p−1)², f_p(k≥2) = 0 (r sade kesir)."""
    r = F(r)
    c = 1.0
    for p, k in carpanlar(r.numerator).items():
        if k >= 2:
            return 0.0
        c *= p / (p - 1) ** 2
    for p, k in carpanlar(r.denominator).items():
        c *= float(p) ** (-k)
    return c


def g_rg(om):
    return np.exp(om / 2.0) / om


# ============================ KALEM listeleri ============================
def fx(*a):
    return [F(x) if not isinstance(x, tuple) else F(*x) for x in a]


KOR = {"pencere": (-2.12, -1.58),
       "ana": fx(5, (11, 2), 6, (13, 2), 7, (15, 2), 8),
       "ceyrek": fx((31, 6), (16, 3), (17, 3), (35, 6), (37, 6), (19, 3), (20, 3), (41, 6),
                    (43, 6), (22, 3), (23, 3), (47, 6), (49, 6)) + fx((29, 6), (14, 3), (25, 3))}
SIKINTI_DOGRU = fx((9, 2), (17, 2), 9)   # KALEM: modelden düşer; yalnız sentetik DOĞRUDA
KAL = {"pencere": (-1.30, -0.55),
       "ana": fx(2, (5, 2), 3, (7, 2)),
       "ceyrek": fx((11, 6), (13, 6), (7, 3), (8, 3), (17, 6), (19, 6), (10, 3), (11, 3))
       + fx((5, 3), (23, 6))}
Z_AILE = fx(5, 6, 7, 8)
H_AILE = fx((11, 2), (13, 2), (15, 2))
ANA6 = {"+log2": np.log(2), "+log3": np.log(3), "+log6": np.log(6), "+log(3/2)": np.log(1.5),
        "-log2": -np.log(2), "-log3": -np.log(3)}
M2_REF = {("+log3", "+log2"): 0.310, ("+log6", "+log2"): 0.602, ("-log2", "+log2"): 0.240,
          ("-log3", "+log3"): 0.466}
IKINCIL_PENCERE = (0.15, 2.35)
EK_PENCERE = (-0.55, -0.15)      # KALEM DIŞI EK (kabul: KAYIT) — −log(3/2), −log(4/3)
for _x in KOR["ana"] + KAL["ana"]:
    assert abs(c_bk(1 / _x) - 1 / float(_x)) < 1e-12, _x
for _x in KOR["ceyrek"] + KAL["ceyrek"]:
    assert abs(c_bk(1 / _x) - 0.25 / float(_x)) < 1e-12, _x
# ḡ_cal: kalibrasyon çeyrek tepelerinde 0.25/x ağırlıklı g(L + Δω_i), Δω_i = −log x_i
G_CAL = float(sum(0.25 / float(x) * g_rg(L_G - np.log(float(x))) for x in KAL["ceyrek"]) /
              sum(0.25 / float(x) for x in KAL["ceyrek"]))


def k_model(x, model):
    """Kör çeyrek genlik çarpanı: M_BK 1; M_Rg g(L + Δω_i)/ḡ_cal."""
    return 1.0 if model == "BK" else float(g_rg(L_G - np.log(float(x))) / G_CAL)


def ad_x(x):
    return f"x={x}"


def kat192():
    K = json.load(open(S192 / "ONKAYIT_192.json"))["B"]["katalog"]
    return {ad: {"r": F(int(v["a"]), int(v["b"])), "delta_omega": float(v["delta_omega"]),
                 "liste": v["liste"]} for ad, v in K.items()}


def s192():
    """192 sinyal ölçeği: ters-varyans ağ. ort. {2·P(−log2), 3·P(−log3)} (HUKUM_192 A2_P;
    görülmüş bant, zaten kayıtlı)."""
    A = json.load(open(S192 / "HUKUM_192.json"))["A2_P"]
    v = [(2 * A["-log2"]["P"], 2 * A["-log2"]["se"]), (3 * A["-log3"]["P"], 3 * A["-log3"]["se"])]
    w = [1 / e ** 2 for _, e in v]
    return float(sum(wi * p for wi, (p, _) in zip(w, v)) / sum(w))


def bilesenler_yan(pencere, esik=0.05, mu0=()):
    """Pencere içindeki (merkez ± yarım dilim) r = a/b uyduları: BK c ≥ esik + verilen μ=0."""
    lo, hi = np.exp(pencere[0] - 0.5 * DW), np.exp(pencere[1] + 0.5 * DW)
    out = {}
    for b in range(1, 41):
        for a in range(1, int(hi * b) + 2):
            r = F(a, b)
            if r.denominator != b or not (lo <= float(r) <= hi) or r in out:
                continue
            c = c_bk(r)
            if c >= esik:
                out[r] = c
    for r in mu0:
        if lo <= float(r) <= hi:
            out.setdefault(r, c_bk(r))
    return dict(sorted(out.items(), key=lambda t: float(t[0])))


# ============================ çizgi biçimi önbelleği ============================
class Bicim:
    def __init__(self, mid, kb, Lb, kap, J):
        self.a = (mid, kb, Lb, kap, J)
        self.on = {}

    def __call__(self, log_r):
        k = round(float(log_r), 12)
        if k not in self.on:
            self.on[k] = cizgi_bicimi(*self.a, float(log_r))
        return self.on[k]


# ============================ profil uyumu ============================
def uyum(prof, H, J, pencere, serbest, ortak=(), ortak_c=()):
    """KALEM PROFİL UYUMU (doğrusal EKK, pencere = dilim MERKEZİ ∈ [lo, hi]):
    κ_Σ(j) = Σ_k A_k h(j; r_k) [+ q_f Σ_i c_i h(j; r_i)] + a + bΔω_j; her replika kendi h[v]'siyle.
    Pencerede Σh(tam) < AGIRLIK_ESIK olan serbest bileşen TANIMSIZ (sütun düşer)."""
    m = (J * DW >= pencere[0] - EPS) & (J * DW <= pencere[1] + EPS)
    assert np.all(np.isfinite(prof[:, m])), "pencerede tanımsız profil dilimi"
    adlar, fn, dusen = [], [], []
    for ad, lr in serbest:
        if H(lr)[0, m].sum() < AGIRLIK_ESIK:
            dusen.append(ad)
            continue
        adlar.append(ad)
        fn.append(lambda v, lr=lr: H(lr)[v, m])
    if len(ortak):
        adlar.append("q_f")
        fn.append(lambda v: sum(ci * H(lr)[v, m] for lr, ci in zip(ortak, ortak_c)))
    adlar += ["a", "b"]
    nv = prof.shape[0]
    th = np.zeros((nv, len(adlar)))
    model = np.zeros((nv, int(m.sum())))
    Xs = []
    for v in range(nv):
        X = np.stack([f(v) for f in fn] + [np.ones(int(m.sum())), J[m] * DW], 1)
        t, *_ = np.linalg.lstsq(X, prof[v, m], rcond=None)
        th[v], model[v] = t, X @ t
        Xs.append(X)
    sv = np.linalg.svd(Xs[0], compute_uv=False)
    return {"adlar": adlar, "th": th, "m": m, "model": model, "X": Xs,
            "artik": prof[0, m] - model[0], "dusen": dusen,
            "rank": int((sv > sv[0] * 1e-10).sum()), "kosul": float(sv[0] / sv[-1])}


def par(fit, ad):
    return fit["th"][:, fit["adlar"].index(ad)] if ad in fit["adlar"] else None


def uyum_kal(prof, H, J):
    return uyum(prof, H, J, KAL["pencere"],
                [(ad_x(x), -np.log(float(x))) for x in KAL["ana"]],
                [-np.log(float(x)) for x in KAL["ceyrek"]],
                [0.25 / float(x) for x in KAL["ceyrek"]])


def uyum_kor(prof, H, J, qf_cal, model):
    """KÖR: serbest ana A_k + a + bΔω; çeyrek SABİT = q_f^cal^{(v)}·Σ_i 0.25/x_i·k_i·h_i^{(v)}."""
    sabit = np.array([qf_cal[v] * sum(0.25 / float(x) * k_model(x, model) *
                                      H(-np.log(float(x)))[v] for x in KOR["ceyrek"])
                      for v in range(prof.shape[0])])
    f = uyum(prof - sabit, H, J, KOR["pencere"],
             [(ad_x(x), -np.log(float(x))) for x in KOR["ana"]])
    f["sabit"] = sabit[:, f["m"]]
    f["model_toplam"] = f["model"] + f["sabit"]
    f["artik"] = prof[0, f["m"]] - f["model_toplam"][0]
    return f


def nls(xs, A, se, p0=None):
    """A_x = C·x^{−β}, ağırlıklı doğrusal-olmayan EKK (ağırlık 1/se²)."""
    from scipy.optimize import least_squares
    x = np.array([float(v) for v in xs])
    if p0 is None:
        if np.all(A > 0):
            sl, ic = np.polyfit(np.log(x), np.log(A), 1)
            p0 = [np.exp(ic), -sl]
        else:
            p0 = [float(np.mean(A * x)), 1.0]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = least_squares(lambda p: (p[0] * x ** (-p[1]) - A) / se, p0, method="lm")
    return r.x


def ozet(v):
    v = np.asarray(v, float)
    return {"deger": float(v[0]), "se": jk(v[1:])}


# ============================ türetilmiş nicelikler ============================
def turet(kor, s, w_a):
    A = {x: par(kor, ad_x(x)) for x in KOR["ana"]}
    R = {x: A[x] * float(x) / s for x in KOR["ana"]}
    seR = {x: jk(R[x][1:]) for x in R}
    wZ = {x: 1.0 / seR[x] ** 2 for x in Z_AILE}
    RZ = sum(wZ[x] * R[x] for x in Z_AILE) / sum(wZ.values())
    sZ = jk(RZ[1:])
    gbar = sum(w_a[a] * g_rg(L_G - np.log(a)) for a in w_a) / sum(w_a.values())
    Rg = {x: float(g_rg(L_G - np.log(float(x))) / gbar) for x in KOR["ana"]}
    RgZ = float(sum(wZ[x] * Rg[x] for x in Z_AILE) / sum(wZ.values()))
    seA = {x: jk(A[x][1:]) for x in KOR["ana"]}
    AZ = np.array([A[x] for x in Z_AILE])
    sez = np.array([seA[x] for x in Z_AILE])
    p_t = nls(Z_AILE, AZ[:, 0], sez)
    beta = np.array([p_t[1]] + [nls(Z_AILE, AZ[:, v], sez, p_t)[1] for v in range(1, NJACK + 1)])
    wH = {x: 1.0 / seR[x] ** 2 for x in H_AILE}
    RH = sum(wH[x] * R[x] for x in H_AILE) / sum(wH.values())
    rho = RH / RZ
    RgH = float(sum(wH[x] * Rg[x] for x in H_AILE) / sum(wH.values()))
    nneg = int(sum(A[x][0] < 0 for x in H_AILE))
    x3 = Z_AILE[:3]
    A3 = np.array([A[x] for x in x3])
    se3 = np.array([seA[x] for x in x3])
    q_t = nls(x3, A3[:, 0], se3)
    q_v = [q_t] + [nls(x3, A3[:, v], se3, q_t) for v in range(1, NJACK + 1)]
    psi = np.array([A[F(8)][v] / (q_v[v][0] * 8.0 ** (-q_v[v][1])) for v in range(NJACK + 1)])
    chi = {ad: {"BK": float(sum((R[x][0] - 1) ** 2 / seR[x] ** 2 for x in xs)),
                "Rg": float(sum((R[x][0] - Rg[x]) ** 2 / seR[x] ** 2 for x in xs)), "n": len(xs)}
           for ad, xs in (("Z", Z_AILE), ("ana7", KOR["ana"]))}
    return {"A": {ad_x(x): {**ozet(A[x]), "c_BK": 1 / float(x)} for x in KOR["ana"]},
            "R": {ad_x(x): {**ozet(R[x]), "R_g": Rg[x]} for x in R},
            "R_Z": {"deger": float(RZ[0]), "sigma_Z": sZ,
                    "agirlik": {ad_x(x): wZ[x] for x in wZ}},
            "R_g_Z": RgZ, "g_bar": float(gbar), "beta": ozet(beta), "C_Z": float(p_t[0]),
            "rho": ozet(rho), "R_H": ozet(RH), "H_negatif_A": nneg, "psi": ozet(psi),
            "R_g_H": RgH, "rho_g": RgH / RgZ,
            "A8_hat": float(q_t[0] * 8.0 ** (-q_t[1])), "KAYIT_chi2": chi,
            "uyum": {"parametre": {a: ozet(kor["th"][:, i]) for i, a in enumerate(kor["adlar"])},
                     "rank": kor["rank"], "kosul": kor["kosul"], "n_dilim": int(kor["m"].sum())},
            "_RZ": RZ}


def hesap(prof, H, J, tam=True, f=None):
    kal = uyum_kal(prof, H, J)
    qf = par(kal, "q_f")
    oran = {a: float(a) * par(kal, ad_x(F(a))) for a in (2, 3)}
    se_o = {a: jk(oran[a][1:]) for a in oran}
    w_a = {a: 1.0 / se_o[a] ** 2 for a in oran}
    s = sum(w_a[a] * oran[a] for a in oran) / sum(w_a.values())
    fits = {mo: uyum_kor(prof, H, J, qf, mo) for mo in MODELLER}
    out = {"s": ozet(s), "s_agirlik": {str(a): w_a[a] for a in w_a},
           "s_oranlari": {f"A{a}·{a}": ozet(oran[a]) for a in oran},
           "q_f_cal": ozet(qf), "g_cal": G_CAL,
           "kalibrasyon_uyum": {"parametre": {a: ozet(kal["th"][:, i])
                                              for i, a in enumerate(kal["adlar"])},
                                "rank": kal["rank"], "kosul": kal["kosul"],
                                "n_dilim": int(kal["m"].sum())},
           "modeller": {mo: turet(fits[mo], s, w_a) for mo in MODELLER},
           "_fit": {"kal": kal, **fits}, "_s": s}
    RZb = out["modeller"]["BK"]["_RZ"]
    out["KAYIT_qf"] = {"q_f_cal_bolu_s": ozet(qf / s), "q_f_cal_bolu_s_bolu_RZ_BK":
                       ozet(qf / s / RZb)}
    fZ = 1.0 if f is None else f["Z"]
    sB, sR = out["modeller"]["BK"]["R_Z"]["sigma_Z"], out["modeller"]["Rg"]["R_Z"]["sigma_Z"]
    # M4 (KALEM): σ_eff = f_Z·σ_jk ≤ 0.10, M_BK VE M_Rg
    out["M4"] = {"gecti": bool(fZ * sB <= 0.10 and fZ * sR <= 0.10), "sigma_BK": sB,
                 "sigma_Rg": sR, "f_Z": fZ, "sigma_eff_BK": fZ * sB, "sigma_eff_Rg": fZ * sR,
                 "esik": 0.10}
    if not tam:
        return out
    # ---- M2 ----
    P = {ad: parlaklik(prof, J, h) for ad, h in ANA6.items()}
    m2_yan = {ad: {"P": P[ad][0], "se": P[ad][2], "P_bolu_se": P[ad][0] / P[ad][2],
                   "yanar": bool(P[ad][0] > 3 * P[ad][2]), "dilim_j": P[ad][3].tolist()}
              for ad in P}
    m2_or = {}
    for (u, d), ref in M2_REF.items():
        r = P[u][0] / P[d][0]
        m2_or[f"{u}/{d}"] = {"oran": r, "se": jk(P[u][1] / P[d][1]), "ref_192": ref,
                             "goreli_sapma": r / ref - 1, "icinde": bool(abs(r / ref - 1) <= 0.30)}
    M2 = bool(all(v["yanar"] for v in m2_yan.values()) and all(v["icinde"] for v in m2_or.values()))
    # ---- M5 ----
    mk = kal["m"]
    rms = float(np.sqrt(np.mean(kal["artik"] ** 2)))
    std = float(np.std(prof[0, mk]))
    jw = J[mk]
    sel = np.abs(jw * DW + np.log(2)) <= 0.06 + EPS
    korr = float(np.corrcoef(prof[0, mk][sel], kal["model"][0][sel])[0, 1])
    out["kapilar"] = {
        "M2": {"gecti": M2, "yanma": m2_yan, "oranlar": m2_or},
        "M4": out["M4"],
        "M5": {"gecti": bool(rms <= 0.35 * std and korr >= 0.85), "artik_rms": rms,
               "pencere_std": std, "oran": rms / std, "esik_oran": 0.35,
               "mlog2_korelasyon": korr, "esik_korr": 0.85, "mlog2_dilim_j": jw[sel].tolist()}}
    return out


# ============================ hüküm kuralları (KALEM metni birebir) ============================
def hukum_a(mB, mR, fZ):
    """H-197a öz-tutarlılık, σ_eff = f_Z·σ_jk (ÜÇÜNCÜ MAKİNE RAPORU): BK öz-tutarlı ⇔
    |R̄_Z^BK − 1| ≤ 2σ_eff^BK; R-g öz-tutarlı ⇔ |R̄_Z^Rg − R̄^g| ≤ 2σ_eff^Rg (R̄^g M_Rg
    ağırlıklarıyla). β karardan ÇIKTI (KAYIT)."""
    eB = abs(mB["R_Z"]["deger"] - 1)
    eR = abs(mR["R_Z"]["deger"] - mR["R_g_Z"])
    sB, sR = fZ * mB["R_Z"]["sigma_Z"], fZ * mR["R_Z"]["sigma_Z"]
    bk, rg = eB <= 2 * sB, eR <= 2 * sR
    if bk and not rg:
        return "BK TUTAR"
    if rg and not bk:
        return "R-g KAZANIR (BK'nın çıplak eşlemesi ÖLÜR)"
    if eB > 3 * sB and eR > 3 * sR:
        return "İKİSİ DE ÖLÜR"
    return "belirsiz (KAYIT)"


def hukum_b(m, frho):
    """H-197b; σ_ρ = f_ρ·σ_jk (ÜÇÜNCÜ MAKİNE RAPORU)."""
    rho, s_rho, nneg = m["rho"]["deger"], frho * m["rho"]["se"], m["H_negatif_A"]
    if rho < 0.50 or nneg >= 2:          # (4) ÖLÜM önce; "üç hedeften ikisinde" = en az iki
        return "ÖLDÜ"
    if abs(rho - 1) <= max(0.25, 2 * s_rho):
        return "TUTAR"
    return "belirsiz (KAYIT)"


def hukum_c(m):
    psi = m["psi"]["deger"]
    if psi < 0.30:
        return "ÖLDÜ"
    if 0.60 <= psi <= 1.60:
        return "TUTAR"
    return "belirsiz (KAYIT)"


def hukumler(o, m6_gecti, f, m1m3=True, m2_atla=False):
    mB, mR = o["modeller"]["BK"], o["modeller"]["Rg"]
    ham = {"H-197a": hukum_a(mB, mR, f["Z"]), "H-197b": hukum_b(mB, f["rho"]),
           "H-197c": hukum_c(mB)}
    k = o["kapilar"]
    # (3) kapı önceliği: M1/M3 > M2 > M6 > M5 > M4. M2 gerçek-veri makine kapısıdır;
    # sentetik provada ATLANIR (m2_atla, loglanır — BK-doğru sentetikte kalması beklenir).
    if not m1m3:
        etiket = "GEÇERSİZ (M1/M3)"
    elif not k["M2"]["gecti"] and not m2_atla:
        etiket = "SINAV GEÇERSİZ (M2; BK ölümü değil)"
    elif not m6_gecti:
        etiket = "yalnız KAYIT (hüküm-anı M6 kaldı)"
    elif not k["M5"]["gecti"]:
        etiket = "yalnız KAYIT (M5: çizgi-biçimi modeli GEÇERSİZ)"
    elif not k["M4"]["gecti"]:
        etiket = "belirsiz (M4: güç yetersiz)"
    else:
        etiket = None
    son = {h: (v if etiket is None else f"{etiket} — hesaplanan: {v}") for h, v in ham.items()}
    return {"hukum": son, "kural_ciktisi": ham, "kapi_etiketi": etiket,
            "KAYIT_Rg_modeli": {"H-197b": hukum_b(mR, f["rho"]), "H-197c": hukum_c(mR)},
            "KAYIT_beta": {"BK": mB["beta"], "Rg": mR["beta"],
                           "BK_0.60-1.30_icinde": bool(0.60 <= mB["beta"]["deger"] <= 1.30)},
            "M2_atlandi_sentetik": bool(m2_atla), "f": f}


# ============================ M6 (iki-doğrulu güç sınavı) ============================
def jk_gurultu(t, se_j, rng):
    """e_g ~ N(0, 8·se_j²) (8 grup); tam = t_0 + ort_g e_g; replika i = t_i + ort_{g≠i} e_g."""
    e = rng.normal(0.0, np.sqrt(NJACK) * se_j, (NJACK, t.shape[1]))
    p = np.empty_like(t)
    p[0] = t[0] + e.mean(0)
    for i in range(NJACK):
        p[1 + i] = t[1 + i] + (e.sum(0) - e[i]) / (NJACK - 1)
    return p


def dogru_profil(H, J, S, a0, b0, dogru):
    """Sentetik doğru: genlik = S·c(r)·k(Δω); BK: k = 1; R-g: k = g(L + Δω)/ḡ_s,
    ḡ_s = (g(L − log2) + g(L − log3))/2 (s ölçeği iki doğruda da ≈ S)."""
    gs = 0.5 * (g_rg(L_G - np.log(2)) + g_rg(L_G - np.log(3)))
    kd = (lambda x: 1.0) if dogru == "BK" else (lambda x: float(g_rg(L_G - np.log(float(x))) / gs))
    t = np.zeros((NJACK + 1, len(J))) + a0 + b0 * J * DW
    th = {"kor": {}, "kal": {}}
    for x in KOR["ana"] + SIKINTI_DOGRU:
        t += S * kd(x) / float(x) * H(-np.log(float(x)))
        if x in KOR["ana"]:
            th["kor"][ad_x(x)] = S * kd(x) / float(x)
    for x in KAL["ana"]:
        t += S * kd(x) / float(x) * H(-np.log(float(x)))
        th["kal"][ad_x(x)] = S * kd(x) / float(x)
    for x in KOR["ceyrek"] + KAL["ceyrek"]:
        t += S * kd(x) * 0.25 / float(x) * H(-np.log(float(x)))
    # q_f "doğrusu": BK S; R-g'de kalibrasyon modeli (sabit 0.25/x oranları) yanlış belirli —
    # tanım: S·(0.25/x ağırlıklı k ortalaması) = S·ḡ_cal/ḡ_s
    th["kal"]["q_f"] = S * sum(0.25 / float(x) * kd(x) for x in KAL["ceyrek"]) / \
        sum(0.25 / float(x) for x in KAL["ceyrek"])
    for w in ("kor", "kal"):
        th[w]["a"], th[w]["b"] = a0, b0
    return t, th


# f kalibrasyonu ve güç replikaları AYRI tohum kümeleriyle (kesişimsiz; replika başına tohum)
K_KAL = {"BK": list(range(1970000, 1970200)), "Rg": list(range(1971000, 1971200))}
K_GUC = {"BK": list(range(1972000, 1972200)), "Rg": list(range(1973000, 1973200))}
assert not (set(K_KAL["BK"]) | set(K_KAL["Rg"])) & (set(K_GUC["BK"]) | set(K_GUC["Rg"]))


def f_kalibrasyon(H, J, se_j, S, a0, b0):
    """ÜÇÜNCÜ MAKİNE RAPORU: z = (R̄_Z − doğru)/σ_jk — BK-doğruda M_BK (doğru 1), R-g-doğruda
    M_Rg (doğru R̄^g, replikanın M_Rg ağırlıklarıyla); z_ρ aynı eşleşmeyle (BK: ρ doğru 1;
    R-g: M_Rg ρ, doğru ρ^g = R̄^g_H/R̄^g_Z). f = iki doğrunun z standart sapmalarının BÜYÜĞÜ
    (ddof=1), f_Z ve f_ρ ayrı. Tohumlar K_KAL (replika başına)."""
    out = {}
    for dogru in MODELLER:
        t, _ = dogru_profil(H, J, S, a0, b0, dogru)
        zZ, zr = [], []
        for tohum in K_KAL[dogru]:
            o = hesap(jk_gurultu(t, se_j, np.random.default_rng(tohum)), H, J, tam=False)
            m = o["modeller"][dogru]
            hZ = 1.0 if dogru == "BK" else m["R_g_Z"]
            hr = 1.0 if dogru == "BK" else m["rho_g"]
            zZ.append((m["R_Z"]["deger"] - hZ) / m["R_Z"]["sigma_Z"])
            zr.append((m["rho"]["deger"] - hr) / m["rho"]["se"])
        zZ, zr = np.array(zZ), np.array(zr)
        out[dogru] = {"model": f"M_{dogru}", "n": len(zZ),
                      "z_Z_ort": float(zZ.mean()), "z_Z_sd": float(zZ.std(ddof=1)),
                      "z_rho_ort": float(zr.mean()), "z_rho_sd": float(zr.std(ddof=1))}
    return {"f_Z": max(out[d]["z_Z_sd"] for d in MODELLER),
            "f_rho": max(out[d]["z_rho_sd"] for d in MODELLER),
            "dogrular": out, "tohumlar": {d: [K_KAL[d][0], K_KAL[d][-1]] for d in MODELLER},
            "S": S, "a0": a0, "b0": b0, "gurultu_se_medyan": float(np.nanmedian(se_j))}


def m6(H, J, se_j, S, a0, b0, f, tohumlar=None):
    """KALEM M6: (i) gürültüsüz A_k, q_f %5; (ii) KARARA GİREN niceliklerde kapsama ≥ %88 —
    R̄_Z ve ρ |Â − A| ≤ 2σ_eff (f_Z, f_ρ; doğru-eşleşik model: BK-doğru M_BK/1, R-g-doğru
    M_Rg/R̄^g, ρ^g), kalibrasyon çapaları A(−log2), A(−log3), q_f^cal ve kör ana A_x
    (doğru-eşleşik model) |Â − A| ≤ 2.36σ_jk; (iii) doğru hipotez
    kazanma ≥ %80, YANLIŞ ≤ %5 (her iki doğru), kazanma = H-197a σ_eff = f_Z·σ_jk ile (kapılar
    uygulanmadan). Güç replikaları K_GUC tohumlarıyla (f kalibrasyonundan AYRI)."""
    tohumlar = K_GUC if tohumlar is None else tohumlar
    sonuc, gecti = {}, True
    for dogru in MODELLER:
        nrep = len(tohumlar[dogru])
        t, th = dogru_profil(H, J, S, a0, b0, dogru)
        o0 = hesap(t, H, J, tam=False, f=f)
        f_kor, f_kal = o0["_fit"][dogru], o0["_fit"]["kal"]
        gh = {}
        for w, ff in (("kor", f_kor), ("kal", f_kal)):
            for i, a in enumerate(ff["adlar"]):
                if a in th[w] and a not in ("a", "b"):
                    gh[f"{w}:{a}"] = float(abs(ff["th"][0, i] / th[w][a] - 1))
        kap = {f"{w}:{a}": 0 for w in th for a in th[w]}
        say = {"BK TUTAR": 0, "R-g KAZANIR": 0, "İKİSİ DE ÖLÜR": 0, "belirsiz": 0}
        sB, sR, m4, tani = [], [], 0, []
        kk = {"R_Z (2σ_eff)": 0, "rho (2σ_eff)": 0, "kal:x=2": 0, "kal:x=3": 0, "kal:q_f": 0,
              **{f"kor:{ad_x(x)}": 0 for x in KOR["ana"]}}
        for tohum in tohumlar[dogru]:
            o = hesap(jk_gurultu(t, se_j, np.random.default_rng(tohum)), H, J, tam=False, f=f)
            for w, ff in (("kor", o["_fit"][dogru]), ("kal", o["_fit"]["kal"])):
                for i, a in enumerate(ff["adlar"]):
                    if a in th[w]:
                        kap[f"{w}:{a}"] += int(abs(ff["th"][0, i] - th[w][a]) <=
                                               T7 * jk(ff["th"][1:, i]))
            mB_, mR_ = o["modeller"]["BK"], o["modeller"]["Rg"]
            # (ii) KARARA GİREN nicelikler (DÖRDÜNCÜ MAKİNE RAPORU): doğru-eşleşik model
            md = o["modeller"][dogru]
            hZ = 1.0 if dogru == "BK" else md["R_g_Z"]
            hr = 1.0 if dogru == "BK" else md["rho_g"]
            kk["R_Z (2σ_eff)"] += int(abs(md["R_Z"]["deger"] - hZ) <= 2 * f["Z"] * md["R_Z"]["sigma_Z"])
            kk["rho (2σ_eff)"] += int(abs(md["rho"]["deger"] - hr) <= 2 * f["rho"] * md["rho"]["se"])
            for w, ff, a in ([("kal", o["_fit"]["kal"], a) for a in ("x=2", "x=3", "q_f")] +
                             [("kor", o["_fit"][dogru], ad_x(x)) for x in KOR["ana"]]):
                i = ff["adlar"].index(a)
                kk[f"{w}:{a}"] += int(abs(ff["th"][0, i] - th[w][a]) <= T7 * jk(ff["th"][1:, i]))
            hk = hukum_a(mB_, mR_, f["Z"])
            say[next(k for k in say if hk.startswith(k))] += 1
            tani.append(((mB_["R_Z"]["deger"] - 1) / mB_["R_Z"]["sigma_Z"],
                         (mR_["R_Z"]["deger"] - mR_["R_g_Z"]) / mR_["R_Z"]["sigma_Z"],
                         mB_["beta"]["deger"], mB_["R_Z"]["deger"], mR_["R_Z"]["deger"],
                         hukum_b(mB_, f["rho"])))
            sB.append(o["M4"]["sigma_eff_BK"])
            sR.append(o["M4"]["sigma_eff_Rg"])
            m4 += int(o["M4"]["gecti"])
        tn = np.array([x[:5] for x in tani], float)
        hb = [x[5] for x in tani]
        kap_tum = {a: v / nrep for a, v in kap.items()}           # KAYIT: tüm parametreler
        kap = {a: v / nrep for a, v in kk.items()}                  # (ii) karara girenler
        oran = {k: v / nrep for k, v in say.items()}
        dogru_k = "BK TUTAR" if dogru == "BK" else "R-g KAZANIR"
        yanlis_k = "R-g KAZANIR" if dogru == "BK" else "BK TUTAR"
        g_i = max(gh.values()) <= 0.05
        g_ii = min(kap.values()) >= 0.88
        g_iii = oran[dogru_k] >= 0.80 and oran[yanlis_k] <= 0.05
        gecti &= g_i and g_ii and g_iii
        sonuc[dogru] = {"i_gurultusuz_goreli_hata": gh, "i_gecti": bool(g_i),
                        "ii_kapsama_karara_giren": kap, "ii_gecti": bool(g_ii),
                        "KAYIT_tum_parametre_kapsama_2.36sigma_jk": kap_tum,
                        "iii_oranlar": oran, "iii_dogru_kazanma": oran[dogru_k],
                        "iii_yanlis_kazanma": oran[yanlis_k], "iii_gecti": bool(g_iii),
                        "sigma_eff_BK_medyan": float(np.median(sB)),
                        "sigma_eff_Rg_medyan": float(np.median(sR)), "M4_gecme": m4 / nrep,
                        "KAYIT_H197b_oranlari": {h: hb.count(h) / nrep for h in set(hb)},
                        "KAYIT_tani": {
                            "z_BK_jk_ort_sd": [float(np.mean(tn[:, 0])), float(np.std(tn[:, 0]))],
                            "z_Rg_jk_ort_sd": [float(np.mean(tn[:, 1])), float(np.std(tn[:, 1]))],
                            "beta_BK_ort_sd": [float(np.mean(tn[:, 2])), float(np.std(tn[:, 2]))],
                            "R_Z_BK_ort_sd": [float(np.mean(tn[:, 3])), float(np.std(tn[:, 3]))],
                            "R_Z_Rg_ort_sd": [float(np.mean(tn[:, 4])), float(np.std(tn[:, 4]))]},
                        "tohumlar": [tohumlar[dogru][0], tohumlar[dogru][-1]],
                        "gecti": bool(g_i and g_ii and g_iii)}
    return {"gecti": bool(gecti), "nrep": {d: len(tohumlar[d]) for d in MODELLER}, "S": S,
            "a0": a0, "b0": b0, "f": f, "gurultu_se_medyan": float(np.nanmedian(se_j)),
            "g_cal": G_CAL, "dogrular": sonuc}


# ============================ İKİNCİL (+ yan) ve EK (− yan) ============================
def ikincil(prof, H, J, pencere, mu0, kat, isaret):
    bil = bilesenler_yan(pencere, 0.05, mu0)
    f = uyum(prof, H, J, pencere, [(str(r), float(np.log(float(r)))) for r in bil])
    ref = F(2) if isaret > 0 else None
    Aref = par(f, str(ref)) if ref is not None else None
    konum = {}
    for ad, k in kat.items():
        if np.sign(k["delta_omega"]) != isaret or not (
                pencere[0] - EPS <= k["delta_omega"] <= pencere[1] + EPS):
            continue
        r = k["r"]
        Ar = par(f, str(r))
        if Ar is None:
            konum[ad] = {"durum": "tanımsız/bileşen yok", "c_BK": c_bk(r)}
            continue
        d = {"A": ozet(Ar), "A_bolu_se": float(Ar[0] / jk(Ar[1:])) if jk(Ar[1:]) > 0 else None,
             "c_BK": c_bk(r), "liste_192": k["liste"]}
        if Aref is not None and c_bk(r) > 0:
            d["R_plus"] = ozet((Ar / c_bk(r)) / (Aref / c_bk(ref)))
        konum[ad] = d
    return {"pencere": list(pencere), "n_bilesen": len(bil),
            "bilesenler": {str(r): c for r, c in bil.items()},
            "tanimsiz": f["dusen"], "rank": f["rank"], "kosul": f["kosul"],
            "konumlar_192": konum}, f


# ============================ figür ============================
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#a3a29c"
C_BKM, C_BK, C_RG = "#2a78d6", "#eb6834", "#1baf7a"


def ciz(yol, prof, J, o, H, fik, baslik):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK2,
                         "xtick.color": INK2, "ytick.color": INK2, "axes.titlesize": 10})
    se = np.array([jk(prof[1:, j]) for j in range(len(J))])
    c = J * DW
    fig, ax = plt.subplots(2, 2, figsize=(14, 8.4))
    fk = o["_fit"]
    s = o["_s"][0]

    def veri(a, m):
        a.fill_between(c[m], (prof[0] - se)[m], (prof[0] + se)[m], color=MUTED, alpha=0.35,
                       lw=0, label="κ_Σ ± jk se")
        a.plot(c[m], prof[0, m], color=INK, lw=1.2, marker="o", ms=3, label="κ_Σ (ölçülen)")

    def etiket(a, hedefler):
        dik = len(hedefler) > 8
        for ad, h in hedefler:
            a.axvline(h, color=MUTED, lw=0.7, zorder=0)
            a.text(h, 0.99 if dik else 1.0, ad, transform=a.get_xaxis_transform(),
                   ha="center", va="top" if dik else "bottom", rotation=90 if dik else 0,
                   fontsize=6 if dik else 7, color=INK2)
        a.axhline(0, color=MUTED, lw=0.6)
        a.set_xlabel("Δω (blok-yerel)")
        a.set_ylabel("κ_Σ")

    a = ax[0, 0]
    m = fk["BK"]["m"]
    veri(a, m)
    a.plot(c[m], fk["BK"]["model_toplam"][0], color=C_BKM, lw=2, label="M_BK uyumu")
    a.plot(c[m], fk["Rg"]["model_toplam"][0], color=C_RG, lw=1.6, ls=":", label="M_Rg uyumu")
    thB = fk["BK"]["th"][0]
    ong = thB[-2] + thB[-1] * c[m] + fk["BK"]["sabit"][0]
    for x in KOR["ana"]:
        ong = ong + s / float(x) * H(-np.log(float(x)))[0, m]
    a.plot(c[m], ong, color=C_BK, lw=1.6, ls="--", label="BK öngörüsü (A = s/x)")
    etiket(a, [(str(x), -np.log(float(x))) for x in KOR["ana"]])
    a.set_title("KÖR BANT [−2.12, −1.58] — iki model", pad=14)
    a.legend(fontsize=7, frameon=False, loc="lower left")
    a = ax[0, 1]
    m = fk["kal"]["m"]
    veri(a, m)
    a.plot(c[m], fk["kal"]["model"][0], color=C_BKM, lw=2, label="kalibrasyon uyumu")
    etiket(a, [(str(x), -np.log(float(x))) for x in KAL["ana"]])
    a.set_title("KALİBRASYON [−1.30, −0.55] (görülmüş)", pad=14)
    a.legend(fontsize=7, frameon=False, loc="lower left")
    a = ax[1, 0]
    m = fik["m"]
    veri(a, m)
    a.plot(c[m], fik["model"][0], color=C_BKM, lw=2, label="profil uyumu")
    etiket(a, [(k, v["delta_omega"]) for k, v in kat192().items()
               if IKINCIL_PENCERE[0] <= v["delta_omega"] <= IKINCIL_PENCERE[1]])
    a.set_title(f"İKİNCİL + yan {list(IKINCIL_PENCERE)} (KAYIT, hüküm dışı)", pad=14)
    a.legend(fontsize=7, frameon=False, loc="lower left")
    a = ax[1, 1]
    xs = np.array([float(x) for x in KOR["ana"]])
    mB, mR = o["modeller"]["BK"], o["modeller"]["Rg"]
    for mo, off, mk, fc in ((mB, -0.05, "o", INK), (mR, 0.05, "s", "white")):
        a.errorbar(xs + off, [mo["R"][ad_x(x)]["deger"] for x in KOR["ana"]],
                   yerr=[mo["R"][ad_x(x)]["se"] for x in KOR["ana"]], fmt=mk, color=INK,
                   mfc=fc, ms=5, capsize=3,
                   label=f"R_x — {'M_BK' if mo is mB else 'M_Rg'} (±jk se)")
    xx = np.linspace(4.8, 8.3, 100)
    a.axhline(1.0, color=C_BK, lw=1.6, ls="--", label="BK: R = 1")
    a.plot(xx, g_rg(L_G - np.log(xx)) / mR["g_bar"], color=C_RG, lw=1.6, ls=":",
           label="R-g: R^g(x)")
    a.set_xlabel("x (Δω = −log x)")
    a.set_ylabel("R_x = A_x/(s/x)")
    a.set_title(f"R̄_Z: M_BK {mB['R_Z']['deger']:.3f}±{mB['R_Z']['sigma_Z']:.3f} | M_Rg "
                f"{mR['R_Z']['deger']:.3f}±{mR['R_Z']['sigma_Z']:.3f} (R̄^g {mR['R_g_Z']:.3f})")
    a.legend(fontsize=7, frameon=False)
    fig.suptitle(baslik, fontsize=11, color=INK)
    fig.tight_layout()
    fig.savefig(yol, dpi=150)
    plt.close(fig)


# ============================ ana akış ============================
def temiz(d):
    if isinstance(d, dict):
        return {k: temiz(v) for k, v in d.items() if not str(k).startswith("_")}
    return d


def calistir(prof_duz, prof_ham, H, J, m1m3_ok, f, m6_se=None, cikti=None, png=None,
             baslik="", m2_atla=False):
    o = hesap(prof_duz, H, J, f=f)
    se_j = np.array([jk(prof_duz[1:, j]) for j in range(len(J))]) if m6_se is None else m6_se
    thB = o["_fit"]["BK"]["th"][0]
    M6 = m6(H, J, se_j, s192(), float(thB[-2]), float(thB[-1]), f)
    hk = hukumler(o, M6["gecti"], f, m1m3_ok, m2_atla)
    kat = kat192()
    ik, fik = ikincil(prof_duz, H, J, IKINCIL_PENCERE,
                      [k["r"] for k in kat.values() if k["liste"] == "sonmeli"], kat, +1)
    ek, _ = ikincil(prof_duz, H, J, EK_PENCERE, [], kat, -1)
    oh = hesap(prof_ham, H, J, f=f)
    hk_ham = hukumler(oh, M6["gecti"], f, m1m3_ok, m2_atla)
    sonuc = {"birincil": "K_duz (S1b)", "hukum": hk["hukum"], "kapi_etiketi": hk["kapi_etiketi"],
             "kural_ciktisi": hk["kural_ciktisi"], "f": f,
             "M2_atlandi_sentetik": bool(m2_atla),
             "KAYIT_Rg_modeli_hukumleri": hk["KAYIT_Rg_modeli"], "KAYIT_beta": hk["KAYIT_beta"],
             "kapilar": {**o["kapilar"], "M6_hukum_ani": M6},
             "nicelikler": temiz(o),
             "KAYIT": {"K_ham": {"hukum": hk_ham["hukum"], "nicelikler": temiz(oh)},
                       "IKINCIL_arti": ik,
                       "EK_eksi_KALEM_DISI": {**ek, "not": "kabul edilen EK (KAYIT): −log(3/2), "
                                                           "−log(4/3) için − yan pencere"}}}
    if png:
        ciz(png, prof_duz, J, o, H, fik, baslik)
    if cikti:
        json.dump(sonuc, open(cikti, "w"), indent=1, ensure_ascii=False, default=float)
    return sonuc


def kinematik32():
    E = np.load(ZD / "eta_dusuk_t0.4_c4000.npz")
    mid = np.asarray(E["mid"], float)
    kb, Lb = h197.bloklar(len(mid), mid, h197.NBLOK)
    return mid, kb, Lb


def m1_gurultu():
    """Ölçüm-öncesi gürültü ölçeği: M1 ω-parçaları (görülmüş bant, j ≥ −51; 8 blok, 190
    düzeni), HAVUZ' sütunları (τ < 0.74), mix tam + 8 loo → 8-blok κ_Σ profili jk se(j)."""
    b188 = h197.b188
    G = np.load(ZD / "K1_gercek_dusuk.npz")
    OZ = np.load(ZD / "OZ_gercek_dusuk.npz")
    P = np.load(ZD / "G1_proj_gercek_dusuk.npz")
    tau, aq = np.asarray(G["tau"], float), np.asarray(G["aq"], float)
    win = (tau >= 0.45) & (tau < 0.86)
    hv = tau[win] < h197.HAVUZ_197[1]
    mixes = [b188.karisim(P, OZ, aq, v - 1)[win][hv] for v in range(NJACK + 1)]
    ONK = json.load(open(S190 / "ONKAYIT_190.json"))
    L, Lb = float(ONK["pencere"]["L"]), np.array(ONK["bloklar"]["L_b"])
    parca = sorted((S197 / "m1" / "omega").glob("b*_c*.npz"))
    assert parca, "M1 ω-parçaları yok (önce 197b m1)"
    jmin = min(int(np.load(p)["j"].min()) for p in parca)
    jmax = max(int(np.load(p)["j"].max()) for p in parca)
    assert jmin >= h197.J_GOR
    J = np.arange(jmin, jmax + 1)
    T = np.full((NJACK + 1, NJACK, len(J)), np.nan)
    nb = np.zeros(NJACK)
    for p in parca:
        d = np.load(p)
        b = int(d["blok"])
        C = 2.0 * (d["re"] + 1j * d["im"])[:, hv] / int(d["nb"])
        for v, mx in enumerate(mixes):
            T[v, b, d["j"] - jmin] = -(C @ np.conj(mx)).real / np.sum(np.abs(mx) ** 2)
        nb[b] = int(d["nb"])
    c = J * DW
    kap = np.array([(c - DW / 2 >= 0.86 * L - Lb[b] - EPS) & (c + DW / 2 <= 1.30 * L - Lb[b] + EPS)
                    for b in range(NJACK)]) & np.all(np.isfinite(T), 0)
    rep = np.full((NJACK, len(J)), np.nan)
    with np.errstate(invalid="ignore", divide="ignore"):
        for i in range(NJACK):
            kk = kap.copy()
            kk[i] = False
            rep[i] = np.nansum(np.where(kk, T[1 + i] * nb[:, None], 0), 0) / (nb[:, None] * kk).sum(0)
    se = np.array([jk(rep[:, j]) if np.all(np.isfinite(rep[:, j])) and kap[:, j].sum() >= 4
                   else np.nan for j in range(len(J))])
    return J, se


def se_esle(Jh, J_m1, se_m1):
    """M1 se(j)'yi 32-blok ızgarasına taşı: ortak j → aynı değer; yoksa medyan (kör bant dahil)."""
    med = float(np.nanmedian(se_m1))
    out = np.full(len(Jh), med)
    for k, j in enumerate(Jh):
        i = j - J_m1[0]
        if 0 <= i < len(J_m1) and np.isfinite(se_m1[i]):
            out[k] = se_m1[i]
    return out, med


def hazirlik():
    mid, kb, Lb = kinematik32()
    jlo, jhi = int(round(h197.DW_ARALIK[0] / DW)), int(round(h197.DW_ARALIK[1] / DW))
    J = np.arange(jlo, jhi + 1)
    H = Bicim(mid, kb, Lb, np.ones((len(Lb), len(J)), bool), J)
    J1, se1 = m1_gurultu()
    se_j, med = se_esle(J, J1, se1)
    return J, H, se_j, med, int(np.isfinite(se1).sum())


F_DOSYA = S197 / "F_KALIBRASYON.json"
# DÖRDÜNCÜ MAKİNE RAPORU: f SABİT (f_kal, K_KAL tohumları, 2×200; yeniden kalibre EDİLMEZ)
F_SABIT = {"Z": 1.5650650491670195, "rho": 1.3736133274633446}


def f_oku():
    d = json.load(open(F_DOSYA))
    if d["f_Z"] != F_SABIT["Z"] or d["f_rho"] != F_SABIT["rho"]:
        raise SystemExit("F_KALIBRASYON değerleri SABİT f'lerle uyumsuz")
    return dict(F_SABIT)


def m6_yazdir(M6):
    for d, r in M6["dogrular"].items():
        print(f"  [{d}-doğru] (i) gürültüsüz maks göreli hata {max(r['i_gurultusuz_goreli_hata'].values()):.2e}"
              f" → {r['i_gecti']}; (ii) 2.36σ kapsama min {min(r['ii_kapsama_karara_giren'].values()):.3f}"
              f" → {r['ii_gecti']}; (iii) doğru kazanma {r['iii_dogru_kazanma']:.3f}, yanlış "
              f"{r['iii_yanlis_kazanma']:.3f} → {r['iii_gecti']}", flush=True)
        print(f"      oranlar {r['iii_oranlar']}; σ_eff^BK medyan {r['sigma_eff_BK_medyan']:.3f}, "
              f"σ_eff^Rg medyan {r['sigma_eff_Rg_medyan']:.3f}; M4 geçme {r['M4_gecme']:.2f}; "
              f"H-197b {r['KAYIT_H197b_oranlari']}", flush=True)
        print(f"      kapsama: { {a: round(v, 3) for a, v in r['ii_kapsama_karara_giren'].items()} }",
              flush=True)
    print(f"  M6 → {'GEÇTİ' if M6['gecti'] else 'KALDI'}", flush=True)


def mod_f_kal():
    if F_DOSYA.exists():
        raise SystemExit("f SABİT (DÖRDÜNCÜ MAKİNE RAPORU) — yeniden kalibre EDİLMEZ")
    J, H, se_j, med, n1 = hazirlik()
    S = s192()
    print(f"  gürültü (M1, görülmüş): medyan jk se = {med:.3e}; S = {S:.4f}; tohumlar K_KAL "
          f"{ {d: [K_KAL[d][0], K_KAL[d][-1]] for d in MODELLER} }", flush=True)
    fk = f_kalibrasyon(H, J, se_j, S, 2 * med, med)
    for d, r in fk["dogrular"].items():
        print(f"  [{d}-doğru, {r['model']}] z_Z ort {r['z_Z_ort']:+.3f} sd {r['z_Z_sd']:.3f}; "
              f"z_ρ ort {r['z_rho_ort']:+.3f} sd {r['z_rho_sd']:.3f} (n={r['n']})", flush=True)
    print(f"  f_Z = {fk['f_Z']:.4f}   f_ρ = {fk['f_rho']:.4f}", flush=True)
    fk.update({"kod_sha256": sha(Path(__file__).resolve()),
               "zaman": time.strftime("%Y-%m-%d %H:%M:%S %z")})
    json.dump(fk, open(F_DOSYA, "w"), indent=1, ensure_ascii=False, default=float)
    print(f"  -> {F_DOSYA}", flush=True)


def mod_m6_on():
    J, H, se_j, med, n1 = hazirlik()
    S = s192()
    f = f_oku()
    print(f"  gürültü (M1, görülmüş, 8 blok, HAVUZ'): medyan jk se = {med:.3e} ({n1} dilim; kör "
          f"pencerede medyan); S = {S:.4f}; ḡ_cal = {G_CAL:.3f}; f = {f}; güç tohumları K_GUC",
          flush=True)
    M6 = m6(H, J, se_j, S, 2 * med, med, f)
    m6_yazdir(M6)
    out = {"tur": "ÖLÇÜM-ÖNCESİ M6 (iki-doğrulu; σ_eff = f·σ_jk; gürültü M1, sinyal 192)",
           "M6": M6, "F_KALIBRASYON_sha256": sha(F_DOSYA),
           "kod_sha256": sha(Path(__file__).resolve()),
           "zaman": time.strftime("%Y-%m-%d %H:%M:%S %z")}
    json.dump(out, open(S197 / "M6_on_deneme.json", "w"), indent=1, ensure_ascii=False,
              default=float)
    print(f"  -> {S197/'M6_on_deneme.json'}", flush=True)


def mod_sentetik():
    """Uçtan uca prova: BK-doğru sentetik κ_Σ (tüm c ≥ 0.02 uyduları, iki yan) + jk gürültü.
    M2 ATLANIR (gerçek-veri makine kapısı; BK-doğru sentetikte kalması beklenir)."""
    J, H, se_j, med, _ = hazirlik()
    rng = np.random.default_rng(1974000)
    S = s192()
    f = f_oku()
    bil = {**bilesenler_yan((-2.3, -0.02), 0.02), **bilesenler_yan((0.02, 2.45), 0.02)}
    t = np.zeros((NJACK + 1, len(J))) + 2 * med + med * J * DW
    for r, cc in bil.items():
        t += S * cc * H(np.log(float(r)))
    d = S197 / "sentetik_test"
    d.mkdir(exist_ok=True)
    son = calistir(jk_gurultu(t, se_j, rng), jk_gurultu(t, se_j, rng), H, J, True, f,
                   m6_se=se_j, cikti=d / "HUKUM_197_sentetik.json",
                   png=d / "197_ayna_yasasi_sentetik.png", m2_atla=True,
                   baslik="197 — SENTETİK PROVA (BK-doğru; gerçek veri DEĞİL)")
    n = son["nicelikler"]["modeller"]
    k = son["kapilar"]
    print(f"  sentetik: {len(bil)} BK uydusu (c ≥ 0.02), S = {S:.4f}, f = {f}", flush=True)
    print(f"  M2 ATLANDI (sentetik prova; gerçek-veri kapısı) — bilgi: M2 {k['M2']['gecti']}, "
          f"oranlar { {a: round(v['oran'], 3) for a, v in k['M2']['oranlar'].items()} }",
          flush=True)
    print(f"  hüküm: {son['hukum']}", flush=True)
    for mo in MODELLER:
        print(f"  {mo}: R̄_Z = {n[mo]['R_Z']['deger']:.3f} ± {n[mo]['R_Z']['sigma_Z']:.3f} (jk); "
              f"R̄^g = {n[mo]['R_g_Z']:.3f}; β = {n[mo]['beta']['deger']:.2f}; "
              f"ρ = {n[mo]['rho']['deger']:.2f}; ψ = {n[mo]['psi']['deger']:.2f}", flush=True)
    print(f"  kapılar: M4 {k['M4']['gecti']} M5 {k['M5']['gecti']} (oran {k['M5']['oran']:.2f}, "
          f"korr {k['M5']['mlog2_korelasyon']:.2f}) M6 {k['M6_hukum_ani']['gecti']}", flush=True)
    print(f"  -> {d}", flush=True)


def mod_hukum():
    o = json.load(open(S197 / "ONKAYIT_197.json"))
    if sha(QM / "197_configs" / "197a_onkayit.py") != o["sha256"]:
        raise SystemExit("ON-KAYIT SHA UYUMSUZ")
    if sha(KALEM) != o["kalem_sha256"]:
        raise SystemExit("KALEM ön-kayıttan sonra DEĞİŞMİŞ")
    for rel, h in o["girdi_sha256"].items():
        if sha(QM / rel) != h:
            raise SystemExit(f"GİRDİ SHA UYUMSUZ: {rel}")
    on = json.load(open(S197 / "M6_on_deneme.json"))
    if not on["M6"]["gecti"]:
        raise SystemExit("ölçüm-öncesi M6 KALMIŞ — KALEM: sınav yapılmaz (tasarım yeniden)")
    f = {"Z": float(o["f"]["f_Z"]), "rho": float(o["f"]["f_rho"])}      # ön-kayıtta SABİT
    assert f == f_oku(), "ön-kayıt f'leri F_KALIBRASYON ile uyumsuz"
    K1 = json.load(open(S197 / "K1_harita_197.json"))
    assert K1["sha_onkayit"] == o["sha256"], "harita başka ön-kayıtla üretilmiş"
    m1 = json.load(open(S197 / "m1" / "M1_sonuc.json"))
    m1m3 = bool(m1["gecti"] and K1["M3"]["kesisim"] == 0 and max(K1["M_kesik"]) <= 1e-10)
    Hm = np.load(S197 / "harita_omega32.npz")
    assert str(Hm["birincil"]) == "duz"
    tam = Hm["tam"]
    J = Hm["J"][tam]
    pd_ = Hm["kappa_sigma_duz"][:, tam]
    ph_ = Hm["kappa_sigma_ham"][:, tam]
    mid, kb, Lb = kinematik32()
    assert np.array_equal(kb, Hm["kb"]) and np.array_equal(Lb, Hm["L_b"])
    H = Bicim(mid, kb, Lb, Hm["kapsama"][:, tam], J)
    son = calistir(pd_, ph_, H, J, m1m3, f, png=S197 / "197_ayna_yasasi.png",
                   baslik="197 — AYNA YASASI: kör bant profil sınavı (BİRİNCİL K_düz)")
    son.update({"sha_onkayit": o["sha256"], "kod_sha256": sha(Path(__file__).resolve()),
                "M1": {"gecti": m1["gecti"], "maks_bagil": m1["maks_bagil_eleman"]},
                "M3": K1["M3"], "M_kesik": K1["M_kesik"], "M6_on": on["M6"]["gecti"],
                "zaman": time.strftime("%Y-%m-%d %H:%M:%S %z")})
    json.dump(son, open(S197 / "HUKUM_197.json", "w"), indent=1, ensure_ascii=False, default=float)
    print(json.dumps({"hukum": son["hukum"], "kapi_etiketi": son["kapi_etiketi"]},
                     ensure_ascii=False, indent=1))
    print(f"-> {S197/'HUKUM_197.json'}, {S197/'197_ayna_yasasi.png'}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("mod", choices=["hukum", "f_kal", "m6_on", "sentetik"])
    a = ap.parse_args()
    print("=" * 78)
    print(f"197c / {a.mod.upper()}  {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 78, flush=True)
    {"hukum": mod_hukum, "f_kal": mod_f_kal, "m6_on": mod_m6_on,
     "sentetik": mod_sentetik}[a.mod]()
