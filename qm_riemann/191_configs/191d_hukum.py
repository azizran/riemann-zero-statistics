# -*- coding: utf-8 -*-
"""
191d — K3: HÜKÜM (H-191a/b/c) + KİMLİK KAYDI (4 derinlik)
===========================================================
ONKAYIT_191'de donan tanımlar ve eşiklerle (kurtarma yok):
  H-191a  zarf sürer (yön): f_1.30 > f_1.20 + 2se (HAVUZ, ortak-blok loo)
          VE 8 bantta r_1.20 ≤ r_1.30 + 2σ (184 konv., bağımsız yayılım).
          ÖLÜM: f_1.30 ≤ f_1.20 + 2se.
  H-191b  harita-doğrusal yasa (nicel): f_1.30 ∈ [0.58, 0.64].
  H-191c  σ_ε yakınsaması: σ_ε(1.30) < σ_ε(1.20) − 2se VE ≥ 0.40921 − 2se
          (aşma yok); nicel bant [0.4140, 0.4160].
  K3      dört derinlikte (1.00,1.10,1.20,1.30) r_D(τ)=1−c·τ^α (184c AYNEN);
          f_D, g_D (σ-kapanış payı) vs Δz_D doğrusallığı (orijinden-geçen +
          serbest fit); yasa tutarsa f_∞ = eğim×0.2529 ± se; 189 A/B ile
          yan yana.
Girdi: 189/zincir_{Hkeskin,Hderin110,Hderin120,Hderin130}.{json,npz},
       189/zincir_gercek.npz, 189/muhur_Hkeskin.json, 189/HUKUM_189.json,
       191/ONKAYIT_191.json, 191/ZINCIR_191_ozet.json,
       184/K1_gercek.npz, 189/K1_<ad>.npz.
Çıktı: 191/HUKUM_191.json (+ ekran defteri).
"""
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S184, S189, S191 = SCR / "184", SCR / "189", SCR / "191"
NJACK = 8
DER = {"1.00": "Hkeskin", "1.10": "Hderin110", "1.20": "Hderin120",
      "1.30": "Hderin130"}
BANT = ["0.45-0.50", "0.50-0.55", "0.55-0.60", "0.60-0.65",
        "0.65-0.70", "0.70-0.75", "0.75-0.80", "0.80-0.86"]
IH = 9          # 185-defterinde HAVUZ satırı (0..7 bant, 8 KUYRUK, 9 HAVUZ)
IH7 = 8         # 187c defterinde HAVUZ satırı

spec = importlib.util.spec_from_file_location(
    "b184c", QM / "184_configs" / "184c_K2_yarisma.py")
b184c = importlib.util.module_from_spec(spec)
spec.loader.exec_module(b184c)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def jk(reps):
    reps = np.asarray(reps, float)
    return np.sqrt((NJACK - 1) / NJACK * np.sum((reps - reps.mean(0)) ** 2, 0))


def onkayit():
    o = json.load(open(S191 / "ONKAYIT_191.json"))
    s = sha(QM / "191_configs" / "191a_onkayit.py")
    if s != o["sha256"]:
        raise SystemExit(f"ON-KAYIT SHA UYUMSUZ: {s} != {o['sha256']}")
    if sha(QM / "189_configs" / "189c_zincir.py") != o["zincir"]["sha_189c"]:
        raise SystemExit("189c K0'dan beri DEĞİŞMİŞ")
    if sha(QM / "189_configs" / "189d_hukum.py") != o["zincir"]["sha_189d"]:
        raise SystemExit("189d K0'dan beri DEĞİŞMİŞ")
    return o


def yukle_zincir(ad):
    js = json.load(open(S189 / f"zincir_{ad}.json"))
    nz = np.load(S189 / f"zincir_{ad}.npz", allow_pickle=True)
    return js, nz


def fit184c(g, tw):
    """184c makinesi AYNEN (189d.fit184c ile bit-bit aynı çağrı sırası):
    bant_defteri + H_Z2 fit; ek: loo-jk yeniden fit."""
    sat = b184c.bant_defteri(g, tw, "x")[:8]
    tau = np.array([s["tau"] for s in sat])
    r = np.array([s["r"] for s in sat])
    se = np.array([s["se_r"] for s in sat])
    p0, bnd = [0.3, 3.0], ([0.0, 0.1], [10.0, 20.0])
    popt, pcov, chi2, err = b184c.fit(b184c.H_Z2, tau, r, se, p0, bnd)
    perr = np.sqrt(np.diag(pcov)) if pcov is not None else [np.nan, np.nan]
    ah_g, jk_g = b184c.band_ahat(g)
    ah_t, jk_t = b184c.band_ahat(tw)
    KEN = b184c.KENAR
    reps = []
    for k in range(NJACK):
        rk = []
        for i in range(8):
            msk = (g["tau"] > KEN[i]) & (g["tau"] <= KEN[i + 1])
            rk.append(jk_g[k][msk].sum() / jk_t[k][msk].sum())
        pk, _, _, _ = b184c.fit(b184c.H_Z2, tau, np.array(rk), se, p0, bnd)
        reps.append(pk if pk is not None else [np.nan, np.nan])
    reps = np.array(reps)
    sinir = bool(popt is not None and (popt[0] <= 1e-9 or popt[0] >= 10 - 1e-9
                                       or popt[1] <= 0.1 + 1e-9
                                       or popt[1] >= 20 - 1e-9))
    return dict(tau=tau.tolist(), r=r.tolist(), se_r=se.tolist(),
                c=float(popt[0]), alfa=float(popt[1]),
                perr=[float(x) for x in perr], jk_se=jk(reps).tolist(),
                chi2_dof=float(chi2 / 6), sinirda=sinir, hata=err)


def fit_serbest(tau, r, se):
    """ÖN-KAYITSIZ EK (hüküm dışı): işaret-serbest 1 − c·τ^α (c<0 izinli)."""
    from scipy.optimize import curve_fit
    f = lambda t, c, a: 1.0 - c * t ** a
    try:
        popt, pcov = curve_fit(f, tau, r, p0=[0.1, 1.3], sigma=se,
                               absolute_sigma=True,
                               bounds=([-10.0, 0.1], [10.0, 20.0]), maxfev=20000)
        chi2 = float(np.sum(((r - f(tau, *popt)) / se) ** 2))
        return dict(c=float(popt[0]), alfa=float(popt[1]),
                    perr=np.sqrt(np.diag(pcov)).tolist(), chi2_dof=chi2 / 6)
    except Exception as e:
        return dict(hata=str(e))


def dogrusal_orijin(x, y, se):
    """Orijinden geçen tek-parametre ağırlıklı en-küçük-kare: y = m·x."""
    x, y, se = np.asarray(x, float), np.asarray(y, float), np.asarray(se, float)
    w = 1.0 / se ** 2
    m = np.sum(w * x * y) / np.sum(w * x * x)
    se_m = np.sqrt(1.0 / np.sum(w * x * x))
    return float(m), float(se_m)


def dogrusal_serbest(x, y, se):
    """Serbest (kesişimli) ağırlıklı en-küçük-kare: y = a + m·x."""
    x, y, se = np.asarray(x, float), np.asarray(y, float), np.asarray(se, float)
    W = np.diag(1.0 / se ** 2)
    A = np.vstack([np.ones_like(x), x]).T
    cov = np.linalg.inv(A.T @ W @ A)
    beta = cov @ (A.T @ W @ y)
    return float(beta[1]), float(np.sqrt(cov[1, 1])), float(beta[0]), \
        float(np.sqrt(cov[0, 0]))


if __name__ == "__main__":
    ONK = onkayit()
    print("=" * 78)
    print(f"191d / HÜKÜM (H-191a/b/c) + KİMLİK KAYDI  "
          f"[on-kayit {ONK['zaman']} sha {ONK['sha256'][:12]}]")
    print("=" * 78, flush=True)
    out = {"sha_onkayit": ONK["sha256"], "damga_onkayit": ONK["zaman"]}

    muhur = json.load(open(S189 / "muhur_Hkeskin.json"))
    out["makine_muhru"] = muhur
    print(f"makine mührü (Hkeskin bit-bit, 189'dan referans): "
          f"{'TUTTU' if muhur['HEPSI'] else 'TUTMADI'}")
    if not muhur["HEPSI"]:
        raise SystemExit("MAKİNE MÜHRÜ TUTMADI — hüküm verilmez")

    insa130 = json.load(open(S189 / "insa_kapi_Hderin130.json"))
    print(f"inşa kapıları (Hderin130): {'GEÇTİ' if insa130['HEPSI'] else 'KALDI'}"
          f"  maks|F|={insa130['maks_F']:.3e}  sıralı={insa130['sirali']}  "
          f"L={insa130['L_8hane']}")
    out["insa_kapisi_Hderin130"] = insa130

    Z = {D: yukle_zincir(DER[D]) for D in DER}
    js = {D: Z[D][0] for D in Z}
    nz = {D: Z[D][1] for D in Z}
    ZG = np.load(S189 / "zincir_gercek.npz")

    # ---------------- bant tablosu ----------------
    print("\nBANT TABLOSU (185-defter; r = w_g/w_HD ± σ_r [184 konv.]):")
    for D in ("1.00", "1.10", "1.20", "1.30"):
        j = js[D]
        print(f"  D={D} ({DER[D]}): L={j['L']:.9f}  N={j['N']}")
        for k in list(range(8)) + [IH]:
            et = j["etiket185"][k]
            print(f"   {et:>11} w_HD={j['w_HD'][k]:.4f}±{j['se_w_HD'][k]:.4f} "
                  f"w^öz={j['w_oz_HD'][k]:.4f}±{j['se_w_oz_HD'][k]:.4f} "
                  f"m={j['m_HD'][k]:.4f}±{j['se_m_HD'][k]:.4f} "
                  f"r={j['r'][k]:.4f}±{j['sig_r'][k]:.4f}")

    r = {D: np.array(js[D]["r"]) for D in Z}
    s = {D: np.array(js[D]["sig_r"]) for D in Z}
    rr = {D: nz[D]["d185_reps"][:, :, 0] for D in Z}       # (8, 10)

    f = {}
    f_se = {}
    for D in ("1.10", "1.20", "1.30"):
        f[D] = (r[D] - r["1.00"]) / (1 - r["1.00"])
        reps_f = (rr[D] - rr["1.00"]) / (1 - rr["1.00"])
        f_se[D] = jk(reps_f)

    # ---------------- H-191a ----------------
    d_f_reps = ((rr["1.30"] - rr["1.00"]) / (1 - rr["1.00"])
               - (rr["1.20"] - rr["1.00"]) / (1 - rr["1.00"]))
    d_f_se = jk(d_f_reps)
    sD = np.sqrt(s["1.20"] ** 2 + s["1.30"] ** 2)
    c_havuz = bool(f["1.30"][IH] > f["1.20"][IH] + 2 * d_f_se[IH])
    band_c = [bool(r["1.20"][k] <= r["1.30"][k] + 2 * sD[k]) for k in range(8)]
    band_gecen = int(sum(band_c))
    kosul = c_havuz and band_gecen == 8
    olum = bool(f["1.30"][IH] <= f["1.20"][IH] + 2 * d_f_se[IH])
    hk_a = "ÖLDÜ" if olum else ("MÜHÜR" if kosul else "KAYIT")
    ha = dict(
        f110_havuz=float(f["1.10"][IH]), f110_se=float(f_se["1.10"][IH]),
        f120_havuz=float(f["1.20"][IH]), f120_se=float(f_se["1.20"][IH]),
        f130_havuz=float(f["1.30"][IH]), f130_se=float(f_se["1.30"][IH]),
        d_f_130_120=float(f["1.30"][IH] - f["1.20"][IH]),
        d_f_130_120_se=float(d_f_se[IH]),
        havuz_yon_kosulu=c_havuz, bant_kosulu=[bool(x) for x in band_c],
        bant_gecen=band_gecen, kosul=kosul, olum=olum, hukum=hk_a,
        r_havuz={D: float(r[D][IH]) for D in r},
        sig_r_havuz={D: float(s[D][IH]) for D in s})
    out["H-191a"] = ha
    print(f"\nH-191a: HAVUZ r = {r['1.00'][IH]:.4f} → {r['1.10'][IH]:.4f} → "
          f"{r['1.20'][IH]:.4f} → {r['1.30'][IH]:.4f}±{s['1.30'][IH]:.4f}")
    print(f"  f_1.10={f['1.10'][IH]:.4f}±{f_se['1.10'][IH]:.4f}  "
          f"f_1.20={f['1.20'][IH]:.4f}±{f_se['1.20'][IH]:.4f}  "
          f"f_1.30={f['1.30'][IH]:.4f}±{f_se['1.30'][IH]:.4f}")
    print(f"  Δf(1.30−1.20)={ha['d_f_130_120']:+.4f}±{ha['d_f_130_120_se']:.4f}"
          f"  (>2se yön: {c_havuz});  8-bant r_1.20≤r_1.30+2σ: "
          f"{band_gecen}/8  -> {hk_a}")
    for k in range(8):
        print(f"   {BANT[k]}: r 1.00={r['1.00'][k]:.4f} 1.10={r['1.10'][k]:.4f}"
              f" 1.20={r['1.20'][k]:.4f} 1.30={r['1.30'][k]:.4f}  "
              f"f130={f['1.30'][k]:.3f}±{f_se['1.30'][k]:.3f}  "
              f"c={'E' if band_c[k] else 'H'}")

    # ---------------- H-191b ----------------
    bant_lo, bant_hi = 0.58, 0.64
    f130 = float(f["1.30"][IH])
    hb = dict(f130=f130, f130_se=float(f_se["1.30"][IH]), bant=[bant_lo, bant_hi],
             tuttu=bool(bant_lo <= f130 <= bant_hi),
             hukum="MÜHÜR" if bant_lo <= f130 <= bant_hi else "ÖLDÜ")
    out["H-191b"] = hb
    print(f"\nH-191b: f_1.30={f130:.4f}±{hb['f130_se']:.4f}  bant "
          f"[{bant_lo},{bant_hi}]  -> {hb['hukum']}")

    # ---------------- H-191c ----------------
    sg = {D: nz[D]["sigma_eps"].item() for D in Z}
    sg_se = {D: nz[D]["sigma_eps_se"].item() for D in Z}
    sg_reps = {D: nz[D]["sigma_eps_reps"] for D in Z}
    sgg = float(ZG["sigma_eps"]); sgg_se = float(ZG["sigma_eps_se"])
    sgg_reps = ZG["sigma_eps_reps"]
    d_120_130 = sg["1.20"] - sg["1.30"]
    d_120_130_se = float(jk(sg_reps["1.20"] - sg_reps["1.30"]))
    d_gercek = sg["1.30"] - sgg
    d_gercek_se = float(jk(sg_reps["1.30"] - sgg_reps))
    asma_yok = bool(sg["1.30"] >= sgg - 2 * d_gercek_se)
    azaliyor = bool(sg["1.30"] < sg["1.20"] - 2 * d_120_130_se)
    bant_lo_s, bant_hi_s = 0.4140, 0.4160
    bant_ic = bool(bant_lo_s <= sg["1.30"] <= bant_hi_s)
    kosul_c = azaliyor and asma_yok and bant_ic
    olum_c = bool(sg["1.30"] >= sg["1.20"] - 2 * d_120_130_se
                 or sg["1.30"] < sgg - 2 * d_gercek_se or not bant_ic)
    hk_c = "ÖLDÜ" if olum_c else ("MÜHÜR" if kosul_c else "KAYIT")
    hc = dict(sigma_eps={D: sg[D] for D in sg}, se=sg_se, gercek=[sgg, sgg_se],
             fark_120_130=d_120_130, se_fark_120_130=d_120_130_se,
             fark_130_gercek=d_gercek, se_fark_130_gercek=d_gercek_se,
             azaliyor=azaliyor, asma_yok=asma_yok, bant=[bant_lo_s, bant_hi_s],
             bant_ic=bant_ic,
             kapanan_130=(sg["1.00"] - sg["1.30"]) / (sg["1.00"] - sgg),
             hukum=hk_c)
    out["H-191c"] = hc
    print(f"\nH-191c: σ_ε Hk={sg['1.00']:.5f} → 110={sg['1.10']:.5f} → "
          f"120={sg['1.20']:.5f} → 130={sg['1.30']:.5f}±{sg_se['1.30']:.5f} "
          f"(gerçek {sgg:.5f}±{sgg_se:.5f})")
    print(f"  Δ(120−130)={d_120_130:+.5f}±{d_120_130_se:.5f} (azalıyor: "
          f"{azaliyor})  130−gerçek={d_gercek:+.5f}±{d_gercek_se:.5f} "
          f"(aşma yok: {asma_yok})  bant[{bant_lo_s},{bant_hi_s}] içinde: "
          f"{bant_ic}  kapanan pay={hc['kapanan_130']:.3f}  -> {hk_c}")

    # ---------------- K3 — KİMLİK KAYDI ----------------
    k3 = {"fit": {}, "fit_serbest_EK": {}}
    g = b184c.yukle("gercek")
    for D in DER:
        d = np.load(S189 / f"K1_{DER[D]}.npz")
        tw = {k: d[k] for k in d.files}
        fr = fit184c(g, tw)
        k3["fit"][D] = fr
        k3["fit_serbest_EK"][D] = fit_serbest(np.array(fr["tau"]),
                                             np.array(fr["r"]),
                                             np.array(fr["se_r"]))
    print("\nK3 — r_D(τ) = 1 − c·τ^α (dört derinlik, 184c makinesi AYNEN):")
    for D, fr in k3["fit"].items():
        ek = k3["fit_serbest_EK"][D]
        print(f"   D={D}: c={fr['c']:.4f}±{fr['perr'][0]:.4f}/±{fr['jk_se'][0]:.4f} "
              f"α={fr['alfa']:.3f}±{fr['perr'][1]:.3f}/±{fr['jk_se'][1]:.3f} "
              f"χ²/dof={fr['chi2_dof']:.3f}{'  [SINIRDA]' if fr['sinirda'] else ''}"
              f"   | EK serbest: c={ek.get('c', float('nan')):+.4f} "
              f"α={ek.get('alfa', float('nan')):.3f}")

    # --------- f, g vs Δz_D doğrusallığı (4 nokta; D=1.00 orijin) ---------
    # Δz_D KAYNAKTAN (K0'da zaten doğrulandı) — burada doğrudan ONK'tan al
    dzK = ONK["yeniden_hesap_kaynaktan"]["Delta_z_D"]
    dz = {"1.00": 0.0, "1.10": dzK["1.10"], "1.20": dzK["1.20"],
         "1.30": dzK["1.30"]}
    dz_inf = dzK["inf"]

    g_D = {}
    g_D_se = {}
    for D in ("1.10", "1.20", "1.30"):
        g_D[D] = (sg["1.00"] - sg[D]) / (sg["1.00"] - sgg)
        reps_g = ((sg_reps["1.00"] - sg_reps[D]) / (sg["1.00"] - sgg))
        g_D_se[D] = float(jk(reps_g))

    xs = np.array([dz["1.00"], dz["1.10"], dz["1.20"], dz["1.30"]])
    ys_f = np.array([0.0, f["1.10"][IH], f["1.20"][IH], f["1.30"][IH]])
    se_f = np.array([1e-6, f_se["1.10"][IH], f_se["1.20"][IH], f_se["1.30"][IH]])
    ys_g = np.array([0.0, g_D["1.10"], g_D["1.20"], g_D["1.30"]])
    se_g = np.array([1e-6, g_D_se["1.10"], g_D_se["1.20"], g_D_se["1.30"]])

    # orijin NOKTASI (D=1.00, Δz=0, f=g=0) bir ÖLÇÜM değil TANIM gereği
    # (f_1.00≡0, g_1.00≡0) — "serbest" (kesişimli) fit yalnız ÜÇ ÖLÇÜLEN
    # noktayla (1.10,1.20,1.30) yapılır: kesişimin 0'a ne kadar yakın
    # çıktığı, orijin-noktasını fite ZORLAMADAN, bağımsız bir sınamadır.
    m_f_o, se_m_f_o = dogrusal_orijin(xs[1:], ys_f[1:], se_f[1:])
    m_f_s, se_m_f_s, a_f_s, se_a_f_s = dogrusal_serbest(xs[1:], ys_f[1:], se_f[1:])
    m_g_o, se_m_g_o = dogrusal_orijin(xs[1:], ys_g[1:], se_g[1:])
    m_g_s, se_m_g_s, a_g_s, se_a_g_s = dogrusal_serbest(xs[1:], ys_g[1:], se_g[1:])

    f_inf = m_f_o * dz_inf
    f_inf_se = se_m_f_o * dz_inf
    g_inf = m_g_o * dz_inf
    g_inf_se = se_m_g_o * dz_inf

    k3["dogrusallik"] = dict(
        Delta_z={"1.00": dz["1.00"], "1.10": dz["1.10"], "1.20": dz["1.20"],
                "1.30": dz["1.30"], "inf": dz_inf},
        f={"1.10": float(f["1.10"][IH]), "1.20": float(f["1.20"][IH]),
          "1.30": float(f["1.30"][IH])},
        f_se={"1.10": float(f_se["1.10"][IH]), "1.20": float(f_se["1.20"][IH]),
             "1.30": float(f_se["1.30"][IH])},
        g={"1.10": float(g_D["1.10"]), "1.20": float(g_D["1.20"]),
          "1.30": float(g_D["1.30"])},
        g_se={"1.10": g_D_se["1.10"], "1.20": g_D_se["1.20"],
             "1.30": g_D_se["1.30"]},
        f_egim_orijin=m_f_o, f_egim_orijin_se=se_m_f_o,
        f_egim_serbest=m_f_s, f_egim_serbest_se=se_m_f_s,
        f_kesisim_serbest=a_f_s, f_kesisim_serbest_se=se_a_f_s,
        g_egim_orijin=m_g_o, g_egim_orijin_se=se_m_g_o,
        g_egim_serbest=m_g_s, g_egim_serbest_se=se_m_g_s,
        g_kesisim_serbest=a_g_s, g_kesisim_serbest_se=se_a_g_s,
        f_inf=f_inf, f_inf_se=f_inf_se, kalinti_inf=1 - f_inf,
        kalinti_inf_se=f_inf_se, g_inf=g_inf, g_inf_se=g_inf_se)
    print(f"\nDoğrusallık f_D vs Δz_D (4 nokta): orijinden eğim="
          f"{m_f_o:.4f}±{se_m_f_o:.4f}  serbest eğim={m_f_s:.4f}±"
          f"{se_m_f_s:.4f} kesişim={a_f_s:+.4f}±{se_a_f_s:.4f}")
    print(f"  f_∞ = eğim(orijin) × Δz_∞({dz_inf:.4f}) = {f_inf:.4f}±"
          f"{f_inf_se:.4f}  ⇒ kalıntı payı = {1-f_inf:.4f}±{f_inf_se:.4f}")
    print(f"Doğrusallık g_D vs Δz_D (4 nokta): orijinden eğim="
          f"{m_g_o:.4f}±{se_m_g_o:.4f}  serbest eğim={m_g_s:.4f}±"
          f"{se_m_g_s:.4f} kesişim={a_g_s:+.4f}±{se_a_g_s:.4f}")
    print(f"  g_∞ = {g_inf:.4f}±{g_inf_se:.4f}")

    # ---------------- 189 A/B modelleriyle yan yana ----------------
    h189 = json.load(open(S189 / "HUKUM_189.json"))
    A_Pder = h189["K3"]["D_sonsuz"]["A_P_der"][8]      # HAVUZ
    A_Pder_se = h189["K3"]["D_sonsuz"]["A_P_der_se"][8]
    B_Pder = h189["K3"]["D_sonsuz"]["B_P_der"][8]
    B_Pder_se = h189["K3"]["D_sonsuz"]["B_P_der_se"][8]
    k3["189_model_karsilastirma"] = dict(
        A_P_der_havuz=A_Pder, A_P_der_se_havuz=A_Pder_se,
        B_P_der_havuz=B_Pder, B_P_der_se_havuz=B_Pder_se,
        f_inf_191=f_inf, f_inf_191_se=f_inf_se)
    print(f"\n189 model A (P_der HAVUZ)={A_Pder:.3f}±{A_Pder_se:.3f}  "
          f"189 model B (P_der HAVUZ)={B_Pder:.3f}±{B_Pder_se:.3f}  "
          f"191 ampirik yasa f_∞={f_inf:.3f}±{f_inf_se:.3f}")

    out["K3"] = k3
    json.dump(out, open(S191 / "HUKUM_191.json", "w"), indent=1,
              ensure_ascii=False, default=float)
    print(f"\n-> {S191/'HUKUM_191.json'}  BİTTİ", flush=True)
