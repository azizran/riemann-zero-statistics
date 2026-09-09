# -*- coding: utf-8 -*-
"""
184c — K1 bant defteri + K2 HİPOTEZ YARIŞI
===========================================
K0/ONKAYIT_184'te donan bant ızgarası, jackknife ve üç hipotez formuyla:
  • bant başına w_g, w_ikiz(twin), r(τ)=â_g/â_twin, D=w_g−w_twin + jackknife se
  • r(τ)'ye H-Z1 (erfc), H-Z2 (güç), H-Z3 (Gauss-DW) fit + ölüm eşikleri
Birincil ikiz = Hkeskin (sadakatli keskin, τ_c→∞). Ek: HA4 (erfc 0.68),
keskin (152 kaba — ön-kayıt w≈1 kontrolünü GEÇEMEYEN; kayıt için).

jackknife: c_q^(−k)=2((RE−RE_k)+i(IM−IM_k))/(N−n_k); â bandı yeniden kurulur;
se²=(K−1)/K Σ(x_(−k)−x̄)². Oranlar bağımsız gazlarda hata yayılımıyla.
Kullanım: 184c_K2_yarisma.py
"""
import json
import math
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S184 = SCR / "184"
ONK = json.load(open(S184 / "ONKAYIT_184.json"))
KENAR = ONK["kenar"]
KUYRUK = ONK["kuyruk_tau"]
ESIK = ONK["esik"]


def yukle(ad):
    d = np.load(S184 / f"K1_{ad}.npz")
    return {k: d[k] for k in d.files}


def band_ahat(D):
    """Her gaz için: tam â ve 8 jackknife-blok leave-one-out â (nline)."""
    re = D["re_blok"]; im = D["im_blok"]; nb = D["nb"]; N = int(D["N"])
    re_t = re.sum(0); im_t = im.sum(0)
    c = 2.0 * (re_t + 1j * im_t) / N
    ahat = np.abs(c)
    K = re.shape[0]
    ah_jk = np.empty((K, len(ahat)))
    for k in range(K):
        ck = 2.0 * ((re_t - re[k]) + 1j * (im_t - im[k])) / (N - nb[k])
        ah_jk[k] = np.abs(ck)
    return ahat, ah_jk


def agg(ahat, ah_jk, msk):
    """Bir bandın Σâ ve jackknife se'si."""
    s = float(ahat[msk].sum())
    sk = ah_jk[:, msk].sum(1)
    K = len(sk)
    se = math.sqrt((K - 1) / K * np.sum((sk - sk.mean()) ** 2))
    return s, se


def bant_defteri(g, tw, ad_tw):
    """g (gerçek) ve tw (twin) sözlükleri; bant başına w/r/D + se."""
    tau = g["tau"]; aeff = g["aq_eff"]
    ah_g, jk_g = band_ahat(g)
    ah_t, jk_t = band_ahat(tw)
    satir = []
    edges = KENAR + ["KUYRUK"]
    for i in range(len(KENAR) - 1):
        lo, hi = KENAR[i], KENAR[i + 1]
        msk = (tau > lo) & (tau <= hi)
        satir.append(_satir(lo, hi, msk, tau, aeff, ah_g, jk_g, ah_t, jk_t))
    # kuyruk toplu bandı τ>KUYRUK
    msk = tau > KUYRUK
    satir.append(_satir(KUYRUK, KENAR[-1], msk, tau, aeff, ah_g, jk_g,
                        ah_t, jk_t, etiket="KUYRUK>0.70"))
    return satir


def _satir(lo, hi, msk, tau, aeff, ah_g, jk_g, ah_t, jk_t, etiket=None):
    n = int(msk.sum())
    tb = float(np.average(tau[msk], weights=aeff[msk])) if n else float("nan")
    Sg, seg = agg(ah_g, jk_g, msk)
    St, set_ = agg(ah_t, jk_t, msk)
    Se = float(aeff[msk].sum())
    wg = Sg / Se; wt = St / Se
    se_wg = seg / Se; se_wt = set_ / Se
    r = Sg / St
    se_r = r * math.sqrt((seg / Sg) ** 2 + (set_ / St) ** 2)
    D = wg - wt
    se_D = math.sqrt(se_wg ** 2 + se_wt ** 2)
    return dict(lo=lo, hi=hi, etiket=etiket or f"{lo:.2f}-{hi:.2f}", n=n,
                tau=tb, w_g=wg, se_w_g=se_wg, w_tw=wt, se_w_tw=se_wt,
                r=r, se_r=se_r, D=D, se_D=se_D)


# ---- hipotez formları ----
def H_Z1(tau, tau_c, delta):
    from scipy.special import erfc
    return 0.5 * erfc((tau - tau_c) / delta)


def H_Z2(tau, c, alpha):
    return np.clip(1.0 - c * tau ** alpha, 0.0, None)


def H_Z3(tau, sigma):
    return np.exp(-(2 * np.pi * tau) ** 2 * sigma ** 2 / 2.0)


def fit(fn, tau, y, se, p0, bounds):
    from scipy.optimize import curve_fit
    try:
        popt, pcov = curve_fit(fn, tau, y, p0=p0, sigma=se,
                               absolute_sigma=True, bounds=bounds, maxfev=20000)
    except Exception as e:
        return None, None, float("inf"), str(e)
    res = (y - fn(tau, *popt)) / se
    chi2 = float(np.sum(res ** 2))
    return popt, pcov, chi2, None


def kos():
    g = yukle("gercek")
    twins = {}
    for ad in ("Hkeskin", "HA4", "ikiz"):
        p = S184 / f"K1_{ad}.npz"
        if p.exists():
            twins[ad] = yukle(ad)
    out = dict(bant={}, fit={})
    for ad, tw in twins.items():
        sat = bant_defteri(g, tw, ad)
        out["bant"][ad] = sat
    # K2: r(τ) fit — BİRİNCİL twin Hkeskin
    birincil = "Hkeskin" if "Hkeskin" in twins else list(twins)[0]
    sat = out["bant"][birincil]
    fit_sat = [s for s in sat if s["etiket"].startswith(("0.4", "0.5", "0.6",
                                                          "0.7", "0.8"))
               and not s["etiket"].startswith("KUYRUK")]
    tau = np.array([s["tau"] for s in fit_sat])
    r = np.array([s["r"] for s in fit_sat])
    se = np.array([s["se_r"] for s in fit_sat])
    ndof_base = len(tau)
    yaris = {}
    for nm, fn, p0, bnds, npar in (
        ("H-Z1_erfc", H_Z1, [0.75, 0.15], ([0.4, 0.01], [3.0, 2.0]), 2),
        ("H-Z2_guc", H_Z2, [0.3, 3.0], ([0.0, 0.1], [10.0, 20.0]), 2),
        ("H-Z3_gaussDW", H_Z3, [0.05], ([1e-3], [2.0]), 1),
    ):
        popt, pcov, chi2, err = fit(fn, tau, r, se, p0, bnds)
        dof = max(ndof_base - npar, 1)
        c2d = chi2 / dof if np.isfinite(chi2) else float("inf")
        # kuyruk sınavı (τ>0.70 bantları)
        kmask = tau > KUYRUK
        if popt is not None and kmask.any():
            resk = (r[kmask] - fn(tau[kmask], *popt)) / se[kmask]
            kz = float(abs(resk.mean()) / (1.0 / math.sqrt(len(resk))) if len(resk) else 0)
            kz = float(abs(resk.mean()) * math.sqrt(len(resk)))
        else:
            kz = float("nan")
        yasiyor = c2d <= ESIK["survival_chi2_dof"]
        kuyruk_ok = (not np.isfinite(kz)) or (kz <= ESIK["kuyruk_z"])
        yaris[nm] = dict(par=None if popt is None else [float(x) for x in popt],
                         perr=None if pcov is None else
                         [float(x) for x in np.sqrt(np.diag(pcov))],
                         chi2=chi2, dof=dof, chi2_dof=c2d, kuyruk_z=kz,
                         yasiyor=bool(yasiyor), kuyruk_gecti=bool(kuyruk_ok),
                         hata=err)
    out["fit"] = dict(birincil_twin=birincil,
                      tau=[float(x) for x in tau],
                      r=[float(x) for x in r], se=[float(x) for x in se],
                      yaris=yaris)
    (S184 / "K2_sonuc.json").write_text(json.dumps(out, indent=1, default=float))
    _yaz(out, twins)
    return out


def _yaz(out, twins):
    print("=" * 78)
    print("184c — K1 BANT DEFTERİ + K2 YARIŞI")
    print("=" * 78)
    for ad in out["bant"]:
        print(f"\n### İKİZ = {ad}   (w_g/w_tw/r/D, genlik-ağırlıklı, ±jk se)")
        print(f"{'bant':>12} {'n':>5} {'τ':>6} {'w_g':>14} {'w_tw':>14} "
              f"{'r=g/tw':>14} {'D=g−tw':>14}")
        for s in out["bant"][ad]:
            print(f"{s['etiket']:>12} {s['n']:>5} {s['tau']:>6.3f} "
                  f"{s['w_g']:>7.3f}±{s['se_w_g']:.3f} "
                  f"{s['w_tw']:>7.3f}±{s['se_w_tw']:.3f} "
                  f"{s['r']:>7.3f}±{s['se_r']:.3f} "
                  f"{s['D']:>+7.3f}±{s['se_D']:.3f}")
    f = out["fit"]
    print(f"\n### K2 — r(τ) fit (birincil ikiz {f['birincil_twin']}, "
          f"n_bant={len(f['tau'])})")
    for nm, y in f["yaris"].items():
        par = y["par"]
        pstr = ", ".join(f"{v:.4f}" for v in par) if par else "—"
        print(f"  {nm:14s} par=[{pstr}]  χ²/dof={y['chi2_dof']:.3f}  "
              f"kuyruk_z={y['kuyruk_z']:.2f}  "
              f"{'YAŞIYOR' if y['yasiyor'] else 'ÖLDÜ'}"
              f"{'' if y['kuyruk_gecti'] else '  [KUYRUKTA ÖLÜR]'}")


if __name__ == "__main__":
    kos()
