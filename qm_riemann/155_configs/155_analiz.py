"""
155 — τ₀ ÇIKARIMI ve HÜKÜM TABLOSU
==================================
Girdi: 155_kos.py'nin JSON'ları (çizgi-bazında kayıtlı).

TEMEL FİKİR: φ_Γ ilkel gözlenebilirdir; R ondan türetilir. Bu yüzden
τ₀ ÖNCE φ'den, sonra (karşılaştırma için) R'den çıkarılır.

Her şey ÇİZGİ-BAZINDA kayıttan yeniden hesaplanır:
    Γ(bin) = Σ_i (cr_on,i − cr_off,i) / Σ_i (pow_on,i − pow_off,i)
Böylece
  * bant ızgarası ölçümden SONRA değiştirilebilir (bant-genişliği
    sistematiği ÖLÇÜLÜR, varsayılmaz),
  * jackknife (grup sil) BÜTÜN BANTLARDA eşzamanlı yapılabilir →
    τ₀'ın kendi hatası doğrudan replikalardan gelir (154'ün MC
    yayılımından daha az varsayımlı),
  * apsis olarak bant ortası (τ̄) yerine GÜÇ-AĞIRLIKLI τ_eff
    kullanılabilir — sıfır-geçişi ölçümünde bu ~0.001'lik gerçek
    bir kaymadır.
"""
import json
import sys
from pathlib import Path

import numpy as np

SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/155")
NJACK = 8


# ---------------------------------------------------------------- veri
def yukle(veri, taban, bset="ince"):
    p = SCR / f"tau0_{veri}_t{taban}_{bset}.json"
    if not p.exists():
        return None
    return json.loads(p.read_text())


def cizgiler(j, tmin=-1, tmax=9):
    """JSON'daki bütün bantların çizgilerini havuzla (bantlar τ'da ayrık)."""
    L = []
    for b in j["bantlar"]:
        if not b.get("olculdu"):
            continue
        for c in b["cizgi"]:
            if tmin <= c["tau"] <= tmax:
                L.append(c)
    L.sort(key=lambda c: c["tau"])
    return L


def bin_phi(L, kenar, sil=None):
    """Çizgi havuzunu kenarlara böl; her bin için (τ_eff, φ, |Γ|, n).

    sil: silinecek jackknife grubu (None = tam örneklem).
    """
    out = []
    for lo, hi in zip(kenar[:-1], kenar[1:]):
        S = [c for c in L if lo < c["tau"] <= hi
             and (sil is None or c["grup"] != sil)]
        if len(S) < 2:
            out.append(None)
            continue
        u = np.array([c["pow_on"] - c["pow_off"] for c in S])
        cr = np.array([complex(c["cr_on_re"] - c["cr_off_re"],
                               c["cr_on_im"] - c["cr_off_im"]) for c in S])
        tv = np.array([c["tau"] for c in S])
        gp = np.array([c["gp"] for c in S])
        d = float(u.sum())
        if abs(d) < 1e-300:
            out.append(None)
            continue
        G = complex(cr.sum() / d)
        out.append(dict(lo=lo, hi=hi, tbar=0.5 * (lo + hi),
                        teff=float((tv * u).sum() / d), n=len(S),
                        phi=float(np.angle(G)), absG=float(abs(G)),
                        rho=d / float(gp.sum())))
    return out


# ---------------------------------------------------- τ₀ tahmin edicisi
def fit_tau0(t, p, w, derece=1):
    """φ = Σ c_k τ^k ağırlıklı fit; kökü τ₀ (0.35–0.65 aralığında)."""
    t = np.asarray(t, float); p = np.asarray(p, float); w = np.asarray(w, float)
    X = np.vstack([t**k for k in range(derece + 1)]).T
    W = np.diag(w)
    beta = np.linalg.solve(X.T @ W @ X, X.T @ W @ p)
    kok = np.roots(beta[::-1])
    kok = [x.real for x in kok if abs(x.imag) < 1e-9 and 0.35 < x.real < 0.65]
    art = p - X @ beta
    chi2 = float(np.sum(w * art**2))
    if not kok:
        return float("nan"), beta, chi2, art
    # sıfıra en yakın kök (birden çok kök varsa: veri aralığına en yakın)
    kok = sorted(kok, key=lambda x: abs(x - t.mean()))
    return float(kok[0]), beta, chi2, art


def tau0_jk(L, kenar, apsis="teff", derece=1, njack=NJACK, tmin=None,
            tmax=None):
    """Tam örneklem τ₀ + grup-sil jackknife hatası (bütün bantlar birlikte).

    Ağırlık: her binin φ hatası jackknife'tan; fit ağırlığı 1/σ_φ².
    """
    tam = bin_phi(L, kenar)
    rep = [bin_phi(L, kenar, sil=k) for k in range(njack)]
    idx = [i for i, b in enumerate(tam)
           if b is not None and all(r[i] is not None for r in rep)]
    if tmin is not None:
        idx = [i for i in idx if tam[i]["tbar"] >= tmin - 1e-9]
    if tmax is not None:
        idx = [i for i in idx if tam[i]["tbar"] <= tmax + 1e-9]
    if len(idx) < derece + 1:      # derece+1 = tam interpolasyon (154 usulü)
        return None
    sph = []
    for i in idx:
        v = np.array([r[i]["phi"] for r in rep])
        sph.append(np.sqrt((len(v) - 1) / len(v) * np.sum((v - v.mean())**2)))
    sph = np.maximum(np.array(sph), 1e-4)
    w = 1.0 / sph**2

    def apsisler(bl):
        return np.array([bl[i][apsis if apsis in bl[i] else "tbar"]
                         for i in idx])

    t0, beta, chi2, art = fit_tau0(apsisler(tam),
                                   np.array([tam[i]["phi"] for i in idx]),
                                   w, derece)
    t0r = []
    for r in rep:
        v, _, _, _ = fit_tau0(apsisler(r),
                              np.array([r[i]["phi"] for i in idx]), w, derece)
        if np.isfinite(v):
            t0r.append(v)
    t0r = np.array(t0r)
    s = (float(np.sqrt((len(t0r) - 1) / len(t0r) * np.sum((t0r - t0r.mean())**2)))
         if len(t0r) >= 4 else float("nan"))
    dof = max(len(idx) - (derece + 1), 1)
    # χ² > dof ise jackknife (yalnız çizgi-seçimi) hatayı EKSİK tahmin
    # ediyordur; saçılmayla şişirilmiş hatayı da ver.
    s_sac = s * float(np.sqrt(max(chi2 / dof, 1.0)))
    return dict(tau0=t0, s_jk=s, s_sac=s_sac, chi2=float(chi2), dof=dof,
                n=len(idx), beta=[float(x) for x in beta],
                egim=float(beta[1]) if len(beta) > 1 else float("nan"),
                artik=[float(x) for x in art],
                tau=[float(x) for x in apsisler(tam)],
                phi=[float(tam[i]["phi"]) for i in idx],
                sphi=[float(x) for x in sph])


def tau0_R(j, derece=1):
    """R'nin KENDİ sıfırı — ölçülen R bantlarına parabol/doğru (154 usulü)."""
    B = [b for b in j["bantlar"] if b.get("olculdu")
         and b["artik"] < 0.02 and np.isfinite(b.get("sR_jk", np.nan))]
    B = sorted(B, key=lambda b: b["tau"])[:max(derece + 2, 4)]
    if len(B) < derece + 2 or not (min(b["R"] for b in B) < 0 <
                                   max(b["R"] for b in B)):
        return None
    t = np.array([b["tau_eff"] for b in B])
    R = np.array([b["R"] for b in B])
    w = 1.0 / np.array([max(b["sR_jk"], 1e-3) for b in B])**2
    Rjk = np.array([[b["R_jk"][k] for b in B] for k in range(NJACK)])
    t0, beta, chi2, art = fit_tau0(t, R, w, derece)
    rep = []
    for k in range(NJACK):
        if not np.all(np.isfinite(Rjk[k])):
            continue
        v, _, _, _ = fit_tau0(t, Rjk[k], w, derece)
        if np.isfinite(v):
            rep.append(v)
    rep = np.array(rep)
    s = (float(np.sqrt((len(rep) - 1) / len(rep)
                       * np.sum((rep - rep.mean())**2)))
         if len(rep) >= 4 else float("nan"))
    return dict(tau0=t0, s_jk=s, n=len(B),
                tau=[float(x) for x in t], R=[float(x) for x in R])
