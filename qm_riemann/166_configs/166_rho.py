"""
166 — ρ(τ_bant) SPEKTROSKOPİSİ: 144'ün YÖNTEMİ, ÜÇ GAZDA
========================================================
KALEM'in H-K2 adayı `kalib(τ_bant) = ρ(τ_bant)` (bant çizgisinin
tarak-sağkalımı). Bu script ρ'yu 144'ün BİREBİR estimatörüyle ölçer:

    ρ(bant) = Σ_çizgi [ |ĥ(w)|² − |ĥ(w+gap/2)|² ] / Σ_çizgi b_q²
    ĥ(W) = 2⟨η e^{−iW m_n}⟩ ,  b_q = 2 a_q sin(πτ_q)

İKİ KONVANSİYON (146'nın DERSİ — ρ TABAN-BAĞIMLIDIR, gizlenmez):
  • `t0.52/c720`  — 144/145/146'nın standart konvansiyonu (ρ = 0.38/0.27/
    0.16/0.08 orada ölçüldü). Yalnız lo ≥ 0.52 bantları ölçülebilir.
  • `t0.40/c4000` — 165'in KENDİ konvansiyonu (kalibrasyon orada doğdu).
    146 E2: aynı bantta ρ ~1.5× daha büyük ve eğim ~2.6 (6.0 değil).
Yarışta ikisi de koşar; hangisinin kullanıldığı her satırda yazılıdır.

İKİ BANT KÜMESİ:
  • `b144` — 144'ün bantları (0.525,0.60](0.60,0.70](0.70,0.78](0.78,0.84],
    tohum 11, 260 örnek: `son`/t0.52 için 0.381/0.266/0.159/0.080
    DOĞRULAMA hedefidir (V-R1).
  • `b165` — 163/165'in `bant_adaylari`'sı (IZGARA_T1, tohum 21, 220
    örnek, gap < 2.5·dres): kalib ile ρ AYNI ÇİZGİLER üzerinde ölçülür.
    Yarışın hakemi bu kümedir.

Kullanım:  166_rho.py <gaz>            (gaz ∈ son | Hkeskin | HA4)
Çıktı:     scratchpad/166/RHO_<gaz>.json
"""
import importlib
import json
import sys
import time
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
for _p in ("165_configs", "163_configs", "160_configs", "159_configs",
           "155_configs", "154_configs"):
    sys.path.insert(0, str(QM / _p))
K165 = importlib.import_module("165_cekirdek")
C163 = importlib.import_module("163_cekirdek")
K155 = importlib.import_module("155_cekirdek")
KOS155 = importlib.import_module("155_kos")
C154 = importlib.import_module("154_cekirdek")

TWO_PI = 2 * np.pi
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/166")
BANT144 = [(0.525, 0.60), (0.60, 0.70), (0.70, 0.78), (0.78, 0.84)]


class Shim:
    """bant_adaylari'nın istediği asgari arayüz (L, mid) — 165 ile aynı."""

    def __init__(self, L, mid):
        self.L, self.mid = L, mid


def proj(eta, mid, W, blok=25000, fblok=256):
    """ĥ(W) = 2⟨η e^{−iW m}⟩ — bloklu."""
    N = len(mid)
    W = np.asarray(W, float)
    acc = np.zeros(len(W), dtype=complex)
    for f0 in range(0, len(W), fblok):
        fs = slice(f0, min(f0 + fblok, len(W)))
        Wc = W[fs]
        for a in range(0, N, blok):
            sl = slice(a, min(a + blok, N))
            arg = np.outer(mid[sl], Wc)
            acc[fs] += eta[sl] @ np.cos(arg) - 1j * (eta[sl] @ np.sin(arg))
            del arg
    return 2 * acc / N


def cizgiler_144(L, mid, bantlar, tohum=11, nmax=260):
    """144'ün aday seçimi BİREBİR (allq = pk(exp(0.85·L)), gap<2.5·dres)."""
    allq = C154.pk(int(np.exp(0.85 * L)))
    allw = np.array([np.log(q) for q in allq])
    dres = TWO_PI / (mid[-1] - mid[0])
    rng = np.random.default_rng(tohum)
    qm = C154.pk_m(int(np.exp(0.85 * L)))
    out = []
    for lo, hi in bantlar:
        cand = [(q, w) for q, w in zip(allq, allw) if lo < w / L <= hi]
        if len(cand) > nmax:
            idx = rng.choice(len(cand), nmax, replace=False)
            cand = [cand[i] for i in idx]
        Ls, kul = [], 0
        for q, w in cand:
            j = np.searchsorted(allw, w)
            koms = [allw[k] for k in (j - 1, j + 1)
                    if 0 <= k < len(allw) and abs(allw[k] - w) > 1e-12]
            gap = min(abs(w - k) for k in koms)
            if gap < 2.5 * dres:
                continue
            aq = 1.0 / (np.pi * qm[q] * np.sqrt(q))
            Ls.append(dict(q=int(q), w=float(w), tau=float(w / L),
                           gap=float(gap), grup=kul % 8,
                           gp=float((2 * aq * np.sin(np.pi * w / L)) ** 2)))
            kul += 1
        out.append(dict(lo=lo, hi=hi, tau=round(0.5 * (lo + hi), 4),
                        kul=kul, cizgi=Ls))
    return out


def rho_bantlar(eta, mid, ban, njack=8):
    """Her bant için ρ = Σ(on−off)/Σ gp + 8-grup jackknife hatası."""
    W = []
    for b in ban:
        for r in b["cizgi"]:
            W.append(r["w"])
            W.append(r["w"] + r["gap"] / 2)
    if not W:
        return []
    H = proj(eta, mid, np.array(W))
    P = np.abs(H) ** 2
    out, k = [], 0
    for b in ban:
        n = len(b["cizgi"])
        on = P[k:k + 2 * n:2]
        off = P[k + 1:k + 2 * n:2]
        gp = np.array([r["gp"] for r in b["cizgi"]])
        grp = np.array([r["grup"] for r in b["cizgi"]])
        tau = np.array([r["tau"] for r in b["cizgi"]])
        k += 2 * n
        if n == 0 or gp.sum() == 0:
            out.append(dict(tau=b["tau"], lo=b["lo"], hi=b["hi"], kul=0,
                            rho=float("nan"), srho=float("nan")))
            continue
        rho = float((on - off).sum() / gp.sum())
        jk = []
        for g in range(njack):
            m = grp != g
            if not m.any() or gp[m].sum() == 0:
                continue
            jk.append(float((on[m] - off[m]).sum() / gp[m].sum()))
        jk = np.array(jk)
        s = (float(np.sqrt((len(jk) - 1) / len(jk) * np.sum((jk - jk.mean()) ** 2)))
             if len(jk) >= 4 else float("nan"))
        out.append(dict(tau=b["tau"], lo=b["lo"], hi=b["hi"], kul=int(n),
                        rho=rho, srho=s,
                        tau_agir=float((tau * gp).sum() / gp.sum()),
                        on=float(on.sum()), off=float(off.sum()),
                        onp=float(gp.sum())))
    return out


def kos(veri):
    t0 = time.time()
    z = KOS155.veri_yukle(veri)
    mid = 0.5 * (z[:-1] + z[1:])
    L = float(np.log(mid / TWO_PI).mean())
    print(f"=== 166-RHO {veri}  N={len(mid)} L={L:.5f} ===", flush=True)
    ban144 = cizgiler_144(L, mid, BANT144)
    ban165 = C163.bant_adaylari(Shim(L, mid), K165.IZGARA_T1)
    sonuc = {}
    for taban, cap in ((0.52, 720), (0.40, 4000)):
        Ç = K155.eta_onbellek(z, veri, taban, cap)
        eta = Ç["eta"]
        assert abs(Ç["L"] - L) < 1e-9, (Ç["L"], L)
        etiket = f"t{taban}_c{cap}"
        print(f"  [{etiket}] nq={Ç['nq']} σ_η²={Ç['s_eta']:.5f} "
              f"c₁={Ç['c1']:+.5f}", flush=True)
        for ad, ban in (("b144", ban144), ("b165", ban165)):
            # taban altındaki bantlar ÖLÇÜLEMEZ (regresyon çizgiyi sildi)
            kul = [b for b in ban if b["lo"] >= taban - 1e-9]
            tm = time.time()
            R = rho_bantlar(eta, mid, kul)
            sonuc[f"{etiket}_{ad}"] = R
            print(f"   {ad}: " + "  ".join(
                f"({b['lo']:.3f},{b['hi']:.2f}] n={b['kul']:3d} "
                f"ρ={b['rho']:+.4f}±{b['srho']:.4f}" for b in R)
                + f"   [{time.time()-tm:.0f}s]", flush=True)
    out = dict(veri=veri, L=L, N=int(len(mid)), kosum=sonuc,
               sure_s=time.time() - t0)
    SCR.mkdir(parents=True, exist_ok=True)
    p = SCR / f"RHO_{veri}.json"
    p.write_text(json.dumps(out, indent=1))
    print(f"-> {p}  ({(time.time()-t0)/60:.1f} dk)", flush=True)


if __name__ == "__main__":
    kos(sys.argv[1])
