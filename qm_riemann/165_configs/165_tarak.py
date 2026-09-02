"""
165 — DÜZGÜN TARAK KONTROLÜ: Ç3'ü Ç4'ten ayıran ölçü
=====================================================
Aynı çizgi genlikleriyle ama sitelerin DÜZGÜN olduğu bir tarakta
(s_n → s̄ + ḡ(n−n̄)) κ(Δ) Dirichlet çekirdeğidir: Δ = 0'da 1,
Δ = k·dres'te ÖZDEŞ 0, aradaki tepe yalnız |Δ| ≲ dres. Yani

    düzgün tarak toplamı  =  Ç1 + Ç2 + Ç3(|Δ| ≲ dres)
    gerçek tarak toplamı  =  yukarıdakiler + **Ç4** (tarak-rezonansı:
                             |Δ| ≈ ω_r'de κ'nın 𝒢-tepeleri)

163 §6d bu kontrolü ölçmüş ve `s2*`'ın ×10–×40 çöktüğünü görmüştü;
burada aynı kontrol 165'in TAM dörtlü toplamına uygulanır ve Ç4'ün
payını doğrudan verir.

Kullanım:  165_tarak.py <veri> [taban] [kaynak] [tau_c]
"""
import importlib
import json
import sys
import time
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
sys.path.insert(0, str(QM / "165_configs"))
sys.path.insert(0, str(QM / "163_configs"))
K = importlib.import_module("165_cekirdek")
C163 = importlib.import_module("163_cekirdek")

TWO_PI = 2 * np.pi
SCR = K.SCR


def kos(veri, taban=0.40, kaynak="olculen", tau_c=0.95):
    t0 = time.time()
    Y = K.gaz(veri, taban, 4000)
    Mo = K.Model165(veri, taban, 4000, 0.95, kaynak, Y=Y, tau_c=tau_c)
    Mo.sec(kaynak, tau_c)
    gbar = TWO_PI / Y.L
    N = len(Mo.s)
    su = Mo.s[0] + gbar * np.arange(N)
    print(f"=== 165-TARAK {veri} kaynak={kaynak} τ_c={tau_c} ===", flush=True)
    print(f"  gerçek s: [{Mo.s[0]:.1f},{Mo.s[-1]:.1f}] ; düzgün s: "
          f"[{su[0]:.1f},{su[-1]:.1f}]  ḡ={gbar:.5f}", flush=True)

    w = Mo.M["w"][Mo.msk]
    amp_h, amp_y = Mo.hp[Mo.msk], Mo.y[Mo.msk]
    ban = C163.bant_adaylari(Y, [b for b in K.IZGARA_T1 if b[0] >= 0.52 - 1e-9])
    Wall, bas = [], []
    for b in ban:
        bas.append(len(Wall))
        for r in b["cizgi"]:
            Wall.append(r["w"])
            Wall.append(r["w"] + r["gap"] / 2)
    Wall = np.array(Wall)

    J = {}
    for ad, sit in (("gercek", Mo.s), ("duzgun", su)):
        tm = time.time()
        E = K.sentez(sit, w, amp_h)
        X = K.sentez(sit, w, amp_y)
        X = X - X.mean()
        G = [E, E * X, E * X * X]
        J[ad] = [x / 2.0 for x in K.tayf_s(sit, G, Wall)]
        print(f"  [{ad}] Var(E)={np.var(E):.5f} Var(X)={np.var(X):.5f} "
              f"({time.time()-tm:.0f}s)", flush=True)

    olc = [C163.olc_cizgi(Y, float(x), kmax=2) for x in Wall]
    qidx = {int(q): i for i, q in enumerate(Mo.M["q"])}
    cikti = []
    print("   τ_eff    A²s2 ölç   GERÇEK tarak  DÜZGÜN tarak   Ç4 = G−D"
          "     Ç4/G    Ç1        u0(G)  u0(D)")
    for bi, b in enumerate(ban):
        L = b["cizgi"]
        pN = np.array([olc[bas[bi] + 2 * i]["pow"] for i in range(len(L))])
        pO = np.array([olc[bas[bi] + 2 * i + 1]["pow"] for i in range(len(L))])
        Ao = np.array([olc[bas[bi] + 2 * i]["A"] for i in range(len(L))])
        Af = np.array([olc[bas[bi] + 2 * i + 1]["A"] for i in range(len(L))])
        tv = np.array([r["tau"] for r in L])
        u = pN - pO
        payda = float(u.sum())

        def agg(ad, k):
            vN = np.array([K.s_den_J(olc[bas[bi] + 2 * i]["h"],
                                     J[ad][k][bas[bi] + 2 * i],
                                     olc[bas[bi] + 2 * i]["rho_ort"])[0]
                           for i in range(len(L))])
            vO = np.array([K.s_den_J(olc[bas[bi] + 2 * i + 1]["h"],
                                     J[ad][k][bas[bi] + 2 * i + 1],
                                     olc[bas[bi] + 2 * i + 1]["rho_ort"])[0]
                           for i in range(len(L))])
            return float((pN * Ao ** 2 * vN - pO * Af ** 2 * vO).sum() / payda)

        def agg_u0(ad):
            vN = np.array([K.s_den_J(olc[bas[bi] + 2 * i]["h"],
                                     J[ad][0][bas[bi] + 2 * i],
                                     olc[bas[bi] + 2 * i]["rho_ort"])[1]
                           for i in range(len(L))])
            vO = np.array([K.s_den_J(olc[bas[bi] + 2 * i + 1]["h"],
                                     J[ad][0][bas[bi] + 2 * i + 1],
                                     olc[bas[bi] + 2 * i + 1]["rho_ort"])[1]
                           for i in range(len(L))])
            return float((pN * vN - pO * vO).sum() / payda)

        m = float((pN * Ao ** 2 * np.array(
            [olc[bas[bi] + 2 * i]["s2"] for i in range(len(L))])
            - pO * Af ** 2 * np.array(
            [olc[bas[bi] + 2 * i + 1]["s2"] for i in range(len(L))])).sum()
            / payda)
        g = agg("gercek", 2)
        d = agg("duzgun", 2)
        c1 = float((pN * Ao ** 2 * np.array(
            [K.s_den_J(olc[bas[bi] + 2 * i]["h"],
                       Mo.kanal1(qidx[int(L[i]['q'])])[0],
                       olc[bas[bi] + 2 * i]["rho_ort"])[0]
             for i in range(len(L))])).sum() / payda)
        rec = dict(tau=b["tau"], tau_eff=float((tv * u).sum() / payda),
                   olc=m, gercek=g, duzgun=d, C4=g - d, C1=c1,
                   u0g=agg_u0("gercek"), u0d=agg_u0("duzgun"))
        cikti.append(rec)
        print(f"   {rec['tau_eff']:.4f} {m:+.6f} {g:+.6f}   {d:+.6f}   "
              f"{g-d:+.6f}  {(g-d)/g if g else float('nan'):+6.3f} "
              f"{c1:+.6f} {rec['u0g']:+.3f} {rec['u0d']:+.3f}", flush=True)

    out = dict(veri=veri, kaynak=kaynak, tau_c=tau_c, bant=cikti,
               sure_s=time.time() - t0)
    SCR.mkdir(parents=True, exist_ok=True)
    p = SCR / f"TAR_{veri}_{kaynak}_tc{tau_c}.json"
    p.write_text(json.dumps(out, indent=1))
    print(f"-> {p}  ({(time.time()-t0)/60:.1f} dk)", flush=True)


if __name__ == "__main__":
    kos(sys.argv[1],
        float(sys.argv[2]) if len(sys.argv) > 2 else 0.40,
        sys.argv[3] if len(sys.argv) > 3 else "olculen",
        float(sys.argv[4]) if len(sys.argv) > 4 else 0.95)
