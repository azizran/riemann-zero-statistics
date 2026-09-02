"""
165 — T3: b'NİN İŞARETİ HANGİ KANALDAN / HANGİ τ_q BÖLGESİNDEN?
================================================================
163 §6b, X̃'yi τ-dilimlerine ayırıp ÖLÇÜLEN `A²s2`'nin hücre matrisini
çıkarmıştı (işaret τ ≤ 0.4 çift-bloğunda doğuyordu). 165 aynı matrisi
İKİ KEZ kurar:

  (M) ÖLÇÜLEN:  A²⟨σ X^{(j)}X^{(k)}⟩/⟨ρ⟩ , η GERÇEK  (163 §6b ile aynı;
      kapanış Σ_{jk} = ölçülen A²s2 ÖZDEŞ — dilimler + artık TAM ayrışım)
  (P) BELİRLENİMLİ: aynı matris, η yerine model alanı E_mod (yani
      dörtlü toplamın kendisi); kapanış Σ_{jk} = öngörülen A²s2.

Ayrıca Ç1 kanalının τ_q yoğunluğu KAPALI biçimde bilindiğinden
(F(τ) = Re[h'_q ȳ_q], §T2) Ç1'in dilim dağılımı doğrudan verilir.

Kullanım:  165_T3.py <veri> [taban] [kaynak] [tau_c]
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
KOS = importlib.import_module("165_kos")

TWO_PI = 2 * np.pi
SCR = K.SCR
DILIM = [(0.0, 0.2), (0.2, 0.4), (0.4, 0.55), (0.55, 0.7), (0.7, 0.85),
         (0.85, 0.95)]


def kos(veri, taban=0.40, kaynak="olculen", tau_c=0.95, bantlar=None):
    t0 = time.time()
    bantlar = bantlar or [b for b in K.IZGARA_T1 if b[0] >= 0.52 - 1e-9]
    Y = K.gaz(veri, taban, 4000)
    Mo = K.Model165(veri, taban, 4000, 0.95, kaynak, Y=Y, tau_c=tau_c)
    Mo.sec(kaynak, tau_c).alanlar(kmax=2)
    print(f"=== 165-T3 {veri} kaynak={kaynak} τ_c={tau_c} ===", flush=True)
    print(f"  g_E={Mo.gE:.4f} g_X={Mo.gX:.4f} g={Mo.gcal:.4f}", flush=True)

    # --- X dilimleri --------------------------------------------------
    tau, w = Mo.M["tau"], Mo.M["w"]
    Xs = []
    for lo, hi in DILIM:
        m = (tau > lo) & (tau <= hi) & Mo.msk
        v = K.sentez(Mo.s, w[m], Mo.y[m]) if m.any() else np.zeros(len(Mo.s))
        Xs.append(v - v.mean())
    Xres = Y.Xtil0 - sum(Xs)
    ad = [f"{lo}-{hi}" for lo, hi in DILIM] + ["artık"]
    XM = Xs + [Xres]                      # ÖLÇÜLEN ayrışım (7 bileşen)
    XP = Xs                               # ÖNGÖRÜ ayrışımı (6 bileşen)
    print("  Var: " + "  ".join(f"{a}={np.var(x):.5f}" for a, x in
                                zip(ad, XM)), flush=True)

    E_olc = Y.e1 - Y.e1.mean()
    alan, etik = [], []
    for i in range(len(XM)):
        for j in range(i, len(XM)):
            alan.append(E_olc * XM[i] * XM[j])
            etik.append(("M", i, j))
    for i in range(len(XP)):
        for j in range(i, len(XP)):
            alan.append(Mo.E * XP[i] * XP[j])
            etik.append(("P", i, j))

    ban = C163.bant_adaylari(Y, bantlar)
    Wall, bas = [], []
    for b in ban:
        bas.append(len(Wall))
        for r in b["cizgi"]:
            Wall.append(r["w"])
            Wall.append(r["w"] + r["gap"] / 2)
    Wall = np.array(Wall)
    print(f"  {len(ban)} bant, {len(Wall)} frekans, {len(alan)} alan",
          flush=True)
    tm = time.time()
    J = [x / 2.0 for x in K.tayf_s(Mo.s, alan, Wall, fblok=128)]
    print(f"  tayf {time.time()-tm:.0f}s", flush=True)
    olc = [C163.olc_cizgi(Y, float(x), kmax=2) for x in Wall]

    cikti = []
    for bi, b in enumerate(ban):
        L = b["cizgi"]
        pN = np.array([olc[bas[bi] + 2 * li]["pow"] for li in range(len(L))])
        pO = np.array([olc[bas[bi] + 2 * li + 1]["pow"]
                       for li in range(len(L))])
        Ao = np.array([olc[bas[bi] + 2 * li]["A"] for li in range(len(L))])
        Af = np.array([olc[bas[bi] + 2 * li + 1]["A"] for li in range(len(L))])
        u = pN - pO
        payda = float(u.sum())
        tv = np.array([r["tau"] for r in L])
        MM = np.zeros((len(XM), len(XM)))
        PP = np.zeros((len(XP), len(XP)))
        for t, (tip, i, j) in enumerate(etik):
            vN = np.array([K.s_den_J(olc[bas[bi] + 2 * li]["h"],
                                     J[t][bas[bi] + 2 * li],
                                     olc[bas[bi] + 2 * li]["rho_ort"])[0]
                           for li in range(len(L))])
            vO = np.array([K.s_den_J(olc[bas[bi] + 2 * li + 1]["h"],
                                     J[t][bas[bi] + 2 * li + 1],
                                     olc[bas[bi] + 2 * li + 1]["rho_ort"])[0]
                           for li in range(len(L))])
            val = float((pN * Ao ** 2 * vN - pO * Af ** 2 * vO).sum() / payda)
            c = 1.0 if i == j else 2.0
            if tip == "M":
                MM[i, j] = c * val
            else:
                PP[i, j] = c * val
        rec = dict(tau=b["tau"], tau_eff=float((tv * u).sum() / payda),
                   olc_s2=float((pN * Ao ** 2 * np.array(
                       [olc[bas[bi] + 2 * li]["s2"] for li in range(len(L))])
                       - pO * Af ** 2 * np.array(
                       [olc[bas[bi] + 2 * li + 1]["s2"]
                        for li in range(len(L))])).sum() / payda),
                   M=MM.tolist(), P=PP.tolist(),
                   M_top=float(MM.sum()), P_top=float(PP.sum()), ad=ad)
        cikti.append(rec)
        print(f"  τe={rec['tau_eff']:.4f}  ÖLÇÜLEN A²s2={rec['olc_s2']:+.6f}"
              f"  matris kapanış M={rec['M_top']:+.6f} "
              f"(fark {rec['M_top']-rec['olc_s2']:+.1e})  P={rec['P_top']:+.6f}",
              flush=True)
        for tip, A in (("M", MM), ("P", PP)):
            print(f"    [{tip}] " + " | ".join(
                f"{ad[i]}×{ad[j]}:{A[i, j]:+.4f}"
                for i in range(A.shape[0]) for j in range(i, A.shape[0])
                if abs(A[i, j]) > 0.004), flush=True)

    out = dict(veri=veri, kaynak=kaynak, tau_c=tau_c, dilim=ad,
               gE=Mo.gE, gX=Mo.gX, gcal=Mo.gcal,
               varX=[float(np.var(x)) for x in XM], bant=cikti,
               sure_s=time.time() - t0)
    SCR.mkdir(parents=True, exist_ok=True)
    p = SCR / f"T3_{veri}_{kaynak}_tc{tau_c}.json"
    p.write_text(json.dumps(out, indent=1))
    print(f"-> {p}  ({(time.time()-t0)/60:.1f} dk)", flush=True)


if __name__ == "__main__":
    kos(sys.argv[1],
        float(sys.argv[2]) if len(sys.argv) > 2 else 0.40,
        sys.argv[3] if len(sys.argv) > 3 else "olculen",
        float(sys.argv[4]) if len(sys.argv) > 4 else 0.95)
