"""
169 — K2(b) TAMAMLAYICI: BÖLÜNMÜŞ-MERDİVEN ile GERÇEK ÇOK-BACAKLI KIRPMA
========================================================================
169_k2.py'nin (X1X2) ve (E+X1X2) sütunları KULLANILAMAZ: aynı alanı iki
kez kırpmak `clip(X)·clip(X) = σ_X²` sabitini verir (rezonans yok olur).
Burada bacaklar GERÇEKTEN ayrılır: model çizgileri iki AYRIK yarıya
bölünür (τ'ya göre sıralı, bir atlamalı) ve

    X = X_a + X_b ,  X² ⊃ 2 X_a X_b        (ÇAPRAZ terim: iki BAĞIMSIZ bacak)

Sınav nesnesi  G = E · X_a · X_b  (Ç4 ile aynı yapı: üç alan bacağı +
taşıyıcı). Varyantlar `clip(F) = σ_F sgn(F)`:

    000 tam | 100 E | 010 X_a | 011 X_a+X_b | 110 E+X_a | 111 üçü de

ρ(varyant) = Pu2(varyant)/Pu2(tam), 166_T1'in bant birleştirmesiyle.
Gauss beklentisi bacak başına √(2/π) = 0.79788.

AYRICA: arcsine ARA-DEĞER eğrisi — her bant için model alanı YALNIZ o
bandın τ penceresindeki çizgilerden kurulur (E_b), ölçülen η ile
r = korr ve ⟨sgn·sgn⟩ ölçülür; Gauss çiftinde (2/π)arcsin r beklenir.

Kullanım: 169_k2b.py <gaz> [taban] [tau_c]
Çıktı:    scratchpad/169/K2b_<gaz>.json
"""
import importlib
import json
import sys
import time
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
for _p in ("167_configs", "166_configs", "165_configs", "163_configs",
           "160_configs", "159_configs", "155_configs", "154_configs"):
    sys.path.insert(0, str(QM / _p))
K = importlib.import_module("165_cekirdek")
C163 = importlib.import_module("163_cekirdek")
ORT = importlib.import_module("167_ortak")

TWO_PI = 2 * np.pi
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/169")
LO_MIN, LO_MAX = 0.52, 0.68
KMAX = 3
K.PENCERE.update(ORT.pencere_dict())


def main(veri="Hkeskin", taban=0.40, tau_c=0.95):
    t0 = time.time()
    SCR.mkdir(parents=True, exist_ok=True)
    Y = K.gaz(veri, taban, 4000)
    Mo = K.Model165(veri, taban, 4000, 0.95, "olculen", Y=Y, tau_c=tau_c)
    Mo.sec("olculen", tau_c).alanlar(kmax=KMAX)
    print(f"=== 169-K2b {veri} taban={taban} τ_c={tau_c} ===", flush=True)
    idx = np.nonzero(Mo.msk)[0]
    o = np.argsort(Mo.M["tau"][idx])
    ia, ib = idx[o[0::2]], idx[o[1::2]]
    w = Mo.M["w"]
    Xa = K.sentez(Mo.s, w[ia], Mo.y[ia])
    Xb = K.sentez(Mo.s, w[ib], Mo.y[ib])
    Xa -= Xa.mean()
    Xb -= Xb.mean()
    E = Mo.E
    print(f"  bölünmüş merdiven: |A|={len(ia)} |B|={len(ib)}  "
          f"σ_Xa={np.std(Xa):.5f} σ_Xb={np.std(Xb):.5f} "
          f"korr(Xa,Xb)={np.corrcoef(Xa,Xb)[0,1]:+.4f}", flush=True)

    def clip(F):
        return float(np.std(F)) * np.sign(F)

    Ec, Xac, Xbc = clip(E), clip(Xa), clip(Xb)
    for et, F, Fc in (("E", E, Ec), ("Xa", Xa, Xac), ("Xb", Xb, Xbc)):
        print(f"  bacak {et}: ⟨clip·F⟩/⟨F²⟩ = "
              f"{float(np.dot(Fc,F)/np.dot(F,F)):.5f}   "
              f"(Gauss √(2/π)={np.sqrt(2/np.pi):.5f})")

    VAR = {"000_tam": (E, Xa, Xb), "100_E": (Ec, Xa, Xb),
           "010_Xa": (E, Xac, Xb), "001_Xb": (E, Xa, Xbc),
           "011_XaXb": (E, Xac, Xbc), "110_EXa": (Ec, Xac, Xb),
           "111_hepsi": (Ec, Xac, Xbc)}
    NB = {"000_tam": 0, "100_E": 1, "010_Xa": 1, "001_Xb": 1,
          "011_XaXb": 2, "110_EXa": 2, "111_hepsi": 3}

    ban = C163.bant_adaylari(Y, K.IZGARA_T1)
    hedef = [b for b in ban if LO_MIN - 1e-9 <= b["lo"] <= LO_MAX + 1e-9]
    Wall, bas = [], []
    for b in hedef:
        bas.append(len(Wall))
        for r in b["cizgi"]:
            Wall.append(r["w"])
            Wall.append(r["w"] + r["gap"] / 2)
    Wall = np.array(Wall)
    tm = time.time()
    olc = [C163.olc_cizgi(Y, float(x), kmax=KMAX) for x in Wall]
    print(f"  ölçüm {time.time()-tm:.0f}s ({len(Wall)} frekans)", flush=True)
    tm = time.time()
    JJ = {ad: [j / 2.0 for j in K.tayf_s(Mo.s, [a * b1 * b2], Wall)][0]
          for ad, (a, b1, b2) in VAR.items()}
    print(f"  {len(VAR)} varyantın öngörüsü {time.time()-tm:.0f}s", flush=True)

    ad_sirali = [a for a in VAR if a != "000_tam"]
    print("\n  bant τ_eff | " + " ".join(f"ρ({a.split('_')[1]})"
                                         for a in ad_sirali))
    tab = []
    for bi, b in enumerate(hedef):
        rec = {}
        pN, pO, Ao, Af, tv = [], [], [], [], []
        U = {ad: ([], []) for ad in VAR}
        for li, r in enumerate(b["cizgi"]):
            k0 = bas[bi] + 2 * li
            on, of = olc[k0], olc[k0 + 1]
            pN.append(on["pow"])
            pO.append(of["pow"])
            Ao.append(on["A"])
            Af.append(of["A"])
            tv.append(r["tau"])
            for ad in VAR:
                U[ad][0].append(K.s_den_J(on["h"], JJ[ad][k0],
                                          on["rho_ort"])[1])
                U[ad][1].append(K.s_den_J(of["h"], JJ[ad][k0 + 1],
                                          of["rho_ort"])[1])
        pN, pO = np.array(pN), np.array(pO)
        Ao, Af, tv = np.array(Ao), np.array(Af), np.array(tv)
        dd = float((pN - pO).sum())
        tef = float((tv * (pN - pO)).sum() / dd)
        P = {ad: float((pN * Ao ** 2 * np.array(U[ad][0])
                        - pO * Af ** 2 * np.array(U[ad][1])).sum() / dd)
             for ad in VAR}
        rec = dict(tau_eff=tef, lo=b["lo"],
                   oran={ad: P[ad] / P["000_tam"] for ad in ad_sirali})
        tab.append(rec)
        print(f"  {tef:.4f}    | " +
              " ".join(f"{P[a]/P['000_tam']:.4f}   " for a in ad_sirali))
    g1 = np.sqrt(2 / np.pi)
    ort = {a: float(np.mean([r["oran"][a] for r in tab])) for a in ad_sirali}
    print("  ORTALAMA  | " + " ".join(f"{ort[a]:.4f}   " for a in ad_sirali))
    print(f"\n  Gauss beklentisi: 1 bacak {g1:.4f} | 2 bacak {g1**2:.4f} "
          f"(=2/π) | 3 bacak {g1**3:.4f} | 4 bacak {g1**4:.4f} (=4/π²)")
    print("  ölçülen bacak-sayısı ortalaması:")
    for nb in (1, 2, 3):
        v = [ort[a] for a in ad_sirali if NB[a] == nb]
        print(f"    {nb} bacak: {np.mean(v):.4f}   (Gauss {g1**nb:.4f}, "
              f"fark {100*(np.mean(v)/g1**nb-1):+.1f}%)")
    p1 = np.mean([ort[a] for a in ad_sirali if NB[a] == 1])
    print(f"    ÇARPIMSAL uzatma 4 bacak: {p1**4:.4f}   "
          f"(4/π² = {4/np.pi**2:.4f}, fark {100*(p1**4/(4/np.pi**2)-1):+.1f}%)")

    # ---------- arcsine ARA-DEĞER eğrisi: bant-sınırlı model alanları ----
    print("\n  --- arcsine ARA-DEĞER EĞRİSİ (bant-sınırlı model alanı ↔ η) ---")
    print("   bant lo  n_çizgi   r        ⟨sgn·sgn⟩   (2/π)arcsin r   fark%")
    e1 = Y.e1 - Y.e1.mean()
    se = np.sign(e1)
    ars = []
    for b in hedef:
        m = Mo.msk & (Mo.M["tau"] >= b["lo"]) & (Mo.M["tau"] < b["hi"])
        if m.sum() < 2:
            continue
        Eb = K.sentez(Mo.s, w[m], Mo.hp[m])
        Eb -= Eb.mean()
        r = float(np.corrcoef(Eb, e1)[0, 1])
        ss = float(np.mean(np.sign(Eb) * se))
        pr = float((2 / np.pi) * np.arcsin(r))
        ars.append(dict(lo=b["lo"], n=int(m.sum()), r=r, ss=ss, pred=pr))
        print(f"   {b['lo']:.2f}     {int(m.sum()):5d}   {r:+.5f}   {ss:+.5f}"
              f"     {pr:+.5f}      {100*(ss/pr-1):+6.2f}")

    out = dict(veri=veri, bant=tab, ortalama=ort, arcsine=ars,
               nleg=NB, sure_s=time.time() - t0)
    p = SCR / f"K2b_{veri}.json"
    p.write_text(json.dumps(out, indent=1, default=float))
    print(f"\n-> {p}  ({(time.time()-t0)/60:.1f} dk)", flush=True)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "Hkeskin",
         float(sys.argv[2]) if len(sys.argv) > 2 else 0.40,
         float(sys.argv[3]) if len(sys.argv) > 3 else 0.95)
