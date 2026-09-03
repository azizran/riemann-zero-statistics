"""
169 — K3(a): β = 0.2175 NEREDEN GELİYOR? (BOŞ ÇİZGİ deneyi)
============================================================
168 §A2.7'nin ÖLÇÜLEN yasası:  θ/θ₀ = 1 − β φ ,  φ = BOŞ ÇİZGİ KESRİ.
Kalemle türetilemedi (168 §6.1'in borcu). Burada TÜRETİM YERİNE yasanın
ADRESİ sınanır — tek ve keskin bir hipotezle:

  **HİPOTEZ B (boş-çizgi seyrelmesi).** Kesim gazında model hâlâ τ≤0.95'e
  kadar 8981 çizgi uydurur; bunların φ kadarı gazın merdiveninde YOKTUR
  ("boş çizgi") ve gürültüye uyar. Eğer θ'nın φ-bağımlılığı BU çizgilerden
  geliyorsa, model alanından boş çizgileri ÇIKARINCA φ-bağımlılığı da
  KAYBOLMALIDIR:      θ_gerçek-çizgi (K070) / θ(Hkeskin) → 1.

Ölçülenler (her gaz için, ortak pencerede):
  ρ_E = Var(E_boş)/Var(E) , ρ_X = Var(X_boş)/Var(X)      seyrelme payı
  Π   = Pu2(E_g X_g²)/Pu2(E X²)                          öngörünün payı
  g_cal^g = g_E^g (g_X^g)²   (yalnız gerçek çizgilerle regresyon)
  θ^g = (Mu2/Pu2^g)/g_cal^g   ↔   θ = (Mu2/Pu2)/g_cal

Kullanım: 169_k3_bos.py <gaz> [taban] [tau_c]
Çıktı:    scratchpad/169/BOS_<gaz>.json
"""
import importlib
import json
import math
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


def main(veri="K070", taban=0.40, tau_c=0.95):
    t0 = time.time()
    SCR.mkdir(parents=True, exist_ok=True)
    kun = ORT.KUNYE.get(veri, dict(lam=1.0, tau_ust=1.00, pen=None))
    tau_ust = kun.get("tau_ust") or 1.00
    pen = kun.get("pen")
    Y = K.gaz(veri, taban, 4000)
    Mo = K.Model165(veri, taban, 4000, 0.95, "olculen", Y=Y, tau_c=tau_c)
    Mo.sec("olculen", tau_c).alanlar(kmax=KMAX)
    tau = Mo.M["tau"]
    if pen is not None:
        wq = np.array([0.5 * math.erfc((t - pen[0]) / pen[1]) for t in tau])
    else:
        wq = (tau <= tau_ust + 1e-12).astype(float)
    ger = Mo.msk & (wq > 0.5)                # "GERÇEK" çizgi
    bos = Mo.msk & (wq <= 0.5)               # "BOŞ" çizgi
    phi = float(1 - wq[Mo.msk].mean())
    print(f"=== 169-K3(a) BOŞ-ÇİZGİ {veri} (τ_ust={tau_ust} pen={pen}) ===")
    print(f"  çizgi: toplam {int(Mo.msk.sum())}  gerçek {int(ger.sum())}  "
          f"boş {int(bos.sum())}   φ(ağırlıklı) = {phi:.4f}", flush=True)

    w = Mo.M["w"]
    Eg = K.sentez(Mo.s, w[ger], Mo.hp[ger])
    Xg = K.sentez(Mo.s, w[ger], Mo.y[ger])
    Xg -= Xg.mean()
    E, X = Mo.E, Mo.X
    rE = 1 - float(np.var(Eg) / np.var(E))
    rX = 1 - float(np.var(Xg) / np.var(X))
    e1 = Y.e1 - Y.e1.mean()
    x1 = Y.Xtil0
    gEg = float(np.dot(Eg, e1) / np.dot(Eg, Eg))
    gXg = float(np.dot(Xg, x1) / np.dot(Xg, Xg))
    gcg = gEg * gXg * gXg
    A = Mo.artik
    print(f"  Var payı: boş çizgiler E'nin %{100*rE:.1f}'ini, X'in "
          f"%{100*rX:.1f}'ini taşıyor")
    print(f"  g_E: {A['gE']:.4f} → {gEg:.4f} (gerçek)   g_X: {A['gX']:.4f} → "
          f"{gXg:.4f}   g_cal: {A['gcal']:.4f} → {gcg:.4f}", flush=True)

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
    Jt = [j / 2.0 for j in K.tayf_s(Mo.s, [E * X * X, Eg * Xg * Xg], Wall)]
    print(f"  öngörü {time.time()-tm:.0f}s", flush=True)

    print("\n  bant τ_eff |  Mu2       Pu2(tam)   Pu2(gerçek)   Π=oran   "
          "KALİB    KALİB^g")
    tab = []
    for bi, b in enumerate(hedef):
        pN, pO, Ao, Af, tv = [], [], [], [], []
        M, P0, P1 = ([], []), ([], []), ([], [])
        for li, r in enumerate(b["cizgi"]):
            k0 = bas[bi] + 2 * li
            for j, (o, kk) in enumerate(((olc[k0], k0), (olc[k0 + 1], k0 + 1))):
                M[j].append(o["u2"])
                P0[j].append(K.s_den_J(o["h"], Jt[0][kk], o["rho_ort"])[1])
                P1[j].append(K.s_den_J(o["h"], Jt[1][kk], o["rho_ort"])[1])
            pN.append(olc[k0]["pow"])
            pO.append(olc[k0 + 1]["pow"])
            Ao.append(olc[k0]["A"])
            Af.append(olc[k0 + 1]["A"])
            tv.append(r["tau"])
        pN, pO = np.array(pN), np.array(pO)
        Ao, Af, tv = np.array(Ao), np.array(Af), np.array(tv)
        dd = float((pN - pO).sum())
        tef = float((tv * (pN - pO)).sum() / dd)

        def AG(V):
            return float((pN * Ao ** 2 * np.array(V[0])
                          - pO * Af ** 2 * np.array(V[1])).sum() / dd)

        mu2, p0, p1 = AG(M), AG(P0), AG(P1)
        tab.append(dict(tau_eff=tef, lo=b["lo"], Mu2=mu2, Pu2=p0, Pu2g=p1,
                        Pi=p1 / p0, KALIB=mu2 / p0, KALIBg=mu2 / p1))
        print(f"  {tef:.4f}    | {mu2:+.5f}  {p0:+.5f}   {p1:+.5f}   "
              f"{p1/p0:.4f}   {mu2/p0:.4f}  {mu2/p1:.4f}")

    kal = float(np.exp(np.mean(np.log([r["KALIB"] for r in tab]))))
    kalg = float(np.exp(np.mean(np.log([r["KALIBg"] for r in tab]))))
    th = kal / A["gcal"]
    thg = kalg / gcg
    print(f"\n  KALİB(tam) = {kal:.4f}   θ = {th:.4f}")
    print(f"  KALİB(gerçek çizgi) = {kalg:.4f}   θ^g = {thg:.4f}")
    print(f"  Hkeskin referansı: θ₀ = 0.8884  (168 §A2.6, 169 K1.1'de "
          f"yeniden üretildi)")
    print(f"  ⇒ θ/θ₀ = {th/0.8884:.4f}   (yasa: 1−0.2175φ = "
          f"{1-0.2175*phi:.4f})")
    print(f"  ⇒ θ^g/θ₀ = {thg/0.8884:.4f}   "
          f"[HİPOTEZ B: bu 1.0000 olmalıydı]")
    out = dict(veri=veri, phi=phi, rE=rE, rX=rX, gE=A["gE"], gX=A["gX"],
               gcal=A["gcal"], gEg=gEg, gXg=gXg, gcalg=gcg,
               KALIB=kal, KALIBg=kalg, theta=th, theta_g=thg,
               nger=int(ger.sum()), nbos=int(bos.sum()), bant=tab,
               sure_s=time.time() - t0)
    p = SCR / f"BOS_{veri}.json"
    p.write_text(json.dumps(out, indent=1, default=float))
    print(f"\n-> {p}  ({(time.time()-t0)/60:.1f} dk)", flush=True)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "K070",
         float(sys.argv[2]) if len(sys.argv) > 2 else 0.40,
         float(sys.argv[3]) if len(sys.argv) > 3 else 0.95)
