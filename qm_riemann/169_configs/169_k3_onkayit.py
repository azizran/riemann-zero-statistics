"""
169 — K3(b) ÖN KAYIT: L070 için c ÖNGÖRÜSÜ (ÖLÇÜMDEN ÖNCE)
==========================================================
Bu script `c`'yi HESAPLAMAZ. Yalnız yeni gazın (λ=0.70) DW çarpanı
W_X(τ)'yu — ki bu gazın X̃ MARJİNALİNİN bir özelliğidir, korelatörün
değil — 167_olcum.py ile BİREBİR aynı tanımla ölçer:

    W_X(τ_q) = ⟨e^{−2πiτ_q X̃0}⟩       (166_bacak.karakteristik)
    bant değeri = aday çizgilerde gp-ağırlıklı ortalama

ve iki RAKİP öngörüyü sayısallaştırır:

  P_A  (168 §A3'ün ölçülmüş yasası: KALİB_u2 λ-DEĞİŞMEZ)
       c(L070, bant) = c(Hkeskin, bant) · W_X(Hkeskin,bant)/W_X(L070,bant)
  P_B  (H-C1: c EVRENSEL doymuş sabit = 4/π²)
       c(L070, bant) = c(Hkeskin, bant)      (yani ≈ 4/π², λ'dan bağımsız)

Kullanım: 169_k3_onkayit.py [gaz=L070]
Çıktı:    scratchpad/169/ONKAYIT_<gaz>.json
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
B166 = importlib.import_module("166_bacak")
ORT = importlib.import_module("167_ortak")

TWO_PI = 2 * np.pi
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
LO_MIN, LO_MAX = 0.52, 0.68
C_HIP = 4.0 / np.pi ** 2


def main(veri="L070", taban=0.40):
    t0 = time.time()
    K.PENCERE.update(ORT.pencere_dict())
    Y = K.gaz(veri, taban, 4000)
    print(f"=== 169 K3(b) ÖN KAYIT — {veri} (λ={ORT.lam_of(veri)}) ===")
    print(f"  N={len(Y.m0)} L={Y.L:.5f} σ_ds={np.std(Y.ds):.5f} "
          f"σ_X̃={np.std(Y.Xtil0):.5f}", flush=True)
    ban = C163.bant_adaylari(Y, K.IZGARA_T1)
    hedef = [b for b in ban if LO_MIN - 1e-9 <= b["lo"] <= LO_MAX + 1e-9]
    tauc = np.array([r["tau"] for b in hedef for r in b["cizgi"]])
    WX = B166.karakteristik(Y.Xtil0, -TWO_PI * tauc).real
    Wa = B166.karakteristik(Y.ds, np.pi * tauc).real

    ref = json.load(open(SCR / "169/K1.json"))
    D0 = json.load(open(SCR / "167/C_Hkeskin.json"))
    B0 = {round(b["lo"], 4): b for b in D0["bant"] if b.get("olculdu")}

    print("\n  bant lo   n_çizgi   W_X(L070)  W_amp(L070) | W_X(Hk)  "
          "c(Hk)   | P_A öngörü   P_B öngörü")
    i, rows = 0, []
    for b in hedef:
        nb = len(b["cizgi"])
        sl = slice(i, i + nb)
        gp = np.array([r["gp"] for r in b["cizgi"]])
        gw = gp / gp.sum()
        wx = float((WX[sl] * gw).sum())
        wa = float((Wa[sl] * gw).sum())
        i += nb
        b0 = B0[round(b["lo"], 4)]
        c0 = b0["KALIB_u2"] / b0["W_X"]
        pA = c0 * b0["W_X"] / wx
        rows.append(dict(lo=b["lo"], tau=b["tau"], n=nb, WX=wx, Wamp=wa,
                         WX0=b0["W_X"], c0=c0, P_A=pA, P_B=c0))
        print(f"  {b['lo']:.2f}      {nb:4d}    {wx:.4f}     {wa:.4f}    | "
              f"{b0['W_X']:.4f}  {c0:.4f}  |  {pA:.4f}      {c0:.4f}")

    gA = float(np.exp(np.mean(np.log([r["P_A"] for r in rows]))))
    gB = float(np.exp(np.mean(np.log([r["P_B"] for r in rows]))))
    print(f"\n  ** ÖN KAYIT (ölçümden ÖNCE yazıldı) **")
    print(f"     P_A  (KALİB λ-değişmez, 168 §A3)     : c(L070) = {gA:.4f}")
    print(f"     P_B  (H-C1: c evrensel = 4/π²)       : c(L070) = {gB:.4f}"
          f"   [4/π² = {C_HIP:.4f}]")
    print(f"     ayrım = {100*(gB/gA-1):+.1f}%  ;  beklenen ölçüm hatası "
          f"≈ ±0.003 ⇒ ayrım ≈ {(gB-gA)/0.003:.0f}σ")
    p = SCR / f"169/ONKAYIT_{veri}.json"
    p.write_text(json.dumps(dict(veri=veri, bant=rows, P_A=gA, P_B=gB,
                                 C_HIP=C_HIP, sigX=float(np.std(Y.Xtil0)),
                                 sigds=float(np.std(Y.ds)),
                                 sure_s=time.time() - t0), indent=1))
    print(f"\n-> {p}  ({(time.time()-t0)/60:.1f} dk)")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "L070")
