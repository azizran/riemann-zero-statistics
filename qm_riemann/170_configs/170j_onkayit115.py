"""
170 — EK KAPI (K2'), ÖN KAYIT: λ = 1.15 (L115) — λ EKSENİNİN DERİN UCU
======================================================================
`c`'yi HESAPLAMAZ, korelatöre (J₂/KALİB_u2) BAKMAZ. Bant seçimi ve
`W_X` tanımı `170c_onkayit60.bant_c` ile AYNI nesneden gelir (import),
ölçüm zinciri parçaları `165/163/166`'dan import edilir.

SORU: λ ≤ 0.70'te `c` bir TABANA oturdu (0.3690 / 0.3689) ve λ = 1.00'de
0.4035. λ'yı 1.00'ın ÖTESİNE (1.15) taşırsak ne olur?

 Q1  H-D1 (T-yasası)      : σ_Ĉ büyür ⇒ Π_b T(σ_b) küçülür ⇒ c DÜŞER
 Q2  H-C1                 : c = 4/π² (λ'dan bağımsız)
 Q3  "derin limit ALTTAN" : c λ ile artıp 4/π²'ye yaklaşır (üstünü
                            geçmez)  ⇒  c(1.15) ∈ (0.4035, 0.4053]
 Q4  λ-kuadratik (4 nokta, EMPİRİK) : bant bant ikinci derece uyum
 Q5  W_X^ν üyesi (4 λ gazı, EMPİRİK): c = A·W_X^{ν−1}
 Q6  KALİB λ-değişmez (168 §A3, ölü): c = c(Hk)·W_X(Hk)/W_X(1.15)

**Q1 ile Q3/Q4/Q5 İŞARETTE ayrışır** (Q1 aşağı, ötekiler yukarı).

ÖN-MÜHÜR: yalnız gazın marjinallerine bakılır; `C_L115.json` ölçümü bu
betikten SONRA (`170k_olcum115.py`) alınır.

Çıktı: scratchpad/170/ONKAYIT_L115.json
"""
import importlib
import json
import sys
import time
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
for _p in ("170_configs", "169_configs", "167_configs", "166_configs",
           "165_configs", "163_configs", "160_configs", "159_configs",
           "155_configs", "154_configs"):
    sys.path.insert(0, str(QM / _p))
K = importlib.import_module("165_cekirdek")
C163 = importlib.import_module("163_cekirdek")
B166 = importlib.import_module("166_bacak")
ORT = importlib.import_module("167_ortak")
TB = importlib.import_module("170b_T_yasasi")
OK60 = importlib.import_module("170c_onkayit60")     # bant_c + tanımlar

TWO_PI = 2 * np.pi
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
LO_MIN, LO_MAX = 0.52, 0.68
C_HIP = 4.0 / np.pi ** 2
LAM4 = (("Hkeskin", 1.00), ("L085", 0.85), ("L070", 0.70), ("L060", 0.60))


def main(veri="L115", lam=1.15, taban=0.40):
    t0 = time.time()
    K.PENCERE.update(ORT.pencere_dict())
    Y = K.gaz(veri, taban, 4000)
    print(f"=== 170 K2' ÖN KAYIT — {veri} (λ={lam}) ===")
    print(f"  N={len(Y.m0)} L={Y.L:.5f} σ_ds={np.std(Y.ds):.5f} "
          f"σ_X̃={np.std(Y.Xtil0):.5f}", flush=True)
    Mo = K.Model165(veri, taban, 4000, 0.95, "olculen", Y=Y, tau_c=0.95)
    Mo.sec("olculen", 0.95)
    m = Mo.msk
    tl = Mo.M["tau"][m]
    wh, wy = np.abs(Mo.hp[m]), np.abs(Mo.y[m])
    tE = float((tl * wh).sum() / wh.sum())
    tX = float((tl * wy).sum() / wy.sum())
    sC = Mo.sigC
    print(f"  σ_Ĉ={sC:.5f}   ⟨τ⟩_E={tE:.4f}   ⟨τ⟩_X={tX:.4f}", flush=True)

    ban = C163.bant_adaylari(Y, K.IZGARA_T1)
    hedef = [b for b in ban if LO_MIN - 1e-9 <= b["lo"] <= LO_MAX + 1e-9]
    tauc = np.array([r["tau"] for b in hedef for r in b["cizgi"]])
    WX = B166.karakteristik(Y.Xtil0, -TWO_PI * tauc).real

    REF = {a: OK60.bant_c(a)[0] for a, _ in LAM4}
    KT = json.load(open(SCR / "170/K1_T.json"))["gaz"]
    OK6 = json.load(open(SCR / "170/ONKAYIT_L060.json"))

    def PiT(tauQ, te, tx, sc):
        return float(np.prod([TB.T_sigma(TWO_PI * t * sc)
                              for t in (tauQ, te, tx, tx)]))

    # Q5 uyumu: dört λ gazının BÜTÜN bantları, log KALİB ↔ log W_X
    lx, ly = [], []
    for a, _ in LAM4:
        for lo, (kal, wx, c) in REF[a].items():
            lx.append(np.log(wx))
            ly.append(np.log(kal))
    p5 = np.polyfit(lx, ly, 1)
    nu, lnA = float(p5[0]), float(p5[1])
    print(f"  [Q5 uyumu] log KALİB = {nu:+.4f} log W_X + {lnA:.4f}  "
          f"(artık rms {np.std(np.array(ly)-np.polyval(p5,lx)):.4f})")

    lams = np.array([l for _, l in LAM4])
    V = np.vstack([np.ones(4), lams, lams ** 2]).T
    print("\n  bant lo  n  W_X(L115) | c(Hk)  c(.85) c(.70) c(.60) | "
          "Q1     Q3     Q4     Q5     Q6")
    i, rows = 0, []
    for b in hedef:
        nb = len(b["cizgi"])
        sl = slice(i, i + nb)
        gp = np.array([r["gp"] for r in b["cizgi"]])
        wx = float((WX[sl] * (gp / gp.sum())).sum())
        i += nb
        lo, tQ = round(b["lo"], 4), b["tau"]
        cs = np.array([REF[a][lo][2] for a, _ in LAM4])
        q1 = cs[0] * PiT(tQ, tE, tX, sC) / PiT(
            tQ, KT["Hkeskin"]["tau_E"], KT["Hkeskin"]["tau_X"],
            KT["Hkeskin"]["sigC"])
        q3 = C_HIP
        co, *_ = np.linalg.lstsq(V, cs, rcond=None)
        q4 = float(co[0] + co[1] * 1.15 + co[2] * 1.3225)
        q5 = float(np.exp(lnA) * wx ** (nu - 1.0))
        q6 = cs[0] * REF["Hkeskin"][lo][1] / wx
        rows.append(dict(lo=lo, n=nb, WX=wx, cs=list(cs), Q1=q1, Q3=q3,
                         Q4=q4, Q5=q5, Q6=q6))
        print(f"  {lo:.2f}   {nb:4d}  {wx:.4f}  | " +
              " ".join(f"{x:.4f}" for x in cs) +
              f" | {q1:.4f} {q3:.4f} {q4:.4f} {q5:.4f} {q6:.4f}")

    G = {k: float(np.exp(np.mean(np.log([r[k] for r in rows]))))
         for k in ("Q1", "Q3", "Q4", "Q5", "Q6")}
    G["Q2"] = C_HIP
    print("\n  " + "*" * 72)
    print("  ** ÖN KAYIT — λ = 1.15 (ÖLÇÜMDEN ÖNCE YAZILDI) **")
    print(f"     Q1  H-D1 (T-yasası, oranlı)      : c(L115) = {G['Q1']:.4f}"
          "   ← TEK 'AŞAĞI' diyen")
    print(f"     Q2  H-C1 (c = 4/π²)              : c(L115) = {G['Q2']:.4f}")
    print(f"     Q3  derin limite ALTTAN yaklaşma : c(L115) ∈ "
          f"(0.4035, {C_HIP:.4f}]")
    print(f"     Q4  λ-kuadratik (4 nokta, EMPİRİK): c(L115) = {G['Q4']:.4f}")
    print(f"     Q5  W_X^ν üyesi (ν={nu:.3f}, EMPİRİK): c(L115) = "
          f"{G['Q5']:.4f}")
    print(f"     Q6  KALİB λ-değişmez (168 §A3)    : c(L115) = {G['Q6']:.4f}")
    print(f"     çapa: c(1.00) = 0.4035 ± 0.0018")
    print("  " + "*" * 72)

    out = dict(veri=veri, lam=lam, sigds=float(np.std(Y.ds)),
               sigX=float(np.std(Y.Xtil0)), sigC=sC, tau_E=tE, tau_X=tX,
               nu=nu, A=float(np.exp(lnA)), bant=rows, ongoru=G,
               C_HIP=C_HIP, sure_s=time.time() - t0)
    p = SCR / "170/ONKAYIT_L115.json"
    p.write_text(json.dumps(out, indent=1, default=float))
    print(f"\n-> {p}  ({(time.time()-t0)/60:.1f} dk)")


if __name__ == "__main__":
    main()
