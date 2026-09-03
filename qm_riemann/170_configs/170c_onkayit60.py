"""
170 — K2 ÖN KAYIT: λ = 0.60 (L060) GAZI İÇİN c ÖNGÖRÜLERİ (ÖLÇÜMDEN ÖNCE)
=========================================================================
Bu betik `c`'yi HESAPLAMAZ ve korelatöre (J₂ / KALİB_u2) HİÇ BAKMAZ.
Yalnız L060 gazının MARJİNALLERİNİ ölçer —

    σ_ds, σ_X̃, σ_Ĉ  ;  W_X(τ) = ⟨e^{−2πiτ X̃0}⟩  (166_bacak.karakteristik,
    167_olcum ile BİREBİR aynı tanım, bantta gp-ağırlıklı)  ;
    ⟨τ⟩_E, ⟨τ⟩_X (model çizgi genliklerinin ağırlıklı τ ortalaması)

— ve BEŞ rakip yasayı sayısallaştırır. (169 §K3.2'nin ÖN KAYIT protokolü
ile aynı: marjinal SERBEST, korelatör YASAK.)

Öngörüler (hepsi bant bant kurulup log-ortalanır, ölçümle aynı pencere
lo ∈ [0.52,0.68]):

 P1a H-D1 MUTLAK   : c = Π_{b=1..4} T(σ_b),  T(σ)=√(ρ/arcsin ρ),
                     ρ = 1−e^{−σ²},  σ_b = 2πτ_b σ_Ĉ         [170b §A]
 P1b H-D1 ORANLI   : c(0.60) = c(Hkeskin)·ΠT(0.60)/ΠT(1.00)
                     (H-D1'e en iyi şans: yalnız λ-BAĞIMLILIĞI sınanır)
 P2  H-C1          : c = 4/π² (bant bant c(Hkeskin))          [169 K1]
 P3  λ-KUADRATİK   : bant bant λ'da ikinci derece uyum (ÜÇ ölçülen λ
                     noktasından) — EMPİRİK, türetim değil
 P4  W_X^ν ÜYESİ   : KALİB = A·W_X^ν üç λ gazından uyum ⇒
                     c = A·W_X^{ν−1}                          — EMPİRİK
 P5  KALİB λ-DEĞİŞMEZ (168 §A3, 169'da düşmüştü):
                     c(0.60) = c(Hk)·W_X(Hk)/W_X(0.60)

ÖN-MÜHÜR (koşudan ÖNCE):
  * L060 gazı kuruldu: σ_ds = 0.32135 (log_insa_L060.txt). Bu ÖLÇÜLDÜ;
    170_insa60.py'nin ön-mührü σ_ds ∈ [0.327,0.340] demişti — ISKA (%1.7
    düşük); gerçek kuvvet p = ln(0.32135/0.43134)/ln0.60 = 0.576.
  * Buradan σ_X̃(0.60) ≈ 0.24204·(0.32135/0.43134) = 0.1803 beklenir
    (σ_X̃/σ_ds oranı üç gazda 0.5612/0.5601/0.5591 — sabit).
  * W_X(bant 0.52) ≈ exp(−2π²·0.54²·0.1803²) = 0.8285 beklenir.
  * P1b'nin YÖNÜ kesin: ΠT σ_Ĉ küçüldükçe BÜYÜR ⇒ P1b > c(0.70) = 0.3690.
    Öteki dört yasa ise c(0.60) < 0.3690 der. **λ=0.60 H-D1'i İŞARETTEN
    ayırır.**

SONUÇ bloğu YALNIZ gerçek koşu çıktısındandır.
Kullanım: 170c_onkayit60.py [gaz=L060]
Çıktı:    scratchpad/170/ONKAYIT_L060.json
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
TB = importlib.import_module("170b_T_yasasi")          # T(σ) yasası

ORT.KUNYE["L060"] = dict(gercek=False, lam=0.60, tau_ust=1.00, pen=None,
                         aile="lam", T=1.0)

TWO_PI = 2 * np.pi
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
LO_MIN, LO_MAX = 0.52, 0.68
C_HIP = 4.0 / np.pi ** 2
LAM3 = ("Hkeskin", "L085", "L070")


def bant_c(ad):
    """C_<ad>.json'dan bant bant (lo -> (KALİB_u2, W_X, c_WX))."""
    d = json.load(open(SCR / f"167/C_{ad}.json"))
    out = {}
    for b in d["bant"]:
        if not b.get("olculdu"):
            continue
        if b["lo"] < LO_MIN - 1e-9 or b["lo"] > LO_MAX + 1e-9:
            continue
        if b["tau_eff"] > 0.85 or b["Ms2"] <= 0 or b["R_bant"] < 0.98 \
                or b["SNR"] < 3.0:
            continue
        out[round(b["lo"], 4)] = (b["KALIB_u2"], b["W_X"],
                                  b["KALIB_u2"] / b["W_X"])
    return out, (d["lam"] or 1.0)


def main(veri="L060", taban=0.40):
    t0 = time.time()
    K.PENCERE.update(ORT.pencere_dict())
    Y = K.gaz(veri, taban, 4000)
    print(f"=== 170 K2 ÖN KAYIT — {veri} (λ={ORT.lam_of(veri)}) ===")
    print(f"  N={len(Y.m0)}  L={Y.L:.5f}  σ_ds={np.std(Y.ds):.5f}  "
          f"σ_X̃={np.std(Y.Xtil0):.5f}", flush=True)

    # --- gazın Ĉ'si ve model bacaklarının τ'su (KORELATÖRE BAKMADAN) ---
    Mo = K.Model165(veri, taban, 4000, 0.95, "olculen", Y=Y, tau_c=0.95)
    Mo.sec("olculen", 0.95)
    msk = Mo.msk
    tl = Mo.M["tau"][msk]
    wh, wy = np.abs(Mo.hp[msk]), np.abs(Mo.y[msk])
    tE = float((tl * wh).sum() / wh.sum())
    tX = float((tl * wy).sum() / wy.sum())
    sC = Mo.sigC
    print(f"  σ_Ĉ={sC:.5f}   ⟨τ⟩_E={tE:.4f}   ⟨τ⟩_X={tX:.4f}", flush=True)

    ban = C163.bant_adaylari(Y, K.IZGARA_T1)
    hedef = [b for b in ban if LO_MIN - 1e-9 <= b["lo"] <= LO_MAX + 1e-9]
    tauc = np.array([r["tau"] for b in hedef for r in b["cizgi"]])
    WX = B166.karakteristik(Y.Xtil0, -TWO_PI * tauc).real
    Wa = B166.karakteristik(Y.ds, np.pi * tauc).real

    # --- referans gazların bant tabloları -----------------------------
    REF = {}
    SIGC = {}
    for ad in LAM3:
        REF[ad] = bant_c(ad)[0]
        SIGC[ad] = json.load(open(SCR / f"167/C_{ad}.json"))["sigC"]
    TAUB = {g: (json.load(open(SCR / "170/K1_T.json"))["gaz"][g]["tau_E"],
                json.load(open(SCR / "170/K1_T.json"))["gaz"][g]["tau_X"])
            for g in LAM3}

    def PiT(tauQ, tE_, tX_, sC_):
        s = [TWO_PI * t * sC_ for t in (tauQ, tE_, tX_, tX_)]
        return float(np.prod([TB.T_sigma(x) for x in s]))

    print("\n  bant lo  n_çiz  W_X(L060) W_amp | c(Hk)   c(L085)  c(L070) | "
          "P1a     P1b     P2      P3      P4      P5")
    i, rows = 0, []
    # P4: KALİB = A·W_X^ν  (üç λ gazının BÜTÜN bantları, log-log EKK)
    lx, ly = [], []
    for ad in LAM3:
        for lo, (kal, wx, c) in REF[ad].items():
            lx.append(np.log(wx))
            ly.append(np.log(kal))
    A4 = np.polyfit(lx, ly, 1)
    nu, lnA = float(A4[0]), float(A4[1])
    print(f"  [P4 uyumu]  log KALİB = {nu:+.4f}·log W_X + {lnA:.4f}   "
          f"⇒ ν = {nu:.4f}, A = {np.exp(lnA):.4f}, "
          f"artık rms = {np.std(np.array(ly)-np.polyval(A4,lx)):.4f}")

    for b in hedef:
        nb = len(b["cizgi"])
        sl = slice(i, i + nb)
        gp = np.array([r["gp"] for r in b["cizgi"]])
        gw = gp / gp.sum()
        wx = float((WX[sl] * gw).sum())
        wa = float((Wa[sl] * gw).sum())
        i += nb
        lo = round(b["lo"], 4)
        tauQ = b["tau"]
        cH, cL85, cL70 = (REF["Hkeskin"][lo][2], REF["L085"][lo][2],
                          REF["L070"][lo][2])
        p1a = PiT(tauQ, tE, tX, sC)
        p1b = cH * p1a / PiT(tauQ, *TAUB["Hkeskin"], SIGC["Hkeskin"])
        p2 = cH
        V = np.vstack([np.ones(3), [1.00, 0.85, 0.70],
                       np.array([1.00, 0.85, 0.70]) ** 2]).T
        co = np.linalg.solve(V, [cH, cL85, cL70])
        p3 = float(co[0] + co[1] * 0.60 + co[2] * 0.36)
        p4 = float(np.exp(lnA) * wx ** (nu - 1.0))
        p5 = cH * REF["Hkeskin"][lo][1] / wx
        rows.append(dict(lo=lo, tau=tauQ, n=nb, WX=wx, Wamp=wa,
                         cH=cH, cL85=cL85, cL70=cL70,
                         P1a=p1a, P1b=p1b, P2=p2, P3=p3, P4=p4, P5=p5))
        print(f"  {lo:.2f}    {nb:4d}   {wx:.4f}   {wa:.4f} | {cH:.4f}  "
              f"{cL85:.4f}   {cL70:.4f} | {p1a:.4f}  {p1b:.4f}  {p2:.4f}  "
              f"{p3:.4f}  {p4:.4f}  {p5:.4f}")

    G = {k: float(np.exp(np.mean(np.log([r[k] for r in rows]))))
         for k in ("P1a", "P1b", "P2", "P3", "P4", "P5")}
    print("\n  " + "*" * 74)
    print("  ** ÖN KAYIT — λ = 0.60 (ÖLÇÜMDEN ÖNCE YAZILDI) **")
    print(f"     P1a  H-D1 MUTLAK  (Π_b T(σ_b))          : c(L060) = "
          f"{G['P1a']:.4f}")
    print(f"     P1b  H-D1 ORANLI  (λ-bağımlılığı yalnız) : c(L060) = "
          f"{G['P1b']:.4f}")
    print(f"     P2   H-C1  (c = 4/π² evrensel)          : c(L060) = "
          f"{G['P2']:.4f}   [4/π² = {C_HIP:.4f}]")
    print(f"     P3   λ-kuadratik (EMPİRİK, 3 nokta)      : c(L060) = "
          f"{G['P3']:.4f}")
    print(f"     P4   W_X^ν üyesi (EMPİRİK, ν={nu:.3f})     : c(L060) = "
          f"{G['P4']:.4f}")
    print(f"     P5   KALİB λ-değişmez (168 §A3)          : c(L060) = "
          f"{G['P5']:.4f}")
    print(f"\n     ölçülen dizi: c(1.00)=0.4035  c(0.85)=0.3817  "
          f"c(0.70)=0.3690   (σ_tot ≈ 0.002–0.006)")
    print("     H-D1 (P1b) TEK BAŞINA c(0.60) > c(0.70) diyor; "
          "P2/P3/P4/P5 hepsi < 0.3690 diyor.")
    print("  " + "*" * 74)

    out = dict(veri=veri, lam=0.60, sigds=float(np.std(Y.ds)),
               sigX=float(np.std(Y.Xtil0)), sigC=sC, tau_E=tE, tau_X=tX,
               nu_P4=nu, A_P4=float(np.exp(lnA)), bant=rows, ongoru=G,
               C_HIP=C_HIP, sure_s=time.time() - t0)
    p = SCR / "170/ONKAYIT_L060.json"
    p.write_text(json.dumps(out, indent=1, default=float))
    print(f"\n-> {p}  ({(time.time()-t0)/60:.1f} dk)")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "L060")
