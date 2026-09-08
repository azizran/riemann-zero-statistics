# -*- coding: utf-8 -*-
"""
182g — K2: E-BACAĞININ ANATOMİSİ (ÇİZGİ / BANT AYRIŞIMI + SURROGATE'LAR)
=========================================================================
ÖN-MÜHÜR: bütün tanımlar ve ölüm/yaşama koşulları `182/ONKAYIT_182.json`
içinde, bu betik koşmadan ÖNCE donduruldu. Burada hiçbir eşik seçilmez.

MAKİNE KOPYALANMAZ: gaz/model/alan zinciri `165_cekirdek`ten AYNEN
kurulur — `K.gaz`, `K.Model165(...,'olculen',tau_c=0.95)`, `.alanlar()`,
bölünmüş merdiven (τ'ya göre sıralı, bir atlamalı: ia=o[0::2],
ib=o[1::2]), `K.sentez`, `K.tayf_s`, `C163.bant_adaylari` — yani
`169_k2b.main`'in ilk yarısının BİREBİR aynısı. ρ değerleri YENİDEN
hesaplanmaz; diskteki `169/K2b_<gaz>.json`'dan okunur.

ÖLÇÜLENLER
----------
(1) ÇİZGİ çarpanı  κ_F := ⟨clip(F)·F⟩/⟨F²⟩ = E|F|/σ_F   (F = E, Xa, Xb)
    + fazla basıklık kurt(F).  Gauss: κ = √(2/π), kurt = 0.
(2) BANT çarpanı   β_E := ρ(E)/κ_E   ve log-dilinde PAY_ÇİZGİ.
(3) τ eğimleri: ρ(E), ρ(Xa), ρ(Xb) ~ τ_eff doğrusal eğimleri (5 bant).
(4) BAĞIMSIZ KESTİRİMCİ (H-G3):
        ρ_J(var) := Σ_k Re[J_var(W_k)·conj(J_000(W_k))] / Σ_k |J_000|²
    — yalnız MODEL tarafı; ölçülen gazın hiçbir niceliği girmez.
(5) SURROGATE'LAR (H-G1 / H-G2):
    faz-rastgeleleştirme: amp_q → |amp_q|·e^{iψ_q}, ψ ~ U(0,2π) bağımsız.
    Güç tayfı BİT-BİT korunur, alan Gauss'laşır, bacaklar bağımsızlaşır.
      ρ°  : ÜÇ bacak da surrogate  → referansın kendi varsayımları altında
            sınavı (V1∧V2∧V3∧V4).
      ρ*  : YALNIZ E surrogate     → geri kazanım
            Rec := (ρ*(E) − ρ_J(E)) / (√(2/π) − ρ_J(E)).
(6) korr(X_a, X_b) — tam bağımsızlıkta 0 beklenir.
(7) R_η ve ÖZDEŞLİK Kov(η_artık, η)/Var(η) = 1 − R_η (176/K1_*.json).

Kullanım: 182g_anatomi.py VF1 VF2 ... (varsayılan VF1..VF8)
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

SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S182, S176, S169 = SCR / "182", SCR / "176", SCR / "169"
LO_MIN, LO_MAX = 0.52, 0.68           # 169_k2b ile AYNI
KMAX = 3
K.PENCERE.update(ORT.pencere_dict())
G1 = np.sqrt(2 / np.pi)


def kap(F):
    """κ_F = ⟨clip(F)·F⟩/⟨F²⟩, clip(F) = σ_F sgn(F) — 169_k2b'nin sayısı."""
    Fc = float(np.std(F)) * np.sign(F)
    return float(np.dot(Fc, F) / np.dot(F, F))


def kurt(F):
    d = F - F.mean()
    s = d.std()
    return float((d ** 4).mean() / s ** 4 - 3.0)


def egim(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    return float(np.polyfit(x, y, 1)[0])


def kos(ad, tohum):
    t0 = time.time()
    print("=" * 78, flush=True)
    print(f"182g — E-BACAĞI ANATOMİSİ: {ad}  (169_k2b zinciri aynen)",
          flush=True)
    print("=" * 78, flush=True)

    Y = K.gaz(ad, 0.40, 4000)
    Mo = K.Model165(ad, 0.40, 4000, 0.95, "olculen", Y=Y, tau_c=0.95)
    Mo.sec("olculen", 0.95).alanlar(kmax=KMAX)
    idx = np.nonzero(Mo.msk)[0]
    o = np.argsort(Mo.M["tau"][idx])
    ia, ib = idx[o[0::2]], idx[o[1::2]]
    w = Mo.M["w"]
    Xa = K.sentez(Mo.s, w[ia], Mo.y[ia]); Xa -= Xa.mean()
    Xb = K.sentez(Mo.s, w[ib], Mo.y[ib]); Xb -= Xb.mean()
    E = Mo.E
    rXaXb = float(np.corrcoef(Xa, Xb)[0, 1])
    print(f"  bölünmüş merdiven: |A|={len(ia)} |B|={len(ib)}  "
          f"korr(Xa,Xb)={rXaXb:+.5f}", flush=True)

    # ---------------- (1) ÇİZGİ çarpanları + basıklık ------------------
    kE, kXa, kXb = kap(E), kap(Xa), kap(Xb)
    uE, uXa, uXb = kurt(E), kurt(Xa), kurt(Xb)
    print(f"  κ_E ={kE:.5f}  κ_Xa={kXa:.5f}  κ_Xb={kXb:.5f}   "
          f"(Gauss {G1:.5f})", flush=True)
    print(f"  kurt: E {uE:+.5f}  Xa {uXa:+.5f}  Xb {uXb:+.5f}   (Gauss 0)",
          flush=True)

    # ---------------- SURROGATE alanları -------------------------------
    rng = np.random.default_rng(1000 + int(tohum))

    def sur(W, amp):
        ps = rng.uniform(0, 2 * np.pi, len(amp))
        f = K.sentez(Mo.s, W, np.abs(amp) * np.exp(1j * ps))
        return f - f.mean()

    Es = sur(w[Mo.msk], Mo.hp[Mo.msk])
    Xas = sur(w[ia], Mo.y[ia])
    Xbs = sur(w[ib], Mo.y[ib])
    print(f"  surrogate: κ_E°={kap(Es):.5f} κ_Xa°={kap(Xas):.5f} "
          f"κ_Xb°={kap(Xbs):.5f}  kurt_E°={kurt(Es):+.5f}  "
          f"korr(Xa°,Xb°)={float(np.corrcoef(Xas,Xbs)[0,1]):+.5f}", flush=True)

    # ---------------- W kümesi (169_k2b ile AYNI) ----------------------
    ban = C163.bant_adaylari(Y, K.IZGARA_T1)
    hedef = [b for b in ban if LO_MIN - 1e-9 <= b["lo"] <= LO_MAX + 1e-9]
    Wall = []
    for b in hedef:
        for r in b["cizgi"]:
            Wall.append(r["w"])
            Wall.append(r["w"] + r["gap"] / 2)
    Wall = np.array(Wall)

    def clip(F):
        return float(np.std(F)) * np.sign(F)

    Ec, Xac, Xbc = clip(E), clip(Xa), clip(Xb)
    Esc, Xasc, Xbsc = clip(Es), clip(Xas), clip(Xbs)
    URUN = {
        "000_tam":   E * Xa * Xb,
        "100_E":     Ec * Xa * Xb,
        "010_Xa":    E * Xac * Xb,
        "111_hepsi": Ec * Xac * Xbc,
        "S000":      Es * Xas * Xbs,      # ° : üçü de surrogate
        "S100":      Esc * Xas * Xbs,
        "S111":      Esc * Xasc * Xbsc,
        "M000":      Es * Xa * Xb,        # * : yalnız E surrogate
        "M100":      Esc * Xa * Xb,
    }
    tm = time.time()
    adlar = list(URUN)
    JL = K.tayf_s(Mo.s, [URUN[a] for a in adlar], Wall)
    J = {a: j / 2.0 for a, j in zip(adlar, JL)}
    print(f"  {len(adlar)} çarpım alanının tayfı {time.time()-tm:.0f}s "
          f"({len(Wall)} frekans)", flush=True)

    def rho_J(var, ref="000_tam"):
        n = float(np.sum(np.real(J[var] * np.conj(J[ref]))))
        d = float(np.sum(np.abs(J[ref]) ** 2))
        return n / d

    rJ_E = rho_J("100_E")
    rJ_Xa = rho_J("010_Xa")
    rJ_3 = rho_J("111_hepsi")
    rJo_E = rho_J("S100", "S000")          # ρ°(E)
    rJo_3 = rho_J("S111", "S000")          # ρ°(hepsi)
    rJm_E = rho_J("M100", "M000")          # ρ*(E)
    print(f"  ρ_J(E)={rJ_E:.5f}  ρ_J(Xa)={rJ_Xa:.5f}  ρ_J(hepsi)={rJ_3:.5f}",
          flush=True)
    print(f"  ρ°(E)={rJo_E:.5f}  ρ°(hepsi)={rJo_3:.5f}   "
          f"(Gauss {G1:.5f} / {G1**3:.5f})", flush=True)
    print(f"  ρ*(E)={rJm_E:.5f}   Rec=(ρ*−ρ_J)/(Gauss−ρ_J) = "
          f"{(rJm_E-rJ_E)/(G1-rJ_E):+.4f}", flush=True)

    # ---------------- diskteki ρ (169_k2b) -----------------------------
    d = json.load(open(S169 / f"K2b_{ad}.json"))
    rho = d["ortalama"]
    tau = [b["tau_eff"] for b in d["bant"]]
    sE = egim(tau, [b["oran"]["100_E"] for b in d["bant"]])
    sXa = egim(tau, [b["oran"]["010_Xa"] for b in d["bant"]])
    sXb = egim(tau, [b["oran"]["001_Xb"] for b in d["bant"]])
    bE = rho["100_E"] / kE
    payC = (np.log(kE) - np.log(G1)) / (np.log(rho["100_E"]) - np.log(G1))
    print(f"  ρ(E)={rho['100_E']:.5f} (disk)   β_E=ρ/κ={bE:.5f}   "
          f"PAY_ÇİZGİ={payC:+.4f}", flush=True)
    print(f"  τ eğimleri: s_E={sE:+.4f}  s_Xa={sXa:+.4f}  s_Xb={sXb:+.4f}",
          flush=True)

    # ---------------- R_η (176/K1_*.json, varsa) -----------------------
    Ret = None
    p = S176 / f"K1_{ad}.json"
    if p.exists():
        k1 = json.load(open(p))["eta"]
        Ret = dict(R=k1["R"], Var=k1["Var"], P=k1["P"],
                   kov_art_x=k1["kov_art_x"], ham_art_x=k1["ham_art_x"],
                   Var_cizgi=k1["Var_cizgi"], Var_artik=k1["Var_artik"],
                   ic_bant=k1["ic_bant"], capraz=k1["capraz"])
        print(f"  R_η={k1['R']:.5f}   1−R_η={1-k1['R']:+.5f}   "
              f"Kov(artık,η)/Var(η)={k1['kov_art_x']/k1['Var']:+.5f}   "
              f"Var(çizgi)/P={k1['Var_cizgi']/k1['P']:.5f}   "
              f"çapraz/iç={k1['capraz']/k1['ic_bant']:+.5f}", flush=True)
    else:
        print(f"  (K1_{ad}.json yok — R_η bu tohumda ölçülmedi)", flush=True)

    rec = dict(ad=ad, tohum=int(tohum), korr_XaXb=rXaXb,
               kappa=dict(E=kE, Xa=kXa, Xb=kXb),
               kappa_sur=dict(E=kap(Es), Xa=kap(Xas), Xb=kap(Xbs)),
               kurt=dict(E=uE, Xa=uXa, Xb=uXb, Es=kurt(Es)),
               korr_XaXb_sur=float(np.corrcoef(Xas, Xbs)[0, 1]),
               rho_disk=rho, rho_J=dict(E=rJ_E, Xa=rJ_Xa, hepsi=rJ_3),
               rho_sur=dict(E=rJo_E, hepsi=rJo_3), rho_Esur=rJm_E,
               Rec=(rJm_E - rJ_E) / (G1 - rJ_E),
               beta_E=bE, pay_cizgi=payC,
               egim=dict(E=sE, Xa=sXa, Xb=sXb), tau=tau,
               bant_rhoE=[b["oran"]["100_E"] for b in d["bant"]],
               R_eta=Ret, sur_tohum=1000 + int(tohum),
               nW=int(len(Wall)), sure_s=time.time() - t0)
    S182.mkdir(parents=True, exist_ok=True)
    q = S182 / f"ANATOMI_{ad}.json"
    q.write_text(json.dumps(rec, indent=1, ensure_ascii=False, default=float))
    print(f"  -> {q}  ({(time.time()-t0)/60:.1f} dk)", flush=True)
    return rec


if __name__ == "__main__":
    gz = sys.argv[1:] or [f"VF{i}" for i in range(1, 9)]
    for g in gz:
        kos(g, int(g.replace("VF", "")))
