"""
169 — K2 (MEKANİZMA): sarılma testi + işaret-aktarımı
=====================================================
Ölçüm/öngörü zinciri KOPYALANMAZ: `165_cekirdek.Model165`, `163_cekirdek.
olc_cizgi/bant_adaylari`, `166_T1.bant_agg` aynen import edilir.

K2(a) SARILMA. 165 §7 / 168 §A1.5'in geometrisi:
    s_n = s̄ + ḡ(n−n̄) + ḡ·Ĉ_n ,  ḡ = 2π/L ,  Ĉ = kübik-trendsiz Σ X̃0
    bir ν frekansının taşıdığı SİTE-SEĞİRME fazı:  ϑ_n(ν) = ν ḡ Ĉ_n = 2π τ_ν Ĉ_n
Sarılma ölçütleri (N = 3e5 örnek):
    R(τ)  = |⟨e^{2πiτĈ}⟩|          (= 166_bacak.karakteristik = W_pos(τ))
    V(τ)  = Kuiper istatistiği (mod 2π; düzgünlük H0)
    √N·V  düzgünlükte ≈ 1.62 (medyan); √N·V > 2.00 → %1'de RED
KONTRAST: düşük-τ (q=2, τ=0.058) DÜZGÜN OLMAMALI; Ç4'ün dört-frekans
toplamı (τ_Q+τ₁+τ₂+τ₃ ≈ 2–3) TAM düzgün olmalı.

K2(b) İŞARET-AKTARIMI. Tanımlar AÇIK:
  * model alanları (165 §1): E = Σ_q Re[h'_q e^{iω_q s}], X = Σ_q Re[y_q e^{iω_q s}]
  * VARYANS-EŞLEŞMİŞ işaret kırpması:  clip(F) ≡ σ_F · sgn(F)
    (Gauss alanı için ⟨clip(F)·F⟩/⟨F²⟩ = √(2/π) = 0.79788 — "bacak başına
     Gauss aktarımı"; iki bacak ⇒ 2/π, DÖRT bacak ⇒ (2/π)² = 4/π².)
  * öngörü nesnesi 165/167 ile birebir: J₂(W) = ⟨G₂ e^{−iW s}⟩,
    u₂ = Re[h̄_Q J₂]/(2⟨ρ⟩), bant birleştirmesi 166_T1.bant_agg
  * varyantlar:  G₂ = E·X·X (tam) → clip'lenmiş bacak sayısı 1,2,3
  * aktarım oranı  ρ_var ≡ Pu2(varyant)/Pu2(tam)   (bant bant)
  * AYRICA arcsine ailesi: her bacak için r = korr(model, ölçülen) ve
    ⟨sgn·sgn⟩; Gauss çiftinde ⟨sgn x sgn y⟩ = (2/π) arcsin r.

Kullanım: 169_k2.py <gaz> [taban] [tau_c]
Çıktı:    scratchpad/169/K2_<gaz>.json
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
T166 = importlib.import_module("166_T1")
B166 = importlib.import_module("166_bacak")
ORT = importlib.import_module("167_ortak")

TWO_PI = 2 * np.pi
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/169")
LO_MIN, LO_MAX = 0.52, 0.68
KMAX = 3
K.PENCERE.update(ORT.pencere_dict())


def kuiper(th):
    """Kuiper V = D+ + D− ; mod 2π açılar için düzgünlük testi."""
    u = np.sort((np.asarray(th) % TWO_PI) / TWO_PI)
    n = len(u)
    i = np.arange(1, n + 1)
    return float(np.max(i / n - u) + np.max(u - (i - 1) / n))


def main(veri="Hkeskin", taban=0.40, tau_c=0.95):
    t0 = time.time()
    SCR.mkdir(parents=True, exist_ok=True)
    Y = K.gaz(veri, taban, 4000)
    Mo = K.Model165(veri, taban, 4000, 0.95, "olculen", Y=Y, tau_c=tau_c)
    Mo.sec("olculen", tau_c).alanlar(kmax=KMAX)
    A = Mo.artik
    print(f"=== 169-K2 {veri} taban={taban} τ_c={tau_c} ===", flush=True)
    print(f"  N={len(Y.m0)} L={Y.L:.5f} T={Mo.T:.1f} σ_Ĉ={Mo.sigC:.5f} "
          f"g_E={A['gE']:.4f} g_X={A['gX']:.4f} çizgi={A['nline']}", flush=True)

    # ---------------- K2(a) SARILMA ------------------------------------
    C = np.cumsum(Y.Xtil0)
    n = np.arange(len(C), dtype=float)
    n = (n - n.mean()) / (n[-1] if len(n) > 1 else 1.0)
    V = np.vstack([np.ones_like(n), n, n * n, n ** 3]).T
    cc, *_ = np.linalg.lstsq(V, C, rcond=None)
    Chat = C - V @ cc
    NN = len(Chat)
    sq = np.sqrt(NN)
    print(f"\n  Ĉ: σ={np.std(Chat):.5f}  N={NN}  (√N·V düzgünlükte ≈1.62, "
          f"%1 red eşiği ≈2.00)", flush=True)

    # model alanlarındaki genlik-ağırlıklı ortalama τ (Ç4 bacakları)
    msk = Mo.msk
    tl = Mo.M["tau"][msk]
    wh = np.abs(Mo.hp[msk])
    wy = np.abs(Mo.y[msk])
    tE = float((tl * wh).sum() / wh.sum())
    tX = float((tl * wy).sum() / wy.sum())
    print(f"  model bacaklarının genlik-ağırlıklı τ'su: ⟨τ⟩_E={tE:.4f}  "
          f"⟨τ⟩_X={tX:.4f}", flush=True)

    ban = C163.bant_adaylari(Y, K.IZGARA_T1)
    hedef = [b for b in ban if LO_MIN - 1e-9 <= b["lo"] <= LO_MAX + 1e-9]
    sar = []
    print("\n  --- K2(a) SARILMA: ϑ = 2πτĈ mod 2π ---")
    print("   etiket                          τ        R(τ)=|⟨e^{iϑ}⟩|   "
          "√N·V(Kuiper)   düzgün mü?")

    def satir(et, tau):
        th = TWO_PI * tau * Chat
        R = float(abs(np.mean(np.exp(1j * th))))
        Vk = kuiper(th) * sq
        duz = "DÜZGÜN" if Vk < 2.00 else "DEĞİL"
        print(f"   {et:30s} {tau:7.4f}   {R:10.5f}     {Vk:8.2f}      {duz}")
        sar.append(dict(etiket=et, tau=float(tau), R=R, kuiperV=Vk,
                        duzgun=(Vk < 2.00)))

    satir("KONTROL q=2 (tek çizgi)", np.log(2) / Y.L)
    satir("KONTROL q=7 (tek çizgi)", np.log(7) / Y.L)
    satir("KONTROL q=1747 (tek çizgi)", np.log(1747) / Y.L)
    satir("merdiven ucu τ=0.95 (tek)", 0.95)
    for b in hedef:
        satir(f"bant τ_Q={b['tau']:.3f} (TEK taşıyıcı)", b["tau"])
    for b in hedef:
        satir(f"bant τ_Q={b['tau']:.3f} Ç4 (4 frekans)",
              b["tau"] + tE + 2 * tX)
    satir("Ç4 alt sınır (τ_Q+3·0.5)", hedef[0]["tau"] + 1.5)

    # ---------------- K2(b) İŞARET-AKTARIMI ----------------------------
    print("\n  --- K2(b) İŞARET-AKTARIMI ---")
    E, X = Mo.E, Mo.X
    sE, sX = float(np.std(E)), float(np.std(X))
    Ec, Xc = sE * np.sign(E), sX * np.sign(X)
    # bacak başına Gauss aktarımı (doğrudan ölçüm)
    aE = float(np.dot(Ec, E) / np.dot(E, E))
    aX = float(np.dot(Xc, X) / np.dot(X, X))
    print(f"   bacak başına ⟨clip(F)F⟩/⟨F²⟩:  E: {aE:.5f}   X: {aX:.5f}   "
          f"(Gauss: √(2/π) = {np.sqrt(2/np.pi):.5f})")
    # arcsine ailesi: model ↔ ÖLÇÜLEN alan çiftleri
    e1 = Y.e1 - Y.e1.mean()
    x1 = Y.Xtil0 - Y.Xtil0.mean()
    ars = {}
    for et, (u, v) in (("E", (E, e1)), ("X", (X, x1))):
        r = float(np.corrcoef(u, v)[0, 1])
        ss = float(np.mean(np.sign(u) * np.sign(v)))
        pred = (2 / np.pi) * np.arcsin(r)
        ars[et] = dict(r=r, sgnsgn=ss, arcsin=float(pred), oran=ss / r)
        print(f"   arcsine bacak {et}: r={r:.5f}  ⟨sgn·sgn⟩={ss:.5f}  "
              f"(2/π)arcsin r={pred:.5f}  fark={100*(ss/pred-1):+.2f}%  "
              f"⟨sgn sgn⟩/r={ss/r:.5f}")

    VAR = {                        # ad -> (E alanı, X₁, X₂), kırpılan bacak
        "tam":       (E, X, X, 0),
        "E":         (Ec, X, X, 1),
        "X1":        (E, Xc, X, 1),
        "X1X2":      (E, Xc, Xc, 2),
        "E+X1":      (Ec, Xc, X, 2),
        "E+X1X2":    (Ec, Xc, Xc, 3),
    }
    Wall, bas = [], []
    for b in hedef:
        bas.append(len(Wall))
        for r in b["cizgi"]:
            Wall.append(r["w"])
            Wall.append(r["w"] + r["gap"] / 2)
    Wall = np.array(Wall)
    print(f"   {len(hedef)} bant, {len(Wall)//2} çizgi ({len(Wall)} frekans)",
          flush=True)
    tm = time.time()
    olc = [C163.olc_cizgi(Y, float(w), kmax=KMAX) for w in Wall]
    print(f"   ölçüm {time.time()-tm:.0f}s", flush=True)
    tm = time.time()
    JJ = {}
    for ad, (a, b1, b2, nk) in VAR.items():
        JJ[ad] = [j / 2.0 for j in K.tayf_s(Mo.s, [a * b1 * b2], Wall)][0]
    print(f"   {len(VAR)} varyantın öngörüsü {time.time()-tm:.0f}s", flush=True)

    print("\n   bant τ_eff |  " + "  ".join(f"ρ({a})" for a in VAR if a != "tam")
          + "   [ρ = Pu2(varyant)/Pu2(tam)]")
    tab = []
    for bi, b in enumerate(hedef):
        L = []
        for li, r in enumerate(b["cizgi"]):
            d = dict(tau=r["tau"], gp=r["gp"], grup=r["grup"])
            k0 = bas[bi] + 2 * li
            for et, off in (("on", 0), ("off", 1)):
                o = olc[k0 + off]
                d[f"A_{et}"] = o["A"]
                d[f"pow_{et}"] = o["pow"]
                for kk in range(KMAX + 1):
                    d[f"s{kk}_{et}"] = o[f"s{kk}"]
                    d[f"u{kk}_{et}"] = o[f"u{kk}"]
                for ad in VAR:
                    sp, up = K.s_den_J(o["h"], JJ[ad][k0 + off], o["rho_ort"])
                    d[f"s{ad}_{et}"] = sp
                    d[f"u{ad}_{et}"] = up
                # bant_agg'in beklediği alan adları (yalnız k=2 kullanılacak)
                sp, up = K.s_den_J(o["h"], JJ["tam"][k0 + off], o["rho_ort"])
                for kk in range(KMAX + 1):
                    d[f"s{kk}p_{et}"] = sp
                    d[f"u{kk}p_{et}"] = up
            L.append(d)
        pN = np.array([r["pow_on"] for r in L])
        pO = np.array([r["pow_off"] for r in L])
        Ao = np.array([r["A_on"] for r in L])
        Af = np.array([r["A_off"] for r in L])
        dd = float((pN - pO).sum())
        tef = float((np.array([r["tau"] for r in L]) * (pN - pO)).sum() / dd)

        def AG(key):
            vN = np.array([r[f"{key}_on"] for r in L])
            vO = np.array([r[f"{key}_off"] for r in L])
            return float((pN * Ao ** 2 * vN - pO * Af ** 2 * vO).sum() / dd)

        P = {ad: AG(f"u{ad}") for ad in VAR}
        Pu2 = float((pN * Ao ** 2 * np.array([r["u2_on"] for r in L])
                     - pO * Af ** 2 * np.array([r["u2_off"] for r in L])).sum()
                    / dd)
        row = dict(tau_eff=tef, lo=b["lo"], Mu2=Pu2, Pu2_tam=P["tam"],
                   KALIB=Pu2 / P["tam"],
                   oran={ad: P[ad] / P["tam"] for ad in VAR})
        tab.append(row)
        print(f"   {tef:.4f}    | " +
              "  ".join(f"{P[a]/P['tam']:.4f}" for a in VAR if a != "tam") +
              f"   (KALİB_u2={Pu2/P['tam']:.4f})")

    print("\n   BEKLENEN (Gauss kırpma, bacak başına √(2/π)):")
    g1 = np.sqrt(2 / np.pi)
    print(f"     1 bacak: {g1:.4f}   2 bacak: {g1**2:.4f} (=2/π)   "
          f"3 bacak: {g1**3:.4f}   4 bacak: {g1**4:.4f} (=4/π²)")

    out = dict(veri=veri, taban=taban, tau_c=tau_c, N=NN,
               sigChat=float(np.std(Chat)), tau_E=tE, tau_X=tX,
               sarilma=sar, arcsine=ars, aE=aE, aX=aX,
               bant=tab, sure_s=time.time() - t0)
    p = SCR / f"K2_{veri}.json"
    p.write_text(json.dumps(out, indent=1))
    print(f"\n-> {p}  ({(time.time()-t0)/60:.1f} dk)", flush=True)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "Hkeskin",
         float(sys.argv[2]) if len(sys.argv) > 2 else 0.40,
         float(sys.argv[3]) if len(sys.argv) > 3 else 0.95)
