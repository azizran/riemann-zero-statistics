"""
166 — HİPOTEZ YARIŞI (T2) + KAPANIŞ (T3) + ŞİŞME (T4)
=====================================================
Girdi: 166_T1 (kalib tabloları + jackknife), 166_rho (ρ bantları),
166_bacak (sönüm/sızıntı çarpanları). Hiçbir sayı elle yazılmaz.

HEDEF BÜYÜKLÜKLER
    KALİB_u2 = A²u2_ölç / A²u2_öng   — 165'in KULLANDIĞI kalibrasyon
                                       (yalnız İKİNCİ momentler)
    KALİB_s2 = A²s2_ölç / A²s2_öng   — üçüncü momentin İSTEDİĞİ kalibrasyon
Yasalar KALİB_u2'ye (ikinci moment) sabitlenir, KALİB_s2 üzerinde
SINANIR: üçüncü momentten hiçbir girdi alınmaz.

ÖLÇEK SINIFLARI
    S0  — SIFIR ölçek: yasa doğrudan kalib'dir.
    S1  — TEK ölçek: üç gaz ve bütün bantlar için TEK global c.
    S1g — gaz-başına ölçek (ZAYIF; yalnız tanı için, hükme girmez).

HAKEMLER
    (i)  bant-şekli: log-artık rms (%), maks sapma, χ²/dof (jackknife)
    (ii) ÇAPRAZ-GAZ: her gazın tek başına isteyeceği ölçek c_g; yayılım
         maks/min — S1 iddiası ancak bu 1'e yakınsa geçerlidir
    (iii) artık yapısı: log-artığın τ ile eğimi (± jackknife)
"""
import importlib
import json
import sys
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/166")
GAZ = ("son", "Hkeskin", "HA4")
KAZANANLAR = ("K1o", "K1b", "K1d")
LO_MIN = 0.52          # görevin τ kısıtı: hüküm lo ≥ 0.52 bantlarından


def yukle():
    D = {}
    for g in GAZ:
        T1 = json.load(open(SCR / f"T1_{g}_olculen_tc0.95.json"))
        RH = json.load(open(SCR / f"RHO_{g}.json"))
        BC = json.load(open(SCR / f"BACAK_{g}.json"))
        r52 = {round(b["lo"], 4): b for b in RH["kosum"]["t0.52_c720_b165"]}
        r40 = {round(b["lo"], 4): b for b in RH["kosum"]["t0.4_c4000_b165"]}
        bc = {round(b["lo"], 4): b for b in BC["bant"]}
        rec = []
        for b in T1["bant"]:
            if not b.get("olculdu"):
                continue
            lo = round(b["lo"], 4)
            d = dict(gaz=g, lo=lo, hi=b["hi"], tau_eff=b["tau_eff"],
                     Ms2=b["Ms2"], Ps2=b["Ps2"], Mu2=b["Mu2"], Pu2=b["Pu2"],
                     K_s2=b["KALIB_s2"], K_u2=b["KALIB_u2"],
                     sK_s2=b.get("sKALIB_s2", np.nan),
                     sK_u2=b.get("sKALIB_u2", np.nan),
                     sMs2=b.get("sMs2", np.nan), sPs2=b.get("sPs2", np.nan))
            d["rho52"] = r52[lo]["rho"] if lo in r52 else np.nan
            d["srho52"] = r52[lo]["srho"] if lo in r52 else np.nan
            d["rho40"] = r40[lo]["rho"] if lo in r40 else np.nan
            d.update({k: bc[lo][k] for k in
                      ("lamE", "lamX", "absE", "absX", "W_amp", "W_X",
                       "W_pos", "R_h", "R_y", "tau_gp")})
            rec.append(d)
        # Λ (merdiven-ortalamalı sızıntı) DOĞRU AĞIRLIKLA yeniden:
        # JSON'daki alan b²-ağırlıklıdır ve τ ≤ taban'da h'=0 olduğundan
        # patlıyor. Kimlik: Var(E_mod)/⟨Eη⟩ = Σ|h'|²Reλ_E / Σ|h'|² .
        z = np.load(SCR / f"BACAK_{g}.npz")
        Lam = dict(BC["Lam"])
        for et, am, lm in (("E", z["hp"], z["lamE"]), ("X", z["y"], z["lamX"])):
            wgt = np.abs(am) ** 2 * z["msk"]
            fin = np.isfinite(lm.real) & (wgt > 0)
            Lam[f"Lam{et}"] = float(np.sum(lm.real[fin] * wgt[fin])
                                    / np.sum(wgt[fin]))
        D[g] = dict(T1=T1, bant=rec, Lam=Lam, sigC=BC["sigC"],
                    sigds=BC["sigds"], sigX=BC["sigX"])
    return D


def saglikli(D):
    """Hükme giren bantlar: lo ≥ 0.52 ve 160'ın sağlık kuralı (τ_eff < 0.85,
    ölçülen A²s2 > 0). 165'in 20 bandını birebir verir."""
    out = []
    for g in GAZ:
        for b in D[g]["bant"]:
            if b["lo"] < LO_MIN - 1e-9:
                continue
            if b["tau_eff"] > 0.85 or b["Ms2"] <= 0:
                continue                     # HA4'ün ‡ bandı burada eleniyor
            out.append(b)
    return out


# ----------------------------------------------------------------------
def yasalar(D):
    """(ad, sinif, fn) — fn(bant, gaz_kaydi) → yasanın ÇIPLAK değeri."""
    Y = []
    A = Y.append
    A(("K0  sabit (bant bağımsız)", "S1", lambda b, g: 1.0))
    A(("K0  g = g_E g_X²", "S0", lambda b, g: g["Lam"]["gcal"]))
    # --- K2: tarak-sağkalımı ρ ---------------------------------------
    A(("K2a ρ(τ) [0.52-taban]", "S0", lambda b, g: b["rho52"]))
    A(("K2b c·ρ(τ) [0.52-taban]", "S1", lambda b, g: b["rho52"]))
    A(("K2c ρ(τ) [0.40-taban]", "S0", lambda b, g: b["rho40"]))
    A(("K2d c·ρ(τ) [0.40-taban]", "S1", lambda b, g: b["rho40"]))
    A(("K2e c·√ρ [0.52-taban]", "S1",
       lambda b, g: np.sqrt(max(b["rho52"], 1e-9))))
    A(("K2f c·√ρ [0.40-taban]", "S1",
       lambda b, g: np.sqrt(max(b["rho40"], 1e-9))))
    # --- K1: bacak sönümleri ------------------------------------------
    A(("K1a W_pos(τ)=⟨e^{2πiτĈ}⟩", "S0", lambda b, g: b["W_pos"]))
    A(("K1b c·W_pos(τ)", "S1", lambda b, g: b["W_pos"]))
    A(("K1c c·W_X(τ)=⟨e^{−2πiτX̃}⟩", "S1", lambda b, g: b["W_X"]))
    A(("K1d c·W_amp(τ)=⟨cos πτ ds⟩", "S1", lambda b, g: b["W_amp"]))
    A(("K1e c·W_pos·W_amp", "S1", lambda b, g: b["W_pos"] * b["W_amp"]))
    A(("K1f c·W_pos²", "S1", lambda b, g: b["W_pos"] ** 2))
    A(("K1g g·W_pos (sıfır ölçek)", "S0",
       lambda b, g: g["Lam"]["gcal"] * b["W_pos"]))
    A(("K1o c·W_X·W_amp (X-bacağın tamamı)", "S1",
       lambda b, g: b["W_X"] * b["W_amp"]))
    A(("K1p c·√(W_X·W_amp)", "S1",
       lambda b, g: np.sqrt(b["W_X"] * b["W_amp"])))
    A(("K1q c·W_X^{3/2}", "S1", lambda b, g: b["W_X"] ** 1.5))
    A(("K1r c·exp(−2π²σ_Ĉ²τ²)  [Gauss]", "S1",
       lambda b, g: np.exp(-2 * np.pi ** 2 * g["sigC"] ** 2
                           * b["tau_eff"] ** 2)))
    A(("K1s c·exp(−2π²σ_X̃²τ²)  [Gauss]", "S1",
       lambda b, g: np.exp(-2 * np.pi ** 2 * g["sigX"] ** 2
                           * b["tau_eff"] ** 2)))
    A(("K1t c·exp(−½π²σ_ds²τ²) [Gauss]", "S1",
       lambda b, g: np.exp(-0.5 * np.pi ** 2 * g["sigds"] ** 2
                           * b["tau_eff"] ** 2)))
    A(("K1u c·W_pos·W_X", "S1", lambda b, g: b["W_pos"] * b["W_X"]))
    # --- Gram sızıntısı (bacakların ÖLÇÜLEN sönümü) -------------------
    A(("K1h 1/λ_X(ω_Q)", "S0", lambda b, g: 1.0 / b["lamX"]))
    A(("K1i c/λ_X(ω_Q)", "S1", lambda b, g: 1.0 / b["lamX"]))
    A(("K1j 1/λ_E(ω_Q)", "S0", lambda b, g: 1.0 / b["lamE"]))
    A(("K1k 1/(λ_E λ_X²)(ω_Q)", "S0",
       lambda b, g: 1.0 / (b["lamE"] * b["lamX"] ** 2)))
    A(("K1l c/(λ_E λ_X²)(ω_Q)", "S1",
       lambda b, g: 1.0 / (b["lamE"] * b["lamX"] ** 2)))
    A(("K1m 1/(λ_X·Λ_E Λ_X)  [Λ=merdiven ort.]", "S0",
       lambda b, g: 1.0 / (b["lamX"] * g["Lam"]["LamE"] * g["Lam"]["LamX"])))
    A(("K1n c/(λ_X·λ_E)", "S1", lambda b, g: 1.0 / (b["lamX"] * b["lamE"])))
    # --- K3: çarpımlar -------------------------------------------------
    A(("K3a c·ρ52·W_pos", "S1", lambda b, g: b["rho52"] * b["W_pos"]))
    A(("K3b c·ρ40·W_pos", "S1", lambda b, g: b["rho40"] * b["W_pos"]))
    A(("K3c c·ρ52/λ_X", "S1", lambda b, g: b["rho52"] / b["lamX"]))
    A(("K3d c·ρ40/λ_X", "S1", lambda b, g: b["rho40"] / b["lamX"]))
    A(("K3e c·√ρ40·W_pos", "S1",
       lambda b, g: np.sqrt(max(b["rho40"], 1e-9)) * b["W_pos"]))
    A(("K3f c·W_pos/λ_X", "S1", lambda b, g: b["W_pos"] / b["lamX"]))
    return Y


def degerlendir(ad, sinif, fn, B, D, hedef="K_u2"):
    v = np.array([fn(b, D[b["gaz"]]) for b in B])
    k = np.array([b[hedef] for b in B])
    s = np.array([b["s" + hedef] for b in B])
    ok = np.isfinite(v) & np.isfinite(k) & (v > 0) & (k > 0)
    if ok.sum() < 5:
        return None
    lr = np.log(k[ok]) - np.log(v[ok])          # log-artık
    c = float(np.exp(lr.mean())) if sinif != "S0" else 1.0
    pred = c * v[ok]
    res = np.log(k[ok] / pred)
    gz = np.array([b["gaz"] for b in B])[ok]
    cg = {g: float(np.exp(lr[gz == g].mean())) for g in GAZ if (gz == g).any()}
    yay = max(cg.values()) / min(cg.values())
    tau = np.array([b["tau_eff"] for b in B])[ok]
    egim = float(np.polyfit(tau, res, 1)[0])
    chi = float(np.sum(((k[ok] - pred) / s[ok]) ** 2))
    dof = ok.sum() - (0 if sinif == "S0" else 1)
    return dict(ad=ad, sinif=sinif, n=int(ok.sum()), c=c,
                rms=float(np.sqrt(np.mean(res ** 2))),
                maks=float(np.max(np.abs(res))), chi2=chi, dof=int(dof),
                cg=cg, yayilim=float(yay), egim=egim, v=v, ok=ok)


def t3_say(v, B, c=1.0, hedef=0.25):
    """Yasayı ÇIPLAK öngörüye çarp: (c·yasa·Ps2)/Ms2 ∈ [1−h, 1+h]?"""
    o, sat = 0, []
    for i, b in enumerate(B):
        r = c * v[i] * b["Ps2"] / b["Ms2"] if b["Ms2"] else np.nan
        iyi = np.isfinite(r) and abs(r - 1.0) <= hedef
        o += bool(iyi)
        sat.append((b["gaz"], b["tau_eff"], float(r), bool(iyi)))
    return o, sat


def main():
    D = yukle()
    B = saglikli(D)
    print(f"=== 166 YARIŞ === hüküm bantları: {len(B)} "
          f"({', '.join(f'{g}:{sum(1 for b in B if b[chr(103)+chr(97)+chr(122)]==g)}' for g in GAZ)})")
    for g in GAZ:
        L = D[g]["Lam"]
        print(f"  {g:9s} σ_Ĉ={D[g]['sigC']:.4f} σ_ds={D[g]['sigds']:.4f} "
              f"σ_X̃={D[g]['sigX']:.4f} g_E={L['gE']:.4f} g_X={L['gX']:.4f} "
              f"g={L['gcal']:.4f} Λ_E={L['LamE']:.4f} Λ_X={L['LamX']:.4f}")

    print("\n--- T1 HAM TABLO (hüküm bantları) ---")
    print("gaz       τ_eff   KALİB_s2 ±jk        KALİB_u2 ±jk     "
          "ρ52     ρ40    W_pos   λ_E    λ_X")
    for b in B:
        print(f"{b['gaz']:9s} {b['tau_eff']:.4f} {b['K_s2']:+.4f}±"
              f"{b['sK_s2']:.4f}  {b['K_u2']:.4f}±{b['sK_u2']:.4f}  "
              f"{b['rho52']:+.4f} {b['rho40']:+.4f} {b['W_pos']:.4f} "
              f"{b['lamE']:.4f} {b['lamX']:.4f}")

    Bg = [b for g in GAZ for b in D[g]["bant"]
          if b["Ms2"] > 0 and b["tau_eff"] < 0.85]      # 27 bant (lo ≥ 0.44)
    Bp = [b for b in B if 0.53 < b["tau_eff"] < 0.71]   # 165'in sağlıklı
    for hedef, BB, et in (("K_u2", B, "20 hüküm bandı"),
                          ("K_u2", Bg, "27 bant — GENİŞ şekil tanısı"),
                          ("K_s2", B, "20 hüküm bandı"),
                          ("K_s2", Bp, "15 bant — τ_eff∈(0.53,0.71)")):
        print(f"\n=== YARIŞ — hedef {hedef} [{et}] ===")
        print(f"{'yasa':38s} sın  n   c       rms%   maks%   χ²/dof   "
              f"çapraz-gaz(maks/min)  eğim")
        R = []
        for ad, sn, fn in yasalar(D):
            r = degerlendir(ad, sn, fn, BB, D, hedef)
            if r:
                R.append(r)
        R.sort(key=lambda r: r["rms"])
        for r in R:
            print(f"{r['ad']:38s} {r['sinif']:3s} {r['n']:2d} {r['c']:7.4f} "
                  f"{100*r['rms']:6.2f} {100*r['maks']:6.2f} "
                  f"{r['chi2']/max(r['dof'],1):9.1f}  {r['yayilim']:6.3f}  "
                  f"[{' '.join(f'{k}:{v:.3f}' for k, v in r['cg'].items())}]"
                  f"  {r['egim']:+.3f}")


    # ---------------- T3 -------------------------------------------
    print("\n=== T3 — KAPANIŞ: yasa ile ±%25 sayımı (u2-sabitlemesiz) ===")
    print("  (ölçek c YALNIZ KALİB_u2'ye, yani İKİNCİ momentlere oturtuldu)")
    ref = []
    for ad, sn, fn in yasalar(D):
        r = degerlendir(ad, sn, fn, B, D, "K_u2")
        if not r:
            continue
        v = np.array([fn(b, D[b["gaz"]]) for b in B])
        o, sat = t3_say(v, B, r["c"])
        gsay = {}
        for g, t, x, iy in sat:
            gsay.setdefault(g, [0, 0])
            gsay[g][1] += 1
            gsay[g][0] += iy
        print(f"  {ad:38s} {sn:3s} c={r['c']:7.4f}  **{o}/{len(B)}**   "
              + " ".join(f"{g}:{a}/{b}" for g, (a, b) in gsay.items()))
        if ad.split()[0] in KAZANANLAR:
            ref.append((ad, r, v, sat))
    # 165'in kendi u2 kalibrasyonu (referans)
    v0 = np.array([b["K_u2"] for b in B])
    o0, sat0 = t3_say(v0, B, 1.0)
    print(f"  {'[165] u2-normalizasyonu (referans)':38s} —   c=1.0000  "
          f"**{o0}/{len(B)}**")

    for ad, r, v, sat in ref:
        print(f"\n  --- {ad} ile bant bant (c={r['c']:.4f}) ---")
        print("  gaz       τ_eff   kalib_yasa  öngörü_kal   ölçüm     oran  ±%25")
        for i, b in enumerate(B):
            kl = r["c"] * v[i]
            print(f"  {b['gaz']:9s} {b['tau_eff']:.4f}  {kl:8.4f}    "
                  f"{kl*b['Ps2']:+.6f}  {b['Ms2']:+.6f}  {sat[i][2]:6.3f}  "
                  f"{'✓' if sat[i][3] else ''}")

    # ---------------- T4 -------------------------------------------
    print("\n=== T4 — GÖSTERİM ŞİŞMESİ ===")
    print("  gaz       Var(E_mod)/Var(η)  Var(X_mod)/Var(X̃)   1/g_E   1/g_X"
          "   Λ_E    Λ_X   Σ|h'|²/2Var(η)  1/W̄_pos  1/W̄_pos²")
    for g in GAZ:
        L = D[g]["Lam"]
        print(f"  {g:9s} {L['varE_mod']/L['varE_olc']:14.3f} "
              f"{L['varX_mod']/L['varX_olc']:17.3f} "
              f"{1/L['gE']:9.3f} {1/L['gX']:6.3f} {L['LamE']:6.3f} "
              f"{L['LamX']:6.3f} {L['sum_hp2']/L['varE_olc']:12.3f} "
              f"{1/L['W_pos_bar']:9.3f} {1/L['W_pos_bar']**2:9.3f}")


if __name__ == "__main__":
    main()
