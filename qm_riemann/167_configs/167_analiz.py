"""
167 — ABLASYON HARİTASI (T2) + HASSAS TABAN (T1) + KAPANIŞ (T4)
===============================================================
Girdi: scratchpad/167/C_<gaz>.json (167_olcum'un çıktıları).
Hiçbir sayı elle yazılmaz.

ÖLÇÜLEN BÜYÜKLÜK
    KALİB_u2(bant) = A²u2_ölç / A²u2_öng          (166'nın hedefi)
    c_M(bant)      = KALİB_u2 / M(τ_bant)          M ∈ {1, W_amp, W_X,
                                                        W_amp·W_X, W_pos}
    c(gaz)         = hüküm bantlarında log-ortalama; hata = bant-içi
                     jackknife'ların log-ortalamaya taşınmış birleşimi
                     (+ bantlar arası saçılım AYRI raporlanır)

HÜKÜM BANTLARI: lo ≥ 0.52, τ_eff < 0.85, A²s2_ölç > 0 (160'ın sağlık
kuralı) VE sadakat R_bant ≥ 0.98 (164 ölçütü) VE bant SNR ≥ 4.
"""
import importlib
import json
import sys
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
sys.path.insert(0, str(QM / "167_configs"))
ORT = importlib.import_module("167_ortak")

SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/167")
LO_MIN, SNR_MIN, R_MIN = 0.52, 3.0, 0.98   # SNR eşiği 164 §3b'nin ölçütü

UYELER = (("K0  sabit", lambda b: 1.0),
          ("K1d W_amp", lambda b: b["W_amp"]),
          ("K1c W_X", lambda b: b["W_X"]),
          ("K1o W_amp·W_X", lambda b: b["W_amp"] * b["W_X"]),
          ("K1p √(W_amp W_X)", lambda b: np.sqrt(b["W_amp"] * b["W_X"])),
          ("K1b W_pos", lambda b: b["W_pos"]))


def yukle(adlar=None):
    D = {}
    for p in sorted(SCR.glob("C_*.json")):
        ad = p.stem[2:]
        if adlar and ad not in adlar:
            continue
        D[ad] = json.load(open(p))
    return D


def hukum(d, gevsek=False, lo_max=None):
    """Hüküm bantları. lo_max verilirse ORTAK PENCERE (ablasyon eksenlerinde
    gazlar arası adil karşılaştırma; K070/E060 ve T/4 dilimleri yüksek
    bantları kaybeder, c(τ) ise τ ile yükselir — aynı bant kümesi şart)."""
    out = []
    for b in d["bant"]:
        if not b.get("olculdu"):
            continue
        if b["lo"] < LO_MIN - 1e-9 or b["tau_eff"] > 0.85 or b["Ms2"] <= 0:
            continue
        if lo_max is not None and b["lo"] > lo_max + 1e-9:
            continue
        if not gevsek and (b["R_bant"] < R_MIN or b["SNR"] < SNR_MIN):
            continue
        out.append(b)
    return out


def c_gaz(B, uye=lambda b: b["W_amp"] * b["W_X"]):
    """log-ortalama c ve iki hata: (i) jackknife'ların birleşimi,
    (ii) bantlar arası saçılım (sistematik)."""
    if not B:
        return dict(c=np.nan, n=0)
    v = np.array([b["KALIB_u2"] / uye(b) for b in B])
    s = np.array([b.get("sKALIB_u2", np.nan) / b["KALIB_u2"] for b in B])
    lg = np.log(v)
    c = float(np.exp(lg.mean()))
    n = len(B)
    jk = float(c * np.sqrt(np.nansum(s ** 2)) / n) if np.isfinite(s).all() \
        else float("nan")
    sac = float(c * np.std(lg, ddof=1) / np.sqrt(n)) if n > 1 else 0.0
    te = np.array([b["tau_eff"] for b in B])
    eg = float(np.polyfit(te, lg - lg.mean(), 1)[0]) if n > 2 else float("nan")
    return dict(c=c, n=n, jk=jk, sac=sac, rms=float(np.std(lg, ddof=1))
                if n > 1 else 0.0, egim=eg,
                tmin=float(te.min()), tmax=float(te.max()))


def main(adlar=None):
    D = yukle(adlar)
    if not D:
        raise SystemExit("C_*.json yok")
    print(f"=== 167 ABLASYON HARİTASI === {len(D)} gaz: "
          f"{', '.join(sorted(D))}\n")

    # ---------------- gaz künyeleri ----------------
    print("gaz       aile      λ     τ_ust  erfc        N      T        "
          "dres      σ_ds    σ_X̃    σ_Ĉ    g_E    g_X    g=g_Eg_X²")
    for ad in sorted(D, key=lambda a: (ORT.KUNYE.get(a, {}).get("aile", "z"), a)):
        d = D[ad]
        A = d["artik"]
        k = d.get("kunye", {})
        print(f"{ad:9s} {k.get('aile','?'):9s} {str(d['lam']):5s} "
              f"{str(d['tau_ust']):5s}  {str(d['pen']):11s} {d['N']:7d} "
              f"{d['T']:8.0f} {d['dres']:.3e} {d['sigds']:.4f} "
              f"{d['sigX']:.4f} {d['sigC']:.4f} {A['gE']:.4f} {A['gX']:.4f} "
              f"{A['gE']*A['gX']**2:.4f}")

    # ---------------- bant tabloları ----------------
    for ad in sorted(D):
        d = D[ad]
        B = hukum(d)
        Bg = hukum(d, gevsek=True)
        print(f"\n--- {ad}: {len(B)}/{len(Bg)} hüküm bandı "
              f"(R≥{R_MIN}, SNR≥{SNR_MIN}) ---")
        print("  τ_eff   KALİB_u2±jk       W_amp   W_X    W_pos   "
              "c_o=K/(WaWx)±jk   R_bant   SNR")
        for b in Bg:
            m = "*" if b in B else " "
            print(f" {m}{b['tau_eff']:.4f} {b['KALIB_u2']:.4f}±"
                  f"{b.get('sKALIB_u2', float('nan')):.4f}   {b['W_amp']:.4f} "
                  f"{b['W_X']:.4f} {b['W_pos']:.4f}  {b['c_u2']:.4f}±"
                  f"{b.get('sc_u2', float('nan')):.4f}    {b['R_bant']:.4f} "
                  f"{b['SNR']:8.1f}")

    # ---------------- c(gaz) — bütün üyeler ----------------
    print("\n=== c(gaz), YASA ÜYESİ BAZINDA ===")
    print(f"{'gaz':9s} {'n':>3s} " + "".join(f"{u[0][:14]:>16s}"
                                             for u in UYELER))
    R = {}
    for ad in sorted(D, key=lambda a: (ORT.KUNYE.get(a, {}).get("aile", "z"), a)):
        B = hukum(D[ad])
        row = {}
        for nm, fn in UYELER:
            row[nm] = c_gaz(B, fn)
        R[ad] = row
        print(f"{ad:9s} {len(B):3d} " + "".join(
            f"{row[nm]['c']:9.4f}±{row[nm]['sac']:5.4f}" for nm, _ in UYELER))

    # ---------------- ABLASYON EKSENLERİ ----------------
    # ORTAK PENCERE: lo ∈ [0.52, 0.68] (dört bant). c(τ) τ ile yükseldiği
    # için eksenler ancak AYNI bant kümesinde karşılaştırılabilir; K070 ve
    # T/4 dilimleri yüksek bantları kaybeder.
    RO = {ad: {nm: c_gaz(hukum(D[ad], lo_max=0.68), fn)
               for nm, fn in UYELER} for ad in D}
    for et, RR in (("KENDİ hüküm bantları", R),
                   ("ORTAK PENCERE lo ≤ 0.68", RO)):
        print(f"\n############ EKSENLER — {et} ############")
        print("\n=== EKSEN (a) GENLİK ÖLÇEĞİ λ ===")
        eks_yaz(D, RR, ["L070", "L085", "Hkeskin", "L115"],
                lambda ad: D[ad]["lam"], "λ")
        print("\n=== EKSEN (b) PENCERE T ===")
        for kok, gr in (("Hkeskin", ["HkT4a", "HkT4b", "HkT4c", "HkT4d",
                                     "HkT2a", "HkT2b", "Hkeskin"]),
                        ("son", ["snT4a", "snT4d", "snT2a", "snT2b", "son"])):
            eks_yaz(D, RR, gr,
                    lambda ad: f"T={D[ad]['T']:.0f}/N={D[ad]['N']}", "T")
        print("\n=== EKSEN (c) MERDİVEN KESİMİ ===")
        eks_yaz(D, RR, ["K070", "K090", "Hkeskin", "E060", "HA4"],
                lambda ad: str(D[ad]["pen"] or f"keskin≤{D[ad]['tau_ust']}"),
                "kesim")

    # ---------------- ÜÇ TABAN GAZI ----------------
    print("\n=== T1 — HASSAS TABAN (üç gaz) ===")
    for nm, _ in UYELER:
        vs = [R[a][nm]["c"] for a in ("son", "Hkeskin", "HA4") if a in R]
        if len(vs) == 3 and all(np.isfinite(vs)):
            print(f"  {nm:18s} son={vs[0]:.4f} Hkeskin={vs[1]:.4f} "
                  f"HA4={vs[2]:.4f}  ortak={np.exp(np.mean(np.log(vs))):.4f} "
                  f"yayılım={max(vs)/min(vs):.3f}")
    # ---------------- BANT ŞEKLİNİN TEK SAYISI: γ_eff -----------------
    # kalib(τ) ≈ c₀·exp(−γ τ²)  (iki parametreli, hüküm bantlarında)
    # DW adayları:  γ_amp = ½π²σ_ds² , γ_X = 2π²σ_X̃² , γ_pos = 2π²σ_Ĉ²
    print("\n=== BANT ŞEKLİ: γ_eff  (kalib ∝ e^{−γτ²}) ↔ DW adayları ===")
    print("  gaz        n   c₀      **γ_eff**  γ_amp   γ_X     γ_pos  "
          "γ_amp+γ_X   γ_eff/γ_X  artık rms%   α_amp  α_X   α_pos  (kalib ∝ W^α)")
    GAM = {}
    for ad in sorted(D, key=lambda a: (ORT.KUNYE.get(a, {}).get("aile", "z"), a)):
        d = D[ad]
        B = hukum(d)
        if len(B) < 3:
            continue
        x = np.array([b["tau_eff"] ** 2 for b in B])
        y = np.log([b["KALIB_u2"] for b in B])
        p = np.polyfit(x, y, 1)
        res = y - np.polyval(p, x)
        ga = 0.5 * np.pi ** 2 * d["sigds"] ** 2
        gx = 2 * np.pi ** 2 * d["sigX"] ** 2
        gp = 2 * np.pi ** 2 * d["sigChat"] ** 2
        alf = {}
        for k in ("W_amp", "W_X", "W_pos"):
            lw = np.log([b[k] for b in B])
            alf[k] = float(np.polyfit(lw, y, 1)[0])
        GAM[ad] = dict(c0=float(np.exp(p[1])), gam=float(-p[0]), ga=float(ga),
                       gx=float(gx), gp=float(gp), alfa=alf,
                       rms=float(np.std(res, ddof=1)))
        print(f"  {ad:10s} {len(B):2d} {np.exp(p[1]):7.4f} {-p[0]:9.4f} "
              f"{ga:7.4f} {gx:7.4f} {gp:7.4f} {ga+gx:9.4f} "
              f"{-p[0]/gx:10.4f} {100*np.std(res, ddof=1):9.2f}   "
              f"{alf['W_amp']:5.3f} {alf['W_X']:5.3f} {alf['W_pos']:5.3f}")

    # ---------------- GRAM SIZINTISI HİPOTEZİ -------------------------
    # Model alanı E ≈ Λ_E·(η'nın çizgi izdüşümü) ⇒ kalib ≈ θ/(Λ_E Λ_X²).
    # Λ − 1, tarağın YOĞUNLUK dalgalanmasıyla (u ∝ λ) doğrusal olmalı.
    print("\n=== GRAM SIZINTISI: Λ = 1/g ve kalib/(1/Λ_EΛ_X²) ===")
    print("  gaz        λ     Λ_E=1/g_E  Λ_X=1/g_X  1/(Λ_EΛ_X²)  ⟨KALİB_u2⟩  "
          "θ=oran   (Λ_E−1)/λ (Λ_X−1)/λ  ⟨W_ampW_X⟩  c=θ/(ΛΛ²⟨W⟩)")
    for ad in sorted(D, key=lambda a: (ORT.KUNYE.get(a, {}).get("aile", "z"), a)):
        d = D[ad]
        B = hukum(d)
        if not B:
            continue
        A = d["artik"]
        LE, LX = 1 / A["gE"], 1 / A["gX"]
        gc = 1.0 / (LE * LX ** 2)
        kb = float(np.exp(np.mean(np.log([b["KALIB_u2"] for b in B]))))
        wb = float(np.exp(np.mean(np.log([b["W_amp"] * b["W_X"] for b in B]))))
        lam = d["lam"] or 1.0
        print(f"  {ad:10s} {lam:.2f}  {LE:9.4f}  {LX:9.4f}  {gc:10.4f}  "
              f"{kb:10.4f}  {kb/gc:7.4f}  {(LE-1)/lam:8.4f} {(LX-1)/lam:8.4f}"
              f"  {wb:10.4f}  {gc/wb:10.4f}")

    # ---------------- T4 — KAPANIŞ (165'in T1 tablosu) ----------------
    taban3 = [a for a in ("son", "Hkeskin", "HA4") if a in D]
    if len(taban3) == 3:
        print("\n=== T4 — KAPANIŞ: 165'in T1 tablosu, TAM parametresiz "
              "(20 hüküm bandı = 160'ın sağlık kuralı) ===")
        B20 = [(a, b) for a in taban3 for b in hukum(D[a], gevsek=True)]
        print(f"  bant sayısı: {len(B20)} "
              + " ".join(f"{a}:{sum(1 for x,_ in B20 if x==a)}"
                         for a in taban3))
        for nm, fn in UYELER:
            cc = c_gaz([b for _, b in B20], fn)["c"]
            o, w = 0, 0
            say = {a: [0, 0] for a in taban3}
            for a, b in B20:
                r = cc * fn(b) * b["Ps2"] / b["Ms2"] if b["Ms2"] else np.nan
                iy = np.isfinite(r) and abs(r - 1.0) <= 0.25
                o += bool(iy)
                say[a][1] += 1
                say[a][0] += bool(iy)
                if 0.53 < b["tau_eff"] < 0.71:
                    w += bool(iy)
            npen = sum(1 for _, b in B20 if 0.53 < b["tau_eff"] < 0.71)
            print(f"  {nm:18s} c={cc:.4f}  ±%25: **{o}/{len(B20)}**  "
                  f"(sağlıklı pencere {w}/{npen})  " +
                  " ".join(f"{a}:{say[a][0]}/{say[a][1]}" for a in taban3))
        # 165'in kendi bant-bant kalibrasyonu (referans)
        o0 = sum(1 for _, b in B20
                 if b["Ms2"] and abs(b["KALIB_u2"] * b["Ps2"] / b["Ms2"] - 1)
                 <= 0.25)
        print(f"  {'[165] u2 bant-bant':18s} c=  —     ±%25: **{o0}/{len(B20)}**")

    json.dump({a: {k: {kk: (float(vv) if isinstance(vv, (int, float)) else vv)
                       for kk, vv in v.items()} for k, v in r.items()}
               for a, r in R.items()},
              open(SCR / "ozet_c.json", "w"), indent=1)
    print(f"\n-> {SCR/'ozet_c.json'}")


def eks_yaz(D, R, grup, keyfn, etiket):
    grup = [g for g in grup if g in R and np.isfinite(R[g][UYELER[0][0]]["c"])]
    if not grup:
        print("  (veri yok)")
        return
    print(f"  {'gaz':9s} {etiket:>20s} {'n':>3s} " +
          "".join(f"{nm[:12]:>14s}" for nm, _ in UYELER))
    rad = next((a for a in ("Hkeskin", "son") if a in grup), grup[0])
    ref = R[rad]
    for ad in grup:
        r = R[ad]
        print(f"  {ad:9s} {str(keyfn(ad)):>20s} {r[UYELER[0][0]]['n']:3d} " +
              "".join(f"{r[nm]['c']:14.4f}" for nm, _ in UYELER))
    print(f"  --- oran / {rad} ---")
    for ad in grup:
        r = R[ad]
        print(f"  {ad:9s} {str(keyfn(ad)):>20s} {'':>3s} " +
              "".join(f"{r[nm]['c']/ref[nm]['c']:14.4f}"
                      for nm, _ in UYELER))


if __name__ == "__main__":
    main(sys.argv[1:] or None)
