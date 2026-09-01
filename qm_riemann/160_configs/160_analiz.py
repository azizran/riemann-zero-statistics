"""
160 — ANALİZ
B1  δ'nın TAM muhasebesi, bant bant: δ vs A·S vs K1 vs K1+K2
B2  Kesmenin kapalı formu: kümülant açılımı (Aκ1, −A³κ3/6) ve sönüm
B3  σ'nın kimliği: kuadratür regresyonu + teşhis korelasyonları
B4  σ kanalının serisi: s0, A·s1, −A²s2/2, −A³s3/6 vs ölçülen K2
B5  a ve b'nin HASSAS türetimi (158'in 4 pencere × 4 taban ensemble'ı)
B6  Taban dayanıklılığı (0.40 ↔ 0.52)
"""
import importlib
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
C = importlib.import_module("160_cekirdek")
SCR = C.SCR
TWO_PI = 2 * np.pi
FOUR_PI = 4 * np.pi

GAZ = ["son", "orta", "keskin", "A4", "P1"]
ETI = {"son": "gerçek-son", "orta": "gerçek-orta", "keskin": "keskin",
       "A4": "A4", "P1": "P1(GUE)"}

# ölçülen mekanizma eğrileri: ad -> (bant anahtarı, jackknife anahtarı)
EGRI = {"φ": ("phi", "sPhi_jk"),
        "δ": ("delta", "sDelta_jk"),
        "K1+K2": ("Kfull", "sKf_jk"),
        "K1": ("K1", "sK1_jk"),
        "M1(dsΔ)": ("M1", "sM1_jk"),
        "A·S": (("AS", "S_tam"), "sAS_jk")}


def yukle(veri, taban, g):
    p = SCR / f"S_{veri}_t{taban}_{g}.json"
    return json.load(open(p)) if p.exists() else None


def bantlar(d, lo=None, hi=None):
    B = [b for b in d["bantlar"] if b.get("olculdu")]
    if lo is not None:
        B = [b for b in B if b["tau"] >= lo - 1e-9]
    if hi is not None:
        B = [b for b in B if b["tau"] <= hi + 1e-9]
    return B


def saglik(B):
    """159 §7'nin kuralı, birebir: |Γ|≤1, |φ|≤2.8, |τ_eff−τ̄|≤yarım bant;
    ilk bozulmadan SONRAKİ bantlar da düşer."""
    ok, kir = [], True
    for b in B:
        yw = 0.5 * (b["hi"] - b["lo"])
        iyi = (abs(b["absG"]) <= 1.0 and abs(b["phi"]) <= 2.8
               and abs(b["tau_eff"] - b["tau"]) <= yw + 1e-9)
        if kir and iyi:
            ok.append(True)
        else:
            kir = False
            ok.append(False)
    return ok


def al(b, key):
    return b["AS"][key[1]] if isinstance(key, tuple) else b[key]


# --------------------------------------------------------------- fit
def _poli(x, y, s, deg):
    x = np.asarray(x, float); y = np.asarray(y, float)
    s = np.asarray(s, float) if s is not None else None
    if s is None or not np.all(np.isfinite(s)) or np.any(s <= 0):
        w = np.ones_like(x)
    else:
        w = 1.0 / s
    return np.polyfit(x, y, deg, w=w), w


def kok(x, y, s, deg):
    c, _ = _poli(x, y, s, deg)
    r = np.roots(c)
    r = r[np.abs(r.imag) < 1e-9].real
    return float(r[np.argmin(np.abs(r - np.mean(x)))]) if len(r) else float("nan")


def egim(x, y, s, deg, t0):
    """(dy/dx|t0 , ½·d²y/dx²|t0). deg=2'de ikincisi = c₀ (158'in b'si)."""
    c, w = _poli(x, y, s, deg)
    d1 = np.polyder(c, 1); d2 = np.polyder(c, 2)
    rez = np.asarray(y) - np.polyval(c, x)
    dof = max(len(x) - (deg + 1), 1)
    chi = float(np.sum((rez * w) ** 2) / dof)
    return float(np.polyval(d1, t0)), float(0.5 * np.polyval(d2, t0)), chi


PENCERE = {"W-A": dict(lo=0.43, hi=0.61, deg=2),
           "W-D": dict(lo=0.43, hi=0.59, deg=2),
           "W-B": dict(merkez=True, deg=2),
           "W-C": dict(lo=0.43, hi=0.61, deg=3)}


def fit_seti(d, pen):
    """Bir (gaz,taban,pencere) için: τ₀(φ), a(φ), b(φ) ve her mekanizma
    eğrisinin τ₀(φ)'daki eğimi/eğriliği. 158'in konvansiyonu."""
    P = PENCERE[pen]
    if P.get("merkez"):
        B0 = bantlar(d, 0.43, 0.61)
        if len(B0) < 4:
            return None
        x0 = [b["tau_eff"] for b in B0]
        t00 = kok(x0, [b["phi"] for b in B0], [b["sPhi_jk"] for b in B0], 2)
        B = [b for b in bantlar(d)
             if -0.075 - 1e-9 <= b["tau_eff"] - t00 <= 0.105 + 1e-9]
    else:
        B = bantlar(d, P["lo"], P["hi"])
    if len(B) < P["deg"] + 2:
        return None
    deg = P["deg"]
    x = [b["tau_eff"] for b in B]
    sp = [b["sPhi_jk"] for b in B]
    yp = [b["phi"] for b in B]
    t0 = kok(x, yp, sp, deg)
    if not np.isfinite(t0):
        return None
    a_m, b_m, chi = egim(x, yp, sp, deg, t0)
    out = dict(n=len(B), t0=t0, a=a_m, b=b_m, chi=chi, deg=deg)
    for ad, (ky, sky) in EGRI.items():
        y = [al(b, ky) for b in B]
        s = [b[sky] for b in B]
        g1, g2, _ = egim(x, y, s, deg, t0)
        out[f"d_{ad}"] = g1
        out[f"b_{ad}"] = g2
        g1w, g2w, _ = egim(x, y, sp, deg, t0)     # ağırlık-eşleşmiş sürüm
        out[f"dw_{ad}"] = g1w
        out[f"bw_{ad}"] = g2w
    return out


# --------------------------------------------------------------- B1
def B1(veri, taban, g="t1"):
    d = yukle(veri, taban, g)
    if d is None:
        print(f"  [{veri} t{taban} {g}] YOK")
        return
    B = bantlar(d)
    sg = saglik(B)
    print()
    print("=" * 100)
    print(f"B1 — δ'NIN TAM MUHASEBESİ  ({ETI.get(veri,veri)}, taban {taban}, "
          f"ızgara {g}, {len(B)} bant)")
    print("=" * 100)
    print("  τ_eff  A·σX  |Γ|    δ_meas    A·S     K1      K2    K1+K2  |  "
          "δ/(A·S)  δ/K1  K1/(A·S) |  kapanış K1 %   kapanış K1+K2 %   sağ")
    R = []
    for b, s in zip(B, sg):
        AS = b["AS"]["S_tam"]
        k1, k2, kf = b["K1"], b["K2"], b["Kfull"]
        dl = b["delta"]
        A = b["A_eff"]
        asx = A * np.sqrt(max(b["KUM"]["A2k2"] / A ** 2, 0.0))
        e_as = 100 * (AS - dl) / abs(dl) if dl else np.nan
        e_k1 = 100 * (k1 - dl) / abs(dl) if dl else np.nan
        e_kf = 100 * (kf - dl) / abs(dl) if dl else np.nan
        R.append((b["tau_eff"], dl, AS, k1, k2, kf, e_as, e_k1, e_kf, s))
        print(f"  {b['tau_eff']:.4f} {asx:.3f} {b['absG']:.3f} {dl:+8.4f} "
              f"{AS:+8.4f} {k1:+8.4f} {k2:+7.4f} {kf:+8.4f}  | "
              f"{dl/AS if AS else np.nan:+7.3f} {dl/k1 if k1 else np.nan:+6.3f} "
              f"{k1/AS if AS else np.nan:+7.3f}  | {e_k1:+10.2f}      "
              f"{e_kf:+12.4f}      {'✓' if s else '✗'}")
    a = np.array([(r[6], r[7], r[8]) for r in R if r[9]], float)
    if len(a):
        print(f"\n  SAĞLIKLI {len(a)}/{len(B)} bantta δ'ya göre BAĞIL HATA (%):")
        for i, ad in enumerate(("A·S  (KALEM, 1. mertebe)",
                                "K1   (reel ağırlık, TAM kar. fonk.)",
                                "K1+K2 (kuadratür dahil = TAM)")):
            v = a[:, i]
            v = v[np.isfinite(v)]
            ic10 = int(np.sum(np.abs(v) <= 10))
            print(f"    {ad:36s} ort={np.mean(v):+9.3f}  "
                  f"medyan={np.median(v):+9.3f}  menzil=[{v.min():+.3f},"
                  f"{v.max():+.3f}]   ±%10 içinde {ic10}/{len(v)}")
    return B


# --------------------------------------------------------------- B2
def B2(veri="son", taban=0.4, g="t1"):
    d = yukle(veri, taban, g)
    if d is None:
        return
    B = bantlar(d)
    print()
    print("-" * 100)
    print(f"B2 — KESMENİN KAPALI FORMU: kümülant açılımı "
          f"({ETI.get(veri,veri)}, taban {taban}, {g})")
    print("-" * 100)
    print("  τ_eff   A·κ1(=A·S)  −A³κ3/6    toplam     K1(ölçülen)  "
          "artık(K1−toplam)  | log|E1|  −A²κ2/2  +A⁴κ4/24  toplam")
    for b in B:
        K = b["KUM"]
        t1 = K["A1k1"]; t3 = -K["A3k3"] / 6.0
        lg = np.log(max(b["absE1"], 1e-300))
        d2 = -K["A2k2"] / 2.0; d4 = K["A4k4"] / 24.0
        print(f"  {b['tau_eff']:.4f} {t1:+11.5f} {t3:+10.5f} {t1+t3:+10.5f} "
              f"{b['K1']:+12.5f} {b['K1']-(t1+t3):+16.5f}  | "
              f"{lg:+8.4f} {d2:+8.4f} {d4:+9.4f} {d2+d4:+8.4f}")
    print("\n  A·σ_X (açılım parametresi) ve mertebe oranları:")
    print("  τ_eff   A       κ1        κ2       κ3        A·σ_X   "
          "|A³κ3/6| / |A·κ1|   K1/(A·S)")
    for b in B:
        K = b["KUM"]; A = b["A_eff"]
        k1 = K["A1k1"] / A; k2 = K["A2k2"] / A ** 2; k3 = K["A3k3"] / A ** 3
        r = abs(K["A3k3"] / 6.0) / abs(K["A1k1"]) if K["A1k1"] else np.nan
        AS = b["AS"]["S_tam"]
        print(f"  {b['tau_eff']:.4f} {A:.3f} {k1:+9.5f} {k2:8.5f} {k3:+9.5f} "
              f"  {A*np.sqrt(max(k2,0)):.3f}   {r:14.3f}   "
              f"{b['K1']/AS if AS else np.nan:+9.3f}")


# --------------------------------------------------------------- B3
def B3(gazlar=("son", "keskin", "A4", "P1"), taban=0.4, g="t1"):
    print()
    print("=" * 100)
    print("B3 — σ'NIN KİMLİĞİ: kuadratür regresyonu σ ≈ α·ρ + β·Δρ  "
          "(öngörü α=tan(A/2), β=1/sinA)")
    print("=" * 100)
    print("  gaz         τ_eff   A      α_öl     α_öngörü   β_öl     "
          "β_öngörü    R²     | korel(σ,ρ) korel(σ,Δρ) korel(σ,X̃) "
          "korel(σ,dX̃/dn) korel(ρ,X̃) korel(ρ,dX̃/dn)")
    for v in gazlar:
        d = yukle(v, taban, g)
        if d is None:
            continue
        B = bantlar(d)
        sg = saglik(B)
        for b, s in zip(B, sg):
            if not s or "TANI" not in b:
                continue
            T = b["TANI"]
            print(f"  {ETI.get(v,v):11s} {b['tau_eff']:.4f} {b['A_eff']:.3f} "
                  f"{T['q_al']:+8.4f} {T['q_al_p']:+9.4f} {T['q_be']:+9.4f} "
                  f"{T['q_be_p']:+9.4f} {T['q_R2']:8.5f} | "
                  f"{T['c_sr']:+10.4f} {T['c_sdr']:+11.4f} {T['c_sx']:+11.5f} "
                  f"{T['c_sdx']:+14.5f} {T['c_rx']:+11.5f} {T['c_rdx']:+14.5f}")


# --------------------------------------------------------------- B4
def B4(veri="son", taban=0.4, g="t1"):
    d = yukle(veri, taban, g)
    if d is None:
        return
    B = bantlar(d)
    print()
    print("-" * 100)
    print(f"B4 — σ KANALININ SERİSİ ({ETI.get(veri,veri)}, taban {taban}, {g})")
    print("    Es = s0 + iA·s1 − A²s2/2 − iA³s3/6 + … ;  "
          "K2 ≈ Re[Es·conj(E1)]/|E1|²")
    print("-" * 100)
    print("  τ_eff   s0(=⟨σ⟩/⟨ρ⟩)  A·s1      −A²s2/2   −A³s3/6   Re[Es]   "
          "Im[Es]  | K2(ölçülen)  Re[Es conj E1]/|E1|²  K2(2.mert. seri)  "
          "seri/ölçülen")
    for b in B:
        K = b["KUM"]
        E1 = complex(b["E1_re"], b["E1_im"]); Es = complex(b["Es_re"], b["Es_im"])
        lin = float((Es * np.conj(E1)).real / abs(E1) ** 2)
        Es2 = complex(K["s0"] - K["A2s2"] / 2, K["A1s1"] - K["A3s3"] / 6)
        ser = float((Es2 * np.conj(E1)).real / abs(E1) ** 2)
        print(f"  {b['tau_eff']:.4f} {K['s0']:+12.2e} {K['A1s1']:+9.5f} "
              f"{-K['A2s2']/2:+10.5f} {-K['A3s3']/6:+9.5f} {Es.real:+9.5f} "
              f"{Es.imag:+8.5f} | {b['K2']:+11.5f} {lin:+21.5f} "
              f"{ser:+17.5f} {ser/b['K2'] if b['K2'] else np.nan:+12.3f}")


# --------------------------------------------------------------- B5
def B5():
    print()
    print("=" * 100)
    print("B5 — a ve b'nin HASSAS TÜRETİMİ")
    print("=" * 100)
    print("  Fit konvansiyonu 158/157: apsis τ_eff, ağırlık 1/σ_jk, "
          "τ₀ = φ'nin kökü (x ort.'na en yakın),")
    print("  a = dy/dτ|τ₀, b = ½·d²y/dτ²|τ₀ (kuadratikte c₀). "
          "Pencereler W-A/W-D/W-B/W-C = 158'inkiler.")

    # --- 5a: tek fit tablosu (taban 0.40, W-A) ---
    print("\n  5a. TEK FİT (g158, taban 0.40, W-A) — 159 §6c/6d'nin karşılığı")
    print("  gaz          n  τ₀(φ)   a(φ)   | 4π+dδ/dτ 4π+d(K1+K2)/dτ  "
          "4π+dK1/dτ  4π+dM1/dτ  4π+dA·S/dτ |  b(φ)   b(δ)  b(K1+K2) b(K1) "
          " b(M1)  b(A·S)")
    tek = {}
    for v in GAZ:
        d = yukle(v, 0.4, "g158")
        if d is None:
            continue
        f = fit_seti(d, "W-A")
        if f is None:
            continue
        tek[v] = f
        print(f"  {ETI.get(v,v):11s} {f['n']:2d} {f['t0']:.4f} {f['a']:7.3f} | "
              f"{FOUR_PI+f['d_δ']:8.3f} {FOUR_PI+f['d_K1+K2']:14.3f} "
              f"{FOUR_PI+f['d_K1']:10.3f} {FOUR_PI+f['d_M1(dsΔ)']:10.3f} "
              f"{FOUR_PI+f['d_A·S']:11.3f} | {f['b']:+6.2f} {f['b_δ']:+6.2f} "
              f"{f['b_K1+K2']:+8.2f} {f['b_K1']:+6.2f} "
              f"{f['b_M1(dsΔ)']:+6.2f} {f['b_A·S']:+7.2f}")

    print("\n  hata % (a):")
    print("  gaz          a(φ)  | 4π+dδ/dτ  %   4π+dK1/dτ  %   4π+dA·S/dτ  %"
          "   çıplak 4π  %")
    ks = list(tek)
    for v in ks:
        f = tek[v]
        e = lambda p: 100 * (p - f["a"]) / f["a"]
        p1 = FOUR_PI + f["d_δ"]; p2 = FOUR_PI + f["d_K1"]
        p3 = FOUR_PI + f["d_A·S"]
        print(f"  {ETI.get(v,v):11s} {f['a']:6.3f} | {p1:8.3f} {e(p1):+6.1f} "
              f"{p2:10.3f} {e(p2):+6.1f} {p3:11.3f} {e(p3):+6.1f} "
              f"{FOUR_PI:10.3f} {e(FOUR_PI):+6.1f}")
    am = np.array([tek[v]["a"] for v in ks])
    for nm, key in (("4π+dδ/dτ", "d_δ"), ("4π+dK1/dτ", "d_K1"),
                    ("4π+dM1/dτ", "d_M1(dsΔ)"), ("4π+dA·S/dτ", "d_A·S")):
        pr = np.array([FOUR_PI + tek[v][key] for v in ks])
        sr = (np.argsort(np.argsort(am)) == np.argsort(np.argsort(pr))).all()
        print(f"    korel(a_meas, {nm:12s}) = {np.corrcoef(am, pr)[0,1]:+.4f}  "
              f"sıralama {'AYNI' if sr else 'FARKLI'}")

    print("\n  b'nin KAPANIŞI (|b_pred|/|b(φ)|, işaret ✓/✗):")
    print("  gaz          b(φ)   | b(A·S)  kapanış% iş | b(K1)  kapanış% iş "
          "| b(δ)=b(K1+K2)  kapanış% iş")
    for v in ks:
        f = tek[v]
        o = []
        for key in ("b_A·S", "b_K1", "b_δ"):
            bp = f[key]
            o.append((bp, 100 * abs(bp) / abs(f["b"]),
                      "✓" if np.sign(bp) == np.sign(f["b"]) else "✗"))
        print(f"  {ETI.get(v,v):11s} {f['b']:+6.2f} | "
              + " | ".join(f"{a:+7.2f} {b:7.1f} {c}" for a, b, c in o))

    # --- 5b: gerçek gazın 32-fit ensemble'ı ---
    print("\n  5b. GERÇEK GAZ ENSEMBLE — 158'in 32 fiti "
          "(2 pencere[son,orta] × 4 taban × 4 fit-penceresi)")
    kol = ["a", "4π+dδ/dτ", "4π+d(K1+K2)/dτ", "4π+dK1/dτ", "4π+dM1/dτ",
           "4π+dA·S/dτ"]
    key = {"a": None, "4π+dδ/dτ": "d_δ", "4π+d(K1+K2)/dτ": "d_K1+K2",
           "4π+dK1/dτ": "d_K1", "4π+dM1/dτ": "d_M1(dsΔ)",
           "4π+dA·S/dτ": "d_A·S"}
    acc = {k: [] for k in kol}
    accw = {k: [] for k in kol}
    accb = {k: [] for k in ("b", "b(δ)", "b(K1)", "b(M1)", "b(A·S)")}
    accbw = {k: [] for k in ("b", "b(δ)", "b(K1)", "b(A·S)")}
    n = 0
    print("  veri taban pen   n  τ₀(φ)    a(φ)   4π+dδ/dτ  4π+dK1/dτ  "
          "4π+dA·S/dτ    b(φ)   b(δ)   b(K1)  b(A·S)")
    for v in ("son", "orta"):
        for t in (0.28, 0.34, 0.4, 0.46):
            d = yukle(v, t, "g158")
            if d is None:
                print(f"  {v:5s} {t:<5}  --- KOŞU YOK ---")
                continue
            for pen in ("W-A", "W-D", "W-B", "W-C"):
                f = fit_seti(d, pen)
                if f is None:
                    continue
                n += 1
                acc["a"].append(f["a"]); accw["a"].append(f["a"])
                for k in kol[1:]:
                    acc[k].append(FOUR_PI + f[key[k]])
                    accw[k].append(FOUR_PI + f["dw_" + key[k][2:]])
                accb["b"].append(f["b"]); accb["b(δ)"].append(f["b_δ"])
                accb["b(K1)"].append(f["b_K1"])
                accb["b(M1)"].append(f["b_M1(dsΔ)"])
                accb["b(A·S)"].append(f["b_A·S"])
                accbw["b"].append(f["b"]); accbw["b(δ)"].append(f["bw_δ"])
                accbw["b(K1)"].append(f["bw_K1"])
                accbw["b(A·S)"].append(f["bw_A·S"])
                print(f"  {v:5s} {t:<5} {pen} {f['n']:3d} {f['t0']:.4f} "
                      f"{f['a']:7.3f} {FOUR_PI+f['d_δ']:9.3f} "
                      f"{FOUR_PI+f['d_K1']:10.3f} {FOUR_PI+f['d_A·S']:11.3f} "
                      f"{f['b']:+7.2f} {f['b_δ']:+6.2f} {f['b_K1']:+6.2f} "
                      f"{f['b_A·S']:+7.2f}")
    if n:
        print(f"\n  → {n} fit.  REFERANS (158): a = 10.713 ± 0.059, "
              f"b = −7.22 ± 0.79, σ_konv = 0.075")
        print("  nicelik              ortalama    sd     menzil            "
              "fark(10.713)  kaç σ(0.059)")
        for k in kol:
            a = np.array(acc[k], float)
            print(f"  {k:20s} {a.mean():8.3f} {a.std(ddof=1):7.3f}  "
                  f"[{a.min():7.3f},{a.max():7.3f}]  {a.mean()-10.713:+9.3f}"
                  f"  {(a.mean()-10.713)/0.059:+9.2f}")
        print("\n  (ağırlık-EŞLEŞMİŞ sürüm: bütün eğriler σ_φ ile tartıldı — "
              "fit ağırlığı sistematiğini yalıtır)")
        for k in kol:
            a = np.array(accw[k], float)
            print(f"  {k:20s} {a.mean():8.3f} {a.std(ddof=1):7.3f}  "
                  f"[{a.min():7.3f},{a.max():7.3f}]  {a.mean()-10.713:+9.3f}"
                  f"  {(a.mean()-10.713)/0.059:+9.2f}")
        bm = np.mean(accb["b"])
        print("\n  b ensemble (32 fit):     ortalama    sd    |ort|/sd   "
              "negatif/n   kapanış%(aynı ensemble b(φ))  kapanış%(158 W-A −7.22)")
        for k in ("b", "b(δ)", "b(K1)", "b(M1)", "b(A·S)"):
            a = np.array(accb[k], float)
            print(f"  {k:20s} {a.mean():+9.3f} {a.std(ddof=1):7.3f} "
                  f"{abs(a.mean())/a.std(ddof=1):9.2f}  {int((a<0).sum())}/{len(a)}"
                  f"     {100*abs(a.mean())/abs(bm):19.1f} {100*abs(a.mean())/7.22:22.1f}")
        bmw = np.mean(accbw["b"])
        print("  (ağırlık-EŞLEŞMİŞ b: bütün eğriler σ_φ ile tartıldı — "
              "b(φ) ≡ b(δ) kimliğini yalıtır)")
        for k in ("b", "b(δ)", "b(K1)", "b(A·S)"):
            a = np.array(accbw[k], float)
            print(f"  {k:20s} {a.mean():+9.3f} {a.std(ddof=1):7.3f} "
                  f"{abs(a.mean())/a.std(ddof=1):9.2f}  {int((a<0).sum())}/{len(a)}"
                  f"     {100*abs(a.mean())/abs(bmw):19.1f}")
        # 158'in "8 fit, yalnız W-A" referansı (a = 10.717±0.067, b = −7.22±0.79)
        print("\n  5c. 158'in İKİNCİ referansı: yalnız W-A, 8 fit "
              "(son+orta × 4 taban) — a = 10.717 ± 0.067, b = −7.2229 ± 0.7909")
        A8 = {k: [] for k in kol}
        B8 = {k: [] for k in ("b", "b(δ)", "b(K1)", "b(A·S)")}
        for v in ("son", "orta"):
            for t in (0.28, 0.34, 0.4, 0.46):
                d = yukle(v, t, "g158")
                if d is None:
                    continue
                f = fit_seti(d, "W-A")
                if f is None:
                    continue
                A8["a"].append(f["a"])
                for k in kol[1:]:
                    A8[k].append(FOUR_PI + f[key[k]])
                B8["b"].append(f["b"]); B8["b(δ)"].append(f["b_δ"])
                B8["b(K1)"].append(f["b_K1"]); B8["b(A·S)"].append(f["b_A·S"])
        for k in kol:
            a = np.array(A8[k], float)
            if not len(a):
                continue
            print(f"  {k:20s} {a.mean():8.3f} ± {a.std(ddof=1):.3f}   "
                  f"fark(10.717) = {a.mean()-10.717:+.3f} = "
                  f"{(a.mean()-10.717)/0.067:+.2f} σ(0.067)")
        bm8 = np.mean(B8["b"])
        for k in ("b", "b(δ)", "b(K1)", "b(A·S)"):
            a = np.array(B8[k], float)
            print(f"  {k:20s} {a.mean():+8.3f} ± {a.std(ddof=1):.3f}   "
                  f"kapanış (|ort|/|b(φ)|) = {100*abs(a.mean())/abs(bm8):.1f}%")

        # 5d. gaz gaz konvansiyon bütçesi (taban 0.40, 4 pencere)
        print("\n  5d. GAZ GAZ (taban 0.40, 4 pencere) — ortalama ± pencere sd")
        print("  gaz          a(φ)          4π+dδ/dτ       4π+dK1/dτ     "
              " 4π+dA·S/dτ   |  b(φ)          b(δ)          b(K1)        "
              " b(A·S)")
        for v in GAZ:
            d = yukle(v, 0.4, "g158")
            if d is None:
                continue
            F = [fit_seti(d, p) for p in ("W-A", "W-D", "W-B", "W-C")]
            F = [f for f in F if f]
            if not F:
                continue
            def ms(g):
                a = np.array([g(f) for f in F], float)
                return f"{a.mean():7.3f}±{a.std(ddof=1):.3f}"
            print(f"  {ETI.get(v,v):11s} {ms(lambda f: f['a'])} "
                  f"{ms(lambda f: FOUR_PI+f['d_δ'])} "
                  f"{ms(lambda f: FOUR_PI+f['d_K1'])} "
                  f"{ms(lambda f: FOUR_PI+f['d_A·S'])} | "
                  f"{ms(lambda f: f['b'])} {ms(lambda f: f['b_δ'])} "
                  f"{ms(lambda f: f['b_K1'])} {ms(lambda f: f['b_A·S'])}")


# --------------------------------------------------------------- B6
def B6():
    print()
    print("=" * 100)
    print("B6 — TABAN DAYANIKLILIĞI (0.40 ↔ 0.52), örtüşen bantlar")
    print("=" * 100)
    for v, g in (("son", "t1"), ("son", "g158"), ("A4", "g158")):
        a = yukle(v, 0.4, g); b = yukle(v, 0.52, g)
        if a is None or b is None:
            continue
        A = {x["tau"]: x for x in bantlar(a)}
        Bd = {x["tau"]: x for x in bantlar(b)}
        ort = sorted(set(A) & set(Bd))
        if not ort:
            continue
        print(f"\n  {ETI.get(v,v)} / {g} — {len(ort)} örtüşen bant")
        print("  τ̄     δ(.40)   δ(.52)  K1(.40)  K1(.52)  K2(.40)  K2(.52)  "
              "δ/K1(.40) δ/K1(.52)  δ/A·S(.40) δ/A·S(.52)")
        for t in ort:
            x, y = A[t], Bd[t]
            print(f"  {t:.2f} {x['delta']:+8.4f} {y['delta']:+8.4f} "
                  f"{x['K1']:+8.4f} {y['K1']:+8.4f} {x['K2']:+8.4f} "
                  f"{y['K2']:+8.4f} {x['delta']/x['K1']:+9.3f} "
                  f"{y['delta']/y['K1']:+9.3f} "
                  f"{x['delta']/x['AS']['S_tam']:+10.3f} "
                  f"{y['delta']/y['AS']['S_tam']:+10.3f}")


def B7():
    """Muhasebenin ÖZETİ: iki çarpan, K3(dsΔ) ve kapanış seviyeleri."""
    print()
    print("=" * 100)
    print("B7 — ÖZET: iki çarpan (kesme × kuadratür), K3 ve kapanış seviyesi")
    print("=" * 100)
    print("  gaz/taban        sağlıklı  |  K1/(A·S) kesme çarpanı   "
          "δ/K1 kuadratür çarpanı   δ/(A·S) çarpım | K2/δ payı % | "
          "maks|δ−(K1+K2)|/|δ|  maks|δ−M0|/|δ|")
    for v, t, g in (("son", 0.4, "t1"), ("son", 0.52, "t1"),
                    ("keskin", 0.4, "t1"), ("A4", 0.4, "t1"),
                    ("P1", 0.4, "t1"), ("son", 0.4, "g158"),
                    ("orta", 0.4, "g158")):
        d = yukle(v, t, g)
        if d is None:
            continue
        B = bantlar(d)
        sg = saglik(B)
        H = [b for b, s in zip(B, sg) if s]
        if not H:
            continue
        # sıfır geçişinden uzak bantlar: |δ| ve |A·S| bandın en büyüğünün %20'si üstü
        mx = max(abs(b["delta"]) for b in H)
        U = [b for b in H if abs(b["delta"]) > 0.2 * mx
             and abs(b["AS"]["S_tam"]) > 0.05 * mx and abs(b["K1"]) > 0.05 * mx]
        if not U:
            U = H
        f1 = np.array([b["K1"] / b["AS"]["S_tam"] for b in U])
        f2 = np.array([b["delta"] / b["K1"] for b in U])
        f3 = np.array([b["delta"] / b["AS"]["S_tam"] for b in U])
        p2 = np.array([100 * b["K2"] / b["delta"] for b in U])
        e3 = max(abs(b["delta"] - b["Kfull"]) / abs(b["delta"]) for b in H)
        e4 = max(abs(b["delta"] - b["M0"]) / abs(b["delta"]) for b in H)
        print(f"  {ETI.get(v,v):11s} t{t} {g:5s} {len(H):2d}/{len(B):2d} | "
              f"{np.median(f1):+7.3f} [{f1.min():+.2f},{f1.max():+.2f}]   "
              f"{np.median(f2):+7.3f} [{f2.min():+.2f},{f2.max():+.2f}]   "
              f"{np.median(f3):+7.3f} | {np.median(p2):+8.1f} | "
              f"{e3:.2e}          {e4:.2e}")
    print("\n  (çarpan özetleri δ'nın sıfır geçişinden UZAK bantlarda — "
          "|δ| bandın maksimumunun %20'si üstü; oran orada tanımlı.)")


if __name__ == "__main__":
    for v in ("son", "keskin", "A4", "P1"):
        B1(v, 0.4, "t1")
    B1("son", 0.52, "t1")
    B2("son", 0.4, "t1")
    B2("A4", 0.4, "t1")
    B3()
    B4("son", 0.4, "t1")
    B4("A4", 0.4, "t1")
    B5()
    B6()
    B7()
