"""
158 — ANALİZ: a = dφ_Γ/dτ|_{τ₀} gaz gaz, taban taban
====================================================
Girdi : scratchpad/158/F_<veri>_t<taban>.json      (158_kos.py)
        scratchpad/157/R_<veri>_t<taban>_ince.json (157'nin gerçek koşuları,
                                                    yalnız DENETİM için)
Çıktı : scratchpad/158/analiz_cikti.txt

FİT KONVANSİYONU — 157'nin bolum_E'siyle BİREBİR:
  apsis  x = τ_eff (bant ORTASI değil; 155'in F2'si)
  ordinat y = φ_Γ = arg Γ_rot
  ağırlık w = 1/σ_φ(jackknife)   [np.polyfit(..., w=1/σ)]
  derece 2 (kuadratik):  φ = a(τ−τ₀) + b(τ−τ₀)²
      τ₀ = köklerden x ortalamasına en yakın olanı
      a  = c₁ + 2c₀τ₀      (τ₀'DAKİ TÜREV — 155'in ham c₁'i DEĞİL)
      b  = c₀
  ΔAIC = (χ²_kuad + 2·3) − (χ²_dogr + 2·2)

PENCERELER
  W-A  τ̄ ∈ [0.43, 0.61]                 ← 157'nin penceresi (BİRİNCİL)
  W-B  τ_eff − τ₀ ∈ [−0.075, +0.105]     ← SIFIR-MERKEZLİ (gazın kendi τ₀'ına
                                            göre; a'nın pencere yerleşiminden
                                            gelip gelmediğini ayırır)
  W-C  W-A üstünde KÜBİK fit             ← eğrilik modelinin a'ya sızması

HATA
  istatistik : 8 gruplu jackknife — her replika için (teff_jk, phi_jk) ile
               fit tekrarlanır (155/157'nin hata mimarisi)
  sistematik : taban ekseni boyunca yayılım (yarı-menzil ve sd)
"""
import json
import sys
from pathlib import Path

import numpy as np

SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/158")
SCR157 = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
              "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/157")
SCR155 = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
              "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/155")

TABAN = (0.28, 0.34, 0.40, 0.46, 0.52)
# korelasyon merdiveni sırasıyla (155 §4: τ₀ = 0.4885 → 0.509+)
MERDIVEN = ("keskin", "A4", "J14", "J26", "N5", "N5z", "P1", "son", "orta")
GERCEK = ("son", "orta")
OUT = []

# 157'nin gerçek-gaz referansı (157_Rlad_yarisi_RAPOR.md §5c / (iv))
REF_A, REF_EA = 10.759, 0.114
REF_B, REF_EB = -7.121, 0.776


def yaz(*a):
    s = " ".join(str(x) for x in a)
    OUT.append(s)
    print(s, flush=True)


def yukle(kok=SCR, kalip="F_{v}_t{t}.json"):
    D = {}
    for v in MERDIVEN + ("P0",):
        for t in TABAN:
            for isim in (kalip.format(v=v, t=t), kalip.format(v=v, t=f"{t:g}")):
                p = kok / isim
                if p.exists():
                    D[(v, t)] = json.load(open(p))
                    break
    return D


def bant(d):
    return [b for b in d["bantlar"] if b.get("olculdu")]


# ------------------------------------------------------------------- fit
def _fit(xx, yy, ee, derece):
    """Ağırlıklı polinom fit → (τ₀, a, b, χ², n, katsayılar)."""
    m = np.isfinite(xx) & np.isfinite(yy)
    xx, yy, ee = xx[m], yy[m], ee[m]
    if len(xx) < derece + 2:
        return None
    ee = np.where(np.isfinite(ee) & (ee > 0), ee, np.nanmedian(ee))
    c = np.polyfit(xx, yy, derece, w=1.0 / ee)
    r = np.roots(c)
    r = np.array([z.real for z in r if abs(z.imag) < 1e-9])
    if len(r) == 0:
        return None
    t0 = float(r[np.argmin(np.abs(r - xx.mean()))])
    d1 = np.polyder(c)
    a = float(np.polyval(d1, t0))
    b = float(0.5 * np.polyval(np.polyder(d1), t0))   # ½φ''(τ₀)
    chi2 = float(np.sum(((yy - np.polyval(c, xx)) / ee)**2))
    return dict(t0=t0, a=a, b=b, chi2=chi2, n=len(xx), c=c)


def secim(bs, pencere, t0=None):
    """Bantları pencereye göre seç; (x, y, e) döndür."""
    x = np.array([b["tau_eff"] for b in bs], float)
    y = np.array([b["phi"] for b in bs], float)
    e = np.array([b["sPhi_jk"] for b in bs], float)
    tb = np.array([b["tau"] for b in bs], float)
    if pencere == "A":
        m = (tb >= 0.43 - 1e-9) & (tb <= 0.61 + 1e-9)
    elif pencere == "D":          # daha dar üst sınır (pencere sağlamlığı)
        m = (tb >= 0.43 - 1e-9) & (tb <= 0.59 + 1e-9)
    elif pencere == "B":
        m = (x - t0 >= -0.075) & (x - t0 <= 0.105)
    else:
        raise ValueError(pencere)
    return x[m], y[m], e[m], m


def olc(bs, pencere="A", derece=2, njack=8):
    """Bir koşudan (τ₀, a, b) + jackknife hataları."""
    if pencere == "B":
        # önce W-A ile kaba τ₀, sonra sıfır-merkezli pencere
        ilk = olc(bs, "A", derece=2, njack=0)
        if ilk is None:
            return None
        t0k = ilk["t0"]
        for _ in range(3):
            x, y, e, m = secim(bs, "B", t0k)
            f = _fit(x, y, e, derece)
            if f is None:
                return None
            if abs(f["t0"] - t0k) < 1e-6:
                break
            t0k = f["t0"]
    else:
        x, y, e, m = secim(bs, pencere, None)
        f = _fit(x, y, e, derece)
        if f is None:
            return None
    # doğrusal karşılaştırma (aynı bantlar)
    fl = _fit(x, y, e, 1)
    f["daic"] = (f["chi2"] + 2 * (derece + 1)) - (fl["chi2"] + 4) \
        if fl else float("nan")
    f["chi2dof_dog"] = fl["chi2"] / max(len(x) - 2, 1) if fl else float("nan")
    f["a_dog"] = fl["a"] if fl else float("nan")
    f["t0_dog"] = fl["t0"] if fl else float("nan")
    f["chi2dof"] = f["chi2"] / max(f["n"] - derece - 1, 1)
    # ---- pencere τ₀'ı KAPSIYOR mu? (sentetiklerde yüksek tabanlarda
    #      τ₀ ≈ 0.489 pencerenin altında kalır → a EKSTRAPOLASYONDUR)
    f["n_neg"] = int(np.sum(x < f["t0"]))
    f["n_poz"] = int(np.sum(x > f["t0"]))
    f["kenar"] = float(f["t0"] - x.min())     # <0 ise τ₀ pencerenin dışında
    f["xmin"], f["xmax"] = float(x.min()), float(x.max())
    # ---- jackknife
    rep = {"t0": [], "a": [], "b": []}
    idx = np.where(m)[0]
    for k in range(njack):
        xk = np.array([bs[i]["teff_jk"][k] for i in idx], float)
        yk = np.array([bs[i]["phi_jk"][k] for i in idx], float)
        ek = np.array([bs[i]["sPhi_jk"] for i in idx], float)
        g = _fit(xk, yk, ek, derece)
        if g is None:
            continue
        for q in rep:
            rep[q].append(g[q])
    for q in rep:
        v = np.array(rep[q], float)
        v = v[np.isfinite(v)]
        f["e_" + q] = (float(np.sqrt((len(v) - 1) / len(v)
                                     * np.sum((v - v.mean())**2)))
                       if len(v) >= 4 else float("nan"))
    return f


def ozet(vals):
    v = np.array([x for x in vals if np.isfinite(x)], float)
    if len(v) == 0:
        return float("nan"), float("nan"), float("nan")
    return float(v.mean()), float(v.std(ddof=1)) if len(v) > 1 else 0.0, \
        float(v.max() - v.min())


# ================================================================== V
def bolum_V(D):
    """DENETİM 1: 157'nin a/b tablosu, 157'nin KENDİ json'larından, bu
    scriptin fitiyle yeniden üretiliyor mu?
       DENETİM 2: 158'in gerçek-gaz φ'si = 157'nin φ'si mi (ortak bantlar)?"""
    yaz("\n" + "=" * 100)
    yaz("V — DENETİM: fitin ve ölçümün 157'ye karşı yeniden üretimi")
    yaz("=" * 100)
    yaz("\nV1: 157'nin §6 tablosu (τ̄≤0.61, kuadratik) — 157'nin json'ları, "
        "bu scriptin fiti")
    yaz("    pen  taban | 157: τ₀      a        b     | burada: τ₀      a"
        "        b     | Δa      Δb")
    ref = {("son", 0.28): (0.4999, 10.790, -5.890),
           ("son", 0.34): (0.5043, 10.689, -6.877),
           ("son", 0.40): (0.5088, 10.615, -7.759),
           ("son", 0.46): (0.5120, 10.724, -7.502),
           ("son", 0.52): (0.5127, 11.013, -7.272),
           ("orta", 0.28): (0.4995, 10.814, -6.495),
           ("orta", 0.34): (0.5047, 10.695, -7.123),
           ("orta", 0.40): (0.5087, 10.656, -7.801),
           ("orta", 0.46): (0.5121, 10.753, -8.338),
           ("orta", 0.52): (0.5124, 10.845, -6.152)}
    D157 = yukle(SCR157, "R_{v}_t{t}_ince.json")
    D157t = yukle(SCR157, "R_{v}_t{t}_tau0.json")
    D157.update({k: v for k, v in D157t.items() if k not in D157})
    da, db = [], []
    for (v, t), (r0, ra, rb) in ref.items():
        if (v, t) not in D157:
            yaz(f"    {v:4s} {t:.2f}  | 157 json YOK")
            continue
        bs = [b for b in bant(D157[(v, t)]) if b["tau"] <= 0.61 + 1e-9]
        f = _fit(np.array([b["tau_eff"] for b in bs]),
                 np.array([b["phi"] for b in bs]),
                 np.array([b["sPhi_jk"] for b in bs]), 2)
        da.append(f["a"] - ra); db.append(f["b"] - rb)
        yaz(f"    {v:4s} {t:.2f}  | {r0:.4f} {ra:7.3f} {rb:+7.3f} | "
            f"      {f['t0']:.4f} {f['a']:7.3f} {f['b']:+7.3f} | "
            f"{f['a']-ra:+.3f}  {f['b']-rb:+.3f}")
    yaz(f"    maks |Δa| = {np.max(np.abs(da)):.4f}   "
        f"maks |Δb| = {np.max(np.abs(db)):.4f}   "
        f"(157'nin tablosu 3 haneye yuvarlı)")

    yaz("\nV2: 158'in gerçek-gaz φ'si = 157'nin φ'si? (ortak τ̄, ortak taban)")
    yaz("    pen  taban  ortak bant  maks|Δφ|      maks|Δτ_eff|   maks|Δσ_φ|")
    for v in GERCEK:
        for t in TABAN:
            if (v, t) not in D or (v, t) not in D157:
                continue
            A = {round(b["tau"], 4): b for b in bant(D[(v, t)])}
            B = {round(b["tau"], 4): b for b in bant(D157[(v, t)])}
            ort = sorted(set(A) & set(B))
            if not ort:
                continue
            dphi = max(abs(A[k]["phi"] - B[k]["phi"]) for k in ort)
            dte = max(abs(A[k]["tau_eff"] - B[k]["tau_eff"]) for k in ort)
            dse = max(abs(A[k]["sPhi_jk"] - B[k]["sPhi_jk"]) for k in ort)
            yaz(f"    {v:4s} {t:.2f}  {len(ort):8d}   {dphi:.3e}   "
                f"{dte:.3e}   {dse:.3e}")

    yaz("\nV3: 158'in sentetik φ'si = 155'in φ'si? (ortak τ̄, ortak taban)")
    yaz("    gaz   taban  ortak bant  maks|Δφ|      maks|Δτ_eff|")
    D155 = yukle(SCR155, "tau0_{v}_t{t}_genis.json")
    for (v, t), d155 in sorted(D155.items()):
        if (v, t) not in D:
            continue
        A = {round(b["tau"], 4): b for b in bant(D[(v, t)])}
        B = {round(b["tau"], 4): b for b in bant(d155)}
        ort = sorted(set(A) & set(B))
        if not ort:
            continue
        dphi = max(abs(A[k]["phi"] - B[k]["phi"]) for k in ort)
        dte = max(abs(A[k]["tau_eff"] - B[k]["tau_eff"]) for k in ort)
        yaz(f"    {v:6s} {t:.2f} {len(ort):8d}   {dphi:.3e}   {dte:.3e}")


# ================================================================== T0
def bolum_T0(D):
    yaz("\n" + "=" * 100)
    yaz("T0 — KÜNYE: her koşunun gaz momentleri ve η artığı")
    yaz("=" * 100)
    yaz("   gaz  taban   L      σΔ²      κ₃Δ      çarp    σ_ds²   σ_η²    "
        "c₁      regresör  bant")
    for v in MERDIVEN + ("P0",):
        for t in TABAN:
            if (v, t) not in D:
                continue
            d = D[(v, t)]
            m = d["mom"]
            yaz(f"  {v:6s} {t:.2f} {d['L']:7.4f} {d['sA2']:.5f} "
                f"{m['k3']:+.5f} {m['carp']:+.3f} {d['s_ds']:7.4f} "
                f"{d['s_eta']:7.4f} {d['c1']:+.5f} {'':8s} "
                f"{len(bant(d)):4d}")


# ================================================================== T1
def bolum_T1(D, taban=0.40):
    yaz("\n" + "=" * 100)
    yaz(f"T1 — φ_Γ(τ_eff) HAM TAYFI  (taban {taban:.2f}, ızgara 0.02, "
        "negatif dal dahil)")
    yaz("=" * 100)
    tl = sorted({round(b["tau"], 4) for v in MERDIVEN
                 if (v, taban) in D for b in bant(D[(v, taban)])})
    yaz("   gaz   " + " ".join(f"{t:7.2f}" for t in tl))
    for v in MERDIVEN:
        if (v, taban) not in D:
            continue
        h = {round(b["tau"], 4): b for b in bant(D[(v, taban)])}
        yaz(f"  {v:6s} " + " ".join(
            (f"{h[t]['phi']:+7.4f}" if t in h else "      —") for t in tl))
    yaz("   (±σ_φ jackknife)")
    for v in MERDIVEN:
        if (v, taban) not in D:
            continue
        h = {round(b["tau"], 4): b for b in bant(D[(v, taban)])}
        yaz(f"  {v:6s} " + " ".join(
            (f"{h[t]['sPhi_jk']:7.4f}" if t in h else "      —") for t in tl))
    yaz("   |Γ|  (patlak denetimi: |Γ|>1.5 ölçüm sayılmaz)")
    for v in MERDIVEN:
        if (v, taban) not in D:
            continue
        h = {round(b["tau"], 4): b for b in bant(D[(v, taban)])}
        yaz(f"  {v:6s} " + " ".join(
            (f"{h[t]['absG']:7.3f}" if t in h else "      —") for t in tl))


# ================================================================== T2
def bolum_T2(D, pencere="A", derece=2, baslik=""):
    yaz("\n" + "=" * 100)
    yaz(f"T2{baslik} — (τ₀, a, b) her gaz × her taban   "
        f"[pencere W-{pencere}, derece {derece}]")
    yaz("=" * 100)
    yaz("   gaz  taban bant neg   τ₀      ±jk    |     a      ±jk   |     b   "
        "   ±jk   | χ²/dof  ΔAIC   hüküm")
    R = {}
    for v in MERDIVEN + ("P0",):
        for t in TABAN:
            if (v, t) not in D:
                continue
            bs = bant(D[(v, t)])
            f = olc(bs, pencere, derece)
            if f is None:
                yaz(f"  {v:6s} {t:.2f}  —    ÖLÇÜLEMEDİ (yeterli bant yok)")
                continue
            R[(v, t)] = f
            uy = "" if f["n_neg"] >= 1 else "  ‡EKSTRAPOLASYON (τ₀ pencerenin " \
                f"{abs(f['kenar']):.3f} altında)"
            yaz(f"  {v:6s} {t:.2f} {f['n']:4d} {f['n_neg']:3d} {f['t0']:.4f} "
                f"{f['e_t0']:.4f} | {f['a']:7.3f} {f['e_a']:5.3f} | "
                f"{f['b']:+7.3f} {f['e_b']:5.3f} | {f['chi2dof']:6.2f} "
                f"{f['daic']:+7.1f}  "
                f"{'EĞRİ' if f['daic'] < -2 else 'ayrılamaz'}{uy}")
    return R


# ================================================================== T3
def bolum_T3(R, etiket="", yalniz_kapsayan=True):
    yaz("\n" + "=" * 100)
    yaz(f"T3{etiket} — HÜKÜM TABLOSU: gaz | τ₀ | a | b   "
        "(taban ekseni üzerinden ortalama ± sd, yayılım)")
    if yalniz_kapsayan:
        yaz("     [YALNIZ τ₀'ı KAPSAYAN pencereler; ekstrapolasyon satırları "
            "(‡) dışarıda]")
    yaz("=" * 100)
    yaz("   gaz    n_taban |  τ₀      sd      yay   |   a       sd    yay% "
        " |   b       sd     yay%")
    sat = []
    for v in MERDIVEN:
        ks = [(v, t) for t in TABAN if (v, t) in R
              and (not yalniz_kapsayan or R[(v, t)]["n_neg"] >= 1)]
        if not ks:
            continue
        t0m, t0s, t0y = ozet([R[k]["t0"] for k in ks])
        am, as_, ay = ozet([R[k]["a"] for k in ks])
        bm, bs_, by = ozet([R[k]["b"] for k in ks])
        sat.append((v, len(ks), t0m, t0s, t0y, am, as_, ay, bm, bs_, by))
        yaz(f"  {v:6s} {len(ks):5d}   | {t0m:.4f} {t0s:.4f} {t0y:.4f} | "
            f"{am:7.3f} {as_:6.3f} {100*ay/abs(am):5.1f} | "
            f"{bm:+7.3f} {bs_:6.3f} {100*by/abs(bm):6.1f}")
    yaz(f"\n  157 gerçek referansı (2 pencere × 5 taban, τ̄≤0.61): "
        f"a = {REF_A:.3f} ± {REF_EA:.3f} · b = {REF_B:.3f} ± {REF_EB:.3f}")
    return sat


# ================================================================== T4
def bolum_T4(sat, D, R):
    yaz("\n" + "=" * 100)
    yaz("T4 — a MERDİVENİ İZLİYOR MU? (τ₀ ile korelasyon + boyutsuz oranlar)")
    yaz("=" * 100)
    yaz("   gaz     τ₀      a       b    |  σΔ²     P(s<0.3)*  |  a·σΔ²   "
        "a/(2π)   a·τ₀")
    P03 = {"son": 0.02420, "keskin": 0.00013, "A4": 0.00013, "N5": 0.02324,
           "N5z": 0.21955, "J14": 0.02196, "J26": 0.07170, "P1": 0.10752}
    xs, ys, bsq = [], [], []
    for (v, n, t0m, t0s, t0y, am, as_, ay, bm, bs_, by) in sat:
        sA2 = None
        for t in TABAN:
            if (v, t) in D:
                sA2 = D[(v, t)]["sA2"]; break
        p = P03.get(v, float("nan"))
        yaz(f"  {v:6s} {t0m:.4f} {am:7.3f} {bm:+7.3f} | {sA2:.5f} "
            f"{p:8.5f}  | {am*sA2:7.4f} {am/(2*np.pi):7.4f} {am*t0m:7.4f}")
        xs.append(t0m); ys.append(am); bsq.append(bm)
    xs = np.array(xs); ys = np.array(ys); bsq = np.array(bsq)
    if len(xs) > 2:
        yaz(f"\n  korel(τ₀, a) = {np.corrcoef(xs, ys)[0,1]:+.3f}   "
            f"korel(τ₀, b) = {np.corrcoef(xs, bsq)[0,1]:+.3f}   (n={len(xs)})")
        A = np.vstack([xs, np.ones_like(xs)]).T
        sl, ic = np.linalg.lstsq(A, ys, rcond=None)[0]
        yaz(f"  a = {sl:+.2f}·τ₀ + {ic:.2f}    "
            f"(τ₀'ın 0.0206'lık merdiveni a'da {abs(sl)*0.0206:.3f} "
            f"= %{100*abs(sl)*0.0206/ys.mean():.1f} demek)")
    yaz("\n  * P(s<0.3): 155 §2 künyesinden okundu (yeniden ölçülmedi).")
    yaz("  ** Boyutsuz oranlar YALNIZ KAYIT. Türetimsiz kapalı-form iddiası "
        "YOK (siren protokolü).")


# ================================================================== T6
def bolum_T6(sat, D, taban=0.40, tref=0.53):
    """a NEYİ İZLİYOR? — aynı koşudan okunan yapısal büyüklüklerle korelasyon.
    Bütün nicelikler taban 0.40'ta, τ̄ = tref bandından (aynı bant, her gaz)."""
    yaz("\n" + "=" * 100)
    yaz(f"T6 — a NEYİ İZLİYOR? aday korelatlar (hepsi taban {taban:.2f}, "
        f"τ̄ = {tref:.2f} bandı)")
    yaz("=" * 100)
    yaz("   gaz      a       τ₀     |  |Γ|     ρ      argM_emp  σΔ²      "
        "σ_ds²    σ_η²")
    kol = {k: [] for k in ("a", "t0", "absG", "rho", "argMe", "sA2", "sds",
                           "seta")}
    isim = []
    for (v, n, t0m, t0s, t0y, am, as_, ay, bm, bs_, by) in sat:
        if (v, taban) not in D:
            continue
        d = D[(v, taban)]
        h = [b for b in bant(d) if abs(b["tau"] - tref) < 1e-6]
        if not h:
            continue
        b0 = h[0]
        isim.append(v)
        for k, val in (("a", am), ("t0", t0m), ("absG", b0["absG"]),
                       ("rho", b0["rho"]), ("argMe", b0["argMe"]),
                       ("sA2", d["sA2"]), ("sds", d["s_ds"]),
                       ("seta", d["s_eta"])):
            kol[k].append(val)
        yaz(f"  {v:6s} {am:7.3f} {t0m:.4f} | {b0['absG']:.4f} "
            f"{b0['rho']:.4f} {b0['argMe']:+.5f}  {d['sA2']:.5f} "
            f"{d['s_ds']:.4f} {d['s_eta']:.4f}")
    if len(isim) >= 4:
        a = np.array(kol["a"])
        yaz("\n   korel(a, ·)  — n = %d gaz (gerçek dahil)" % len(isim))
        for k in ("t0", "absG", "rho", "argMe", "sA2", "sds", "seta"):
            v = np.array(kol[k])
            yaz(f"     {k:8s} {np.corrcoef(a, v)[0,1]:+.3f}")
        yaz("   (KAYIT — nedensellik iddiası yok; korelatlar tek τ bandından "
            "okundu.)")


# ================================================================== T5
def bolum_T5(RA, RB, RC, RD):
    yaz("\n" + "=" * 100)
    yaz("T5 — PENCERE SAĞLAMLIĞI: a, pencere/derece konvansiyonuna ne kadar "
        "asılı?")
    yaz("=" * 100)
    yaz("   gaz    W-A τ̄≤0.61 kuad | W-D τ̄≤0.59 kuad | W-B sıfır-merkezli |"
        " W-C W-A KÜBİK   | a yayılımı")
    yaz("            a       b     | a       b     | a       b     | "
        "a       b     | (dört pencere)")
    for v in MERDIVEN:
        row = []
        for R in (RA, RD, RB, RC):
            ks = [(v, t) for t in TABAN if (v, t) in R
                  and R[(v, t)]["n_neg"] >= 1]
            if not ks:
                row.append((float("nan"), float("nan")))
                continue
            row.append((ozet([R[k]["a"] for k in ks])[0],
                        ozet([R[k]["b"] for k in ks])[0]))
        aa = np.array([a for a, _ in row], float)
        aa = aa[np.isfinite(aa)]
        yay = (aa.max() - aa.min()) if len(aa) else float("nan")
        yaz(f"  {v:6s} " + " | ".join(f"{a:7.3f} {b:+7.3f}" for a, b in row)
            + f" | {yay:.3f}")

    # ---- toplu istatistik: gaz başına BÜTÜN pencere × taban fitleri
    yaz("\n   Bütün fitler (gaz başına 4 pencere × ≤4 taban, ‡ hariç):")
    yaz("     küme                       n   a ort     sd     min      max")
    for ad, gaz in (("gerçek (son+orta)", GERCEK),
                    ("sentetik (7 gaz)",
                     ("keskin", "A4", "J14", "J26", "N5", "N5z", "P1"))):
        aa = np.array([R[(v, t)]["a"] for R in (RA, RD, RB, RC)
                       for v in gaz for t in TABAN
                       if (v, t) in R and R[(v, t)]["n_neg"] >= 1], float)
        yaz(f"     {ad:24s} {len(aa):3d} {aa.mean():8.3f} {aa.std(ddof=1):6.3f}"
            f" {aa.min():8.3f} {aa.max():8.3f}")
    for ad, R, dg in (("W-A", RA, 2),):
        for k in ("a", "b", "t0"):
            v8 = np.array([R[(u, t)][k] for u in GERCEK for t in TABAN
                           if (u, t) in R and R[(u, t)]["n_neg"] >= 1], float)
            yaz(f"     GERÇEK {ad} 8 fit (‡ hariç): {k:3s} = {v8.mean():8.4f} ± "
                f"{v8.std(ddof=1):.4f}   menzil {v8.max()-v8.min():.4f}")
    for k in ("a", "t0"):
        v32 = np.array([R[(u, t)][k] for R in (RA, RD, RB, RC) for u in GERCEK
                        for t in TABAN
                        if (u, t) in R and R[(u, t)]["n_neg"] >= 1], float)
        yaz(f"     GERÇEK 32 fit:     {k:3s} = {v32.mean():8.4f} ± "
            f"{v32.std(ddof=1):.4f}   menzil {v32.max()-v32.min():.4f} "
            f"(%{100*(v32.max()-v32.min())/abs(v32.mean()):.1f})")

    # ---- AYIRICILIK: gerçek − A4 farkı taban taban
    yaz("\n   AYIRICILIK — (gerçek − A4) farkı taban ekseninde ne kadar "
        "kararlı? (W-A)")
    yaz("     taban |   Δτ₀      Δa")
    dt, da = [], []
    for t in TABAN:
        if ("A4", t) not in RA or not all((u, t) in RA for u in GERCEK):
            continue
        if RA[("A4", t)]["n_neg"] < 1 or any(RA[(u, t)]["n_neg"] < 1
                                             for u in GERCEK):
            continue      # ‡ ekstrapolasyon satırı — ayırıcılık ölçülemez
        gt = np.mean([RA[(u, t)]["t0"] for u in GERCEK])
        ga = np.mean([RA[(u, t)]["a"] for u in GERCEK])
        x, y = gt - RA[("A4", t)]["t0"], ga - RA[("A4", t)]["a"]
        dt.append(x); da.append(y)
        yaz(f"      {t:.2f} | {x:+.4f}  {y:+8.3f}")
    # ---- b'nin İŞARET kararlılığı (b sentetiklerde ölçülemiyor; işaret
    #      taşıdığı tek bilgi)
    yaz("\n   b İŞARET KARARLILIĞI (gaz başına 4 pencere × ≤4 taban, ‡ hariç)")
    yaz("     gaz     n_fit  b<0  b>0 | b ortalama | |ort|/sd")
    for v in MERDIVEN:
        bb = np.array([R[(v, t)]["b"] for R in (RA, RD, RB, RC) for t in TABAN
                       if (v, t) in R and R[(v, t)]["n_neg"] >= 1], float)
        if not len(bb):
            continue
        yaz(f"     {v:6s} {len(bb):5d} {int((bb < 0).sum()):4d} "
            f"{int((bb > 0).sum()):4d} | {bb.mean():+10.3f} | "
            f"{abs(bb.mean())/bb.std(ddof=1):8.2f}")

    dt, da = np.array(dt), np.array(da)
    for ad, w in (("Δτ₀", dt), ("Δa ", da)):
        yaz(f"     {ad}: ort {w.mean():+.4f}, taban menzili "
            f"{w.max()-w.min():.4f} = ortalamanın %"
            f"{100*(w.max()-w.min())/abs(w.mean()):.1f}'i; "
            f"en büyük/en küçük = {max(abs(w))/min(abs(w)):.2f} kat")


# ================================================================== T7
def bolum_T7(RA, RD, RB, RC, taban=None):
    """HÜKÜM: (H-M) evrensel mi, (H-K) merdiveni izliyor mu?

    Her gaz için a'nın KONVANSİYON BÜTÇESİ = taban yayılımı ⊕ pencere
    yayılımı; gerçekten ayrılma bu bütçe biriminde ölçülür. Ayrıca
    155'in τ₀ merdiven kesirleri (A4 → N5 → P1 → gerçek) a için
    yeniden hesaplanır.
    """
    yaz("\n" + "=" * 100)
    et = "taban ekseni ortalaması" if taban is None else f"taban {taban:.2f}"
    yaz(f"T7 — HÜKÜM ARİTMETİĞİ  [{et}]")
    yaz("=" * 100)

    def al(R, v):
        ks = [(v, t) for t in TABAN if (v, t) in R and R[(v, t)]["n_neg"] >= 1
              and (taban is None or t == taban)]
        if not ks:
            return None
        return ozet([R[k]["a"] for k in ks]), ozet([R[k]["b"] for k in ks]), \
            ozet([R[k]["t0"] for k in ks])

    yaz("\n  (a) a'nın KONVANSİYON BÜTÇESİ ve gerçekten ayrılma")
    yaz("     gaz      a      σ_taban  σ_pencere  σ_konv  | a−a(gerçek) | "
        "kaç σ_konv")
    ger = np.mean([al(RA, v)[0][0] for v in GERCEK])
    ger_s = np.mean([np.hypot(al(RA, v)[0][1],
                              np.std([al(R, v)[0][0] for R in
                                      (RA, RD, RB, RC)], ddof=1))
                     for v in GERCEK])
    for v in MERDIVEN:
        o = al(RA, v)
        if o is None:
            continue
        am, asd, _ = o[0]
        pen = [al(R, v)[0][0] for R in (RA, RD, RB, RC) if al(R, v)]
        psd = float(np.std(pen, ddof=1)) if len(pen) > 1 else 0.0
        kon = float(np.hypot(asd, psd))
        d = am - ger
        tot = float(np.hypot(kon, ger_s))
        yaz(f"    {v:6s} {am:7.3f}  {asd:6.3f}   {psd:6.3f}   {kon:6.3f}  | "
            f"{d:+8.3f}    | {abs(d)/tot:6.1f}")
    yaz(f"    [gerçek referansı = son/orta ortalaması {ger:.3f}, "
        f"σ_konv {ger_s:.3f}]")

    yaz("\n  (b) MERDİVEN KESİRLERİ: A4 → N5 → P1 → gerçek yolunun ne kadarı?")
    yaz("     nicelik | A4      N5      P1      gerçek | N5 %   P1 %")
    for ad, i in (("τ₀", 2), ("a", 0), ("b", 1)):
        dg = {}
        for v in ("A4", "N5", "P1"):
            o = al(RA, v)
            dg[v] = o[i][0] if o else float("nan")
        dg["ger"] = np.mean([al(RA, v)[i][0] for v in GERCEK])
        yol = dg["ger"] - dg["A4"]
        f5 = 100 * (dg["N5"] - dg["A4"]) / yol
        fp = 100 * (dg["P1"] - dg["A4"]) / yol
        yaz(f"     {ad:6s}  | {dg['A4']:7.4f} {dg['N5']:7.4f} {dg['P1']:7.4f} "
            f"{dg['ger']:7.4f} | {f5:5.1f}  {fp:5.1f}")
    yaz("     [155'in τ₀ merdiveni (taban 0.40, aralık (0.42,0.56]): "
        "N5 %25, P1 %52]")

    yaz("\n  (c) KONTROLLÜ ÇİFT N5 ↔ J14 (momentleri eşleşmiş; N5 korelasyonlu "
        "itme, J14 yapısız titreşim)")
    for ad, i in (("τ₀", 2), ("a", 0), ("b", 1)):
        n5, j14 = al(RA, "N5")[i], al(RA, "J14")[i]
        d = n5[0] - j14[0]
        s = float(np.hypot(n5[1], j14[1]))
        yaz(f"     Δ{ad:3s} = {d:+8.4f}   (taban σ: N5 {n5[1]:.4f}, "
            f"J14 {j14[1]:.4f} → {abs(d)/s if s else float('nan'):5.1f} σ_taban)")
    yaz("\n  (d) DOZ SINAVI: yapısız titreşim (A4→J14→J26) ve korelasyonlu "
        "itme (A4→N5→N5z)")
    for zin in (("A4", "J14", "J26"), ("A4", "N5", "N5z")):
        for ad, i in (("τ₀", 2), ("a", 0)):
            vals = [al(RA, v)[i][0] for v in zin]
            yaz(f"     {ad:3s}: " + " → ".join(f"{v}={x:.4f}"
                                               for v, x in zip(zin, vals)))


if __name__ == "__main__":
    D = yukle()
    yaz(f"158 ANALİZ — {len(D)} koşu yüklendi: "
        f"{sorted(set(k[0] for k in D))}")
    bolum_V(D)
    bolum_T0(D)
    bolum_T1(D, 0.40)
    bolum_T1(D, 0.28)
    RA = bolum_T2(D, "A", 2)
    sat = bolum_T3(RA)
    RB = bolum_T2(D, "B", 2, " (sıfır-merkezli pencere)")
    bolum_T3(RB, " (W-B)")
    RC = bolum_T2(D, "A", 3, " (KÜBİK)")
    bolum_T3(RC, " (W-C kübik)")
    RD = bolum_T2(D, "D", 2, " (dar pencere τ̄≤0.59)")
    bolum_T3(RD, " (W-D)")
    bolum_T5(RA, RB, RC, RD)
    bolum_T4(sat, D, RA)
    bolum_T6(sat, D, 0.40, 0.53)
    bolum_T6(sat, D, 0.46, 0.55)
    bolum_T7(RA, RD, RB, RC)
    bolum_T7(RA, RD, RB, RC, 0.40)
    (SCR / "analiz_cikti.txt").write_text("\n".join(OUT))
    print(f"\n-> {SCR/'analiz_cikti.txt'}")
