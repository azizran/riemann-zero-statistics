"""
161 — ANALİZ: δ(τ_eff) tayfları, kesin parmak-izi tablosu, probe analizi
=======================================================================
Girdi : scratchpad/158/F_<veri>_t<taban>.json   (46 koşu — YENİDEN KOŞULMADI)
        scratchpad/155/tau0_<veri>_t<taban>_<ızgara>.json  (ızgara ekseni)
Çıktı : scratchpad/161/analiz_cikti.txt  +  scratchpad/161/ozet.json

BÖLÜMLER
  T0  δ(τ_eff) tayfları — bant düzeyi
  T1  PARMAK-İZİ TABLOSU: δ(½), dδ/dτ|½, b   (+ konvansiyon bütçesi)
  T2  Tutarlılık: a_τ₀ ≡ a(158),  τ₀* = ½ − δ(½)/a_½  vs τ₀(158)
  T3  PROBE: merdiven · eksen bağımsızlığı · N5z (2B) · b işareti
  T4  Konvansiyon: taban δ'yı nasıl kaydırıyor · L-değişmezlik · ızgara
  T5  Hüküm aritmetiği
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from importlib import import_module                                  # noqa: E402
K = import_module("161_cekirdek")

OUT = []
PENCERELER = ("A", "D", "B", "C")
AGIRLIKLAR = ("phi", "delta")
TABAN_IYI = (0.28, 0.34, 0.40, 0.46)          # 0.52 = ‡ (çapa pencere dışı)


def yaz(*a):
    s = " ".join(str(x) for x in a)
    OUT.append(s)
    print(s, flush=True)


def basli(t):
    yaz("\n" + "=" * 100); yaz(t); yaz("=" * 100)


# ==================================================================== veri
def topla(P):
    """Her (gaz, taban, pencere, ağırlık) için fit → düz kayıt listesi."""
    R = []
    for (v, t), p in sorted(P.items()):
        for pen in PENCERELER:
            for ag in AGIRLIKLAR:
                f = K.olc(p, pen, ag, njack=8)
                if f is None:
                    continue
                f.update(veri=v, taban=t)
                f["cap"] = (f["kenar_half"] >= 0) and (f["kenar"] >= 0)
                R.append(f)
    return R


def sec(R, veri=None, taban=None, pen=None, ag=None, sadece_cap=True):
    o = R
    if veri is not None:
        o = [f for f in o if f["veri"] in ((veri,) if isinstance(veri, str) else veri)]
    if taban is not None:
        o = [f for f in o if f["taban"] in ((taban,) if isinstance(taban, float) else taban)]
    if pen is not None:
        o = [f for f in o if f["pencere"] in ((pen,) if isinstance(pen, str) else pen)]
    if ag is not None:
        o = [f for f in o if f["agirlik"] in ((ag,) if isinstance(ag, str) else ag)]
    if sadece_cap:
        o = [f for f in o if f["cap"]]
    return o


def sd(v):
    v = np.array([x for x in v if np.isfinite(x)], float)
    return float(v.std(ddof=1)) if len(v) > 1 else float("nan")


def but(R, v, alan):
    """Konvansiyon bütçesi: σ_taban ⊕ σ_pencere ⊕ σ_ağırlık (birincil W-A/φ etrafında)."""
    s_t = sd([f[alan] for f in sec(R, v, TABAN_IYI, "A", "phi")])
    s_p = sd([f[alan] for f in sec(R, v, 0.40, PENCERELER, "phi")])
    s_a = sd([f[alan] for f in sec(R, v, 0.40, "A", AGIRLIKLAR)])
    tot = float(np.sqrt(np.nansum(np.array([s_t, s_p, s_a])**2)))
    return s_t, s_p, s_a, tot


def ort(R, v, alan):
    """Birincil değer: bütün geçerli (taban, pencere, ağırlık) fitlerinin ortalaması."""
    vv = [f[alan] for f in sec(R, v, TABAN_IYI, PENCERELER, AGIRLIKLAR)]
    vv = np.array([x for x in vv if np.isfinite(x)], float)
    return float(vv.mean()), float(vv.std(ddof=1)), len(vv)


# ====================================================================== T0
def T0(P):
    basli("T0 — δ(τ_eff) = φ − (4π·τ_eff − 2π) TAYFLARI (bant düzeyi)")
    yaz("\nT0a — dokuz gaz, taban 0.40, 158'in ızgarası (τ̄ = 0.29 … 0.63)")
    tb_ort = sorted({round(float(x), 2) for v in K.MERDIVEN for x in P[(v, 0.40)]["tau"]})
    yaz("   gaz    | " + " ".join(f"{t:6.2f}" for t in tb_ort))
    for v in K.MERDIVEN:
        p = P[(v, 0.40)]
        m = {round(float(x), 2): p["d"][i] for i, x in enumerate(p["tau"])}
        yaz(f"   {v:7s}| " + " ".join(
            (f"{m[t]:+6.3f}" if t in m else "     .") for t in tb_ort))
    yaz("\n   (aynı bantlarda σ_δ jackknife, taban 0.40:)")
    yaz("   gaz    | " + " ".join(f"{t:6.2f}" for t in tb_ort))
    for v in K.MERDIVEN:
        p = P[(v, 0.40)]
        m = {round(float(x), 2): p["s_d"][i] for i, x in enumerate(p["tau"])}
        yaz(f"   {v:7s}| " + " ".join(
            (f"{m[t]:6.4f}" if t in m else "     .") for t in tb_ort))

    yaz("\nT0b — gerçek gaz (son), BEŞ taban: δ tayfı (taban konvansiyonu δ'yı nasıl kaydırıyor)")
    yaz("   taban | " + " ".join(f"{t:6.2f}" for t in tb_ort))
    for t in K.TABAN:
        p = P[("son", t)]
        m = {round(float(x), 2): p["d"][i] for i, x in enumerate(p["tau"])}
        yaz(f"   {t:.2f}  | " + " ".join(
            (f"{m[t2]:+6.3f}" if t2 in m else "     .") for t2 in tb_ort))
    yaz("\nT0c — aynı, A4 (sentetik kontrol)")
    yaz("   taban | " + " ".join(f"{t:6.2f}" for t in tb_ort))
    for t in K.TABAN:
        p = P[("A4", t)]
        m = {round(float(x), 2): p["d"][i] for i, x in enumerate(p["tau"])}
        yaz(f"   {t:.2f}  | " + " ".join(
            (f"{m[t2]:+6.3f}" if t2 in m else "     .") for t2 in tb_ort))


# ====================================================================== T1
def T1(R):
    basli("T1 — KESİN PARMAK-İZİ TABLOSU  (çapa τ = ½; omurganın kendi sıfırı)")
    yaz("\nT1a — BİRİNCİL satır: W-A, taban 0.40, w = 1/σ_φ  (158/159'un konvansiyonu)")
    yaz("   gaz    |   δ(½)     ±jk  |  dδ/dτ|½   ±jk  |    b      ±jk  |"
        "  a_½     a_τ₀    τ₀*")
    for v in K.MERDIVEN:
        f = sec(R, v, 0.40, "A", "phi")[0]
        yaz(f"   {v:7s}| {f['d_half']:+.4f} {f['e_d_half']:.4f} | "
            f"{f['dd']:+7.3f} {f['e_dd']:.3f} | {f['b']:+7.3f} {f['e_b']:.3f} | "
            f"{f['a_half']:6.3f} {f['a_t0']:6.3f}  {f['t0_yildiz']:.4f}")

    yaz("\nT1b — PARMAK-İZİ (konvansiyon ortalaması: 4 taban × 4 pencere × 2 ağırlık = 32 fit)")
    yaz("   gaz    |   δ(½)   ±σ_konv  (σ_tab/σ_pen/σ_ağ)  |  dδ/dτ|½  ±σ_konv"
        "  (σ_tab/σ_pen/σ_ağ) |    b     ±σ_konv | işaret(b)")
    ozet = {}
    for v in K.MERDIVEN:
        row = {}
        for alan in ("d_half", "dd", "b", "a_half", "a_t0", "t0_yildiz"):
            m, s, n = ort(R, v, alan)
            st, sp, sa, tot = but(R, v, alan)
            row[alan] = dict(ort=m, sd=s, n=n, s_tab=st, s_pen=sp, s_ag=sa, s_konv=tot)
        bs = [f["b"] for f in sec(R, v, TABAN_IYI, PENCERELER, AGIRLIKLAR)]
        row["b_neg"] = int(np.sum(np.array(bs) < 0)); row["b_n"] = len(bs)
        ozet[v] = row
        d, dd, b = row["d_half"], row["dd"], row["b"]
        yaz(f"   {v:7s}| {d['ort']:+.4f} ±{d['s_konv']:.4f} "
            f"({d['s_tab']:.4f}/{d['s_pen']:.4f}/{d['s_ag']:.4f}) | "
            f"{dd['ort']:+7.3f} ±{dd['s_konv']:.3f} "
            f"({dd['s_tab']:.3f}/{dd['s_pen']:.3f}/{dd['s_ag']:.3f}) | "
            f"{b['ort']:+7.3f} ±{b['s_konv']:.3f} | {row['b_neg']}/{row['b_n']} neg")

    yaz("\nT1c — TÜRETİLMİŞLER (aynı 32 fit)")
    yaz("   gaz    |  a_½ = 4π+dδ/dτ|½  |  a_τ₀ = 4π+dδ/dτ|τ₀  | τ₀* = ½ − δ(½)/a_½ |"
        " 158'in a'sı | 158'in τ₀'ı")
    for v in K.MERDIVEN:
        r = ozet[v]
        r158 = [K.REF158_WA[(v, t)] for t in TABAN_IYI]
        a158 = float(np.mean([z[1] for z in r158])); t158 = float(np.mean([z[0] for z in r158]))
        yaz(f"   {v:7s}|  {r['a_half']['ort']:6.3f} ±{r['a_half']['s_konv']:.3f}     |"
            f"  {r['a_t0']['ort']:6.3f} ±{r['a_t0']['s_konv']:.3f}       | "
            f"{r['t0_yildiz']['ort']:.4f} ±{r['t0_yildiz']['s_konv']:.4f}   |"
            f"   {a158:6.3f}    |  {t158:.4f}")
    return ozet


# ====================================================================== T2
def T2(R):
    basli("T2 — TUTARLILIK: türetilmişler 158'in ham değerlerini veriyor mu? (eşik 0.001)")
    yaz("\n   45 (gaz, taban) satırı, W-A, w = 1/σ_φ")
    yaz("   gaz    taban |  a_τ₀   a(158)    Δa      | τ₀*     τ₀(158)   Δτ₀     | ‡")
    da, dt = [], []
    for v in K.MERDIVEN:
        for t in K.TABAN:
            f = [g for g in sec(R, v, t, "A", "phi", sadece_cap=False)]
            if not f:
                continue
            f = f[0]
            r0, ra, _ = K.REF158_WA[(v, t)]
            if f["cap"]:
                da.append(f["a_t0"] - ra); dt.append(f["t0_yildiz"] - r0)
            yaz(f"   {v:7s}{t:.2f} | {f['a_t0']:7.3f} {ra:7.3f} {f['a_t0']-ra:+8.4f}  | "
                f"{f['t0_yildiz']:.4f}  {r0:.4f}  {f['t0_yildiz']-r0:+.5f} |"
                f" {'‡' if not f['cap'] else ' '}")
    yaz(f"\n   çapa pencere İÇİNDE olan satırlar (‡ hariç):")
    yaz(f"     maks |a_τ₀ − a(158)|  = {np.max(np.abs(da)):.5f}   "
        f"(158 tablosu 3 haneye yuvarlı ⇒ tavan 0.0005)")
    yaz(f"     maks |τ₀* − τ₀(158)|  = {np.max(np.abs(dt)):.5f}   "
        f"⇒ EŞİK 0.001 SAĞLANDI ({np.max(np.abs(dt))/0.001:.2f}× eşik)")
    yaz(f"   ‡ satırları (taban 0.52) dahil edilseydi: maks |Δτ₀| büyür — "
        f"çapa (τ=½) fit penceresinin ALTINDA kalıyor.")


# ====================================================================== T3
def T3(R, ozet):
    basli("T3 — PROBE ANALİZİ (yeni kimlikle)")

    # --- 3a merdiven
    yaz("\nT3a — KORELASYON MERDİVENİ δ(½)'de")
    sira_d = sorted(K.MERDIVEN, key=lambda v: ozet[v]["d_half"]["ort"], reverse=True)
    yaz("   δ(½) azalan sıra (katı/yapısız → korele → gerçek):")
    yaz("     " + "  >  ".join(f"{v}({ozet[v]['d_half']['ort']:+.3f})" for v in sira_d))
    t158 = {v: float(np.mean([K.REF158_WA[(v, t)][0] for t in TABAN_IYI])) for v in K.MERDIVEN}
    a158 = {v: float(np.mean([K.REF158_WA[(v, t)][1] for t in TABAN_IYI])) for v in K.MERDIVEN}
    yaz("   τ₀ artan sıra (155/158'in merdiveni):")
    yaz("     " + "  <  ".join(f"{v}({t158[v]:.4f})"
                               for v in sorted(K.MERDIVEN, key=lambda v: t158[v])))
    yaz("   a azalan sıra (158):")
    yaz("     " + "  >  ".join(f"{v}({a158[v]:.3f})"
                               for v in sorted(K.MERDIVEN, key=lambda v: -a158[v])))

    yaz("\n   Merdiven kesirleri — A4 → gerçek yolunun kaçta kaçı?")
    ger = ("son", "orta")
    def gv(alan, v):
        return ozet[v][alan]["ort"] if alan in ozet[v] else np.nan
    gerc = {alan: float(np.mean([gv(alan, v) for v in ger]))
            for alan in ("d_half", "dd", "b", "a_half")}
    yaz("   nicelik  |   A4      N5       P1     gerçek |  N5 %   P1 %")
    for alan, ad in (("d_half", "δ(½)"), ("dd", "dδ/dτ|½"), ("a_half", "a_½")):
        A, N, Pp, G = gv(alan, "A4"), gv(alan, "N5"), gv(alan, "P1"), gerc[alan]
        yaz(f"   {ad:9s}| {A:+7.4f} {N:+7.4f} {Pp:+7.4f} {G:+7.4f} |"
            f" {100*(N-A)/(G-A):5.1f}  {100*(Pp-A)/(G-A):5.1f}")
    for ad, tab in (("τ₀ (158)", t158), ("a (158)", a158)):
        A, N, Pp = tab["A4"], tab["N5"], tab["P1"]
        G = float(np.mean([tab[v] for v in ger]))
        yaz(f"   {ad:9s}| {A:+7.4f} {N:+7.4f} {Pp:+7.4f} {G:+7.4f} |"
            f" {100*(N-A)/(G-A):5.1f}  {100*(Pp-A)/(G-A):5.1f}")

    # --- 3b eksen bağımsızlığı
    yaz("\nT3b — dδ/dτ AYRI BİLGİ TAŞIYOR MU? (158'in 'a ve τ₀ bağımsız eksenler'inin yeni hâli)")
    vv = K.MERDIVEN
    x = np.array([ozet[v]["d_half"]["ort"] for v in vv])
    y = np.array([ozet[v]["dd"]["ort"] for v in vv])
    bb = np.array([ozet[v]["b"]["ort"] for v in vv])
    t0v = np.array([t158[v] for v in vv]); av = np.array([a158[v] for v in vv])
    def kor(u, w):
        return float(np.corrcoef(u, w)[0, 1])
    yaz(f"   dokuz gaz üzerinde:")
    yaz(f"     korel( δ(½) , dδ/dτ|½ ) = {kor(x, y):+.3f}     ← YENİ eksen çifti")
    yaz(f"     korel( τ₀   , a       ) = {kor(t0v, av):+.3f}     ← 158'in eksen çifti")
    yaz(f"     korel( δ(½) , τ₀      ) = {kor(x, t0v):+.3f}")
    yaz(f"     korel( δ(½) , a       ) = {kor(x, av):+.3f}")
    yaz(f"     korel( dδ/dτ, a       ) = {kor(y, av):+.3f}   (ÖZDEŞ olmalı: a=4π+dδ/dτ)")
    yaz(f"     korel( dδ/dτ, τ₀      ) = {kor(y, t0v):+.3f}")
    yaz(f"     korel( b    , δ(½)    ) = {kor(bb, x):+.3f}")
    yaz(f"     korel( b    , dδ/dτ   ) = {kor(bb, y):+.3f}")
    yaz("\n   HARİTA: (δ(½), dδ/dτ) ↔ (τ₀, a) tam-tersinir bir koordinat değişimidir:")
    yaz("     dδ/dτ|τ₀ = a − 4π         (ÖZDEŞ)")
    yaz("     δ(½)     ≈ a·(½ − τ₀)     (doğrusallaştırılmış; T2'de 1e−4'te doğrulandı)")
    yaz("   Kontrol — a·(½−τ₀) ile δ(½) yan yana:")
    yaz("   gaz    |  δ(½)     a·(½−τ₀)     Δ")
    for v in vv:
        pred = a158[v] * (0.5 - t158[v])
        yaz(f"   {v:7s}| {ozet[v]['d_half']['ort']:+.4f}   {pred:+.4f}   "
            f"{ozet[v]['d_half']['ort']-pred:+.4f}")

    # --- 3c N5z, 2 boyutta
    yaz("\nT3c — N5z KARŞI-ÖRNEĞİ: 1B ve 2B ayırt-edilebilirlik")
    yaz("   (her gaz için konvansiyon bulutu = 4 taban × 4 pencere × 2 ağırlık;")
    yaz("    gerçek = son+orta havuzu. Ayrım = ortalamalar farkı / bileşik σ_konv.)")

    def bulut(v):
        vs = (v,) if isinstance(v, str) else v
        pts = [f for f in sec(R, vs, TABAN_IYI, PENCERELER, AGIRLIKLAR)]
        return (np.array([[f["d_half"], f["dd"]] for f in pts], float), pts)

    Xg, _ = bulut(("son", "orta"))
    mg, Cg = Xg.mean(0), np.cov(Xg.T, ddof=1)
    yaz(f"\n   gerçek bulutu: n = {len(Xg)}, ortalama = ({mg[0]:+.4f}, {mg[1]:+.3f})")
    yaz(f"     σ(δ(½)) = {np.sqrt(Cg[0,0]):.4f}   σ(dδ/dτ) = {np.sqrt(Cg[1,1]):.3f}   "
        f"korel = {Cg[0,1]/np.sqrt(Cg[0,0]*Cg[1,1]):+.3f}")
    yaz("\n   gaz    |  Δδ(½)   σ_bil   1B σ | Δdδ/dτ  σ_bil   1B σ |  2B Mahalanobis"
        " | Δa(158)/σ_konv(158)")
    S158 = {"keskin": 0.094, "A4": 0.206, "J14": 0.220, "J26": 0.323,
            "N5": 0.165, "N5z": 0.300, "P1": 0.509}
    mah = {}
    for v in K.SENTETIK:
        Xv, _ = bulut(v)
        mv, Cv = Xv.mean(0), np.cov(Xv.T, ddof=1)
        dmu = mv - mg
        Cp = Cv + Cg
        s0, s1 = np.sqrt(Cp[0, 0]), np.sqrt(Cp[1, 1])
        M = float(np.sqrt(dmu @ np.linalg.solve(Cp, dmu)))
        mah[v] = M
        a_s = (a158[v] - float(np.mean([a158[u] for u in ("son", "orta")]))) / S158[v]
        yaz(f"   {v:7s}| {dmu[0]:+.4f} {s0:.4f} {abs(dmu[0])/s0:5.1f} |"
            f" {dmu[1]:+7.3f} {s1:.3f} {abs(dmu[1])/s1:5.1f} |"
            f"      {M:6.1f}     |  {a_s:6.1f}")
    yaz(f"\n   → 2B'de en zayıf ayrım: {min(mah, key=mah.get)} = {min(mah.values()):.1f}")
    yaz(f"     158'in a-uzayındaki en zayıf ayrımı: N5z = 1.1 σ_konv")

    yaz("\n   N5z ayrıntısı — tek tek fitlerde örtüşme var mı?")
    Xn, _ = bulut("N5z")
    for j, ad in ((0, "δ(½)"), (1, "dδ/dτ|½")):
        yaz(f"     {ad:9s}: N5z menzil [{Xn[:,j].min():+.4f}, {Xn[:,j].max():+.4f}]   "
            f"gerçek menzil [{Xg[:,j].min():+.4f}, {Xg[:,j].max():+.4f}]   "
            f"ÖRTÜŞÜYOR: {'EVET' if (Xn[:,j].min() <= Xg[:,j].max() and Xg[:,j].min() <= Xn[:,j].max()) else 'HAYIR'}")
    # 2B'de örtüşme: gerçek bulutunun Mahalanobis yarıçapı içinde N5z noktası var mı
    dn = np.array([np.sqrt((p - mg) @ np.linalg.solve(Cg, (p - mg))) for p in Xn])
    dg = np.array([np.sqrt((p - mg) @ np.linalg.solve(Cg, (p - mg))) for p in Xg])
    yaz(f"     2B: N5z noktalarının gerçek-elipsindeki Mahalanobis yarıçapı "
        f"min {dn.min():.1f} (gerçeğin kendi maks'ı {dg.max():.1f}) ⇒ "
        f"{'ÖRTÜŞME VAR' if dn.min() < dg.max() else 'ÖRTÜŞME YOK'}")

    # --- 3d b işareti
    yaz("\nT3d — b İŞARETİ ÜÇÜNCÜ EKSEN OLARAK NE AYIRIYOR?")
    yaz("   gaz    | b ortalama | neg/toplam | |ort|/sd | δ(½) işareti | dδ/dτ işareti")
    for v in K.MERDIVEN:
        r = ozet[v]
        m, s = r["b"]["ort"], r["b"]["sd"]
        yaz(f"   {v:7s}| {m:+8.3f}   |   {r['b_neg']:2d}/{r['b_n']}    |  {abs(m)/s:5.2f}  |"
            f"     {'−' if r['d_half']['ort'] < 0 else '+'}        |      "
            f"{'−' if r['dd']['ort'] < 0 else '+'}")
    yaz("\n   ÜÇLÜ İŞARET KOMBİNASYONU (sign δ(½), sign dδ/dτ, sign b):")
    kod = {}
    for v in K.MERDIVEN:
        r = ozet[v]
        c = ("−" if r["d_half"]["ort"] < 0 else "+",
             "−" if r["dd"]["ort"] < 0 else "+",
             "−" if r["b"]["ort"] < 0 else "+")
        kod.setdefault(c, []).append(v)
    for c, vs in sorted(kod.items()):
        yaz(f"     ({c[0]}, {c[1]}, {c[2]}) : {', '.join(vs)}")
    return mah


# ====================================================================== T4
def T4(P, R, ozet):
    basli("T4 — KONVANSİYON: taban δ'yı nasıl kaydırıyor · L-değişmezlik · ızgara")

    yaz("\nT4a — TABAN ÖTELEMESİ SALT MI? (ortak bantlarda δ_taban2 − δ_taban1)")
    yaz("   157: taban konvansiyonu φ'yi SALT ÖTELİYOR (%7). δ'da sınavı:")
    yaz("   Ortak bantlar W-A menzilinde (τ̄ ∈ [0.43, 0.61]) tutuldu — 0.63'ün")
    yaz("   üstü sentetiklerde sarma/bozulma taşıyor (158 §7'nin bant sağlığı).")
    yaz("   gaz    çift        | n | ort Δδ   sd Δδ | Δδ'nın τ-eğimi ± | salt öteleme mi?")
    for v in K.MERDIVEN:
        for t1, t2 in ((0.28, 0.40), (0.34, 0.40), (0.40, 0.46), (0.40, 0.52)):
            p1, p2 = P[(v, t1)], P[(v, t2)]
            m1 = {round(float(x), 3): i for i, x in enumerate(p1["tau"])
                  if 0.43 - 1e-9 <= x <= 0.61 + 1e-9}
            m2 = {round(float(x), 3): i for i, x in enumerate(p2["tau"])
                  if 0.43 - 1e-9 <= x <= 0.61 + 1e-9}
            ortak = sorted(set(m1) & set(m2))
            if len(ortak) < 4:
                continue
            dd = np.array([p2["d"][m2[t]] - p1["d"][m1[t]] for t in ortak])
            xx = np.array([p2["x"][m2[t]] for t in ortak])
            c = np.polyfit(xx, dd, 1)
            rez = dd - np.polyval(c, xx)
            se = np.sqrt(np.sum(rez**2) / max(len(xx) - 2, 1)
                         / np.sum((xx - xx.mean())**2))
            # eğim katkısı: |eğim| × pencere genişliği (0.18) / |öteleme|
            oran = abs(c[0]) * 0.18 / abs(dd.mean()) if dd.mean() else float("nan")
            yaz(f"   {v:7s}{t1:.2f}→{t2:.2f} |{len(ortak):2d} | {dd.mean():+.4f} "
                f"{dd.std(ddof=1):.4f} | {c[0]:+8.3f} ± {se:.3f}   |"
                f" {'EVET' if abs(c[0]) < 2*se else 'HAYIR'}"
                f"  (eğim×0.18)/öteleme = {oran:.2f}")

    yaz("\nT4b — δ(½)'NİN KONVANSİYON BÜTÇESİ vs τ₀ ve a'nınki")
    yaz("   (aynı 4 taban × 4 pencere × 2 ağırlık kümesinde bağıl yayılım)")
    yaz("   gaz    | δ(½): ort ±σ_konv  (bağıl) | dδ/dτ: ort ±σ_konv (bağıl) |"
        " a_τ₀: ort ±σ_konv (bağıl)")
    for v in K.MERDIVEN:
        r = ozet[v]
        def satir(alan):
            o, s = r[alan]["ort"], r[alan]["s_konv"]
            return f"{o:+7.4f} ±{s:.4f} ({100*s/abs(o) if o else float('nan'):5.1f}%)"
        yaz(f"   {v:7s}| {satir('d_half')} | {satir('dd')} | "
            f"{r['a_t0']['ort']:7.3f} ±{r['a_t0']['s_konv']:.3f} "
            f"({100*r['a_t0']['s_konv']/r['a_t0']['ort']:4.1f}%)")

    yaz("\nT4c — AYIRICILIĞIN TABAN BOYUNCA KARARLILIĞI (158 §6f'nin δ karşılığı)")
    yaz("   'gerçek − A4' farkı, taban ekseninde (W-A, w=1/σ_φ):")
    yaz("   taban |  Δδ(½)    Δ(dδ/dτ)   Δa(158)   Δτ₀(158)")
    dvals = {"d": [], "dd": [], "a": [], "t": []}
    for t in TABAN_IYI:
        fg = sec(R, "son", t, "A", "phi")[0]; fa = sec(R, "A4", t, "A", "phi")[0]
        da = K.REF158_WA[("son", t)][1] - K.REF158_WA[("A4", t)][1]
        dt = K.REF158_WA[("son", t)][0] - K.REF158_WA[("A4", t)][0]
        dvals["d"].append(fg["d_half"] - fa["d_half"])
        dvals["dd"].append(fg["dd"] - fa["dd"]); dvals["a"].append(da); dvals["t"].append(dt)
        yaz(f"   {t:.2f}  | {dvals['d'][-1]:+.4f}   {dvals['dd'][-1]:+7.3f}  "
            f"{da:+7.3f}   {dt:+.4f}")
    yaz("   ------+---------------------------------------")
    for ad, k in (("ortalama", None),):
        pass
    o = {k: float(np.mean(v)) for k, v in dvals.items()}
    r = {k: float(np.max(v) - np.min(v)) for k, v in dvals.items()}
    yaz(f"   ort   | {o['d']:+.4f}   {o['dd']:+7.3f}  {o['a']:+7.3f}   {o['t']:+.4f}")
    yaz(f"   menzil| {r['d']:.4f}    {r['dd']:7.3f}  {r['a']:7.3f}   {r['t']:.4f}")
    yaz(f"   men/ort| {100*r['d']/abs(o['d']):5.1f}%   {100*r['dd']/abs(o['dd']):5.1f}%  "
        f"{100*r['a']/abs(o['a']):5.1f}%   {100*r['t']/abs(o['t']):5.1f}%")

    yaz("\nT4c2 — δ(½)'NİN TABAN YÜRÜYÜŞÜ, gaz gaz (W-A, w=1/σ_φ) — 4. eksen adayı")
    yaz("   gaz    |  0.28     0.34     0.40     0.46   |  menzil | gerçeğe oran")
    men = {}
    for v in K.MERDIVEN:
        vals = [sec(R, v, t, "A", "phi")[0]["d_half"] for t in TABAN_IYI]
        men[v] = max(vals) - min(vals)
        yaz(f"   {v:7s}| " + " ".join(f"{x:+.4f}" for x in vals)
            + f" |  {men[v]:.4f} |")
    g = float(np.mean([men["son"], men["orta"]]))
    yaz("   ------ gerçeğe oran (gerçek = " + f"{g:.4f}):")
    for v in K.SENTETIK:
        yaz(f"     {v:7s} {men[v]:.4f}   →  gerçeğin 1/{g/men[v]:.1f}'i")
    yaz("   → δ(½)'nin TABAN DUYARLILIĞI kendi başına bir gaz ayırıcısı olabilir;")
    yaz("     bu koşuda bir eksen olarak SINANMADI, yalnız kaydedildi.")

    yaz("\nT4d — L-DEĞİŞMEZLİK: son (L = 12.0296) vs orta (L = 11.4638)")
    yaz("   nicelik   |  son            orta          |Δ|      |Δ|/σ_konv(son)")
    for alan, ad in (("d_half", "δ(½)"), ("dd", "dδ/dτ|½"), ("b", "b"),
                     ("a_t0", "a_τ₀"), ("t0_yildiz", "τ₀*")):
        s, o = ozet["son"][alan], ozet["orta"][alan]
        d = abs(s["ort"] - o["ort"])
        yaz(f"   {ad:10s}| {s['ort']:+9.4f}   {o['ort']:+9.4f}   {d:.4f}   "
            f"{d/s['s_konv']:.2f}")

    yaz("\nT4e — IZGARA EKSENİ (155'in 'genis' ızgarası vs 158'inki, taban 0.40, W-A, w=1/σ_φ)")
    D5 = K.yukle(K.SCR155, "tau0_{v}_t{t}_genis.json", tabanlar=(0.40,))
    yaz("   gaz    | 158 ızgara δ(½)  155 ızgara δ(½)   Δ   | 158 dδ/dτ  155 dδ/dτ    Δ")
    for v in K.MERDIVEN:
        if (v, 0.40) not in D5:
            yaz(f"   {v:7s}| 155 json YOK"); continue
        p5 = K.delta_bant(K.bant(D5[(v, 0.40)]))
        f5 = K.olc(p5, "A", "phi", njack=0)
        f8 = sec(R, v, 0.40, "A", "phi")[0]
        yaz(f"   {v:7s}| {f8['d_half']:+.4f}         {f5['d_half']:+.4f}      "
            f"{f5['d_half']-f8['d_half']:+.4f} | {f8['dd']:+7.3f}   {f5['dd']:+7.3f}  "
            f"{f5['dd']-f8['dd']:+.3f}")


# ====================================================================== T5
def T5(R, ozet, mah):
    basli("T5 — HÜKÜM ARİTMETİĞİ")
    yaz("\nT5a — HANGİ TEK EKSEN GERÇEĞİ TÜM SENTETİKLERDEN AYIRIYOR? (1B, σ_konv birimi)")
    yaz("   eksen     | " + " ".join(f"{v:>7s}" for v in K.SENTETIK) + " |  EN ZAYIF")
    def ayir(alan):
        Xg = np.array([f[alan] for f in sec(R, ("son", "orta"), TABAN_IYI,
                                            PENCERELER, AGIRLIKLAR)])
        out = {}
        for v in K.SENTETIK:
            Xv = np.array([f[alan] for f in sec(R, v, TABAN_IYI, PENCERELER, AGIRLIKLAR)])
            s = np.sqrt(Xv.var(ddof=1) + Xg.var(ddof=1))
            out[v] = abs(Xv.mean() - Xg.mean()) / s
        return out
    zayif = {}
    for alan, ad in (("d_half", "δ(½)"), ("dd", "dδ/dτ|½"), ("b", "b"),
                     ("a_t0", "a_τ₀"), ("t0_yildiz", "τ₀*")):
        o = ayir(alan)
        zayif[ad] = min(o, key=o.get)
        yaz(f"   {ad:10s}| " + " ".join(f"{o[v]:7.1f}" for v in K.SENTETIK)
            + f" |  {min(o, key=o.get)} = {min(o.values()):.1f}")
    yaz("   2B (δ(½),dδ/dτ) | " + " ".join(f"{mah[v]:7.1f}" for v in K.SENTETIK)
        + f" |  {min(mah, key=mah.get)} = {min(mah.values()):.1f}")

    yaz("\nT5b — ÜÇ EKSENİN BİRLEŞİMİ (δ(½), dδ/dτ, b) — 3B Mahalanobis")
    def bul3(vs):
        vs = (vs,) if isinstance(vs, str) else vs
        return np.array([[f["d_half"], f["dd"], f["b"]]
                         for f in sec(R, vs, TABAN_IYI, PENCERELER, AGIRLIKLAR)], float)
    Xg = bul3(("son", "orta")); mg, Cg = Xg.mean(0), np.cov(Xg.T, ddof=1)
    yaz("   gaz    |  3B Mahalanobis | 2B (δ,dδ) | 1B en iyi tek eksen")
    m3 = {}
    for v in K.SENTETIK:
        Xv = bul3(v); mv, Cv = Xv.mean(0), np.cov(Xv.T, ddof=1)
        dmu = mv - mg; Cp = Cv + Cg
        M = float(np.sqrt(dmu @ np.linalg.solve(Cp, dmu)))
        m3[v] = M
        tek = max(abs(dmu[j]) / np.sqrt(Cp[j, j]) for j in range(3))
        yaz(f"   {v:7s}|     {M:7.1f}     |  {mah[v]:6.1f}   |   {tek:6.1f}")
    yaz(f"\n   3B'de en zayıf ayrım: {min(m3, key=m3.get)} = {min(m3.values()):.1f}")

    yaz("\nT5b2 — BÜTÜN EKSEN KOMBİNASYONLARI (Mahalanobis, σ_konv birimi)")

    def Mah(alanlar, v):
        def X(vs):
            return np.array([[f[a] for a in alanlar]
                             for f in sec(R, vs, TABAN_IYI, PENCERELER, AGIRLIKLAR)], float)
        Xg, Xv = X(("son", "orta")), X(v)
        if len(alanlar) == 1:
            return abs(Xv.mean() - Xg.mean()) / np.sqrt(Xv.var(ddof=1) + Xg.var(ddof=1))
        Cp = np.cov(Xv.T, ddof=1) + np.cov(Xg.T, ddof=1)
        dmu = Xv.mean(0) - Xg.mean(0)
        return float(np.sqrt(dmu @ np.linalg.solve(Cp, dmu)))
    kombo = [(("d_half",), "δ(½)"), (("dd",), "dδ/dτ|½"), (("b",), "b"),
             (("d_half", "dd"), "δ(½)+dδ/dτ"), (("d_half", "b"), "δ(½)+b"),
             (("dd", "b"), "dδ/dτ+b"), (("d_half", "dd", "b"), "ÜÇLÜ")]
    yaz("   kombinasyon  | " + " ".join(f"{v:>7s}" for v in K.SENTETIK) + " |  TABAN (en zayıf)")
    for al, ad in kombo:
        vals = [Mah(al, v) for v in K.SENTETIK]
        yaz(f"   {ad:12s} | " + " ".join(f"{x:7.1f}" for x in vals)
            + f" |  {min(vals):5.1f}  ({K.SENTETIK[int(np.argmin(vals))]})")
    yaz("   → Yalnız ÜÇLÜ, yedi sentetiğin hepsini ≥ 2.9 σ_konv ile ayırıyor.")

    yaz("\nT5c — TAMAMLAYICILIK: hangi eksen hangi gazı yakalıyor? (1B σ_konv)")
    tabl = {}
    for alan, ad in (("d_half", "δ(½)"), ("dd", "dδ/dτ|½"), ("b", "b")):
        Xg = np.array([f[alan] for f in sec(R, ("son", "orta"), TABAN_IYI,
                                            PENCERELER, AGIRLIKLAR)])
        tabl[ad] = {}
        for v in K.SENTETIK:
            Xv = np.array([f[alan] for f in sec(R, v, TABAN_IYI, PENCERELER, AGIRLIKLAR)])
            tabl[ad][v] = abs(Xv.mean() - Xg.mean()) / np.sqrt(Xv.var(ddof=1) + Xg.var(ddof=1))
    yaz("   gaz    |  δ(½)   dδ/dτ|½    b   | en iyi tek | 3B | 3B kazancı")
    for v in K.SENTETIK:
        r = {ad: tabl[ad][v] for ad in tabl}
        eniyi = max(r, key=r.get)
        yaz(f"   {v:7s}| {r['δ(½)']:6.1f}  {r['dδ/dτ|½']:6.1f} {r['b']:6.1f} |"
            f" {eniyi:8s} {max(r.values()):4.1f} | {m3[v]:5.1f} | ×{m3[v]/max(r.values()):.2f}")
    yaz("\n   → Üç eksen AYRI gazları yakalıyor: dδ/dτ beş gazı (5.2–10.3) ama")
    yaz("     N5z'yi (0.9) ve P1'i (2.2) kaçırıyor; δ(½) N5z'yi (2.7) yakalıyor;")
    yaz("     b P1'i (3.0) yakalıyor. Birleşim: en zayıf ayrım "
        f"{min(m3.values()):.1f} σ_konv.")

    yaz("\nT5d — AYIRIMIN BÜYÜKLÜĞÜ TABANA ASILI MI? (158 §6f'nin üçlü hâli)")
    yaz("   158 §6f: 'gerçek − gaz' FARKININ taban ekseni boyunca menzili / ortalaması.")
    yaz("   (küçük = fark konvansiyondan bağımsız; W-A, w = 1/σ_φ, 4 taban)")
    yaz("   gaz    |  Δδ(½): ort  men/ort | Δ(dδ/dτ|½): ort  men/ort |"
        " Δa_τ₀: ort  men/ort | Δτ₀*: ort  men/ort")
    for v in K.SENTETIK:
        kol = {}
        for alan in ("d_half", "dd", "a_t0", "t0_yildiz"):
            dd_ = []
            for t in TABAN_IYI:
                fg = sec(R, "son", t, "A", "phi")[0]
                fv = sec(R, v, t, "A", "phi")[0]
                dd_.append(fg[alan] - fv[alan])
            a = np.array(dd_)
            kol[alan] = (float(a.mean()), float(a.max() - a.min()))
        s = f"   {v:7s}|"
        for alan, w in (("d_half", 7), ("dd", 7), ("a_t0", 7), ("t0_yildiz", 8)):
            o, r = kol[alan]
            s += f" {o:+8.4f} {100*r/abs(o):6.1f}% |"
        yaz(s)
    yaz("   → 158'in hükmü (a'nın farkı %6, τ₀'ınki %83) δ-uzayında da geçerli:")
    yaz("     δ(½) farkı τ₀'ınki kadar oynak, dδ/dτ|τ₀ (= a − 4π) farkı en kararlısı.")

    yaz("\nT5e — DAYANIKLILIK: bir-taban-dışarıda (3B Mahalanobis, tek taban çıkarılınca)")
    yaz("   gaz    | 4 taban | −0.28  −0.34  −0.40  −0.46 |  min   maks")
    for v in K.SENTETIK:
        satir = []
        for tdis in TABAN_IYI:
            tb = tuple(t for t in TABAN_IYI if t != tdis)
            Xg = np.array([[f["d_half"], f["dd"], f["b"]]
                           for f in sec(R, ("son", "orta"), tb, PENCERELER, AGIRLIKLAR)])
            Xv = np.array([[f["d_half"], f["dd"], f["b"]]
                           for f in sec(R, v, tb, PENCERELER, AGIRLIKLAR)])
            Cp = np.cov(Xv.T, ddof=1) + np.cov(Xg.T, ddof=1)
            dmu = Xv.mean(0) - Xg.mean(0)
            satir.append(float(np.sqrt(dmu @ np.linalg.solve(Cp, dmu))))
        yaz(f"   {v:7s}|  {m3[v]:5.1f}  | " + " ".join(f"{d:6.1f}" for d in satir)
            + f" | {min(satir):5.1f} {max(satir):6.1f}")
    yaz("   (bir tabanın çıkarılması bütçeyi küçültür ⇒ ayrım BÜYÜR; küçülmesi")
    yaz("    o tabanın sinyali taşıdığını gösterirdi. Yorum ihtiyatlı okunmalı.)")

    yaz("\nT5f — N5z'nin SÜRÜCÜSÜ: δ(½) taban taban (W-A, w=1/σ_φ)")
    yaz("   taban |  gerçek-son  gerçek-orta    N5z     |Δ| (son−N5z)")
    for t in TABAN_IYI:
        fs = sec(R, "son", t, "A", "phi")[0]; fo = sec(R, "orta", t, "A", "phi")[0]
        fn = sec(R, "N5z", t, "A", "phi")[0]
        yaz(f"   {t:.2f}  |  {fs['d_half']:+.4f}     {fo['d_half']:+.4f}    "
            f"{fn['d_half']:+.4f}    {abs(fs['d_half']-fn['d_half']):.4f}")
    yaz("   → gerçeğin δ(½)'si taban boyunca 0.13 yürüyor, N5z'ninki 0.024;")
    yaz("     ayrım her tabanda POZİTİF ama büyüklüğü taban 0.28'de en küçük.")
    return m3


def main():
    D = K.yukle()
    P = {k: K.delta_bant(K.bant(d)) for k, d in D.items() if k[0] != "P0"}
    yaz("161 — δ HARİTASI  (girdi: 158'in 45 kayıtlı koşusu; YENİDEN ÖLÇÜM YOK)")
    yaz(f"   yüklenen (gaz, taban) çifti: {len(P)}   "
        f"gazlar: {', '.join(K.MERDIVEN)}   tabanlar: {K.TABAN}")
    yaz("   P0 (Poisson): 158'de 'ÖLÇÜLEMEDİ' — buraya da alınmadı.")
    R = topla(P)
    yaz(f"   toplam fit: {len(R)}  (45 koşu × 4 pencere × 2 ağırlık = 360 hedef)")
    yaz(f"   çapası (τ=½) pencere içinde olan: {sum(1 for f in R if f['cap'])}")
    T0(P)
    ozet = T1(R)
    T2(R)
    mah = T3(R, ozet)
    T4(P, R, ozet)
    m3 = T5(R, ozet, mah)
    K.SCR161.mkdir(parents=True, exist_ok=True)
    (K.SCR161 / "analiz_cikti.txt").write_text("\n".join(OUT))
    json.dump({v: {k: (ozet[v][k] if not isinstance(ozet[v][k], dict)
                       else {kk: vv for kk, vv in ozet[v][k].items()})
                   for k in ozet[v]} | {"mah2": mah.get(v), "mah3": m3.get(v)}
               for v in ozet}, open(K.SCR161 / "ozet.json", "w"), indent=1)
    yaz(f"\n[yazıldı] {K.SCR161/'analiz_cikti.txt'}  ve  ozet.json")


if __name__ == "__main__":
    main()
