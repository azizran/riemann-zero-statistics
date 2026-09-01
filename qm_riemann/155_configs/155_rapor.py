"""
155 — RAPOR ÜRETİCİ: bütün tablolar tek koşudan
===============================================
Kullanım: 155_rapor.py > scratchpad/155/rapor_cikti.txt
Raporun HER sayısı buradan gelir; elle girilen sayı yoktur.

τ₀ KONVANSİYONU (açıkça yazılır, her bileşeni D'de ölçülür):
  bin 0.02 · apsis = güç-ağırlıklı τ_eff · φ(τ)'ya KUADRATİK fit ·
  aralık (0.42,0.56] · hata = grup-sil jackknife (8 grup, bütün bantlar
  eşzamanlı); χ²>dof ise saçılmayla şişirilmiş sürüm de verilir.
"""
import importlib
import json
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
A = importlib.import_module("155_analiz")
SCR = A.SCR

SIRA = ["son", "orta", "dusuk", "keskin", "A4", "N5", "N5z", "J14", "J26",
        "P1", "P0", "sonS0.4-0.46", "sonS0.4-0.52", "dusukS0.4-0.46"]
ETIKET = {"son": "gerçek L=12.03", "orta": "gerçek L=11.46",
          "dusuk": "gerçek L=10.48", "keskin": "sent. keskin merdiven",
          "A4": "sent. erfc-0.68", "N5": "sent. erfc+itme (N5)",
          "N5z": "sent. erfc+itme ×8 (N5z)", "J14": "sent. erfc+titreşim .141",
          "J26": "sent. erfc+titreşim .260", "P1": "sent. GUE gap + boya",
          "P0": "sent. Poisson gap + boya",
          "sonS0.4-0.46": "SAHTE TABAN 0.40->0.46 (asal olmayan frekans)",
          "sonS0.4-0.52": "SAHTE TABAN 0.40->0.52 (asal olmayan frekans)",
          "dusukS0.4-0.46": "SAHTE TABAN, düşük-L"}


def ke(lo, hi, w):
    n = int(round((hi - lo) / w))
    return [round(lo + w * i, 6) for i in range(n + 1)]


def tum_dosyalar():
    """scratchpad'teki bütün tau0_*.json'ları (veri, taban, bantset) ile."""
    D = {}
    for p in sorted(SCR.glob("tau0_*.json")):
        m = re.match(r"tau0_(.+)_t([0-9.]+)_(\w+)\.json$", p.name)
        if not m:
            continue
        v, t, bs = m.group(1), float(m.group(2)), m.group(3)
        D[(v, t, bs)] = json.loads(p.read_text())
    return D


def havuz(D, v, t):
    """Aynı (veri,taban) için bütün bant setlerinin çizgilerini birleştir.

    Setler farklı ızgaralar; aynı q iki sette de geçebilir → q'ya göre
    tekilleştirilir (aynı q'nun cr/pow değerleri özdeştir, ölçüm
    frekansı yalnız q'ya bağlı).
    """
    H = {}
    for (vv, tt, bs), j in D.items():
        if vv != v or abs(tt - t) > 1e-9:
            continue
        for c in A.cizgiler(j):
            H.setdefault(c["q"], c)
    return sorted(H.values(), key=lambda c: c["tau"])


def t0(H, lo, hi, w=0.02, deg=2, apsis="teff"):
    return A.tau0_jk(H, ke(lo, hi, w), apsis, deg)


def sec(D, v, t, bs=None):
    for b in (bs,) if bs else ("ince", "genis", "kaba"):
        if (v, t, b) in D:
            return D[(v, t, b)]
    return None


def main():
    D = tum_dosyalar()
    ciftler = sorted({(v, t) for (v, t, _) in D},
                     key=lambda x: (SIRA.index(x[0]) if x[0] in SIRA else 99,
                                    x[1]))
    HAV = {c: havuz(D, *c) for c in ciftler}

    print("=" * 110)
    print("A) VERİ KÜNYELERİ   (κ₃Δ = ⟨dsΔ³⟩; çarp = κ₃/σ³)")
    print(f"{'veri':7s} {'taban':>5s} {'L':>8s} {'σΔ²':>8s} {'κ₃Δ':>9s} "
          f"{'çarp':>6s} {'σ_ds²':>7s} {'σ_η²':>7s} {'c₁':>8s} "
          f"{'ρ@0.55':>7s} {'|Γ|@.55':>7s}  açıklama")
    for (v, t) in ciftler:
        j = sec(D, v, t)
        m = j["mom"]
        b55 = [b for b in j["bantlar"] if b.get("olculdu")
               and abs(b["tau"] - 0.55) < 1e-6]
        r55 = b55[0]["rho"] if b55 else float("nan")
        g55 = b55[0]["absG"] if b55 else float("nan")
        print(f"{v:7s} {t:5.2f} {j['L']:8.4f} {m['sA2']:8.5f} "
              f"{m['k3']:+9.5f} {m['carp']:+6.3f} {j['s_ds']:7.4f} "
              f"{j['s_eta']:7.4f} {j['c1']:+8.5f} {r55:7.4f} {g55:7.3f}  "
              f"{ETIKET.get(v,'')}")

    print()
    print("=" * 110)
    print("B) φ(τ) — ORTAK BANT TABLOSU (0.02 ızgara, taban 0.40)")
    tg = [round(0.43 + 0.02 * i, 3) for i in range(7)]
    print(f"{'veri':8s} " + "".join(f"{('τ='+str(x)):>10s}" for x in tg))
    for (v, t) in ciftler:
        if abs(t - 0.40) > 1e-9:
            continue
        j = sec(D, v, t)
        M = {round(b["tau"], 3): b for b in j["bantlar"] if b.get("olculdu")}
        print(f"{v:8s} " + "".join(
            f"{(M[x]['phi'] if x in M else float('nan')):+10.4f}" for x in tg))
    print("  (jackknife ±φ B2'de)")
    print(f"\n{'veri':8s} " + "".join(f"{('±'+str(x)):>10s}" for x in tg))
    for (v, t) in ciftler:
        if abs(t - 0.40) > 1e-9:
            continue
        j = sec(D, v, t)
        M = {round(b["tau"], 3): b for b in j["bantlar"] if b.get("olculdu")}
        print(f"{v:8s} " + "".join(
            f"{(M[x]['sPhi_jk'] if x in M else float('nan')):10.4f}"
            for x in tg))

    print()
    print("=" * 110)
    print("C) φ ŞEKLİ — DOĞRUSAL mı EĞRİLİKLİ mi?  (taban 0.40, (0.42,0.56])")
    print(f"{'veri':8s} {'n':>3s} | {'DOĞRUSAL':>24s} | "
          f"{'KUADRATİK':>32s} | {'ΔAIC':>6s} {'hüküm':>10s}")
    print(f"{'':8s} {'':>3s} | {'τ₀':>8s} {'a':>7s} {'χ²/dof':>7s} | "
          f"{'τ₀':>8s} {'a':>7s} {'b':>8s} {'χ²/dof':>6s} |")
    for (v, t) in ciftler:
        if abs(t - 0.40) > 1e-9:
            continue
        d1 = t0(HAV[(v, t)], 0.42, 0.56, 0.02, 1)
        d2 = t0(HAV[(v, t)], 0.42, 0.56, 0.02, 2)
        if not d1 or not d2:
            continue
        a1, a2 = d1["chi2"] + 4, d2["chi2"] + 6
        print(f"{v:8s} {d1['n']:3d} | {d1['tau0']:8.4f} {d1['egim']:7.3f} "
              f"{d1['chi2']/d1['dof']:7.2f} | {d2['tau0']:8.4f} "
              f"{d2['beta'][1]:7.3f} {d2['beta'][2]:8.2f} "
              f"{d2['chi2']/d2['dof']:6.2f} | {a2-a1:+6.2f} "
              f"{('EĞRİ' if a2 < a1-2 else 'ayrılamaz'):>10s}")

    print()
    print("=" * 110)
    print("D) KONVANSİYON SİSTEMATİĞİ (taban 0.40; ortak aralık (0.46,0.56])")
    print(f"{'veri':8s} {'w=.01':>8s} {'w=.02':>8s} {'w=.025':>8s} "
          f"{'w=.05/d1':>9s} | {'τ̄ apsis':>9s} {'Δapsis':>8s} | "
          f"{'deg1':>8s} {'deg2':>8s} {'deg3':>8s} | {'(.42,.56]':>10s} "
          f"{'(.48,.58]':>10s}")
    for (v, t) in ciftler:
        if abs(t - 0.40) > 1e-9:
            continue
        H = HAV[(v, t)]
        nan = float("nan")

        def g(d):
            return d["tau0"] if d else nan
        row = [g(t0(H, 0.46, 0.56, w, 2)) for w in (0.01, 0.02, 0.025)]
        d5 = t0(H, 0.46, 0.56, 0.05, 1)
        db = t0(H, 0.46, 0.56, 0.02, 2, "tbar")
        dd = [t0(H, 0.46, 0.56, 0.02, k) for k in (1, 2, 3)]
        w1 = t0(H, 0.42, 0.56, 0.02, 2)
        w2 = t0(H, 0.48, 0.58, 0.02, 2)
        print(f"{v:8s} " + "".join(f"{x:8.4f}" for x in row)
              + f"{(d5['tau0'] if d5 else nan):9.4f} | "
              f"{g(db):9.4f} {g(db)-row[1]:+8.4f} | "
              + "".join(f"{(x['tau0'] if x else float('nan')):8.4f}"
                        for x in dd)
              + f" | {(w1['tau0'] if w1 else float('nan')):10.4f} "
              f"{(w2['tau0'] if w2 else float('nan')):10.4f}")

    print()
    print("=" * 110)
    print("E) TABAN (REGRESYON KONVANSİYONU) TARAMASI — EN BÜYÜK SİSTEMATİK")
    print("   Aynı bantlar, aynı fit; DEĞİŞEN tek şey η regresyonundan")
    print("   çıkarılan çizgilerin üst sınırı (τ ≤ taban).")
    print(f"{'veri':8s} {'taban':>5s} {'σ_η²':>7s} {'c₁':>8s} {'ρ@.55':>7s} "
          f"{'τ₀ (0.46,0.56]':>15s} {'±jk':>7s} {'τ₀ (0.42,0.56]':>15s} "
          f"{'τ₀ (0.52,0.60] EKS':>19s}")
    TAB = {}
    for (v, t) in ciftler:
        j = sec(D, v, t)
        H = HAV[(v, t)]
        b55 = [b for b in j["bantlar"] if b.get("olculdu")
               and abs(b["tau"] - 0.55) < 1e-6]
        d46 = t0(H, 0.46, 0.56, 0.02, 2) if t <= 0.46 + 1e-9 else None
        d42 = t0(H, 0.42, 0.56, 0.02, 2) if t <= 0.42 + 1e-9 else None
        d52 = t0(H, 0.52, 0.60, 0.02, 2)
        print(f"{v:8s} {t:5.2f} {j['s_eta']:7.4f} {j['c1']:+8.5f} "
              f"{(b55[0]['rho'] if b55 else float('nan')):7.4f} "
              f"{(d46['tau0'] if d46 else float('nan')):15.4f} "
              f"{(d46['s_jk'] if d46 else float('nan')):7.4f} "
              f"{(d42['tau0'] if d42 else float('nan')):15.4f} "
              f"{(d52['tau0'] if d52 else float('nan')):19.4f}")
        TAB[(v, t)] = dict(seta=j["s_eta"], c1=j["c1"],
                           t46=d46["tau0"] if d46 else None,
                           s46=d46["s_jk"] if d46 else None,
                           t42=d42["tau0"] if d42 else None,
                           t52=d52["tau0"] if d52 else None)

    print("\n   TABAN EĞİLİMİ:  τ₀ = a + b·σ_η²  (aynı veride ≥3 taban varsa)")
    print(f"{'veri':8s} {'ntaban':>6s} {'b (dτ₀/dσ_η²)':>15s} "
          f"{'a = τ₀(σ_η²→0)':>16s} {'yayılım':>9s}")
    EKS = {}
    for v in SIRA:
        pts = [(TAB[(v, t)]["seta"], TAB[(v, t)]["t46"])
               for (vv, t) in ciftler if vv == v and TAB[(v, t)]["t46"]]
        if len(pts) < 3:
            continue
        x = np.array([p[0] for p in pts]); y = np.array([p[1] for p in pts])
        b, a = np.polyfit(x, y, 1)
        EKS[v] = (a, b)
        print(f"{v:8s} {len(pts):6d} {b:15.3f} {a:16.4f} "
              f"{y.max()-y.min():9.4f}")

    print()
    print("=" * 110)
    print("F) HÜKÜM TABLOSU — τ₀ (kuadratik, 0.02, τ_eff)")
    print(f"{'veri':8s} {'L':>7s} | {'taban 0.28':>10s} {'0.34':>8s} "
          f"{'0.40':>8s} {'0.46':>8s} | {'0.40 derin':>10s} "
          f"{'(σ_η²→0)':>9s}  açıklama")
    for v in SIRA:
        ts = [t for (vv, t) in ciftler if vv == v]
        if not ts:
            continue
        j = sec(D, v, ts[0])
        row = []
        for tt in (0.28, 0.34, 0.40, 0.46):
            k = (v, tt)
            row.append(TAB[k]["t46"] if k in TAB and TAB[k]["t46"] else
                       float("nan"))
        derin = (TAB[(v, 0.40)]["t42"] if (v, 0.40) in TAB
                 and TAB[(v, 0.40)]["t42"] else float("nan"))
        print(f"{v:8s} {j['L']:7.3f} | " + "".join(f"{x:{w}.4f}" for x, w in
                                                   zip(row, (10, 8, 8, 8)))
              + f" | {derin:10.4f} "
              f"{(EKS[v][0] if v in EKS else float('nan')):9.4f}  "
              f"{ETIKET.get(v,'')}")

    print()
    print("=" * 110)
    print("F2) 154 İLE MUTABAKAT — 0.5153 nereden 0.51'e indi?")
    print("    154'ün konvansiyonu: 0.03'lük bantlar (0.525'ten), apsis bant")
    print("    ORTASI (τ̄), en alttaki ÜÇ banda parabol, EKSTRAPOLASYON.")
    print(f"{'veri/taban':14s} {'154 usulü (τ̄,0.03,3 bant)':>27s} "
          f"{'aynı ama τ_eff':>16s} {'155 usulü':>11s} {'Δ(apsis)':>9s} "
          f"{'Δ(usul)':>9s}")
    k154 = [0.525, 0.555, 0.585, 0.615]
    for (v, t) in ciftler:
        if v not in ("son", "orta", "dusuk") or t not in (0.46, 0.52):
            continue
        H = HAV[(v, t)]
        a1 = A.tau0_jk(H, k154, "tbar", 2)
        a2 = A.tau0_jk(H, k154, "teff", 2)
        a3 = t0(H, max(0.46, t), 0.60 if t >= 0.5 else 0.56, 0.02, 2)
        nan = float("nan")
        print(f"{v+'/'+str(t):14s} {(a1['tau0'] if a1 else nan):27.4f} "
              f"{(a2['tau0'] if a2 else nan):16.4f} "
              f"{(a3['tau0'] if a3 else nan):11.4f} "
              f"{((a2['tau0']-a1['tau0']) if a1 and a2 else nan):+9.4f} "
              f"{((a3['tau0']-a2['tau0']) if a2 and a3 else nan):+9.4f}")
    print("    [154, taban 0.52 / son için 0.5157; taban 0.46 R-sıfırı için "
          "0.5153 bildirmişti]")

    print()
    print("=" * 110)
    print("G) H-L (τ₀ = ½ + c/L) vs SABİT — üç gerçek pencere, HER TABANDA")
    ozet = {"HL": {}}
    for tt in (0.28, 0.34, 0.40, 0.46, 0.52):
        dat = []
        for v in ("son", "orta", "dusuk"):
            k = (v, tt)
            if k not in TAB:
                continue
            val = TAB[k]["t46"] or TAB[k]["t52"]
            s = TAB[k]["s46"] or 0.001
            if val:
                dat.append((sec(D, v, tt)["L"], val, max(s, 1e-4), v))
        if len(dat) < 3:
            continue
        Lv = np.array([d[0] for d in dat]); tv = np.array([d[1] for d in dat])
        sv = np.array([d[2] for d in dat]); w = 1 / sv**2
        x = 1.0 / Lv
        c0 = float(np.sum(w * tv) / np.sum(w))
        chi_s = float(np.sum(w * (tv - c0)**2))
        c_hl = float(np.sum(w * x * (tv - 0.5)) / np.sum(w * x * x))
        chi_h = float(np.sum(w * (tv - 0.5 - c_hl * x)**2))
        chi18 = float(np.sum(w * (tv - 0.5 - 0.18 * x)**2))
        # H-L'nin ÖNGÖRDÜĞÜ fark (son ↔ dusuk) ile ÖLÇÜLEN fark
        i_s = [d[3] for d in dat].index("son")
        i_d = [d[3] for d in dat].index("dusuk")
        ong = 0.18 * (x[i_d] - x[i_s])
        olc = tv[i_d] - tv[i_s]
        sfark = float(np.hypot(sv[i_d], sv[i_s]))
        print(f"\n  --- taban {tt} ---")
        for (Lx, tx, sx, vv) in dat:
            print(f"    {vv:8s} L={Lx:7.4f}  τ₀={tx:.4f} ± {sx:.4f}   "
                  f"H-L(c=.18) öngörü {0.5+0.18/Lx:.4f}")
        print(f"    ÖLÇÜLEN  τ₀(dusuk)−τ₀(son) = {olc:+.4f} ± {sfark:.4f}"
              f"   |   H-L ÖNGÖRÜSÜ = {ong:+.4f}   "
              f"→ {abs(olc-ong)/max(sfark,1e-9):.1f}σ")
        print(f"    sabit χ²={chi_s:8.2f} | H-L(c serbest={c_hl:+.4f}) "
              f"χ²={chi_h:8.2f} | H-L(c=0.18, 0 par) χ²={chi18:8.2f}")
        ozet["HL"][str(tt)] = dict(L=Lv.tolist(), tau0=tv.tolist(),
                                   sig=sv.tolist(), c0=c0, chi_s=chi_s,
                                   c_hl=c_hl, chi_h=chi_h, chi18=chi18,
                                   olc=olc, ong=ong, sfark=sfark)

    print()
    print("=" * 110)
    print("H) H-N — δ = τ₀ − ½ hangi istatistikle? (taban 0.40, (0.42,0.56])")
    print(f"{'veri':8s} {'σΔ²':>8s} {'σ_ds²':>7s} {'κ₃Δ':>9s} {'çarp':>6s} "
          f"{'τ₀':>8s} {'δ':>8s} | {'γ=δ/σΔ²':>9s} {'β=δ/κ₃':>9s}  açıklama")
    MN = []
    for v in SIRA:
        k = (v, 0.40)
        if k not in TAB or TAB[k]["t42"] is None \
                or not np.isfinite(TAB[k]["t42"]):
            continue                      # P0: |Γ|<0.3, n_eff≈1 → ölçüm yok
        j = sec(D, v, 0.40); m = j["mom"]
        tt = TAB[k]["t42"]; dl = tt - 0.5
        print(f"{v:8s} {m['sA2']:8.5f} {j['s_ds']:7.4f} {m['k3']:+9.5f} "
              f"{m['carp']:+6.3f} {tt:8.4f} {dl:+8.4f} | "
              f"{dl/m['sA2']:9.3f} "
              f"{(dl/m['k3'] if abs(m['k3'])>1e-6 else float('nan')):9.1f}  "
              f"{ETIKET.get(v,'')}")
        if "S0." not in v:      # sahte-taban satırları gerçeğin KOPYASIDIR;
            MN.append(dict(v=v, sA2=m["sA2"], sds2=j["s_ds"], k3=m["k3"],
                           carp=m["carp"], tau0=tt, delta=dl,
                           seta=j["s_eta"], c1=j["c1"]))   # fite girmez
    ozet["HN"] = MN
    if len(MN) >= 4:
        d = np.array([x["delta"] for x in MN])
        print()
        for ad, g in (("σΔ²", np.array([x["sA2"] for x in MN])),
                      ("σ_ds²", np.array([x["sds2"] for x in MN])),
                      ("κ₃Δ", np.array([x["k3"] for x in MN])),
                      ("çarp γ₁", np.array([x["carp"] for x in MN]))):
            kk = float(np.sum(g * d) / np.sum(g * g))
            art = d - kk * g
            print(f"    δ = k·{ad:8s}: k={kk:9.3f}  rms artık="
                  f"{np.sqrt(np.mean(art**2)):.5f}  (δ'nın kendi rms'i="
                  f"{np.sqrt(np.mean(d**2)):.5f})")
        print("    → rms artık δ'nın rms'inden küçük DEĞİLSE o moment "
              "açıklamıyor.")

    print()
    print("=" * 110)
    print("I) R'NİN KENDİ SIFIRI (154 usulü) — φ sıfırıyla farkı")
    print(f"{'veri':8s} {'taban':>5s} {'τ_R0':>8s} {'±jk':>7s} "
          f"{'τ₀(φ)':>8s} {'fark':>8s} {'arg M_emp':>10s}")
    for (v, t) in ciftler:
        j = sec(D, v, t)
        d = A.tau0_R(j, 1)
        if not d:
            continue
        H = HAV[(v, t)]
        dp = t0(H, max(0.42, t + 0.02), 0.56, 0.02, 2)
        b0 = [x for x in j["bantlar"] if x.get("olculdu")]
        print(f"{v:8s} {t:5.2f} {d['tau0']:8.4f} {d['s_jk']:7.4f} "
              f"{(dp['tau0'] if dp else float('nan')):8.4f} "
              f"{(d['tau0']-dp['tau0'] if dp else float('nan')):+8.4f} "
              f"{b0[0]['argMe']:+10.4f}")

    ozet["taban"] = {f"{v}|{t}": TAB[(v, t)] for (v, t) in TAB}
    json.dump(ozet, open(SCR / "tau0_ozet.json", "w"), indent=1, default=float)
    print("\n-> tau0_ozet.json")


if __name__ == "__main__":
    main()
