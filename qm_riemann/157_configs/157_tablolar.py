"""
157 — RAPOR TABLOLARI (markdown). Her sayı JSON'lardan okunur; elle
girilen tek şey 154/155/156'nın referans değerleridir (karşılaştırma
sütunlarında, kaynağı belirtilerek).

Çıktı: scratchpad/157/tablolar.md
"""
import importlib
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
A = importlib.import_module("157_analiz")

M = []


def w(*a):
    M.append(" ".join(str(x) for x in a))


def T0(D):
    w("### T0 — kanal künyeleri (ayrışım ds = drift + lad + eta)\n")
    w("| pencere | taban | regresör | σ_ds² | σ_lad² | σ_η² | c₁(η) | σΔ² | "
      "pay Cov(X_lad,X_tam)/σΔ² | korel(X_lad,X_tam) | ölçek s_lad |")
    w("|---|---|---|---|---|---|---|---|---|---|---|")
    for v in A.VERI:
        for t in A.TABAN:
            k = (v, t, "kaba")
            if k not in D:
                continue
            d = D[k]
            b = A.bant(d)[0]
            w(f"| {v} | {t:.2f} | {d['nq']} | {d['s_ds']:.5f} | "
              f"{d['s_lad']:.5f} | {d['s_eta']:.5f} | {d['c1']:+.5f} | "
              f"{d['sA2']:.5f} | {d['kov']['lad']/d['sA2']:+.4f} | "
              f"{d['korel']['lad']:+.4f} | {d['olcek']['lad']:.4f} |")
    w("")
    w("σΔ² = Var(X_tam) ve X_tam TABANDAN BAĞIMSIZDIR (aynı sütun her")
    w("tabanda aynı) — raporun omurgası olan gözlem.\n")


def T1(D, taban=0.46):
    w(f"### T1 — R tayfı, taban {taban}, 0.03'lük ızgara, apsis τ_eff\n")
    w("| τ̄ | τ_eff | pencere | \\|Γ\\| | φ_Γ | arg M_emp | R_tam | ±jk | "
      "\\|ΔRe\\|_tam | R_lad | ±jk | \\|ΔRe\\|_lad | n_eff(lad) | ρ | not |")
    w("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    taus = sorted({b["tau"] for v in A.VERI if (v, taban, "kaba") in D
                   for b in A.bant(D[(v, taban, "kaba")])})
    for tau in taus:
        for v in A.VERI:
            if (v, taban, "kaba") not in D:
                continue
            bb = [b for b in A.bant(D[(v, taban, "kaba")])
                  if b["tau"] == tau]
            if not bb:
                continue
            b = bb[0]
            kt, kl = b["kanal"]["tam"], b["kanal"]["lad"]
            nt = []
            if b["tau"] > A.TAU_UST:
                nt.append("τ>0.76 fit-dışı")
            if kl["dRe"] >= 0.06:
                nt.append("ΔRe>0.06")
            if kl["artik"] >= 0.02:
                nt.append("faz tutmadı")
            w(f"| {tau:.4f} | {b['tau_eff']:.4f} | {v} | {b['absG']:.3f} | "
              f"{b['phi']:+.4f} | {b['argMe']:+.4f} | {kt['Rham']:+.3f} | "
              f"{kt['sRham_jk']:.3f} | {kt['dRe']:.4f} | "
              f"{kl['Rham']:+.3f} | {kl['sRham_jk']:.3f} | {kl['dRe']:.4f} | "
              f"{kl['neff']:.0f} | {b['rho']:.4f} | {', '.join(nt)} |")
    w("")


def T2(D, v="son"):
    w(f"### T2 — TABAN EKSENİ: R(taban) sabit τ'da (pencere {v})\n")
    taus = [x for x in sorted({b["tau"] for t in A.TABAN
                               if (v, t, "kaba") in D
                               for b in A.bant(D[(v, t, "kaba")])})
            if 0.58 < x <= A.TAU_UST]
    for kan in ("tam", "lad"):
        w(f"**R_{kan} (ham)**\n")
        w("| taban | σ_η² | σ_lad² | pay_lad | " +
          " | ".join(f"τ̄={x:.3f}" for x in taus) + " |")
        w("|---" * (4 + len(taus)) + "|")
        for t in A.TABAN:
            if (v, t, "kaba") not in D:
                continue
            d = D[(v, t, "kaba")]
            row = []
            for x in taus:
                bb = [b for b in A.bant(d) if b["tau"] == x]
                if not bb or bb[0]["kanal"][kan]["artik"] >= 0.02:
                    row.append("—")
                else:
                    f = "" if A.saglam(bb[0], kan) else " ⚠"
                    row.append(f"{bb[0]['kanal'][kan]['Rham']:+.3f}{f}")
            w(f"| {t:.2f} | {d['s_eta']:.4f} | {d['s_lad']:.4f} | "
              f"{d['kov']['lad']/d['sA2']:+.4f} | " + " | ".join(row) + " |")
        for etiket, tbs in (("156 penceresi {0.46,0.52,0.58}", A.TABAN_156),
                            ("TAM eksen {0.28…0.58}", A.TABAN)):
            row = []
            for x in taus:
                vv = []
                for t in tbs:
                    if (v, t, "kaba") not in D:
                        continue
                    bb = [b for b in A.bant(D[(v, t, "kaba")])
                          if b["tau"] == x]
                    if bb and A.saglam(bb[0], kan):
                        vv.append(bb[0]["kanal"][kan]["Rham"])
                row.append(f"**{100*(max(vv)-min(vv))/abs(np.mean(vv)):.1f}%**"
                           f" ({len(vv)})" if len(vv) >= 2 else "—")
            w(f"| _yayılım %_ | _{etiket}_ |  |  | " + " | ".join(row) + " |")
        w("")
        w("⚠ = n_eff < 3000 (156'nın kuralı); yayılım satırlarında bu "
          "hücreler SAYILMADI. Parantez: kaç taban.\n")


def T3(D, kayit):
    w("### T3 — τ₀ (apsis τ_eff, kuadratik fit, ±grup-sil jackknife)\n")
    w("| pencere | taban | σ_η² | τ₀(φ_Γ) | τ₀\\*=(φ−argM_emp) sıfırı | "
      "τ₀(R_tam) | τ₀(R_lad) | Δ(lad−tam) | 155'in τ₀(φ)'si |")
    w("|---|---|---|---|---|---|---|---|---|")
    R155 = {("son", 0.28): 0.4998, ("son", 0.34): 0.5044,
            ("son", 0.40): 0.5091, ("son", 0.46): 0.5122,
            ("son", 0.52): 0.5126, ("orta", 0.28): 0.4996,
            ("orta", 0.34): 0.5049, ("orta", 0.40): 0.5091,
            ("orta", 0.46): 0.5122, ("orta", 0.52): 0.5125}
    for v in A.VERI:
        for t in A.TABAN:
            if (v, t) not in kayit:
                continue
            k = kayit[(v, t)]
            w(f"| {v} | {t:.2f} | {k['seta']:.4f} | "
              f"{k['phi']:.4f} ± {k['ephi']:.4f} | "
              f"{k['dphi']:.4f} ± {k['edphi']:.4f} | "
              f"{k['tam']:.4f} ± {k['etam']:.4f} | "
              f"{k['lad']:.4f} ± {k['elad']:.4f} | {k['lad']-k['tam']:+.4f} | "
              f"{R155.get((v,t),'—')} |")
    w("")
    w("| pencere | nicelik | dτ₀/dσ_η² | kesişim (σ_η²→0) | "
      "yayılım (taban ekseni) |")
    w("|---|---|---|---|---|")
    for v in A.VERI:
        for ad in ("phi", "dphi", "tam", "lad"):
            xs = [kayit[(v, t)]["seta"] for t in A.TABAN if (v, t) in kayit
                  and np.isfinite(kayit[(v, t)][ad])]
            ys = [kayit[(v, t)][ad] for t in A.TABAN if (v, t) in kayit
                  and np.isfinite(kayit[(v, t)][ad])]
            if len(xs) < 3:
                continue
            p = np.polyfit(xs, ys, 1)
            w(f"| {v} | {ad} | {p[0]:+.3f} | {p[1]:.4f} | "
              f"{max(ys)-min(ys):.4f} |")
    w("\n155'in (gerçek, üç pencere) ölçtüğü: eğim −0.189/−0.189/−0.195, "
      "kesişim 0.5191/0.5191/0.5196.\n")


def T4(D, kayit, kanal="lad", taban=0.46, g="kaba", etiket=""):
    P = A.veri_kur(D, kanal, g, taban)
    if len(P) < 6:
        return
    t0m = float(np.nanmean([kayit[(v, taban)][kanal] for v in A.VERI
                            if (v, taban) in kayit
                            and np.isfinite(kayit[(v, taban)][kanal])]))
    w(f"### T4{etiket} — YARIŞ: R_{kanal}, taban {taban}, {g} ızgara, "
      f"n={len(P)} (iki pencere ortak)\n")
    tum = {}
    for mod in ("I", "II", "III", "IV"):
        sg = A.sigma(P, mod)
        AD = A.adaylar(P, sg, t0m)
        res = []
        for ad, (tanim, k, chi2, th0, mfn) in AD.items():
            c2, th = A._opt(chi2, th0)
            res.append(dict(ad=ad, tanim=tanim, k=k, chi2=c2,
                            aic=c2 + 2 * k, th=th, mfn=mfn, chi2f=chi2))
        res.sort(key=lambda r: r["aic"])
        tum[mod] = (res, float(np.mean(sg)))
    sira = [r["ad"] for r in tum["I"][0]]
    w("| aday | tanım | k | " + " | ".join(
        f"AIC [{m}] χ²/dof" for m in ("I", "II", "III", "IV")) + " |")
    w("|---|---|---|---|---|---|---|")
    for ad in sira:
        r0 = [r for r in tum["I"][0] if r["ad"] == ad][0]
        cells = []
        for m in ("I", "II", "III", "IV"):
            r = [q for q in tum[m][0] if q["ad"] == ad][0]
            cells.append(f"{r['aic']:.1f}  ({r['chi2']/max(len(P)-r['k'],1):.2f})")
        w(f"| **{ad}** | {r0['tanim']} | {r0['k']} | " + " | ".join(cells)
          + " |")
    w("")
    for m in ("I", "II", "III", "IV"):
        w(f"* hata modeli **{m}** = {A.HATA_MODU[m][0]}; ⟨σ⟩ = "
          f"{tum[m][1]:.4f}; kazanan: **{tum[m][0][0]['ad']}** "
          f"(ΔAIC = {tum[m][0][1]['aic']-tum[m][0][0]['aic']:.2f} ikinciye)")
    w("")
    w("**En iyi parametreler (hata modeli I, istatistik):**\n")
    w("| aday | parametreler |")
    w("|---|---|")
    for r in tum["I"][0][:6]:
        w(f"| {r['ad']} | " + ", ".join(f"{u:+.4f}" for u in r["th"]) + " |")
    w("")
    ps = [p for p in P if p["v"] == "son"]
    w("**Artık yapıları (R_ölç − model), pencere son:**\n")
    w("| aday | " + " | ".join(f"{p['x']:.3f}" for p in ps) + " |")
    w("|---" * (1 + len(ps)) + "|")
    y = np.array([p["R"] for p in P])
    for r in tum["I"][0][:6]:
        mm = r["mfn"](r["th"])
        vals = [f"{yy-m2:+.3f}" for yy, m2, p in zip(y, mm, P)
                if p["v"] == "son"]
        w(f"| {r['ad']} | " + " | ".join(vals) + " |")
    w("")
    return tum


def T5(D):
    w("### T5 — ilkel gözlenebilir φ_Γ(τ_eff): doğrusal mı, eğri mi?\n")
    w("| pencere | taban | bant | doğrusal τ₀ | a | χ²/dof | kuadratik τ₀ | "
      "a | b | χ²/dof | ΔAIC | hüküm |")
    w("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for v in A.VERI:
        for t in A.TABAN:
            k = (v, t, "ince") if (v, t, "ince") in D else (v, t, "tau0")
            if k not in D:
                continue
            bs = [b for b in A.bant(D[k]) if b["tau"] <= 0.61]
            if len(bs) < 5:
                continue
            xx = np.array([b["tau_eff"] for b in bs])
            yy = np.array([b["phi"] for b in bs])
            ee = np.array([b["sPhi_jk"] for b in bs])
            ee = np.where(np.isfinite(ee) & (ee > 0), ee, np.nanmedian(ee))
            out = []
            for deg in (1, 2):
                c = np.polyfit(xx, yy, deg, w=1.0 / ee)
                r = np.roots(c)
                r = np.array([z.real for z in r if abs(z.imag) < 1e-9])
                t0 = float(r[np.argmin(np.abs(r - xx.mean()))]) \
                    if len(r) else np.nan
                c2 = float(np.sum(((yy - np.polyval(c, xx)) / ee)**2))
                out.append((t0, c, c2))
            d1, d2 = out
            daic = (d2[2] + 6) - (d1[2] + 4)
            a_q = float(d2[1][1] + 2 * d2[1][0] * d2[0])
            w(f"| {v} | {t:.2f} | {len(xx)} | {d1[0]:.4f} | "
              f"{float(d1[1][0]):.3f} | {d1[2]/max(len(xx)-2,1):.2f} | "
              f"{d2[0]:.4f} | {a_q:.3f} | {float(d2[1][0]):+.3f} | "
              f"{d2[2]/max(len(xx)-3,1):.2f} | {daic:+.2f} | "
              f"{'**EĞRİ**' if daic < -2 else 'ayrılamaz'} |")
    w("\n155 (taban 0.40, τ∈(0.42,0.56]): son b = −6.47, orta −6.10, "
      "düşük-L −6.73.\n")


if __name__ == "__main__":
    D = A.yukle()
    kayit = {}
    for v in A.VERI:
        for t in A.TABAN:
            k = (v, t, "ince") if (v, t, "ince") in D else (v, t, "tau0")
            if k not in D:
                continue
            bs = A.bant(D[k])
            ust = 0.57 if len([b for b in bs if b["tau"] <= 0.57]) >= 5 \
                else 0.63
            a, ea, n = A.tau0_hatali(bs, "phi", ust)
            s, es, _ = A.tau0_hatali(bs, "dphi", ust)
            b_, eb, _ = A.tau0_hatali(bs, "tam", ust)
            c_, ec, _ = A.tau0_hatali(bs, "lad", ust)
            kayit[(v, t)] = dict(seta=D[k]["s_eta"], phi=a, ephi=ea,
                                 dphi=s, edphi=es, tam=b_, etam=eb,
                                 lad=c_, elad=ec, n=n, ust=ust)
    T0(D); T1(D, 0.46); T2(D, "son"); T2(D, "orta"); T3(D, kayit)
    T4(D, kayit, "lad", 0.46, "kaba", "a")
    T4(D, kayit, "tam", 0.46, "kaba", "b [KONTROL]")
    T4(D, kayit, "lad", 0.52, "kaba", "c [ROBUSTLUK: taban 0.52]")
    T4(D, kayit, "lad", 0.46, "ince", "d [ROBUSTLUK: 0.02 ızgara]")
    T5(D)
    p = A.SCR / "tablolar.md"
    p.write_text("\n".join(M))
    print("->", p, f"({len(M)} satır)")
