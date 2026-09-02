"""
164 — ANALİZ: eski-inşa ↔ sadakatli-inşa ↔ gerçek HÜKÜM TABLOLARI
=================================================================
Hiçbir yeni ölçüm yapmaz. Okuduğu ham kayıtlar:

  scratchpad/158/F_<veri>_t<taban>.json    ← 158'in koşuları (eski gazlar, gerçek)
  scratchpad/164/F_<veri>_t<taban>.json    ← 164'ün koşuları (sadakatli gazlar)
  scratchpad/164/faz152_<veri>.json        ← 152 konvansiyonunda faz oranı
  scratchpad/164/sadakat_<veri>.json       ← genlik sadakat denetimi
  scratchpad/164/insa_<veri>.json          ← inşa tanıları

Fit çekirdeği 161_cekirdek'ten IMPORT edilir (kopya değil): aynı
`np.polyfit(x, y, derece, w=1/σ)`, aynı pencere seçimi, aynı jackknife,
aynı ½ çapası.

Kullanım:  164_analiz.py [> çıktı.txt]
"""
import importlib
import json
import sys
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
sys.path.insert(0, str(QM / "161_configs"))
C161 = importlib.import_module("161_cekirdek")

SCR = C161.SCR
SCR158, SCR164 = SCR / "158", SCR / "164"
TWO_PI, FOUR_PI = C161.TWO_PI, C161.FOUR_PI

ESKI = ("keskin", "A4")
YENI = ("Skeskin", "SA4", "SA1", "NKkeskin", "Hkeskin", "HA4")
GERCEK = ("son", "orta")
TABAN = (0.28, 0.34, 0.40, 0.46, 0.52)
# 161 §3b'nin konvansiyon bulutu: 4 taban × 4 pencere × 2 ağırlık = 32 fit.
# (taban 0.52'de fit penceresinin 0.43–0.61'i büyük ölçüde taban altında
#  kalıyor; 161 onu bulutun DIŞINDA tutmuştu — burada da aynısı.)
TABAN_KONV = (0.28, 0.34, 0.40, 0.46)
PENCERELER = ("A", "B", "C", "D")
AGIRLIKLAR = ("phi", "delta")


def _yukle_bir(kok, v, t):
    for isim in (f"F_{v}_t{t}.json", f"F_{v}_t{t:g}.json"):
        p = kok / isim
        if p.exists():
            return json.load(open(p))
    return None


def topla():
    D = {}
    for v in ESKI + GERCEK:
        for t in TABAN:
            d = _yukle_bir(SCR158, v, t)
            if d:
                D[(v, t)] = d
    for v in YENI:
        for t in TABAN:
            d = _yukle_bir(SCR164, v, t)
            if d:
                D[(v, t)] = d
    return D


def parmak(D, v, t, pencere="A", agirlik="phi"):
    d = D.get((v, t))
    if d is None:
        return None
    bs = C161.bant(d)
    if len(bs) < 5:
        return None
    P = C161.delta_bant(bs)
    return C161.olc(P, pencere=pencere, agirlik=agirlik, capa=0.5)


def satir(D, v):
    """Birincil satır (taban 0.40, W-A, 1/σ_φ) + konvansiyon ortalaması."""
    o = parmak(D, v, 0.40)
    kv = {k: [] for k in ("d_half", "dd", "b", "a_t0", "t0_yildiz")}
    nfit = 0
    for t in TABAN_KONV:
        for pn in PENCERELER:
            for ag in AGIRLIKLAR:
                r = parmak(D, v, t, pn, ag)
                if r is None:
                    continue
                nfit += 1
                for k in kv:
                    kv[k].append(r[k])
    ozet = {k: C161.ozet(vv) for k, vv in kv.items()}
    return o, ozet, nfit


def tablo_parmak(D, gazlar):
    print("\n### T-A. PARMAK İZİ — birincil satır (taban 0.40, W-A, w=1/σ_φ)")
    print("| gaz | δ(½) | ±jk | dδ/dτ\\|½ | ±jk | b | ±jk | a_½ | a_τ₀ | τ₀* | n_bant |")
    print("|---|---|---|---|---|---|---|---|---|---|---|")
    for v in gazlar:
        o = parmak(D, v, 0.40)
        if o is None:
            print(f"| {v} | — | | | | | | | | | |")
            continue
        print(f"| {v} | {o['d_half']:+.4f} | {o['e_d_half']:.4f} | "
              f"{o['dd']:+.3f} | {o['e_dd']:.3f} | {o['b']:+.3f} | "
              f"{o['e_b']:.3f} | {o['a_half']:.3f} | {o['a_t0']:.3f} | "
              f"{o['t0_yildiz']:.4f} | {o['n']} |")
    print("\n### T-B. PARMAK İZİ — konvansiyon ortalaması "
          "(taban × pencere × ağırlık)")
    print("| gaz | n_fit | δ(½) ± σ_k | dδ/dτ\\|½ ± σ_k | b ± σ_k | "
          "a_τ₀ ± σ_k | τ₀* ± σ_k |")
    print("|---|---|---|---|---|---|---|")
    for v in gazlar:
        o, oz, nf = satir(D, v)
        if nf == 0:
            print(f"| {v} | 0 | — | — | — | — | — |")
            continue
        f = lambda k, p=4: (f"{oz[k][0]:+.{p}f} ± {oz[k][1]:.{p}f}")
        print(f"| {v} | {nf} | {f('d_half')} | {f('dd',3)} | {f('b',3)} | "
              f"{f('a_t0',3)} | {f('t0_yildiz')} |")


def tablo_delta_bant(D, gazlar, taban=0.40):
    print(f"\n### T-C. δ BANT BANT (taban {taban}); "
          f"δ = φ_Γ − (4π·τ_eff − 2π)")
    kol = {}
    for v in gazlar:
        d = D.get((v, taban))
        if d is None:
            continue
        bs = C161.bant(d)
        P = C161.delta_bant(bs)
        kol[v] = {round(float(t), 4): (float(x), float(dd), float(ss),
                                       float(ph), float(ag))
                  for t, x, dd, ss, ph, ag
                  in zip(P["tau"], P["x"], P["d"], P["s_d"], P["phi"],
                         P["absG"])}
    tlar = sorted({t for c in kol.values() for t in c})
    bas = "| τ̄ | " + " | ".join(f"{v}: δ (τ_eff, \\|Γ\\|)" for v in kol) + " |"
    print(bas)
    print("|" + "---|" * (len(kol) + 1))
    for t in tlar:
        hu = []
        for v in kol:
            c = kol[v].get(t)
            hu.append("—" if c is None
                      else f"{c[1]:+.4f} ({c[0]:.4f}, {c[4]:.2f})")
        print(f"| {t:.2f} | " + " | ".join(hu) + " |")
    return kol


def tablo_phi_oran(gazlar, ref="son"):
    """152'nin özgün faz oranı: arg Γ(sentetik)/arg Γ(gerçek), bant bant."""
    D = {}
    for v in gazlar + (ref,):
        p = SCR164 / f"faz152_{v}.json"
        if p.exists():
            D[v] = {round(b["tau"], 4): b for b in json.load(open(p))["bantlar"]
                    if b.get("olculdu")}
    if ref not in D:
        print("\n(152-konvansiyon faz koşusu yok)")
        return {}
    print("\n### T-D. ×1.4 FAZ ORANI — 152'nin ÖZGÜN ÖLÇÜSÜ "
          "(taban 0.52, cap 720)")
    tl = sorted(D[ref])
    print("| gaz | " + " | ".join(f"τ={t:.4f}" for t in tl) +
          " | ort (b1–4) |")
    print("|" + "---|" * (len(tl) + 2))
    for v in (ref,) + gazlar:
        if v not in D:
            continue
        hu = []
        o4 = []
        for i, t in enumerate(tl):
            b = D[v].get(t)
            r = D[ref].get(t)
            if b is None or r is None or abs(r["phi"]) < 1e-12:
                hu.append("—")
                continue
            o = b["phi"] / r["phi"]
            pat = abs(b["absG"]) > 1.5
            hu.append(f"{o:.2f}" + ("‡" if pat else ""))
            if i < 4 and not pat:
                o4.append(o)
        hu.append(f"**{np.mean(o4):.2f}**" if o4 else "—")
        print(f"| {v} | " + " | ".join(hu) + " |")
    print("\n(‡ = |Γ| > 1.5, 152'nin patlak ölçütü; ortalamaya girmiyor.)")
    print("\n| gaz | " + " | ".join(f"φ(τ={t:.4f})" for t in tl) + " |")
    print("|" + "---|" * (len(tl) + 1))
    for v in (ref,) + gazlar:
        if v not in D:
            continue
        print(f"| {v} | " + " | ".join(
            f"{D[v][t]['phi']:+.4f}" if t in D[v] else "—" for t in tl) + " |")
    print("\n| gaz | " + " | ".join(f"\\|Γ\\|(τ={t:.4f})" for t in tl) + " |")
    print("|" + "---|" * (len(tl) + 1))
    for v in (ref,) + gazlar:
        if v not in D:
            continue
        print(f"| {v} | " + " | ".join(
            f"{D[v][t]['absG']:.3f}" if t in D[v] else "—" for t in tl) + " |")
    print("\n| gaz | " + " | ".join(f"R(τ={t:.4f})" for t in tl) +
          " | " + " | ".join(f"n_eff(τ={t:.4f})" for t in tl) + " |")
    print("|" + "---|" * (2 * len(tl) + 1))
    for v in (ref,) + gazlar:
        if v not in D:
            continue
        print(f"| {v} | " + " | ".join(
            f"{D[v][t]['R']:+.2f}" if t in D[v] else "—" for t in tl) +
            " | " + " | ".join(
            f"{D[v][t]['neff']:.0f}" if t in D[v] else "—" for t in tl) + " |")
    return D


def tablo_S3(D, gazlar):
    print("\n### T-E. S3 istatistikleri (taban 0.40 koşusundan) ve n_eff rejimi")
    print("| gaz | L | σ_ds² | σ_η² | c₁ | σΔ² | n_eff (τ̄=0.55) | "
          "n_eff (τ̄=0.61) | ⟨n_eff⟩ τ̄≥0.5 |")
    print("|---|---|---|---|---|---|---|---|---|")
    for v in gazlar:
        d = D.get((v, 0.40))
        if d is None:
            continue
        bs = {round(b["tau"], 4): b for b in C161.bant(d)}
        ne = [b["neff"] for t, b in bs.items() if t >= 0.5]
        print(f"| {v} | {d['L']:.4f} | {d['s_ds']:.4f} | {d['s_eta']:.4f} | "
              f"{d['c1']:+.5f} | {d['sA2']:.4f} | "
              f"{bs.get(0.55,{}).get('neff',float('nan')):.0f} | "
              f"{bs.get(0.61,{}).get('neff',float('nan')):.0f} | "
              f"{np.mean(ne) if ne else float('nan'):.0f} |")
    print("\n[gerçek hedef: σ_ds² = 0.1674, σ_η² = 0.0227, c₁ = −0.01158]")


def tablo_sadakat(gazlar):
    print("\n### T-F. GENLİK SADAKATİ — gerçekleşen/nominal (bant bant)")
    D = {}
    for v in gazlar:
        p = SCR164 / f"sadakat_{v}.json"
        if p.exists():
            D[v] = json.load(open(p))
    if not D:
        print("(sadakat koşusu yok)")
        return {}
    tl = sorted({(b["lo"], b["hi"]) for d in D.values() for b in d["bant"]})
    print("| τ-bant | " + " | ".join(D) + " |")
    print("|" + "---|" * (len(D) + 1))
    for lo, hi in tl:
        hu = []
        for v in D:
            b = next((x for x in D[v]["bant"]
                      if x["lo"] == lo and x["hi"] == hi), None)
            hu.append("—" if b is None else
                      f"{b['R']:.3f}" + (f" (SNR {b['snr']:.1f})"
                                         if b["snr"] < 3 else ""))
        print(f"| {lo:.2f}–{hi:.2f} | " + " | ".join(hu) + " |")
    print("\n| çizgi (163 §1c) | " + " | ".join(D) + " |")
    print("|" + "---|" * (len(D) + 1))
    for q in (2, 3, 5, 7, 11, 101, 1009):
        hu = []
        for v in D:
            r = D[v]["cizgi"].get(str(q))
            if r is None:
                hu.append("—")
                continue
            c = complex(*r["con"])
            hu.append(f"{abs(c)/r['b_pen']:.3f}")
        print(f"| q={q} | " + " | ".join(hu) + " |")
    return D


def main():
    D = topla()
    print("# 164 — ANALİZ ÇIKTISI (ham)")
    print(f"\nyüklenen koşular: {sorted({k[0] for k in D})}")
    for v in sorted({k[0] for k in D}):
        print(f"  {v}: taban {sorted(t for (vv, t) in D if vv == v)}")
    gazlar = tuple(v for v in ("son", "keskin", "NKkeskin", "Skeskin", "Hkeskin",
                    "A4", "SA4", "HA4", "SA1")
                   if any(k[0] == v for k in D))
    tablo_parmak(D, gazlar)
    tablo_delta_bant(D, gazlar)
    tablo_S3(D, gazlar)
    sent = tuple(v for v in ("keskin", "NKkeskin", "Skeskin", "Hkeskin",
                             "A4", "SA4", "HA4", "SA1")
                 if (SCR164 / f"faz152_{v}.json").exists())
    tablo_phi_oran(sent)
    tablo_sadakat(tuple(v for v in
                        ("son", "keskin", "NKkeskin", "Skeskin", "Hkeskin",
                         "A4", "SA4", "HA4", "SA1", "eski_keskin",
                         "eski_A4")
                        if (SCR164 / f"sadakat_{v}.json").exists()))
    print("\n### T-G. İNŞA TANILARI")
    print("| gaz | tip | h | ızgara | ΔG<0 kesri | maks\\|F\\| | "
          "\\|F\\|>1e−8 tekne | sıralı | min Δz | σ_ds² |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    for p in sorted(SCR164.glob("insa_*.json")):
        j = json.load(open(p))
        print(f"| {j['ad']} | {j['tip']} | {j.get('h','—')} | "
              f"{j.get('ng','—')} | {j.get('fr_dG_neg',float('nan')):.4f} | "
              f"{j['maxF']:.2e} | {j.get('nF_asan','—')} | "
              f"{'TAM' if j['sirali'] else 'BOZUK'} | "
              f"{j['min_dz']:.5f} | {j['sigma_ds2']:.4f} |")


if __name__ == "__main__":
    main()
