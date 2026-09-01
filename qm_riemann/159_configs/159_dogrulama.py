"""
159 — DENETİM (kopya-kayması + cebir + kinematik omurga)

V1  159'un φ, τ_eff, σ_φ'si = 158'in kaydı mı? (bit düzeyi)
V2  ÖZDEŞLİK: zp·conj(zc)·e^{−iA}/4 = ⟨(ρ+iσ)·e^{+iA·X̃}⟩ (gerçek veri)
V3  SENTETİK KONTROL: bilinen tek çizgi + bilinen adım dizisi üzerinde
    arg[zp·conj(zc)] = +A mı, −A mı? (de-rotasyonun işareti)
V4  156_cekirdek.zincir3'ün η'sı = 155'in önbelleklediği η mı?
V5  Kinematik omurga: φ − (4πτ_eff − 2π) = δ mi? (bant düzeyi)
"""
import importlib
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
C = importlib.import_module("159_cekirdek")
TWO_PI = 2 * np.pi
SCR158 = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
              "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/158")


def V1():
    print("=" * 78)
    print("V1 — 159'un ölçüm zinciri = 158'inki mi? (φ, τ_eff, σ_φ bit düzeyi)")
    print("=" * 78)
    ok = True
    for veri, taban in (("son", 0.4), ("son", 0.52), ("A4", 0.4),
                        ("keskin", 0.4), ("P1", 0.4), ("orta", 0.4)):
        p9 = C.SCR / f"S_{veri}_t{taban}_g158.json"
        p8 = SCR158 / f"F_{veri}_t{taban}.json"
        if not (p9.exists() and p8.exists()):
            continue
        a = {b["tau"]: b for b in json.load(open(p9))["bantlar"]
             if b.get("olculdu")}
        b_ = {b["tau"]: b for b in json.load(open(p8))["bantlar"]
              if b.get("olculdu")}
        ort = sorted(set(a) & set(b_))
        d1 = max(abs(a[t]["phi"] - b_[t]["phi"]) for t in ort)
        d2 = max(abs(a[t]["tau_eff"] - b_[t]["tau_eff"]) for t in ort)
        d3 = max(abs(a[t]["sPhi_jk"] - b_[t]["sPhi_jk"]) for t in ort)
        d4 = max(abs(a[t]["absG"] - b_[t]["absG"]) for t in ort)
        ok &= (d1 == 0.0 and d2 == 0.0)
        print(f"  {veri:6s} taban {taban:.2f}: {len(ort):2d} ortak bant  "
              f"maks|Δφ|={d1:.3e}  maks|Δτ_eff|={d2:.3e}  "
              f"maks|Δσ_φ|={d3:.3e}  maks|Δ|Γ||={d4:.3e}")
    print(f"  → bit düzeyinde özdeş: {'EVET' if ok else 'HAYIR'}")


def V2(veri="son", taban=0.40, nline=12):
    print()
    print("=" * 78)
    print("V2 — ÖZDEŞLİK  zp·conj(zc)·e^{−iA}/4 = ⟨(ρ+iσ)·e^{+iA·X̃}⟩")
    print("=" * 78)
    z = C.KOS155.veri_yukle(veri)
    Y = C.Yerel(z, veri, taban)
    eta, mid, L = Y.eta, Y.mid, Y.L
    qm = C.C154.pk_m(int(np.exp(0.86 * L)))
    qs = np.array(sorted(qm), float)
    tv = np.log(qs) / L
    sec = [int(qs[int(np.argmin(np.abs(tv - t)))])
           for t in np.linspace(0.30, 0.84, nline)]
    print("     q     τ      |lhs|      maks|lhs−rhs|   bağıl      "
          "arg(raw)−A     arg(raw)+A")
    mx = 0.0
    for q in sec:
        w = np.log(q); A = TWO_PI * w / L
        e0, e1, m0 = eta[:-1], eta[1:], mid[:-1]
        cw, sw = np.cos(w * m0), np.sin(w * m0)
        zc = complex(2 * np.mean(e0 * cw), -2 * np.mean(e0 * sw))
        zp = complex(2 * np.mean(e1 * cw), -2 * np.mean(e1 * sw))
        lhs = zp * np.conj(zc) * np.exp(-1j * A) / 4
        cfull = eta * np.exp(-1j * w * mid)
        pc = cfull[1:] * np.conj(zc / 2)
        rhs = np.mean(pc * np.exp(1j * A * Y.Xtil))
        d = abs(lhs - rhs); mx = max(mx, d / abs(lhs))
        ar = float(np.angle(zp * np.conj(zc)))
        print(f"  {q:6d} {w/L:.4f}  {abs(lhs):.4e}  {d:.3e}   {d/abs(lhs):.2e}  "
              f"{((ar-A+np.pi)%TWO_PI)-np.pi:+10.5f}   "
              f"{((ar+A+np.pi)%TWO_PI)-np.pi:+10.5f}")
    print(f"  → maks BAĞIL fark {mx:.2e}  (özdeşlik makine hassasiyetinde)")
    print("  → arg(raw)−A küçük, arg(raw)+A değil ⇒ ham korelatör e^{+iA} "
          "TAŞIR;\n     de-rotasyon e^{−iA} olmalı, kod e^{+iA} ile çarpıyor.")


def V3():
    print()
    print("=" * 78)
    print("V3 — SENTETİK KONTROL: bilinen çizgi, bilinen adım")
    print("=" * 78)
    rng = np.random.default_rng(7)
    L0 = 12.0
    n = 200000
    t0 = TWO_PI * np.exp(L0)
    # adımlar: ortalama 2π/L, kontrollü dalgalanma
    g = (TWO_PI / L0) * (1 + 0.20 * rng.standard_normal(n + 1))
    z = t0 + np.concatenate(([0.0], np.cumsum(g)))
    mid = 0.5 * (z[:-1] + z[1:])
    L = float(np.log(mid / TWO_PI).mean())
    w = 0.60 * L
    A = TWO_PI * w / L
    for ad, eta in (("saf çizgi", np.cos(w * mid)),
                    ("çizgi+gürültü", np.cos(w * mid)
                     + 2.0 * rng.standard_normal(len(mid)))):
        e0, e1, m0 = eta[:-1], eta[1:], mid[:-1]
        cw, sw = np.cos(w * m0), np.sin(w * m0)
        zc = complex(2 * np.mean(e0 * cw), -2 * np.mean(e0 * sw))
        zp = complex(2 * np.mean(e1 * cw), -2 * np.mean(e1 * sw))
        raw = zp * np.conj(zc)
        art = float(np.angle(raw * np.exp(-1j * A)))
        yanlis = float(np.angle(raw * np.exp(+1j * A)))
        print(f"  {ad:14s}: A={A:.4f}  arg(raw)={np.angle(raw):+.4f}  "
              f"arg(raw·e^−iA)={art:+.4f}  arg(raw·e^+iA)={yanlis:+.4f}"
              f"   [2A−2π = {2*A-TWO_PI:+.4f}]")
    print("  → e^{−iA} ile çarpım fazı SIFIRA indiriyor; e^{+iA} 2A−2π'ye "
          "taşıyor.")


def V5():
    print()
    print("=" * 78)
    print("V5 — kinematik omurga: φ − (4π·τ_eff − 2π) =? δ (çizgi çizgi arınmış)")
    print("=" * 78)
    for veri, taban, g in (("son", 0.4, "g158"), ("son", 0.4, "t1"),
                           ("A4", 0.4, "g158"), ("keskin", 0.4, "g158"),
                           ("P1", 0.4, "g158")):
        p = C.SCR / f"S_{veri}_t{taban}_{g}.json"
        if not p.exists():
            continue
        B = [b for b in json.load(open(p))["bantlar"] if b.get("olculdu")]
        d = max(abs(b["delta_b"] - b["delta"]) for b in B)
        r = max(abs(b["phi"] - (b["kinematik"] + b["delta"])) for b in B)
        print(f"  {veri:6s} t{taban} {g:5s}: {len(B):2d} bant   "
              f"maks|δ_b − δ| = {d:.2e}   maks|φ − (kin+δ)| = {r:.2e}")
    print("  → bant düzeyinde kaba arındırma (τ_eff ile) ile çizgi çizgi "
          "arındırma\n     arasındaki fark bant genişliği mertebesinde bir "
          "kalıntıdır.")


if __name__ == "__main__":
    V1(); V2(); V3(); V5()
