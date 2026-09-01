"""
160 — σ'NIN KİMLİĞİ: derin tanı (çizgi bazında)
==============================================
Üç aday sınanır:
  (a) yerel faz-gradyanı            → corr(σ, dX̃/dn)
  (b) çizginin yerel frekans kayması → SON BONDUN kendi adımı fazörü
      döndürüyor: u_n = u⁰_n·e^{−iA·X̃_n}
  (c) komşu-çizgi kaçağı            → etki komşuluk boşluğuna (gap) ve
      çizgi gücüne asılı olmalı; A'ya değil

TAM KİMLİK (yaklaşımsız):
    u_n ≡ ρ_n + iσ_n = c_{n+1}conj⟨c⟩ = η_{n+1}·|⟨c⟩|·e^{iθ_n}
    θ_n = −(W·m_{n+1} + arg⟨c⟩) = sabit − A(n+1) − A·C_n,  C_n = Σ_{k≤n} X̃_k
⇒ (ρ+iσ)·e^{+iA(n+1)} = η_{n+1}·e^{−iA·C_n}·(e^{−iW m_0}conj⟨c⟩)     [D1]

Yani σ, çizginin küresel taşıyıcıya göre BİRİKTİRDİĞİ fazın kuadratür
izidir. SON adımın payı ayrılabilir: u⁰_n ≡ u_n·e^{+iAX̃_n} (X̃_n'siz fazör)
    ρ = ρ⁰·cos(AX̃) + σ⁰·sin(AX̃)
    σ = σ⁰·cos(AX̃) − ρ⁰·sin(AX̃)                                     [D2]
σ⁰ ve X̃_n bağımsızsa Cov(σ⁰cos, X̃) ≈ 0 ve BÜTÜN σ–adım bağlaşımı
−ρ⁰sin(AX̃) teriminden gelir ⇒ kapalı form
    Cov(σ, X̃)/⟨ρ⟩ ≈ −A·σΔ²                                          [D3]
(Gauss adım için ⟨sin(AX)X⟩ = A s² e^{−A²s²/2} ve ⟨ρ⟩ ≈ ⟨ρ⁰⟩e^{−A²s²/2}.)
"""
import importlib
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
C = importlib.import_module("160_cekirdek")
TWO_PI = 2 * np.pi


def kor(a, b):
    return C._kor(a, b)


def cizgiler(Y, L, lo, hi, nmax=40, tohum=21):
    """159/160 ile AYNI aday seçimi ve gap filtresi (ilk nmax çizgi)."""
    qm = C.C154.pk_m(int(np.exp(0.86 * L)))
    allq = sorted(qm)
    allw = np.array([np.log(q) for q in allq])
    dres = TWO_PI / (Y.mid[-1] - Y.mid[0])
    rng = np.random.default_rng(tohum)
    tumu = [(q, w) for q, w in zip(allq, allw) if lo < w / L <= hi]
    cand = tumu
    if len(cand) > 220:
        idx = rng.choice(len(cand), 220, replace=False)
        cand = [cand[i] for i in idx]
    out = []
    for q, w in cand:
        j = np.searchsorted(allw, w)
        koms = [allw[k] for k in (j - 1, j + 1)
                if 0 <= k < len(allw) and abs(allw[k] - w) > 1e-12]
        gap = min(abs(w - k) for k in koms)
        if gap < 2.5 * dres:
            continue
        out.append((q, w, gap))
        if len(out) >= nmax:
            break
    return out, dres


def tani(veri, taban, lo, hi, nmax=40):
    z = C.KOS155.veri_yukle(veri)
    Y = C.Yerel160(z, veri, taban, 4000)
    L = Y.L
    Ls, dres = cizgiler(Y, L, lo, hi, nmax)
    N = len(Y.m0)
    n1 = np.arange(1, N + 1, dtype=float)
    Xt = Y.Xtil                      # ham X̃ (merkezlenmemiş)
    X0 = Y.Xtil0
    Cn = np.cumsum(Xt)               # C_n = Σ_{k≤n} X̃_k
    sX2 = float(np.var(X0))
    print("=" * 108)
    print(f"σ TANISI — {veri}, taban {taban}, bant ({lo},{hi}], "
          f"{len(Ls)} çizgi   [σΔ² = {sX2:.5f}]")
    print("=" * 108)
    print("  q        τ      A      gap/dres  |  D1 kimlik   | Cov(σ,X̃)/⟨ρ⟩  "
          "  −A·σΔ²   oran | Cov(σ⁰cos,X̃)/⟨ρ⟩  Cov(−ρ⁰sin,X̃)/⟨ρ⟩ | "
          "kor(σ,ρ) kor(σ,Δρ) kor(σ,X̃) kor(σ,dX̃/dn)")
    R = []
    for q, w, gap in Ls:
        Wf = w
        A = TWO_PI * Wf / L
        cw, sw = np.cos(Wf * Y.m0), np.sin(Wf * Y.m0)
        cw1, sw1 = np.cos(Wf * Y.mid[1:]), np.sin(Wf * Y.mid[1:])
        zc = complex(2 * np.mean(Y.e0 * cw), -2 * np.mean(Y.e0 * sw))
        zb = zc / 2.0
        c1r, c1i = Y.e1 * cw1, -Y.e1 * sw1
        rho = c1r * zb.real + c1i * zb.imag
        sig = c1i * zb.real - c1r * zb.imag
        rm = float(rho.mean())
        # --- D1: demodüle fazör kimliği ---------------------------
        ph = np.exp(1j * A * n1)                     # e^{+iA(n+1)}
        v = (rho + 1j * sig) * ph
        K = np.exp(-1j * Wf * Y.mid[0]) * np.conj(zb)
        v_pred = Y.e1 * np.exp(-1j * A * Cn) * K
        d1 = float(np.max(np.abs(v - v_pred)) / np.mean(np.abs(v)))
        # --- D2: son adımın ayrılması ------------------------------
        dr = cw1 * cw + sw1 * sw
        di = sw1 * cw - cw1 * sw
        ca, sa = float(np.cos(A)), float(np.sin(A))
        Cx = dr * ca + di * sa                       # cos(A X̃)
        Sx = di * ca - dr * sa                       # sin(A X̃)
        r0 = rho * Cx - sig * Sx                     # ρ⁰
        s0 = sig * Cx + rho * Sx                     # σ⁰
        def cov(a):
            return float((np.dot(a, X0) / len(a) - a.mean() * X0.mean()) / rm)
        c_all = cov(sig)
        c_old = cov(s0 * Cx)
        c_new = cov(-r0 * Sx)
        dro = np.concatenate((np.diff(rho), np.zeros(1)))
        R.append(dict(q=int(q), tau=w / L, A=A, gap=gap / dres, pow=abs(zc)**2,
                      d1=d1, c_all=c_all, c_old=c_old, c_new=c_new,
                      pred=-A * sX2,
                      k_sr=kor(sig, rho), k_sdr=kor(sig, dro),
                      k_sx=kor(sig, X0), k_sdx=kor(sig, Y.dX),
                      k_srx=kor(sig, rho * X0)))
        r = R[-1]
        print(f"  {q:7d} {r['tau']:.4f} {A:.3f} {r['gap']:9.2f}  | "
              f"{d1:.2e} | {c_all:+13.5f} {r['pred']:+10.5f} "
              f"{c_all/r['pred']:6.3f} | {c_old:+16.5f} {c_new:+18.5f} | "
              f"{r['k_sr']:+8.4f} {r['k_sdr']:+9.4f} {r['k_sx']:+9.5f} "
              f"{r['k_sdx']:+12.5f}")
    if not R:
        return R
    a = lambda k: np.array([x[k] for x in R], float)
    print(f"\n  ÖZET ({len(R)} çizgi):")
    print(f"    D1 kimlik (bağıl) maks              = {a('d1').max():.2e}")
    print(f"    Cov(σ,X̃)/⟨ρ⟩ / (−A·σΔ²)  ort={np.mean(a('c_all')/a('pred')):+.4f} "
          f" sd={np.std(a('c_all')/a('pred')):.4f}  "
          f"menzil=[{np.min(a('c_all')/a('pred')):+.3f},"
          f"{np.max(a('c_all')/a('pred')):+.3f}]")
    pay = np.abs(a('c_old')) / (np.abs(a('c_old')) + np.abs(a('c_new')))
    print(f"    'eski kuadratür' payı |Cov(σ⁰cos,X̃)|/(|·|+|Cov(−ρ⁰sin,X̃)|)"
          f"  ort={pay.mean():.4f}  maks={pay.max():.4f}")
    print(f"    korel(σ,ρ)      ort={a('k_sr').mean():+.5f}  "
          f"maks|·|={np.abs(a('k_sr')).max():.5f}")
    print(f"    korel(σ,Δρ)     ort={a('k_sdr').mean():+.5f}  "
          f"maks|·|={np.abs(a('k_sdr')).max():.5f}")
    print(f"    korel(σ,X̃)      ort={a('k_sx').mean():+.5f}  "
          f"maks|·|={np.abs(a('k_sx')).max():.5f}")
    print(f"    korel(σ,dX̃/dn)  ort={a('k_sdx').mean():+.5f}  "
          f"maks|·|={np.abs(a('k_sdx')).max():.5f}")
    print(f"    korel(σ,ρ·X̃)    ort={a('k_srx').mean():+.5f}  "
          f"maks|·|={np.abs(a('k_srx')).max():.5f}")
    # (c) komşu-çizgi kaçağı sınavı
    y = a('c_all')
    print("\n  KAÇAK SINAVI — Cov(σ,X̃)/⟨ρ⟩ neye asılı?")
    for nm, xv in (("A (taşıyıcı faz)", a('A')), ("gap/dres (komşuluk)",
                                                  a('gap')),
                   ("log çizgi gücü", np.log(np.maximum(a('pow'), 1e-300)))):
        print(f"    korel(Cov(σ,X̃)/⟨ρ⟩, {nm:22s}) = "
              f"{np.corrcoef(y, xv)[0,1]:+.4f}")
    print(f"    artık (kapalı form çıkarıldıktan sonra) sd = "
          f"{np.std(y - a('pred')):.5f}  vs sd(y) = {np.std(y):.5f}")
    return R


if __name__ == "__main__":
    out = {}
    for veri, taban in (("son", 0.40), ("A4", 0.40), ("keskin", 0.40)):
        for lo, hi in ((0.60, 0.64), (0.76, 0.80)):
            k = f"{veri}_t{taban}_{lo}"
            out[k] = tani(veri, taban, lo, hi, nmax=30)
            print()
    json.dump({k: v for k, v in out.items()},
              open(C.SCR / "sigma_tani.json", "w"), indent=1)
    print(f"-> {C.SCR/'sigma_tani.json'}")
