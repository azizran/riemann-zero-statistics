"""
165 — T2: F(τ_q; τ_band) — ÇİZGİ-BAŞINA KATKI YOĞUNLUĞU
========================================================
İki ayrı yoğunluk ölçülür:

(1) **Ç1'in yoğunluğu — KAPALI FORMDA.** §2'de türetildi:
       J_2^{Ç1(B+C)} = (y_Q/2)·Σ_{q} Re[h'_q ȳ_q]
    yani Ç1'in çizgi-başına yoğunluğu F(τ_q) = Re[h'_q ȳ_q]'dır ve
    **τ_band'dan BAĞIMSIZDIR** (bant yalnız ortak taşıyıcı çarpanı
    Im[h̄_Q y_Q]/(4⟨ρ⟩) = −½sin(2πτ_Q) ile girer).
    ÇIPLAK (yarım-gap ailesi, SIFIR serbest katsayı):
       F_çıplak(τ) = b_τ B_τ cos(πτ) = 4a² sin²(πτ)cos²(πτ)
                   = **a² sin²(2πτ)**            (163 §2e ile aynı)
    Yarı-analitik sınav: F_emp/a² ↔ c₁ sin²(2πτ) [1 katsayı] ve
    c₁ sin²(2πτ) + c₂ sin⁴(πτ) [2 katsayı].

(2) **TOPLAMIN yoğunluğu.** J_2 = Σ_q ⟨E X^{(q)} X e^{−iω_Q s}⟩ ⇒
       F_tot(q; Q) = ½[ y_q Ĝ(ω_q−ω_Q) + ȳ_q Ĝ(−ω_q−ω_Q) ],
       Ĝ(ν) = ⟨(E·X) e^{iν s}⟩
    (Ç1 bunun içindedir; fark Ç3'ün yoğunluğudur.) Bant bağımlılığı
    BURADA ölçülür.

Kullanım:  165_T2.py <veri> [taban] [kaynak] [tau_c]
"""
import importlib
import json
import sys
import time
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
sys.path.insert(0, str(QM / "165_configs"))
sys.path.insert(0, str(QM / "163_configs"))
K = importlib.import_module("165_cekirdek")
C163 = importlib.import_module("163_cekirdek")

TWO_PI = 2 * np.pi
SCR = K.SCR


def fit(tau, F, a, msk):
    """F/a² ↔ c₁ sin²(2πτ) (+ c₂ sin⁴(πτ)) — en küçük kareler."""
    y = F[msk] / a[msk] ** 2
    s2 = np.sin(TWO_PI * tau[msk]) ** 2
    s4 = np.sin(np.pi * tau[msk]) ** 4
    w = a[msk] ** 2                       # katkı ağırlıklı (a² F ölçüsü)
    c1 = float(np.sum(w * y * s2) / np.sum(w * s2 * s2))
    A = np.vstack([s2, s4]).T
    Wd = np.diag(w) if len(w) < 4000 else None
    AtA = (A * w[:, None]).T @ A
    Atb = (A * w[:, None]).T @ y
    c12 = np.linalg.solve(AtA, Atb)
    r0 = float(np.sum(w * (y - c1 * s2) ** 2) / np.sum(w * y * y))
    r1 = float(np.sum(w * (y - A @ c12) ** 2) / np.sum(w * y * y))
    return c1, c12.tolist(), r0, r1


def kos(veri, taban=0.40, kaynak="olculen", tau_c=0.95, nband=3):
    t0 = time.time()
    Y = K.gaz(veri, taban, 4000)
    Mo = K.Model165(veri, taban, 4000, 0.95, kaynak, Y=Y, tau_c=tau_c)
    Mo.sec(kaynak, tau_c).alanlar(kmax=2)
    tau, a, w = Mo.M["tau"], Mo.M["a"], Mo.M["w"]
    hp, y = Mo.hp, Mo.y
    F = (hp * np.conj(y)).real
    Fb = a ** 2 * np.sin(TWO_PI * tau) ** 2
    msk = (tau > taban + 1e-12) & Mo.msk
    print(f"=== 165-T2 {veri} kaynak={kaynak} τ_c={tau_c} ===", flush=True)
    print(f"  Σ F_emp = {F[msk].sum():.6f}   Σ F_çıplak = {Fb[msk].sum():.6f}"
          f"   oran = {F[msk].sum()/Fb[msk].sum():.4f}", flush=True)
    c1, c12, r0, r1 = fit(tau, F, a, msk)
    print(f"  F/a² ≈ c₁ sin²(2πτ):  c₁ = {c1:.4f}  (artık {r0:.4f})")
    print(f"  F/a² ≈ c₁ sin²(2πτ) + c₂ sin⁴(πτ):  c = {c12[0]:.4f}, "
          f"{c12[1]:.4f}  (artık {r1:.4f})", flush=True)
    print("   τ-bant   ΣF_emp     ΣF_çıplak   oran    ncizgi")
    kut = []
    for lo in np.arange(0.40, 0.95, 0.05):
        m = msk & (tau > lo) & (tau <= lo + 0.05)
        if not m.any():
            continue
        kut.append(dict(lo=float(lo), F=float(F[m].sum()),
                        Fb=float(Fb[m].sum()), n=int(m.sum())))
        print(f"   {lo:.2f}-{lo+0.05:.2f} {F[m].sum():+.6f} "
              f"{Fb[m].sum():+.6f} {F[m].sum()/Fb[m].sum() if Fb[m].sum() else float('nan'):+7.3f}"
              f" {int(m.sum()):6d}", flush=True)

    # --- (2) TOPLAMIN yoğunluğu (temsili çizgilerde) ------------------
    G = Mo.E * Mo.X
    bant = [b for b in K.IZGARA_T1 if b[0] >= 0.52 - 1e-9]
    ban = C163.bant_adaylari(Y, bant)
    sec = []
    step = max(1, len(ban) // nband)
    for b in ban[::step][:nband]:
        if not b["cizgi"]:
            continue
        sec.append(max(b["cizgi"], key=lambda r: r["gp"]))
    tot = []
    for r in sec:
        WQ = float(r["w"])
        # tayf_s frekans W'de ⟨G e^{−iWs}⟩ verir; bize Ĝ(ν)=⟨G e^{+iνs}⟩
        # gerekiyor ⇒ W = −ν. G REEL olduğundan conj ile de alınabilir.
        nu = np.concatenate((w - WQ, -w - WQ))
        tm = time.time()
        Gh = K.tayf_s(Mo.s, [G], nu, fblok=192)[0] / 2.0
        n = len(w)
        Ftot = 0.5 * (y * np.conj(Gh[:n]) + np.conj(y) * np.conj(Gh[n:]))
        o = C163.olc_cizgi(Y, WQ, kmax=2)
        hQ, rm = o["h"], o["rho_ort"]
        s2F = (np.conj(hQ) * Ftot).imag / (2 * rm)
        print(f"  [toplam yoğunluğu] Q={r['q']} τ_Q={r['tau']:.4f} "
              f"({time.time()-tm:.0f}s)  Σ s2F = {s2F.sum():+.6f} "
              f"(ölçülen s2 = {o['s2']:+.6f})", flush=True)
        kb = []
        for lo in np.arange(0.0, 0.95, 0.10):
            m = Mo.msk & (tau > lo) & (tau <= lo + 0.10)
            if not m.any():
                continue
            kb.append(dict(lo=float(lo), s2=float(s2F[m].sum()),
                           n=int(m.sum())))
        print("     " + "  ".join(f"{d['lo']:.1f}:{d['s2']:+.5f}"
                                  for d in kb), flush=True)
        tot.append(dict(q=int(r["q"]), tau=r["tau"], s2_olc=o["s2"],
                        s2_top=float(s2F.sum()), kutu=kb))

    out = dict(veri=veri, kaynak=kaynak, tau_c=tau_c, c1=c1, c12=c12,
               art1=r0, art2=r1, kutu=kut,
               F=F.tolist(), Fb=Fb.tolist(), tau=tau.tolist(),
               a=a.tolist(), q=Mo.M["q"].tolist(), toplam=tot,
               sure_s=time.time() - t0)
    SCR.mkdir(parents=True, exist_ok=True)
    p = SCR / f"T2_{veri}_{kaynak}_tc{tau_c}.json"
    p.write_text(json.dumps(out, indent=1))
    print(f"-> {p}  ({(time.time()-t0)/60:.1f} dk)", flush=True)


if __name__ == "__main__":
    kos(sys.argv[1],
        float(sys.argv[2]) if len(sys.argv) > 2 else 0.40,
        sys.argv[3] if len(sys.argv) > 3 else "olculen",
        float(sys.argv[4]) if len(sys.argv) > 4 else 0.95)
