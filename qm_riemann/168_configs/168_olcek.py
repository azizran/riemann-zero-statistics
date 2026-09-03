"""
168 — A2(i): ÖLÇEK YASALARININ AYRIŞTIRILMASI (ucuz; yalnız C_*.json)
=====================================================================
Girdi: scratchpad/167/C_<gaz>.json  (yeniden gaz kurulmuyor)

Sorduğu tek soru: 167'nin İKİ ölçek yasası (c(T), c(kesim)) hangi çarpandan
geliyor?  c = KALİB_u2/(W_amp·W_X) ayrışması:

    KALİB_u2 = ölç/öngörü  =  g_cal · ĉ · (W_amp W_X)/(…)
    g_cal    = g_E·g_X²    (165 §7'nin ÇİZGİ-UYUM sızıntısı: g_E =
                            ⟨η E_mod⟩/⟨E_mod²⟩ regresyon katsayısı)
    ĉ        ≡ c / g_cal   (sızıntıdan arındırılmış "saf" ortak kazanç)

Ayrıca A1'in ÖNGÖRÜLERİ sayısal olarak buraya konur:
  * "yerleşim" (occupancy) formu:  c ∝ w_eff = 2π/T   ⇒ c(T/2)/c(T) = 2
  * "aşırı-uyum sızıntısı" formu:  1 + 2·n_çizgi/N    ⇒ T ile küçülür
  * kesim ekseni için merdiven integralleri
        Σ𝒢²(τ_c) ≈ ∫_0^{τc} τ cos²(πτ) dτ ,
        Σ|𝒢|(τ_c) ≈ ∫_0^{τc} |cos(πτ)| dτ          (asal yoğunluğu ile)
"""
import glob
import json
import os

import numpy as np

SCR = ("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
       "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/167")
LO_MIN, LO_MAX, SNR_MIN, R_MIN = 0.52, 0.68, 3.0, 0.98


def yukle():
    D = {}
    for p in sorted(glob.glob(SCR + "/C_*.json")):
        D[os.path.basename(p)[2:-5]] = json.load(open(p))
    return D


def bantlar(d, lo_max=LO_MAX):
    out = []
    for b in d["bant"]:
        if not b.get("olculdu"):
            continue
        if b["lo"] < LO_MIN - 1e-9 or b["lo"] > lo_max + 1e-9:
            continue
        if b["tau_eff"] > 0.85 or b["Ms2"] <= 0:
            continue
        if b["R_bant"] < R_MIN or b["SNR"] < SNR_MIN:
            continue
        out.append(b)
    return out


def logort(v):
    v = np.asarray(v, float)
    return float(np.exp(np.mean(np.log(v))))


# ---------- merdiven integralleri (asal yoğunluğu dn/dω = e^ω/ω) ----------
def I_G2(tc):
    """Σ_r 𝒢_r²  ≈ ∫_0^{tc} τ cos²(πτ) dτ  (𝒢=−τcos(πτ)/(m√q), m=1 baskın,
    asal yoğunluğu e^{τL}/(τL)·L dτ ile 1/q = e^{−τL} tam sadeleşir)."""
    return (tc ** 2 / 4 + tc * np.sin(2 * np.pi * tc) / (4 * np.pi)
            + (np.cos(2 * np.pi * tc) - 1) / (8 * np.pi ** 2))


# NOT: Σ_r |𝒢_r| için kapalı form YOKTUR — yoğunluk × |𝒢| ~ e^{τL/2}|cos πτ|
# olduğundan toplam τ_c ile üstel büyür; her yerde sayısal merdiven kullanılır.


def ladder(L=12.0296, tau_c=0.95, cap=4000):
    """(τ_r, 𝒢_r) — asal kuvvet merdiveni, 165/143 ile aynı a_q = 1/(πm√q)."""
    qmax = int(np.exp(tau_c * L)) + 1
    sieve = np.ones(qmax + 1, bool)
    sieve[:2] = False
    for i in range(2, int(qmax ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = False
    pr = np.nonzero(sieve)[0]
    q, m = [], []
    for p in pr:
        v, k = int(p), 1
        while v <= qmax:
            q.append(v)
            m.append(k)
            v *= int(p)
            k += 1
    q = np.array(q, float)
    m = np.array(m, float)
    tau = np.log(q) / L
    ok = tau <= tau_c + 1e-12
    q, m, tau = q[ok], m[ok], tau[ok]
    a = 1.0 / (np.pi * m * np.sqrt(q))
    G = -np.pi * tau * a * np.cos(np.pi * tau)
    return tau, G, q, m, a


def main():
    D = yukle()
    print("=" * 100)
    print("A2(i) — ÖLÇEK YASALARININ AYRIŞTIRILMASI  (ortak pencere "
          f"lo ∈ [{LO_MIN},{LO_MAX}])")
    print("=" * 100)
    R = {}
    print("\ngaz        nb   c        g_cal    ĉ=c/g_cal  θ=KALİB/g  KALİB "
          "   W_a·W_X   n/N      1+2n/N    T         nline")
    for ad in sorted(D):
        d = D[ad]
        B = bantlar(d)
        if not B:
            continue
        A = d["artik"]
        g = A["gE"] * A["gX"] ** 2
        c = logort([b["c_u2"] for b in B])
        kal = logort([b["KALIB_u2"] for b in B])
        w = logort([b["W_amp"] * b["W_X"] for b in B])
        n, N = A["nline"], d["N"]
        R[ad] = dict(c=c, g=g, ch=c / g, th=kal / g, kal=kal, w=w,
                     n=n, N=N, T=d["T"], nb=len(B), gE=A["gE"], gX=A["gX"],
                     lam=d["lam"] or 1.0, sigC=d["sigChat"],
                     vE=A["varE_mod"] / A["varE_olc"],
                     vX=A["varX_mod"] / A["varX_olc"])
        print(f"{ad:10s} {len(B):2d}  {c:.4f}   {g:.4f}   {c/g:.4f}    "
              f"{kal/g:.4f}    {kal:.4f}   {w:.4f}   {n/N:.4f}   "
              f"{1+2*n/N:.4f}  {d['T']:9.1f} {n:6d}")

    def eksen(baslik, ref, grup, etik):
        print(f"\n--- {baslik} (referans {ref}) ---")
        print(f"  {'gaz':10s} {etik:>14s}  c/c₀    g/g₀    ĉ/ĉ₀   θ/θ₀   "
              "(1+2n/N)₀/(1+2n/N)  A1-yerleşim (2π/T)/(2π/T₀)")
        r0 = R[ref]
        f0 = 1 + 2 * r0["n"] / r0["N"]
        for ad, lab in grup:
            if ad not in R:
                continue
            r = R[ad]
            f = 1 + 2 * r["n"] / r["N"]
            print(f"  {ad:10s} {lab:>14s}  {r['c']/r0['c']:.4f}  "
                  f"{r['g']/r0['g']:.4f}  {r['ch']/r0['ch']:.4f}  "
                  f"{r['th']/r0['th']:.4f}   {f0/f:.4f}"
                  f"               {r0['T']/r['T']:.4f}")

    eksen("EKSEN (b) PENCERE T", "Hkeskin",
          [("HkT2a", "T/2 (ilk)"), ("HkT2b", "T/2 (son)"),
           ("HkT4a", "T/4 (ilk)"), ("HkT4b", "T/4 (2.)")], "pencere")
    eksen("EKSEN (c) MERDİVEN KESİMİ", "Hkeskin",
          [("K090", "keskin≤0.90"), ("K070", "keskin≤0.70"),
           ("E060", "erfc 0.60"), ("HA4", "erfc 0.68"),
           ("son", "gerçek")], "kesim")
    eksen("EKSEN (a) GENLİK ÖLÇEĞİ λ", "Hkeskin",
          [("L070", "λ=0.70"), ("L085", "λ=0.85"), ("L115", "λ=1.15")], "λ")

    # ---------------- A1'in kesim-ekseni merdiven integralleri -----------
    print("\n--- A1: kesim ekseninin merdiven integralleri ---")
    L = D["Hkeskin"]["L"]
    tau, G, q, m, a = ladder(L, 0.95)
    print(f"  merdiven: {len(tau)} çizgi (τ≤0.95, L={L:.4f});  "
          f"Σ𝒢² = {np.sum(G**2):.6f}, Σ|𝒢| = {np.sum(np.abs(G)):.4f}")
    print("  τ_c    Σ𝒢²(sayısal)  ∫τcos²(πτ)dτ   Σ|𝒢|    N_r    "
          "Σ𝒢²/Σ𝒢²(0.95)  Σ|𝒢|/…   N_r/…")
    ref = None
    for tc in (0.60, 0.68, 0.70, 0.80, 0.90, 0.95, 1.00):
        s = tau <= tc + 1e-12
        v = (float(np.sum(G[s] ** 2)), float(np.sum(np.abs(G[s]))),
             int(s.sum()))
        if tc == 0.95:
            ref = v
        print(f"  {tc:.2f}   {v[0]:.6f}      {I_G2(tc):.6f}    {v[1]:.4f}  "
              f"{v[2]:6d}", end="")
        if ref:
            print(f"     {v[0]/ref[0]:.4f}        {v[1]/ref[1]:.4f}   "
                  f"{v[2]/ref[2]:.4f}")
        else:
            print()
    # erfc penceresi (HA4/E060)
    from math import erfc
    for (tc0, dl) in ((0.68, 0.125), (0.60, 0.125)):
        wgt = np.array([0.5 * erfc((t - tc0) / dl) for t in tau])
        print(f"  erfc({tc0},{dl}): Σ(𝒢w)² = "
              f"{np.sum((G*wgt)**2)/ref[0]:.4f}·Σ𝒢²(0.95), "
              f"Σ|𝒢w| = {np.sum(np.abs(G*wgt))/ref[1]:.4f}·Σ|𝒢|(0.95), "
              f"Σw = {wgt.sum()/len(tau):.4f}·N_r")

    # -------- KESİM YASASI: θ = 1 − β·φ  (φ = BOŞ ÇİZGİ KESRİ) ----------
    # φ = gazın merdiveninde OLMAYAN ama modelin (τ_c=0.95) FİT ETTİĞİ
    #     çizgilerin payı:  φ = 1 − Σ_q w_q / N_r    (keskin kesimde
    #     w_q = 1{τ_q ≤ τ_ust}; erfc'de ½erfc((τ−τ_c)/Δ))
    from math import erfc as _erfc
    print("\n--- KESİM YASASI: θ/θ₀ = 1 − β·φ  (φ = boş çizgi kesri) ---")
    print("  gaz        τ_ust/erfc     φ       θ/θ₀ ölç   β=(1−θ/θ₀)/φ   "
          "1−0.218φ öngörü   fark%")
    th0 = R["Hkeskin"]["th"]
    bet = []
    for ad in ("Hkeskin", "son", "L085", "K090", "K070", "E060", "HA4",
               "HkT2a", "HkT2b", "HkT4a", "HkT4b"):
        if ad not in R:
            continue
        d = D[ad]
        pen, tu = d["pen"], d["tau_ust"] or 1.0
        if pen is not None:
            wq = np.array([0.5 * _erfc((t - pen[0]) / pen[1]) for t in tau])
        else:
            wq = (tau <= tu + 1e-12).astype(float)
        phi = float(1 - wq.mean())
        r = R[ad]["th"] / th0
        b = (1 - r) / phi if phi > 1e-6 else float("nan")
        if phi > 0.1:
            bet.append(b)
        pred = 1 - 0.218 * phi
        print(f"  {ad:10s} {str(pen or f'≤{tu:.2f}'):>13s}  {phi:.4f}  "
              f"{r:9.4f}   {b:11.4f}    {pred:12.4f}   "
              f"{100*(pred/r-1):+7.2f}")
    if bet:
        print(f"  ⇒ β = {np.mean(bet):.4f} ± {np.std(bet, ddof=1) if len(bet)>1 else 0:.4f} "
              f"({len(bet)} kesim gazı)")

    json.dump({k: {kk: float(vv) for kk, vv in v.items()}
               for k, v in R.items()},
              open(os.path.join(SCR, "..", "168", "168_olcek.json"), "w"), indent=1)
    print("\n-> ../168/168_olcek.json")


if __name__ == "__main__":
    main()
