"""
163 / T2 — YAKIN-REZONANS KATMANI: rezonans koşulunu pencere çözünürlüğüne
kadar gevşetmek.   Kullanım: 163_T2.py <veri> [taban] [tau_ust]

ÇEKİRDEK FİKİR.  §163_cekirdek'in türetimi

    T_c = (1/16)·Σ_{q₁} [ h₁·Γ(2π(τ₁−τ_Q),k)·Ψ(ω_Q−ω₁)
                        + h̄₁·Γ(−2π(τ₁+τ_Q),k)·Ψ(ω_Q+ω₁) ]

biçimindedir; Ψ(μ) X̃'nin ÇİFT-TAYFIdır:

    Ψ(μ) = Σ_{a,b} Σ_{ε_a,ε_b} x_a^{(ε_a)} x_b^{(ε_b)} · δ( ε_aω_a + ε_bω_b − μ )

TAM ÇARPIMSAL üçlüler δ'nın SIFIR genişlikli halidir (ω_a ± ω_b = μ tam
olarak, yani q_a^{±}q_b^{±} = Q q₁^{−1} çarpımsal özdeşliği). YAKIN-rezonans
δ'yı W genişliğinde bir kutuya açar: |ε_aω_a+ε_bω_b − μ| < W.

Bu, ω-EKSENİNDE BİR IZGARAYA binleyip FFT ile hesaplanır:
    g[j] = Σ_{ω_q ∈ kutu j} x_q
    D = IFFT(|FFT(g)|²)          (fark kanadı Σ x_a x̄_b δ(ω_a−ω_b−μ))
    S = IFFT(FFT(g)²)            (toplam kanadı Σ x_a x_b δ(ω_a+ω_b−μ))
    Ψ(μ) = 2[ D(μ) + S(μ)·1{μ>0} + conj(S(−μ))·1{μ<0} ]
Kutu genişliği W = pencere çözünürlüğünün katıdır: W = j·dres, dres = 2π/T.

W → dres:  "yakın-rezonans dahil" öngörü.
W = 0   :  yalnız TAM çarpımsal üçlüler (163_kos.py'nin T1 öngörüsü).
"""
import importlib
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
K = importlib.import_module("163_cekirdek")
TWO_PI = 2 * np.pi
KMAX = 2


def cift_tayf(w, x, W, wmax, nfft=None):
    """Ψ(μ) — X̃ çizgilerinin W kutulu çift-tayfı. Döner: (mu_grid, Psi)."""
    nb = int(np.ceil(wmax / W)) + 2
    n = 1 << int(np.ceil(np.log2(4 * nb)))
    g = np.zeros(n, dtype=complex)
    idx = np.clip((w / W + 0.5).astype(int), 0, nb - 1)
    np.add.at(g, idx, x)
    G = np.fft.fft(g)
    D = np.fft.ifft(G * np.conj(G))          # D[m] = Σ x_a x̄_b, ω_a−ω_b = mW
    Sm = np.fft.ifft(G * G)                  # S[m] = Σ x_a x_b,  ω_a+ω_b = mW
    return n, D, Sm


def Psi(n, D, Sm, mu, W):
    """Ψ(μ) — vektörel arama."""
    m = np.rint(mu / W).astype(int)
    pos = m >= 0
    out = np.empty(len(mu), dtype=complex)
    im = np.mod(m, n)
    Dv = np.where(pos, D[im], np.conj(D[np.mod(-m, n)]))
    Sv = np.where(pos, Sm[im], np.conj(Sm[np.mod(-m, n)]))
    out = 2 * (Dv + Sv)
    return out


def kos(veri, taban=0.40, tau_ust=0.95, izgara="t1", carpanlar=(0, 1, 2, 5, 10, 25, 60)):
    t0 = time.time()
    Y = K.gaz(veri, taban)
    L = Y.L
    dres = TWO_PI / (Y.mid[-1] - Y.mid[0])
    mer = K.merdiven(L, 1.00)
    sel = mer["tau"] <= tau_ust
    Wl, TAUl = mer["w"][sel], mer["tau"][sel]
    d = np.load(K.SCR / f"sp_{veri}_{tau_ust}.npz")
    h1, x1 = d["h"], d["x"]
    PT = K.PhiTablo(Y)
    print(f"[{veri}] L={L:.4f} dres={dres:.3e} çizgi={len(Wl)}  "
          f"({time.time()-t0:.0f}s)", flush=True)

    # W-bağımsız kısım: her (q₁, Q) çifti için Γ
    bantlar = [b for b in (K.IZGARA_T1 if izgara == "t1" else K.IZGARA158)
               if b[0] >= taban - 1e-9]
    BA = K.bant_adaylari(Y, bantlar)

    # çift-tayflar (kutu genişliği başına bir kez)
    TAB = {}
    for c in carpanlar:
        if c == 0:
            continue
        Wb = c * dres
        TAB[c] = (Wb,) + cift_tayf(Wl, x1, Wb, L * 1.05)
    print(f"  çift-tayflar hazır ({time.time()-t0:.0f}s)", flush=True)

    satir = []
    for bd in BA:
        if bd["kul"] == 0:
            continue
        rec = []
        for cz in bd["cizgi"]:
            r = {}
            for Wf, et in ((cz["w"], "on"), (cz["w"] + cz["gap"] / 2, "off")):
                o = K.olc_cizgi(Y, Wf, KMAX)
                tQ = Wf / L
                msk = np.abs(TAUl - tQ) > 1e-12
                # TAM çarpımsal (W=0): evrensel üçlü
                pr0 = K.ongor_cizgi(o["h"], o["x"], tQ, TAUl[msk], h1[msk],
                                    x1[msk], PT, KMAX)
                r[f"A_{et}"] = o["A"]
                r[f"pow_{et}"] = o["pow"]
                r[f"s2_{et}"] = o["s2"]
                r[f"u1_{et}"] = o["u1"]
                r[f"s2_W0_{et}"] = pr0["s2_pred"]
                r[f"u1_W0_{et}"] = pr0["u1_pred"]
                # pencereli: Ψ ile
                th_f = TWO_PI * (TAUl - tQ)
                th_s = -TWO_PI * (TAUl + tQ)
                Gf2 = PT.Gam(th_f, 2)
                Gs2 = PT.Gam(th_s, 2)
                Gf1 = PT.Gam(th_f, 1)
                Gs1 = PT.Gam(th_s, 1)
                muf = Wf - Wl
                mus = Wf + Wl
                hb = np.conj(o["h"])
                n2 = abs(o["h"]) ** 2
                Ta2 = o["h"] / 2 * PT.mom[2]
                Ta1 = o["h"] / 2 * PT.mom[1]
                for c in carpanlar:
                    if c == 0:
                        continue
                    Wb, nf, D, Sm = TAB[c]
                    Pf = Psi(nf, D, Sm, muf, Wb)
                    Ps = Psi(nf, D, Sm, mus, Wb)
                    J2 = Ta2 + (np.sum(h1 * Gf2 * Pf)
                                + np.sum(np.conj(h1) * Gs2 * Ps)) / 16.0
                    J1 = Ta1 + (np.sum(h1 * Gf1 * Pf)
                                + np.sum(np.conj(h1) * Gs1 * Ps)) / 16.0
                    r[f"s2_W{c}_{et}"] = 2 * float((hb * J2).imag) / n2
                    r[f"u1_W{c}_{et}"] = 2 * float((hb * J1).real) / n2
            r.update(q=cz["q"], tau=cz["tau"], gp=cz["gp"])
            rec.append(r)
        pN = np.array([r["pow_on"] for r in rec])
        pO = np.array([r["pow_off"] for r in rec])
        Aon = np.array([r["A_on"] for r in rec])
        Aof = np.array([r["A_off"] for r in rec])
        tv = np.array([r["tau"] for r in rec])
        u = pN - pO
        payda = float(u.sum())

        def AGG(key, pw):
            von = np.array([r[f"{key}_on"] for r in rec])
            vof = np.array([r[f"{key}_off"] for r in rec])
            return float((pN * Aon ** pw * von - pO * Aof ** pw * vof).sum() / payda)

        dd = dict(tau=bd["tau"], kul=bd["kul"],
                  tau_eff=float((tv * u).sum() / payda),
                  A2s2=AGG("s2", 2), A1u1=AGG("u1", 1),
                  A2s2_W0=AGG("s2_W0", 2), A1u1_W0=AGG("u1_W0", 1))
        for c in carpanlar:
            if c == 0:
                continue
            dd[f"A2s2_W{c}"] = AGG(f"s2_W{c}", 2)
            dd[f"A1u1_W{c}"] = AGG(f"u1_W{c}", 1)
        satir.append(dd)
        wl = "  ".join(f"W{c}={dd[f'A2s2_W{c}']:+.5f}" for c in carpanlar if c)
        print(f"  τ={bd['tau']:.4f} A²s2={dd['A2s2']:+.6f}  W0={dd['A2s2_W0']:+.6f}  "
              f"{wl}   [{time.time()-t0:.0f}s]", flush=True)
    out = dict(veri=veri, taban=taban, dres=dres, carpanlar=list(carpanlar),
               bantlar=satir, sure_s=time.time() - t0)
    p = K.SCR / f"T2_{veri}_t{taban}.json"
    p.write_text(json.dumps(out, indent=1))
    print(f"-> {p}  ({(time.time()-t0)/60:.1f} dk)", flush=True)


if __name__ == "__main__":
    kos(sys.argv[1],
        float(sys.argv[2]) if len(sys.argv) > 2 else 0.40,
        float(sys.argv[3]) if len(sys.argv) > 3 else 0.95)
