"""
163 / T4 — İŞARET MUHASEBESİ + TAŞIYICI TANISI
Kullanım: 163_T4.py <veri> [taban]

İki ölçüm:

(A) τ-BÖLGESİ MATRİSİ.  X̃ kendi merdiven çizgilerinden τ-dilimlerine
    ayrıştırılır:  X̃ = Σ_j X^{(j)} + X^{art}.  Bant-toplu
        A²·⟨σ X^{(j)} X^{(k)}⟩/⟨ρ⟩
    matrisi, ⟨σX̃²⟩'nin (ve dolayısıyla b'nin) İŞARETİNİ hangi çizgi-çifti
    bölgesinin ürettiğini doğrudan gösterir. Çıplak üçlü-toplamın öngörüsü
    "yalnız (Q, q₁) evrensel çifti" olduğundan, matrisin köşegen-dışı
    yoğunluğu doğrudan tarak-aracılı/yakın-rezonans payını ölçer.

(B) DÜZGÜN-TARAK KONTROLÜ.  Aynı çizgi genlik/fazlarıyla alanlar DÜZGÜN
    bir tarakta (m̄_n = m_0 + ḡn) yeniden kurulur; orada ⟨e^{iνm}⟩ = 0
    (ν ≠ 0) olduğundan TARAK-ARACILI katman ÖLÜR. Gerçek tarakla farkı
    G-katmanının büyüklüğüdür.
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
DILIM = [(0.0, 0.20), (0.20, 0.40), (0.40, 0.55), (0.55, 0.70),
         (0.70, 0.85), (0.85, 0.95)]


def sentez(amp, Wv, mm, blok=20000, fblok=500):
    out = np.zeros(len(mm))
    for f0 in range(0, len(Wv), fblok):
        fs = slice(f0, min(f0 + fblok, len(Wv)))
        Wc, A = Wv[fs], amp[fs]
        for s in range(0, len(mm), blok):
            sl = slice(s, min(s + blok, len(mm)))
            arg = np.outer(mm[sl], Wc)
            out[sl] += np.cos(arg) @ A.real - np.sin(arg) @ A.imag
            del arg
    return out


def kos(veri, taban=0.40, tau_ust=0.95, hedef=(0.54, 0.62, 0.70)):
    t0 = time.time()
    Y = K.gaz(veri, taban)
    L = Y.L
    d = np.load(K.SCR / f"sp_{veri}_{tau_ust}.npz")
    Wl, x1, TAUl = d["W"], d["x"], d["tau"]
    m = Y.m0
    XD = []
    for lo, hi in DILIM:
        s = (TAUl > lo) & (TAUl <= hi)
        v = sentez(x1[s], Wl[s], m)
        XD.append(v - v.mean())
        print(f"  dilim τ∈({lo},{hi}]: {int(s.sum())} çizgi, Var={np.var(v):.5f} "
              f"[{time.time()-t0:.0f}s]", flush=True)
    Xart = Y.Xtil0 - sum(XD)
    XD.append(Xart - Xart.mean())
    ad = [f"{lo}-{hi}" for lo, hi in DILIM] + ["artık"]
    nD = len(XD)
    print(f"  artık Var={np.var(Xart):.5f}  (toplam Var(X̃0)={np.var(Y.Xtil0):.5f})",
          flush=True)

    bantlar = [b for b in K.IZGARA_T1 if b[0] >= taban - 1e-9]
    BA = K.bant_adaylari(Y, bantlar)
    cikti = []
    for bd in BA:
        if bd["kul"] == 0 or round(bd["tau"], 4) not in [round(h, 4) for h in hedef]:
            continue
        acc = {}
        rec = []
        for cz in bd["cizgi"]:
            r = {}
            for Wf, et in ((cz["w"], "on"), (cz["w"] + cz["gap"] / 2, "off")):
                cw, sw = np.cos(Wf * m), np.sin(Wf * m)
                zc = complex(2 * np.mean(Y.e0 * cw), -2 * np.mean(Y.e0 * sw))
                cw1, sw1 = np.cos(Wf * Y.mid[1:]), np.sin(Wf * Y.mid[1:])
                c1r, c1i = Y.e1 * cw1, -Y.e1 * sw1
                zb = zc / 2.0
                rho = c1r * zb.real + c1i * zb.imag
                sig = c1i * zb.real - c1r * zb.imag
                rm = float(rho.mean())
                r[f"pow_{et}"] = float(abs(zc) ** 2)
                r[f"A_{et}"] = TWO_PI * Wf / L
                r[f"s2_{et}"] = float(np.mean(sig * Y.Xtil0 ** 2) / rm)
                for j in range(nD):
                    for k in range(j, nD):
                        f = 1.0 if j == k else 2.0
                        r[f"M{j}_{k}_{et}"] = float(f * np.mean(sig * XD[j] * XD[k]) / rm)
            r.update(tau=cz["tau"], q=cz["q"])
            rec.append(r)
        pN = np.array([r["pow_on"] for r in rec])
        pO = np.array([r["pow_off"] for r in rec])
        Aon = np.array([r["A_on"] for r in rec])
        Aof = np.array([r["A_off"] for r in rec])
        payda = float((pN - pO).sum())

        def AGG(key):
            von = np.array([r[f"{key}_on"] for r in rec])
            vof = np.array([r[f"{key}_off"] for r in rec])
            return float((pN * Aon ** 2 * von - pO * Aof ** 2 * vof).sum() / payda)

        M = np.zeros((nD, nD))
        for j in range(nD):
            for k in range(j, nD):
                M[j, k] = AGG(f"M{j}_{k}")
        tot = AGG("s2")
        print(f"\n[{veri}] bant τ={bd['tau']:.2f}  A²s2={tot:+.6f}   "
              f"matris toplamı={M.sum():+.6f}", flush=True)
        print("      " + "".join(f"{a:>10}" for a in ad))
        for j in range(nD):
            print(f"{ad[j]:>6}" + "".join(
                (f"{M[j,k]:+10.5f}" if k >= j else " " * 10) for k in range(nD)),
                flush=True)
        cikti.append(dict(tau=bd["tau"], A2s2=tot, M=M.tolist(), ad=ad))
    p = K.SCR / f"T4_{veri}_t{taban}.json"
    p.write_text(json.dumps(dict(veri=veri, dilim=DILIM, bantlar=cikti,
                                 sure_s=time.time() - t0), indent=1))
    print(f"\n-> {p}  ({(time.time()-t0)/60:.1f} dk)", flush=True)


if __name__ == "__main__":
    kos(sys.argv[1], float(sys.argv[2]) if len(sys.argv) > 2 else 0.40)
