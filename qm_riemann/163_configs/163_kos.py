"""
163 — KOŞU SÜRÜCÜSÜ.  Kullanım: 163_kos.py <veri> [taban] [tau_ust]

Her bant için 160'ın ADAY LİSTESİYLE (tohum 21, gap filtresi, on/off ara-nokta)
hem ÖLÇÜM (u_k, s_k) hem ÜÇLÜ-TOPLAM ÖNGÖRÜSÜ hesaplanır ve 160'ın bant
birleştirmesiyle (Σ(pN·A^p·v_on − pO·A'^p·v_off)/Σ(pN−pO)) toplanır.
Çıktı: scratchpad/163/P_<veri>_t<taban>.json
"""
import importlib
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
K = importlib.import_module("163_cekirdek")

KMAX = 3


def kos(veri, taban=0.40, tau_ust=0.95, izgara="t1"):
    t0 = time.time()
    Y = K.gaz(veri, taban)
    L = Y.L
    mer = K.merdiven(L, 1.00)
    sel = mer["tau"] <= tau_ust
    Wl, TAUl, Ql = mer["w"][sel], mer["tau"][sel], mer["q"][sel]
    p = K.SCR / f"sp_{veri}_{tau_ust}.npz"
    if p.exists():
        d = np.load(p)
        h1, x1 = d["h"], d["x"]
    else:
        h1, x1, _ = K.tayf(Y, Wl)
        np.savez(p, h=h1, x=x1, W=Wl, tau=TAUl, q=Ql)
    print(f"[{veri}] L={L:.4f} sA2={Y.sA2:.5f}  merdiven {len(Wl)} çizgi "
          f"(τ≤{tau_ust})  {time.time()-t0:.0f}s", flush=True)
    PT = K.PhiTablo(Y)
    print(f"  Φ tablosu: hata0={PT.hata0:.1e}  ⟨X0²⟩={PT.mom[2]:.5f} "
          f"⟨X0³⟩={PT.mom[3]:.5f}", flush=True)

    bantlar = [b for b in (K.IZGARA_T1 if izgara == "t1" else K.IZGARA158)
               if b[0] >= taban - 1e-9]
    BA = K.bant_adaylari(Y, bantlar)
    satir = []
    for bd in BA:
        if bd["kul"] == 0:
            satir.append(dict(tau=bd["tau"], olculdu=False))
            continue
        rec = []
        for c in bd["cizgi"]:
            r = {}
            for Wf, et in ((c["w"], "on"), (c["w"] + c["gap"] / 2, "off")):
                o = K.olc_cizgi(Y, Wf, KMAX)
                tQ = Wf / L
                msk = np.abs(TAUl - tQ) > 1e-12
                pr = K.ongor_cizgi(o["h"], o["x"], tQ, TAUl[msk], h1[msk],
                                   x1[msk], PT, KMAX)
                r[f"A_{et}"] = o["A"]
                r[f"pow_{et}"] = o["pow"]
                r[f"absh_{et}"] = abs(o["h"])
                r[f"absx_{et}"] = abs(o["x"])
                r[f"argh_{et}"] = float(np.angle(o["h"]))
                r[f"argx_{et}"] = float(np.angle(o["x"]))
                for k in range(KMAX + 1):
                    for nm in (f"u{k}", f"s{k}"):
                        r[f"{nm}_{et}"] = o[nm]
                    for nm in (f"u{k}_pred", f"s{k}_pred", f"s{k}_fark",
                               f"s{k}_topl", f"s{k}_Tb", f"u{k}_fark"):
                        r[f"{nm}_{et}"] = pr[nm]
            r.update(q=c["q"], tau=c["tau"], grup=c["grup"], gp=c["gp"])
            rec.append(r)
        pN = np.array([r["pow_on"] for r in rec])
        pO = np.array([r["pow_off"] for r in rec])
        Aon = np.array([r["A_on"] for r in rec])
        Aof = np.array([r["A_off"] for r in rec])
        tv = np.array([r["tau"] for r in rec])
        u = pN - pO
        payda = float(u.sum())
        tau_eff = float((tv * u).sum() / payda)

        def AGG(key, pw):
            von = np.array([r[f"{key}_on"] for r in rec])
            vof = np.array([r[f"{key}_off"] for r in rec])
            return float((pN * Aon ** pw * von - pO * Aof ** pw * vof).sum() / payda)

        d = dict(tau=bd["tau"], lo=bd["lo"], hi=bd["hi"], kul=bd["kul"],
                 N=bd["N"], olculdu=True, tau_eff=tau_eff,
                 payda=payda, A_eff=float(2 * np.pi * tau_eff),
                 rho=payda / float(sum(c["gp"] for c in bd["cizgi"])))
        for k in range(KMAX + 1):
            d[f"A{k}u{k}"] = AGG(f"u{k}", k)
            d[f"A{k}s{k}"] = AGG(f"s{k}", k)
            d[f"A{k}u{k}_pred"] = AGG(f"u{k}_pred", k)
            d[f"A{k}s{k}_pred"] = AGG(f"s{k}_pred", k)
            d[f"A{k}s{k}_fark"] = AGG(f"s{k}_fark", k)
            d[f"A{k}s{k}_topl"] = AGG(f"s{k}_topl", k)
            d[f"A{k}s{k}_Tb"] = AGG(f"s{k}_Tb", k)
            # yalnız ON kanadı (ara-nokta çıkarımı olmadan)
            von = np.array([r[f"s{k}_on"] for r in rec])
            vpr = np.array([r[f"s{k}_pred_on"] for r in rec])
            d[f"A{k}s{k}_on"] = float((pN * Aon ** k * von).sum() / payda)
            d[f"A{k}s{k}_on_pred"] = float((pN * Aon ** k * vpr).sum() / payda)
        d["cizgi"] = rec
        satir.append(d)
        print(f"  τ={bd['tau']:.4f} τe={tau_eff:.4f} kul={bd['kul']:3d} "
              f"A²s2={d['A2s2']:+.6f} öngörü={d['A2s2_pred']:+.6f} "
              f"(fark {d['A2s2_fark']:+.6f} topl {d['A2s2_topl']:+.6f})  "
              f"A·S={d['A1u1']:+.5f} öngörü={d['A1u1_pred']:+.5f}  "
              f"[{time.time()-t0:.0f}s]", flush=True)
    out = dict(veri=veri, taban=taban, L=L, sA2=Y.sA2, tau_ust=tau_ust,
               izgara=izgara, sure_s=time.time() - t0, bantlar=satir)
    K.SCR.mkdir(parents=True, exist_ok=True)
    q = K.SCR / f"P_{veri}_t{taban}_{izgara}.json"
    q.write_text(json.dumps(out, indent=1))
    print(f"-> {q}  ({(time.time()-t0)/60:.1f} dk)", flush=True)


if __name__ == "__main__":
    veri = sys.argv[1]
    taban = float(sys.argv[2]) if len(sys.argv) > 2 else 0.40
    tu = float(sys.argv[3]) if len(sys.argv) > 3 else 0.95
    izg = sys.argv[4] if len(sys.argv) > 4 else "t1"
    kos(veri, taban, tu, izg)
