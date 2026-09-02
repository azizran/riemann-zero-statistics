"""
163 — DENETİM.  Kullanım: 163_dogrulama.py

V1  BİT DÜZEYİ: 163'ün ölçüm zinciri 160'ın kayıtlı JSON'unu yeniden
    üretiyor mu? (t1 ızgarası, taban 0.40, üç gaz, 27 bant; A²s2, A·s1,
    A³s3, A·κ1, s0, τ_eff, kullanılan çizgi sayısı)
V2  Φ tablosu: Φ_j(0) = ⟨X̃0^j⟩ (binleme hatası)
V3  ⟨σ⟩/⟨ρ⟩ = s0 gerçekten sıfır mu (kendi-çizgi teriminin s_k'ya katkısının
    sıfır olduğunun ölçülmüş karşılığı)
V4  Öngörü cebrinin normalizasyonu: u0_pred = 1 mi (birinci-mertebe hatası)
"""
import importlib
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
K = importlib.import_module("163_cekirdek")
S160 = K.SCR.parent / "160"

ANAHTAR = [("A2s2", "A2s2"), ("A1s1", "A1s1"), ("A3s3", "A3s3"),
           ("A1k1", "A1u1"), ("s0", "A0s0")]

if __name__ == "__main__":
    print("=" * 74)
    print("V1 — 163'ün ölçümü = 160'ın kayıtlı JSON'u mu? (t1, taban 0.40)")
    tot = 0
    for v in ("son", "keskin", "A4"):
        a = json.loads((S160 / f"S_{v}_t0.4_t1.json").read_text())
        b = json.loads((K.SCR / f"P_{v}_t0.4_t1.json").read_text())
        A = {r["tau"]: r for r in a["bantlar"] if r.get("olculdu")}
        B = {r["tau"]: r for r in b["bantlar"] if r.get("olculdu")}
        ort = sorted(set(A) & set(B))
        d = {}
        for t in ort:
            for k160, k163 in ANAHTAR:
                d.setdefault(k160, []).append(abs(A[t]["KUM"][k160] - B[t][k163]))
            d.setdefault("tau_eff", []).append(abs(A[t]["tau_eff"] - B[t]["tau_eff"]))
            d.setdefault("kul", []).append(abs(A[t]["kul"] - B[t]["kul"]))
            tot += 1
        print(f"  {v}: {len(ort)} ortak bant — maks fark " +
              "  ".join(f"{k}={max(x):.1e}" for k, x in d.items()))
    print(f"  TOPLAM {tot} bant.  (160'ın A²κ₂'si KÜMÜLANTTIR, μ₂−μ₁²; "
          f"163'ün A²u₂'si HAM momenttir — karşılaştırmaya girmez.)")

    print("=" * 74)
    for v in ("son", "keskin", "A4"):
        Y = K.gaz(v, 0.40)
        PT = K.PhiTablo(Y)
        print(f"V2 [{v}] Φ_j(0) − ⟨X̃0^j⟩ maks = {PT.hata0:.2e}   "
              f"⟨X̃0²⟩={PT.mom[2]:.6f} ⟨X̃0³⟩={PT.mom[3]:.6f}")
        d = np.load(K.SCR / f"sp_{v}_0.95.npz")
        tau, x = d["tau"], d["x"]
        mer = K.merdiven(Y.L, 1.00)
        sel = mer["tau"] <= 0.95
        print(f"      Σ_{{τ>0.4}}|x_q|² = {float(np.sum(np.abs(x[tau>0.4])**2)):.6f}"
              f"   çıplak Σ b_qB_qcos(πτ) = "
              f"{float(np.sum((mer['b']*mer['B']*np.cos(np.pi*mer['tau']))[sel][mer['tau'][sel]>0.4])):.6f}")
        for tq in (0.54, 0.62, 0.70):
            i = int(np.argmin(np.abs(tau - tq)))
            o = K.olc_cizgi(Y, np.log(float(d["q"][i])), 2)
            msk = np.abs(tau - tau[i]) > 1e-12
            pr = K.ongor_cizgi(o["h"], o["x"], tau[i], tau[msk], d["h"][msk],
                               d["x"][msk], PT, 2)
            print(f"      V3/V4 τ={tau[i]:.4f}: s0_ölç={o['s0']:+.2e} "
                  f"(≡0 beklenir)   u0_öngörü={pr['u0_pred']:.4f} (≡1 beklenir)")
