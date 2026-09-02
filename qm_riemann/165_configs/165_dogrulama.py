"""
165 — DENETİM.  Kullanım: 165_dogrulama.py

V1  BİT DÜZEYİ (160): 165'in ölçüm zinciri 160'ın kayıtlı JSON'unu
    yeniden üretiyor mu? (t1 ızgarası, taban 0.40; A²s2, A·s1, A³s3,
    A·κ1, s0, τ_eff, çizgi sayısı)
V2  KALEM (163): Ç1'in ÇIPLAK limiti gerçekten 163 §2e'nin kapalı formu
    mu? `kanal1`, h'=b, y=B e^{−iπτ}, h_Q=b_Q, ⟨ρ⟩=b_Q²/4 ile
    −½sin(2πτ_Q)Σa²sin²(2πτ)'yı vermeli.
V3  SADAKAT (164): 165'in s-sitesi tayfı 164 §3g'nin çizgi oranlarını
    yeniden üretiyor mu? (|c_q|/b_q, arg c_q)
V4  Ç1'in KAPALI TOPLAMI = doğrudan sayım (küçük merdivende, terim terim)
V5  MODEL alanı ile ÖZDEŞLİK: ⟨E X² e^{−iWs}⟩ = Σ_terim coef·κ(Δ)
    (165_kanal.py V4'te alt-pencerede bit düzeyinde)
"""
import importlib
import json
import sys
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
sys.path.insert(0, str(QM / "165_configs"))
sys.path.insert(0, str(QM / "163_configs"))
K = importlib.import_module("165_cekirdek")
C163 = importlib.import_module("163_cekirdek")

TWO_PI = 2 * np.pi
S160 = K.SCR.parent / "160"
S164 = K.SCR.parent / "164"
ANAHTAR = [("A2s2", "A2s2"), ("A1s1", "A1s1"), ("A3s3", "A3s3"),
           ("A1k1", "A1u1"), ("s0", "A0s0")]


def v1():
    print("=" * 74)
    print("V1 — 165'in ölçümü = 160'ın kayıtlı JSON'u mu? (t1, taban 0.40)")
    for v in ("son",):
        pa = S160 / f"S_{v}_t0.4_t1.json"
        pb = K.SCR / f"K_{v}_t0.4.json"
        if not (pa.exists() and pb.exists()):
            print(f"  {v}: dosya yok, atlandı")
            continue
        a = json.loads(pa.read_text())
        b = json.loads(pb.read_text())
        kk = sorted(b["kosum"])[0]
        A = {r["tau"]: r for r in a["bantlar"] if r.get("olculdu")}
        B = {r["tau"]: r for r in b["kosum"][kk]["bant"] if r.get("olculdu")}
        ort = sorted(set(A) & set(B))
        d = {}
        for t in ort:
            for k160, k165 in ANAHTAR:
                d.setdefault(k160, []).append(
                    abs(A[t]["KUM"][k160] - B[t][k165]))
            d.setdefault("tau_eff", []).append(
                abs(A[t]["tau_eff"] - B[t]["tau_eff"]))
            d.setdefault("kul", []).append(abs(A[t]["kul"] - B[t]["kul"]))
        print(f"  {v}: {len(ort)} bant | " + "  ".join(
            f"{k}: {max(x):.1e}" for k, x in d.items()))


def v2():
    print("=" * 74)
    print("V2 — Ç1'in ÇIPLAK limiti 163 §2e'nin kapalı formu mu?")

    class Sahte:
        pass
    L = 12.0296
    M = C163.merdiven(L, 0.95)
    tau, a, b, B = M["tau"], M["a"], M["b"], M["B"]
    taban = 0.40
    hp = b * (tau > taban)
    y = B * np.exp(-1j * np.pi * tau)
    Mo = Sahte()
    Mo.hp, Mo.y, Mo.M = hp, y, M
    Mo.qlim = int(np.exp(taban * L))
    print(f"{'q':>8} {'tau_Q':>7} {'s2 Ç1 (kod)':>14} "
          f"{'s2 çıplak (kalem)':>18} {'fark':>10}")
    for tq in (0.54, 0.58, 0.62, 0.66, 0.70, 0.74, 0.78):
        i = int(np.argmin(np.abs(tau - tq)))
        Jc1 = K.Model165.kanal1(Mo, i)[0]
        hQ = complex(b[i])
        rm = abs(hQ) ** 2 / 4
        s2 = K.s_den_J(hQ, Jc1, rm)[0]
        ref = K.c1_capla(tau[i], tau, a, taban, 0.95)
        print(f"{int(M['q'][i]):>8d} {tau[i]:7.4f} {s2:+14.9f} "
              f"{ref:+18.9f} {abs(s2-ref):10.2e}")


def v3(veriler=("Hkeskin", "HA4")):
    print("=" * 74)
    print("V3 — 165'in s-sitesi tayfı ↔ 164 §3g'nin çizgi oranları")
    Q163 = (2, 3, 5, 7, 11, 101, 1009)
    for v in veriler:
        p = S164 / f"sadakat_{v}.json"
        ref = json.loads(p.read_text()) if p.exists() else None
        Mo = K.Model165(v, 0.40, 4000, 0.95, "olculen")
        tau, B = Mo.M["tau"], Mo.M["B"]
        print(f"  {v}: (|y_q|/B_q , arg y_q + πτ_q)  [164: |c_q|/b_q]")
        for q in Q163:
            i = int(np.searchsorted(Mo.M["q"], q))
            r = abs(Mo.y_olc[i]) / abs(B[i])
            ph = float(np.angle(Mo.y_olc[i]) + np.pi * tau[i])
            ph = (ph + np.pi) % TWO_PI - np.pi
            e = ""
            if ref and str(q) in ref.get("cizgi", {}):
                c = ref["cizgi"][str(q)]
                cq = complex(*c["con"])
                e = (f"  [164: |c_q|/b_q={abs(cq)/c['b_pen']:.3f} "
                     f"arg={np.angle(cq):+.4f}]")
            print(f"    q={q:5d} τ={tau[i]:.4f}  |y|/B={r:7.4f}  "
                  f"argy+πτ={ph:+.4f}{e}")


def v4(veri="Hkeskin", tau_s=0.46, nsite=3000):
    print("=" * 74)
    print("V4 — Ç1'in KAPALI toplamı = doğrudan sayım (alt-merdiven)")
    KAN = importlib.import_module("165_kanal")
    Mo = K.Model165(veri, 0.40, 4000, 0.95, "olculen")
    ban = C163.bant_adaylari(Mo.Y, [b for b in K.IZGARA_T1
                                    if b[0] >= 0.52 - 1e-9])
    r = max(ban[0]["cizgi"], key=lambda x: x["gp"])
    iQ = int(np.searchsorted(Mo.M["q"], r["q"]))
    W = float(Mo.M["w"][iQ])
    a, bq, nw, npair = KAN.V4(Mo, iQ, W, tau_s, nsite)
    print(f"  {veri} Q={r['q']} τ_Q={r['tau']:.4f}, alt-merdiven τ≤{tau_s} "
          f"({nw} işaretli çizgi, {npair} çift), {nsite} site:")
    print(f"    SENTEZ = {a:.12e}")
    print(f"    SAYIM  = {bq:.12e}")
    print(f"    bağıl fark = {abs(a-bq)/max(abs(a),1e-300):.3e}")


if __name__ == "__main__":
    v1()
    v2()
    v3()
    v4()
