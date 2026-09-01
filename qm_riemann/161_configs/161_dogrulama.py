"""
161 — DOĞRULAMA: bu raporun okuma zinciri 158/159'la tutarlı mı?

V1  φ fiti  : 161'in fit çekirdeği 158'in W-A tablosunu (τ₀, a, b)
              45 satırda yeniden üretiyor mu?   (158'in json'larından)
V2  ÖZDEŞLİK: aynı ağırlıkla c_φ = c_δ + (0, 4π, −2π) — makine mertebesi?
              ⇒ b(δ) ≡ b(φ) ve a_τ₀ = 4π + dδ/dτ|τ₀ ≡ a(φ)
V3  159 §2c : gerçek gaz, taban 0.40, bant düzeyinde δ tablosu birebir mi?
V4  159 §2d : δ(½), dδ/dτ ve τ₀ = ½ − δ(½)/a — 159'un beş gazlık tablosu
V5  σ_δ vs σ_φ: δ'nın kendi jackknife hatası φ'ninkinden farklı mı?
V6  ÇAPA    : dδ/dτ|½ ile dδ/dτ|τ₀ farkı = 2b(½−τ₀) mü (ölçüldü)?
Çıktı: scratchpad/161/dogrulama_cikti.txt
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from importlib import import_module                                  # noqa: E402
K = import_module("161_cekirdek")

OUT = []


def yaz(*a):
    s = " ".join(str(x) for x in a)
    OUT.append(s)
    print(s, flush=True)


def main():
    D = K.yukle()
    P = {k: K.delta_bant(K.bant(d)) for k, d in D.items() if k[0] != "P0"}

    yaz("=" * 96)
    yaz("V1 — 161'in fit çekirdeği 158'in W-A tablosunu üretiyor mu?")
    yaz("=" * 96)
    yaz("   gaz    taban |   158: τ₀      a        b    |  161: τ₀      a        b"
        "    |   Δτ₀      Δa      Δb")
    dt, da, db = [], [], []
    for (v, t), (r0, ra, rb) in sorted(K.REF158_WA.items()):
        if (v, t) not in P:
            yaz(f"   {v:7s}{t:.2f} | json YOK"); continue
        f = K.olc(P[(v, t)], "A", njack=0)
        if f is None:
            yaz(f"   {v:7s}{t:.2f} | fit YOK"); continue
        dt.append(f["t0_phi"] - r0); da.append(f["a_phi"] - ra); db.append(f["b_phi"] - rb)
        yaz(f"   {v:7s}{t:.2f} |  {r0:.4f} {ra:7.3f} {rb:+8.3f}  |  "
            f"{f['t0_phi']:.4f} {f['a_phi']:7.3f} {f['b_phi']:+8.3f} | "
            f"{f['t0_phi']-r0:+.4f}  {f['a_phi']-ra:+.4f} {f['b_phi']-rb:+.4f}")
    yaz(f"\n   maks |Δτ₀| = {np.max(np.abs(dt)):.5f}   maks |Δa| = {np.max(np.abs(da)):.5f}"
        f"   maks |Δb| = {np.max(np.abs(db)):.5f}    (158'in tablosu 3–4 haneye yuvarlı)")

    yaz("\n" + "=" * 96)
    yaz("V2 — ÖZDEŞLİK: c_φ − c_δ = (0, 4π, −2π)?  ⇒ b(δ)≡b(φ), a_τ₀≡a(φ)")
    yaz("=" * 96)
    mc, mb, ma, n = 0.0, 0.0, 0.0, 0
    for pen in ("A", "D", "B", "C"):
        for (v, t) in sorted(P):
            f = K.olc(P[(v, t)], pen, njack=0)
            if f is None:
                continue
            n += 1
            cd, cp = f["fd"]["c"], f["fphi"]["c"]
            hedef = np.zeros_like(cp); hedef[-2] += K.FOUR_PI; hedef[-1] -= K.TWO_PI
            mc = max(mc, float(np.max(np.abs((cp - cd) - hedef))))
            ma = max(ma, abs(f["a_t0"] - f["a_phi"]))
            if pen != "C":      # kuadratikte b çapadan bağımsız
                mb = max(mb, abs(f["b"] - f["b_phi"]))
    yaz(f"   {n} fit (45 koşu × 4 pencere, w=1/σ_φ; 2'si taban 0.52 + W-B'de kurulamıyor):")
    yaz(f"     maks |Δkatsayı|              = {mc:.3e}")
    yaz(f"     maks |b(δ)|½ − b(φ)|τ₀| (kuadratik) = {mb:.3e}")
    yaz(f"     maks |a_τ₀ − a(φ)|          = {ma:.3e}")
    yaz("   (W-C kübiktir; orada b çapaya bağlıdır ve ½ ile τ₀ okumaları")
    yaz("    farklı sayılardır — o pencere bu satırdan çıkarıldı.)")
    yaz("   ⇒ ÖZDEŞ (aynı ağırlıkla omurga model uzayındadır). b, rotor")
    yaz("     düzeltmesinden ETKİLENMEYEN tek sayıdır: yeniden ölçülmedi, AYNI sayıdır.")

    yaz("\n" + "=" * 96)
    yaz("V3 — 159 §2c bant tablosu (gerçek-son, taban 0.40): δ birebir mi?")
    yaz("=" * 96)
    ref = {0.41: +0.0262, 0.43: +0.0049, 0.45: -0.0198, 0.47: -0.0463,
           0.49: -0.0817, 0.51: -0.1150, 0.53: -0.1592, 0.55: -0.2083,
           0.57: -0.2561, 0.59: -0.3121, 0.61: -0.3887, 0.63: -0.4542}
    p = P[("son", 0.40)]
    yaz("     τ̄     τ_eff     φ       omurga      δ(161)     δ(159)     Δ")
    mx = 0.0
    for i, tb in enumerate(p["tau"]):
        r = ref.get(round(float(tb), 2))
        om = K.omurga(p["x"][i])
        s = (f"   {tb:.2f}  {p['x'][i]:.4f}  {p['phi'][i]:+.4f}  {om:+.4f}  "
             f"{p['d'][i]:+.4f}")
        if r is not None:
            mx = max(mx, abs(p["d"][i] - r)); s += f"   {r:+.4f}   {p['d'][i]-r:+.5f}"
        yaz(s)
    yaz(f"   maks |Δδ| = {mx:.5f}   (159'un tablosu 4 haneye yuvarlı)")

    yaz("\n" + "=" * 96)
    yaz("V4 — 159 §2d: δ(½), dδ/dτ, τ₀ = ½ − δ(½)/a   (taban 0.40, W-A)")
    yaz("=" * 96)
    yaz("   gaz     | δ(½) 161  δ(½) 159    Δ    | dδ/dτ 161  159      Δ    |"
        " τ₀* 161  τ₀ 158     Δ")
    for v in ("son", "orta", "keskin", "A4", "P1"):
        f = K.olc(P[(v, 0.40)], "A", njack=0)
        r1, r2 = K.REF159_DELTA_HALF[v], K.REF159_DDELTA[v]
        r3 = K.REF158_WA[(v, 0.40)][0]
        # 159 dδ/dτ'yı τ₀'da okuyor (a = 4π + dδ/dτ|τ₀)
        dd_t0 = f["a_t0"] - K.FOUR_PI
        yaz(f"   {v:7s} | {f['d_half']:+.4f}   {r1:+.4f}  {f['d_half']-r1:+.4f} | "
            f"{dd_t0:+.3f}   {r2:+.3f}  {dd_t0-r2:+.3f} | "
            f"{f['t0_yildiz']:.4f}  {r3:.4f}  {f['t0_yildiz']-r3:+.5f}")
    yaz("   (δ(½) 161 çapası τ=½; 159 aynı çapayı kullanıyor. dδ/dτ burada")
    yaz("    τ₀ çapasında okundu — 159'un a = 4π + dδ/dτ|τ₀ tanımıyla eşleşsin diye.)")

    yaz("\n" + "=" * 96)
    yaz("V5 — σ_δ (δ'nın kendi jackknife'ı) vs σ_φ")
    yaz("=" * 96)
    yaz("   gaz    taban | bant | medyan σ_φ | medyan σ_δ | medyan σ_δ/σ_φ | menzil")
    for v in K.MERDIVEN:
        for t in (0.40,):
            p = P[(v, t)]
            m = K.secim(p, "A")
            r = p["s_d"][m] / p["s_phi"][m]
            yaz(f"   {v:7s}{t:.2f} |  {int(m.sum()):2d}  | {np.median(p['s_phi'][m]):.5f}    | "
                f"{np.median(p['s_d'][m]):.5f}    |     {np.median(r):.3f}      | "
                f"{r.min():.3f} – {r.max():.3f}")

    yaz("\n" + "=" * 96)
    yaz("V6 — ÇAPA farkı: dδ/dτ|½ − dδ/dτ|τ₀ = 2b(½ − τ₀)?")
    yaz("=" * 96)
    yaz("   gaz    taban |  dδ/dτ|½   dδ/dτ|τ₀    fark    2b(½−τ₀)    Δ")
    mx = 0.0
    for v in K.MERDIVEN:
        for t in (0.40,):
            f = K.olc(P[(v, t)], "A", njack=0)
            fark = f["dd"] - (f["a_t0"] - K.FOUR_PI)
            pred = 2 * f["b"] * (0.5 - f["t0_phi"])
            mx = max(mx, abs(fark - pred))
            yaz(f"   {v:7s}{t:.2f} | {f['dd']:+8.4f}  {f['a_t0']-K.FOUR_PI:+8.4f}  "
                f"{fark:+8.4f}  {pred:+8.4f}  {fark-pred:+.2e}")
    yaz(f"   maks |Δ| = {mx:.3e}  (ÖZDEŞ — kuadratikte tam)")

    (K.SCR161 / "dogrulama_cikti.txt").write_text("\n".join(OUT))
    yaz(f"\n[yazıldı] {K.SCR161/'dogrulama_cikti.txt'}")


if __name__ == "__main__":
    main()
