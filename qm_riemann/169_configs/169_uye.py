"""
169 — ÜYE SINAVI: hangi DW çarpanına bölünmeli? İKİ ÖLÇÜT ÇATIŞIYOR MU?
=======================================================================
168 §A2.8 üyeyi TEK ölçütle seçti: `log(KALİB/üye)`'nin BANT (τ_eff)
eğimi en küçük olan üye. 169 ikinci bir ölçüt ekliyor: λ ekseninde
DEĞİŞMEZLİK (Hkeskin / L085 / L070 — üçü de aynı merdiven, yalnız genlik
ölçeği farklı; fizik aynı olmalı).

  ölçüt (A) τ-eğimi   : |d log c_m / d τ_eff|  (Hkeskin, dokuz bant)
  ölçüt (B) λ-saçılımı: log c_m'in {λ=1.00, 0.85, 0.70} arası sd'si
  ölçüt (C) mutlak    : c_m(Hkeskin) — 4/π² ile karşılaştırılır

Çıktı: ekrana tablo + scratchpad/169/UYE.json
"""
import json
import os

import numpy as np

SCR = ("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
       "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
LO_MIN, LO_MAX, SNR_MIN, R_MIN = 0.52, 0.68, 3.0, 0.98
C_HIP = 4.0 / np.pi ** 2

UYE = {
    "K0 (sabit)":        lambda b: 1.0,
    "W_amp":             lambda b: b["W_amp"],
    "√W_X":              lambda b: np.sqrt(b["W_X"]),
    "W_X":               lambda b: b["W_X"],
    "√(W_amp·W_X)":      lambda b: np.sqrt(b["W_amp"] * b["W_X"]),
    "W_pos":             lambda b: b["W_pos"],
    "W_amp·W_X":         lambda b: b["W_amp"] * b["W_X"],
    "W_X^{3/2}":         lambda b: b["W_X"] ** 1.5,
}


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


def main():
    D = {a: json.load(open(f"{SCR}/167/C_{a}.json"))
         for a in ("Hkeskin", "L085", "L070")}
    B9 = bantlar(D["Hkeskin"], 0.80)
    print("=" * 96)
    print("ÜYE SINAVI — iki ölçüt (A: τ-eğimi, 168'in ölçütü; B: λ-değişmezlik,"
          " 169'un ölçütü)")
    print("=" * 96)
    print(f"  (A) Hkeskin, {len(B9)} bant (lo≤0.80);  (B) ortak pencere "
          f"lo∈[0.52,0.68], λ = 1.00 / 0.85 / 0.70")
    print("\n  üye              τ-eğimi   artık_rms |  c(λ=1.00) c(0.85) "
          " c(0.70)  **λ-sd**  λ-menzil%  | 4/π²'ye %")
    R = {}
    for ad, f in UYE.items():
        lg = np.log([b["KALIB_u2"] / f(b) for b in B9])
        te = np.array([b["tau_eff"] for b in B9])
        p = np.polyfit(te, lg, 1)
        egim = float(p[0])
        rms = float(np.std(lg - np.polyval(p, te)))
        cl = []
        for g in ("Hkeskin", "L085", "L070"):
            Bb = bantlar(D[g])
            cl.append(float(np.exp(np.mean(np.log(
                [b["KALIB_u2"] / f(b) for b in Bb])))))
        sd = float(np.std(np.log(cl), ddof=1))
        men = 100 * (max(cl) / min(cl) - 1)
        R[ad] = dict(egim=egim, rms=rms, c=cl, lam_sd=sd, lam_menzil=men)
        print(f"  {ad:16s} {egim:+7.3f}   {rms:.4f}   | {cl[0]:.4f}   "
              f"{cl[1]:.4f}   {cl[2]:.4f}   {sd:.4f}    {men:6.2f}    | "
              f"{100*(cl[0]/C_HIP-1):+7.2f}")
    print(f"\n  4/π² = {C_HIP:.4f}")
    print("\n  ÇATIŞMA: (A) en küçük |τ-eğimi| olan üye vs (B) en küçük "
          "λ-sd olan üye")
    a_best = min(R, key=lambda k: abs(R[k]["egim"]))
    b_best = min(R, key=lambda k: R[k]["lam_sd"])
    print(f"    (A) seçer: **{a_best}**  (eğim {R[a_best]['egim']:+.3f}, "
          f"λ-menzil {R[a_best]['lam_menzil']:.2f}%)")
    print(f"    (B) seçer: **{b_best}**  (eğim {R[b_best]['egim']:+.3f}, "
          f"λ-menzil {R[b_best]['lam_menzil']:.2f}%)")
    json.dump(R, open(f"{SCR}/169/UYE.json", "w"), indent=1, default=float)
    print(f"\n-> {SCR}/169/UYE.json")


if __name__ == "__main__":
    main()
