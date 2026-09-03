# -*- coding: utf-8 -*-
"""
171a — ÇARPANLAŞMA SINAMASI: KALİB_u2(λ,τ) = A(τ)·M(λ) ?  (3 Eylül 2026, gece-2)
Kaptan kalemi (T1). Veri: 170'in K0.json'u — YENİ KOŞU YOK.

ÖN-MÜHÜR (KALEM_NU_UYE_03EYL2026.md):
 Ö1  λ-serisi 5 gazda (L060..L115) KALİB bant-şekli kolapsı: ilk üç bant ≤%0.5,
     tüm bantlar rms ≤%2, en kötü (son bant) %3-4.
 Ö2  M(λ) (orta banttan) ≈ 1.00 / 1.00 / 1.01 / 1.06 / 1.12  (λ = 1.15→0.60).
 Ö3  ν_band·σ_X̃² yaklaşık sabit (%10-15 içinde) — A(τ) sabit-şekil sonucu.
 Ö4  Gerçek gaz (son) aynı A(τ) üstünde (şekil sapması kolaps bandı içinde);
     farkı yalnız M'de (+%5-6).

KALİB_b yeniden kurulumu: KALIB_b = c_bant_b · W_X_b,
W_X_b = exp(−(2πτ_b)²σ_X̃²/2)  (orta bantta JSON'un W_X'iyle %0.12-0.56 içinde).

SONUÇ (gerçek koşudan, 3 Eyl gece-2):
 Ö1 ✓  KOLAPS MÜHÜRLÜ: rms %0.81, en kötü %2.76 (L115 son bandı); ilk üç bant
       ≤%0.33. A(τ) = 1.1048 / 1.0448 / 1.0000 / 0.9423 / 0.8832.
 Ö2 ✓  M(λ) = 1.0002 / 1.0000 / 1.0169 / 1.0647 / 1.1272 (λ = 1.15→0.60).
 Ö3 ~  ν_band·σ_X̃² = 0.0498–0.0616 (±%13, sistematik eğilimli) — ölçek doğru,
       saf 1/σ² değil (A·M çerçevesinde beklenen ikinci-mertebe düzeltmeler).
 Ö4 ✓✓ GERÇEK GAZ AYNI A(τ) ÜSTÜNDE (sapmalar +0.87/−0.25/0.00/−0.02/−0.56% —
       sentetik kolaps bandının içinde); M(son) = 1.0559 ⇒ +3σ fazlası tamamen
       SEVİYE farkı (ΔM). c-oranı sağlaması: 1.0559·(W_X_Hk/W_X_son) = 1.0219
       ↔ ölçülen 0.4122/0.4035 = 1.0216 ✓.
 ν_λ yeniden üretimi: −0.004 / +0.247 / +0.607 / +1.024 (170 ile birebir).
 M koklaması: M−1 varyans-açığında doğrusal DEĞİL; kare-aday oranı 0.62–0.72
 bandında ama λ>1 tarafında ölüyor (M TEK-TARAFLI doymuş: λ≥1'de düz) —
 kimlik yarışı 171'in işi.
"""
import json, math

Y = "/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/170/K0.json"
d = json.load(open(Y))["gaz"]

LAMBDA_SIRA = ["L115", "Hkeskin", "L085", "L070", "L060"]
def wx(tau, sig): return math.exp(-((2*math.pi*tau)**2) * sig*sig / 2.0)

def kalib_bant(g):
    r = d[g]
    return [c * wx(t, r["sigX"]) for c, t in zip(r["c_bant"], r["tau_eff"])], r

print("== W_X yeniden-kurulum kontrolü (orta bant, JSON'un W_X'i ile) ==")
for g in LAMBDA_SIRA + ["son"]:
    r = d[g]
    w3 = wx(r["tau_eff"][2], r["sigX"])
    print("  %-8s W_X(json) = %.6f   parametrik = %.6f   fark %+.3f%%"
          % (g, r["W_X"], w3, 100*(w3/r["W_X"]-1)))

print()
print("== T1: ŞEKİL KOLAPSI — S_b = KALIB_b / KALIB_orta ==")
S = {}
K3mid = {}
for g in LAMBDA_SIRA + ["son"]:
    kb, r = kalib_bant(g)
    K3mid[g] = kb[2]
    S[g] = [k / kb[2] for k in kb]
    print("  %-8s lam=%-5s  S = %s" % (g, str(r["lam"]) if r["lam"] else "-",
          "  ".join("%.4f" % s for s in S[g])))

A = [sum(S[g][b] for g in LAMBDA_SIRA)/5 for b in range(5)]
print("  A(τ) (5-gaz ort)   :", "  ".join("%.4f" % a for a in A))
print("  bant-bant kolaps sapmaları (gaz - A, %):")
en_kotu = 0.0; rms_top = 0.0; n = 0
for g in LAMBDA_SIRA:
    sap = [100*(S[g][b]/A[b]-1) for b in range(5)]
    for s_ in sap: rms_top += s_*s_; n += 1; en_kotu = max(en_kotu, abs(s_))
    print("    %-8s %s" % (g, "  ".join("%+.2f" % s_ for s_ in sap)))
print("  KOLAPS rms = %.2f%%   en kötü = %.2f%%" % (math.sqrt(rms_top/n), en_kotu))
sap_son = [100*(S["son"][b]/A[b]-1) for b in range(5)]
print("  GERÇEK (son) sapması: %s   (Ö4 hakemi)" % "  ".join("%+.2f" % s_ for s_ in sap_son))

print()
print("== M(λ) — orta-bant seviyesi (Hkeskin = 1) ==")
for g in LAMBDA_SIRA + ["son"]:
    print("  %-8s lam=%-5s  M = %.4f   (sigX = %.5f, sigds = %.5f)"
          % (g, str(d[g]["lam"]) if d[g]["lam"] else "-",
             K3mid[g]/K3mid["Hkeskin"], d[g]["sigX"], d[g]["sigds"]))

print()
print("== ν_band ve 1/σ_X̃² ölçeği (Ö3) ==")
for g in LAMBDA_SIRA + ["son"]:
    kb, r = kalib_bant(g)
    ws = [wx(t, r["sigX"]) for t in r["tau_eff"]]
    x = [math.log(w) for w in ws]; y = [math.log(k) for k in kb]
    mx = sum(x)/5; my = sum(y)/5
    nu = sum((a-mx)*(b-my) for a, b in zip(x, y)) / sum((a-mx)**2 for a in x)
    print("  %-8s ν_band = %.3f   ν_band·σ_X̃² = %.5f" % (g, nu, nu*r["sigX"]**2))

print()
print("== ν_λ sağlaması (ardışık λ çiftlerinden, orta bant) ==")
for g1, g2 in zip(LAMBDA_SIRA[:-1], LAMBDA_SIRA[1:]):
    w1 = wx(d[g1]["tau_eff"][2], d[g1]["sigX"]); w2 = wx(d[g2]["tau_eff"][2], d[g2]["sigX"])
    nu_l = (math.log(K3mid[g2]) - math.log(K3mid[g1])) / (math.log(w2) - math.log(w1))
    print("  %s → %s :  ν_λ = %+.3f" % (g1, g2, nu_l))

print()
print("== M(λ) aday koklamaları (kaba) ==")
s1 = d["Hkeskin"]["sigX"]
for g in LAMBDA_SIRA:
    r = d[g]; M = K3mid[g]/K3mid["Hkeskin"]
    if g != "Hkeskin":
        print("  %-8s M−1 = %+.4f   (σ1²−σ²)/σ1² = %+.4f   oran = %+.3f"
              % (g, M-1, (s1**2 - r["sigX"]**2)/s1**2,
                 (M-1)/((s1**2 - r["sigX"]**2)/s1**2) if abs(s1**2-r["sigX"]**2) > 1e-12 else float('nan')))
