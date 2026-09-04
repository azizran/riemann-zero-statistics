# -*- coding: utf-8 -*-
"""
172a — VADİ AYRIŞTIRMASI: M(λ)'yı çarpanlarına böl (4 Eylül 2026, sabah kalemi)
Veri: 167/C_*.json önbellekleri — YENİ KOŞU YOK. M = (g_E/g₀)·(g_X/g₀)²·(θ/θ₀),
θ_b := KALİB_u2,b / g_cal (tanım gereği özdeşlik; içerik, çarpanların λ-eğrilerinde).

ÖN-MÜHÜR (koşudan önce):
 Ö1  (g_X/g₀)² λ boyunca ~düz (171: g_X λ-değişmez): U'ya katkısı ≤ %1-2.
 Ö2  M'nin yükselen kolunu g_E ve θ birlikte taşır (~eşit paylı, ΔM(son)'un
     %45/%52'si gibi); λ≥1 düzlüğünü hangisinin taşıdığı AÇIK — ölçüm söylesin.
 Ö3  nline bütün λ-gazlarında sabit (8981) ⇒ θ'nin λ-hareketi sayım değil
     ağırlık/rezonans-şiddeti işi.
 Ö4  g_E/g₀, korE/korE₀ ile birlikte hareket eder (çizgi-uyum sızıntı ailesi);
     gerçek gazın (son) payları 171-T3 ile tutarlı çıkar (g_E %45, θ %52).

SONUÇ (gerçek koşudan, 4 Eyl sabahı):
 Ö1 ~ KISMİ ISKA: (gX/g₀)² çekirdekte düz (λ 0.6-1.3: ±%0.7) ama 0.60→0.50
      adımında uyanıyor (Δlog +0.0157 = adımın %17'si; öngörü ≤%1-2 idi).
 Ö2 ✓ Yükselen kolu gE (+0.0417) ile θ (+0.0330) birlikte taşıyor; gerçek gazın
      payları gE %43 / gX² %3 / θ %55 (171-T3 ile uyumlu).
 Ö3 ✓ nline = 8981 bütün gazlarda sabit ⇒ θ'nin λ-hareketi sayım değil AĞIRLIK işi.
 Ö4 ✗/✓ gE ∝ korE^p ÖLDÜ (log-log oran 2.85→−0.22, işaret dönüyor); son'un
      pay dağılımı tutarlı çıktı.
 YENİDEN ÇERÇEVE: θ/θ₀ TÜM merdivende λ'da TEKDÜZE (1.0936→0.9668) — tek işaretli
 "doyum bacağı" adayı; gE tek-taraflı hokey sopası (λ≥1 düz, altta yükseliş);
 M'nin [1.00,1.15] platosu derin yasa değil, θ↓ ile gX²↑'nin yaklaşık iptali.
 KOKLAMA (ön-mühürsüz, etiketli): gE, varE_res/varE_olc oranında DOĞRUSAL —
 λ ailesi + GERÇEK gaz aynı doğrunun üstünde (%0.1-0.3; kaba gE ≈ 0.798−0.373·rE);
 L130 bilinen kırık bölgesinde sapıyor. θ(son) ise her marjinal eğrinin DIŞINDA
 (+%3) — aritmetik fazlanın dar adresi θ olabilir. 172 kapıları bu iki kokudan.
"""
import json, math

KOK = "/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/167/"
SIRA = ["L050", "L060", "L070", "L085", "Hkeskin", "L115", "L130"]
LAM = {"L050": 0.50, "L060": 0.60, "L070": 0.70, "L085": 0.85, "Hkeskin": 1.00, "L115": 1.15, "L130": 1.30}

def oku(g):
    d = json.load(open(KOK + "C_%s.json" % g))
    a = d["artik"]
    orta = [b for b in d["bant"] if abs(b["lo"] - 0.60) < 1e-9][0]
    return {"gE": a["gE"], "gX": a["gX"], "gcal": a["gcal"], "korE": a["korE"],
            "korX": a["korX"], "nline": a["nline"], "sigX": d["sigX"],
            "K": orta["KALIB_u2"], "WX": orta["W_X"], "th": orta["KALIB_u2"] / a["gcal"]}

V = {g: oku(g) for g in SIRA + ["son"]}
H = V["Hkeskin"]

print("== ÇARPAN EĞRİLERİ (orta bant lo=0.60; Hkeskin = 1) ==")
print("gaz       λ     M       gE/g0   (gX/g0)²  θ/θ0    | korE/k0  korX/k0  nline")
for g in SIRA + ["son"]:
    v = V[g]
    M = v["K"] / H["K"]; fE = v["gE"] / H["gE"]; fX = (v["gX"] / H["gX"]) ** 2; ft = v["th"] / H["th"]
    print("%-8s %-5s %.4f  %.4f  %.4f   %.4f  | %.4f   %.4f   %d"
          % (g, ("%.2f" % LAM[g]) if g in LAM else "-", M, fE, fX, ft,
             v["korE"] / H["korE"], v["korX"] / H["korX"], v["nline"]))

print()
print("== LOG-EĞİMLER (ardışık λ çiftleri; dlog/dλ) ==")
print("adım            dlogM    dlog gE   dlog gX²  dlog θ    dlog W_X")
for g1, g2 in zip(SIRA[:-1], SIRA[1:]):
    v1, v2 = V[g1], V[g2]; dl = LAM[g2] - LAM[g1]
    def s(a, b): return (math.log(b) - math.log(a)) / dl
    print("%-5s→%-8s %+7.3f  %+7.3f   %+7.3f   %+7.3f   %+7.3f"
          % (g1, g2, s(v1["K"], v2["K"]), s(v1["gE"], v2["gE"]),
             s(v1["gX"] ** 2, v2["gX"] ** 2), s(v1["th"], v2["th"]), s(v1["WX"], v2["WX"])))

print()
print("== GERÇEK GAZIN ΔM PAYLARI (log-pay, Ö4 sağlaması) ==")
v = V["son"]
dM = math.log(v["K"] / H["K"]); dE = math.log(v["gE"] / H["gE"])
dX = math.log((v["gX"] / H["gX"]) ** 2); dt = math.log(v["th"] / H["th"])
print("ΔlogM = %+.5f  |  gE: %+.5f (%.0f%%)  gX²: %+.5f (%.0f%%)  θ: %+.5f (%.0f%%)"
      % (dM, dE, 100 * dE / dM, dX, 100 * dX / dM, dt, 100 * dt / dM))

print()
print("== TEK-TARAFLILIK ve VADİ SORUSU ==")
print("λ≥1 bölgesi (1.00→1.30): Δlog gE = %+.4f, Δlog θ = %+.4f, Δlog gX² = %+.4f"
      % (math.log(V["L130"]["gE"] / H["gE"]), math.log(V["L130"]["th"] / H["th"]),
         math.log((V["L130"]["gX"] / H["gX"]) ** 2)))
print("λ≤0.65 kolu (0.60→0.50): Δlog gE = %+.4f, Δlog θ = %+.4f, Δlog gX² = %+.4f"
      % (math.log(V["L050"]["gE"] / V["L060"]["gE"]), math.log(V["L050"]["th"] / V["L060"]["th"]),
         math.log((V["L050"]["gX"] / V["L060"]["gX"]) ** 2)))

print()
print("== Ö4: gE — korE bağlantısı (log-log) ==")
for g in SIRA + ["son"]:
    v = V[g]
    le = math.log(v["gE"] / H["gE"]); lk = math.log(v["korE"] / H["korE"])
    print("  %-8s log(gE/g0) = %+.4f   log(korE/k0) = %+.4f   oran = %s"
          % (g, le, lk, ("%+.2f" % (le / lk)) if abs(lk) > 1e-6 else "  —"))
