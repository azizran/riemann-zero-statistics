# -*- coding: utf-8 -*-
"""
175g — κ'NIN DÜZELTİLMİŞ KAPALI YASASI — **ÖN-MÜHÜRSÜZ İKİNCİ TUR**
=====================================================================
DÜRÜSTLÜK: bu betik `175e`'nin ÖLÇÜMÜNDEN SONRA yazıldı ve öyle
etiketlenir (174e'nin aynı konumu gibi). 175a'nın ön-mühürlü K3
maddeleri DEĞİŞMEZ: K3-a öldü, K3-c'nin rakibi kazandı.

175e iki şey ölçtü:
  (i) ASAL çizgilerde `|κ| / |−πA_qτ_q cos(πτ_q)|` medyanı τ = 0.60…0.95
      boyunca **0.90 – 1.01** — yani `W_pos` sönümü YOKTUR;
  (ii) KULE çizgilerinde aynı oran **×2.0 – 2.3** — ve kulelerin medyan
      çokluğu tam **2.0**.

(ii)'nin kaynağı, seviye denkleminin İKİNCİ MERTEBE tersidir:
    N̄(z_n) + S(z_n) = n − ½ ,   N̄(z̄_n) = n − ½ ,  d := z_n − z̄_n
    N̄'d = −S − S'd + O(d³)  ⇒  d = −ḡS + ḡ²S S' + …     (ḡ := 1/N̄')
`S S'` ÇARPIMSAL frekansları doğurur (143'ün G-yasası köşesi):
    S S' = ½ Σ_{c₁c₂} A₁A₂ω₂ [ sin((ω₁+ω₂)t) + sin((ω₁−ω₂)t) ]
Yer değiştirmenin ω_Q bileşeni (ḡ ile bölünmüş, yani ETKİN genlik):

    A_Q^eff = A_Q + (πτ_Q/2)·Σ_{q₁q₂=Q} A₁A₂  −  πτ_Q·Σ_{q₁/q₂=Q} A₁A₂
              ────────────  ─────────────────      ─────────────────
               1. mertebe     ÇARPIMSAL (kule)       FARK kanalı
    (toplamlar SIRALI merdiven çiftleri üzerinde; ḡω_Q = 2πτ_Q)

    ***  κ(ω_Q) ≈ −π τ_Q cos(πτ_Q) · A_Q^eff  ***       (175g yasası)

SIFIR PARAMETRE. Kapalı kehanet: `Q = p²` için ikinci terim birincinin
tam **τ_Q** katıdır (A_p²/A_{p²} = 2/π özdeşliğinden), `Q = p³` için
**1.5 τ_Q**. Yani kulelerin ×2 fazlası bir uydurma değil, ARİTMETİĞİN
kendisidir: **`κ`'nın kule fazlası ÇARPIMSAL kanaldır.**

SINAV: `|κ| / |κ̂₂|` medyanı ASAL ve KULE çizgilerinde AYNI ve 1'e yakın
olmalı. (175e'de asallar 0.90–1.01, kuleler 1.9–2.3 idi.)
"""
import importlib
import math
import json
import sys
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
for _p in ("167_configs", "165_configs", "163_configs", "160_configs",
           "159_configs", "156_configs", "155_configs", "154_configs"):
    sys.path.insert(0, str(QM / _p))
K165 = importlib.import_module("165_cekirdek")
ORT = importlib.import_module("167_ortak")
K165.PENCERE.update(ORT.pencere_dict())

SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S175, S174, S172 = SCR / "175", SCR / "174", SCR / "172"
G1 = json.load(open(S172 / "G1.json"))
GAZ = sys.argv[1:] or ["son", "Hkeskin", "K090", "K070"]
INCE = np.arange(0.40, 0.96, 0.05)

print("=" * 78)
print("175g — κ'nın düzeltilmiş kapalı yasası  (**ÖN-MÜHÜRSÜZ İKİNCİ TUR**)")
print("=" * 78)
OUT = {}
for g in GAZ:
    p = S174 / f"K3_{g}.json"
    if not p.exists():
        print(f"  ({g}: K3 yok)")
        continue
    d = json.load(open(p))
    lam, tau_ust = d["lam"], d["tau_ust"]
    q = np.array(d["q"], dtype=np.int64)
    tau = np.array(d["tau"])
    kap = np.array(d["kap_abs"])
    L = float(np.log(q[-1].astype(float)) / tau[-1])
    M = K165.merdiven(L, 0.95, g)
    qm = {int(x): i for i, x in enumerate(M["q"])}
    idx = np.array([qm[int(x)] for x in q])
    mask = (M["tau"] <= tau_ust + 1e-12).astype(float)
    Aall = lam * M["a"] * mask
    tall = M["tau"]
    pos = {int(x): i for i, x in enumerate(M["q"])}
    mult = M["mult"][idx].astype(int)
    A = Aall[idx]

    Smul = np.zeros(len(q))          # Σ_{q₁q₂=Q} A₁A₂  (sıralı)
    Dfar = np.zeros(len(q))          # Σ_{q₁/q₂=Q} A₁A₂ (sıralı, tek yön)
    nS = nD = 0
    for j, Q in enumerate(q):
        k = int(mult[j])
        P = int(round(float(Q) ** (1.0 / k)))
        while P ** k != Q:
            P += 1 if P ** k < Q else -1
        s = 0.0
        for i in range(1, k):                     # (p^i , p^{k−i})
            a1, a2 = pos.get(P ** i), pos.get(P ** (k - i))
            if a1 is not None and a2 is not None:
                s += Aall[a1] * Aall[a2]
        t = 0.0
        b = 1
        while True:                               # (p^{b+k} , p^b)
            a1, a2 = pos.get(P ** (b + k)), pos.get(P ** b)
            if a1 is None or a2 is None:
                break
            t += Aall[a1] * Aall[a2]
            b += 1
        Smul[j], Dfar[j] = s, t
        nS += s != 0
        nD += t != 0
    A_eff = A + (np.pi * tau / 2.0) * Smul - np.pi * tau * Dfar
    kh0 = np.abs(-np.pi * A * tau * np.cos(np.pi * tau))
    kh2 = np.abs(-np.pi * A_eff * tau * np.cos(np.pi * tau))
    r0 = kap / np.maximum(kh0, 1e-300)
    r2 = kap / np.maximum(kh2, 1e-300)
    kule = mult >= 2
    print("\n  %s:   çizgi %d (asal %d, kule %d) ; Σ_çarpımsal ≠ 0: %d ; "
          "Σ_fark ≠ 0: %d" % (g, len(q), int((~kule).sum()), int(kule.sum()),
                              nS, nD))
    print("     τ bandı      n  n_kule |  ASAL: |κ|/|κ̂₀|  |κ|/|κ̂₂| | "
          " KULE: |κ|/|κ̂₀|  |κ|/|κ̂₂| |  A^eff/A (kule)  öngörü 1+c_kτ")
    sat = []
    # yasa yalnız gazın MERDİVENİNDE gerçekten çizgi olan yerde sınanabilir:
    # keskin kesimde τ > τ_ust'ta A_q ≡ 0, erfc'de zarf çizgiyi silmiştir —
    # orada ölçülen κ bir TABANDIR, bir tepki değil (oran tanımsız).
    pen = ORT.KUNYE.get(g, {}).get("pen")
    wq = (np.array([0.5 * math.erfc((float(t) - pen[0]) / pen[1])
                    for t in tau]) if pen else np.ones_like(tau))
    wq = wq * (tau <= tau_ust + 1e-12)
    var = wq >= 0.02
    for lo, hi in zip(INCE[:-1], INCE[1:]):
        m = (tau > lo) & (tau <= hi) & var
        if m.sum() < 3:
            continue
        ma, mk = m & ~kule, m & kule
        ck = np.where(mult[mk] == 2, 1.0, 1.5) * tau[mk] if mk.sum() else None
        rec = dict(lo=float(lo), hi=float(hi), n=int(m.sum()),
                   nk=int(mk.sum()),
                   a0=float(np.median(r0[ma])), a2=float(np.median(r2[ma])),
                   k0=float(np.median(r0[mk])) if mk.sum() else np.nan,
                   k2=float(np.median(r2[mk])) if mk.sum() else np.nan,
                   aeff=float(np.median((A_eff / np.maximum(A, 1e-300))[mk]))
                   if mk.sum() else np.nan,
                   ong=float(np.median(1 + ck)) if mk.sum() else np.nan)
        sat.append(rec)
        print("     %.2f–%.2f %6d %6d | %14.4f %11.4f | %14s %11s | %11s %11s"
              % (lo, hi, rec["n"], rec["nk"], rec["a0"], rec["a2"],
                 "%.4f" % rec["k0"] if rec["k0"] == rec["k0"] else "—",
                 "%.4f" % rec["k2"] if rec["k2"] == rec["k2"] else "—",
                 "%.4f" % rec["aeff"] if rec["aeff"] == rec["aeff"] else "—",
                 "%.4f" % rec["ong"] if rec["ong"] == rec["ong"] else "—"))
    ust = [r for r in sat if r["lo"] >= 0.60]
    ma0 = float(np.median([r["a0"] for r in ust]))
    mk0 = float(np.median([r["k0"] for r in ust if r["k0"] == r["k0"]]))
    mk2 = float(np.median([r["k2"] for r in ust if r["k2"] == r["k2"]]))
    print("     τ>0.60 medyanları:  ASAL |κ|/|κ̂₀| = %.4f  |  KULE "
          "|κ|/|κ̂₀| = %.4f → |κ|/|κ̂₂| = %.4f   (kule fazlası %s)"
          % (ma0, mk0, mk2,
             "**KAPANDI**" if abs(mk2 - 1) < abs(mk0 - 1) * 0.5 else
             "kısmen kapandı" if abs(mk2 - 1) < abs(mk0 - 1) else "KAPANMADI"))
    OUT[g] = dict(bant=sat, asal0=ma0, kule0=mk0, kule2=mk2,
                  nS=int(nS), nD=int(nD))
json.dump(OUT, open(S175 / "K3G.json", "w"), indent=1, ensure_ascii=False)
print("\n-> %s" % (S175 / "K3G.json"))
