# -*- coding: utf-8 -*-
"""
175e — K3: ÇARPIMSAL κ (H-K3) ve ÖN-KAYITLI RAKİBİ
====================================================
174e ölçtü: `κ`'nın kapalı yasası `κ̂ = −πA_qτ_q cos(πτ_q) W_pos(τ_q)`
τ ≤ 0.50'de ±%20 tutuyor, τ > 0.50'de ölçülen `|κ|` yasayı ×1.9 … ×4.2
aşıyor ve fazla τ ile büyüyor. 174e fazlanın adayını koydu ama saymadı:
**çarpımsal çözümler** (`q₁q₂ = Q`). 175 o borcu ödüyor.

TÜRETİM (175a ön-kaydında, ölçümden önce):
    κ⁽²⁾ = (π²τ_Q²/2)·(S_Q − D_Q)
    S_Q = Σ_{q₁q₂=Q}  A₁A₂cos(πτ₁)cos(πτ₂)      (ÇARPIMSAL)
    D_Q = Σ_{q₁/q₂=Q} A₁A₂cos(πτ₁)cos(πτ₂)      (FARK; iki sıralama)
Merdiven asal kuvvetlerdir ⇒ `q₁q₂ = Q`'nun çözümü yalnız `Q = p^k`,
`k ≥ 2` (KULE); `q₁/q₂ = Q` ise `Q = p^k`'nın hepsinde (asallar dahil)
ama yalnız `p^{b+k} ≤ e^{0.95L}` olan küçük p'lerde.

KOD NOTU (düzeltme): 174e `A_q`'yu `163_cekirdek.merdiven`'den alıyordu —
o merdiven erfc penceresini UYGULAMAZ. Erfc gazları (HA4/E060) için
175e `165_cekirdek.merdiven(L, 0.95, veri)`yi kullanır (pencere uygulanır);
keskin gazlarda ikisi özdeştir (denetim Ö-e0 raporlanır).

ÖN-MÜHÜRLER 175a'da donduruldu (P7):
  K3-a  medyan(|κ|/|κ̂|) ASAL çizgilerde τ>0.60'ta ≥ 1.5 ⇒ H-K3 ÖLÜR.
  K3-b  KULE çizgilerinde κ⁽²⁾ eklenince |1 − oran| KÜÇÜLMELİ.
  K3-c  RAKİP: κ̂₀ = −πA_Qτ_Q cos(πτ_Q) (W_pos YOK). τ>0.55'te medyan
        |κ|/|κ̂₀| ∈ [0.8, 1.25] ⇒ rakip kazanır, H-K3 gereksizdir.
"""
import importlib
import json
import sys
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
for _p in ("167_configs", "165_configs", "163_configs", "160_configs",
           "159_configs", "156_configs", "155_configs", "154_configs"):
    sys.path.insert(0, str(QM / _p))
K165 = importlib.import_module("165_cekirdek")
C163 = importlib.import_module("163_cekirdek")
ORT = importlib.import_module("167_ortak")
K165.PENCERE.update(ORT.pencere_dict())

SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S175, S174 = SCR / "175", SCR / "174"
ONK = json.load(open(S175 / "ONKAYIT_K2.json"))
TWO_PI = 2 * np.pi
GAZ = sys.argv[1:] or ["son", "Hkeskin", "HA4"]
INCE = np.arange(0.40, 0.96, 0.05)

print("=" * 74)
print("175e — K3: ÇARPIMSAL κ   [ön-kayıt %s  sha %s]"
      % (ONK["zaman"], ONK["sha256"][:16]))
print("=" * 74)


def w_pos_egri(veri, tgrid):
    """W_pos(τ) = ⟨e^{2πiτ Ĉ_n}⟩ , Ĉ = cumsum(X̃0) kübik trendsiz
    (174e / 165.koherans / 167_olcum ile AYNI tanım)."""
    Y = K165.gaz(veri, 0.40, 4000)
    C = np.cumsum(Y.Xtil0)
    n = np.arange(len(C), dtype=float)
    n = (n - n.mean()) / n[-1]
    V = np.vstack([np.ones_like(n), n, n * n, n ** 3]).T
    c, *_ = np.linalg.lstsq(V, C, rcond=None)
    Cd = C - V @ c
    return (np.array([np.mean(np.exp(1j * TWO_PI * t * Cd)) for t in tgrid]),
            float(np.std(Cd)), Y)


tg = np.linspace(0.0, 0.96, 241)
OUT = {}
for g in GAZ:
    p = S174 / f"K3_{g}.json"
    if not p.exists():
        print(f"  ({g}: K3 yok, atlanıyor)")
        continue
    d = json.load(open(p))
    lam, tau_ust = d["lam"], d["tau_ust"]
    W, sC, Y = w_pos_egri(g, tg)
    q = np.array(d["q"], dtype=np.int64)
    tau = np.array(d["tau"])
    kap = np.array(d["kap_abs"])
    dphi = np.array(d["dphi"])
    hp = np.array(d["hp_abs"])
    xi = np.array(d["xi"])
    L = float(Y.L)

    # --- merdiven: PENCERELİ (165) ve pencereSİZ (163) ------------------
    M5 = K165.merdiven(L, 0.95, g)          # erfc penceresi uygulanır
    M3 = C163.merdiven(L, 0.95)             # 174e'nin kullandığı (pencereSİZ)
    qm = {int(x): i for i, x in enumerate(M5["q"])}
    idx = np.array([qm[int(x)] for x in q])
    mask = (M5["tau"] <= tau_ust + 1e-12).astype(float)
    A_tam = lam * M5["a"] * mask            # bütün merdiven üzerinde A_q
    A = A_tam[idx]
    A174 = (lam * M3["a"] * mask)[idx]
    e0 = float(np.max(np.abs(A - A174)))
    mult = M5["mult"][idx].astype(int)
    Wq = np.interp(tau, tg, W.real) + 1j * np.interp(tau, tg, W.imag)
    C_ = np.cos(np.pi * tau)
    kap_hat = np.abs(-np.pi * A * tau * C_ * Wq)      # 174e yasası
    kap_hat0 = np.abs(-np.pi * A * tau * C_)          # RAKİP: W_pos yok

    # --- κ⁽²⁾: çarpımsal (S) ve fark (D) kanalları ---------------------
    qall = M5["q"].astype(np.int64)
    Aall, tall, mall = A_tam, M5["tau"], M5["mult"].astype(int)
    pos = {int(x): i for i, x in enumerate(qall)}
    CA = Aall * np.cos(np.pi * tall)          # A_q cos(πτ_q)
    S = np.zeros(len(q))
    D = np.zeros(len(q))
    for j, Q in enumerate(q):
        k, P = int(mall[idx[j]]), None
        # taban asal p:  Q = p^k
        P = int(round(float(Q) ** (1.0 / k))) if k > 1 else int(Q)
        while P ** k != Q:                      # kök yuvarlama düzeltmesi
            P += 1 if P ** k < Q else -1
        # ÇARPIMSAL: (p^i, p^{k−i}), i = 1..k−1  (sıralı çiftler)
        s = 0.0
        for i in range(1, k):
            a1, a2 = pos.get(P ** i), pos.get(P ** (k - i))
            if a1 is not None and a2 is not None:
                s += CA[a1] * CA[a2]
        # FARK: (p^{b+k}, p^b) ve tersi, b ≥ 1
        t = 0.0
        b = 1
        while True:
            hi = P ** (b + k)
            a1, a2 = pos.get(hi), pos.get(P ** b)
            if a1 is None or a2 is None:
                break
            t += 2 * CA[a1] * CA[a2]
            b += 1
        S[j], D[j] = s, t
    kap2 = (np.pi ** 2 * tau ** 2 / 2.0) * (S - D)
    kap1 = -np.pi * A * tau * C_
    # her iki mertebeye AYNI sönüm (W_pos) uygulanır — ek varsayım yok
    kap_hat_M = np.abs((kap1 + kap2) * Wq)
    oran2 = np.abs(kap2) / np.maximum(np.abs(kap1), 1e-300)
    kule = mult >= 2
    print("\n  %s:  σ_Ĉ = %.5f   ort(E) = %+.6f   λ=%s τ_ust=%s"
          % (g, sC, d["ortE"], lam, tau_ust))
    print("     [Ö-e0] maks |A(165,pencereli) − A(163,174e'nin)| = %.2e  %s"
          % (e0, "(erfc gazı: 174e'nin A'sı YANLIŞ olurdu)" if e0 > 1e-12
             else "(özdeş)"))
    print("     çizgi tipi: asal %d, kule %d ; S≠0 olan %d, D≠0 olan %d"
          % (int((~kule).sum()), int(kule.sum()), int((S != 0).sum()),
             int((D != 0).sum())))
    print("     τ bandı      n   n_kule  med|κ|/|κ̂|  med(ASAL)  med(KULE) "
          " med|κ|/|κ̂+κ⁽²⁾|(KULE)  RAKİP med|κ|/|κ̂₀|  |κ⁽²⁾/κ⁽¹⁾|(KULE)")
    sat = []
    for lo, hi in zip(INCE[:-1], INCE[1:]):
        m = (tau > lo) & (tau <= hi)
        if m.sum() < 3:
            continue
        r = kap / np.maximum(kap_hat, 1e-300)
        r0 = kap / np.maximum(kap_hat0, 1e-300)
        rM = kap / np.maximum(kap_hat_M, 1e-300)
        mk, ma = m & kule, m & ~kule
        rec = dict(lo=float(lo), hi=float(hi), n=int(m.sum()),
                   n_kule=int(mk.sum()), med=float(np.median(r[m])),
                   med_asal=float(np.median(r[ma])) if ma.sum() else np.nan,
                   med_kule=float(np.median(r[mk])) if mk.sum() else np.nan,
                   med_kule_M=float(np.median(rM[mk])) if mk.sum() else np.nan,
                   med0=float(np.median(r0[m])),
                   med0_asal=float(np.median(r0[ma])) if ma.sum() else np.nan,
                   o2_kule=float(np.median(oran2[mk])) if mk.sum() else np.nan,
                   o2_asal=float(np.median(oran2[ma])) if ma.sum() else np.nan,
                   o2_max=float(oran2[m].max()),
                   W=float(np.interp(0.5 * (lo + hi), tg, np.abs(W))))
        sat.append(rec)
        print("     %.2f–%.2f %6d %6d  %9.4f  %9.4f  %9s  %14s  %14.4f  %12s"
              % (lo, hi, rec["n"], rec["n_kule"], rec["med"], rec["med_asal"],
                 "%.4f" % rec["med_kule"] if rec["med_kule"] == rec["med_kule"]
                 else "—",
                 "%.4f" % rec["med_kule_M"]
                 if rec["med_kule_M"] == rec["med_kule_M"] else "—",
                 rec["med0"],
                 "%.4f" % rec["o2_kule"] if rec["o2_kule"] == rec["o2_kule"]
                 else "—"))
    OUT[g] = dict(sigma_C=sC, bant=sat, e0=e0,
                  n_asal=int((~kule).sum()), n_kule=int(kule.sum()),
                  nS=int((S != 0).sum()), nD=int((D != 0).sum()))

# ─────────────────────────────────────────────────────────────────────
print("\n" + "=" * 74)
print("ÖN-MÜHÜRLERİN HÜKMÜ (175a-P7)")
print("=" * 74)
for g, d in OUT.items():
    ust = [b for b in d["bant"] if b["lo"] >= 0.60]
    us5 = [b for b in d["bant"] if b["lo"] >= 0.55]
    ma = float(np.median([b["med_asal"] for b in ust]))
    mk = [b["med_kule"] for b in ust if b["med_kule"] == b["med_kule"]]
    mkM = [b["med_kule_M"] for b in ust if b["med_kule_M"] == b["med_kule_M"]]
    m0 = float(np.median([b["med0"] for b in us5]))
    a = "K3-a  medyan(ASAL, τ>0.60) = %.3f  (eşik 1.5) → **%s**" % (
        ma, "H-K3 ÖLDÜ" if ma >= 1.5 else "ayakta")
    iy = (np.median([abs(1 - x) for x in mkM]) <
          np.median([abs(1 - x) for x in mk])) if mk and mkM else False
    b = ("K3-b  |1−oran| KULE: yasa %.3f → +κ⁽²⁾ %.3f  → **%s**"
         % (np.median([abs(1 - x) for x in mk]) if mk else float("nan"),
            np.median([abs(1 - x) for x in mkM]) if mkM else float("nan"),
            "iyileşti" if iy else "İYİLEŞMEDİ ⇒ H-K3 ÖLDÜ"))
    c = ("K3-c  RAKİP medyan |κ|/|κ̂₀| (τ>0.55) = %.3f  (bant [0.8,1.25]) "
         "→ **%s**" % (m0, "RAKİP KAZANDI" if 0.8 <= m0 <= 1.25 else "rakip de ıska"))
    print("  %s:\n    %s\n    %s\n    %s" % (g, a, b, c))
    d["hukum"] = dict(K3a_med_asal=ma, K3a_oldu=bool(ma >= 1.5),
                      K3b_iyilesti=bool(iy), K3c_med0=m0,
                      K3c_kazandi=bool(0.8 <= m0 <= 1.25))
json.dump(OUT, open(S175 / "K3M.json", "w"), indent=1, ensure_ascii=False)
print("\n-> %s" % (S175 / "K3M.json"))
