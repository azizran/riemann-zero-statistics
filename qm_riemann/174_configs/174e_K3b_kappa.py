# -*- coding: utf-8 -*-
"""
174e — K3b: κ'NIN KAPALI YASASI ve GERÇEK−İKİZ FARKININ ANATOMİSİ
==================================================================
**ÖN-MÜHÜRSÜZ İKİNCİ TUR (dürüstlük).** 174d'nin ön-mühürlü çıplak
limiti `κ ≈ −πA_qτ_q` İKİ maddede öldü: (Ö1) ort(E) POZİTİF çıktı,
(Ö3) kaçak birkaç büyük çizgide değil binlerce yüksek-τ çizgisinde.
Ölüm kurtarılmadı; bu betik ÖLÇÜMDEN SONRA yazıldı ve öyle etiketlenir.

174d'nin bant defteri iki şeyi gösterdi:
  (i) ⟨cos Δφ⟩ τ < 0.5'te ≈ −1, τ > 0.5'te ≈ +1 — yani κ'nın işareti
      τ = 0.5'te DÖNÜYOR;
  (ii) |κ| çıplak −πA_qτ_q'nun ≈ 0.2 katı ve τ ile hızla düşüyor.
İkisinin de kapalı karşılığı var ve **yeni parametre gerektirmiyor**:

  * İşaret: ölçüm sitesi `s_n = mid_{n+1}` — SIFIR değil ORTA NOKTA
    tarağıdır. Yarım-gap kayması κ'ya `cos(πτ)` çarpanı takar
    (`B_q = b_q cos(πτ_q)` ile aynı çarpan), ve `cos(πτ)` τ = 0.5'te
    işaret değiştirir.
  * Sönüm: `s_n = s_0 + ḡ(n + Ĉ_n)` ⇒ rezonans genliği tarağın kendi
    faz seğirmesiyle çarpılır; bu, 166/167'nin ZATEN ÖLÇTÜĞÜ
    karakteristik fonksiyondur: `W_pos(τ) = ⟨e^{2πiτ Ĉ_n}⟩`.

  ***  κ(ω_q) ≈ −π A_q τ_q cos(πτ_q) · W_pos(τ_q)  ***      (K3b yasası)

`A_q = λ a_q w_q` gazın İNŞA genliği, `W_pos` ölçülür — **sıfır serbest
parametre**. Sınav: ölçülen κ ile bu yasanın oranı.

Ayrıca gerçek−ikiz farkı üç çarpana ayrıştırılır (bant bant):
    ort(E) = Σ_q |hp_q|·|κ_q|·cos Δφ_q
    genlik payı (|hp|), rezonans payı (|κ|), faz payı (cos Δφ).
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
ORT = importlib.import_module("167_ortak")
K165.PENCERE.update(ORT.pencere_dict())

SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S174 = SCR / "174"
TWO_PI = 2 * np.pi
GAZ = sys.argv[1:] or ["son", "Hkeskin", "L085"]
BANT = [0.40, 0.50, 0.60, 0.70, 0.80, 0.95001]


def w_pos_egri(veri, tgrid):
    """W_pos(τ) = ⟨e^{2πiτ Ĉ_n}⟩ , Ĉ = cumsum(X̃0) kübik trendsiz
    (165.koherans / 167_olcum ile AYNI tanım)."""
    Y = K165.gaz(veri, 0.40, 4000)
    C = np.cumsum(Y.Xtil0)
    n = np.arange(len(C), dtype=float)
    n = (n - n.mean()) / n[-1]
    V = np.vstack([np.ones_like(n), n, n * n, n ** 3]).T
    c, *_ = np.linalg.lstsq(V, C, rcond=None)
    Cd = C - V @ c
    out = np.array([np.mean(np.exp(1j * TWO_PI * t * Cd)) for t in tgrid])
    return out, float(np.std(Cd)), Y


print("=" * 74)
print("174e — K3b: κ'nın kapalı yasası (ÖN-MÜHÜRSÜZ, ikinci tur)")
print("=" * 74)
tg = np.linspace(0.0, 0.96, 241)
D = {}
for g in GAZ:
    p = S174 / f"K3_{g}.json"
    if not p.exists():
        print(f"  ({g} K3'ü yok, atlanıyor)")
        continue
    d = json.load(open(p))
    W, sC, Y = w_pos_egri(g, tg)
    tau = np.array(d["tau"])
    lam = d["lam"]
    q = np.array(d["q"])
    # a_q'yu merdivenden yeniden üret (mult gerekli), A_q = λ a_q w_q
    M = importlib.import_module("163_cekirdek").merdiven(float(Y.L), 0.95)
    qm = {int(x): i for i, x in enumerate(M["q"])}
    idx = np.array([qm[int(x)] for x in q])
    a = M["a"][idx]
    A = lam * a * (tau <= d["tau_ust"] + 1e-12)
    kap_abs = np.array(d["kap_abs"])
    dphi = np.array(d["dphi"])
    hp = np.array(d["hp_abs"])
    xi = np.array(d["xi"])
    Wq = np.interp(tau, tg, W.real) + 1j * np.interp(tau, tg, W.imag)
    kap_hat = -np.pi * A * tau * np.cos(np.pi * tau) * Wq
    D[g] = dict(tau=tau, q=q, A=A, kap_abs=kap_abs, dphi=dphi, hp=hp, xi=xi,
                kap_hat=kap_hat, sC=sC, ortE=d["ortE"], V_O=d["V_O"],
                mu2=d["mu2"], W=W)
    print(f"\n  {g}:  σ_Ĉ(trendsiz) = {sC:.5f}   ort(E) = {d['ortE']:+.6f}"
          f"   μ̂² = {d['mu2']:.5f}")
    print("     τ bandı      n   medyan|κ|  medyan|κ̂|  ORAN(med)  "
          "işaret uyumu  ⟨cosΔφ⟩_ağ   Σξ")
    ince = np.arange(0.40, 0.96, 0.05)
    for lo, hi in zip(ince[:-1], ince[1:]):
        m = (tau > lo) & (tau <= hi)
        if m.sum() < 3:
            continue
        r = float(np.median(kap_abs[m] / np.abs(kap_hat[m])))
        # işaret yasası: sign(κ) = −sign(cos πτ) ⇔ cos Δφ'nin işareti
        bek = -np.sign(np.cos(np.pi * tau[m]))
        uy = float(np.mean(np.sign(np.cos(dphi[m])) == bek))
        cd = float(np.sum(xi[m]) / np.sum(hp[m] * kap_abs[m]))
        print("     %.2f–%.2f %6d  %9.3e %9.3e %9.4f   %8.3f     %+8.4f  "
              "%+.3e" % (lo, hi, int(m.sum()), float(np.median(kap_abs[m])),
                         float(np.median(np.abs(kap_hat[m]))), r, uy, cd,
                         float(np.sum(xi[m]))))

# ---- gerçek − ikiz farkının ayrıştırması ----------------------------
if "son" in D and "Hkeskin" in D:
    print("\n" + "=" * 74)
    print("GERÇEK − İKİZ: DC kaçağının bant-bant farkı")
    print("=" * 74)
    a_, b_ = D["son"], D["Hkeskin"]
    print("  ort(E): son %+.6f   ikiz %+.6f   oran %.4f   (μ̂² oranı %.4f)"
          % (a_["ortE"], b_["ortE"], a_["ortE"] / b_["ortE"],
             a_["mu2"] / b_["mu2"]))
    print("\n   τ bandı    Σξ(son)     Σξ(ikiz)    Δ           Δ payı%   "
          "|hp| oranı  |κ| oranı  cosΔφ oranı")
    tot = a_["ortE"] - b_["ortE"]
    sat = []
    for i in range(len(BANT) - 1):
        m = (a_["tau"] > BANT[i]) & (a_["tau"] <= BANT[i + 1])
        if not m.any():
            continue
        xa, xb = float(a_["xi"][m].sum()), float(b_["xi"][m].sum())
        ha = float(np.sum(a_["hp"][m])) / float(np.sum(b_["hp"][m]))
        ka = float(np.sum(a_["kap_abs"][m])) / float(np.sum(b_["kap_abs"][m]))
        ca = (xa / float(np.sum(a_["hp"][m] * a_["kap_abs"][m]))) / \
             (xb / float(np.sum(b_["hp"][m] * b_["kap_abs"][m])))
        sat.append(dict(lo=BANT[i], hi=BANT[i + 1], n=int(m.sum()),
                        xi_son=xa, xi_ikiz=xb, d=xa - xb,
                        pay=100 * (xa - xb) / tot, hp=ha, kap=ka, cos=ca))
        print("   %.2f–%.2f  %+.4e %+.4e %+.4e  %+7.1f   %8.4f  %8.4f  "
              "%+9.4f" % (BANT[i], BANT[i + 1], xa, xb, xa - xb,
                          100 * (xa - xb) / tot, ha, ka, ca))
    print("   TOPLAM     %+.4e %+.4e %+.4e   100.0"
          % (a_["ortE"], b_["ortE"], tot))
    json.dump(dict(bant=sat, ortE_son=a_["ortE"], ortE_ikiz=b_["ortE"],
                   mu2_oran=a_["mu2"] / b_["mu2"]),
              open(S174 / "K3b.json", "w"), indent=1)
    print("\n-> %s" % (S174 / "K3b.json"))
