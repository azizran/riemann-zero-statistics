# -*- coding: utf-8 -*-
"""
174d — K3: μ̂²_E'NİN ANATOMİSİ (DC KAÇAĞININ ÇİZGİ-ÇİZGİ TAYFI)
================================================================
172-G3'ün en aykırı tek sayısı `μ̂²_E = −%23.8`'di ve şu adres verilmişti:
"model alanının DC kaçağı = sıfır tarağının asal rezonansları". 174-K3
o cümleyi ÖLÇÜYOR.

ÖZDEŞLİK (165_cekirdek.Model165.alanlar'dan, yaklaşım yok):
    E_n = Σ_q Re[hp_q e^{iω_q s_n}]   ⇒   ort(E) = Σ_q Re[hp_q κ(ω_q)] ,
    κ(ν) := ⟨e^{iν s_n}⟩ ,  μ̂²_E = ort(E)²/Var(e1)
Yani DC kaçağı ÇİZGİ-ÇİZGİ ayrışır:  ξ_q := Re[hp_q κ(ω_q)] ,
    ξ_q = |hp_q|·|κ(ω_q)|·cos Δφ_q ,  Δφ_q := arg hp_q + arg κ(ω_q)
üç çarpan: **çizgi genliği**, **tarak rezonansı**, **bağıl faz**.

ÇIPLAK LİMİT (ön-mühür: koşudan önce türetildi):
    dN = (N̄' + S')dt , S(t) = −Σ_q A_q sin(ω_q t) , A_q = λ·a_q·w_q
    ⇒ κ(ω_Q) ≈ −A_Q ω_Q/(2N̄') = −π A_Q τ_Q          (REEL, NEGATİF)
    ve hp_q ≈ b_q = 2A_q sin(πτ_q) (reel, τ<1'de pozitif) ⇒
    **ξ_q^çıplak = −2π A_q² τ_q sin(πτ_q)** ,  Δφ_q^çıplak = π
    ⇒ ort(E) < 0 ve |ort(E)| ∝ λ² (μ̂² ∝ λ², çünkü Var(η) ∝ λ²).
Bu, 172b'nin "μ̂²_E λ ile tekdüze büyür" gözleminin kapalı formudur ve
**(−) işaretin kökenidir**: kaçak, tarağın kendi asal rezonansının
çizgiyle TERS fazda olmasından doğar.

ÖN-MÜHÜR (174d, koşudan önce):
 Ö1 ort(E) < 0 (üç gazda da). Aksi çıkarsa çıplak türetim ölür.
 Ö2 Σ_q ξ_q'dan hesaplanan μ̂² = ort(E)²/Var(e1), 172/G1.json'un
    μ̂²_E'sini ‰5 içinde vermeli (bağımsız yol denetimi).
 Ö3 Kaçak KÜÇÜK q'da (büyük genlikli, düşük τ) yoğunlaşmalı: en büyük
    |ξ_q| taşıyan 10 çizginin payı ≥ %50.
 Ö4 |κ(ω_q)| çıplak −πA_qτ_q ile aynı mertebede (0.3 … 3 kat) olmalı.
 Ö5 GERÇEK−İKİZ FARKI: |ort(E)| gerçekte ikizden KÜÇÜK olmalı
    (μ̂² −%23.8'in kaynağı). Fark, "genlik payı" / "rezonans payı" /
    "faz payı" olarak ayrıştırılır (log-ayrışım, özdeş).
Kullanım: 174d_K3_dc.py <gaz> [<gaz> ...]
"""
import importlib
import json
import sys
import time
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
S174, S172 = SCR / "174", SCR / "172"
G1 = json.load(open(S172 / "G1.json"))
TAU_C = 0.95
BANT = [0.0, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.95001]


def kappa(s, nu, blok=20000, fblok=512):
    """κ(ν) = ⟨e^{iν s_n}⟩ (bellek: fblok×blok)."""
    N = len(s)
    acc = np.zeros(len(nu), dtype=complex)
    one = None
    for f0 in range(0, len(nu), fblok):
        fs = slice(f0, min(f0 + fblok, len(nu)))
        nc = nu[fs]
        for a in range(0, N, blok):
            sl = slice(a, min(a + blok, N))
            if one is None or len(one) != (sl.stop - sl.start):
                one = np.ones(sl.stop - sl.start)
            arg = np.outer(s[sl], nc)
            acc[fs] += one @ np.cos(arg) + 1j * (one @ np.sin(arg))
            del arg
    return acc / N


def kos(veri):
    t0 = time.time()
    print("=" * 72, flush=True)
    print(f"174d / K3 — {veri}", flush=True)
    print("=" * 72, flush=True)
    kun = ORT.KUNYE.get(veri, {})
    lam = kun.get("lam") or 1.0
    tau_ust = kun.get("tau_ust") or 1.00
    Mo = K165.Model165(veri, 0.40, 4000, 0.95, "olculen", tau_c=TAU_C)
    Mo.sec("olculen", TAU_C)
    M, msk = Mo.M, Mo.msk
    w = M["w"][msk]
    tau = M["tau"][msk]
    q = M["q"][msk]
    mult = M["mult"][msk]
    a = M["a"][msk]                        # çıplak a_q (λ ve kesim HARİÇ)
    hp = Mo.hp[msk]
    e1 = Mo.Y.e1 - Mo.Y.e1.mean()
    V_O = float(np.var(e1))
    # gazın GERÇEK merdiven genliği: A_q = λ·a_q·w_q  (keskin kesim)
    A = lam * a * (tau <= tau_ust + 1e-12)
    print(f"  {len(w)} çizgi (τ≤{TAU_C}), λ={lam} τ_ust={tau_ust}, "
          f"Var(e1)={V_O:.6f}", flush=True)

    tm = time.time()
    kap = kappa(Mo.s, w)
    print(f"  [κ(ω_q) {time.time()-tm:.0f}s]", flush=True)

    xi = (np.conj(hp) * 0 + hp * kap).real       # ξ_q = Re[hp_q κ_q]
    ortE = float(xi.sum())
    mu2 = ortE ** 2 / V_O
    mu2_ref = G1[veri]["E"]["mu2"]
    dphi = np.angle(hp) + np.angle(kap)
    dphi = (dphi + np.pi) % (2 * np.pi) - np.pi
    kap_cip = -np.pi * A * tau
    xi_cip = -2 * np.pi * A ** 2 * tau * np.sin(np.pi * tau)
    ok = tau > 0.40                              # η regresyonunun bıraktığı
    print("  ort(E) = %+.6f   μ̂²(ξ toplamı) = %.6f   μ̂²(172/G1) = %.6f"
          "   fark %.2e" % (ortE, mu2, mu2_ref, abs(mu2 / mu2_ref - 1)),
          flush=True)
    print("  ÇIPLAK: Σξ^çıplak(τ>0.40) = %+.6f  (ölçülen %+.6f, oran %.3f)"
          % (xi_cip[ok].sum(), xi[ok].sum(),
             xi[ok].sum() / xi_cip[ok].sum()), flush=True)
    print("  Ö1 ort(E)<0 : %s   Ö2 (‰5): %s"
          % ("✓" if ortE < 0 else "✗ ÖLDÜ",
             "✓" if abs(mu2 / mu2_ref - 1) < 5e-3 else "✗ ÖLDÜ"), flush=True)

    idx = np.argsort(-np.abs(xi))
    top = idx[:10]
    pay = float(np.abs(xi[top]).sum() / np.abs(xi).sum())
    print("  Ö3 en büyük 10 çizginin |ξ| payı = %.1f%%  (%s)"
          % (100 * pay, "✓" if pay >= 0.50 else "✗ ÖLDÜ"), flush=True)
    print("\n   q      τ     |hp|      |κ|      Δφ/π     ξ_q        "
          "ξ^çıplak   |κ|/|κ^çıplak|")
    for i in top:
        print("  %6d %.4f %.6f %.6f %+7.4f %+.3e %+.3e %8.3f"
              % (q[i], tau[i], abs(hp[i]), abs(kap[i]), dphi[i] / np.pi,
                 xi[i], xi_cip[i],
                 abs(kap[i]) / abs(kap_cip[i]) if kap_cip[i] else np.nan),
              flush=True)

    # --- bant defteri --------------------------------------------------
    bant = []
    print("\n  BANT DEFTERİ:  τ      n     Σξ_q      Σ|hp||κ|   ⟨cosΔφ⟩_ağ"
          "   Σξ^çıplak")
    for i in range(len(BANT) - 1):
        m = (tau > BANT[i]) & (tau <= BANT[i + 1]) if i else (tau <= BANT[1])
        if not m.any():
            continue
        mod = float(np.sum(np.abs(hp[m]) * np.abs(kap[m])))
        s_ = float(xi[m].sum())
        rec = dict(lo=BANT[i], hi=BANT[i + 1], n=int(m.sum()), xi=s_,
                   mod=mod, cos=s_ / mod if mod else np.nan,
                   xi_cip=float(xi_cip[m].sum()),
                   hp=float(np.sum(np.abs(hp[m]) ** 2) / 2),
                   kap=float(np.mean(np.abs(kap[m]))))
        bant.append(rec)
        print("   %.2f–%.2f %5d  %+.3e %.3e  %+7.4f   %+.3e"
              % (rec["lo"], rec["hi"], rec["n"], rec["xi"], rec["mod"],
                 rec["cos"], rec["xi_cip"]), flush=True)

    out = dict(veri=veri, lam=lam, tau_ust=tau_ust, nline=int(len(w)),
               V_O=V_O, ortE=ortE, mu2=mu2, mu2_G1=mu2_ref,
               xi_toplam_cip=float(xi_cip[ok].sum()),
               xi_toplam_olc=float(xi[ok].sum()),
               pay_top10=pay, bant=bant,
               q=[int(x) for x in q], tau=[float(x) for x in tau],
               hp_abs=[float(abs(x)) for x in hp],
               kap_abs=[float(abs(x)) for x in kap],
               dphi=[float(x) for x in dphi], xi=[float(x) for x in xi],
               kap_cip=[float(x) for x in kap_cip],
               xi_cip=[float(x) for x in xi_cip],
               sure_s=time.time() - t0)
    S174.mkdir(parents=True, exist_ok=True)
    p = S174 / f"K3_{veri}.json"
    p.write_text(json.dumps(out))
    print(f"\n  -> {p}  ({(time.time()-t0)/60:.1f} dk)", flush=True)
    return out


if __name__ == "__main__":
    for g in (sys.argv[1:] or ["son"]):
        kos(g)
