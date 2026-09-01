"""
162 / T1 — KİMLİK SINAVI:  C_n = cumsum(ds) ↔ −S_merdiven(z_n)
==============================================================
KALEM: ds_k = S(z_k) − S(z_{k+1}) (teleskop) ⇒ C_n = Σ_{k<n} ds_k
= −S(z_n) + sabit. Bu modül o kimliği ölçer ve C'nin künyesini çıkarır:

  (a) korelasyon r, varyans oranı, regresyon eğimi          [C_ds ↔ −S(z)]
  (b) aynı sınav BOND değişkeniyle                          [Ĉ ↔ −S(mid)]
  (c) artığın TAYFI: her bantta |Ĉ(w)|², |−S(w)|², |artık(w)|²
      → hangi τ'larda eksik (kesme τ>0.9 mi, yoksa taban-içi mi?)
  (d) C'nin künyesi: Var(C), Var(ΣX̃) (sürüklenme), C↔X̃ korelasyonu,
      C'nin bant güçleri ve beklenen merdiven gücü Σ(a_q)²/2 ile oranı.

Kullanım:  162_T1.py [veri ...]      (varsayılan: son keskin A4)
"""
import importlib
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
C = importlib.import_module("162_cekirdek")
TWO_PI = 2 * np.pi
TAU_MAX = 0.9                      # merdiven kesmesi (görev)


def _kor(a, b):
    a = a - a.mean(); b = b - b.mean()
    return float(np.dot(a, b) / np.sqrt(np.dot(a, a) * np.dot(b, b)))


def t1(veri, taban=0.40, cap=4000, tau_analiz=1.15, bant_cizgi=140,
       tohum=21):
    out = {"veri": veri, "taban": taban}
    z = C.KOS155.veri_yukle(veri)
    Ç = C.K155.eta_onbellek(z, veri, taban, cap)
    mid, L, ds, eta = Ç["mid"], Ç["L"], Ç["ds"], Ç["eta"]
    Nn = len(ds)
    Xtil = (mid[1:] - mid[:-1]) * L / TWO_PI - 1.0
    dsB = C.K156._bond(ds)
    r = Xtil - dsB
    Cds = np.cumsum(ds)                       # ↔ −S(z_{n+1})
    Chat = np.concatenate(([0.0], np.cumsum(dsB)))     # ↔ −S(mid_n)
    CX = np.cumsum(Xtil)                      # D1'in ham C'si (sürüklenmeli)

    print(f"\n=== T1 / {veri} (taban {taban}) ===")
    print(f"  N={len(z)}  Nn={Nn}  L={L:.5f}  t∈[{z[0]:.0f},{z[-1]:.0f}]")
    print(f"  Var(ds)={np.var(ds):.5f}  Var(η)={np.var(eta):.5f}  "
          f"Var(X̃)={np.var(Xtil):.5f}  Var(dsΔ)={np.var(dsB):.5f}")
    print(f"  Var(C_ds)={np.var(Cds):.5f}  Var(Ĉ)={np.var(Chat):.5f}  "
          f"Var(ΣX̃)={np.var(CX):.4g}  ← ΣX̃'nin devasa varyansı SÜRÜKLENME")
    out.update(N=len(z), Nn=Nn, L=L, var_ds=float(np.var(ds)),
               var_eta=float(np.var(eta)), var_Xtil=float(np.var(Xtil)),
               var_dsB=float(np.var(dsB)), var_Cds=float(np.var(Cds)),
               var_Chat=float(np.var(Chat)), var_CX=float(np.var(CX)),
               r_std=float(np.std(r)), r_min=float(r.min()),
               r_max=float(r.max()))

    # --- merdiven ----------------------------------------------------
    t0 = time.time()
    Sz, mk = C.merdiven(z[1:], L, TAU_MAX)          # S(z_{n+1})
    Sm, _ = C.merdiven(mid, L, TAU_MAX)             # S(mid_n)
    print(f"  merdiven: {mk['nq']} çizgi (τ≤{TAU_MAX}), Σa²/2="
          f"{0.5*np.sum(mk['aq']**2):.5f}  ({time.time()-t0:.0f} s)")
    out["merdiven_nq"] = mk["nq"]
    out["merdiven_var_kuram"] = float(0.5 * np.sum(mk["aq"] ** 2))
    out["merdiven_var_olculen_z"] = float(np.var(Sz))
    out["merdiven_var_olculen_mid"] = float(np.var(Sm))

    for ad, Cv, Sv in (("C_ds ↔ −S(z_{n+1})", Cds, -Sz),
                       ("Ĉ    ↔ −S(mid_n)", Chat, -Sm)):
        rr = _kor(Cv, Sv)
        vr = float(np.var(Cv) / np.var(Sv))
        a = Cv - Cv.mean(); b = Sv - Sv.mean()
        egim = float(np.dot(a, b) / np.dot(b, b))
        art = a - egim * b
        print(f"  {ad}:  r={rr:+.5f}  Var oranı={vr:.4f}  eğim={egim:.4f}  "
              f"artık payı={np.var(art)/np.var(a):.4f}")
        out[ad.split()[0]] = dict(r=rr, var_orani=vr, egim=egim,
                                  artik_payi=float(np.var(art) / np.var(a)))

    # --- artığın tayfı (bant bant, konum uzayında) --------------------
    # Çizgi sayısı τ ile üstel büyür (τ=1.15'te ~4e4). Her 0.05'lik bantta
    # EN ÇOK `bant_cizgi` çizgi rastgele (sabit tohum) örneklenir; bandın
    # toplam gücü ⟨|c|²/2⟩·n ile kestirilir. Bu, tayfın ŞEKLİNİ ölçmek
    # için yeterlidir ve maliyeti sabit tutar.
    qm = C.C154.pk_m(int(np.exp(tau_analiz * L)))
    qs_all = np.array(sorted(qm))
    w_all = np.log(qs_all.astype(float))
    tau_all = w_all / L
    rng = np.random.default_rng(tohum)
    ken = np.arange(0.0, tau_analiz + 1e-9, 0.05)
    sec, agir, bkn, nn = [], [], [], []
    for i in range(len(ken) - 1):
        idx = np.where((tau_all > ken[i]) & (tau_all <= ken[i + 1]))[0]
        if len(idx) == 0:
            continue
        j = idx if len(idx) <= bant_cizgi else rng.choice(idx, bant_cizgi,
                                                          replace=False)
        sec.append(np.sort(j))
        agir.append(len(idx) / len(j))
        ms_ = np.array([qm[int(x)] for x in qs_all[j]], float)
        aq_ = 1.0 / (np.pi * ms_ * np.sqrt(qs_all[j].astype(float)))
        bkn.append(float(np.sum(0.5 * aq_ ** 2) * len(idx) / len(j)))
        nn.append((ken[i], ken[i + 1], len(idx), len(j)))
    w = w_all[np.concatenate(sec)]
    a = Chat - Chat.mean()
    b = (-Sm) - (-Sm).mean()
    art = a - b * (np.dot(a, b) / np.dot(b, b))
    cA, cB, cR = C.cizgi_cikar(mid, [a, b, art], w)
    print(f"\n  ARTIĞIN TAYFI (Ĉ vs −S_merdiven, konum uzayı; {len(w)} "
          f"örneklenmiş çizgi, τ≤{tau_analiz})")
    print("   τ bandı     çizgi(örnek)  |Ĉ|²/2      |S|²/2      |artık|²/2  "
          "artık/Ĉ   Ĉ/beklenen  cos²(πτ)")
    print("   (Ĉ, dsΔ'nın kümülatifidir: bond ortalaması bir tona cos(πτ) "
          "kazancı verir ⇒ beklenen oran cos²(πτ). τ>0.9'da merdiven KESİK "
          "ve τ≈1'de kafes ALIAS'ı devrededir — o satırlar ölçüm değil.)")
    tayf = []
    p = 0
    for k, (lo_, hi_, ntot, nsec) in enumerate(nn):
        sl = slice(p, p + nsec); p += nsec
        g_ = agir[k]
        PA = float(np.sum(np.abs(cA[sl]) ** 2) / 2 * g_)
        PB = float(np.sum(np.abs(cB[sl]) ** 2) / 2 * g_)
        PR = float(np.sum(np.abs(cR[sl]) ** 2) / 2 * g_)
        PE = bkn[k]
        cs2 = float(np.cos(np.pi * 0.5 * (lo_ + hi_)) ** 2)
        print(f"   {lo_:.2f}–{hi_:.2f}  {ntot:6d}({nsec:4d})  "
              f"{PA:.4e}  {PB:.4e}  {PR:.4e}  {PR/PA if PA else np.nan:7.4f}  "
              f"{PA/PE if PE else np.nan:8.4f}   {cs2:.4f}")
        tayf.append(dict(lo=float(lo_), hi=float(hi_), n=int(ntot),
                         nsec=int(nsec), PC=PA, PS=PB, PR=PR, Pbek=PE))
    out["tayf"] = tayf

    # --- C'nin künyesi -----------------------------------------------
    kk = dict(kor_C_X=_kor(Chat[:-1], Xtil), kor_C_dsB=_kor(Chat[:-1], dsB),
              kor_C_eta=_kor(Chat, eta), kor_Cds_X=_kor(Cds[:-1], Xtil),
              kor_r_C=_kor(np.cumsum(r), Chat[1:]))
    print(f"\n  C'nin künyesi: corr(Ĉ,X̃)={kk['kor_C_X']:+.5f}  "
          f"corr(Ĉ,dsΔ)={kk['kor_C_dsB']:+.5f}  corr(Ĉ,η)={kk['kor_C_eta']:+.5f}"
          f"  corr(Σr,Ĉ)={kk['kor_r_C']:+.5f}")
    out["kunye"] = kk
    return out


if __name__ == "__main__":
    veriler = sys.argv[1:] or ["son", "keskin", "A4"]
    C.SCR.mkdir(parents=True, exist_ok=True)
    hepsi = [t1(v) for v in veriler]
    p = C.SCR / "T1.json"
    p.write_text(json.dumps(hepsi, indent=1))
    print(f"\n-> {p}")
