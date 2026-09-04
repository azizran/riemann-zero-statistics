# -*- coding: utf-8 -*-
"""
174b — K1: GİRİŞİM ORANLARI, ÜÇ (VE DAHA ÇOK) GAZDA, 162'NİN MAKİNESİYLE
=========================================================================
KALEM/174 K1: "162'nin girişim-oranı ölçümü AYNI pencere/tabanla ÜÇ gazda
koşulur: gerçek (son), sadakatli ikiz (Hkeskin), kontrol L085. η-kanalı +
Ĉ-kanalı + bant-bant defteri. 162'nin 1.288/1.026'sı YENİDEN ölçülür."

ÖLÇÜM MAKİNESİ KOPYALANMADI: `162_configs/162_cekirdek.Taban162` AYNEN
kurulur (taban 0.40, cap 4000, τ_çizgi 0.86, NITER 0, BLOK 2000). Yani
çizgi evreni `pk_m(e^{0.86L})`, çıkarım `c_q(x) = 2⟨x_n e^{−i w_q m_n}⟩`
ve çizgi/süreklilik ayrımı 162'nin ta kendisidir.

EKLENEN ÜÇ ŞEY (162'de YOK, ön-kayıtta bilinmeyen olarak listelendi):
  1. X̃ (= ΔĈ) KANALI: c_q(X̃0) 165'in sitesinde (s_n = mid_{n+1}) —
     165'in `y_q = 2⟨X̃0_n e^{−iω_q s_n}⟩` tanımıyla aynı hizada.
  2. BANT-BANT DEFTERİ: çizgiler τ bantlarına bölünüp her bandın alanı
     ayrı sentezlenir; ÖZDEŞ defter
         Var(x) = Σ_b P_b + Cov(x_artık, x)   ,  P_b = Σ_{q∈b}|c_q|²/2
     (çünkü Cov(x_çizgi,b , x) = P_b ÖZDEŞ), ve bantlar arası kovaryans
     matrisi Var(x_çizgi) = Σ_b V_b + Σ_{b≠b'} Cov_{bb'}.
  3. ÜÇÜNCÜ MOMENTLER: rastgele fazlı çizgi alanında ÖZDEŞ SIFIR olan
     nesneler — skew ve m3 = ⟨e1·x1²⟩/(σ_e1σ_x1²).

ÖN-MÜHÜR (174b, koşudan önce — 174a'nın ön-kaydıyla tutarlı):
  Ö1  R_η(son) 162'nin 1.288'ini ±%1 içinde yeniden vermeli (aynı makine,
      aynı pencere). Vermezse boru hattı kimliği kırılmıştır ve yazılır.
  Ö2  R_η(Hkeskin) < R_η(son)  [H-G1'in yönü]. Ters çıkarsa H-G1 ÖLÜR.
  Ö3  R_Ĉ üç gazda 1.00-1.05 aralığında ve birbirinden ±%2 içinde.
  Ö4  Bant defteri ÖZDEŞ kapanmalı: |Σ_b P_b + Cov(artık,x) − Var(x)|
      / Var(x) ≤ 1e−12.
  Ö5  Çizgi payları 162'nin gerçek gaz için kaydettiği değerleri
      (η 0.487, Ĉ 0.918) ±0.005 içinde vermeli.

SONUÇ: yalnız gerçek koşudan (scratchpad/174/K1_<gaz>.json + log).
Kullanım: 174b_K1_girisim.py <gaz> [<gaz> ...]
"""
import importlib
import json
import sys
import time
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
for _p in ("162_configs", "160_configs", "159_configs", "156_configs",
           "155_configs", "154_configs"):
    sys.path.insert(0, str(QM / _p))
C162 = importlib.import_module("162_cekirdek")

SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/174")
TABAN = 0.40
CAP = 4000
# bant kenarları: regresyonun sildiği bölge (τ ≤ 0.40) + 0.40..0.86
KENAR = [0.0, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.86001]


def _skew(x):
    x = np.asarray(x, float)
    d = x - x.mean()
    s = d.std()
    return float((d ** 3).mean() / s ** 3) if s > 0 else float("nan")


def _m3(a, b):
    """⟨a·b²⟩/(σ_a σ_b²) — θ'nın ⟨e1 x1² e^{−iWs}⟩ çekirdeğinin W=0 vekili."""
    a = a - a.mean()
    b = b - b.mean()
    sa, sb = a.std(), b.std()
    return float((a * b * b).mean() / (sa * sb * sb)) if sa * sb > 0 else np.nan


def kanal_defteri(ad, x, c, w, m, tau, kenar, blok=C162.BLOK):
    """Bir kanalın girişim oranı + bant-bant defteri (fit yok)."""
    x0 = x - x.mean()
    V = float(np.var(x))
    P = float(np.sum(np.abs(c) ** 2) / 2.0)
    # bant alanları
    bnt = []
    alan = []
    for i in range(len(kenar) - 1):
        msk = (tau > kenar[i] - (1e-12 if i == 0 else 0.0)) & (tau <= kenar[i + 1])
        if i == 0:
            msk = tau <= kenar[1]
        nb = int(msk.sum())
        if nb == 0:
            bnt.append(dict(lo=kenar[i], hi=kenar[i + 1], n=0))
            continue
        f = C162.cizgi_kur(m, [c[msk]], w[msk], blok=blok)[0]
        alan.append(f)
        Pb = float(np.sum(np.abs(c[msk]) ** 2) / 2.0)
        bnt.append(dict(lo=kenar[i], hi=kenar[i + 1], n=nb, P=Pb,
                        V=float(np.var(f)),
                        ham_x=float(np.mean(f * x)),      # ≡ P_b (özdeşlik)
                        kov_x=float(np.mean((f - f.mean()) * x0)),
                        oran=Pb / float(np.var(f)) if np.var(f) > 0 else np.nan))
    A = np.array(alan)                       # (nb_dolu, N)
    A0 = A - A.mean(axis=1, keepdims=True)
    G = (A0 @ A0.T) / A.shape[1]             # bantlar arası kovaryans
    xciz = A.sum(axis=0)
    art = x - xciz
    # ÖZDEŞ defter (ham momentlerle):  ⟨x²⟩ = Σ_b P_b + ⟨x·x_artık⟩
    defter = float(np.sum([b["P"] for b in bnt if b.get("n")]) +
                   np.mean(art * x) - np.mean(x * x))
    return dict(
        ad=ad, Var=V, P=P, R=P / V, Var_cizgi=float(np.var(xciz)),
        Var_artik=float(np.var(art)),
        pay=float(1.0 - np.var(art) / V),
        kapanis=float((P + np.var(art)) / V),
        kov_ciz_x=float(np.mean((xciz - xciz.mean()) * x0)),
        kov_art_x=float(np.mean((art - art.mean()) * x0)),
        ham_ciz_x=float(np.mean(xciz * x)), ham_art_x=float(np.mean(art * x)),
        ort=float(np.mean(x)), ort_ciz=float(np.mean(xciz)),
        defter_kalinti=defter, defter_bagil=abs(defter) / float(np.mean(x * x)),
        ic_bant=float(np.trace(G)), capraz=float(G.sum() - np.trace(G)),
        G=[[float(v) for v in r] for r in G], bant=bnt,
        skew=_skew(x), skew_cizgi=_skew(xciz),
    ), xciz, art


def kos(veri):
    t0 = time.time()
    print("=" * 72, flush=True)
    print(f"174b / K1 — {veri}  (162 makinesi: taban={TABAN} cap={CAP} "
          f"τ_çizgi={C162.TAU_CIZGI} NITER={C162.NITER})", flush=True)
    print("=" * 72, flush=True)
    T = C162.Taban162(veri, taban=TABAN, cap=CAP)
    Y = T.Y
    print(f"  [Taban162 kuruldu  {time.time()-t0:.0f}s]", flush=True)
    tau = T.w / Y.L
    m = Y.mid

    out = dict(veri=veri, N=int(len(Y.m0)), Nn=int(T.Nn), L=float(Y.L),
               nline=int(len(T.q)), taban=TABAN, cap=CAP,
               tau_cizgi=C162.TAU_CIZGI, m_kimlik=T.m_kimlik,
               tau_min=float(tau.min()), tau_max=float(tau.max()),
               kenar=KENAR)

    # --- η kanalı (162'nin ta kendisi) --------------------------------
    tm = time.time()
    d_eta, eta_ciz, eta_art = kanal_defteri("eta", Y.eta, T.c_eta, T.w, m,
                                            tau, KENAR)
    d_eta["Var_e1"] = float(np.var(Y.e1))
    d_eta["R_e1"] = d_eta["P"] / d_eta["Var_e1"]
    print(f"  η   : Σ|c|²/2 = {d_eta['P']:.6f}  Var = {d_eta['Var']:.6f}  "
          f"**R_η = {d_eta['R']:.4f}**  (Var(e1) ile {d_eta['R_e1']:.4f})  "
          f"[{time.time()-tm:.0f}s]", flush=True)
    print(f"        çizgi payı {d_eta['pay']:.4f}  kapanış "
          f"{d_eta['kapanis']:.4f}  defter kalıntısı "
          f"{d_eta['defter_bagil']:.2e}", flush=True)

    # --- Ĉ kanalı ------------------------------------------------------
    tm = time.time()
    d_C, C_ciz, C_art = kanal_defteri("Chat", T.Chat, T.c_C, T.w, m, tau,
                                      KENAR)
    print(f"  Ĉ   : Σ|c|²/2 = {d_C['P']:.6f}  Var = {d_C['Var']:.6f}  "
          f"**R_Ĉ = {d_C['R']:.4f}**  çizgi payı {d_C['pay']:.4f}  "
          f"[{time.time()-tm:.0f}s]", flush=True)

    # --- X̃ (=ΔĈ) kanalı: 165'in sitesi s_n = mid[1:] -------------------
    tm = time.time()
    s = np.ascontiguousarray(m[1:])
    x1 = Y.Xtil0
    c_X = C162.cizgi_cikar(s, [x1], T.w)[0]
    print(f"  [X̃ çıkarımı {time.time()-tm:.0f}s]", flush=True)
    tm = time.time()
    d_X, X_ciz, X_art = kanal_defteri("Xtil", x1, c_X, T.w, s, tau, KENAR)
    print(f"  X̃   : Σ|c|²/2 = {d_X['P']:.6f}  Var = {d_X['Var']:.6f}  "
          f"**R_X = {d_X['R']:.4f}**  çizgi payı {d_X['pay']:.4f}  "
          f"[{time.time()-tm:.0f}s]", flush=True)

    # --- üçüncü momentler (rastgele fazda ÖZDEŞ SIFIR) ------------------
    e1 = Y.e1 - Y.e1.mean()
    out["m3"] = dict(
        olc=_m3(e1, x1),                       # ⟨e1 x1²⟩ / (σ σ²)
        cizgi=_m3(eta_ciz[1:], X_ciz),         # çizgi alanlarının eşleniği
        capraz=_m3(e1, X_ciz), capraz2=_m3(eta_ciz[1:], x1),
        skew_e1=_skew(e1), skew_x1=_skew(x1),
        skew_eta_ciz=d_eta["skew_cizgi"], skew_X_ciz=d_X["skew_cizgi"],
        skew_Chat=d_C["skew"], skew_C_ciz=d_C["skew_cizgi"],
        skew_ds=_skew(Y.ds),
    )
    print("  m3  : ⟨e1x1²⟩/(σσ²) ölç = %+.5f   çizgi = %+.5f   "
          "çapraz = %+.5f / %+.5f" % (out["m3"]["olc"], out["m3"]["cizgi"],
                                      out["m3"]["capraz"],
                                      out["m3"]["capraz2"]), flush=True)
    print("        skew: e1 %+.4f  x1 %+.4f  η_çiz %+.4f  X_çiz %+.4f  "
          "Ĉ %+.4f  ds %+.4f"
          % (out["m3"]["skew_e1"], out["m3"]["skew_x1"],
             out["m3"]["skew_eta_ciz"], out["m3"]["skew_X_ciz"],
             out["m3"]["skew_Chat"], out["m3"]["skew_ds"]), flush=True)

    out["eta"], out["Chat"], out["Xtil"] = d_eta, d_C, d_X
    out["sigma"] = dict(ds=float(np.std(Y.ds)), X=float(np.std(x1)),
                        Chat=float(np.std(T.Chat)), eta=float(np.std(Y.eta)))
    out["sure_s"] = time.time() - t0
    SCR.mkdir(parents=True, exist_ok=True)
    p = SCR / f"K1_{veri}.json"
    p.write_text(json.dumps(out, indent=1))
    np.savez_compressed(SCR / f"K1c_{veri}.npz", q=T.q, w=T.w, tau=tau,
                        c_eta=T.c_eta, c_C=T.c_C, c_X=c_X)
    print(f"  -> {p}   ({out['sure_s']/60:.1f} dk)", flush=True)

    print("\n  BANT DEFTERİ (η):  lo–hi   n   P_b      V_b     P/V    "
          "Kov(f_b,η)", flush=True)
    for b in d_eta["bant"]:
        if not b.get("n"):
            continue
        print("    %.2f–%.2f %5d  %8.6f %8.6f %6.3f  %+9.6f"
              % (b["lo"], b["hi"], b["n"], b["P"], b["V"], b["oran"],
                 b["kov_x"]), flush=True)
    print("    iç-bant Σ V_b = %.6f   bantlar-arası Σ Kov = %+.6f   "
          "Var(η_çizgi) = %.6f"
          % (d_eta["ic_bant"], d_eta["capraz"], d_eta["Var_cizgi"]),
          flush=True)
    return out


if __name__ == "__main__":
    for g in (sys.argv[1:] or ["son"]):
        kos(g)
