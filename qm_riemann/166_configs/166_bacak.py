"""
166 — BACAK YAPISI ve SÖNÜM ÇARPANLARI (H-K1'in TÜRETİMİ)
==========================================================
Görev: "165 çekirdeğindeki bacak yapısından sönüm çarpımını TÜRET —
hangi bacak hangi A ile e^{−A²σ²/2}-tipi faktör taşıyor".

──────────────────────── TÜRETİM (kalem) ────────────────────────
Merdiven:  ds_n = Σ_q 2a_q sin(ω_q g_n/2)·cos(ω_q m_n)
           m_n = (z_n+z_{n+1})/2 ,  g_n = z_{n+1}−z_n = ḡ(1+ds_n)
165'in sitesi  s_n ≡ m_{n+1}.  Üç alan, üç FARKLI faz referansı:

(E bacağı)  η_{n+1} ⊃ 2a_q sin(πτ_q(1+ds_{n+1}))·cos(ω_q s_n)
    ⇒ h'_q = 2⟨η_{n+1}e^{−iω_q s_n}⟩ = b_q·⟨cos(πτ_q ds)⟩ + O(⟨sin⟩)
    **A_E = πτ_q , dalgalanma ds  ⇒  W_amp(τ) = ⟨cos(πτ·ds)⟩
      ≈ exp(−½π²τ²σ_ds²)**  (YARIM-GAP genlik modülasyonu; 165 §6b'nin
      sin⁴(πτ) üyesinin kaynağı)

(X bacağı)  X̃_n = ½(ds_n+ds_{n+1}); q-bileşeni s-sitesinde
    ½b_q[e^{−iω_q(m_{n+1}−m_n)} + 1]·e^{iω_q s_n},
    ω_q(m_{n+1}−m_n) = ω_q ḡ(1+X̃_n) = 2πτ_q(1+X̃_n)
    ⇒ y_q = ½ b_q W_amp(τ_q)·[ W_X(τ_q)e^{−2πiτ_q} + 1 ]
    **A_X = 2πτ_q , dalgalanma X̃  ⇒  W_X(τ) = ⟨e^{−2πiτ·X̃}⟩
      ≈ exp(−2π²τ²σ_X̃²)**   (SİTE-KAYMASI sönümü)
    ÇIPLAK limit W_amp=W_X=1: y_q = b_q cos(πτ_q)e^{−iπτ_q} = B_q e^{−iπτ_q}
    — 165 §2a AYNEN. W_X < 1 iken arg y_q + πτ_q ARTIK SIFIR DEĞİLDİR
    (165 V3'ün τ=0.5'te 1.6 rad'a fırlayan artığının kalem karşılığı).

(TAŞIYICI)  e^{−iω_Q s_n}: hem ölçümde hem öngörüde AYNI sitede ⇒
    doğrudan bir sönüm çarpanı YOKTUR. Taşıyıcının τ_Q bağımlılığı
    yalnız GRAM SIZINTISI'ndan gelebilir (aşağıda λ).

(GRAM SIZINTISI — dördüncü, ölçülen çarpan)  Merdiven frekansları sonlu
    tarakta DİK DEĞİLDİR. Model alanı E = Σ_q Re[h'_q e^{iω_q s}] bir
    çizgi frekansına yansıtılınca komşularının gücünü de toplar:
        λ_E(ω) ≡ 2⟨E e^{−iωs}⟩ / h'_ω ,  λ_X(ω) ≡ 2⟨X e^{−iωs}⟩ / y_ω
    λ ≥ 1 ve çizgi yoğunluğuyla (yani τ ile) BÜYÜR. 163'ün V4'ü tam
    olarak `u0_pred = Re[h̄_Q·2⟨E e^{−iω_Q s}⟩]/|h_Q|² ≈ Re λ_E(ω_Q)`
    çarpanını 1'e çeker; 165'in `u2`-normalizasyonu aynı ailenin ikinci
    momentli üyesidir.
    ⇒ ÜÇ BACAKLI ÖNGÖRÜNÜN SIZINTI ÇARPIMI:  **λ_E(ω₁)·λ_X(ω₂)·λ_X(ω₃)**
    Baskın yapılandırmada (Ç1-B ailesi ve Ç4) bir X bacağı TAŞIYICIDA
    oturur (ω₂ = ω_Q), diğer ikisi merdiven boyunca toplanır ⇒
        kalib(τ_bant) ≈ 1 / [ λ_X(ω_Q) · Λ ] ,  Λ = bant-BAĞIMSIZ
    Bu, sıfır-parametreli bir BANT ŞEKLİ öngörüsüdür.

Bu script yukarıdaki BEŞ çarpanı da ÖLÇER (yaklaşımsız karakteristik
fonksiyonlarla, Gauss limiti yalnız karşılaştırma için yazılır).

Kullanım:  166_bacak.py <gaz> [taban] [tau_c]
Çıktı:     scratchpad/166/BACAK_<gaz>.npz (+ .json özeti)
"""
import importlib
import json
import sys
import time
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
for _p in ("166_configs", "165_configs", "163_configs", "160_configs",
           "159_configs", "155_configs", "154_configs"):
    sys.path.insert(0, str(QM / _p))
K = importlib.import_module("165_cekirdek")
C163 = importlib.import_module("163_cekirdek")

TWO_PI = 2 * np.pi
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/166")


def karakteristik(x, A):
    """⟨e^{iA x}⟩ — TAM (Gauss varsayımı YOK). A dizi, x dizi."""
    A = np.asarray(A, float)
    out = np.empty(len(A), dtype=complex)
    for i, a in enumerate(A):
        out[i] = complex(np.mean(np.cos(a * x)), np.mean(np.sin(a * x)))
    return out


def kos(veri, taban=0.40, tau_c=0.95):
    t0 = time.time()
    Y = K.gaz(veri, taban, 4000)
    Mo = K.Model165(veri, taban, 4000, 0.95, "olculen", Y=Y, tau_c=tau_c)
    Mo.sec("olculen", tau_c).alanlar(kmax=2)
    M = Mo.M
    tau, w = M["tau"], M["w"]
    msk = Mo.msk
    print(f"=== 166-BACAK {veri} τ_c={tau_c} ===")
    print(f"  σ_Ĉ={Mo.sigC:.5f} σ_ds={Mo.sigds:.5f} "
          f"σ_X̃={np.std(Y.Xtil0):.5f} çizgi={int(msk.sum())}", flush=True)

    # --- (1) üç karakteristik fonksiyon, TAM ---------------------------
    C = np.cumsum(Y.Xtil0)
    n = np.arange(len(C), dtype=float)
    n = (n - n.mean()) / (n[-1] if len(n) > 1 else 1.0)
    V = np.vstack([np.ones_like(n), n, n * n, n ** 3]).T
    c, *_ = np.linalg.lstsq(V, C, rcond=None)
    Chat = C - V @ c                      # TRENDSİZ birikmiş kayma
    ds, Xt = Y.ds, Y.Xtil0
    tg = np.arange(0.0, 1.0001, 0.005)    # τ ızgarası (rapor eğrileri)
    W_amp_g = karakteristik(ds, np.pi * tg).real          # A = πτ  (ds)
    W_X_g = karakteristik(Xt, -TWO_PI * tg)               # A = 2πτ (X̃)
    W_pos_g = karakteristik(Chat, TWO_PI * tg)            # A = 2πτ (Ĉ)

    # --- (2) çizgi-bazlı ölçülen/çıplak oranları -----------------------
    W_amp = karakteristik(ds, np.pi * tau).real
    W_X = karakteristik(Xt, -TWO_PI * tau)
    W_pos = karakteristik(Chat, TWO_PI * tau)
    b, B = M["b"], M["B"]
    hp, y = Mo.hp_olc, Mo.y_olc
    with np.errstate(divide="ignore", invalid="ignore"):
        R_h = np.where(b != 0, hp.real / b, np.nan)       # ↔ W_amp?
        R_y = np.where(B != 0, np.abs(y) / np.abs(B), np.nan)
    # kalemin y öngörüsü: ½ b W_amp (W_X e^{−2πiτ} + 1)
    y_kalem = 0.5 * b * W_amp * (W_X * np.exp(-2j * np.pi * tau) + 1.0)

    # --- (3) GRAM SIZINTISI λ_E, λ_X (merdivenin TAMAMINDA) -----------
    tm = time.time()
    hE, hX = K.tayf_s(Mo.s, [Mo.E, Mo.X], w)
    print(f"  λ tayfı {time.time()-tm:.0f}s", flush=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        lamE = np.where(np.abs(hp) > 0, hE / hp, np.nan)
        lamX = np.where(np.abs(y) > 0, hX / y, np.nan)

    # --- (4) τ kutularında özet ---------------------------------------
    kenar = np.arange(0.0, 1.0001, 0.04)
    ozet = []
    for i in range(len(kenar) - 1):
        m = (tau > kenar[i]) & (tau <= kenar[i + 1]) & msk
        if m.sum() < 3:
            continue
        gw = b[m] ** 2                     # çizgi gücü ağırlığı
        gw = gw / gw.sum() if gw.sum() else np.ones(m.sum()) / m.sum()
        ozet.append(dict(
            lo=float(kenar[i]), hi=float(kenar[i + 1]), n=int(m.sum()),
            tau=float((tau[m] * gw).sum()),
            R_h=float(np.nansum(R_h[m] * gw)),
            W_amp=float((W_amp[m] * gw).sum()),
            R_y=float(np.nansum(R_y[m] * gw)),
            R_y_kalem=float(np.nansum(
                (np.abs(y_kalem[m]) / np.abs(B[m])) * gw)),
            lamE=float(np.nansum(lamE[m].real * gw)),
            lamX=float(np.nansum(lamX[m].real * gw)),
            absE=float(np.nansum(np.abs(lamE[m]) * gw)),
            absX=float(np.nansum(np.abs(lamX[m]) * gw)),
            W_X=float((W_X[m].real * gw).sum()),
            W_pos=float((W_pos[m].real * gw).sum())))
    print("\n   τ    n    R_h=|h'|/b  W_amp   R_y=|y|/B  kalem   "
          "λ_E     λ_X    W_X    W_pos")
    for o in ozet:
        print(f"  {o['tau']:.3f} {o['n']:4d}  {o['R_h']:+8.4f} "
              f"{o['W_amp']:7.4f}  {o['R_y']:8.4f} {o['R_y_kalem']:7.4f} "
              f"{o['lamE']:7.4f} {o['lamX']:7.4f} {o['W_X']:6.3f} "
              f"{o['W_pos']:6.3f}", flush=True)

    # --- (5) bant frekanslarında λ (yarışın girdisi) -------------------
    ban = C163.bant_adaylari(Y, K.IZGARA_T1)
    bantlar = []
    qpos = {int(q): i for i, q in enumerate(M["q"])}
    for bd in ban:
        if not bd["cizgi"]:
            continue
        idx = np.array([qpos[int(r["q"])] for r in bd["cizgi"]])
        gp = np.array([r["gp"] for r in bd["cizgi"]])
        gw = gp / gp.sum()
        bantlar.append(dict(
            tau=bd["tau"], lo=bd["lo"], hi=bd["hi"], kul=bd["kul"],
            tau_gp=float((tau[idx] * gw).sum()),
            lamE=float(np.nansum(lamE[idx].real * gw)),
            lamX=float(np.nansum(lamX[idx].real * gw)),
            absE=float(np.nansum(np.abs(lamE[idx]) * gw)),
            absX=float(np.nansum(np.abs(lamX[idx]) * gw)),
            W_amp=float((W_amp[idx] * gw).sum()),
            W_X=float((W_X[idx].real * gw).sum()),
            W_pos=float((W_pos[idx].real * gw).sum()),
            R_h=float(np.nansum(R_h[idx] * gw)),
            R_y=float(np.nansum(R_y[idx] * gw))))

    # --- (6) merdiven-ortalamalı sızıntılar (Λ adayları) ---------------
    gw = (b ** 2 * msk)
    gw = gw / gw.sum()
    Lam = dict(
        LamE=float(np.nansum(np.where(np.isfinite(lamE.real), lamE.real, 0)
                             * gw)),
        LamX=float(np.nansum(np.where(np.isfinite(lamX.real), lamX.real, 0)
                             * gw)),
        W_amp_bar=float((W_amp * gw).sum()),
        W_X_bar=float((W_X.real * gw).sum()),
        W_pos_bar=float((W_pos.real * gw).sum()),
        varE_mod=Mo.artik["varE_mod"], varE_olc=Mo.artik["varE_olc"],
        varX_mod=Mo.artik["varX_mod"], varX_olc=Mo.artik["varX_olc"],
        gE=Mo.artik["gE"], gX=Mo.artik["gX"], gcal=Mo.artik["gcal"],
        EX_ic=float(np.dot(Mo.E, Y.e1 - Y.e1.mean()) / len(Mo.E)),
        sum_hp2=float(0.5 * np.sum(np.abs(hp[msk]) ** 2)),
        sum_b2=float(0.5 * np.sum(b[msk] ** 2)),
        sum_y2=float(0.5 * np.sum(np.abs(y[msk]) ** 2)),
        sum_B2=float(0.5 * np.sum(B[msk] ** 2)))
    print("\n  Λ (merdiven ortalamaları): " + " ".join(
        f"{k}={v:.5f}" for k, v in Lam.items()), flush=True)

    SCR.mkdir(parents=True, exist_ok=True)
    np.savez(SCR / f"BACAK_{veri}.npz", tau=tau, w=w, b=b, B=B, msk=msk,
             hp=hp, y=y, hE=hE, hX=hX, lamE=lamE, lamX=lamX,
             W_amp=W_amp, W_X=W_X, W_pos=W_pos, y_kalem=y_kalem,
             tg=tg, W_amp_g=W_amp_g, W_X_g=W_X_g, W_pos_g=W_pos_g)
    out = dict(veri=veri, taban=taban, tau_c=tau_c, sigC=Mo.sigC,
               sigds=Mo.sigds, sigX=float(np.std(Y.Xtil0)),
               sigChat=float(np.std(Chat)), ozet=ozet, bant=bantlar,
               Lam=Lam, sure_s=time.time() - t0)
    p = SCR / f"BACAK_{veri}.json"
    p.write_text(json.dumps(out, indent=1))
    print(f"-> {p}  ({(time.time()-t0)/60:.1f} dk)", flush=True)


if __name__ == "__main__":
    kos(sys.argv[1],
        float(sys.argv[2]) if len(sys.argv) > 2 else 0.40,
        float(sys.argv[3]) if len(sys.argv) > 3 else 0.95)
