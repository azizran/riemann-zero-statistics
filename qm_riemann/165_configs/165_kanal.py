"""
165 — Ç3'ÜN |Δ| ÇÖZÜNÜRLÜĞÜ ve PENCERE ÇEKİRDEĞİ κ'NIN TANISI
==============================================================
`165_cekirdek` toplamı ALAN ÖZDEŞLİĞİ ile (budamasız) veriyor; bu modül
aynı toplamı ALT-MERDİVENDE terim terim sayarak

  (i) ÖZDEŞLİĞİ doğrular (V4: sayım = sentez, alt-pencerede bit düzeyi),
  (ii) Ç3'ün |Δ| dağılımını çıkarır (pencere-fazlı yakın-rezonanslar
       gerçekten nerede yaşıyor?),
  (iii) κ(ν)'yi ölçer: sinc zarfı, pencere-merkezi fazı ν·s̄, ve
        tarağın kendi rezonans tepeleri κ(±ω_q) ↔ 143'ün 𝒢'si.

κ'nın AYRIŞTIRILMASI (hesabın kilidi):
    κ(Δ) = ⟨e^{iΔ s_n}⟩ = e^{iΔ s̄} · κ̃(Δ),  κ̃(Δ) = ⟨e^{iΔ(s_n−s̄)}⟩
`arg κ` s̄ ≈ 1.05e6 yüzünden Δ ~ 1e−6'da tam tur döner (163'ün kutu
penceresinin gördüğü şey buydu); κ̃ ise dres = 2π/T ölçeğinde YUMUŞAKTIR.
Terim başına e^{iΔs̄} TAM uygulanır, κ̃ ise dres/20 aralıklı tablodan
doğrusal interpolasyonla okunur (ikinci-mertebe interpolasyon hatası
(1/20)²/8 ≈ 3e−04 bağıl).

Kullanım:  165_kanal.py <veri> [tau_s] [taban] [kaynak]
"""
import importlib
import json
import sys
import time
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
sys.path.insert(0, str(QM / "165_configs"))
sys.path.insert(0, str(QM / "163_configs"))
K = importlib.import_module("165_cekirdek")
C163 = importlib.import_module("163_cekirdek")

TWO_PI = 2 * np.pi
SCR = K.SCR


# ---------------------------------------------------------------------
def kappa_tablo(s, sbar, dmax, adim):
    """κ̃(Δ) = ⟨e^{iΔ(s−s̄)}⟩ tablosu, Δ ∈ [−dmax, dmax]."""
    n = int(np.ceil(dmax / adim))
    grid = np.arange(-n, n + 1) * adim
    ds = s - sbar
    tab = np.zeros(len(grid), dtype=complex)
    blok = max(500, int(2e7 / max(len(grid), 1)))
    for a in range(0, len(ds), blok):
        sl = slice(a, min(a + blok, len(ds)))
        arg = np.outer(ds[sl], grid)
        tab += np.cos(arg).sum(0) + 1j * np.sin(arg).sum(0)
        del arg
    return grid, tab / len(ds)


def kap(D, grid, tab, sbar, adim):
    """κ(Δ) — e^{iΔs̄} TAM, κ̃ tablodan interpolasyon."""
    f = (D - grid[0]) / adim
    i = np.clip(f.astype(np.int64), 0, len(grid) - 2)
    w = f - i
    kt = tab[i] * (1 - w) + tab[i + 1] * w
    return np.exp(1j * D * sbar) * kt


# ---------------------------------------------------------------------
def kanal_tara(Mo, iQ, W, tau_s, dcuts, adim_bol=20, alt_site=None):
    """Alt-merdivende (τ ≤ tau_s) |Δ| ≤ dcut terimlerini TAM sayar."""
    M = Mo.M
    tau, w = M["tau"], M["w"]
    sub = np.where(tau <= tau_s + 1e-12)[0]
    n2 = len(sub)
    # η çizgileri: regresyon τ ≤ taban'ı sildi ⇒ genlik ~0, yine de al
    We = np.concatenate((w[sub], -w[sub]))
    He = np.concatenate((Mo.hp[sub], np.conj(Mo.hp[sub])))
    Wx = We
    Ax = np.concatenate((Mo.y[sub], np.conj(Mo.y[sub])))
    # çift tayfı
    Fp = (Wx[:, None] + Wx[None, :]).ravel()
    Vp = (Ax[:, None] * Ax[None, :]).ravel()
    o = np.argsort(Fp)
    Fp, Vp = Fp[o], Vp[o]

    s = Mo.s if alt_site is None else Mo.s[alt_site]
    sbar = float(s.mean())
    dmax = max(dcuts) * 1.001
    adim = Mo.dres / adim_bol
    grid, tab = kappa_tablo(s, sbar, dmax, adim)

    D_all, C_all = [], []
    dcut = max(dcuts)
    for i1 in range(len(We)):
        lo = W - dcut - We[i1]
        hi = W + dcut - We[i1]
        a = int(np.searchsorted(Fp, lo))
        b = int(np.searchsorted(Fp, hi))
        if b <= a:
            continue
        D = We[i1] + Fp[a:b] - W
        C = He[i1] * Vp[a:b] / 8.0
        D_all.append(D)
        C_all.append(C)
    if not D_all:
        return dict(n=0)
    D = np.concatenate(D_all)
    C = np.concatenate(C_all)
    kv = kap(D, grid, tab, sbar, adim)
    terim = C * kv
    ad = np.abs(D)
    out = dict(n=int(len(D)), n_sifir=int((ad < 1e-12).sum()),
               J_sifir=complex(terim[ad < 1e-12].sum()))
    kum = []
    for c in dcuts:
        m = ad <= c
        kum.append(dict(dcut=float(c), dcut_dres=float(c / Mo.dres),
                        n=int(m.sum()), J=complex(terim[m].sum())))
    out["kum"] = kum
    return out


# ---------------------------------------------------------------------
def kappa_tani(Mo, nq=400):
    """κ'nın sinc zarfı ve tarak rezonans tepeleri."""
    s, sbar = Mo.s, Mo.sbar
    dres = Mo.dres
    # NOT: κ̃ düzgün tarakta Dirichlet çekirdeğidir ve Δ = k·dres'te
    # ÖZDEŞ SIFIRDIR; bu yüzden YARIM-tamsayı katlarında da örneklenir.
    D = np.array([j * dres / 8 for j in range(0, 81)])
    _, tab = kappa_tablo(s, sbar, D[-1] * 1.001, dres / 8)
    kt = tab[len(tab) // 2:]
    N = len(s)
    gbar = TWO_PI / Mo.Y.L
    Tn = N * gbar
    x = D * Tn / 2
    sinc = np.where(x == 0, 1.0, np.sin(x) / np.where(x == 0, 1.0, x))
    # tarak tepeleri: κ(ω_q)
    q = Mo.M["q"][:nq]
    wq = Mo.M["w"][:nq]
    kq = K.tayf_s(s, [np.ones(len(s))], wq)[0] / 2.0
    G = -np.pi * Mo.M["tau"][:nq] * Mo.M["a"][:nq] * np.cos(
        np.pi * Mo.M["tau"][:nq])
    return dict(D=D.tolist(), kt_abs=np.abs(kt).tolist(),
                kt_re=kt.real.tolist(), sinc=sinc.tolist(),
                q=q.tolist(), kq_abs=np.abs(kq).tolist(),
                kq_re=kq.real.tolist(), kq_im=kq.imag.tolist(),
                G=G.tolist(), Tn=float(Tn), T=float(Mo.T))


# ---------------------------------------------------------------------
def V4(Mo, iQ, W, tau_s=0.46, nsite=2000):
    """SAYIM = SENTEZ özdeşliği (alt-pencerede, budamasız)."""
    M = Mo.M
    tau, w = M["tau"], M["w"]
    sub = np.where(tau <= tau_s + 1e-12)[0]
    We = np.concatenate((w[sub], -w[sub]))
    He = np.concatenate((Mo.hp[sub], np.conj(Mo.hp[sub])))
    Ax = np.concatenate((Mo.y[sub], np.conj(Mo.y[sub])))
    s = Mo.s[:nsite]
    # sentez
    E = K.sentez(s, w[sub], Mo.hp[sub])
    X = K.sentez(s, w[sub], Mo.y[sub])
    J_sentez = complex(np.mean(E * X * X * np.exp(-1j * W * s)))
    # sayım (TAM, budamasız)
    Fp = (We[:, None] + We[None, :]).ravel()
    Vp = (Ax[:, None] * Ax[None, :]).ravel()
    J_sayim = 0j
    for i1 in range(len(We)):
        D = We[i1] + Fp - W
        kv = np.exp(1j * np.outer(D, s)).mean(1)
        J_sayim += complex(np.sum(He[i1] * Vp / 8.0 * kv))
    return J_sentez, J_sayim, len(We), len(Fp)


# ---------------------------------------------------------------------
if __name__ == "__main__":
    veri = sys.argv[1]
    tau_s = float(sys.argv[2]) if len(sys.argv) > 2 else 0.60
    taban = float(sys.argv[3]) if len(sys.argv) > 3 else 0.40
    kaynak = sys.argv[4] if len(sys.argv) > 4 else "olculen"
    t0 = time.time()
    Mo = K.Model165(veri, taban, 4000, 0.95, kaynak)
    Y = Mo.Y
    print(f"=== 165-KANAL {veri} tau_s={tau_s} kaynak={kaynak} ===")
    print(f"  N={len(Y.m0)} dres={Mo.dres:.4e} s̄={Mo.sbar:.1f} T={Mo.T:.1f}",
          flush=True)

    # --- κ tanısı ----------------------------------------------------
    kt = kappa_tani(Mo)
    print("  κ̃(Δ) ↔ sinc(ΔT_N/2):", flush=True)
    for j in (0, 2, 4, 6, 8, 12, 16, 20, 28, 40, 60, 80):
        print(f"    Δ/dres={kt['D'][j]/Mo.dres:6.2f}  |κ̃|={kt['kt_abs'][j]:.5f} "
              f"  sinc={kt['sinc'][j]:+.5f}   Re κ̃={kt['kt_re'][j]:+.5f}")
    print("  tarak tepeleri κ(ω_q) ↔ 𝒢 = −πτ a cos(πτ):", flush=True)
    for j in range(8):
        print(f"    q={kt['q'][j]:5d} |κ(ω_q)|={kt['kq_abs'][j]:.5f} "
              f"Re={kt['kq_re'][j]:+.5f} Im={kt['kq_im'][j]:+.5f}  "
              f"𝒢={kt['G'][j]:+.5f}")

    # --- bant çizgileri ----------------------------------------------
    ban = C163.bant_adaylari(Y, K.IZGARA_T1)
    qidx = {int(q): i for i, q in enumerate(Mo.M["q"])}
    dcuts = [Mo.dres * x for x in (0.0, 0.5, 1, 2, 3, 5, 10, 20, 50, 100,
                                   200)]
    dcuts[0] = 1e-13
    kayit = []
    for b in ban:
        if b["lo"] < 0.52 - 1e-9 or not b["cizgi"]:
            continue
        # bandın en güçlü 6 çizgisi (gp ağırlığı)
        Ls = sorted(b["cizgi"], key=lambda r: -r["gp"])[:6]
        for r in Ls:
            iQ = qidx[int(r["q"])]
            W = float(r["w"])
            o = C163.olc_cizgi(Y, W, kmax=2)
            hQ, rm = o["h"], o["rho_ort"]
            res = kanal_tara(Mo, iQ, W, tau_s, dcuts)
            Jc1 = Mo.kanal1(iQ)[0]
            # alt-merdiven sentez toplamı (referans)
            sub = np.where(Mo.M["tau"] <= tau_s + 1e-12)[0]
            E = K.sentez(Mo.s, Mo.M["w"][sub], Mo.hp[sub])
            X = K.sentez(Mo.s, Mo.M["w"][sub], Mo.y[sub])
            Jsub = complex(np.mean(E * X * X * np.exp(-1j * W * Mo.s)))
            rec = dict(q=int(r["q"]), tau=r["tau"], bant=b["tau"],
                       s2_olc=o["s2"],
                       s2_sub=K.s_den_J(hQ, Jsub, rm)[0],
                       s2_C1=K.s_den_J(hQ, Jc1, rm)[0],
                       s2_sifir=K.s_den_J(hQ, res["J_sifir"], rm)[0],
                       n_sifir=res["n_sifir"], n=res["n"],
                       kum=[dict(dres=c["dcut_dres"], n=c["n"],
                                 s2=K.s_den_J(hQ, c["J"], rm)[0])
                            for c in res["kum"]])
            kayit.append(rec)
            print(f"  q={rec['q']:6d} τ={rec['tau']:.4f} bant={b['tau']:.2f} "
                  f"s2ölç={rec['s2_olc']:+.5f} s2altmerd={rec['s2_sub']:+.5f} "
                  f"Ç1={rec['s2_C1']:+.5f} Δ≡0={rec['s2_sifir']:+.5f} "
                  f"(n0={rec['n_sifir']})", flush=True)
            print("     kümülatif |Δ|≤c·dres: " + "  ".join(
                f"{c['dres']:.0f}:{c['s2']:+.5f}" for c in rec["kum"]),
                flush=True)
            break       # bant başına bir çizgi (süre)

    # --- V4 ----------------------------------------------------------
    r0 = kayit[0] if kayit else None
    if r0:
        iQ = qidx[int(r0["q"])]
        W = float(Mo.M["w"][iQ])
        a, bq, nw, npair = V4(Mo, iQ, W, 0.44, 1000)
        print(f"  V4 (τ_s=0.44, {nw} işaretli çizgi, {npair} çift, 1000 site): "
              f"sentez={a:.12e} sayım={bq:.12e} "
              f"bağıl fark={abs(a-bq)/max(abs(a),1e-300):.2e}", flush=True)

    out = dict(veri=veri, tau_s=tau_s, kaynak=kaynak, dres=Mo.dres,
               kappa=kt, kayit=kayit,
               V4=dict(sentez=[a.real, a.imag], sayim=[bq.real, bq.imag],
                       bagil=abs(a - bq) / max(abs(a), 1e-300)) if r0 else None,
               sure_s=time.time() - t0)
    SCR.mkdir(parents=True, exist_ok=True)
    p = SCR / f"KAN_{veri}_ts{tau_s}_{kaynak}.json"
    p.write_text(json.dumps(out, indent=1, default=str))
    print(f"-> {p}  ({(time.time()-t0)/60:.1f} dk)", flush=True)
