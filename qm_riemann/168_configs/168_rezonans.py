"""
168 — A2(ii): Ç4 REZONANS-İNTEGRALİNİN SAYISAL YÜZLEŞMESİ (yalın)
=================================================================
167_rezonans.py'nin ★ kapalı formunu KULLANIR (yeniden türetmez), ama
GERÇEK-tarak J_2'yi (en pahalı parça, bekçiyi patlatan) HİÇ hesaplamaz:
c zaten ölçülmüş (C_*.json). Burada yalnız A1'in ORANI ölçülür.

★   Ç4(W) = Σ_r 𝒢_r [ Ĵ_D(W−ω_r) + Ĵ_D(W+ω_r) ] ,  𝒢_r = −πτ_r a_r cos(πτ_r)

A1'in iddiası:  model bu toplamı ÇIPLAK AYRIK olarak, tam ω_r noktalarında
değerlendirir.  Fiziksel karşılık ise profil-integralidir:

    c_pred = [ n_Δ(ω_r) · ∫κ_tepe ] / [ çıplak ayrık değerlendirme ]
           = ⟨Ĵ_D⟩_yerel(K·dres) / Ĵ_D(tam ω_r)

ÖLÇÜLENLER
  (P1) κ'nın TEPE PROFİLİ (gerçek tarak, ω_r çevresinde ince ızgara)
       → w_eff = ∫κ dΔ / κ(ω_r)  ve  ∫κ_tepe / (𝒢_r·2π/T)
  (P2) S_exact = Σ_r 𝒢_r[Ĵ_D(W∓ω_r)]  ↔  S_smear(K) (θ-ızgarasında K·dres
       genişlikli kutu ortalaması ≡ pencereyi T/K'ya kısaltmak)
       → oran(K) = A1'in c_pred'i, bant-birleştirilmiş
  (P3) ÇAPRAZ SINAV: merdiven noktaları ÖZEL mi?  ω_r → ω_r + δ_r
       (δ_r ∈ ±Δ rastgele) ile aynı toplam.

Kullanım: 168_rezonans.py <gaz> [nbant] [taban] [tau_c]
"""
import importlib
import json
import sys
import time
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
for _p in ("167_configs", "165_configs", "163_configs", "160_configs",
           "159_configs", "155_configs", "154_configs"):
    sys.path.insert(0, str(QM / _p))
K = importlib.import_module("165_cekirdek")
C163 = importlib.import_module("163_cekirdek")
ORT = importlib.import_module("167_ortak")

TWO_PI = 2 * np.pi
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/168")
K.PENCERE.update(ORT.pencere_dict())
OVS = 8
KSMEAR = (1, 2, 4, 8, 32)          # kutu YARI-genişliği, dres biriminde


class DTFT:
    """(1/N)Σ f_n e^{−iθn}; M = 2^k ≥ OVS·N ızgarada TAM, arada doğrusal."""

    def __init__(self, f=None, F=None, N=None):
        if F is None:
            self.N = len(f)
            self.M = 1 << int(np.ceil(np.log2(OVS * self.N)))
            self.F = np.fft.fft(np.asarray(f, float), self.M) / self.N
        else:
            self.F, self.N, self.M = F, N, len(F)

    def kutu(self, half_bins):
        """θ'da ±half_bins kutu ortalaması (≡ pencereyi kısaltmak)."""
        h = int(half_bins)
        if h < 1:
            return self
        F = self.F
        c = np.concatenate(([0], np.cumsum(np.concatenate([F, F[:2 * h + 1]]))))
        idx = np.arange(self.M)
        sm = (c[idx + 2 * h + 1] - c[idx]) / (2 * h + 1)
        return DTFT(F=np.roll(sm, h), N=self.N)

    def __call__(self, th):
        x = np.mod(np.asarray(th, float), TWO_PI) * self.M / TWO_PI
        i = np.floor(x).astype(np.int64)
        w = x - i
        return self.F[i % self.M] * (1 - w) + self.F[(i + 1) % self.M] * w


def kos(veri, nbant=3, taban=0.40, tau_c=0.95):
    t0 = time.time()
    SCR.mkdir(parents=True, exist_ok=True)
    Y = K.gaz(veri, taban, 4000)
    Mo = K.Model165(veri, taban, 4000, 0.95, "olculen", Y=Y, tau_c=tau_c)
    Mo.sec("olculen", tau_c)
    M = Mo.M
    gbar = TWO_PI / Y.L
    N = len(Mo.s)
    su = Mo.s[0] + gbar * np.arange(N)
    dres = Mo.dres
    print(f"=== 168-REZONANS {veri} taban={taban} τ_c={tau_c} ===", flush=True)
    print(f"  N={N} L={Y.L:.5f} ḡ={gbar:.6f} dres={dres:.4e} T={Mo.T:.1f} "
          f"çizgi={int(Mo.msk.sum())}  ({time.time()-t0:.0f}s)", flush=True)

    w = M["w"][Mo.msk]
    tau_r = M["tau"][Mo.msk]
    a_r = M["a"][Mo.msk]
    q_r = M["q"][Mo.msk]
    G = -np.pi * tau_r * a_r * np.cos(np.pi * tau_r)

    # ---------- (P1) κ'nın tepe profili, GERÇEK tarak ------------------
    tm = time.time()
    sel = [0, 1, 2, 4, 8, 20, 60, 300, 1500, 6000]
    sel = [i for i in sel if i < len(w)]
    JJ = np.arange(-30, 31)
    step = dres / 6.0
    nu = np.concatenate([w[i] + JJ * step for i in sel])
    kp = K.kappa(Mo.s, nu)
    P1 = []
    print(f"\n  (P1) κ TEPE PROFİLİ (gerçek tarak, ±5·dres, adım dres/6) "
          f"[{time.time()-tm:.0f}s]", flush=True)
    print("      q     τ_r      𝒢_r      Reκ(ω_r)   ∫Reκ dΔ      w_eff/dres"
          "   ∫Reκ/(𝒢·dres)   κ(±dres)/κ(0)")
    for k, i in enumerate(sel):
        pr = kp[k * len(JJ):(k + 1) * len(JJ)].real
        c0 = pr[len(JJ) // 2]
        integ = float(pr.sum() * step)
        i_p = len(JJ) // 2 + 6
        i_m = len(JJ) // 2 - 6
        P1.append(dict(q=int(q_r[i]), tau=float(tau_r[i]), G=float(G[i]),
                       k0=float(c0), integ=integ,
                       weff=integ / c0 / dres if c0 else float("nan"),
                       prof=[float(x) for x in pr]))
        print(f"  {int(q_r[i]):7d} {tau_r[i]:.4f} {G[i]:+.6f} {c0:+.6f} "
              f"{integ:+.3e}  {integ/c0/dres if c0 else float('nan'):9.3f}"
              f"    {integ/(G[i]*dres):9.3f}      "
              f"{0.5*(pr[i_p]+pr[i_m])/c0 if c0 else float('nan'):+.4f}",
              flush=True)

    # ---------- düzgün tarak alanları + DTFT ---------------------------
    tm = time.time()
    Ed = K.sentez(su, w, Mo.hp[Mo.msk])
    Xd = K.sentez(su, w, Mo.y[Mo.msk])
    Xd = Xd - Xd.mean()
    f2 = Ed * Xd * Xd
    D2 = DTFT(f2)
    DS = {k: D2.kutu(k * OVS) for k in KSMEAR}
    print(f"\n  düzgün tarak + DTFT + {len(KSMEAR)} kutu "
          f"({time.time()-tm:.0f}s) Var(Ed)={np.var(Ed):.5f} "
          f"Var(Xd)={np.var(Xd):.5f} M={D2.M}", flush=True)

    def JD(D, W):
        W = np.asarray(W, float)
        return np.exp(-1j * W * su[0]) * D(W * gbar)

    # ---------- (P2)+(P3) bant bant ------------------------------------
    ban = C163.bant_adaylari(Y, K.IZGARA_T1)
    ban = [b for b in ban if b["lo"] >= 0.52 - 1e-9 and b["cizgi"]][:nbant]
    Wall = []
    for b in ban:
        for r in b["cizgi"]:
            Wall.append(r["w"])
            Wall.append(r["w"] + r["gap"] / 2)
    Wall = np.array(Wall)
    tm = time.time()
    olc = [C163.olc_cizgi(Y, float(x), kmax=2) for x in Wall]
    print(f"  {len(ban)} bant, {len(Wall)} frekans; olc_cizgi "
          f"{time.time()-tm:.0f}s", flush=True)

    rng = np.random.default_rng(7)
    dlt = rng.uniform(-3 * dres, 3 * dres, size=len(w))     # (P3) kaydırma
    dlt2 = rng.uniform(-30 * dres, 30 * dres, size=len(w))

    tm = time.time()
    ETIK = ["exact"] + [f"K{k}" for k in KSMEAR] + ["shift3", "shift30"]
    Sv = {e: np.zeros(len(Wall), dtype=complex) for e in ETIK}
    for i, Wv in enumerate(Wall):
        Sv["exact"][i] = np.sum(G * (JD(D2, Wv - w) + JD(D2, Wv + w)))
        for k in KSMEAR:
            Sv[f"K{k}"][i] = np.sum(G * (JD(DS[k], Wv - w)
                                         + JD(DS[k], Wv + w)))
        Sv["shift3"][i] = np.sum(G * (JD(D2, Wv - w - dlt)
                                      + JD(D2, Wv + w + dlt)))
        Sv["shift30"][i] = np.sum(G * (JD(D2, Wv - w - dlt2)
                                       + JD(D2, Wv + w + dlt2)))
    print(f"  rezonans toplamları {time.time()-tm:.0f}s", flush=True)

    cikti = []
    print("\n  (P2/P3) BANT BANT — A²s2 biriminde Ç4 kapalı formu")
    print("   τ_eff   exact      " + "".join(f"K{k:<9d}" for k in KSMEAR)
          + "shift3     shift30")
    print("           (oran)     " + "".join(f"{'':10s}" for k in KSMEAR))
    for bi, b in enumerate(ban):
        Ls = b["cizgi"]
        i0 = sum(len(x["cizgi"]) for x in ban[:bi]) * 2
        pN = np.array([olc[i0 + 2 * i]["pow"] for i in range(len(Ls))])
        pO = np.array([olc[i0 + 2 * i + 1]["pow"] for i in range(len(Ls))])
        Ao = np.array([olc[i0 + 2 * i]["A"] for i in range(len(Ls))])
        Af = np.array([olc[i0 + 2 * i + 1]["A"] for i in range(len(Ls))])
        tv = np.array([r["tau"] for r in Ls])
        payda = float((pN - pO).sum())

        def agg(Jv):
            vN = np.array([K.s_den_J(olc[i0 + 2 * i]["h"], Jv[i0 + 2 * i],
                                     olc[i0 + 2 * i]["rho_ort"])[0]
                           for i in range(len(Ls))])
            vO = np.array([K.s_den_J(olc[i0 + 2 * i + 1]["h"],
                                     Jv[i0 + 2 * i + 1],
                                     olc[i0 + 2 * i + 1]["rho_ort"])[0]
                           for i in range(len(Ls))])
            return float((pN * Ao ** 2 * vN - pO * Af ** 2 * vO).sum() / payda)

        vals = {e: agg(Sv[e]) for e in ETIK}
        te = float((tv * (pN - pO)).sum() / payda)
        cikti.append(dict(tau_eff=te, **vals))
        e0 = vals["exact"]
        print(f"   {te:.4f} {e0:+.6f}  " + "  ".join(
            f"{vals[f'K{k}']:+.6f}" for k in KSMEAR)
            + f"  {vals['shift3']:+.6f}  {vals['shift30']:+.6f}")
        print("           oran:      " + "  ".join(
            f"{vals[f'K{k}']/e0:+.4f}    " for k in KSMEAR)
            + f"  {vals['shift3']/e0:+.4f}     {vals['shift30']/e0:+.4f}")

    out = dict(veri=veri, taban=taban, tau_c=tau_c, N=N, L=float(Y.L),
               T=Mo.T, dres=dres, ovs=OVS, ksmear=list(KSMEAR),
               P1=P1, P2=cikti, sure_s=time.time() - t0)
    p = SCR / f"REZ168_{veri}.json"
    p.write_text(json.dumps(out, indent=1))
    print(f"\n-> {p}  ({(time.time()-t0)/60:.1f} dk)", flush=True)


if __name__ == "__main__":
    kos(sys.argv[1],
        int(sys.argv[2]) if len(sys.argv) > 2 else 3,
        float(sys.argv[3]) if len(sys.argv) > 3 else 0.40,
        float(sys.argv[4]) if len(sys.argv) > 4 else 0.95)
