"""
167 — T3: Ç4'ÜN REZONANS İNTEGRALİ (kapalı form) ve SAYISAL SINAVI
==================================================================
TÜRETİM (kalem)
---------------
165'in tek-site çekirdeğinde  J_2(W) = ⟨E X² e^{−iWs}⟩  ve
κ(ν) = ⟨e^{iνs}⟩ olduğundan bütün üçlü toplam κ-ağırlıklıdır. Tarağın
KENDİSİ bir merdivendir: sayma fonksiyonu N̄(t)+S(t) = n ⇒ yoğunluk

    dN/dt = N̄'(t)·[1 + u(t)] ,   u(t) = 2 Σ_r 𝒢_r cos(ω_r t) ,
    𝒢_r = −π τ_r a_r cos(π τ_r)          (143'ün yasası; cos(πτ) yarım-gap
                                          orta-nokta form faktöründen)

Dolayısıyla HER pencere ortalaması ikiye ayrılır (`s` üzerinden ölçmek =
düzgün tarak + yoğunluk dalgalanması):

    ⟨f(s)⟩ = (1/T)∫f(t)[1+u(t)]dt
           = ⟨f⟩_düzgün + Σ_r 𝒢_r [ ⟨f e^{+iω_r t}⟩_düz + ⟨f e^{−iω_r t}⟩_düz ]

f = E X² e^{−iWt} koyunca **Ç4'ÜN KAPALI FORMU**:

    ★  J_2(W) = Ĵ_D(W) + Σ_r 𝒢_r [ Ĵ_D(W−ω_r) + Ĵ_D(W+ω_r) ] + O(u²)
       Ç4(W)  = Σ_r 𝒢_r [ Ĵ_D(W−ω_r) + Ĵ_D(W+ω_r) ]

Ĵ_D = DÜZGÜN tarak öngörüsü (165 §4b'nin kontrolü). Yani Ç4, düzgün-tarak
tayfının merdiven frekanslarında ÖRNEKLENMESİ ve 𝒢 ile tartılmasıdır —
"rezonans integrali" tam olarak budur. Tepe YÜKSEKLİĞİ 𝒢_r, tepe
PROFİLİ Ĵ_D'nin kendi Dirichlet çekirdeği, YOĞUNLUK merdivenin
d(sayı)/dω ≈ e^ω/ω'sidir.

HESAP HİLESİ (★'yı ucuzlatan): düzgün tarakta s_n = s_0 + ḡn olduğundan
    Ĵ_D(W) = ⟨f e^{−iWs}⟩_düz = e^{−iWs_0} · F(Wḡ) ,
    F(θ) = (1/N)Σ_n f_n e^{−iθn}            (N-noktalı DTFT)
F, f'nin M = 16N'e sıfır-doldurulmuş FFT'sinden TAM okunur; böylece
2·8981 kaydırılmış frekans **bedava** olur. (θ = Wḡ mod 2π ⇒ düzgün
tarakta ALIAS gerçektir: W ve W ± L aynı değeri verir. Bu bir kusur
değil, düzgün tarağın kendi fiziğidir.)

Kullanım: 167_rezonans.py <gaz> [taban] [tau_c] [nbant]
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
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/167")
K.PENCERE.update(ORT.pencere_dict())
OVS = 8           # FFT sıfır-doldurma katsayısı


class DTFT:
    """f_n dizisinin (1/N)Σ f_n e^{−iθn} DTFT'si; M = OVS·N ızgarada TAM,
    arasında doğrusal interpolasyon (bağıl hata ~ (1/OVS)²/8)."""

    def __init__(self, f, ovs=OVS):
        self.N = len(f)
        self.M = 1 << int(np.ceil(np.log2(ovs * self.N)))   # 2'nin kuvveti
        self.F = np.fft.fft(np.asarray(f, float), self.M) / self.N

    def __call__(self, th):
        x = np.mod(np.asarray(th, float), TWO_PI) * self.M / TWO_PI
        i = np.floor(x).astype(np.int64)
        w = x - i
        i0 = i % self.M
        i1 = (i + 1) % self.M
        return self.F[i0] * (1 - w) + self.F[i1] * w


def kos(veri, taban=0.40, tau_c=0.95, nbant=99):
    t0 = time.time()
    Y = K.gaz(veri, taban, 4000)
    Mo = K.Model165(veri, taban, 4000, 0.95, "olculen", Y=Y, tau_c=tau_c)
    Mo.sec("olculen", tau_c)
    M = Mo.M
    gbar = TWO_PI / Y.L
    N = len(Mo.s)
    su = Mo.s[0] + gbar * np.arange(N)
    print(f"=== 167-REZONANS {veri} taban={taban} τ_c={tau_c} ===")
    print(f"  N={N} L={Y.L:.5f} ḡ={gbar:.6f} dres={Mo.dres:.4e} "
          f"çizgi={int(Mo.msk.sum())}", flush=True)

    w = M["w"][Mo.msk]
    tau_r = M["tau"][Mo.msk]
    a_r = M["a"][Mo.msk]
    G_formul = -np.pi * tau_r * a_r * np.cos(np.pi * tau_r)

    # --- κ ölçümü (gerçek tarak) merdiven frekanslarında ---------------
    tm = time.time()
    kap_olc = K.tayf_s(Mo.s, [np.ones(N)], w)[0] / 2.0
    print(f"  κ(ω_r) ölçüldü ({time.time()-tm:.0f}s); "
          f"|Im|maks={np.abs(kap_olc.imag).max():.2e}", flush=True)

    # --- düzgün tarakta alanlar + DTFT'ler -----------------------------
    tm = time.time()
    Ed = K.sentez(su, w, Mo.hp[Mo.msk])
    Xd = K.sentez(su, w, Mo.y[Mo.msk])
    Xd = Xd - Xd.mean()
    f2 = Ed * Xd * Xd
    D1 = DTFT(np.ones(N))            # κ_D  (düzgün tarağın çekirdeği)
    D2 = DTFT(f2)                    # Ĵ_D
    print(f"  düzgün tarak + DTFT ({time.time()-tm:.0f}s) "
          f"Var(Ed)={np.var(Ed):.5f} Var(Xd)={np.var(Xd):.5f}", flush=True)

    def JD(W):
        W = np.asarray(W, float)
        return np.exp(-1j * W * su[0]) * D2(W * gbar)

    def KD(nu):
        nu = np.asarray(nu, float)
        return np.exp(1j * nu * su[0]) * np.conj(D1(nu * gbar))

    # --- (R1) κ'nın kapalı formu ---------------------------------------
    #     κ(Δ) = κ_D(Δ) + Σ_r 𝒢_r[κ_D(Δ−ω_r) + κ_D(Δ+ω_r)]
    r1 = []
    for j in range(min(12, len(w))):
        d = float(w[j])
        pk = complex(KD(np.array([d]))[0]
                     + np.sum(G_formul * (KD(d - w) + KD(d + w))))
        r1.append(dict(q=int(M["q"][Mo.msk][j]), tau=float(tau_r[j]),
                       kap_olc=[kap_olc[j].real, kap_olc[j].imag],
                       G=float(G_formul[j]), kap_kapali=[pk.real, pk.imag]))
    print("\n  (R1) κ(ω_r): ölçülen ↔ 𝒢 (143) ↔ kapalı form (★'nın κ'sı)")
    print("     q     τ_r    Re κ_ölç     𝒢         fark%    Re κ_kapalı  fark%")
    for r in r1:
        ko, kk = r["kap_olc"][0], r["kap_kapali"][0]
        print(f"  {r['q']:6d} {r['tau']:.4f} {ko:+.6f} {r['G']:+.6f} "
              f"{100*(r['G']-ko)/abs(ko):+7.2f}  {kk:+.6f} "
              f"{100*(kk-ko)/abs(ko):+7.2f}", flush=True)

    # --- (R1b) GRAM ŞİŞMESİNİN aynı rezonans integralinden gelmesi -----
    #   Λ = ⟨f⟩_gerçek/⟨f⟩_düzgün  ,  f = X_mod² ya da E_mod²
    #   ★ ⇒ ⟨f⟩_gerçek − ⟨f⟩_düz = Σ_r 𝒢_r[F̂_D(−ω_r) + F̂_D(+ω_r)]
    E = K.sentez(Mo.s, w, Mo.hp[Mo.msk])
    X = K.sentez(Mo.s, w, Mo.y[Mo.msk])
    X = X - X.mean()
    r1b = []
    print("\n  (R1b) Gram şişmesi ★'dan: ⟨f⟩_G − ⟨f⟩_D  ↔  Σ_r 𝒢_r[F̂_D(∓ω_r)]")
    print("     f          ⟨f⟩_G      ⟨f⟩_D      fark      kapalı(𝒢)   oran"
          "    kapalı(κ_ölç)  oran     Λ_ölç")
    for ad, fg, fd in (("X_mod²", X * X, (Xd - Xd.mean()) ** 2),
                       ("E_mod²", E * E, Ed * Ed),
                       ("E·X²", E * X * X, Ed * (Xd - Xd.mean()) ** 2)):
        Dk = DTFT(fd)
        g_ = float(np.mean(fg))
        d_ = float(np.mean(fd))
        sh = (np.exp(1j * w * su[0]) * np.conj(Dk(w * gbar))
              + np.exp(-1j * w * su[0]) * Dk(w * gbar))
        cG = float(np.sum(G_formul * sh).real)
        cK = float(np.sum(kap_olc.real * sh).real)
        r1b.append(dict(f=ad, G=g_, D=d_, fark=g_ - d_, kap_G=cG, kap_K=cK))
        print(f"    {ad:8s} {g_:+.6f} {d_:+.6f} {g_-d_:+.6f}  {cG:+.6f} "
              f"{cG/(g_-d_) if g_ != d_ else float('nan'):7.3f}  {cK:+.6f} "
              f"{cK/(g_-d_) if g_ != d_ else float('nan'):7.3f}  "
              f"{g_/d_ if d_ else float('nan'):7.4f}", flush=True)

    # --- (R2) Ç4'ün kapalı formu, bant bant ----------------------------
    ban = C163.bant_adaylari(Y, K.IZGARA_T1)
    ban = [b for b in ban if b["lo"] >= 0.52 - 1e-9 and b["cizgi"]][:nbant]
    Wall = []
    for b in ban:
        for r in b["cizgi"]:
            Wall.append(r["w"])
            Wall.append(r["w"] + r["gap"] / 2)
    Wall = np.array(Wall)
    tm = time.time()
    E = K.sentez(Mo.s, w, Mo.hp[Mo.msk])
    X = K.sentez(Mo.s, w, Mo.y[Mo.msk])
    X = X - X.mean()
    Jg = K.tayf_s(Mo.s, [E * X * X], Wall)[0] / 2.0        # GERÇEK tarak
    print(f"\n  gerçek tarak J_2 ({time.time()-tm:.0f}s)", flush=True)

    olc = [C163.olc_cizgi(Y, float(x), kmax=2) for x in Wall]
    cikti = []
    print("  (R2) Ç4 = J_gerçek − Ĵ_D  ↔  Σ_r 𝒢_r[Ĵ_D(W∓ω_r)]  (bant bant)")
    print("   τ_eff   A²s2(G)    A²s2(D)    Ç4_ölç     Ç4_kapalı(𝒢)  oran  "
          " Ç4_kapalı(κ_ölç)  oran   Ç4/G")
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

        Jd = JD(Wall)
        Jc_G = np.empty(len(Wall), dtype=complex)
        Jc_K = np.empty(len(Wall), dtype=complex)
        TRCUT = (0.20, 0.35, 0.50, 0.65, 0.80, 0.95)
        Jc_c = {x: np.zeros(len(Wall), dtype=complex) for x in TRCUT}
        for i in range(i0, i0 + 2 * len(Ls)):
            Wv = Wall[i]
            sh = JD(Wv - w) + JD(Wv + w)
            Jc_G[i] = np.sum(G_formul * sh)
            Jc_K[i] = np.sum(kap_olc.real * sh)
            for x in TRCUT:
                Jc_c[x][i] = np.sum(np.where(tau_r <= x, G_formul, 0.0) * sh)
        g = agg(Jg)
        d = agg(Jd)
        cG = agg(Jc_G)
        cK = agg(Jc_K)
        rec = dict(tau=b["tau"], lo=b["lo"], hi=b["hi"],
                   tau_eff=float((tv * (pN - pO)).sum() / payda),
                   G=g, D=d, C4=g - d, C4_kapali_G=cG, C4_kapali_K=cK,
                   kum_tr={str(x): agg(Jc_c[x]) for x in TRCUT})
        cikti.append(rec)
        print(f"   {rec['tau_eff']:.4f} {g:+.6f} {d:+.6f} {g-d:+.6f}  "
              f"{cG:+.6f}  {cG/(g-d) if g != d else float('nan'):6.3f}  "
              f"{cK:+.6f}  {cK/(g-d) if g != d else float('nan'):6.3f}  "
              f"{(g-d)/g if g else float('nan'):6.3f}", flush=True)
        print("      Ç4_kapalı(τ_r ≤ x) kümülatif: " + "  ".join(
            f"{x}:{rec['kum_tr'][str(x)]/cG if cG else float('nan'):.3f}"
            for x in TRCUT), flush=True)

    out = dict(veri=veri, taban=taban, tau_c=tau_c, N=N, L=float(Y.L),
               dres=Mo.dres, ovs=OVS, R1=r1, R2=cikti,
               sure_s=time.time() - t0)
    SCR.mkdir(parents=True, exist_ok=True)
    p = SCR / f"REZ_{veri}.json"
    p.write_text(json.dumps(out, indent=1))
    print(f"-> {p}  ({(time.time()-t0)/60:.1f} dk)", flush=True)


if __name__ == "__main__":
    kos(sys.argv[1],
        float(sys.argv[2]) if len(sys.argv) > 2 else 0.40,
        float(sys.argv[3]) if len(sys.argv) > 3 else 0.95,
        int(sys.argv[4]) if len(sys.argv) > 4 else 99)
