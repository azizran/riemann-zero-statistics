"""
170 — K1: T(σ) AKTARIM FONKSİYONUNUN TÜRETİMİ + BACAK-BAŞINA σ_b ÖLÇÜMÜ
=======================================================================
Ölçüm/model parçası KOPYALANMAZ: `165_cekirdek.Model165/sentez/gaz`,
`163_cekirdek.bant_adaylari`, `169_k2.kuiper` AYNEN import edilir.

──────────────────────────────────────────────────────────────────────
A. T(σ)'NİN TÜRETİMİ  (kalemle; §A3 sayısal olarak DOĞRULAR)
──────────────────────────────────────────────────────────────────────
A1. TANIM. Bir bacağın "aktarımı", o bacağın alanı F yerine doyurulmuş
    karşılığı g(F) konduğunda korelatörün kazandığı çarpandır. Korelatör
    F'ye DOĞRUSAL olduğu ve öngörü varyans-eşleştirildiği için aktarım,
    g(F) ile F'nin KORELASYON KATSAYISIDIR:

        T ≡ ⟨g(F)·F⟩ / (σ_{g(F)} σ_F)                                  (T1)

    TAM doyumda g = sgn:  ⟨sgn(F)F⟩ = ⟨|F|⟩ = σ√(2/π), σ_{sgn} = 1 ⇒
        **T_∞ = √(2/π) = 0.797885**  — 169 §K2.3'ün bacak-başı sayısı.

A2. KISMİ DOYUM AİLESİ. Yumuşak doyum çekirdeği erf'tir (sgn'in Gauss
    yumuşatması):  g_α(F) = erf(F/(α σ_F √2)),  α → 0 ⇒ sgn.
    F ~ N(0,σ²) için (Stein + Gauss integralleri):

        ⟨F g_α⟩ = σ √(2/π) · √ρ ,        Var(g_α) = (2/π) arcsin ρ ,
        ρ ≡ 1/(1+α²) ∈ (0,1]

    ⇒  **T(ρ) = √( ρ / arcsin ρ )**                                    (T2)
        T(ρ→0) = 1        (doyum yok, doğrusal)
        T(ρ=1) = √(2/π)   (TAM doyum)                    ⇒ A1 ile TUTARLI

A3. İKİ BACAK ⇒ ARCSINE AİLESİ (169'un ara-değer eğrisi). İki bacak da
    kırpıldığında Gauss çiftinin işaret korelasyonu ⟨sgn x sgn y⟩ =
    (2/π)arcsin r, doğrusal karşılığı r ⇒ iki-bacak aktarımı

        **T₂(r) = (2/π) arcsin(r) / r**                                (T3)
        T₂(r→0) = 2/π = T_∞² ,  T₂(1) = 1

    (T2) ve (T3) AYNI arcsine ailesidir: her ikisi de ρ ↔ arcsin ρ
    oranıdır; (T3) 169 §K2.3'ün ölçtüğü nesnedir ve §B'de o beş nokta
    yeniden üretilir.

A4. σ ↔ ρ KÖPRÜSÜ (sarılmış-Gauss). Bacağın taşıdığı birikmiş faz
    ϑ_b = 2π τ_b Ĉ, ϑ_b ~ N(0, σ_b²). Sarılmış Gauss'un KOHERENT
    (temel-harmonik) genliği |⟨e^{iϑ}⟩| = e^{−σ_b²/2}; koherent GÜÇ payı
    e^{−σ_b²}. Geri kalan güç sarılmıştır ⇒ DOYMUŞ pay

        **ρ(σ_b) = 1 − e^{−σ_b²}**                                     (T4)

    ⇒ T(σ_b) = √( (1−e^{−σ_b²}) / arcsin(1−e^{−σ_b²}) ) ;
       σ_b→∞ ⇒ √(2/π) (KALEM'in şartı), σ_b→0 ⇒ 1.

A5. **H-D1'İN İŞARETİ (türetimin kendi öngörüsü, ölçümden ÖNCE).**
    T(σ) σ'da MONOTON AZALAN ve her yerde ≥ √(2/π)'dir. λ genlikleri
    küçültür ⇒ σ_b küçülür ⇒ T BÜYÜR ⇒ **c(λ) λ azaldıkça ARTMALIDIR
    ve 4/π²'nin ÜSTÜNDE kalmalıdır.** Ölçülen (169): c = 0.4035 →
    0.3817 → 0.3690, yani AZALIYOR ve 4/π²'nin ALTINA iniyor.
    Bu bir uyum sorunu değil, İŞARET sorunudur; K2 onu ön-kayıtla
    resmîleştirir.

──────────────────────────────────────────────────────────────────────
ÖN-MÜHÜR (koşudan ÖNCE yazılan sayısal tahminler)
──────────────────────────────────────────────────────────────────────
  * (T2)'nin Monte-Carlo doğrulaması: |MC − formül| < 0.003 (N = 4e6).
  * (T3), 169 §K2.3'ün beş bandını (%0.8–3.5) DÖRT HANEDE yeniden üretir.
  * σ_b = 2πτ_bσ_Ĉ:  Hkeskin'de ⟨τ⟩_E = 0.7054, ⟨τ⟩_X = 0.7166,
    σ_Ĉ = 0.27768 ⇒ σ_E ≈ 1.231, σ_X ≈ 1.250 rad; ρ ≈ 0.78 ⇒ T ≈ 0.934.
    ⇒ ΠT (4 bacak) ≈ 0.76 — ölçülen c = 0.4035'in NEREDEYSE İKİ KATI.
    Yani H-D1 λ=1.00'de zaten büyüklük olarak da ıskalamalıdır.
  * Tek bacak MOD 2π DÜZGÜN DEĞİLDİR (169 §K2.1: tek çizgide √N·V=78–500);
    yalnız dört-frekans toplamı düzgündür (σ_Σ = 2π·2.68·0.2777 ≈ 4.68 rad).
  * Alan düzeyi T_ölç = ⟨|F|⟩/σ_F: Hkeskin E ≈ 0.828, X ≈ 0.821 (169'un
    ölçtüğü sayılar); Edgeworth T = √(2/π)(1 − γ₂/24) bunları %1 içinde
    vermeli ve γ₂ < 0 (basıklık eksiği) çıkmalı.

SONUÇ bloğu YALNIZ gerçek koşu çıktısındandır.

Kullanım: 170b_T_yasasi.py [gaz1,gaz2,...] [çıktı_eki]
Çıktı:    scratchpad/170/K1_T<çıktı_eki>.json
"""
import importlib
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
for _p in ("169_configs", "167_configs", "166_configs", "165_configs",
           "163_configs", "160_configs", "159_configs", "155_configs",
           "154_configs"):
    sys.path.insert(0, str(QM / _p))
K = importlib.import_module("165_cekirdek")
C163 = importlib.import_module("163_cekirdek")
ORT = importlib.import_module("167_ortak")
KU = importlib.import_module("169_k2")            # kuiper() buradan

SCR = ("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
       "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/170")
TWO_PI = 2 * np.pi
LO_MIN, LO_MAX = 0.52, 0.68
T_INF = np.sqrt(2 / np.pi)
K.PENCERE.update(ORT.pencere_dict())


# ---------- (T2)/(T3)/(T4): kapalı formlar --------------------------
def T_rho(rho):
    """(T2)  T(ρ) = √(ρ/arcsin ρ)."""
    rho = np.clip(np.asarray(rho, float), 1e-12, 1.0)
    return np.sqrt(rho / np.arcsin(rho))


def T2_r(r):
    """(T3)  iki-bacak arcsine aktarımı."""
    r = np.asarray(r, float)
    return (2 / np.pi) * np.arcsin(r) / r


def rho_sigma(sig):
    """(T4)  sarılmış-Gauss doymuş güç payı."""
    return 1.0 - np.exp(-np.asarray(sig, float) ** 2)


def T_sigma(sig):
    return T_rho(rho_sigma(sig))


# ---------- A3: (T2)'nin Monte-Carlo doğrulaması --------------------
def dogrula_T2(n=4_000_000, tohum=1701):
    from scipy.special import erf
    rng = np.random.default_rng(tohum)
    F = rng.standard_normal(n)
    print("\n  --- A3: (T2)'nin MONTE-CARLO DOĞRULAMASI "
          f"(N = {n:,}, erf yumuşatması) ---")
    print("    α       ρ=1/(1+α²)   T_MC        T(ρ) formül    fark")
    out = []
    for al in (0.0, 0.25, 0.5, 1.0, 2.0, 4.0):
        g = np.sign(F) if al == 0 else erf(F / (al * np.sqrt(2)))
        Tmc = float(np.mean(g * F) / (np.std(g) * np.std(F)))
        rho = 1.0 / (1.0 + al ** 2)
        Tf = float(T_rho(rho))
        print(f"   {al:5.2f}   {rho:8.5f}   {Tmc:.6f}    {Tf:.6f}    "
              f"{Tmc-Tf:+.6f}")
        out.append(dict(alpha=al, rho=rho, T_mc=Tmc, T_form=Tf))
    return out


# ---------- B: 169'un arcsine ara-değer noktaları -------------------
def arcsine_169():
    """169 §K2.3'ün ÖLÇÜLEN beş noktası (K2b_Hkeskin.json'dan okunur)."""
    p = Path(SCR).parent / "169" / "K2b_Hkeskin.json"
    d = json.loads(p.read_text())
    print("\n  --- B: (T3) ↔ 169 §K2.3'ün ARCSINE ARA-DEĞER NOKTALARI ---")
    print("    bant lo  n_çizgi     r      ⟨sgn·sgn⟩/r (ÖLÇÜLEN)   "
          "T₂(r) = (2/π)arcsin r / r    fark%")
    out = []
    for a in d["arcsine"]:
        r, ss = a["r"], a["ss"]
        olc = ss / r
        pr = float(T2_r(r))
        print(f"    {a['lo']:.2f}    {a['n']:5d}   {r:+.5f}      "
              f"{olc:.5f}                 {pr:.5f}            "
              f"{100*(olc/pr-1):+6.2f}")
        out.append(dict(lo=a["lo"], n=a["n"], r=r, T2_olc=olc, T2_form=pr,
                        fark=100 * (olc / pr - 1)))
    return out


# ---------- C: bacak-başına σ_b ve alan-düzeyi aktarım --------------
def gaz_bacaklari(veri, taban=0.40, tau_c=0.95):
    t0 = time.time()
    Y = K.gaz(veri, taban, 4000)
    Mo = K.Model165(veri, taban, 4000, 0.95, "olculen", Y=Y, tau_c=tau_c)
    Mo.sec("olculen", tau_c).alanlar(kmax=2)
    msk = Mo.msk
    tl = Mo.M["tau"][msk]
    wh, wy = np.abs(Mo.hp[msk]), np.abs(Mo.y[msk])
    tE = float((tl * wh).sum() / wh.sum())
    tX = float((tl * wy).sum() / wy.sum())
    # katılım oranı (participation ratio) — alan şeklinin tek sayısı
    def npar(a):
        a2 = np.abs(a) ** 2
        return float(a2.sum() ** 2 / (a2 ** 2).sum())
    npE, npX = npar(Mo.hp[msk]), npar(Mo.y[msk])

    ban = C163.bant_adaylari(Y, K.IZGARA_T1)
    hed = [b for b in ban if LO_MIN - 1e-9 <= b["lo"] <= LO_MAX + 1e-9]
    tQ = float(np.mean([b["tau"] for b in hed]))
    sC = Mo.sigC

    # Ĉ (kübik-trendsiz) — sarılma testleri için
    C = np.cumsum(Y.Xtil0)
    n = np.arange(len(C), dtype=float)
    n = (n - n.mean()) / (n[-1] if len(n) > 1 else 1.0)
    V = np.vstack([np.ones_like(n), n, n * n, n ** 3]).T
    cc, *_ = np.linalg.lstsq(V, C, rcond=None)
    Chat = C - V @ cc
    sq = np.sqrt(len(Chat))

    def sarilma(tau):
        th = TWO_PI * tau * Chat
        return (float(abs(np.mean(np.exp(1j * th)))),
                float(KU.kuiper(th) * sq))

    bac = []
    for et, tau in (("taşıyıcı h_Q", tQ), ("E (h'_1)", tE),
                    ("X (y_2)", tX), ("X (y_3)", tX)):
        sb = TWO_PI * tau * sC
        R, V4 = sarilma(tau)
        bac.append(dict(bacak=et, tau=tau, sigma_b=sb, rho=float(rho_sigma(sb)),
                        T=float(T_sigma(sb)), R=R, kuiperV=V4))
    tS = tQ + tE + 2 * tX
    RS, VS = sarilma(tS)

    # ALAN DÜZEYİ: T_ölç = ⟨|F|⟩/σ_F, basıklık, Edgeworth
    alan = {}
    for et, F in (("E", Mo.E), ("X", Mo.X)):
        F = F - F.mean()
        s = float(np.std(F))
        Tm = float(np.mean(np.abs(F)) / s)
        g2 = float(np.mean((F / s) ** 4) - 3.0)
        alan[et] = dict(T_olc=Tm, gamma2=g2,
                        T_edge=float(T_INF * (1 - g2 / 24)),
                        n_p_kur=float(-1.5 / g2) if g2 < 0 else float("nan"),
                        sigma=s)
    alan["E"]["n_p_amp"] = npE
    alan["X"]["n_p_amp"] = npX

    lam_ = (ORT.KUNYE[veri].get("lam") or 1.0) if veri in ORT.KUNYE else 1.0
    r = dict(veri=veri, lam=lam_, sigC=sC, sigds=Mo.sigds,
             sigX=float(np.std(Y.Xtil0)), tau_Q=tQ, tau_E=tE, tau_X=tX,
             bacak=bac, tau_toplam=tS, R_toplam=RS, kuiperV_toplam=VS,
             sigma_toplam=TWO_PI * tS * sC, alan=alan,
             PiT=float(np.prod([b["T"] for b in bac])),
             sure_s=time.time() - t0)
    return r


def main(gazlar, ek=""):
    os.makedirs(SCR, exist_ok=True)
    print("=" * 116)
    print("K1 — T(σ) AKTARIM YASASI: TÜRETİM + BACAK-BAŞINA σ_b ÖLÇÜMÜ")
    print("=" * 116)
    print(f"  T_∞ = √(2/π) = {T_INF:.6f}   4/π² = {4/np.pi**2:.6f}")
    mc = dogrula_T2()
    ars = arcsine_169()

    print("\n  --- C: BACAK-BAŞINA σ_b = 2πτ_b σ_Ĉ  ve  T(σ_b) ---")
    R = {}
    for g in gazlar:
        r = gaz_bacaklari(g)
        R[g] = r
        print(f"\n  [{g}]  λ={r['lam']}  σ_Ĉ={r['sigC']:.5f}  "
              f"σ_ds={r['sigds']:.5f}  σ_X̃={r['sigX']:.5f}   "
              f"({r['sure_s']:.0f}s)")
        print("    bacak          τ_b      σ_b(rad)   ρ=1−e^{−σ²}   T(σ_b)   "
              "R(τ_b)     √N·V (mod 2π düzgün mü?)")
        for b in r["bacak"]:
            duz = "DÜZGÜN" if b["kuiperV"] < 2.0 else "DEĞİL"
            print(f"    {b['bacak']:13s} {b['tau']:.4f}   {b['sigma_b']:7.4f}"
                  f"    {b['rho']:.5f}     {b['T']:.5f}  {b['R']:.5f}  "
                  f"{b['kuiperV']:8.2f}  {duz}")
        duzS = "DÜZGÜN" if r["kuiperV_toplam"] < 2.0 else "DEĞİL"
        print(f"    Ç4 TOPLAM     {r['tau_toplam']:.4f}   "
              f"{r['sigma_toplam']:7.4f}    "
              f"{1-np.exp(-r['sigma_toplam']**2):.5f}     "
              f"{float(T_sigma(r['sigma_toplam'])):.5f}  {r['R_toplam']:.5f}  "
              f"{r['kuiperV_toplam']:8.2f}  {duzS}")
        print(f"    ⇒ H-D1'in ÇARPIMI  Π_b T(σ_b) = {r['PiT']:.5f}"
              f"      [4/π² = {4/np.pi**2:.5f}]")
        print("    ALAN DÜZEYİ (varyans-eşleşmiş kırpma):")
        for et in ("E", "X"):
            a = r["alan"][et]
            print(f"      {et}: T_ölç=⟨|F|⟩/σ = {a['T_olc']:.5f}   "
                  f"γ₂ = {a['gamma2']:+.4f}   Edgeworth √(2/π)(1−γ₂/24) = "
                  f"{a['T_edge']:.5f}  (fark {100*(a['T_olc']/a['T_edge']-1):+.2f}%)"
                  f"   n_p(genlik)={a['n_p_amp']:.2f}")

    # λ ekseninde H-D1'in öngördüğü YÖN
    print("\n" + "=" * 116)
    print("K1(d) — H-D1'İN λ EKSENİNDEKİ YÖNÜ (türetimden, ölçüme bakmadan)")
    print("=" * 116)
    print(f"  {'gaz':9s} {'λ':>5s} {'σ_Ĉ':>8s} {'σ_E':>7s} {'σ_X':>7s} "
          f"{'Π_b T(σ_b)':>11s}  {'Π T / (4/π²)':>12s}")
    for g in gazlar:
        r = R[g]
        sE = [b for b in r["bacak"] if b["bacak"].startswith("E")][0]["sigma_b"]
        sX = [b for b in r["bacak"] if b["bacak"].startswith("X")][0]["sigma_b"]
        print(f"  {g:9s} {r['lam']:5.2f} {r['sigC']:8.5f} {sE:7.4f} {sX:7.4f} "
              f"{r['PiT']:11.5f}  {r['PiT']/(4/np.pi**2):12.4f}")

    json.dump(dict(T_inf=float(T_INF), mc=mc, arcsine=ars, gaz=R),
              open(SCR + f"/K1_T{ek}.json", "w"), indent=1, default=float)
    print(f"\n-> {SCR}/K1_T{ek}.json")


if __name__ == "__main__":
    gz = (sys.argv[1].split(",") if len(sys.argv) > 1
          else ["Hkeskin", "L085", "L070", "son"])
    main(gz, sys.argv[2] if len(sys.argv) > 2 else "")
