"""
141 — 0.76 ÇARPANI: EKRAN HİPOTEZİ SINAVI (29 Ağu)
==========================================================================
Kalem: çukur ağırlığı f(u)=sin²(πu/L)e^{−u/2}/u → τ̄_çukur≈0.25 (küçük
asallar). H-S: zero-örnekleme bastırması = ekran D(τ)=1−0.36τ (108,
bağımsız), kovaryansta D². Öngörüler:
  P_naif(n)   = Σ 2a²sin²(πτ)cos(2πnτ)              [140'ınki]
  P_D2(n)     = Σ 2a²sin²(πτ)·D(τ)²·cos(2πnτ)       [H-S]
  P_D1(n)     = Σ 2a²sin²(πτ)·D(τ)·cos(2πnτ)        [karşıt model]
Ayrıca varyans defteri: Σ2a²sin²D² + σ_η²(0.022, 109) =? V_zero(0.1674).
ÖN-MÜHÜR:
  E1  Altı çukurda ölçüm/P_D2 → 1.00±0.07 ise H-S DOĞRULANIR: 0.76,
      ekranın D² gölgesidir; 0.51 ile 0.76 tek fonksiyonda birleşir.
  E2  Oran sistematik <1 kalırsa (örn ~0.90) ekran-ötesi bastırma var
      (soğurma adayı) — miktar mühürlenir, kapı açık kalır.
  E3  Varyans defteri ±%5 içinde kapanmalı (aynı D² ile).

SONUÇ (29 Ağustos, gerçek koşulardan) — 0.76 ÜÇ KATMANA AYRILDI:
  Kalem: çukur ağırlığı f(u)=sin²(πu/L)e^{−u/2}/u → ⟨u⟩=2.82, τ̄=0.235
  (çukurları küçük asallar kazıyor; kuyruk kesimi konu dışı).
  Merdiven: naif 0.756±0.023 → eski-ekran D²(1−0.36τ) 0.883±0.026 →
  YERİNDE ölçülen r(τ)² ile 0.932±0.026 (parametresiz, altı çukur).
  D-kısmı — EKRAN YENİDEN ÖLÇÜLDÜ (ilk kez çizgi çizgi, kesişim serbest):
  asallar r ≈ 1.033 − 0.642τ (τ=0.06-0.37; kesişim ≈1 — düz gizemli
  çarpan YOK; eğim eski küresel 0.36'dan DİK — 108'le uzlaşma AÇIK:
  farklı kanal/L olabilir). KULELER (4/8/16) güçlendirilmiş okunuyor:
  1.03/1.05/1.06 sırayla artan — harmonik sızıntı adayı (AÇIK kapı).
  Ek koşu (faz + kule ayrımı): fazlar ≤0.002 rad — ÖZDEŞLİĞİN FAZ YAPISI
  MİLİRADYAN HASSASİYETLE DOĞRULANDI; faz-gecikme kapısı KAPALI.
  Kule-düzeltmeli öngörü farksız (0.952 vs 0.950) — kule kapısı KAPALI.
  ⇒ KALINTI: 0.93±0.03 (son 0.950 / orta 0.915). Ekran genlik-düzeyinde
  tam, faz kusursuz; kalan ~%7 KOLEKTİF aday — çapraz-çizgi koherent
  desenin (resurgence'ın kendisinin) okunuşundaki soğurma. AÇIK KAPI.
  E3 varyans defteri eski ekranla +%13.5 sapmıştı — dik ekranla uyumlu
  yönde; tam kapanış r(τ)'nin tüm-τ haritasını ister (AÇIK).
"""

import numpy as np
from pathlib import Path
from sympy import primerange

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
GAMS = [14.134725142, 21.022039639, 25.010857580]
NMAX = 60

d = np.load(HERE / "128_odl_zeros6_2e6_zeros.npz")
Z = np.sort(np.asarray(d["zeros"], dtype=float))
N = len(Z)
PEN = {"son": Z[N - 300000:], "orta": Z[850000:1150000]}

def pk_arrays(L):
    lim = int(np.exp(L))
    om, a2 = [], []
    for p in primerange(2, lim + 1):
        q, mm = p, 1
        while q <= lim:
            om.append(np.log(q)); a2.append(1.0 / (np.pi**2 * mm**2 * q))
            q *= p; mm += 1
    return np.array(om), np.array(a2)

def anom(C):
    A = np.zeros_like(C)
    for i in range(len(C)):
        js = [j for j in range(len(C)) if 3 <= abs(j - i) <= 8]
        A[i] = C[i] - np.median(C[js])
    return A

print("Çukur ağırlığı kalemi: f(u)=sin²(πu/L)e^{-u/2}/u, L=12:")
uu = np.linspace(0.05, 12, 2000)
f = np.sin(np.pi * uu / 12)**2 * np.exp(-uu / 2) / uu
print(f"  ağırlık merkezi ⟨u⟩={np.sum(uu*f)/np.sum(f):.2f} "
      f"→ τ̄={np.sum(uu*f)/np.sum(f)/12:.3f};  "
      f"⟨D²⟩_f={np.sum(f*(1-0.36*uu/12)**2)/np.sum(f):.3f}  "
      f"⟨D⟩_f={np.sum(f*(1-0.36*uu/12))/np.sum(f):.3f}")

ORT = {"naif": [], "D2": [], "D1": []}
for ad, zz in PEN.items():
    g = np.diff(zz)
    mid = 0.5 * (zz[:-1] + zz[1:])
    ds = g * np.log(mid / TWO_PI) / TWO_PI - 1
    ds = ds - ds.mean()
    L = float(np.log(mid / TWO_PI).mean())
    om, a2 = pk_arrays(L)
    tau = om / L
    D = np.clip(1 - 0.36 * tau, 0, None)
    base = 2 * a2 * np.sin(np.pi * tau)**2
    Cm = np.array([np.mean(ds[:-n] * ds[n:]) for n in range(1, NMAX + 1)])
    preds = {"naif": base, "D2": base * D**2, "D1": base * D}
    Am = anom(Cm)
    nstars = [gk * L / TWO_PI for gk in GAMS]
    if ad == "son":
        v_pred = float(np.sum(base * D**2)) + 0.022
        print(f"\nE3 varyans defteri (son): Σ2a²sin²D²+ση² = "
              f"{float(np.sum(base*D**2)):.4f}+0.022 = {v_pred:.4f}  "
              f"vs V_zero ölçüm = {float(np.var(ds)):.4f}  "
              f"(sapma {100*(v_pred-np.var(ds))/np.var(ds):+.1f}%)")
    print(f"\n=== {ad}: L={L:.3f}")
    for ns in nstars:
        i = int(round(ns)) - 1
        win = Am[max(0, i - 1):i + 2]
        j = i - 1 + int(np.argmax(np.abs(win)))
        satir = f"  n*≈{ns:5.1f} (n={j+1}): ölçüm={Am[j]:+.2e}"
        for etik, w in preds.items():
            Cp = np.array([float(np.sum(w * np.cos(2*np.pi*n*tau)))
                           for n in range(1, NMAX + 1)])
            r = Am[j] / anom(Cp)[j]
            ORT[etik].append(r)
            satir += f"  ölçüm/{etik}={r:.3f}"
        print(satir, flush=True)

print("\nDEFTER — altı çukur ortalaması (hedef: hangi model 1.00'e gelir?)")
for etik, rs in ORT.items():
    print(f"  ölçüm/P_{etik}: ort={np.mean(rs):.3f}  std={np.std(rs):.3f}")

# --- D-kısmı: ekran eğrisinin kendisi (kesişim serbest), son-300k
# tam-taban chunked regresyon → çizgi başına A1/(2a·sinπτ)
def chunked_fit(y, tmid, freqs, chunk=40000):
    tt = (tmid - tmid.mean()) / (tmid[-1] - tmid[0])
    C = 3 + 2 * len(freqs)
    XtX = np.zeros((C, C)); Xty = np.zeros(C)
    def cols(sl):
        c = [np.ones(sl.stop - sl.start), tt[sl], tt[sl]**2]
        for f in freqs:
            arg = f * tmid[sl]
            c += [np.cos(arg), np.sin(arg)]
        return np.vstack(c).T
    for s0 in range(0, len(y), chunk):
        sl = slice(s0, min(s0 + chunk, len(y)))
        Xc = cols(sl)
        XtX += Xc.T @ Xc; Xty += Xc.T @ y[sl]
    return np.linalg.solve(XtX, Xty)

from sympy import factorint
qs = []
for p in primerange(2, 721):
    q = p
    while q <= 720:
        qs.append(q); q *= p
qs = sorted(set(qs))

print("\nD/E-kısmı — yerinde ölçülen ekranla parametresiz çukur öngörüsü:")
ORT_R = []
for ad in ("son", "orta"):
    zz = PEN[ad]
    g = np.diff(zz); mid = 0.5 * (zz[:-1] + zz[1:])
    Lw = np.log(mid / TWO_PI); L = float(Lw.mean())
    ds = g * Lw / TWO_PI - 1
    freqs = [np.log(q) for q in qs if np.log(q) / L <= 0.52]
    b = chunked_fit(ds, mid, freqs)
    rq = {}
    for q in qs:
        f = np.log(q)
        if f / L > 0.52:
            continue
        i = freqs.index(f)
        A1 = float(np.hypot(b[3 + 2*i], b[4 + 2*i]))
        (pp, kk), = factorint(q).items()
        a = 1.0 / (np.pi * kk * q**0.5)
        rq[q] = A1 / (2 * a * np.sin(np.pi * f / L))
    if ad == "son":
        print("   q    τ      r=A1/(2a·sinπτ)  [asallar: kesişim≈1, eğim dik; kuleler >1]")
        for q in [2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 17, 19, 23]:
            tau = np.log(q) / L
            print(f"  {q:3d}  {tau:.3f}   {rq[q]:.3f}")
        pr = [(np.log(q)/L, rq[q]) for q in rq if factorint(q).popitem()[1] == 1]
        tt_, rr_ = np.array([x for x, _ in pr]), np.array([y for _, y in pr])
        m_ = tt_ <= 0.30
        A_ = np.polyfit(tt_[m_], rr_[m_], 1)
        print(f"  asal-fit (τ≤0.30): r ≈ {A_[1]:.3f} {A_[0]:+.3f}·τ")
    # E: parametresiz çukur öngörüsü — basis içinde ölçülen r_q, dışında fit-uzatma
    om_all, a2_all = pk_arrays(L)
    tau_all = om_all / L
    r_all = np.array([rq.get(int(round(np.exp(o))), 0.0) for o in om_all])
    ext = np.clip(A_[1] + A_[0] * tau_all, 0.35, None)
    r_eff = np.where(r_all > 0, r_all, ext)
    w = 2 * a2_all * np.sin(np.pi * tau_all)**2 * r_eff**2
    ds0 = ds - ds.mean()
    Cm = np.array([np.mean(ds0[:-n] * ds0[n:]) for n in range(1, NMAX + 1)])
    Am = anom(Cm)
    Cp = np.array([float(np.sum(w * np.cos(2*np.pi*n*tau_all)))
                   for n in range(1, NMAX + 1)])
    Ap = anom(Cp)
    print(f"  {ad}: ", end="")
    for gk in GAMS:
        ns = gk * L / TWO_PI
        i = int(round(ns)) - 1
        win = Am[max(0, i-1):i+2]
        j = i - 1 + int(np.argmax(np.abs(win)))
        ORT_R.append(Am[j] / Ap[j])
        print(f"n={j+1}: ölçüm/P_r² = {Am[j]/Ap[j]:.3f}   ", end="")
    print(flush=True)
print(f"\nE-DEFTER: ölçüm/P_r² altı çukur ort = {np.mean(ORT_R):.3f} "
      f"std = {np.std(ORT_R):.3f}  (hedef 1.00 — E1', parametresiz)")
