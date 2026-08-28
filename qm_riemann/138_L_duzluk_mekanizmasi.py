"""
138 — −2'NİN L-DÜZLÜK MEKANİZMASI: PNT + τ-YALNIZ REGÜLATÖR (28 Ağu)
==========================================================================
137 ayrışımı mühürledi: κ_toplam(−2.00) = çekirdek(−3.05) + adyabatik
(+1.05). SORU: bu toplamı her L'de −2'de düz tutan ne?

KALEM (mekanizma adayı — iki L-iptali üst üste):
  (i)  PNT öz-benzerliği: kuyruk gücü τ-koordinatında
       dΣA² ≈ (2sinπτ)²·dτ/(π²τ)  [Σ_{q≤x}Λ(q)²/q ~ ½log²x'ten;
       log q = τL koyunca L düşer] — ölçü dτ/τ, L-BAĞIMSIZ.
  (ii) Regülatör: sönüm argümanı ω·σ_t = τL·σ_u·2π/L = 2πσ_u·τ —
       L kendini İPTAL eder; kesici yalnız τ'nun fonksiyonu.
  Ölçü de kesici de yalnız τ ⇒ çekirdek κ_k(τ) L-değişmez; κ_ad(τ)
  yerel kinematik ⇒ toplam yasa −2cosπτ L-değişmez. Kalıntı: σ_u(L)
  yavaş büyümesi → ±%3 (136'nın gördüğü). Düzlük Gauss-DW'ye özgü
  DEĞİL — her τ-yalnız regülatörde çıkar (değer ayrı, düzlük ayrı).

SINAV (üç zeros6 penceresi, L≈10.4 / 11.3 / 12.0):
  B1  κ_ad(τ*) ortak ızgarada (τ*=0.10/0.16/0.24) üç L'de ölçülür.
  B2  κ_arit(τ_p) her pencerede; fark(τ_p) = κ_arit − κ_ad(interp).
  B3  Ölçü kontrolü: τ-kutularında ΣA² / (PNT süreklisi) oranı, L'den
      bağımsız →1 mi?
ÖN-MÜHÜR:
  M1  κ_ad(τ*) ve fark(τ) eğrileri üç L'de ÇAKIŞIR (±%5) → düzlük =
      PNT öz-benzerliği + τ-yalnız regülatör; mekanizma MÜHÜRLENİR.
  M2  Çakışmaz (parça başına L-kayması, toplamda telafi) → konspirasyon;
      mekanizma başka yerde aranır.
Not: p=2 (τ≈0.06) farkı, κ_ad ızgarasının altına düşer — interp kenara
kıstırılır; fark karşılaştırması p≥3 üstünden okunur.

SONUÇ (28 Ağustos) — M1 DOĞRULANDI, MEKANİZMA MÜHÜRLENDİ:
  B3  Ölçü kontrolü: ΣA²/PNT oranı = 0.97-1.00, İKİ L'DE AYNI —
      kuyruk gücü τ-koordinatında evrensel dτ/τ ölçüsüne oturuyor. ✓
  B1  κ_ad ortak ızgara (erken/orta/son): 0.10→ +0.97/+0.94/+0.92;
      0.16→ +1.06/+0.99/+1.07; 0.24→ +1.34/+1.36/+1.34 — ÇAKIŞIK ±%4. ✓
  B2  fark(τ) üç pencerede TEK EĞRİ: τ=0.06'da ≈−2.87'den τ=0.25'te
      ≈−3.25'e yükselen evrensel çekirdek fonksiyonu K(τ); τ-eşleşmeli
      çapraz okumalar ±%2 (örn. −3.250@0.229 ↔ −3.245@0.224). ✓
  ⇒ −2'nin L-düzlüğü KONSPİRASYON DEĞİL: iki parça da ayrı ayrı
  L-değişmez τ-fonksiyonu. Mekanizma: (i) PNT öz-benzerliği (ölçü
  dτ/τ, L düşer) + (ii) τ-yalnız regülatör (ωσ_t = 2πσ_u·τ, L iptal).
  Düzlük herhangi bir τ-yalnız kesicide çıkar — Gauss-DW'ye özgü değil.
  BONUS: K(τ) ince yapısı artık ölçülü evrensel bir eğri (−2.87→−3.25,
  τ=0.06-0.25) — 132c'nin δ(τ) notunun gerçek-veri karşılığı; kapalı
  formun L=12 inceltmesi için hedef eğri bu.
"""

import numpy as np
from pathlib import Path
from sympy import primerange

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi

def pk_list(lim):
    out = []
    for p in primerange(2, int(lim) + 1):
        q = p
        while q <= lim:
            out.append(q); q *= p
    return sorted(set(out))

def chunked_fit(y, tmid, freqs, chunk=40000):
    tt = (tmid - tmid.mean()) / (tmid[-1] - tmid[0])
    C = 3 + 2 * len(freqs)
    XtX = np.zeros((C, C)); Xty = np.zeros(C)
    n = len(y)
    def cols(sl):
        c = [np.ones(sl.stop - sl.start), tt[sl], tt[sl]**2]
        for f in freqs:
            arg = f * tmid[sl]
            c += [np.cos(arg), np.sin(arg)]
        return np.vstack(c).T
    for s0 in range(0, n, chunk):
        sl = slice(s0, min(s0 + chunk, n))
        Xc = cols(sl)
        XtX += Xc.T @ Xc; Xty += Xc.T @ y[sl]
    b = np.linalg.solve(XtX, Xty)
    fit = np.empty(n)
    for s0 in range(0, n, chunk):
        sl = slice(s0, min(s0 + chunk, n))
        fit[sl] = cols(sl) @ b
    return b, fit

def coef(b, freqs, f):
    i = freqs.index(f)
    return b[3 + 2 * i], b[3 + 2 * i + 1]

def rnn_at(z, hedefler, ekstra=()):
    g = np.diff(z)
    m = 0.5 * (z[:-1] + z[1:])
    Lw = np.log(m / TWO_PI)
    L = float(Lw.mean())
    ds = g * Lw / TWO_PI - 1
    qs = pk_list(min(int(np.exp(0.52 * L)), 720))
    freqs = [np.log(q) for q in qs] + list(ekstra)
    b1, f1 = chunked_fit(ds, m, freqs)
    eta = ds - f1
    ee = eta[:-1] * eta[1:]
    c1 = float(ee.mean())
    mm = 0.5 * (m[:-1] + m[1:])
    b3, _ = chunked_fit(ee - c1, mm, freqs)
    out = {}
    for f in hedefler:
        tau = f / L
        cg, sg = coef(b1, freqs, f)
        A1 = np.hypot(cg, sg); ph = np.arctan2(sg, cg)
        Rn = (coef(b3, freqs, f)[0] * np.cos(ph) +
              coef(b3, freqs, f)[1] * np.sin(ph)) / (2 * c1 * A1)
        out[f] = (tau, Rn, Rn / np.cos(np.pi * tau))
    return out, L

d = np.load(HERE / "128_odl_zeros6_2e6_zeros.npz")
Z = np.sort(np.asarray(d["zeros"], dtype=float))
N = len(Z)
PENCERELER = {
    "erken":  Z[200000:500000],
    "orta":   Z[850000:1150000],
    "son":    Z[N - 300000:],
}
PR = [2, 3, 5, 7, 11, 13, 17]
TSTAR = [0.10, 0.16, 0.24]
LINES720 = [np.log(q) for q in pk_list(720)]

# B3 — ölçü kontrolü (kalem, veri gerektirmez): τ-kutularında
# Σ(2a sinπτ)² / [4sin²πτ/(π²τ)·Δτ] oranı, iki L'de
print("B3 — PNT ölçü kontrolü (oran → 1, L'den bağımsız mı?)")
QM = pk_list(1500000)
QOM = np.array([np.log(q) for q in QM])
QP = np.array([float(np.log(next(iter(__import__('sympy').factorint(q)))))
               for q in QM])
for Lc in (10.4, 12.0):
    tau = QOM / Lc
    a = QP / (np.pi * np.sqrt(np.exp(QOM)) * QOM)
    A2 = (2 * a * np.sin(np.pi * tau))**2
    print(f"  L={Lc}: ", end="")
    for t0 in (0.6, 0.8, 1.0):
        m = (tau >= t0) & (tau < t0 + 0.1)
        xs = np.linspace(t0, t0 + 0.1, 201)
        ys = 4 * np.sin(np.pi * xs)**2 / (np.pi**2 * xs)
        pnt = float((ys[:-1] + ys[1:]).sum() / 2 * (xs[1] - xs[0]))
        print(f"τ∈[{t0},{t0+0.1}): {A2[m].sum()/pnt:.3f}  ", end="")
    print(flush=True)

SONUC = {}
for ad, zz in PENCERELER.items():
    L0 = float(np.log(zz.mean() / TWO_PI))
    gbar = TWO_PI / L0
    print(f"\n=== pencere '{ad}': L={L0:.2f} ===", flush=True)
    out, L = rnn_at(zz, [np.log(p) for p in PR])
    kads = []
    for ts in TSTAR:
        om = ts * L
        while min(abs(om - l) for l in LINES720) < 0.012:
            om += 0.013
        U = 0.1 / (2 * np.sin(np.pi * om / L))
        zp = np.sort(zz + U * gbar * np.cos(om * zz))
        outp, _ = rnn_at(zp, [om], ekstra=[om])
        tau, Rn, kap = outp[om]
        kads.append((tau, kap))
        print(f"  boyalı τ*={tau:.3f}  κ_ad={kap:+.3f}", flush=True)
    ka_t = [t for t, k in kads]; ka_v = [k for t, k in kads]
    satirlar = []
    for p in PR:
        tau, Rn, kap = out[np.log(p)]
        kad = float(np.interp(tau, ka_t, ka_v))
        satirlar.append((p, tau, kap, kad, kap - kad))
        print(f"  p={p:>2} τ={tau:.3f}  κ_arit={kap:+.3f}  "
              f"fark={kap-kad:+.3f}", flush=True)
    SONUC[ad] = (L, kads, satirlar)

print("\nDEFTER — çakışma sınavı (M1: sütunlar L'den bağımsız mı?)")
print("κ_ad ortak ızgarada:")
print(f"  {'τ*':>5} " + " ".join(f"{ad:>8}" for ad in SONUC))
for i, ts in enumerate(TSTAR):
    print(f"  {ts:>5.2f} " + " ".join(
        f"{SONUC[ad][1][i][1]:>+8.3f}" for ad in SONUC))
print("fark = κ_arit − κ_ad (τ_p pencereye göre kayar; p etiketiyle):")
print(f"  {'p':>3} " + " ".join(f"{ad:>14}" for ad in SONUC))
for j, p in enumerate(PR):
    hcr = []
    for ad in SONUC:
        _, tau, kap, kad, frk = SONUC[ad][2][j]
        hcr.append(f"{frk:>+7.3f}@{tau:.3f}")
    print(f"  {p:>3} " + " ".join(f"{h:>14}" for h in hcr))
