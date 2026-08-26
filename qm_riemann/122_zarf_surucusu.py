"""
122 — KARE-TERS ZARF SÜRÜCÜSÜ: "2× KİNEMATİK" ZİNCİR HAKEMİ (26 Ağu)
==========================================================================
Kalem zinciri: ν (yoğunluk-artığı) faz-hızı gürültüsüdür ve zarfla
ters-kare sürülür: Var(ν|E) ∝ E^{-4}; E² gap-dalgasına kinematik (+2)
biner (Not 1'in dalga hali) → R_ρ = +4 → R_g = 2−4 = −2.
"Kare bir kat daha ikiler."

ÖN-MÜHÜRLER:
  Z1  NOKTASAL (dalgasız) test: Var(η | a_u-desili) güçlü azalan;
      log-log eğim p ∈ [−2, −4] bandı (a_u = zarf vekili: gap-içi
      maks|Z|, unfold). Kontrol/dürüstlük: Var(η | g-desili) de basılır
      (birinci-moment g↔a bağının heteroskedastisite karışanı).
  Z2  R_w YENİDEN-NORMU: R_w' = P[ξ²]/(2σ_ξ²·A1_w) (genlik-kanalının
      KENDİ dalga genliğiyle). Zincir R_w'nun 1/τ-düşüşünü payda
      hatasına yorar → R_w' yaklaşık DÜZ ve O(+2) olmalı.
  Z3  Çapraz: Var(ξ | a_u) da basılır (genlik gürültüsünün kendi
      zarf-koşullaması — resmin iç tutarlılığı).
Veri: 41 (gaps+amps+tmid, iki pencere); η ve ξ tam-taban artıkları.

SONUÇ (26 Ağustos):
  Z1 — YÖN İSABET, ÜS RET: Var(η|a) tekdüze 6× düşüyor (0.054→0.009,
    d10'da hafif dönüş) AMA log-log eğim −0.60 (öngörü −2..−4 bandı
    DEĞİL). Kare-ters zarf zinciri üs düzeyinde ÇÜRÜDÜ; genlik-
    koşullaması gerçek ama yumuşak. (Kontrol: Var(η|ds-desil) de
    benzer düşüyor — paylaşılan yapı; iki-değişkenli koşullama
    [Var(η|a,g)] gelecek inceltme.)
  Z2 — BÜYÜK İSABET: doğru paydayla (A1_w) R_w' = +1.21→+0.72 —
    vahşi +6.3→+0.9 tamamen payda hatasıymış. GENLİK KANALI KENDİ
    DALGASINA ADYABATİK (+1) NEFES ALIYOR. Yeni keskin yapı:
    anomali GAP/YOĞUNLUK SEKTÖRÜNE HAS — |Z| sektörü öz-eş noktada
    (+1), gap sektörü −2'de. "Neden −2" sorusu sektör-seçici hale
    geldi: iki sektör aynı dalgaları taşırken biri denge-gibi, öteki
    ayna-noktasında.
  Z3 — genlik gürültüsü çarpımsal (Var(ξ|a) 0.09→0.27 yükselir) —
    doğal, resimle tutarlı.
  DURUM: −2'nin türetimi açık ama iki kez daralmış: (i) evrensel
  sabit (120), (ii) sektör-seçici (122). Sıradaki kalem: gap-sektörünü
  ayna-noktasına iten kısıt — sayım-özdeşliğinin kendisi mi (gap'ler
  ∫ρ=1 kısıtını TAM taşır, |Z| taşımaz — kısıtlı sektör anomalik)?
"""

import numpy as np
from pathlib import Path
from sympy import primerange

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
A_CG = (np.e**2 - 5) / 2
B0, B1 = 2.7580, -0.0543

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

d41 = np.load(HERE / "41_bigT_windows.npz")
PR = [2, 3, 5, 7, 11, 13, 17]
print("Z1/Z3 — desil tabloları ve Z2 — R_w yeniden-normu:")
RWp = {}
for k in ["120k", "200k"]:
    gz, am, tm = d41[f"gaps_{k}"], d41[f"amps_{k}"], d41[f"tmid_{k}"]
    Lw = np.log(tm / TWO_PI)
    L = float(Lw.mean())
    ds = gz * Lw / TWO_PI - 1
    a_u = am / np.sqrt(A_CG * Lw + B0 + B1 / Lw)
    a_u = a_u / np.sqrt((a_u**2).mean())
    qs = pk_list(min(np.exp(0.52 * L), 720))
    freqs = [np.log(q) for q in qs]
    bg, fg = chunked_fit(ds, tm, freqs)
    ba, fa = chunked_fit(a_u, tm, freqs)
    eta = ds - fg
    xi = a_u - fa
    # Z1: Var(η | a_u-desili) + kontrol Var(η | g-desili)
    d10 = np.quantile(a_u, np.linspace(0, 1, 11))
    print(f"\n[{k}] L={L:.2f}  Var(η|a_u-desil) [kontrol: Var(η|ds-desil)]:")
    xs, ys = [], []
    g10 = np.quantile(ds, np.linspace(0, 1, 11))
    for i in range(10):
        m1 = (a_u >= d10[i]) & (a_u < d10[i + 1])
        m2 = (ds >= g10[i]) & (ds < g10[i + 1])
        va = float((eta[m1]**2).mean())
        vg = float((eta[m2]**2).mean())
        am_ = float(a_u[m1].mean())
        xs.append(np.log(am_)); ys.append(np.log(va))
        print(f"  d{i+1}: ⟨a⟩={am_:5.2f}  Var(η|a)={va:.5f}   "
              f"[Var(η|ds-desil)={vg:.5f}]")
    p_fit = np.polyfit(xs[1:-1], ys[1:-1], 1)[0]
    print(f"  → log-log eğim (iç desiller): p = {p_fit:+.2f}  "
          f"(öngörü bandı −2..−4)")
    # Z3: Var(ξ | a_u)
    v_lo = float((xi[a_u < d10[3]]**2).mean())
    v_hi = float((xi[a_u > d10[7]]**2).mean())
    print(f"  Z3: Var(ξ|a düşük-3desil)={v_lo:.4f}  yüksek-3desil={v_hi:.4f}")
    # Z2: R_w', A1_w ile
    sx2 = float((xi**2).mean())
    b2, _ = chunked_fit(xi**2 - sx2, tm, freqs)
    for p in PR:
        f = np.log(p)
        cg, sg = coef(bg, freqs, f)
        A1g = np.hypot(cg, sg); ph = np.arctan2(sg, cg)
        ca, sa = coef(ba, freqs, f)
        A1w = np.hypot(ca, sa)
        cw, sw = coef(b2, freqs, f)
        P2 = cw * np.cos(ph) + sw * np.sin(ph)
        RWp.setdefault(p, []).append((P2 / (2 * sx2 * A1g),
                                      P2 / (2 * sx2 * A1w), A1w / A1g))
print(f"\nZ2 tablosu (iki pencere ort.):")
print(f"{'p':>3} {'R_w(eski)':>10} {'R_w(A1_w)':>10} {'A1_w/A1_g':>10}")
for p in PR:
    arr = np.array(RWp[p]).mean(axis=0)
    print(f"{p:>3} {arr[0]:>+10.3f} {arr[1]:>+10.3f} {arr[2]:>10.3f}")
