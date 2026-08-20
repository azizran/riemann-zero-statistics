"""
102b — KANAT ÖNGÖRÜSÜ: 87-T3'ÜN "TDS TAŞIYICILARI" GERÇEK Mİ? (20 Ağustos)
==========================================================================
SEFERİN KALBİ. 87-T3, giydirmeyi taşıyan çiftlerin (Δω, Σω)
frekanslarında |Ĝ| = 0.002–0.009 ölçtü; "p^k çizgilerinde DEĞİL, ama
yerel karanlık tabanın (5e-4) 4–15 KATI, ω arttıkça tarağa doğru
yükseliyor" dedi ve kimliği TERMAL DİFÜZ SAÇILMA (TDS) koydu. 101j o
tabanın TAPERSİZ ölçümde saf pencere sızıntısı olduğunu gösterdi.

HİPOTEZ (test edilen): çizgi-dışı frekanslardaki TÜM Ĝ içeriği,
komşu çizgilerin PENCERE ÇEKİRDEĞİ KANATLARI (+ üst mertebe yan
bantları) ile açıklanır; ayrı bir "fonon difüz" bileşeni gerekmez.

BİRİM: I(ω) = |Σ w e^{iω t}|²/Σw²  (Poisson→1, w=1 iken n|Ĝ|²,
bin-ortalaması Montgomery rampası α). |Ĝ| = √(I/n).

PARAMETRESİZ MODELLER (hiçbir uydurma yok):
  M-K  KANAT ÇEKİRDEĞİ (birinci mertebe, tutarsız toplam):
       I_kanat(ω) = I_düz(ω) + Σ_q |A_q|²·[I_düz(ω−log q) + I_düz(ω+log q)]
       I_düz = RvM pürüzsüz ızgaranın ÖLÇÜLEN I'sı (= n_eff|K|², yani
       pencerenin kendi çekirdeği: DC tepesi + tarak, ölçümle geliyor);
       A_q = Λ(q)/(L√q)·cos(πτ)·e^{−ω²σ_t²/2}  (89'un MUTLAK benek
       yasası — kanal ölçümü yok, serbest parametre yok).
  M-L  LABORATUVAR ÖRGÜSÜ (tüm mertebeler): pürüzsüz sıfır ızgarasına
       açık formülün yerdeğiştirmesi u(t) = (2/L(t))Σ_q Λ(q)/(√q log q)
       ·sin(t log q) uygulanır (q ≤ Q), ORTA NOKTALARI alınır ve AYNI
       boru hattından geçirilir. Gürültü yok, ζ yok, sadece asallar +
       pencere optiği. İkinci/üçüncü mertebe yan bantları (log(q₁q₂)…)
       otomatik içerir.
  M-S  SIZINTI TABANI: RvM pürüzsüz ızgara tek başına (DC + tarak
       kanadı) — "hiç aritmetik yokken alet ne gösterir".

ÖLÇÜMLER
  B1  87'nin T2/T3'ü birebir yeniden koşulur (iki pencere); en büyük 8
      taşıyıcı çiftin Δω, Σω frekansları ÇIKARILIR (87'nin gerçek
      örneklem noktaları).
  B2  O frekanslarda: ölçülen I (tapersiz/Hann), M-S, M-K, M-L; oran
      tablosu. 87'nin "yerel taban"ı (±0.04 medyanı) da aynı şekilde.
  B3  ÇİZGİ MESAFESİ DENETİMİ: 87 "çizgide değil" hükmünü q≤128'lik
      EKSİK bir listeyle kurmuştu. Tam asal-kuvvet listesiyle (q ≤ e^{11})
      en yakın çizgi mesafesi, pencere çözünürlüğü 2π/T birimlerinde.
  B4  ORMAN TESTİ: taşıyıcı frekansı çevresinde ±0.02'lik İNCE tarama
      (adım ≈ çözünürlük/4), Hann taperli. Sürekli bir difüz omuz mu,
      yoksa ÇÖZÜLMÜŞ ÇİZGİ ORMANI mı? Aynı tarama M-L'de de yapılır.
  B5  SICAKLIK KAPISI: laboratuvar örgüsünün σ_t'si gerçek pencereninki
      ile uyuşuyor mu (yani doğru "ısıda" mı karşılaştırıyoruz)?

ÖN-MÜHÜR (ölçümden ÖNCE yazıldı):
  H1 87'nin taşıyıcı frekansları, TAM çizgi listesiyle bakıldığında
     çizgilerden ~çözünürlüğün 10–100 katı uzakta çıkacak (yani "çizgide
     değil" hükmü teknik olarak doğru), AMA çevrelerindeki çizgi
     ORMANI o kadar sık olacak ki (Δlog p ≈ log p/p ~ 0.01) kanat
     toplamı kaçınılmaz.
  H2 Tapersiz ölçüm ile M-S (saf sızıntı) 87'nin "yerel taban"ını
     (5e-4) büyük ölçüde açıklayacak → "4–15×" oranının PAYDASI alet.
  H3 Taşıyıcı düzeyi (0.005 civarı, I≈1) taperle ÇÖKMEYECEK: ω≈5–9'da
     α = ω/L ≈ 0.5–0.9, yani Montgomery rampası zaten bu düzeyde; içerik
     GERÇEK. Ama kimliği "fonon TDS" değil, sık ÇİZGİ ORMANI olacak
     (B4 çözecek) — yani 87'nin GİYDİRME mekanizması ayakta, TAŞIYICI
     KİMLİĞİ değişecek.
  H4 M-L, taşıyıcı frekanslarındaki düzeyi mertebe içinde (0.3–3×)
     yeniden üretecek; M-K (yalnız birinci mertebe kanat) daha düşük
     kalabilir (üst mertebe yan bantları eksik).
  RİSK: H3 yanlışsa (taperle çöküyorsa) 87'nin T3 hükmü tamamen düşer;
  H3 doğru ama B4 sürekli omuz gösterirse TDS kimliği AYAKTA kalır.

==========================================================================
SONUÇ (20 Ağustos, koşu ~1.5 dk) — H1 ✓, H2 RET, H3 ✓✓, H4 ✓(M-L)/RET(M-K)
KANAT HİPOTEZİ REDDEDİLDİ; İÇERİK GERÇEK VE TAMAMEN ARİTMETİK — AMA
SIFIR ÖRGÜSÜNÜN DEĞİL, ORTA-NOKTA ÖRNEKLEMESİNİN ÜRÜNÜ.
==========================================================================
B1 ✓ 87 birebir yeniden üretildi. Taşıyıcılar: Δω = 3.80–4.76
  (α = 0.37–0.46) ve Σω = 6.14–7.02 (α = 0.59–0.68). Ölçülen |Ĝ|
  (tapersiz) = 0.0016–0.0073 — 87'nin "0.002–0.009"u ✓.

B3 ✓ (H1): 87 "çizgide değil (uzaklık 0.03–0.44)" derken q ≤ 128'lik
  EKSİK bir liste kullanmıştı. Tam listeyle (7327 asal-kuvvet, log q ≤
  11.2) uzaklıklar 0.0018–0.049, yani ÇÖZÜNÜRLÜĞÜN (2π/T = 2.6e-4)
  7–199 KATI. Çizgi-üstü değiller (doğru), ama çevrelerindeki çizgi
  ormanı ω≈6.5'te Δlog p ≈ 0.010 aralıkla oturuyor.

B2 — SEVİYELER (I = n|Ĝ|² birimi, 200k / 120k medyanları):
  ölçülen (tapersiz)          0.367 / 0.604      [rampa α ≈ 0.63/0.64]
  ölçülen (Hann)              0.241 / 0.340
  taperli/tapersiz            0.67  / 0.62   → İÇERİK GERÇEK (H3 ✓✓)
  ölçülen / M-S (saf sızıntı) 3.4e4 / 7.0e4 (ham); 1.3e20 / 4.9e20 (Hann)
  ölçülen / M-K (1. mertebe kanat) 2291× / 3588×  → KANAT HİPOTEZİ RET
  ölçülen / M-L (laboratuvar örgüsü) 1.06 / 1.07 (ham); 0.96 / 1.12 (Hann)
  ★ M-L: yalnızca AÇIK FORMÜLDEN kurulmuş, ζ verisi görmemiş, gürültüsüz
  bir örgü, 87'nin taşıyıcı düzeyini MEDYANDA %6 içinde veriyor.
  (Nokta-nokta saçılma 0.48–4.7×, geometrik ~2×: kimlik "düzey ve
  mekanizma" düzeyinde, nokta-nokta değil — dürüst kayıt.)
  B5 SICAKLIK KAPISI: σ_lab/σ_gerçek = 1.216 (200k), 0.998 (120k) ✓.

  87'NİN "YEREL KARANLIK TABAN"I (±0.04 medyanı, Δω ≈ 4.3):
  ölçülen |Ĝ| = 5.5e-4 (87: 5e-4 ✓), I = 0.0121 ham / 0.0075 Hann;
  M-S saf sızıntı I ≈ 1e-5 → TABAN SIZINTI DEĞİL, GERÇEK (H2 RET).
  Dolayısıyla 87'nin "taban üstü 4–15×" oranı AYAKTA (Δω=4.309'da
  0.00218/0.00055 = 4.0×). 101j'nin "87'nin tabanı sızıntı düzeyinde"
  şüphesi ω ≈ 4–7 bandı için YANLIŞ çıktı; yalnız ω < 1 için doğru
  (102a).

B4 — ★ ORMAN TESTİ (ω = 7.017 ± 0.02, adım = çözünürlük/4, Hann):
  ızgara      medyan I    ort I    maks    çizgi-üstü  çizgi-dışı
  gerçek ORTA  7.0e-02   3.8e-01   3.37      1.17        0.220
  SIFIRLAR     9.1e-06   7.2e-01  11.10      4.20        4.3e-04
  M-L lab      2.6e-01   6.0e-01   4.73      1.37        0.442
  M-S düz      5.8e-21   8.5e-21   4.5e-20      —           —
  SIFIR tayfı: ÇÖZÜLMÜŞ ÇİZGİ ORMANI + GERÇEK BOŞLUK (maks/medyan =
  1.2e6; çizgi-üstü/çizgi-dışı ≈ 10⁴). ORTA-NOKTA tayfı: aynı orman,
  ama boşluklar 0.22 düzeyinde DOLDURULMUŞ (kontrast yalnız 5.3×).

★ IZGARA AYRIMI (aynı taşıyıcı frekanslarında, Hann, medyan I):
  ORTA NOKTALAR t_n = z_n + g_n/2 : 0.241 (200k) / 0.340 (120k)
  SIFIRLAR      z_n               : 3.7e-09      / 7.9e-10
  oran                             : 9.1e+07×     / 2.6e+08×
  Yani 87'nin "termal difüz saçılma" dediği içeriğin ~%100'ü ORTA-NOKTA
  İNŞASININ ürünüdür (102a'nın örnekleme-fazı eğrisiyle aynı mekanizma:
  Ĝ_orta = (1/n)Σ e^{iωz_n}e^{iωg_n/2}, yarım-gap kayması doğrusal
  olmayan yan bantlar üretir).

GENEL HÜKÜM:
  1) 87'nin GİYDİRME mekanizması AYAKTA: Gram=Ĝ kimliği kesin, transfer
     cebiri kesin, taşıyıcı düzeyleri GERÇEK (taperi geçiyor, sızıntının
     10⁴–10²⁰ katı üstünde, laboratuvar örgüsüyle %6 içinde).
  2) TAŞIYICININ KİMLİĞİ DEĞİŞİYOR: "sıfır kristalinin termal difüz
     saçılması" DEĞİL; "asal çizgilerinin, orta-nokta örneklemesinin
     dalgalanan yarım-gap kaymasıyla ürettiği DOĞRUSAL-OLMAYAN YAN
     BANTLARI". Sıfır örgüsünün kendi tayfı aynı frekanslarda ATOMİK
     ve aralarda BOŞ.
  3) Görevin sorduğu "kanat toplamı" (birinci mertebe pencere kanadı)
     açıklamıyor: 2291–3588× eksik. Artık, üst mertebelerdedir.

ARTEFAKT ŞÜPHELERİM (gizlemiyorum):
  • M-K tutarsız (rms) bir zarf kestirimidir; koherent toplam ±birkaç
    kat oynatabilir ama 3 kademeyi kapatamaz.
  • M-L'nin medyan uyumu (1.06) nokta-nokta uyumdan (0.48–4.7) çok daha
    iyi; "seviye kimliği" evet, "frekans-frekans kimliği" HAYIR.
  • M-L örgüsü Q = e^L'de kesiliyor; kesim noktası σ_lab'ı ve dolayısıyla
    yüksek-ω düzeyini etkiler (200k'da σ %22 fazla).
  • B4 orman testinde M-L'nin çizgi-üstü/çizgi-dışı kontrastı 120k'da
    0.41× (ters işaret) — laboratuvar örgüsü ormanın İNCE yapısını
    güvenilir üretmiyor, yalnız düzeyini.
"""

import numpy as np
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sympy import primerange

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
A_CG = (np.e**2 - 5) / 2
B0, B1 = 2.7580, -0.0543
P11 = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]

def prime_powers(qmax):
    out = []
    for p in primerange(2, int(qmax) + 1):
        q = p
        while q <= qmax:
            out.append((q, float(np.log(p))))     # (q, Λ(q))
            q *= p
    return sorted(out)

def hann(t):
    return 0.5 * (1 - np.cos(TWO_PI * (t - t[0]) / (t[-1] - t[0])))

def I_omega(t, omegas, taper=False, chunk=256):
    omegas = np.asarray(omegas, float)
    w = hann(t) if taper else np.ones_like(t)
    nrm = float((w**2).sum())
    tc = t - 0.5 * (t[0] + t[-1])
    out = np.empty(len(omegas))
    for s0 in range(0, len(omegas), chunk):
        ob = omegas[s0:s0 + chunk]
        S = (np.exp(1j * np.outer(ob, tc)) * w[None, :]).sum(axis=1)
        out[s0:s0 + chunk] = np.abs(S)**2 / nrm
    return out

def n_eff(t, taper):
    w = hann(t) if taper else np.ones_like(t)
    return float(w.sum()**2 / (w**2).sum())

def Ghat(t, omegas, chunk=256):
    omegas = np.asarray(omegas, float)
    tc = t - 0.5 * (t[0] + t[-1])
    out = np.empty(len(omegas), dtype=complex)
    for s0 in range(0, len(omegas), chunk):
        ob = omegas[s0:s0 + chunk]
        out[s0:s0 + chunk] = np.exp(1j * np.outer(ob, tc)).sum(axis=1)
    return out / len(t)

def rvm_N(t):
    x = t / TWO_PI
    return x * np.log(x / np.e) + 7 / 8

def smooth_zeros(z0, m):
    """N(z_k) = N(z_0) + k olan pürüzsüz sıfır ızgarası, k=0..m-1."""
    kk = np.arange(m, dtype=float)
    zs = z0 + kk * TWO_PI / np.log(z0 / TWO_PI)
    for _ in range(8):
        fdel = rvm_N(zs) - rvm_N(z0) - kk
        zs = zs - fdel / (np.log(zs / TWO_PI) / TWO_PI)
    return zs

def u_explicit(t, Q, chunk=120):
    """u(t) = (2/log(t/2π)) Σ_{q≤Q} Λ(q)/(√q log q) · sin(t log q)."""
    PQ = prime_powers(Q)
    acc = np.zeros_like(t)
    for s0 in range(0, len(PQ), chunk):
        blk = PQ[s0:s0 + chunk]
        qs = np.array([b[0] for b in blk], float)
        lam = np.array([b[1] for b in blk], float)
        lq = np.log(qs)
        c = lam / (np.sqrt(qs) * lq)
        acc += (np.sin(np.outer(t, lq)) * c[None, :]).sum(axis=1)
    return 2.0 * acc / np.log(t / TWO_PI), len(PQ)

def unfold(gaps, amps, tmid):
    Lw = np.log(tmid / TWO_PI)
    g_u = gaps * Lw / TWO_PI
    a_u = amps / np.sqrt(A_CG * Lw + B0 + B1 / Lw)
    a_u /= np.sqrt((a_u**2).mean())
    return np.log(a_u), g_u, float(Lw.mean())

# ============================ VERİ ============================
d41 = np.load(HERE / "41_bigT_windows.npz")
K41 = sorted({x.split("_")[1] for x in d41.files}, key=lambda s: int(s[:-1]))
LINES_FULL = [(q, lam) for q, lam in prime_powers(int(np.exp(11.2)))]
LOGQ_FULL = np.array([np.log(q) for q, _ in LINES_FULL])
LAM_FULL = np.array([lam for _, lam in LINES_FULL])
print(f"[102b] tam çizgi listesi: {len(LINES_FULL)} asal-kuvvet, "
      f"log q ≤ {LOGQ_FULL.max():.3f}", flush=True)

OUT = {}
for kkey in [K41[1], K41[0]]:
    gaps, amps, tmid = d41[f"gaps_{kkey}"], d41[f"amps_{kkey}"], d41[f"tmid_{kkey}"]
    ya, g_u, L = unfold(gaps, amps, tmid)
    n = len(ya)
    res_T = TWO_PI / (tmid[-1] - tmid[0])
    band = [p for p in primerange(int(np.exp(0.505*L)), int(np.exp(0.55*L)))][:10]
    kucuk = [p for p in P11 if np.log(p) / L < 0.45]
    print(f"\n{'='*74}\nPENCERE {kkey}: L={L:.4f}  n={n}  çözünürlük 2π/T="
          f"{res_T:.3e}\n  band asalları {band}\n  küçük küme {kucuk}",
          flush=True)

    # ---------- B1: 87'nin T2/T3'ü birebir ----------
    def cols(ps):
        c = []
        for p in ps:
            arg = tmid * np.log(p)
            c += [np.cos(arg), np.sin(arg)]
        return np.vstack(c).T
    Xb = np.hstack([np.ones((n, 1)), cols(band)])
    Xs = cols(kucuk)
    yc = ya - ya.mean()
    bf, *_ = np.linalg.lstsq(np.hstack([Xb, Xs]), yc, rcond=None)
    b_s = bf[Xb.shape[1]:]
    T = np.linalg.solve(Xb.T @ Xb, Xb.T @ Xs)
    pairs = []
    for i, p in enumerate(band):
        for j, q in enumerate(kucuk):
            contrib = (T[1+2*i, 2*j] * b_s[2*j]
                       + T[1+2*i, 2*j+1] * b_s[2*j+1]) / p**-0.5
            pairs.append((abs(contrib), contrib, p, q))
    pairs.sort(reverse=True)
    TOP = pairs[:8]
    freqs, etiket = [], []
    for _, c_, p, q in TOP:
        freqs += [np.log(p) - np.log(q), np.log(p) + np.log(q)]
        etiket += [f"Δ({p},{q})", f"Σ({p},{q})"]
    freqs = np.array(freqs)

    # ---------- ızgaralar ve modeller ----------
    z = np.empty(n + 1)
    z[0] = tmid[0] - gaps[0] / 2
    z[1:] = z[0] + np.cumsum(gaps)
    zs = smooth_zeros(z[0], n + 1)
    ts_sm = 0.5 * (zs[:-1] + zs[1:])                 # M-S: pürüzsüz ortalar
    sig_real = float((tmid - ts_sm).std())

    Qlab = int(np.exp(L))
    ul, nlines = u_explicit(zs, Qlab)
    z_lab = zs + ul
    t_lab = 0.5 * (z_lab[:-1] + z_lab[1:])           # M-L: laboratuvar
    sig_lab = float((t_lab - ts_sm).std())
    print(f"  [B5] SICAKLIK KAPISI: σ_t(gerçek)={sig_real:.4f}  "
          f"σ_t(lab, Q={Qlab}, {nlines} çizgi)={sig_lab:.4f}  "
          f"oran {sig_lab/sig_real:.3f}", flush=True)

    rng = np.random.default_rng(1022)
    g_sh = rng.permutation(gaps)
    t_sh = tmid[0] + np.cumsum(g_sh) - g_sh / 2

    # ---------- M-K: kanat çekirdeği modeli ----------
    def I_wing(om_list, taper):
        """I_kanat(ω) = I_düz(ω) + Σ_q |A_q|²[I_düz(ω−logq)+I_düz(ω+logq)]"""
        om_list = np.asarray(om_list, float)
        base = I_omega(ts_sm, om_list, taper=taper)
        out = base.copy()
        for i, om in enumerate(om_list):
            m = np.abs(LOGQ_FULL - om) < 0.5
            lq = LOGQ_FULL[m]; lam = LAM_FULL[m]
            qv = np.exp(lq)
            tau = lq / L
            A = lam / (L * np.sqrt(qv)) * np.cos(np.pi * tau) \
                * np.exp(-lq**2 * sig_real**2 / 2)
            dd = np.concatenate([om - lq, om + lq])
            Id = I_omega(ts_sm, dd, taper=taper)
            out[i] += float((A**2 * (Id[:len(lq)] + Id[len(lq):])).sum())
        return out, base

    print(f"\n  [B3] ÇİZGİ MESAFESİ (tam liste) + [B2] SEVİYELER", flush=True)
    Ir_h = I_omega(tmid, freqs, taper=False)
    Ir_t = I_omega(tmid, freqs, taper=True)
    Iz_h = I_omega(z, freqs, taper=False)          # SIFIRLARIN kendisi
    Iz_t = I_omega(z, freqs, taper=True)
    Izl_t = I_omega(z_lab, freqs, taper=True)      # laboratuvar SIFIRLARI
    Is_h = I_omega(ts_sm, freqs, taper=False)
    Is_t = I_omega(ts_sm, freqs, taper=True)
    Il_h = I_omega(t_lab, freqs, taper=False)
    Il_t = I_omega(t_lab, freqs, taper=True)
    Ish_h = I_omega(t_sh, freqs, taper=False)
    Ik_h, _ = I_wing(freqs, False)
    Ik_t, _ = I_wing(freqs, True)

    print(f"{'çift':>11} {'ω':>7} {'α=ω/L':>7} {'d_çizgi':>9} {'/(2π/T)':>8} "
          f"{'|Ĝ|ham':>8} {'I ham':>9} {'I Hann':>9} {'M-S ham':>9} "
          f"{'M-K ham':>9} {'M-L ham':>9} {'M-L Hann':>9} {'rampa α':>8}")
    TAB = []
    for i, om in enumerate(freqs):
        dmin = float(np.min(np.abs(LOGQ_FULL - om)))
        Gh = np.sqrt(Ir_h[i] / n)
        print(f"{etiket[i]:>11} {om:>7.3f} {om/L:>7.3f} {dmin:>9.5f} "
              f"{dmin/res_T:>8.1f} {Gh:>8.5f} {Ir_h[i]:>9.4f} {Ir_t[i]:>9.4f} "
              f"{Is_h[i]:>9.4f} {Ik_h[i]:>9.4f} {Il_h[i]:>9.4f} "
              f"{Il_t[i]:>9.4f} {om/L:>8.4f}")
        TAB.append([om, om/L, dmin, dmin/res_T, Ir_h[i], Ir_t[i], Is_h[i],
                    Is_t[i], Ik_h[i], Ik_t[i], Il_h[i], Il_t[i]])
    TAB = np.array(TAB)
    print(f"  ORANLAR (medyan): ölçülen/M-S ham = "
          f"{np.median(Ir_h/Is_h):.2f}×   ölçülen/M-K ham = "
          f"{np.median(Ir_h/Ik_h):.2f}×   ölçülen/M-L ham = "
          f"{np.median(Ir_h/Il_h):.2f}×")
    print(f"                    Hann: /M-S = {np.median(Ir_t/Is_t):.3e}×  "
          f"/M-K = {np.median(Ir_t/Ik_t):.2f}×  /M-L = "
          f"{np.median(Ir_t/Il_t):.2f}×   ölçülen/rampa = "
          f"{np.median(Ir_t/(freqs/L)):.2f}×")
    print(f"  taperli/tapersiz (ölçülen) = {np.median(Ir_t/Ir_h):.3f}  "
          f"(1'e yakınsa içerik GERÇEK, ≪1 ise sızıntıydı)")
    print(f"  karıştırılmış-gap vekili (ham): medyan I = {np.median(Ish_h):.4f}")
    print(f"\n  ★ IZGARA AYRIMI (aynı ω'larda, medyan I):")
    print(f"    ORTA NOKTALAR t_n=z_n+g_n/2 : ham {np.median(Ir_h):.4e}  "
          f"Hann {np.median(Ir_t):.4e}")
    print(f"    SIFIRLAR      z_n           : ham {np.median(Iz_h):.4e}  "
          f"Hann {np.median(Iz_t):.4e}")
    print(f"    lab SIFIRLARI (M-L, z)      : Hann {np.median(Izl_t):.4e}")
    print(f"    orta/sıfır oranı (Hann) = {np.median(Ir_t/Iz_t):.3e}×  "
          f"→ taşıyıcı içeriği hangi nesnede yaşıyor?")

    # ---------- 87'nin "yerel karanlık taban"ı ----------
    dw0 = np.log(TOP[0][2]) - np.log(TOP[0][3])
    scan = dw0 + np.linspace(-0.04, 0.04, 81)
    for isim, tt in [("gerçek", tmid), ("M-S pürüzsüz", ts_sm),
                     ("M-L lab", t_lab)]:
        sh = I_omega(tt, scan, taper=False)
        st = I_omega(tt, scan, taper=True)
        print(f"  87'nin ±0.04 taban taraması @Δω={dw0:.3f} [{isim:>12}]: "
              f"ham medyan |Ĝ|={np.sqrt(np.median(sh)/n):.5f} (I={np.median(sh):.4f}) "
              f"| Hann I={np.median(st):.4f}", flush=True)

    # ---------- B4: ORMAN TESTİ ----------
    om0 = freqs[np.argmax(TAB[:, 4])]
    fine = om0 + np.arange(-0.02, 0.02, res_T / 4)
    If_r = I_omega(tmid, fine, taper=True)
    If_l = I_omega(t_lab, fine, taper=True)
    If_s = I_omega(ts_sm, fine, taper=True)
    If_z = I_omega(z, fine, taper=True)
    near = LOGQ_FULL[(LOGQ_FULL > fine[0]) & (LOGQ_FULL < fine[-1])]
    # zirveleri çizgilere hizala: her çizginin ±2 çözünürlük komşuluğu
    onl = np.zeros(len(fine), bool)
    for l in near:
        onl |= np.abs(fine - l) < 2 * res_T
    print(f"\n  [B4] ORMAN TESTİ @ω={om0:.4f} (±0.02, adım {res_T/4:.2e}, "
          f"{len(fine)} nokta, {len(near)} asal-kuvvet çizgisi):")
    for isim, Iv in [("gerçek", If_r), ("SIFIRLAR", If_z), ("M-L lab", If_l),
                     ("M-S düz", If_s)]:
        print(f"    {isim:>9}: medyan I={np.median(Iv):.4e}  ort={Iv.mean():.4e}  "
              f"maks={Iv.max():.4e}  maks/medyan={Iv.max()/max(np.median(Iv),1e-30):.1f}  "
              f"çizgi-üstü ort={Iv[onl].mean() if onl.any() else np.nan:.4e}  "
              f"çizgi-dışı ort={Iv[~onl].mean():.4e}")
    if onl.any():
        print(f"    ★ çizgi-üstü/çizgi-dışı kontrast: gerçek "
              f"{If_r[onl].mean()/If_r[~onl].mean():.2f}×  "
              f"lab {If_l[onl].mean()/If_l[~onl].mean():.2f}×")
    OUT[kkey] = dict(TAB=TAB, freqs=freqs, etiket=etiket, L=L, n=n,
                     fine=fine, If_r=If_r, If_l=If_l, If_s=If_s, If_z=If_z,
                     om0=om0, sig_real=sig_real, sig_lab=sig_lab,
                     res_T=res_T, Iz_h=Iz_h, Iz_t=Iz_t, Izl_t=Izl_t)

# ---------------- figür ----------------
k0 = K41[1]
D = OUT[k0]
fig, ax = plt.subplots(1, 2, figsize=(14, 5.0))
x = np.arange(len(D["freqs"]))
ax[0].semilogy(x, D["TAB"][:, 4], "o-", c="firebrick", ms=5, label="ölçülen (tapersiz)")
ax[0].semilogy(x, D["TAB"][:, 5], "^-", c="darkorange", ms=5, label="ölçülen (Hann)")
ax[0].semilogy(x, D["TAB"][:, 6], "s-", c="0.5", ms=4, label="M-S sızıntı (pürüzsüz)")
ax[0].semilogy(x, D["TAB"][:, 8], "v-", c="steelblue", ms=4, label="M-K kanat modeli")
ax[0].semilogy(x, D["TAB"][:, 10], "d-", c="seagreen", ms=4, label="M-L laboratuvar örgüsü")
ax[0].semilogy(x, D["TAB"][:, 1], "k:", lw=1.3, label="Montgomery rampası α")
ax[0].set_xticks(x); ax[0].set_xticklabels(D["etiket"], rotation=90, fontsize=6)
ax[0].set_ylabel("I(ω) = n|Ĝ|²")
ax[0].set_title("87-T3 taşıyıcı frekansları: ölçüm vs parametresiz modeller")
ax[0].legend(fontsize=8); ax[0].grid(alpha=0.25, which="both")
ax[1].semilogy(D["fine"] - D["om0"], D["If_r"], "-", lw=0.8, c="firebrick",
               label="gerçek ızgara (Hann)")
ax[1].semilogy(D["fine"] - D["om0"], D["If_l"], "-", lw=0.8, c="seagreen",
               label="laboratuvar örgüsü (Hann)")
ax[1].semilogy(D["fine"] - D["om0"], D["If_s"], "-", lw=0.8, c="0.6",
               label="pürüzsüz (sızıntı)")
ax[1].set_xlabel(f"ω − {D['om0']:.4f}"); ax[1].set_ylabel("I(ω)")
ax[1].set_title("ORMAN TESTİ: sürekli omuz mu, çözülmüş çizgiler mi?")
ax[1].legend(fontsize=8); ax[1].grid(alpha=0.25, which="both")
fig.tight_layout()
fig.savefig(HERE / "102b_kanat.png", dpi=125)
np.savez(HERE / "102b_kanat.npz",
         **{f"{k}_{kk}": v for k, D_ in OUT.items()
            for kk, v in D_.items() if not isinstance(v, list)})
print("\n102b_kanat.png + 102b_kanat.npz yazıldı.", flush=True)
