"""
119 — GENLİK-KANALI NEFESİ R_w: KORUNUM-ÇERÇEVESİ SINAVI (26 Ağustos)
==========================================================================
Kalem defteri Ö2: |Z| genlik-gürültüsünün varyansı, aritmetik gap-
dalgasının fazında nasıl salınır? İlk kez ölçülüyor. İki hipotez
(ön-mühür; ikisi de bilgilendirici):
  H-korunum: kanallar-arası bütçe → gap-gürültüsü kısılan (gerilen)
     bölgede genlik-gürültüsü BÜYÜR: R_w işareti R_p'nin TERSİ (+).
  H-uniform: gaz gerilen bölgede HER kanalda sessizleşir: R_w ≈ R_p (−).
Makine: 41 (gaps+amps+tmid, iki pencere); ds tam-taban fiti → çizgi
fazları; a_u = amps/√(A·L+B0+B1/L) (rms-normlu), ξ = a_u'nun tam-taban
artığı; ξ²-dalgası gap-dalgası yönüne izdüşürülür:
  R_w(p) = P[ξ²]_p / (2σ_ξ²·A1_gap(p)).
Plasebo: +0.037 kaydırılmış sahte çizgiler.

SONUÇ (26 Ağustos) — H-KORUNUM EZİCİ KAZANDI: ÇAPRAZ-KANAL NEFES TAKASI:
  R_w = +6.26/+3.74/+2.32/+1.76/+1.26/+1.12/+0.94 (p=2..17) —
  GÜÇLÜ POZİTİF; R_p aynı koşuda −0.69..−0.92; plasebo 0.05-0.18.
  Aritmetik dalganın GERDİĞİ bölgede gap-gürültüsü kısılırken
  GENLİK-GÜRÜLTÜSÜ KÜKRÜYOR — iki kanal, gazın kendi dalgasının
  fazına kilitli biçimde dalgalanma TAKASI yapıyor. 68'in korunum
  yasası (√w+(β/2)v=1) ve Not 1'in gap-genlik korelasyonu, GÜRÜLTÜ
  DÜZEYİNDE de doğrulandı: kayıpsız sistemin ikinci-moment sureti.
  Yapı notu: R_w ~ hiperbolik düşüş (R_w·τ ≈ 0.44→0.26 yavaş iniş) —
  biçim yasası açık. Denge-dışı ilkenin adı netleşiyor:
  "SABİT DALGALANMA BÜTÇESİ: koheran aritmetik + iki kanalın gürültüsü,
  faz-çözümlü bir korunumla bağlı." Türetim (hangi büyüklük korunuyor —
  |Z|²'nin mi, enerji-akısının mı?) sıradaki kalem ödevi; Ö1
  (adalarda R) hâlâ mühürlü ayrıştırıcı.
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
SON = {}
for k in ["120k", "200k"]:
    gz, am, tm = d41[f"gaps_{k}"], d41[f"amps_{k}"], d41[f"tmid_{k}"]
    Lw = np.log(tm / TWO_PI)
    L = float(Lw.mean())
    ds = gz * Lw / TWO_PI - 1
    a_u = am / np.sqrt(A_CG * Lw + B0 + B1 / Lw)
    a_u = a_u / np.sqrt((a_u**2).mean())
    qs = pk_list(min(np.exp(0.52 * L), 720))
    freqs = [np.log(q) for q in qs] + [np.log(p) + 0.037 for p in PR]
    bg, fg = chunked_fit(ds, tm, freqs)
    ba, fa = chunked_fit(a_u, tm, freqs)
    xi = a_u - fa
    sx2 = float((xi**2).mean())
    b2, _ = chunked_fit(xi**2 - sx2, tm, freqs)
    # gap-artığı nefesi (çapraz, aynı koşuda): η²-dalgası
    eta = ds - fg
    se2 = float((eta**2).mean())
    b3, _ = chunked_fit(eta**2 - se2, tm, freqs)
    print(f"[{k}] L={L:.2f}  σ_ξ={np.sqrt(sx2):.4f}  σ_η={np.sqrt(se2):.4f}",
          flush=True)
    for p in PR:
        f = np.log(p)
        cg, sg = coef(bg, freqs, f)
        A1 = np.hypot(cg, sg); ph = np.arctan2(sg, cg)
        cw, sw = coef(b2, freqs, f)
        Rw = (cw * np.cos(ph) + sw * np.sin(ph)) / (2 * sx2 * A1)
        ce, se_ = coef(b3, freqs, f)
        Rp = (ce * np.cos(ph) + se_ * np.sin(ph)) / (2 * se2 * A1)
        fpl = np.log(p) + 0.037
        cpl, spl = coef(b2, freqs, fpl)
        plw = np.hypot(cpl, spl) / (2 * sx2 * A1)
        SON.setdefault(p, []).append((Rw, Rp, plw))

print(f"\n{'p':>3} {'R_w':>7} {'R_p':>7} {'plasebo_w':>10}   "
      f"(H-korunum: R_w>0; H-uniform: R_w≈R_p<0)")
for p in PR:
    arr = np.array(SON[p])
    a = arr.mean(axis=0)
    print(f"{p:>3} {a[0]:>+7.3f} {a[1]:>+7.3f} {a[2]:>10.3f}")
