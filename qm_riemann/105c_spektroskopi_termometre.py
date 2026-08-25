"""
105c — YENİ ADALARDA SPEKTROSKOPİ + TARAK-TERMOMETRESİ (25 Ağustos 2026)
==========================================================================
99'un kampanya kalıbı üç yeni adaya (chi5e, chi8e, chi8o) uygulanır:
benek genlikleri + fazlar + MUTLAK BENEK YASASI (89) + ada-içi 3
pencerede tarak-termometresi (koro-yasası, 97).

MUTLAK YASA (89, parametresiz):
    |Ĝ(log q)| = Λ(q)/(L·√q) · cos(πτ) · DW,   τ = log q / L,
    DW = exp(−ω²σ_t²/2),  σ_t = σ_u·2π/L  (tarak-termometresinden)
FAZ YASASI (96/97/99):  faz = 180° + arg χ(q)

ÖN-MÜHÜRLÜ ÖNGÖRÜLER (105 kampanya mühründen, ölçümden ÖNCE):
  P1  GAMMA EVRENSELLİĞİ: çift adalar (a=0) her şeyde tek adalar gibi
      davranır — mutlak yasa oranları ~1 (%1-5 bandında), log-ısınma,
      hiperuniform sertifika. Gamma faktörü HİÇBİR yasaya girmemeli.
  P2  SPEKTROSKOPİ (fazlar = 180° + arg χ, reel karakter → 0° veya 180°):
      chi5e: 2→0°, 3→0°, 4→180°, 7→0°, 9→180°;  ölü 5, 25
      chi8e: 3→0°, 5→0°, 7→180°, 9→180°;        ölü 2, 4, 8
      chi8o: 3→180°, 5→0°, 7→0°;                ölü 2, 4, 8
  P4  (kısmî — tam hüküm 105d ile) KORO/β-SOĞUKLUĞU: koro-düzeltmeli
      tam sıcaklıklar σ_tam² = σ_u² + Σ_{ölü q} U_q²/2,
      U_q = 2Λ(q)/(L√q log q)  (97 yöntemi). chi8e/chi8o, β ile AYNI
      artık-soğukluğu gösterirse "2'li ailesizlik" açıklaması güçlenir.

VERİ: 105b_{ada}_zeros.npz (sertifikalı). Kıyas için ESKİ dört ada
101f_{ada}_zeros.npz — aynı sertifika ayağında (99'un tablosu
99_*_zeros.npz üstündeydi; fark ayrıca raporlanır).

==========================================================================
SONUÇ (25 Ağustos) — P2 ✓✓✓  P1 ✓✓  P4 ✓✓✓
==========================================================================
P2 SPEKTROSKOPİ ✓✓✓ — ÜÇÜ DE DERECENİN ONDA BİRİ İÇİNDE:
  chi5e maks sapma 0.11° (2→0.00, 3→0.02, 4→179.97, 7→359.99, 9→179.89)
  chi8e maks sapma 0.07° (3→0.05, 5→359.93, 7→180.02, 9→180.00)
  chi8o maks sapma 0.04° (3→179.99, 5→359.99, 7→0.04)
  Mezarlıklar: chi5e {5,25} = 0.0000/0.0001; chi8e/chi8o {2,4,8} =
  0.0000-0.0001 — canlı medyanın binde 1-2'si (plasebo dibi).
  chi8e ile chi8o AYNI |Ĝ| (0.0593/0.0618/0.0583/...) ama fazları tam
  arg χ kadar ayrı: aynı örgü, farklı imza.

P1 MUTLAK YASA ✓✓ — asal oranları medyanı: chi5e 1.009 (a=0!),
  chi8e 0.986 (a=0), chi8o 0.985 (a=1), chi3 0.997, β 0.976,
  chi5 1.009, chi7 1.013. chi5e (a=0) ile chi5 (a=1) AYNI q=5,
  aynı L: 1.009 vs 1.009 — gamma faktörü yasaya GİRMİYOR.
  Kuvvet açığı da paritesiz: chi8e/chi8o (1−ε)/k sütunları 0.03/0.12/
  0.12/0.31 ile BİREBİR aynı. Log-ısınma üç yeni adada da net.

P4 KORO/β-SOĞUKLUĞU ✓✓✓ — ÖLÜ 2-AİLESİ AÇIKLAMASI DOĞRULANDI:
  σ_tam²(L) = A + B·L fitiyle L=9'da:
    canlı-2 adaları: chi3 0.0582, chi5 0.0630, chi7 0.0641,
                     chi5e 0.0629  → ortalama 0.0620
    ölü-2 adaları:   β 0.0496 (−20.0%), chi8e 0.0500 (−19.4%),
                     chi8o 0.0500 (−19.4%)
  Üç ada, iki farklı iletken (4 ve 8), iki farklı parite — AYNI artık
  soğukluk. β'nın soğukluğu iletken-4'e veya imprimitifliğe özgü
  DEĞİL; 2-ailesinin ölü olmasının sonucu. Koro yasası ("her ada
  kendi korosu kadar ısınır") artık koroyu TAM açıklamıyor: susturulan
  2-ailesinin geri eklenen varyansı yetmiyor, ~%20'lik bir artık var
  ve o artık DA yalnız çizgi envanterine bağlı.
  Not: chi5e (a=0) ve chi5 (a=1) 0.0629 vs 0.0630 — binde bir.

P1 (3. ayak) HİPERUNİFORMLUK ✓ — sayı varyansı Σ²(n) yedi adada da
  n'de LOGARİTMİK (Poisson'un 8-124 katı bastırılmış), GUE Dyson-Mehta
  formunun 0.65-1.55 katı:
    n=       5      10      20      50   (GUE: 0.384/0.454/0.525/0.617)
    chi5e 0.594   0.477   0.541   0.576   (a=0)
    chi8e 0.520   0.459   0.454   0.500   (a=0)
    chi8o 0.529   0.459   0.449   0.498   (a=1)
    chi5  0.539   0.443   0.500   0.518   (a=1, chi5e'nin ikizi)
    chi3  0.504   0.466   0.524   0.472
    beta  0.456   0.402   0.412   0.431
    chi7  0.377   0.346   0.354   0.402
  chi8e ile chi8o farkı ≤0.01 (parite görünmüyor); chi5e ile chi5
  farkı ≤0.06. Kaba ama bağımsız bir hiperuniformluk sertifikası
  (102a'nın taperli/vekilli titiz ölçümü DEĞİL — o düzeyde iddia
  edilmiyor; burada yalnızca a=0'ın a=1'den ayrılmadığı gösteriliyor).
"""

import numpy as np
from pathlib import Path

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
exec(open(HERE / "98_L_motoru.py").read().split('CHI4 = ')[0])

# ---- karakter tabloları ----
CHI3 = {0: 0, 1: 1, 2: -1}
CHI4 = {0: 0, 1: 1, 2: 0, 3: -1}
CHI5 = {0: 0, 1: 1, 2: 1j, 3: -1j, 4: -1}
z6 = np.exp(1j * np.pi / 3)
CHI7 = {0: 0, 1: 1, 3: z6, 2: z6**2, 6: z6**3, 4: z6**4, 5: z6**5}
CHI5E = {0: 0, 1: 1, 2: -1, 3: -1, 4: 1}
CHI8E = {0: 0, 1: 1, 2: 0, 3: -1, 4: 0, 5: -1, 6: 0, 7: 1}
CHI8O = {0: 0, 1: 1, 2: 0, 3: 1, 4: 0, 5: -1, 6: 0, 7: -1}

PK = [(2, 2, 1), (3, 3, 1), (4, 2, 2), (5, 5, 1), (7, 7, 1), (8, 2, 3),
      (9, 3, 2), (11, 11, 1), (13, 13, 1), (16, 2, 4), (17, 17, 1),
      (19, 19, 1), (23, 23, 1), (25, 5, 2), (27, 3, 3), (29, 29, 1),
      (31, 31, 1), (49, 7, 2)]

ADALAR = [("chi5e", 5, CHI5E, 0, "105b_chi5e_zeros.npz", "YENİ"),
          ("chi8e", 8, CHI8E, 0, "105b_chi8e_zeros.npz", "YENİ"),
          ("chi8o", 8, CHI8O, 1, "105b_chi8o_zeros.npz", "YENİ"),
          ("chi3",  3, CHI3,  1, "101f_chi3_zeros.npz",  "eski"),
          ("beta",  4, CHI4,  1, "101f_beta_zeros.npz",  "eski"),
          ("chi5",  5, CHI5,  1, "101f_chi5_zeros.npz",  "eski"),
          ("chi7",  7, CHI7,  1, "101f_chi7_zeros.npz",  "eski")]


def Ghat(t, omegas, chunk=30000):
    out = np.zeros(len(omegas), dtype=complex)
    for s0 in range(0, len(t), chunk):
        tt = t[s0:s0 + chunk]
        out += np.exp(1j * np.outer(omegas, tt)).sum(axis=1)
    return out / len(t)


def olu_aile(q, tab):
    """Susturulmuş Euler ailesi: χ(p)=0 olan p'lerin tüm kuvvetleri."""
    olu = []
    for p in [2, 3, 5, 7, 11, 13]:
        if tab[p % q] == 0:
            k = 1
            while p**k <= 512:
                olu.append((p**k, np.log(p)))
                k += 1
    return olu


SUM = {}
print("=" * 78, flush=True)
print("105c — SPEKTROSKOPİ + TARAK-TERMOMETRESİ", flush=True)
print("=" * 78, flush=True)

# --- KOD KAPISI: koro formülü 97'nin yayımlanmış sayılarını üretiyor mu ---
print("KOD KAPISI — koro formülü vs 97'nin yayımlanmış tam sıcaklıkları:",
      flush=True)


def _kor(L, p, kmax=9):
    return sum((2 * np.log(p) / (L * np.sqrt(p**k) * np.log(p**k)))**2 / 2
               for k in range(1, kmax + 1))


_kapi = True
for _et, _su, _L, _p, _ref in [("χ₃", 0.196, 6.92, 3, 0.0535),
                               ("β", 0.161, 7.20, 2, 0.0481),
                               ("χ₅", 0.220, 7.30, 5, 0.0562)]:
    _v = _su**2 + _kor(_L, _p)
    _d = 100 * (_v / _ref - 1)
    if abs(_d) > 2:
        _kapi = False
    print(f"  {_et:>3}: hesap {_v:.4f} vs 97 {_ref:.4f} ({_d:+.2f}%)",
          flush=True)
print(f"  KAPI: {'✓ GEÇTİ' if _kapi else '✗ RET'} "
      f"(97'nin L değerleri yuvarlanmış; %2 tolerans)\n", flush=True)

for etiket, q, tab, a, dosya, tur in ADALAR:
    p = HERE / dosya
    if not p.exists():
        print(f"[{etiket}] {dosya} YOK — atlanıyor", flush=True); continue
    zz = np.load(p)["zeros"]
    gaps = np.diff(zz)
    mids = 0.5 * (zz[:-1] + zz[1:])
    n = len(mids)
    Leff = float(np.log(q * mids / TWO_PI).mean())

    def Nsm(t, _q=q):
        x = t / TWO_PI
        return x * np.log(_q * x / np.e)

    print(f"\n===== {etiket} ({tur}; q={q}, a={a}): n={n}, "
          f"L_eff={Leff:.3f}, t∈[{zz[0]:.0f},{zz[-1]:.0f}] =====", flush=True)

    # --- ada-içi 3 pencere: tarak-termometresi ---
    pens = []
    edges = np.exp(np.linspace(np.log(zz[0] + 1), np.log(zz[-1]), 4))
    for i in range(3):
        m = (mids >= edges[i]) & (mids < edges[i + 1])
        mm = mids[m]
        Lw = float(np.log(q * mm / TWO_PI).mean())
        xw = Nsm(mm) - Nsm(mm[0])
        comb = abs(np.exp(2j * np.pi * xw).mean())
        su = np.sqrt(-2 * np.log(comb)) / TWO_PI
        pens.append((Lw, su, int(m.sum())))
        print(f"  pencere {i+1}: L={Lw:.2f}  σ_u = {su:.4f}  (n={m.sum()})",
              flush=True)
    x_all = Nsm(mids) - Nsm(zz[0])
    comb_all = abs(np.exp(2j * np.pi * x_all).mean())
    sig_u = np.sqrt(-2 * np.log(comb_all)) / TWO_PI
    sig_t = sig_u * TWO_PI / np.log(q * mids.mean() / TWO_PI)

    # --- benekler: genlik + faz + mutlak yasa ---
    print(f"  {'q':>3} {'sınıf':>6} {'ölçüm':>8} {'öngörü':>8} {'oran':>6} "
          f"{'faz':>8} {'ö-faz':>6} {'Δfaz':>6}", flush=True)
    rows = []
    for qq, pp, kk in PK:
        om = np.log(qq)
        tau = om / Leff
        G = Ghat(mids, np.array([om]))[0]
        chi = tab[qq % q]
        ph = (np.degrees(np.angle(G)) + 360) % 360
        lam = np.log(pp)
        if chi == 0:
            print(f"  {qq:>3} {'ÖLÜ':>6} {abs(G):>8.4f} {'—':>8} {'—':>6} "
                  f"{ph:>7.1f}° {'—':>6} {'—':>6}", flush=True)
            rows.append((qq, kk, tau, abs(G), np.nan, ph, np.nan, np.nan))
        else:
            dw = np.exp(-om**2 * sig_t**2 / 2)
            pred = lam / (Leff * np.sqrt(qq)) * np.cos(np.pi * tau) * dw
            pph = (180 + np.degrees(np.angle(complex(chi)))) % 360
            dph = (ph - pph + 180) % 360 - 180
            print(f"  {qq:>3} {qq % q:>6} {abs(G):>8.4f} {pred:>8.4f} "
                  f"{abs(G)/pred:>6.3f} {ph:>7.1f}° {pph:>5.0f}° "
                  f"{dph:>+6.1f}", flush=True)
            rows.append((qq, kk, tau, abs(G), abs(G) / pred, ph, pph, dph))
    SUM[etiket] = dict(L=Leff, sig_u=sig_u, sig_t=sig_t, pens=pens,
                       rows=rows, n=n, q=q, a=a, tab=tab, tur=tur)

# ================== P1 / P2 HÜKÜMLERİ ==================
print("\n" + "=" * 78, flush=True)
print("P2 — SPEKTROSKOPİ HÜKMÜ (faz = 180° + arg χ; reel karakter → 0/180°)",
      flush=True)
ONGORU = {"chi5e": {2: 0, 3: 0, 4: 180, 7: 0, 9: 180},
          "chi8e": {3: 0, 5: 0, 7: 180, 9: 180},
          "chi8o": {3: 180, 5: 0, 7: 0}}
OLU_ONG = {"chi5e": [5, 25], "chi8e": [2, 4, 8], "chi8o": [2, 4, 8]}
for et, ong in ONGORU.items():
    if et not in SUM:
        continue
    d = SUM[et]
    rr = {r[0]: r for r in d["rows"]}
    sap = []
    for qq, pfaz in ong.items():
        r = rr[qq]
        dd = (r[5] - pfaz + 180) % 360 - 180
        sap.append(abs(dd))
        print(f"  {et} q={qq:>2}: ölçüm {r[5]:7.2f}°  öngörü {pfaz:>3}°  "
              f"Δ = {dd:+6.2f}°  |Ĝ| = {r[3]:.4f}", flush=True)
    olu_g = [rr[qq][3] for qq in OLU_ONG[et]]
    canli_g = [r[3] for r in d["rows"] if not np.isnan(r[4])]
    print(f"  {et}: ölü benekler {OLU_ONG[et]} → "
          + " ".join(f"{g:.4f}" for g in olu_g)
          + f"  | canlı medyan {np.median(canli_g):.4f} "
          f"(oran {np.median(olu_g)/np.median(canli_g):.3f})", flush=True)
    print(f"  → {et} MAKS FAZ SAPMASI = {max(sap):.2f}°", flush=True)

print("\n" + "=" * 78, flush=True)
print("P1 — MUTLAK YASA ORANLARI (asallar; 1.00 = parametresiz yasa)",
      flush=True)
print(f"{'ada':>6} {'a':>2} {'L':>6}  asal oranları (q: oran)", flush=True)
for et in SUM:
    d = SUM[et]
    asal = [(r[0], r[4]) for r in d["rows"]
            if r[1] == 1 and not np.isnan(r[4])]
    med = np.median([x[1] for x in asal])
    print(f"{et:>6} {d['a']:>2} {d['L']:>6.2f}  "
          + " ".join(f"{q}:{o:.3f}" for q, o in asal[:9])
          + f"   | medyan {med:.3f}", flush=True)
print("\nKUVVET AÇIĞI (k≥2 örgüler): (1−ε)/k vs τ", flush=True)
print(f"{'ada':>6}  " + "  ".join(f"{q:>10}" for q in [4, 8, 9, 25, 27, 49]),
      flush=True)
for et in SUM:
    d = SUM[et]
    rr = {r[0]: r for r in d["rows"]}
    cells = []
    for qq in [4, 8, 9, 25, 27, 49]:
        r = rr.get(qq)
        if r is None or np.isnan(r[4]):
            cells.append(f"{'ölü':>10}")
        else:
            cells.append(f"{r[4]:.3f}/{(1-r[4])/r[1]:.3f}"[:10].rjust(10))
    print(f"{et:>6}  " + "  ".join(cells), flush=True)

# ================== P4 — KORO TERMOMETRESİ ==================
print("\n" + "=" * 78, flush=True)
print("P4 — KORO-DÜZELTMELİ TAM SICAKLIK  σ_tam² = σ_u² + Σ_ölü U_q²/2",
      flush=True)
print(f"{'ada':>6} {'a':>2} {'L':>6} {'σ_u':>7} {'σ_u²':>8} {'koro':>8} "
      f"{'σ_tam²':>8} {'ölü aile':>12}", flush=True)
TAM = {}
for et in SUM:
    d = SUM[et]
    olu = olu_aile(d["q"], d["tab"])
    L = d["L"]
    kor = sum((2 * lam / (L * np.sqrt(qq) * np.log(qq)))**2 / 2
              for qq, lam in olu)
    s2 = d["sig_u"]**2
    TAM[et] = (L, d["sig_u"], s2, kor, s2 + kor)
    aile = ",".join(str(qq) for qq, _ in olu[:4]) if olu else "—"
    print(f"{et:>6} {d['a']:>2} {L:>6.2f} {d['sig_u']:>7.4f} {s2:>8.4f} "
          f"{kor:>8.4f} {s2+kor:>8.4f} {aile:>12}", flush=True)

print("\nPENCERE-BAZLI (ada-içi ısınma + koro düzeltmesi):", flush=True)
print(f"{'ada':>6} " + " ".join(f"{'L':>6} {'σ_u':>6} {'σ_tam²':>7}"
                                for _ in range(3)), flush=True)
PW = {}
for et in SUM:
    d = SUM[et]
    olu = olu_aile(d["q"], d["tab"])
    cells, ws = [], []
    for (Lw, su, nn) in d["pens"]:
        kor = sum((2 * lam / (Lw * np.sqrt(qq) * np.log(qq)))**2 / 2
                  for qq, lam in olu)
        cells.append(f"{Lw:>6.2f} {su:>6.4f} {su**2+kor:>7.4f}")
        ws.append((Lw, su, su**2 + kor))
    PW[et] = ws
    print(f"{et:>6} " + " ".join(cells), flush=True)

# --- β vs mod-8 karşılaştırması: aynı ölü aile (2'nin kuvvetleri) ---
print("\nP4 TEŞHİS TABLOSU — '2'li ailesizlik' mi, β-özgü mü?", flush=True)
print("  (σ_tam² ~ L'de lineer ısınıyor; her adayı ortak L ızgarasında",
      flush=True)
print("   lineer fitle kıyasla: σ_tam²(L) = A + B·L)", flush=True)
print(f"{'ada':>6} {'a':>2} {'ölü aile':>10} {'A':>8} {'B':>8} "
      f"{'σ_tam²@L=9':>11} {'σ_tam²@L=10':>12}", flush=True)
FIT = {}
for et in PW:
    ws = np.array(PW[et])
    if len(ws) < 2:
        continue
    B, A = np.polyfit(ws[:, 0], ws[:, 2], 1)
    d = SUM[et]
    olu = olu_aile(d["q"], d["tab"])
    aile = str(olu[0][0]) + "-ailesi" if olu else "—"
    FIT[et] = (A, B)
    print(f"{et:>6} {d['a']:>2} {aile:>10} {A:>8.4f} {B:>8.5f} "
          f"{A+9*B:>11.4f} {A+10*B:>12.4f}", flush=True)

ref = [et for et in FIT if et in ("chi3", "chi5", "chi7", "chi5e")]
if ref and "beta" in FIT:
    r9 = np.mean([FIT[e][0] + 9 * FIT[e][1] for e in ref])
    print(f"\n  'canlı-2' adalarının ortalaması @L=9: {r9:.4f}", flush=True)
    for et in ("beta", "chi8e", "chi8o"):
        if et in FIT:
            v = FIT[et][0] + 9 * FIT[et][1]
            print(f"  {et:>6}: {v:.4f}  → artık soğukluk "
                  f"{100*(v/r9-1):+.1f}%", flush=True)

# ============ P1 (3. ayak) — HİPERUNİFORMLUK: SAYI VARYANSI ============
print("\n" + "=" * 78, flush=True)
print("P1 (3. ayak) — HİPERUNİFORMLUK: Σ²(n) = Var[N(kutu)] ", flush=True)
print("  kutu = n ortalama boşluk; Poisson: Σ² = n;", flush=True)
print("  GUE (Dyson-Mehta): Σ² ≈ (1/π²)(ln 2πn + γ + 1) − 1/8", flush=True)
GAM = 0.5772156649


def sigma2(zz, q, nlist=(5, 10, 20, 50)):
    x = (zz / TWO_PI) * np.log(q * zz / TWO_PI / np.e)   # açılmış sayım
    x = x - x[0]
    out = []
    for n in nlist:
        kutu = np.arange(x[0], x[-1] - n, n)
        cnt = np.diff(np.searchsorted(x, np.concatenate([kutu, [kutu[-1] + n]])))
        out.append(float(np.var(cnt[:-1])))
    return out


NL = (5, 10, 20, 50)
print(f"\n{'ada':>6} {'a':>2} " + " ".join(f"{'n='+str(n):>18}" for n in NL),
      flush=True)
print(f"{'GUE':>6} {'—':>2} " + " ".join(
    f"{(np.log(2*np.pi*n)+GAM+1)/np.pi**2 - 0.125:>18.3f}" for n in NL),
    flush=True)
for et in SUM:
    d = SUM[et]
    zz = np.load(HERE / dict((x[0], x[4]) for x in ADALAR)[et])["zeros"]
    s2v = sigma2(zz, d["q"], NL)
    gue = [(np.log(2 * np.pi * n) + GAM + 1) / np.pi**2 - 0.125 for n in NL]
    print(f"{et:>6} {d['a']:>2} " + " ".join(
        f"{v:6.3f} ({v/g:4.2f}G {n/v:4.0f}P)"
        for v, g, n in zip(s2v, gue, NL)), flush=True)
print("  (parantez: GUE'ye oran, Poisson'a göre kaç kat bastırılmış)",
      flush=True)

np.savez(HERE / "105c_spektroskopi.npz",
         **{f"rows_{et}": np.array([[r[0], r[1], r[2], r[3], r[4], r[5],
                                     r[6] if not np.isnan(r[6]) else np.nan,
                                     r[7] if not np.isnan(r[7]) else np.nan]
                                    for r in SUM[et]["rows"]])
            for et in SUM},
         **{f"pens_{et}": np.array(SUM[et]["pens"]) for et in SUM},
         **{f"meta_{et}": np.array([SUM[et]["L"], SUM[et]["sig_u"],
                                    SUM[et]["sig_t"], SUM[et]["n"],
                                    SUM[et]["q"], SUM[et]["a"]])
            for et in SUM})
print("\nBİTTİ — 105c_spektroskopi.npz yazıldı", flush=True)
