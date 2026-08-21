"""
87b — 87'NİN KAPALI DEVRESİ: TAPERLİ YENİDEN DENETİM (20 Ağustos)
==========================================================================
GÖREV (kaptan): 101j, 87'nin T3 hükmünün ("taşıyıcılar yerel karanlık
tabanın 4–15 katı → TDS") TAPERSİZ Ĝ ile kurulduğunu ve o "taban"ın
dikdörtgen-pencere sızıntısı düzeyinde olduğunu gösterdi. 87'nin üç
testi tabana ne kadar yaslanıyor, taperle ne kalıyor?

DURUM NOTU (dürüstlük): Bu denetim yazılırken 102b (kanat öngörüsü)
zaten koşulmuş ve mühürlenmişti; ön-mühür tahminleri 102b'nin
bulgularını BİLEREK yazıldı. 102b'nin kapsamadığı boşluklar burada:
  - 87-T1 hiç yeniden koşulmadı (102b yalnız T2/T3'ü koştu) ve
    kimliğin TAPERLİ biçimi hiç test edilmedi;
  - 87-T2'nin kapanış artığı (b_çıplak − b_tam − T·b_s) hiç ölçülmedi;
  - T3 için frekans-frekans "taşıyıcı / taperli yerel taban" tablosu
    yok (102b medyan düzeyinde ve M-S/M-K/M-L modellerine karşı baktı;
    çizgi-DIŞLANMIŞ yerel taban hiç ölçülmedi).

ÖLÇÜMLER (kurulum 87 ile birebir: 41_bigT_windows, band [0.505,0.55]L,
P11 küçük küme, iki pencere):
  D1  T1 KİMLİK: 87'nin Gram=Ĝ doğrulaması birebir + AYNI kimliğin
      Hann-AĞIRLIKLI biçimi: (1/Σw)Σ w·cos(ω_a t)cos(ω_b t) =
      ½Re[Ĝ_w(ω_a−ω_b)+Ĝ_w(ω_a+ω_b)], Ĝ_w = Σ w e^{iωt}/Σw.
      Kimlik tabana/kestirimciye yaslanıyor mu?
  D2  T2 TRANSFER: 87 birebir + kapanış artığı
      max|b_çıplak − (b_tam + T·b_s)| (kesin cebir makinede kapanmalı).
      Kod denetimi: T2'ye herhangi bir spektral kestirimci giriyor mu?
  D3  T3 TAŞIYICI/TABAN: 87'nin en büyük 8 çiftinin 16 (Δω, Σω)
      frekansında |Ĝ| tapersiz VE Hann. Yerel taban: ω±0.04, adım
      2·(2π/T), TAM p^k listesiyle (log q ≤ 11.2) çizgilerin ±4
      çözünürlüğü ve taşıyıcının ±6 çözünürlüğü DIŞLANMIŞ medyan.
      Sızıntı kontrolü: aynı tarama pürüzsüz orta ızgarada (M-S).
      Kimlik kontrolü: aynı 16 frekansta sıfır ızgarası z_n (Hann).

ÖN-MÜHÜR (koşudan ÖNCE yazıldı; 102b bilgisi dahil):
  P1  T1 kimliği tapersiz VE Hann-ağırlıklı biçimde hata < 1e-6
      (Gram ölçeği ~0.006) ile geçecek. Kimlikte taban diye bir terim
      YOK; hüküm tabandan bağımsız çıkacak.
  P2  T2 kapanış artığı lstsq hassasiyetinde (≲1e-8) kapanacak; T2'ye
      Ĝ de taban da girmiyor → 87-T2 tabandan bağımsız.
  P3  Taşıyıcılar Hann ile ÇÖKMEYECEK (102b: I oranı 0.62–0.67).
      Taperli çizgi-dışlamalı taban 87'nin 5e-4'ünün altında ama M-S
      sızıntısının ÇOK üstünde (gerçek içerik) kalacak; taşıyıcı/taban
      oranı 87'nin 4–15× bandıyla aynı mertebede çıkacak.
  P4  Aynı frekanslarda sıfır ızgarası ≥10⁷ kat karanlık (102b) →
      "sıfır kristalinin TDS'si" kimliği RET, taşıyıcı içeriği
      ORTA-NOKTA örneklemesinin yan bandı olarak yeniden adlandırılır;
      87'nin giydirme mekanizması ve düzeyleri AYAKTA kalır.

==========================================================================
SONUÇ (20 Ağustos, koşu 8 s) — P1 ✓, P2 ✓, P3 ✓ (zayıflatılmış), P4 ✓
T1/T2 TABANDAN BAĞIMSIZ VE AYAKTA. T3 TAŞIYICILARI TAPERİ GEÇİYOR VE
TABANIN ÜSTÜNDE DURUYOR — AMA FREKANS-EŞLENİK TABANLA ORAN 3–7×,
87'NİN "4–15×"İ DEĞİL; VE KİMLİK SIFIR KRİSTALİNİN DEĞİL.
==========================================================================
D1 T1 KİMLİK ✓ (P1): tapersiz hata 4.7e-08 / 3.3e-08 (87'nin ~5e-8'i ✓),
  Hann-AĞIRLIKLI aynı kimlik 5.0e-08 / 3.0e-08 (Gram ölçeği ~0.006).
  Kimlik her iki kestirimcide kesin; tabanla hiçbir teması yok.
D2 T2 TRANSFER ✓ (P2): kapanış artığı max|b_ç−(b_t+T·b_s)| = 6.7e-15 /
  9.9e-15; B̂ farkı +0.429 / +0.450 (87 birebir ✓). Kod denetimi: T2'de
  Ĝ de taban da geçmiyor — yalnız X'X cebiri. Tabandan bağımsız.
D3 T3 (16 frekans × 2 pencere; |Ĝ| birimi, Ĝ_w = Σw e^{iωt}/Σw):
  • İÇERİK GERÇEK: Hann/tapersiz |Ĝ|² medyanı 0.666 / 0.621; taşıyıcı
    |Ĝ| taperle pratik olarak değişmiyor (ör. 2.18e-3 → 2.18e-3).
    102b çaprazı birebir (I medyanları 0.367→0.241 / 0.604→0.340 ✓).
  • TABAN GERÇEK: taperli, çizgi-dışlamalı yerel taban 4.2e-4–1.6e-3
    (medyan 8.3e-4 / 7.8e-4); M-S pürüzsüz-örgü sızıntısı aynı taramada
    ~2e-13–5e-13 → taban sızıntının 2.7e9 / 4.5e9 KATI. 101j'nin "taban
    alet-sınırlı" şüphesi ω≈4–7 bandında KESİN YANLIŞ (ω<1 için 102a
    haklı çıkarmıştı; burada değil).
  • ORAN DÜZELİYOR: taşıyıcı/taban (Hann, frekans-eşlenik) medyan 4.0×
    (3.0–5.1) / 5.0× (3.1–7.2). Tapersiz frekans-eşlenik de 3.5–7.3.
    87'nin "15×" ucu, yüksek-ω taşıyıcılarını (Σω, |Ĝ|≈5–9e-3) tek
    düşük-ω tabanına (Δω'daki 5e-4) bölmesinin eseriydi: Σω civarının
    kendi tabanı ~1e-3, oran oraya da ~4–7×.
  • KİMLİK SIFIRLARIN DEĞİL (P4): aynı frekanslarda sıfır ızgarası z_n
    (Hann) medyan |Ĝ_z| = 3.7e-7 / 1.7e-7 → orta-nokta/sıfır oranı
    9.3e3× / 1.6e4× (I'da ~9e7 / 2.6e8, 102b ✓). Taşıyıcı içeriği
    orta-nokta örneklemesinin (t_mid = z + g/2) malı.
GENEL HÜKÜM (görevin üç sorusu):
  1) Taşıyıcılar taperli tabanın üstünde DURUYOR; çökme yok. 87-T3'ün
     sayısal iskeleti düzeltilmiş tabanla hayatta: oran bandı 3–7×.
  2) TDS ATFI YİNE DE DÜZELTİLMELİ: "sıfır kristalinin termal difüz
     saçılması" değil; sıfır tayfı aynı frekanslarda ~10⁴ kat (|Ĝ|)
     karanlık. Doğru cümle: "bükülmeyi taşıyan şey, orta-nokta
     ızgarasının asal-çizgi yan bantlarından oluşan difüz alanıdır"
     (102b kimliği; burada frekans-frekans çapraz doğrulandı).
  3) T1 (Gram=Ĝ) ve T2 (transfer) tabana hiç yaslanmıyor; taper
     ikisini de değiştirmiyor. 87'nin bu iki hükmü OLDUĞU GİBİ AYAKTA.
DÜRÜST KAYITLAR: (i) ön-mühür 102b bilinerek yazıldı — P3/P4 bağımsız
  öngörü değil çapraz doğrulama; (ii) taban medyanı dışlama yarıçapına
  (çizgi ±4, taşıyıcı ±6 çözünürlük) ve 87'nin ±0.04 pencere seçimine
  duyarlı — oranlar ±%30 oynayabilir, 3–7× bandı sağlam; (iii) sıfır
  ızgarası z_n gaps'ten kümülatif kuruldu (102b ile aynı); mutlak sıfır
  konumu değil gap dizisi kullanılıyor, ω≈4–7'de bu yeterli.
"""

import numpy as np
from pathlib import Path
from sympy import primerange
import time

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
            out.append(float(np.log(q)))
            q *= p
    return np.array(sorted(out))

def hann(t):
    return 0.5 * (1 - np.cos(TWO_PI * (t - t[0]) / (t[-1] - t[0])))

def unfold(gaps, amps, tmid):
    Lw = np.log(tmid / TWO_PI)
    g_u = gaps * Lw / TWO_PI
    a_u = amps / np.sqrt(A_CG * Lw + B0 + B1 / Lw)
    a_u /= np.sqrt((a_u**2).mean())
    return np.log(a_u), g_u, float(Lw.mean())

def Ghat_iki(t, omegas, w, center=True, chunk=192):
    """Tek geçişte iki kestirimci: (|Ĝ_1|, Ĝ_1, |Ĝ_w|, Ĝ_w).
    Ĝ_1 = Σ e^{iωt}/n (87'nin tapersizi), Ĝ_w = Σ w e^{iωt}/Σw (Hann)."""
    omegas = np.asarray(omegas, float)
    tc = t - 0.5 * (t[0] + t[-1]) if center else t
    sw = float(w.sum())
    G1 = np.empty(len(omegas), complex)
    Gw = np.empty(len(omegas), complex)
    for s0 in range(0, len(omegas), chunk):
        ob = omegas[s0:s0 + chunk]
        E = np.exp(1j * np.outer(ob, tc))
        G1[s0:s0 + chunk] = E.sum(axis=1) / len(t)
        Gw[s0:s0 + chunk] = (E * w[None, :]).sum(axis=1) / sw
    return G1, Gw

def rvm_N(t):
    x = t / TWO_PI
    return x * np.log(x / np.e) + 7 / 8

def smooth_zeros(z0, m):
    kk = np.arange(m, dtype=float)
    zs = z0 + kk * TWO_PI / np.log(z0 / TWO_PI)
    for _ in range(8):
        fdel = rvm_N(zs) - rvm_N(z0) - kk
        zs = zs - fdel / (np.log(zs / TWO_PI) / TWO_PI)
    return zs

t0 = time.time()
LOGQ_FULL = prime_powers(int(np.exp(11.2)))
print(f"[87b] tam çizgi listesi: {len(LOGQ_FULL)} asal-kuvvet, "
      f"log q ≤ {LOGQ_FULL.max():.3f}", flush=True)

d41 = np.load(HERE / "41_bigT_windows.npz")
K41 = sorted({x.split("_")[1] for x in d41.files}, key=lambda s: int(s[:-1]))
OUT = {}

for kkey in [K41[1], K41[0]]:
    gaps, amps, tmid = d41[f"gaps_{kkey}"], d41[f"amps_{kkey}"], d41[f"tmid_{kkey}"]
    ya, g_u, L = unfold(gaps, amps, tmid)
    n = len(ya)
    res_T = TWO_PI / (tmid[-1] - tmid[0])
    band = [p for p in primerange(int(np.exp(0.505*L)), int(np.exp(0.55*L)))][:10]
    kucuk = [p for p in P11 if np.log(p) / L < 0.45]
    w_h = hann(tmid)
    n_eff_h = w_h.sum()**2 / (w_h**2).sum()
    print(f"\n{'='*74}\nPENCERE {kkey}: L={L:.4f}  n={n}  2π/T={res_T:.3e}  "
          f"n_eff(Hann)={n_eff_h:.0f}", flush=True)

    def cols(ps):
        c = []
        for p in ps:
            arg = tmid * np.log(p)
            c += [np.cos(arg), np.sin(arg)]
        return np.vstack(c).T

    Xb = np.hstack([np.ones((n, 1)), cols(band)])
    Xs = cols(kucuk)
    yc = ya - ya.mean()

    # ---------------- D1: T1 kimliği, tapersiz + Hann-ağırlıklı ----------------
    need = set()
    for p in band:
        for q in kucuk:
            need.add(round(np.log(p) - np.log(q), 10))
            need.add(round(np.log(p) + np.log(q), 10))
    need = sorted(need)
    # kimlik cols(t) ile aynı fazda kurulmalı: MERKEZLEMESİZ Ĝ
    G1n, Gwn = Ghat_iki(tmid, np.array(need), w_h, center=False)
    Gv1 = dict(zip(need, G1n))
    Gvw = dict(zip(need, Gwn))

    def gram_err(weights, Gv):
        sw = float(weights.sum())
        G_dir = (Xb[:, 1:].T @ (weights[:, None] * Xs)) / sw
        G_rec = np.zeros_like(G_dir)
        for i, p in enumerate(band):
            for j, q in enumerate(kucuk):
                d_ = Gv[round(np.log(p) - np.log(q), 10)]
                s_ = Gv[round(np.log(p) + np.log(q), 10)]
                G_rec[2*i, 2*j]     = 0.5 * (d_.real + s_.real)
                G_rec[2*i+1, 2*j+1] = 0.5 * (d_.real - s_.real)
                G_rec[2*i, 2*j+1]   = 0.5 * (s_.imag - d_.imag)
                G_rec[2*i+1, 2*j]   = 0.5 * (d_.imag + s_.imag)
        return np.abs(G_dir - G_rec).max(), np.abs(G_dir).max()

    e1, s1 = gram_err(np.ones(n), Gv1)
    ew, sw_ = gram_err(w_h, Gvw)
    print(f"D1 T1 KİMLİK: tapersiz hata {e1:.2e} (ölçek {s1:.4f}) | "
          f"Hann-ağırlıklı hata {ew:.2e} (ölçek {sw_:.4f}) → "
          f"{'✓ kimlik her iki kestirimcide' if max(e1, ew) < 1e-6 else 'KONTROL!'}",
          flush=True)

    # ---------------- D2: T2 transferi + kapanış artığı ----------------
    bb_bare, *_ = np.linalg.lstsq(Xb, yc, rcond=None)
    bf, *_ = np.linalg.lstsq(np.hstack([Xb, Xs]), yc, rcond=None)
    b_true_b, b_s = bf[:Xb.shape[1]], bf[Xb.shape[1]:]
    T = np.linalg.solve(Xb.T @ Xb, Xb.T @ Xs)
    art = np.abs(bb_bare - (b_true_b + T @ b_s)).max()
    def Bhat(vec):
        return np.array([vec[1 + 2*i] / p**-0.5
                         for i, p in enumerate(band)]).mean()
    print(f"D2 T2 TRANSFER: B̂ çıplak {Bhat(bb_bare):+.4f} | tam "
          f"{Bhat(b_true_b):+.4f} | kapanış artığı max|b_ç−(b_t+T·b_s)| = "
          f"{art:.2e} → {'✓ kesin cebir' if art < 1e-8 else 'KONTROL!'}",
          flush=True)
    print("   kod denetimi: T2'de Ĝ/taban terimi YOK — yalnız X'X cebiri "
          "(tabandan bağımsız).", flush=True)

    # ---------------- D3: T3 taşıyıcı vs taperli yerel taban ----------------
    pairs = []
    for i, p in enumerate(band):
        for j, q in enumerate(kucuk):
            contrib = (T[1+2*i, 2*j] * b_s[2*j]
                       + T[1+2*i, 2*j+1] * b_s[2*j+1]) / p**-0.5
            pairs.append((abs(contrib), contrib, p, q))
    pairs.sort(reverse=True)
    freqs, etiket = [], []
    for _, c_, p, q in pairs[:8]:
        freqs += [np.log(p) - np.log(q), np.log(p) + np.log(q)]
        etiket += [f"Δ({p},{q})", f"Σ({p},{q})"]
    freqs = np.array(freqs)

    # ızgaralar: sıfırlar ve pürüzsüz ortalar (sızıntı kontrolü)
    z = np.empty(n + 1)
    z[0] = tmid[0] - gaps[0] / 2
    z[1:] = z[0] + np.cumsum(gaps)
    zs = smooth_zeros(z[0], n + 1)
    ts_sm = 0.5 * (zs[:-1] + zs[1:])
    w_z = hann(z)
    w_sm = hann(ts_sm)

    # taşıyıcı değerleri
    Gc1, Gcw = Ghat_iki(tmid, freqs, w_h)
    Gz1, Gzw = Ghat_iki(z, freqs, w_z)
    Gs1, Gsw = Ghat_iki(ts_sm, freqs, w_sm)

    # taban taramaları (çizgi-dışlamalı medyan)
    step = 2 * res_T
    off = np.arange(-0.04, 0.04 + step / 2, step)
    tab = np.zeros((len(freqs), 4))          # ham, hann, M-S ham, M-S hann
    for i, om in enumerate(freqs):
        scan = om + off
        m = np.ones(len(scan), bool)
        for l in LOGQ_FULL[(LOGQ_FULL > scan[0] - 0.01) & (LOGQ_FULL < scan[-1] + 0.01)]:
            m &= np.abs(scan - l) > 4 * res_T
        m &= np.abs(scan - om) > 6 * res_T
        S1, Sw = Ghat_iki(tmid, scan[m], w_h)
        M1, Mw = Ghat_iki(ts_sm, scan[m], w_sm)
        tab[i] = [np.median(np.abs(S1)), np.median(np.abs(Sw)),
                  np.median(np.abs(M1)), np.median(np.abs(Mw))]
        print(f"   taban {i+1:2d}/{len(freqs)} ω={om:.3f}  "
              f"({m.sum()}/{len(scan)} nokta)  "
              f"Hann taban={tab[i,1]:.2e}  M-S={tab[i,3]:.2e}", flush=True)

    print(f"\nD3 TABLO — |Ĝ| birimi (Ĝ_w = Σw e^{{iωt}}/Σw):")
    print(f"{'çift':>11} {'ω':>7} {'|Ĝ|ham':>9} {'taban_h':>9} {'oran_h':>7} "
          f"{'|Ĝ|Hann':>9} {'tabanHan':>9} {'oranHan':>8} {'M-S_Han':>9} "
          f"{'|Ĝ_z|Han':>9}")
    oranlar = []
    for i, om in enumerate(freqs):
        oh = np.abs(Gc1[i]) / tab[i, 0]
        ot = np.abs(Gcw[i]) / tab[i, 1]
        oranlar.append(ot)
        print(f"{etiket[i]:>11} {om:>7.3f} {np.abs(Gc1[i]):>9.2e} "
              f"{tab[i,0]:>9.2e} {oh:>7.1f} {np.abs(Gcw[i]):>9.2e} "
              f"{tab[i,1]:>9.2e} {ot:>8.1f} {tab[i,3]:>9.2e} "
              f"{np.abs(Gzw[i]):>9.2e}")
    oranlar = np.array(oranlar)
    I_ham = n * np.abs(Gc1)**2
    I_han = n_eff_h * np.abs(Gcw)**2
    print(f"\nD3 ÖZET [{kkey}]:")
    print(f"  Hann/tapersiz |Ĝ|² oranı (medyan)     : "
          f"{np.median(np.abs(Gcw)**2 * n_eff_h / (np.abs(Gc1)**2 * n)):.3f}  "
          f"(102b çaprazı: I medyan ham {np.median(I_ham):.3f}, "
          f"Hann {np.median(I_han):.3f})")
    print(f"  taşıyıcı |Ĝ| Hann (medyan)            : "
          f"{np.median(np.abs(Gcw)):.2e}  (aralık "
          f"{np.abs(Gcw).min():.1e}–{np.abs(Gcw).max():.1e})")
    print(f"  taperli çizgi-dışlamalı taban (medyan): {np.median(tab[:,1]):.2e}  "
          f"(87 tapersiz: 5e-4)")
    print(f"  taban / M-S sızıntı (medyan)          : "
          f"{np.median(tab[:,1] / tab[:,3]):.1e}×  (aralık "
          f"{(tab[:,1]/tab[:,3]).min():.1e}–{(tab[:,1]/tab[:,3]).max():.1e})")
    print(f"  ★ TAŞIYICI/TABAN Hann (medyan)        : {np.median(oranlar):.1f}×  "
          f"(aralık {oranlar.min():.1f}–{oranlar.max():.1f}; 87: 4–15×)")
    print(f"  ★ sıfır ızgarası |Ĝ_z| Hann (medyan)  : "
          f"{np.median(np.abs(Gzw)):.2e}  → orta/sıfır "
          f"{np.median(np.abs(Gcw) / np.abs(Gzw)):.0f}×  "
          f"(I oranı frekans başına {np.median(n_eff_h * np.abs(Gcw)**2 / (n_eff_h * np.abs(Gzw)**2)):.1e}×)",
          flush=True)
    OUT[kkey] = dict(freqs=freqs, Gc1=np.abs(Gc1), Gcw=np.abs(Gcw),
                     Gzw=np.abs(Gzw), tab=tab, oranlar=oranlar, L=L, n=n,
                     res_T=res_T, n_eff_h=n_eff_h)

np.savez(HERE / "87b_denetim.npz",
         **{f"{k}_{kk}": v for k, D_ in OUT.items()
            for kk, v in D_.items() if not isinstance(v, list)})
print(f"\n87b_denetim.npz yazıldı. Toplam {time.time()-t0:.0f} s.", flush=True)
