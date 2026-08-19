"""
87 — GİYDİRMENİN KAPALI DEVRESİ: BÜKÜLME NEREDE YAŞIYOR? (20 Ağustos, gece)
==========================================================================
81 giydirmeyi ölçtü (kat-dibi çıplak okuma, küçük asallar eklenince
−0.354 → −0.126); 86 Ĝ(ω)'yı ölçtü. Bu script halkayı kapatır:

  KİMLİK: iki kolonun Gram'ı Ĝ'den kesin olarak gelir —
    ⟨cos(ω_a t) cos(ω_b t)⟩ = ½ Re[Ĝ(ω_a−ω_b) + Ĝ(ω_a+ω_b)]
    (sin/cos karışımları benzer, Im'lerle). T1 bunu sayısal doğrular.
  AYRIŞTIRMA: atlanmış-değişken transferi kesin cebirdir —
    b_çıplak = b_tam + (Xb'Xb)^{-1} Xb'Xs · b_s.
    T2 transferi atlanan çizgi çizgi ayrıştırır: giydirmeyi kim taşıyor?
  TAŞIYICI KİMLİĞİ: T3 en büyük taşıyıcıların (Δω, Σω) frekanslarında
    Ĝ'yi tarar — çizgiye mi, tarağa mı yakın, yoksa karanlık-alan
    tabanının üstünde yerel tepe mi? (86-T4 gerilimi: taban 2e-4,
    kaba gereksinim ~6e-3 — bükülme düz tabanda YAŞAYAMAZ.)

Kurulum 81-L1 rungunun aynısı: band asalları [0.505, 0.55]·L, küçük
küme P11 (τ<0.45); KONTROLSÜZ saf-dalga devresi (temiz Ĝ hikâyesi).

SONUÇ (iki pencere): T1 kimlik ✓ (hata ~5e-8, Gram ölçeği 0.006);
T2 transfer cebiri birebir kapanıyor (+0.429 / +0.450) ve katkılar
TEK taşıyıcıda değil — 11 atlanan çizginin HEPSİNDE dağılmış, hepsi
pozitif, q ile yavaşça azalan (KOLEKTİF transfer). T3 taşıyıcı kimliği:
çift başına |Ĝ| 0.002-0.009, p^k çizgilerinde DEĞİL, ama yerel karanlık
tabanın (5e-4) 4-15 katı ve ω arttıkça tarağa (ω=L) doğru yükseliyor →
bükülmeyi taşıyan şey KRİSTALİN TERMAL DİFÜZ SAÇILMASI (TDS): hiperuniform
karanlık ile Bragg tarağı arasındaki termal omuz. Kaba 1B termal-kristal
kestirimi |Ĝ|_difüz ~ √((1−e^{−ω²σ²})/2n) ≈ 0.0035 @ω=6.5 — ölçülen
0.005-0.009 ile aynı mertebe. Açık: kontrollü-rung kontrol-aracılı payı;
TDS'nin sıcaklık-ölçeklemesi (pencereler arası); tam difüz-model fiti.
"""

import numpy as np
from pathlib import Path
from sympy import primerange

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
A_CG = (np.e**2 - 5) / 2
B0, B1 = 2.7580, -0.0543
P11 = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]

def unfold(gaps, amps, tmid):
    Lw = np.log(tmid / TWO_PI)
    g_u = gaps * Lw / TWO_PI
    a_u = amps / np.sqrt(A_CG * Lw + B0 + B1 / Lw)
    a_u /= np.sqrt((a_u**2).mean())
    return np.log(a_u), g_u, float(Lw.mean())

def Ghat(t, omegas, chunk=40000):
    out = np.zeros(len(omegas), dtype=complex)
    for s0 in range(0, len(t), chunk):
        tt = t[s0:s0 + chunk]
        out += np.exp(1j * np.outer(omegas, tt)).sum(axis=1)
    return out / len(t)

d41 = np.load(HERE / "41_bigT_windows.npz")
K41 = sorted({x.split("_")[1] for x in d41.files}, key=lambda s: int(s[:-1]))

for kkey in [K41[1], K41[0]]:
    gaps, amps, tmid = d41[f"gaps_{kkey}"], d41[f"amps_{kkey}"], d41[f"tmid_{kkey}"]
    ya, g_u, L = unfold(gaps, amps, tmid)
    band = [p for p in primerange(int(np.exp(0.505*L)), int(np.exp(0.55*L)))][:10]
    kucuk = [p for p in P11 if np.log(p) / L < 0.45]
    n = len(ya)
    print(f"\n================ PENCERE L = {L:.2f} (n={n}) ================")

    def cols(ps):
        c = []
        for p in ps:
            arg = tmid * np.log(p)
            c += [np.cos(arg), np.sin(arg)]
        return np.vstack(c).T

    Xb = np.hstack([np.ones((n, 1)), cols(band)])
    Xs = cols(kucuk)
    yc = ya - ya.mean()

    # ---- T1: Gram kimliği (Ĝ'den yeniden kurulum)
    G_dir = (Xb[:, 1:].T @ Xs) / n
    oms_b = np.repeat([np.log(p) for p in band], 2)
    oms_s = np.repeat([np.log(q) for q in kucuk], 2)
    need = set()
    for ob in [np.log(p) for p in band]:
        for os_ in [np.log(q) for q in kucuk]:
            need.add(round(ob - os_, 10)); need.add(round(ob + os_, 10))
    need = sorted(need)
    Gv = dict(zip(need, Ghat(tmid, np.array(need))))
    G_rec = np.zeros_like(G_dir)
    for i, p in enumerate(band):
        for j, q in enumerate(kucuk):
            d_ = Gv[round(np.log(p) - np.log(q), 10)]
            s_ = Gv[round(np.log(p) + np.log(q), 10)]
            G_rec[2*i, 2*j]     = 0.5 * (d_.real + s_.real)   # cos·cos
            G_rec[2*i, 2*j+1]   = 0.5 * (s_.imag - d_.imag)   # cos·sin
            G_rec[2*i+1, 2*j]   = 0.5 * (d_.imag + s_.imag)   # sin·cos
            G_rec[2*i+1, 2*j+1] = 0.5 * (d_.real - s_.real)   # sin·sin
    err = np.abs(G_dir - G_rec).max()
    print(f"T1 KİMLİK: maks |Gram_direkt − Gram_Ĝ| = {err:.2e}  "
          f"(Gram ölçeği {np.abs(G_dir).max():.4f}) → {'✓ kimlik' if err < 1e-6 else 'kontrol!'}")

    # ---- T2: transfer ayrıştırması (kontrolsüz devre)
    bb_bare, *_ = np.linalg.lstsq(Xb, yc, rcond=None)
    Xf = np.hstack([Xb, Xs])
    bf, *_ = np.linalg.lstsq(Xf, yc, rcond=None)
    b_true_b, b_s = bf[:Xb.shape[1]], bf[Xb.shape[1]:]
    T = np.linalg.solve(Xb.T @ Xb, Xb.T @ Xs)      # transfer matrisi
    bias = T @ b_s
    def Bhat(vec):
        ws = np.array([vec[1 + 2*i] / p**-0.5 for i, p in enumerate(band)])
        return ws.mean()
    print(f"T2 TRANSFER: B̂ çıplak = {Bhat(bb_bare):+.4f} | tam = {Bhat(b_true_b):+.4f} "
          f"| cebirsel fark = {Bhat(bb_bare) - Bhat(b_true_b):+.4f} "
          f"(≡ bias {Bhat(np.concatenate([[0], (T @ b_s)[1:]])):+.4f})")
    katki = []
    for j, q in enumerate(kucuk):
        vj = np.zeros_like(b_s); vj[2*j] = b_s[2*j]; vj[2*j+1] = b_s[2*j+1]
        katki.append((q, Bhat(np.concatenate([[0], (T @ vj)[1:]]))))
    katki.sort(key=lambda x: -abs(x[1]))
    print("  atlanan-çizgi katkıları (B̂-birimi):")
    for q, k in katki:
        print(f"    q={q:>2}: {k:+.4f}")

    # ---- T3: en büyük taşıyıcı çiftlerin frekans kimliği
    print("T3 TAŞIYICILAR (en büyük 8 |T·b| çifti):")
    pairs = []
    for i, p in enumerate(band):
        for j, q in enumerate(kucuk):
            contrib = (T[1+2*i, 2*j] * b_s[2*j] + T[1+2*i, 2*j+1] * b_s[2*j+1]) / p**-0.5
            pairs.append((abs(contrib), contrib, p, q))
    pairs.sort(reverse=True)
    lines = [np.log(x) for x in [2,3,4,5,7,8,9,11,13,16,25,27,32,49,121,125,128]]
    for _, contrib, p, q in pairs[:8]:
        dw = np.log(p) - np.log(q); sw = np.log(p) + np.log(q)
        Gd, Gs = Gv[round(dw, 10)], Gv[round(sw, 10)]
        yak_d = min(abs(dw - l) for l in lines)
        yak_s = min(abs(sw - l) for l in lines + [L])
        print(f"    ({p:>3},{q:>2}): katkı {contrib:+.4f}  Δω={dw:.3f} |Ĝ|={abs(Gd):.4f} "
          f"(çizgiye {yak_d:.3f})  Σω={sw:.3f} |Ĝ|={abs(Gs):.4f} (çizgi/taraka {yak_s:.3f})")
    # yerel taban kıyası: taşıyıcı frekansları ± kaydırılmış
    top = pairs[0]
    dw = np.log(top[2]) - np.log(top[3])
    scan = dw + np.linspace(-0.04, 0.04, 81)
    Gs_ = np.abs(Ghat(tmid, scan))
    print(f"    en büyük taşıyıcının Δω civarı: |Ĝ(Δω)|={Gs_[40]:.4f}, "
          f"±0.04 medyan taban {np.median(Gs_):.4f}, maks {Gs_.max():.4f}")
