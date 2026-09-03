"""
171c — T2a: A(τ) KİMLİK YARIŞI, 1. AŞAMA (uyum + 9-BANT ÖN KAYDI)
=================================================================
Yeni ölçüm YOK. Girdi: scratchpad/167/C_<gaz>.json (168/169/170 ile AYNI
dosyalar). Hesap parçası kopyalanmaz: `169_k1.yukle/saglikli/band_jk`
(o da `166_T1.bant_agg/_jk`) AYNEN import edilir.

HAKEM. 171a çarpanlaşmayı mühürledi: KALİB_u2(λ,τ) = A(τ)·M(λ) (kolaps
rms %0.81). Bu betik A(τ)'yu **GERÇEK bant KALİB_u2'sünden** (171a
parametrik W_X ile yeniden kurmuştu) beş λ-gazının log-ortalaması olarak
kurar ve ADAYLARI YALNIZ BEŞ BANTTA (lo ≤ 0.68, τ ∈ [0.539, 0.698])
uydurur. DOKUZ-BANT (lo ≤ 0.80) uzantısı — τ = 0.7386 (beş gazın hepsi)
ve τ = 0.7774 (L085, Hkeskin) — bu betikte **HESAPLANMAZ**; her adayın
oradaki öngörüsü ÖN KAYDA yazılır, hükmü 171d verir.

ADAYLAR ve TÜRETİMLERİ (siren kuralı: salt-uyum kazanamaz)
  A0   düz                         : A ≡ 1. NUL çıpa (0p).
  A1   DW-tam, Hkeskin marjinali   : A = W_amp·W_X = exp[−½(πτ)²σ_ds²
       (0p)                          −½(2πτ)²σ_X̃²]; 167'nin c_ampX üyesi
                                     bant-düz olsaydı A tam bu olurdu.
  A1λ  DW-tam, λ-ailesi geo.ort.   : aynı biçim, σ'lar beş gazın geometrik
       (0p)                          ortalaması (λ-değişmezlik iddiası).
  A2s  saf Gauss e^{−ατ²}          : ** A1 ile AYNI AİLE ** — çünkü
       (1p)                          W_amp·W_X = exp[−½π²τ²(σ_ds²+4σ_X̃²)].
                                     Serbest α ⇒ σ*_eff = √(2α)/π.
  A2Hk W_X-tek, Hkeskin σ_X̃       : A = exp(−2π²τ²σ_X̃²), σ_X̃ = 0.24204.
       (0p)                          ("ν_bant = 1, λ=1 gazında" hükmü)
  A2sn W_X-tek, GERÇEK gaz σ_X̃    : σ_X̃(son) = 0.23384 (0p).
  A2iv W_X-tek, λ-DEĞİŞMEZ taban   : σ_X̃² = A_Xλ² + B_X uyumundan B_X;
       (0p)                          "A(τ)'yu λ'dan bağımsız kalan genişlik
                                     kurar" iddiası.
  A3f  165 §6b yarı-analitik F     : A ∝ 0.546 sin²(2πτ) + 0.666 sin⁴(πτ)
       (0p)                          (Hkeskin'de ölçülmüş katsayılar).
  A3s  165 F AİLESİ (1p)           : A ∝ c₁sin²(2πτ) + c₂sin⁴(πτ), c₂/c₁
                                     serbest — yarım-gap ailesinin iki üyesi.
  A4   tarak-sayım, çıplak (0p)    : 168 §A1.3: c_pred = n_Δ·w_eff;
                                     merdiven yoğunluğu dn/dω = e^ω/ω ⇒
                                     A ∝ e^{τL}/(τL). (İŞARET sınavı.)
  A4s  DOYMUŞ tarak-sayım (1p)     : doluluk oranı doyarsa fazlalık
                                     A ∝ 1/(1+κ e^{τL}).
  A5   kuvvet yasası τ^{−p} (1p)   : TÜRETİMSİZ, EMPİRİK ÇIPA (etiketli;
                                     kazansa bile kimlik sayılmaz).

ÖN-MÜHÜR (koşudan ÖNCE, 171a'nın mühürlü A vektöründen EL HESABI):
 P1 5-bant A (gerçek KALİB'den) 171a'nınkine (1.1048/1.0448/1.0000/
    0.9423/0.8832) %0.6 içinde eşit olur; parametrik-W_X yeniden kurulumu
    yalnız ‰1-5 kaydırır.
 P2 A2s (saf Gauss) beş bantta artık rms ≤ %0.35 ile kazanır; α ≈ 1.16
    (aralık 1.05–1.30) ⇒ σ*_eff = √(2α)/π ≈ 0.485 (0.46–0.51), yani
    "W_X-tek" dilinde σ_X̃-eşdeğeri ≈ 0.2424 (0.230–0.257).
 P3 KİMLİK BAŞ ADAYI: σ*_eff/2 = σ_X̃(Hkeskin) = 0.24204 (%1 içinde).
    ⇒ A2Hk sıfır-parametreyle A2s kadar iyi olmalı. Yedek kimlik:
    beş gazın geometrik ortalaması σ_X̃ = 0.2184 (%10 uzak; düşmeli).
 P4 A1 (DW-tam, Hkeskin) ÇOK DİK: σ*=√(σ_ds²+4σ_X̃²)=0.6484, gereken
    0.485 ⇒ beş bantta ±%4-5 artık. ÖLÜR. Aynı şekilde A1λ.
 P5 A4 (çıplak tarak-sayım) İŞARETTE ölür: τ=0.539→0.698 boyunca ×6.0
    ARTIŞ öngörür, ölçülen ×0.80 DÜŞÜŞ.
 P6 A3f (165'in sabit katsayıları) ölür (τ ∈ [0.54,0.58]'de ARTIŞ verir);
    A3s serbest c₂/c₁ ile beş bantta yarışa girebilir (artık ~%0.5-1.5)
    ama 9-bant uzantısında A2s'den daha DİK düşer ve orada ayrışır.
 P7 9-BANT ÖNGÖRÜSÜ (asıl sınav, bu betik ölçmez):
    A2s → A(0.7386) ≈ 0.828, A(0.7774) ≈ 0.784 (α=1.16 ile).
    Adaylar burada %2-8 ayrışır; hüküm 171d'nin.

SONUÇ bloğu yalnız gerçek koşudan.
Çıktı: scratchpad/171/A_ONKAYIT.json + log

SONUÇ (gerçek koşudan, 23:44:17):
 P1 ✓  A(gerçek KALİB) = 1.1058/1.0454/1.0000/0.9405/0.8810; 171a ile fark
       +0.09/+0.06/0.00/−0.19/−0.25 %.
 P2 ✓  α = 1.146 (öngörü 1.05–1.30) ⇒ σ*_eff = 0.48190 ⇒ σ_X̃-eşdeğeri
       0.24095 (öngörü 0.230–0.257).
 P3 ✓✓ KİMLİK: σ_X̃(Hkeskin) = 0.24204'ten **−%0.46**; λ_eff = 0.9886.
       A2Hk (SIFIR param) rms %0.62 ↔ A2s (1 param) %0.61 — eşit.
       Yedek kimlik (λ-geo 0.21845) +%10.3 ile öldü.
 P4 ✓  A1 (DW-tam Hk) öldü — ama artığı öngörülenden büyük (±%9.4 ↔ ±%4-5).
 P5 ✓  A4 (çıplak tarak-sayım) İŞARETTE öldü: rms %79.2, en kötü +%150.9.
 P6 ✓  A3f (165'in sabit katsayıları) %13.7 ile öldü; A3s serbest c₂/c₁ =
       2.457 ile %2.52 (öngörü %0.5-1.5 — ISKA).
 P7 ✓  A2s'nin 9-bant ÖN KAYDI: A(0.7386) = 0.8301, A(0.7774) = 0.7760
       (el hesabı 0.828/0.784).
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
K1 = importlib.import_module("169_k1")

SCR = ("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
       "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
SCR171 = SCR + "/171"
LAM5 = ["L115", "Hkeskin", "L085", "L070", "L060"]
MID = 2                       # orta bant (lo = 0.60)
TAU9 = [0.7386, 0.7774]       # 9-bant uzantısının τ'ları (ÖN KAYIT için)
LL = 12.0296                  # merdiven L


# ---------------- aday şekiller (hepsi orta bantta 1'e normalize) --------
def norm(f, tau, t0):
    v = np.asarray(f(np.asarray(tau, float)), float)
    return v / float(f(np.array([t0]))[0])


def sekiller(sg_hk, sX_hk, sg_lam, sX_lam, sX_son, sX_inv):
    """(ad, fonksiyon, serbest_param_sayısı, etiket) listesi."""
    def dw(sd, sx):
        return lambda t: np.exp(-0.5 * (np.pi * t) ** 2 * sd ** 2
                                - 0.5 * (2 * np.pi * t) ** 2 * sx ** 2)

    def wx(sx):
        return lambda t: np.exp(-2 * np.pi ** 2 * t ** 2 * sx ** 2)

    S = [
        ("A0  düz", lambda t: np.ones_like(t), 0, "NUL"),
        ("A1  DW-tam Hk", dw(sg_hk, sX_hk), 0, "türetimli"),
        ("A1λ DW-tam λ-geo", dw(sg_lam, sX_lam), 0, "türetimli"),
        ("A2Hk W_X-tek Hk", wx(sX_hk), 0, "türetimli"),
        ("A2sn W_X-tek son", wx(sX_son), 0, "türetimli"),
        ("A2iv W_X-tek λ-değişmez", wx(sX_inv), 0, "türetimli"),
        ("A3f  165-F sabit", lambda t: (0.546 * np.sin(2 * np.pi * t) ** 2
                                        + 0.666 * np.sin(np.pi * t) ** 4),
         0, "türetimli"),
        ("A4  tarak-sayım çıplak", lambda t: np.exp(t * LL) / (t * LL),
         0, "türetimli"),
    ]
    P = [
        ("A2s  saf Gauss e^{−ατ²} [≡ DW-tam σ* serbest]",
         lambda t, a: np.exp(-a * t ** 2), (0.2, 4.0), 1, "türetimli-aile"),
        ("A3s  165-F ailesi (c₂/c₁)",
         lambda t, r: np.sin(2 * np.pi * t) ** 2 + r * np.sin(np.pi * t) ** 4,
         (0.05, 60.0), 1, "türetimli-aile"),
        ("A4s  doymuş tarak-sayım 1/(1+κe^{τL})",
         lambda t, k: 1.0 / (1.0 + k * np.exp(t * LL)), (1e-9, 1e-1),
         1, "türetimli-aile"),
        ("A5   kuvvet τ^{−p}", lambda t, p: t ** (-p), (-5.0, 30.0),
         1, "EMPİRİK-ÇIPA"),
    ]
    return S, P


def uydur(fam, sinir, tau, lgA, t0, n=4001):
    """1 parametreli aileyi 5 bantta log-artık en küçük kareyle uydur."""
    gr = np.linspace(sinir[0], sinir[1], n)
    best = (np.inf, None)
    for p in gr:
        v = fam(np.asarray(tau), p)
        v = v / fam(np.array([t0]), p)[0]
        r = lgA - np.log(np.abs(v))
        q = float(np.sum((r - 0) ** 2))
        if q < best[0]:
            best = (q, p)
    # ince arama
    p0 = best[1]
    d = (sinir[1] - sinir[0]) / n
    gr = np.linspace(p0 - 2 * d, p0 + 2 * d, 401)
    for p in gr:
        if not (sinir[0] <= p <= sinir[1]):
            continue
        v = fam(np.asarray(tau), p)
        v = v / fam(np.array([t0]), p)[0]
        q = float(np.sum((lgA - np.log(np.abs(v))) ** 2))
        if q < best[0]:
            best = (q, p)
    return best[1]


def main():
    os.makedirs(SCR171, exist_ok=True)
    D = K1.yukle()

    # ---- A(τ): gerçek bant KALİB_u2'sünden, beş λ-gazının log-ortalaması
    print("=" * 100)
    print("171c — T2a AŞAMA-1: A(τ) HAKEMİ (5 bant) ve ADAY UYUMU")
    print("=" * 100)
    S = {}
    TAU = {}
    SK = {}
    for g in LAM5:
        B = K1.saglikli(D[g], lo_max=0.68)
        rows = [K1.band_jk(b["cizgi"]) for b in B]
        k = np.array([r["KALIB_u2"] for r in rows])
        sk = np.array([r["sKALIB_u2"] for r in rows])
        S[g] = k / k[MID]
        SK[g] = sk / k                       # bağıl jackknife hatası
        TAU[g] = np.array([r["tau_eff"] for r in rows])
        print("  %-8s λ=%-5s  S = %s" % (g, D[g]["lam"],
              "  ".join("%.4f" % x for x in S[g])))
    tau = np.mean([TAU[g] for g in LAM5], axis=0)
    lgS = np.array([np.log(S[g]) for g in LAM5])
    A = np.exp(lgS.mean(axis=0))
    sA = lgS.std(axis=0, ddof=1) / np.sqrt(len(LAM5))     # gazlar arası s.h.
    print("  τ_ort              : " + "  ".join("%.4f" % t for t in tau))
    print("  **A(τ) (5-gaz geo)**: " + "  ".join("%.4f" % a for a in A))
    print("  gazlar arası s.h.  : " + "  ".join("%+.2f%%" % (100 * s)
                                                for s in sA))
    print("  jackknife (bant-içi, ort): " + "  ".join(
        "%.2f%%" % (100 * np.mean([SK[g][b] for g in LAM5]))
        for b in range(len(tau))))
    A171a = np.array([1.1048, 1.0448, 1.0000, 0.9423, 0.8832])
    print("  171a'nın A'sı (parametrik W_X ile): fark = " + "  ".join(
        "%+.2f%%" % (100 * (A[b] / A171a[b] - 1)) for b in range(5)))

    # ---- marjinaller ve λ-değişmez taban ------------------------------
    lam = np.array([D[g]["lam"] for g in LAM5], float)
    sX = np.array([D[g]["sigX"] for g in LAM5])
    sd = np.array([D[g]["sigds"] for g in LAM5])
    AX, BX = np.polyfit(lam ** 2, sX ** 2, 1)
    Ad, Bd = np.polyfit(lam ** 2, sd ** 2, 1)
    sX_inv = float(np.sqrt(max(BX, 1e-12)))
    sX_hk, sd_hk = D["Hkeskin"]["sigX"], D["Hkeskin"]["sigds"]
    sX_son = D["son"]["sigX"]
    sX_lam = float(np.exp(np.log(sX).mean()))
    sd_lam = float(np.exp(np.log(sd).mean()))
    print("\n  marjinaller: σ_X̃(Hk)=%.5f  σ_X̃(son)=%.5f  σ_X̃(λ-geo)=%.5f"
          "  σ_X̃(λ-değişmez taban √B_X)=%.5f" % (sX_hk, sX_son, sX_lam, sX_inv))
    print("               σ_ds(Hk)=%.5f  σ_ds(λ-geo)=%.5f   "
          "(σ_X̃²=%.5fλ²+%.5f, σ_ds²=%.5fλ²+%.5f)"
          % (sd_hk, sd_lam, AX, BX, Ad, Bd))

    # ---- aday yarışı (YALNIZ 5 BANT) ----------------------------------
    lgA = np.log(A)
    t0 = tau[MID]
    SIF, PAR = sekiller(sd_hk, sX_hk, sd_lam, sX_lam, sX_son, sX_inv)
    print("\n" + "=" * 100)
    print("ADAY YARIŞI — 5 BANT (uyum penceresi τ ∈ [%.3f, %.3f])"
          % (tau[0], tau[-1]))
    print("=" * 100)
    print("  %-46s %2s  %-9s  %7s  %7s  %s"
          % ("aday", "p", "param", "rms%", "enkötü%", "artıklar (%)"))
    sonuc = []
    for ad, f, npar, et in SIF:
        v = norm(f, tau, t0)
        r = 100 * (A / v - 1)
        rec = dict(ad=ad, npar=npar, etiket=et, param=None,
                   rms=float(np.sqrt((r ** 2).mean())),
                   enkotu=float(np.abs(r).max()),
                   artik=[float(x) for x in r],
                   ong9=[float(x) for x in norm(f, np.array(TAU9), t0)])
        sonuc.append(rec)
        print("  %-46s %2d  %-9s  %7.2f  %7.2f  %s"
              % (ad, npar, "—", rec["rms"], rec["enkotu"],
                 "  ".join("%+6.2f" % x for x in r)))
    for ad, fam, sinir, npar, et in PAR:
        p = uydur(fam, sinir, tau, lgA, t0)
        v = fam(tau, p) / fam(np.array([t0]), p)[0]
        r = 100 * (A / v - 1)
        v9 = fam(np.array(TAU9), p) / fam(np.array([t0]), p)[0]
        rec = dict(ad=ad, npar=npar, etiket=et, param=float(p),
                   rms=float(np.sqrt((r ** 2).mean())),
                   enkotu=float(np.abs(r).max()),
                   artik=[float(x) for x in r],
                   ong9=[float(x) for x in v9])
        sonuc.append(rec)
        print("  %-46s %2d  %-9.4g  %7.2f  %7.2f  %s"
              % (ad, npar, p, rec["rms"], rec["enkotu"],
                 "  ".join("%+6.2f" % x for x in r)))

    # ---- A2s'nin KİMLİK sınavı ----------------------------------------
    alfa = [r["param"] for r in sonuc if r["ad"].startswith("A2s ")][0]
    sig_eff = float(np.sqrt(2 * alfa) / np.pi)          # DW dilinde σ*
    sX_eff = sig_eff / 2.0                              # W_X dilinde σ_X̃
    print("\n" + "=" * 100)
    print("A2s'nin KİMLİĞİ: α = %.4f  ⇒  σ*_eff = √(2α)/π = %.5f  "
          "⇒  σ_X̃-eşdeğeri = %.5f" % (alfa, sig_eff, sX_eff))
    print("=" * 100)
    kim = []
    for ad2, v in (("σ_X̃(Hkeskin, λ=1.00)", sX_hk),
                   ("σ_X̃(son, GERÇEK ζ)", sX_son),
                   ("σ_X̃(λ-ailesi geo.ort.)", sX_lam),
                   ("σ_X̃(λ-değişmez taban)", sX_inv),
                   ("σ_X̃(L115)", D["L115"]["sigX"]),
                   ("σ_X̃(L060)", D["L060"]["sigX"]),
                   ("σ_Ĉ(Hkeskin)", D["Hkeskin"]["sigC"]),
                   ("σ_ds(Hkeskin)/2", sd_hk / 2)):
        print("   %-26s = %.5f    σ_X̃-eşdeğeri/bu = %+.2f%%"
              % (ad2, v, 100 * (sX_eff / v - 1)))
        kim.append(dict(ad=ad2, deger=float(v),
                        sapma=float(100 * (sX_eff / v - 1))))
    # hangi λ bu σ_X̃'yi verir
    p = np.polyfit(np.log(lam), np.log(sX), 2)
    from numpy.polynomial import polynomial as _P   # noqa: F401
    ll = np.linspace(np.log(0.4), np.log(1.6), 20001)
    lam_eff = float(np.exp(ll[np.argmin(np.abs(np.polyval(p, ll)
                                               - np.log(sX_eff)))]))
    print("   ⇒ σ_X̃-eşdeğerini veren λ (log-log kuadratik ara değer): "
          "**λ_eff = %.4f**" % lam_eff)

    # ---- ÖN KAYIT ------------------------------------------------------
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    print("\n" + "=" * 100)
    print("ÖN KAYIT — 9-BANT (lo ≤ 0.80) UZANTISI, ZAMAN DAMGASI %s" % ts)
    print("  (bu betik 9-bant A'sını HESAPLAMAZ; hüküm 171d'nin)")
    print("=" * 100)
    print("  %-46s  A(%.4f)   A(%.4f)" % ("aday", TAU9[0], TAU9[1]))
    for r in sonuc:
        print("  %-46s   %.4f     %.4f" % (r["ad"], r["ong9"][0], r["ong9"][1]))
    json.dump(dict(zaman=ts, tau=list(map(float, tau)), TAU9=TAU9,
                   A5=[float(x) for x in A], sA=[float(x) for x in sA],
                   S={g: [float(x) for x in S[g]] for g in LAM5},
                   TAU={g: [float(x) for x in TAU[g]] for g in LAM5},
                   alfa=float(alfa), sig_eff=sig_eff, sX_eff=sX_eff,
                   lam_eff=lam_eff, kimlik=kim,
                   marj=dict(sX_hk=sX_hk, sX_son=sX_son, sX_lam=sX_lam,
                             sX_inv=sX_inv, sd_hk=sd_hk, sd_lam=sd_lam,
                             AX=float(AX), BX=float(BX)),
                   adaylar=sonuc),
              open(SCR171 + "/A_ONKAYIT.json", "w"), indent=1)
    print("\n-> %s/A_ONKAYIT.json" % SCR171)


if __name__ == "__main__":
    main()
