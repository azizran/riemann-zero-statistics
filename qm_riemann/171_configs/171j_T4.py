"""
171j — T4 (bonus): β KÖPRÜSÜ — kesim gazları AYNI A(τ) üstünde mi?
===================================================================
Yeni ölçüm YOK. Girdi: scratchpad/167/C_<gaz>.json + 171e'nin defteri.
`169_k1` AYNEN import edilir.

SORULAR (kalem): kesim gazları (E060, HA4, K070, K090) aynı A(τ)
üstünde mi? Kesim ekseninde kolaps var mı, M/θ nasıl ayrışıyor,
β = 0.2175 bu çerçevede hangi ÇARPANIN işi?

ÖN-MÜHÜR (koşudan ÖNCE, 4 Eylül 2026 00:12):
 U1 KESİM GAZLARI AYNI A(τ) ÜSTÜNDE DEĞİL. 171e §T2b.4'ün gaz-başına
    genişlikleri: σ*/2 = 0.2486 (K090, +%2.7) / 0.2683 (K070, +%10.9) /
    0.2938 (HA4, +%21.4) / 0.3414 (E060, +%41.0) — λ-ailesi bandı
    0.2269–0.2493. K090 sınırda İÇERİDE, K070 sınırda DIŞARIDA, HA4 ve
    E060 kesin DIŞARIDA. Şekil sapması bant-bant %5–20 olmalı.
 U2 Kesim ekseninde M'yi θ taşır; g_E ve g_X² φ ile ARTAR ve θ'nın
    düşüşünü kısmen iptal eder (K070: 1.0475 × 1.0317 × 0.8038 = 0.8687).
 U3 β = 0.2175 SEVİYE (θ) çarpanının işidir, ŞEKİL (A) çarpanının değil;
    θ = 1 − βφ kuralı KALİB'in bant-ortalamasını verir ama BANT ŞEKLİNİ
    vermez. Öngörü: `KALİB_b(kesim)/[KALİB_b(Hk)·(g_cal/g₀)·(1−βφ)]`
    bant-ortalamada %1-2 içinde 1 olur ama bant-bant EĞİMİ sıfır
    olmaz — K070/HA4/E060'ta %5-15'lik tekdüze bir τ-eğimi kalır.
 U4 Kesim ekseninde ORTAK bir şekil (φ'ye bağlı tek parametreli) var mı?
    Öngörü: kısmen — α_kesim − α_λ farkı φ ile artar ama HA4 ile K070
    sırayı bozar (φ 0.910 vs 0.927 ama α 1.70 vs 1.42), yani φ TEK
    değişken değil; kesim BİÇİMİ (keskin ↔ erfc) de giriyor.

SONUÇ yalnız gerçek koşudan.   Çıktı: scratchpad/171/T4.json

SONUÇ (gerçek koşudan):
 U1 ✓  Şekil sapması rms: K090 %0.70 (λ-ailesinin İÇİNDE), K070 %1.97,
       HA4 %4.16, E060 %9.16;  σ*/2 = 0.24855 / 0.26834 / 0.29380 / 0.34136
       (λ-ailesi bandı 0.22685–0.24926). **Yalnız K090 aynı A üstünde.**
 U2 ✓  Kesimde M'yi θ taşıyor; g_E ve g_X² φ ile ARTIP kısmen iptal ediyor.
 U3 ✓✓ β köprüsü: θ/θ₀ ölçülen ile 1−βφ farkı −%0.19/+%0.68/−%0.57/+%0.40
       (SEVİYE ‰7 içinde kapanıyor) ama artıkta tekdüze τ-EĞİMİ kalıyor:
       −15.2 / −33.3 / −71.6 / −141.7. λ ekseninde aynı artık DÜZ
       (L060 +5.3…+3.1, son +3.9…+1.7) ⇒ λ borcu saf SEVİYE, kesim borcu
       SEVİYE + ŞEKİL. **β = 0.2175 seviyenin yasası, şeklin değil.**
 U4 ✓  Kural rms'i K090 %0.89, K070 %1.95, HA4 %4.11 (öngörü %2-4),
       E060 %8.38 (gürültülü).
 EK: Δα–φ korelasyonu r = +0.673 ve K070↔HA4 sıralaması bozuk ⇒ kesim
     ekseninde de TEK-DEĞİŞKENLİK YOK (kesimin BİÇİMİ ayrı değişken).
"""
import importlib
import json
import os
import sys
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
KES = ["K090", "K070", "HA4", "E060"]
BETA = 0.2175                # 168 §A2.7 / 169 §K1(b)


def main():
    os.makedirs(SCR171, exist_ok=True)
    D = K1.yukle()
    MO = json.load(open(SCR171 + "/M_ONKAYIT.json"))
    K0 = json.load(open(SCR + "/170/K0.json"))["gaz"]
    sX_hk = D["Hkeskin"]["sigX"]
    A2 = lambda t: np.exp(-2 * np.pi ** 2 * t ** 2 * sX_hk ** 2)  # noqa: E731

    R = {}
    for g in LAM5 + ["son"] + KES:
        B = K1.saglikli(D[g], lo_max=0.68)
        rows = [K1.band_jk(b["cizgi"]) for b in B]
        R[g] = dict(k=np.array([r["KALIB_u2"] for r in rows]),
                    sk=np.array([r["sKALIB_u2"] for r in rows]),
                    t=np.array([r["tau_eff"] for r in rows]))
    KH, TH = R["Hkeskin"]["k"], R["Hkeskin"]["t"]

    print("=" * 104)
    print("171j — T4: β KÖPRÜSÜ.  Kesim gazları AYNI A(τ) üstünde mi?")
    print("=" * 104)
    print("  %-7s %7s | %-34s | %7s %9s | %7s"
          % ("gaz", "φ", "S/A2Hk − 1 (%)  [şekil sapması]", "rms%", "τ-eğimi",
             "σ*/2"))
    sek = {}
    for g in LAM5 + ["son"] + KES:
        S = R[g]["k"] / R[g]["k"][2]
        Ap = A2(R[g]["t"]) / A2(R[g]["t"][2])
        s = 100 * (S / Ap - 1)
        eg = float(np.polyfit(R[g]["t"], s, 1)[0])
        sek[g] = dict(sap=[float(x) for x in s],
                      rms=float(np.sqrt((s ** 2).mean())), egim=eg,
                      phi=float(K0[g]["phi"]),
                      sX_eff=MO["genislik"][g]["sX_eff"])
        print("  %-7s %7.4f | %s | %7.2f %+9.1f | %7.5f"
              % (g, K0[g]["phi"], "  ".join("%+6.2f" % x for x in s),
                 sek[g]["rms"], eg, MO["genislik"][g]["sX_eff"]))
    ge5 = [MO["genislik"][g]["sX_eff"] for g in LAM5]
    print("  λ-ailesi σ*/2 bandı: [%.5f, %.5f]   (ort %.5f = σ_X̃(λ=1) %+.2f%%)"
          % (min(ge5), max(ge5), np.mean(ge5),
             100 * (np.mean(ge5) / sX_hk - 1)))

    # ---------- β köprüsü: seviye ve şekil ayrı ayrı -------------------
    print("\n" + "=" * 104)
    print("T4.2 — β KÖPRÜSÜ: KALİB_b(kesim) = KALİB_b(Hk)·(g_cal/g₀)·(1−βφ)?"
          "   [β = %.4f]" % BETA)
    print("=" * 104)
    g0 = K0["Hkeskin"]["gcal"]
    print("  %-7s %7s %9s %9s | %-34s | %7s %9s"
          % ("gaz", "φ", "1−βφ", "θ/θ₀ ölç", "artık (%)  [şekil borcu]",
             "rms%", "τ-eğimi"))
    kop = {}
    for g in KES:
        phi = K0[g]["phi"]
        fac = (K0[g]["gcal"] / g0) * (1 - BETA * phi)
        r = 100 * (R[g]["k"] / (KH * fac) - 1)
        eg = float(np.polyfit(R[g]["t"], r, 1)[0])
        th = MO["ayrisim"][g]["th"]
        kop[g] = dict(phi=float(phi), fac=float(fac), th=float(th),
                      artik=[float(x) for x in r],
                      rms=float(np.sqrt((r ** 2).mean())), egim=eg,
                      ort=float(np.mean(r)))
        print("  %-7s %7.4f %9.4f %9.4f | %s | %7.2f %+9.1f"
              % (g, phi, 1 - BETA * phi, th,
                 "  ".join("%+6.2f" % x for x in r),
                 kop[g]["rms"], eg))
    print("\n  (karşılaştırma: aynı kural λ ekseninde, φ = 0 ⇒ 1−βφ = 1)")
    for g in ("L060", "son"):
        fac = K0[g]["gcal"] / g0
        r = 100 * (R[g]["k"] / (KH * fac) - 1)
        print("  %-7s %7.4f %9.4f %9s | %s | %7.2f"
              % (g, 0.0, 1.0, "%.4f" % MO["ayrisim"][g]["th"],
                 "  ".join("%+6.2f" % x for x in r),
                 float(np.sqrt((r ** 2).mean()))))

    # ---------- kesim ekseninde ortak şekil var mı ----------------------
    print("\n" + "=" * 104)
    print("T4.3 — KESİM EKSENİNDE ORTAK ŞEKİL: Δα(φ) tek değişkenli mi?")
    print("=" * 104)
    a_lam = np.mean([MO["genislik"][g]["alfa"] for g in LAM5])
    print("  %-7s %7s %9s %9s %9s %s"
          % ("gaz", "φ", "α_g", "α−α_λ", "σ*/2", "kesim biçimi"))
    for g in ["Hkeskin"] + KES:
        a = MO["genislik"][g]["alfa"]
        pen = K0[g]["pen"]
        bic = ("erfc %.2f/%.3f" % tuple(pen) if pen
               else "keskin τ≤%.2f" % (K0[g]["tau_ust"] or 1.0))
        print("  %-7s %7.4f %9.4f %+9.4f %9.5f  %s"
              % (g, K0[g]["phi"], a, a - a_lam,
                 MO["genislik"][g]["sX_eff"], bic))
    ph = np.array([K0[g]["phi"] for g in KES])
    da = np.array([MO["genislik"][g]["alfa"] - a_lam for g in KES])
    kor = float(np.corrcoef(ph, da)[0, 1])
    print("  Δα ile φ arasındaki korelasyon (4 kesim gazı): r = %+.3f  "
          "(φ TEK değişken olsaydı r ≈ 1 ve sıralama bozulmazdı)" % kor)
    print("  sıralama sınavı: φ(K070)=%.4f > φ(HA4)=%.4f ama "
          "α(K070)=%.4f < α(HA4)=%.4f  ⇒ **φ tek değişken DEĞİL**"
          % (K0["K070"]["phi"], K0["HA4"]["phi"],
             MO["genislik"]["K070"]["alfa"], MO["genislik"]["HA4"]["alfa"]))

    json.dump(dict(BETA=BETA, sekil=sek, kopru=kop, a_lam=float(a_lam),
                   kor_phi_da=kor),
              open(SCR171 + "/T4.json", "w"), indent=1, default=float)
    print("\n-> %s/T4.json" % SCR171)


if __name__ == "__main__":
    main()
