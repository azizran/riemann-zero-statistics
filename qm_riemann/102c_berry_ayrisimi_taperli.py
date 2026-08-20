"""
102c — BERRY AYRIŞIMININ TAPERLİ DENETİMİ (88-M3) (20 Ağustos)
==========================================================================
DENETİM işi. 88-M3, F(α)=n|Ĝ|² ince ızgarasında iki sayı verdi:
  (a) ilk çizgiden ÖNCE F̄=0.002, GUE rampasının (0.036) 18 KAT ALTI;
  (b) log2-log3 bininde çizgiler DAHİL F̄=0.102 vs rampa 0.089 —
      "rampa beneklerden geri geliyor" (Berry 1988);
  (c) çizgi-arası süreklilik her yerde 10-20 kat altta.
Üçü de TAPERSİZ ölçüldü; 101j o rejimde tabanın pencere sızıntısı
olduğunu gösterdi. 102a ayrıca ORTA NOKTA ızgarası ile SIFIR ızgarasının
karanlık alanlarının 8-10 kademe farklı olduğunu buldu.

BİRİM: I(ω) = |Σ w e^{iωt}|²/Σw².  Bu, taper'dan bağımsız olarak
form-faktörün doğru normalizasyonudur (w=1 iken n|Ĝ|² = 88'in F'i;
bin-ortalaması Montgomery rampası α'ya oturur — türetim 102a'da).

ÖLÇÜMLER
  C1  88-M3'ün ince ızgarası (adım (2π/span)/2.2, ω∈[0.10,4.4)) DÖRT
      ayarda: {orta noktalar, sıfırlar} × {tapersiz, Hann}. Ayrıca RvM
      pürüzsüz ızgara (saf sızıntı tabanı).
  C2  BİN TABLOSU: F̄(hepsi), F̄(çizgi-dışı, 88'in ±3 adım maskesi),
      F̄(çizgi-dışı, ±12 adım = Hann ana lobu için güvenli maske),
      rampa α, ve MUTLAK YASA ÖNGÖRÜSÜ
        F̄_öngörü = (L/W)·Σ_{q∈bin} |A_q|²,
        A_q(orta) = Λ(q)/(L√q)·cos(πτ)·DW,  A_q(sıfır) = Λ(q)/(L√q)·DW.
      (89'un mutlak benek yasası; serbest parametre yok.)
  C3  "18×" DENETİMİ: ilk çizgi öncesi bin (log2'nin altı — açık
      formülde ATOM YOK) her ayarda kaç kat karanlık?
  C4  log2-log3 BİNİ: "rampa dönüşü" gerçek mi? Ayrıştırma —
      binin F̄'sinin ne kadarı çizgi-üstü noktalardan geliyor, ne kadarı
      çizgi-arasından? Çizgi-arası düzey, saf sızıntının (RvM pürüzsüz)
      ve 102a'nın orta-nokta gap-difüz tabanının neresinde?

ÖN-MÜHÜR (ölçümden ÖNCE yazıldı):
  C-S1 (a)'daki 18× ALET-SINIRLI çıkacak: taperli ölçümde ilk-çizgi-
       öncesi bin en az 4-6 kademe daha aşağı inecek (orta noktalarda
       gap-difüz tabanı yüzünden ~1e-6, sıfırlarda ~1e-15).
  C-S2 (b)'deki rampa dönüşü GERÇEK ve KANAT DEĞİL: taperle bin
       ortalaması ~%10 içinde korunacak, çünkü ağırlığı çizgi-ÜSTÜ
       noktalar taşıyor; mutlak yasa öngörüsü (≈0.10) ölçümü
       tutturacak. Berry ayrışımı AYAKTA KALACAK.
  C-S3 (c)'deki "10-20 kat" çizgi-arası süreklilik sayısı DÜŞECEK:
       gerçek çizgi-arası düzey rampanın 10-20 katı değil, 10³-10¹⁴
       katı altında olacak (ızgaraya göre).
  RİSK: C-S2 yanlışsa (taperle bin ortalaması çökerse) 88-M3'ün Berry
  hükmü de düşer — o zaman rampa "geri gelmiyor" demektir.

==========================================================================
SONUÇ (20 Ağustos, koşu ~3 dk) — C-S1 ✓✓, C-S2 ✓✓ (BERRY AYAKTA), C-S3 ✓
"18×" TAMAMEN ALET; "RAMPA DÖNÜŞÜ" TAM ANLAMIYLA GERÇEK VE MUTLAK
YASAYLA %2 İÇİNDE ÖNGÖRÜLÜYOR.
==========================================================================
İÇ KAPI ✓ (I normalizasyonu doğru mu?): bin ortalamaları taperden
  BAĞIMSIZ çıktı — ORTA ham 0.1015 vs ORTA Hann 0.1013 (%0.2);
  sekiz binin hepsinde ≤%0.6. Yani I(ω)=|Σw e^{iωt}|²/Σw² gerçekten
  taper-değişmez form faktörüdür (Parseval kapısı geçildi).

C2 — BİN TABLOSU (F̄; rampa α):
  bin           α      ORTA ham  ORTA Hann  öngörü | SIFIR Hann  öngörü
  [0.10,0.65) 0.036     0.0020    5.6e-06   0.0000 |  1.7e-14   0.0000
  [0.65,1.20) 0.089     0.1015    0.1013    0.0997 |  0.1127    0.1086
  [1.20,1.75) 0.142     0.0846    0.0845    0.0803 |  0.1119    0.1016
  [1.75,2.30) 0.195     0.0808    0.0807    0.0743 |  0.1289    0.1101
  [2.30,2.85) 0.248     0.1233    0.1229    0.1033 |  0.2685    0.2061
  [2.85,3.40) 0.301     0.0758    0.0757    0.0570 |  0.2495    0.1691
  [3.40,3.95) 0.354     0.0614    0.0617    0.0370 |  0.3167    0.1866
  [3.95,4.50) 0.408     0.0475    0.0473    0.0173 |  0.4056    0.2047
  ★ SIFIR ızgarasının bin ortalaması Montgomery rampasını α ≳ 0.2'de
  %±25 içinde İZLİYOR (0.79–1.26); ORTA-NOKTA ızgarası yüksek α'da
  rampanın GERİSİNE düşüyor — 88/89'un "yüksek α açığı" tam olarak
  orta-nokta cos(πτ)+DW sönümüdür, fizik değil ÖRNEKLEME.

C3 — ★ "18× KARANLIK" DENETİMİ (ilk çizgi öncesi bin [0.10,0.65);
  log2'nin ALTI, açık formülde ATOM YOK; rampa 0.0362):
    88'in ölçümü (ORTA, tapersiz)      F̄ = 2.02e-03 → 17.9× altı
    RvM PÜRÜZSÜZ, tapersiz (SAF SIZINTI) F̄ = 2.10e-03 → 17.2× altı
    ORTA noktalar, Hann                F̄ = 5.58e-06 → 6.5e+03× altı
    SIFIRLAR, Hann                     F̄ = 1.72e-14 → 2.1e+12× altı
    RvM pürüzsüz, Hann (saf sızıntı)   F̄ = 1.49e-14 → 2.4e+12× altı
  HÜKÜM: "18×" %100 ALET-SINIRLI. Hiç aritmetik içermeyen pürüzsüz
  kontrol ızgarası AYNI SAYIYI (17.2×) veriyor. Dürüst hâli:
  orta noktalarda ≥6.5e3×, SIFIRLARDA ≥2.1e12× (bu son sayı da
  sızıntı/veri sınırlı — gerçek boşluk daha da derin olabilir).

C4 — ★ log2-log3 BİNİ: RAMPA DÖNÜŞÜ GERÇEK Mİ, KANAT MI? GERÇEK.
  bindeki çizgiler: yalnız q = 2 ve 3.   rampa α = 0.0892
    ORTA, tapersiz  F̄=0.1015  (çizgi-üstü payı %97.8)
    ORTA, Hann      F̄=0.1013  (çizgi-üstü payı %99.91)
    SIFIR, Hann     F̄=0.1127  (çizgi-üstü payı %100.00)
  MUTLAK YASA ÖNGÖRÜSÜ (89'un A_q'su, serbest parametre YOK):
    orta 0.0997 (ölçüm/öngörü = 1.016)   sıfır 0.1086 (1.038)
  Yani bin ortalaması (a) taperden etkilenmiyor, (b) tamamı iki
  ÇİZGİDEN geliyor, (c) mutlak benek yasası tarafından %2–4 içinde
  öngörülüyor. BERRY AYRIŞIMI SADECE AYAKTA DEĞİL, GÜÇLENDİ.
  KANAT SORUSU — bu bindeki çizgi-DIŞI düzey:
    ORTA (Hann)   9.66e-05  = rampanın 9.2e+02 katı altı
    SIFIR (Hann)  3.55e-07  = rampanın 2.5e+05 katı altı
    RvM pürüzsüz (Hann) 1.22e-18 (saf pencere kanadı — ihmal edilebilir)
    88'in tapersiz çizgi-dışı ölçümü (±3 adım maskesi) 0.0086
      → rampanın 10.4 katı altı (88: "10–20 kat" ✓ yeniden üretildi),
      ama aynı maskeyle saf sızıntı 0.0002 ve TAPERLİ gerçek değer
      9.7e-05: 88'in sayısı ~90× ŞİŞİK (log2/log3 ana loblarının
      sinc kuyruğu, ±3 adımlık maske yetersiz).

HÜKÜM ÖZETİ: 88-M3'ün ÜÇ sayısından biri (rampa dönüşü 0.102 vs 0.089)
GERÇEK ve şimdi mutlak yasayla türetilmiş; ikisi ("18× karanlık",
"çizgi-arası 10–20 kat altta") ALET-SINIRLI ve dürüst hâlleri
sırasıyla ≥6.5e3× / ≥2.1e12× ve ≥9.2e2× / ≥2.5e5×.

ARTEFAKT ŞÜPHELERİM (gizlemiyorum):
  • σ_t bu scriptte SIFIRLARIN pürüzsüz ızgaradan sapmasıyla ölçüldü
    (0.198); 89 orta-nokta sapmasını (0.155) kullanıyordu. DW
    konvansiyonu yüksek-α öngörülerini ~1.3× oynatır; log2-log3
    bininde DW ≈ 1 olduğu için oradaki %2 uyum bu seçimden BAĞIMSIZ.
  • ±12 adımlık çizgi-dışı maskesi Hann ana lobu (±2.2 adım yarı-genişlik
    ×2) için güvenlidir; ±3 adım (88'in maskesi) DEĞİLDİR.
  • Yüksek-α binlerinde mutlak yasa öngörüsü ölçümün ~yarısı: kuvvet
    açığı + ikinci mertebe yan bantları + DW konvansiyonu karışıyor;
    bu binler bu işin hükmüne girmiyor.
"""

import numpy as np
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sympy import primerange

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi

def prime_powers(qmax):
    out = []
    for p in primerange(2, int(qmax) + 1):
        q = p
        while q <= qmax:
            out.append((q, float(np.log(p))))
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

def rvm_N(t):
    x = t / TWO_PI
    return x * np.log(x / np.e) + 7 / 8

def smooth_zeros(z0, m):
    kk = np.arange(m, dtype=float)
    zs = z0 + kk * TWO_PI / np.log(z0 / TWO_PI)
    for _ in range(8):
        zs = zs - (rvm_N(zs) - rvm_N(z0) - kk) / (np.log(zs / TWO_PI) / TWO_PI)
    return zs

# ------------------------- veri -------------------------
d41 = np.load(HERE / "41_bigT_windows.npz")
K41 = sorted({x.split("_")[1] for x in d41.files}, key=lambda s: int(s[:-1]))
gaps, amps, tmid = (d41[f"gaps_{K41[1]}"], d41[f"amps_{K41[1]}"],
                    d41[f"tmid_{K41[1]}"])          # 88-M3'ün penceresi
L = float(np.log(tmid / TWO_PI).mean())
n = len(tmid)
z = np.empty(n + 1); z[0] = tmid[0] - gaps[0] / 2
z[1:] = z[0] + np.cumsum(gaps)
zs = smooth_zeros(z[0], n + 1)
ts_sm = 0.5 * (zs[:-1] + zs[1:])
sig_t = float((z - zs).std())
print(f"[102c] L={L:.4f}  n={n}  σ_t(sıfırlar)={sig_t:.4f}", flush=True)

span = tmid[-1] - tmid[0]
step = TWO_PI / span / 2.2
om = np.arange(0.10, 4.4, step)
print(f"[102c] ince ızgara: adım {step:.5e}, {len(om)} nokta (88 ile aynı)",
      flush=True)

PQ = prime_powers(120)
LQ = np.array([np.log(q) for q, _ in PQ])
LAM = np.array([lam for _, lam in PQ])

GR = {}
for isim, t, tp in [("orta ham", tmid, False), ("orta Hann", tmid, True),
                    ("sıfır ham", z, False), ("sıfır Hann", z, True),
                    ("düz Hann", ts_sm, True), ("düz ham", ts_sm, False)]:
    GR[isim] = I_omega(t, om, taper=tp)
    print(f"  [{isim}] bitti", flush=True)

def mask_online(k):
    m = np.zeros(len(om), bool)
    for l in LQ:
        m |= np.abs(om - l) < k * step
    return m
ON3, ON12 = mask_online(3), mask_online(12)
print(f"[102c] çizgi maskesi: ±3 adım {ON3.sum()} nokta, "
      f"±12 adım {ON12.sum()} nokta", flush=True)

W = 0.55
bins = np.arange(0.10, 4.4, W)
print("\n=== C2) BİN TABLOSU (F̄ = ⟨I⟩; rampa = α = ω/L) ===")
print(f"{'bin':>13} {'α':>7} {'rampa':>7} | {'ORTA ham':>9} {'ORTA Hann':>10} "
      f"{'öngörü':>8} | {'SIFIR ham':>10} {'SIFIR Hann':>11} {'öngörü':>8} | "
      f"{'çizgi-dışı(Hann,±12)':>21}")
ROWS = []
for lo in bins:
    m = (om >= lo) & (om < lo + W)
    if m.sum() < 50:
        continue
    a_mid = (lo + W / 2) / L
    qin = (LQ >= lo) & (LQ < lo + W)
    dw = np.exp(-LQ[qin]**2 * sig_t**2 / 2)
    A_mid = LAM[qin] / (L * np.sqrt(np.exp(LQ[qin]))) * \
        np.cos(np.pi * LQ[qin] / L) * dw
    A_zer = LAM[qin] / (L * np.sqrt(np.exp(LQ[qin]))) * dw
    P_mid = float((A_mid**2).sum()) * L / W
    P_zer = float((A_zer**2).sum()) * L / W
    row = dict(lo=lo, alpha=a_mid, nq=int(qin.sum()),
               om=float(GR["orta ham"][m].mean()),
               ot=float(GR["orta Hann"][m].mean()),
               pm=P_mid,
               zm=float(GR["sıfır ham"][m].mean()),
               zt=float(GR["sıfır Hann"][m].mean()),
               pz=P_zer,
               ooff=float(GR["orta Hann"][m & ~ON12].mean()),
               zoff=float(GR["sıfır Hann"][m & ~ON12].mean()),
               ooff3=float(GR["orta ham"][m & ~ON3].mean()),
               doff=float(GR["düz Hann"][m & ~ON12].mean()),
               dham=float(GR["düz ham"][m & ~ON3].mean()))
    ROWS.append(row)
    print(f"[{lo:.2f},{lo+W:.2f}) {a_mid:>7.3f} {a_mid:>7.3f} | "
          f"{row['om']:>9.4f} {row['ot']:>10.4f} {P_mid:>8.4f} | "
          f"{row['zm']:>10.4f} {row['zt']:>11.4f} {P_zer:>8.4f} | "
          f"orta {row['ooff']:>9.2e} sıfır {row['zoff']:>9.2e}")

# --------------- C3: "18×" denetimi ---------------
r0 = ROWS[0]
print(f"\n=== C3) '18× KARANLIK' DENETİMİ — ilk çizgi öncesi bin "
      f"[{r0['lo']:.2f},{r0['lo']+W:.2f}), rampa {r0['alpha']:.4f} ===")
print(f"  (bu bin log2={np.log(2):.4f}'nin ALTINDA: açık formülde ATOM YOK)")
for et, v in [("88'in ölçümü (orta, tapersiz)", r0["om"]),
              ("orta noktalar, Hann", r0["ot"]),
              ("SIFIRLAR, tapersiz", r0["zm"]),
              ("SIFIRLAR, Hann", r0["zt"]),
              ("RvM pürüzsüz, tapersiz (saf sızıntı)",
               float(GR["düz ham"][(om >= r0['lo']) & (om < r0['lo']+W)].mean())),
              ("RvM pürüzsüz, Hann (saf sızıntı)",
               float(GR["düz Hann"][(om >= r0['lo']) & (om < r0['lo']+W)].mean()))]:
    print(f"    {et:>38}: F̄={v:.4e}  → rampanın {r0['alpha']/v:.3e} KATI ALTINDA")

# --------------- C4: log2-log3 bini ---------------
r1 = ROWS[1]
lo = r1["lo"]; m = (om >= lo) & (om < lo + W)
print(f"\n=== C4) log2-log3 BİNİ [{lo:.2f},{lo+W:.2f}) — RAMPA DÖNÜŞÜ ===")
print(f"  bindeki çizgiler: {[q for q,_ in PQ if lo <= np.log(q) < lo+W]}"
      f"  rampa α={r1['alpha']:.4f}")
print(f"  88 dedi: F̄(hepsi)=0.102 vs rampa 0.089")
for et, key in [("orta, tapersiz", "orta ham"), ("orta, Hann", "orta Hann"),
                ("sıfır, tapersiz", "sıfır ham"), ("sıfır, Hann", "sıfır Hann")]:
    tot = float(GR[key][m].mean())
    on = float(GR[key][m & ON12].sum()) / m.sum()
    off = float(GR[key][m & ~ON12].sum()) / m.sum()
    print(f"    {et:>16}: F̄={tot:.4f}  (çizgi-üstü payı {on:.4f} = "
          f"{100*on/tot:.2f}%, çizgi-dışı payı {off:.3e} = {100*off/tot:.3f}%)")
print(f"  MUTLAK YASA ÖNGÖRÜSÜ (kanal ölçümü yok): orta {r1['pm']:.4f}  "
      f"sıfır {r1['pz']:.4f}   [rampa {r1['alpha']:.4f}]")
print(f"  KANAT SORUSU — bu bindeki ÇİZGİ-DIŞI düzey nedir?")
print(f"    orta noktalar (Hann)  : {r1['ooff']:.3e}  "
      f"= rampanın {r1['alpha']/r1['ooff']:.2e} katı altı")
print(f"    sıfırlar     (Hann)  : {r1['zoff']:.3e}  "
      f"= rampanın {r1['alpha']/r1['zoff']:.2e} katı altı")
print(f"    RvM pürüzsüz (Hann)  : {r1['doff']:.3e}  (saf pencere kanadı)")
print(f"    88'in tapersiz çizgi-dışı ölçümü (orta, ±3 adım maskesi): "
      f"{r1['ooff3']:.4f}; aynı maskeyle saf sızıntı {r1['dham']:.4f}")

# --------------- figür ---------------
R = ROWS
al = np.array([r["alpha"] for r in R])
fig, ax = plt.subplots(1, 2, figsize=(14, 5.0))
ax[0].plot(al, [r["om"] for r in R], "o-", c="firebrick", label="orta, tapersiz (88)")
ax[0].plot(al, [r["ot"] for r in R], "^-", c="darkorange", label="orta, Hann")
ax[0].plot(al, [r["zt"] for r in R], "s-", c="seagreen", label="sıfırlar, Hann")
ax[0].plot(al, [r["pz"] for r in R], "k--", lw=1.1, label="mutlak yasa (sıfır)")
ax[0].plot(al, al, "k:", lw=1.2, label="GUE rampası F=α")
ax[0].set_xlabel(r"$\alpha=\omega/L$"); ax[0].set_ylabel(r"$\bar F$")
ax[0].set_title("Berry ayrışımı: bin ortalaması (çizgiler dahil)")
ax[0].legend(fontsize=8); ax[0].grid(alpha=0.25)
ax[1].semilogy(al, [max(r["ooff"], 1e-22) for r in R], "^-", c="darkorange",
               label="çizgi-dışı: orta noktalar (Hann)")
ax[1].semilogy(al, [max(r["zoff"], 1e-22) for r in R], "s-", c="seagreen",
               label="çizgi-dışı: SIFIRLAR (Hann)")
ax[1].semilogy(al, [max(r["ooff3"], 1e-22) for r in R], "o-", c="firebrick",
               label="çizgi-dışı: 88'in tapersiz ölçümü")
ax[1].semilogy(al, [max(r["dham"], 1e-22) for r in R], ":", c="0.5",
               label="saf sızıntı (pürüzsüz, tapersiz)")
ax[1].semilogy(al, al, "k:", lw=1.2, label="rampa α")
ax[1].set_xlabel(r"$\alpha$"); ax[1].set_ylabel(r"$\bar F$ (çizgi-dışı)")
ax[1].set_title("Çizgi-arası süreklilik: '10-20 kat' mı, 10³–10¹⁴ kat mı?")
ax[1].legend(fontsize=8); ax[1].grid(alpha=0.25, which="both")
fig.tight_layout()
fig.savefig(HERE / "102c_berry.png", dpi=125)
np.savez(HERE / "102c_berry.npz", om=om, L=L, n=n, sig_t=sig_t,
         **{k.replace(" ", "_"): v for k, v in GR.items()},
         alpha=al, rows=np.array([[r["lo"], r["alpha"], r["om"], r["ot"],
                                   r["pm"], r["zm"], r["zt"], r["pz"],
                                   r["ooff"], r["zoff"], r["ooff3"],
                                   r["doff"], r["dham"]] for r in R]))
print("\n102c_berry.png + 102c_berry.npz yazıldı.", flush=True)
