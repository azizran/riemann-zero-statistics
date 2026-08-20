"""
101i — VURUŞ (MOIRÉ) TESTİ: DONUK BÖLGEDEKİ ARTIK D NEREDEN GELİYOR? (20 Ağu)
==========================================================================
101h, eşik-altı "yumuşak öncü"nün bir kısmının band kaçağı olduğunu
gösterdi; ama kaçaksız bir ARTIK kaldı (τ=0.050, ±0.005L: ζ .092,
χ₃ .092, χ₅ .065, χ₇ .155; vekil .009-.013). Hipotez: bu artık, üst
modların VURUŞU — iki asal-kuvvet çizgisinin fark frekansındaki moiré.

VURUŞ ENVANTERİ  (u_q = Λ(q)/(√q·log q); genlik çarpımı u_q1·u_q2)
  (q1,q2)   ω=log(q1/q2)   u·u      β'da (χ₄(2)=0 → 2,4,8,16 ÖLÜ)
  (3,2)      0.4055        0.408    ÖLÜ      ← en keskin ayrıştırıcı
  (7,5)      0.3365        0.169    canlı
  (4,3)      0.2877        0.144    ÖLÜ
  (5,4)      0.2231        0.112    ÖLÜ
  (5,3)      0.5108        0.258    canlı
  (7,3)      0.8473        0.218    canlı (β'nın donuk bölgesi log3'e
                                     dek uzanır; ζ'da orası çözülmüş)

YÖNTEM — SABİT-ω ÖLÇÜMÜ (τ değil): her pencere kendi τ_w = ω/L_w'siyle
katkı verir, kinematik c-faktörü PENCERE BAŞINA hesaplanır
(c_w = 2 sin(πτ_w)/(2πτ_w)), sonra pencereler-arası havuz:
   D(ω) = |Σ_w Σ_ω' (G/c_w)·r*| / Σ_w Σ_ω' |r|²
Vekil taban aynı ağırlıklarla, 40 permütasyon. Band ±0.02 (ω BİRİMİNDE,
±0.005L değil — vuruşlar sık, dar tutuldu). Çizgi-dışı filtre (>0.01)
korundu. Makine 101d/101h'den kopya; yalnız τ→ω çevrimi eklendi.

VERİ DİSİPLİNİ: ζ = 41_bigT_windows ilk 6 pencere; β = 101f sertifikalı
sıfırlar; her ikisinde de 101d düzlük segmentasyonu.

ÖN-MÜHÜR (ölçümden ÖNCE yazıldı):
  M1  ζ'da ω=0.4055'te D, 0.360 ve 0.455 kontrollerinden BELİRGİN
      yüksek.
  M2  β'da ω=0.4055'te HİÇBİR ÖZELLİK YOK (vuruş ölü) — en keskin
      öngörü; β'nın 0.4055'i kendi kontrollerinden farksız olmalı.
  M3  ω=0.5108 (5,3) ve 0.3365 (7,5)'te İKİ AİLEDE DE zayıf yükselti.
  M4  β'da ω=0.8473 (7,3)'te kontrollerine (0.78, 0.92) göre zayıf
      yükselti.
EK (ön-mühür sonrası eklendi, açıkça etiketli): (a) ölü vuruşlar
  0.2877 (4,3) ve 0.2231 (5,4) + kontrol 0.25 — M2 ile aynı mantıkta
  iki ek ayrıştırıcı; (b) YOĞUN ω-TARAMASI (ζ ve β, 0.30-0.56 adım
  0.01; β ayrıca 0.75-0.95), çünkü üç noktanın karşılaştırması değil
  TEPE PROFİLİ kanıttır; (c) ÇÖZÜNÜRLÜK KAPISI (aşağıda) — null'un
  anlamlı olması için ölçümün gerçek bir çizgiyi görebildiği
  gösterilmeli.

==========================================================================
SONUÇ (20 Ağustos, koşu 180 s) — M1 RET, M2 boş-İSABET, M3 RET, M4 RET
VURUŞ (MOIRÉ) HİPOTEZİ REDDEDİLDİ. Yerine geçen: ÇİZGİ OMZU.
==========================================================================
NOKTA ÖLÇÜMLERİ (band ±0.02; τ_w = ω/⟨L⟩; hepsi donuk bölgede):
   ω       ζ:D   ζ:vek | β:D   β:vek   ne
  0.2231  0.004  0.001 | 0.007 0.004   vuruş (5,4), β ölü
  0.2500  0.006  0.001 | 0.008 0.004   KONTROL
  0.2877  0.007  0.001 | 0.010 0.005   vuruş (4,3), β ölü
  0.3365  0.010  0.002 | 0.012 0.005   vuruş (7,5), ikisinde canlı
  0.3600  0.009  0.002 | 0.014 0.006   KONTROL
  0.4055  0.015  0.002 | 0.015 0.007   vuruş (3,2) — EN GÜÇLÜ, β ölü
  0.4550  0.019  0.003 | 0.017 0.007   KONTROL
  0.5108  0.028  0.004 | 0.024 0.008   vuruş (5,3), ikisinde canlı
  0.7800  0.254  0.010 | 0.076 0.012   KONTROL (β)
  0.8473  0.190  0.009 | 0.123 0.018   vuruş (7,3), β canlı
  0.9200  0.205  0.010 | 0.184 0.024   KONTROL (β)
VURUŞ/KONTROL oranları: ζ 0.4055 → 1.05×; β 0.4055 → 0.99×;
  β 0.8473 → 0.95×; ζ 0.8473 → 0.83×. HİÇBİRİ YÜKSELTİ DEĞİL.

YOĞUN TARAMA (0.30-0.56, adım 0.01): D(ω) her iki ailede de PÜRÜZSÜZ
ve TEKDÜZE ARTAN — hiçbir vuruş frekansında yerel tümsek yok.
  ζ: 0.008 (ω=0.30) → 0.050 (ω=0.56), maks tarama UCUNDA (0.560)
  β: 0.010 (ω=0.30) → 0.028 (ω=0.56), maks yine ucunda (0.550)
  vuruşların z-skorları (tarama ortalamasına göre): ζ (7,5) −0.96,
  (3,2) −0.52, (5,3) +0.90; β (7,5) −1.04, (3,2) −0.40, (5,3) +1.02.
  z'ler vuruş GENLİĞİYLE değil yalnızca ω KONUMUYLA sıralanıyor —
  en güçlü vuruş (3,2), en zayıflardan (5,3)'ün altında kalıyor.

YOĞUN TARAMA (0.75-0.95): 0.8473'te tümsek yok. Bunun yerine
  ζ: 0.356 (0.75) → 0.187 (0.87) → 0.230 (0.95) — log2 çizgisinden
     uzaklaşırken DÜŞEN, log3'e yaklaşırken yeniden YÜKSELEN bir VADİ.
  β: 0.063 (0.75) → 0.226 (0.95) — tekdüze artan (β'nın log2'si ölü,
     ilk çizgisi log3; vadi yok, sadece log3'e yaklaşma).

ÇÖZÜNÜRLÜK KAPISI (band ±0.004, çizgi filtresi KAPALI) — null'un
anlamlı olduğunun kanıtı, ve kendi başına yeni bir sonuç:
  ζ, ω = log2 = 0.6931 (CANLI çizgi): D 0.353 (ω=0.655) → 1.005
     (çizgide) → 0.426 (ω=0.735). KONTRAST 1.93×, keskin tepe.
     → Ölçüm dar bir çizgiyi ±0.004 bantla RAHATLIKLA görüyor.
  β, AYNI ω (ÖLÜ çizgi, χ₄(2)=0): D 0.038-0.063 arası düz, tepe YOK.
     KONTRAST 0.97×. → Ölü benek D-estimatöründe de ÖLÜ. 101d'nin
     hükmü artık TEK BİR FREKANSTA doğrudan görülüyor (101d bunu
     τ-taramasıyla dolaylı göstermişti).
EŞLEŞTİRİLMİŞ ÇÖZÜNÜRLÜKTE VURUŞ TARAMASI (aynı ±0.004 bant):
  (3,2) 0.3930-0.4180: ζ 0.012→0.015 tekdüze, std 0.0012, tepe yok;
        β 0.014→0.017 tekdüze, tepe yok.
  KONTROL 0.4430-0.4680: ζ 0.017-0.020 — VURUŞ BÖLGESİNDEN YÜKSEK.
  (5,3) 0.4985-0.5235: ζ 0.027→0.033 tekdüze, tepe yok.
  Yani çizgiyi 1.93× kontrastla gören alet, vuruşlarda %1 tümsek bile
  görmüyor. Null gerçek bir null, duyarlılık eksiği değil.

HÜKÜMLER:
  M1 RET — ζ'da 0.4055 kontrollerinden yüksek DEĞİL (1.05×; dar
     bantta kontroller daha da yüksek çıkıyor).
  M2 boş-İSABET — β'da 0.4055'te hiçbir özellik yok; ama ζ'da da yok.
     Öngörü "doğrulandı" ama AYRIŞTIRICI DEĞİL: iki aile de düz.
  M3 RET — 0.5108 ve 0.3365'te yerel yükselti yok; oradaki D pürüzsüz
     eğilimin üstünde değil.
  M4 RET — β'da 0.8473 kontrollerinin ALTINDA (0.95×).

YERİNE GEÇEN RESİM: donuk bölgedeki artık D, VURUŞ moiré'si değil;
İLK SAĞ KALAN ÇİZGİNİN OMZU. D(ω) çizgiler arasında pürüzsüz bir vadi
çizer — çizgiden uzaklaştıkça düşer, sonraki çizgiye yaklaştıkça
yükselir. 101h'de "band kaçağı" olarak bulunan şeyin sürekli hâli bu.
Omuz biçiminin TDS profili olup olmadığı → 101j.
"""

import numpy as np
from pathlib import Path

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
rng = np.random.default_rng(1011)

PKS = [2,3,4,5,7,8,9,11,13,16,17,19,23,25,27,29,31,32,37,41,43,47,49,53,
       59,61,64,67,71,73,79,81,83,89,97,101,103,107,109,113,121,125,127,128]
LINES_ALL = [np.log(qq) for qq in PKS]

# ---- 101d makinesi (kopya) ----

def pencere_hazirla(z, qeff):
    gaps = np.diff(z)
    mids = 0.5 * (z[:-1] + z[1:])
    Lw = float(np.log(qeff * mids / TWO_PI).mean())
    ds = gaps * np.log(qeff * mids / TWO_PI) / TWO_PI - 1
    tt = (mids - mids.mean()) / (mids[-1] - mids[0])
    P = np.vstack([np.ones_like(tt), tt, tt**2, tt**3]).T
    ds = ds - P @ np.linalg.lstsq(P, ds, rcond=None)[0]
    return (z, mids, ds, Lw)

def duz_segmentler(zz, sayim, min_n=3000, win=80, esik=0.5, pad=160):
    d = np.arange(len(zz)) - (sayim(zz) - sayim(zz[0]))
    from numpy.lib.stride_tricks import sliding_window_view
    med = np.median(sliding_window_view(d, win), axis=1)
    adim = med[win + 1:] - med[:-(win + 1)]
    bad = np.zeros(len(zz), bool)
    for i in np.where(np.abs(adim) > esik)[0] + win // 2:
        bad[max(0, i - pad):i + 2 * pad] = True
    seg, kes = [], 0
    s0 = 0
    for i in range(1, len(zz) + 1):
        if i == len(zz) or bad[i] != bad[i - 1]:
            if not bad[s0] and i - s0 >= min_n:
                seg.append(zz[s0:i])
            if i < len(zz) and bad[i]:
                kes += 1
            s0 = i
    return seg, kes

# ---- SABİT-ω varyantı (yeni; c-faktörü pencere başına) ----

def egri_omega(WIN, omegas, band=0.02, nfrek=160, cizgi_filtre=True):
    """Sabit-ω D(ω): her pencere kendi τ_w = ω/L_w ve c_w'siyle katkı."""
    Dv, Vv = [], []
    for om in omegas:
        num = 0.0 + 0j; den = 0.0
        Gs, rs = [], []
        for (zz, tm, ds, Lw) in WIN:
            tau_w = om / Lw
            kap = 2 * np.pi * tau_w
            c = 2 * np.sin(kap / 2) / kap
            oms = om + np.linspace(-band, band, nfrek)
            if cizgi_filtre:
                oms = np.array([o for o in oms
                                if min(abs(o - l) for l in LINES_ALL) > 0.01])
            for s0 in range(0, len(oms), 40):
                ob = oms[s0:s0 + 40]
                rr = np.exp(1j * np.outer(ob, zz)).sum(axis=1)
                GG = (np.exp(1j * np.outer(ob, tm)) * ds[None, :]).sum(axis=1)
                GG = GG / c                      # pencere başına c
                num += (GG * np.conj(rr)).sum()
                den += (np.abs(rr)**2).sum()
                Gs.append(GG); rs.append(rr)
        G_all, r_all = np.concatenate(Gs), np.concatenate(rs)
        fl = []
        for _ in range(40):
            perm = rng.permutation(len(r_all))
            fl.append(abs((G_all * np.conj(r_all[perm])).sum()) / den)
        Dv.append(abs(num) / den); Vv.append(float(np.mean(fl)))
    return Dv, Vv

# ------------------------- pencereler -------------------------

def rvm_sayim(t):
    x = t / TWO_PI
    return x * np.log(x / np.e)

d41 = np.load(HERE / "41_bigT_windows.npz")
K41 = sorted({x.split("_")[1] for x in d41.files}, key=lambda s: int(s[:-1]))
WINZ = []
for k in K41[:6]:
    gz, tm = d41[f"gaps_{k}"], d41[f"tmid_{k}"]
    zz = np.empty(len(gz) + 1)
    zz[0] = tm[0] - gz[0] / 2
    zz[1:] = zz[0] + np.cumsum(gz)
    seg, _ = duz_segmentler(zz, rvm_sayim)
    WINZ += [pencere_hazirla(s, 1.0) for s in seg]

exec(open(HERE / "98_L_motoru.py").read().split('CHI4 = ')[0])
HERE = Path(__file__).resolve().parent
CHI4 = {0: 0, 1: 1, 2: 0, 3: -1}
zb = np.load(HERE / "101f_beta_zeros.npz")["zeros"]
lo = np.exp(np.log(zb[0] + 1) + 0.25 * (np.log(zb[-1]) - np.log(zb[0] + 1)))
zb = zb[zb >= lo]
Mb = Lmotor(4, CHI4, 1)
segb, _ = duz_segmentler(zb, lambda t: Mb.theta(t) / np.pi)
kenar = np.exp(np.linspace(np.log(zb[0]), np.log(zb[-1] * 1.0001), 4))
WINB = []
for s in segb:
    for i in range(3):
        p = s[(s >= kenar[i]) & (s < kenar[i + 1])]
        if len(p) >= 3000:
            WINB.append(pencere_hazirla(p, 4.0))

PEN = {"zeta": WINZ, "beta": WINB}
LORT = {}
for a in PEN:
    Ls = np.array([w[3] for w in PEN[a]])
    ns = np.array([len(w[0]) for w in PEN[a]], float)
    LORT[a] = float(np.average(Ls, weights=ns))
    print(f"[{a}] {len(PEN[a])} pencere  ⟨L⟩_n={LORT[a]:.3f}  "
          f"L∈[{Ls.min():.2f},{Ls.max():.2f}]", flush=True)
TILK = {"zeta": np.log(2) / LORT["zeta"], "beta": np.log(3) / LORT["beta"]}
print(f"  τ_ilk: ζ {TILK['zeta']:.5f} (log2)  β {TILK['beta']:.5f} (log3)",
      flush=True)

# ------------------------- nokta ölçümleri -------------------------
NOKTA = [(0.2231, "vuruş (5,4)  u·u=0.112  β:ÖLÜ   [EK]"),
         (0.2500, "KONTROL                          [EK]"),
         (0.2877, "vuruş (4,3)  u·u=0.144  β:ÖLÜ   [EK]"),
         (0.3365, "vuruş (7,5)  u·u=0.169  β:canlı"),
         (0.3600, "KONTROL"),
         (0.4055, "vuruş (3,2)  u·u=0.408  β:ÖLÜ"),
         (0.4550, "KONTROL"),
         (0.5108, "vuruş (5,3)  u·u=0.258  β:canlı"),
         (0.7800, "KONTROL (β)"),
         (0.8473, "vuruş (7,3)  u·u=0.218  β:canlı"),
         (0.9200, "KONTROL (β)")]
OMS = [o for o, _ in NOKTA]

RES = {}
for a in PEN:
    RES[a] = egri_omega(PEN[a], OMS)
    print(f"  ...{a} nokta ölçümleri bitti", flush=True)

print("\n=== SABİT-ω NOKTA ÖLÇÜMLERİ (band ±0.02, 40-perm vekil) ===")
print(f"{'ω':>7} {'ζ: D':>7} {'vekil':>7} {'τ_w':>13} | "
      f"{'β: D':>7} {'vekil':>7} {'τ_w':>13}  açıklama")
for i, (om, aciklama) in enumerate(NOKTA):
    tz = om / LORT["zeta"]; tb = om / LORT["beta"]
    dz = "ÇÖZÜLMÜŞ" if tz > TILK["zeta"] else "donuk"
    db = "ÇÖZÜLMÜŞ" if tb > TILK["beta"] else "donuk"
    print(f"{om:>7.4f} {RES['zeta'][0][i]:>7.3f} {RES['zeta'][1][i]:>7.3f} "
          f"{tz:>6.4f}/{dz:<7}| "
          f"{RES['beta'][0][i]:>7.3f} {RES['beta'][1][i]:>7.3f} "
          f"{tb:>6.4f}/{db:<7} {aciklama}")

print("\n=== VURUŞ / KONTROL ORANLARI ===")
def dg(a, om):
    return RES[a][0][OMS.index(om)]
print(f"  ζ  0.4055 / ort(0.360, 0.455) = {dg('zeta',0.4055):.3f} / "
      f"{(dg('zeta',0.36)+dg('zeta',0.455))/2:.3f} = "
      f"{dg('zeta',0.4055)/((dg('zeta',0.36)+dg('zeta',0.455))/2):.2f}")
print(f"  β  0.4055 / ort(0.360, 0.455) = {dg('beta',0.4055):.3f} / "
      f"{(dg('beta',0.36)+dg('beta',0.455))/2:.3f} = "
      f"{dg('beta',0.4055)/((dg('beta',0.36)+dg('beta',0.455))/2):.2f}")
print(f"  β  0.8473 / ort(0.780, 0.920) = {dg('beta',0.8473):.3f} / "
      f"{(dg('beta',0.78)+dg('beta',0.92))/2:.3f} = "
      f"{dg('beta',0.8473)/((dg('beta',0.78)+dg('beta',0.92))/2):.2f}")
print(f"  ζ  0.8473 / ort(0.780, 0.920) = {dg('zeta',0.8473):.3f} / "
      f"{(dg('zeta',0.78)+dg('zeta',0.92))/2:.3f} = "
      f"{dg('zeta',0.8473)/((dg('zeta',0.78)+dg('zeta',0.92))/2):.2f}")

# ------------------------- yoğun tarama [EK] -------------------------
SCAN1 = np.round(np.arange(0.30, 0.5601, 0.01), 4)
SCAN2 = np.round(np.arange(0.75, 0.9501, 0.01), 4)
TAR = {}
for a in PEN:
    TAR[(a, 1)] = egri_omega(PEN[a], list(SCAN1))
    print(f"  ...{a} tarama-1 bitti", flush=True)
for a in PEN:
    TAR[(a, 2)] = egri_omega(PEN[a], list(SCAN2))
    print(f"  ...{a} tarama-2 bitti", flush=True)

VUR1 = {0.2231: "(5,4)", 0.2877: "(4,3)", 0.3365: "(7,5)",
        0.4055: "(3,2)", 0.5108: "(5,3)"}
VUR2 = {0.8473: "(7,3)", 0.8109: "(9,4)?", 0.7885: "(11,5)"}

print("\n=== YOĞUN ω-TARAMASI [EK]: 0.30 - 0.56 (adım 0.01) ===")
print(f"{'ω':>7} {'ζ:D':>7} {'ζ:vek':>7} {'β:D':>7} {'β:vek':>7}  vuruş")
for i, om in enumerate(SCAN1):
    et = "".join(f" ←{v} ω={k}" for k, v in VUR1.items() if abs(k - om) < 0.006)
    print(f"{om:>7.3f} {TAR[('zeta',1)][0][i]:>7.3f} "
          f"{TAR[('zeta',1)][1][i]:>7.3f} {TAR[('beta',1)][0][i]:>7.3f} "
          f"{TAR[('beta',1)][1][i]:>7.3f} {et}")

print("\n=== YOĞUN ω-TARAMASI [EK]: 0.75 - 0.95 (adım 0.01) ===")
print(f"{'ω':>7} {'ζ:D':>7} {'ζ:vek':>7} {'β:D':>7} {'β:vek':>7}  vuruş")
for i, om in enumerate(SCAN2):
    et = "".join(f" ←{v} ω={k}" for k, v in VUR2.items() if abs(k - om) < 0.006)
    print(f"{om:>7.3f} {TAR[('zeta',2)][0][i]:>7.3f} "
          f"{TAR[('zeta',2)][1][i]:>7.3f} {TAR[('beta',2)][0][i]:>7.3f} "
          f"{TAR[('beta',2)][1][i]:>7.3f} {et}")

print("\n=== TARAMA-1 TEPE İSTATİSTİĞİ (0.30-0.56) ===")
for a in PEN:
    d = np.array(TAR[(a, 1)][0])
    print(f"  {a}: ort {d.mean():.3f}  std {d.std():.3f}  "
          f"maks {d.max():.3f} @ ω={SCAN1[d.argmax()]:.3f}  "
          f"min {d.min():.3f} @ ω={SCAN1[d.argmin()]:.3f}")
    for om, et in VUR1.items():
        if SCAN1[0] <= om <= SCAN1[-1]:
            v = np.interp(om, SCAN1, d)
            print(f"      {et} ω={om:.4f}: D={v:.3f}  "
                  f"z=(D−ort)/std = {(v-d.mean())/d.std():+.2f}")

# ---- ÇÖZÜNÜRLÜK KAPISI (ön-mühür sonrası; NULL'un anlamlı olması için) ----
# Vuruş bulunamadı. Ama ölçüm DAR bir çizgiyi görebiliyor mu? Band ±0.02'de
# 160 frekans üzerine ortalama alınıyor; keskin bir mod seyrelmiş olabilir.
# POZİTİF KONTROL: ζ için ω=log2 GERÇEK ve CANLI bir çizgi. Aynı makine
# dar bantla (±0.004, çizgi filtresi KAPALI) oradan geçirilir: tepe
# görüyorsak duyarlılık var, null anlamlıdır. β için aynı ω ÖLÜ çizgidir
# (χ₄(2)=0) — yerleşik negatif kontrol.
BK = 0.004
LINE_SCAN = np.round(np.arange(0.655, 0.73501, 0.005), 4)
LKON = {}
for a in PEN:
    LKON[a] = egri_omega(PEN[a], list(LINE_SCAN), band=BK, cizgi_filtre=False)
    print(f"  ...{a} çizgi-kontrol taraması bitti", flush=True)

print(f"\n=== ÇÖZÜNÜRLÜK KAPISI [EK]: log2 ÇİZGİSİNDEN GEÇİŞ "
      f"(band ±{BK}, çizgi filtresi KAPALI) ===")
print("  ζ: log2 CANLI çizgi (pozitif kontrol) | β: log2 ÖLÜ (χ₄(2)=0)")
print(f"{'ω':>8} {'ζ:D':>8} {'ζ:vek':>7} {'β:D':>8} {'β:vek':>7}")
for i, om in enumerate(LINE_SCAN):
    et = "  ←← log2 = 0.6931" if abs(om - np.log(2)) < 0.003 else ""
    print(f"{om:>8.4f} {LKON['zeta'][0][i]:>8.3f} {LKON['zeta'][1][i]:>7.3f} "
          f"{LKON['beta'][0][i]:>8.3f} {LKON['beta'][1][i]:>7.3f}{et}")
for a in PEN:
    d = np.array(LKON[a][0])
    ic = np.abs(LINE_SCAN - np.log(2)) <= 0.006
    dis = np.abs(LINE_SCAN - np.log(2)) >= 0.020
    print(f"  {a}: çizgide {d[ic].mean():.3f}  uzakta {d[dis].mean():.3f}  "
          f"KONTRAST = {d[ic].mean()/max(d[dis].mean(),1e-9):.2f}×")

# ---- EŞLEŞTİRİLMİŞ ÇÖZÜNÜRLÜKTE VURUŞ TARAMASI ----
VS = {"(3,2) 0.4055": np.round(np.arange(0.3930, 0.41801, 0.0025), 4),
      "KONTROL 0.455": np.round(np.arange(0.4430, 0.46801, 0.0025), 4),
      "(5,3) 0.5108": np.round(np.arange(0.4985, 0.52351, 0.0025), 4)}
VR = {}
for ad, gg in VS.items():
    for a in PEN:
        VR[(ad, a)] = egri_omega(PEN[a], list(gg), band=BK)
    print(f"  ...dar-band vuruş taraması '{ad}' bitti", flush=True)

print(f"\n=== DAR-BAND VURUŞ TARAMASI [EK] (band ±{BK}, adım 0.0025) ===")
for ad, gg in VS.items():
    print(f"  --- {ad} ---")
    print(f"{'ω':>8} {'ζ:D':>8} {'ζ:vek':>7} {'β:D':>8} {'β:vek':>7}")
    for i, om in enumerate(gg):
        print(f"{om:>8.4f} {VR[(ad,'zeta')][0][i]:>8.3f} "
              f"{VR[(ad,'zeta')][1][i]:>7.3f} {VR[(ad,'beta')][0][i]:>8.3f} "
              f"{VR[(ad,'beta')][1][i]:>7.3f}")
    for a in PEN:
        d = np.array(VR[(ad, a)][0])
        print(f"    {a}: ort {d.mean():.4f} std {d.std():.4f} "
              f"maks {d.max():.4f} @ {gg[d.argmax()]:.4f}")

np.savez(HERE / "101i_vurus.npz",
         nokta_om=np.array(OMS),
         line_scan=LINE_SCAN,
         **{f"line_{a}_{k}": np.array(LKON[a][k]) for a in PEN for k in (0, 1)},
         **{f"vs_{ad.split()[0]}_{a}_{k}": np.array(VR[(ad, a)][k])
            for ad in VS for a in PEN for k in (0, 1)},
         **{f"vsgrid_{ad.split()[0]}": gg for ad, gg in VS.items()},
         **{f"nokta_{a}_{k}": np.array(RES[a][k])
            for a in PEN for k in (0, 1)},
         scan1=SCAN1, scan2=SCAN2,
         **{f"scan{j}_{a}_{k}": np.array(TAR[(a, j)][k])
            for a in PEN for j in (1, 2) for k in (0, 1)})
print("\n101i_vurus.npz yazıldı.")
