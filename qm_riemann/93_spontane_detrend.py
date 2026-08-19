"""
93 — SPONTANE YANITIN DÜŞÜK-τ SORGUSU: DETREND + VEKİL TABAN (21 Ağustos)
==========================================================================
92-T4'ün τ=0.1 noktası (D_spont = 0.604) şüpheliydi. Üç kontrol:
  (a) YEREL katlama: ds = g·L(t)/2π − 1 (global ort. yerine) + kübik
      polinom detrend — yoğunluk-sürüklenmesi tamamen dışarı.
  (b) VEKİL TABAN: G(ω_i) × ρ(ω_π(i)) karıştırılmış eşleme — estimatörün
      |Σ|-pozitif-yanlılık tabanı.
  (c) İNCE τ-IZGARASI: ilk-çizgi eşiği τ₂ = log2/L ≈ 0.067 civarında
      yapı var mı — Berry kuazi-rastgele bölgesi (form faktörü rampanın
      18 kat altında) gerçek uzun-dalga katılığı verebilir.
Karar çerçevesi: detrend sonrası D→~1 ise artefakttı; düşük kalır ve
taban temizse GERÇEK: zeta'nın uzun-dalga spontane modları CUE'dan
katı (aritmetik-yoksul bölgenin rijitliği).

SONUÇ (21 Ağustos, gece — kapanış):
  1. ARTEFAKT DEĞİL: detrend hiçbir şeyi değiştirmedi (D_detrend =
     D_92tarzı üç ondalıkta, her τ'da); vekil taban 0.01-0.05 (temiz).
  2. YENİ YAPI — DONMA→CUE GEÇİŞİ ÇÖZÜLDÜ: D_spont(τ) = 0.041 (τ=0.04!)
     → 0.240 (0.06) → 0.492 (0.08) → 0.662 (0.10) → 0.788 (0.125) →
     0.925 (0.15) → 0.964 (0.20) → CUE'ya oturuyor (0.3+). İlk asal
     çizgisi eşiğinin (τ₂ = log2/L ≈ 0.067) ALTINDA gaz fiilen DONMUŞ:
     konum içeriği var, gap yanıtı ~SIFIR (katı kolektif öteleme —
     Berry'nin kuazi-rastgele, rampa-altı bölgesi). Çizgiler devreye
     girdikçe süreklilik CUE-termal davranışa geçiyor.
  3. İKİ-DAL YANIT RESMİ TAMAM: asal (koheran) dalı D 0.91→0.48 iner;
     spontane dal ~0→1.18 çıkar; İKİ EĞRİ τ≈0.15-0.2'DE KESİŞİR.
     Uzun dalgada asallar gap'leri spontane modlardan ÇOK daha etkin
     sıkıştırır (adyabatik denge-deformasyonu); kısa dalgada tersine
     perdelenir. "Perde neden asal korosuna iner" sorusu inceldi:
     sorunun yarısı "spontane modlar uzun dalgada neden fonon-gibi
     katı öteler" çıktı.
"""

import numpy as np
from pathlib import Path

HERE = Path(__file__).resolve().parent
TWO_PI = 2 * np.pi
rng = np.random.default_rng(93)

d41 = np.load(HERE / "41_bigT_windows.npz")
K41 = sorted({x.split("_")[1] for x in d41.files}, key=lambda s: int(s[:-1]))
lines_all = [np.log(q) for q in
    [2,3,4,5,7,8,9,11,13,16,17,19,23,25,27,29,31,32,37,41,43,47,49,53,59,
     61,64,67,71,73,79,81,83,89,97,101,103,107,109,113,121,125,127,128]]
CUE_REF = [(0.031, 1.001), (0.051, 1.003), (0.102, 1.010), (0.199, 1.038),
           (0.301, 1.081), (0.398, 1.132), (0.449, 1.163), (0.500, 1.183)]
cue_t, cue_d = np.array(CUE_REF).T

WIN = []
for k in K41[:6]:
    gz, tm = d41[f"gaps_{k}"], d41[f"tmid_{k}"]
    Lz = float(np.log(tm / TWO_PI).mean())
    zz = np.empty(len(gz) + 1)
    zz[0] = tm[0] - gz[0] / 2
    zz[1:] = zz[0] + np.cumsum(gz)
    ds_loc = gz * np.log(tm / TWO_PI) / TWO_PI - 1          # yerel katlama
    tt = (tm - tm.mean()) / (tm[-1] - tm[0])
    P = np.vstack([np.ones_like(tt), tt, tt**2, tt**3]).T
    ds_det = ds_loc - P @ np.linalg.lstsq(P, ds_loc, rcond=None)[0]
    ds_glob = gz / gz.mean() - 1                            # 92'nin tanımı
    WIN.append((zz, tm, ds_det, ds_glob, Lz))

print(f"{'τ0':>6} {'D_detrend':>9} {'D_92tarzı':>9} {'vekil':>6} {'CUE':>6} {'n_frek':>6}")
for tau0 in [0.04, 0.06, 0.08, 0.10, 0.125, 0.15, 0.20, 0.30, 0.40, 0.50]:
    num_d = num_g = 0.0 + 0j
    den = 0.0
    Gs_all, rhos_all = [], []
    ntot = 0
    for zz, tm, ds_det, ds_glob, Lz in WIN:
        oms = tau0 * Lz + np.linspace(-0.02 * Lz, 0.02 * Lz, 160)
        oms = np.array([o for o in oms
                        if min(abs(o - l) for l in lines_all) > 0.01])
        for s0 in range(0, len(oms), 40):
            ob = oms[s0:s0+40]
            rr = np.exp(1j * np.outer(ob, zz)).sum(axis=1)
            E = np.exp(1j * np.outer(ob, tm))
            Gd = (E * ds_det[None, :]).sum(axis=1)
            Gg = (E * ds_glob[None, :]).sum(axis=1)
            num_d += (Gd * np.conj(rr)).sum()
            num_g += (Gg * np.conj(rr)).sum()
            den += (np.abs(rr)**2).sum()
            Gs_all.append(Gd); rhos_all.append(rr)
        ntot += len(oms)
    kap = 2 * np.pi * tau0
    c = 2 * np.sin(kap / 2) / kap
    G_all = np.concatenate(Gs_all)
    r_all = np.concatenate(rhos_all)
    fl = []
    for _ in range(40):
        perm = rng.permutation(len(r_all))
        fl.append(abs((G_all * np.conj(r_all[perm])).sum()) / (c * den))
    cue = np.interp(tau0, cue_t, cue_d)
    print(f"{tau0:>6.3f} {abs(num_d)/(c*den):>9.3f} {abs(num_g)/(c*den):>9.3f} "
          f"{np.mean(fl):>6.3f} {cue:>6.3f} {ntot:>6}")
print("\n(τ₂ = log2/L ≈ 0.062-0.070 pencere aralığı — ilk asal çizgisi eşiği)")
