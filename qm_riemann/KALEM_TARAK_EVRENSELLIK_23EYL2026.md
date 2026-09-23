# KALEM — 190: TARAĞIN EVRENSELLİĞİ — örneklem-dışı sınav (düşük pencere)
(23 Eylül 2026 — 188'in tarak resmi için ilk gerçek öngörü sınavı)

## Durum

188: pencere-ötesi iptal bant-bağımsız bir τ'-tarağı; tepeler Bragg (τ'≈1)
ve asal uydularında; (1.05,1.10] tepesinin %96'sı 1+τ_2, 1+τ_3'te; τ'≈1.15
tepesi L+log 6 gibi davranıyor; tepeler blok-blok yerel L ile kayıyor. H-188a
kaptanın |Ĝ| tasarım hatasıyla öldü; İMZALI öngörücü S_Re ikincil kayıtta
|corr| = 0.948. Bu kalem resmi YENİ bir pencerede, öngörüyle sınar.

## Kalem cebiri — neden Δω

Uydu koşulu ω' = L_yerel ± log n (Bragg n=1; asal uydular n=p; ikinci-mertebe
n=6=2·3). τ' = ω'/L koordinatında konumlar L'ye bağlıdır; ama
**Δω := ω' − L_yerel** koordinatında tepeler YÜKSEKLİKTEN BAĞIMSIZ olmalı:
Δω = 0, ±log 2 (0.693), ±log 3 (1.099), +log 6 (1.792)… Son pencere (L=12.03)
bunu dolaylı gösterdi (188f: tepeler L_b + sabit). Düşük pencere (L_ort = 10.484,
L ∈ [10.008, 10.836]) örneklem-dışı sınavdır.

**Rakip okuma (τ'-sabit):** tepeler son penceredeki τ' konumlarında kalır ⇒
düşük pencerede +log 2 tepesi Δω = 0.0576·10.484 = 0.604'e düşer (0.693 değil;
fark 0.089), +log 3 tepesi 0.957'ye (1.099 değil), 1.15 tepesi 1.562'ye
(1.792 değil). Ayrım blok-içi Δω çözünürlüğünde rahat.

**Veri-öncesi ayırt edilebilirlik:** pencere-içi L kayması ±0.0395 (τ') →
TÜM-PENCERE τ' dilimlerinde tepeler ±0.04 bulanır (konum sınavı yapılamaz);
bu yüzden konum sınavı 8 BLOKTA, her blok kendi L_b'siyle Δω ekseninde yapılır
(blok-içi L yayılımı ~0.10 → Δω bulanıklığı ~0.1·τ' ≲ 0.11 ω-birimi; tepeler
arası uzaklık ≥ 0.4 — ayrım var; son pencerede blok-içi yayılım çok daha dar).

## Hipotezler (ön-kayıtla donar; ölümler kurtarmasız)

- **H-190a — Δω EVRENSELLİĞİ (birincil):** düşük pencerede, 8 bloğun HAVUZ
  çekirdek profili Δω ekseninde (ω-dilim genişliği 0.025) tepelerini
  Δω = 0, +0.693, +1.099, +1.792'de verir: her biri için blok-medyan tepe konumu
  hedefin ±0.05 içinde VE τ'-sabit rakip konumdan (0.604 / 0.957 / 1.562) daha
  yakın. ÖLÜM: +log 2 tepesi (en güçlü, son'da %76) ±0.05 dışında ya da rakibe
  daha yakın.
- **H-190b — İMZALI ÖNGÖRÜCÜ (artık birincil):** tüm-pencere τ' dilimlerinde
  (188 ızgarası AYNEN: (0.86,1.30], 0.005) corr(Re v, S_Re) ≤ −0.80
  (S_Re = Σ_{q'∈s} a'·πτ'·cos(πτ')·Re Ĝ_mid(ω_{q'}), MUTLAK m_n). ÖLÜM:
  |corr| < 0.60.
- **H-190c — YAPISIZ SIFIR:** |corr(Re v, W)| ≥ |corr(Re v, S_Re)| ise tarak
  resmi ölür (W = Σ a').
- **KAYIT:** düşük pencerenin ζ_g (τ_c = 0.86; modül + açı; 8 bant) — iptalin
  ~1/3 büyüklüğü ve 180°'si daha düşük yükseklikte de sürüyor mu (eşiksiz).
  Son pencere Δω profili (188 dilim dosyalarından, blok-blok) yan yana.

## Kapılar

- **K0 — ZİNCİR (gerçek deniz, düşük pencere):** 155 eta_onbellek('dusuk'),
  184b (kos(etiket, 'dusuk')), 185b, 186b (ρ≡1), 187c — yeni etiketle,
  mevcut 'son' dosyalarının ÜZERİNE YAZMADAN (sarmalayıcı; içerik değişikliği
  yok). Kontrol: aynı sarmalayıcı 'son' için 187 ζ_g = 0.3287∠179.96°'yi
  vermeli (makine mührü).
- **K0b — ÖN-KAYIT (sha+damga, haritadan ÖNCE):** hedef ve rakip Δω konumları,
  ω-dilim ızgarası, blok tanımı (8 eşit blok, L_b = blok ortalaması), eşikler.
- **K1 — HARİTA:** 188b makinesi AYNEN (seri bit-bit mührü), düşük pencere,
  τ'≤1.30 (~66k çizgi). Ayrıca blok-blok K(Δω) (ω-dilim 0.025).
- **K2 — HÜKÜM + figür** (sol: Δω profili düşük vs son, hedef/rakip çizgileri;
  sağ: Re v ile S_Re).

Ölümler kurtarmasız; TEK DALGA; nohup + ≤30 sn yoklama; ajan git'e DOKUNMAZ;
sonuç ORTAK TEFTİŞE. Veri: qm_riemann/scratchpad/ (köprülü).
