# KALEM — 197: AYNA YASASI — Bogomolny–Keating'in ayna-yanı uydularının kör profil sınavı
(25 Eylül 2026 — kullanıcı onaylı. TÜRETİM TEFTİŞİNDEN GEÇTİ (Sonnet): cebir p=2..11
için 1e-15 ile doğrulandı; işaret sözleşmesi görülmüş veriyle zaten ayrışmış; makine
sapmaları sağlam (32 blokta blok-içi L yayılımı 0.017-0.041, mix çizgi başına ⇒ havuz alt
kümesi geçerli, jk kenarları birebir). Teftişin KRİTİK bulgusu (ilk taslakta R-g'nin
öngördüğü medyan R = 0.68, TUTAR alt sınırı 0.70'in hemen altında ⇒ sınav ayırt etmiyordu)
ve iki UYARI (pencere sızıntısı; kalibrasyon şemaları farklı) bu sürümde giderildi:
profil düzeyi uyum + BK/R-g ayrım kuralı + güç kapısı.
MAKİNE RAPORU SONRASI (ölçümden önce, kör veri görülmeden): M1 bit-bit geçti (0.0 fark,
70 650 karşılaştırma); iki düzeltme: (1) taranan pencere-içi çizgiler karışım serisinde
zaten var ⇒ dilimin kendi katkısı K'ya pozitif sapma sokuyordu (kör hedeflerde 22-32/32 blok,
kalibrasyonda 0-4/32) ⇒ BİRİNCİL K = K_düz (dilim j'nin katkısı mix'ten çıkarılır; 190'ın
pencere-ötesi düzenine denk); (2) HAVUZ' sınırı 0.62 → 0.74 (104 → ~320 çizgi; kör bant
tüm bloklarda yine kapsanır).
İKİNCİ MAKİNE RAPORU (sentetik; kör veri görülmeden): kör pencerede serbest taban ile sıkışık
tarak (7 ana + 13-16 çeyrek tepe) dejenere ⇒ serbest q_f ile σ_Z ≈ 0.35 (M4 hep kalır);
çeyreksiz model BK-doğruyken R̄_Z'yi 0.73'e çeker (tehlikeli). KARAR: çeyrek aile ölçeği
KALİBRASYON penceresinden taşınır (orada tepeler seyrek, q_f iyi belirli) ve karar İKİ MODELİN
ÖZ-TUTARLILIĞI ile verilir (M_BK / M_Rg); M6 iki-doğrulu güç sınavına çevrildi. Havuz 341
çizgi, 24 yeni τ'-dilimi (703 çizgi), [−2.2, +2.4]'te tüm dilimlerde B_j = 32; M1 (K_ham ve
K_düz) 190 ile bit-bit.
ÜÇÜNCÜ MAKİNE RAPORU (M6 ön-deneme, 2×200 sentetik replika; kör veri görülmeden): M6 KALDI —
BK-doğruda doğru kazanma 0.715 < 0.80 (yanlış kazanma 0.005 / 0.000 — ayırıcılık sağlam);
nedenler: 8-grup jk σ gerçek saçılımı düşük tahmin ediyor; β penceresi σ_β ≈ 0.24 için fazla
dar. KARAR: (1) tüm karar σ'ları σ_eff = f·σ_jk, f (R̄_Z için f_Z, ρ için f_ρ) = M6 KALİBRASYON replikalarında (güç
replikalarından AYRI tohumlar) z = (tahmin − doğru)/σ_jk'nın standart sapması, iki doğrunun
büyüğü; f ÖN-KAYDA yazılır; (2) β karar koşulundan çıkarıldı (KAYIT).
DÖRDÜNCÜ MAKİNE RAPORU (M6, ayrı tohumlarla; kör veri görülmeden): f_Z = 1.565, f_ρ = 1.374
(K_KAL tohumları, 2×200); güç (K_GUC, 2×200): doğru kazanır 0.930 / 0.860, yanlış kazanır
0.000 / 0.005 ✓; (ii) tek başarısızlık R-g-doğruda KARARA GİRMEYEN kalibrasyon bileşeni
A(7/2) kapsaması 0.875 < 0.88 (binom sd ≈ 0.023; R-g altında kalibrasyon modeli bilerek
yanlış belirli). KARAR: (ii) yalnız KARARA GİREN niceliklerde, kararda kullanılan σ ile
sınanır (aşağıda M6).)

## Durum

192 kataloğu + 196 ŞERHİ: BK uydu katsayıları c(r) AYNA YANINDA tutuyor (−log2 0.96,
−log3 1.05, çıplak c ile, κ(+log2) normalizasyonu), saf asallar sönük, KARIŞIK oranlarda
ıskalıyor (+log(5/3) −9.6σ …). Nokta okumaları 8-blok çizgi-biçimiyle bulanık: blok içi
L yayılımı 0.07-0.15 ⇒ her uydu ~0.1 genişliğinde bir kutuya yayılır; komşular karışır.

SORU: BK'nın ayna-yanı öngörüsü — −log x zarfı TAM 1/x; tam sayı ve buçuklu aileler AYNI
zarf üstünde; asal kuvvetleri aynada bastırılmaz — hiç görüntülenmemiş/hesaplanmamış bir
Δω bandında, dar çizgi-biçimiyle ve profil düzeyinde tutuyor mu?

## Kalem cebiri (kaptan, veri-öncesi; teftişte doğrulandı)

1. BK (Bogomolny 2007, arXiv:0708.4223, eş. 7.5-7.6): R₂^off(ε) ∝ e^{iLε} f(ε) + c.c.,
   f(ε) = |ζ(1+iε)|² Φ_off(ε) = Π_p g_p(ε log p),
   g_p(θ) = |1 − e^{−iθ}/p|^{−2} · [1 − (1 − e^{iθ})²/(p−1)²].
2. g_p'nin Fourier katsayıları (e^{ikθ}): ilk çarpan A_k = p^{−|k|}/(1 − p^{−2});
   ikinci çarpan B_0 = 1 − 1/(p−1)², B_1 = 2/(p−1)², B_2 = −1/(p−1)² (e^{0, iθ, 2iθ}).
   ĝ_p(k) = Σ_j B_j A_{k−j}:
   - k ≤ 0: ĝ_p(k)/ĝ_p(0) = p^{k};
   - k = 1: ĝ_p(1)/ĝ_p(0) = p/(p−1)²;
   - k ≥ 2: ĝ_p(k) ∝ p^{−k}(B_0 + pB_1 + p²B_2) = p^{−k}[1 − (p−1)²/(p−1)²] = 0.
3. Uydu Δω = log r'nin göreli genliği c(r) = Π_p f_p(k_p), r = Π p^{k_p} (sade kesir);
   f_p(k<0) = p^{k}, f_p(1) = p/(p−1)², f_p(k≥2) = 0. İşaret sözleşmesi: f'nin +log r
   Bohr bileşeni ölçüm ekseninde Δω = +log r'ye düşer (görülmüş veride +log2 en parlak
   ↔ c(2) = 2; ters sözleşme orada bastırma öngörürdü).
4. AYNA YANI (r = 1/x, Δω = −log x), kesirler sade:
   - x = n tam sayı: c = 1/n (asal kuvvetleri DAHİL: c(1/8) = 1/8);
   - x = m/2, m tek: c(2/m) = 2/m = 1/x;
   - x = b/3 (3∤b): c = 0.75/b = 0.25/x; x = b/6 (gcd(b,6)=1): c = 1.5/b = 0.25/x;
     x = b/5 (5∤b): 0.0625/x; x = b/4 (b TEK): 0 (pay 4 ⇒ f_2(2) = 0).
   ⇒ Ayna yanında baskın tarak x ∈ ½ℕ'de, zarf TAM 1/x; aralarda çeyrek zarflı
   (0.25/x) küçük tepeler.
5. EŞLEME VARSAYIMI (sınanan): ölçülen parlaklık κ ∝ c(r) ÇIPLAK. RAKİP R-g:
   κ ∝ c(r)·g(ω'), g(ω') = e^{ω'/2}/ω' (dilimdeki çizgi sayısı × a_q' ∝ √q'/log q'),
   ω' = L + Δω, L = 10.48393.

## Kör bant ve görülmüşlük beyanı (teftişle düzeltildi)

- 190b haritası Δω ∈ [−1.775, +3.55] HESAPLANDI (bloğa göre). GÖRÜNTÜLENEN/ARANAN en alt
  Δω: 190d figürü ARALIK = (−1.3, 2.3); 192e (−1.35, 2.35); 190c/190e arama pencereleri
  ≥ −1.249; 192 kataloğu ≥ −1.099; log'larda Δω < −1.6 bandına ait κ yok (yalnız
  sayaçlar). Δω < −1.35 hiç görüntülenmedi; Δω < −1.775 HİÇ HESAPLANMADI. (Statik denetim
  kayıt dışı etkileşimli incelemeyi dışlayamaz; kaptanın bilgisi dahilinde yok.)
- SINAV BANDI: Δω ∈ [−2.12, −1.58] (x ∈ [4.86, 8.33]).
  Tam sayı hedefleri x = 5, 6, 7, 8 (Δω = −1.609, −1.792, −1.946, −2.079);
  buçuklu hedefler x = 5.5, 6.5, 7.5 (Δω = −1.705, −1.872, −2.015).

## Ölçüm tanımı (190b/188b makinesi; ÜÇ BEYANLI SAPMA, gerisi AYNEN)

- S1 HAVUZ' = pencere çizgileri τ ∈ [0.45, 0.74) (190'ın HAVUZ'unun alt kümesi; mix =
  188b.karisim AYNEN — çizgi başına hesaplanır, havuz filtresi sonradan).
- S1b K_düz (BİRİNCİL): her ω-dilimi j için mix, j dilimindeki çizgilerin kesik seriye
  katkısı ÇIKARILMIŞ seriden alınır (190'da τ' > 0.86 çizgileri seride hiç yoktu ⇒ orada
  K_düz ≡ K; M1 bunu sınar). K_ham (çıkarmasız) KAYIT.
- S2 Tarama: τ' ∈ [0.74, 0.86) için YENİ τ'-dilimleri (0.005 ızgarası,
  188b.dilimleri_kur AYNEN) + mevcut (0.86, 1.30] dilimleri. Δω blok-yerel:
  j = floor((log q' − L_b)/0.025 + 0.5).
- S3 Bloklar: 32 eşit bitişik blok (kenar = linspace(0, N, 33).astype(int); L_b = blok
  içi mean log(m_n/2π)). Jackknife: 8 grup × 4 ardışık blok (grup kenarları 190'ın 8-blok
  kenarlarıyla birebir); mix loo 188b AYNEN.
- Profil: κ_Σ(j) = Σ_{b ∈ B_j} κ_b(j)·n_b / Σ_{b∈B_j} n_b, B_j = j dilimini TAM kapsayan
  bloklar (192 AYNEN).
- ÇİZGİ BİÇİMİ (parametresiz, veriden): Δω = log r'deki bir uydu için
  h(j; r) = Σ_{b∈B_j} n_b·h_b(j; r) / Σ_{b∈B_j} n_b,
  h_b(j; r) = (1/n_b)·#{n ∈ b : floor((L_n − L_b + log r)/0.025 + 0.5) = j},
  L_n = log(m_n/2π). (Uydunun blok b'deki konumu, gap n'in yerel L'sine göre kayar.)
- PROFİL UYUMU (doğrusal en küçük kareler, jk-loo ile se):
  κ_Σ(j) = Σ_k A_k·h(j; r_k) + q_f·Σ_{çeyrek} (0.25/x_i)·h(j; r_i) + a + b·Δω_j.
  (Ayna yanında BK tepeleri İKİ zarf düzeyinde: ana 1/x (x ∈ ½ℕ) ve çeyrek 0.25/x
  (x = b/3, b/6); diğer aileler c < 0.03 — kapalı biçimle tarandı.)
  - KÖR uyum penceresi Δω ∈ [−2.12, −1.58]: serbest A_k — 7 ana hedef x = 5, 11/2, 6,
    13/2, 7, 15/2, 8; doğrusal taban a + bΔω (serbest); çeyrek aile x = 31/6, 16/3,
    17/3, 35/6, 37/6, 19/3, 20/3, 41/6, 43/6, 22/3, 23/3, 47/6, 49/6 (+ kenar: 29/6,
    14/3, 25/3) SABİT genlikle, iki modelde ayrı:
      M_BK: genlik_i = q_f^cal · 0.25/x_i;
      M_Rg: genlik_i = q_f^cal · 0.25/x_i · g(L + Δω_i)/ḡ_cal,
    q_f^cal = kalibrasyon penceresinde serbest uydurulan q_f; ḡ_cal = kalibrasyon
    penceresindeki çeyrek tepelerde g'nin (0.25/x ağırlıklı) ortalaması. Kenardan taşan
    ana tepeler (9/2, 17/2, 9) kör pencerede ağırlık taşımıyor (çizgi genişliği ≈ 0.01) ⇒
    modelden düşer.
  - KALİBRASYON uyum penceresi Δω ∈ [−1.30, −0.55] (görülmüş bant): serbest A_k — ana
    x = 2, 5/2, 3, 7/2; çeyrek aile x = 11/6, 13/6, 7/3, 8/3, 17/6, 19/6, 10/3, 11/3
    (+ kenar: 5/3, 23/6) ortak q_f ile; taban a + bΔω.
  - s = ters-varyans ağırlıklı ortalama {A(−log2)/(1/2), A(−log3)/(1/3)} (ağırlıklar
    jk se'lerinden). R_x = A_x/(s/x).
- RAKİP R-g'nin öngörüsü AYNI ağırlıklarla: R^g_x = g(L − log x) / ḡ,
  ḡ = Σ_a w_a g(L − log x_a)/Σ_a w_a (a ∈ {2, 3}, w_a = s tahminindeki ağırlıklar).
  (Eşit ağırlıkta R^g = 0.754, 0.702, 0.662, 0.629 (x = 5..8); ağırlık −log2'ye kaydıkça
  daha da küçülür.)

## Kapılar (K0 — ölçümden ÖNCE tanımlı, hepsi zorunlu)

- M1 (kod yolu): S1-S3 kapalıyken yeni kod 190'ın `harita_omega_dusuk.npz` K'sını
  ≤ 1e-10 bağıl farkla üretir (tam harita ya da Δω ∈ [−1.3, 3.5]'e yayılmış ≥ 12 ω-dilimi,
  50 K sütununun tümü).
- M2 (kontrol bandı, yeni konfig): 192'nin ana-6'sı (+log2, +log3, +log6, +log(3/2),
  −log2, −log3) yanar (P = Σ_{|c−h|≤0.03} κ_Σ > 3se) ve oranlar +log3/+log2, +log6/+log2,
  −log2/+log2, −log3/+log3 192 A2 değerlerinden (0.310, 0.602, 0.240, 0.466) ±%30 içinde.
  BAŞARISIZSA SINAV GEÇERSİZ (BK ölümü değil).
- M3 (öz-terim): tarama çizgileri ∩ HAVUZ' = ∅ (assert).
- M4 (güç): tam sayı ailesinin ağırlıklı ortalaması R̄_Z'nin σ_eff'i ≤ 0.10 (M_BK ve M_Rg) (BK 1.0 ile
  R-g ~0.68 arasını ≥ 3σ ayırmak için). Değilse tüm kör hükümler "belirsiz".
- M6 (uyum makinesi + güç, sentetik; ölçümden ÖNCE M1-tabanlı gürültü ölçeğiyle, hükümde
  ölçülen jk gürültüsüyle yeniden): sentetik profiller (aynı h(j; r), taban, 192 sinyal
  ölçeği) İKİ doğruyla üretilir — BK-doğru ve R-g-doğru — her biri ≥ 100 replika.
  Geçme: (i) gürültüsüz geri kazanım A_k, q_f %5 içinde; (ii) KARARA GİREN niceliklerde —
  R̄_Z ve ρ (σ_eff = f_Z·σ_jk, f_ρ·σ_jk; |Â − A| ≤ 2σ_eff) ile kalibrasyon çapaları
  A(−log2), A(−log3), q_f^cal ve kör ana genlikler A_x (x = 5..8, 11/2, 13/2, 15/2)
  (|Â − A| ≤ 2.36σ_jk, t₇) — kapsama ≥ %88; (iii) doğru hipotezin kazanma oranı ≥ %80 ve YANLIŞ
  hipotezin kazanma oranı ≤ %5 (her iki doğru için; f kalibrasyon replikalarından AYRI
  ≥ 200'er güç replikasında, σ_eff ile). Ölçüm-öncesi M6 kalırsa sınav
  YAPILMAZ (tasarım yeniden); hüküm-anı M6 kalırsa hüküm yalnız KAYIT.
- M5 (çizgi biçimi): kalibrasyon penceresinde uyum artığı rms'si ≤ 0.35 × κ_Σ'nin
  pencere-içi std'si VE −log2 tepesinin ölçülen/model profil korelasyonu (|Δω+log2| ≤ 0.06
  dilimleri) ≥ 0.85. Başarısızsa çizgi-biçimi modeli GEÇERSİZ ⇒ hüküm yalnız KAYIT.

## Hipotezler (eşikler DONMUŞ; kurtarma yok)

R̄_Z = x = 5, 6, 7, 8 için R_x'in ters-varyans ağırlıklı ortalaması (ağırlık 1/se_x²,
se'ler 8-grup jk; σ_Z = R̄_Z'nin jk se'si); R̄^g aynı ağırlıklarla R^g_x ortalaması.
β: A_x = C·x^{−β} ağırlıklı doğrusal-olmayan en küçük kareler (ağırlık 1/se_x², x = 5..8;
se jk). BK 1, R-g ≈ 1.4-1.5. Uyum penceresi üyeliği: dilim MERKEZİ pencere içindeyse.

- H-197a (zarf düzeyi, BK vs R-g; iki modelin öz-tutarlılığı): R̄_Z^BK, σ^BK ← M_BK uyumu;
  R̄_Z^Rg, σ^Rg ← M_Rg uyumu (σ'lar σ_eff = f·σ_jk). BK öz-tutarlı ⇔ |R̄_Z^BK − 1| ≤
  2σ^BK; R-g öz-tutarlı ⇔ |R̄_Z^Rg − R̄^g| ≤ 2σ^Rg.
  - BK TUTAR: BK öz-tutarlı, R-g DEĞİL. (β^BK ve se'si KAYIT.)
  - R-g KAZANIR (BK'nın çıplak eşlemesi ÖLÜR): R-g öz-tutarlı, BK DEĞİL.
  - İKİSİ DE ÖLÜR: |R̄_Z^BK − 1| > 3σ^BK VE |R̄_Z^Rg − R̄^g| > 3σ^Rg.
  - Diğer (ikisi de öz-tutarlı ya da ara durum): belirsiz (KAYIT).
  H-197b ve H-197c M_BK uyumundan hesaplanır; M_Rg ile aynı hükümler KAYIT.
- H-197b (buçuklu aile aynı zarfta; kalibrasyondan bağımsız): ρ = R̄_H/R̄_Z, R̄_H = x = 5.5,
  6.5, 7.5 ağırlıklı ortalaması; se jk. TUTAR: |ρ − 1| ≤ max(0.25, 2σ_ρ) (σ_ρ = f·σ_jk). ÖLÜM: ρ < 0.50
  (a = 2 ailesi bastırılmış, μ-benzeri koherens) ya da üç hedeften ikisinde A < 0.
- H-197c (asal kuvveti aynada bastırılmaz; kalibrasyondan bağımsız): ψ = A_8/Â_8,
  Â_8 = x = 5, 6, 7'ye uydurulan A ∝ x^{−β} zarfının x = 8'deki değeri. TUTAR:
  ψ ∈ [0.60, 1.60]. ÖLÜM: ψ < 0.30 (μ(8) = 0 sönmesinin aynaya uygulanması).

## KAYIT kalemleri

- q_f / R̄_Z (BK: 1 — çeyrek ailelerin zarfı ana zarfın ¼'ü).
- χ²(BK) vs χ²(R-g) kör hedeflerde; β ve se'si.
- İKİNCİL (görülmüş bant, HÜKÜM DIŞI): 32 bloklu dar çizgi-biçimiyle 192'nin 23 konumu aynı
  profil uyumuyla yeniden okunur (+ yan için pencere Δω ∈ [0.15, 2.35], tepeler BK c ≥ 0.05
  + μ=0 konumları serbest). 196 şerhinin (ii) okuması: karışık oranlar (+log(5/3),
  +log(5/2), −log(4/3), +log(5/4)) ve μ=0 çukurları BK'ya yaklaşıyor mu? Yaklaşmıyorsa
  (i) okuması (çekirdek ≠ R₂ Bohr katsayıları) güçlenir.

## Ölçüm notları

- Girdi: 190k0 'dusuk' zinciri (scratchpad/190/zincir_dusuk), zeros6 Z[200000:500000],
  L = 10.48393, N = 299 999. 190'ın τ' > 0.86 dilim serileri (scratchpad/190/dilim_*.npz)
  yeniden kullanılır; yeni dilimler ayrı dizine (scratchpad/197).
- Hesap: yeni τ'-dilimleri (q' ∈ (~2340, 8250)); 32 bloklu ω-dilim izdüşümü
  Δω ∈ [−2.2, +2.4]; tam koşu ~dakikalar (makine raporu).
- Ajanlar git'e dokunmaz; K0 (ONKAYIT_197.json: KALEM sha256 + girdi sha256'ları + damga)
  push edilmeden ölçüm başlamaz.
