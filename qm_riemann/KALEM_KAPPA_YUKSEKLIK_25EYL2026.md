# KALEM — 198: κ_p'NİN YÜKSEKLİK BAĞIMLILIĞI — paydaki asal sönmesi sonlu-yükseklik düzeltmesi mi, sabit yapı mı?
(25 Eylül 2026 — kullanıcı onaylı; TÜRETİM TEFTİŞİNDEN GEÇTİ (Sonnet) — iki KRİTİK bulgu bu
sürümde kapatıldı: (1) W_son yarı-körlüğü nicel olarak gerekçelendirildi (görülen 192 oranları
H_S/H_C'yi ayırmıyordu) ve tasarım tutarlı hâle getirildi; (2) ölçüm geometrisi karıştırıcısı
(W_alt'ta 32-blok içi L yayılımı 0.026-0.217, düşükte 0.017-0.041, sonda ~0.005 — genişleyen
çizgi biçimi H_S yönünde sahte sönme üretebilir) → eşit-ΔL bloklar (≤ 0.03), pencere başına
sentetik enjeksiyon yanlılık kapısı (M8) ve havuz değişmezliği kapısı (M7). L_alt gerçek
sıfırlardan 9.3426. Eşikler M6'dan sonra donar.
MAKİNE RAPORU (ölçümden önce; kör pencerelerde hiçbir κ/K/P hesaplanmadı): W_alt N = 179 999,
L_W = 9.342576, L ∈ [7.963, 10.008], zincir 198k0 (8 sn), HAVUZ' 166 çizgi; W_son zinciri tam
örneklem (N = 299 999). Eşit-ΔL bloklar zincirin 8 jk grubu İÇİNDE (mix loo tutarlılığı):
W_alt 72 blok (yayılım 0.024-0.029), W_düşük 35, W_son 36. M1 GEÇTİ (s ve 19 İKİNCİL genliği
bağıl fark 0). M7 GEÇTİ (χ²₄ = 5.79, p = 0.216; havuz τ ∈ [0.45, 0.74) kalır). M6: f_γ = 1.242,
f_Λ = 1.165; H-198a doğru 0.930 / 0.935, yanlış 0.010 / 0.010 ✓; H-198b gücü 0.33 / 0.40 < 0.80
⇒ H-198b KAYIT (kural gereği). M8 (ilk tanım) KALDI: tek tek κ_p yanlılığı 0.005'i aşıyor (en
kötü W_düşük κ_5 −0.022; kaynak s'deki ORTAK-MOD kalibrasyon yanlılığı — çeyrek oranları sabit
0.25/x, doğruda κ_3, κ_2κ_3 ile sönük; b/5 ailesi modelde yok); karar niceliği γ̂'nın yanlılığı
+0.013 (0.054σ_γ). KAPTAN KARARLARI: (i) W_alt kapsaması KISMİ (192 kuralı: dilim j'ye yalnız
onu TAM kapsayan bloklar; çizgi biçimi aynı B_j ile; 13/72 düşük-L blok alt uçta kapsamıyor,
B_j ≥ 59/72); (ii) M8 KARAR NİCELİĞİ ile: |bias_γ| ≤ 0.25σ_γ (GEÇİYOR); tek tek κ_p
yanlılıkları KAYIT ve mutlak κ değerleri ~0.02 model yanlılığı taşır (197'nin veri-sonrası
κ_p ≈ p^{−1/3} okumasına da not düşülür).)

## Durum

197 (kör, mühürlü): BK'nın ayna yanı (negatif üsler, f_p(k<0) = p^k) TAM tutuyor
(R̄ = 1.023 ± 0.023). VERİ-SONRASI (197 İKİNCİL, düşük pencere L_d = 10.484): pozitif üslü
(k_p = +1, paydaki) asallar BK'ya göre asal-başına bir çarpanla sönük:
κ_2 = 0.795 ± 0.015, κ_3 = 0.720 ± 0.012, κ_5 = 0.608 ± 0.018, κ_7 = 0.530 ± 0.020 (jk),
κ_p(W) := A_W(+log p) / (s_W · c_BK(p)), s_W = aynı penceredeki ayna-yanı ölçeği.
Görünüş κ_p ≈ p^{−1/3}; çarpımsallık yaklaşık (log6 1.14×, log10 1.21×).

SORU: κ_p yükseklikle 1'e mi gidiyor (BK asimptotik bir sanı; sonlu L düzeltmesi), yoksa
L'den bağımsız mı (ölçtüğümüz çekirdeğin ya da sıfırların gerçek bir yapısı)?

## Kalem cebiri (kaptan, veri-öncesi)

Genel biçim: 1 − κ_p(L) = B_p · (L/L_d)^{−γ}, γ ortak.
- H_S (sonlu yükseklik): γ = 1 (BK'nın ilk düzeltmesi O(1/L) — tipik alt-mertebe terimler
  log p/L ya da 1/L ile ölçeklenir).
- H_C (sabit yapı): γ = 0.
Öngörüler (B_p = 1 − κ_p^d, düşük pencereden; pencere L_W K0'da kesinleşir):
| pencere | L_W (yaklaşık) | H_S: κ_2, κ_3, κ_5, κ_7 | H_C |
|---|---|---|---|
| W_alt | 9.3426 | 0.770, 0.686, 0.560, 0.473 | 0.795, 0.720, 0.608, 0.530 |
| W_düşük | 10.484 | (referans) | (referans) |
| W_son | 12.030 | 0.821, 0.756, 0.659, 0.591 | 0.795, 0.720, 0.608, 0.530 |

## Pencereler ve görülmüşlük beyanı

- W_alt (YENİ, KÖR): zeros6 Z[20 000 : 200 000] (~1.8×10⁵ aralık; L ≈ 7.9-10.0; kesin
  sınırlar ve L_W = Σ_b n_b L_b / Σ n_b K0'da). Bu pencerede tarak uyduları HİÇ
  hesaplanmadı (K0'da 155/184-197 pencere listeleriyle doğrulanır).
- W_düşük (190/197 zinciri, Z[200 000 : 500 000]): REFERANS; 197'nin İKİNCİL sayıları
  aynı boru hattıyla yeniden üretilir (M1).
- W_son (190'ın 'son' zinciri, L_son = 12.030): YARI-KÖR — 192'de 8-blok A2 oranları ve B
  kataloğu son sütunu görüldü; 192 KEŞİF(c) "göreli parlaklıklar kabaca yükseklikten
  bağımsız" (fark −2.3σ…+0.8σ) dedi. NİCEL GEREKÇE: o oranlar H_S/H_C'yi AYIRMIYORDU —
  −log2/+log2 ∝ 1/κ_2: H_S 0.232, H_C 0.240, görülen 0.240 ± 0.014 (0.6σ); +log3/+log2 ∝
  κ_3/κ_2: H_S 0.315, H_C 0.310, görülen 0.314 ± 0.008 (ikisine de uyumlu). Bu yüzden W_son
  H-198a'ya GİRER; tamamen kör karar ayrıca H-198b (yalnız W_alt) ile verilir. 32-blok
  profil uyumu ve κ_p tanımıyla hiç ölçülmedi.

## Ölçüm tanımı (197 makinesi AYNEN; pencereye özgü parametreler)

- Her pencerede: HAVUZ' = τ ∈ [0.45, 0.74) (τ = log q / L_W); K_düz BİRİNCİL; BLOKLAR
  EŞİT-ΔL: blok sayısı n_B = max(32, ceil(ΔL_W/0.03)), her blok eşit L genişliğinde (≤ 0.03;
  W_düşük ~32 × 0.026, W_alt ~69 × 0.030, W_son 32 × ~0.005); jk 8 bitişik grup (gap
  sayısınca olabildiğince eşit); blok-yerel Δω; tarama τ' ≥ 0.74, üst sınır her blokta Δω ∈ [−1.40, +2.40]
  TAM kapsanacak şekilde (W_alt'ta τ' üst sınırı gerekirse 1.30'un üstüne çıkar).
- 197c fonksiyonları AYNEN: kalibrasyon penceresi [−1.30, −0.55] (s_W, q_f^cal), İKİNCİL +
  yan penceresi [0.15, 2.35] (34 bileşen) ⇒ κ_p(W) = A_W(+log p)/(s_W·c_BK(p)),
  p ∈ {2, 3, 5, 7}; se = 8-grup jk.
- γ̂: 3 pencere × 4 asal üzerinde ağırlıklı doğrusal-olmayan en küçük kareler,
  1 − κ_p(W) = B_p (L_W/L_d)^{−γ}, serbest B_p (4) ve ortak γ; ağırlık 1/se². σ_γ: pencere
  başına jk loo'larının birleşimi (pencereler bağımsız) × f_γ (M6 kalibrasyonu).
- Λ_alt (kör bacak): Λ = Σ_p w_p (1 − κ_p^alt)/(1 − κ_p^d) / Σ w_p (ters-varyans).
  H_S: Λ = L_d/L_alt; H_C: Λ = 1. σ_Λ × f_Λ.

## Kapılar (K0 — ölçümden ÖNCE, hepsi zorunlu)

- M1 (kod yolu): W_düşük'te yeni boru hattı 197'nin HUKUM_197.json İKİNCİL genliklerini ve
  s'yi ≤ 1e-10 bağıl farkla üretir.
- M2 (her pencere): 192 ana-6'sı yanar (> 3se); ayna ölçeği tutarlı: A(−log2)·2 /
  A(−log3)·3 ∈ [0.85, 1.15].
- M3: tarama ∩ HAVUZ' = ∅ (her pencere).
- M6 (güç, sentetik; ölçümden önce): H_S-doğru ve H_C-doğru sentetik κ_p üçlüleri (197'nin
  B_p'leri, pencere başına gürültü = 197 düşük-pencere jk se'si × √(N_d/N_W)), her doğru
  için ≥ 200 replika; f_γ, f_Λ AYRI tohumlarla kalibre. Geçme: doğru hipotez ≥ %80,
  yanlış ≤ %5 (H-198a için). Kalırsa sınav YAPILMAZ (pencere tasarımı yeniden).
- M4 (güç, gerçek): σ_eff(γ) ≤ 0.35; değilse H-198a belirsiz.
- M5' (uyum kalitesi, HER pencere): 197 M5 AYNEN (artık rms ≤ 0.35 × pencere std; −log2
  profil korelasyonu ≥ 0.85) kalibrasyon ve + yan pencerelerinde; ayna R̄_W (tam sayılar
  x = 2, 3 için A·x / s_W) ∈ [0.9, 1.1]. Kalan pencere karardan çıkar (rapor edilir).
- M7 (havuz değişmezliği, W_düşük — görülmüş veri): κ_p iki havuz yarısıyla (τ ∈ [0.45,
  0.60) ve [0.60, 0.74)) ayrı ölçülür; 4 asal üzerinden χ² p-değeri > 0.01. KALIRSA tüm
  pencerelerde ORTAK MUTLAK HAVUZ q ∈ [224, 992] (üç pencerenin ortak aralığı) kullanılır;
  M6/M8 bu havuzla yeniden koşulur.
- M8 (geometri yanlılığı, HER pencere, sentetik enjeksiyon): pencerenin kendi blok yapısı
  ve çizgi biçimleri h_W(j; r) ile H_C-doğru (κ_p sabit) sentetik profiller + ölçülen
  gürültü düzeyi; KARAR NİCELİĞİ γ̂'nın yanlılığı |bias_γ| ≤ 0.25 σ_γ (makine raporu sonrası
  tanım; tek tek κ_p yanlılıkları KAYIT). Kalırsa sınav YAPILMAZ.

## Hipotezler (eşikler teftiş + M6 sonrası donar; kurtarma yok)

- H-198a (γ, üç pencere): H_S TUTAR ⇔ |γ̂ − 1| ≤ 2σ_eff VE |γ̂| > 2σ_eff.
  H_C TUTAR ⇔ |γ̂| ≤ 2σ_eff VE |γ̂ − 1| > 2σ_eff. İKİSİ DE ÖLÜR ⇔ |γ̂| > 3σ_eff VE
  |γ̂ − 1| > 3σ_eff. Diğer: belirsiz.
- H-198b (TAMAMEN KÖR bacak, W_alt): aynı kural Λ_alt için (H_S: L_d/L_alt = 1.122;
  H_C: 1). M6'da H-198b'nin gücü < %80 çıkarsa H-198b yalnız KAYIT.

## KAYIT

- Pencere başına κ_p, p^{−1/3} ile kıyas; bileşik/karışık oranların κ-çarpımına oranı
  (+log6, +log10, +log(3/2), +log(5/2), +log(5/4), +log(5/3), −log(3/2), −log(4/3));
  buçuklu aile ve q_f/s; ayna R̄ (1 beklenir); μ=0 konumlarının negatif artıkları (L ile
  değişiyor mu?).
- Λ_son (H_S: 0.871; H_C: 1) — yarı-kör bacak.
- Çapasız oran Λ_son/Λ_alt = (1−κ_p^son)/(1−κ_p^alt) ağırlıklı (H_S: L_alt/L_son = 0.777;
  H_C: 1) — düşük pencereyi (desenin fark edildiği yer) tamamen dışarıda bırakır; seçim
  etkisi (ortalamaya dönüş) kontrolü.

## Ölçüm notları

- W_alt için zincir (190k0 muadili) kurulmalı: 155 PENCERE / 184-188 boru hattı; süre ve
  yeniden-kur gereksinimi makine raporunda.
- Ajanlar git'e dokunmaz; K0 (ONKAYIT_198.json: sha256 + damga) push edilmeden ölçüm
  başlamaz. W_alt ve W_son'da ölçüm-öncesi hiçbir κ/K/P değeri görülmez.
