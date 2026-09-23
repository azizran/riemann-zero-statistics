# KALEM — 192: UYDULARIN TEORİSİ — tarak uydusu kalıntı sınıfını okur (μ(n)/φ(n) seçim kuralı)
(23 Eylül 2026 gecesi — derin ikiz 1.30 inşa edilirken; kullanıcı onaylı A+B)

## Durum

188-190: pencere-ötesi iptal, orta-nokta örgüsünün Bragg tarağı ve
uydularıyla taşınır; uydular Δω = ω' − L_yerel ekseninde log 2, log 3, log 6'da
ve yükseklikten bağımsız (190 MÜHÜR). KONUMLAR anlaşıldı; PARLAKLIKLAR değil
(son pencerede 1+log2 0.056 ≫ 1+log3 0.015; artı/eksi asimetrisi ~4×).

## Kalem cebiri (kaptan, veri-öncesi türetim)

1. Seviye koşulu N(t_n) = n − ½, N = N̄ + S, S(t) = −Σ_q a_q sin(ω_q t),
   a_q = Λ(q)/(π√q log q). 2πN̄(t) = tL(t) − t + 7π/4, L(t) = log(t/2π).
2. e^{iω t_n} = −e^{i(ω t_n − 2πN̄(t_n))}·e^{−2πiS(t_n)}; Jacobi–Anger:
   e^{−2πiS} = Π_q Σ_k J_k(2π a_q) e^{ik ω_q t}.
3. {k_q} terimi için faz φ(t) = (ω' + Σ k_q ω_q) t − 2πN̄(t); DURAĞAN nokta
   L(t*) = ω' − log n, n := Π q^{−k_q} (rasyonel olabilir) ⇒ uydular
   Δω = ω' − L = log n (190'ın evrenselliği buradan).
4. Durağan fazdaki değer: φ* = t* − 7π/4, t* = 2π e^{L(t*)} = 2π q'/n. n = a/b
   (en sade) için **φ* ≡ 2π q' b/a − 7π/4 (mod 2π) — faz q' mod a'ya bağlı.**
5. Çizgiler (asal kuvvetler) a'ya asal kalıntı sınıflarına eşit dağılır ⇒
   uydunun koherens çarpanı = sınıflar üzerinde e^{2πi r/a} ortalaması =
   c_a(1)/φ(a) = **μ(a)/φ(a)** (Ramanujan toplamı).
6. Genlik (birinci mertebe): tek asal p için |J_1(2πa_p)| ≈ π a_p = 1/√p;
   uydudaki çizgi sayısı ∝ n, çizgi ağırlığı a' ∝ 1/√n ⇒ artı-taraf kare-siz
   tamsayı n için parlaklık ∝ |μ(n)|/φ(n); işaret (−1)^{ω(n)}·μ(n) = +1 ⇒
   HEPSİ aynı yönde (iptal) — gözlenen 180° ile tutarlı. Eksi taraf
   (n = 1/m, a = 1): koherens 1 (her m için), parlaklık ∝ 1/m.
   Akrabalık notu: Landau–Gonek formülü (Σ x^{iγ} ↔ Λ(x)/√x) ailesinden; x ~ t/2π
   rejimindeki bu uydu yapısının literatürdeki yeri KONTROL EDİLECEK
   (yenilik iddiası yok).

## Öngörüler (ön-kayıt)

**B — SEÇİM KURALI KATALOĞU (Δω ∈ [−1.3, 2.3]; n = a/b sade):**
- YANMALI (μ(a) ≠ 0): +log 2 (0.693), +log 3 (1.099), +log 5 (1.609),
  +log 6 (1.792), +log 7 (1.946), +log 10 (2.303, kenar), +log(3/2) (0.405),
  +log(5/3) (0.511), +log(5/2) (0.916), +log(5/4) (0.223), −log 2 (−0.693),
  −log 3 (−1.099), −log(3/2) (−0.405; a=2), −log(4/3) (−0.288; a=3).
- SÖNMELİ (μ(a) = 0 — birinci mertebe değil, TÜM mertebelerde koherens sıfır):
  +log 4 (1.386), +log 8 (2.079), +log 9 (2.197), +log(4/3) (0.288),
  +log(8/3) (0.981), +log(8/5) (0.470), +log(9/2) (1.504), +log(9/4) (0.811),
  +log(9/5) (0.588).

**A1 — KALINTI-SINIFI DÖNMESİ (en keskin; hiçbir verisine bakılmadı):** bir
uydunun çizgilerini q' mod a sınıflarına ayırınca her sınıfın karmaşık katkısı
K_r, toplamın yönüne göre (2π r b/a) kadar dönük durur:
- +log 3: r = 1, 2 sınıfları birbirine göre **120°±15°** ayrık; |K_r|/|K_top| ∈ [0.7, 1.4]
  (her sınıf tek başına toplam kadar güçlü; toplam −½ koherensten).
- +log 5: dört sınıf ardışık **72°±15°** aralıklı; |K_r|/|K_top| ∈ [0.6, 1.5].
- +log 6: r = 1, 5 sınıfları **120°±15°** ayrık (±60° toplam yönünden).
- **Kontrol, −log 3:** r = 1, 2 sınıfları AYNI yönde (fark < 15°), |K_r|/|K_top| ∈ [0.35, 0.65].
- **Kontrol, +log 2:** tüm çizgiler tek sınıf (tek), dönme yok.

**A2 — PARLAKLIK ORANLARI (ikincil; son penceresi görülmüştü → sınav DÜŞÜK
pencerede):** artı taraf, log 2'ye göre: log3 0.5, log5 0.25, log6 0.5, log7 0.17,
log10 0.25 (bant ×/÷ 2); eksi/artı: −log2/+log2 ≈ 0.5, −log3/+log3 ≈ 0.67 (bant ×/÷ 2).
Not: son penceresinde görülen oranlar (0.27, 0.08, 0.62, 0.03) teoriden daha hızlı
düşüyor — bilinen gerilim; A2 bunu düşük pencerede ölçer.

## Hipotezler ve ölüm koşulları

- **H-192a (seçim kuralı, B):** "sönmeli" listesinin HEPSİ söner (|κ| < 2 jk-se
  ya da yerel taban içinde) VE "yanmalı" listesinin ana 6'sı (±log2, ±log3,
  +log6, +log(3/2)) yanar (> 3 jk-se). ÖLÜM: sönmeli listeden herhangi biri
  > 4 jk-se yanarsa (seçim kuralı yanlış).
- **H-192b (kalıntı dönmesi, A1):** +log 3 ve −log 3 koşulları ikisi birden
  tutar. ÖLÜM: +log 3'te sınıflar aynı yönde (< 45°) ya da −log 3'te > 45° ayrık.
- **H-192c (parlaklık, A2):** ikincil, KAYIT niteliğinde bant sınavı.

## Kapılar ve İŞLEMCİ KURALI

- **İŞLEMCİ:** derin ikiz 1.30 inşası 7 çekirdeği kullanıyor. Bu kalemin bütün
  koşuları `nice -n 19` ve TEK süreçle; ağır tam-harita yeniden koşusu YOK.
- **K0 — ÖN-KAYIT:** bu kalemin commit'i (push'lu hash + damga) + ajan ONKAYIT_192.json.
- **K1 — B KATALOĞU:** mevcut Δω profillerinden (190: harita_omega_dusuk,
  profiller_190; son: 188 dosyaları) — yeni koşu gerekmez; blok jackknife se.
- **K2 — A1 SINIF AYRIŞTIRMASI:** son penceresinde (188 makinesi AYNEN, seri
  fonksiyonu bit-bit mühürlü), yalnız uydu pencerelerindeki çizgiler
  (+log3: τ'∈[1.078,1.104]; −log3: τ'∈[0.896,0.922]; +log5: τ'∈[1.121,1.146];
  +log6: τ'∈[1.137,1.161]; +log2: τ'∈[1.045,1.070]) q' mod a'ya göre ayrılıp
  karmaşık K_r (HAVUZ) hesaplanır; alt-örneklem her 3. nokta.
- **K3 — A2 + HÜKÜM + figür.**

Ölümler kurtarmasız; ajan git'e DOKUNMAZ; sonuç ORTAK TEFTİŞE (sabah).
