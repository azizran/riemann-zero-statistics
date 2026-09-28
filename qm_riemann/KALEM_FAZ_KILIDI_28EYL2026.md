# KALEM — 201: FAZ KİLİDİ — sıfırlar asal dalgalarını hangi fazda yakalıyor, ve bu asal izinin sönmesini açıklıyor mu?
(28 Eylül 2026 — TÜRETİM TEFTİŞİNDEN GEÇTİ (Sonnet): açıklama-gücü sınavı KARAR olmaktan çıkarıldı (harmonik kesimine göre işaret değiştiriyor, asallar arası kovaryans yok); çift kilit öngörüsü düzeltildi (yerel GUE itmesi ⇒ 4×; adaylar 1/2/4); sınıflama kuralının yapay kesinlik tuzağı giderildi; asallar arası jackknife kovaryansı; ε sağlamlığı. Kaynak fikir: kaptanın "tepede mi çukurda mı?"
sorusu (27 Eyl) ve `problem_nesnesi.svg`'deki C tutamağı. Eşikler M6'dan sonra donar.)

## Durum

200-A/B/C (mühürlü): κ₃ asal izi koşullanan sıfır sayısıyla sönüyor (L≈11: +0.283 → +0.150 → +0.086),
sönme yükseklikle kapanıyor (L = 22'de +0.28 / +0.184 / +0.147); κ₂ izi her basamakta ≈ a_k.
Standart hibrit resimde (Gonek–Hughes–Keating) Euler kısmı sıfır koşullamasından BAĞIMSIZ — o
zaman sönme olmamalıydı. Aday mekanizma: sıfırlar asal fazlarını kilitliyor.

## Kuram (Landau 1911; Gonek 1993 düzgün hâli; türetim teftişi yapılacak)

Σ_{0<γ≤T} x^{iγ} = −(T/2π)·Λ(x)/√x + (küçük terimler) ⇒ bir L-penceresinde sıfırlar üzerinden
ortalama: **E[cos(m·γ log p)] = −log p / (p^{m/2} · L)**, E[sin(·)] ≈ 0.
Sezgi: sıfır yoğunluğu n(t) ≈ (L/2π)[1 − (2/L)Σ_p (log p/√p) cos(t log p) + …] — asal dalgası
tepedeyken sıfır azalır. 2 için L = 11'de: E cos φ = −0.042, E cos 2φ = −0.030 ⇒ en olası faz ≈ ±111°
(`faz_kilidi_ongoru.svg`). Yakın çift (iki sıfır, aynı yükseklikte): yoğunluklar çarpılır ⇒
yakın çift için birinci-mertebe öngörü aşağıda (adaylar 1 / 2 / 4 kat).

**Kilidin asal izine etkisi (KAYIT, karar YOK):** bağımsız Euler çarpanları modeliyle çeviri
kararsız çıktı (teftiş): L = 10.59, çift kilitte Δκ₃ harmonik kesimine göre İŞARET değiştiriyor
(m ≤ 2: +0.026; m ≤ 7: −0.009; m → ∞: −0.009); asallar arası kovaryans modelde yok (ortak γ ile
bütün asalların eğimleri üst üste binince "yoğunluk" negatife düşüyor). Sağlam olan tek cümle:
**basit bağımsız-Euler çevirisi, gözlenen κ₃ sönmesini (b1−b0 ≈ −0.13) açıklayacak büyüklükte
değil (tek kilitte ters işaretli +0.005) ve κ₂'de gözlenmeyen bir düşüş (−0.078) öngörüyor.**
"Kilit sönmeyi açıklıyor mu?" sorusu bu kalemde KARARA bağlanmaz; ölçülen kilit yalnızca KAYIT olarak
bu çeviriye beslenir.

## Veri ve gözlenebilirler (yalnız sıfır KONUMLARI; hiç hesaplanmadı)

Pencereler: zeros6 W1–W3 (L 9.35 / 10.59 / 11.66) + LMFDB C1–C4 (14.19 / 16.58 / 18.88 / 22.31),
200-A/B/C ile aynı sıfır kümeleri. Asallar p ∈ {2, 3, 5, 7, 11, 13}; harmonikler m = 1, 2.
 - b=1: bütün sıfırlar γ_n, φ = γ_n log p (mod 2π).
 - b=2: yakın çiftler (δ̃ < 0.2), φ = m_n log p, m_n orta nokta.
 - b=0 (denetim): 10⁶ düzgün rastgele t ⇒ c_m ≈ 0 olmalı (boru hattı sağlaması).
 İstatistikler: ĉ_m(p) = ortalama cos(mφ), ŝ_m(p) = ortalama sin(mφ); SE: 64 t-bloklu jackknife.
 Kilit oranı: R_b = ĉ_1 / c_1^{LG} (asallar ve pencereler üzerinde ters-varyans havuzu; m=2 ayrıca).

## Hipotezler ve karar kuralları

Önce çift için doğru birinci-mertebe öngörü (teftiş sonrası): asal dalgaları sıfır yoğunluğunu
yavaşça modüle eder (2 için dalga boyu 2π/log 2 ≈ 9 ≈ 16 ortalama aralık). Yerel olarak GUE
istatistiği yerel yoğunluk n(t) ile geçerliyse, KÜÇÜK aralıklı çiftlerin yoğunluğu
R₂(t; δ) ≈ n(t)²·[1 − sinc²(π δ n(t))] ≈ n(t)⁴·π²δ²/3 ⇒ orta nokta kilidi tek sıfırınkinin 4 katı.
Üç aday: **R₂ = 1** (TEK: çift, tek sıfır gibi), **R₂ = 2** (ÇİFT: yoğunlukların çarpımı, itme
yoğunluğa duyarsız), **R₂ = 4** (DÖRT: yerel GUE itmesi yerel yoğunlukla ölçeklenir).
(Not: çiftler δ̃ = δ·L/2π ile DÜZGÜN açılımla seçiliyor; yerel yoğunlukla değil — R₂ = 4
öngörüsünün dayanağı bu.)

 - **H-201-1 (b=1, Landau–Gonek sınavı):** R₁ = ĉ₁ / c₁^{LG} (p ∈ {2,…,13}; pencere başına GLS,
   asallar arası tam jackknife kovaryansıyla; pencereler arası ters-varyans). |R₁ − 1| ≤ 0.15 ve
   ≤ 3σ ⇒ TUTAR. Ayrıca m = 2 için R₁^{(2)} = ĉ₂/c₂^{LG} (aynı kural, KAYIT) ve 1/L ölçeklemesi:
   7 pencerede R₁'in log L'ye eğimi (|eğim| ≤ 0.1 ⇒ ölçekleme TUTAR).
 - **H-201-2 (b=2, yeni; orta nokta, δ̃ < 0.2):** R̂₂ ± σ (aynı GLS). Sınıf c ∈ {1, 2, 4}:
   |R̂₂ − c| ≤ 0.25·c VE diğer her sınıfın ±%25 bandından ≥ 3σ uzak ⇒ KESİN: c. Hiçbir sınıfın
   ±%25 bandında değilse ⇒ ARA DEĞER (R̂₂ raporlanır; tamsayı olmayan bir kilit kendi başına
   bulgu). Aksi ⇒ BELİRSİZ.
 - **İKİNCİL (KAYIT):** ε ∈ {0.1, 0.3} için R̂₂; b=2'de tek tek sıfırların (orta nokta yerine)
   fazları; ŝ_m ≈ 0 kontrolü; ölçülen kilidin bağımsız-Euler çevirisi (m ≤ 2 ve m ≤ 20) ve
   gözlenen κ₂/κ₃ basamak farklarıyla karşılaştırması.

## Kapılar

**M6 (sentetik güç, 28 Eyl, geçti):** LG ile modüle edilmiş yoğunluktan inceltilmiş sentetik süreç
(L = 11, altı asal, m ≤ 2) — hat R₁ = 0.932 ± 0.004 (m=1), 0.952 (m=2) geri buluyor; orta noktaları
n(t)^k ile örneklenen çiftler: k = 1 → 0.993, k = 2 → 1.850, k = 4 → 3.242. Üç sınıf ayırt ediliyor;
büyük modülasyonda doğrusal olmama nedeniyle ~%7 (k=1) ile ~%19 (k=4) aşağı kayma var — sınıf
bantları (±%25) içinde. (`201_configs/201_olcum.py --sentetik`, `201_analiz.py --sentetik`.)



M1: b=0 denetimi (|ĉ_m| ≤ 3 SE); M6: sentetik güç (LG kilitli fazlarla örneklenmiş sahte "sıfırlar"
⇒ R₁ ≈ 1 geri kazanılıyor mu; çift kilit ⇒ R₂ ≈ 2); M9: körlük (ölçümden önce faz istatistiği yok).

## Süreç

KALEM → Sonnet türetim teftişi → makine (basit; konum verisi) → ONKAYIT_201 → ölçüm → ortak teftiş.
