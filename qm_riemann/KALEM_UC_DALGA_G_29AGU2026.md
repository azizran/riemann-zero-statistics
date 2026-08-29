# KALEM OTURUMU — ÜÇ-DALGA MAKİNESİNİN YASASI (29 Ağu 2026)

Hedef: 142'nin Gram bağlaşımını (çarpımsal çiftlerde 9.6×) aracı-dalga
genliğinden türetmek. Sonuç: iki terimli, sıfır parametreli kapalı form;
19 çiftte 19 isabet (143).

## 1. Yer değiştirme alanı

Kesin sayımdan z_n = z̄_n − S(z_n)·ḡ_t ve S(z) = −Σ_c a_c sin(ω_c z)
⇒ tekne kayması d(z) = +ḡ_t Σ_c a_c sin(ω_c z). Bağ orta-noktası iki
komşunun ortalaması olduğundan yarım-gap kosinüsü belirir:

  δ(m̄) = ḡ_t Σ_c a_c cos(πτ_c) sin(ω_c m̄)

## 2. Gram açılımı ve hayatta kalma koşulu

G_cc(q₁,q₂) = 2⟨cos ω₁m cos ω₂m⟩, m = m̄+δ. Birinci mertebede türev
terimi ½Σω·sin(Σm) + ½ω₃·sin(ω₃m) verir (ω₃=ω₂−ω₁, Σ=ω₁+ω₂).
Katı-örgü ortalamasında yalnız δ'nın AYNI frekanstaki bileşenleri
hayatta kalır: fark frekansı ω₃ ⇔ q₃=q₂/q₁ örgüde (çarpımsal çift!);
toplam frekansı Σ ⇔ q₁q₂ asal-kuvvet (yalnız KULE çiftleri!).

## 3. Yasa (sıfır parametre; a_Q = 1/(π m_Q √Q), τ_Q = log Q/L)

  G_cc = −π[ τ₃ a₃ cos(πτ₃) + (τ₁+τ₂) a_Σ cos(πτ_Σ)·1{q₁q₂ örgüde} ]
  G_ss = −π[ τ₃ a₃ cos(πτ₃) − (τ₁+τ₂) a_Σ cos(πτ_Σ)·1{...} ]
  G_cs = 0

Kehanetler: aracı kuralı (eş-asal çiftte G yalnız q₃'e bağlı); kule
çiftleri Σ-terimiyle güçlü ve merdivende yukarı zayıflar; kanal
yarılması yalnız kulede, işareti cos(πτ_Σ) ile döner; çapraz kanal
özdeş sıfır — fazların hiç kıpırdamamasının mekaniği.

## 4. Sınav (143, son-300k) — HEPSİ TUTTU

- İşaret 19/19 negatif ✓. Değerler: oran ölçüm/öngörü = 0.91–0.98
  (ort ~0.95; tekdüze ~%5 açık — ikinci mertebe adayı, kayıtta).
- Aracı kuralı: (5,10)(7,14)(11,22)(13,26) → −0.0384/−0.0382/−0.0382/
  −0.0382 (±%1!) ✓
- Kule merdiveni: |G(2,4)|>|G(4,8)|>|G(8,16)| = 562>442>393 (öngörü
  575>464>416) ✓
- Kanal yarılması: (2,4) cc/ss = −562/−236 (öngörü −575/−226!) ✓;
  yarılmanın işaret DÖNÜŞÜ (25,125)'te tam öngörüldüğü yerde ✓
- G_cs ≤ 1e-4, 19/19 ✓ (gürültü tabanının altında)
- Kayıt: eş-asal çiftlerde küçük cc/ss asimetrisi (~0.0025, ss yüksek)
  — öngörü eşit diyor; ikinci-mertebe not.

## 5. Anlamı ve sıradaki

142'nin bütün fenomenolojisi (taban-akışı, çapraz-güç, kule-güçlenmesi,
soğurma) artık TÜRETİLMİŞ bir mikro-yasadan akıyor. Sıradaki kalem:
bu G'yi topluca ilerletip (Gram'ın tersinden) ρ_tail(τ) profilini,
C(n) çukur genliğini ve nihayet K(τ) çekirdeğini KOLEKTİF olarak
hesaplamak — makinenin kendisi artık elimizde.
