# Keşif Seferi — Kuantum Kaos / RMT Literatüründe Bulgularımızın Akrabaları
**Tarih:** 29 Ağustos 2026
**Kapsam:** Riemann/L-fonksiyon sıfır aralık istatistikleri ↔ periyodik yörünge teorisi, rastgele matris teorisi, aritmetik düzeltmeler
**Kural:** Aşağıdaki kaynakların hepsi gerçekten açılıp okundu. Okunamayan/yalnız künyesi doğrulananlar §5'te ayrı işaretlendi. Uydurma referans yok.

---

## 1. Yönetici Özeti

1. **Bulgu 1'in ikinci momenti LİTERATÜRDE VAR — hem de birebir.** `ds_n = Σ_q 2a_q sin(ω_q g_n/2) cos(ω_q m_n)` özdeşliğinin varyansını alınca prime-power başına `2a_q² sin²(πτ_q)` çıkar; bu **Berry (1988)** sayım varyansı formülünün δ=1 terimiyle terim terim aynıdır (§3.1'de cebir gösterildi).
2. Ancak Bulgu 1'in **örnek-yolu (sample-path) düzeyi** — tek tek `n` için 0.975 korelasyon — literatürde bulunamadı. Literatür ikinci moment düzeyinde kalıyor.
3. **Bulgu 2'nin (Nefes Yasası) "nesne"si akraba: VAR.** İki-nokta istatistiğinin her asal-kuvvete `cos(2πxτ_q)` fazıyla ve `log²p/p^m` ağırlığıyla ayrışması **Berry–Keating (1999) denk. (4.20)**'de aynen var. Yarım-aralık konumunda (x=1/2) bu tam olarak **cos(πτ)** verir — bizim `R_nn = −2cos(πτ)` bağ kanalıyla form düzeyinde örtüşüyor.
4. Fakat **koşullu-varyans / tepki okuması (R_nn) GÖRÜNMÜYOR.** Literatür asal terimleri korelasyon fonksiyonuna *katkı* olarak yazıyor; "tek bir asalın fazına koşullanmış artık varyans" diye bir nesne bulunamadı.
5. **Bulgu 3'ün ayrışımı (adyabatik +1 vs. öz-çekirdek −3) GÖRÜNMÜYOR** — ama üretecek mekanizma biliniyor: köşegen-dışı (off-diagonal) asal-asal korelasyonları (Hardy–Littlewood / Bogomolny–Keating / Sieber–Richter).
6. Kritik gözlem: köşegen yaklaşımında asal fazları bağımsızdır ⇒ `R_nn ≡ 0`. Dolayısıyla **`R_nn ≠ 0` doğrudan K_off'un (köşegen-dışı) bir ölçüsüdür.** Berry–Keating (1999, s. 248) bunu açıkça yazıyor: asal logaritmaları ikili olarak korelasyonsuz olsaydı, K_off rastgele fazların ortalaması olurdu ve sıfır çıkardı. Bu, bulgumuzun teorik yerini kesin biçimde çiviler.
7. **Bulgu 4'ün (L-değişmezlik) mekanizması AKRABA — VAR.** Berry–Keating (2.14): Riemann dinamiği *ölçekleyicidir*, tüm yörüngelerin kararsızlık üsteli λ_p = 1'dir; ayrıca (6.7)–(6.8): `H = XP` dilatasyon üretir, yani **zaman ötelemesi = ölçek değişimi**. Bu yüzden doğal değişken τ = T/T_H = log q / L'dir ve istatistikler t'den bağımsız hale gelir. `dτ/τ` ölçüsü ise **Hannay–Ozorio de Almeida** toplam kuralının ta kendisidir.
7b. **8 L-ailesindeki evrenselliğin literatürde hazır gerekçesi var.** Berry–Keating (6.9): her Dirichlet L-fonksiyonu `XP`'nin *farklı bir öz-eşlenik genişlemesidir* (χ karakterleriyle); **dinamik hepsinde aynıdır.** Yalnız dinamiğe bağlı bir istatistiğin ailelerde aynı çıkması beklenendir.
8. **Bulgu 5 (koherans / keskin-kesim) için tam bir literatür emsali VAR:** Berry–Keating (1999) denk. (5.24), Dirichlet serisinin keskin kesimini **Erfc yumuşatmasıyla** (Gauss'un integrali, saf Gauss değil) değiştiriyor. Yani "ne keskin ne Gauss, arada Erfc" ailesi zaten hazır alettir.
9. **Ölçü aralığımız tam olarak literatürün "evrensellik kırılma bölgesi"dir.** Berry–Keating: en kısa yörünge T₀ = log 2 olduğu için evrensellik ~L/log 2 ortalama aralıktan sonra bozulur. τ ∈ [0.06, 0.25] penceresi (L≈11) tam olarak q = 2'den q ≈ 13–16'ya kadardır — yani maksimum evrensel-olmayan pencere.
10. **Bonus / uyarı:** Programın `N_eff = L + c` bulgusu, Keating–Snaith'in yoğunluk eşleme kuralı `N = log(T/2π) = L` ile aynı ailededir; ama literatürde aralık dağılımına kalibre edilmiş **rakip** bir tanım var: Bogomolny–Bohigas–Leboeuf–Monastra, `N_eff = log(E/2π)/√(12Λ)`, Λ = 1.57314… (asal toplamlarıyla tanımlı). İki konvansiyon çakışmıyor; çapraz kontrol edilmeli.

---

## 2. Nesne–Nesne Eşleme Tablosu

| # | Bizim nesne | Literatür nesnesi | Kaynak (doğrulanmış) | Hüküm |
|---|---|---|---|---|
| A | Asal dalgası: `a_q = 1/(π m p^{m/2})`, `ω_q = log q` | `N_fl(t) = −(1/π) Σ_p Σ_m (p^{−m/2}/m) sin(t m log p)` | Berry–Keating SIAM Rev. 41 (1999) **denk. (2.6)** | **BİLİNEN** (birebir aynı nesne) |
| B | Göreli frekans `τ = log q / L` | Ölçekli zaman `τ = T/T_H`; Riemann'da `T_j = m log p`, `⟨d⟩ = L/2π` ⇒ `τ = m log p / L` | Berry–Keating **(2.14), (4.3), (2.21)** | **BİLİNEN** |
| C | `ds_n = S(z_n) − S(z_{n+1})` = asal dalgaları toplamı (korel. 0.975) | Sayım fonksiyonu ayrışımı `N = ⟨N⟩ + N_fl`, `N_fl = (1/π) Im log ζ` | Berry–Keating **(2.2)–(2.4), (2.6)** | Özdeşliğin kendisi **BİLİNEN**; 0.975'lik *nokta-bazlı* koherans iddiası **YENİ GÖRÜNÜYOR** |
| D | C'nin varyansı: `Σ_q 2a_q² sin²(πτ_q)` | Berry sayım varyansı (evrensel-olmayan rejim): `(1/π²) Σ_n [Λ²(n)/(n log²n)]·(1 − cos(2πδ log n / log T))` | Berry, *Nonlinearity* **1** (1988) 399; formül metni: arXiv **2211.14918**, Conj. 1.4.1 | **BİLİNEN** — δ=1'de terim terim özdeş (cebir §3.1) |
| E | `R_nn`'in `cos(πτ)` şekli | Evrensel-olmayan çift korelasyon düzeltmesi `R_c¹(x) ∝ Σ_{p,m} (log²p/p^m) cos(2πxτ_{p^m})`; **x = 1/2 ⇒ cos(πτ)** | Berry–Keating **(4.20)**; kapalı form **(4.23)** | **AKRABA** (form birebir; katsayı/işaret türetilmedi) |
| F | `R_nn` = tek asalın fazına koşullu artık varyans (tepki fonksiyonu) | — | (aranan, bulunamayan) | **YENİ GÖRÜNÜYOR** |
| G | `−2 = (+1.05 adyabatik) + (−3.05 öz-çekirdek)` ayrışımı | Form faktörünün `K = K_diag + K_off` ayrışımı; `K_diag = |τ|` (HOdA), `K_off = Θ(|τ|−1)(1−|τ|)` | Berry–Keating **(4.8)–(4.12)**; Hannay–Ozorio de Almeida (1984) | **AKRABA** (aynı ayrışım felsefesi, farklı gözlemlenebilir) |
| H | Öz-korelasyon çekirdeği `K(τ)` (−2.87 → −3.25) | Köşegen-dışı asal-asal korelasyonları (Hardy–Littlewood tekil serisi `C(k)`); yörünge-eylem korelasyonları | Berry–Keating **(4.13)–(4.16)**; Argaman ve ark., PRL **71** (1993) 4326; Bogomolny–Keating, *Nonlinearity* **8** (1995) 1115 | **AKRABA** — eğrinin kendisi **YENİ GÖRÜNÜYOR** |
| I | "Denge gazı bunu üretemez" | RMT'de aritmetik bilgi yoktur: "rastgele matrisler Riemann sıfırlarının konumları hakkında bilgi içermez" | Snaith, *Riemann zeros and RMT* (2009), §4 & Fig. 3 | **BİLİNEN** (aynı gözlem, aynı gerekçe) |
| J | L-değişmezlik mekanizması (PNT öz-benzerliği) | Riemann dinamiği ölçekleyici; λ_p = 1 (hepsi eşit kararsız); yörünge yoğunluğu `ρ(T) ~ e^{λT}/T` ⇒ PNT | Berry–Keating **(2.14), (3.1), (3.2)** | **BİLİNEN / AKRABA** |
| K | `dτ/τ` kuyruk ölçüsü | HOdA toplam kuralı: `lim (1/T) Σ_j A_j² δ(T − T_j) = 1` | Hannay–Ozorio de Almeida, *J. Phys. A* **17** (1984) 3429; Berry–Keating **(3.2)** | **BİLİNEN** |
| L | Gauss Debye–Waller sönümü fazla agresif; gerçek çekirdek keskin-kesime yakın | Riemann–Siegel keskin kesim `n* = ⌊√(t/2π)⌋` → **Erfc** ile yumuşatılmış kesim + optimizasyon parametresi K | Berry–Keating **(5.3), (5.24)** | **AKRABA — hazır alet** |
| M | Ölçüm penceresi τ ∈ [0.06, 0.25] | Evrensellik `~log(t/2π)/log 2` aralıktan sonra bozulur (T₀ = log 2 en kısa yörünge) | Berry–Keating §3 (s. 245) | **BİLİNEN** — penceremiz tam kırılma bölgesi |
| N | `N_eff = L + c` | (i) Keating–Snaith yoğunluk eşlemesi `N = log(T/2π)`; (ii) rakip: `N_eff = log(E/2π)/√(12Λ)`, Λ=1.57314 (asal toplamı) | Snaith (2009) **denk. (3.2)**; Bogomolny–Bohigas–Leboeuf–Monastra, *J. Phys. A* **39** (2006) 10743 | (i) ile **AKRABA**, (ii) ile **ÇELİŞKİ ADAYI** — çapraz kontrol gerekli |
| O | 8 L-fonksiyon ailesinde evrensellik | Ratios conjecture'ın L-fonksiyon ailelerine genellemesi; alt-mertebe terimler | Conrey–Snaith, arXiv **math/0509480**, *Proc. LMS* **93** (2007) 594 | **AKRABA** |
| P | Asal-dalga alanının öz-kovaryansı (bizim `S(z)` alanı) | Asal zeta fonksiyonu `P(½+it) = Σ_p p^{−½−it}` kritik doğruda asimptotik normal; öz-kovaryansı `log|ζ|` ile yaklaşık, ve **Riemann sıfır yüksekliklerine eşit ayrımlarda belirgin negatif** | Chavez & Allawala, *J. Stat. Mech.* (2021) 073206, arXiv **2102.02280** | **AKRABA** — aynı alanın kovaryansı, farklı gözlemlenebilir |
| R | `R_nn`'in 8 L-fonksiyon ailesinde aynı çıkması | `XP` dilatasyon simetrisi `X→KX, P→P/K`; **her L-fonksiyonu `XP`'nin farklı bir öz-eşlenik genişlemesidir** (karakterler χ(n), çarpımsal grup mod k) — dinamik aynı, yalnız genişleme değişir | Berry–Keating **(6.7)–(6.9)** ve izleyen tartışma | **BİLİNEN GEREKÇE** — aile-evrenselliğinin literatürdeki açıklaması |

---

## 3. Ana Kaynak Özetleri

### 3.1 M. V. Berry, "Semiclassical formula for the number variance of the Riemann zeros", *Nonlinearity* **1** (1988) 399–407
**Ne yapıyor.** Sayım varyansını yarı-klasik iz formülünden hesaplıyor: kısa mesafede GUE, uzun mesafede asallara kilitli, doyuma giden bir formül.
**Formülü (doğrulandı, arXiv 2211.14918 Conj. 1.4.1'den birebir):** evrensel-olmayan rejimde (δ ≫ log T)

```
∫₀ᵀ [S(t + 2πδ/log T) − S(t)]² dt
    = (T/π²) [ Σ_{n≤T} (Λ²(n) / (n log²n)) · (1 − cos(2πδ · log n / log T)) + 1 ] + o(T)
```

**Bizimkiyle tam ilişki — FORMÜL DÜZEYİNDE KARŞILAŞTIRILABİLİR, ve eşleşiyor.**
`n = p^m` için `Λ²(n)/(n log²n) = 1/(m² p^m)`, dolayısıyla `(1/π²)·Λ²/(n log²n) = 1/(π² m² p^m) = a_q²`.
Ve `log n / log T = τ_q`. Yani Berry'nin q-terimi:
`a_q² (1 − cos(2πδτ_q)) = 2 a_q² sin²(πδτ_q)`.
Bizim özdeşliğimizin q-bileşeni `2a_q sin(ω_q g_n/2) cos(ω_q m_n)`; `ω_q g/2 = log q · (2π/L)/2 = πτ_q` ve `⟨cos²⟩ = 1/2` ile varyansı:
`4a_q² sin²(πτ_q) · (1/2) = 2 a_q² sin²(πτ_q)`.
**δ = 1'de birebir aynı.** Yani Bulgu 1'in varyansı Berry 1988'dir; katkımız varyanstan değil, *örnek-yolu koheransından* ve *koşullu* istatistiklerden geliyor.
**Uyarı (dürüstlük payı):** Berry sabit bir `h` penceresinde `S(t+h) − S(t)` alıyor; biz "tam bir aralık" lag'inde alıyoruz. δ=1'de eşleşme sezgisel-kesin, ama koşullama farkı bir düzeltme doğurabilir; bu ayrıca hesaplanmalı.

### 3.2 M. V. Berry & J. P. Keating, "The Riemann Zeros and Eigenvalue Asymptotics", *SIAM Review* **41** (1999) 236–266
**Programımızın en yakın ve en kullanışlı kaynağı — tam metin okundu.**
- (2.6): `N_fl(t) = −(1/π) Σ_p Σ_m (exp(−½ m log p)/m) sin(t m log p)` — bizim asal dalgalarımız, birebir.
- (2.14): Kuantum ↔ Riemann sözlüğü. **Kararlılıklar: `½λ_p T_p ↔ ½ log p ⇒ λ_p = 1`.** Yani Riemann dinamiği homojen kararsız ve *ölçekleyici*; ℏ görünmüyor. **Bu bizim L-değişmezliğimizin literatürdeki gerekçesidir.**
- §3, s. 245: Evrensellik `log(t/2π)/log 2` ortalama aralıktan sonra kırılır, çünkü en kısa yörünge `T₀ = log 2`'dir. Ölçüm penceremiz (τ = 0.06 ≈ log2/L) tam bu eşikten başlıyor.
- (4.8)–(4.12): `K_diag(τ) = |τ|` (HOdA'dan), `K_off(τ) = Θ(|τ|−1)(1−|τ|)`. Ayrışım felsefesi bizimkiyle aynı.
- **(4.20) — bizim R_nn'in en yakın akrabası:**
  `R_c¹(x) = [1/(2(π⟨d⟩)²)] Σ_{m,p} (log²p / p^m) cos{ x m log p / ⟨d⟩ } − 2∫₀^{τ*} dτ cos(2πxτ)`
  `⟨d⟩ = L/2π` olduğundan faz `x m log p/⟨d⟩ = 2π x τ_q`. **x = 1/2 (bağ/yarım-aralık konumu) ⇒ cos(πτ).**
- (4.23): Kapalı form. İçinde `−∂²_ξ Re log ζ(1 − iξ)` ve `−Re Σ_p log²p/(p e^{iξ log p} − 1)²` var. `ζ(1−iξ)` terimi **"resurgence"**: yüksek sıfırların çift korelasyonunda *alçak Riemann sıfırları* rezonans olarak beliriyor. Bu, "gaz kendi aritmetiğini tanır" ifadesinin literatürdeki en sert hali.
- Şekil 3: `n = 10¹²` civarında sayım varyansı; `t₁⟨d⟩, t₂⟨d⟩, t₃⟨d⟩` konumlarında resurgence tepeleri. Not: **"teori verideki küçük, hızlı salınımları yakalayamıyor."** (Sonra (4.27)–(4.28) ile kısmen düzeltiliyor.) Bulgu 5'imizin ilgi alanı burası.
- (5.22): Riemann–Siegel yeniden toplama, `T* = πℏ⟨d⟩ = T_H/2` yani **τ = 1/2**'de bölünüyor. Penceremiz tamamen τ < 1/2 tarafında.
- **(5.24): Keskin kesim yerine Erfc'li yumuşak kesim** + optimizasyon parametresi K; doğruluk `exp(−t²)` (Riemann–Siegel'in `exp(−πt)`'sinden iyi). Bulgu 5 için hazır alet.
- **§6 (s. 260–263) — iki kritik ek:**
  (i) **(6.7)–(6.8):** `H = XP` dilatasyon üretir; `X → KX, P → P/K` dönüşümü "`log K` kadar zaman evrimine" karşılık gelir. **Zaman ötelemesi = ölçek değişimi.** Bu, log-frekansların doğallığının ve L-değişmezliğin operatör düzeyindeki kaynağıdır.
  (ii) **(6.9) ve izleyen paragraf:** Her Dirichlet L-fonksiyonu, `XP`'nin *farklı bir öz-eşlenik genişlemesine* karşılık gelir (çarpımsal grup mod k'nin χ(n) karakterleriyle). **Dinamik hepsinde aynıdır.** Dolayısıyla yalnız dinamiğe bağlı (karaktere bağlı olmayan) bir istatistiğin 8 ailede de aynı çıkması *beklenen* şeydir — `R_nn`'in aile-evrenselliğinin literatürdeki hazır gerekçesi budur.
  (iii) **Uyarı (madde e, s. 260):** Riemann dinamiğinde her ilkel yörünge *kendi sembolüyle* (asalıyla) etiketlenir; sonlu alfabeli sembolik dinamik yoktur. Bu, Sieber–Richter/Müller kombinatoriğinin doğrudan aktarılamamasının teknik sebebidir (bkz. §4-D).

### 3.3 E. B. Bogomolny & J. P. Keating, "Random matrix theory and the Riemann zeros I / II", *Nonlinearity* **8** (1995) 1115–1131 ve **9** (1996) 911–935
**Ne yapıyor.** Zeta'nın açık formülünü (Gutzwiller iz formülünün muadili) ve Hardy–Littlewood asal-çifti varsayımını kullanarak 3- ve 4-nokta (sonra tüm n-nokta) korelasyon fonksiyonlarının asimptotik olarak GUE'ye eşit olduğunu gösteriyor; **ve alt-mertebe (aritmetik) düzeltmeleri veriyor.**
**Bizimkiyle ilişki.** Bizim `R_nn`'imizin *var olabilmesi* için gereken tek mekanizma budur: köşegen yaklaşımında farklı asalların fazları bağımsızdır, dolayısıyla bir asalın fazına koşullamak diğerlerinin varyansını değiştiremez ⇒ `R_nn ≡ 0`. Ölçtüğümüz `R_nn = −2cos(πτ) ≠ 0`, **doğrudan K_off'un (Hardy–Littlewood asal-asal korelasyonlarının) bir ölçümüdür.** Bu, bulgunun teorik adresini kesinleştirir.
**Formül düzeyinde karşılaştırma mümkün mü?** Evet: iki-nokta sonucunun açık hali Snaith derlemesinde Teorem 4.3 (denk. 4.19) olarak veriliyor; içinde `(ζ'/ζ)'(1+ir)`, Euler çarpımı `A(η) = Π_p [(1−p^{−1−η})(1−2/p+p^{−1−η})/(1−1/p)²]` ve **`B(η) = Σ_p (log p/(p^{1+η}−1))²`** var. `B` bizim `Σ_q a_q²`'nin doğrudan akrabası.
*Not: I. makalenin yalnız özeti okunabildi (IOPscience tam metin kapalı); formüller Berry–Keating (1999) ve Snaith (2009) üzerinden doğrulandı.*

### 3.4 N. Argaman, F.-M. Dittes, E. Doron, J. P. Keating, A. Yu. Kitaev, M. Sieber, U. Smilansky, PRL **71** (1993) 4326–4329
**Ne yapıyor.** Klasik periyodik yörünge *eylemlerinin* iki-nokta korelasyonlarını inceliyor. İz formülü kesin ve spektral istatistik RMT ise, eylemler arasında **önemsiz-olmayan, evrensel biçimli korelasyonlar** bulunmak zorundadır. Bunu asal sayıların çift korelasyonu problemiyle örnekliyor ve üç kaotik sistemde sayısal doğrulama veriyor.
**Bizimkiyle ilişki — YÖN SORUSUNA CEVAP: TERS YÖN.**
Argaman ve ark.: *kuantum evrenselliği varsay ⇒ klasik/aritmetik korelasyonları çıkar.*
Biz: *aritmetik dalgaların fazına koşullu artık varyansı ölç ⇒ gazın kendi aritmetiğine tepkisini çıkar.*
Yani aynı nesnenin (asal-asal / yörünge-yörünge korelasyonları) iki zıt çıkarım yönü. Onların sonucu bizim ölçümümüzün *nedeni*; bizimki onların varsaydığı evrenselliğin *ölçülmüş imzası*. Kavramsal olarak çelişki yok, tamamlayıcılık var.
*Not: Yalnız özet doğrulandı (PRL tam metin kapalı, arXiv kopyası yok).*

### 3.5 M. Sieber & K. Richter, *Phys. Scr.* **T90** (2001) 128; ve S. Müller, S. Heusler, P. Braun, F. Haake, A. Altland, PRL **93** (2004) 014103
**Ne yapıyor.** Sieber–Richter: küçük açıyla bir kez kendi kendini kesen uzun yörüngeler ("encounter") ve eş yörüngeleri, form faktörüne RMT ile uyuşan `τ²` katkısını verir. Müller ve ark.: bunu tüm mertebelere genelleyip `K(τ) = τ + Σ_{n≥2} K_n τⁿ` açılımını kuruyor; üniter sınıfta `K_n = 0` (n>1), ortogonal sınıfta `K_orth = 2τ + Σ_{n≥2} [(−2)^{n−1}/(n−1)] τⁿ`. Bileşenler: ergodiklik + hiperboliklik + kombinatorik (permütasyon grubu), artı HOdA toplam kuralı.
**Bizimkiyle ilişki.** Bu, "τ'nun kuvvetlerinde bir çekirdek eğrisi" üretmenin literatürdeki *tek* sistematik makinesidir; `K(τ)` eğrimiz için doğal aday alet. **Ama kritik bir uyarı içeriyor:** Müller ve ark. açıkça yazıyor — "Hecke simetrili dinamiklerdeki gibi **güçlü eylem dejenerasyonları hariç tutulmalıdır**" (ref. Bogomolny–Schmit 2004). Riemann durumunda `log p` frekansları çarpımsal yapıdan ötürü tam olarak böyle dejenerasyonlar taşır. Yani **bizim sistemimiz, evrensellik ispatının dışarıda bıraktığı sınıftadır** — bulgumuzun "hiçbir denge topluluğu üretemiyor" tespitiyle uyumlu.

### 3.6 J. H. Hannay & A. M. Ozorio de Almeida, *J. Phys. A* **17** (1984) 3429–3440
**Ne yapıyor.** Periyodik yörüngelerin faz uzayında (doğal ağırlıkla) düzgün yoğun olduğu ilkesinden hareketle klasik toplam kuralını kuruyor: `lim_{T→∞} (1/T) Σ_j A_j² δ(T − T_j) = 1` (Berry–Keating denk. 3.2). Riemann tarafında bu, `A_j = −log p/p^{m/2}` ve `ρ(T) ~ e^T/T` ile **asal sayı teoreminin ta kendisidir** (Berry–Keating (3.1)).
**Bizimkiyle ilişki.** Bulgu 4'teki `dτ/τ` kuyruk ölçüsü budur. `A_j²·ρ(T) dT = dT/T` ve `τ = T/T_H` ⇒ `dτ/τ`. **L-değişmezliğin yarısı (ölçü tarafı) burada hazır.** Diğer yarısı (yalnız-τ regülatör, `ω σ_t = 2π σ_u τ`) Berry–Keating'in "ölçekleyici dinamik / λ_p = 1" gözleminden çıkar.

### 3.7 Snaith, "Riemann zeros and random matrix theory" (derleme, 2009) + Conrey–Farmer–Zirnbauer ratios
**Ne yapıyor.** Teorem 4.3 (denk. 4.19), ratios conjecture'dan türetilen **ölçeklenmemiş** sıfırların iki-nokta korelasyonunu veriyor; içinde `(ζ'/ζ)'(1+ir)`, `A(ir)` (Euler çarpımı, denk. 4.20) ve `B(ir) = Σ_p (log p/(p^{1+η}−1))²` (denk. 4.21) var. Şekil 3: ilk 100 000 sıfırdan hesaplanan ham iki-nokta istatistiğinde, **alçak Riemann sıfırlarının konumlarında (14.13, 21.02, 25.01, 30.42, 32.93, 37.59) çukurlar**; bunlar `(ζ'/ζ)'(1+ir)` teriminden geliyor. Metnin kendi ifadesiyle: rastgele matrisler bu bilgiyi içeremez.
**Bizimkiyle ilişki.** (i) `B(η)` bizim `Σ_q a_q²` toplamımızın akrabası. (ii) Şekil 3 fenomeni, "gaz kendi spektrumunu tanır"ın literatürdeki en görünür deneysel imzası. (iii) Denk. (3.2): `N = log(T/2π)` — Keating–Snaith yoğunluk eşlemesi; `N_eff = L + c` bulgumuzun ait olduğu konvansiyon.

### 3.8 Bogomolny, Bohigas, Leboeuf, Monastra, *J. Phys. A* **39** (2006) 10743 (arXiv math/0602270)
**Ne yapıyor.** Riemann sıfırlarının **en yakın komşu aralık dağılımının** GUE'den sonlu-E sapmalarını, sonlu boyutlu üniter matrislerin sapmalarıyla açıklıyor:
`N_eff = log(E/2π)/√(12Λ)`, `Λ ≡ γ₀² + 2γ₁ + c₀ = 1.57314…`, `c₀ = Σ_p (log p)⁴ Σ_{r≥1} (r−1)r²/p^r`, ayrıca `Q = Σ_p log³p/(p−1)²`.
**Bizimkiyle ilişki — DİKKAT.** Bu `N_eff` ile Keating–Snaith'in `N = L`'si aynı şey değil (biri yoğunluğa, diğeri aralık dağılımına kalibre). Programın `N_eff = L + c` bulgusu ikincisiyle *doğrudan çelişmez* ama **aynı veriden iki farklı efektif boyut çıkarmak, hangi istatistiğin kalibre edildiğine bağlıdır**. Bu iki tanımın uzlaştırılması ayrı bir çapraz kontrol maddesi olmalı. Λ ve Q'nun asal toplamları olması ayrıca ilginç: sonlu-boyut düzeltmesinin katsayısı bile aritmetiktir.

### 3.9 Bogomolny, Georgeot, Giannoni, Schmit, "Arithmetical chaos", *Phys. Rep.* **291** (1997) 219–326
**Ne yapıyor.** Aritmetik gruplardan üretilen sabit negatif eğrilikli yüzeylerde, klasik kaosa rağmen kuantum istatistikler **Poisson**'a yakındır. Sebep: periyodik yörünge uzunluklarının **üstel dejenerasyonu**.
**Bizimkiyle ilişki.** Bu, "aritmetik yapı ⇒ RMT evrenselliğinin bozulması" tezinin en temiz örneğidir. Riemann durumu bu ailenin sınırındadır (uzunluklar `m log p`; çarpımsal dejenerasyon `log(p q) = log p + log q` üzerinden gelir). "Denge topluluğu bunu üretemez" tespitimizin kavramsal komşusu.

### 3.10 G. Chavez & A. Allawala, "Prime zeta function statistics and Riemann zero-difference repulsion", *J. Stat. Mech.* (2021) 073206 (arXiv 2102.02280)
**Ne yapıyor.** Asal zeta fonksiyonunun `P(½+it) = Σ_p p^{−½−it}` kritik doğruda asimptotik olarak normal dağıldığını, öz-kovaryans fonksiyonunun `log|ζ|` ile yakından yaklaşıldığını gösteriyor. Kritik nokta: bu kovaryans, **Riemann sıfırlarının sanal kısımlarına yaklaşık eşit ayrımlarda belirgin biçimde negatiftir**; bundan da sıfır-farklarının bu değerlerden kaçınması ("repulsion") çıkıyor.
**Bizimkiyle ilişki.** `P(½+it)`, bizim asal-dalga alanımızın (m=1 kısmının) ta kendisidir. Yani literatürde **bu alanın öz-kovaryansı zaten çalışılmış** ve negatif bir yapı taşıdığı bulunmuş. Bizim `R_nn`'imiz bu alanın *fazına koşullu artık varyansı* — bir mertebe daha yukarısı. **Ortak tema: asal-dalga alanının kendi kendisiyle korelasyonu negatiftir.** Bizim `−2` ve `K(τ) ≈ −3` işaretlerimizin aynı fizikten gelip gelmediği doğrudan test edilebilir: onların kovaryans eğrisini bizim `τ` değişkenimize çevirip `K(τ)` ile karşılaştırın. **Alet çantasına eklenecek düşük maliyetli bir madde.**
*Not: Yalnız özet ve künye okundu; tam metin okunmadı.*

---

## 4. `K(τ)` Türetimi İçin Alet Çantası

Hangi teknik, öz-korelasyon çekirdeği `K(τ)`'yu türetmekte doğrudan kullanılabilir:

**A. Bogomolny–Keating asal-tarafı (BİRİNCİ TERCİH).**
*Dayanak (doğrulanmış):* Berry–Keating (1999), s. 248 — asal logaritmaları ikili korelasyonsuz olsaydı `K_off` rastgele fazların ortalaması olur ve sıfır çıkardı; `K_off ≠ 0` olması Hardy–Littlewood asal-çift korelasyonlarının doğrudan sonucudur.
Hedef nesne bir *koşullu ikinci moment*tir; bu, dördüncü mertebeden bir asal toplamı gerektirir: `⟨ (Σ_{q'≠q} …)² · cos(ω_q m) ⟩`. Bunun köşegen kısmı `R_nn = 0` verir; sıfırdan farklı kısım tamamen `Σ_{q₁,q₂} ⟨e^{i(ω_{q₁} − ω_{q₂} ± ω_q)t}⟩` üçlü rezonanslarından gelir. Bu tam olarak Bogomolny–Keating'in Hardy–Littlewood'la beslediği köşegen-dışı makinedir; onların 3- ve 4-nokta hesabı (Nonlinearity 8 (1995) 1115) doğrudan bu üçlü/dörtlü asal rezonanslarını yönetir. **Somut ilk adım:** `log q₁ − log q₂ = ± log q` koşulunu sağlayan asal-kuvvet üçlüleri (yani `q₁ = q₂ q` çarpımsal ilişkileri) üzerinden toplamı yazmak. Riemann'da bu koşul *tam olarak* sağlanabilir (Hecke tipi dejenerasyon) — jenerik kaotik sistemde sağlanamaz. **Bu, `K(τ)`'nun neden var olduğunun ve neden denge gazlarında olmadığının aday açıklamasıdır.**

**B. Berry–Keating (4.20)/(4.23) kapalı formu.**
`R_c¹(x)`'i `x = 1/2` (bağ konumu) ve `x = 1` (ardışık aralık) civarında açmak, `cos(πτ)` ve `cos(2πτ)` katsayılarını doğrudan verir. `−2` önçarpanını ve işareti bu formülden okumaya çalışmak en ucuz testtir. (4.23)'ün `−∂²_ξ Re log ζ(1−iξ)` terimi, `K(τ)`'nun τ ile yavaş kayışını (−2.87 → −3.25) üretebilecek tek analitik yapıdır.

**C. Köşegen yaklaşımı + HOdA (taban çizgisi).**
`K_diag = |τ|` ve `Σ_j A_j² δ(T−T_j) → T` ile "hiç korelasyon yoksa ne olurdu" taban çizgisini kurun. Bu, `+1.05` adyabatik kinematik terimin literatür karşılığını sabitler ve `K(τ)`'yu artık olarak tanımlar. **Ölçülü ayrışımınızı literatür diline çevirmenin en temiz yolu budur.**

**D. Encounter açılımı (Sieber–Richter / Müller–Heusler–Braun–Haake).**
`τ`'nun kuvvetlerinde sistematik seri üretmenin tek genel makinesi. **Ama olduğu gibi uygulanamaz:** Müller ve ark. güçlü eylem dejenerasyonlu (Hecke) sistemleri açıkça dışlıyor. Buradan iki yol var: (i) kombinatoriği Riemann'ın çarpımsal dejenerasyonlarını *içerecek* şekilde yeniden saymak (dışlanan durumu geri koymak); (ii) `K_n = 0` iptalinin nasıl bozulduğunu ölçmek — bozulma miktarı doğrudan `K(τ)` olabilir. **(ii) düşük maliyetli ve yüksek bilgi getirili bir test.**

**E. Erfc-yumuşatılmış kesim (Bulgu 5 için).**
Berry–Keating (5.24): keskin kesim yerine `½Erfc{(log n − θ'(t))·√(t/(2(K² − iθ''(t))))}`. Gauss Debye–Waller ile keskin kesim arasındaki tek-parametreli (K) aile budur ve *zaten optimize edilmiştir*. **Somut adım:** kendi çekirdeğinizi bu aileye fit edip K'yı ölçün; K → ∞ keskin kesime, küçük K Gauss'a gider. Bulgu 5'i tek bir sayıya indirger.

**F. Ölçekleyici dinamik argümanı (Bulgu 4 için).**
Berry–Keating (2.14): λ_p = 1, ℏ yok. Bu, herhangi bir yarı-klasik ifadenin yalnız `τ = T/T_H` üzerinden yazılabileceğini garanti eder. `dτ/τ` ölçüsüyle birleştirince L-değişmezlik **türetilmiş** olur, gözlem olmaktan çıkar. Bu argümanı yazıya dökmek, arXiv notu için en ucuz "teorem" adayıdır.

**G. Ratios conjecture (L-fonksiyon aileleri için).**
8 ailede evrenselliği teorik zemine oturtmak isterseniz, Conrey–Snaith'in reçetesi aynı yapıyı aile-bazında verir; `A(η)`, `B(η)` analogları her aile için hesaplanabilir. `R_nn`'in aileden bağımsızlığı, bu Euler çarpımlarının aileye bağlı kısmının `R_nn`'e girmediği anlamına gelir — test edilebilir bir tahmin.

---

## 5. Kaynakça

### 5.1 Açılıp okunmuş, içeriği doğrulanmış (tam metin veya geniş bölüm)

1. **M. V. Berry & J. P. Keating**, "The Riemann Zeros and Eigenvalue Asymptotics", *SIAM Review* **41**(2) (1999) 236–266. PII S0036144598347497.
   Okunan: s. 236–259 (Bölüm 1–5 tamamı). PDF: `https://empslocal.ex.ac.uk/people/staff/mrwatkin/zeta/berry-keating1.pdf`
2. **N. C. Snaith**, "Riemann zeros and random matrix theory" (derleme, 8 Aralık 2009). Okunan: s. 1–14 (tamamı + kaynakça).
   `https://people.maths.bris.ac.uk/~mancs/papers/SnaithRiemann.pdf`
3. **S. Müller, S. Heusler, P. Braun, F. Haake, A. Altland**, "Semiclassical Foundation of Universality in Quantum Chaos", *Phys. Rev. Lett.* **93** (2004) 014103. arXiv: **nlin/0401021**. Okunan: makalenin tamamı (4 sayfa) + kaynakçası.
   `https://arxiv.org/pdf/nlin/0401021`
4. **E. Bogomolny**, "Quantum and Arithmetical Chaos" (ders notları), arXiv **nlin/0312061** (2003). Okunan: giriş + Bölüm 1 başı (bölüm yapısı ve iddiaları doğrulandı).
   `https://arxiv.org/pdf/nlin/0312061`
5. **M. V. Berry**, "Semiclassical formula for the number variance of the Riemann zeros", *Nonlinearity* **1** (1988) 399–407.
   *Orijinal makaleye erişilemedi; formülün birebir metni şu kaynaktan doğrulandı (Conjecture 1.4.1, hem evrensel hem evrensel-olmayan rejim):* arXiv **2211.14918**, "On the number variance of zeta zeros and a conjecture of Berry", `https://arxiv.org/html/2211.14918`
6. **E. Bogomolny, O. Bohigas, P. Leboeuf, A. G. Monastra**, "On the spacing distribution of the Riemann zeros: corrections to the asymptotic result", *J. Phys. A: Math. Gen.* **39** (2006) 10743–10754. arXiv **math/0602270**. Okunan: özet + Λ, Q, N_eff tanımları (ar5iv HTML).
   `https://arxiv.org/abs/math/0602270`
7. **J. B. Conrey & N. C. Snaith**, "Applications of the L-functions ratios conjectures", *Proc. London Math. Soc.* **93**(3) (2007) 594–646. arXiv **math/0509480**. Okunan: özet + künye.
8. **S. M. Nishigaki**, "Distributions of consecutive level spacings of circular unitary ensemble and their ratio: finite-size corrections and Riemann ζ zeros", arXiv **2507.10193**. Okunan: özet. *(Ardışık aralıkların ORTAK dağılımı ve oranı; CUE'de sonlu-N düzeltmesi oran için O(N⁻⁴), ortak dağılımlar için O(N⁻²); Riemann sıfırlarında sapma `(log(T/2π))⁻³`.)*
9. **F. Bornemann, P. J. Forrester, A. Mays**, "Finite size effects for spacing distributions in random matrix theory: circular ensembles and Riemann zeros", arXiv **1608.04638**. Okunan: özet.
10. **G. Chavez & A. Allawala**, "Prime zeta function statistics and Riemann zero-difference repulsion", *J. Stat. Mech.* (2021) 073206. arXiv **2102.02280** (v1: 3 Şub 2021, v5: 27 Eki 2021). Okunan: özet + künye (yazarlar ve dergi doğrulandı).

### 5.2 Künyesi ve özeti doğrulanmış, tam metni okunamamış

11. **N. Argaman, F.-M. Dittes, E. Doron, J. P. Keating, A. Yu. Kitaev, M. Sieber, U. Smilansky**, "Correlations in the actions of periodic orbits derived from quantum chaos", *Phys. Rev. Lett.* **71**(26) (1993) 4326–4329. DOI **10.1103/PhysRevLett.71.4326**. *(Özet iki bağımsız kaynaktan doğrulandı: APS ve Weizmann kurumsal deposu. Tam metin kapalı, arXiv kopyası yok.)*
12. **E. B. Bogomolny & J. P. Keating**, "Random matrix theory and the Riemann zeros I: three- and four-point correlations", *Nonlinearity* **8** (1995) 1115–1131. DOI **10.1088/0951-7715/8/6/013**. *(Özet IOPscience'tan okundu; formülleri Berry–Keating 1999 ve Snaith 2009 üzerinden dolaylı doğrulandı.)*
13. **E. B. Bogomolny & J. P. Keating**, "Random matrix theory and the Riemann zeros II: n-point correlations", *Nonlinearity* **9** (1996) 911–935. *(Künye Snaith 2009 kaynakça [8]'den doğrulandı.)*
14. **E. B. Bogomolny & J. P. Keating**, "Gutzwiller's Trace Formula and Spectral Statistics: Beyond the Diagonal Approximation", *Phys. Rev. Lett.* **77**(8) (1996) 1472–1475. *(Künye APS ve Snaith kaynakça [7]'den doğrulandı; APS tam metin 403 verdi. İçerik özeti arama sonucundan: köşegen-dışı katkılar köşegen terimlerle ilişkilendirilerek hesaplanıyor.)*
15. **M. Sieber & K. Richter**, "Correlations between periodic orbits and their rôle in spectral statistics", *Physica Scripta* **T90** (2001) 128. DOI **10.1238/Physica.Topical.090a00128**. *(Özet IOPscience'tan okundu.)*
16. **M. Sieber**, "Leading off-diagonal approximation for the spectral form factor for uniformly hyperbolic systems", *J. Phys. A* **35** (2002) L613. arXiv **nlin/0209016**. *(Künye Müller ve ark. kaynakça [11] + arama sonucundan doğrulandı.)*
17. **J. H. Hannay & A. M. Ozorio de Almeida**, "Periodic orbits and a correlation function for the semiclassical density of states", *J. Phys. A* **17** (1984) 3429–3440. DOI **10.1088/0305-4470/17/18/013**. *(Künye IOPscience + Müller kaynakça [10] + Berry–Keating kaynakça [25]'ten üç kez doğrulandı; içeriği Berry–Keating §3'te ayrıntılı aktarılıyor ve o metin okundu.)*
18. **M. V. Berry**, "Semiclassical theory of spectral rigidity", *Proc. R. Soc. London A* **400** (1985) 229. *(Künye Müller ve ark. kaynakça [9]'dan doğrulandı; "köşegen yaklaşımı" olarak Müller metninde açıkça atfediliyor.)*
19. **E. B. Bogomolny, B. Georgeot, M.-J. Giannoni, C. Schmit**, "Arithmetical chaos", *Physics Reports* **291**(5–6) (1997) 219–326. DOI **10.1016/S0370-1573(97)00016-1**. *(Künye + özet arama sonucundan; ScienceDirect 403 verdi.)*
20. **E. Bogomolny & C. Schmit**, "Multiplicities of periodic orbit lengths for non-arithmetic models", *J. Phys. A* **37** (2004) 4501–4526. *(Künye Müller ve ark. kaynakça [13] + arama sonucundan doğrulandı.)*
21. **A. M. Odlyzko**, "On the distribution of spacings between zeros of the zeta function", *Mathematics of Computation* **48** (1987) 273–308. *(Künye + içerik özeti doğrulandı; ilk 10⁵ ve 10¹²+1…10¹²+10⁵ sıfırları, ±10⁻⁸ doğrulukla.)*
22. **M. V. Berry & J. P. Keating**, "A new asymptotic representation for ζ(½+it) and quantum spectral determinants", *Proc. R. Soc. London A* **437** (1992) 151–173. *(Künye doğrulandı; içeriği Berry–Keating 1999 §5'te aktarılıyor ve o metin okundu.)*
23. **M. V. Berry & J. P. Keating**, "H = xp and the Riemann zeros", in *Supersymmetry and Trace Formulae: Chaos and Disorder*, ed. J. P. Keating, D. E. Khmelnitskii, I. V. Lerner (Plenum/Kluwer, 1999), s. 355–367. *(Künye Bristol kurumsal deposu + arama sonucundan doğrulandı. Tam metin okunmadı; H=xp spekülasyonu Berry–Keating 1999 §6'da da geçiyor ve o metin okundu.)*
24. **Y. Y. Atas, E. Bogomolny, O. Giraud, G. Roux**, "Distribution of the Ratio of Consecutive Level Spacings in Random Matrix Ensembles", *Phys. Rev. Lett.* **110** (2013) 084101. arXiv **1212.5611**. *(Künye + özet doğrulandı; PDF sunucusu erişilemedi. Riemann zeta sıfırları örneği içeriyor.)*

### 5.3 DOĞRULANAMADI — kullanmayın

- **"On Statistics of the Riemann Zeros Differences"**, arXiv **1402.0865**: yalnız arama sonucunda göründü, açılmadı.
- **Odlyzko, "The 10²⁰-th zero of the Riemann zeta function and 175 million of its neighbors"**: Berry–Keating ve Snaith'te atıf var, orijinal belge açılmadı. (Snaith Şekil 1–3 bu veriyi kullanıyor.)
- **Keating, "Periodic orbit resummation and the quantization of chaos"** (Proc. R. Soc. A, 1992 civarı): aranan başlık **bulunamadı/doğrulanamadı**. Resummation içeriği Berry–Keating 1999 §5 üzerinden karşılandı; ayrı bir Keating makalesine atıf yapılmamalı.
- **Berry–Keating 1999 denk. (5.24)'ün birincil kaynağı**: SIAM metninde "[55], [4]'teki fikri genelleştirerek" deniyor; [55] ve [4]'ün tam künyeleri okunan sayfalarda görünmedi. **(5.24)'ü Berry–Keating SIAM Rev. 41 (1999) denk. (5.24) olarak atfedin**, başka bir makaleye değil.
- **GUE'de ardışık aralıkların kovaryansı / korelasyon katsayısının sayısal taban değeri**: aranan spesifik sayı (ör. ρ ≈ −0.27) **hiçbir kaynakta doğrulanamadı.** Modern literatür bunun yerine *oran* istatistiğini (Atas ve ark.) kullanıyor. Rapor bu sayıyı iddia etmiyor; kendi taban değerinizi kendiniz üretmeniz gerekiyor.

---

## 6. Sıradaki Somut Adımlar (öneri)

1. **En ucuz test:** Berry–Keating (4.20)'yi `x = 1/2`'de açıp `cos(πτ)` katsayısını çıkarın; `−2` ile karşılaştırın. Bir öğleden sonralık iş, ve tutarsa bulgunun yarısı türetilmiş olur.
2. **İkinci test:** Bulgu 1'in varyansını Berry 1988 (δ=1) ile sayısal olarak karşılaştırın — §3.1'deki cebir doğruysa artık sıfır olmalı. Bu, veri hattının bağımsız bir doğrulamasıdır.
3. **Çelişki kontrolü:** `N_eff = L + c` bulgunuzu Bogomolny–Bohigas–Leboeuf–Monastra `N_eff = L/√(12Λ)` ile aynı istatistik üzerinde karşılaştırın. Hangi konvansiyonda olduğunuzu netleştirin.
4. **Bulgu 5'i tek sayıya indirin:** Berry–Keating (5.24)'ün Erfc ailesine fit edip K parametresini ölçün.
5. **arXiv notu için en güçlü çerçeve:** "R_nn ≠ 0, köşegen-dışı asal-asal korelasyonlarının doğrudan, faz-çözünürlüklü bir ölçümüdür" — çünkü köşegen yaklaşımı zorunlu olarak R_nn = 0 verir. Bu tek cümle, bulguyu 30 yıllık K_off literatürünün tam ortasına yerleştirir.

---

## 7. EK — Yenilik Denetimi (ikinci tekne, 29 Ağu)

Beş iddiamız için hedefli tarama; hükümler:

1. **Tek asalın fazına koşullu aralık-varyansı (R_nn'in tepki okuması): YOK.**
   Hiçbir çalışmada bir spektral varyans tek bir asal dalgasının fazına
   koşullanmamış. En yakın komşu Chavez–Allawala (2102.02280): koşullama
   TERS yönde (sıfırda olmaya koşullu asal-alan istatistiği).
2. **Öz-yanıt vs dış-sürücü asimetrisi:** dış/adyabatik yarı AKRABA
   (Forrester lineer yanıt + mükemmel perdeleme, cond-mat/9411019;
   Stillinger–Lovett kuralları). **NEGATİF ÖZ-YANIT: YOK** — log-gaz/RMT/zeta
   literatüründe hiçbir iz yok. En derin bulgumuz açık arazi.
3. **τ değişkeni ve τ-yalnız istatistik: VAR.** Bohigas–Leboeuf–Sánchez
   (nlin/0012049, Found. Phys. 31 (2001) 489), denk. (32)/(44) doğrulandı:
   C(n) = (2/π²)Σ_{p,r} [sin²(πrτ_p)/(r²p^r)] cos(2πnrτ_p), τ_p = log p/L.
   Bu, dalga resminin KÖŞEGEN lag-n öz-kovaryansının ta kendisi (bizim
   2a_q²sin²(πτ)cos(2πnτ) ile birebir — denetimde doğruladım).
   "dτ/τ ölçüsü ⇒ L-değişmezlik" adlandırması ise literatürde bulunamadı.
4. **Koherans/DW-sönüm formalizmi: VAR** (Gauss yumuşatma = yörünge
   genliğinde Gauss sönümü; erfc keskin kesimden iyi — BK 1992).
5. **Nokta-bazlı özdeşlik (ds_n = dalgaların toplamı, fazlı hali): YOK —
   yazılmamış.** Varyans düzeyi tamamen VAR (Berry 1988 + BLS).

**EN TEHLİKELİ ÖNCEL: BLS nlin/0012049** — τ_p, sin²(πτ) asal toplamı ve
Odlyzko verisiyle doğrulanmış asal-asal ayrıştırma zaten orada ("~12 asalla
yeniden üretiliyor"). Madde 3/5'te yenilik iddiasından ÖNCE tam metin okunmalı.

**NET SINIR (arXiv çerçevesi için):** ikinci-moment/köşegen düzeyinin tamamı
bilinen (Berry 88 = aynı-nokta; BLS 01 = lag-n). Bizim açık arazimiz:
(a) faz-çözünürlüklü KOŞULLU istatistikler (köşegen yaklaşımında özdeş
sıfır olan her şey: R_nn, ayrışım, K(τ)); (b) negatif öz-yanıt;
(c) örnek-yolu koheransı (0.975) ve fazlı özdeşliğin kendisi.

Doğrulanamayanlar (kullanma): Caselle–Magnea adyabatik pasajı; BK 1992
erfc cümlesinin birincil metni; Blaschke–Brack denklemi (yalnız arama
metni); Berry 1988'in kendi türetiminde dτ/τ olup olmadığı.
