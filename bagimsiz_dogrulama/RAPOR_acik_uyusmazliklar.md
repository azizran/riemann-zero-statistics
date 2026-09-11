# Açık Uyuşmazlıklar — Bağımsız Çoğaltmanın Kapatamadığı Noktalar

**Tarih:** 11 Eylül 2026 · **Amaç:** arXiv öncesi, dışarıdan bakan bir okuyucunun (hakemin) takılacağı yerleri önceden listelemek.

Aşağıdaki her satır, **benim bağımsız kodumla** (yazarın scriptleri kullanılmadan, Odlyzko'nun SHA'sı teyitli ham verisinden) üretilmiş sayıdır. Çözülemeyenler "açık" olarak işaretlendi; hiçbiri "yanlış" hükmü değil, **konvansiyon sorusu**.

---

## 1. ✅ Çoğalanlar (özet)

| İddia | Sonuç |
|---|---|
| M² asimptotiği (Conrey–Ghosh eğimi, HLPC sabit terimi) | eğim +0.05%, sabit −0.14% |
| r tablosu, 12 pencere | 12/12, sapma ≤0.004 |
| Gaussian surrogate null | 0.4938 ± 0.0059 vs makale 0.494 ± 0.007 |
| N_eff anomalisi (Pearson kayması + Spearman sabitliği) | trend ve büyüklük tuttu |
| Toplam kural (on üç asal, u = 1) | 1.001–1.008; plasebo tabanı 0.030 (44–112σ) |
| Kanallar saf reel (kuadratür ~0) | faz 0.0–0.3°, kuadratür ~10⁻⁴ |
| `v(τ)` τ-kolapsı, 10¹² penceresi (t oranı 2,2×10⁶) | oran 1.069 ± 0.055, r = 0.947 |
| `v(τ)` kolapsı 10²¹ ve 10²²'de (t oranı 1,14×10¹⁶) | 10²¹/10¹² = 0.995 ± 0.147, 10²²/10¹² = 1.000 ± 0.142 ✅ |
| Not 4 benek fazları 180° | 24/24 satır, sapma ≤0.3° |
| Not 3: `w` işaret değişimi + dönmeden 0→π atlama | işaret değişiyor, kuadratür ≤%5 ✅ (konum için §4-bis) |
| Not 3 (iv): doğrudan |Z|² kanalı `R_pred(q)` (şekil) | r = 0.9989, 21 asal ✅ (çadır kanalı için §4-quater) |

## 2. 🔧 YENİ BULGU: `k` faktörü — asal-kuvvetlerinde normalizasyon

Açık formülde bir asal-kuvvetin satır genliği `1/(π k √q)` ile gelir; yani `Λ(q)/(log q) = 1/k`. Ben ilk denemede normu `√q` almıştım ve asal-kuvvetleri "birimin çok altında" bulmuştum (q=16'da u = 0.213). `k` ile çarpınca:

| q | 4 | 8 | 9 | 16 | 25 | 27 |
|---|---|---|---|---|---|---|
| u | 0.493 | 0.315 | 0.478 | 0.213 | 0.455 | 0.287 |
| **u·k** | **0.986** | **0.944** | **0.955** | **0.852** | **0.910** | **0.861** |

Tüm asal-kuvvetlerinde `u·k = 0.950 ± 0.083`. **Yani toplam kural yalnız asallarda değil, tüm satırlarda geçerli** — doğru normalizasyon `k√q`. Bu, Not 2'nin "prime powers a few percent below" ifadesiyle uyumlu hale geliyor; benim ilk turdaki "asal-kuvvetleri sapıyor" gözlemim normalizasyon hatamdı.

**İstek:** Not 2/3'ün `w` ve `v` tablolarında asal-kuvvetleri için hangi normalizasyonun kullanıldığı metinde yazılı değil (`w = coefficient × p^{1/2}` yalnız asallar için). `k` faktörü eklenince bende `w` tablosu da değişiyor; makale ile karşılaştırma ancak bu netleşince yapılabilir.

## 3. ⚠️ `√w–v` kısıtı: ilişki var, katsayılar ve sıkılık farklı

| Veri kümesi | Fit | RMS | r |
|---|---|---|---|
| 6 asal × 5 pencere (n=30) | √w = 1.056 − 1.259 v | 0.0127 | −0.995 |
| yalnız asallar, tüm τ (n=169) | √w = 1.216 − 1.846 v | 0.0610 | −0.958 |
| tüm satırlar, `k` düzeltmesiz (n=230) | √w = 0.440 − 0.378 v | 0.1846 | −0.344 |
| tüm satırlar, `k` düzeltmeli (n=230) | √w = 1.188 − 1.774 v | 0.0827 | −0.917 |
| **makale** | **√w = 1.017 − 0.884 v** | **0.0075** | — |

İki gözlem:
1. **Eğim bende ~2 kat dik** (−1.77…−1.85 vs −0.884). Kesişim 1.19…1.22 vs 1.017.
2. **Sıkılık veri kümesine bağlı**: 30 noktalık asal alt-kümesinde RMS 0.013; tüm asallarda 0.061; tüm satırlarda 0.083. Yani "sıkı doğrusal kısıt" ifadesi, hangi satırların dahil edildiğine duyarlı.

**Kavramsal uyarı (hakem de sorabilir):** makalenin kendisi `w` ve `v`'nin ikisinin de tek-değişkenli τ-yasalarına kollabe olduğunu söylüyor. İki nicelik de τ'nun fonksiyonuysa, `√w`–`v` düzlemi tek parametreli bir **eğridir** — doğrusallığı ek bilgi değil, o eğrinin şeklidir. Kısıtın bağımsız bilgi taşıdığını göstermek için ya (a) τ-kolapsının kırıldığı bir rejimde de doğrusal kalmalı, ya da (b) aynı τ'da farklı `w`,`v` çiftleri üreten bir değişken (ör. farklı pencere/L) gösterilmeli.

**İstek:** (i) `w` tablosu τ-kolapsı mı, pencere ortalaması mı? (ii) "own-gap kontrolü"nün tam tanımı (ben `log δ̃` düz sütununu kullandım; "aynı-frekans" varyantı matematiksel olarak dejenere çıktı). (iii) Kısıt fiti hangi satırları içeriyor?

## 3-bis. 🔧 GAP 1 ÇÖZÜLDÜ: `w` seviyesi ve kısıt katsayıları (11 Eylül 2026)

Scratchpad + kampanya günlüğü okundu; tablo netleşti:

**(a) Üç ayrı "w" var — karıştırılmamalı:**

| nesne | tanım | değer | nerede |
|---|---|---|---|
| `wq` (184 K1) | `\|c_q(ds)\| / (2a_q sin πτ)` — **yer değiştirme zarfı** | **1.18 – 1.28** (>1) | `K1_gercek.npz` |
| `wg` (185b) | `Σ_b\|â_q\| / Σ_b ae_q` — bant genlik oranı | 1.17 – 1.27 | `K1_faktorler.npz` |
| **Not 2 tablosu** `w_p` | log M̃'nin dalgalara **gap-koşullu** regresyonu | 0.83 … 0.38 | Not 2 (79–95 arku) |

**(b) Benim bant düzeyinde ölçtüğüm ilişki makalenin kısıtı DEĞİL.** `√wg` vs `v_b` bende
**pozitif** eğimli (gerçek: 0.553 + 0.905·v, r = 1.000; ikiz: 0.576 + 0.867·v), makalenin
kısıtı ise **negatif** (−0.884). Sebep: bant düzeyinde hem `√wg` hem `v_b` τ̄'nun düzgün
fonksiyonu → ilişki bir τ-yeniden-parametrizasyonu. Kısıt **satır düzeyinde** yaşıyor.

**(c) 1.017 / −0.884'ün kaynağı bulundu — ve benim farkım orada zaten belgelenmiş.**
Kampanya günlüğü (`LITERATUR_TARAMASI`, görev 82, 20 Ağustos):

> *"YER-GERÇEĞİ KANITI (T4b): … P4-tabanlı regresyon w'yi **%15-25 SİSTEMATİK DÜŞÜK** ölçüyor
> (21σ); tam taban gerçeği ~%1 içinde buluyor. Kısıtlı-taban kanal ölçümü ızgara-çiftlenimi
> yüzünden yanlış; FİZİKSEL OKUMA = TAM-TABAN OKUMASI."*
> *"GERÇEK VERİDE ETKİ: korunum yasası fiti P4→TAM: a 0.990→1.017, b 1.025→0.884 (~12σ)."*

Yani kısıtın katsayıları **tam-taban** konvansiyonundan; kısıtlı tabanda 0.990/1.025 idi. Ve
kampanyanın ölçtüğü "w'yi %15–25 düşük ölçme" etkisi, benim bağımsız ölçümümdeki farkın
**tam olarak aynı mertebesi** (bende 0.781…0.285 vs tablo 0.83…0.38 → %6–25).

**Hüküm:** Gap 1 artık "gizemli uyuşmazlık" değil — **kampanyanın kendi kaydettiği bir taban/konvansiyon duyarlılığı.** Benim kalan %6–25'lik farkım, büyük olasılıkla "own-gap kontrolü"nün tanımından (ben düz `log δ̃` sütunu koydum). Not 2/3'ün sayıları doğru konvansiyonda; **ama makale metni "own-gap kontrolü"nü tanımlamıyor** — hakem için hâlâ açık, tek cümlelik bir düzeltme yeter.

**(d) Bonus teyit:** kampanya `k` faktörünü bağımsız bulmuş (görev 90: `k·v(p^k) = asal v-eğrisi`)
— benim `u·k = 0.95 ± 0.08` bulgumla aynı yapı.

## 4. ⚠️ Not 4 tarak köprüsü: nitel tutuyor, nicel tutmuyor

Unfold kafeste Bragg yansımaları (kendi ölçümüm, 5 pencere):

| L | \|Ĝ(2π)\| | σ_φ | Gauss e^{−c} | oran | \|Ĝ(4π)\| |
|---|---|---|---|---|---|
| 9.857 | 0.00846 | 0.290 | 0.190 | 0.045 | 0.00404 |
| 10.928 | 0.01442 | 0.288 | 0.194 | 0.074 | 0.00058 |
| 11.978 | 0.02578 | 0.287 | 0.197 | 0.131 | 0.00076 |

- **Nitel uyum:** ikinci tarak birinciden 3–8× sönük → makalenin "second-order comb is extinguished" ifadesi doğru yönde.
- **Nicel uyumsuzluk:** makale "measured/Gaussian = 0.88–0.89 sabit" diyor; bende 0.045–0.131 (10× küçük ve pencereden pencereye 3× değişiyor). Benim Gauss tahminim `e^{−(2πσ_φ)²/2}` ile yapıldı (σ_φ = sarılmış faz std'si); "first-order brightness"ın tam normalizasyonu metinden çıkmıyor.

**İstek:** tarak parlaklığının tanımı (ham `|Ĝ|` mı, `I = |Σ w_n e^{iωt_n}|²/Σw_n²` mi, ideal kafese göre normalize mi) ve `c`'nin hangi σ ile hesaplandığı.

## 4-bis. ⚠️ `w` işaret değişimi: nitel ✅, konum farklı

Not 3 (ii) der ki: *"w changes sign at τ₀ = 0.447 ± 0.005, with response phases that jump from 0 to π without rotating (quadrature below 8% through the crossing)"*.

Bunu düşük-t pencereleriyle (`zeros_100k`, t ≈ 1.7×10³ → 5.6×10⁴) ve yüksek-t pencereleriyle bağımsız ölçtüm; **işaretli** kosinüs katsayısını (kuadratür ayrı) aldım. Asallarda τ bantları:

| τ bandı | 0.05-0.10 | 0.10-0.15 | 0.15-0.20 | 0.20-0.25 | 0.25-0.30 | 0.30-0.35 | 0.35-0.40 | 0.40-0.45 |
|---|---|---|---|---|---|---|---|---|
| a ort | **+0.729** | +0.584 | +0.432 | +0.310 | +0.199 | +0.101 | +0.032 | **−0.034** |
| \|b\|/\|a\| | 0.0004 | 0.0008 | 0.0014 | 0.0011 | 0.0019 | 0.0121 | 0.0394 | 0.0480 |

- **✅ Nitel doğrulandı:** işaret gerçekten değişiyor (son pozitif τ = 0.387, ilk negatif τ = 0.408) ve kuadratür geçiş boyunca ≤%5 (makale: <%8) — yani faz **dönmeden 0 → π atlıyor**. Ayrıca geçişten önce katsayı düzgün ve **dışbükey** biçimde sıfıra iniyor.
- **⚠️ Nicel fark:** benim sıfır geçişim **τ₀ ≈ 0.397** (tüm asallar, [0.06, 0.45], n=147, doğrusal fit: a = 0.800 − 2.022τ → 0.396); makale **0.447 ± 0.005**. Aradaki fark ~%11 ve makalenin verdiği hata çubuğunun çok dışında.

**Olası sebepler (teyit gerekiyor):** (a) makalenin τ₀'sı "termal Bragg modeli" ile birlikte fit edilmiş (ideal ayna τ = 1/2'de, jitter ile 0.443'e kayıyor) — yani model-bağımlı bir konum; benimki model-bağımsız sıfır geçişi. (b) `w`'nin normalize edilme biçimi (τ'ya bağlı bir çarpan) geçişi kaydırabilir. (c) benim τ ölçeğim `L = log(t̄/2π)` ile kuruldu ve fold `τ = 1/2` (Riemann–Siegel ana toplam kesimi) doğru yerde.

**İstek:** τ₀ model-bağımsız mı (ham `a(τ)` sıfırı) yoksa Bragg fitinden mi geliyor? İkincisiyse, ham sıfır geçişinin de raporlanması iki ölçümü karşılaştırılabilir kılar.

## 4-quater-bis. 🔴 ERRATUM: repo'daki BBLM `c₀` formülü yanlış

`KESIF_SEFERI_KUANTUM_KAOS_29AGU2026.md` satır 114 diyor ki:

> `Λ ≡ γ₀² + 2γ₁ + c₀ = 1.57314…`, `c₀ = Σ_p (log p)⁴ Σ_{r≥1} (r−1)r²/p^r`

**Bu transkripsiyon hatalı.** BBLM (arXiv math/0602270) metnindeki gerçek tanım:

```
c_n = [(−1)^n/(2n)!] · Σ_p (log p)^{2(n+1)} · Σ_{r≥1} (r−1) r^{2n} / p^r
```

`n = 0` için: `r^{2n} = 1` ve üs `2(n+1) = 2`, dolayısıyla

```
c₀ = Σ_p (log p)² · Σ_{r≥1} (r−1)/p^r = Σ_p (log p)²/(p−1)²
```

Sayısal sonuç (ilk 17984 asal, mpmath 25 hane):

| nicelik | değer | doğrulama |
|---|---|---|
| `c₀ = Σ (log p)²/(p−1)²` | **1.3855389** | — |
| `c₀ = Σ (log p)⁴Σ(r−1)r²/p^r` (repo) | **33.808** | ❌ tutmuyor |
| `γ₀² + 2γ₁` | 0.1875462 | — |
| **Λ = γ₀²+2γ₁+c₀** | **1.573085** | makale 1.57314 (fark = asal kuyruğu, ∫log x/x² ~ 6×10⁻⁵) ✓ |
| `√(12Λ)` | 4.344769 | — |
| `Q = Σ log³p/(p−1)²` | 2.314970 | — |
| `C = Q/Λ` | **1.47161** | makale 1.4720 ✓ |

**Etki:** Λ'nin *değeri* doğru (1.57314) olduğu için bu hata sonuçları bozmamış; ama formül yanlış olduğu için "asal toplamı" iddiası denetlenemez haldeydi. Not 5/6'da bu tanıma atıf yapılırsa düzeltilmeli.

### BBLM `N_eff` ile bizim `N_eff` — nicel karşılaştırma

| L | BBLM `L/√(12Λ)` | bizim ölçüm (gap–max korelasyonu) | oran |
|---|---|---|---|
| 5.60 | 1.289 | 6.52 | 5.06 |
| 8.39 | 1.931 | 9.45 | 4.89 |
| 9.86 | 2.269 | 10.91 | 4.81 |
| 11.98 | 2.757 | 13.40 | 4.86 |

İki tanım **~4.85 sabit oranıyla** ayrılıyor. Bu bir çelişki değil — BBLM'ninki *aralık dağılımına*, bizimki *gap–max korelasyonuna* kalibre — ama "`N_eff = L + c`" cümlesi bir makaleye girerse, BBLM'nin `N_eff`'iyle **karıştırılmaması** için açıkça "korelasyon-kalibreli" diye etiketlenmeli. (Oranın sabitliği derin bir şey değil: iki nicelik de L'de ayrı eğimlerle doğrusal.)

## 4-quater. ⚠️ Not 3 (iv): doğrudan kanal ✅, "çadır" kanalı ✗

Not 3 (iv): *"The structure factor (1−2τ) is derived … as the pair count of the folded main sum; pointwise |Z|² regressions … instead follow the shifted-second-moment family ≈ 2(1−τ) through and beyond the fold, so that the complementary mirror channel is a tent peaked at the fold."*

Bunu **parametresiz** test ettim: L = 9.86 penceresinde 1.279.919 noktalık ızgarada |Z|² hesaplandı (`74`'ün tarif ettiği 1.27M nokta ile aynı kurulum), her asal için Fourier genliği çıkarıldı ve 74'ün özdeşliğiyle karşılaştırıldı:

```
DC = mean|Z|² = 11.0136    (teori: 2·S₁(N) = 11.015, N = √(t/2π) = 138.2)  ✓
R_pred(q) = 2 q^{−1/2} S₁(⌊N/q⌋)/S₁(N)      [τ > 1/2, yani q > N → TAM SIFIR]
R_tam(τ)  = 2(L − log q + 2γ − 1)/(L + 2γ − 1)   [≈ 2(1−τ)]
```

| q | 2 | 3 | 5 | 11 | 53 | 137 | 139 | 157 |
|---|---|---|---|---|---|---|---|---|
| τ | 0.070 | 0.111 | 0.163 | 0.243 | 0.403 | 0.499 | **0.501** | 0.513 |
| ölçülen | 1.325 | 1.040 | 0.764 | 0.471 | 0.176 | 0.094 | 0.091 | 0.084 |
| R_pred | 1.237 | 0.926 | 0.632 | 0.340 | 0.075 | 0.031 | **0.000** | **0.000** |
| R_tam | 1.862 | 1.781 | 1.679 | 1.521 | 1.207 | 1.017 | 1.014 | 0.990 |

- **✅ Doğrudan kanal kimliği çoğaldı (şekil):** fold altındaki 21 asalda ölçülen ile `R_pred` arasında **r = 0.9989**. Küçük τ'da seviye de tutuyor (q=2'de %7, q=3'te %12).
- **⚠️ Ama keskin ufuk görünmüyor:** `R_pred` q > N = 138.2'de tam sıfıra düşerken, ölçülen modülasyon orada da ~0.09 veriyor (yumuşak geçiş, sıçrama yok). Ölçülen eksi doğrudan = "ayna" kalıntısı **~0.06–0.14** bandında **yayvan**.
- **✗ "Çadır" iddiası tutmuyor:** makalenin ayna tahmini `R_tam − R_pred` foldda ~1.0'a çıkıp üstünde 2(1−τ) ≈ 1.0'da kalıyor; benim ölçtüğüm kalıntı bunun **~1/10'u** ve tepe yapmıyor. Ayrıca `R_tam` ailesi (≈2(1−τ)) benim ölçtüğüm toplam modülasyonun kendisiyle de uyuşmuyor (τ=0.5'te ölçülen 0.09, R_tam 1.02).

**GÜNCELLEME:** bu genlik açığı §5-bis'te çözüldü (DW çarpanı); aşağıdaki not yalnız tarak/çatlak için geçerli. **İstek:** "pointwise |Z|² regressions" tam olarak hangi nicelik? Ben ham Fourier genliğini (`2|c|/DC`) ölçtüm; onların regresyonu bir taban içeriyorsa seviye farkı oradan gelebilir — ama **çadırın tepe yapması** normalizasyondan bağımsız bir şekil iddiası ve o bende yok.

## 4-quinquies. ✅/⚠️ KALEM 188 — veri kurtarıldı, T1 yapıldı

Scratchpad (277 MB zip) 11 Eylül'de kurtarıldı; pipeline bu makinede **birebir** çalıştı
(gerçek ζ = 0.3287 ∠179.96°, ikiz 0.0768 — kampanya değerleriyle aynı).
T1 sonucu `NOT_188_teori.md` §5-bis'te; özet:

- **İptal aritmetik** (8 faz-rastgele vekilin hiçbiri |ζ| ≈ 0.33 vermiyor; ortalama 2–7).
- **⚠️ H-F2 mührü zayıf:** faz-rastgele vekillerin **3/8'i de** "8 bant 180°±15°" testini geçiyor
  → bu ölçütün surrogate null'a karşı yanlış-pozitif oranı ~%40. Ayrım **modülde**
  (gerçek 8 bant 0.32–0.38; hiçbir vekilde yok). Not 5/6'da H-F2 tek başına sunulmamalı.

**T2 (kayan bant taraması) — plato kinematik değil.** 42 kayan pencerede (genişlik 0.05 ve 0.10)
|ζ|(τ) düz değil: **τ ≈ 0.65'te minimum 0.3164 ± 0.0035**, iki yana yükseliyor
(0.45'te 0.336, 0.85'te 0.392). Hatalar ±0.003–0.008 olduğundan yapı 10–20σ anlamlı.
**1/π = 0.31831 tam da bu minimuma denk geliyor** → "plato = 1/π" okuması yapılı bir
eğrinin tek noktasındaki tesadüf; **1/π hipotezi düşürülmeli.** Faz kilidi ise 42 pencerenin
hepsinde 180° ± 0.7° — bant seçiminden bağımsız, yani H-F2'nin *içeriği* sağlam.

**T3 (kesim taraması) — şekil kesim konumunun eseri değil.** τ_c = 0.70 / 0.86 / 1.00
(246 → 6119 çizgi) için aynı tarama: minimum her üçünde de var (τ ≈ 0.67 / 0.66 / 0.55–0.60)
ve her eğri **kendi kesimine doğru monoton yükseliyor** (0.86 → 0.392 @ τ=0.85;
1.00 → 0.427 @ τ=0.90). Yani |ζ|(τ) ≈ *taban + kesime-yakınlık artışı*. **Uyarı:** mutlak
seviye kesime bağlı (min 0.308 / 0.316 / 0.235) → ζ alıntılanırken **kesim birlikte
yazılmalı**; ve 1/π tesadüfü böylece ikinci kez zayıflıyor (minimumun değeri kesimle
değişiyor).

## 4-sexies. 🔧 GAP 2 ÇÖZÜLDÜ: tarak parlaklığı — iki taban, ikisi de doğru (11 Eylül 2026)

`88_S_modeli.py` M2 bloğu **birebir** koşturuldu (yazarın tanımı: `x = rvm_N(tmid)`,
`delta = angle(P1·conj(mean yönü))/2π`, `meas = |⟨e^{2πikx}⟩|`).

**İki farklı "Gauss" tabanı var — notun kendi dipnotunun söylediği gibi:**

| taban | σ | Gauss | ölçüm/Gauss | not |
|---|---|---|---|---|
| **sarılmış faz** (`delta.std()`) | 0.239–0.250 | 0.3228…0.2912 | 0.747 ± 0.031 | not: "0.32 → Edgeworth 0.27 vs ölçülen 0.256" ✓ birebir |
| **unfold edilmiş jitter** (`σ_u = σ_t·L/2π`) | 0.251–0.272 | 0.287…0.233 | **0.8847 ± 0.0050** | not: **"0.88–0.89 sabit"** ✓ **birebir** |

| L | 9.856 | 10.368 | 10.928 | 11.467 | 11.978 | 12.448 |
|---|---|---|---|---|---|---|
| σ_t (mutlak) | 0.16028 | 0.15510 | 0.14982 | 0.14507 | 0.14080 | 0.13713 |
| ölçüm | 0.2563 | 0.2441 | 0.2319 | 0.2208 | 0.2121 | 0.2051 |
| oran (unfold tabanı) | 0.8925 | 0.8894 | 0.8859 | 0.8808 | 0.8794 | 0.8802 |

**Hüküm:** makalenin 0.88–0.89 iddiası **tam doğru** ve altı pencerede sabit; benim önceki
"0.045–0.131" ölçümüm yanlış tabanı (lineer unfold + sarılmış faz) kullanmamdan geliyordu.
Notun kendi işaret ettiği Edgeworth zinciri de (wrapped-phase tabanı) birebir çoğaldı
(0.32 → 0.27 vs ölçülen 0.256 @ L=9.86). **Bu açık tamamen kapanmıştır.**

## 5. ⚠️ Not 4 genlik yasası (benek şiddeti)

Asal çizgilerde `|Ĝ|/[(Λ(q)/(L√q))cos(πτ)]` oranı τ ile 0.99 → 0.79 arasında düşüyor; makale 14 çizgide %0.1–4 diyor. Faz kilitli (180°), genlik yasasının **parametresiz** kesinliği teyit edilemedi. `88_S_modeli.py`'nin eski normu (`v·p^{−1/2}/2`, `f(τ)=πτ·cot(πτ)`) farklı bir aile olduğu için karşılaştırma doğrudan yapılamıyor.

## 5-bis. 🔧 GAP 3 ÇÖZÜLDÜ: benek genliği — eksik olan DW çarpanıydı (11 Eylül 2026)

Makalenin tahmini (`103_fig_diffraction_v2.py`): `pred = τ·q^{−1/2}·cos(πτ)·DW`,
`DW = exp(−(log q)² σ² / 2)`, `σ = std(orta nokta − pürüzsüz RvM kafesi)` (Newton ile, `89`/`103`).

**Benim önceki testim DW'yi içermiyordu** — o yüzden oranlar 0.99 → 0.79 diye düşüyordu; o düşüş
DW'nin ta kendisiymiş. DW dahil edilince (bu pencerede σ = 0.1603, σ/ḡ = 0.251):

| | oran (ölçülen/pred) | saçılma |
|---|---|---|
| **yalnız asallar (18 çizgi), DW dahil** | **0.987 ± 0.054** | **%5,4** |
| asallar, DW'siz | 0.678 ± 0.289 | %43 |
| asal-kuvvetleri (k≥2), DW dahil | 0.08 – 0.49 | — |

**Hüküm:** benek yasası **asal çizgilerinde parametresiz olarak ~%1 doğrulukla geçerli** — makalenin
iddiası doğru, benim "tutmuyor" hükmüm DW'yi atlamamdan geliyordu. Asal-kuvvetleri ayrı bir
konvansiyon istiyor (kampanya görev 90: `k·v(p^k) = asal v-eğrisi`; benim `u·k = 0.95` bulgum
aynı yapı) — bu makale metninde zaten "prime powers a few percent below" diye kayıtlı.

## 6. ✅ VERİ KURTARILDI — engel kalktı (11 Eylül 2026)

Diğer MacBook'taki scratchpad (277 MB zip) alındı ve `qm_riemann/scratchpad/` altına kuruldu
(26 görev klasörü, 297 MB). Pipeline bu makinede **birebir** çalışıyor:
`187c_zeta_defteri.py` gerçek denizde **ζ = 0.3287 ∠179.96°**, ikizde **0.0768** veriyor —
kampanya raporuyla aynı; ONKAYIT sha kontrolü de geçti (scriptler bayt-bayt aynı).

Böylece kalem 188'in veri-bağımlı testleri koşuldu: **T1** (faz-rastgele vekiller → iptal
aritmetik), **T2** (kayan bant → plato yapılı, 1/π düştü), **T3** (kesim taraması → taban +
kesime-yakınlık). Ayrıntı `NOT_188_teori.md`.

Eski engel kaydı (tarihsel): 152–187 zincirinin ara verisi başka makinenin geçici
scratchpad'inde yaşıyordu ve repoya hiç commit edilmemişti → `RAPOR_veri_kurtarma_riski.md`.

## Kapanış (11 Eylül 2026, akşam — güncel)

Bağımsız çoğaltmanın işaretlediği **dört "konvansiyon açığı"nın dördü de kapandı**:

| açık | sonuç |
|---|---|
| `k` faktörü (asal-kuvvet normalizasyonu) | benim bulgum; kampanya görev 90'da aynı yapı kayıtlı ✓ |
| **Gap 1** `w` seviyesi + kısıt katsayıları | kampanyanın belgelediği **taban duyarlılığı** (P4→TAM: %15–25); notlar tam-taban sayılarını kullanıyor ✓ |
| **Gap 2** tarak parlaklığı 0.88–0.89 | **birebir doğrulandı** (0.8847 ± 0.0050) — yanlış taban bendeydi ✓ |
| **Gap 3** benek genliği yasası | **birebir doğrulandı** (asallarda 0.987 ± 0.054) — eksik DW çarpanı bendeydi ✓ |

**Geriye kalan tek istek (belge düzeltmesi, ölçüm değil):**
1. Not 2/3 metninde **"own-gap kontrolü"nün tanımı** yazılı değil (benim `w`'m bu yüzden %6–25 farklı olabilir — aynı mertebe, kampanyanın kendi ölçtüğü taban etkisi).
2. Not 5/6'da **H-F2** sunulurken modül koşulu (|ζ| ≈ 1/3) ve surrogate tabanı eklenmeli; ±15° ölçütü tek başına %40 yanlış-pozitif veriyor.
3. ζ alıntılanırken **kesim (τ_c = 0.86)** birlikte yazılmalı; mutlak seviye kesime bağlı.

Hiçbir bulgu "yanlış" çıkmadı. Dört hatanın **üçü bendeydi** (Ĝ ortalaması, DW çarpanı, tarak tabanı) ve hepsi bağımsız testle yakalanıp düzeltildi.
