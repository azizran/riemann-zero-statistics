# κ_p(L) sonlu-yükseklik sönmesi — literatür teftişi

Tarih: 25 Eylül 2026. Konu: KALEM 198'in ölçtüğü κ_p(L) < 1 (pozitif üslü/paydaki asal uyduların
BK'ya göre sönük olması, negatif üslü/ayna-yanının TAM tutması) fenomenini açıklayabilecek
literatür — Conrey–Snaith ratios sanısı, Bogomolny–Bohigas–Leboeuf–Monastra (BBLM) sonlu-E
düzeltmesi, Landau–Gonek ve genellemeleri (özellikle Durkan–Hughes–Pearce-Crump 2026).

KURAL: her formül birincil kaynaktan (arXiv PDF, çoğunlukla bu oturumda pypdf ile bizzat
metin-çıkarımı yapılarak) doğrulandı. Doğrulanamayan her şey **DOĞRULANAMADI** diye
işaretlendi. Bazı künyeler önceki oturumlarda (24 Eylül 2026, `LITERATUR_BOGOMOLNY_KEATING_
24EYL2026.md` ve `LITERATUR_UYDU_KLASIK_MI_24EYL2026.md`) tam-metinden doğrulanmıştı; bu
dosyada o doğrulamalar **aktarılıyor** (yeniden okunmadı) ve öyle belirtiliyor.

## 0. Erişim durumu / birincil kaynak tablosu

| # | Kaynak | Bu oturumda erişim |
|---|---|---|
| 1 | E. Bogomolny, "Riemann zeta function and quantum chaos", arXiv:0708.4223 | TAM METİN (yeniden okundu, `scratchpad/bogomolny_0708_4223.txt`, pypdf) |
| 2 | J.B. Conrey, N.C. Snaith, "Applications of the L-functions ratios conjectures", arXiv:math/0509480 | TAM METİN (yeniden okundu, `scratchpad/conrey_snaith_math0509480.txt`, pypdf) |
| 3 | E. Bogomolny, O. Bohigas, P. Leboeuf, A.G. Monastra, "On the spacing distribution of the Riemann zeros: corrections to the asymptotic result", arXiv:math/0602270 | TAM METİN (bu oturumda ilk kez PDF indirilip pypdf ile okundu, `scratchpad/bblm_0602270.txt`, 9 sayfa) |
| 4 | B. Durkan, C. Hughes, A. Pearce-Crump, "Generalisations of the Landau–Gonek Theorem and applications to mean values of zeta", arXiv:2601.18025 | TAM METİN (bu oturumda ilk kez pypdf ile okundu, `scratchpad/durkan_hughes_pearcecrump_2601_18025.txt`, 19 sayfa; önceki oturumda yalnız özet okunabilmişti) |
| 5 | E.B. Bogomolny, J.P. Keating, "Gutzwiller's Trace Formula and Spectral Statistics: Beyond the Diagonal Approximation", Phys. Rev. Lett. 77 (1996), 1472–1475 | ERİŞİLEMEDİ (APS paywall) — künye 3 bağımsız birincil kaynağın (Bogomolny 2007 ref. 7, BBLM ref. [7], Conrey–Snaith ref. [3]) kendi kaynakçalarından çapraz doğrulandı |
| 6 | E.B. Bogomolny, J.P. Keating, "Random matrix theory and the Riemann zeros I / II", Nonlinearity 8 (1995) 1115–1131 / 9 (1996) 911–935 | ERİŞİLEMEDİ (IOP paywall/bot-koruması) — künye aynı şekilde çapraz doğrulandı (not: Bogomolny 2007'nin kendi kaynakçası Part II'yi "Nonlinearity 9, 1995" diye yanlış yazmış; BBLM'nin kaynakçası da aynı hatayı tekrarlıyor — CrossRef'e göre doğrusu 1996'dır, bkz. önceki oturum notu) |
| 7 | "Bogomolny–Keating 2013" (görevde anılan) | **BULUNAMADI.** arXiv API'de `au:Bogomolny_E AND au:Keating_J` sorgusu 0 sonuç döndürdü (bu oturumda çalıştırıldı) — Bogomolny–Keating ikilisinin arXiv'de hiç ortak yazarlı makalesi yok (1995/96 makaleleri arXiv-öncesi). 2013 tarihli, konuyla en yakın makale **B. Rodgers, "Macroscopic pair correlation of the Riemann zeroes for smooth test functions"**, Q. J. Math. 64 (2013), 1197–1219, arXiv:1203.3275 — ama bu Rodgers'ın TEK yazarlı makalesi, Bogomolny–Keating'in değil. Aşağıda §1'de bu karışıklık ayrıca ele alınıyor. |
| 8 | E. Landau, "Über die Nullstellen der Zetafunction", Math. Ann. 71(4) (1912), 548–564 | Birincil metne erişilmedi; künye kaynak 4'ün (DHPC) kendi kaynakçasından ([19]) birebir doğrulandı — DHPC'nin metni Teoremi "Landau, 1911" diye etiketliyor ama kaynakça yılı 1912 (aşağıda not edildi) |
| 9 | S.M. Gonek, "A formula of Landau and mean values of ζ(s)", Topics in analytic number theory (Austin, Tex., 1982), Univ. Texas Press, 1985, 92–97; ve "An explicit formula of Landau and its applications to the theory of the zeta-function", Contemp. Math. 143 (1993), 395–413 | Birincil metne erişilmedi; her iki künye VE Gonek'in "Theorem 2" ifadesinin TAM formülü (hata terimleriyle) kaynak 4'ün (DHPC) TAM METNİNDEN (Bölüm 1, Teorem 1/2, Sonuç 2.1) birebir aktarıldı — bu, önceki oturumdaki ikincil (Ford–Soundararajan–Zaharescu üzerinden) doğrulamadan daha güçlü bir doğrulama |
| 10 | A. Fujii, "On a theorem of Landau" I/II (1989/1990), "On a conjecture of Shanks" (1994), "On the distribution of values of the derivative..." (2012) | Birincil metne erişilmedi; 4 künye de kaynak 4'ün (DHPC) kendi kaynakçasından ([6],[7],[8],[9]) birebir doğrulandı |

---

## SORU 1 — Conrey–Snaith (math/0509480): BK'nın ötesinde 1/log T mertebesinde, pozitif-frekanslı düzeltme var mı? BBLM (math/0602270) ve "1996/2013" ne diyor?

### 1.1 Conrey–Snaith'in kendi formülü (Teorem 4.1, tam metinden)

Kaynak 2, §4 "Pair-correlation" (metin satır ~1414–1830). Formül, çift-korelasyon toplamı
S(f) = Σ_{0<γ,γ'≤T} f(γ−γ') için, Cauchy integral yöntemiyle (dört kontur parçası I1–I4)
türetiliyor. Sonuç (Teorem 4.1, eş. 4.25–4.27, DOĞRUDAN ALINTI):

```
Σ_{γ,γ'≤T} f(γ−γ') = (1/(2π)²) ∫₀ᵀ ( 2π f(0) log(t/2π)
    + ∫_{-T}^{T} f(r) ( log²(t/2π) + 2[ (ζ'/ζ)'(1+ir)
        + (t/2π)^{-ir} ζ(1-ir)ζ(1+ir) A(ir) − B(ir) ] ) dr ) dt
    + O(T^{1/2+ε})

A(η) = Π_p (1 − p^{-1-η})(1 − 2/p + p^{-1-η}) / (1 − 1/p)²          (4.26)
B(η) = Σ_p [ log p / (p^{1+η} − 1) ]²                                (4.27)
```

Metin kendi ifadesiyle (satır ~1817, DOĞRUDAN ALINTI, Türkçeye çevrilmeden):
**"We believe that this formula, originally found by Bogomolny and Keating [3], is very
accurate, indeed, down to a square root error term. It includes all of the lower order terms
that arise from arithmetical considerations..."**

Bu, SORU 1'in ilk yarısını doğrudan cevaplıyor: **HAYIR** — Conrey–Snaith'in ratios-sanısından
türettiği formül, BK'nın (1996) orijinal formülünün ÜSTÜNE 1/log T mertebesinde YENİ bir terim
eklemiyor; metnin kendi ifadesiyle bu AYNI formül ("originally found by Bogomolny and Keating"),
yalnız farklı bir yöntemle (ratios sanısı) yeniden türetilmiş ve hata O(T^{1/2+ε}) — yani
1/log T'den çok daha küçük — diye iddia ediliyor.

**Kendi kontrolümüz (bu oturumda yapıldı, metinde YAZILI DEĞİL):** A(η)'nın, Bogomolny 2007'nin
Φ^(off)(ε) = Π_p[1−(1−p^{iε})²/(p−1)²] (eş. 7.6) ile cebirsel olarak özdeş olup olmadığını
kontrol ettik. p=2, ε=1 için sayısal hesap: Φ_off,2(1) ≈ 1.35505 + 0.29497i,
A_2(i·1) ≈ 1.3552 − 0.2948i. Yani **A(iε) = conj(Φ^(off)(ε)) = Φ^(off)(−ε)** — aynı fonksiyon,
yalnız ε işareti/karmaşık-eşlenik kuralı farklı (bu, r vs. r' etiketleme kuralından kaynaklanan
bir konvansiyon farkı, yeni bir arithmetik içerik değil). Bu, metnin kendi sözlü iddiasını
(aynı formül) cebirsel olarak teyit ediyor.

**Sonuç (1.1):** Conrey–Snaith'in ratios-sanısı yönteminde κ_p(L)'yi 1'e YAKINSATACAK ya da
ondan UZAKLAŞTIRACAK, BK'nın Φ_off'unun ÖTESİNDE hiçbir ek terim YOK. A(η), B(η) tamamen BK'nın
Φ^(off), Φ^(diag)'ıyla aynı — parametresiz bir κ_p düzeltmesi bu formülden ÇIKMAZ.

### 1.2 BBLM (math/0602270): N_eff, α — ama YANLIŞ REJİMDE

Kaynak 3, tam metin (9 sayfa). Özet (birebir): *"...these deviations are the same as those of
unitary random matrices of finite dimension N_eff = log(E/2π)/√(12Λ), where Λ = 1.57314...
is a well defined constant."* Yayın: J. Phys. A: Math. Gen. 39 (2006), 10743–10754; gönderim
13 Şubat 2006.

**Türetimin iskeleti (eş. 9–19, tam metinden):** BBLM, BK'nın AYNI r₂(ε) formülünü (onların
eş. 9–12, Bogomolny'nin (6.1)–(7.6)'sıyla birebir aynı — Φ^(diag), Φ^(off)) alıp, **ε → 0
(küçük ε, yani ω≈0 civarı — bizim uydu bölgemiz ω≈L'nin TAM TERSİ)** limitinde Taylor açıyor:

```
r_2^(off)(ε) = (1/4π²)[1/ε² + (γ0²+2γ1+c0) + iQε + O(ε²)] e^{i2πρε} + c.c.     (BBLM eş. 15)
Q = Σ_p log³p/(p−1)²  = δ (Bogomolny 2007'deki AYNI sabit, ≈ 2.3157)
Λ = γ0²+2γ1+c0 = β (Bogomolny 2007'deki AYNI sabit, ≈ 1.57314)   [c0=Σ_p log²p/(p−1)² —
   bu oturumda cebirsel olarak Σ_r(r−1)/p^r = 1/(p−1)² özdeşliğiyle Bogomolny'nin β'sına
   eşitliği doğrulandı]

N_eff = πρ/√(3Λ) = log(E/2π)/√(12Λ)                                            (BBLM eş. 19)
α = 1 + Q/(2πρΛ) = 1 + C/log(E/2π),  C = Q/Λ ≈ 1.4720                          (BBLM eş. 18)
```

Bu SAYILAR Bogomolny 2007'nin §8'iyle (eş. 8.6) tam örtüşüyor (bu oturumda çapraz kontrol
edildi: N_eff = πd̄(E)/√(3β), α=1+δ/(β log(E/2π)) — aynı formül, aynı sabitler).

**KRİTİK NOKTA (metinden açıkça okunuyor, satır ~204–206, DOĞRUDAN ALINTI):**
*"We are interested in the corrections to the asymptotic behavior of R2(s) in the limit when
E→∞. In this limit ρ→∞ and the argument of r2 in Eq.(13) becomes small (keeping s finite).
Therefore, one can expand r2(ǫ) for ǫ≪1."*

Yani N_eff/α düzeltmesi, **ε≪1 (unfolded s = O(1), yani en-yakın-komşu ARALIK istatistiği)**
rejiminde, r2(ε)'nin ε=0 civarındaki Taylor açılımından geliyor. Bizim κ_p'miz ise **ε≈L
(çok BÜYÜK ε, uydunun KENDİSİ)** bölgesinde Φ_off'un genliğini ölçüyor — bu, aynı formülün
TAMAMEN FARKLI bir argüman bölgesi. BBLM/Bogomolny'nin metninde, Φ_off'un ε≈L bölgesindeki
genliğine (yani bizim κ_p'ye karşılık gelen niceliğe) yönelik AYRI bir sonlu-E düzeltmesi
**YOK** — orada formül (7.5)/(12) "kapalı form" (resummed Hardy–Littlewood tekil serisi)
olarak sunuluyor ve iddia edilen hata O(1/√E) (satır: *"The very optimistic error in these
formulae is supposed to be of the order of 1/√E"*, Bogomolny 2007 §7) — bu iddia KANIT değil,
"çok iyimser" (very optimistic) diye nitelenen bir varsayım, ve Odlyzko karşılaştırması
(Bogomolny 2007 Şek. 3–5) E ~ 10^23'üncü sıfır civarında (L≈44) yapılmış — **bizim ölçüm
aralığımızın (L≈9–12) ÇOK ÜSTÜNDE**; yani bu "1/√E" iddiasının bizim L aralığımızdaki
doğruluğu birincil kaynakta hiç test edilmemiş.

**Sonuç (1.2):** BBLM'nin N_eff/α'sı gerçek bir "sonlu-yükseklik düzeltmesi" ama başka bir
nesne için (en-yakın-komşu aralık dağılımının GUE-çekirdek zarfı, ε≪1); bunu κ_p(L)'ye
(ε≈L bölgesi) mekanik olarak taşımak literatürde YAPILMIYOR ve YAPILAMAZ (farklı Taylor açılım
rejimi). §4'te (SORU 4) bunun ne kadar spekülatif bir ekstrapolasyon olacağını ayrıca tartışıyoruz.

### 1.3 "1996/2013" — DOĞRULANAMADI

Görevde anılan "Bogomolny–Keating'in 1996/2013 çalışmaları" ifadesindeki 2013 tarihli ortak
yazarlı bir makale **bulunamadı** (bkz. §0, arXiv API sorgusu 0 sonuç). En yakın 2013-tarihli
ilişkili yayın — Rodgers (2013, Q.J.Math, arXiv:1203.3275, TEK yazarlı) — önceki oturumda
(`LITERATUR_UYDU_KLASIK_MI_24EYL2026.md`, madde 7) tam metinden okunmuştu: Montgomery/BK'nın
çift-korelasyon iddiasını yalnızca |α|≤1−ε (bizim ω<L−ε karşılığı) için RİGOROUS kanıtlıyor;
kendi ifadesiyle α'nın bu aralığın ÖTESİNE (bizim uydu bölgemiz, α≳1) genişletilmesi
*"will almost certainly require a breakthrough"*. Bu, κ_p(L) sorusunun matematiksel
literatürün KANITLANMAMIŞ ucunda kaldığını bağımsız olarak teyit ediyor. **DOĞRULANAMADI**:
görevdeki "2013" atfının Rodgers'a mı yoksa gerçekten var olmayan bir BK ortak makalesine mi
işaret ettiği.

---

## SORU 2 — Landau–Gonek ve genellemeleri (özellikle Durkan–Hughes–Pearce-Crump arXiv:2601.18025)

### 2.1 Klasik zincir — tam metinden (DHPC'nin kendi Bölüm 1'inden, birebir aktarılan teoremler)

```
Teorem 1 (Landau, 1911 — DHPC'nin ifadesiyle). Sabit X>1 için
  Σ_{0<γ≤T} X^ρ = −(T/2π) Λ(X) + O(log T),  T→∞.

Teorem 2 (Gonek, 1985 [DHPC bibliyografyasında [11,12] — 1985 duyuru + 1993 tam makale]).
  X,T>1 için DÜZGÜN (uniform):
  Σ_{0<γ≤T} X^ρ = −(T/2π)Λ(X)
      + O( X log(2XT) loglog(3X) )
      + O( log X · min(T, X/⟨X⟩) )
      + O( log(2T) · min(T, 1/logX) )
  burada ⟨X⟩, X'e en yakın (X'ten farklı) asal kuvvete olan uzaklık.

Sonuç 2.1 (Gonek, 1993 — fonksiyonel denklemden, X>1 için X^{-ρ} toplamı, yani
  x=1/X<1 durumunun karşılığı):
  Σ_{0<γ≤T} X^{-ρ} = −(T/2π)·Λ(X)/X
      + O( log(2XT) loglog(3X) )
      + O( log X · min(T/X, 1/⟨X⟩) )
      + O( log(2T) · min(T/X, 1/(X log X)) )
```

**Ana terimlerin simetrisi (bizim türetimimiz, ρ=1/2+iγ ile):** x^{iγ}=x^{ρ-1/2}, dolayısıyla
Σx^{iγ} = x^{-1/2}ΣX^ρ (x=X>1 durumu, Teorem 2 kullanılır) VEYA X^{1/2}·ΣX^{-ρ}
(x=1/X<1 durumu, Sonuç 2.1 kullanılır, X=1/x). Her iki yoldan da **ana terim aynı**:
−T·Λ(X)/(2π√X) — görev metnindeki "−T Λ(x)/(2π√x)" ifadesiyle birebir örtüşüyor
(x=X yerine konursa). Bu simetri BK'nın "ayna" (mirror) yapısının klasik köküdür.

**Hata terimlerinin ASİMETRİSİ (birincil kaynaktan doğrudan okunuyor, görev sorusunun
yanıtı):** Teorem 2'nin (X>1, bizim "paydaki/pozitif üs" yönümüz) baş hata terimi
**O(X log(2XT) loglog(3X))** — X ile ÇARPIMSAL büyüyor. Sonuç 2.1'in (X>1 ama X^{-ρ}
toplamı, bizim "ayna/negatif üs" yönümüz) karşılık gelen terimi ise
**O(log(2XT) loglog(3X))** — **X'siz**, yani çok daha küçük. Bu, iki yönün hata
yapısının GERÇEKTEN ASİMETRİK olduğunu birincil kaynaktan doğruluyor — görev sorusunun
"pozitif ve negatif log x için asimetri var mı?" sorusuna kısmi bir EVET.

**Ama bu asimetri κ_p(L)'yi NİCEL olarak öngörmüyor — üç nedenle (dürüstlük notu):**
1. Bu hata terimleri X SABİT, T→∞ asimptotiğinde; bizim rejimimizde X=p^k SABİT ve L=log(T/2π)
   sadece 9–12 (T~10⁴–10⁶) — asimptotik rejimde miyiz, belirsiz; ayrıca hata terimleri
   ÜST SINIR (upper bound), işaretli/kesin bir düzeltme DEĞİL — yani "asimetri var" demek
   "κ_p bu yüzden <1 çıkar, κ_{1/p} =1 kalır" demek DEĞİL, sadece potansiyel farklı hata
   büyüklüğü.
2. Landau–Gonek TEK bir X (tek asal kuvvet) için Σx^{iγ} (BİR-NOKTA toplamı) veriyor; bizim
   κ_p, BK'nın Φ_off'undaki (İKİ-NOKTA korelasyon, Λ(n1)Λ(n2) çifti) Euler-çarpanı genliği —
   FARKLI bir istatistik. İkisini birbirine bağlayan bir teorem bu 3 kaynakta (Landau, Gonek,
   DHPC) YOK.
3. DHPC'nin KENDİ YENİ SONUCU (Teorem 5, aşağıda) Σx^{iγ}'nin kendisini DEĞİL, χ(ρ)X^ρ
   (fonksiyonel-denklem-burgulu, FARKLI bir toplam) genişletiyor — bu, Σx^{iγ}'ye ikincil
   terim EKLEMİYOR, ayrı bir teorem.

### 2.2 DHPC'nin (2601.18025) GERÇEK yeni katkısı — Σx^{iγ}'nin kendisi DEĞİL

Önceki oturumun notu (`LITERATUR_UYDU_KLASIK_MI_24EYL2026.md`, madde 12) bu makaleyi
"x≍T rejiminde üniform Landau–Gonek güncellemesi... MUTLAKA okunmalı" diye işaretlemişti —
tam metin şimdi okundu ve bu beklenti **KISMEN YANLIŞ ÇIKTI**: makalenin yeni sonucu
(Teorem 5, tam metinden, DOĞRUDAN ALINTI biçiminde yukarıda §durkan bölümünde verildi)

```
S(X,T) = Σ_{T<γ≤2T} χ(ρ) X^ρ
```

toplamı için — burada χ(s), ζ(s)=χ(s)ζ(1−s) fonksiyonel denklem çarpanı — Σx^{iγ}'nin KENDİSİ
DEĞİL. Teorem 5, X'in T/(2π) ve T/π eşiklerine göre ÜÇ rejime ayrıldığını gösteriyor (eş. 2.1,
tam metinden doğrulandı) ve DHPC bunu Shanks sanısının (Σζ'(ρ) gerçel/pozitif) ve genelleşmiş
Shanks sanısının (Σζ^(ν)(ρ)) yeni ispatları için kullanıyor (§5–6). **Σ_{0<γ≤T} x^{iγ}'nin
kendisine (Landau/Gonek'in konusu) DHPC hiçbir yeni ikincil terim EKLEMİYOR** — yalnız
Teorem 1/2/Sonuç 2.1'i (Landau, Gonek) olduğu gibi kaynak olarak kullanıyor (§5, eş. 5.3–5.5
türetiminde). **Sonuç: önceki oturumun "bu makale bizim ihtiyacımız olan güncelleme olabilir"
beklentisi bu tam okumayla DÜZELTİLDİ — değil.**

### 2.3 Fujii'nin ikincil terimleri — hâlâ DOĞRULANAMADI (içerik)

DHPC'nin kendi metni (Bölüm 1, satır ~69) Fujii'yi şöyle anıyor: *"Gonek's result was further
generalised by Fujii [6, 7] who found lower order terms in the expansion."* — bu, Fujii'nin
TAM OLARAK Teorem 2'nin (Gonek) hata terimini daha ince ikincil terimlere açtığını gösteriyor,
ama bu ikincil terimlerin KENDİSİ (işareti, log p/log T mertebesi, x>1/x<1 asimetrisi) DHPC'nin
metninde AÇIK YAZILI DEĞİL (yalnız "lower order terms" diye anılıyor, formül verilmiyor).
Fujii'nin 1989/1990 makalelerinin kendisine erişilemedi (Project Euclid ücretli, önceki
oturumda da böyleydi) — **DOĞRULANAMADI**: bu ikincil terimlerin tam biçimi.

---

## SORU 3 — μ=0 (x=p^k, k≥2) konumlarında küçük negatif ikincil terim öngörülüyor mu?

**Kısa cevap: bu üç kaynakta (Bogomolny 2007/BK, Conrey–Snaith, DHPC/Landau–Gonek) böyle bir
öngörü YOK / bulunamadı — DOĞRULANAMADI.**

Gerekçe (birincil kaynaktan): BK'nın Φ^(off)(ε)=Π_p[1−(1−p^{iε})²/(p−1)²] (Bogomolny 2007 eş.
7.6) Euler çarpanının p asalı için k-inci uyduya (r=p^k) karşılık gelen katsayısı f_p(k≥2)=0
sonucu, Φ_off'un p^{iε} kuvvet serisine AÇILIMINDAN gelir — bu, Hardy–Littlewood tekil serisinin
(eş. 6.8) q,p toplamının RESUMMASYONUNDAN çıkan CEBİRSEL bir sonuç (Bogomolny 2007 §6–7'nin
TAMAMI, eş. 6.5→7.6 türetim zinciri), yani zaten "sıfırıncı mertebe" (leading, exact closed
form) BK tahmininin İÇİNDE, ayrı bir "ikincil düzeltme terimi" DEĞİL.

Landau–Gonek tarafında ise durum FARKLI bir istatistikle ilgili: Λ(p^k)=log p, **k'DAN
BAĞIMSIZ** (von Mangoldt fonksiyonunun tanımı gereği — herhangi bir k≥1 için Λ(p^k)=log p).
Yani Landau–Gonek'in (Teorem 1/2) Σx^{iγ} ANA TERİMİ, k=1 (p^1) ile k≥2 (p^k) arasında HİÇBİR
FARK göstermiyor — bu, BK'nın f_p(k≥2)=0 (İKİ-NOKTA korelasyon, farklı nesne) sonucuyla
DOĞRUDAN karşılaştırılabilecek bir şey değil, çünkü Landau–Gonek BİR-NOKTA toplamı (tek γ
üzerinden), BK ise İKİ-NOKTA (γ,γ' çifti üzerinden). Bu iki nesneyi birbirine bağlayan, "BİR-NOKTA
toplamın k-bağımsızlığından İKİ-NOKTA korelasyonun k≥2'de küçük bir kalıntıya sahip olması
gerektiğini" çıkaran bir teorem/argüman bu üç kaynakta **bulunamadı**.

Bizim ölçtüğümüz μ=0 konumlarındaki küçük negatif artık (ölçeğin ~%2'si, L'den bağımsız) —
bu niceliğin kaynağı olabilecek iki aday (ikisi de metinde AÇIK YAZILI DEĞİL, bizim
spekülasyonumuz):
(a) R_2^(diag)(ε) (Bogomolny 2007 eş. 6.2–6.4) düzgün (smooth) bir zarf olarak ε=k log p
    noktalarında da sıfırdan farklı bir değer taşıyabilir — ama bu bir "çizgi" değil sürekli
    bir fon, doğrudan bizim ayrık-uydu ölçümümüze nasıl karışacağı belirsiz;
(b) Φ_off'un p^{iε} açılımının HARDY–LITTLEWOOD tekil serisinin kendisinin (eş. 6.8) sonlu-q
    kesilmesinden (resummasyonun T'ye bağlı yakınsaklık hızından) gelen bir artık.
Her ikisi de **DOĞRULANAMADI** — bu üç birincil kaynakta açıkça ele alınmıyor.

---

## SORU 4 — Genel hüküm: parametresiz κ_p(L) literatürden çıkarılabilir mi?

**HAYIR — literatürden parametresiz bir κ_p(L) formülü çıkarılamıyor.** Gerekçe, üç bölümün
sentezi:

1. BK'nın Φ_off(ε) formülü (Bogomolny 2007 eş. 7.5–7.6 / BBLM eş. 12 / Conrey–Snaith'in A(η)'sı
   — üçü de bu oturumda CEBİRSEL/SAYISAL olarak aynı nesne diye doğrulandı) **kapalı form**
   (resummed Hardy–Littlewood) olarak sunuluyor; L'ye (yani E/T'ye) bağlı, satelit genliğini
   1'den küçültecek AYRI bir çarpan formülün içinde YOK — formül "exact, hata O(1/√E)" diye
   iddia ediliyor (ama bu iddia "very optimistic" diye nitelenmiş bir varsayım, kanıtlı değil,
   ve bizim L≈9–12 aralığında hiç test edilmemiş — bkz. §1.2).
2. Var olan TEK somut "sonlu-yükseklik/sonlu-boyut" düzeltmesi — BBLM'nin N_eff, α'sı — BAŞKA
   bir rejimde (ε≪1, en-yakın-komşu aralık istatistiği) türetilmiş; bunu ε≈L (uydu) bölgesine
   TAŞIYAN bir teorem yok.
3. Landau–Gonek/DHPC tarafı BİR-NOKTA istatistiği; İKİ-NOKTA (BK/κ_p) istatistiğine BAĞLAYAN
   bir köprü teoremi yok (§2.1, madde 2).

**Eksik parça, açıkça:** Φ_off(ε)'un (ya da eşdeğer A(η)'nın) türetimi — Bogomolny 2007 §6–7,
eş. 6.5→7.6 — sonsuz-T (T→∞) limitinde YAPILIYOR: "n1,n2→∞" toplamları sürekli hale getirilip
(eş. 6.6→7.1, "Changing the sum over n to the integral") kapalı forma resumme ediliyor. Bu
adımda (ayrık toplamdan sürekli integrale geçiş) SİSTEMATİK olarak atılan, T'ye (dolayısıyla
L'ye) bağlı bir kalan terim OLASI — ama bu kalan terimin p ve L'ye göre nasıl ölçeklendiği
(bizim ölçtüğümüz 1−κ_p ∝ L^{−γ}, γ=1.36±0.21 gibi) hiçbir okunan kaynakta AÇIKÇA
hesaplanmamış/yazılmamış. Bu, tam olarak Rodgers (2013)'ın *"to rigorously extend the range...
will almost certainly require a breakthrough"* dediği boşluk (§1.3) — κ_p(L)'nin L'ye göre
davranışı, matematiksel literatürün bugün (2026 itibarıyla, okunan kaynaklara göre) eli
ULAŞMADIĞI bir soru.

**Salt sağlama amaçlı, literatür-DIŞI bir karşılaştırma (KESİNLİKLE bir öngörü değil, sadece
bir büyüklük-mertebesi kontrolü — bu paragrafın hiçbir cümlesi birincil kaynaktan alıntı
DEĞİLDİR):** BBLM'nin CUE_N açılımında ilk düzeltme 1/N² mertebesinde (eş. 8.7, Bogomolny
2007: p^(CUE_N)(s)=p₀(s)+p₁(s)/N_eff²+O(N_eff⁻⁴)) — bu "tipik RMT sonlu-boyut düzeltmesi
γ=2'dir" sezgisini çağırıyor (N_eff∝L olduğundan 1/N_eff²∝1/L², yani γ=2). Bu, bizim ölçtüğümüz
γ=1.36±0.21 ile GERİLİMLİ (2 için z-skoru ≈ (2−1.36)/0.21 ≈ 3.0, yani ~3σ uzakta) — ama bu
karşılaştırmanın kendisi de sağlam değil: κ_p'nin BBLM'nin N_eff'iyle AYNI mekanizmadan geldiğine
dair hiçbir teorem yok (§1.2, madde 3'te gerekçelendirilen nedenle), dolayısıyla bu ne γ=1'i
(H_S) ne γ=0'ı (H_C) ne de γ=2'yi teorik olarak destekliyor/çürütüyor — yalnızca "eğer κ_p
BBLM'inkiyle aynı jenerik RMT mekanizmasından gelseydi γ≈2 beklenirdi, ama gelmiyor gibi
görünüyor" diye okunabilecek, zayıf ve spekülatif bir gözlem.

**Sayısal öngörü tablosu (L=9.34, 10.48, 12.03; p=2,3,5,7):** Okunan literatürden çıkan TEK
somut, parametresiz nicel öngörü, BK'nın kendisinin sıfırıncı-mertebe (L-BAĞIMSIZ) iddiasıdır:

| p | 2 | 3 | 5 | 7 |
|---|---|---|---|---|
| Literatürün κ_p(L) öngörüsü (BK, tüm L için) | 1 | 1 | 1 | 1 |
| Ölçülen κ_p, L=9.34 | 0.784 | 0.663 | 0.583 | 0.521 |
| Ölçülen κ_p, L=10.48 | 0.789/0.795* | 0.713/0.720* | 0.604/0.608* | 0.538/0.530* |
| Ölçülen κ_p, L=12.03 | 0.837 | 0.770 | 0.712 | 0.616 |

*İki değer görevin bağlam metninde ve KALEM_KAPPA_YUKSEKLIK_25EYL2026.md'de hafifçe farklı
(0.789 vs 0.795 vb.) — bu dosyanın konusu değil, kaynak kalem notuna bakılmalı.

Literatürden L'ye bağlı BAŞKA bir sayı (κ_p(L)≠1 öngören) çıkarılamadı — bu satır boş
bırakılıyor, uydurulmadı.

---

## Kaynakça (bu dosyada geçen, doğrulanmış künyeler)

1. E. Bogomolny, "Riemann zeta function and quantum chaos", arXiv:0708.4223 [nlin.CD] (30 Ağu
   2007); Prog. Theor. Phys. Suppl. 166 (2007), 19–44.
2. J.B. Conrey, N.C. Snaith, "Applications of the L-functions ratios conjectures",
   arXiv:math/0509480 [math.NT] (21 Eyl 2005, v2 20 Tem 2007); Proc. London Math. Soc. 94:3
   (2007), 594–646.
3. E. Bogomolny, O. Bohigas, P. Leboeuf, A.G. Monastra, "On the spacing distribution of the
   Riemann zeros: corrections to the asymptotic result", arXiv:math/0602270 [math.NT]
   (13 Şub 2006); J. Phys. A: Math. Gen. 39 (2006), 10743–10754.
4. B. Durkan, C. Hughes, A. Pearce-Crump, "Generalisations of the Landau–Gonek Theorem and
   applications to mean values of zeta", arXiv:2601.18025 [math.NT] (25 Oca 2026).
5. E.B. Bogomolny, J.P. Keating, "Gutzwiller's Trace Formula and Spectral Statistics: Beyond
   the Diagonal Approximation", Phys. Rev. Lett. 77 (1996), 1472–1475. DOI:
   10.1103/PhysRevLett.77.1472.
6. E.B. Bogomolny, J.P. Keating, "Random matrix theory and the Riemann zeros I: three- and
   four-point correlations", Nonlinearity 8 (1995), 1115–1131; "...II: n-point correlations",
   Nonlinearity 9 (1996), 911–935.
7. B. Rodgers, "Macroscopic pair correlation of the Riemann zeroes for smooth test functions",
   arXiv:1203.3275 [math.NT] (2012); Q. J. Math. 64:4 (2013), 1197–1219.
8. E. Landau, "Über die Nullstellen der Zetafunction", Math. Ann. 71:4 (1912), 548–564
   (DHPC'nin kendi metni Teoremi "Landau, 1911" diye etiketliyor).
9. S.M. Gonek, "A formula of Landau and mean values of ζ(s)", in Topics in Analytic Number
   Theory (Austin, Tex., 1982), Univ. Texas Press, Austin, TX, 1985, 92–97.
10. S.M. Gonek, "An explicit formula of Landau and its applications to the theory of the
    zeta-function", in A Tribute to Emil Grosswald: Number Theory and Related Analysis,
    Contemp. Math. 143, Amer. Math. Soc., Providence, RI, 1993, 395–413.
11. A. Fujii, "On a theorem of Landau", Proc. Japan Acad. Ser. A Math. Sci. 65:2 (1989), 51–54.
12. A. Fujii, "On a theorem of Landau, II", Proc. Japan Acad. Ser. A Math. Sci. 66:9 (1990),
    291–296.
13. A. Fujii, "On a conjecture of Shanks", Proc. Japan Acad. Ser. A Math. Sci. 70:4 (1994),
    109–114.
14. A. Fujii, "On the distribution of values of the derivative of the Riemann zeta function at
    its zeros. I", Tr. Mat. Inst. Steklova 276 (2012), 51–76.

(Künye 5–6, 8–14: bu oturumda birincil metne erişilemedi; bibliyografik ayrıntılar kaynak 1,
3, 4'ün KENDİ tam-metin kaynakçalarından — yani birincil kaynaklardan — birebir aktarıldı.)
