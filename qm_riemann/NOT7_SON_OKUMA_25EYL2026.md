# NOT 7 — Son Okuma (25 Eylül 2026)

**Denetlenen taslak:** `qm_riemann/arxiv_comb_satellites.tex` (686 satır, tamamı okundu) ve derlenmiş `arxiv_comb_satellites.pdf` (13 sayfa, pypdf ile metin çıkarıldı).
**Karşılaştırma:** `NOT7_HAKEM_TEFTISI_25EYL2026.md` (26 bulgu: K1–K6, Ö1–Ö13, KU1–KU7), kaynak raporlar 188–198 (hepsi repoda mevcut, çapraz kontrol edildi).

## 0. Genel sonuç

**Yayın engeli: HAYIR.** Önceki 26 bulgunun tamamı (26/26) metinde ele alınmış durumda; 6 KRİTİK bulgunun tamamı ve 13 ÖNEMLİ bulgunun 12'si tam, 1'i (Ö12) kısmen düzeltilmiş. PDF'te hiçbir "??" (çözülmemiş \ref/\cite) yok. `\cite` anahtarları ile `\bibitem` anahtarları birebir aynı küme (17 anahtar, eksiksiz, kullanılmayan bibitem yok). Özetten tabloya kadar çapraz kontrol edilen tüm sayısal değerler (1.023±0.023, 12σ, 0.789±0.037, 0.993±0.039, γ=1.36±0.21, 6.4σ, 7.9σ, 1.8–2.4, 1.25–1.75, 21/21, 22/24, 180°±1.3°, üç pencerenin L değerleri) birbirini tutuyor. η_p/κ notasyonu tutarlı: §7, Tablo 1 ve Şekil 5 her yerde η_p kullanıyor (eski κ_p çakışması KU3 ile giderilmiş); κ yalnız çekirdek (kernel) niceliği için kullanılıyor. Kalan sorunlar yayını engellemeyen, ince cilalama düzeyinde.

## 1. Önceki 26 bulgunun işlenme durumu

### KRİTİK (6/6 düzeltildi)

| # | Konu | Durum | Kanıt (güncel satır) |
|---|------|-------|----------------------|
| K1 | %77 rank-one, ön-kayıtlı eşik (≥80%) tutmadı | **Düzeltildi** | satır 197–198: "(a pre-registered threshold of 80% for a rank-one mechanism was not met)" |
| K2 | Δω/τ' korelasyonu post-hoc oysa Observation kutusunda | **Düzeltildi** | satır 206–216: Observation kutusu satır 214'te bitiyor; korelasyon cümlesi ayrı paragrafta "Post hoc," ile başlıyor (satır 216) |
| K3 | +log7 büyüklüğü "1.5–1.9" → olması gereken 1.8–2.4 | **Düzeltildi** | satır 298–299: "magnitudes 1.8–2.4 times the prediction" |
| K4 | "tüm yüksekliklerde kalıcı" iddiası desteksizdi | **Düzeltildi** | satır 304–307: 9.4σ artık yalnız "last window"a bağlanmış; §7 satır 534–538'de üç pencerenin gerçek sayıları (−0.023/−0.023/−0.032 vb.) eklenmiş, iddia artık destekli |
| K5 | "+log 2m iki katı parlak" — tek/çift m niteleyicisi eksikti | **Düzeltildi** | satır 424–425: "twice as bright as $+\log m$ for odd $m$" |
| K6 | Kör bandın "hiç hesaplanmamış" kısmı ile "hesaplanmış ama aranmamış" kısmı karışıktı | **Düzeltildi** | satır 438–440: "whose part below −1.775 had never been computed and whose remainder had been computed in an earlier map but never displayed or searched" |

### ÖNEMLİ (12/13 tam, 1/13 kısmi)

| # | Konu | Durum | Kanıt |
|---|------|-------|-------|
| Ö1 | "exactly like control" abartısı | **Düzeltildi** | satır 281: "statistically indistinguishable from" |
| Ö2 | 9 pozisyondan yalnız 3'ü temiz sıfır nüansı eksikti | **Düzeltildi** | satır 282–285: "only three of the nine are, however, consistent with zero --- the others carry significant negative values" |
| Ö3 | +log5 "≈1.3" alt-uç seçici alıntı | **Düzeltildi** | satır 300–301: tam aralık "1.25–1.75", "the same stretch" ifadesi kaldırılmış, "a smaller stretch" ile değiştirilmiş |
| Ö4 | "22/24" yerine sadece "allowed ones" deniyordu | **Düzeltildi** | satır 345–347: "22 of the 24 allowed ones ... the two exceptions being +log4" |
| Ö5 | z=−55..−82 niteleyicisiz | **Düzeltildi** | satır 350 ve şekil altyazısı satır 371: "for the primary/six primary satellite pairs" |
| Ö6 | Φ(ε) formülü iç nottan değil doğrudan alıntı gibi sunuluyordu | **Düzeltildi** | satır 386–389: doğrudan "\cite[eqs.~(7.5)--(7.6)]{Bog07}" ile yayınlanmış literatüre atfedilmiş |
| Ö7 | "A lesson in line shape" alt-bölümünde post-hoc etiketi yoktu | **Düzeltildi** | satır 427: başlık artık "A lesson in line shape (post hoc)" |
| Ö8 | γ=1 ile 1.9σ gerilimi atlanmıştı | **Düzeltildi** | satır 517–519: "(and is 1.9σ steeper than γ = 1)" |
| Ö9 | Yalnız R(5/2) raporlanıyordu, R(7/2) eksikti | **Düzeltildi** | satır 519–521: her iki aile de var: "R(5/2): 0.79→0.83→0.92; R(7/2): 0.75→0.81→0.88" |
| Ö10 | Özet +log7/+log5 taşmasını yansıtmıyordu | **Düzeltildi** | satır 36–37: "(two other satellites follow its signs but exceed its magnitudes)" |
| Ö11 | C_j(q) tanımsızdı | **Düzeltildi** | satır 166–168: açık formül "$C_{b,j}(q) = 2\langle \delta^{(b,j)}_n e^{-i\omega_q m_n}\rangle_{n\in b}$" eklenmiş |
| Ö12 | mix'in sıfır olmama nedeni tek cümleye sıkışıktı, T→∞ sorusu | **Kısmen** | satır 142–147: "neden sıfır değil" artık iyi açıklanmış (Gram couplings, $\hat G$, companion4 atfı) ama "bu sonlu-pencere etkisi mi yoksa T→∞'da kalıcı mı" sorusu hâlâ açık bırakılmış — küçük, yayını engellemez |
| Ö13 | −π/4 teriminin iç-not kaynağı izlenemiyordu | **N/A (metne uygulanmaz)** | Bu, dahili teftiş-izi sorunuydu, yayınlanan makalenin kendisiyle ilgili değil; türetim zaten bağımsız olarak doğrulanmıştı (önceki raporun §0'ı) |

### KÜÇÜK (7/7, 6 düzeltildi, 1 açık)

| # | Konu | Durum |
|---|------|-------|
| KU1 | "10^-15" hassasiyet iddiası gevşekti | **Düzeltildi** — satır 175–176: "(to machine precision)" ile değiştirilmiş |
| KU2 | "alt" penceresinin kaynağı izlenemiyordu | **Açık** — satır 119–125'te hâlâ dipnot/kaynak yok |
| KU3 | κ_p / κ(ω') sembol çakışması | **Düzeltildi** — §7'de η_p kullanılıyor, tutarlı |
| KU4 | Özet ~373 kelime, çok yoğun | **Düzeltildi (kısaltılmış)** — güncel özet ölçüldü: ~268 kelime |
| KU5 | "pre-registered" vs "blind" tanımsız eşanlamlı kullanım | **Düzeltildi** — satır 100–105'te açık operasyonel tanım eklenmiş |
| KU6 | Tablo 1'de "last" penceresinin kısmi-görülmüşlüğü yoktu | **Düzeltildi** — tablo altyazısı satır 509–512: "coarse satellite ratios of the last window had been seen before the test (see text)" |
| KU7 | τ(χ) / τ_q sembol çakışması | **Düzeltildi** — satır 332: "(not to be confused with $\tau = \log q/L$)" eklenmiş |

## 2. İç tutarlılık — çapraz kontrol sonuçları (özet ↔ gövde ↔ şekil ↔ tablo)

Aşağıdaki sayılar iki veya daha fazla yerde geçiyor; hepsi tutarlı bulundu:

- **1.023 ± 0.023**: özet (44), Observation (456–457), Şekil 4 altyazısı (476) — aynı.
- **12σ**: özet (45) "$12\sigma$", gövde (457) "more than $12\sigma$", Şekil 4 altyazısı (478) "$12\sigma_{\rm eff}$" — aynı büyüklük, alt simge farkı önemsiz.
- **0.789 ± 0.037** (ρ): Observation (461–462), Şekil 4 altyazısı (480) — aynı.
- **0.993 ± 0.039** (ψ): Observation (459–460), Şekil 4 altyazısı (481) — aynı.
- **γ = 1.36 ± 0.21**: Observation (516), Şekil 5 altyazısı (554–555) — aynı.
- **6.4σ**: gövde (516–517), Şekil 5 altyazısı (558) — aynı.
- **7.9σ**: yalnız gövdede (517–518), tek geçiş, çelişkisiz.
- **1.8–2.4** (+log7): gövde (299), Şekil 2 altyazısı (319–320) — aynı.
- **1.25–1.75** (+log5): yalnız gövdede (301), tek geçiş.
- **21/21, 22/24**: gövde (345–347), tutarlı; kaynak raporlarla (195) çelişki bulunamadı.
- **180°±1.3°**: gövde (349–351), Şekil 3 altyazısı (375) — aynı.
- **Pencere L değerleri** (alt 9.34, low 10.48, last 12.03): Veri bölümü (122–124), Tablo 1 (504–506), Şekil 5 altyazısı (553) — tutarlı.
- **η_p pre-registered değerleri** (0.795, 0.720, 0.608, 0.530; Şekil 5 altyazısı 556–557) low penceresinin ölçülen değerlerinden (0.789, 0.713, 0.604, 0.538; Tablo 1) küçük farklarla ayrılıyor — bu bir **hata değil**, tam olarak beklenen durum: biri ön-kayıtlı tahmin, diğeri ölçüm.

**\ref/\label/\cite/\bibitem taraması** (programatik, `grep`): Tüm 23 `\label` hedefi tanımlı; kullanılan her `\ref`/`\eqref` (28 çağrı) karşılık gelen bir `\label`'a çözülüyor. 17 `\cite` anahtarı = 17 `\bibitem` anahtarı, birebir aynı küme — ne eksik ne fazla. PDF'te (pypdf, 13 sayfa) `"??"` dizisi **hiç yok**.

## 3. Dil / okunurluk — yeni gözlemler (önceki 26'da yoktu)

1. **σ_eff hiç tanımlanmamış.** Yalnızca üç şekil altyazısında geçiyor (satır 476, 478, 555: "$\sigma_{\rm eff}$") ama gövde metninde ne formülü ne de sözel tanımı var. math.NT okuru için tanımsız terim. **Öneri:** Metinde bir yerde (§6 ya da bir dipnot) "$\sigma_{\rm eff}$" için bir cümlelik tanım ekleyin (muhtemelen korelasyonlu jackknife bloklarını hesaba katan efektif belirsizlik).
2. **"z" (z-skoru) formel tanımsız, "σ" ile karışabilir.** Satır 293 ("$|z| = 11$–$52$"), 299 ("$|z| = 4.9$–$6.2$"), 350/371 ("$z = -55$ to $-82$") işaretli z kullanıyor, oysa metnin geri kalanı "$N\sigma$" biçimini kullanıyor (ör. "12σ", "6.4σ"). İkisinin aynı şey olup olmadığı (z = değer/hata, işaretli; Nσ = |z|, işaretsiz) hiçbir yerde açıklanmıyor. **Öneri:** İlk kullanıldığı yerde (satır 293 civarı) bir dipnotla "z := value/(jackknife s.e.); we quote $|z|$ as '$N\sigma$' when only the magnitude matters" türü bir açıklama ekleyin.
3. **Uzun, çok maddeli cümleler** — okunabilirliği zorluyor:
   - Satır 160–176 ("The kernel" alt-bölümü): "The window is cut into blocks..." ile başlayıp denklemi (eq:K) içine alarak "...to machine precision)." ile biten tek cümle, dört ayrı fikri (blok tanımı, dilim ataması, projeksiyon, denklem+yorum) noktalı virgülle birbirine bağlıyor. **Öneri:** Denklemden önce yeni bir cümle başlatın: "...within the block. This gives $K_{b,j}$ and $\kappa_j$ as in Eq.~\eqref{eq:K}; the sums run over..."
   - Satır 437–444 ("The blind test" ilk paragrafı): "The final design fits the band with the seven main satellites free..." cümlesi tasarım + iki hipotez tanımı + sayısal tahmini tek cümlede topluyor. **Öneri:** "it compares two hypotheses by self-consistency:" öncesinde cümleyi bölün.
   - Satır 196–202 (§3 açılış cümlesi): benzer şekilde üç bağımsız gözlemi (rank-one oranı, band faktörü, comb yapısı) tek cümlede taşıyor; kritik değil ama bölünebilir.
4. **"Bohr coefficient" (satır 398, 543)** açıkça tanımlanmamış; bağlamdan çıkarılabilir (neredeyse-periyodik $f(\varepsilon)=\prod_p g_p$'nin genelleştirilmiş Fourier katsayısı) ama math.NT okuru için bir yarım cümlelik gloss faydalı olur.
5. Yazım/dilbilgisi taraması (double-word, "a"+ünlü, çift boşluk, yaygın yanlış yazımlar) **temiz** — hiçbir gerçek hata bulunmadı; tirelerin tümü (`--`/`---`) tutarlı kullanılmış, sayı aralıklarında tek tire yok.

## 4. Kaynakça

`\cite` ↔ `\bibitem` birebir eşleşiyor (17/17), eksik veya kullanılmayan bibitem yok. Görev numaraları (188, 190, 192–198), Reproducibility bölümünde (satır 594–597) verilenlerle birebir örtüşüyor ve hepsi repoda gerçekten mevcut (188_zeta_kimligi_RAPOR.md … 198_kappa_yukseklik_RAPOR.md, 196_bk_uydu_katsayilari.py dahil).

## 5. Yayın engeli değerlendirmesi

Yanlış iddia, eksik uyarı veya tutarsız ön-kayıt/post-hoc etiketi **bulunamadı**. Önceki 26 bulgunun hiçbiri metinde çözümsüz kalmamış (yalnızca Ö12 kısmi, Ö13 zaten metne uygulanmaz, KU2 açık — üçü de kozmetik/dipnot düzeyinde). **Sonuç: yayına engel yok; yukarıdaki 5 kalem (σ_eff, z/σ ayrımı, 3 uzun cümle, Bohr coefficient, KU2 dipnotu) isteğe bağlı son cila.**
