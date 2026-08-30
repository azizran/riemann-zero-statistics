# BLS TAM METİN OKUMASI — Bohigas, Leboeuf, Sánchez (2001)

**Makale:** O. Bohigas, P. Lebœuf, M. J. Sánchez, "Spectral spacing correlations
for chaotic and disordered systems", *Foundations of Physics* **31** (2001)
489–517; arXiv:nlin/0012049. Gutzwiller Festschrift özel sayısı için yazılmış;
metin içi tarih 1 Temmuz 1999.

**Okunan kaynak:** ar5iv HTML tam metni
(`https://ar5iv.labs.arxiv.org/html/nlin/0012049`), LaTeXML dönüşümü. Baştan
sona okundu: Bölüm I–IV, denklem (1)–(53), Tablo I, Şekil 1–6 altyazıları,
25 kaynaklık kaynakça. **Erişilemeyen tek şey: şekillerin kendisi (grafik
görüntüleri) ve dergi sayfa numaraları.** Bu yüzden aşağıda hiçbir yerde sayfa
numarası verilmiyor — referanslar denklem / bölüm / şekil numarasıyla. Bunlar
arXiv sürümünün numaralandırmasıdır; dergi sürümüyle aynı olması beklenir ama
alıntı yapmadan önce PDF'ten teyit edilmeli.

**Okuma tarihi:** 30 Ağustos 2026.

---

## 1. MAKALENİN TAM ÖZETİ — ne yapıyor, nasıl

### 1.1 Sorduğu soru

Kuantum kaotik ve difüzif sistemlerin spektral dalgalanmaları için standart
nesne iki-nokta fonksiyonu `R₂(ε)` ya da sayım varyansı `Σ²(L)`. BLS bunun
yerine **ardışık seviye aralıklarının öz-kovaryansını** çalışıyor:

> `C(n) = ⟨(s_m − ⟨s_m⟩)(s_{m+n} − ⟨s_{m+n}⟩)⟩ = ⟨s_m s_{m+n}⟩ − 1 = I(n) − 1`
> — denk. (1), `s_m = x_m − x_{m−1}`, `x` katlanmamış (unfolded) seviyeler.

Motivasyon (Bölüm I, açık biçimde yazılmış): `C(n)`, `R₂`'nin **ayrık** bir
sürümüdür; ayrıklaştırma ortalama-aralık ölçeğinde bir yumuşatma yapar ve o
ölçekteki yapıları güçlü biçimde bastırır. Zaman dilinde bu, Heisenberg
zamanında bir kesim demektir. Sonuç: **yarı-klasik yaklaşımlar `C(n)` için
`R₂` için olduğundan çok daha iyi çalışır.** Bu, makalenin bütün stratejik
zeminidir.

### 1.2 Türetme yolu (Bölüm II.A) — sayım varyansından aralık kovaryansına

Zincir kısa ve tamamen açık:

1. `σ²(n)` = `n` ardışık aralıktan oluşan `S = Σ s_i` uzunluğunun varyansı;
   denk. (2) ile `I(j)`'ye bağlı, denk. (3)–(4) ile ters çevrilmiş:
   `C(n) = ½[σ²(n+1) − 2σ²(n) + σ²(n−1)]`, `n ≥ 2`.
2. **Köprü varsayımı** — denk. (5): `Σ²(L=n) − σ²(n) ≈ 1/6`. Kaynağı
   French–Mello–Pandey (kaynakça 7). Metin diyor: prensipte büyük `n` için
   geçerli, ama sayısal hesaplar ve analitik tahminler **küçük `n` için de
   `~0.01/n²` hatayla** tuttuğunu gösteriyor. Bölüm I'de ayrıca uyarı var:
   bu bağıntı Gauss ansambllerine, kaotik sistemlere (evrensel rejimin
   ötesinde bile) ve integrallenebilir spektrumlara uygulanır, **ama Poisson
   spektrumlarına uygulanmaz** — spektral katılık (rigidity) kavramına bağlı.
3. Denk. (6): `C(n) ≈ ½[Σ²(n+1) − 2Σ²(n) + Σ²(n−1)]`, `n ≥ 2`. **Denk. (3)'ün
   yapısı gereği (5)'teki sabitin (1/6) tam değeri `n ≥ 2` için önemsiz.**
   Yani `C(n)`, `Σ²`'nin **tamsayılarda alınmış ayrık eğriliğidir**.
4. Karşılaştırma: sürekli eğrilik `Y₂(L) = −½ ∂²Σ²/∂L²` (denk. 8). Yani
   `C(n) ↔ −Y₂` benzeşimi; makalenin bir amacı farklarını saptamak.

**İlk büyük yapısal sonuç (Bölüm II.A sonu, denk. 15–20):** `Σ²_rm` açılımındaki
alt-baskın salınımlı terimler (`cos(2πL)/(2πL)²` vb.) `C(n)`'de **tamsayı `L`'de
değerlendirildiği için donuyor** (0 ya da 1). Dolayısıyla `Y₂`'de ortalama-aralık
ölçeğinde salınım varken (denk. 20: `Y₂ = 1/2π²L² − cos(2πL)/2π²L²`, β=2),
**`C(n)` evrensel rejimde monoton, salınımsızdır.** RMT sonucu:

> `C⁰_rm(n) = (1/βπ²) log(1 − 1/n²)` — denk. (14)
> `C_rm(n) = C⁰_rm(n) + (1/βπ²)(λ_β/n⁴ + α_β/n⁶) + O(1/n⁸)` — denk. (16)
> `λ₁ = 3/2π²`, `α₁ = 15/6π² − 105/6π⁴` — denk. (17)
> `λ₂ = −3/2π²`, `α₂ = −15/6π² + 135/6π⁴` — denk. (18)
> `β=4` için `β=1`'e Mehta–Dyson eşlemesi: denk. (19).

Denk. (16) `n ≈ 2`'ye kadar iyi. Uyarı (denk. 21): Gauss ansambllerinde geçerli
toplam kuralı `Σ_j C(j) + σ²(1)/2 = 0` bu yaklaşımda **ihlal ediliyor**, çünkü
`j = 1, 2` terimleri kaba (özellikle `C(1)`, denk. (7)'den).

### 1.3 Difüzif sistemler (Bölüm II.B)

Kısa bölüm. Difüzif form faktörü `K_nu(τ) = (2/β)(τ_c/4π)^{d/2} τ^{1−d/2}`
(denk. 22, Altshuler–Shklovskii); `K = K_rm + K_nu` toplanabilirliği varsayılıyor.
Sonuç denk. (25):

> `C_dif(n) = [a_d d(d/2−1)/βπ²] (τ_c/4π)^{d/2} n^{−(2−d/2)}`,
> `a_d = π, π²/2, 4π²/3` (d=1,2,3).

`d=2`'de tam olarak sıfır; `d=1` negatif, `d=3` pozitif; her `d ≠ 2` için
`n^{−(2−d/2)}` ile söndüğünden `C_rm`'in `n^{−2}`'sinden **yavaş** ve kuyruğa
hâkim. `τ_c = 1/g` (boyutsuz iletkenlik) olduğundan bu bir sonlu-iletkenlik
düzeltmesi.

### 1.4 Kaotik sistemler (Bölüm II.C) — **köşegen (diagonal) yaklaşım**

Berry'nin (Nonlinearity **1** (1988) 399, kaynakça 8) yarı-klasik sayım
varyansından yola çıkıyor:

> `Σ²_po(n) = (4/βπ²) Σ_{rτ_p < τ_*} sin²(πnrτ_p)/(r²|det(M_p^r − 1)|)` — denk. (26)
> `Σ²_* = (2/βπ²)[Ci(2πnτ_*) − log(2πnτ_*) − γ]` — denk. (27)

Hannay–Ozorio de Almeida toplam kuralı sayesinde `C_po + C_*` toplamında `τ_*`
sonsuza uzatılabiliyor (kesim bağımsızlığı) ve `C_*`'ın limiti `C⁰_rm`'i tam
olarak götürüyor. Kalan:

> **`C(n) = (4/βπ²) Σ_{p,r} [sin²(πrτ_p)/(r²|det(M_p^r − 1)|)] cos(2πnrτ_p)
> + (1/βπ²)(λ_β/n⁴ + α_β/n⁶) + O(1/n⁸)`** — denk. (31)
> ve pratikte sadece ilk terim: **denk. (32)**.

**Varsayımlar, açıkça:** Bu **köşegen yaklaşımdır** — metin bunu her yerde
"diagonal approximation" diye adlandırıyor. Bölüm IV'ün ilk paragrafı köşegen-dışı
katkıların yerini de söylüyor: *"The off-diagonal semiclassical contributions are
here expected to reproduce the smooth higher-order corrections written in
Eq.(31)"* — yani denk. (31)'in `1/n⁴`, `1/n⁶` pürüzsüz kuyruğu. `R₂`'de
köşegen-dışıya bağlanan **analitik-olmayan salınımlı yapı `C(n)`'de yoktur.**
Bunun gerekçesi §1.2'deki tamsayı-donması. Bölüm IV: *"contrary to the
autocovariances of spacings, the diagonal approximation is a poorer approximation
for `R₂`."*

Sayısal reçete (dipnot 2): toplamı `rτ_p < τ_*`'de kes, kalan için `−C_*` ekle
(denk. 29).

### 1.5 Klasik zeta fonksiyonları ve rezonans formülü (Bölüm II.C, ikinci yarı)

Makalenin "ana sonuç" diye işaretlediği kısım. `sin²`'i seriye açıp
`τ_p^{2k} = r τ_p^{2k} − (r−1) τ_p^{2k}` ayrıştırmasıyla (denk. 33) toplam iki
klasik dinamik zeta fonksiyonuna dönüştürülüyor:

> `Z(s) = Π_p Π_{r≥1} exp(−e^{srT_p}/(r|det(M_p^r−1)|))` — denk. (35)
> `F(s) = Π_p Π_{r≥2} exp(−(r−1)e^{srT_p}/(r²|det(M_p^r−1)|))` — denk. (38)
> yaklaşık hali: denk. (40)
>
> **`C(n) = (2/βπ²) Σ_{k≥1} (1/(2k)!) Re ∂^{2k}/∂n^{2k} log[Z(in/ħρ̄)/F(in/ħρ̄)]`
> — denk. (39)**

`Z(s)`'in sıfırları Perron–Frobenius evrim operatörünün tayfı `γ_μ`; `F(s)`
makalenin kendi ifadesiyle **daha önce çalışılmamış** ve fiziksel yorumu açık
problem (Bölüm IV'te tekrar söyleniyor).

`Z, F`'yi tekil noktalarından çarpanlara ayırınca:

> `C(n) = (1/βπ²) Σ_sp sign(sp) log|1 − 1/(n + i n_sp)²|` — denk. (43),
> `n_μ = ħρ̄ γ_μ`.

Her tekil nokta `C(n)`'de bir tepe üretiyor; şekli denk. (42) ile veriliyor:
merkez `n ≈ n_μI`, yükseklik `H = (2βπ²)^{−1} log[1 + 2(n_μR² + 1)/n_μR⁴]`,
genişlik `W = 2n_μR`. Yani **hayali kısım konumu, reel kısım yükseklik ve
genişliği belirliyor.** İşaretler: `Z`'nin sıfırları ile `F`'nin kutupları aynı
işaret; `F`'nin sıfırları ile `Z`'nin kutupları ters işaret.

`T_H → ∞` limitinde `γ₀ = n₀ = 0` (ergodik sıfır) dışındaki bütün rezonanslar
sonsuza itiliyor → sadece evrensel denk. (14) kalıyor. Sonlu `T_H`'de
evrensel-olmayanlar katkı veriyor.

### 1.6 Riemann uygulaması (Bölüm III) — makalenin ağırlık merkezi

**Sözlük (Tablo I), tam olarak okunduğu gibi:**

| Nesne | Riemann karşılığı |
|---|---|
| periyodik yörünge etiketi | asal sayılar `p` |
| `ħ` | `→ 1` |
| simetri sınıfı | `β → 2` |
| asimptotik yoğunluk | `ρ̄ → log(t/2π)/(2π)` |
| Heisenberg zamanı | `T_H = hρ̄ → log(t/2π)` |
| eylem | `S_p → t log p` |
| periyot | `T_p → log p` |
| **ölçekli periyot** | **`τ_p → log p / log(t/2π)`** |
| Lyapunov üsteli | `λ_p → 1` |
| kararlılık çarpanı | `|det(M_p^r − 1)| → p^r` |
| dinamik zeta | `Z(s) → ζ^{−1}(1 − s)` |
| | `F(s) → Π_{r≥2}[ζ(r − rs)]^{(r−1)/r²}` |
| `Z`'nin "ergodik" sıfırı | `γ₀ →` `ζ(s)`'in kutbu |
| diğer sıfır/kutuplar | `→` `ζ(s)`'in sıfırları |

> ⚠️ **Kaynakta tutarsızlık:** Tablo I'deki `F(s)` üsteli `+(r−1)/r²`, oysa
> gövde metnindeki **denk. (46)** `F(s) = Π_{r≥2}[ζ(r − rs)]^{−(r−1)/r²}`
> (eksi üstel). Denk. (40) + `Π_p(1 − p^{rs−r}) = ζ^{−1}(r − rs)` zincirinden
> **denk. (46)'nın doğru, Tablo I girdisinin dizgi hatası olduğu** çıkıyor.
> Alıntı yaparken denk. (46)'yı kullanın.

**Ana formül:**

> **`C(n) = (2/π²) Σ_{p,r} [sin²(πrτ_p)/(r² p^r)] cos(2πnrτ_p)`,
> `τ_p = log p / log(t/2π)` — denk. (44)**

Bu bizim işaretlediğimiz nesnenin ta kendisidir; **birebir doğrulandı.**
(β=2 ⇒ `4/βπ² = 2/π²`; `|det| → p^r`.)

**Zeta diline çeviri:** `Z(s) = Π_p(1 − p^{s−1}) = ζ^{−1}(1 − s)` (denk. 45) —
yani Riemann durumunda **klasik dinamik zeta fonksiyonu, ötelenmiş Riemann
zeta'nın tersidir.** `F(s)` de ζ cinsinden (denk. 46). Bunları denk. (39)'a
koyunca denk. (47), sonra Hadamard çarpanlamasıyla (denk. 48, Titchmarsh) ve
**RH varsayımıyla** (`κ = 1/2 + i t_μ`) rezonans formülü:

> **denk. (49)** — üç grup terim (kutup / kritik sıfırlar / trivial sıfırlar),
> dış toplam tekrarlar `r` üzerinde, `r=1` teriminin işareti `r ≥ 2`'ninkine
> ters, global çarpan `(δ_{r,1} + 1 − r)/r²`.
>
> Rezonans konumları — **denk. (50)**:
> `n_{0r} = ρ̄(1 − 1/r)`
> **`n_{μr} = ρ̄[1 − (1/2 − i t_μ)/r]`**
> `n_{mr} = ρ̄[1 + 2(m+1)/r]`

**`r = 1` okuması (makalenin manşeti):**
- `n_{01} = 0` → evrensel terim denk. (14).
- `n_{μ1} = ρ̄(1/2 + i t_μ)` → **her kritik sıfır `C(n)`'de `n = ρ̄ t_μ`
  konumunda NEGATİF bir korelasyon tepesi üretiyor.** Hepsinin reel kısmı aynı
  (`ρ̄/2`) olduğu için **yükseklik ve genişlik sabit**.
- Odlyzko penceresi için `ρ̄ = 3.89533` ⇒ tepeler `n = 55.1, 81.9, 97.4, …`
  (`μ = 1,2,3`), `H = (4π²)^{−1} log[1 + 32(ρ̄²/4 + 1)/ρ̄⁴] = 0.013`, `W = ρ̄`.
- Trivial sıfırların katkısı küçük (hayali eksenden uzaklar).

**`r ≥ 2` okuması (resurgence):** Her kritik sıfır, ana (negatif) tepesinin
yanında `n = ρ̄t_μ/2, ρ̄t_μ/3, …` konumlarında **daha küçük POZİTİF alt-rezonanslar**
üretiyor; yükseklik `H = (r−1)(4π²r²)^{−1} log[1 + 2((Re n_{μr})² + 1)/(Re n_{μr})⁴]`.
Verilen somut örnek: `n ≈ 18.4` tepesi, ilk Riemann sıfırının (`t₁ ≈ 14.13`)
üçüncü mertebe alt-rezonansı `ρ̄t₁/3`.

**Rezonans örtüşme eşiği:** `n_c ≈ log(t/2π)·exp(2π)` — Odlyzko'nun `t`'si için
`≈ 13100`. Bunun ötesinde tepeler örtüşüyor ve `C(n)` düzensiz salınıyor.

**Lehmer olgusu:** İki sıfır çok yakınsa katkıları koherent toplanıp **iki katı**
tepe veriyor. Somut: `t = 7702/ρ̄ = 1977.24` civarındaki neredeyse-dejenere çift,
`n = 7702`'de `H ≈ 0.026` tepesi (Şekil 6, `Z(t) = e^{iθ(t)}ζ(1/2+it)`, denk. 51).

**"Neden Riemann'a özel?"** — Bölüm III sonu, üç sebep sıralıyor: (i) yoğunluk
için iz formülü **kesin**; (ii) yörünge periyotları `T_p = log p` **enerjiden
bağımsız**; (iii) monodromi matrisinin **tek genişleyen özdeğeri** var. Metnin
kendi ifadesiyle bu "genel bir dinamik sistem özelliğinden çok bir komplo gibi
görünüyor" ve daha fazla çalışma hak ediyor.

### 1.7 Odlyzko verisiyle karşılaştırma — **tam kapsam ve hassasiyet**

Bu, görevin özel sorusu; dikkatli okundu.

**Veri:** Şekil 2 altyazısı: Odlyzko'nun, `10¹²`-inci sıfır civarında
`t = 267653395648.8475`'ten başlayan ve **50 000 sıfır** içeren aralıktan
hesapladığı değerler. (`ρ̄ = 3.89533` bu pencereye karşılık geliyor.)

**Verinin gerçekten göründüğü yer sadece İKİ şekil:**

| Şekil | `n` aralığı | İçerik |
|---|---|---|
| **Şekil 2** | `n ≤ 20` | üçgen = RMT (denk. 16), **kare = Odlyzko verisi**, daire = teori (denk. 44) |
| **Şekil 3** | `9980 ≤ n ≤ 10000` | **kare = veri**, daire = teori |
| Şekil 4 (a,b,c) | üç geniş aralık | **yalnız teori** (daire, denk. 44) |
| Şekil 5 | Şekil 4(a)'nın sol kısmı, büyütülmüş | **yalnız teori**: daire = denk. (44), kare = denk. (49) `r=1`, elmas = denk. (49) `r≤3` |
| Şekil 6 | — | `Z(t)`, Lehmer çifti |

> **Kritik okuma:** Rezonans yapısı (alçak sıfır tepeleri, alt-rezonanslar,
> Lehmer çifti) **hiçbir şekilde nokta-nokta veriyle karşılaştırılmıyor.**
> Şekil 4 ve 5 tamamen teorik eğriler. Metin yalnızca sözel olarak "Eq.(44),
> which follows very closely the numerical results of Odlyzko" diyor.
> Veriyle doğrudan yüzleşme sadece iki adet 20 noktalık pencerede
> (`n ≤ 20` ve `9980 ≤ n ≤ 10000`).

**Hassasiyet:** Sayısal bir uyum ölçütü (χ², RMS, hata çubuğu) **verilmiyor.**
Verilen tek nicel çerçeve genlik mertebeleri:
- küçük `n` (Şekil 2): evrensel-olmayan dalgalanmaların tipik genliği `~10⁻³`;
  RMT'den farklar **zaten `n ≈ 3`'te** görünür hale geliyor.
- büyük `n` (Şekil 3): salınım genliği `~10⁻²`.
- gürültü tabanı: `s_m`'ler korelasyonsuz varsayılırsa sonlu-boyut istatistiksel
  dalgalanması `~10⁻⁴` (Odlyzko'ya atfen).
- Şekil 3 ayrıca ardışık noktalar arasında **işaret korelasyonu** gösteriyor.

Metnin kendi yargısı (Bölüm III ve IV): *"The agreement is good, although some
small deviations remain"* ve *"It is unclear whether the remaining differences
are due to the approximations made in our derivation or to numerical
inaccuracies. Probably to both. This deserves further study."*

**"~12 asal yeter" ifadesinin tam bağlamı:** Bölüm III, Şekil 4 tartışması.
`50 < n < 500` aralığında `C(n)` küçük pozitif bir korelasyon üstünde 3–4
aralık genişliğinde büyük anti-korelasyon tepeleriyle "noktalanmış" görünüyor.
Cümle: **denk. (44)'te yaklaşık 12 asal üzerinden (tekrarlarsız) toplayarak
Şekil 4(a)'daki ana rezonansların NİTEL özellikleri yeniden üretiliyor.**
Yani: (i) *nitel* bir ifade, nicel değil; (ii) `r = 1` (tekrarsız); (iii)
yalnız `50 < n < 500` penceresi için; (iv) veriye değil, denk. (44)'ün kendi
tam toplamına karşı bir sadeleştirme iddiası.

### 1.8 Sonuç bölümünün ekstra kazancı: `R₂` için aynı makine

Bölüm IV, `C(n)` makinesini iki-nokta fonksiyonuna taşıyor:

> `R₂(ε) − 1 = (4/β) Σ_{p,r} [τ_p²/|det(M_p^r−1)|] cos(2πεrτ_p) + (R₂^RMT(ε) − 1)
> + 1/(βπ²ε²)` — denk. (52)
> `R₂^diag(ε) − 1 = (1/βπ²) Re ∂²/∂ε² log[Z(iε/ħρ̄)/F(iε/ħρ̄)]` — denk. (53)

**`C(n)` ile `R₂` arasındaki tek yapısal fark, `sin²` çarpanıdır** — metnin
kendi ifadesiyle *"the presence of the square of the sine function in the
latter, which is due to the discrete nature of the curvature `C(n)`"*.
`R₂`'de `τ_p²`, `C(n)`'de `sin²(πrτ_p)/r²`.

`R₂^diag ↔ Z(s)` bağıntısının Agam–Altshuler–Andreev (PRL **75** (1995) 4389,
kaynakça 19) tarafından bulunduğu; BLS'in katkısının **ikinci zeta fonksiyonu
`F(s)`'i eklemek** olduğu açıkça söyleniyor. Ayrıca alçak sıfırların `R₂`'de
negatif rezonans olarak belirdiğinin **Berry–Keating (SIAM Review 41 (1999) 236,
kaynakça 20) tarafından zaten fark edildiği** iki ayrı yerde kabul ediliyor.

---

## 2. NESNE-NESNE SINIR TABLOSU

Hüküm dili: **VAR** = makalede aynı nesne var, yenilik iddiası ölü.
**AKRABA** = form/felsefe örtüşüyor, nesne farklı. **YOK** = makalede izi yok.

| # | Bizim nesne | BLS'te durumu | Referans | HÜKÜM |
|---|---|---|---|---|
| **C(n)** | `τ_p = log p/log(t/2π)` ve `C(n) = (2/π²)Σ_{p,r}[sin²(πrτ_p)/(r²p^r)]cos(2πnrτ_p)` | **Birebir aynı formül.** Tablo I'deki `τ_p` tanımı bizimkiyle harfi harfine aynı. | **denk. (44)** (Riemann), genel hali **denk. (32)**; sözlük **Tablo I** | **VAR — tam** |
| **(a)** | Tek asalın FAZINA KOŞULLU artık-varyans/kovaryans modülasyonu (`R_nn` tepki okuması) | **Hiçbir iz yok.** Tam metinde "conditional", "response", "susceptibility", "moment", "phase" kelimelerinin **hiçbiri geçmiyor** (kelime taraması yapıldı). Bütün nesneler koşulsuz ikinci momentler. Dahası: türetme baştan sona **köşegendir** ⇒ BLS çerçevesinde asal fazları bağımsızdır ⇒ `R_nn ≡ 0`; formül bu nesneyi **ilke olarak üretemez.** | — (yokluk); köşegenlik: **denk. (32)** ve **Bölüm IV, 1. paragraf** | **YOK** |
| **(b)** | Nokta-bazlı (örnek-yolu) özdeşlik `ds_n = Σ 2a_q sin(ω_q g/2) cos(ω_q m_n)` | **Özdeşliğin kendisi YOK** — makale tek bir `s_m`'i asallar cinsinden hiç yazmıyor, yalnız `⟨s_m s_{m+n}⟩` ile çalışıyor. **AMA bu özdeşliğin İKİNCİ MOMENTİ tam olarak denk. (44)'tür** (cebir aşağıda §3.4'te açık yazıldı; `2a_q² sin²(πτ_q)cos(2πnτ_q)` ⇔ `(2/π²)sin²(πrτ_p)cos(2πnrτ_p)/(r²p^r)`, `q = p^r`). | ikinci moment: **denk. (44)**; türetme yolu **denk. (2)–(6) + (26)** | **Özdeşlik YOK / ikinci momenti VAR** |
| **(c)** | Çift-korelasyonda alçak-sıfır rezonans ÇUKURLARININ gap-gecikme dilinde konumu `n* = γL/2π` | **BLS'in MANŞET SONUCU.** `ρ̄ = log(t/2π)/(2π) = L/2π` olduğundan onların `n = ρ̄ t_μ`'sü **bizim `n* = γL/2π`'mizin aynısıdır.** Türetilmiş (denk. 49–50), sayısallaştırılmış (55.1, 81.9, 97.4), şekli/yüksekliği/genişliği verilmiş (`H = 0.013`, `W = ρ̄`), örtüşme eşiği hesaplanmış (`n_c ≈ 13100`), üstüne **alt-rezonanslar** `ρ̄t_μ/r` ve **Lehmer ikilenmesi** eklenmiş. | **denk. (49), (50), (42)**; **Şekil 4, 5**; `R₂` versiyonunun önceliği Berry–Keating'e veriliyor (Bölüm III sonu + Bölüm IV) | **VAR — tam (teori)** |
| **(c′)** | Aynı çukurların **veriden ölçülmesi** | **YOK.** Şekil 4 ve 5 tamamen teorik. Odlyzko verisi yalnız `n ≤ 20` (Şekil 2) ve `9980 ≤ n ≤ 10000` (Şekil 3) pencerelerinde, toplam ~40 noktada görünüyor — **rezonans tepelerinin bulunduğu `n = 55–500` bandında hiç veri gösterilmiyor.** Nicel uyum ölçütü de hiç verilmemiş. | **Şekil 2–5 altyazıları**; Bölüm IV'ün "unclear whether the remaining differences…" cümlesi | **YOK — açık arazi** |
| **(d)** | Bond-adımı kompleks transferi / dispersiyonu `Θ(τ) = 2πτ + φ(τ)` | **Yok.** "transfer" kelimesi metinde geçmiyor; komşu aralıklar arası bir transfer katsayısı, kompleks kazanç ya da faz gecikmesi nesnesi hiç kurulmuyor. En yakın **görüntü** aldatıcı: denk. (42)'deki `n_μ` komplekstir — ama o, **klasik rezonansın gecikme düzlemindeki konumudur** (reel kısım = tepe genişliği), iki gap değişkeni arasındaki bir transfer fonksiyonu değil. Ayrıca BLS **tam olarak `n = 1`'de en zayıf**: `C(1)` ayrı bir denklemden (denk. 7) geliyor ve toplam kuralını (denk. 21) ihlal edecek kadar kaba olduğu **kendileri tarafından yazılıyor.** | yokluk; `n=1` zaafı: **denk. (7), (21)** ve onları izleyen paragraf | **YOK** (uyarı: yanlış-akraba denk. 42) |
| **(e)** | Negatif öz-yanıt (`−2` yasası) | **Yok.** Ne öz-yanıt ne koşullu bir ikinci moment var. **Karışma tehlikesi yüksek:** makalede iki ayrı "negatiflik" var — (i) RMT `C(n)`'in negatif ve monoton olması (`−1/βπ²n²`, denk. 14); (ii) alçak sıfır rezonanslarının **negatif** tepeler vermesi (işaret, denk. 49'da kutup-vs-sıfır yapısından geliyor). İkisi de **koşulsuz kovaryans işaretidir**, tepki katsayısı değil. | denk. (14); denk. (49)–(50) işaret tartışması | **YOK** (ama açıkça ayrıştırılmalı) |

### 2.1 Tabloya ek — iki tehlikeli sınır çizgisi

**(i) `sin²` vs `sin`.** BLS'te `sin(πrτ_p)` **yalnız karesi içinde** görünür.
Dolayısıyla `C(n)` — ve genel olarak herhangi bir koşulsuz ikinci moment —
`sin(πτ)`'nın **işaretine kördür.** Bizim yarım-gap fazı `sin(ω_q g/2)`
genlik düzeyinde lineerdir. Bu, (a)/(e)'nin ikinci momentlerden bağımsız
olduğunun **formül düzeyinde kanıtıdır** ve savunmada en sağlam noktamız.

**(ii) BLS'in "köşegen yeterli" iddiası bizi vurmuyor.** Metin, köşegen
yaklaşımın `C(n)` için çok iyi olduğunu ve köşegen-dışının yalnız pürüzsüz
`1/n⁴`, `1/n⁶` düzeltmeleri ürettiğini söylüyor. **Ama bu ifade koşulsuz
`C(n)` hakkındadır.** Koşullama farklı bir izdüşümdür: köşegen kısım
koşullamada **tam olarak sıfırlanır**, dolayısıyla koşullu sinyalin tamamı
köşegen-dışıdır ve BLS'in "küçüktür" tahmini o kanala uygulanmaz. Bu ayrımı
notta açıkça yazmak gerekir; yazılmazsa hakem "BLS köşegenin yettiğini
göstermiş" diye reddedebilir.

---

## 3. ALET NOTLARI — bizim işimize yarayacaklar

### 3.1 Sayım varyansı → aralık kovaryansı köprüsü (en değerli alet)

`C(n) ≈ ½Δ²Σ²(n)` (denk. 6), `Σ²(L=n) − σ²(n) ≈ 1/6` (denk. 5) üzerinden.

Bize kazandırdığı: **literatürdeki HER `Σ²(L)` ya da `K(τ)` öngörüsü, tek satırlık
ayrık ikinci farkla doğrudan gap-gecikme diline çevrilebilir.** Bizim bütün
gözlemlenebilirlerimiz gap tarafında, literatürün çoğu `Σ²`/`K` tarafında;
bu köprü ikisini parasız birleştiriyor.

Kullanım kuralları (hepsi metinden):
- `n ≥ 2` için sabitin (1/6) değeri **önemsiz** — Δ² onu yok ediyor.
- `n = 1` ayrı formül gerektiriyor (denk. 7) ve **güvenilmez**.
- Küçük `n` hatası `~0.01/n²`.
- **Poisson spektrumlarına uygulanamaz** (spektral katılık gerekiyor).

### 3.2 Ayrıklaştırmanın "salınım öldürme" teoremi

`Σ²` açılımındaki `cos(2πL)`, `sin(2πL)` terimleri tamsayı `L`'de donuyor ⇒
`C(n)` evrensel rejimde monoton, `Y₂`'nin ortalama-aralık ölçekli salınımı
yok (denk. 15 vs denk. 20).

Bizim için iki sonucu var:
1. Gap-gecikme gözlemlenebilirleri **ortalama-aralık ölçeğindeki salınımlardan
   yapısal olarak arındırılmış** gelir — `R₂` ile uğraşırken çekilen dert bizde
   yok.
2. Ters yönde uyarı: bu yüzden gap-gecikme dilinde ölçtüğümüz *herhangi bir*
   `2π`-periyotlu salınımın kaynağı **`Σ²`'nin düz açılımı olamaz**; başka bir
   yerden gelmek zorundadır. Sinyalin adresini daraltan ücretsiz bir filtre.

### 3.3 Sonlu-boyut / pencere düzeltmeleri

- **RMT kuyruğu:** `C_rm(n) = C⁰_rm(n) + (1/βπ²)(λ_β/n⁴ + α_β/n⁶) + O(1/n⁸)`,
  β=2 için `λ₂ = −3/2π²`, `α₂ = −15/6π² + 135/6π⁴` (denk. 16, 18). `n ≈ 2`'ye
  kadar iyi. **Hazır taban çizgisi** — kendi RMT null'umuzu üretirken doğrudan
  kullanılabilir.
- **Toplam kuralı tanısı:** `Σ_j C(j) + σ²(1)/2 = 0` (denk. 21, Pandey; Gauss
  ansambllerinde geçerli). BLS'in kendi yaklaşımı bunu ihlal ediyor.
  Bizim sayısal `C(n)` boru hattımız için **ucuz ve keskin bir sağlama**.
- **Toplam kesme reçetesi:** `rτ_p < τ_*`'de kes, `−C_*` ekle (denk. 29, dipnot 2).
- **Pencere/örtüşme sınırı:** `n_c ≈ log(t/2π)·e^{2π}`. Odlyzko `t`'si için
  `≈ 13100`; bizim `t`'lerimiz için yeniden hesaplanmalı. Gecikmeyi bunun
  ötesine iterek "yapı yok" demek anlamsız — orada zaten örtüşme rejimi var.
- **Rezonans şekli:** denk. (42) — merkez `Im n_μ`, genişlik `2 Re n_μ`,
  yükseklik `Re n_μ`'den. Riemann'da bütün kritik-sıfır rezonanslarının
  `Re n_{μ1} = ρ̄/2` olması ⇒ **sabit yükseklik/genişlik**. Ölçtüğümüz herhangi
  bir tepe ailesi için hazır bir kalıp.

### 3.4 Unfolding konvansiyonu — **bizimkiyle fark YOK**

- Unfolding: `x_i = ∫^{t_i} ρ̄(t)dt`, pencere `Δt << t` (Bölüm II.A).
- `ρ̄ → log(t/2π)/(2π)`, `T_H → log(t/2π)` (Tablo I).
- **`τ_p = log p / log(t/2π)`** — yani BLS'in paydası **`log(t/2π)`**, `log T`
  değil. Bizim `L = log(t/2π)` ile **birebir aynı.** Bu soruda risk yok.
- Tekrarlar: BLS `τ_p` asala göre tanımlıyor, kuvvetler `rτ_p` olarak giriyor.
  Bizim `τ_q = log q/L` (`q = p^m`) tanımı `m τ_p`'ye eşit — **aynı değişken.**
- Ölçek çevrimi: `n` katlanmış (unfolded) gecikme; `t` birimine geçiş için
  ortalama aralık `2π/L`. Dolayısıyla `ω_q · (m_{n+k} − m_n) = 2πkτ_q`.
- Normalizasyon uyarısı: BLS'in `C(n)`'i **birim ortalamalı** `s`'lerin
  kovaryansıdır (korelasyon katsayısı değil). Sayı karşılaştırırken
  `σ²`'ye bölünüp bölünmediğine dikkat.

**Cebirsel doğrulama (bizim özdeşliğin ikinci momenti = denk. 44).** Elle
yapıldı, geçti:
`ds_n = Σ_q 2a_q sin(ω_q g/2) cos(ω_q m_n)`, `g ≈ 2π/L`, `a_q = 1/(πm p^{m/2})`
⇒ `sin(ω_q g/2) = sin(πτ_q)`, `⟨cos(ω_q m_n)cos(ω_q m_{n+k})⟩ = ½cos(2πkτ_q)`
⇒ `Cov = Σ_q 2a_q² sin²(πτ_q) cos(2πkτ_q) = (2/π²) Σ_{p,m} sin²(πmτ_p) cos(2πkmτ_p)/(m²p^m)`
= **denk. (44)**, `r ↔ m`.

Aynı `sin²` çarpanı BLS'te bambaşka bir yerden geliyor: `Σ²_po`'daki
`sin²(πnrτ_p)` üzerine ayrık ikinci fark uygulanınca
`Δ²[cos(2πnrτ)] = −4sin²(πrτ)cos(2πnrτ)` özdeşliğiyle doğuyor. **Yani BLS'te
`sin²(πrτ_p)` bir AYRIKLAŞTIRMA artifaktı, bizde bir YARIM-GAP FAZIDIR.**
Aynı çarpanın iki farklı okuması — bu, ucuz ve dürüst bir yorumsal katkıdır.

### 3.5 `R₂` tarafı için hazır formül

denk. (52): `R₂(ε) − 1` köşegen kısmı `Σ_{p,r} τ_p²/|det| cos(2πεrτ_p)`.
Riemann sözlüğüyle `(4/2)(log²p/L²)/p^r = 2log²p/(L²p^r)` — bu, önceki
keşif raporunda Berry–Keating (4.20) diye işaretlediğimiz yapının **BLS
notasyonundaki tam karşılığıdır.** İki kaynağı birbirine bağlamak için
kullanışlı; `R_nn`'in `cos(πτ)` şeklinin literatürdeki adresi olarak
**BK (4.20) ∪ BLS denk. (52)** diye yazılabilir.

---

## 4. ARXİV NOTUMUZ İÇİN KONUMLANDIRMA CÜMLELERİ

### 4.1 Atıfla teslim ettiklerimiz (yenilik iddia ETMEYECEĞİZ)

Bunlar tartışmasız BLS'indir; notun ilgili yerlerinde **öne çıkarılarak**
teslim edilmeli, dipnotta değil:

1. **`τ_p = log p/log(t/2π)` değişkeni ve gap-gecikme öz-kovaryansının
   asal-toplamı.** Taslak cümle:
   > "The lag-`n` autocovariance of consecutive gaps was expressed as a sum
   > over primes by Bohigas, Lebœuf and Sánchez [BLS], Eq. (44), in terms of
   > the rescaled period `τ_p = log p/log(t/2π)`; our second moment reproduces
   > their formula identically, and we adopt their normalization throughout."

2. **Alçak Riemann sıfırlarının `n = ρ̄ t_μ` konumunda rezonans olarak
   belirmesi** (ve alt-rezonanslar, Lehmer ikilenmesi, örtüşme eşiği).
   Taslak cümle:
   > "That the low-lying zeros reappear as resonances in the correlations of
   > zeros high on the critical line is not new: the effect was noticed in the
   > two-point function by Berry and Keating, and was derived in the
   > gap-lag language — with explicit positions `n = ρ̄ t_μ`, `ρ̄ = L/2π`,
   > constant peak height and width, sub-resonances at `ρ̄ t_μ/r`, and a
   > resonance-overlap threshold `n_c ≍ L e^{2π}` — by BLS, Eqs. (49)–(50)."

   ⚠️ **Bu, notun en tehlikeli yeridir.** `n* = γL/2π` bizim değil, BLS'indir.
   İddiayı buradan çekmek zorundayız.

3. **Aralık-kovaryansının sayım varyansından türetilmesi** (denk. 2–7) ve
   **RMT taban çizgisi** (denk. 14, 16–18). Alet olarak atıfla kullanılır.

4. **Köşegen yaklaşımın `C(n)` için `R₂`'den iyi çalışması** ve ayrıklaştırmanın
   Heisenberg-ölçekli yapıyı öldürmesi (Bölüm II.A, IV).

### 4.2 Yeni iddia edebileceklerimiz (ve iddianın tam ifadesi)

**Y1 — Koşullu / faz-çözünürlüklü istatistikler (`R_nn`).**
> "Every prime-resolved statement in the literature we are aware of — Berry's
> number variance, BLS's Eq. (44), Berry–Keating's `R_c¹` — is an *unconditional
> second moment*, in which the prime phase enters only through `sin²(πrτ_p)`
> and `cos(2πnrτ_p)`. Such expressions are, by construction, blind to the sign
> of `sin(πτ)`. The object measured here is different in kind: the residual
> variance of the gap sequence *conditioned on the phase of a single prime
> wave*. In the diagonal approximation — which is the approximation underlying
> Eq. (44) — distinct prime phases are independent, so this quantity vanishes
> identically. A nonzero measurement is therefore a direct, phase-resolved
> probe of the off-diagonal (prime–prime) correlations."

Ekle (savunma): BLS'in "köşegen yeterli" ifadesi koşulsuz `C(n)` içindir;
koşullu kanalda köşegen kısım tam olarak sıfırlandığı için o tahmin bu kanala
taşınmaz.

**Y2 — Örnek-yolu (nokta-bazlı) özdeşlik ve ölçülen koherans.**
> "BLS derive Eq. (44) from Berry's number variance by taking a discrete second
> difference; the factor `sin²(πrτ_p)` there is an artifact of the
> discretization. We obtain the same second moment from a *sample-path*
> identity for the individual gap, `ds_n = Σ_q 2a_q sin(ω_q g/2) cos(ω_q m_n)`,
> in which the same factor appears as a *half-gap phase* and, crucially, at
> linear order in the amplitude. The identity itself — and the pointwise
> coherence it achieves gap by gap — does not appear in the literature; what is
> known is its second moment."

**Y3 — Bond-adımı kompleks transferi / dispersiyonu.**
> "BLS's framework is weakest exactly at `n = 1`: `C(1)` requires a separate
> equation (their Eq. (7)) and is inaccurate enough that the Gaussian-ensemble
> sum rule (their Eq. (21)) is violated — a limitation they state explicitly.
> The adjacent-gap channel is therefore open ground. The object we measure
> there is not a covariance at all but a *complex transfer coefficient* with a
> nontrivial phase, `Θ(τ) = 2πτ + φ(τ)`."

**Y4 — Negatif öz-yanıt (`−2`).**
Notta **açık bir ayrım paragrafı** şart, yoksa BLS'in negatif tepeleriyle
karıştırılır:
> "Two distinct negativities must not be conflated. In BLS, `C(n)` is negative
> in the universal regime (`−1/βπ²n²`) and each critical zero contributes a
> *negative* resonance peak; both are signs of an unconditional covariance,
> fixed by the pole-versus-zero structure of Eq. (49). The quantity reported
> here is a *response coefficient* of a conditional second moment, and its
> negativity is a statement about how the gas reacts to its own prime field —
> a quantity that is identically zero in the diagonal theory that produces
> those peaks."

**Y5 — Veri tarafı: rezonans çukurlarının doğrudan ölçümü.**
Bu, tam metin okumasının **beklenmedik kazancı**:
> "Although BLS derived the resonance structure in full, their comparison with
> Odlyzko's data is confined to two windows, `n ≤ 20` and
> `9980 ≤ n ≤ 10000` (their Figs. 2 and 3); Figs. 4 and 5, which display the
> resonances at `n = ρ̄ t_μ`, are entirely theoretical, and no quantitative
> measure of agreement is reported anywhere in the paper. A direct measurement
> of the low-zero resonances in the band `50 ≲ n ≲ 500`, with an error budget,
> therefore appears not to have been carried out."

⚠️ Bu cümleyi kullanmadan önce **1987–2026 arası literatürde bu ölçümün
yapılmadığı ayrıca taranmalı.** BLS'te yapılmamış olması, 25 yılda başkasınca
yapılmadığı anlamına gelmez. (Bu okuma yalnız BLS'i kapsıyordu.)

### 4.3 Notun anlatı sırası önerisi

1. BLS denk. (44) ile başla, **bizim özdeşliğimizin ikinci momenti olarak**
   teslim et (güven inşa eder, hakemi silahsızlandırır).
2. `sin²` → `sin` geçişini yap: "ikinci moment işarete kör" argümanı.
3. Koşullu kanalı tanımla, köşegende özdeş sıfır olduğunu göster.
4. Ölçümü ver.
5. `n = ρ̄t_μ` rezonanslarını **BLS'e atıfla** anlat, kendi ölçümümüzü onların
   teorik eğrisine karşı konumlandır.

---

## 5. OKUMANIN SINIRLARI (dürüstlük kaydı)

- **Şekillerin görselleri okunmadı** — yalnız altyazıları. Şekil 2/3'teki
  veri-teori uyumunun göz kararı kalitesi hakkında hüküm veremem; metnin
  "agreement is good" ifadesine dayanıyorum.
- **Dergi sayfa numaraları yok** (ar5iv HTML'inde pagination yok). Bu raporda
  hiçbir sayfa numarası verilmedi; hepsi denklem/bölüm/şekil referansı.
  Dergi sürümüne atıf yapacaksak sayfa numaraları PDF'ten teyit edilmeli.
- **Denk. (46) vs Tablo I üstel işareti tutarsız** (§1.6). Denk. (46)'nın
  doğru olduğu türetmeden çıkarıldı, ama PDF'ten teyit edilmesi iyi olur.
- ar5iv/LaTeXML dönüşümünde birkaç `\hbox{...}` kaçağı görüldü (`τ ≲ τ_e`
  gibi eşitsizlik işaretleri). Bunlar okumayı etkilemedi ama **bu rapordan
  formül kopyalayıp doğrudan TeX'e yapıştırmayın**; denklemleri PDF'ten
  yeniden alın.
- **Denk. (50) numarasına dikkat.** Rezonans konumları (`n_{0r}`, `n_{μr}`,
  `n_{mr}`) üç satırlık bir grup ve HTML'de "(50)" etiketini taşıyor; ama
  gövde metni bunlara iki yerde **"Eq.(III)"** diye atıf yapıyor — ar5iv'in
  çözemediği bir çapraz-referans. Dergi sürümünde numara farklı olabilir;
  atıf yapmadan önce PDF'ten teyit edin.
- Dönüşümün eksiksizliği **denetlendi**: numaralı denklem bloklarının
  1–53 arası tamamı mevcut, hiç eksik yok; bölüm zinciri (Özet → I → II.1/2/3
  → III → IV → Teşekkür → 25 kaynak → Tablo I → Şekil 1–6 altyazıları) tam.
  §2'deki "kelime metinde geçmiyor" hükümleri kelime taramasıyla doğrulandı
  (`phase`, `moment`, `response`, `conditional`, `transfer`, `susceptibility`
  = 0 kez; kontrol olarak `resonance` = 29, `diagonal` = 12).
- Kaynakça 25 künye olarak okundu; hiçbirinin içeriği bu oturumda **açılmadı**.
