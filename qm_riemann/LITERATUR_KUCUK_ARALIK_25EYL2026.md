# Küçük-Aralık (δ_n, M_n) Ortak Yasası — Literatür + Veri Envanteri

**Tarih:** 25 Eylül 2026
**Bağlam:** Not 1 (`arxiv_gap_amplitude.tex`) Riemann sıfırları için (δ_n, M_n) ortak
yasasını ilk kez sayısal olarak ölçtü ve "Açık problemler" bölümünde (1) maddesi
olarak *gap-koşullu* E[M_n² | δ_n] türetimini açıkça bir açık soru olarak
bıraktı. Not 1'in bir okurunun önerisi ve kaptanın veri-öncesi fikri (küçük δ'da
M_n ≈ |Z''|δ_n²/8; CUE'de iki özdeğer θ'da birleşirken kalan N−2'nin ölçüsü
|Λ_{N−2}(θ)|⁴ ile eğilir ⇒ koşullu M/δ² yasası = |Λ|⁴-eğilmiş Keating–Snaith
dağılımı) bu belgede sınanıyor.

**Yöntem notu:** Tüm atıflar arXiv abs sayfası, Crossref API veya dergi
sayfasından doğrudan doğrulandı (WebFetch/WebSearch ile bu oturumda). Yalnız
ikincil kaynaktan (başka bir makalenin atıfı, ResearchGate özeti) teyit
edilebilenler "DOĞRULANAMADI (yalnız ikincil kaynak)" olarak işaretlendi.
KÖRLÜK KURALINA UYULDU: bu belgenin hiçbir yerinde M_n/δ_n² dağılımı,
momentleri veya korelasyonu hesaplanmadı/yazdırılmadı/incelenmedi. Veri
bölümünde yalnız dosya yapısı (anahtar/boyut/t-aralığı) ve unfolded-δ
SAYIMLARI (δ̃ eşiklerinin altında kaç aralık) verilmiştir.

---

## 1. Literatürde bu (veya çok yakın bir) sonuç var mı?

### Kısa cevap
**Hayır — tam bu sonucu (küçük-δ koşullu M/δ² yasasının |Λ_{N−2}|⁴-eğilmiş
Keating–Snaith dağılımı olduğu, ve zeta tarafında aritmetik çarpanlı hâli)
literatürde bulamadım.** Bu, hem Not 1'in kendi "Açık problemler" listesiyle
(madde 1 ve 3), hem de projenin 16 Ağustos 2026 tarihli önceki literatür
taramasıyla (`LITERATUR_TARAMASI_16AGU2026.md`, satır 33-34, 54, 673:
"Bizim ölçtüğümüz şey ... Conrey–Ghosh'un gap-koşullu rafinesidir ...
Lehmer literatürü nitel kalıyor") bağımsız olarak örtüşüyor. Aşağıda her
aday kaynak tek tek değerlendiriliyor.

### 1.1 Lehmer çiftleri ve Z″/ξ″ ilişkisi (nitel önceki)

**Csordas, G.; Smith, W.; Varga, R. S.**, *Lehmer pairs of zeros, the
de Bruijn–Newman constant Λ, and the Riemann Hypothesis*, Constructive
Approximation **10** (1994), 107–129. DOI: `10.1007/BF01205170`.
DOĞRULANDI (Crossref).
— Lehmer çifti tanımını ve Λ ile bağını verir: küçük aralık ⇒ ξ (ya da H₀)
için bir alt-sınır/itme koşulu. **Sonuç:** nitel — "küçük aralık küçük ara
tepe zorlar" folklor ilişkisinin kaynağı (Not 1'in `Leh56` atfıyla aynı
aile). M/δ² için nicel bir dağılım/momet vermez.

**Stopple, J.**, *Lehmer pairs revisited*, arXiv:1508.05870 (2015);
Experimental Mathematics'te kabul edildi. DOĞRULANDI (arXiv abs).
— "Güçlü Lehmer çifti" tanımını Ξ(t)'nin pre-Schwarzian türevi üzerinden
verir (Teorem 2-3: ζ′(ρ) ve komşu ζ′ sıfırları cinsinden). t=10⁶'da 114.661
çiftten 855 güçlü Lehmer çifti bulunmuş. **Sonuç:** yine nitel/eşitsizlik
düzeyinde; M_n'nin kendisini veya dağılımını hesaplamıyor, doğrudan Z″ veya
M/δ² momenti içermiyor.

**Rodgers, B.; Tao, T.**, *The de Bruijn–Newman constant is non-negative*,
Forum of Mathematics, Pi **8** (2020), e6. DOI: `10.1017/fmp.2020.6`.
DOĞRULANDI (Crossref).
— Newman sanısını kanıtlıyor (Λ≥0); yerel sıfır-itme argümanları kullanıyor
ama küçük-aralık koşullu M momentiyle ilgisi yok. **Sonuç:** İLGİSİZ (yalnız
tarihsel/motivasyonel bağlam için faydalı).

### 1.2 "Small gaps" + ζ′'nin kritik çizgiye yakınlığı

**Farmer, D. W.; Gonek, S. M.; Lee, Y.**, *Pair correlation of the zeros of
the derivative of the Riemann ξ-function*, J. London Math. Soc. (2) **90**
(2014), no. 1, 241–269. DOI: `10.1112/jlms/jdu026`. DOĞRULANDI (Crossref).
— Ana sonuç: ζ sıfırları arasında pozitif oranda küçük aralık ⟺ ζ′(s)
sıfırlarının pozitif oranda kritik çizgiye çok yakın olması (Farmer–Ki
sanısının güçlü biçimde çözümü). **Sonuç:** yapısal/nitel eşdeğerlik —
küçük δ ile "yerel düzlük" arasındaki bağı teyit eder ama M_n'nin kendisini
veya bir M/δ² yasasını türetmez.

**"Farmer–Gonek–Lee–Yang" (görevde anılan 4 yazarlı biçim): DOĞRULANAMADI.**
arXiv/Crossref taramasında bu 4 yazarlı bir makale bulunamadı; en yakın
gerçek eşleşme yukarıdaki 3-yazarlı Farmer–Gonek–Lee (2014) makalesidir.
Muhtemelen bir karıştırma (belki Farmer–Gonek–Hughes, J. Reine Angew. Math.
**609** (2007), 215–236, "The maximum size of L-functions" ile, ya da
başka bir kombinasyonla) — bu ikincisi de küçük-aralık koşullu M momenti
içermiyor, yalnız ikincil kaynaktan (arama sonucu özetinden) teyit edildi,
birincil metne bakılmadı; gerekirse ayrıca doğrulanmalı.

### 1.3 C. Hughes'in 2001 tezi ve "ortak momentler" ailesi — EN YAKIN AKRABA, AMA FARKLI SORU

**Hughes, C. P.**, PhD tezi, *On the characteristic polynomial of a random
unitary matrix and the Riemann zeta function*, University of Bristol, 2001.
DOĞRULANAMADI (yalnız ikincil kaynaklardan — ResearchGate, IMRN makalesinin
girişi — teyit edildi; tez arXiv'de değil, Bristol kurumsal deposuna
birincil erişim bu oturumda denenmedi). Tezin içeriği, onu kanıtlayan
makalelerin girişlerinden güvenilir biçimde yeniden kurulabiliyor:

> **Hughes'in sanısı** (gerçek üs h,s için): E[|Λ_N(θ)|^{2s}|Λ_N′(θ)|^{2h}]
> asimptotiğinin genel gerçek s,h için var olduğu ve kapalı formu — tamsayı
> s,h için tez kendisi kanıtlamış, F(s,h) katsayısını kombinatorik bir
> toplamla vermiş.

Bu sanı şurada kanıtlandı:
**Assiotis, T.; Keating, J. P.; Warren, J.**, *On the joint moments of the
characteristic polynomials of random unitary matrices*, Int. Math. Res.
Not. **2022** (2022), no. 18, 14564–14603. DOI: `10.1093/imrn/rnab336`.
DOĞRULANDI (Crossref). arXiv:2005.13961.

İlgili aile: **Dehaye, P.-O.**, *Joint moments of derivatives of
characteristic polynomials*, Algebra & Number Theory **2** (2008), no. 1,
31–68. DOI: `10.2140/ant.2008.2.31`. DOĞRULANDI (Crossref). — tamsayı
üslerde |Λ|^{2k}|Λ′|^{2h} ortak momentlerini rasyonel fonksiyon olarak
verir; **makalenin kendisi zeta ve Hardy'nin Z-fonksiyonu için benzer
sanılar önerdiğini** ikincil kaynaklardan görüyoruz (birincil metin bu
oturumda tam okunmadı — teyit gerekiyor).
Ayrıca **Conrey, J. B.; Rubinstein, M. O.; Snaith, N. C.**, *Moments of the
derivative of characteristic polynomials with an application to the
Riemann zeta function*, Comm. Math. Phys. **267** (2006), no. 3, 611–629 —
karakteristik polinom ile türevinin karışık momentini Painlevé III'e
bağlıyor. DOĞRULANAMADI (yalnız ikincil kaynak — başlık/dergi/cilt arama
sonuçlarından toplandı, Crossref ile bire bir teyit edilmedi).

**Neden bu farklı bir soru:** Hughes/Dehaye/CRS06/AKW22 ailesi, karakteristik
polinomun ve **kendi türevinin AYNI NOKTADA** (tek bir θ'da, basit bir
sıfırın yakınında) ortak momentini soruyor — zeta tarafındaki analoğu
Z(γ_n), Z′(γ_n) momentleridir (Hughes–Keating–O'Connell 2000, Proc. R. Soc.
A **456** (2000), 2611–2627 — DOĞRULANAMADI, ikincil kaynak). Kaptanın
sorusu ise **İKİ AYRI özdeğerin (θ₁,θ₂) BİRLEŞMESİ** durumunda geri kalan
N−2 özdeğerin polinomunun θ'daki değeridir — yani gap'in KENDİSİ sıfıra
giderken kalan polinomun davranışı. Matematiksel olarak ilişkili (bir
sıfırın "ikizlenmesi" limitinde |Λ_N(θ)| ~ (θ−θ₁)(θ−θ₂)|Λ_{N−2}(θ)| olur ve
CRS06/Dehaye tipi türev-momentleri bu limitte kaptanın |Λ_{N−2}|⁴ tiltine
bağlanabilir) ama **literatürde bu köprü açıkça kurulmuş görünmüyor.**

### 1.4 CUE'yi noktalara koşullamak / eğmek (RMT tarafının teknik akrabası)

**Charlier, C.; Claeys, T.**, *Thinning and conditioning of the Circular
Unitary Ensemble*, arXiv:1604.08399 (2016). DOĞRULANDI (arXiv abs).
— CUE'yi bir yay üzerinde özdeğer OLMAYACAK şekilde koşullandırıyor (boşluk
olasılıkları, Toeplitz determinantı / Fisher–Hartwig tekillikleri ile).
Teknik olarak "koşullama + tiltleme" aynı ailede (FH tekillikleri = eğilmiş
ölçü) ama **soru tam tersi**: burada bir bölgede sıfır özdeğer olasılığı
soruluyor, kaptanın sorduğu "iki özdeğer TAM O NOKTADA birleşsin, geri kalanı
nasıl dağılır" değil. **Sonuç:** komşu teknik altyapı (FH/Toeplitz), doğrudan
sonuç değil.

Kaptanın Vandermonde/Coulomb-gaz argümanı (iki koşullama noktası θ→aynı yere
gelince kalan N−2'nin ölçüsü ∏|θ_i−θ|⁴ ile eğilir) **kendi başına RMT'de
temel/standart bir gözlemdir** (determinantal süreçlerin Palm dağılımının
yine determinantal olması ve koşullama noktası başına |Λ|² kazanılması genel
bir olgu — bkz. Forrester'ın *Log-Gases and Random Matrices* kitabındaki
Coulomb-gaz formalizmi, ikincil kaynak, birincil sayfa numarası bu oturumda
doğrulanmadı). Standart olması, sonucun bilinmediği anlamına gelmiyor — ama
bunu **açıkça küçük-δ zeta/CUE gap-max yasasına uygulayan bir kaynak
bulunamadı.**

### 1.5 CG85 / HLPC24 / PC24 — Not 1'in zaten kullandığı üçlü (yeniden teyit + "koşullu mu?" kontrolü)

**Conrey, J. B.; Ghosh, A.**, *A mean value theorem for the Riemann
zeta-function at its relative extrema on the critical line*, J. London
Math. Soc. (2) **32** (1985), 193–202. DOI: `10.1112/jlms/s2-32.2.193`.
DOĞRULANDI (Crossref) — Not 1'in `CG85` atfıyla birebir örtüşüyor.
**Yalnız marjinal** (koşulsuz) mean(M_n²) ~ (e²−5)/2·log T veriyor.

**Hughes, C.; Lugmayer, S.; Pearce-Crump, A.**, *The second moment of the
Riemann zeta function at its local extrema*, arXiv:2411.05573 (2024); J.
London Math. Soc. (2) **112** (2025), e70250. DOĞRULANDI (arXiv abs+PDF
metadata okundu). **Özet/PDF metaverisinde "small gaps", "conditioned on
the gap", "tilted", "characteristic polynomial" ifadeleri YOK** — bu makale
de yalnız marjinal (δ'ya göre integrallenmiş) alt-mertebe terimleri veriyor,
Not 1'in zaten kullandığı biçimde.

**Pearce-Crump, A.**, *Moments of the Riemann zeta function at its local
extrema*, arXiv:2411.05568 (2024); Mathematika **71** (2025), e70035.
DOĞRULANDI (arXiv abs+PDF metadata okundu). Aynı şekilde: gap-koşullu ifade
**YOK**, marjinal ilk moment (türevlerle) veriyor.

**Sonuç:** Not 1'in zaten dayandığı bu iki en-güncel makale bile, kendi
metinlerinde küçük-aralık koşullu bir analiz içermiyor — bu, sorunun güncel
literatürün (2024-2025) sınırında bile açık kaldığının bağımsız bir işareti.

### 1.6 Özet tablo

| Aday | Durum | Verdiği | Kaptanın sorusuyla ilişki |
|---|---|---|---|
| Csordas–Smith–Varga 1994 | DOĞRULANDI | Lehmer çifti ↔ Λ alt sınırı | nitel, aynı folklor ailesi |
| Stopple 2015 | DOĞRULANDI | güçlü Lehmer çifti tanımı (Ξ pre-Schwarzian) | nitel |
| Rodgers–Tao 2020 | DOĞRULANDI | Λ≥0 kanıtı | ilgisiz (bağlam) |
| Farmer–Gonek–Lee 2014 | DOĞRULANDI | küçük-δ ⟺ ζ′ çizgiye yakın | yapısal eşdeğerlik, nicel değil |
| "Farmer–Gonek–Lee–Yang" | DOĞRULANAMADI | — | muhtemel karıştırma |
| Hughes tezi 2001 | DOĞRULANAMADI (ikincil) | Λ,Λ′ ortak moment sanısı | **en yakın akraba, farklı soru** (aynı nokta, iki nokta değil) |
| Assiotis–Keating–Warren 2022 | DOĞRULANDI | Hughes sanısını kanıtlıyor | aynı — tek nokta |
| Dehaye 2008 | DOĞRULANDI | tamsayı Λ,Λ′ ortak momentleri; zeta/Z için sanı öneriyor (ikincil kaynak) | akraba, teyit edilmeli |
| Conrey–Rubinstein–Snaith 2006 | DOĞRULANAMADI (ikincil) | Λ,Λ′ karışık momenti ↔ Painlevé III | akraba, tek nokta |
| Charlier–Claeys 2016 | DOĞRULANDI | CUE'yi bir yayda koşullama (boşluk olasılığı) | komşu teknik, ters soru |
| Conrey–Ghosh 1985 | DOĞRULANDI (Not 1'de zaten var) | marjinal mean(M²) | koşulsuz |
| Hughes–Lugmayer–Pearce-Crump 2024 | DOĞRULANDI | marjinal alt-mertebe terimler | **küçük-aralık YOK** (kontrol edildi) |
| Pearce-Crump 2024 | DOĞRULANDI | marjinal ilk moment+türevler | **küçük-aralık YOK** (kontrol edildi) |
| Fyodorov–Gnutzmann–Keating 2018 | Not 1'de zaten var | global max BÜYÜK aralıkta oturur | ters rejim (büyük δ, global max) |
| Durkan–Hughes–Pearce-Crump 2026 (arXiv:2601.18025) | DOĞRULANDI (proje içi önceki oturumdan + bu oturumda arXiv abs) | Landau–Gonek genellemesi, zeta ortalamaları | ilgisiz görünüyor (küçük-aralık dili yok) |

---

## 2. Birim/normalizasyon

- **Z(t) ↔ |ζ(1/2+it)|:** Tanım gereği |Z(t)| = |ζ(1/2+it)| (Z(t) =
  e^{iθ(t)}ζ(1/2+it), θ = Riemann–Siegel theta fonksiyonu) — standart, ders
  kitabı düzeyinde bir kimlik (Titchmarsh, *The Theory of the Riemann
  Zeta-Function*, 2. baskı, Bölüm 4.17; bu oturumda birincil sayfaya
  bakılmadı ama kimlik tartışmasız).
- **Unfolding (δ̃):** Not 1'in kendi tanımı — δ̃_n = δ_n·L/2π, L=log(t/2π) —
  aralık envanterinde de aynen bu tanım kullanıldı (bkz. §3). mean(δ̃)=1
  hedeflenir.
- **CUE θ birimleri:** özfazlar θ_j∈[0,2π), ortalama aralık 2π/N; unfold
  edilmiş CUE aralığı (N/2π)·(gerçek aralık) — standart Mehta/Keating–Snaith
  kuralı.
- **M ve δ'yı karşılaştırılabilir kılan normalizasyon:** Not 1 zaten bunu
  çözmüş durumda — M̃_n = M_n/√(AL+2.758−0.054/L), A=(e²−5)/2, ile
  mean(M̃²)=1 sağlanıyor (§3 of Not 1). Bu, kaptanın M/δ² fikri için de doğru
  başlangıç noktası: CUE tarafında |Λ_N(θ)|² 'nin ortalaması M_N(1)=N
  (Keating–Snaith momentinden, k=1) olduğundan, N=L eşlemesiyle
  mean(|Λ|²)~L uyumlu — Not 1'in mean(M̃²)=1 normalizasyonu bu yüzden RMT
  tarafıyla tutarlı bir taban seçimi.
- **N ↔ L eşlemesi:** standart taban N=L (yoğunluk eşleştirme, KS00). Not 1
  BUNUN YETERSİZ olduğunu, gözlemlenebilir-bağımlı bir N_eff=L+c gerektiğini
  BULDU (Pearson: +0.75→+1.30 sürükleniyor; Spearman: sabit +1.99-2.14).
  **Küçük-δ tilt fikrini sınarken bare N=L değil, Not 1'in ölçtüğü N_eff
  kullanılmalı** — aksi halde tilt modelinin başarısızlığı N_eff
  yanlışlığından mı yoksa tilt fikrinin kendisinden mi geldiği
  karıştırılabilir (confound riski, aşağıda §4'te tekrar).

---

## 3. Veri envanteri (qm_riemann/)

### 3.1 Zaten hesaplanmış (δ_n, M_n) çiftleri — dosya dosya

Tüm M_n değerleri **iki motordan biriyle** hesaplandı (Not 1, §2):

- **Motor A** (mpmath, ince-hassasiyetli kök-bulma tipi): her aralıkta
  `scipy.optimize.minimize_scalar(..., method='bounded')` ile
  `-|mpmath.siegelz(t)|` minimize edilir (xatol 0.02–0.05, dps=10–12).
  Kaynak: `33_T_buyuk.py`, `36_T100k.py`.
- **Motor B** (vektörize Riemann–Siegel, float64): her aralıkta 24 noktalı
  iç ızgara + üç-noktalı parabolik rafinasyon ile tepe konumu ve değeri
  bulunur (bkz. `53_odlyzko_motor.py` satır ~110-127, `41_bigT_scan.py`,
  `55_kucuk_tau_scan.py`). Ek önlem: şüpheli (iç min|Z|<0.1) aralıklarda
  10× ince yeniden-tarama (`41_bigT_scan.py` docstring) — bu, tam da küçük
  aralıkların kaçırılmaması için konmuş bir güvenlik önlemi.

| Dosya | Anahtarlar (shape) | t-aralığı | L-aralığı | n (aralık) | Motor |
|---|---|---|---|---|---|
| `33_max_amps_T30000.npz` | intervals,max_amps,t_mid (35672,) | [17.6, 3.0e4] | [1.03,8.47] | 35 672 | A |
| `36_T100k.npz` | intervals,max_amps,t_mid (99999,) | [17.6, 7.49e4] | [1.03,9.39] | 99 999 | A |
| `41_bigT_windows.npz` | gaps_/amps_/tmid_{120k,200k,350k,600k,1000k,1600k} (~40000 her biri) | 1.07e5–1.61e6 | 9.75–12.45 | 6×~40 000 | B |
| `53_odlyzko_amps.npz` | gaps,max_amps,t_mid (9999,), L (skaler) | t≈2.6765×10¹¹ (T0=267653395647 ofsetli) | 24.475 (sabit) | 9 999 | B (çapa açılımlı) |
| `55_win_1e+08.npz` | gaps,amps,tmid (29999,) | ~1e8 | 16.583 | 29 999 | B |
| `55_win_1e+09.npz` | gaps,amps,tmid (29999,) | ~1e9 | 18.885 | 29 999 | B |
| `55_win_1e+10.npz` | gaps,amps,tmid (24998,) | ~1e10 | 21.188 | 24 998 | B |
| `55_win_1e+11.npz` | gaps,amps,tmid (20000,) | ~1e11 | 23.491 | 20 000 | B |

**`84_yasa_ciftleri.npz` — DİKKAT, bu ham (δ,M) çifti DEĞİL.** Anahtarlar
`tr` (44,4), `oos` (22,4), `odl` (6,5) — yukarıdaki ham pencerelerden
(`41`, `55_win_*`, `53`) `83_buyuk_yeniden_olcum.py`'daki `kanal_tam()`
fonksiyonuyla türetilmiş (τ, w, s, v) tipi parametreler (farklı bir
"yasa"/kanal-uyum projesinin ürünü, bkz. `NOT7_*`/`LITERATUR_KAPPA_*`
dosyaları). **Küçük-δ sefer için doğrudan kullanılamaz**; ham çiftler
yukarıdaki kaynak pencerelerden yeniden türetilmeli.

### 3.2 Yalnız sıfır tabloları (M_n YOK — bir motordan geçirilmesi gerekir)

| Dosya | Anahtarlar | n (sıfır) | L-aralığı | Not |
|---|---|---|---|---|
| `128_odl_zeros1_1e5_zeros.npz` | zeros,taban,taban_str,n,L,L_min,L_max | 100 000 | [0.81,9.39] | `36`'nın kaynak tablosuyla aynı (δ-sayımları birebir örtüşüyor) |
| `128_odl_zeros6_2e6_zeros.npz` | aynı şema | 2 001 052 | [0.81,12.10] | **En büyük sürekli blok**, M_n hesaplanmamış — 41'in 6 penceresi bunun küçük bir alt-örneği |
| `128_odl_t1e12_1e4_zeros.npz` | aynı şema | 10 000 | 24.475 (sabit) | Odlyzko'nun ~10¹²'inci sıfır civarı tablosu (dosya adındaki "t1e12" YÜKSEKLİK değil, YAKLAŞIK SIFIR İNDEKSİ — taban≈2.68×10¹¹) |
| `128_odl_t1e21_1e4_zeros.npz` | aynı şema | 10 000 | 44.580 (sabit) | taban≈1.44×10²⁰ (yine indeks≈10²¹, yükseklik değil) |
| `128_odl_t1e22_1e4_zeros.npz` | aynı şema | 10 000 | 46.832 (sabit) | taban≈1.37×10²¹ |

`128_odlyzko_edinme.py` dört-kapılı bir doğrulama uyguluyor (satır sayısı,
uç değerler, kesin artanlık, sayım yoğunluğu N′(t)≈L/2π ±%5, ds varyansı
GUE'ye yakın) — veri kalitesi kontrolü mevcut ve kod içinde belgeli.

### 3.3 Küçük-aralık SAYIMLARI (yalnız δ̃, M'ye dokunulmadı)

Her dosya için δ̃_n = δ_n·L(t)/2π hesaplanıp eşik-altı sayıldı (L, ham sıfır
tablolarında noktasal; sabit-t pencerelerinde dosyanın kendi L değeri):

| Kaynak | n | δ̃<0.1 | δ̃<0.2 | δ̃<0.3 |
|---|---|---|---|---|
| 33 (T≤3e4) | 35 672 | 25 (0.070%) | 214 (0.600%) | 705 (1.976%) |
| 36 (T≤7.49e4) | 99 999 | 79 (0.079%) | 613 (0.613%) | 2109 (2.109%) |
| 41[120k] | 39 992 | 39 (0.098%) | 281 (0.703%) | 899 (2.248%) |
| 41[200k] | 39 997 | 38 (0.095%) | 286 (0.715%) | 900 (2.250%) |
| 41[350k] | 39 998 | 32 (0.080%) | 301 (0.753%) | 981 (2.453%) |
| 41[600k] | 39 998 | 42 (0.105%) | 307 (0.768%) | 1002 (2.505%) |
| 41[1000k] | 39 999 | 42 (0.105%) | 299 (0.748%) | 935 (2.338%) |
| 41[1600k] | 40 000 | 30 (0.075%) | 295 (0.738%) | 950 (2.375%) |
| 53 (L=24.475) | 9 999 | 6 (0.060%) | 83 (0.830%) | 262 (2.620%) |
| 55[1e8] | 29 999 | 32 (0.107%) | 247 (0.823%) | 775 (2.583%) |
| 55[1e9] | 29 999 | 32 (0.107%) | 258 (0.860%) | 790 (2.633%) |
| 55[1e10] | 24 998 | 31 (0.124%) | 202 (0.808%) | 674 (2.696%) |
| 55[1e11] | 20 000 | 13 (0.065%) | 177 (0.885%) | 551 (2.755%) |
| **128_zeros6 (L≤12.10, ham, M yok)** | 2 001 051 | **1824 (0.091%)** | 14 279 (0.714%) | 47 458 (2.372%) |
| 128_t1e12 (L=24.475, ham) | 9 999 | 6 (0.060%) | 83 (0.830%) | 262 (2.620%) |
| 128_t1e21 (L=44.580, ham) | 9 999 | 10 (0.100%) | 83 (0.830%) | 257 (2.570%) |
| 128_t1e22 (L=46.832, ham) | 9 999 | 7 (0.070%) | 79 (0.790%) | 266 (2.660%) |

(33⊂36 ve zeros1=36'nın kaynak tablosu, birbirinin tekrarı — çakışan
sayılmamalı. 128_t1e12 = 53 ile aynı kaynak tablo.)

**Toplam bağımsız, ŞU AN M_n'si hesaplanmış küçük-gap havuzu** (33/36
tekrarını 36 olarak sayarak, 128_t1e12'yi 53 olarak sayarak):
δ̃<0.1 için 36+41(6 pencere)+53+55(4 pencere) ≈ 79+265+6+108 ≈ **458 olay**;
δ̃<0.2 için ≈ **2778 olay**; δ̃<0.3 için ≈ **9636 olay**.
Buna karşılık `128_odl_zeros6_2e6_zeros.npz` TEK BAŞINA (M_n hesaplanırsa)
δ̃<0.1'de 1824 — yani mevcut havuzun ~4 katı — sağlayabilir; zero-bulma işi
zaten bitmiş, yalnız M_n motoru (Motor B, 24-nokta+parabolik) bu tablo
üzerinde çalıştırılmamış.

### 3.4 Hassasiyet (küçük δ'da M çok küçük olur — yeterli mi?)

- Motor B'nin nokta-bazlı mutlak hatası Not 1'e göre ≤7.6×10⁻⁶ (10⁵≤t≤1.6e6
  aralığında, mpmath'e karşı doğrulanmış). Bu, tipik M~O(1–10) için
  fazlasıyla yeterli, ama **M/δ² oranı küçük δ'da hatayı δ⁻² ile büyütür** —
  δ̃~0.05 (L~10) civarında gerçek δ~0.03 olur; M de orantılı küçülüp
  ~10⁻³–10⁻² mertebesine inebilir (Lehmer-tipi aşırı uçlarda çok daha
  küçük). 7.6×10⁻⁶ mutlak hata bu ölçekte hâlâ küçük ama **M/δ² gibi bir
  türetilmiş büyüklük için hata yayılımı bu oturumda ayrıca
  karakterize edilmedi** — sefer öncesi kontrol listesine eklenmeli (bkz.
  §4).
- Motor A (mpmath, dps=10–12) için ayrı bir sistematik hata karakterizasyonu
  script içinde yok (yalnız Motor B'ye karşı örtüşme-bölgesi karşılaştırması
  var, Not 1 §2: "%0.1 (medyan: exact)").
- **Z″ hesaplanabilirliği:** bu oturumda doğrudan test edildi —
  `mpmath.diff(mpmath.siegelz, t, 2)` (sayısal türev, dps=30) çalışıyor ve
  ~0.077 sn/çağrı sürüyor (tek çekirdek, `.venv` içinde ölçüldü). Mevcut hiçbir
  script (27-60, 84, 128 taraması) şu an Z′ veya Z″'yi hesaplamıyor —
  yalnız `mpmath.siegelz` (Z'nin kendisi) kullanılıyor. Yani Z″/ξ″ tabanlı
  bir yaklaşım (kaptanın M_n≈|Z″|δ_n²/8 önerisi) **sıfırdan
  uygulanmalı**, ama yöntem (mp.diff) doğrulandı ve makul hızda.

---

## 4. Ön-kayıtlı sefer iskeleti (öneri — karar kaptanın)

**Hipotezler (rakip modeller):**
1. **H_tilt (kaptanın fikri):** δ̃→0 limitinde, (M/δ²)|δ̃ küçük dağılımı,
   |Λ_{N−2}(θ)|⁴-eğilmiş Keating–Snaith yasasına yakınsar; momentler
   E[(M/δ²)^k | küçük δ̃] ∝ M_{N−2}(k+2)/M_{N−2}(2), N=N_eff(L)
   (Not 1'in ÖLÇTÜĞÜ N_eff, bare L değil — confound'u önlemek için).
2. **H_untilted (rakip 1):** eğilmemiş sıradan Keating–Snaith |Λ_N|
   dağılımı (tilt gereksiz/gözlenemez).
3. **H_gauss (rakip 2, Not 1'in Null 1'inin küçük-δ dilimine özelleşmiş
   hâli):** faz-rastgeleştirilmiş Gauss yüzey null'u, yalnız küçük-δ
   diliminde.
4. **H_arith (rakip 3):** H_tilt + Keating–Snaith aritmetik çarpanı a(k)
   (kaydırılmış moment) — zeta'ya özgü düzeltme var mı yok mu.
5. **H_CUE_conditioned_lit:** Charlier–Claeys tipi FH/Toeplitz koşullamanın
   (eğer birincil metinden doğrulanırsa) verdiği kapalı form, varsa.

**Test edilecek istatistikler:** δ̃ için küçülen bantlarda (ör.
[0,0.05),[0.05,0.10),...) koşullu M/δ² momentleri (k=1,2,3), model
tahminleriyle moment-oranı testi; CUE tarafında aynı analiz Monte Carlo ile
(ucuz, kesin) birincil doğrulama zemini olarak kullanılmalı — zeta
tarafındaki örneklem küçük (§3.3).

**Kaba güç tahmini:** mevcut havuzda δ̃<0.1 için ~458, δ̃<0.2 için ~2778 olay
var (§3.3); `128_odl_zeros6` üzerinde M_n hesaplanırsa δ̃<0.1 tek başına
~1824'e çıkar. k=1,2 moment tahmini için birkaç yüz-bin olay (δ̃<0.3
bandında ~9600 mevcut, zeros6 eklenince ~57 000) muhtemelen yeterli;
k=3+ momentler veya asimptotik δ̃→0 ekstrapolasyonu için CUE Monte Carlo'nun
ana kanıt kaynağı olması, zeta verisinin doğrulama/kalibrasyon rolü
üstlenmesi mantıklı görünüyor.

**Sefer öncesi kontrol listesi (önerilir):**
- (a) M/δ² için hata yayılımı karakterizasyonu (Motor B'nin 7.6e-6 mutlak
  hatasının küçük-δ bandında M/δ²'ye etkisi).
- (b) N_eff(L) confound'unu izole etmek için hem bare N=L hem Not 1'in
  ölçtüğü N_eff ile ayrı ayrı fit.
- (c) `128_odl_zeros6_2e6_zeros.npz` üzerinde Motor B'nin (24-nokta+
  parabolik, gerekirse şüpheli-aralık ince-tarama) çalıştırılması —
  mevcut en büyük, kullanılmayan ham havuz.
- (d) Z″ yolunu (mp.diff, doğrulandı) küçük bir alt-örnekte Motor B'nin
  M_n'siyle çapraz doğrulamak (M_n ≈ |Z″|δ_n²/8 yaklaşıklığının kendisini,
  düşük δ̃'de, KÖRLÜK KURALINI ihlal etmeden — yalnız yöntem/kod
  doğrulaması, tam istatistik değil).

---

## Kaynakça (bu belgede DOĞRULANDI işaretli olanlar)

- Csordas, Smith, Varga (1994), Constructive Approximation 10, 107–129. DOI 10.1007/BF01205170.
- Stopple (2015), arXiv:1508.05870, Experimental Mathematics.
- Rodgers, Tao (2020), Forum of Math. Pi 8, e6. DOI 10.1017/fmp.2020.6.
- Farmer, Gonek, Lee (2014), JLMS (2) 90, 241–269. DOI 10.1112/jlms/jdu026.
- Assiotis, Keating, Warren (2022), IMRN 2022(18), 14564–14603. DOI 10.1093/imrn/rnab336. arXiv:2005.13961.
- Dehaye (2008), Algebra & Number Theory 2(1), 31–68. DOI 10.2140/ant.2008.2.31.
- Charlier, Claeys (2016), arXiv:1604.08399.
- Conrey, Ghosh (1985), JLMS (2) 32, 193–202. DOI 10.1112/jlms/s2-32.2.193.
- Hughes, Lugmayer, Pearce-Crump (2024/2025), arXiv:2411.05573, JLMS (2) 112, e70250.
- Pearce-Crump (2024/2025), arXiv:2411.05568, Mathematika 71, e70035.
- Keating, Snaith (2000), Comm. Math. Phys. 214, 57–89. DOI 10.1007/s002200000261.
- Durkan, Hughes, Pearce-Crump (2026), arXiv:2601.18025.
