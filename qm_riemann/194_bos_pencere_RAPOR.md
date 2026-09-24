# 194 — BOŞ PENCERE TABANI: 193'ün germesi ve 192'nin μ=0 çukurları tek bir eksi taban mı?

**Soru (KALEM_BOS_PENCERE_TABANI_24EYL2026, commit 9e74e4a — ön-kayıtlı):** 193,
kalıntı-sınıfı payı desenini (cos yasası) doğru yönde ama ortak bir çarpanla
GERİLMİŞ buldu (+log7/+log5'te ölçülen ≈1.3–2.3× öngörü). 192, μ(a)=0
konumlarında anlamlı negatif çukurlar gördü. **Hipotez:** her δ=0.03 blok-yerel
Δω penceresinde, uydudan ve kalıntı sınıfından BAĞIMSIZ, sabit bir eksi taban
κ_boş ≈ −0.0025 var; germe ve çukurlar bunun iki görünümü. **Rakip:** komşu
uydu sızıntısı — bu doğruysa uydulardan UZAK ("boş") pencerelerde taban ≈ 0
çıkar.

Tek dalga koşuldu. Ölümler kurtarılmadı. Ön-kayıtlı kapılar (makine mührü
dahil) atlanmadı. Git'e dokunulmadı. Sonuç ORTAK TEFTİŞE sunulur, commit
sonra.

---

## Ölçüm tanımı (193b_sinif.py AYNEN; yalnız pencere merkezleri değişti)

188'in seri fonksiyonu, blok-yerel Δω seçimi (δ=0.03, her blok kendi yerel
L_b'siyle), sınıf ayrıştırması, HAVUZ κ = −Re K_HAVUZ tanımı ve 8-blok
jackknife (loo) — 193b_sinif.py'den `importlib` ile AYNEN yüklendi (dosya
düzenlenmedi), yalnız hedef merkezleri (ve boş pencerelerde sınıf modülü: a=1
birincil "sınıfsız" κ, ayrıca mod 3/5/10) değişti [194b_olcum.py].

---

## K0 — DONMUŞ ÖN-KAYIT (yalnız aritmetik, veri okunmadı)  [194a_onkayit.py]

`scratchpad/194/ONKAYIT_194.json`:
- **sha256 = f0ab5e767771fd8b7ee46719b8e0ffffd3502983d67d6f2b8bc12eec0b01cbbe**
- **damga Thu Sep 24 13:46:26 +03 2026**

**BOŞ PENCERE arama:** a,b≤12 sade kesirlerin log(a/b) konumları (hem μ≠0 hem
μ=0) taranarak 6'nın altında çıktı (45 katalog noktası → yalnız **2** boş
pencere) → kalemin talimatıyla sınır a,b≤10'a indirildi (31 katalog noktası →
**4** boş pencere). **a,b≤10'da da 4 < 6** — kalemde bundan öte bir gevşetme
adımı YOK; eşik gevşetilmedi, mevcut 4 pencereyle devam edildi (bu sapma
burada kaydedilir) [194a_onkayit.py].

BOŞ PENCERELER (N=4; merkez = boşluğun ortası, her iki komşudan eşit uzak):

| ad | Δω merkez | sol komşu (kataloğ) | sağ komşu (kataloğ) | uzaklık (her iki yöne) |
|---|---|---|---|---|
| bos1 | 1.31953 | 7/2 = 1.25276 (μ=−1) | 4/1 = 1.38629 (μ=0) | 0.06677 |
| bos2 | 1.70060 | 5/1 = 1.60944 (μ=−1) | 6/1 = 1.79176 (μ=1) | 0.09116 |
| bos3 | 1.86883 | 6/1 = 1.79176 (μ=1) | 7/1 = 1.94591 (μ=−1) | 0.07708 |
| bos4 | 2.01268 | 7/1 = 1.94591 (μ=−1) | 8/1 = 2.07944 (μ=0) | 0.06677 |

Birbirinden ayrım: [0.3811, 0.1682, 0.1438] ≥ 0.06 — tamamı sağlandı
[194a_onkayit.py].

MU=0 PENCERELERİ (kalemdeki liste AYNEN; 8'i Δω∈[0.30,2.30] içinde, 1'i
(+log(4/3)=0.28768) aralık dışında → KAYIT, ölçülmedi):

| ad | a/b | Δω merkez | durum |
|---|---|---|---|
| +log(8/5) | 8/5 | 0.47000 | aralık içi |
| +log(9/5) | 9/5 | 0.58779 | aralık içi |
| +log(9/4) | 9/4 | 0.81093 | aralık içi |
| +log(8/3) | 8/3 | 0.98083 | aralık içi |
| +log4 | 4/1 | 1.38629 | aralık içi |
| +log(9/2) | 9/2 | 1.50408 | aralık içi |
| +log8 | 8/1 | 2.07944 | aralık içi |
| +log9 | 9/1 | 2.19722 | aralık içi |
| +log(4/3) | 4/3 | 0.28768 | **KAYIT — Δω<0.30, aralık dışı** |

H-194a..e eşikleri kalemden AYNEN donduruldu [194a_onkayit.py].

---

## MAKİNE MÜHÜRÜ  [194b_olcum.py]

+log10 hedefi, 193'ün KENDİ hedef tanımıyla (a193.HEDEFLER, aynı a=10,
sınıflar {1,3,7,9}, aynı log_n) bu betiğin kendi çağrı zinciriyle (193b'den
AYNEN alınan `blok_yerel_D`, `K_hav`, 8-blok jackknife) yeniden ölçüldü:

| mühür | 194 (yeniden ölçüm) | 193 (arşiv, K1_sinif.json) | \|Δ\| |
|---|---|---|---|
| +log10 κ_top | 0.0140115331 | 0.0140115331 | **0.0** (bit-bit) |

Bit-bit eşleşme — 194'ün kinematik/pencere/asal-kuvvet/jackknife kurulumu
193'le TIPATIP aynı [194b_olcum.py].

---

## K1 — ÖLÇÜMLER (blok-yerel Δω, HAVUZ κ)  [194b_olcum.py]

**BOŞ PENCERELER** (a=1, sınıfsız — tüm seçilen çizgiler tek "havuz" serisi):

| pencere | Δω | κ (HAVUZ) ± jk se | çizgi sayısı |
|---|---|---|---|
| bos1 | 1.31953 | **+0.00074 ± 0.00028** | 22 572 |
| bos2 | 1.70060 | **−0.00116 ± 0.00028** | 32 073 |
| bos3 | 1.86883 | **−0.00087 ± 0.00042** | 37 647 |
| bos4 | 2.01268 | **+0.00225 ± 0.00040** | 42 849 |

**κ̄_boş (havuzlanmış — pencereler arası, 8-blok jackknife replikaları önce
pencereler üzerinden ortalanıp sonra jk() uygulanarak) = +0.00024 ± 0.00017**
(+1.4σ; N=4 pencere) [194c_hukum.py].

**MU=0 PENCERELERİ** (a=1, aralık-içi 8 tanesi):

| pencere | Δω | κ (HAVUZ) ± jk se | çizgi sayısı |
|---|---|---|---|
| +log(8/5) | 0.47000 | −0.00193 ± 0.00028 | 10 294 |
| +log(9/5) | 0.58779 | −0.00104 ± 0.00027 | 11 521 |
| +log(9/4) | 0.81093 | −0.00149 ± 0.00023 | 14 067 |
| +log(8/3) | 0.98083 | −0.00147 ± 0.00031 | 16 565 |
| +log4 | 1.38629 | −0.00290 ± 0.00038 | 24 031 |
| +log(9/2) | 1.50408 | −0.00100 ± 0.00031 | 26 844 |
| +log8 | 2.07944 | −0.00155 ± 0.00023 | 45 643 |
| +log9 | 2.19722 | −0.00206 ± 0.00029 | 50 989 |

**κ̄_μ0 (havuzlanmış) = −0.00168 ± 0.00012** (N=8 pencere) — sekizinin
SEKİZİ de negatif; boş pencerelerin dörttte ikisi pozitif [194b/194c].

**Boş pencerelerde mod 3/5/10 sınıf payları (havuzlanmış s_r, H-194c
için):**

| modül | sınıf | s_r (havuz) ± se | öngörü 1/φ(a) | ±[0.10/φ(a)+2σ] bant |
|---|---|---|---|---|
| mod3 (φ=2) | r1 | 0.3029 ± 0.0568 | 0.5000 | **HAYIR** |
| | r2 | 0.6971 ± 0.0568 | 0.5000 | **HAYIR** |
| mod5 (φ=4) | r1 | 0.2576 ± 0.1594 | 0.2500 | EVET |
| | r2 | 0.1950 ± 0.1904 | 0.2500 | EVET |
| | r3 | 0.1475 ± 0.1020 | 0.2500 | EVET |
| | r4 | 0.3999 ± 0.1430 | 0.2500 | EVET |
| mod10 (φ=4) | r1 | 0.2576 ± 0.1594 | 0.2500 | EVET |
| | r3 | 0.1475 ± 0.1020 | 0.2500 | EVET |
| | r7 | 0.1950 ± 0.1904 | 0.2500 | EVET |
| | r9 | 0.4000 ± 0.1430 | 0.2500 | EVET |

[194c_hukum.py]

---

## HÜKÜM (eşikler K0'da donmuş; kurtarma yok)  [194c_hukum.py]

| hipotez | hüküm | dayanak |
|---|---|---|
| **H-194a** (EKSİ TABAN VAR) | **ÖLDÜ** | κ̄_boş=+0.00024±0.00017 (+1.4σ, POZİTİF); bant [−0.0035,−0.0015] dışında; ölüm koşulu (\|κ̄_boş\|<2σ VEYA pozitif) sağlandı — κ̄_boş zaten pozitif |
| **H-194b** (SABİTLİK, KAYIT) | **SABİT** | eğim=+0.00092±0.00055 (+1.7σ, 2σ içinde) — Δω ile anlamlı eğilim yok |
| **H-194c** (SINIFTAN BAĞIMSIZ, eşiksiz) | **KISMEN** | mod5/mod10 bantta (geniş se ile); mod3'te r1=0.30/r2=0.70 bant dışı (~3.5σ) — boş pencerelerde de hafif mod-3 asimetrisi var, SINANMADI (bkz. aşağı) |
| **H-194d** (ÇUKUR=TABAN) | **ÖLDÜ** | κ̄_μ0−κ̄_boş = −0.00192, σ_comb=0.00020 → **−9.4σ** — ölüm koşulu (μ0 ≥3σ daha derin) fazlasıyla sağlandı; çukurlar taban DEĞİL, ayrı bir olgu |
| **H-194e** (GERME KAPANIŞI) | **KAYIT (kısmi)** | κ̄_boş≈0 olduğu için düzeltme ihmal edilebilir düzeyde (+log10 payları ~0.01 kayıyor); +log10'un 4 sınıfından yalnız r7 ±0.10 bandında, ölüm eşiği ±0.15 hiçbirinde aşılmadı; +log7/+log5 hemen hemen değişmeden ±0.10 dışında kalıyor |

**TAM-ÖRNEKLEM/MAKİNE MÜHÜRÜ:** TUTTU (|Δ|=0.0, bit-bit) [194b_olcum.py].

---

## Dürüst okuma

- **Ana hipotez (H-194a) temiz öldü.** Boş pencerelerin ortalaması
  (+0.00024±0.00017) kalemin öngördüğü −0.0025'in tam tersi yönde: hafifçe
  POZİTİF, sıfırdan sadece 1.4σ uzakta. Dört pencerenin ikisi pozitif
  (bos1=+0.00074, bos4=+0.00225), ikisi negatif (bos2=−0.00116,
  bos3=−0.00087) — sistematik bir eksi taban değil, penceden pencereye
  işaret değiştiren bir gürültü deseni. Bu tam olarak **rakip hipotezin
  (komşu uydu sızıntısı) öngördüğü şey**: uydulardan uzak pencerelerde taban
  ≈ 0.
- **H-194d, çukurların tabanla açıklanamayacağını kesin gösteriyor.**
  κ̄_μ0 = −0.00168±0.00012, κ̄_boş'tan **−9.4σ** uzakta — μ=0 pencerelerin
  SEKİZİNİN SEKİZİ de negatif (κ̄_boş'un aksine), çok daha güçlü ve tutarlı
  bir etki. 192'nin μ=0 çukurları, 194'ün ölçtüğü "boş pencere" gürültüsüyle
  AYNI OLGU DEĞİL — kendi başına, daha derin, daha sistematik bir şey.
- **H-194e mekanik olarak neredeyse hiçbir şey değiştirmiyor**, çünkü
  düzeltme terimi κ̄_boş/φ(a) ≈ 0.0001 mertebesinde — 193'ün ham paylarının
  yanında ihmal edilebilir. Bu yüzden 193'ün asıl bulmacası (+log7/+log5'in
  neden ~1.3–2.3× büyük çıktığı) 194 ile ÇÖZÜLMEDİ; germe olgusu hâlâ
  açıklanmamış duruyor.
- **Beklenmedik yan-gözlem (SINANMADI):** boş pencerelerde bile mod-3 payı
  simetrik değil (r1=0.30/r2=0.70, ~3.5σ). δ=0.03 blok-yerel seçiminin
  kendisi mi (örn. q' aralığının belirli mod-3 sınıflarını sistematik
  favoriye alan bir kenar etkisi), yoksa gerçek bir ikincil yapı mı —
  bilinmiyor; H-194c eşiksiz olduğu için hükme girmedi, ama not düşülür.
- **N=4 boş pencere küçük bir örneklem** — kalemin kendi ≥6 eşiği bile
  a,b≤10'da karşılanamadı (a,b≤12'de yalnız 2, a,b≤10'da 4 çıktı; ötesine
  geçilmedi). κ̄_boş'un işareti/büyüklüğü konusundaki güven bu küçük N ile
  sınırlı — ama H-194a'nın ölüm koşulu (pozitif YA DA <2σ) zaten net
  sağlandığı için ek pencere sayısı hükmü değiştirmez (pozitif olmak tek
  başına yeterli ölüm nedeni).

---

## TÜRETİLEN vs ÖLÇÜLEN

**Türetilen (veri-öncesi, kalem):**
- Boş pencere listesi: Δω∈[0.30,2.30], a,b≤12(→10) sade kesir kataloğundan
  ≥0.06 uzak boşluk-ortaları — saf aritmetik, veriye bakılmadan [194a].
- μ=0 pencere listesi kalemden AYNEN.
- H-194a..e eşikleri kalemden AYNEN.
- H-194e düzeltme formülü s_r' = (κ_r−κ̄_boş/φ(a))/(κ_top−κ̄_boş) — kalemin
  "veri-öncesi düzeltme hesabı" örneğiyle aynı formül.

**Ölçülen:**
- 4 boş + 8 μ=0 pencerenin κ'sı (HAVUZ, 8-blok jk se); +log10 makine
  mührü (bit-bit); boş pencerelerde mod 3/5/10 sınıf payları.
- κ̄_boş, κ̄_μ0 (havuzlanmış, pencereler-arası pooled jackknife); Δω eğilimi;
  193'ün düzeltilmiş payları (delta-yöntemi se yayılımıyla, κ_r/κ_top 193'ten
  sabit alınıp κ̄_boş ile bağımsız varsayılarak kombine edildi — kovaryans
  hesaba katılmadı, kalem bunu istemiyordu).

**Türetilmemiş/açık kalanlar:**
- 193'ün germe olgusunun (+log7/+log5 ~1.3–2.3× büyük) asıl nedeni — H-194e
  bunu ÇÖZMEDİ, çünkü κ̄_boş ≈ 0 çıktı.
- 192'nin μ=0 çukurlarının asıl nedeni — H-194d bunun "boş pencere tabanı"
  OLMADIĞINI gösterdi ama alternatif bir mekanizma sunmuyor.
- Boş pencerelerde görülen mod-3 asimetrisi (r1=0.30/r2=0.70) — SINANMADI,
  yeni bir açık kalem konusu olabilir.
- α=1.30 bağı; bu kalemde dokunulmadı.

---

## MANŞET (aday cümle)

> **Boş pencere tabanı hipotezi ÖLDÜ: uydulardan ≥0.06 uzak 4 pencerenin
> κ ortalaması +0.00024±0.00017 — kalemin öngördüğü −0.0025'lik eksi tabanın
> tam tersi yönde, sıfırla uyumlu** (H-194a). **μ=0 pencerelerinin ortalaması
> (−0.00168±0.00012) bu "boş taban"dan −9.4σ uzakta** (H-194d ÖLDÜ) —
> 192'nin çukurları ayrı bir olgu, 194'ün ölçtüğü pencere-gürültüsüyle
> açıklanamıyor. κ̄_boş≈0 olduğu için 193'ün ham payları üzerinde düzeltme
> etkisi ihmal edilebilir düzeyde kaldı; germe bulmacası (+log7/+log5'in
> ~1.3–2.3× büyük çıkması) ÇÖZÜLMEDEN açık kalıyor. Rakip hipotez (komşu uydu
> sızıntısı) bu turda kazandı: uydulardan uzak pencerelerde taban ≈ 0
> çıktığı gözlendi, tıpkı öngördüğü gibi.

---

Teslim:
- Bu rapor.
- `194_configs/`:
  - 194a_onkayit.py (K0 — yalnız aritmetik)
  - 194b_olcum.py (K1 — makine mührü + ölçüm)
  - 194c_hukum.py (K2 — hüküm)
  - 194d_figur.py (figür)
- `194_bos_pencere.png`: sol Δω ekseninde boş (gri) ve μ=0 (kırmızı)
  pencerelerin κ±se, κ̄_boş bandı (amber) ve κ̄_μ0 çizgisi (kırmızı kesikli);
  sağ 193 payları ham (gri) vs düzeltilmiş (elmas, yeşil=±0.10 bantta) vs
  cos yasası öngörüsü (siyah çizgi, ±0.10 gölge).
- `scratchpad/194/`:
  - ONKAYIT_194.json (boş+μ0 pencere listeleri, eşikler dahil)
  - log_194b.txt
  - muhur_194.json
  - proj_bos_*.npz (4 pencere × 4 modül), proj_mu0_*.npz (8 pencere)
  - K1_bos_194.json
  - HUKUM_194.json
