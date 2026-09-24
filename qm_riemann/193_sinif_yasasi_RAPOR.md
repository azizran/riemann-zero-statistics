# 193 — KALINTI-SINIFI YASASI: cos(2πrb/a)/μ(a) kör sınavı

**Soru (KALEM_SINIF_YASASI_24EYL2026):** 192, en keskin kör öngörüyü ("sınıflar
120° döner") öldürdü, ama ön-kayıtsız keşfinde (192f, (a)) veri-sonrası şunu
gördü: sınıflar dönmüyor, toplamın ekseninde **cos(2π r b/a)** işaretleriyle
diziliyor — +log3, −log3, +log6'da desen tuttu, +log5'te işaret deseni tuttu ama
büyüklük 1.5-2.6σ sapıyordu. Bu kalem o gözlemi türetip (kaptan hatası
düzeltilerek: φ0 ≡ 0, teftiş (c)) **KÖR** olarak sınıyor: en keskin, daha önce
hiç görülmemiş iki hedef **+log10** (birincil) ve **+log7** (ikincil), artı
doğrulayıcı/kontrol hedefleri +log5, +log3, −log3.

Her sayının yanında onu üreten betik adı vardır. Tek dalga koşuldu. Ölümler
kurtarılmadı. Ön-kayıtlı kapılar (makine mührü dahil) atlanmadı. Git'e
dokunulmadı. Sonuç ORTAK TEFTİŞE sunulur, commit sonra.

---

## Kalem cebiri (teftiş düzeltmeleriyle; kaynak KALEM_SINIF_YASASI_24EYL2026)

1. Durağan faz: Δω = log(a/b) uydusunda, durağan noktası blok içinde olan her
   çizgi q' için Ĝ_mid(ω') katkısı ∝ A_n·e^{i(2π q'b/a + φ0)}, A_n gerçel.
2. **φ0 türetildi:** 2πN̄(t) = tL − t + 7π/4 ⇒ durağan değer t*−7π/4; faz
   eğriliği φ''=−1/t<0 ⇒ durağan-faz çarpanı e^{−iπ/4}; orta noktalarda N(m_n)
   tam sayı ⇒ ek işaret yok. Toplam sabit −7π/4−π/4=−2π ⇒ **φ0 ≡ 0**
   (192'nin kaleminde −π/4 eksikti — kaptan hatası).
3. Çekirdek gerçel kısmı duyar (S_Re; corr −0.95/−0.85, EŞİTLİK DEĞİL, ~%90-95
   açıklayan vekil) ⇒ sınıf r (q' ≡ r mod a) katkısı ∝ cos(2π r b/a).
4. **PAY YASASI:** s_r := κ_r/κ_top = cos(2π r b/a) / μ(a) (payda Ramanujan
   toplamı Σ_{r coprime}cos(2πrb/a) = μ(a); çizgiler sınıflara eşit dağılır).

Öngörü tablosu (K0'da doğrulandı; ölçüme bakmadan aritmetikle yeniden üretildi,
KALEM'in yuvarlanmış sayılarıyla < 5e-4 fark) [193a]:

| uydu (Δω) | a,b | sınıf payları s_r (tam) | statü |
|---|---|---|---|
| **+log10** (2.3026) | 10,1 | r1:+0.8090 r9:+0.8090 r3:−0.3090 r7:−0.3090 | **KÖR, birincil** |
| **+log7** (1.9459) | 7,1 | r1:−0.6235 r6:−0.6235 r2:+0.2225 r5:+0.2225 r3:+0.9010 r4:+0.9010 | **KÖR, ikincil** |
| +log5 (1.6094) | 5,1 | r1:−0.3090 r4:−0.3090 r2:+0.8090 r3:+0.8090 | doğrulayıcı KAYIT |
| +log3 (1.0986) | 3,1 | r1:+0.5000 r2:+0.5000 | kontrol |
| −log3 (−1.0986) | 1 (a_teori) | keyfi mod-3: r1,r2 → ½,½ (sıfır-yapı) | kontrol |

---

## Ölçüm tanımı (K0'da donmuş; 188b/190b/192c makinesi AYNEN, yalnız SEÇİM yeni)

Son penceresi (L_son = 12.029593242). 192c'nin tüm pencerede ORTAK (blok-bağımsız)
τ'-penceresi yerine, **her blok kendi yerel L_b'siyle** Δω_b(q') = log q' − L_b
hesaplar; uydu n için |Δω_b(q') − log n| < **δ = 0.03** olan çizgiler O BLOĞUN
KENDİ noktalarıyla (190b'nin ω-dilim / tek-blok izdüşüm deseni, `izdusum_blok`)
sınıflandırılıp seri + izdüşüm hesaplanır; q' mod a sınıflarına ayrılır (a'yı
bölen çizgiler ayrı "bölünen" sınıfı, **payadan HARİÇ**, ayrı raporlanır).
8 bloğun katkıları sonra 188b.c_proj/K_matris'e AYNEN verilir (blok-loo
destekli) — yani blok-yerel (RE,IM) dizisi, biçimce 188b.izdusum'un ürettiği
(8, sınıf, pencere) diziyle özdeş; geri kalan makine (κ = −Re K_HAVUZ, 8-blok
jackknife, K(b,s) = Σ c^{(s)}·conj(karışım)/Σ|karışım|², karışım = c^kesik −
c^öz, 186/187 defterleri) AYNEN 192c/188b [193b].

8 eşit blok (kenar = linspace(0,N,9), N=299999), L_b = blok içi ortalama
log(m_n/2π) [190a/192a AYNEN]:

| blok | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|---|---|
| L_b | 11.9632 | 11.9829 | 12.0023 | 12.0212 | 12.0398 | 12.0580 | 12.0759 | 12.0934 |

Alt-örneklem her 3. nokta (100 000; 188/192 ile aynı); +log10 r=1 sınıfında TAM
örneklem kontrolü (N=299 999).

---

## K0 — DONMUŞ ÖN-KAYIT  [193a_onkayit.py]

`scratchpad/193/ONKAYIT_193.json`:
- **sha256 = ff0852669f9f579d…** (tam: ff0852669f9f579dc79dc55035631b12011e9e9343b047dcd7cc5c1538fef67b)
- **damga Thu Sep 24 11:29:56 +03 2026**

Donanlar: öngörü tablosu (yukarıdaki, pay yasasıyla kalemin yuvarlanmış
sayılarına < 5e-4 uyumu doğrulanarak), δ=0.03, blok tanımı + L_b, sınıf
tanımları (a, b, sınıflar, "bölünen payadan hariç"), H-193a/b eşikleri ve ölüm
koşulları, KAYIT tanımları, 11 girdi dosyasının sha256'sı, **KOMŞU RAPORU**
(yalnız aritmetik, profillere/verilere BAKILMADAN).

**Evren kontrolü** [193a]: her hedef için gereken τ' üst sınırı (8 blok, en
kötü durum) τ'≤1.30 evreninin İÇİNDE: +log10 maks 1.1992, +log7 maks 1.1696,
+log5 maks 1.1416, +log3 maks 1.0991, −log3 maks 0.9165 — beşi de evren içi.

**KOMŞU RAPORU** (a,b ≤ 20 sade rasyoneller, |log(a/b) − hedef| < 0.10;
2δ=0.06 içi = "kirlenme riskli") [193a]:

| hedef | komşu (a/b) | Δfark | riskli mi |
|---|---|---|---|
| +log10 | 19/2 | −0.0513 | **EVET** |
| | 11/1 | +0.0953 | hayır |
| +log7 | 20/3 | −0.0488 | **EVET** |
| | 15/2 | +0.0690 | hayır |
| | 13/2 | −0.0741 | hayır |
| +log5 | 19/4 | −0.0513 | **EVET** |
| | 16/3, 14/3, 11/2 | +0.065/−0.069/+0.095 | hayır |
| +log3 | 20/7, 19/6, 17/6 | −0.049/+0.054/−0.057 | **EVET (üçü de)** |
| | 16/5, 14/5, 13/4, 11/4 | +0.065…−0.087 | hayır |
| −log3 | 7/20, 6/19, 6/17 | +0.049/−0.054/+0.057 | **EVET (üçü de)** |
| | 5/16, 5/14, 4/13, 4/11 | −0.065…+0.087 | hayır |

Beş hedefin BEŞİ de en az bir "kirlenme riskli" komşu taşıyor (δ=0.03'ün 2×'i,
0.06, ızgarada seyrek değil). Hüküm yine de ön-kayıtlı literal kuralla verildi
(etiketle); +log10'un tek riskli komşusu 19/2 (μ=−1) yalnızca −0.051 uzakta.

---

## MAKİNE MÜHÜRÜ  [193b_sinif.py]

188b'nin seri fonksiyonunun 193'ün kendi çağrı zincirinde 188 ile bit-bit aynı
olduğu, 192c ile TIPATIP aynı sınavla doğrulandı:

| mühür | sonuç |
|---|---|
| 188 `dilim_44.npz` (τ'∈(1.080,1.085], 2070 çizgi) aynı liste + seri_ve_G ile | dds, Gre, Gim, q, a, τ **BİT-BİT** |
| tek-dilim K (50 bant) vs 188 `harita_K_gercek.npz` K[:,44] | maks\|Δ\| = 1.02e-17 (tam), 1.30e-17 (loo) |
| doğrusallık: dilim-44 çizgilerinin mod-3 sınıf serileri toplamı = dilim serisi | maks\|Δ\| = 9.71e-17 (maks\|seri\| 0.099) |

**Tam-örneklem kontrolü (+log10 r1, blok-yerel, N=299999)** [193b]: κ_alt =
0.012596±0.000303, κ_tam = 0.012344±0.000171, |Δ| = 2.53e-04 = **0.83 se_alt**
— alt-örneklem yansız.

---

## K1 — SINIF ÖLÇÜMLERİ (blok-yerel Δω)  [193b_sinif.py]

s_r = κ_r/κ_top; κ_top = Σ_{r∈sınıflar} κ_r (**bölünen hariç**); z = s_r/se_s.
Bölünen sınıf her beş hedefte de neredeyse tamamen boş (0-3 çizgi, κ≈0) —
prime-KUVVET (2^k, 3^k, 5^k, 7^k …) çizgileri bu ölçekte (q'~10⁶) çok seyrek,
narrow δ-penceresine rastlaması olağan değil.

**+log10 (birincil, KÖR)** — τ'∈(1.1834, 1.1992], κ_top = 0.01401±0.00049:

| sınıf | çizgi | κ_r | s_r (ölçülen) | öngörü | z |
|---|---|---|---|---|---|
| r1 | 14 098 | +0.01260±0.00030 | **+0.8990±0.0231** | +0.8090 | **+38.9** |
| r9 | 13 927 | +0.01264±0.00034 | **+0.9024±0.0172** | +0.8090 | **+52.4** |
| r3 | 14 072 | −0.00582±0.00045 | **−0.4156±0.0345** | −0.3090 | **−12.0** |
| r7 | 14 014 | −0.00541±0.00045 | **−0.3858±0.0353** | −0.3090 | **−10.9** |
| bölünen | 0 | 0 | — (payda dışı) | — | — |

**+log7 (ikincil, KÖR)** — τ'∈(1.1537, 1.1696], κ_top = 0.00249±0.00044
(4σ eşik = 0.00175; erişilebilir):

| sınıf | çizgi | κ_r | s_r (ölçülen) | öngörü | z |
|---|---|---|---|---|---|
| r1 | 6 795 | −0.00359±0.00020 | −1.4391±0.2938 | −0.6235 | −4.90 |
| r6 | 6 798 | −0.00379±0.00016 | −1.5190±0.2955 | −0.6235 | −5.14 |
| r3 | 6 750 | +0.00425±0.00016 | +1.7054±0.2898 | +0.9010 | +5.89 |
| r4 | 6 726 | +0.00413±0.00012 | +1.6571±0.2695 | +0.9010 | +6.15 |
| r2 (yan) | 6 738 | +0.00083±0.00028 | +0.3319±0.0961 | +0.2225 | +3.45 |
| r5 (yan) | 6 671 | +0.00066±0.00023 | +0.2637±0.0953 | +0.2225 | +2.77 |
| bölünen | 0 | 0 | — | — | — |

**+log5 (doğrulayıcı, KAYIT)** — τ'∈(1.1258, 1.1416], κ_top = 0.00602±0.00036:

| sınıf | çizgi | κ_r | s_r | öngörü | z |
|---|---|---|---|---|---|
| r1 | 7 475 | −0.00319±0.00032 | −0.5307±0.0715 | −0.3090 | −7.42 |
| r4 | 7 375 | −0.00327±0.00032 | −0.5441±0.0558 | −0.3090 | −9.76 |
| r2 | 7 333 | +0.00641±0.00036 | +1.0661±0.0547 | +0.8090 | +19.5 |
| r3 | 7 409 | +0.00607±0.00027 | +1.0087±0.0603 | +0.8090 | +16.7 |
| bölünen | 0 | 0 | — | — | — |

**+log3 (kontrol, ½/½±0.10)** — τ'∈(1.0833, 1.0991], κ_top = 0.01926±0.00054:

| sınıf | çizgi | κ_r | s_r | bantta (½±0.1) |
|---|---|---|---|---|
| r1 | 9 238 | +0.00889±0.00036 | +0.4618±0.0181 | **EVET** |
| r2 | 9 217 | +0.01037±0.00051 | +0.5382±0.0181 | **EVET** |
| bölünen | 3 | ≈0 | — | — |

**−log3 (kontrol, sıfır-yapı, ½/½±0.10)** — τ'∈(0.9007, 0.9165],
κ_top = 0.00940±0.00034:

| sınıf | çizgi | κ_r | s_r | bantta (½±0.1) |
|---|---|---|---|---|
| r1 | 1 250 | +0.00469±0.00025 | +0.4989±0.0213 | **EVET** |
| r2 | 1 248 | +0.00471±0.00028 | +0.5011±0.0213 | **EVET** |
| bölünen | 3 | ≈0 | — | — |

---

## HÜKÜM (eşikler K0'da donmuş; kurtarma yok)  [193c_hukum.py]

| hipotez/kapı | hüküm | dayanak |
|---|---|---|
| **K0** ön-kayıt | **GEÇTİ** | sha ff0852669f9f…, 11:29:56; ölçüm (193b) 11:34'te başladı [193a] |
| makine mührü | **TUTTU** | dilim_44 bit-bit; K 1.02e-17; doğrusallık 9.71e-17 [193b] |
| tam-örneklem kontrolü | **TUTTU** (KAYIT) | 0.83 se_alt [193b] |
| **H-193a** +log10 (birincil, KÖR) | **MÜHÜR** | r1,r9 pozitif ≥2σ (z=+38.9,+52.4); r3,r7 negatif ≥2σ (z=−12.0,−10.9); dördü de öngörünün ±0.25 içinde (fark 0.077-0.107); ölüm koşulu (yanlış yön ≥2σ) sağlanmadı [193b,193c] |
| **H-193b** +log7 (ikincil, KÖR) | **KAYIT (kısmi)** | yön dördünde de doğru ve ≥2σ (r1,r6 negatif z=−4.9/−5.1; r3,r4 pozitif z=+5.9/+6.2) — ölüm koşulu sağlanmadı; ANCAK nicel ±0.30 bandı dördünde de AŞILDI (ölçülen ≈1.5-1.9× öngörü büyüklüğünde); κ_top = 5.7σ ≥ 4σ eşiği ⇒ erişilemedi DEĞİL [193b,193c] |
| +log5 (doğrulayıcı KAYIT) | işaret deseni 4/4 doğru, büyüklük ~1.3-1.8× öngörü | eşiksiz kayıt [193b,193c] |
| +log3 (kontrol) | **bantta** (½±0.10) | r1=0.462, r2=0.538, ikisi de \|s−0.5\|≤0.10 [193b,193c] |
| −log3 (kontrol, sıfır-yapı) | **bantta** (½±0.10) | r1=0.499, r2=0.501 — 0.5'e 0.001 farkla; sıfır-yapı (a_teori=1) temiz doğrulandı [193b,193c] |

**Dürüst okuma:**
- **+log10 kör sınavı temiz ve keskin geçti.** Dört sınıfın hepsi doğru yönde,
  hepsi öngörünün ±0.11 içinde, z-değerleri 11 ile 52 arasında — 192'nin
  ön-kayıtsız keşfinin (192f, (a)) KÖR doğrulaması budur. 192'de "dönme yok,
  cos deseni" yalnız +log3/−log3/+log6'da (a=2,3 gibi küçük/simetrik
  modüllerde) görülmüştü; +log10 (a=10, dört farklı sınıf, hem +hem− işaret)
  ilk kez görülen, en ayırt edici desendi ve tuttu.
- **+log7 yön olarak da tutuyor ama büyüklük sistematik olarak büyük.**
  Σ_r s_r ≡ 1 özdeşliği (κ_top = Σκ_r tanımından) hem veride hem öngörüde
  sağlanıyor (kontrol: −1.439+0.332+1.705+1.657+0.264−1.519 = 1.000), yani bu
  bir global genlik kayması değil, sınıflar arası bir yeniden-dağılım: ölçülen
  desen öngörülenden daha "keskin" (uç sınıflar daha ağır, orta sınıflar
  görece az farklı). Aynı yön-doğru/büyüklük-fazla deseni +log5'te de var.
  Hem +log7 hem +log5'in komşu raporunda "kirlenme riskli" bir rakip var
  (20/3, 19/4) — ama +log10'un da var (19/2) ve o kusursuz tuttu, o yüzden
  kirlenme tek başına açıklayıcı değil. Kalemin S_Re yaklaşıklığı zaten
  "%90-95 açıklayan vekil, eşitlik değil" olarak işaretlenmişti (madde 3) —
  büyüklük fazlalığı bu yaklaşıklığın kalıntısı olabilir. **SINANMADI.**
- **Kontroller (+log3, −log3) mükemmel tuttu**, özellikle −log3'ün 0.499/0.501
  sonucu: gerçek yapı yokluğunun (a_teori=1) keyfi mod-3 bölmesinde tam
  simetrik göründüğünü gösteriyor — sıfır-yapı sınavı beklendiği gibi "düz".
- **Bölünen sınıf her beş hedefte de anlamsız (0-3 çizgi, κ≈0)** — pay
  yasasının "payadan hariç" tanımı bu ölçekte pratikte otomatik, çünkü
  q'~e^{13-14} civarında 2^k/3^k/5^k/7^k gibi tekil kuvvetler δ=0.03
  penceresine neredeyse hiç düşmüyor.

---

## TÜRETİLEN vs ÖLÇÜLEN

**Türetilen (veri-öncesi, kalem; 192f'nin ön-kayıtsız keşfinden türetildi):**
- φ0 ≡ 0 (teftiş düzeltmesiyle, madde 2).
- Pay yasası s_r = cos(2πrb/a)/μ(a) (madde 3-4).
- +log10, +log7 öngörü tabloları — 192'de HİÇ görülmemiş, tamamen kör hedefler.

**Körlük durumu:** +log10 ve +log7 kesinlikle kördü — 192'nin ön-kayıtsız
keşfi yalnız +log3/−log3/+log5/+log6'ya bakmıştı (192f). a=10 ve a=7
sınıflandırması bu kalemde ilk kez yapıldı, hiçbir profile bakılmadan K0
donduruldu (komşu raporu dahil, saf aritmetik). +log5/+log3/−log3 doğrulayıcı/
kontrol statüsündeydi (192'de görülmüştü), buradaki katkıları yeni-hedeflerin
tutarlılığını sağlamlaştırmak.

**Ölçülen:** beş uydunun sınıf κ_r, s_r, z değerleri; bölünen sınıfların
neredeyse-boşluğu; +log10/+log7'nin karşılaştırmalı büyüklük davranışı; 8-blok
jackknife se'leri; tam-örneklem/makine mührü kontrolleri.

**Türetilmemiş ve açık kalanlar:**
- +log7/+log5'in sistematik büyüklük fazlalığı (yön doğru, ~1.3-2.3×) —
  S_Re vekilinin %90-95 sınırından mı, yoksa a=5,7'ye özgü bir ikinci
  mertebe etkiden mi? SINANMADI.
- Kirlenme riskli komşuların (2δ içi) gerçek etkisi ölçülmedi — yalnız
  etiketlendi (ön-kayıt kuralı böyle).
- α=1.30 bağı; bu kalemde dokunulmadı.

---

## MANŞET (aday cümle)

> **Kalıntı-sınıfı yasası kör sınavı GEÇTİ: +log10'un dört sınıfı da
> cos(2πr/10) işaretiyle, öngörünün ±0.11 içinde ve z=11-52 arasında dizildi**
> (H-193a MÜHÜR). +log7'nin dört sınıfı da doğru YÖNDE ve anlamlı (z=4.9-6.2)
> ama büyüklükçe öngörünün ~1.5-1.9 katı — nicel bant tutmadı (H-193b KAYIT,
> ölüm yok). Kontroller temiz: +log3 (0.46/0.54) ve özellikle −log3
> (0.499/0.501, sıfır-yapı) ½/½'ye ±0.1 bandında oturdu. Bölünen sınıf (a'yı
> bölen asal-kuvvetler) bu ölçekte pratikte boş çıktı. 192'nin ön-kayıtsız
> keşfi (cos deseni, dönme yok) böylece en keskin kör biçimiyle doğrulandı;
> açık kalan tek soru büyüklüğün neden a=7,5'te fazla çıktığıdır.

---

Teslim:
- Bu rapor.
- `193_configs/`:
  - 193a_onkayit
  - 193b_sinif
  - 193c_hukum
  - 193d_figur
- `193_sinif_yasasi.png`: beş uydu için sınıf payları s_r (çubuk, ± jk se) ve
  KALEM öngörüsü (kesikli çizgi); renk: yeşil = yön+nicel OK, turuncu = yalnız
  yön OK, kırmızı = yön HAYIR/bant dışı, mavi = yan-KAYIT (hükme girmez).
- `scratchpad/193/`:
  - ONKAYIT_193.json (komşu raporu dahil)
  - log_193b.txt
  - muhur_193.json
  - proj_plog10.npz, proj_plog7.npz, proj_plog5.npz, proj_plog3.npz, proj_mlog3.npz
  - tam_plog10_r1.npz
  - K1_sinif.json
  - HUKUM_193.json
