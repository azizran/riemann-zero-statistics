# 195 — UYDULARIN KARAKTERİ: L-fonksiyonu adalarında tarak iptali

**Soru (KALEM_UYDU_KARAKTERI_24EYL2026, commit 50df21e):** Zeta'da (188-193)
pencere-ötesi iptal, orta-nokta örgüsünün Bragg tarağı ve uydularıyla taşınıyor.
Uydular Δω = log(a/b) konumunda duruyor, çekirdek koherensi μ(a)/φ(a) ile
geliyor. Dirichlet L(s,χ) adalarında (Not 5 sıfırları, 118b: β=χ₄, χ₃, χ₅ₑ, χ₈ₑ,
χ₈ₒ) bu mekanizma nasıl değişir? Kalemin veri-öncesi resmi (R3) üç şey söylüyor:
- k'yı bölen bir asal içeren uydular SÖNER.
- İzinli uydular her iki paritede zeta ile AYNI işaretlidir.
- Çekirdek büyüklüğü √k/φ(k) ile ölçeklenir.

İki rakip resim var. R1 iletkeni unutur: adalarda iptal hiç olmaz. R2 pariteyi
unutur: tek adalarda faz 90° olur, gerçel kısım sıfırlanır.

Her sayının yanında onu üreten betik adı vardır. Tek dalga koşuldu. Ölümler
kurtarılmadı. Ön-kayıtlı kapılar atlanmadı. Git'e dokunulmadı. Sonuç ORTAK
TEFTİŞE sunulur, commit sonra.

**İşaret referansı:** hipotezler "K̃_ζ ile AYNI işaret" diye ifade edildi
(teftiş (g)). Ölçülen K̃_ζ yanan uydularda NEGATİF çıktı (∠180°). Yani kalemin B
paragrafındaki "K̃ < 0" ifadesi de doğruymuş [195k0c].

---

## Ölçüm tanımı (K0d'de donmuş; kaynak `195o_ortak.py`)

- **Seri:** çizgi ağırlığı a_q^χ = χ(q)Λ(q)/(π√q log q); χ(q)=0 olan çizgi yok.
  Seri s(n) = Σ 2a_{q'}^χ sin(ω'g_n/2) cos(ω'm_n), 188b.seri_ve_G'nin dds yolu
  AYNEN. İzdüşüm c_{q,b} = 2⟨s·e^{−iω_q m}⟩_b, 190b.izdusum_blok AYNEN.
- **Öz-terim:** c^öz_{q,b} = 2⟨2a_q^χ sin(ω_q g/2) cos(ω_q m) e^{−iω_q m}⟩_b,
  185b.c_oz ile aynı nicelik.
- **Bloklar:** adalarda üst bölge L_χ = log(kt/2π) ≥ 8.5 alındı. Bloklar
  8.5-çapalı, eşit ΔL = 0.05 genişliklidir (üstteki kısmi blok ≥ 0.0125 ise
  dahil). L_b blok içi ortalama L_χ'dir.

  | ada | blok | N (orta nokta) |
  |---|---|---|
  | β | 38 | 65 334 |
  | χ₃ | 34 | 68 034 |
  | χ₅ₑ | 41 | 65 546 |
  | χ₈ₑ | 48 | 61 482 |
  | χ₈ₒ | 48 | 61 482 |

  [195a]
- **Çizgi seçimi:** her blokta pencere çizgileri τ = log q/L_b ∈ [0.45, 0.86),
  χ(q) ≠ 0. Uydu çizgileri |log q' − L_b − log n| < 0.03, χ(q') ≠ 0.
- **Çekirdek:** K̃ = Σ_b N_b·num_b / Σ_b N_b·den_b. Burada
  num_b = Σ_q c^{(uydu)}_{q,b} conj(c^öz_{q,b}) ve den_b = Σ_q |c^öz_{q,b}|².
  Blok-loo jackknife. İşaretli istatistik Re K̃'dir (z = Re/se_jk); Im ve açı
  KAYIT.
- **A (güç):** Ĝ_b(Δω) = mean e^{i(L_b+Δω)m_n}, MUTLAK m_n ile, Δω ∈ [−2.0, 2.6]
  aralığında 0.001 adımla hesaplandı.
  - Pencere |Δω − log n| < 0.03.
  - Halka 0.06–0.25; a,b ≤ 6 rasyonel komşuların ±0.035 çevresi halkadan
    dışlandı.
  - D_b = pencere − halka ortalaması. D = N_b-ağırlıklı ortalama, z = D/se_jk.
- **L-eşli bant (B3):** L ∈ [10.0, 10.4), ızgara blokları j=30..37. Ada verisi
  10.4'ten önce bitiyorsa üst sınır U_χ = L_max olur ve zeta referansı AYNI
  kısmi blokla alınır: β U = 10.368, χ₃ U = 10.176 [195a].
- **B'den dışlanan:** −log5. Uydu çizgileri her ada bloğunda τ' < 0.86 bölgesine,
  yani pencerenin içine düşüyor; bu öz-terim kirliliği demek. Dışlama ölçümden
  önce, aritmetikle donduruldu. −log5 A'da ölçüldü [195a].
- **Zeta blokları:**
  - son: 193 blokları; seri ve izdüşüm ALT=3 (193 AYNEN), öz-terim tam blok.
  - düşük: oran kapısı için 190 blokları, B3 referansı ve A kontrolü için
    adalarla AYNI ΔL ızgarası.

---

## K0a — C_χ MUTLAK KALİBRASYONU  [195k0ab_kapilar.py, 195k0a2_kusur.py]

Beş dosya da t ≈ 200'den başlıyor. Yani ilk sıfır GERÇEK ilk sıfır DEĞİL.
İndeks kaydırması j0 (z_1'den önceki sıfır sayısı) üç yoldan belirlendi:
- (i) ortalama N̄_χ(m_i) − i yuvarlaması,
- (ii) tek nokta N̄_χ(z_1) − ½,
- (iii) BAĞIMSIZ sayım: mpmath Hardy Z_χ(t) işaret değişimleri, (0, m_1)
  aralığında adım 0.02, yakın-kaçış incelemesiyle. max|Im Z|/|Z| ≤ 4.1e-19.

| ada | z_1 | j0 (i) | j0 (ii) | j0 (iii, bağımsız) | kusur | düzeltilmiş mean[N̄−n] TÜM / ÜST | kapı |
|---|---|---|---|---|---|---|---|
| β | 200.170 | 122 | 122 | 122 | yok | −0.00001 / −0.00001 | GEÇTİ |
| χ₃ | 201.589 | **116** | 114 | 114 | **t≈246.30, 246.42 eksik** | −0.00000 / +0.00000 | GEÇTİ |
| χ₅ₑ | 200.201 | 129 | 129 | 129 | yok | −0.00001 / +0.00000 | GEÇTİ |
| χ₈ₑ | 200.054 | 144 | 144 | 144 | yok | −0.00000 / −0.00001 | GEÇTİ |
| χ₈ₒ | 201.010 | **147** | 145 | 145 | **t≈242.55, 242.69 eksik** | +0.00000 / +0.00001 | GEÇTİ |
| ζ son | — | n0 = 1 701 053 (sıfır dosyası indeksi) | | | — | −0.00000 | GEÇTİ |
| ζ düşük | — | n0 = 200 001 | | | — | −0.00000 | GEÇTİ |

**Veri kusuru (açıkça raporlanıyor).** χ₃ ve χ₈ₒ'da (i) ile (ii)/(iii) arasında
2'lik TAMSAYI fark çıktı. Kaynağı tek bir büyük boşluk:
- χ₃'te orta-nokta 32'de 3.947'lik boşluk (≈ 3 ortalama aralık).
- χ₈ₒ'da orta-nokta 36'da 3.109'luk boşluk.

mpmath her birinin içinde dosyada OLMAYAN birer yakın sıfır çifti buldu:
χ₃ 246.303/246.415 ve χ₈ₒ 242.545/242.691 [195k0a2]. İki kusur da L ≈ 4.8'de,
yani üst bölgenin (L ≥ 8.5) çok altında. Dosyanın geri kalanında başka tamsayı
sıçrama yok.

Kusur düzeltilince, yani sıçramadan sonra n'ye 2 eklenince, kalibrasyon
|r| ≤ 1e-5 veriyor. **LİTERAL sayılar da KAYIT:** tek-nokta j0 ile, düzeltmesiz
tüm dosya ortalaması χ₃'te +1.9992, χ₈ₒ'da +1.9989 çıkıyor. Bu tamsayı kayması
C_χ ile ilgisiz; orta-nokta fazını değiştirmez.

A ve B ölçümleri n'yi KULLANMAZ; yalnız mutlak m_n ve g_n kullanılır. Kusurlar
ölçüm bölgesinin dışında. Kapının amacı olan C_χ yarım-tamsayı hatası yok.
**Yanlış parite kontrolü** (C → −C) her adada ±0.250 verdi; kapı bu farkı
ayırıyor [195k0ab]. **Hüküm: K0a GEÇTİ.**

*Süreç notu:* 195k0a2 ön-kayıttan ÖNCE üç kez koşuldu. Birinci koşuda aynı
boşluk iki kez sayılmıştı. İkincide χ₈ₒ'nın sağ-kenar sıfırı (dosya 244.609,
mpmath 244.599) 0.01 payında yanlışlıkla "eksik" sayıldı. Eşleşme "aralıktaki
mpmath değişimi − aralıktaki dosya sıfırı" olarak düzeltildi. Hiçbir ada
ölçümüne bakılmadı.

---

## K0b — KARAKTER KİMLİĞİ  [195k0ab_kapilar.py]

118b, 118a'yı exec ediyor (metinde doğrulandı). χ tabloları 118a'dan okundu.
χ₅ₑ/χ₈ₑ/χ₈ₒ tabloları ilk üretici 105b ile de karşılaştırıldı. Her adada
şunlar assert edildi: değerler, gerçellik, tam-çarpımsallık, parite (χ(k−1)),
İLKELLİK (her öz bölen d|k için m ≡ 1 (d) üzerinde χ ≢ 1), npz 'q','a' ve Gauss
toplamı.

| ada | k | parite | χ (1..k−1) | τ(χ) | 105b |
|---|---|---|---|---|---|
| β (χ₄) | 4 | tek | (1, 0, −1) | 2i = i√4 | — |
| χ₃ | 3 | tek | (1, −1) | i√3 | — |
| χ₅ₑ | 5 | çift | (1, −1, −1, 1) = (·/5) | √5 | ✓ |
| χ₈ₑ | 8 | çift | (1,·,−1,·,−1,·,1) | √8 | ✓ |
| χ₈ₒ | 8 | tek | (1,·,1,·,−1,·,−1) | i√8 | ✓ |

χ₈ₑ = (1,−1,−1,1) ve χ₈ₒ = (1,1,−1,−1) (1,3,5,7'de), kalemle aynı. İkisi de
ilkel; mod 4'ten indüklenmiş değil, çünkü χ(5) = −1. **Hüküm: K0b GEÇTİ**
[195k0ab].

---

## K0c — ZETA KAPISI (yeni kod, χ ≡ 1, k = 1)  [195k0c_zeta.py]

| mühür | sonuç |
|---|---|
| M1: 195o.seri vs 188 `dilim_44.npz` (son, alt-örneklem, 2070 çizgi) | dds, q, a, τ **BİT-BİT** (maks\|Δ\| = 0.0) |
| M2a: 195o.izdusum_blok vs 190b.izdusum_blok | **BİT-BİT** |
| M2b: 195o.c_oz_blok vs 185b.c_oz (düşük zinciri, 1083 çizgi, tüm pencere) | göreli maks fark **1.5e-15** |
| eski yolun 193'ü yeniden üretimi (κ_top + bölünen) | +log3: 0.01925655 (fark 2.5e-15), +log10: 0.01401153 (fark 4.3e-15) |

**Oran kapısı (±%10).** ρ(n) = K(+log n)/K(+log2). Eski değer mix-normalize κ'dan
(193b.blok_yerel_D + 188b AYNEN), yeni değer öz-normalize Re K̃'dan geliyor:

| pencere | uydu | ρ_eski | ρ_yeni | sapma | |
|---|---|---|---|---|---|
| son | +log3 | 0.3224 | 0.3207 | −0.5% | TUTTU |
| son | +log6 | 0.6235 | 0.6204 | −0.5% | TUTTU |
| son | +log10 | 0.2346 | 0.2337 | −0.4% | TUTTU |
| düşük | +log3 | 0.3085 | 0.3088 | +0.1% | TUTTU |
| düşük | +log6 | 0.6057 | 0.6023 | −0.6% | TUTTU |
| düşük | +log10 | 0.2174 | 0.2141 | −1.5% | TUTTU |

**İşaret referansı:** Re K̃_ζ (son). Düşük pencerede de aynı işaret var; tek
istisna −log5, o da düşük pencerede pencere-içi olduğu için zaten dışlandı
[195k0c].

| uydu | +log2 | −log2 | +log3 | −log3 | +log4 | +log5 | +log6 | +log10 | +log(3/2) |
|---|---|---|---|---|---|---|---|---|---|
| K̃_ζ son | −0.05174(69) | −0.01307(51) | −0.01660(41) | −0.00813(30) | **+0.00248(40)** | −0.00525(35) | −0.03210(56) | −0.01209(41) | −0.00754(25) |
| K̃_ζ düşük (ΔL ızgarası) | −0.06664(98) | −0.01737(33) | −0.02088(27) | −0.01056(26) | **+0.00317(14)** | −0.00642(26) | −0.04049(44) | −0.01435(36) | −0.00952(26) |

Yanan uyduların hepsinde ∠K̃ = 180° ± 3.4° çıktı. Açı istisnası +log4
(μ(4) = 0): orada pozitif bir çukur var, 192'nin +log4 çukurunun (κ < 0) aynısı.
**Hüküm: K0c GEÇTİ** [195k0c].

---

## K0d — DONMUŞ ÖN-KAYIT  [195a_onkayit.py]

`scratchpad/195/ONKAYIT_195.json`:
- **sha256 = f1e1fa8254b4656d…** (tam:
  f1e1fa8254b4656d9c9bb720069cfa4795242f68d32750733f5d9cb9dadaa234)
- **damga Thu Sep 24 14:46:48 +0300 2026**

Sıralama dosya zamanlarından okunur:
1. K0ab 14:38:32
2. K0c 14:43:50
3. K0a_kusur 14:46:16
4. **ÖN-KAYIT 14:46:48**
5. A 14:47:39
6. B 14:49:18

Betik şu durumlarda yazmayı reddeder: K0 kapılarından biri kalmışsa, bir kusur
üst bölgedeyse, ya da ONKAYIT zaten varsa. İlk denemesi bir taşma hatasıyla
dosya YAZMADAN düştü (ΠJ₀ döngüsü; düzeltildi).

Donanlar:
- Kalemin öngörü tablosu ve asal-çarpan kuralından türeyen 5 ada × 10 uydu
  yasak/izinli matrisi. Kalem listeleriyle assert edildi.
- Tüm tanımlar, ada blok listeleri (kinematik), A ızgarası, pencere ve halka.
  Halkada 143–380 nokta var.
- B uydu listesi ve −log5 dışlaması.
- H-195A/B1/B2/B3 karar kuralları; R1/R2/R3 imzaları.
- K0c işaret referansı ve L-eşli zeta bant referansları. Zeta referansı bu
  noktada donduruldu, adalardan ÖNCE.
- K0a/b/c özetleri ve 21 girdi dosyasının sha256'sı.
- **Ajan veri-öncesi notu (HÜKME GİRMEZ).** Ada/zeta oranının kalemdeki
  √k/φ(k)'ye iki ek çarpan alabileceği yazıldı:
  - (i) L-eşli kıyasta durağan-faz noktası t*_χ = t*_ζ/k. Blok Δt = t·ΔL
    olduğundan blok-ortalamalı Ĝ çizgi başına √k büyür, oran k/φ(k) olur.
  - (ii) χ(q)=0 çizgilerinin Jacobi–Anger J₀(2πa_q) çarpanları adada 1 olur.
    Bu, 1/Π_{p|k} J₀ çarpanı getirir: ΠJ₀ = 0.514 (k=4,8), 0.671 (k=3),
    0.801 (k=5).
  - Alternatif sayılar: (k/φ(k))/ΠJ₀ = 3.888 (β, χ₈), 2.236 (χ₃), 1.561 (χ₅ₑ).

  B3 kalemdeki gibi √k/φ(k) ile sınandı; alternatif yalnız KAYIT.

---

## A — YAPI ÇARPANI GÜCÜ  [195b_guc.py]

z = D/se_jk (blok-loo). **Y** = yasak (n'nin bir asal çarpanı k'yı böler).
+log(3/2) KAYIT. Ada satırları üst bölgeden, zeta satırları kontrol amaçlı.

| küme | +log2 | −log2 | +log3 | −log3 | +log4 | +log5 | −log5 | +log6 | +log10 | +log(3/2) |
|---|---|---|---|---|---|---|---|---|---|---|
| β (k=4) | −11.3 **Y** | −11.8 **Y** | **+11.0** | **+10.8** | −7.3 **Y** | +11.3 | +11.2 | −11.6 **Y** | −9.8 **Y** | −9.0 Y |
| χ₃ (k=3) | **+11.5** | **+11.2** | −11.6 **Y** | −12.0 **Y** | −12.5 | +12.1 | +11.8 | −12.3 **Y** | +11.6 | −8.4 Y |
| χ₅ₑ (k=5) | **+10.2** | **+9.8** | **+10.3** | **+9.9** | −10.9 | −8.7 **Y** | −10.8 **Y** | +10.4 | −7.2 **Y** | +11.0 |
| χ₈ₑ (k=8) | −9.9 **Y** | −10.3 **Y** | **+10.1** | **+9.9** | −6.9 **Y** | +10.3 | +10.3 | −9.7 **Y** | −9.0 **Y** | −8.4 Y |
| χ₈ₒ (k=8) | −10.8 **Y** | −10.5 **Y** | **+10.1** | **+9.9** | −7.4 **Y** | +10.3 | +10.2 | −9.6 **Y** | −8.6 **Y** | −8.5 Y |
| ζ son | +126.5 | +85.0 | +98.1 | +55.9 | **−33.0** | +36.5 | +32.4 | +78.7 | +50.7 | +38.9 |
| ζ düşük | +15.2 | +14.4 | +15.0 | +13.7 | **−15.2** | +17.0 | +12.5 | +15.5 | +15.6 | +15.6 |

Gürültü-tabanı biriminde (N_b·D) değerler şöyle [195b]:
- İzinli-güçlüler: +6 … +30.
- Yasaklar: −0.16 … −0.80.
- Profilde yasak konumlarda HİÇ tepe yok. Örneğin β'nın log2 penceresi 0.15 taban
  birimi, komşuları 0.2–0.6; izinli log3 tepesi 25.5.

**Dürüst okuma (A):**
- **Seçim kuralı güçte tuttu.** A listesinde 21 yasak (ada × uydu) var; 21'i de
  sönük (+log(3/2) KAYIT ile 25/25). 24 izinliden 22'si +9.8σ ile +12.1σ
  arasında yanıyor. χ₃'te log10, χ₅ₑ'de log(3/2) ve log6 izinli; bunlar da
  yanıyor. Yanmayan iki izinli, χ₃ ve χ₅ₑ'deki +log4: μ(4)=0 ve zeta'da da güç
  yok (aşağıda KAYIT).
- **Yasakların z'si NEGATİF (−7 … −12σ), sıfır değil.** Halka tabanı, uydusuz bir
  pencerenin tabanından sistematik olarak yüksek. Nedenler: halkaya düşen a,b > 6
  zayıf rasyoneller ve rasyonel konumlardaki yerel bir çukur. Zeta'nın +log4'ü
  de −33σ veriyor. Ön-kayıtlı kural tek yönlü (< 2σ) olduğu için hüküm
  etkilenmiyor. Ama bu z, sıfır-merkezli kalibre bir null testi DEĞİL; teftiş
  notu.
- **KAYIT — kalemin A notu veriyle çelişiyor.** Kalem "μ kuralı burada YOK — ör.
  zeta'da log4 güçte görünür" diyordu. Ölçümde zeta'da log4 gücü YOK: son −33σ,
  düşük −15σ, profil merkezinde 0.28 taban birimi. Karakterce izinli olduğu
  χ₃ (−12.5) ve χ₅ₑ (−10.9) adalarında da yok. Güç de karesiz-benzeri bir kurala
  uyuyor olabilir. **SINANMADI.**

---

## B — ÇEKİRDEK K̃  [195c_cekirdek.py]

### B-i. Hüküm çiftleri (üst bölge; Re K̃ ± jk se; z)

| rol | ada (parite) | uydu | Re K̃ | z_Re | ∠K̃ | Im z | ζ işareti |
|---|---|---|---|---|---|---|---|
| izinli-birincil | χ₃ (tek) | +log2 | **−0.1720 ± 0.0023** | −73.8 | +179.9° | +0.5 | − (aynı) |
| izinli-birincil | χ₅ₑ (çift) | +log2 | **−0.1131 ± 0.0018** | −62.0 | −180.0° | −0.1 | − (aynı) |
| izinli-birincil | χ₅ₑ (çift) | +log3 | **−0.03717 ± 0.00067** | −55.2 | +179.8° | +0.4 | − (aynı) |
| izinli-birincil | β (tek) | +log3 | **−0.1023 ± 0.0013** | −81.9 | −179.9° | −0.2 | − (aynı) |
| B2 birincil | χ₈ₑ (çift) | +log3 | −0.09450 ± 0.00154 | −61.6 | −179.9° | −0.4 | − (aynı) |
| B2 birincil | χ₈ₒ (tek) | +log3 | −0.09539 ± 0.00160 | −59.5 | +179.7° | +0.8 | − (aynı) |
| yasak-birincil | χ₃ | +log3 | **+0.00503 ± 0.00042** | **+11.9** | +0.9° | +0.3 | ZIT |
| yasak-birincil | χ₅ₑ | +log5 | **+0.00420 ± 0.00040** | **+10.6** | +2.7° | +0.9 | ZIT |
| yasak-birincil | β | +log2 | **+0.00594 ± 0.00053** | **+11.3** | −0.5° | −0.2 | ZIT |

### B-ii. Tam tablo (KAYIT; üst bölge Re K̃ ± se, parantezde z_Re; **Y** = yasak)

| uydu | β (k=4, tek) | χ₃ (k=3, tek) | χ₅ₑ (k=5, çift) | χ₈ₑ (k=8, çift) | χ₈ₒ (k=8, tek) |
|---|---|---|---|---|---|
| +log2 | **Y** +0.00594±53 (+11.3) | −0.17200±233 (−73.8) | −0.11311±182 (−62.0) | **Y** +0.00479±55 (+8.7) | **Y** +0.00517±51 (+10.2) |
| −log2 | **Y** +0.00629±46 (+13.8) | −0.04861±78 (−62.0) | −0.03094±77 (−39.9) | **Y** +0.00614±64 (+9.6) | **Y** +0.00502±48 (+10.4) |
| +log3 | −0.10227±125 (−81.9) | **Y** +0.00503±42 (+11.9) | −0.03717±67 (−55.2) | −0.09450±154 (−61.6) | −0.09539±160 (−59.5) |
| −log3 | −0.05417±87 (−62.2) | **Y** +0.00605±43 (+14.2) | −0.01875±73 (−25.7) | −0.05064±95 (−53.5) | −0.05075±83 (−61.3) |
| +log4 | **Y** +0.00378±50 (+7.6) | +0.00548±43 (+12.9) | +0.00447±47 (+9.5) | **Y** +0.00374±50 (+7.4) | **Y** +0.00378±56 (+6.8) |
| +log5 | −0.03633±70 (−51.8) | −0.01885±61 (−30.9) | **Y** +0.00420±40 (+10.6) | −0.03317±69 (−48.1) | −0.03288±87 (−37.7) |
| +log6 | **Y** +0.00473±44 (+10.8) | **Y** +0.00534±43 (+12.4) | −0.06874±100 (−68.7) | **Y** +0.00287±45 (+6.4) | **Y** +0.00399±48 (+8.4) |
| +log10 | **Y** +0.00220±54 (+4.0) | −0.03846±60 (−64.2) | **Y** +0.00255±45 (+5.7) | **Y** +0.00144±50 (+2.9) | **Y** +0.00226±48 (+4.7) |
| +log(3/2) | **Y** +0.00505±55 (+9.1) | **Y** +0.00467±51 (+9.2) | −0.01647±66 (−25.0) | **Y** +0.00520±48 (+10.9) | **Y** +0.00474±52 (+9.1) |
| −log5 | dışlandı (pencere-içi) | dışlandı | dışlandı | dışlandı | dışlandı |

(se son hanelerle yazıldı: "±53" = ±0.00053.) Bütün izinli (μ≠0) uydular
∠180° ± 1.3° ve |Im z| ≤ 2.1 veriyor; yani pariteden bağımsız, gerçel ve zeta
işaretli. Bütün yasaklar (−log5 hariç) ∠0° ± 7.2° ve +0.0014…+0.0063
aralığında [195c].

### B-iii. L-eşli bant oranları  R = Re K̃_χ,bant / Re K̃_ζ,düşük,bant  [195c, 195d]

| ada | uydu | K̃_χ,bant | K̃_ζ,bant | R ± se | kalem √k/φ(k) | bant ±%35 | R / [(k/φ(k))/ΠJ₀] (KAYIT) |
|---|---|---|---|---|---|---|---|
| β | +log3 | −0.09529±0.00101 | −0.02202±0.00025 | **4.328 ± 0.067** | 1.000 | [0.650, 1.350] — DIŞINDA (44σ) | 1.113 |
| χ₃ | +log2 | −0.15736±0.00242 | −0.07249±0.00073 | **2.171 ± 0.040** | 0.866 | [0.563, 1.169] — DIŞINDA (25σ) | 0.971 |
| χ₅ₑ | +log2 | −0.10709±0.00275 | −0.07065±0.00060 | **1.516 ± 0.041** | 0.559 | [0.363, 0.755] — DIŞINDA (19σ) | 0.971 |
| χ₅ₑ | +log3 | −0.03563±0.00087 | −0.02178±0.00037 | **1.636 ± 0.048** | 0.559 | [0.363, 0.755] — DIŞINDA (18σ) | 1.048 |
| χ₈ₑ | +log3 | −0.09394±0.00190 | −0.02178±0.00037 | **4.314 ± 0.114** | 0.707 | [0.460, 0.955] — DIŞINDA (30σ) | 1.110 |
| χ₈ₒ | +log3 | −0.09808±0.00137 | −0.02178±0.00037 | **4.504 ± 0.099** | 0.707 | [0.460, 0.955] — DIŞINDA (36σ) | 1.158 |

Bant içeriği şöyle:

| ada | bant blokları | N |
|---|---|---|
| β | 8 | 24 990 |
| χ₃ | 4 | 14 241 |
| χ₅ₑ | 8 | 22 129 |
| χ₈ₑ | 8 | 13 830 |
| χ₈ₒ | 8 | 13 829 |

χ₃'ün jk se'si 4 bloktan geldiği için kaba. Altı çiftin altısı da ölçülebilir:
ada ve ζ tarafında |z| ≥ 3 [195d].

**KAYIT (diğer izinli uydular, aynı bant):**

| ada | uydu → R |
|---|---|
| χ₅ₑ | −log2 1.567, −log3 1.495, +log6 1.598, +log(3/2) 1.600 |
| χ₃ | −log2 2.306, +log5 2.741, +log10 2.641 |
| β | −log3 4.600, +log5 5.679 |
| χ₈ₑ | −log3 4.610, +log5 5.513 |
| χ₈ₒ | −log3 4.577, +log5 5.879 |

Oran ada başına neredeyse sabit. χ₅ₑ'nin altı izinli uydusu 1.50–1.64 aralığına
düşüyor, ajan notu 1.56 [195d].

---

## HÜKÜM (eşikler K0d'de donmuş; kurtarma yok)  [195d_hukum.py]

| hipotez/kapı | hüküm | dayanak |
|---|---|---|
| **K0a** C_χ kalibrasyonu | **GEÇTİ** | düzeltilmiş \|r\| ≤ 1e-5, 5 ada + 2 zeta. χ₃/χ₈ₒ'da L≈4.8'de birer eksik yakın çift (üst bölge dışı); literal düzeltmesiz r = +1.999 KAYIT [195k0ab, 195k0a2] |
| **K0b** karakter kimliği | **GEÇTİ** | tablolar, ilkellik, parite, τ(χ), 105b eşleşmesi [195k0ab] |
| **K0c** zeta kapısı | **GEÇTİ** | seri bit-bit; 193 2.5e-15/4.3e-15; oran sapmaları ≤ %1.5 [195k0c] |
| **K0d** ön-kayıt | **GEÇTİ** | sha f1e1fa82…, 14:46:48; A 14:47:39'da, B 14:49:18'de [195a] |
| **H-195A** güç seçimi | **MÜHÜR** | yasakların hepsi z < 2 (maks −6.9); izinli-güçlüler z > 3 (min +9.8); ölüm (yasak ≥ 4σ) yok. UYARI: yasak z'leri negatif, null kalibre değil (yukarıda) [195b, 195d] |
| **H-195B1** varlık + işaret | **KAYIT** | 4 izinli-birincil ζ ile AYNI işaretli, z −55…−82 → ölüm yok. AMA 3 yasak-birincil \|z\| < 2 DEĞİL: +10.6…+11.9σ, zeta'ya ZIT işaretli küçük artık (+0.0042…+0.0059, izinlilerin ~%3-5'i) → MÜHÜR değil [195c, 195d] |
| **H-195B2** parite | **MÜHÜR** | tek (β, χ₃, χ₈ₒ) ve çift (χ₅ₑ, χ₈ₑ) birincillerin hepsi ζ işaretli, \|z\| ≥ 55; ρ₈ = χ₈ₒ/χ₈ₑ = **1.009 ± 0.024** ∈ [0.60, 1.40]; R2 ölüm deseni yok; tek adalarda \|Im z\| ≤ 0.8 [195d] |
| **H-195B3** L-eşli ölçek | **ÖLDÜ** | 6/6 ölçülebilir çift bant DIŞINDA, 18-44 se uzakta; ölçülen R = 1.52–4.50, öngörü 0.56–1.00 [195d] |
| **Rakip resimler** | **R3 (işaret/parite/seçim) kazandı; R3'ün ölçeği öldü** | R1 ölü (izinliler −55…−82σ); R2 ölü (tek adalarda faz 180°, Re ≠ 0); R3'ün √k/φ(k) ölçeği ölü [195d] |

**Dürüst okuma:**
- **Kalemin nitel resmi (R3) tam tuttu.** k'yı bölen asallı uydular hem güçte hem
  çekirdekte söndü; izinliler her iki paritede zeta'nın işaretinde yandı.
  χ₈ₑ/χ₈ₒ çifti izinli uydularda %1 içinde aynı (+log3, −log3, +log5). Karakter işaretlerinin sadeleşmesi
  (kalem madde 5) ve Gauss-toplamı fazının durağan-faz fazını götürmesi
  (tek: i·e^{−iπ/2} = 1) veriyle tutarlı.
- **"Sönme = sıfır" TUTMADI; sönme "zıt işaretli küçük artık".** Yasak konumların
  hepsi +0.0014…+0.0063 veriyor (−log5 hariç). Bu, izinlilerin ~%2-6'sı büyüklüğünde ve zıt
  işaretli. Aynı işaret ve düzey zeta'nın μ(a)=0 konumu +log4'te de var
  (+0.0025/+0.0032; 192'nin log4 çukuru). Karakterce izinli ama μ(4)=0 olan
  χ₃/χ₅ₑ +log4'te de görülüyor (+0.0055/+0.0045). Yani "yasak" uydu, zeta'nın
  μ=0 uydusuyla AYNI davranıyor. Ön-kayıtlı kural "|z| < 2" dediği için B1
  MÜHÜR değil, KAYIT.
- **Ölçek kalemden 2.7-6.4 kat büyük.** B3'ün ölümü kesin (18-44 se).

---

## ÖN-KAYITSIZ KEŞİF (hükümler görüldükten SONRA; HÜKÜM DIŞI)  [195f_kesif.py]

**Soru:** Yasak konumlardaki pozitif artık evrensel bir "boş konum" tabanı mı?

Aynı makineyle 18 "boş" Δω konumunda K̃ ölçüldü (a,b ≤ 10 rasyonellerinden
≥ 0.045 uzak). Bulgular:
- Güçlü uyduların omzundaki konumlar (≤ 0.05 uzak) güçlü NEGATİF çıktı, −0.01 …
  −0.027. Nedeni ada blok yayılımı 0.05 + δ 0.03 = 0.055 omuz genişliği.
- Omuzdan uzak 9 konumun ortalaması:

  | küme | uzak-boş ort. | ort. z |
  |---|---|---|
  | β | +0.00209 | +4.9 |
  | χ₃ | +0.00181 | +4.1 |
  | χ₅ₑ | +0.00152 | +3.0 |
  | χ₈ₑ | +0.00215 | +4.8 |
  | χ₈ₒ | +0.00147 | +2.9 |
  | ζ düşük | −0.00003 | ≈0 |
  | ζ son | −0.00012 | ≈0 |

- Karşılaştırma: ada yasak ortalamaları +0.0034…+0.0053.

**Okuma (SINANMADI):** Adalarda uydusuz konumlarda da zayıf pozitif bir taban
var (~+0.002). Zeta'da bu taban yok. Yasak konumlar bu tabanın yaklaşık 2 katında
duruyor, yani zeta'nın μ=0 çukuru düzeyinde. Yasak artığının bir kısmı "ada
tabanı", bir kısmı "μ=0 benzeri çukur" olabilir. Ada tabanının kaynağı açık:
seyrek χ≠0 çizgi kümesi mi, yoksa ada kinematiğinin kısa t'si (L-eşli
kıyasta t_χ = t_ζ/k) mi? **SINANMADI.**

---

## TÜRETİLEN vs ÖLÇÜLEN

**Türetilen (kalem, veri-öncesi):**
- Seçim kuralı: n'nin bir asal çarpanı k'yı bölüyorsa uydu yok.
- Karakter işaretlerinin sadeleşmesi.
- Faz sabiti: τ(χ)·e^{−i(2πC_χ+π/4)} her iki paritede gerçel ve zeta işaretli.
- Ölçek √k/φ(k).
- A notu: "zeta'da log4 güçte görünür".

**Ölçülen:**
- Seçim kuralı güçte (21/21 yasak sönük, 22/24 izinli yanıyor; yanmayan ikisi
  +log4) ve çekirdekte işaret düzeyinde tuttu. Yasaklar
  izinlilerin ~%3-5'i düzeyinde zıt işaretli artık bırakıyor (türetilmemişti).
- Faz, beş adanın hepsinde 180° ± 1.3° çıktı. Tek/çift farkı yok: ρ₈ = 1.009.
- Ölçek TUTMADI: R = 1.52 (χ₅ₑ) … 4.50 (χ₈ₒ).
- A notu TUTMADI: zeta'da ve χ₃/χ₅ₑ'de log4 gücü yok.

**Ajan veri-öncesi notu (ön-kayıtta donmuş, HÜKME GİRMEZ):**
(k/φ(k))/Π_{p|k}J₀(2πa_{p^j}) = 3.888 / 2.236 / 1.561. Ölçülen R'nin bu sayıya
oranı altı birincil çiftte 0.97–1.16:
- χ₃ 0.971
- χ₅ₑ 0.971 ve 1.048
- β 1.113
- χ₈ₑ 1.110
- χ₈ₒ 1.158

k=4,8 adalarında (2'nin kuvvetleri kaldırılınca) sistematik +%11-16 kalıyor.
Notun iki bileşeni:
- √k: L-eşli kıyasta durağan-faz genliği √(2πt*)/Δt.
- 1/ΠJ₀: χ(q)=0 çizgilerinin Jacobi–Anger arka plan çarpanı.

Bu bir ÖNGÖRÜ olarak sınanmadı; ön-kayıtta not olarak duruyordu. Kaptanın ortak
teftişte değerlendirmesi için. Sayısal uyum dikkat çekici ama hipotez değil.

**Türetilmemiş ve açık kalanlar:**
- Yasak artığının ve ada tabanının kaynağı.
- A'da yasak z'lerinin negatifliği, yani halka tabanının kalibrasyonu.
- log4 gücünün yokluğu.
- k=4,8'de alternatif sayıdan kalan %11-16.
- α=1.30 bağı; bu kalemde dokunulmadı.

---

## MANŞET (aday cümle)

> **Tarak iptali L-fonksiyonu adalarında da var ve karakteri tam kalemin dediği
> gibi seçiyor.** k'yı bölen asallı uydular hem yapı çarpanında (21/21 yasak
> sönük; izinliler +9.8…+12.1σ) hem çekirdekte sönüyor. İzinli uydular beş adanın
> beşinde, tek ve çift paritede, zeta ile AYNI işaretli ve gerçel: ∠180° ± 1.3°,
> z = −55…−82. χ₈ₑ ile χ₈ₒ ayırt edilemiyor (1.009 ± 0.024). R1 (iletkensiz) ve
> R2 (paritesiz) öldü; H-195A ve H-195B2 MÜHÜR. **Ama iki şey kalemden
> sapıyor.** (1) "Sönme" sıfır değil: yasak uydular izinlilerin ~%3-5'i
> büyüklüğünde, zıt işaretli bir artık bırakıyor. Bu, zeta'nın μ=0 çukuruyla
> (+log4) aynı düzeyde; H-195B1 KAYIT. (2) L-eşli ölçek √k/φ(k) ÖLDÜ: ölçülen
> oranlar 1.5-4.5, öngörü 0.56-1.0, fark 18-44σ (H-195B3). Ön-kayıtta HÜKME
> GİRMEZ diye donmuş ajan notu (k/φ(k))/ΠJ₀ ise altı çiftte %3-16 içinde.

---

Teslim:
- Bu rapor.
- `195_configs/`:
  - 195o_ortak (genelleştirilmiş makine)
  - 195k0ab_kapilar
  - 195k0a2_kusur
  - 195k0c_zeta
  - 195a_onkayit
  - 195b_guc
  - 195c_cekirdek
  - 195d_hukum
  - 195e_figur
  - 195f_kesif (hüküm dışı)
- `195_uydu_karakteri.png`:
  - Üst sıra: her ada için N_b·|Ĝ_b(Δω)|² profili (log ölçek, gri = ζ düşük);
    yasak konumlar turuncu kesikli, izinliler mavi noktalı.
  - Alt sıra: L-eşli bantta K̃_χ/K̃_ζ çubukları (izinli mavi, yasak turuncu,
    ± jk se). Kalem √k/φ(k) siyah çizgi ve ±%35 şeridi; ajan notu gri kesikli
    (hükme girmez).
- `scratchpad/195/`:
  - ONKAYIT_195.json
  - K0ab.json, K0a_kusur.json
  - K0c.json, K0c_bloklar.npz
  - A_guc.json, A_profiller.npz
  - B_cekirdek.json, B_bloklar.npz
  - HUKUM_195.json
  - K_kesif_195.json
  - log_195k0ab, log_195k0a2, log_195k0c, log_195a, log_195b, log_195c,
    log_195d, log_195f (.txt)
