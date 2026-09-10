# 187 — PENCERE-ÖTESİ FAZ-İPTAL DEFTERİ (α=1.30'un güncel adresi)

**Soru (KALEM_FAZ_IPTAL_10EYL2026):** 186'nın bonus teşhisi — gerçek deniz,
pencere-içi (τ'≤0.86) merdiven karışımının ~üçte birini pencere-ötesi
içerikle FAZ-UYUMLU iptal ediyor (m^kesik/m_ölç = 1.46→1.58; ikizde
1.04→1.10) — bu kalemde DEFTERE bağlanır: iptali hangi τ'-katmanı taşıyor,
ne kadar derine iniyor, açısı gerçekten 180° mi, ve defter HA4'ün %26'sını
kapatıyor mu?

Her sayının yanında onu üreten betik adı vardır. Tek dalga; ölümler
kurtarmasız; git'e dokunulmadı; sonuç KAPTAN+KULLANICI ortak teftişine
sunulur, commit sonra.

---

## Ölçüm tanımı (K0'da donmuş; 184-186 makinesi AYNEN)

Katman i = τ' ∈ (ızgara[i−1], ızgara[i]] dilimindeki TÜM asal-kuvvetler;
ızgara GERÇEK [0.86, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20], kontrol
İKİZ+HA4 [0.86, 0.90, 1.00, 1.10]. Nominal a_q = Λ(q)/(π√q·log q)
(makine kontrolü: 184 evrenine maks sapma 2.8e-17). Katman-artım serisi
`dds_i(n) = Σ_{q∈i} 2a_q sin(ω_q g_n/2) cos(ω_q m_n)` (ρ≡1); izdüşüm
c^{Δi}_q = 2⟨dds_i e^{−iω_q m}⟩; `c^ya(τ_c) = c^kesik + Σ_{i≤τ_c} c^{Δi}`
(c^kesik = 186 G1_proj re_bir/im_bir AYNEN, yeniden hesap yok);
`m^ya_b(τ_c) = [Σ_b |c^ya| − Σ_b â^öz]/Σ_b ae`. İptal katsayısı:
Δ_q = c^ölç_q − c^kesik_q, karışım^kesik_q = c^kesik_q − c^öz_q,
**ζ_b = Σ_b Δ_q·conj(karışım_q)/Σ_b |karışım_q|²** (bant-içi kompleks
en-küçük-kare oranı; modül = iptal payı, açı = faz sınavı). Bantlar +
8-blok loo-jackknife 184-186 AYNEN; **HAVUZ** = τ∈[0.45,0.86) birleşik
maske (H-F1/H-F3 hükümleri bu sütunda). H-F4: M^pred_HA4 = m^env(τ_c=0.86;
**HA4-kinematiği**)/m_Hk(ölçülü) — erfc-ağırlıklı pencere-içi karışım,
iptalsiz; 186e bu koşuyu yalnız yabancı (Hk/gerçek) kinematiklerde yapmıştı.

---

## K0 — DONMUŞ ÖN-KAYIT  [187a_onkayit.py]

`187/ONKAYIT_187.json` yazıldı: **sha256 = 42a91c3a…**, damga
**Thu Sep 10 15:37:26 +03 2026** (tüm ölçümlerden önce). Donan: katman
ızgaraları; kestirimciler (katman-artım izdüşümü, ζ ve açısı, bant/jk);
H-F1 çatalının iki okuması; H-F2 açı penceresi 180°±15°; H-F3/H-F4
eşikleri; K3 tanımları (ι birincil = −Δm/m^kesik; C rekonstrüksiyonu;
biçim KAYIT makinesi); İKİNCİL sonda tanımı ((1.20,1.30], 100k
alt-örneklem, ana hükümlere girmez).

**Makine mühürleri:** (1) yeniden kurulan taban serisinin momentleri 186
G1_proj(ρ≡1) ile bağıl fark ≤ 5.3e-16 [187b]; (2) m^ya(0.86), 186'nın
g_bir m_pred defterini hane hane verdi (maks|Δ| = 1.1e-16) [187b];
(3) ölçülü M_HA4 defteri 186e'yi hane hane verdi (maks|Δ| = 2.2e-16) [187f];
(4) ikiz-kinematik taban momenti 186 ile birebir (fark 0.0) [187d].

---

## K1 — KATMAN DEFTERİ (gerçek kinematik; ana koşu)  [187b_katman_defteri.py]

Yeni çizgiler: q ≤ e^{1.20·L} = 1.86e6, katman başına 1826 / 3730 / 6461 /
11198 / 19494 / 33971 / 59365 (toplam 136 045; asal-kuvvetler dahil).
Katman serileri kontrol-noktalı kaydedildi (katman_1..7_gercek.npy).

**Yakınsama defteri** m^ya(τ_c; bant) [hedef m_ölç; jk se tablo dosyada]:

| bant | m_ölç | 0.86 | 0.90 | 0.95 | 1.00 | 1.05 | 1.10 | 1.15 | 1.20 | m^ya(1.20)/m_ölç | kapanan pay |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.45-0.50 | 0.3383 | 0.5003 | 0.4904 | 0.4723 | 0.4548 | 0.4435 | 0.4062 | 0.3937 | 0.3784 | 1.118 | 0.753±.015 |
| 0.50-0.55 | 0.3949 | 0.5804 | 0.5691 | 0.5491 | 0.5299 | 0.5166 | 0.4753 | 0.4614 | 0.4449 | 1.126 | 0.731±.011 |
| 0.55-0.60 | 0.4489 | 0.6564 | 0.6465 | 0.6264 | 0.6049 | 0.5898 | 0.5441 | 0.5265 | 0.5084 | 1.133 | 0.713±.007 |
| 0.60-0.65 | 0.5034 | 0.7351 | 0.7239 | 0.7040 | 0.6793 | 0.6620 | 0.6100 | 0.5905 | 0.5706 | 1.133 | 0.710±.008 |
| 0.65-0.70 | 0.5477 | 0.7992 | 0.7896 | 0.7693 | 0.7437 | 0.7285 | 0.6723 | 0.6533 | 0.6298 | 1.150 | 0.674±.012 |
| 0.70-0.75 | 0.5837 | 0.8620 | 0.8511 | 0.8303 | 0.8053 | 0.7889 | 0.7236 | 0.7018 | 0.6764 | 1.159 | 0.667±.011 |
| 0.75-0.80 | 0.6033 | 0.9105 | 0.9045 | 0.8823 | 0.8535 | 0.8343 | 0.7637 | 0.7418 | 0.7134 | 1.182 | 0.642±.008 |
| 0.80-0.86 | 0.5941 | 0.9410 | 0.9417 | 0.9236 | 0.8984 | 0.8814 | 0.7965 | 0.7694 | 0.7368 | 1.240 | 0.588±.007 |
| **HAVUZ** | 0.5294 | 0.7944 | 0.7868 | 0.7668 | 0.7425 | 0.7263 | 0.6648 | 0.6443 | 0.6200 | **1.171** | **0.658±.003** |

**Katman artımları** (HAVUZ; ±jk se): −0.0076±.0008 / −0.0200±.0009 /
−0.0244±.0051 / −0.0162±.0055 / **−0.0615±.0015** / −0.0205±.0052 /
−0.0242±.0053 — hepsi NEGATİF (iptal yönü), son iki katman hâlâ ~2σ'nın
üstünde. **Korelasyon yan sütunu:** korr(ds^ya(τ_c), ds) = 0.9526 (0.86,
186 AYNEN) → 0.9592 (1.00) → **0.9674** (1.20) — beklendiği gibi yükseliyor.

**H-F1 HÜKMÜ:** yaklaşma HAVUZ'da MONOTON (EVET) ve τ_c=1.20'ye dek açığın
**%65.8±0.3'ü kapandı** (≥%50: EVET) → **H-F1 birincil beklenti MÜHÜR.**
Çatal okumaların İKİSİ DE ateşlemedi: (a) SIĞ İPTAL değil — art arda iki
katmanda |Δm|<2se görülmedi, artımlar 1.20'de hâlâ sürüyor; (b) DERİN
KUYRUK değil — kalan açık %34.2 ≤ %50. Ölçülen karakter: **ORTA-DERİN,
yayvan iptal** — üçte ikisi (0.86,1.20]'de ödeniyor, kalan üçte biri
τ'>1.20 kuyruğunda. (Özdeşlik gereği tam-derinlik limiti m_ölç'tür; buradaki
içerik LİMİT değil, YOLUN YAPISIDIR.)

---

## K2 — ζ DEFTERİ: iptal payı ve açısı  [187c_zeta_defteri.py]

| bant | \|ζ\| gerçek | açı gerçek | \|ζ\| ikiz | açı ikiz | genlik-okuma (gerçek) |
|---|---|---|---|---|---|
| 0.45-0.50 | 0.3284±.0032 | +179.55°±0.48 | 0.0892±.0080 | +177.96°±1.85 | 0.3238±.0033 |
| 0.50-0.55 | 0.3230±.0032 | −179.46°±0.65 | 0.0869±.0083 | +178.68°±1.98 | 0.3195±.0033 |
| 0.55-0.60 | 0.3184±.0060 | −179.75°±0.48 | 0.0796±.0074 | −176.49°±2.41 | 0.3161±.0061 |
| 0.60-0.65 | 0.3182±.0030 | +179.33°±0.34 | 0.0796±.0067 | +178.00°±1.54 | 0.3152±.0027 |
| 0.65-0.70 | 0.3184±.0039 | +179.51°±0.55 | 0.0736±.0062 | −179.90°±2.53 | 0.3147±.0039 |
| 0.70-0.75 | 0.3274±.0041 | −179.94°±0.45 | 0.0695±.0046 | −179.66°±2.50 | 0.3229±.0041 |
| 0.75-0.80 | 0.3456±.0052 | −179.31°±0.33 | 0.0685±.0087 | −177.99°±2.10 | 0.3374±.0045 |
| 0.80-0.86 | 0.3811±.0072 | +179.93°±0.38 | 0.0591±.0086 | +175.92°±1.38 | 0.3687±.0045 |
| HAVUZ | 0.3287±.0023 | +179.96°±0.07 | 0.0768±.0064 | +179.76°±0.34 | 0.3335±.0026 |

- **H-F2 MÜHÜR — SAF YIKICI:** gerçek denizde 8/8 bant 180°±15° penceresinin
  içinde; dahası açı her bantta **180°±0.7°** — pencereden 20 kat dar.
  KALEM'in beklediği modül bölgesi (0.31-0.37) birebir tutmuş: |ζ| = 0.318
  → 0.381, τ ile hafif YÜKSELİYOR (iptal payı kuyruğa doğru artar).
- İkizde iptal küçük (|ζ| = 0.089 → 0.059, τ ile DÜŞÜYOR) ama o da 180°'de:
  ikizin kendi sığ kuyruğu (0.86,1.00] az ama uyumlu iptal ediyor.
- Genlik-okuması (1 − m_ölç/m^kesik) ζ modülüyle %1-3 içinde örtüşür:
  iki lehçe aynı defteri okuyor.

---

## K1-KONTROL — İKİZ SAĞIRLIĞI (H-F3) + HA4 makine mührü  [187d_ikiz_sagirlik.py]

**Makine mührü TUTTU:** ikizin kendi merdiveni τ≤1.00 — τ'≤1.00 katmanları
eklenince m^ya_Hk(1.00) ölçülü m_Hk'yi **bağıl fark ≤ 2e-4** ile geri verdi
(8 bant + HAVUZ; ikiz inşasının 132c-özdeşlik kapanışı; inşa-kesim farkı
L_inşa=12.0304'ün ~7 sınır çizgisi bu düzeyin altında kalır).

**H-F3 ÖLDÜ (kurtarmasız):** ikiz kinematiğinde (1.00,1.10] katman artımı
HAVUZ'da **Δm = −0.0803 ± 0.0046** → |Δm|/2se = 8.6 (eşik ≤1); bant bant
−0.045→−0.106, **0/8 bant geçti**. İkiz o çizgileri hiç "bilmeden"
kurulduğu hâlde katman-sondasına güçlü, İPTAL-YÖNLÜ tepki veriyor; büyüklük
gerçeğinkiyle karşılaştırılabilir (gerçek (1.00,1.10] toplamı −0.078).
**Ön-kayıtlı sonuç uygulanır: 186'nın "iptal gerçeğe özgü" okuması BU
(katman-sondası) biçimiyle düşer.** Ölçülen ayrım şudur (kurtarma değil,
K2'nin bağımsız defteri): ikizin KENDİ pencere-ötesi içeriği azdır ve az
iptal eder (|ζ|_Hk ≈ 0.06-0.09 vs gerçek 0.33) — denizler arası fark
İÇERİK farkıdır; katman-sondası ise kinematikten bağımsız olarak iptal-yönlü
tepki üreten GENEL bir kanalı da ölçer (deniz "sağır" değildir). İki cümle
birlikte 187'nin ana dersidir: **katman defteri yol yapısını verir ama
deniz-ayrıştırıcı değildir; deniz-ayrıştırıcı olan ζ'dir.**

**HA4 (env-ağırlıklı kontrol ızgarası, kendi kinematiği):** erfc-kuyruk
katman artımları HAVUZ'da |Δm| ≤ 4e-5 (env pencere-ötesinde fiilen ölü —
"iptalsiz" varsayımı doğrulandı); makine mührü m^ya,env(1.10) vs ölçülü
m_HA4 bağıl fark ≤ 2e-4 (erfc-inşanın kapanışı).

---

## K3 — YAPI: ι profili, C(τ) rekonstrüksiyonu, biçim KAYITLARI  [187e_yapi.py]

**ι(τ') iptal-yoğunluk profili** (HAVUZ; birincil −Δm/m^kesik; ±jk se):

| katman | (0.86,0.90] | (0.90,0.95] | (0.95,1.00] | (1.00,1.05] | (1.05,1.10] | (1.10,1.15] | (1.15,1.20] |
|---|---|---|---|---|---|---|---|
| ι | .0095±.0010 | .0252±.0011 | .0307±.0064 | .0204±.0068 | **.0774±.0019** | .0258±.0066 | .0305±.0067 |
| ι/Δτ' | 0.24 | 0.50 | 0.61 | 0.41 | **1.55** | 0.52 | 0.61 |

İptal **kesim-hemen-ötesinde yoğun DEĞİL** — profil yayvan (0.4-0.6
taban düzeyi) + (1.05,1.10]'da tabanın ~3 katı bir **TEPE** (bant bant
görünür; ikizin kontrol defterinde de (1.00,1.10] dilimi en güçlü —
tepenin adresi kinematik-genel kanalla örtüşüyor, τ'≈1 "hat-yoğunluğu
frekansı" bölgesi; kayıt). Koherent yan sütun: **her katmanın ζ_i açısı
ayrı ayrı ~180°** (7/7); Σζ_i(≤1.20) = 0.228 ∠180.0° → ζ(HAVUZ)=0.329'un
**%69'u koherent kapsandı** (K1'in %65.8 genlik-kapanışıyla tutarlı).

**C(τ) rekonstrüksiyonu:** C = m^kesik/m_ölç − 1 aynı makinede yeniden
0.479→0.584 (186'nın 1.46→1.58 çarpanının −1'i, hane hane). 1.20-derinlik
rekonstrüksiyonu C_rekon = m^kesik/m^ya(1.20) − 1 = 0.322→0.277:
**tutarlılık mührü 1.20 derinliğinde TUTMADI** (z −18…−52) — bu, ayrı bir
çelişki değil, H-F1'in ölçtüğü %34 kalan kuyruğun aynı cümlesidir: C'nin
tamamı (0.86,1.20]'den ödenmiyor. C'nin τ-büyümesi (0.46→0.58) kesilmiş
defterde ancak kısmen (0.32→0.28, ters-yönlü) yeniden kuruluyor: hedef
bağımlılığın taşıyıcısı ağırlıkla τ'>1.20 kuyruğudur (KAYIT).

**Biçim KAYITLARI (türetim DEĞİL, kayıt):** C(τ) = k·τ^β için en iyi
(k, β) = (0.526±0.015, 0.196±0.048), fit-χ²/dof = 14.3 — **kötü** (C
profili U-biçimli; tek güç yasasına uymuyor). ι(τ') için üstel
A·e^{−(τ'−0.86)/λ}: (A, λ) = (0.025, 2.00), χ²/dof = 218; güç
A·(τ'−0.86)^{−p}: p = 0.05, χ²/dof = 222 — **ikisi de kötü**: profil ne
üstel ne güç; "yayvan taban + (1.05,1.10] tepesi" biçimsiz kaydıdır.
Açık cümle: **bu bölümde hiçbir şey türetilmedi; ölçülen profiller ve
başarısız biçim denemeleri kaydedildi.**

---

## K4 — HA4 %26 KAPANIŞI (H-F4)  [187f_HA4.py]

Birincil defter aritmetiği: M^pred = m^env(τ_c=0.86; HA4-kin)/m_Hk(ölçülü):

| bant τ̄ | M_HA4 ölçülü (186e AYNEN) | M^pred (iptalsiz, 0.86) | z | M^pred (tam-env, 1.10) |
|---|---|---|---|---|
| 0.475 | 1.2640±.0063 | 1.2646±.0080 | +0.1 | 1.2640 |
| 0.526 | 1.2547±.0068 | 1.2552±.0083 | +0.1 | 1.2546 |
| 0.576 | 1.2348±.0085 | 1.2351±.0082 | +0.0 | 1.2346 |
| 0.626 | 1.1930±.0104 | 1.1931±.0065 | +0.0 | 1.1928 |
| 0.676 | 1.1299±.0099 | 1.1297±.0061 | −0.0 | 1.1296 |
| 0.725 | 1.0698±.0084 | 1.0695±.0049 | −0.0 | 1.0696 |
| 0.775 | 1.0322±.0146 | 1.0319±.0051 | −0.0 | 1.0321 |
| 0.830 | 1.0190±.0229 | 1.0186±.0076 | −0.0 | 1.0189 |

**bant-χ²/dof = 0.00 ≤ 2 → H-F4 MÜHÜR: 186-K4 çatlağı KAPANDI.**
%26'nın defteri: HA4'ün merdiveni erfc'yle pencere içinde ölür → iptal
edecek pencere-ötesi içeriği YOKTUR → karışımı "iptalsiz pencere-içi erfc
karışımı"dır; Hkeskin'in paydası ise kendi kuyruğuyla iptallidir. Oran
1.264→1.019'u parametresiz verir. **Şerh (türetilen vs ölçülen):** HA4'ün
ds'i inşa gereği erfc-merdivenin kendisi olduğundan bu kapanış
denizler-arası bir öngörü değil, DEFTER-ARİTMETİĞİ kapanışıdır (doğru
kinematikte muhasebe). 186e'nin ölümünün (χ²/dof=350) sebebi burada
görünür kılındı: aynı erfc-merdiven YABANCI kinematikte (Hk/gerçek)
koşulmuştu; %26 uyuşmazlığın tamamı kinematik seçimindeydi, mekanizmada
değil.

---

## İKİNCİL — DERİN SONDA (1.20,1.30]; 100k alt-örneklem  [187g_derin_sonda.py]

*(İKİNCİL damgalı; ayrı se; ana hükümlere girmez.)* (1.20,1.30] katmanı
(286 009 çizgi, q≤6.2e6) her-3.-nokta alt-örnekleminde (100 000 nokta;
m^ya(1.20) ve (1.30) AYNI alt-örneklemde izdüşürüldü): HAVUZ
**Δm = −0.0238 ± 0.0017** (bantlar −0.012…−0.030, hepsi negatif) — iptal
τ'>1.20'de SÜRÜYOR; 0.05-genişlik başına hız (−0.012) 1.10-1.20
katmanlarının (−0.021…−0.024) altında: kuyruk yavaşlayarak ödüyor. Şerh:
alt-örneklemin |·|-defteri seviye olarak yukarı-yanlıdır (m^ya_alt(1.20)
HAVUZ = 0.687 vs tam-örneklem 0.620 — daha az noktayla izdüşüm gürültü
tabanı büyür); İKİNCİL'in anlamlı niceliği seviye değil, aynı-örneklem
ARTIMIDIR.

---

## HÜKÜM (eşikler K0'da donmuş; kurtarma yok)

| hipotez/kapı | hüküm | dayanak |
|---|---|---|
| **H-F1** (yakınsama; monoton + ≥%50 @1.20) | **MÜHÜR** | monoton EVET; kapanan %65.8±0.3 [187b] |
| — çatal (a) SIĞ İPTAL | ateşlemedi | art arda iki \|Δm\|<2se yok; artımlar sürüyor [187b] |
| — çatal (b) DERİN KUYRUK | ateşlemedi | kalan açık %34.2 ≤ %50 [187b] |
| **H-F2** (faz-uyum 180°±15°) | **MÜHÜR — SAF YIKICI** | 8/8 bant; açı 180°±0.7° (pencereden 20 kat dar) [187c] |
| **H-F3** (ikiz sağırlığı) | **ÖLDÜ** | Δm(1.00,1.10] = −0.080±0.005, 8.6σ; 0/8 bant; katman-sondası kinematik-genel kanal ölçüyor [187d] |
| — ikiz makine mührü (τ'≤1.00 → m_Hk) | **TUTTU** | bağıl fark ≤2e-4, 8/8 [187d] |
| **H-F4** (HA4 %26 kapanışı) | **MÜHÜR — KAPANDI** | χ²/dof = 0.00 ≤ 2; 1.264→1.019 hane hane [187f] |
| **K3** (yapı) | **KAYIT** | yayvan iptal + (1.05,1.10] tepesi; C_rekon@1.20 mührü tutmadı (=%34 kuyruk); biçimler (güç/üstel) uymadı [187e] |
| makine mühürleri (4 adet) | **TUTTU** | ≤5.3e-16 / 1.1e-16 / 2.2e-16 / 2e-4 [187b,187d,187f] |

---

## MANŞET (aday cümle)

> **İptalin açısı ölçüldü: gerçek deniz, pencere-içi merdiven karışımının
> %33'ünü pencere-ötesi içerikle tam 180°'de (±0.7°; H-F2 penceresinden
> 20 kat dar) yok ediyor — "faz-uyumlu iptal" artık sıfat değil, sayı:
> ζ = 0.329∠180.0° (ikizde 0.077∠179.8° — içerik farkı 4 kat).** Katman
> defteri iptalin üçte ikisini (0.86,1.20] diliminde buldu (kapanan pay
> %65.8±0.3, monoton; H-F1 mühür; iki çatal okuma da ateşlemedi), profil
> kesim-dibinde yoğun değil YAYVAN, (1.05,1.10]'da tabanın ~3 katı bir
> tepeyle; kalan üçte bir τ'>1.20 kuyruğunda. Aynı defter HA4'ün %26
> bilmecesini parametresiz kapattı (H-F4: χ²/dof=0.00): HA4'ün karışımı
> iptalsiz pencere-içi erfc karışımıdır — 186'daki 350'lik χ²'nin tamamı
> yanlış kinematik seçimiydi. Bir ölüm dürüst kaydedildi: H-F3 ikiz-sağırlık
> kontrolü ÖLDÜ — ikiz kinematiği de derin katman-sondasına iptal-yönlü
> tepki veriyor (−0.080±0.005); katman-sondası deniz-ayrıştırıcı değildir,
> deniz-ayrıştırıcı olan ζ'dir. α=1.30 sorusunun adresi daralıyor: iptal
> payı |ζ(τ)| = 0.318→0.381 kendi τ-profiliyle yeni İLKEL nesnedir.

Teslim: bu rapor + `187_configs/` (187a-187h) + `187_faz_iptal.png`
(sol: yakınsama defteri + ι profili; sağ: ζ açı defteri + M_HA4
öngörü-vs-ölçüm) + `scratchpad/187/` (ONKAYIT_187.json,
katman_*_{gercek,Hkeskin,HA4env}.npy, katmanproj_*.npz,
K1_katman_defteri.json, K2_zeta.json, K1k_ikiz_sagirlik.json,
K3_yapi.json, K4_HA4.json, IKINCIL_derin_sonda.json, 187b-187g logları).
