# KALEM — 188: ζ(τ)'NİN KİMLİĞİ — İPTAL ÇEKİRDEĞİNİN HARİTASI (tarak köprüsü)
(23 Eylül 2026 — 187'nin açık kalemi; iki haftalık aradan sonra kaldığımız yer)

## Durum

187: gerçek deniz pencere-içi (τ≤0.86) karışımın ~1/3'ünü pencere-ötesi
çizgilerle 180°'de iptal ediyor: ζ = 0.329∠180.0°±0.07 (τ_c=0.86); bant
profili |ζ_g| = 0.328, 0.323, 0.318, 0.318, 0.318, 0.327, 0.346, 0.381
(ortada taban, kesime doğru yükseliş); ikizde 0.089→0.059. Katman defteri:
iptal yayvan + (1.05,1.10]'da ~3× tepe; ikiz kinematiği de derin sondaya
iptal-yönlü tepki veriyor (H-F3 öldü: kinematik-genel kanal). Dış girdi
(bağımsız paket T2/T3, örtüşen pencereler — yalnız NİTEL): minimum τ≈0.6-0.65,
yükseliş kesimi izliyor. Soru: iptali HANGİ çekirdek taşıyor; ζ neden ~1/3,
neden kesime doğru yükseliyor?

## Kalem cebiri

(1) Teleskop özdeşliği: 2 sin(ω g_n/2) cos(ω m_n) = sin(ω z_{n+1}) − sin(ω z_n).
Her pencere-ötesi çizgi q', f(t)=a_{q'} sin(ω_{q'} t)'nin sıfırlardaki ayrık
türevi olarak girer.

(2) Çekirdek ds'i taşır: g_n = ḡ_n(1+ds_n) ⇒ sin(ω' g_n/2) = sin(πτ'(1+ds_n)).
ds_n pencere çizgisi q'yu içerdiğinden, q''nün q-izdüşümüne katkısı birinci
mertebede ∝ a_q × πτ' cos(πτ') × a_{q'} × Ĝ_mid(ω_{q'}) — yani q'NUN KENDİ
katsayısına ORANTILI (açı 0° ya da 180°; τ'≈1'de cos(πτ')≈−1 → YIKICI), ve
ağırlığı orta-nokta örgüsünün YAPI ÇARPANI Ĝ_mid(ω') = ⟨e^{iω' m_n}⟩ (mutlak
m_n). Ĝ_mid büyük olduğu yerler: Bragg tarağı ω'≈L_yerel (τ'≈1.00, pencere
içi L kayması ile ±0.006 yayvan) ve asal UYDULARI ω' = L ± log p:
  τ' = 1 ± τ_p:  1.0576 / 0.9424 (p=2), 1.0913 / 0.9087 (p=3),
                 1.1338 / 0.8662 (p=5), 1.1618 / 0.8382 (p=7).
Bu resim H-F3'ün ölümünü açıklar (her sadakatli örgünün tarağı var → sonda
her kinematikte iptal-yönlü) ve ζ farkını İÇERİĞE bağlar (gerçek deniz bu
çizgileri taşır, ikiz τ'>1.00'da taşımaz). (Taylor burada YALNIZ yorumdur —
πτσ_ε~1 dersi; bütün sayılar kesin beklenen-değerle.)

(3) VERİ-ÖNCESİ İPUCU (187 defterinden, ölçüm değil): jackknife se'si
alışılmadık küçük katmanlar (0.86,0.90] ±.0008, (0.90,0.95] ±.0009,
(1.05,1.10] ±.0015 — diğerleri ±.005 — tam da 1−τ_5, 1−τ_3/1−τ_2 ve
1+τ_2/1+τ_3 uydularını içeren katmanlar. AMA iz KISMİ: 1+τ_5 (1.10,1.15]
ve 1+τ_7 (1.15,1.20] katmanlarının se'si büyük (±.0052/±.0053). Yani ipucu
yalnız 2 ve 3'ün uydularında tutuyor; sınanacak, kanıt sayılmaz.

## Ölçülecek nesne — ÇEKİRDEK HARİTASI K(τ_b, τ')

Pencere-ötesi dilimler τ' ∈ (0.86, 1.30], genişlik 0.005 (88 dilim; ~420k
asal-kuvvet). İnce bantlar τ_b genişlik 0.01 (41 bant) + 184'ün 8 bandı.
K(b,s) = Σ_{q∈b} c_q^{(s)}·conj(karışım_q) / Σ_{q∈b} |karışım_q|²  (KARMAŞIK;
187c'nin ζ tanımının dilim-dilim açılımı; Σ_s K(b,s) = ζ'nin (0.86,1.30]
payı). Ayrıca her dilim için S(s) = Σ_{q'∈s} a_{q'}·|Ĝ_mid(ω_{q'})| ve
W(s) = Σ_{q'∈s} a_{q'} (yapısız ağırlık).

## Hipotezler (ön-kayıtla donar; ölümler kurtarmasız)

- **H-188a — AYRIK GERİ-BESLEME (tarak köprüsü):** K ≈ u(τ_b)·v(τ') (rang-1);
  |v(τ')| S(τ')'yi izler; tepeler Bragg + asal uyduları konumunda.
  Kapılar: (i) karmaşık SVD rang-1 varyans payı ≥ 0.80; (ii) corr(|v|, S) ≥
  0.80 VE corr(|v|,S) − corr(|v|,W) ≥ 0.15; (iii) 187'nin (1.05,1.10] tepesinin
  ≥ %60'ı 1+τ_2, 1+τ_3 uydu pencerelerinde (±0.006). Ölüm: (i) < 0.60 veya
  (ii)'nin ilk koşulu < 0.50.
- **H-188b — KÖŞEGEN / KESİME-YAKINLIK:** rang-1 artığı köşegen sırtlar
  taşır (τ' = 2 − τ_b katlanması ve/veya τ' − τ_b küçük kenar-sızıntısı) ve
  |ζ|'nin kesime doğru yükselişinin ≥ %60'ını taşır. Ölçüt: τ_b ≥ 0.70
  bantlarında köşegen şeridi (±0.02) kontrastı ≥ 3× köşegen-dışı ortalama.
- **H-188c — YAYVAN TABAN (sıfır-hipotezi):** yapı çarpanı hiçbir şey
  eklemiyor: corr(|v|,W) ≥ corr(|v|,S). Tutarsa tarak köprüsü ölür.
- **H-188d — KİNEMATİK-GENEL / İÇERİK AYRIMI (H-F3'ün açıklaması):** aynı
  harita İKİZ kinematiğinde (τ' ≤ 1.20, dilim 0.01) şekilce gerçeğinkine
  benzer: dilim-profilleri corr ≥ 0.80. Tutarsa: sonda kinematik-genel,
  ζ farkı içerik — 187'nin iki cümlesi tek mekanizmada birleşir.
- **K3 — YENİDEN KURULUM:** u(τ_b) ve rang-1 artığından 8-bant |ζ^{≤1.30}(τ_b)|
  profilinin şekli (orta taban + kesime yükseliş) yeniden kurulur; hangi
  bileşenin hangi özelliği taşıdığı bant-bant raporlanır. "Türetilen vs
  ölçülen" ayrımı açık cümleyle. α=1.30 bağı: ζ-makası → M → r zinciri
  (187'nin kapanış aritmetiği) kazanan çekirdekle ileri-hesaplanabiliyorsa
  kayıt; değilse dürüstçe "açık".

## Kapılar

- **K0 — YENİDEN KURULUM KAPISI (veri yalnız Mac mini'de kaldı):** gerçek deniz
  `128_odl_zeros6_2e6_zeros.npz` son-300k'dan; Hkeskin `164_configs/164_insa.py`
  ile (maks|F|≤1e-8, sıralılık TAM, L=12.02959324). Zincir 184b→185b→186b
  (ρ≡1)→187c; GEÇİŞ ŞARTI: 184 w-bantları (gerçek 1.170…1.256; Hkeskin
  1.239…1.406) ±se içinde VE ζ_g = 0.3287±0.0023 ∠179.96°, ζ_Hk = 0.0768±0.0064.
  Tutmazsa DUR, raporla (yeni ölçüm yok).
- **K0b — ÖN-KAYIT (sha+damga, K0'dan sonra, haritadan ÖNCE):** dilim/bant
  ızgaraları, K/S/W tanımları, uydu konumları, H-188a..d eşikleri, K3 kuralları.
- **K1 — HARİTA (gerçek):** 88 dilim; dilim serileri alt-örneklemde (her 3.
  nokta, 100k — K doğrusal, alt-örneklem yansız; 187g şerhi: modül değil
  karmaşık toplam kullan); iki dilimde tam-örneklem kontrolü. Dilim başına
  kontrol noktası.
- **K2 — İKİZ HARİTASI** (H-188d) + Ĝ_mid ölçümü (her pencere-ötesi çizgide).
- **K3 — ANALİZ:** SVD, korelasyonlar, köşegen kontrastı, yeniden kurulum.

Ölümler kurtarmasız; TEK DALGA; koşular nohup + ≤30 sn yoklama (5 dk ön-plan
YASAK); ajan git'e DOKUNMAZ; sonuç ORTAK TEFTİŞE gelir, commit sonra.
Bağımsız paketten alınan dersler: Ĝ'de MUTLAK t (ortalama çıkarma yok);
açı daima modülle birlikte; ζ alıntısında τ_c.
