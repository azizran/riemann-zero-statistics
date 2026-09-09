# 186 — M(τ) İLETİMİNİN TÜRETİMİ (zarf kendi karışımını doğurur mu?)

**Soru (KALEM_M_ILETIMI_09EYL2026):** 185 zarfın taşıyıcısını buldu:
karışım-iletim profili M(τ) = m_g/m_Hk = 0.796±.004 → 0.721±.021
(m = w_ölçülü − w^öz). "Neden α=1.30?" sorusu M'nin türetimine indirgendi.
Üç güzergâh: **G1** öz-tutarlılık (yeniden-ağırlıklı merdiven — zarf kendi
karışımını doğurur mu?), **G2** v-kanalının yüksek-τ İLK ölçümü + Not-2/3
kısıt köprüsü, **G3** biçim defteri (kısıt, türetim değil).

Her sayının yanında onu üreten betik adı vardır. Tek dalga; ölümler
kurtarmasız; git'e dokunulmadı; sonuç KAPTAN+KULLANICI ortak teftişine
sunulur, commit sonra.

---

## Ölçüm tanımı (K0'da donmuş; 184/185 makinesi AYNEN)

m/M tanımları 185'ten AYNEN (w-lehçesi, toplam-oranı; ölçülü â'lar 184 K1
npz'lerinden, â^öz'ler 185 OZ önbelleklerinden). **G1:** yeniden-ağırlıklı
merdiven `ds^ya_n = Σ_q' 2a_q'·ρ(τ_q')·sin(ω_q' g_n/2)·cos(ω_q' m_n)`,
ρ(τ')=1−0.149τ'^1.30 (donmuş), GERÇEK kinematikle (g kesin inversiyonla;
evren = 184'ün 3425 çizgisi, τ≤0.86); m^pred = izdüşüm − ρ_q·c_q^öz (kesin);
M^pred = m^pred/m_Hk. Donmuş yan-lehçe: aynı koşular Hkeskin kinematiğinde
(ikiz-çekirdek okuması). **G2:** v̂_q = |Ĝ_q|/(ḡ·ae_q), Ĝ_q=⟨(g−ḡ)e^{−iω_q m}⟩,
iki denizde özdeş; birincil lehçe ORAN V=v_g/v_Hk; köprü
M^pred = [(1.017−0.884v_g)/(1.017−0.884v_Hk)]². Birincil σ(M) = defter
konvansiyonu (bağımsız yayılım); ölüm eşikleri: G1/G2 ve K3-zincir
bant-χ²/dof ≤ 2. Makine mührü: ρ≡1 kontrolü 185-B ters-yönünü (r_pred>1)
yeniden vermeli.

---

## K0 — DONMUŞ ÖN-KAYIT  [186a_onkayit.py]

`186/ONKAYIT_186.json` yazıldı: **sha256 = 23a8130a…**, damga
**Wed Sep 9 16:15:54 +03 2026** (tüm ölçümlerden önce). Donan: m/M
tanımları; G1 kestirimcisi + ρ≡1 makine-mührü + yan-lehçe; v̂
normalizasyonu; G3 üç biçim; bant ızgarası `[.45 … .86]` + kuyruk τ>0.70
(184-185 AYNEN); 8-blok loo-jackknife; ölüm eşikleri; K3 zincir + fit
makinesi (185d AYNEN); K4'ün kazanan-bağımlı iki varyantı.

**Ölçülü-M yeniden üretimi (makine kontrolü):** 186b'nin defteri 185f'yi
hane hane verdi (0.7959±.0044 → 0.7206±.0212). [186b_g1_oztutarlilik.py]

---

## K1 — G1 ÖZ-TUTARLILIK SINAVI  [186b_g1_oztutarlilik.py]

**MAKİNE MÜHRÜ TUTTU (iki lehçede de):** ρ≡1 kontrol zinciri 8 bandın
hepsinde r_pred>1 — birincil (gerçek-kin.) 1.075→1.140, yan (Hk-kin.)
1.048→1.080; 185-B'nin ters-yön ölümü yeniden üretildi. Cebirsel köprü:
w-lehçesinde M^pred≡1 zinciri 185-B defterini hane hane verir
(1.0145 / 1.0242 / 1.0344 / 1.0568 ≈ 185-B'nin 1.0145 / 1.0242 / 1.0344 /
1.0565'i) → ölümler güvenilir makinede ölçülmüştür.

**G1 bant defteri** (birincil: gerçek kinematik, ρ=184 yasası; ±jk se;
σ = defter-seM):

| τ̄ | M ölçülü | M^pred (G1) | z | M^pred (ρ≡1 kontrol) | M^pred (yan, Hk-çekirdek) |
|---|---|---|---|---|---|
| 0.475 | 0.7959±.0044 | 1.1283±.0073 | +76.1 | 1.1771 | 1.0521 |
| 0.526 | 0.7893±.0044 | 1.1117±.0091 | +73.6 | 1.1598 | 1.0484 |
| 0.576 | 0.7816±.0072 | 1.0956±.0096 | +43.6 | 1.1429 | 1.0402 |
| 0.626 | 0.7761±.0077 | 1.0862±.0070 | +40.3 | 1.1333 | 1.0383 |
| 0.676 | 0.7688±.0079 | 1.0742±.0075 | +38.6 | 1.1219 | 1.0297 |
| 0.725 | 0.7582±.0078 | 1.0705±.0050 | +39.9 | 1.1198 | 1.0212 |
| 0.775 | 0.7466±.0134 | 1.0733±.0091 | +24.3 | 1.1268 | 1.0108 |
| 0.830 | 0.7206±.0212 | 1.0793±.0091 | +16.9 | 1.1414 | 0.9754 |

**χ²/dof: birincil = 2336, ρ≡1 = 3095, yan (Hk-çekirdek, ρ=yasa) = 1517**
(eşik 2) → **G1 ÖLDÜ, kurtarmasız.** Zarfla yeniden-ağırlıklama, karışımı
yalnız ~4-6 puan bastırıyor (1.18→1.13); ölçülü bastırma 20-28 puan.
"Zarf kendi karışımını doğurur" sabit-nokta hipotezi yanlış: karışım
bastırması, çizgi-genliklerinin zarf-katlanmasından ÇOK daha güçlü.

**BONUS TEŞHİS — pencere-dışı YIKICI girişim (yeni bulgu):** kesik
(pencere-içi, τ'≤0.86) nominal merdivenin karışımı, gerçek denizde ölçülü
karışımın **1.46–1.58 KATI** (bantlara göre m^ya(ρ≡1)/m_ölç = 1.479→1.584);
ikizde yalnız **1.04–1.10 katı**. Yani gerçek deniz, pencere-dışı
(τ'>0.86) içerikle karışımın ~üçte birini FAZ-UYUMLU biçimde İPTAL ediyor;
ikizde bu iptal neredeyse yok. M(τ)<1'in mekanizma adresi genlik-katlanması
değil, **çizgiler-arası faz-örgüsünün pencere-ötesi iptali**. (Sızıntı
defteri: korr(ds^ya,ds)=0.953–0.957 gerçek, 0.964 ikiz; artık-varyans
payı 0.10–0.13.) [186b]

---

## K2 — G2: v-KANALININ YÜKSEK-τ İLK ÖLÇÜMÜ + KISIT KÖPRÜSÜ  [186c_g2_vkanal.py]

**v(τ) profilleri (bağımsız bulgu; τ∈[0.45,0.86]'da İLK ölçüm; ±jk se):**

| τ̄ | v_g | v_Hk | v_HA4 | V = v_g/v_Hk |
|---|---|---|---|---|
| 0.475 | 0.5849±.0011 | 0.6193±.0004 | 0.6683±.0012 | 0.9445±.0019 |
| 0.526 | 0.5973±.0008 | 0.6388±.0008 | 0.6841±.0014 | 0.9351±.0013 |
| 0.576 | 0.6085±.0015 | 0.6577±.0014 | 0.6794±.0015 | 0.9253±.0018 |
| 0.626 | 0.6206±.0012 | 0.6769±.0020 | 0.6446±.0018 | 0.9169±.0013 |
| 0.676 | 0.6285±.0012 | 0.6910±.0018 | 0.5787±.0017 | 0.9096±.0013 |
| 0.725 | 0.6344±.0025 | 0.7032±.0021 | 0.5125±.0019 | 0.9021±.0022 |
| 0.775 | 0.6352±.0021 | 0.7075±.0027 | 0.4636±.0021 | 0.8978±.0018 |
| 0.830 | 0.6280±.0039 | 0.7032±.0053 | 0.4358±.0034 | 0.8931±.0016 |

- v_g yüksek-τ'da **~0.628–0.635'te doyuyor**; en düşük bantta 0.5849 —
  Not 2'nin ~0.55-0.60 doyum aralığının içinde (konvansiyon şerhiyle).
  v_Hk 0.703'e tırmanıyor; **v_HA4, erfc kesimini birebir izleyerek
  0.668→0.436'ya çöküyor** (v-kanalında erfc-çöküşünün ilk görüntüsü).
- **YAPISAL BULGU — V ≡ r dejenerasyonu:** ORAN lehçesinde V=v_g/v_Hk,
  184'ün r(τ) defteriyle 4 haneye kadar ÖZDEŞ (0.9445→0.8931). Sebep
  cebirsel: g_n, ds'in KESİN inversiyonu olduğundan |Ĝ|-tabanlı (modül)
  her aralık-serisi kestirimcisi w-kanalıyla dejeneredir — bu lehçede
  yüksek-τ'da BAĞIMSIZ bir v-kanalı yok; iki kanal teke katlanıyor.
  (Bağımsız v için faz/kuadratur istatistiği gerekir — 187+ adayı; bu
  ölümün kapsam şerhi.)
- **Kısıt köprüsü:** M^pred = 1.134→1.365, ölçülü M'nin ters yönünde;
  **χ²/dof = 4187 → G2 ÖLDÜ, kurtarmasız.** Yan-okuma köprü-vs-r de ölü
  (χ²/dof=8058). Kısıt-seviye kontrolü: √w = 1.08–1.19'a karşı
  1.017−0.884v̂ = 0.39–0.50 (fark rms 0.64/0.74 iki denizde) — Not-2/3
  kısıtı **bu konvansiyonda** yüksek-τ'ya UZANMIYOR (uzanma hipotezi
  sınandı ve öldü; ölümün bir payı v̂-konvansiyonunda olabilir, şerh).

---

## K3 — KAPANIŞ: zincir + türetilmiş (c,α) + sabit-nokta  [186d_kapanis.py]

Kazanan kuralı (K0): eşiği geçen yok → zincir, χ²'si küçük olan G1'e
YALNIZ TEŞHİS olarak takıldı (hüküm ÖLDÜ; kurtarma yok).

| τ̄ | r ölçülü (184) | r_pred (zincir, G1) | z |
|---|---|---|---|
| 0.475 | 0.9445±.0017 | 1.0585±.0022 | +68.4 |
| 0.626 | 0.9169±.0040 | 1.0655±.0028 | +37.6 |
| 0.725 | 0.9021±.0044 | 1.0730±.0023 | +39.3 |
| 0.830 | 0.8931±.0135 | 1.1033±.0048 | +15.6 |

**Zincir bant-χ²/dof = 1987 → KAPANMADI.** Güç-yasası fitleri (dof=6):

| eğri | c | α | fit-χ²/dof |
|---|---|---|---|
| ölçülü (kontrol) | **+0.1492** | **1.304** | 0.34 — 184'ün (0.149, 1.30, 0.34)'ü AYNEN |
| zincir r_pred (G1) | −0.0884 | 0.571 | 0.95 |
| yan (Hk-çekirdek) | −0.0599 | 0.800 | 0.16 |

Türetilmiş (c,α) hedefin işaretçe karşısında — 185-H-W2'nin ölümü, karışım
kanalının içinden bir kez daha doğrulandı.

**SABİT-NOKTA KALEM-SAYISALI** (ρ_c=1−c·τ^1.30 ailesi; ρ=yasa ve ρ≡1
koşularından doğrusallaştırma; yeni koşu yok): harita güçlü kontraktif
(|dr_pred/dρ| = 0.275–0.311 ≪ 1 → yineleme yakınsar) ama sabit nokta
HER bantta **r\* = 1.107→1.203, c\* = −0.281…−0.209** — yani öz-tutarlı
çözüm bir ANTİ-ZARF (genlik büyümesi), ölçülü zarfın (r<1) tarafında bile
değil. **Güç-yasası zarfı, karışım-katlanması dinamiğinden kendini
ÖZ-ÜRETMİYOR; α=1.30 bu denklemin çözümü değil.** [186d]

---

## K4 (bonus, örneklem-dışı) — HA4: karışım-eşitsizliği  [186e_HA4.py]

Ölçülü defter (m_HA4 = w_HA4 − w^öz,erfc; M_HA4 = m_HA4/m_Hk ± defter-σ):
**M_HA4 = 1.2640±.0068 → 1.0190±.0207** — 185'in "~%8-10 eşitsizliği"
(w-düzeyinde r_HA4≈1.079) karışım-lehçesinde ilk kez bant bant ölçüldü:
orta-τ'da **%26'ya varıyor**, kuyrukta %2'ye iniyor.

G1-env teşhisi (birincil, Hk-kinematik): M^pred = 1.072→0.829,
**χ²/dof = 350 → ÇATLAK** (dürüst). AMA eşitsizliğin YÖNÜ (M_HA4>1, ilk
4 bant) parametresiz yakalandı: erfc ile kuyruğu ölmüş çizgi kümesinin
karışımı, tam merdiveninkinden BÜYÜK — pencere-dışı yıkıcı iptalin (K1
bonus teşhisi) bağımsız teyidi: kuyruk çizgileri karışımı İPTAL eden
taraftadır, kaldırınca karışım artar. Yan sütun (gerçek-kinematik):
χ²/dof = 127; kuyruk bantlarında z = −1.1 / −2.5 — kuyruk-ucu neredeyse
yakalanıyor. G2-köprü yan sütunu: χ²/dof = 5739 (ölü). Büyüklük (%26)
hiçbir varyantla açıklanamadı → 185-K4 çatlağı AÇIK kalır.

---

## G3 — BİÇİM DEFTERİ (kısıt; türetim değil)  [186f_bicim.py]

| biçim | parametreler (±jk se) | fit-χ²/dof (dof=6) |
|---|---|---|
| doğrusal M=a+b·τ | a=0.8722±.0089, b=−0.1580±.0107 | **0.23** |
| güç M=1−k·τ^β | k=0.2774±.0034, **β=0.424±.034** | 0.32 |
| v-afin M=a+b·V | a=−0.0871±.0417, b=+0.9368±.0471 | 0.49 |

Üçü de yaşıyor (pencere dar, biçimler yakın-dejenere; doğrusal/güç ayırt
edilemez). Kayıt: 1−M ≈ 0.277·τ^0.42 — KALEM'in "√τ kokusu" ön-bakışına
yakın ama daha sığ (β=0.42±0.03; √τ=0.5 buradan 2.2σ). M'nin (k,β) =
(0.277, 0.42) çifti, zarfın (0.149, 1.30)'undan hem ölçekçe hem üsçe AYRI:
**M(τ) kendi üssünü taşıyan ayrı bir ilkel nesnedir.** v-afin bilgi notu:
V≡r olduğundan bu satır "M, r'de afin" demektir (M ≈ 0.94·r − 0.09).

---

## HÜKÜM (ölüm eşikleri K0'da donmuş; kurtarma yok)

| hipotez/kapı | hüküm | dayanak |
|---|---|---|
| **G1** (öz-tutarlılık: zarf karışımını doğurur) | **ÖLDÜ** | χ²/dof=2336 ≫ 2; öngörü zarfın ters yönünde [186b] |
| — makine mührü (ρ≡1 → 185-B ters-yönü) | **TUTTU** | r_pred(ρ≡1)>1 tüm bantlar, iki lehçede; M≡1 zinciri 185-B'yi hane hane verir [186b] |
| **G2** (v-köprüsü; kısıtın yüksek-τ uzantısı) | **ÖLDÜ** | χ²/dof=4187; kısıt seviyede de uzanmıyor (rms 0.64/0.74); V≡r dejenerasyon şerhi [186c] |
| **K3 zinciri** | **KAPANMADI** | χ²/dof=1987; türetilmiş (c,α)=(−0.088, 0.57) hedefin işaretçe karşısında [186d] |
| — sabit-nokta analizi | **ANTİ-ZARF** | r\*=1.11→1.20; harita kontraktif ama yanlış merkez; α=1.30 öz-üretilmiyor [186d] |
| K4 (HA4, bonus) | **ÇATLAK** (yön yakalandı) | χ²/dof=350; M_HA4>1 yönü parametresiz doğru; büyüklük (%26) açık [186e] |
| G3 (biçim) | **KAYIT** | 1−M=0.277·τ^0.42 (β=0.42±0.03); doğrusalla dejenere [186f] |

---

## MANŞET (aday cümle)

> **Zarf kendi karışımını DOĞURMUYOR: M(τ) karışım-iletimi, zarfın
> parametresiz katlanmasından da (G1: χ²/dof=2336), Not-2/3 kısıt
> köprüsünden de (G2: χ²/dof=4187) türetilemedi — makine mührü iki ölümde
> de tuttu (ρ≡1, 185-B ters-yönünü hane hane yeniden üretti).**
> Öz-tutarlılık haritasının sabit noktası bir ANTİ-ZARF (r\*≈1.11→1.20):
> α=1.30, karışım-katlanması denkleminin çözümü değil. İki yeni adres
> bulundu: (1) gerçek deniz, pencere-içi merdiven karışımının ~üçte birini
> pencere-ötesi (τ>0.86) çizgilerle FAZ-UYUMLU biçimde iptal ediyor
> (ikizde bu iptal yok denecek kadar az; HA4'ün M>1 eşitsizliğinin yönü de
> aynı mekanizmayla parametresiz yakalandı) — bastırmanın mekanizması
> genlik değil, çizgiler-arası faz-örgüsü; (2) v-kanalının yüksek-τ İLK
> ölçümünde ORAN lehçesi V=v_g/v_Hk, r(τ) ile 4 haneye kadar özdeş çıktı:
> aralık serisinin modül-kestirimcisi w-kanalıyla cebirsel dejeneredir —
> bu pencerede "iki kanal" tektir. M(τ) = 1 − 0.277·τ^0.42 kendi (sığ)
> üssünü taşıyan, zarfın (0.149, 1.30)'undan ayrı bir İLKEL nesne olarak
> kalıyor; "neden 1.30?" sorusunun adresi artık pencere-ötesi faz-iptal
> defteridir.

Teslim: bu rapor + `186_configs/` (186a–186g) + `186_M_turetim.png`
(sol: ölçülü M vs G1/G2 öngörüleri; sağ: zincir r_pred + sabit-nokta) +
`scratchpad/186/` (ONKAYIT_186.json, G1_proj_{gercek,Hkeskin}.npz,
G1_sonuc.json, G2_sonuc.json, K3_sonuc.json, K4_HA4.json, G3_bicim.json,
186b–186f logları).
