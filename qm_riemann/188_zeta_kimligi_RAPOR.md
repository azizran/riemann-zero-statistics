# 188 — ζ(τ)'NİN KİMLİĞİ: İPTAL ÇEKİRDEĞİNİN HARİTASI (tarak köprüsü)

**Soru (KALEM_ZETA_KIMLIGI_23EYL2026):** 187, gerçek denizin pencere-içi
(τ≤0.86) merdiven karışımının ~1/3'ünü pencere-ötesi çizgilerle tam 180°'de
iptal ettiğini ölçtü: ζ = 0.329∠180.0° (τ_c=0.86), bant profili orta taban +
kesime yükseliş, katman defterinde (1.05,1.10]'da ~3× tepe; ikiz kinematiği de
derin sondaya iptal-yönlü tepki veriyor (H-F3 öldü). Bu kalem iptali **hangi
çekirdeğin** taşıdığını K(τ_b, τ') haritasıyla sorar: ayrık tarak geri-beslemesi
mi (Bragg + asal uyduları; H-188a), köşegen/kesime-yakınlık mı (H-188b), yayvan
taban mı (H-188c, sıfır), ve harita kinematik-genel mi (H-188d)?

Her sayının yanında onu üreten betik adı vardır. Tek dalga; ölümler
kurtarmasız; git'e dokunulmadı; sonuç ORTAK TEFTİŞE sunulur, commit sonra.

---

## Ölçüm tanımı (K0b'de donmuş; 187b/187c makinesi AYNEN)

Dilim s = τ' ∈ (e_s, e_{s+1}], 88 dilim × 0.005, (0.86, 1.30]; dilimdeki TÜM
asal-kuvvetler, a_q = Λ(q)/(π√q·log q) (187b.asal_kuvvetler). Dilim serisi
`dds_s(n) = Σ_{q'∈s} 2a' sin(ω' g_n/2) cos(ω' m_n)` (187b.katman_serisi, ρ≡1)
ALT-ÖRNEKLEMDE (her 3. nokta, 100 000); izdüşüm c^{(s)}_q = 2⟨dds_s e^{−iω_q m}⟩.
**K(b,s) = Σ_{q∈b} c^{(s)}_q·conj(karışım_q) / Σ_{q∈b}|karışım_q|²**, karışım =
c^kesik − c^öz (187c AYNEN, tam örneklem) — Σ_s K(b,s) = ζ'nin (0.86,1.30]
payı. **Ĝ_mid(ω) = mean_n e^{iω m_n}, MUTLAK m_n**; S(s) = Σ a'|Ĝ_mid(ω')|,
W(s) = Σ a'. İptal-yönlü bileşen κ(s) = −Re K_HAVUZ(s). Bantlar: 41 ince
(0.01) + 184'ün 8 bandı + HAVUZ [0.45,0.86). Jackknife: 8-blok loo (184-187
AYNEN; alt-örneklem noktası orijinal indisinin bloğuna düşer).

---

## K0 — YENİDEN KURULUM KAPISI (scratchpad temizdi; veri zinciri sıfırdan)

**Yapılan (içerik değişikliği YOK):**
- `155/eta_son_t0.4_c4000.npz`: 155'in KENDİ üreticisiyle
  (155_kos.veri_yukle('son') + 155_cekirdek.eta_onbellek(z,'son',0.4,4000);
  sürücü `scratchpad/188/k0_eta_son.py`): N=299 999, L=12.029593242,
  var(ds)=0.16744; tutarlılık maks|g − (ds+1)·2π/log(mid/2π)| = 2.2e-16.
- `z_Hkeskin.npy`: `164_insa.py Hkeskin` (nohup, 10.6 dk): **maks|F| =
  1.863e-9 ≤ 1e-8, sıralılık TAM** (min Δz = 0.1067, sıra bozan 0),
  σ_ds² = 0.1861, **L = 12.02959324** (184'teki inşa ile aynı: 1.9e-9 / TAM /
  0.186).
- ONKAYIT_184…187 JSON'ları silinmişti: `184a…187a_onkayit.py` betikleri
  DEĞİŞTİRİLMEDEN yeniden koşuldu → **sha256'lar raporlardakilerle aynı**
  (1e0de637… / d49a52a4… / 23a8130a… / 42a91c3a…); damgalar yeni
  (Wed Sep 23 12:14:10 +03 2026). Betiklere dokunulmadı.
- 184b yalnız 'gercek' gazı için ince bir sürücüyle koşuldu
  (`scratchpad/188/k0_184b_gercek.py`, 184b.kos() AYNEN): ONKAYIT_184'ün
  'ikiz'=keskin gazı 152'nin silinmiş z_keskin'ini ister ve 184-187 zincirinin
  ikizi zaten 184b2 ile Hkeskin'dir. Zincir: 184b(gercek) → 184b2(Hkeskin) →
  185b → 186b → 187c (`scratchpad/188/k0_zincir.sh`).

**Yeniden üretilen sayılar (hepsi HANE HANE):**

| nicelik | 184-187 raporu | 188 yeniden kurulum | betik |
|---|---|---|---|
| w_gerçek 8 bant | 1.170±.001 … 1.256±.013 | 1.170±.001, 1.195, 1.217, 1.241, 1.257, 1.269, 1.270, 1.256±.013 | 184b→185b |
| w_Hkeskin 8 bant | 1.239±.002 … 1.406±.016 | 1.239±.002, 1.278, 1.315, 1.354, 1.382, 1.406, 1.415, 1.406±.016 | 184b2→185b |
| medyan w (g / Hk) | 1.25 / 1.389 | 1.2516 / 1.3890 | 184b / 184b2 |
| H-W1a yapısal kapanış | ~1e-15 | 2.2e-16 | 185b |
| 186 G1 χ²/dof; r_pred(ρ≡1)>1 | 2336; EVET | 2336.45; EVET | 186b |
| pencere-içi pay m^kesik/m_ölç | 1.46→1.58 | 1.479→1.584 | 186b |
| **ζ_g(HAVUZ)** | **0.3287±0.0023 ∠+179.96°±0.07** | **0.3287±0.0023 ∠+179.96°±0.07** | 187c |
| **ζ_Hk(HAVUZ)** | **0.0768±0.0064 ∠+179.76°** | **0.0768±0.0064 ∠+179.76°±0.34** | 187c |
| 8-bant \|ζ_g\| | 0.3284 … 0.3811 | 0.3284, 0.3230, 0.3184, 0.3182, 0.3184, 0.3274, 0.3456, 0.3811 (se'ler dahil aynı) | 187c |

w-bantları ve ζ tablosu (gerçek + ikiz, 9 satır, modül/açı/se/genlik-okuma)
187 K2 tablosuyla birebir aynı. **K0 GEÇTİ — sapma yok.**

---

## K0b — DONMUŞ ÖN-KAYIT  [188a_onkayit.py]

`188/ONKAYIT_188.json`: **sha256 = 0a73fa37cece39c2…**, damga **Wed Sep 23
12:27:41 +0300 2026** — K0'ın son adımından (187c, 12:27:17) SONRA, ilk harita
dilim dosyasından ÖNCE. (Betik K0 koşarken yazıldı, K0 geçtikten sonra
damgalandı.) Donan: dilim/ikiz/ince-bant ızgaraları, seri/izdüşüm/K/Ĝ/S/W
tanımları, uydu konumları 1±τ_p (p=2,3,5,7) ±0.006 ve dilim üyelik kuralı
(merkez-içinde; yarı-örtüşme kuralıyla 8 uyduda da aynı üyelik — denetlendi),
κ tanımı, H-188a..d eşikleri + hüküm kuralları (H-188c tutarsa H-188a ölür),
K3 ayrışım kuralları, İKİNCİL kayıtlar (gürültü-düzeltmeli rang-1 payı;
imzalı öngörücü S_Re = Σ a' πτ' cos(πτ') Re Ĝ ile corr(Re v, S_Re)).

---

## K1 — HARİTA (gerçek kinematik)  [188b_harita.py]

q ≤ e^{1.30·L} = 6 190 158: 425 479 asal-kuvvet; pencere-ötesi (0.86,1.30]
**422 054 çizgi**; τ'≤1.20 kısmı **136 045** (187 katman toplamıyla aynı).
88 dilim 6 süreçle 387 s; kontrol noktaları `188/dilim_<i>.npz`.

**Makine mührü:** 188'in seri fonksiyonu (Ĝ eklentili) 187b.katman_serisi ile
**bit-bit aynı** (maks|Δ| = 0.0) [188b; ikiz koşusunda da 0.0, 188c].

**Tam-örneklem kontrolü** (K_alt vs K_tam, karmaşık) [188b]:

| dilim | K_alt(HAVUZ) | K_tam(HAVUZ) | \|Δ\|/se_alt | 8-bant \|Δ\|/se | 41-bant karmaşık korr |
|---|---|---|---|---|---|
| (1.050,1.055] | −0.01679−0.00032i | −0.01612−0.00009i | 0.07 | 0.06–0.25 | 0.984 |
| (0.900,0.905] | −0.00216−0.00000i | −0.00191+0.00003i | 0.15 | 0.09–0.80 | 0.799 |

Alt-örneklem yansız (tüm farklar <1 se). **K1 tutarlılık:** HAVUZ
Σ_{τ'≤1.20} K = **0.2261±0.0026 ∠−179.82°±0.27** — 187 K3'ün Σζ_i(≤1.20) =
0.228∠180°'siyle uyumlu (fark 0.002 < 1 se). Tüm derinlik: Σ_{τ'≤1.30} K =
**0.2566±0.0032 ∠+179.92°±0.24** = ζ(0.3287)'nin **%78.1'i** [188b].

**Uyarı (jackknife yorumu):** tek-dilim K'larının jk se'si tepelerde çok büyük
(ör. (1.050,1.055] HAVUZ se = 0.0100, değer 0.0168) ama tepeyi bütün içeren
toplamlarda küçük. Sebep gürültü DEĞİL, blok-blok yerel L kaymasıdır (bkz.
KEŞİF (c)): her jk bloğunda tepe kendi L_b'sinin öngördüğü dilime kayıyor.
Bu yüzden İKİNCİL "gürültü-düzeltmeli rang-1 payı" TANIMSIZ çıktı (1.42 > 1;
hücre-jk varyansı gürültüyü değil L-kaymasını ölçüyor) [188d].

---

## K2 — İKİZ HARİTASI (Hkeskin; τ'≤1.20, dilim 0.01)  [188c_ikiz_harita.py]

34 dilim, 136 045 çizgi, aynı makine (mühür 0.0). HAVUZ Σ_{τ'≤1.20} K_Hk =
**0.2430±0.0028 ∠−179.88°±0.38** [188c]. Derinlik-eşli karşılaştırma
(|Σ K_HAVUZ| ± jk) [188f]:

| derinlik | gerçek | ikiz |
|---|---|---|
| τ' ≤ 0.90 | 0.0132±0.0010 | 0.0126±0.0010 |
| τ' ≤ 1.00 | 0.0758±0.0067 | **0.0772±0.0065** |
| τ' ≤ 1.10 | 0.1717±0.0028 | 0.1795±0.0030 |
| τ' ≤ 1.20 | 0.2261±0.0026 | 0.2430±0.0028 |

İkizin kendi içeriği τ'≤1.00'de biter; haritası o derinlikte **0.0772±0.0065**
veriyor — 187c'nin ölçtüğü ζ_Hk = **0.0768±0.0064** ile birebir (alt-örneklem
içinde). Aynı derinlikte gerçek deniz 0.0758 veriyor.

---

## K3 — ANALİZ  [188d_analiz.py]

**(i) Rang-1 (41×88 karmaşık K, ağırlıksız SVD):** rang-1 varyans payı
**0.773 ± 0.024** (ikinci tekil değer 0.025) [188d]. Eşik ≥0.80 TUTMADI
(1.1σ altında), ölüm eşiği 0.60'ın üstünde. u(τ_b) açısı her ince bantta
**180°±5.1°**; |u| 0.0446 (τ_b=0.665) … 0.0642 (0.845). v(τ') faz-sabitli
(Σv = 5.117, reel); |v|>0.1 olan 14 dilimde açı ±3.4° içinde. 88 dilimin
72'sinde κ>0 (Σκ⁺ = 0.268, Σκ⁻ = −0.012). Harita görsel olarak DİKEY
ŞERİTLİ: iptal τ'-tarağında, bant-bağımsız.

**(ii) Korelasyonlar (88 dilim)** [188d]:

| nicelik | değer |
|---|---|
| corr(\|v\|, S) | **0.127 ± 0.040** |
| corr(\|v\|, W) | −0.108 ± 0.018 |
| fark | 0.235 ± 0.023 |
| İKİNCİL corr(Re v, S_Re) | **−0.948 ± 0.006** |

S = Σ a'|Ĝ| izleyici DEĞİL: |Ĝ| çoğu çizgide 1/√N'lik tutarsız tabanda
(S/W ≈ 0.0023-0.0051) ve bu taban W ile birlikte (0.86→1.30'da W 8.8×)
yükselir; |v| yükselmez.
Ön-kayıtlı İKİNCİL imzalı öngörücü S_Re (Re Ĝ, πτ'cos πτ' ağırlığı; tutarsız
fazlar toplamda söner) v'yi |corr| = 0.95 ile izliyor (işaret beklenen:
cos πτ' < 0, v pozitif fazlı). Bu hükmü DEĞİŞTİRMEZ (bkz. HÜKÜM).

**(iii) (1.05,1.10] tepesinin kimliği** [188d; pencere dökümü 188f(f)]:
(1.05,1.10] κ-toplamı 0.0735±0.0022'nin **%96.4 ± 3.1'i** 1+τ_2 (1.0576) ve
1+τ_3 (1.0913) pencerelerinde (sıfır beklentisi: dilim payı %50). Döküm:
1+τ_2 **0.0560±0.0008** (tepenin %76'sı), 1+τ_3 0.0149±0.0030 (%20).
Ön-kayıtlı tüm pencereler (HAVUZ Σκ ± jk): Bragg 0.0236±0.0042; 1+τ_2
0.0560; 1−τ_2 0.0134±0.0008; 1+τ_3 0.0149; 1−τ_3 0.0064±0.0014; 1+τ_5
0.0047±0.0014; 1−τ_5 0.0043±0.0009; 1+τ_7 0.0019±0.0012 (1−τ_7 aralık dışı).
Pencereler (dilimlerin %21.6'sı) (0.86,1.30] iptalinin **%48.8'ini** taşıyor.
**5 ve 7'nin uyduları ZAYIF**; en güçlü pencere-dışı tepe **τ'≈1.15**
(1.1425–1.1525: Σκ = 0.0349±0.0008 — Bragg'den büyük), 1+τ_5 ve 1+τ_7
pencerelerinin arasında.

**(iv) Köşegen şerit (rang-1 artığı R, τ_b ≥ 0.70, ±0.02)** [188d]:
kontrast = **0.839 ± 0.060** (eşik ≥3; 128 şerit hücresi; katlanma şeridi
τ'=2−τ_b tek başına 0.817, yakın şerit τ'−τ_b 1.503). **Köşegen sırt YOK.**
Yükseliş payı: ζ^{≤1.30}'un kesime yükselişi Δz = 0.029±0.013 (B_min =
0.65-0.70); artığın payı **−0.92 ± 0.89** (artık yükselişe KARŞI çalışıyor).

**(v) H-188d (κ-profili, 34 dilim)** [188d]: corr(κ_gerçek, κ_ikiz) =
**0.9988 ± 0.0003** (İKİNCİL |K|-profilleri 0.9988±0.0003). İkiz profili
tepelerde %5-15 daha güçlü, şekil aynı.

**(vi) K3 — 8-bant ζ^{≤1.30} yeniden kurulumu** (z = −Re ζ; ± jk) [188d]:

| bant | \|ζ\| tam (187c) | z = ζ^{≤1.30} | açı | z_rang1 | z_artık | z_r1/\|ζ\|tam | τ'>1.30 kuyruk payı |
|---|---|---|---|---|---|---|---|
| 0.45-0.50 | 0.3284±.0032 | 0.2744±.0093 | +179.3° | 0.2577±.0093 | +0.0167±.0085 | 0.785 | 0.164 |
| 0.50-0.55 | 0.3230±.0032 | 0.2671±.0056 | −179.7° | 0.2415±.0070 | +0.0256±.0062 | 0.747 | 0.173 |
| 0.55-0.60 | 0.3184±.0060 | 0.2576±.0084 | +179.8° | 0.2387±.0029 | +0.0188±.0067 | 0.750 | 0.191 |
| 0.60-0.65 | 0.3182±.0030 | 0.2576±.0082 | +178.1° | 0.2412±.0073 | +0.0164±.0075 | 0.758 | 0.190 |
| 0.65-0.70 | 0.3184±.0039 | 0.2362±.0033 | −177.7° | 0.2410±.0032 | −0.0048±.0037 | 0.757 | 0.258 |
| 0.70-0.75 | 0.3274±.0041 | 0.2366±.0078 | +177.4° | 0.2485±.0054 | −0.0119±.0076 | 0.759 | 0.277 |
| 0.75-0.80 | 0.3456±.0052 | 0.2678±.0073 | −178.0° | 0.2745±.0098 | −0.0067±.0070 | 0.794 | 0.225 |
| 0.80-0.86 | 0.3811±.0072 | 0.2653±.0116 | −178.8° | 0.2969±.0096 | −0.0316±.0066 | 0.779 | 0.304 |

(ince→8 bant toplama kontrolü maks|Δ| = 1.1e-16 [188d]; son iki sütun 188f(d);
kuyruk payı (|ζ|tam − z)/|ζ|tam, tam- ve alt-örneklem farkı olduğundan jk se
**erişilemedi**.)

Bant-bant okuma: **orta taban VE kesime yükseliş RANG-1 bileşenindedir** —
z_r1 profili tam-derinlik |ζ| profilini ~0.76 sabit ölçekle izliyor
(z_r1/|ζ| = 0.747-0.794; corr(z_r1, |ζ|tam) = **0.976**) [188f(d)]; yani
şekli τ'-tarağı değil bant çarpanı u(τ_b) taşıyor. Artık alt bantlarda +
(0.017-0.026), kesimde − (−0.032): yükselişi TAŞIMIYOR, törpülüyor
(0.65-0.70 → 0.80-0.86: z +0.029 = z_r1 +0.056 + z_artık −0.027). Harita
toplamı z tek başına tam profili izlemiyor (corr 0.351): tam-derinlik
yükselişinin (0.063) ancak ~yarısı τ'≤1.30'da; τ'>1.30 kuyruk payı kesime
doğru 0.16 → 0.30 büyüyor.

**α=1.30 bağı: AÇIK.** Ön-kayıtlı kazanan çekirdek yok (H-188a öldü);
ζ-makası → M → r zinciri bu kalemde ileri-hesaplanmadı.

---

## ÖN-KAYITSIZ KEŞİF (harita görüldükten SONRA; hüküm DIŞI, KAYIT)  [188f_kesif.py]

(a) **187'nin jk-se ipucu çözüldü.** Harita-ζ_i 187 katman ızgarasında
(HAVUZ): 0.0132±.0010 / 0.0290±.0020 / 0.0336±.0064 / 0.0224±.0066 /
**0.0735±.0022** / 0.0247±.0060 / 0.0297±.0060 (+ (1.20,1.25] 0.0098±.0008,
(1.25,1.30] 0.0207±.0009) — 187'nin ι profiliyle aynı şekil ((1.05,1.10] ~3×).
Büyük se'li katmanlar, sınırında bir tepeyi İKİYE BÖLEN katmanlardır:
Bragg (1.00) (0.95,1.00]|(1.00,1.05]'i, τ'≈1.15 tepesi (1.10,1.15]|(1.15,1.20]'yi
böler; blok-blok L-kayması tepeyi sınırın iki yanına taşıyınca se büyür.
Küçük se'li katmanlar tepelerini bütün içerir. KALEM'in "1+τ_5/1+τ_7 iz
KISMİ" gözlemi böyle açıklanır: o katmanlarda 5/7 uyduları zayıf, büyük se
1.15 tepesinin bölünmesinden.

(b) **Derinlik-eşli ζ iki denizde aynı** (K2 tablosu): τ'≤1.00'de gerçek
0.0758 vs ikiz 0.0772 = ζ_Hk. ζ farkı (0.33 vs 0.077) çekirdek farkı değil,
içeriğin DERİNLİĞİ: ikizin içeriği 1.00'de biter, gerçeğinki sürer.

(c) **Tepeler yerel L ile kayıyor (tarak kimliği doğrudan):** yalnız-blok
K_HAVUZ tepe konumları 8 blokta L_b/L (Bragg), (L_b+log 2)/L (1+τ_2) ve
(L_b+log 6)/L (τ'≈1.15) öngörülerini maks|Δ| = 0.0022 / 0.0025 / 0.0024
(dilim 0.005'in yarısı) ile izliyor; corr(öngörü, tepe) = 0.945 / 0.926 /
0.909. τ'≈1.15 tepesi bir uydu gibi davranıyor (L + sabit, sabit ≈ log 6 =
log 2 + log 3, ±0.06); 0.005 çözünürlükte küçük rasyonellerin yoğunluğu
yüzünden etiket TEKİL DEĞİL — sınanmadı. Diğer tepeler: simetrik çift
0.9625/1.0375 (1∓0.035), 0.9125/1.0875, 1.1925, 1.2225, 1.2825.

(e) Naif birinci-mertebe bant çarpanı U_T = Σ c^ölç conj(karışım)/Σ|karışım|²
u'yu izlemiyor: corr(|u|,|U_T|) = −0.444 (41 bant) — bant çarpanının kimliği
AÇIK.

---

## HÜKÜM (eşikler K0b'de donmuş; kurtarma yok)

| hipotez/kapı | hüküm | dayanak |
|---|---|---|
| **K0** yeniden kurulum | **GEÇTİ** | w-bantları, ζ_g=0.3287±0.0023∠179.96°, ζ_Hk=0.0768±0.0064, 8-bant \|ζ_g\| hane hane [184b-187c] |
| makine mührü (seri bit-bit) | **TUTTU** | maks\|Δ\| = 0.0 [188b, 188c] |
| tam-örneklem kontrolü | **TUTTU** | \|Δ\| = 0.07 / 0.15 se_alt [188b] |
| K1 tutarlılık (Σ_{≤1.20} K vs 0.228) | **TUTTU** | 0.2261±0.0026∠−179.8° [188b] |
| **H-188a** tarak köprüsü | **ÖLDÜ** | (i) 0.773±0.024 < 0.80 (ölüm eşiği üstü); (ii) corr(\|v\|,S) = 0.127 < 0.50 → ÖLÜM; (iii) 0.964±0.031 ≥ 0.60 geçti [188d] |
| **H-188b** köşegen/kesime-yakınlık | **ÖLDÜ** | kontrast 0.839±0.060 < 3; artık payı −0.92±0.89 [188d] |
| **H-188c** yayvan taban (sıfır) | **ÖLDÜ** | corr(\|v\|,W) = −0.108 < corr(\|v\|,S) = 0.127 (fark 0.235±0.023) [188d] |
| **H-188d** kinematik-genel | **MÜHÜR** | corr(κ_gerçek, κ_ikiz) = 0.9988±0.0003 ≥ 0.80 [188d] |
| **K3** yeniden kurulum | **KAYIT** | taban + yükseliş rang-1 bant çarpanında (corr(z_r1,\|ζ\|) = 0.976); artık yükselişe karşı; yükselişin ~yarısı τ'>1.30; α=1.30 bağı AÇIK [188d, 188f] |
| İKİNCİL kayıtlar | KAYIT | corr(Re v, S_Re) = −0.948±0.006; gürültü-düzeltmeli pay TANIMSIZ (1.42) [188d] |

**H-188a ölümünün dürüst okuması:** ölen, ön-kayıtlı ÖLÇÜTTÜR — "|v| modül-S'yi
izler" ve "rang-1 ≥ 0.80". Tarak resminin konum kısmı ölçümde duruyor ((iii)
%96; tepeler Bragg + 1±τ_2, 1±τ_3 + yerel L ile kayan τ'≈1.15), ama bu kalemin
hükmünü değiştirmez; imzalı-Ĝ öngörücüsü (İKİNCİL, |corr| 0.95) yeni bir
ön-kayıt gerektirir.

---

## TÜRETİLEN vs ÖLÇÜLEN

Bu kalemde **hiçbir şey türetilmedi.** Uydu/Bragg KONUMLARI (1±log p/L, L_b/L)
kinematik öngörüdür ve KALEM'de ölçümden önce yazılmıştı; bütün GENLİKLER,
paylar, korelasyonlar ve K3 ayrışımı ÖLÇÜLENDİR (harita, kesin beklenen-değer
izdüşümü; Taylor yalnız yorum). ζ'nin ~1/3 büyüklüğü ve u(τ_b)'nin bant
profili türetilmedi: açık.

---

## MANŞET (aday cümle)

> **İptalin çekirdeği haritalandı: pencere-ötesi iptal bant-bağımsız bir
> τ'-TARAĞIDIR — K(τ_b,τ') %77 rang-1, bant çarpanı her bantta 180°±5°, tepeler
> Bragg (τ'=1.00) ve uydularda; 187'nin (1.05,1.10] tepesinin %96'sı 2 ve 3'ün
> uydularında (1+τ_2 tek başına %76), tepeler yerel L ile blok-blok kayıyor.
> Çekirdek kinematik-geneldir (ikiz profili corr 0.999, H-188d MÜHÜR) ve
> derinlik-eşli ζ iki denizde aynıdır (τ'≤1.00: 0.076 vs 0.077 = ζ_Hk):
> 0.33 ile 0.077 arasındaki fark çekirdek değil, içeriğin derinliğidir.** Üç
> ön-kayıtlı ölçüt öldü: |Ĝ|-ağırlıklı yapı çarpanı v'yi izlemiyor (corr 0.13;
> H-188a ÖLDÜ — imzalı Re Ĝ öngörücüsü ikincil kayıtta 0.95 ile izliyor), köşegen
> sırt yok (kontrast 0.84; H-188b ÖLDÜ), yayvan sıfır da öldü (H-188c). ζ
> profilinin orta tabanı ve kesime yükselişi TARAKTA değil bant çarpanı
> u(τ_b)'dedir (corr 0.976); yükselişin yarısı τ'>1.30 kuyruğundadır. α=1.30
> bağı açık.

---

Teslim: bu rapor + `188_configs/` (188a_onkayit, 188b_harita, 188c_ikiz_harita,
188d_analiz, 188e_figur, 188f_kesif) + `188_zeta_kimligi.png` (sol üst |K|
haritası + uydu/Bragg/katlanma çizgileri; sağ üst |v| vs S, W (+ İKİNCİL
−S_Re) normalize; alt 8-bant |ζ| ölçülen vs z, z_rang1, z_artık) +
`scratchpad/188/` (ONKAYIT_188.json, dilim_0..87.npz, ikiz_dilim_0..33.npz,
harita_proj_{gercek,Hkeskin}.npz, kontrol_tam_gercek.npz,
harita_K_{gercek,Hkeskin}.npz, K1_harita_{gercek,Hkeskin}.json,
K3_analiz.{json,npz}, K_kesif.json, K0 sürücüleri k0_eta_son.py /
k0_184b_gercek.py / k0_zincir.sh, tüm loglar) + yeniden kurulan zincir
önbellekleri (`scratchpad/{155,164,184,185,186,187}/`).
