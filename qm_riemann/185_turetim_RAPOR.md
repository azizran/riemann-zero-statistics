# 185 — ZARFIN TÜRETİMİ (neden 1 − 0.149·τ^1.30?)

**Soru (KALEM_W_TURETIMI_09EYL2026):** 184 zarfın yüzünü ölçtü —
r(τ) = â_gerçek/â_Hkeskin = 1 − 0.149·τ^1.30 (χ²/dof=0.34). Bu biçim,
132c'nin KESİN özdeşliğinden, öz-izdüşüm ayrışımıyla PARAMETRESİZ
türetilebilir mi? Zarfı hangi faktör taşıyor?

Her sayının yanında onu üreten betik adı vardır. Tek dalga; ölümler
kurtarmasız; git'e dokunulmadı.

---

## Ölçüm tanımı (KALEM cebiri; ONKAYIT_185 ile donmuş)

Özdeşlik kesin: `ds_n = Σ_Q 2a_Q sin(ω_Q g_n/2)·cos(ω_Q m_n)`. Tek-çizgi
ÖZ-İZDÜŞÜMÜ (Taylor YOK, kesin beklenen-değer; a_q NOMİNAL):

```
c_q^öz = (4a_q/N) Σ_n sin(ω_q g_n/2)·cos(ω_q m_n)·e^{−iω_q m_n},   â_q^öz = |c_q^öz|
```

Üç-faktörlü KESİN ayrışım (çizgi düzeyinde teleskopik):

```
r(τ) = [â_g/â_g^öz] × [â_g^öz/â_Hk^öz] × [â_Hk^öz/â_Hk]
        (F_ANOMALİ)     (F_KİNEMATİK)      (F_KOMŞU⁻¹)
```

Kinematik kaynaklar: gerçek gazın g_n'si eta_son önbelleğinden KESİN
inversiyonla `g=(ds+1)·2π/log(mid/2π)` (kalıntı 2.3e-10; 184b veri yolu
yalnız mid+ds taşır); Hkeskin/HA4 için z-dosyadan (184b2 yolu). Ölçülü
â'lar 184 K1 npz'lerinden AYNEN (yeniden ölçüm yok). Bant/jackknife
184'ünkiyle özdeş; bant değeri Σâ/Σae, r = Σâ_g/Σâ_Hk; birincil σ = 184
defter konvansiyonu (w_g, w_Hk jk-se'lerinin bağımsız yayılımı — 184'ün
yayımlanmış ±.002…±.013'ünü hane hane yeniden üretir).

---

## K0 — DONMUŞ ÖN-KAYIT  [185a_onkayit.py]

`185/ONKAYIT_185.json` yazıldı: **sha256 = d49a52a4…**, damga
**Wed Sep 9 14:20:47 +03 2026** (tüm ölçümlerden önce). Donan: â^öz/Ĝ/üç
faktör tanımları; bant ızgarası `[.45 .50 .55 .60 .65 .70 .75 .80 .86]` +
kuyruk τ>0.70 (184 AYNEN); 8-blok loo-jackknife; H-W0..W3 formları ve ölüm
eşikleri; H-W2 ileri-hesap adayları (A çarpımsal-kinematik, B toplamsal
ortak-karışım) ve mekanizma-eşleme kuralı; K4'ün birincil biçimi (toplamsal
öz-değişim, erfc 0.68/0.125 donmuş) ölçüm-öncesi gerekçesiyle.

---

## K1 — ÖZ-MUHASEBE: üç faktörün bant defteri  [185b_oz_muhasebe.py]

σ_ε (aralık-saçılımı): gerçek **0.4092**, Hkeskin **0.4313** — πτσ_ε ~ 1,
Taylor'un neden yasak olduğunun kaydı. Toplam-oranı lehçesi, ±jk se:

| bant τ | τ̄ | r=g/Hk | F_ANOMALİ | F_KİNEMATİK | F_KOMŞU⁻¹ |
|---|---|---|---|---|---|
| 0.45–0.50 | 0.475 | 0.9445±.0017 | 1.407±.001 | 1.0221±.0005 | 0.6568±.0012 |
| 0.50–0.55 | 0.526 | 0.9351±.0019 | 1.494±.001 | 1.0290±.0007 | 0.6083±.0010 |
| 0.55–0.60 | 0.576 | 0.9253±.0033 | 1.584±.004 | 1.0366±.0008 | 0.5634±.0015 |
| 0.60–0.65 | 0.626 | 0.9169±.0040 | 1.682±.003 | 1.0465±.0010 | 0.5208±.0019 |
| 0.65–0.70 | 0.676 | 0.9096±.0043 | 1.772±.003 | 1.0593±.0013 | 0.4845±.0017 |
| 0.70–0.75 | 0.725 | 0.9021±.0044 | 1.852±.007 | 1.0760±.0015 | 0.4527±.0015 |
| 0.75–0.80 | 0.775 | 0.8978±.0081 | 1.904±.007 | 1.0991±.0018 | 0.4289±.0021 |
| 0.80–0.86 | 0.830 | 0.8931±.0135 | 1.897±.012 | 1.1374±.0020 | 0.4138±.0035 |
| **KUYRUK>0.70** | 0.782 | 0.8972±.0078 | 1.886±.005 | 1.1062±.0018 | 0.4301±.0020 |

- **H-W1a (yapısal kapanış):** toplam-oranı lehçesinde ΠF ≡ r; maks fark
  **2.2e-16** — makine sızdırmıyor. [185b]
- **H-W1b (çizgi-lehçesi sınavı):** ae-ağırlıklı çizgi-faktör
  ortalamalarının çarpımı r̂_b, 184 defterine karşı **χ²/dof = 0.197**
  (maks |z| = 0.94, τ=0.83 bandında; eşik 2) → üç faktör 184'ün 8 bandını
  parametresiz yeniden ÜRETİYOR. [185b]
- **TUTARLILIK — F_KOMŞU = â_Hk/â_Hk^öz:** 1'e yakın DEĞİL:
  **1.522 → 2.417** (bant 1→8). İkizin ds'i TAM nominal merdiven olduğu
  hâlde ölçülü çizgi genliği, öz-teriminin 1.5–2.4 katı: aynı apsisteki
  öteki çizgilerin doğrusal-olmayan karışımı (komşu payı) BÜYÜKTÜR ve
  τ ile büyür. w>1 mühürlerinin (184: w_Hk=1.39) kaynağı da budur —
  öz-terim tek başına w^öz<1 verir (aşağıda). [185b]
- **w^öz profilleri:** w_g^öz = 0.831→0.662, w_Hk^öz = 0.813→0.582 —
  kesin beklenen-değer, geniş aralık-saçılımıyla nominalin ALTINA iner;
  ölçülü w'lerin 1'i aşması tümüyle karışım katkısıdır. [185b]

**K1 HÜKMÜ:** öz-muhasebe kapanıyor (H-W1 MÜHÜR) ama beklenen resmin tersi:
zarf tek-çizgi teriminde değil — hem gerçekte hem ikizde ölçülü genliğin
yarısından fazlası komşu-karışım kanalında yaşıyor.

---

## K2 — ANATOMİ: zarfı kim taşıyor?  [185c_anatomi.py]

**Log-paylar** (pay_i = log F_i/log r; Σ=1 kesin), kuyruk τ>0.70:
F_ANOMALİ **−5.85**, F_KİNEMATİK **−0.93**, F_KOMŞU⁻¹ **+7.78**;
çift F_ANOM×F_KOMŞU⁻¹ payı **+1.93**. Yani:

- **F_KİNEMATİK zarfın TERSİ yönünde** (1.022→1.137): gerçeğin aralık
  saçılımı ikizden DAR (Δσ_ε² = −0.0186) ⇒ öz-terimi ikizden BÜYÜK; bu
  faktör tek başına olsaydı r kuyrukta 1.14 olurdu, 0.89 değil.
- Zarfı taşıyan, **F_ANOMALİ×F_KOMŞU⁻¹ çifti** (kuyrukta log-payın
  ~%193'ü; kinematik −%93'ü geri alıyor).

**Taylor-teşhisi (yalnız yorum):** ΔT1 (uyumlu geri-besleme farkı,
Ĝ-tabanlı) ihmal edilebilir (≤0.0008); T2 = −(πτ)²Δσ_ε²/2 = +0.021→+0.063,
log F_KİN'in düşük bantlarda ~%95'ini, kuyrukta ancak ~%49'unu açıklıyor
(0.063 vs 0.129) — πτσ_ε~1'de Taylor'un çöktüğünün, kesin hesabın neden
şart olduğunun kaydı. [185c]

**Ĝ(τ) profili (köprü kaydı):** Re Ĝ_g = 0.0108→0.0007 (τ ile hızla düşer),
Hkeskin hep ~%8-10 üstünde; 2ω kopyaları ≤0.0004. Öz-izdüşümün DC/2ω
anatomisi: 2ω payı ≤%0.3 — â^öz pratikte ⟨sin(ω_q g/2)⟩ istatistiğidir. [185c]

**KARIŞIM-ORANI M(τ) (yorum türevi; hizalı-faz yaklaşımı, K3-B'nin bant
kapanışıyla doğrulanır):** m = w_ölçülü − w^öz kanalında
M = m_g/m_Hk: [185f_figur.py]

| τ̄ | 0.475 | 0.526 | 0.576 | 0.626 | 0.676 | 0.725 | 0.775 | 0.830 |
|---|---|---|---|---|---|---|---|---|
| M | 0.796±.004 | 0.789±.004 | 0.782±.007 | 0.776±.008 | 0.769±.008 | 0.758±.008 | 0.747±.013 | 0.721±.021 |

Gerçek gaz, ikizin uyumlu karışım genliğinin yalnız **%80→%72**'sini
iletiyor — zarfın gerçek taşıyıcısının doğrudan ölçümü; Not 2 w-kanalının
(w₂=0.91 @ τ=0.028; v-kanalı ~0.55-0.60'a doyar) τ∈[0.45,0.86]'daki
keşfedilmemiş devamı budur.

**H-W0 — pencere içi L kayması:** L_min=11.953, L_maks=12.102,
ΔL_rms=0.043 (%0.36); δτ_rms(kuyruk)=0.0028. Zarfa net birinci-mertebe
etki |dr/dτ|·δτ = 5.6e-4 → (1−r)'nin **%0.54'ü < %1 → KAPANDI**. (Yan not:
tek-faktör F_KİN profil-yayılması sınırı 1.8e-3, %1.70 — r'ye ortak-mod
τ-etiketi üzerinden girer, net etki üstteki sayıdır.) [185c]

---

## K3 — KAPANIŞ: ileri-hesaplar ve türetilmiş (c,α)  [185d_kapanis.py]

Mekanizma-eşleme (K0'da donmuş): F1/F3 çifti baskın → birincil = **B**.
İki parametresiz ileri-hesap (ölçülü â_g HİÇ girmez), 184 defter-σ'sıyla:

| bant τ̄ | r ölçülü | r_pred^A (kinematik) | z_A | r_pred^B (ortak-karışım) | z_B |
|---|---|---|---|---|---|
| 0.475 | 0.9445±.0017 | 1.0221±.0005 | +46.5 | 1.0145±.0003 | +42.0 |
| 0.626 | 0.9169±.0040 | 1.0465±.0010 | +32.8 | 1.0242±.0005 | +27.1 |
| 0.725 | 0.9021±.0044 | 1.0760±.0015 | +39.9 | 1.0344±.0006 | +30.4 |
| 0.830 | 0.8931±.0135 | 1.1374±.0020 | +18.1 | 1.0565±.0007 | +12.1 |

**χ²/dof: A = 1321, B = 922** (dof=8; eşik 2) → **H-W2 ÖLDÜ,
kurtarmasız.** İki aday da zarfın TERS yönünde (r_pred>1): tek-çizgi
kinematiği de, "komşu-karışım ortak-mod + öz-terim değişimi" de zarfı
vermiyor. Zarf, karışım genliğinin KENDİSİNDEKİ bastırmada (M(τ), K2) —
kalem dilinde F_ANOMALİ dalı: bu dalda parametresiz türetim yok, doğrudan
ölçülmüş köprü var. [185d]

**Güç-yasası fitleri** (1−c·τ^α; özet amaçlı):

| eğri | c | α | fit-χ²/dof |
|---|---|---|---|
| ölçülü (kontrol) | **+0.1492** | **1.304** | 0.34 — 184'ün (0.149, 1.30, 0.34)'ü AYNEN |
| r_pred^A | −0.2115 | 3.115 | 0.67 |
| r_pred^B | −0.0701 | 2.161 | 0.15 |

Türetilmiş (c,α), hedef (0.149, 1.30)'un işaretçe bile karşısında —
türetim başarısızlığının dürüst özeti. (Kontrol satırı, 185 muhasebesinin
184 makinesiyle bit-bit uyumunun mührüdür.) [185d]

**H-W3 — script-100 şekil sınavı:** ölçek bilgisi bu okuma listesinde YOK →
**erişilemedi**; yalnız en-iyi-ölçek şekil sınavı: 1−c₃·τ^3.3 için en iyi
c₃=0.3824 → **χ²/dof = 70.1** (dof=7). α=3.3 ≠ 1.30 ve şekil, serbest
ölçekle bile ölü → mekanik ε(q,k) yasası zarfın taşıyıcısı DEĞİL
(dürüst kayıt). [185d]

---

## K4 (bonus, örneklem-dışı) — HA4 SINAVI  [185e_HA4.py]

Aynı öz-muhasebe, erfc-ikiz HA4'e (σ_ε=0.4343; erfc 0.68/0.125 donmuş).
**Birincil K4-b (toplamsal öz-değişim):** c_pred = (c_Hk − c_Hk^öz) +
env·c_HA4^öz:

| bant τ̄ | w_HA4 ölçülü | w_HA4 pred | z |
|---|---|---|---|
| 0.475 | 1.3365±.0010 | 1.2243±.0024 | −43.6 |
| 0.626 | 1.2891±.0031 | 1.1639±.0054 | −20.2 |
| 0.725 | 1.0249±.0023 | 0.9714±.0052 | −9.3 |
| 0.775 | 0.9272±.0058 | 0.9021±.0086 | **−2.4** |
| 0.830 | 0.8716±.0092 | 0.8613±.0146 | **−0.6** |

bant-χ²/dof = **648 > 2 → ÇATLAK** (bonus sınav; hüküm dürüst). AMA:
kuyruk-çöküşünün kendisi (0.87'ye düşüş) iki kuyruk bandında **−2.4σ ve
−0.6σ ile parametresiz yakalandı** — "ölü merdiven + hayatta kalan ortak
karışım" mekanizması doğru; kırılan, orta-τ'da ikizler-arası
karışım-eşitliği (~%8-10: HA4'ün karışımı Hkeskin'inkinden büyük; ölçülü
r_HA4 düşük-τ'da 1.079 olması da bunu bağımsız söylüyor). [185e]

**Yan K4-a (çarpımsal taşıma):** kuyrukta z −295'e kadar, χ²/dof=30503 —
oran-modeli ölü merdivende sonlu genlik üretemez; komşu payının TOPLAMSAL
olduğunun nicel gömütü. [185e]

---

## HÜKÜM (ölüm eşikleri K0'da donmuş; kurtarma yok)

| hipotez | hüküm | dayanak |
|---|---|---|
| **H-W0** (L-akışı) | **MÜHÜR/KAPANDI** | zarfa net katkı %0.54 < %1 [185c] |
| **H-W1** (öz-muhasebe) | **MÜHÜR** | yapısal 2.2e-16; çizgi-lehçesi χ²/dof=0.197 ≤ 2 [185b] |
| **H-W2** (baskın-terim türetimi) | **ÖLDÜ** | birincil ileri-hesap χ²/dof=922 ≫ 2; her iki aday zarfın ters yönünde [185d] |
| **H-W3** (script-100 şekli) | **ÖLDÜ** (ölçek: erişilemedi) | en-iyi-ölçekle bile χ²/dof=70; α=3.3≠1.30 [185d] |
| K4 (HA4, bonus) | **ÇATLAK** (kuyruk-çöküşü yakalandı) | χ²/dof=648; kuyruk bantları −2.4σ/−0.6σ; çarpımsal −295σ [185e] |

H-W2'nin ölümü kalem'in öngördüğü ikinci dala düşer: baskın olan
F_ANOMALİ(×F_KOMŞU⁻¹) kanalıdır ve bu dalda türetim değil, DOĞRUDAN ÖLÇÜM
vardır — M(τ)=0.80→0.72 (Not 2 w-kanalının yüksek-τ devamı; Not 6'ya
köprü). "1.30 üssü nereden geliyor" sorusu AÇIK kalır; artık adresi
bellidir: karışım-genliği iletim profili M(τ).

---

## MANŞET (aday cümle)

> **Zarf (1−0.149·τ^1.30) özdeşlik-makinesinin İÇİNDE eksiksiz
> muhasebeleşiyor (H-W1: χ²/dof=0.20) ama tek-çizgi kinematiğinden
> TÜRETİLEMİYOR (H-W2: her iki parametresiz ileri-hesap zarfın ters
> yönünde, χ²/dof≥922 → öldü).** Kesin öz-izdüşüm ayrışımı zarfın
> taşıyıcısını buldu: kinematik değil (o, tersine, zarfı %14'e kadar
> GERİ sarıyor — gerçeğin aralık-saçılımı ikizden dar), komşu-çizgi
> KARIŞIM genliğinin bastırılması: gerçek gaz, sadakatli-keskin ikizin
> uyumlu karışım genliğinin yalnız **%80→%72**'sini iletiyor
> (M(τ), τ∈[0.45,0.86]) — Not 2'nin w-kanalının yüksek-τ'daki
> keşfedilmemiş devamının ilk doğrudan ölçümü. Aynı makine, HA4'ün
> kuyruk-çöküşünü (0.87) toplamsal öz-değişimle parametresiz yakaladı
> (kuyrukta −2.4σ/−0.6σ; çarpımsal taşıma −295σ ile gömüldü): komşu payı
> TOPLAMSALDIR, ve zarfın "neden 1.30?" sorusunun adresi artık
> karışım-iletim profili M(τ)'dir.

Teslim: bu rapor + `185_configs/` (185a–185f) + `185_turetim.png`
(sol: üç faktör + M(τ); sağ: r ölçülü vs iki ileri-hesap) +
`scratchpad/185/` (ONKAYIT_185.json, OZ_*.npz, K1_faktorler.npz,
K2_anatomi.json, K2_G_profil.npz, K3_kapanis.json, K4_HA4.json).
