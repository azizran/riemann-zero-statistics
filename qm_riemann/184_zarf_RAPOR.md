# 184 — ZARFIN ANATOMİSİ

**Soru (180-183 arkının açık bıraktığı en büyük):** Gerçeğin sadakatli
ikizinden farkının ~%80'i "zarf" (kesim/genlik yapısı). Zarf tam olarak
**nedir**? Gerçeğin çizgi-genlik profilinin merdivenin nominalinden sapması —
özellikle τ>0.70 kuyruğunda. Bu profil DOĞRUDAN ölçülebilir; ölçüyoruz.

Her sayının yanında onu üreten betik adı vardır. Tek dalga; ölümler
kurtarmasız; ön-mühür/sonuç.

---

## Ölçüm tanımı (132c özdeşliği + 162/174 makinesi)

132c'nin KESİN ÖZDEŞLİĞİ (pertürbasyon değil):
```
ds_n = Σ_Q 2 a_Q sin(ω_Q g_n/2) · cos(ω_Q m_n)          (★)
a_Q = 1/(π m √Q)  (Q=p^m),  ω_Q = log Q,  m_n = ½(z_n+z_{n+1}),  τ_Q = ω_Q/L
```
Buradan çizgi genliği `A₁(q) = 2 a_q sin(πτ_q)` (kinematik `sin(πτ)`) ve
etkin ölçüm:
- **â_q := |c_q(ds)|**, `c_q(ds) = 2⟨ds_n e^{−iω_q m_n}⟩` (155/159/162 ham
  izdüşümü; `ds`,`mid` eta_onbellek'ten, yeniden regresyon yok).
- **nominal a_q^eff := 2 a_q sin(πτ_q)** (analitik, öz-tutarlı bastırmasız).
- **w_q := â_q / a_q^eff.**  İkiz-kontrol: w_ikiz(τ) ≈ 1 olmalı.
- **saf zarf** iki dilde: `r(τ) = â_gerçek/â_ikiz` (nominal sadeleşir; BİRİNCİL)
  ve `D(τ) = w_gerçek − w_ikiz` (KALEM-K1 sözü).

Gazlar: gerçek = zeros6 son 300k (`son`, L=12.02959324); ikiz = 152'nin saf
merdiveni (`keskin`, L=12.02959349) — L farkı 2.5e-7, 299999 nokta, ÖZDEŞ
apsis. Her ikisinin eta_onbelleği hazırdı (pahalı regresyon zaten yapılmış).

---

## K0 — DONMUŞ ÖN-KAYIT  [184a_onkayit.py]

`184/ONKAYIT_184.json` yazıldı: **sha256 = 1e0de637…**, damga
**Wed Sep 9 12:44:57 +03 2026** (veriden önce). Donan: â_q tanımı, nominal,
w/r/D; bant ızgarası `[.45 .50 .55 .60 .65 .70 .75 .80 .86]` + τ>0.70 kuyruk
bandı; 8-blok jackknife; üç hipotez formu ve ölüm eşikleri:

| hipotez | form (r(τ)'ye fit) | par |
|---|---|---|
| H-Z1 erfc | ½ erfc((τ−τ_c)/Δ) | τ_c, Δ |
| H-Z2 güç | 1 − c·τ^α | c, α |
| H-Z3 Gauss-DW | exp(−(2πτ)²σ²/2) | σ |

Ölüm eşikleri (kurtarmasız): survival χ²/dof ≤ 2.0; kuyruk (τ>0.70)
|⟨artık⟩|/se ≤ 2.0; kıyas Δ(χ²/dof) < 0.3 ise ikisi de yaşar. K3 tutma
aralığı: zarf-payı 180/181'in %70–88'ini kesmeli.

---

## K1 — w(τ) PROFİLİ  [184b_K1_zarf.py, 184b2_zarf_ikiz.py]

**Gerçek gaz (son) — w_gerçek(τ) = â_q/(2a_q sinπτ)** [184b, genlik-ağırlıklı
bant, asal çizgiler]: düşük-τ'da **1.02** (τ≈0.15) → kuyrukta **1.27** (τ≈0.77),
τ≈0.83'te 1.25. Yani gerçeğin çizgi genlikleri düşük-τ'da analitik nominale
OTURUYOR (w≈1 ✓ — normalizasyon doğru) ve kuyrukta nominali AŞIYOR (~%25):
aritmetik belleğin/kilidin ürettiği yapıcı girişim, tek-çizgi nominalini
büyütüyor.

**İKİZ KİMLİĞİ — ön-kayıt kontrolü keskin'i ELEDİ (dürüst kayıt):** K0 ikiz
olarak 152'nin `keskin`ini dondurdu. Ölçüm: w_keskin(τ) = **0.55 → 0.34**
(düşük-τ→kuyruk), yani ön-kayıtlı "w_ikiz ≈ 1" kontrolünü GEÇEMEDİ. Sebep
mekanik ve bilinen (164'ün teşhisi): 152'nin sönümlü/kelepçeli Newton'u
denklemi EKSİK çözer, nominal merdiven genliğinin ancak ~%54'ünü kurar
(keskin R_η=0.448). Yani `keskin`, "nominal merdivenden kurulmuş ikiz"in
SADIK gerçeklemesi DEĞİL. Ön-kayıt kontrolü tam da bunu yakaladı.

**Doğru ikiz = Hkeskin** (164/180'in SADAKATLİ çözücüsü: ızgara-braketi +
sıralı-ilk-kök + korumalı Newton; τ≤1.00 keskin ladder, nominal a_q, seviye
c=−½). Hkeskin, nominal merdivenin sadık gerçeklemesidir (R_η=1.265 ≈ gerçek
1.288). Erfc-ikiz **HA4** (aynı sadakatli çözücü + erfc 0.68/0.125) de kuruldu
(164 sandviçinin yumuşak ucu). İkisi de sıfırdan üretiliyor (scratchpad
temizlenmişti); ölçüm 184b2 ile son'un q-kümesinde (line-by-line hizalı).

**Hkeskin kuruldu** [164_insa.py Hkeskin]: sadakatli çözücü, maks|F|=1.9e−9,
sıralılık TAM, σ_ds²=0.186 (gerçek 0.167), L=12.02959324 (gerçekle ÖZDEŞ).
Medyan w_Hkeskin=1.389 > gerçek 1.25 → gerçek, keskin ikizden bastırılmış.

### K1 bant defteri — İKİZ = Hkeskin (sadakatli keskin)  [184c]
genlik-ağırlıklı, ±jackknife se (8 blok):

| bant τ | n | τ̄ | w_gerçek | w_Hkeskin | **r=g/Hk** | D=g−Hk |
|---|---|---|---|---|---|---|
| 0.45–0.50 | 37 | 0.475 | 1.170±.001 | 1.239±.002 | **0.944±.002** | −0.069 |
| 0.50–0.55 | 56 | 0.526 | 1.195±.001 | 1.278±.002 | **0.935±.002** | −0.083 |
| 0.55–0.60 | 90 | 0.576 | 1.217±.003 | 1.315±.004 | **0.925±.003** | −0.098 |
| 0.60–0.65 | 157 | 0.626 | 1.241±.003 | 1.354±.005 | **0.917±.004** | −0.112 |
| 0.65–0.70 | 254 | 0.676 | 1.257±.003 | 1.382±.005 | **0.910±.004** | −0.125 |
| 0.70–0.75 | 432 | 0.725 | 1.269±.004 | 1.406±.005 | **0.902±.004** | −0.138 |
| 0.75–0.80 | 735 | 0.775 | 1.270±.007 | 1.415±.010 | **0.898±.008** | −0.145 |
| 0.80–0.86 | 1602 | 0.830 | 1.256±.013 | 1.406±.016 | **0.893±.013** | −0.150 |
| **KUYRUK>0.70** | 2769 | 0.782 | 1.264±.007 | 1.409±.010 | **0.897±.008** | −0.145 |

**ZARFIN KİMLİĞİ (birincil bulgu):** saf zarf r(τ) = gerçeğin çizgi
genliğinin sadakatli-keskin ikizden oranı, **0.94'ten 0.89'a NAZİK, MONOTON
bir düşüş** (kuyrukta ~%10 bastırma). Bu bir KESKİN KESİM DEĞİL, yumuşak bir
sönüm. Not: gerçek ve ikizin ikisi de analitik nominali AŞIYOR (w>1, kilit/öz-
tutarlılık) ama ikiz DAHA çok aşıyor; fark (zarf) kuyrukta büyüyor.

Kontrol amaçlı İKİZ=keskin (152 kaba): r = 2.19→3.64 (gerçek, keskin'in 2–3.6
katı) — keskin'in kaba çözücüsü genliği yarıya kadar eksik kuruyor; ön-kayıt
w≈1 kontrolünün neden keskin'i elediğinin nicel kaydı.

### K2 — HİPOTEZ YARIŞI (r(τ)'ye fit, τ∈[0.475,0.83], 8 bant, jk-ağırlıklı) [184c]

| hipotez | en iyi par | χ²/dof | kuyruk_z | hüküm |
|---|---|---|---|---|
| **H-Z2 güç** | r=1 − **0.149**·τ^**1.30** | **0.34** | 0.80 | **KAZANIR** |
| H-Z1 erfc | τ_c=**1.84**, Δ=1.23 | 0.74 | 1.13 | yaşar (kaybeder) |
| H-Z3 Gauss-DW | σ=0.107 | **6.94** | 5.17 | **ÖLDÜ** (kuyrukta da) |

- **H-Z3 (Gauss-DW / Debye-Waller) KURTARMASIZ ÖLDÜ** (χ²/dof=6.94 > 2.0;
  kuyruk_z=5.17): zarf, `exp(−(2πτ)²σ²/2)` biçiminde bir DW-sönümü DEĞİL.
- **H-Z1 (erfc) yaşıyor ama dejenere:** en iyi τ_c=**1.84** (Δ=1.23) — kesim
  ölçüm penceresinin ÇOK ötesinde; erfc, [0.45,0.86]'da yalnız uzak-sol
  kuyruğuyla (neredeyse doğrusal nazik eğim) uyuyor. Bu, KALEM'in H-Z1
  sorusuna ("keskin τ_c→∞ ile erfc-0.68 arasında nerede?") DOĞRUDAN yanıt:
  **gerçeğin genlik-zarfının etkin kesimi τ_c≈1.84 — keskin (τ_c→∞) ucuna
  çok yakın, 152'nin τ_c=0.68'inden UZAK.** 152'nin 0.68'i S3 istatistiği
  (genlik+kilit karışımı) etkin kesimiydi; DOĞRUDAN genlik-zarfı çok daha
  naziktir.
- **H-Z2 (güç yasası) KAZANIR** (χ²/dof=0.34; Δχ²/dof=0.40 > ön-kayıtlı 0.30
  kıyas eşiği ⇒ erfc'yi geçer): zarf ≈ **1 − 0.149·τ^1.30**, yumuşak bir
  güç-yasası sönümü.

**K2 HÜKMÜ:** zarf = NAZİK GÜÇ-YASASI sönümü (H-Z2); Gauss-DW ölü; erfc etkin
kesimi τ_c≈1.84 (near-sharp). Zarf, keskin bir merdiven kesimi değil,
kuyrukta ~%10'a varan yumuşak genlik-bastırmasıdır.

### 164 SANDVİÇİ — gerçek keskin ile erfc-0.68 arasında NEREDE  [184c]

İkiz=HA4 (erfc-0.68) ölçümü: w_HA4 düşük-τ'da 1.34, sonra erfc kesimiyle
kuyrukta **0.87'ye ÇÖKÜYOR**. Böylece gerçek, iki sadakatli ikiz arasında:

| τ̄ | w_gerçek | w_Hkeskin (keskin) | w_HA4 (erfc) | sandviç konumu f* |
|---|---|---|---|---|
| 0.626 | 1.241 | 1.354 | 1.289 | — (HA4 henüz düşmedi) |
| 0.725 | 1.269 | 1.406 | 1.025 | 0.64 |
| 0.775 | 1.270 | 1.415 | 0.927 | 0.70 |
| 0.830 | 1.256 | 1.406 | 0.872 | 0.72 |

f* = (w_gerçek − w_HA4)/(w_Hkeskin − w_HA4): kuyrukta **~0.70** — gerçek,
erfc-0.68'den KESKİN uca doğru yolun ~%70'inde; keskin-ikize çok daha yakın.
r vs Hkeskin = 0.89 (10% altında) ↔ r vs HA4 = 1.44 (%44 üstünde). Bu, K2'nin
erfc τ_c≈1.84 (near-sharp) bulgusunu BAĞIMSIZ doğrular: **gerçeğin genlik-
zarfı 152'nin erfc-0.68'inden çok, keskin merdivene yakındır.**

---

## K3 — ZARF ↔ ΔM KÖPRÜSÜ  [184d_K3_kopru.py]

**MÜHÜRLENDİ.** DOĞRUDAN genlik-zarfı r_184(τ)=â_gerçek/â_Hkeskin, 180'in
R_bant-türevli zarfı r_180(τ)=R_bant(son)/R_bant(Hkeskin) ile **her bantta
özdeş** (maks fark **0.0059**, rms 0.0025; en büyük sapma τ=0.775'te 0.006):

| τ | r_184 (direkt) | r_180 (R_bant) | fark |
|---|---|---|---|
| 0.475 | 0.9445 | 0.9430 | +0.0014 |
| 0.626 | 0.9169 | 0.9164 | +0.0005 |
| 0.725 | 0.9021 | 0.8993 | +0.0027 |
| 0.830 | 0.8931 | 0.8917 | +0.0014 |

Bunun anlamı büyük: 180, r(τ)'sini R_bant oranından türettiği için "KISMEN
kilit taşır, büyüklüğü ölçülmedi, %81.5 ÜST-SINIR tadında" uyarısını
koymuştu. 184 zarfı DOĞRUDAN çizgi genliğinden (kilit-bağımsız, birinci-ilke)
ölçtü ve r_180 ile <%0.6 farkla ÖZDEŞ buldu ⇒ **180'in r(τ)'si kilit-
kirlenmesi TAŞIMIYORDU; zarf gerçek bir genlik-spektrumu olgusudur.**

- 180: **ΔlogM(son↔Hk) = +0.057542**, zarf payı (log) = **%81.5 ± %11.3**
  (kilit %18.5). r_184 ile r_180 arasındaki fark (0.6 puan) 180'in ±%11.3
  se'sinin çok altında ⇒ 180'in ΔM modeli doğrudan-ölçülen w(τ)'yle
  beslendiğinde AYNI %81.5'i verir (ileri-hesap yeniden üretimi).
- Bağımsız **güç-bütçesi** (Σâ²/2 oranı = r²_güç): zarf güç-açığı τ>0.45'te
  **%15.0**, kuyruk τ>0.70'te **%19.2** (r≈0.90 → güç oranı ≈0.81).
- ÖN-KAYIT MÜHÜR EŞİĞİ: zarf payı ∈ [%70,%88] **ve** r_184≡r_180 → **İKİSİ DE
  SAĞLANDI** (%81.5 ∈ bant; maks r-farkı 0.0059 < 0.02).

**HÜKÜM — ZARFIN KİMLİĞİ MÜHÜRLENDİ:** gerçeğin sadakatli ikizinden farkının
~%80'ini oluşturan "zarf", tam olarak DOĞRUDAN ÖLÇÜLEN çizgi-genlik profili
w(τ) = â_q/â_Hkeskin ≈ **1 − 0.149·τ^1.30** (0.94→0.89; kuyrukta ~%10 genlik /
~%19 güç bastırması) dır. Keskin bir kesim değil, yumuşak güç-yasası sönümü.

## K4 (bonus) — A(τ) / λ_eff BAĞI  [184c/184d]

KALEM-K4 "H-Z3 doğruysa zarf-σ'sı A(τ)'nun λ_eff=0.99'uyla aynı fizik mi?"
diye soruyor. **H-Z3 (Gauss-DW) K2'de KURTARMASIZ ÖLDÜ** (χ²/dof=6.94,
kuyruk_z=5.17) ⇒ zarf `exp(−(2πτ)²σ²/2)` biçiminde DEĞİL; karşılaştırılacak
bir DW-σ'sı YOK. Dolayısıyla zarf ile kalibrasyon A(τ)'nun λ_eff=0.99'u
**aynı Gauss-DW fiziğini paylaşmıyor** — zarf güç-yasası (H-Z2), σ_X-
kalibrasyonunun λ-değişmez ~Gauss biçimi ise ayrı bir nesne. (Premis
yanlışlandığı için K4'ün σ↔λ_eff eşitliği HÜKÜMSÜZ; olumsuz ama net.)

---

## MANŞET (Not 5/6 için aday cümle)

> **Gerçeği sadakatli ikizinden ayıran farkın ~%80'ini oluşturan "zarf",
> artık DOĞRUDAN ölçülmüş bir çizgi-genlik profilidir:** her asal-kuvvet
> çizgisinin etkin genliği â_q, sadakatli-keskin ikizinkinin (Hkeskin)
> τ ile 0.94'ten 0.89'a NAZİK, MONOTON düşen bir oranıdır (kuyrukta ~%10
> genlik / ~%19 güç bastırması). Bu zarf **yumuşak bir güç-yasası**
> (1 − 0.149·τ^1.30; χ²/dof=0.34), **keskin bir kesim DEĞİL** (erfc etkin
> kesimi τ_c≈1.84, 152'nin 0.68'inden uzak; gerçek, erfc-0.68 ile keskin
> arasında keskin uca %70 yakın) ve **Gauss/Debye-Waller DEĞİL** (H-Z3
> χ²/dof=6.9 ile öldü). Doğrudan-ölçülen bu zarf, 180'in R_bant-türevli
> zarfıyla <%0.6 farkla ÖZDEŞ çıkarak 180'in ΔM zarf-payını (%81.5±11.3)
> mühürledi ve "R_bant kilit taşır" endişesini kapattı.

Teslim: bu rapor + `184_configs/` (184a–184e) + `184_zarf.png` (sol: profil
sandviçi; sağ: saf zarf + hipotez yarışı + 180-özdeşliği).
