# KALEM — 186: M(τ) İLETİMİNİN TÜRETİMİ (zarf kendi karışımını doğurur mu?)
(9 Eylül 2026 — 185'in bıraktığı adres: karışım-iletim profili M(τ))

## Durum

185: zarfın taşıyıcısı tek-çizgi kinematiği DEĞİL; komşu-çizgi karışım
genliğinin bastırılması: M(τ) = m_g/m_Hk = 0.796±.004 → 0.721±.021
(m = w_ölçülü − w^öz). "Neden α=1.30?" sorusu M'nin türetimine indirgendi.
Bu kalem M'yi üç güzergâhtan kovalar. NOT: sonuç raporu commit'ten önce
KAPTAN+KULLANICI ortak teftişine sunulur (yeni usul).

## Kalem cebiri — üç güzergâh

**G1 — ÖZ-TUTARLILIK (yeniden-ağırlıklı merdiven):** karışım m_q, DİĞER
çizgilerin q'ya taşan katkısıdır. Hipotez: gerçek denizde her çizgi q'
kendi zarfını (r(τ')) taşıdığından, q'daki karışım da bu bastırılmış
genliklerin katlanmasıdır — zarf kendi karışımını doğurur (sabit-nokta).
Sınama parametresiz: yeniden-ağırlıklı merdiven
  ds^ya_n = Σ_q' 2a_q'·ρ(τ_q')·sin(ω_q' g_n/2)·cos(ω_q' m_n),
ρ(τ') = 184'ün ölçülü yasası (1−0.149τ'^1.30), GERÇEK kinematikle; izdüşüm
− yeniden-ağırlıklı öz-terim = m^pred; M_pred = m^pred/m_Hk. Kontrol:
ρ≡1 aynı makinede 185-B'nin ölümünü (r_pred>1) yeniden vermeli.
Dürüstlük şerhi: bu bir DIŞ-İLKE türetimi değil, sabit-nokta tutarlılığı
— çekirdek (hangi çizgi kime taşar) bağımsız kinematik, ρ ölçülü zarf.
Kapanırsa α=1.30 sorusu sabit-nokta denkleminin çözümüne iner (rapora
doğrusallaştırılmış sabit-nokta kalem-bölümü: güç-yasası öz-üretiyor mu?).
Varsayım kaydı: çizgiler faz-sadık, yalnız genlik yeniden-ağırlıklı —
kırılırsa kilidin karışımdaki rolü öğrenilir (o da bulgudur).

**G2 — v-KANALININ YÜKSEK-τ İLK ÖLÇÜMÜ + KISIT KÖPRÜSÜ:** Not 2/3'ün
iki-kanal kısıtı √w = 1.017 − 0.884·v yalnız τ<0.47'de ölçüldü (v doyumu
~0.55-0.60 @ τ≈0.4). v(τ)'yi τ∈[0.45,0.86]'da İLK KEZ ölç (aralık
serisinin çizgi-başına iletimi; her iki denizde AYNI kestirimci —
konvansiyon riskine karşı ORAN lehçesi birincil: v_g/v_Hk). Köprü:
kısıt her iki denizde geçerliyse M(τ) = w_g/w_Hk = [(1.017−0.884v_g)/
(1.017−0.884v_Hk)]² — parametresiz M_pred. (Kısıtın yüksek-τ'ya uzanması
kendisi de sınanan bir yasadır; ölürse ölümü raporlanır.)

**G3 — BİÇİM DEFTERİ (türetim değil, kısıt):** M(τ)'ye aday biçimler
(doğrusal; 1−k·τ^β; v-afin) jackknife'la — güzergâhları sınırlamak için
dürüst ampirik kayıt. (Ön-bakış: 1−M ≈ τ^0.56-ish — √τ kokusu; kayda.)

## Kapılar (186)

- **K0 — KURAL-ÖNCE ÖN-KAYIT (sha+damga):** m/M tanımları 185'ten AYNEN
  (OZ_*.npz önbellekleri); yeniden-ağırlıklı merdiven kestirimcisi; v̂
  kestirimcisi (iki denizde özdeş); bant ızgarası 184-185 AYNEN; 8-blok
  jackknife; ölüm eşikleri: G1/G2 M_pred bant-χ²/dof ≤ 2; kapanış eşiği.
- **K1 — G1 SINAVI:** M_pred(öz-tutarlılık) vs ölçülü M, 8 bant; ρ≡1
  kontrolü (185-B ölümünü yeniden üretmeli — makine mührü).
- **K2 — G2 SINAVI:** v_g(τ), v_Hk(τ) ilk yüksek-τ ölçümü (+ jackknife);
  kısıt-köprüsü M_pred vs ölçülü M. v-profilinin kendisi bağımsız bulgu
  olarak raporlanır (Not 2'nin açık sorusuna veri).
- **K3 — KAPANIŞ:** kazanan güzergâhın M_pred'i zincire takılır:
  r_pred(τ) = (w_g^öz + M_pred·m_Hk)/(w_Hk^öz + m_Hk) → 184 bantlarına
  χ²/dof + türetilmiş (c,α) vs (0.149, 1.30). G1 kazanırsa sabit-nokta
  analizi (kalem-bölümü) rapora.
- **K4 (bonus) — HA4:** kazanan güzergâh HA4'ün karışım genliğini de
  vermeli (185'in bulduğu ~%8-10 karışım-eşitsizliği dahil — açıklarsa
  185-K4 çatlağı da kapanır).

Ölümler kurtarmasız; TEK DALGA; koşular nohup+kısa yoklama (5 dk ön-plan
YASAK); ajan git'e DOKUNMAZ; sonuç ORTAK TEFTİŞE sunulur, commit sonra.

## Veri/makine

scratchpad/155: eta_son (ds, mid), z_Hkeskin.npy, z_HA4.npy;
scratchpad/185: OZ_{gercek,Hkeskin,HA4}.npz, K1_faktorler.npz (m, M
bantları burada); scratchpad/184: K1_*.npz, ONKAYIT_184.json.
Yeniden-ağırlıklı merdiven O(N×Q) ≈ 1e9 — parça-parça vektörle, dakikalar.
