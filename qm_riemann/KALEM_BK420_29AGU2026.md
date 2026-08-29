# KALEM OTURUMU — BK (4.20), x = ½ açılımı (29 Ağu 2026)

Hedef: Berry–Keating (1999, SIAM Rev 41) denk. (4.20)'nin x=½ değerini
kendi dilimize çevirmek; R_nn = −2cos(πτ) yasasının ŞEKLİ ile KATSAYISININ
literatür sınırını çizmek.

## 1. Çeviri: (4.20) bizim birimlerde

BK (4.20), unfolded çift-korelasyonun asal (evrensel-olmayan) düzeltmesi:

  R_c¹(x) = [1/(2(π⟨d⟩)²)] Σ_{p,m} (log²p/p^m) cos(x·m·log p/⟨d⟩) − düz

⟨d⟩ = L/2π ve a_q = 1/(π·m·p^{m/2}) ⇒ log²p/p^m = π²ω_q²a_q²,
x·ω_q/⟨d⟩ = 2πx·τ_q. Terim terim:

  **R_c¹(x) = Σ_q 2π² τ_q² a_q² · cos(2πx·τ_q) − düz**

Ağırlık 2π²τ²a_q² — köşegen (bağımsız-faz) katman. Bu, kesin özdeşliğin
ikinci momentinin ta kendisi: aynı cebir Berry-88'i (rapor §3.1) ve
BLS (2001) denk. (44)'ü verir. Yani (4.20) = bizim dalga-resminin
köşegen izdüşümü. ✓

## 2. x = ½ sözlüğü: iki türetim, tek kinematik

- **Onların yolu:** x=½ (yarım-aralık ayrımı) ⇒ cos(2π·½·τ) = **cos(πτ)**.
- **Bizim yolumuz (134/135):** bağ orta-noktası mm'den okuma; iki gap
  merkezi mm ∓ ḡ/2'de ⇒ faz farkı ω_q·ḡ/2 = **πτ_q** ⇒ cos(πτ).

Aynı kinematik: **yarım-gap fazı.** R_nn'in cos(πτ) ŞEKLİ artık çift
çapalı — bizim çarpanlara ayırma + BK (4.20). Şekil türetilmiş sayılır.

## 3. Katsayı muhasebesi (oturumun teorem cümlesi)

Bizim ölçü nesnemiz KOŞULLU: η_nη_{n+1}'in, q-dalgasının fazına kilitli
modülasyonu / (2c₁A1_q). Köşegen dünyada (fazlar bağımsız; BK99 s.248):

  R_nn^köşegen ≡ 0.

Üstelik tam-taban regresyonumuz birinci-mertebe dalga içeriğini (köşegen
katmanın taşıdığı her şeyi) η'dan ZATEN söküyor. Dolayısıyla:

  **ŞEKİL köşegen-kinematik (yarım-gap fazı; 4.20 türetir);
  KATSAYI (−2) tamamen KÖŞEGEN-DIŞI (4.20 İLKE OLARAK üretemez).**

−2'ye giden iki makine kalıyor: (i) BK95 Hardy–Littlewood köşegen-dışı
hesabı (q₁ = q·q₂ çarpımsal rezonansları — NOT: 133 statik çift-vuru
toplamını ≈0 ölçtü; canlı mekanizma konum-fazı modülasyonu, 134/135),
(ii) bizim 135 öz-tutarlı çekirdeği. Ayrışım defteri:

  R_nn(τ) = 0 (köşegen) + κ_ad(τ)·cosπτ (adyabatik kinematik, ölçülü
  +0.92..+1.35) + K(τ)·cosπτ (köşegen-dışı çekirdek, ölçülü −2.87..−3.25)

## 4. Sıradaki kalem hedefi: (4.23) resurgence terimi

(4.23) kapalı formunda −∂²_ξ Re log ζ(1−iξ) var; log ζ(1−iξ) =
Σ_q (p^{−m}/m) e^{iξω_q} ⇒ bu terim aynı asal merdiveninin toplanmış
hali + ζ'nin 1-doğrusu yakınındaki yapısı ("resurgence": alçak sıfırlar
yüksek sıfırların korelasyonunda rezonans). K(τ)'nun yavaş kayışını
(−2.87→−3.25) üretebilecek tek analitik yapı adayı. → sonraki oturum.

## 5. Oturumun ölçüm ayağı (139)

Rapor testi #2: Berry-88'i zeros6'da SABİT-ARALIK nesnesiyle sına —
V(ḡ) = Var[S(t+ḡ)−S(t)] (t düzgün ızgarada) ↔ Σ_q 2a_q²sin²(πτ_q) + 1/π²
(τ≤1 keskin kesim). Karşıt nesne: Var(ds_n) (SIFIRLARDA örneklenmiş).
Beklenti: sabit-aralık ~0.34'ü vurur (formül doğrulanır, veri hattı
bağımsız onay); sıfırlarda ~0.17 (koşullama/örnekleme farkı — Berry'nin
kendi uyarısı). Fark, ekran ailemizin (D_fiz) frekans-toplamlı hali;
τ-çözünürlüğü ayrı iş.
