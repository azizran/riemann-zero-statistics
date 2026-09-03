# 167 — Ç4 ORTAK-KAZANCI c: ABLASYON HARİTASI (3 Eylül 2026)

Durum: KISMİ TESLİM — ölçümler ve harita tamam (T1, T2-kesim,
T2-pencere, T4); λ-ablasyonu (L085/L115) İNŞA EDİLDİ ama ölçülmedi;
T3 analitik rezonans-integrali YAPILMADI. Ölçümler Opus tayfasının
(167_configs zinciri, iki bekçi-kesintisi); derleme ve hüküm Fable
(167_analiz.py yeniden koşusu, tüm C_*.json üstünde). Ham çıktılar
scratchpad/167.

## H1 HÜKMÜ: c SAF SABİT DEĞİL — sonlu-boyut/merdiven FONKSİYONELİ

1. **Kesim ekseni:** c(kesim)/c(keskin≤1.0), tüm üyelerde tutarlı:
   keskin≤0.7 → 0.86-0.90; keskin≤0.9 → 0.90-0.93; erfc-0.68 →
   0.87-0.89. c, merdiven derinleştikçe +%10-14 BÜYÜYOR.
2. **Pencere ekseni:** c₀(T) = 0.397 (T) → 0.344-0.367 (T/2) →
   0.283-0.317 (T/4). c, pencere kısaldıkça KÜÇÜLÜYOR (−%10 @T/2,
   −%25 @T/4) — κ-tepe genişliği/rezonans-sayımı etkisi.
3. **Sonuç:** iki eksene de bağımlılık ⇒ **saf geometrik sabit
   adayları (1/√π, 9/16) ÖLÜR.** c, (pencere × merdiven-yoğunluğu)
   rezonans-sayım fonksiyonelidir; T3 türetimi bu iki ölçülü ölçek
   yasasını üretmek zorunda (gelecek kalem için iki sert kısıt).
4. λ (genlik) ekseni: ölçülmedi (açık).

## T1 — hassas taban (üç gaz, jackknife)

W_X üyesinde ÇARPICI: son = 0.4053, Hkeskin = 0.4051 (%0.05!) —
gerçek gaz bu fonksiyonelde tam-derinlikli keskin merdiven gibi;
HA4 = 0.3345 (etkin kesim sığ). 166'nın gerçek↔sentetik %0.06
taşınması bağımsız doğrulandı.

## T4 — kapanış (tek küresel c, 20 hüküm bandı)

W_amp·W_X: c = 0.5647 → **13/20 (13/15 sağlıklı)** — 165'in bant-bant
u2'sini (12/20) TEK sayıyla geçiyor; 166'nın hükmü doğrulandı.
Diğer üyeler 11-12/20.

## Bant-şekli ve Gram notları

γ_eff (kalib ∝ e^{−γτ²}) gazlar arası 0.75-2.21 — tek DW üyesi
evrensel değil (166'nın "üyeler ayırt edilemez" kaydıyla tutarlı).
Gram-özdeşliği denemesi c = θ/(Λ_EΛ_X²⟨W⟩) sabitlenmiyor (0.51-0.76)
— kapanmadı, kayıt.

## Açık borçlar

λ-ablasyonu ölçümü; T3 analitik rezonans-integrali (artık iki ölçülü
kısıtla: c(T) ve c(kesim) ölçek yasaları); HA4−%17.5'in kesim-ekseni
okuması (kısmen açıklandı: HA4 ≈ K070 sınıfı).
