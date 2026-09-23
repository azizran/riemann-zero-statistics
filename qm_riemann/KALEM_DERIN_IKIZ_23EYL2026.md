# KALEM — 189: DERİN İKİZ — ZARFIN KİMLİĞİ SINAVI
(23 Eylül 2026 — 188 ortak teftişinin açtığı soru; kullanıcı onaylı yön)

## Durum ve soru

188: iptal çekirdeği kinematik-genel (corr 0.999) ve DERİNLİK-EŞLİ ζ iki denizde
aynı (τ'≤1.00: gerçek 0.0758, ikiz 0.0772). Yani ikizle gerçeğin ζ farkı
(0.077 vs 0.33), ikizin merdiveninin kurgu gereği τ=1'de (q = t/2π) bitmesinden.
4. demirin manşeti "gerçeğin ikiz-fazlasının ~%80'i nazik bir zarf
w(τ) = 1 − 0.149·τ^1.30" idi. **Soru: bu zarf, ikizin merdiven kesiminin parmak
izi mi?** Sınav: ikizi DAHA DERİN merdivenle (τ≤1.10, τ≤1.20) aynı sadakatli
çözücüyle yeniden kur; zarf öngörülen hızla kapanıyor mu?

## ÖN-MÜHÜR — veri-öncesi sayısal öngörü (kaptan, 23 Eyl, 188 haritasından)

Kinematiği SABİT tutan aritmetik (187 kapanışı): m^kesik_Hk = m_Hk/(1−z_1.00),
m_HD = m^kesik_Hk·(1 − z_D), r_D = w_g/(w^öz_Hk + m_HD); z_D = −Re Σ_{τ'≤D} K_Hk
(188 ikiz haritası). Kaynak: scratchpad/185/K1_faktorler.npz + 188/harita_K_Hkeskin.npz.

| bant | r_1.00 (184) | z_1.00 | z_1.10 | z_1.20 | r_1.10 öng. | r_1.20 öng. | kapanan % (1.10 / 1.20) |
|---|---|---|---|---|---|---|---|
| 0.45-0.50 | 0.9445 | 0.0920 | 0.1884 | 0.2480 | 0.9802 | 1.0037 | 64 / 107 |
| 0.50-0.55 | 0.9351 | 0.0841 | 0.1792 | 0.2421 | 0.9747 | 1.0028 | 61 / 104 |
| 0.55-0.60 | 0.9253 | 0.0821 | 0.1808 | 0.2418 | 0.9708 | 1.0013 | 61 / 102 |
| 0.60-0.65 | 0.9169 | 0.0803 | 0.1788 | 0.2403 | 0.9666 | 1.0003 | 60 / 100 |
| 0.65-0.70 | 0.9096 | 0.0705 | 0.1694 | 0.2314 | 0.9624 | 0.9987 | 58 / 99 |
| 0.70-0.75 | 0.9021 | 0.0678 | 0.1752 | 0.2323 | 0.9628 | 0.9985 | 62 / 99 |
| 0.75-0.80 | 0.8978 | 0.0710 | 0.1886 | 0.2619 | 0.9677 | 1.0171 | 68 / 117 |
| 0.80-0.86 | 0.8931 | 0.0637 | 0.1810 | 0.2629 | 0.9638 | 1.0203 | 66 / 119 |

Okuma: sabit kinematikle zarf 1.10'da ~%60, 1.20'de ~%100 kapanıyor. AMA 1.20'deki
"tam kapanış" kısmen tesadüftür: kalan derinlik açığını (gerçek τ'>1.20'yi de
duyar) kinematik çarpan F_KİN = w^öz_g/w^öz_Hk > 1 (1.02→1.14; gerçeğin aralık
saçılımı ikizden dar: σ_ε 0.4092 vs 0.4313) telafi ediyor. Derin ikizin kinematiği
de değişecek — bu yüzden birincil öngörü YÖN ve KABA HIZ, sayısal kıyas ikincil.

**Mekanizma öngörüsü (kalem cebiri, 188):** τ'≈1 çizgilerinin çekirdeği
sin(πτ'(1+ds)) ≈ −π·ds — aralık sapmasına karşı NEGATİF GERİ-BESLEME. Derin
merdiven bu geri-beslemeyi ekler ⇒ derin ikizin aralık saçılımı σ_ε AZALMALI,
gerçeğin 0.4092'sine doğru.

## Hipotezler (ön-kayıtla donar; ölümler kurtarmasız)

- **H-189a — ZARF = MERDİVEN DERİNLİĞİ (birincil):** HAVUZ'da ve ≥7/8 bantta
  r_1.00 < r_1.10 ve r_1.10 ≤ r_1.20 + 2σ (zarf derinlikle kapanır); HAVUZ
  kapanan payı f_1.10 = (r_1.10 − r_1.00)/(1 − r_1.00) ≥ 0.30. ÖLÜM: f_1.10 < 0.15
  ya da r HAVUZ'da azalırsa.
- **H-189b — SABİT-KİNEMATİK ARİTMETİĞİ (ikincil, nicel):** r_1.10 ve r_1.20
  bant başına yukarıdaki öngörünün max(2σ, 0.01) içinde. Tutmazsa sapma
  kinematik değişimle (w^öz_HD) açıklanıyor mu — ayrışım tablosu (KAYIT).
- **H-189c — GERİ-BESLEME KİNEMATİĞİ:** σ_ε(Hderin110) < σ_ε(Hkeskin) = 0.4313 ve
  σ_ε(Hderin120) < σ_ε(Hderin110); ikisi de gerçeğin 0.4092'sine doğru (aşarsa
  kayıt). ÖLÜM: σ_ε derinlikle artarsa.
- **H-189d — HARİTA GEÇERLİLİĞİ (kapı-niteliğinde):** derin ikizin KENDİ 187c-ζ'si
  harita öngörüsüne uyar: HAVUZ |ζ_HD − z_D|/z_D ≤ %10 (z_1.10 ≈ 0.18,
  z_1.20 ≈ 0.24; bant bazında da raporlanır).
- **K3 — ZARFIN KİMLİĞİ (KAYIT):** üç derinlikten (1.00, 1.10, 1.20) r_D(τ)
  profilleri; her birine 1 − c·τ^α → (c, α) derinlikle nasıl sönüyor; D→∞
  kestirimi (kinematik ve derinlik payı ayrı); 4. demirin %80 zarfının ne kadarı
  derinlik, ne kadarı gerçek kinematik fark. "Türetilen vs ölçülen" açık cümleyle.

## Kapılar

- **K0 — ÖN-KAYIT (sha+damga, inşadan ÖNCE):** yukarıdaki öngörü tablosu AYNEN,
  hipotez eşikleri, ikiz adları (Hderin110, Hderin120), ölçüm zinciri.
- **K1 — İNŞA:** 164_insa.py'yi DÜZENLEMEDEN, sarmalayıcıyla `merdiven(...,
  tau_ust=1.10 / 1.20)`; diğer her şey Hkeskin'in aynısı (sadakatli, c=−½, h=0.015).
  Kapılar: maks|F| ≤ 1e-8, sıralılık TAM, L = 12.02959324. IZGARA SINAVI: derin
  merdivende en kısa dalga boyu küçülür — 5000 sıfırlık dilimde h=0.015 vs
  h=0.0075 ilk-kök farkı ≤ 5/5000 (aşarsa h=0.0075 ile inşa, raporla). Önce 1.10
  (~30 dk), sonra 1.20 (~95 dk); nohup + ≤30 sn yoklama.
- **K2 — ZİNCİR:** her derin ikiz için 184b2 → 185b → 186b → 187c (sarmalayıcı
  sürücülerle, içerik değişikliği yok); w-bantları, r_D, w^öz, σ_ε, ζ_HD.
- **K3 — HÜKÜM + KİMLİK KAYDI.**

Ölümler kurtarmasız; TEK DALGA; ajan git'e DOKUNMAZ; sonuç ORTAK TEFTİŞE.
Veri: `qm_riemann/scratchpad/` (köprü: eski geçici yol → kalıcı klasör).
