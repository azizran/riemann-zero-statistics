# KALEM — 187: PENCERE-ÖTESİ FAZ-İPTAL DEFTERİ (α=1.30'un güncel adresi)
(10 Eylül 2026 — 186'nın bonus teşhisinden doğan kalem; ortak-teftiş usulüyle)

## Durum

186: gerçek deniz, pencere-içi (τ'≤0.86) merdiven karışımının ~üçte birini
pencere-ötesi içerikle FAZ-UYUMLU iptal ediyor (m^kesik/m_ölç = 1.46→1.58;
ikizde 1.04→1.10); HA4 yön-teyidi (kuyruğu ölünce karışım ARTAR). M(τ) =
1−0.277·τ^0.42 ayrı ilkel. Bu kalem iptalin DEFTERİNİ çıkarır: kim taşıyor
(hangi τ'-katmanı), ne kadar derine iniyor, açısı gerçekten 180° mi, ve
defter HA4'ün %26'sını kapatıp M'nin yapısını veriyor mu.

## Ölçülecek nesneler

- **KATMAN DEFTERİ:** nominal merdiven GERÇEK kinematikte katman katman
  derinleştirilir: τ_c ızgarası [0.86, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15,
  1.20] (τ_c=1.20 → q≤1.8e6, ~135k çizgi — hesap sınırı). Her katman
  artımı ayrı seri olarak kaydedilir (kontrol noktası!); bant-izdüşümüyle
  Δm(katman; bant) defteri. m^ya(τ_c)'nin m_ölç'e yaklaşması = yakınsama
  defteri; korr(ds^ya(τ_c), ds) yan sütunu (186'da 0.86'da 0.953-0.957 —
  yükselmeli). Not: özdeşlik gereği tam-derinlik limiti m_ölç'ün kendisi;
  türetim içeriği LİMİT değil, YOLUN YAPISIDIR (hangi katman ne taşıyor).
- **İPTAL KATSAYISI ve AÇISI:** Δ_n = ds_n − ds^kesik_n (gerçeğin
  pencere-ötesi içeriği); her bantta kompleks oran
  ζ(bant) = izdüşüm(Δ)/izdüşüm(karışım^kesik) — modülü iptal payı
  (beklenen ~0.31-0.37), AÇISI faz-uyum sınavı. Jackknife'la açı hatası.
- **İKİZ SAĞIRLIĞI (kontrol):** aynı katman defteri Hkeskin kinematiğinde
  (kontrol ızgarası [0.90, 1.00, 1.10]): ikiz o çizgileri hiç "bilmeden"
  kuruldu (kendi merdiveni τ≤1.00) — τ'>1.00 katmanları ikiz kinematikte
  ≈0 vermeli. HA4 kinematiği de (K4 için, aynı kontrol ızgarası).

## Hipotezler (ön-kayıtla donar; ölümler kurtarmasız)

- **H-F1 (YAKINSAMA ÇATALI):** gerçek kinematikte m^ya(τ_c), τ_c ile
  m_ölç'e MONOTON yaklaşır ve τ_c=1.20'ye dek açığın ≥ %50'si kapanır.
  İki okuma ÖN-KAYITLI: (a) katman artımları art arda iki katmanda
  |Δm| < 2·se → "SIĞ İPTAL — yakınsadı" mührü; (b) 1.20'de açık hâlâ
  > %50 → "DERİN KUYRUK" bulgusu (bu da ölçümdür, kurtarma değil).
- **H-F2 (FAZ-UYUM):** ζ açısı ön-kayıtlı pencerede 180°±15° (bant
  başına) → "saf yıkıcı" mührü; dışındaysa döndürülmüş faz-örgüsü KAYDI
  (ölüm değil, karakter tayini — açı ne çıkarsa deftere).
- **H-F3 (İKİZ SAĞIRLIĞI):** ikiz kinematikte τ'>1.00 katman artımları
  |Δm| ≤ 2·se. Bu kontrol ÖLÜRSE (ikiz de derin çizgilere kulak veriyorsa)
  186'nın "iptal gerçeğe özgü" okuması düşer — açıkça raporlanır.
- **H-F4 (HA4 %26 KAPANIŞI):** defter aritmetiği (HA4 = erfc-ağırlıklı
  pencere-içi karışım, iptalsiz; Hkeskin = tam pencere-içi + kendi sığ
  kuyruğu) M_HA4(bant) = 1.264→1.019'u versin: bant-χ²/dof ≤ 2 →
  186-K4 çatlağı KAPANIR. Geçemezse ölüm dürüst.
- **K3 (yapı):** iptal-yoğunluğu ι(τ'; bant) profili — iptal kesim
  hemen-ötesinde mi yoğun (0.86-1.1) yoksa yayvan mı; hedef-bağımlılığı
  C(τ) = m^kesik/m_ölç − 1'in τ ile büyümesi (1.46→1.58) katman
  defterinden yeniden kuruluyor mu (M_rekon tutarlılık mührü); C(τ)'ye
  biçim defteri (türetim değil, kayıt). "Türetilen vs ölçülen" ayrımı
  raporda açık cümleyle.

## Kapılar (187)

- **K0 — KURAL-ÖNCE ÖN-KAYIT (sha+damga):** katman ızgaraları, kestirimci
  tanımları (katman-artımı izdüşümü, ζ ve açısı, bant/8-blok jackknife
  184-186 AYNEN), H-F1 çatalının iki okuması, H-F2 açı penceresi,
  H-F3/H-F4 eşikleri — sayı doğmadan donar.
- **K1 — KATMAN DEFTERİ (gerçek):** ana koşu; katman başına kontrol
  noktası; yakınsama + korr yan sütunu.
- **K2 — ζ DEFTERİ:** iptal payı ve açısı bant bant; ikizde aynı ölçüm
  (küçük/gürültülü beklenir).
- **K3 — YAPI/REKONSTRÜKSİYON:** ι(τ') profili; C(τ) yeniden-kurulumu;
  biçim kaydı.
- **K4 — HA4:** defter aritmetiğiyle M_HA4 öngörüsü vs 186 ölçüsü.

## Hesap notları

Ana yük: gerçek kinematikte τ_c=1.20'ye merdiven ≈ 135k çizgi × 300k
nokta ≈ 8e10 trig — parça-parça vektörle, KATMAN BAŞINA kaydet (yarıda
kesilirse kaldığı katmandan devam), nohup + ≤30 sn yoklama; ~20-40 dk
beklenir. İkiz/HA4 kontrol ızgaraları daha ucuz. İsteğe bağlı ikincil:
gerçek için (1.20, 1.40] katmanları 100k alt-örneklemle (ikincil damgalı,
ayrı se'yle). Veri: scratchpad/155 (eta_son, z_Hkeskin, z_HA4),
scratchpad/186 (G1_proj_*.npz, ONKAYIT_186), scratchpad/184-185 (K1/OZ).
Makine: 186b'nin merdiven+izdüşüm çekirdeği AYNEN.

Ölümler kurtarmasız; TEK DALGA; 5 dk ön-plan YASAK; ajan git'e DOKUNMAZ;
sonuç ORTAK TEFTİŞE gelir, commit sonra.
