# KALEM — 185: ZARFIN TÜRETİMİ (neden 1 − 0.149·τ^1.30?)
(9 Eylül 2026 — dördüncü demirin açık defterindeki en tatlı kalem)

## Durum

184 zarfın YÜZÜNÜ ölçtü: r(τ) = â_gerçek/â_Hkeskin = 1 − 0.149·τ^1.30
(χ²/dof=0.34; 8 bant, τ∈[0.475,0.83]; kuyrukta ~%10 genlik bastırması).
Ama "neden bu biçim?" açık. Bu kalem türetim kalemidir: zarfı ölçmek değil,
KESİN ÖZDEŞLİKTEN hesaplayıp ölçülmüş bantlara parametresiz vurmak.

## Kalem cebiri (türetimin omurgası)

Özdeşlik kesin: ds_n = Σ_Q 2a_Q sin(ω_Q g_n/2)·cos(ω_Q m_n). İzdüşüm
â_q = |2⟨ds·e^{−iω_q m}⟩|. Tek-çizgi ÖZ-İZDÜŞÜMÜ tanımla (Taylor'suz, kesin):

  â_q^öz := |2⟨ 2a_q sin(ω_q g_n/2) · cos(ω_q m_n) · e^{−iω_q m_n} ⟩|

— yani "çizgi q'nun kendi terimi, GERÇEK aralık/orta-nokta kinematiğiyle,
tam beklenen-değer". Bununla üç-faktörlü KESİN ayrışım:

  r(τ) = [â_g/â_g^öz] × [â_g^öz/â_Hk^öz] × [â_Hk^öz/â_Hk]
          (F_ANOMALİ)     (F_KİNEMATİK)      (F_KOMŞU⁻¹)

- **F_KİNEMATİK:** aynı τ'da gerçek ile ikizin aralık/orta-nokta
  istatistiği farkı. Yorum için Taylor (yalnız yorum; πτσ_ε ~ 1 olduğundan
  hesap DAİMA kesin beklenen-değerle):
  - T2 (titreşim): −(πτ)²Δσ_ε²/2 — τ² biçimli, aralık-saçılımı farkından.
  - T1 (uyumlu geri-besleme): (πτ)·cot(πτ) × Ĝ_q — aralık serisinin çizgi-q
    frekansındaki UYUMLU içeriği (Ĝ_q = ⟨(g_n−ḡ)e^{−iω_q m_n}⟩ tabanlı).
    Dikkat: iskelet, κ çizgi-yasasının kinematiğiyle AYNI (143/165:
    κ = −πAτ·cos(πτ); orta-nokta kinematiği cos(πτ)). τ>½'de cot(πτ)<0 —
    tam bizim pencerede işaret değiştiriyor; zarfın penceresi tesadüf değilse
    buradan görünecek.
- **F_ANOMALİ:** gerçeğin ds'inin nominal tek-çizgi teriminden GERÇEK sapması
  — Not 2'nin w-kanalının (genlik iletimi) yüksek-τ devamının bizim dildeki
  karşılığı. Not 2 w'yu τ<0.47'de ölçtü (w₂=0.91 @ τ=0.028; doğrusal kılavuz
  0.982−2.598τ; v-kanalı τ≈0.4'te 0.55-0.60'a doyuyor). Bizim pencere
  τ∈[0.45,0.86] = o eğrilerin KEŞFEDİLMEMİŞ devamı. Köprü buradan kurulur.
- **F_KOMŞU:** komşu-çizgi payı (aynı apsiste diğer çizgilerin izdüşüme
  taşması). Gerçek ve ikiz apsisleri ÖZDEŞ ⇒ oranda büyük ölçüde sadeleşmeli.

Beklenen resim (önceki mühürlerle tutarlılık): w>1 her iki denizde
(gerçek 1.25, Hkeskin 1.39) ⇒ T1-tipi uyumlu güçlendirme BÜYÜK ve ortak
(kilit/H-G2'nin yeniden-kilitlemesi); zarf = bu büyük ortak terimin
gerçek−ikiz FARKI + F_ANOMALİ. Hangisi baskın — kalemin ana sorusu bu.

## Hipotezler (ön-kayıtla donar; ölümler kurtarmasız)

- **H-W0 (L-akışı sıfır-adayı):** pencere içi L kayması (±%3) zarfı
  açıklayamaz; katkı sınırı hesaplanır, <%1 beklenir → kapatılır.
- **H-W1 (ÖZ-MUHASEBE):** üç faktör, 184'ün 8 bandını PARAMETRESİZ yeniden
  üretir. Ölüm: bant-χ²/dof > 2 → zarf özdeşlik-makinesinin dışında
  (bu ölüm bile manşet olur; kurtarma yok).
- **H-W2 (BASKIN TERİM + KÖPRÜ):** zarfı taşıyan faktör teşhis edilir.
  T1/F_KİNEMATİK baskınsa: Ĝ(τ) bağımsız ölçülür (yalnız aralık+orta-nokta
  serisinden) ve kinematik çarpanla r_pred(τ) İLERİ-HESAP edilir — serbest
  sabit YOK. F_ANOMALİ baskınsa: Not 2 w-kanalının yüksek-τ devamı doğrudan
  ölçülmüş olur (Not 6'ya köprü cümlesi). Ölüm: baskın-faktör ileri-hesabı
  bant hatalarının 2σ dışında.
- **H-W3 (kuvvet-açığı):** script-100'ün mekanik ε(q,k) yasası (f≈τ^3.3)
  ileri-hesap; α=3.3 ≠ 1.30 ve ölçek küçük beklenir → dürüst kayıt/ölüm.

## Kapılar (185)

- **K0 — KURAL-ÖNCE ÖN-KAYIT (sha+damga):** â^öz/Ĝ/faktör tanımları, bant
  ızgarası (184'ünkü AYNEN), jackknife (8 blok), ölüm eşikleri, H-W0..W3
  formları — sayılar doğmadan donar.
- **K1 — ÖZ-MUHASEBE:** â^öz gerçek+Hkeskin (kesin beklenen-değer; yalnız
  z-dizileri ve eta önbelleğinden ds); üç faktörün bant defteri; H-W1 hükmü.
- **K2 — ANATOMİ:** faktörlerin τ-profilleri; baskın terimin kimliği;
  Taylor-yorumu (T1/T2) yalnız teşhis için.
- **K3 — KAPANIŞ:** baskın terimin ileri-hesabı 1−0.149·τ^1.30'u bant
  düzeyinde veriyor mu (χ²/dof ≤ 2, parametresiz)? Ek: ileri-hesap eğrisine
  aynı güç-yasası fiti → "türetilmiş (c, α)" 184'ün (0.149, 1.30)'uyla yan
  yana raporlanır.
- **K4 (bonus, örneklem-dışı) — HA4 SINAVI:** AYNI öz-muhasebe, erfc-ikiz
  HA4'ün ölçülü profilini (kuyrukta 0.87'ye çöküş, K1_HA4.npz) parametresiz
  yeniden üretmeli. İkinci deniz, sıfır yeni ayar.

Ölümler kurtarmasız; TEK DALGA; ölçüm akışları ELLE/nohup+kısa yoklama
(183 dersi); ajan git'e DOKUNMAZ.

## Veri/makine (hepsi hazır; pahalı regresyon YOK)

scratchpad/155: eta_son (ds, mid), z_Hkeskin.npy, z_HA4.npy;
scratchpad/184: K1_gercek/K1_Hkeskin/K1_HA4.npz (ölçülü â bantları),
ONKAYIT_184.json (bant ızgarası). Makine: 184b_K1_zarf.py (izdüşüm),
184b2_zarf_ikiz.py (ikiz ds'i nominal merdivenden kurma). Tüm yeni hesaplar
vektörel ve saniyeler mertebesinde.
