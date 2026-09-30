# 203 — Oran sanısından log|ζ(½+it)|'nın sonlu-yükseklik kümülantları: b0 bilmecesi çözüldü
(28 Eylül 2026 gece. Tahminler serbest parametresiz; veriyle karşılaştırmadan ÖNCE commit'lendi:
κ₂ 182190f, κ₃ bf214f0. Veri: mühürlü 200-A/C b0 penceleri. Kod: `203_configs/`.)

## Yöntem (koşullu: CFKRS/Conrey–Farmer–Zirnbauer oran sanısı)
log ζ(s) = −∫_0^∞ (ζ′/ζ)(s+a) da; κ₂ = ½∫∫E[X X̄], κ₃ = −(3/4)∫∫∫E[X X X̄] (X = ζ′/ζ kaydırılmış).
 - Çift ortalaması: G(x) = (ζ′/ζ)′(1+x) − Σ_p log²p/(p^{1+x}−1)² + e^{−Lx}ζ(1+x)ζ(1−x)A(x).
 - Üçlü ortalama (kendi türetimimiz): kimlik terimi (bağımsız-faz Euler modeli) + iki takas terimi
   F(x_i)·K(α_j−α_i, x_i); F(x) = e^{−Lx}ζ(1+x)ζ(1−x)A(x),
   K(u,x) = −Σ_p log p·q/(1−q)·(1−1/p)(1−p^{−x})/(1−2/p+p^{−1−x}), q = p^{−1−u} (analitik devamıyla).
 - Kimlik terimi TAM OLARAK Fazzari–Gerspach'ın c_P'sini verir (= κ₃^{arit} = 0.233653); takas terimleri
   L → ∞'da −π²/4'e (FG'nin c_Z'si) gider (L = 60: −2.2072, L = 120: −2.2186; limit −2.2337, (log L)/L ile).
 - Tuzak (kayıt): takas terimlerini ayırıp esas-değer almak koşullu yakınsak ve tam −π²/4 kaybettiriyor;
   birleşik çekirdek Φ(x₁,x₂)·min(x₁,x₂) ile CUE'de κ₃(N) N = 1, 2, 5, 10 için 5 basamak kesin.
 - Denetimler: ızgara 16→24 düğüm 6·10⁻⁴; kesme Wc 0.5↔0.75 6·10⁻⁵.
 - Teftiş (Sonnet, `203_configs/203_TEFTIS_sonnet.md`): kurulum, işaretler, 3B→2B indirgeme, K'nın analitik
   devamı (asal başına tam özdeşlik), kimlik → c_P: GEÇTİ. Açık kalan tek kalem (takas teriminin yerel
   Euler çarpanı E_θ ve ∂α₂ log A₁) tarafımızca doğrudan doğrulandı: θ üzerinden sayısal integral 1e-12,
   türev özdeşliği 1e-10.

## Sonuç (7 pencere, L = 9.35 … 22.31; serbest parametre YOK)
| | χ²/7 | en büyük sapma |
|---|---|---|
| **κ₂: oran sanısı** | **6.4** | 1.6σ |
| **κ₃: oran sanısı** | **2.5** | 0.9σ |
| κ₃: yalnız FG limiti (−π²/4 + c_P) | 273.8 | 8.4σ |
| κ₃: eski model CUE(N=L) + c_P (a_k) | 58.2 | 3.6σ |

κ₃ ayrıntı: gözlenen / oran — 9.35: −2.089 / −2.095 · 10.59: −2.124 / −2.109 · 11.66: −2.124 / −2.118 ·
14.19: −2.133 / −2.136 · 16.58: −2.133 / −2.149 · 18.88: −2.153 / −2.158 · 22.31: −2.151 / −2.169.

## Okuma
 - "b0 platosu" (0.28 vs 0.234) tamamen bir sonlu-yükseklik etkisi ve oran sanısı onu parametresiz,
   hata payı içinde veriyor. a_k/CUE(N=L) modeli sonlu-L terimini yanlış veriyordu (χ² 58).
 - ζ'nin sonlu-yükseklik terimi rastgele matrisinkinden farklı; fark asal–sıfır etkileşiminden
   (takas terimindeki aritmetik çarpan K) geliyor.
 - Statü: KOŞULLU TÜRETİM (oran sanısı + limit/integral değiş tokuşları). Titiz bir koşullu teorem için
   oran sanısının kaydırma aralığında düzgünlüğü ve hata terimleri gerekir. Literatürde bu hesabın
   (özellikle κ₃ çekirdeği K) yapıldığını görmedik — öncelik taraması yapılmalı.
 - Sıradaki doğal adım: aynı makineyi sıfırlar üzerindeki ortalamalara taşımak (b1: log|ζ′(ρ)|, b2: yakın
   çift) ⇒ κ₃ merdiveninin "sönmesi" de açıklanabilir mi?

---
## b1: log|ζ′(ρ)|'nun varyansı sıfırlar üzerinde (29 Eylül 2026, gece)
(Tahmin karşılaştırmadan ÖNCE commit'lendi: a056eb0. Kod `203_configs/203_oran_b1.py`.)

**Yöntem.** Sıfır yoğunluğu ağırlığı: Σ_γ f(γ) = ∫ f(t)(1/2π)[L + 2 Re (ζ′/ζ)(½+it)] dt ⇒
E_sıfır[M] = E_t[M] + (E_t[M X₀] + E_t[M X̄₀])/L. log ζ′(ρ) = −∫_0^∞[X(a) − 1_{a<1}/a] da. Böylece b1 varyansı yalnız
b0 için türetilmiş İKİ ve ÜÇ nokta formüllerinden çıkar (yeni sanı yok). CUE modunda aynı kod kesin Palm değerlerini
(κ₂ = Σ_{j<N}[ψ′(j+2) − ½ψ′(j+1)]) N = 2, 5, 10.59, 11.66, 22.3 için 6 basamak veriyor.
**Sayısal tuzak (kayıt):** a → 0'da ~1/a mertebeli terimler götürüşür; asal toplamı kesme hataları büyüyüp ilk koşuda
ızgaraya duyarlı çöp üretti (12 vs 16 düğüm: 0.571 vs 0.606). Bu sayılar commit'lenmedi, veriyle karşılaştırılmadı.
Çözüm: a_min = 0.01/L + doğrusal şerit (gerçek iç integral a → 0'da doğrusal sıfıra gider); sonra 12/16/20 düğüm
5·10⁻⁶ içinde. Kesme Wc 0.6↔0.75: L = 9.35'te 6·10⁻⁴ (en kötü durum), yükseldikçe ihmal edilebilir.

| L | gözlenen κ₂ (b1) | oran sanısı | sapma | eski a_k modeli |
|---|---|---|---|---|
| 9.35 | 0.4672(22) | 0.4691 | −0.8σ | 0.4527 (+6.6σ) |
| 10.59 | 0.5200(14) | 0.5190 | +0.7σ | 0.5035 (+11.8σ) |
| 11.66 | 0.5586(11) | 0.5584 | +0.2σ | 0.5437 (+13.6σ) |
| 14.19 | 0.64107(55) | 0.64099 | +0.2σ | 0.6277 (+24.4σ) |
| 16.58 | 0.70628(54) | 0.70783 | −2.9σ | 0.6958 (+19.4σ) |
| 18.88 | 0.76478(61) | 0.76496 | −0.3σ | 0.7539 (+17.9σ) |
| 22.31 | 0.83829(82) | 0.83925 | −1.2σ | 0.8294 (+10.9σ) |

**χ²/7: oran sanısı 10.9 (p ≈ 0.14), eski a_k modeli 1780.** Parametresiz; hata payları b0'dakinin ~5–8 katı dar.
Tek büyük sapma C2 (−2.9σ, 0.0015); oran sanısının kendi hata terimi O(T^{−1/2+ε}) bu yüksekliklerde ~10⁻⁴ mertebesinde.

**Okuma.** Sıfırların üzerindeki ortalama da aynı makineyle, aynı formüllerle, binde birkaç hassasiyetle tutuyor. b1'deki
"a_k'dan +0.015 fazlalık" (200-A) tamamen sonlu-yükseklik etkisiymiş. Sıradaki: b1 κ₃ (dört-nokta ortalamalar gerekir;
tasarım notu hafızada) ve b2.

**Teftiş (sonradan; protokolde karşılaştırmadan önce olmalıydı — sıra kaçırıldı, dürüstçe kayıt):** Sonnet,
`203_configs/203b1_TEFTIS_sonnet.md`. Madde 1–4 GEÇTİ, 5 düzeltmeyle GEÇTİ. Öne çıkanlar: (i) yoğunluk ağırlığı gerçek
Odlyzko sıfırları üzerinden Poisson toplamıyla 5·10⁻⁷ doğrulandı; (ii) BAĞIMSIZ SINAMA: E_ρ[X(a)X̄(b)]-tipi integrand
L = 14.19'da 8·10⁵ LMFDB sıfırından doğrudan hesaplandı — 15 hücrenin hepsi |çekme| ≤ 0.5 (CUE %25 sapıyor);
(iii) CUE Palm κ₂ N = 2…22.3 için −8·10⁻⁸; (iv) Wc etkisi ζ için ~+0.9·10⁻⁴ (L = 9.35); (v) DÜZELTME: W1–W3 geniş L aralıkları
(Var L = 0.26/0.08/0.09; W1'in sıfır-ağırlıklı L̄ = 9.3427) ve κ₂(L) içbükey ⇒ pencere-ortalamalı tahminler −7.6·10⁻⁴ /
−1.4·10⁻⁴ / −1.7·10⁻⁴ kayar; χ² 10.93 → 10.69 (p = 0.15). Aynı düzeltme b0 W pencerelerine de uygulanmalı (b0 hata payları
bunun 5–25 katı; etkisi ihmal edilebilir ama bir sonraki sürümde yapılacak). Kanıtlanmamış kalan: oran sanısının c → 0'da
düzgünlüğü ve a < 0.1/L bölgesi (sıfırlarla sınanamıyor).

---
## b1: log|ζ′(ρ)|'nun ÜÇÜNCÜ kümülantı (29 Eylül 2026 akşam)
(Tahmin karşılaştırmadan ÖNCE commit'lendi: e5f2580. Kod: `203_configs/203_b1k3_zeta.py` (+ `203_cue4.py`, `203_b1k3_cue.py`,
`203_yerel4*.py`, `203_hizli.py`).)

**Yöntem.** κ₃ = ¼[κ(Y,Y,Y) + 3κ(Y,Y,Ȳ)], bağlantılı sıfır-kümülantları yoğunluk ağırlığıyla; üç-nokta (T) ve dört-nokta
(Q: 3+1, Q2: 2+2, çift takas dahil) oran ortalamaları. Hesap: κ₃ = κ₃^CUE(L) (kesin Palm, analitik) + ∭Δ, Δ = I_ζ − I_CUE
yalnız a,b,c ≥ s0 bölgesinde; s0-levhaları doğrusal uzatmayla.
**Doğrulama merdiveni:** (1) CUE dört-nokta formülleri Weyl integraline karşı N = 1, 2, 3'te 8 basamak; CUE b1 κ₃ Palm değeri
2–5·10⁻⁷. (2) Asal başına 4-nokta aritmetik parçalar θ-integraliyle 10⁻¹², türevliler sayısal türevle 10⁻⁸. (3) Hızlı tablolar
kesin sürümlere karşı (K maks 2.6·10⁻⁶; p = 2 kesin). (4) BAĞIMSIZ: bağlantılı üç-nokta kümülantları 1.5·10⁶ gerçek LMFDB
sıfırından (Hadamard toplamı, K→∞ Richardson) — 14/14 |çekme| ≤ 1.1, CUE %9–26 sapıyor (`sifir_sinamasi_b1k3/RAPOR.md`).
(5) Duyarlılık (L = 11.66): s0 0.005/0.0075/0.01 → 5·10⁻⁵; ızgara 4→6 → 10⁻⁴; çift-takas ara değeri → 5·10⁻⁵ (9.35'te 9·10⁻⁵).
Seçimler (s0 = 0.0075, n = 6) iç ölçütle, karşılaştırmadan önce.
**Yolda yakalanan hatalar (kayıt):** çift takas Y_p'de sıfırlanan paydaların (1−1/p)⁴ çarpanı eksikti; f_ij'nin büyük-p
asimptotiği yanlış varsayılmıştı (~−22 sahte kayma); ızgara s0'dan küçük kırılma noktasıyla başlıyordu (levha ile üst üste
binme). Hepsi karşılaştırmadan önce, iç sınamalarla bulundu. A_DS'nin Euler çarpımı yalnız x < ~0.45'te yakınsar; dışında
tekil limitleri (sayısal doğrulandı) koruyan ara değer kullanıldı — etkisi ≤ 10⁻⁴.
**Dürüstlük notu:** W3'ün gözlenen değeri (−0.0509) 200-A raporundan önceden biliniyordu; hiçbir seçim ona göre yapılmadı
(tahmin W3'te −1.9σ). Biçimsel Sonnet türetim teftişi bu adım için henüz yapılmadı (gerçek-sıfır sınaması yapıldı).

| L | gözlenen κ₃ (b1) | oran sanısı | sapma | eski a_k | CUE(N=L) |
|---|---|---|---|---|---|
| 9.35 | −0.0465(27) | −0.0488 | +0.8σ | +0.043 | −0.190 |
| 10.59 | −0.0519(21) | −0.0492 | −1.3σ | +0.036 | −0.198 |
| 11.66 | −0.0509(10) | −0.0490 | −1.9σ | +0.031 | −0.203 |
| 14.19 | −0.0487(13) | −0.0480 | −0.6σ | +0.022 | −0.212 |
| 16.58 | −0.0440(12) | −0.0467 | +2.2σ | +0.015 | −0.218 |
| 18.88 | −0.0445(15) | −0.0454 | +0.6σ | +0.011 | −0.223 |
| 22.31 | −0.0446(19) | −0.0437 | −0.5σ | +0.005 | −0.229 |

**χ²/7: oran sanısı 11.6 (p ≈ 0.11) · eski a_k modeli 16 826 · CUE 90 185.**

**Okuma.** Merdivenin ilk iki basamağı (b0: κ₂, κ₃; b1: κ₂, κ₃) sonlu yükseklikte, serbest parametresiz, oran sanısından
çıkıyor. 200-A'nın "κ₃ asal izi b1'de sönüyor" (+0.150 vs a_k +0.234) bulgusu tamamen bir sonlu-yükseklik etkisi; a_k'yı
CUE(N=L)'ye eklemek bu yükseklikte yanlış tahmin veriyor (χ² 16 826). Sıradaki: b2 (yakın çiftler).

**b1 κ₃ teftişi (Sonnet, karşılaştırmadan SONRA; `203_configs/203b1k3_TEFTIS_sonnet.md`):** reçete defter tutması (asal başına
bağımsız θ-karesel 10⁻¹⁶), CUE (Weyl N=1–4, 8–9 basamak; Palm 1e-7), montaj (CUE girdileriyle 6e-14): GEÇTİ. HATA: f_jk kuyruğu
log²p için k = 3 formülüyle hesaplanmış (10× büyük) ⇒ −3.1e-4 (L=11.66), −4.6e-4 (L=9.35). Sistematik bütçe "~1e-4" iyimserdi:
A_DS ara değeri +3e-4 (9.35) / +7e-5 (11.66); p<1000 çarpım yanlılığı +1.8e-4 / +1.0e-4; levha gerekçesi yalnız simetrik Δ için
geçerli (kenar ≤ 5e-5). Tüm düzeltmelerle tek tam koşu (L=11.66): −0.049262 (commit'li −0.049007; fark −2.6e-4). Pencere-karışım
düzeltmeleri ≤ 2e-4. Düzeltilmiş χ² 11.7/7 — sonuç değişmiyor (hata payları 1–2.7e-3). Hata kodda düzeltildi; 7 pencere
teftiş-sonrası yeniden koşuluyor (kayıt amaçlı; commit'li tahminler değiştirilmez).

**Teftiş-sonrası düzeltilmiş b1 κ₃ (f_jk kuyruk hatası giderildi; kayıt amaçlı — commit'li tahminler e5f2580'de değişmeden durur):**
L = 9.35 / 10.59 / 11.66 / 14.19 / 16.58 / 18.88 / 22.31 ⇒ −0.04924 / −0.04952 / −0.04932 / −0.04819 / −0.04686 / −0.04559 / −0.04385;
χ²/7 = 10.9 (commit'li tahminle 11.6). Sonuç değişmiyor.
