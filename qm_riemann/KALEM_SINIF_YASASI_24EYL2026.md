# KALEM — 193: KALINTI-SINIFI YASASI — temiz, kör sınav
(24 Eylül 2026 — 192'nin düzeltilmiş cos okumasının kör sınavı; kullanıcı onaylı;
TÜRETİM TEFTİŞİNDEN GEÇTİ — Sonnet, "düzeltmeyle mühürlenebilir"; düzeltmeler işlendi)

## Durum

192: μ(a)/φ(a) seçim kuralı ayakta; "sınıflar 120° döner" öngörüsü öldü (kaptan
hatası — çekirdek yalnız Re Ĝ duyar). Veri-sonrası okuma: sınıflar aynı eksende
cos(2πrb/a) ağırlıklarıyla. Bu kalem cos yasasını KÖR ve KİRLİLİKSİZ sınar.

## Kalem cebiri (teftiş düzeltmeleriyle)

1. (192) Durağan faz: Δω = log(a/b) uydusunda, durağan noktası blok içinde olan
   her çizgi q' için Ĝ_mid(ω') katkısı ∝ A_n · e^{i(2π q' b/a + φ0)}, A_n gerçel.
2. **φ0 TÜRETİLDİ (teftiş (c) düzeltmesi):** 2πN̄(t) = tL − t + 7π/4 ⇒ durağan
   değer t* − 7π/4; faz eğriliği φ'' = −1/t < 0 ⇒ durağan-faz çarpanı e^{−iπ/4};
   orta noktalarda N(m_n) tam sayı ⇒ ek işaret yok. Toplam sabit −7π/4 − π/4 =
   −2π ⇒ **φ0 ≡ 0.** (192 kaleminde −π/4 yazılmamıştı — kaptan eksiği.) +log3
   ve +log6'nın ½/½ görünümü bununla tutarlı; artık uydurma değil türetim.
3. Çekirdek gerçel kısmı duyar: κ(s) ≈ Σ a'·πτ'cos(πτ')·Re Ĝ_mid (S_Re; corr
   −0.95/−0.85 — EŞİTLİK DEĞİL, ~%90-95 açıklayan vekil; teftiş (b)). ⇒ sınıf r
   (q' ≡ r mod a) katkısı ∝ cos(2π r b/a).
4. **PAY YASASI:** s_r := κ_r/κ_top = cos(2π r b/a) / μ(a) (payda Ramanujan
   toplamı; çizgiler sınıflara eşit dağılır — Dirichlet).

## Öngörüler

| uydu (Δω) | a | sınıf payları s_r | statü |
|---|---|---|---|
| **+log 10 (2.303)** | 10 | r∈{1,9}: **+0.809**; r∈{3,7}: **−0.309** | **KÖR, birincil** |
| **+log 7 (1.946)** | 7 | r∈{1,6}: **−0.623**; r∈{2,5}: **+0.223**; r∈{3,4}: **+0.901** | **KÖR, ikincil** |
| +log 5 (1.609) | 5 | r∈{1,4}: −0.309; r∈{2,3}: +0.809 | doğrulayıcı KAYIT (192 keşfinde görülmüştü — teftiş (f)) |
| +log 3 (1.099) | 3 | ½, ½ | kontrol (192 tekrarı) |
| −log 3 (−1.099) | 1 | tek koherent sınıf; keyfi mod-3 dilimleri ½, ½ çıkmalı (sıfır-yapı kontrolü) | kontrol |

+log(5/2), +log(5/3) ÇIKARILDI (teftiş (e): μ=0 çukurlu +log(8/3), +log(9/5)'e 0.065/0.077).

## Ölçüm (kirliliksiz)

Son penceresi (L=12.03), 188 makinesi AYNEN (seri bit-bit mühür). 8 blok; her
blokta yerel L_b ile Δω = ω_{q'} − L_b; uydu seçimi |Δω − log n| < **δ = 0.03**.
ÖLÇÜMDEN ÖNCE: her hedef için a, b ≤ 20 rasyonel kataloğunda |Δω farkı| < 0.10
olan komşuların listesi (μ(a) ve 192 çukur bilgisiyle) raporlanır; 2δ içinde
komşusu olan hedef "kirlenme riskli" etiketlenir (hüküm yine verilir, etiketle).
Sınıf başına κ_r (−Re K katkısı, HAVUZ), blok jackknife; s_r = κ_r/κ_top ± se.
a'yı bölen çizgiler ayrı "bölünen" sınıfı (paydan hariç, raporlanır).

## Hipotezler (ön-kayıtla donar; ölümler kurtarmasız)

- **H-193a (birincil, +log10):** r∈{3,7} payları negatif, r∈{1,9} pozitif (her biri
  ≥ 2σ); nicel: her pay öngörünün ±0.25 içinde (±0.20 + vekil payı 0.05).
  ÖLÜM: r∈{3,7}'den biri ≥ 2σ pozitif ya da r∈{1,9}'dan biri ≥ 2σ negatif.
- **H-193b (ikincil, +log7):** r1, r6 negatif; r3, r4 pozitif (≥ 2σ); nicel ±0.30.
  ÖLÜM: r1 ya da r6 ≥ 2σ pozitif. Güç yetersizse (κ_top < 4σ) "erişilemedi".
- **KAYIT:** +log5 (doğrulayıcı), +log3 ve −log3 kontrolleri (½/½ ± 0.1).

Ölümler kurtarmasız; ajan git'e DOKUNMAZ; sonuç ORTAK TEFTİŞE.
