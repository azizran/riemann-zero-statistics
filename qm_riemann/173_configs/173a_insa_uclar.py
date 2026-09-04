"""
173a — HAKEM GAZLARI: λ = 0.40 ve λ = 1.45 İNŞASI
==================================================
Hiçbir çözüm/inşa parçası KOPYALANMAZ: `167_insa.main` AYNEN çağrılır
(170_insa60.py / 171b_insa_uclar.py ile birebir aynı yol); tek yapılan
`167_insa.KONFIG`'e λ = 0.40, λ = 1.45 (ve YEDEK λ = 1.40) satırlarını
EKLEMEKTİR.

    S(t) = −Σ_q (λ a_q w_q) sin(ω_q t),  a_q = 1/(π m √q)
    keskin kesim τ ≤ 1.00 (Hkeskin/L050…L130 ile AYNI merdiven);  c = −½;
    300000 sıfır, h = 0.015, 15450 çizgi.

NEDEN BU İKİ UÇ (173 kalemi): 172 §G4.4, `172/G4.json`'da İKİ BAĞIMSIZ
YOL bu iki λ için M ve c mühürledi ve λ = 0.40'ta **%6.7 ayrışıyorlar**.
Tek ölçüm hakemdir. Ayrıca dokuz-noktalı merdiven, 172'nin iki tümseğini
(Q_E ↔ λ_c, α ↔ λ*) ve 171'in ν(λ) resmini uçlarda sınar.

══════════════════════════════════════════════════════════════════════
ÖN-MÜHÜR — İNŞA MARJİNALLERİ (koşudan ÖNCE yazıldı, 4 Eylül 2026 13:44)
══════════════════════════════════════════════════════════════════════
Girdi: yalnız ölçülmüş yedi λ gazı (önbellekten, hiçbiri yeniden
koşulmadı):

  λ        0.50      0.60      0.70      0.85      1.00      1.15      1.30
  σ_ds   0.27903   0.32135   0.35586   0.39765   0.43134   0.45958   0.48890
  σ_X̃    0.15576   0.17951   0.19897   0.22272   0.24204   0.25839   0.27713
  σ_Ĉ    0.18131   0.20730   0.22892   0.25565   0.27768   0.29653   0.31923
  minΔz  0.18913   0.16397   0.14485   0.12343   0.10675   0.09339   0.08348

(A) MERDİVENDEN ÖZDEŞ GELEN (uyum YOK, λ ile TAM doğrusal):
    rms S′ = 1.8940·λ  ⇒  **0.75760** (λ=0.40),  **2.74630** (λ=1.45)
    Σ a_q ω_q = 259.9·λ ⇒ 103.96 / 376.86 ;  rms S = 0.3825·λ ⇒ 0.1530 / 0.5546
    N̄′ = 1.9147 (λ'dan bağımsız).  nline = 15450, Σa/Σa_ham = λ.

(B) UZATMA ve NEDEN BANT GENİŞ TUTULUYOR — DÜRÜSTLÜK NOTU.
    171 §T2c.0: "kuvvet-yasası uzatması λ > 1'de ÜÇÜNCÜ kez aşağıdan
    ıskaladı". Bunu bu betikte SAYIYLA ölçtüm (bir-adım-dışarı sınavı,
    log-log kuadratik; ölçülmüş yedi noktanın kendi üstünde):
        1.30'u öngör (0.50…1.15 ile): σ_ds −1.87%, σ_X̃ −2.48%, σ_Ĉ −2.62%
        1.30'u öngör (0.50…1.00 ile): σ_ds −2.80%, σ_X̃ −3.41%, σ_Ĉ −3.36%
        0.50'yi öngör (0.60…1.30 ile): σ_ds +1.97%, σ_X̃ +2.37%, σ_Ĉ +2.26%
    Yani uzatma ÜST uçta AŞAĞIDAN, ALT uçta YUKARIDAN ıskalıyor —
    log-log'da eğrilik kuadratiğin taşıdığından fazla. Bu yüzden aşağıdaki
    bantlar hem GENİŞ hem de **ÇARPIK** (üst uçta yukarı, alt uçta aşağı
    kaydırılmış) yazılmıştır. Kurtarma değil: bant koşudan önce yazıldı.

    ham kuadratik uzatma        →  çarpıklık düzeltmesi  →  ÖN-KAYIT BANDI
    σ_ds(0.40) = 0.23378           −2.0%  ⇒ 0.2292        **[0.224, 0.238]**
    σ_ds(1.45) = 0.50583           +2.0%  ⇒ 0.5160        **[0.503, 0.527]**
    σ_X̃(0.40) = 0.13096           ⇒ 0.1273               **[0.1240, 0.1330]**
    σ_X̃(1.45) = 0.28731           ⇒ 0.2913               **[0.2845, 0.2990]**
    σ_Ĉ(0.40) = 0.15430           ⇒ 0.1478               **[0.1440, 0.1545]**
    σ_Ĉ(1.45) = 0.33184           ⇒ 0.3336               **[0.3265, 0.3420]**
    (σ_X̃/σ_ds ve σ_Ĉ/σ_ds oranları logλ'da doğrusal uzatıldı:
     0.5554 / 0.6449 (λ=0.40) ve 0.5651 / 0.6471 (λ=1.45).)

    min Δz: min_dz ≈ λ^(−0.857) — bir-adım-dışarı sınavı burada ±%0.6
    (uzatma İYİ). ÖN-KAYIT: **0.223 ∈ [0.214, 0.233]** (λ=0.40) ve
    **0.0749 ∈ [0.069, 0.081]** (λ=1.45).

    ΔG<0 kesri (F′ işaret değiştirme): Gauss kestirimi P(Z < −N̄′/rms S′)
    0.00575 (0.40) ve 0.24284 (1.45); gözlem/Gauss oranı λ ile düşüyor
    (0.49 @0.50, 0.96 @0.60, 1.23 @1.00, 1.19 @1.15, 1.15 @1.30).
    ÖN-KAYIT: **0.002 ∈ [0.0005, 0.005]** (0.40, oran ≈0.2–0.9 belirsiz)
    ve **0.272 ∈ [0.258, 0.288]** (1.45, oran ≈1.12).

(C) İNŞA KAPILARI — HEPSİ GEÇİLMELİ, KURTARMA YOK:
    (K1) ilk-kök hücreleri benzersiz **300000 / 300000**
    (K2) maks|F| ≤ **1e−8**
    (K3) sıralılık **TAM**
    Biri düşerse gaz REDDEDİLİR, ölçüm YAPILMAZ, gerekçesiyle raporlanır
    ve λ = 1.45 yerine **λ = 1.40**'a çekilinir (aşağıda ön-kayıtlı).

(D) RİSK — λ = 1.45. min Δz ≈ 0.075, ızgara adımı h = 0.015 ⇒ en dar
    tekne yalnız **~5 hücre**. F′ zamanın ~%27'sinde negatif. Hücre
    benzersizliği kırılırsa (K1) gaz düşer. Ölçek kestirimi: min_dz = h
    ancak λ ≈ 10'da olur, yani K1'in gerçek riski DÜŞÜK; asıl risk
    Newton'un braket içinde daha çok ikiye-bölme istemesidir (L130'da 6
    adım gerekti) — bu bir kapı değil, yalnız süre.

(E) YEDEK ÖN-KAYIT — λ = 1.40 (yalnız λ=1.45 kapılardan düşerse kurulur):
    rms S′ = 2.65160 ; σ_ds ≈ **0.5100** [0.497, 0.521] ;
    σ_X̃ ≈ **0.2880** [0.2810, 0.2950] ; σ_Ĉ ≈ **0.3300** [0.3230, 0.3380] ;
    min Δz ≈ **0.0776** [0.072, 0.084] ; ΔG<0 ≈ **0.263** [0.250, 0.278].

SONUÇ bloğu YALNIZ gerçek koşu çıktısındandır
(scratchpad/173/log_insa_*.txt).

Kullanım: 173a_insa_uclar.py <L040|L145|L140> [h] [nz] [nwork]

══════════════════════════════════════════════════════════════════════
SONUÇ (yalnız gerçek koşudan; scratchpad/173/log_insa_*.txt)
══════════════════════════════════════════════════════════════════════
 L040 (9.5 dk) — ÜÇ KAPI DA GEÇTİ (300000/300000, maks|F| = 1.863e-9,
   sıralılık TAM, 0 ikiye-bölme).
   rms S' = 0.7576 ✓;  σ_ds = 0.22728 ✓ (bandın alt kenarında);
   σ_X̃ = 0.12673 ✓;  σ_Ĉ = 0.15011 ✓;  min Δz = 0.223681 ✓ (merkez!);
   ΔG<0 = 0.0002 — ÖN-KAYIT [0.0005, 0.005] ⇒ **ISKA** (Gauss/gözlem
   oranı λ=0.50'de 0.49 idi, 0.40'ta 0.035'e düştü).
   NOT: ÇARPIKLIK DÜZELTMESİ ÇALIŞTI — ham kuadratik üç σ'yı da
   +%2.7…+%3.2 yukarıdan ıskalıyordu; düzeltilmiş merkez yalnız
   +%0.8/+%0.4/−%1.5 ıskaladı.
 L145 (9.8 dk) — ÜÇ KAPI DA GEÇTİ (300000/300000, 1.863e-9, TAM);
   36 ikiye-bölme adımı gerekti (yazılı risk yalnız SÜRE olarak
   gerçekleşti).
   rms S' = 2.7463 ✓;  min Δz = 0.075756 ✓;  ΔG<0 = 0.2717 ✓;
   **σ_ds = 0.5444 — ÖN-KAYIT [0.503, 0.527] ⇒ ISKA (+%3.3 ÜSTÜNDE)**;
   σ_X̃ = 0.31782 (ISKA +%6.3);  σ_Ĉ = 0.37422 (ISKA +%9.4).
   Yerel kuvvet üsteli p = dlogσ_ds/dlogλ 1.30→1.45 adımında 0.504'ten
   **0.985**'e fırlıyor: bu bir uzatma hatası değil, REJİM DEĞİŞİMİ.
   Marjinallerin İÇ oranları da kırılıyor (σ_X̃/σ_ds 0.5567→0.5837,
   σ_Ĉ/σ_ds 0.6530→0.6873) — yedi gazın yayılımının dışında.
 L140 (9.4 dk, YEDEK — λ=1.45'in İNŞA değil ÖLÇÜM kapısı düştüğü için;
   ayrı ön-mühür scratchpad/173/ONMUHUR_L140_Rbant.txt, 14:14:49):
   ÜÇ KAPI DA GEÇTİ (300000/300000, 1.863e-9, TAM, 30 ikiye-bölme);
   rms S' = 2.6516 ✓;  σ_ds = 0.52231 (ISKA +%0.25 — ön-mühür bu ıskayı
   ÖNCEDEN uyarmıştı);  min Δz = 0.078138 ✓;  ΔG<0 = 0.2654 ✓.
   λ = 1.40'ta 169'un sağlık filtresi BEŞ BANTTA DA GEÇİYOR
   (R_bant 0.9829-1.0370) ⇒ SADAKAT SINIRI 1.40 ile 1.45 ARASINDADIR.
"""
import importlib
import sys
from pathlib import Path

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
sys.path.insert(0, str(QM / "167_configs"))
sys.path.insert(0, str(QM / "164_configs"))
I167 = importlib.import_module("167_insa")
I164 = importlib.import_module("164_insa")

I167.KONFIG["L040"] = dict(lam=0.40)
I167.KONFIG["L145"] = dict(lam=1.45)
I167.KONFIG["L140"] = dict(lam=1.40)      # yalnız YEDEK

if __name__ == "__main__":
    ad = sys.argv[1]
    if ad not in ("L040", "L145", "L140"):
        raise SystemExit("kullanım: 173a_insa_uclar.py <L040|L145|L140>")
    I167.main(ad,
              float(sys.argv[2]) if len(sys.argv) > 2 else I164.HIZGARA,
              int(sys.argv[3]) if len(sys.argv) > 3 else I164.NZERO,
              int(sys.argv[4]) if len(sys.argv) > 4 else I164.NWORK)
