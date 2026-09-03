"""
171b — YENİ GAZLAR: λ = 0.50 ve λ = 1.30 (T2c'nin İKİ UCU)
==========================================================
Hiçbir çözüm/inşa parçası KOPYALANMAZ: `167_insa.main` AYNEN çağrılır
(170_insa60.py ile birebir aynı yol); tek yaptığımız `167_insa.KONFIG`'e
λ=0.50 ve λ=1.30 satırlarını EKLEMEK.

    S(t) = −Σ_q (λ a_q w_q) sin(ω_q t),  a_q = 1/(π m √q)
    keskin kesim τ ≤ 1.00 (Hkeskin/L060…L115 ile AYNI merdiven);  c = −½.

NEDEN BU İKİ UÇ (171 kalemi, T2c):
  λ=1.30 — "M(λ) λ≥1'de DÜZ" hükmünün ÖRNEKLEM-DIŞI sınavı (M(1.15)=1.0002).
  λ=0.50 — M'nin yükselen kolunun taban-ötesi ucu (M(0.60)=1.1272).

ÖN-MÜHÜR — MARJİNALLER (koşudan ÖNCE yazıldı, 3 Eylül 2026 23:36).
  Ölçülmüş λ-serisi:
     λ      0.60      0.70      0.85      1.00      1.15
     σ_ds   0.32135   0.35586   0.39765   0.43134   0.45958
     σ_X̃    0.17951   0.19897   0.22272   0.24204   0.25839
     σ_Ĉ    0.20730   0.22892   0.25565   0.27768   0.29653
  Yerel kuvvet üsteli p = dlogσ_ds/dlogλ = 0.662/0.572/0.500/0.454
  (λ_orta = 0.648/0.771/0.922/1.072) ⇒ p(logλ) = 0.4742 − 0.4130·logλ.
  İKİ BAĞIMSIZ ARAÇ (log-log kuadratik uyum ‰1 artıkla; yerel kuvvet yasası):
     σ_ds(0.50) = 0.28198 / 0.28167     σ_ds(1.30) = 0.48173 / 0.48216
  σ_X̃/σ_ds ve σ_Ĉ/σ_ds oranları λ boyunca neredeyse sabit
  (0.5586→0.5622 ve 0.6429→0.6452); logλ'da doğrusal uzatma:
     σ_X̃/σ_ds(0.50) = 0.5574   (1.30) = 0.5627
     σ_Ĉ/σ_ds(0.50) = 0.6439   (1.30) = 0.6442

  ÖNGÖRÜ (koşudan ÖNCE):
     σ_ds(0.50) ∈ [0.279, 0.285]   σ_X̃(0.50) ∈ [0.1555, 0.1590]
                                    σ_Ĉ(0.50) ∈ [0.1795, 0.1835]
     σ_ds(1.30) ∈ [0.477, 0.487]   σ_X̃(1.30) ∈ [0.2685, 0.2742]
                                    σ_Ĉ(1.30) ∈ [0.3073, 0.3138]
     rms S' = 1.894·λ  ⇒  0.947 (λ=0.50) ve 2.462 (λ=1.30);  N̄' = 1.9147.
     ΔG<0 kesri (F' işaret değiştirme): L060 0.0442, L115 0.2252 ölçüldü;
     Gauss kestirimi P(Z < −N̄'/rmsS') ile 0.96–1.19 çarpanı içinde ⇒
        öngörü: λ=0.50 → 0.019 ± 0.005   λ=1.30 → 0.26 ± 0.03
     min Δz: L060 0.1640, L115 0.0934 ⇒ öngörü λ=0.50 → ~0.21,
        λ=1.30 → ~0.055 (0.03–0.08).
  RİSK (açıkça yazılıyor): λ=1.30'da merdiven eğimi zamanın ~%26'sında
  negatif; SIRALILIK KIRILABİLİR. Kırılırsa gaz REDDEDİLİR ve raporda
  kurtarmasız yazılır (ölçüm yapılmaz).

SONUÇ bloğu YALNIZ gerçek koşu çıktısındandır (scratchpad/171/log_insa_*.txt).

Kullanım: 171b_insa_uclar.py <L050|L130> [h] [nz] [nwork]

SONUÇ (gerçek koşudan, 4 Eylül 2026 gece):
 L050 (9.3 dk): rms S' = 0.9470 ✓;  σ_ds = 0.27903 (aralığın ALT sınırında ✓),
   σ_X̃ = 0.15576 ✓, σ_Ĉ = 0.18131 ✓;  ΔG<0 kesri = 0.0106 (öngörü
   0.019 ± 0.005 — ISKA);  min Δz = 0.18913 (öngörü ~0.21 — ISKA);
   sıralılık TAM, maks|F| = 1.86e−9, ikiye-bölme 0, hücreler 300000/300000.
 L130 (9.4 dk): rms S' = 2.4622 ✓;  σ_ds = 0.48890 (aralığın ÜSTÜNDE — ISKA),
   σ_X̃ = 0.27713 (+%1.1 ISKA), σ_Ĉ = 0.31923 (+%1.7 ISKA);
   ΔG<0 kesri = 0.2511 ✓ (öngörü 0.26 ± 0.03);  min Δz = 0.08348 (ISKA);
   SIRALILIK TAM — risk gerçekleşmedi; 6 ikiye-bölme adımı gerekti.
 Kuvvet-yasası uzatması λ > 1'de ÜÇÜNCÜ kez aşağıdan ıskaladı (L115'te
 %1.1, L130'da %1.4) — artık bir yasa gibi okunmalı.
"""
import importlib
import sys
from pathlib import Path

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
sys.path.insert(0, str(QM / "167_configs"))
sys.path.insert(0, str(QM / "164_configs"))
I167 = importlib.import_module("167_insa")
I164 = importlib.import_module("164_insa")

I167.KONFIG["L050"] = dict(lam=0.50)
I167.KONFIG["L130"] = dict(lam=1.30)

if __name__ == "__main__":
    ad = sys.argv[1]
    if ad not in ("L050", "L130"):
        raise SystemExit("kullanım: 171b_insa_uclar.py <L050|L130>")
    I167.main(ad,
              float(sys.argv[2]) if len(sys.argv) > 2 else I164.HIZGARA,
              int(sys.argv[3]) if len(sys.argv) > 3 else I164.NZERO,
              int(sys.argv[4]) if len(sys.argv) > 4 else I164.NWORK)
