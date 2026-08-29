# KALEM OTURUMU — BK (4.23): RESURGENCE (29 Ağu 2026)

## 1. Motoru ele geçirmek: yeniden-toplama özdeşliği (kendi türetimimiz)

Köşegen merdivenin yoğunluk-ağırlıklı hali F(ξ) = Σ_{p,m} log²p·p^{−m(1−iξ)}.
Geometrik seri cebiriyle (Σ m x^m = x/(1−x)², Σ x^m = x/(1−x)):

  **F(ξ) = (ζ'/ζ)'(1−iξ) − Σ_p log²p/(p^{1−iξ}−1)²**

Birinci terim = −∂²_ξ log ζ(1−iξ) — BK (4.23)'ün çekirdeği. ✓ (140'ta
sayısal doğrulama: Re s=1.05'te mutlak-yakınsak iki taraf.)

## 2. Resurgence yapısı

(ζ'/ζ)'(s) = −Σ_ρ 1/(s−ρ)² + düz (Hadamard). s = 1−iξ'de her sıfır
ρ = ½+iγ için s−ρ̄ = ½ − i(ξ−γ):

  Re[−1/(s−ρ̄)²] = −(¼−u²)/(¼+u²)²,  u = ξ−γ

⇒ her alçak sıfır, ξ = γ noktasında **derinlik 4, genişlik ½** bir çukur
kazar. Log-düzeyde (H = Σ p^{−ms}/m², C(n)'nin motoru): tek sıfırın izi
½·log(¼+u²) — çukur derinliği log 2 ≈ 0.69. "Alçak sıfırlar, yüksek
sıfırların korelasyonunda rezonans" — Snaith Şek. 3'ün çukurları,
Lu–Sridhar'ın "resurgence"ı; kaynak bu terim.

## 3. Kapı haritası: resurgence BİZDE nereye bağlanır?

t-birimi ayrım ξ ↔ gap-birimi gecikme n: ξ_n = n·ḡ = 2πn/L.
Rezonanslar ξ = γ_k ⇒ **n*_k = γ_k·L/2π**:

  L=12.03 (son-300k): n* = 27.1, 40.3, 47.9
  L=11.47 (orta-300k): n* = 25.8, 38.4, 45.7   ← L ile KAYAR (imza!)

- **Bağ ölçeği (n=1, bizim K(τ)):** ξ = ḡ ≈ 0.52 — rezonanslardan çok
  uzak; orada (4.23)'ün içeriği s=1 KUTBU (PNT) + m≥2 düzeltmeleri.
  DÜRÜST HÜKÜM: K(τ_p) kayışının (−2.87→−3.25) kaynağı resurgence
  DEĞİL gibi; kayış adayları tam-merdiven/m≥2 inceltmeleri + 135'in
  atılan τ_p-terimleri (K3 listesine işlendi).
- **Uzun gecikmeler (n ≈ 26–48):** resurgence burada yaşar ve
  ÖLÇÜLEBİLİR bir kehanet verir: C(n) = ⟨ds₀ds_n⟩'de n*'larda yerel
  anomali. Köşegen tam-merdiven bunu taşır (Fourier: Δξ~½'lik yapı
  log p ≲ L ister — elimizdeki kuyruk yeter); PÜRÜZSÜZ (PNT-sürekli)
  merdiven taşımaz. Yani anomali, merdivenin aritmetik inceliğinin —
  gazın kendi alçak sıfırlarını tanımasının — gap-dilindeki imzası.

## 4. Oturum kararı

140: C(n) resurgence avı, iki pencerede, ön-mühürlü (konum + işaret +
L-kayması; genlik kayda — 139'un koşullama dersi gereği ekran çarpanı
serbest). Tutarsa: (4.23) kapısı ölçümle bağlanır, "gaz alçak
sıfırlarını da tanır" sözlük adayı olur; tutmazsa: koşullama farkının
uzun-gecikme davranışı öğrenilir — iki durumda da defter kazanır.

## 5. SONUÇ (aynı gün, 140 koşusundan)

Kehanet tuttu — üç kapı da: altı çukur iki pencerede tam γ_k·L/2π'de
(son: 27/40/48, orta: 26/38/46), her biri ~4σ, işaret ve omuz şekli
tam-merdiven öngörüsüyle; pürüzsüz merdiven kör. Genlik oranı tekdüze
0.76±0.03 (koşullama çarpanının uzun-gecikme hali — açık iş). Özdeşlik
motoru mutlak-yakınsak bölgede 7e-6 ile doğrulandı (1.05'teki %15 saf
kesme kuyruğu). Bağ-ölçeği hükmü değişmedi: K(τ_p) kayışı resurgence
işi değil; kayış adayları m≥2/tam-merdiven inceltmesi + 135'in
τ_p-terimleri. Sözlüğe aday: "Gaz alçak sıfırlarını tanır".
