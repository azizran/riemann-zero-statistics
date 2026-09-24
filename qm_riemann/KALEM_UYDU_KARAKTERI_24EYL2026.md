# KALEM — 195: UYDULARIN KARAKTERİ — L-fonksiyonu adalarında tarak iptali
(24 Eylül 2026 — büyük sefer; kullanıcı onaylı. TÜRETİM TEFTİŞİNDEN GEÇTİ (Sonnet):
adımlar 2-5 sağlam; karakter toplamı 56 kombinasyonda sayısal doğrulandı; iki
uyarı kapı olarak eklendi (C_χ mutlak kalibrasyonu, χ₈ ilkellik tuzağı))

## Durum

Zeta'da (188-193): pencere-ötesi iptal orta-nokta örgüsünün Bragg tarağı ve
uydularıyla taşınır; uydu Δω = log(a/b); durağan faz çizgi q'nün kalıntı sınıfını
okur; çekirdek koherensi μ(a)/φ(a); sınıf payları cos(2πrb/a)/μ(a) (193 kör mühür).
Soru: Dirichlet L(s,χ) adalarında (Not 5 sıfırları, 118b: β=χ₄, χ₃, χ₅ₑ, χ₈ₑ, χ₈ₒ
gerçel karakterler; 65-80k sertifikalı sıfır, t ≤ 4-5.5×10⁴) aynı mekanizma nasıl
değişir?

## Kalem cebiri (kaptan, veri-öncesi)

1. İlkel χ mod k, parite κ ∈ {0 (çift), 1 (tek)}. Sayma: N_χ(t) = N̄_χ(t) + S_χ(t),
   2πN̄_χ(t) = t·L_χ(t) − t + 2πC_χ, **L_χ(t) = log(k t/2π)**, C_χ = κ/4 − 1/8
   (çift −1/8, tek +1/8; Γ((s+κ)/2) Stirling'inden). S_χ(t) = −Σ_q χ(q) a_q sin(ω_q t)
   (gerçel χ), a_q = Λ(q)/(π√q log q). Orta noktalar: N(m_n) = n.
2. Jacobi–Anger: e^{−2πiS_χ} = Π_q Σ_k J_k(2π χ(q) a_q) e^{ikω_q t}; J_k(−x) = (−1)^k J_k(x),
   χ(q) = 0 ⇒ yalnız k = 0. Uydu n = a/b (sade): genlik ∝ χ(a)χ(b)·A_n^ζ
   (A_n^ζ zeta genliği); **k'yı bölen asal içeren n için uydu YOK.**
3. Durağan nokta: L_χ(t*) = ω' − log n ⇒ **t* = 2π q' b/(k a)** (iletken girer!);
   faz t* − 2πC_χ; durağan-faz çarpanı e^{−iπ/4} (φ'' < 0).
4. Çekirdek çizgi toplamı (çizgi ağırlığı χ(q')a'): Σ_{q'} χ(q') e(q'b/(ka)).
   gcd(a,k)=1 (2. maddeden) ⇒ ÇKT: = [χ̄(bā)τ(χ)]·[c_a(b k̄)] = χ(a)χ̄(b)·τ(χ)·μ(a).
   Sınıf ortalaması: τ(χ)μ(a)/(φ(k)φ(a)).
5. Birleşim (gerçel χ): genlik χ(a)χ(b) × toplam χ(a)χ̄(b) = χ(a)²|χ(b)|² = 1 ⇒
   **karakter işaretleri SADELEŞİR.** Kalan faz: τ(χ) (çift: √k; tek: i√k) ×
   e^{−i(2πC_χ + π/4)} (çift: e^{0}; tek: e^{−iπ/2}) ⇒ **her iki paritede GERÇEL
   ve zeta ile AYNI işaret.** Genlik çarpanı |τ|/φ(k) = **√k/φ(k)**.
6. Sonuç: adada uydu çekirdeği = zeta çekirdeği × √k/φ(k), k'yı bölen asal içeren
   uydular hariç. (φ(a) iptali açıkça: zeta'nın sınıf ortalaması μ(a)/φ(a); adanınki
   τ(χ)μ(a)/(φ(k)φ(a)) ⇒ oran τ(χ)/φ(k); fazla birlikte gerçel √k/φ(k).)

Üç rakip resim (hepsi ayırt edilir): (R1) iletkeni unutan: t* = 2πq'b/a ⇒ Σχ(q')
e(q'b/a) = 0 ⇒ adalarda iptal YOK; (R2) pariteyi unutan: tek karakterlerde faz
90° ⇒ gerçel kısım 0 ⇒ tek adalarda iptal YOK, çiftlerde var; (R3) kaptanın
resmi: her iki paritede var, √k/φ(k) ölçekli.

## Öngörüler

| ada | k | parite | √k/φ(k) | sönen uydular | yanan (güçlü) uydular |
|---|---|---|---|---|---|
| β (χ₄) | 4 | tek | 1.000 | log2, log6, log10, log(3/2) | ±log3, ±log5 |
| χ₃ | 3 | tek | 0.866 | log3, log6, log(3/2) | ±log2, ±log5, log10 |
| χ₅ₑ | 5 | çift | 0.559 | log5, log10 | ±log2, ±log3, log6 |
| χ₈ₑ | 8 | çift | 0.707 | log2, log6, log10 | ±log3, ±log5 |
| χ₈ₒ | 8 | tek | 0.707 | log2, log6, log10 | ±log3, ±log5 |

**A — yapı çarpanı gücü (ucuz, yalnız sıfırlar):** |Ĝ_mid(ω)|² bloğa-yerel
Δω = ω − L_χ,b ekseninde; uydu gücü, n'nin bütün asal çarpanları k'ya asal
ise VAR, değilse YOK (Bessel genliği; μ kuralı burada YOK — ör. zeta'da log4
güçte görünür). Kontrol: zeta (son, düşük) aynı ölçümle ±log2, ±log3, log4, log6.

**B — çekirdek (asıl sınav):** K̃_sat = Σ_q Σ_{q'∈uydu} c_q^{(q')}·conj(c_q^{öz}) /
Σ_q |c_q^{öz}|² (pencere çizgileri τ∈[0.45,0.86], χ ağırlıklı; öz-terime göre
normalize — karışım tanımından bağımsız). Zeta'da iptal ⇒ K̃ < 0. Öngörü: izinli
uydularda K̃_χ < 0 (her iki parite); yasak uydularda K̃_χ ≈ 0; L-eşli kıyasta
K̃_χ / K̃_ζ = √k/φ(k).

## Kapılar (K0 — ölçümden ÖNCE, hepsi zorunlu)

- **K0a C_χ MUTLAK KALİBRASYONU (teftiş (a)):** her ada ve zeta için
  ortalama[N̄_χ(m_n) − n] ≈ 0 (|·| < 0.05) — değilse (ör. ±0.5) orta-nokta işareti
  döner; DUR ve raporla.
- **K0b KARAKTER KİMLİĞİ (teftiş (i)):** sıfırları üreten Not-5 kodundan (118b/105b)
  her adanın karakter tablosu okunur; ilkel olduğu ve değerleri assert edilir
  (χ₈ₑ = (1,−1,−1,1), χ₈ₒ = (1,1,−1,−1) sırasıyla 1,3,5,7'de; χ₅ₑ = (·/5);
  χ₃, χ₄ standart). Uyuşmazsa DUR.
- **K0c ZETA KAPISI:** genelleştirilmiş kod (χ ≡ 1, k = 1) son ve düşük
  pencerelerde uydu ORANLARINI 188/190/193 ile ±%10 yeniden üretir; K̃_ζ'nin işareti
  bu kodla ölçülür ve REFERANS alınır (öngörüler işarete göre değil, "K̃_ζ ile
  aynı işaret" olarak ifade edilir — teftiş (g)).
- **K0d ÖN-KAYIT:** sha + damga (ONKAYIT_195.json).

## Hipotezler

- H-195A: güç uyduları asal-çarpan kuralına uyar (yasak olanlar < 2σ, izinli
  güçlüler > 3σ); ÖLÜM: yasak bir uydu ≥ 4σ.
- H-195B1 (varlık + işaret): χ₃'te log2, χ₅ₑ'de log2 ve log3, β'da log3: K̃, K̃_ζ ile
  AYNI işaretli, ≥ 3σ;
  yasak uydular (χ₃ log3, χ₅ₑ log5, β log2) |K̃| < 2σ. ÖLÜM: izinli güçlü uydu
  < 2σ (R1 lehine) ya da K̃_ζ'ye ZIT işaretli ≥ 3σ.
- H-195B2 (parite): tek adalar (β, χ₃, χ₈ₒ) ve çiftler (χ₅ₑ, χ₈ₑ) aynı işaret;
  χ₈ₑ ile χ₈ₒ aynı büyüklük (±%40). ÖLÜM: tek adalarda izinli uydular sıfır, çiftlerde
  var (R2).
- H-195B3 (ölçek, ikincil): L-eşli bantta (L ∈ [10.0, 10.4]) K̃_χ/K̃_ζ = √k/φ(k) ± %35.

## Ölçüm notları

L_χ(t) = log(k t/2π) ile blok-yerel Δω; adalarda L hızla değişir (t 200 → 5×10⁴) —
bloklar eşit ΔL genişlikli, yalnız üst bölge (L ≥ 8.5). L-eşli kıyas için zeta
referansı düşük pencere (190; L 10.0-10.84). Genelleştirilmiş kod ÖNCE zeta'da
kapılanır: son penceresinde uydu oranları 188/190/193 ile ±%10.
