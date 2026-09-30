# 203-b1 teftisi (log|zeta'(rho)| varyansi, oran sanisi + sifir-yogunlugu agirligi)

Kapsam: 203_oran_b1.py (+k3,k2). Hicbir dosya degistirilmedi; k2_b1'in tam zeta kosusu tekrarlanmadi.
Yalniz CUE kosulari (ucuz), zeta modunda tekil I(a,b) ve ic-integral J(a) degerleri, ve gercek sifirlar.

## Ozet

| # | Madde | Hukum |
|---|---|---|
| 1 | Yogunluk agirligi (Poisson, isaret, carpan 2, c->0) | PASS |
| 2 | Moment formulleri, T3 atamalari, gozlenebilir | PASS (gercek sifirlarla dogrudan dogrulandi) |
| 3 | Integral, simetri, serit, CUE-kesin | PASS (kucuk kusurla: Wc kesme yapayi ~c/a) |
| 4 | Wc=0.75 kesme etkisi | PASS: L=9.35'te +1.2e-4 (CUE), zeta ~ +0.9e-4 |
| 5 | Eksik terim / veriyle adil karsilastirma | PASS + 1 duzeltme notu (pencere genisligi/L=9.35 vs 9.3427; 0.3 sigma) |

## 1. Agirlik
Elle: Re zeta'/zeta(1/2+c+it) = pi*Sum_g P_c(t-g) - (1/2)log(t/2pi) + O(1/t^2)
=> Sum_g P_c(t-g) = (1/2pi)[L + 2 Re X(c)], artik isaret +, carpan 2. CUE analogu: |q|<1, Re[(1+q)/(1-q)] = (1-|q|^2)/|1-q|^2
=> Sum_n (1-r^2)/|1-q_n|^2 = N + 2 Re Z'/Z (Z=det(1-e^{-s}U)), ayni isaret.
Sayisal: t=60000, c=0.05/0.2, Odlyzko zeros1 (100k sifir) ile Poisson toplami vs (1/2pi)[L+2Re zeta'/zeta]
farki 5e-7..2e-6 (yanlis-isaret varyanti 2.77 vs 0.146: tamamen ayirt edilir). PASS.
E_t[X0]=0 => normalizasyon tam L. c->0: formuller c=0'da analitik (diger kaymalar >0 iken); ratios sanisinin
"Re c >> 1/L" duzgunlugu resmen kanitli degil ama asagidaki dogrudan sifir testleri bunu ampirik dogruluyor.
Ek ince test (bagimsiz): G(x) - L/x + L^2/2 -> 0 (x->0) zeta modunda (L=9.35, 11.66, 22.31; kalan ∝ x, ~1e-3 -> 0).
Bu tam olarak "Re zeta''/(2 zeta')(rho) = -theta'(gamma) = -L/2 deterministik" kimliginin E_z[X(a)] = G(a)/L
uzerinden bir sonucu: hem 1/L agirligi hem carpan 2 hem aritmetik sabitler dogrulaniyor.

## 2. Momentler
E_t[XX]=E_t[XXX]=0 (eslenik yok, yalniz n=1 terimi yasar) dogru. E_t[X(a)Xbar(b)Xbar(0)] = conj(E_t[X(b)X(0)Xbar(a)])
= T3(b,0;a) (T3 gercek: t->-t simetrisi). Kod I(a,b) = 1/2[Cov(X(a),Xbar(b)) + Cov(X(a),X(b))] = Cov(Re X(a), Re X(b)) (dogru).
kappa2 = Var(Re Y) = 1/2[Cov(Y,Ybar)+Re Cov(Y,Y)], -1_{a<1}/a sabiti kovaryansa girmez. log|Z'| = log|zeta'(rho)| (Z' ve zeta' modulu esit);
per-sifir -log L_n yalniz ortalamayi kaydirir; L-bagimli ortalamanin varyans katkisi Var_L(mu)= 9e-6/1.5e-6/1.2e-6 (W1/W2/W3) ihmal edilir.

BAGIMSIZ DOGRUDAN TEST (ratios sanisindan bagimsiz): Re X_rho(a) = Sum_{g'} a/(a^2+(g-g')^2) - L/2 (Poisson, kesin).
Yani I(a,b) = gercek sifirlar uzerinde Cov(S_a,S_b). Karsilastirma:
 * LMFDB C1 penceresinden 8e5 sifir, L=14.192: a,b in {0.1,0.3,1,3,6}/L, 15 hucre: |pull| <= 0.5, rms 0.27.
   Ornekler (data vs formul): (3,3)/L: 1.09118(137) vs 1.09166; (6,6): 0.43334(15) vs 0.43329; (1,1): 1.5634(110) vs 1.5681;
   (3,6): 0.63112(45) vs 0.63121. CUE(N=L) ayni hucrelerde %20-25 uzak (L=9.35'te (3,3): veri 0.4219, zeta-formul 0.4224, CUE 0.5266).
 * Odlyzko zeros1, L=9.35, a,b >= 1/L: (3,3): 0.42189(82) vs 0.42235; (5.6,5.6): 0.18434(16) vs 0.18454. (kucuk a: agir kuyruk, orneklem kovaryansi
   nadir yakin ciftlerle tanimli; jackknife SE hafife alir -> bu bolgede ampirik test bilgi vermez.)
=> G, T3 (kimlik + iki takas), aritmetik carpanlar, 1/L agirligi ve c=0 sinir degeri gercek sifirlarda %0.1-0.5 dogru.

## 3. Integral
* kappa2 = 2 Int da Int_a^inf db I(a,b): I simetrik (G(a+b), T3(a,0,b)+T3(b,0,a), gg, T3(a,b,0) hepsi (a<->b) simetrik) => dogru.
* CUE-kesin: k2_b1(N,cue=True,n=12,Wc=50) vs Sum_{j=1}^{N-1}[psi'(j+2) - psi'(j+1)/2] (tam sayi N) / analitik devam (Sum_j[f(j)-f(j+N-1)], tam olmayan N):
  N=2: 0.0724669768 vs 0.0724670334 (-5.7e-8); N=3: -6.9e-8; N=5: -7.6e-8; N=10: -8.0e-8; N=9.35, 10.59, 11.66, 22.3: -8.0e-8..-8.2e-8. (sistematik -8e-8: serit/kesme.)
  Ayrica bu sayilar formulu elle turetilen Palm degeriyle (KS momentleri, |Lambda|^2 agirligi) ayni.
* Ic integral J(a), L=11.66 zeta modunda kodun kendisi (Wc=0.75, n=12):
  a=0.0025/L: 8.253e-3 (J/a=38.5) | 0.005/L: 8.889e-3 (20.7) | 0.01/L: 1.3383e-2 (15.60) | 0.02/L: 2.4030e-2 (14.01) |
  0.05/L: 5.580e-2 (13.01) | 0.1/L: 0.10404 (12.13) | 0.3/L: 0.25069 (9.74).
  Yani KODUN J(a)'si a -> 0'da dogrusal-sifira gitmiyor (J/a artiyor). Neden: yapay iki kaynaklidan biri BASKIN: Wc kesmesi. b > Wc'de atilan terim
  F(a+b)K(-b,a+b), K'nin icinde -zp(u+x) = -zp(a) ~ 1/a kutbu tasiyor => J'ye c/a yapayi giriyor. Kanit: (J_Wc0.75 - J_Wc0.6)*a sabit (zeta: -4.5e-6, -4.5e-6, -4.4e-6, -4.2e-6
  a=0.005..0.05/L; CUE'de Wc=0.75'in Wc=50'ye gore yapayi *a = 2.27e-6, 2.25e-6, 2.21e-6, 2.10e-6). CUE'de mpmath (80 hane) gercek J: J/a = 16.86, 16.60, 16.31, 15.58, 14.59, 11.74
  (a=0.0025,0.01,0.02,0.05,0.1,0.3 /L) => GERCEK ic integral dogrusal-sifir, J'(0)=~17 (CUE). Zeta'da yapayi c(0.75)~0.72*2.27e-6=1.6e-6 alinirsa J_gercek/a ~ 13.5 (0.01-0.02/L), 12.9 (0.05/L),
  12.1 (0.1/L): duz ve dogrusal. Cift-duyarlik gurultusu ayrica: a=b civari I'da ~2.5e-10/a^2 (cift-duyarlik + asal kesmesi); a<0.005/L'de I(a,b) cop, ama katkisi ihmal.
  Serit: kod 2*J_code(amin)*amin/2 = 1.15e-5; dogru dogrusal 2*13.4*amin^2/2 = 1.0e-5 (fark 1.5e-6). a<0.1/L bolgesinin kappa2'ye toplam katkisi ~1e-3 (kappa2'nin %0.2'si) ise
  bu bolgedeki %20'lik hata <=2e-4, gozlenen (12/16/20 dugum: 5e-6) yakinsama ve CUE-kesin (8e-8) cok daha iyi.
  UYARI: serit 'a<amin' icin dogru; ama [amin, 0.1/L] araligindaki c/a yapayi log(a_hi/amin) ile buyuyor (bkz. madde 4).

## 4. Wc kesmesi
CUE, Wc=0.75 vs Wc=50 (fark = kesme sistematigi): N=9.35: +1.19e-4; 10.59: +4.4e-5; 14.19: +2.7e-6; 22.31: +7e-10.
Wc=0.6: 9.1e-4 / 2.3e-4 / 1.5e-5; Wc=1.0: 6.7e-6 / 1.5e-6.  Zeta'da Wc 0.6->0.75: 5.7e-4 (L=9.35), CUE'de 7.9e-4 => zeta orani ~0.72 => zeta Wc=0.75 hatasi ~ +0.9e-4 (L=9.35).
Etki tumuyle yukaridaki c/a yapayi: 2c*ln(a_hi/amin) ~ 2*1.3e-5*4.6 (c(N=9.35) ~ 2.27e-6*e^{0.75*2.3}); amin'e bagimlilik (CUE, N=9.35): amin=0.003/L: +1.5e-4; 0.01/L: +1.2e-4; 0.03/L: +0.9e-4; 0.1/L: +0.1e-4 (serit+yapay).
Sistematik olcek: <=1.5e-4 (L=9.35), <=5e-5 (L=10.6), 3e-6 (14.2), yok (L>=16.6). SE'ler 2.2e-3, 1.4e-3, 5.5e-4 => <0.1 sigma. PASS.

## 5. Eksik terim / karsilastirma
Isaret/carpan/terim hatasi bulunamadi; D uzerindeki asal kuyrugu (P>1e6, ~2e-4/L) kappa2'de ~1e-7.
Karsilastirmada bulunan (raporda YOK) ince nokta: W1-W3 pencereleri L'de genis (200a_tahmin.json: W1 t=18.7k..138k, Var(L)=0.26; W2 0.079; W3 0.094).
Tahmin tek L'de (9.35/10.59/11.66) hesaplanmis; oysa
 - W1 sifir-agirlikli Lbar1 = 9.3427 (9.35 degil): tahmin +2.9e-4 fazla;
 - kappa2(L) konkav (kappa2'' ~ -0.0035..-0.004): pencere-ortalamali tahmin = E_z[kappa2(L_n)] + Var_L(mu): W1: -4.6e-4, W2: -1.2e-4, W3: -1.2e-4.
Toplam duzeltme (yaklasik; kappa2(L)'ye a+b logL+c/L uyumu 5e-5 icinde): W1 0.46906 -> 0.46832 (sapma -0.85 sigma -> -0.51 sigma), W2 0.51901 -> 0.51888 (+0.71 -> +0.81),
W3 0.55843 -> 0.55826 (+0.15 -> +0.31). chi2/7: 10.93 -> 10.69 (p=0.153). Sonuc degismiyor. C pencereleri dar (Lbar1 = L kullanilan; 14.194254 vs 14.194224).
SE: 128-blok jackknife x kalibrasyon f (200c_analiz); x = log|Z'| ustel kuyruklu (p(delta)~delta^2 => log delta'da e^{3x}) => orneklem varyansi icin adil. C2 -2.9sigma (0.7063 vs 0.7078;
tek-yonlu p~0.002, 7 test icin ~0.03) formul hatasi belirtisi degil: kompanent kontrolleri (madde 2) C1'de %0.1-0.5; tavsiye: C2'nin 1e6 sifiri icin kappa2 blok-duyarliligi (M8c) rapora eklensin.
Uyari (yorum): "oran sanisi" burada onerme olarak degil, 'ampirik dogrulanmis kimlik + integral' olarak da desteklenebilir: I(a,b) = Cov(Re X(a), Re X(b)) zaten dogrudan sifirlardan hesaplanabiliyor.
Kalan acik: (i) ratios sanisinin c->0'da duzgunlugu teorik kanitsiz; (ii) a<0.1/L bolgesi ampirik olarak test edilemedi (agir kuyruk), yalniz CUE-kesin ve kucuk katki argumani.
