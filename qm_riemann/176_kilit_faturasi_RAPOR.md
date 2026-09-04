# 176 — KİLİDİN FATURASI: KARIŞTIRMA SINAVI
### (4 Eylül 2026, Opus tayfası — KALEM_KILIT_FATURASI_04EYL2026.md)

**Soru.** 175, gerçek ζ gazının `ΔM = +%5.92`'lik fazlasını E kanalında
ikiye ayırdı: **kesim %39 + kilit %61** (üç bağımsız kestirici — `f_b`,
`g_E`, `μ̂²_E` — %2.7 içinde aynı sayıyı verdi). Ama bu bir **ATIF**tı:
"kilit" adı, kesim ekseninin taşıyamadığı artığa verildi. 176 o atfı
**nedensel** olarak sınıyor: gerçek gazın merdiveninin **zarfı birebir
korunup yalnız çizgi FAZLARI karıştırılırsa** ve gaz sadakatli zincirle
sıfırdan yeniden kurulursa, defterin hangi satırı çöker?

Betikler (`176_configs/`) — `176g` dışında hepsi **ön-mühürlü**
(hiçbir 176 ölçümü var olmadan önce yazıldı); `176g` figürdür ve
ölçümden sonra yazılmıştır. `176f`, koşulmadan önce **F7'yi F3'e de
uygulayacak** biçimde düzeltildi (§7a; düzeltme hükmü zayıflattı).
Git'e **dokunulmadı**:

| betik | kapı | koşu | ham çıktı |
|---|---|---|---|
| `176a_onkayit.py` | **K2 ön-kaydı** (F0–F9 donduruldu) | 0.1 s | `176/ONKAYIT_K2.json` |
| `176_vekil_cekirdek.py` | faz taşıyan alan değerlendiricisi (çözücü YOK) | — | — |
| `176b_vekil_insa.py` | **K1 (inşa)** — VF1/VF2 vekil gazları | ~10 dk/gaz | `176/z_VF*.npy`, `176/insa_VF*.json` |
| `176c_olcum.py` | **K1 (ölçüm)** — 167 zinciri + çarpan defteri | ~6 dk/gaz | `167/C_VF*.json`, `176/G_VF*.json` |
| `176d_K1.py` | **K1 (girişim)** — 162/174b makinesi | ~2.5 dk/gaz | `176/K1_VF*.json` |
| `176e_K3dc.py` | **K3** — DC kaçağı defteri (174d makinesi) | ~1.5 dk/gaz | `176/K3_VF*.json` |
| `176f_K2_yuzlesme.py` | **K2 + K3** — dondurulmuş kural + ayrışım tablosu | 2 s | `176/K2.json` |
| `176g_figur.py` | figür | 25 s | `176_kilit_faturasi.png` |

Önbellek kökü `/private/tmp/claude-501/.../scratchpad/`; girdi
`128_odl_zeros6_2e6_zeros.npz` (gerçek sıfırlar), `172/G1.json`,
`174/K1_*.json`, `174/K3_*.json`; çıktı `176/`.

---

> ## TEK CÜMLELİK HÜKÜM
>
> **KİLİT ÖLÇÜLDÜ, MÜHÜR ÖLDÜ.** Zarfı `Hkeskin`'inkiyle **bit-bit
> aynı** (sha256 özdeş), yalnız çizgi fazları karıştırılmış ve
> sadakatli zincirle sıfırdan kurulmuş iki gaz, asal fazlarının
> kilidini bu programda **ilk kez nedensel olarak** görünür kıldı:
> `R_η` **1.26504 → 0.8085**, `Q_E` **−%68**, `μ̂²_E` **−%86**,
> `m3_çizgi` **−%77** — hepsi λ ailesinin de kesim ailesinin de çok
> dışında. **Ama ön-kayıtlı mühür (`F2`) kesin biçimde ÖLDÜ:** kilit
> sökülünce `g_E` düşmedi, **0.5837 → 0.7368 (+%26)** yükseldi;
> `M` de **+%64**. Sebep ölçüldü ve cebridir: `g_E ≡ 1 − Q_E/ρ_E` bir
> ORANDIR ve kilit **hem payı (−%68) hem paydayı (−%49)** taşır.
> ⇒ **`g_E` tekdüze bir kilit ölçeri değildir**, ve 175'in "%61 kilit"
> atfı bu para biriminde **sınanamaz**.
> Ayakta kalan tek fatura **θ**'dır: kilit sökülünce θ **0.8933 →
> 0.8494** düşüyor — gerçeğin fazlasının **%160'ı** — ve ödeme oranı
> `ω_θ = +0.716 ∈ (0,1]`. Ön-kayıtlı **beklenti (b) YANLIŞ ÇIKTI**
> (T-4'ün "Ç1 kilide kördür" savı öldü); ama iki tohumun saçılımı
> literal eşiğe ulaşmaya yetmedi ⇒ **H-F1b HÜKÜMSÜZ**, üçüncü tohuma
> borç. **H-F2 de öldü** (`kilit^{HF2} = 0.2499`, hedef 0.61) ve
> ölürken bir **sahte isabet** ortaya çıkardı: τ>0.70 açığının %61.5'i
> **genlik**, yalnız %25'i faz uyumudur.

---

## 0. ÖN-KAYIT (K2) — 4 Eylül 2026, 18:53:39 +03

`176a_onkayit.py` **hiçbir 176 ölçümü koşmadan önce** koştu; kendi
sha256'sını

```
6173943f0f4ef0488dced150f9e5a767e13c6bb24af2837527d2bd3149690244
```

`176/ONKAYIT_K2.json`'a yazdı ve dosya bir daha yazılmaz (varlık
kontrolü). *(Betik metni 175 bittikten sonra, 4 Eylül 16:35'te yazılmış
ve kotayla kesilen ilk başlatmada koşamamıştı; bu koşuda değiştirilmeden
koşuldu ve zaman damgası/sha256 o koşununkidir. Ondan sonra tek bir
karakteri değişmedi.)*

### 0.1 DÜRÜSTLÜK BEYANI (ön-kayıtta aynen duruyor)

1. **KÖRLÜK İDDİASI YOKTUR.** 174/175'in bütün sayıları (ΔM = +%5.92,
   kesim %39 / kilit %61, θ'nın %51.5'i, τ>0.70 faz uyumu %78, genlik
   %77, `Kov(çizgi,artık) = −0.0435`) bu tayfa tarafından okunmuştur.
   Ön-kayıt edilen şey **KURALDIR**.
2. **ÖN-KAYIT ANINDA ÖLÇÜLMEMİŞ OLAN:** faz-karıştırılmış vekil gazın
   **hiçbir niceliği**. Böyle bir gaz bu programda hiç kurulmamıştır.
   162'nin vekilleri **konum uzayında**, kurulmuş gazın ALANLARININ
   çizgi katsayılarını döndürür (`162_cekirdek.Taban162.vekil_yerel`:
   `rot = exp(i·U(0,2π))`, `c_eta·rot`, `c_C·rot` + süreklilikte
   Prichard–Theiler); 176'nınki **MERDİVENİN** fazını karıştırıp gazı
   sıfırdan çözer. İkisi farklı nesnelerdir ve rapor bunu böyle yazar.
3. **VEKİLİN ZARFI GERÇEĞİN ZARFIDIR.** Gerçek gazın nominal merdiveni
   `Hkeskin`'inkidir (λ = 1, keskin τ ≤ 1.00). Vekil bu `a_q` dizisini
   **bit-bit** korur, yalnız her çizgiye rastgele bir `φ_q` takar.
4. **ÖN-KAYIT SONRASI HİÇBİR EŞİK DEĞİŞTİRİLMEZ.** Ölümler kurtarmasız.

### 0.2 TÜRETİMLER (ölçümden önce, kalemle)

**T-1 VEKİL MERDİVENİ.**
`S_v(t) = −Σ_q A_q sin(ω_q t + φ_q)`, `φ_q ~ U(0,2π)` bağımsız,
`|A_q|` birebir `Hkeskin`'inki. Güç tayfı **özdeş** (rms S, rms S',
`Σa_qω_q` hepsi aynı); değişen tek şey çizgiler-arası **bağıl fazdır**.

**T-2 `R_η`'NİN KİLİT İÇERİĞİ.** 174-K1d'nin özdeşliği
`P ≡ Kov(η_çizgi, η)`, `R_η = P/Var(η)`,
`Var(η) = Var(η_çiz) + 2Kov(η_çiz,η_art) + Var(η_art)`;
`R_η > 1 ⟺ Kov(η_çizgi, η_artık) < 0` (Hkeskin: −0.05574). Fazlar
rastgele olunca bu kovaryans yapısal olarak kaybolur ⇒
`R_η(vekil) ≈ P/(P + Var_artık) = 0.6949` (bir NOKTA değil MERKEZ).

**T-3 DC KAÇAĞI KİLİDE KÖRDÜR (ters yönlü öngörü).**
175g'nin yasası `κ(ω_Q) = −π A_Q τ_Q cos(πτ_Q)`; faz taşıyıcısıyla
`κ(ω_Q) ≈ −π τ_Q A_Q cos(πτ_Q)·e^{−iφ_Q}` ve ölçülen `hp_Q ∝ b_Q e^{+iφ_Q}`
olduğundan `Δφ_Q := arg hp_Q + arg κ_Q` **φ_Q'dan bağımsızdır** ⇒ DC
kaçağı bir **öz-rezonanstır** ve faz karıştırmasına kördür.

**T-4 θ VE Ç1 KANALI.** 165'in Ç1 kanalı (`kanal1`, çift-iptal) `q₂=q₃`
eşleşmesiyle çalışır ve taşıyıcıdan başka hiçbir çizginin fazını
taşımaz ⇒ **Ç1 de kilide kördür**. Üçüncü momentin baskın kanalı Ç1 ise
**θ çökmez**. Ön-kayıtlı beklenti: **(b) θ çökmez** ⇒ H-F1b ölür.

**T-5 ÜÇÜNCÜ MOMENTLERİN RICE SIFIRI.** Rastgele fazlı çizgi alanının
bütün üçüncü momentleri özdeş sıfırdır ⇒ `m3_çizgi` ve `skew(η_çizgi)`
çökmelidir. `skew(ds)` için **ön-kayıt yoktur** (162-V4: konum uzayı
vekilinde büyümüştü).

### 0.3 DONDURULMUŞ KURAL (ölçüm buna göre yargılanır)

| # | kural | eşik / ölüm |
|---|---|---|
| **F0** | inşa + ölçüm kapıları | ilk-kök 300000/300000; maks\|F\| ≤ 1e−8; sıralılık TAM; **zarf farkı = 0.0 (TAM sıfır)**; `R_bant` ≥ 0.98 (lo ∈ [0.52,0.68]) |
| **F1** | `R_η(vekil)` | merkez **0.6949**, bant [0.60, 1.05]; **gerek** `Λ_R := R_η(Hk) − R_η(vekil) ≥ 0.0059398`; **ölüm** `R_η(vekil) ≥ 1.15` |
| **F2** | `g_E(vekil)` | bant [0.25, 0.57]; **gerek** `Λ_E := log g_E(Hk) − log g_E(vekil) ≥ 0.0149735` (= 0.61 × 0.0245468); **ölüm** `g_E(vekil) ≥ 0.583701` |
| **F3** | `θ(vekil)` | **(a)** ≤ **0.8656170** ⇒ H-F1b YAŞAR; **(b)** ≥ **0.8933338** ⇒ H-F1b ÖLÜR; ara: KISMİ. *Ön-kayıtlı beklenti: (b)* |
| **F4** | `Δlog M(vekil←Hk)` | bant [−1.00, −0.06]; **gerek** `Λ_M ≥ 0.0575416`; **ölüm** ≥ 0 |
| **F5** | Rice sıfırı | `\|m3_çizgi\| ≤ 0.040` (Hk: 0.0802) ve `\|skew(η_çizgi)\| ≤ 0.130` (Hk: 0.2501); `skew(ds)` için ön-kayıt YOK |
| **F6** | T-3 (ters yönlü) | `ort(E)(vekil) > 0` ve `\|ort(E)/0.0621390 − 1\| ≤ 0.30`; `⟨cosΔφ⟩ ≥ 0.90` (0.50–0.80 bantları) |
| **F7** | tohum saçılımı | ≥ 2 tohum; `\|y₁−y₂\| > \|y_ort − eşik\|` ise madde **HÜKÜMSÜZ** |
| **F8** | H-F2 kaba kalem | `kilit^{HF2} := Δ_faz/(Δ_amp+Δ_res+Δ_faz)`, τ>0.70; `\|·−0.61\| ≤ 0.10` ⇒ YAŞAR |
| **F9** | nihai ayrışım tablosu | `ΔlogM = Δlog g_E + 2Δlog g_X + Δlog θ`; ödeme oranı `ω_j := kilit_j / Δlog_j(Hk←vekil)`; `ω_j ∈ (0,1]` ⇒ fatura ödenebilir |

**MÜHÜR KURALI (dondurulmuş).**
F1 ✓ ∧ F2 ✓ ∧ F3(a) ∧ F4 ✓ ⇒ **MÜHÜR-TAM**;
F1 ✓ ∧ F2 ✓ ∧ F3(b) ⇒ **MÜHÜR-E** (cümle yalnız E kanalında mühürlenir,
θ "üçüncü eksen" olarak keskin bir negatifle açık kalır);
F2 ✗ ⇒ **ÖLÜM** (kilit faturası ödemiyor).

Mühürlenecek cümle:
> **"Gerçeğin üçüncü-moment fazlası, kesim şekli (%39, E) düşüldükten
> sonra asal fazlarının kilidinin faturasıdır."**

### 0.4 ÇAPALAR (ön-kayıt anında önbellekten okundu)

| nicelik | `Hkeskin` (ikiz) | `son` (gerçek) | `HA4` (erfc ikiz) |
|---|---|---|---|
| `R_η` (162 makinesi) | 1.26504 | **1.28829** | 1.29627 |
| `g_E` | 0.583701 | **0.598206** | 0.621589 |
| `g_X` | 0.704216 | 0.704733 | 0.716767 |
| `θ` | 0.893334 | **0.921938** | 0.723573 |
| `M = KALİB_u2` | 0.258590 | **0.273907** | 0.231073 |

`Δlog(son←Hk)`: `g_E` +0.0245468, `2g_X` +0.0014770, `θ` +0.0315178,
**`M` +0.0575416**.
`Δlog(erfc←Hk)`: `g_E` +0.0628900, `2g_X` +0.0353577, `θ` −0.2107590,
`M` −0.1125113.  `f(R_η) = +0.74447`.

---

## 1. K1 (İNŞA) — VEKİL GAZLAR ve ZARF-KORUNUMUNUN KANITI (`176b`)

### 1a. Ne yapıldı (ve neyin kopyalanmadığı)

`176_vekil_cekirdek.py` yalnızca **alan değerlendiricisini** taşır:

```
S_v(z)  = − Σ_q A_q sin(ω_q z + φ_q)
S_v'(z) = − Σ_q A_q ω_q cos(ω_q z + φ_q)
```

blok yapısı (`blok = 800` çizgi, `npt = 10000` nokta) `164_insa.S_ve_Sp`
ile birebir; tek fark `arg`'a `φ_q`'nun eklenmesidir. **Çözücüde tek
satır kopyalanmadı**: `164_insa.coz_sadakatli` (ızgara braketi + SIRALI
İLK-KÖK + korumalı Newton, `h = 0.015`, `nz = 300000`, `c = −½`) aynen
çağrıldı; yalnız onun kullandığı `S_par` / `SSp_par` isimleri koşudan
önce faz taşıyan sürümlerle değiştirildi.

### 1b. ZARF KORUNUMU — kanıt (`176b`, ön-mühürlü kapı G1)

| kanıt | ölçülen |
|---|---|
| `maks \|A_q(vekil) − A_q(Hkeskin)\|` | **0.0e+00** (TAM sıfır, makine hassasiyeti DEĞİL) |
| `maks \|ω_q(vekil) − ω_q(Hkeskin)\|` | **0.0e+00** |
| `sha256` (genlik dizisinin ham baytları), vekil | `45af6fa572bdd99c3a39c1ff75ccb8693038f733bf32d4ef7f921cae9e43c267` |
| `sha256`, `Hkeskin` | **aynı dize** |
| `Σ\|A\|` | 27.815371002 (iki gazda da) |
| `ΣA²` | 0.292663602 |
| `rms S` | 0.3825334 (175 §0.4'ün `Hkeskin` satırı: 0.38253) |
| `rms S'` | 1.8939673 |
| `Σ A_qω_q` | 259.8615 |
| çizgi sayısı | 15450 (τ ≤ 1.00) |

Ayrıca **değerlendiricinin kendisi** sınandı (kapı G2): `φ ≡ 0`
konduğunda faz taşıyan değerlendirici `164_insa.S_ve_Sp` ile
**maks\|ΔS\| = 0.0e+00** ve **maks\|ΔS'\| = 0.0e+00** veriyor (401
nokta). Yani vekil ile ikiz arasındaki *tek* fark `φ_q` dizisidir.

> **Bağımsız üçüncü kanıt:** ızgara aralığı `Σ|A|`'dan türetildiği için
> vekilin ızgarası `Hkeskin`'inkiyle **aynı nokta sayısını** verdi:
> `ng = 10 449 281` (164'ün `insa_Hkeskin.json`'undaki sayı: **10 449 281**).

### 1c. İNŞA KAPILARI (F0) — ölçülen

| kapı | eşik | `VF1` (tohum 1) | `VF2` (tohum 2) |
|---|---|---|---|
| ilk-kök hücre (benzersiz) | 300000 / 300000 | **300000 / 300000** ✓ | **300000 / 300000** ✓ |
| maks \|F\| | ≤ 1e−8 | **1.863e−09** (aşan tekne **0**) ✓ | **1.863e−09** (aşan tekne **0**) ✓ |
| sıralılık | TAM | **TAM**, min Δz = 0.091861 ✓ | **TAM**, min Δz = 0.093265 ✓ |
| zarf farkı | = 0.0 | **0.0e+00** ✓ | **0.0e+00** ✓ |
| φ≡0 sağlaması | = 0.0 | **0.0e+00 / 0.0e+00** ✓ | **0.0e+00 / 0.0e+00** ✓ |
| ızgara nokta sayısı | — | **10 449 281** (= Hkeskin) | **10 449 281** (= Hkeskin) |
| ikiye-bölme adımı | — | 26 (Hkeskin: 0) | 23 |
| `σ_ds` (z'den) | — | 0.44159 | 0.44273 |
| `L` (z'den) | — | 12.029593224 | 12.029593224 (Hkeskin: 12.029593242) |
| ΔG<0 kesri | — | 0.15605 | 0.15577 (Hkeskin: 0.19154) |
| süre | — | 10.0 dk | 10.6 dk |

> **İNŞA KAPILARININ HEPSİ (G1–G5) İKİ TOHUMDA DA GEÇTİ.** Vekiller
> `Hkeskin`'in **tam olarak aynı zarfına**, aynı ızgarasına ve aynı
> sadakat ölçütüne sahiptir; aralarındaki tek fark `φ_q` dizisidir.
> *(Ölçüm tarafındaki `R_bant` kapısı ayrıca §2c'de yargılanır — ve
> ISKALADI; oradaki sayı ve nedeni orada yazılıdır.)*

---

## 2. K1 (ÖLÇÜM) — VEKİLİN DEFTERİ (`176c`, `176d`)

`167_olcum.kos(ad, 0.40, 0.95, kule=0, duz=0)` — λ ve kesim ailesinin
gazlarıyla **aynı çağrı**. Vekilin künyesi `167_ortak.KUNYE`'de yoktur;
`kos` bu durumda λ = 1.00, τ_ust = 1.00, pen = None varsayar, yani
`b_nom = 2a_q sin(πτ_q)` **`Hkeskin`'inkiyle aynı formülden** çıkar.
`m`-kimliği iki vekilde de **0.0e+00**, `nline = 3425` (162 makinesi),
defter kalıntısı **1.27e−13 / 2.83e−14** (eşik 1e−12) ⇒ boru hattı
kimliği kırılmadı.

### 2a. ÖLÇÜLEN DEFTER (`176f`)

| gaz | `R_η` | `R_Ĉ` | `R_X` | `g_E` | `g_X` | `θ` | `M` |
|---|---|---|---|---|---|---|---|
| **`son`** (gerçek) | 1.28829 | 1.02599 | 1.09022 | 0.598206 | 0.704733 | **0.921938** | 0.273907 |
| **`Hkeskin`** (ikiz) | 1.26504 | 0.99086 | 1.07256 | 0.583701 | 0.704213 | 0.893334 | 0.258590 |
| `HA4` (erfc ikiz) | 1.29627 | 1.02270 | 1.06599 | 0.621589 | 0.716774 | 0.723573 | 0.231073 |
| **`VF1`** (tohum 1) | **0.80989** | 0.81215 | 0.80952 | **0.736939** | 0.813474 | **0.866671** | 0.422643 |
| **`VF2`** (tohum 2) | **0.80710** | 0.69177 | 0.78847 | **0.736692** | 0.833331 | **0.832136** | 0.425711 |

### 2b. ALTTAKİ NİCELİKLER — kilidin hangi kaba ne yaptığı

*(`π_E`, `ρ_E`, `Q_E`, `μ̂²_E`, `Var(E_mod)/Var(η)`: `176c` → `176/G_VF*.json`
ve `172/G1.json`; `Kov(η_çiz, η_art) ≡ P_η − Var(η_çizgi)`: `176d` →
`176/K1_VF*.json` ve `174/K1_*.json`.)*

| gaz | `π_E` | `ρ_E` | `Q_E` | `μ̂²_E` | `Var(E_mod)/Var(η)` | `Kov(η_çiz, η_art)` |
|---|---|---|---|---|---|---|
| `son` | 1.30385 | 2.17960 | 0.87575 | 0.040498 | 2.139 | −0.043527 |
| **`Hkeskin`** | 1.28375 | 2.19932 | 0.91558 | 0.056843 | 2.143 | −0.055735 |
| **`VF1`** | **0.83321** | **1.13064** | **0.29743** | **0.007998** | **1.123** | **−0.023755** |
| **`VF2`** | **0.83126** | **1.12837** | **0.29711** | **0.007858** | **1.121** | **−0.023716** |
| değişim (Hk → vekil) | **−%35** | **−%49** | **−%68** | **−%86** | 2.14 → 1.12 | **%43'e düştü** |

> **BULGU (K1-a) — KİLİT GERÇEKTİR ve DEVASADIR.** Faz karıştırması
> ikizin girişim oranını **1.26504 → 0.8085**'e indiriyor: bir λ ya da
> kesim değişikliğinin bütün ailesi (1.174 … 1.296) **bu aralığın
> yüzlerce σ dışındadır**. Aynı çöküş üç kanalda birden oluyor
> (`R_Ĉ` 0.99 → 0.81/0.69, `R_X` 1.07 → 0.81/0.79); projeksiyon kusuru
> `Q_E` **%68**, DC kaçağı `μ̂²_E` **%86** düşüyor. Yani "kilit" bir
> etiket değil, ölçülen bir madde.

### 2c. F0'IN `R_bant` KAPISI **ISKALADI** — ve ıskanın kendisi bir ölçümdür

Ön-kayıt "ölçüm tarafında `R_bant ≥ 0.98` (lo ∈ [0.52, 0.68])" diyordu.
Ölçülen:

| gaz | `R_bant` (9 bant, τ_eff 0.46 → 0.78) | hüküm penceresinde min |
|---|---|---|
| `son` | 1.1576 1.1781 1.1986 1.2150 1.2320 1.2485 1.2553 1.2578 1.2500 | 1.1986 |
| `Hkeskin` | 1.2237 1.2558 1.2857 1.3157 1.3424 1.3719 1.3851 1.4038 1.4019 | 1.2857 |
| `HA4` | 1.3254 1.3827 1.4527 1.5522 1.7204 2.0149 2.5853 3.8074 6.4831 | 1.4527 |
| **`VF1`** | 1.0063 0.9986 **0.9911 0.9772 0.9681 0.9402 0.9125** 0.8796 0.8435 | **0.9125** |
| **`VF2`** | 1.0102 0.9992 **0.9834 0.9757 0.9570 0.9337 0.8925** 0.8729 0.8312 | **0.8925** |

> **KAPI TUTMADI (kurtarma yok).** Ama ıskanın YÖNÜ, kapının koruduğu
> yönün tersidir: kapı "çizgiler yok olmuş olmasın" diye konmuştu;
> ölçüm ise **vekilin çizgilerinin NOMİNAL genlikte olduğunu**
> (`R_bant` ≈ 0.84–1.01) ve asıl **kilitli gazların nominali AŞTIĞINI**
> (`Hkeskin` %22–40, `son` %16–26, `HA4` %33–548) gösteriyor.
> `R_bant` bir sadakat ölçüsüdür: `√(Σ|c_ds|²_net / Σb_nom²)`. Faz
> karıştırılınca **1'e oturuyor.** Yani ön-kayıtlı kapı, üretmeye
> çalıştığı nesnenin ta kendisini eliyor.
> **Aday açıklama (iddia değil):** kilitli gazda artık (süreklilik)
> çizgilerle faz-uyumludur (`Kov(η_çiz, η_art) = −0.0557`), bu yüzden
> çizgi izdüşümü fazladan güç toplar; faz karıştırılınca o kovaryansın
> **%57'si kayboluyor** (−0.0557 → −0.0237) ve izdüşüm nominale iniyor.
> Sınav 177'ye borç.
>
> **Sonuç olarak:** inşa kapıları (G1–G5) tam geçtiği ve `R_bant`
> ıskası bir *inşa* kusuru değil bir *ölçüm bulgusu* olduğu için,
> aşağıdaki F-hükümleri **yazılmıştır**; ama hepsi bu ıskanın
> altında okunmalıdır ve hiçbir mühür bu kapıya dayanarak iddia
> edilmemiştir.

---

## 3. K2 — DONDURULMUŞ KURALIN YÜZLEŞMESİ (`176f`)

### 3a. F1 — `R_η(vekil)`: **✓ GEÇTİ (ezici)**

```
R_η(vekil) = 0.80989 / 0.80710  ⇒ ort 0.80850   tohum saçılımı 0.00280
bant [0.60, 1.05] ✓ İÇİNDE      ölüm eşiği 1.15: ✓ ölmedi
Λ_R = R_η(Hk) − R_η(vekil) = 1.26504 − 0.80850 = +0.4565467
gerek ≥ +0.0059398        ⇒  gereğin  76.9 KATI
```

Ön-kayıtlı merkez 0.6949, ölçülen 0.8085 (merkez %16 ıskaladı ama
bandın içi). Merkezin ıskası açıklanabilir ve ölçülüdür: T-2
`Kov(η_çiz, η_art) → 0` varsayıyordu; ölçüm **−0.0557 → −0.0237**
diyor (sıfır değil, yarısı).

### 3b. F2 — `g_E(vekil)`: **✗ ÖLDÜ (kesin)**

```
g_E(vekil) = 0.736939 / 0.736692 ⇒ ort 0.736815  (tohum saçılımı 0.000248)
ölüm eşiği: g_E(vekil) ≥ g_E(Hkeskin) = 0.583701   ⇒  0.7368 ≥ 0.5837  ✗ ÖLDÜ
bant [0.25, 0.57]: DIŞINDA (üstünde)
Λ_E = log g_E(Hk) − log g_E(vekil) = −0.2329484   gerek ≥ +0.0149735
```

**Ölüm kesindir** (tohum saçılımı 0.00025, eşiğe uzaklık 0.153) ve
mekanizması ölçülmüştür. `g_E ≡ 1 − Q_E/ρ_E` özdeşliğiyle:

| | `Q_E` | `ρ_E` | `g_E = 1 − Q/ρ` |
|---|---|---|---|
| `Hkeskin` | 0.91558 | 2.19932 | 0.583701 |
| vekil (ort) | 0.29727 | 1.12950 | 0.736816 |
| değişim | **−%68** | **−%49** | **+%26** |

> **HÜKÜM (F2).** Projeksiyon kusuru `Q_E` — `g_E`'nin ölçtüğü şeyin ta
> kendisi — kilidin malıdır ve faz karıştırılınca **%68 çöker**. Ama
> model gücü `ρ_E` de kilidin malıdır ve **%49 çöker**. `g_E` bir
> ORANDIR; payı da paydası da kilit taşıdığı için, kilit sökülünce
> `g_E` **yükselir**. ⇒ **`g_E` tekdüze bir kilit ölçeri DEĞİLDİR**, ve
> ön-kayıt onu öyle varsaymıştı. Fatura `g_E` para biriminde
> **ödenemez** (§4'te `ω_E = −0.064 < 0`).

### 3c. F3 — `θ(vekil)`: **(b) KESİN OLARAK DIŞLANDI, (a) TOHUM GÜRÜLTÜSÜNDE**

```
θ(vekil) = 0.866671 / 0.832136  ⇒ ort 0.849404   tohum saçılımı 0.034535
(a) ≤ 0.8656170 (H-F1b YAŞAR):  ort 0.8494 ≤ eşik  → dal A,
      ama eşiğe uzaklık 0.016213 < saçılım 0.034535 ⇒ **HÜKÜMSÜZ (F7)**
(b) ≥ 0.8933338 (H-F1b ÖLÜR):   eşiğe uzaklık 0.043930 > saçılım ⇒ KESİN
```

> **HÜKÜM (F3).** **ÖN-KAYITLI BEKLENTİ (b) YANLIŞ ÇIKTI.** T-4'ün
> savı — "üçüncü momentin baskın kanalı Ç1'dir, Ç1 kilide kördür,
> dolayısıyla θ çökmez" — **ölmüştür**: θ ikizin 0.893334'ünden
> 0.849404'e **düştü** ve düşüş, gerçek gazın ikizden fazlasının
> (`Δlog θ = +0.031518`) **%160'ı** kadar (`Δlog θ(Hk←vekil) =
> +0.050633`). Yani **θ'nın fazlası kilit para biriminde ÖDENEBİLİR**
> (§4: `ω_θ = +0.716`). Ama (a) dalının literal eşiğine tohum
> gürültüsü içinde ulaşılamadı (tohumlar 0.8667 ve 0.8321; ikisi arası
> fark eşiğe uzaklığın iki katı) ⇒ **H-F1b ne yaşadı ne öldü;
> HÜKÜMSÜZ.** Üçüncü bir tohum bu maddeyi tek başına karara bağlar.

### 3d. F4 — `Δlog M(vekil←Hk)`: **✗ ÖLDÜ**

```
M(vekil) = 0.422643 / 0.425711    M(Hkeskin) = 0.258590
Δlog M(vekil←Hk) = +0.4912825 / +0.4985151 ⇒ ort +0.4948988 (saçılım 0.0072)
ölüm: Δlog M ≥ 0  ⇒  ✗ ÖLDÜ ;  bant [−1.00, −0.06] DIŞINDA
```

Sebep §3b'nin aynısıdır: `M = g_E g_X² θ` ve `g_E` (+%26) ile `g_X`
(+%17) yükselişi, θ'nın (−%5) düşüşünü ezici biçimde yeniyor.

### 3e. F5 — RICE SIFIRI: **✓✓ İKİ MADDE DE GEÇTİ**

| gaz | `m3` | **`m3_çizgi`** | **`skew(η_çizgi)`** | `skew(e1)` | `skew(x1)` | `skew(ds)` |
|---|---|---|---|---|---|---|
| `son` | +0.27675 | −0.06575 | −0.22159 | −0.15129 | +0.28174 | +0.46806 |
| `Hkeskin` | +0.24133 | −0.08022 | −0.25010 | −0.19100 | +0.22884 | +0.31850 |
| `HA4` | +0.28313 | −0.03277 | −0.14736 | −0.25369 | +0.29694 | +0.58796 |
| **`VF1`** | +0.18359 | **−0.02227** | **−0.11275** | −0.38024 | +0.16110 | +0.45381 |
| **`VF2`** | +0.18755 | **−0.01533** | **−0.11024** | −0.37831 | +0.21531 | +0.48887 |

`\|m3_çizgi\|` ort **0.01880** ≤ 0.040 ✓ (saçılım 0.0069);
`\|skew(η_çizgi)\|` ort **0.11149** ≤ 0.130 ✓ (saçılım 0.0025).
İkisi de kesin (saçılım ≪ eşiğe uzaklık).

> **HÜKÜM (F5).** T-5 doğrulandı: **çizgi alanının üçüncü momentleri
> faz karıştırmasıyla çöküyor** (`m3_çizgi` −%77, `skew(η_çiz)` −%55).
> Buna karşılık **`skew(ds)` çökmedi, BÜYÜDÜ** (0.3185 → 0.4713 ort) —
> ön-kayıt bunu zaten "sınanmayacak" ilan etmişti (162-V4) ve ölçüm
> onu doğruladı: öz-tutarlılığın doğrusal-olmayanlığı **fazdan bağımsız
> çarpıklık üretir**. `skew(e1)` de büyüdü (−0.191 → −0.379).

### 3f. F6 — T-3 (DC kaçağı kilide kör mü?): **✗ T-3 ÖLDÜ, ama yarısı ayakta**

| gaz | `ort(E)` | `μ̂²_E` | `⟨cosΔφ⟩` (0.40–0.50 / 0.50–0.60 / 0.60–0.70 / 0.70–0.80 / 0.80–0.95) |
|---|---|---|---|
| `son` | +0.046900 | 0.040498 | −0.849 / +1.000 / +0.999 / +0.988 / +0.604 |
| `Hkeskin` | +0.062139 | 0.056843 | −0.912 / +1.000 / +0.999 / +0.992 / +0.772 |
| `HA4` | +0.031399 | 0.015656 | −0.749 / +1.000 / +0.998 / +0.941 / −0.191 |
| **`VF1`** | **+0.027005** | 0.007998 | −0.935 / **+0.934 / +0.958 / +0.915** / +0.451 |
| **`VF2`** | **+0.026771** | 0.007858 | −0.934 / **+0.946 / +0.961 / +0.919** / +0.460 |

```
ort(E) > 0                    ✓  (T-3'ün işaret öngörüsü tuttu)
⟨cosΔφ⟩ ≥ 0.90 (0.50–0.80)    ✓  (min 0.9151; T-3'ün FAZ öngörüsü tuttu)
|ort(E)/0.0621390 − 1| ≤ 0.30 ✗  (ölçülen 0.5673 — %57 düşüş)
⇒ F6 ✗  ⇒  T-3 ÖLDÜ
```

> **HÜKÜM (F6).** Ön-kayıt bu ölümün anlamını önceden yazmıştı:
> *"Tutmazsa T-3 ÖLÜR (ve o zaman kilit DC kaçağını da taşıyor demektir
> — H-F1 için LEHTE bir sürpriz)."* Ölçüm tam olarak budur.
> T-3'ün **faz** yarısı ayakta: `Δφ_Q = arg hp_Q + arg κ_Q` gerçekten
> `φ_Q`'dan bağımsızdır (`⟨cosΔφ⟩` 0.50–0.80'de 0.91–0.96'da kalıyor,
> ikizinki 0.99–1.00). T-3'ün **genlik** yarısı öldü: `ort(E)` %57,
> `μ̂²_E` %86 düşüyor — çünkü hem `Σ|hp|` hem `|κ|` kilidin malıdır.
> Ve en üst bantta (τ > 0.80) faz uyumu da çöküyor: 0.772 (ikiz)
> → 0.604 (gerçek) → **0.455 (vekil)** — **gerçek gaz tam ikizle vekil
> arasındadır** ve vekile ikizden daha yakındır.
> *(174d'nin kendi ön-mühürleri Ö1 ve Ö3 vekilde de öldü —
> `ort(E) > 0`, en büyük 10 çizginin payı %6.7 — gerçek gazlardaki gibi;
> kayda geçer.)*

### 3g. F8 — H-F2 (kilit payının ilk-ilke öngörüsü): **✗ ÖLDÜ**

Dondurulmuş özdeş ayrışım (`ξ_b ≡ amp_b × res_b × faz_b`), τ > 0.70,
8325 çizgi:

| çarpan | `son` | `Hkeskin` | `Δlog` | açığın payı |
|---|---|---|---|---|
| `amp = Σ\|hp\|` | 7.016578 | 9.024681 | **−0.251688** | **%61.5** |
| `res = Σ\|hp\|\|κ\|/Σ\|hp\|` | 4.515989e−03 | 4.772317e−03 | −0.055208 | %13.5 |
| `faz = ξ/Σ\|hp\|\|κ\|` | 0.792881 | 0.878215 | **−0.102218** | **%25.0** |
| **`Σξ(τ>0.70)`** | **+2.512385e−02** | **+3.782353e−02** | **−0.409114** | %100 |

```
kilit^HF2 = Δ_faz / (Δ_amp + Δ_res + Δ_faz) = 0.2499
hedef 0.61 ± 0.10   ⇒  |0.2499 − 0.61| = 0.3601   ⇒  ✗ H-F2 ÖLDÜ
```

> **HÜKÜM (F8) — ve ölümün içindeki uyarı.** Kaba kalem **çok yanlış**
> değil, **tam tersi** yanlıştır: %61 gerçekten defterde vardır, ama
> **FAZ payında değil GENLİK payındadır** (`Δ_amp` payı **%61.5**,
> hedef 0.61'in 0.005'i içinde). Yani "τ>0.70 faz uyumu %78, genlik
> %77" metriklerinden `0.77 × 0.78` gibi bir kombinasyon kurmak
> ölçülen %61'i **rastlantısal olarak** yeniden üretir ve bu bir
> **sahte isabettir**: gerçek ayrışımda faz yalnızca **%25** taşır.
> Bant bant: τ 0.70–0.80'de kilit payı **0.0123** (faz neredeyse hiç
> katkı vermiyor), τ 0.80–0.95'te **0.4314**. Yani H-F2'nin tek bir
> sayısı yoktur.
> ⇒ **175'in "%61 kilit" atfı, ölçülü kilit metriklerinden
> parametresiz TÜRETİLEMEZ.** Aday öldü.

---

## 4. K3 — NİHAİ AYRIŞIM TABLOSU (`176f`)

**Dürüstlük notu (kayda geçiyor).** F9'un ön-kayıttaki *literal* metni
`kesim_j := f_j·Δ_j(erfc←Hk)`, `f_j := Δ_j(son←Hk)/Δ_j(erfc←Hk)`
**dejeneredir** (o zaman `kesim_j ≡ Δ_j(son←Hk)` ve `kilit_j ≡ 0`).
Uygulanan okuma, **F1 ve F2'nin SAYISAL eşiklerinde zaten yazılı
olandır**:
`kesim_j = f_j·Δ_j(son←Hk)`, `kilit_j = (1−f_j)·Δ_j(son←Hk)`
(gerçekten: `0.0149735 = 0.61 × 0.0245468` ve
`0.0059398 = 0.2555 × 0.0232450`). **Hiçbir ölüm eşiği
değiştirilmemiştir.**

Özdeşlik `Δlog M ≡ Δlog g_E + 2Δlog g_X + Δlog θ` (kalıntı **4.65e−16**;
vekil sütununda **2.2e−16 / 5.6e−17**):

| çarpan | `Δ(son←Hk)` | `Δ(erfc←Hk)` | `f_j` | **KESİM_j** | **KİLİT_j** | `Δ(Hk←vekil)` = TOPLAM kilit içeriği | **ödeme oranı ω_j** |
|---|---|---|---|---|---|---|---|
| `g_E` | +0.024547 | +0.062890 | **+0.3903** | +0.009581 | +0.014966 | **−0.232948** (−0.23312 / −0.23278) | **−0.064** ✗ |
| `g_X²` (2·g_X) | +0.001477 | +0.035358 | +0.0418 | +0.000062 | +0.001415 | **−0.312583** (−0.28847 / −0.33670) | **−0.005** ✗ |
| **`θ`** | +0.031518 | −0.210759 | **−0.1495** | −0.004713 | **+0.036231** | **+0.050633** (+0.03030 / +0.07096) | **+0.716** ✓ |
| **`M`** | **+0.057542** | −0.112511 | −0.5114 | −0.029429 | +0.086970 | −0.494899 | −0.176 ✗ |
| `R_η` (doğrusal, 162) | +0.023245 | +0.031223 | **+0.7445** | +0.017305 | +0.005940 | **+0.456547** | **+0.013** ✓ |

### 4a. ΔM = KESİM + KİLİT — E-payı ve θ-payı ayrı satır (parametresiz)

| | **KESİM** | **KİLİT** | toplam |
|---|---|---|---|
| **E kanalı** (`g_E`, `g_X²`) | +0.009643 | +0.016381 | +0.026024 |
| **θ kanalı** | −0.004713 | +0.036231 | +0.031518 |
| **TOPLAM (M)** | **+0.004929** | **+0.052612** | **+0.057542** |
| ΔM'nin yüzdesi | **+%8.6** | **+%91.4** | %100 |

E kanalının kendi içinde: kesim **%37.1**, kilit **%62.9** — 175'in
`g_E`'den okuduğu **%39 / %61** ile tutarlı (fark `2g_X`'in katkısı).

### 4b. FATURA ÖDENEBİLİR Mİ? (ön-kayıtlı ölçüt `ω_j ∈ (0,1]`)

> * **`θ`: ✓ ÖDENEBİLİR.** `ω_θ = +0.716`. Kilidin θ'da taşıdığı toplam
>   (+0.0506) gerçeğin fazlasına atfedilen kilit payından (+0.0362)
>   **büyüktür**, ve işareti doğrudur. *(Tohum saçılımı büyük:
>   ω_θ tohum başına **1.196 / 0.511**; üçüncü bir tohum gerekir.)*
> * **`R_η`: ✓ ÖDENEBİLİR, ezici bollukla.** `ω_R = +0.013` — kilidin
>   toplam içeriği faturanın **77 katı**.
> * **`g_E`, `g_X²`, `M`: ✗ ÖDENEMEZ — İŞARET TERS.** `ω < 0`.
>   Kilit sökülünce bu çarpanlar gerçeğin gittiği yönün **tersine**
>   gidiyor. Ön-kayıtlı ölçütün kendi diliyle: **kilit bu çarpanları
>   taşıyamaz.**

---

## 5. MÜHÜR — ve ÖLÜM

```
F0  ✗ (R_bant kapısı ıskaladı; inşa kapıları G1–G5 TAM)
F1  ✓ (Λ_R = +0.4565, gereğin 76.9 katı)
F2  ✗ ÖLDÜ (kesin)
F3  dal (b) KESİN olarak dışlandı; dal (a) tohum gürültüsünde ⇒ HÜKÜMSÜZ
F4  ✗ ÖLDÜ
F5  ✓✓
F6  ✗ T-3 ÖLDÜ (faz yarısı ayakta)
F8  ✗ H-F2 ÖLDÜ
```

Ön-kayıtlı mühür kuralı: **F2 ✗ ⇒ ÖLÜM.**

> ## MÜHÜR VURULMADI — **ÖLÜM**
>
> Mühürlenecek cümle —
> *"Gerçeğin üçüncü-moment fazlası, kesim şekli (%39, E) düşüldükten
> sonra asal fazlarının kilidinin faturasıdır"* —
> **MÜHÜRLENMEDİ.** Ön-kayıtlı ölçüt (`F2`) açık ve kesin biçimde
> tutmadı: faz kilidi sökülünce `g_E` düşmedi, **%26 yükseldi**;
> `M` de **%64 yükseldi**. Kurtarma yapılmadı.

---

## 6. ÖLÜMÜN BİLGİSİ (ölçülmüş, iddia değil)

1. **KİLİT VAR ve DEVASA.** Faz karıştırması `R_η`'yı 1.265 → 0.808,
   `Q_E`'yi %68, `μ̂²_E`'yi %86, `m3_çizgi`'yi %77 siliyor. Bu çöküşler
   λ ailesinin (1.174–1.296) ve kesim ailesinin (1.221–1.296) hiçbir
   üyesinin yaklaşamadığı bir yerdedir. **"Kilit" bir etiket değildir.**
2. **AMA ΔM DEFTERİNİN PARA BİRİMİ KİLİDİ ÖLÇMÜYOR.** `g_E ≡ 1 − Q_E/ρ_E`
   bir orandır; kilit hem `Q_E`'yi (−%68) hem `ρ_E`'yi (−%49) siler ve
   oran **yükselir**. Aynı şey `M` için de geçerlidir. ⇒ 175'in
   "E kanalında %61 kilit" atfı, **176'nın nedensel sınavıyla
   doğrulanamaz**; çünkü sınavın ölçtüğü yön ile defterin para biriminin
   yönü **ters**.
3. **TEK İSTİSNA θ'DIR — ve orada fatura ödeniyor.** Kilit sökülünce θ
   ikizin altına iniyor (0.8933 → 0.8494) ve düşüş, gerçeğin fazlasının
   %160'ı kadar; `ω_θ = 0.716 ∈ (0,1]`. **ΔM'nin %51.5'lik θ payı için
   "kilidin faturası" okuması ayakta kalan tek okumadır** — ama tohum
   saçılımı (ω_θ: 1.196 / 0.511) onu tek başına mühürleyemiyor.
4. **H-F2 ÖLDÜ ve bir SAHTE İSABET uyarısı bıraktı.** τ>0.70 açığının
   **%61.5'i genlik**, %13.5'i κ, yalnız **%25'i faz uyumudur**.
   "Genlik %77 × faz %78" tipi bir kaba kalem 0.61'i **rastlantıyla**
   yeniden üretir; gerçek ayrışım bunun tersini söylüyor.
5. **YENİ, ÖLÇÜLMÜŞ BİR OLGU: `R_bant` fazlası kilidin malıdır.**
   Nominal çizgi genliğine göre sadakat oranı kilitli gazlarda
   **1.16–1.40** (erfc-ikizde 6.48'e kadar), faz karıştırılınca
   **0.84–1.01**. Yani kurulmuş gazların "nominali aşan" çizgi gücü bir
   çözücü kusuru değil, **fazların kilidinin ürünüdür**. Bu, 167/169'un
   `R_bant ≥ 0.98` filtresinin ne ölçtüğünü yeniden tanımlar.
6. **T-3'ün YARISI AYAKTA.** DC kaçağının **bağıl fazı** kilide kördür
   (`⟨cosΔφ⟩` 0.50–0.80'de 0.91–0.96), **genliği** değildir. Ve en üst
   bantta (τ>0.80) faz uyumu sıralaması **ikiz 0.772 > gerçek 0.604 >
   vekil 0.455**: gerçek gaz, ikizden çok vekile yakındır.

---

## 7. DENETİM

| | sınav | sonuç |
|---|---|---|
| **D1** | zarf birebir mi? | `maks\|A_q(vekil) − A_q(Hkeskin)\| = **0.0e+00**`, aynı `sha256`, `Σ\|A\|`/`ΣA²`/`rms S`/`rms S'`/`Σa_qω_q` özdeş (`176b`) |
| **D2** | değerlendirici doğru mu? | `φ ≡ 0`'da `164_insa.S_ve_Sp` ile maks fark **0.0e+00** (S ve S') (`176b`) |
| **D3** | çözücü aynı mı? | `164_insa.coz_sadakatli` **aynen çağrıldı**; ızgara nokta sayısı `Hkeskin`'inkiyle **birebir aynı** (10 449 281) |
| **D4** | inşa kapıları | ilk-kök **300000/300000**, maks\|F\| **1.863e−09** (aşan 0), sıralılık **TAM** — iki tohumda da |
| **D5** | 162 makinesi bit düzeyinde mi? | `m`-kimliği **0.0e+00**, `nline = 3425`, `L = 12.0296` (`176d`) |
| **D6** | bant defteri özdeş mi? | bağıl kalıntı **1.27e−13 / 2.83e−14** (eşik 1e−12) |
| **D7** | DC kaçağı özdeşliği | `Σξ`'den `μ̂²`, `176c`'nin defteriyle **9.3e−14 / 5.5e−13** |
| **D8** | ΔM defteri özdeş mi? | `Δlog M − (Δlog g_E + 2Δlog g_X + Δlog θ)` = **4.65e−16**; vekil sütununda **2.2e−16 / 5.6e−17** |
| **D9** | ön-kayıt sonradan değişti mi? | `176/ONKAYIT_K2.json` bir kez yazıldı; betik varlık kontrolüyle üzerine yazmayı reddediyor; sha256 §0'da |
| **D10** | git | **dokunulmadı** |

### 7a. Dürüstlük notları (kayda geçiyor)

* **`176a` bu koşuda yazılmadı.** Betik, kotayla kesilen ilk başlatmada
  (4 Eylül 16:35) yazılmış, ama **koşamamıştı** (`176/` dizini boştu).
  Bu koşuda **tek karakteri değiştirilmeden** koşuldu (18:53:39) ve
  raporun sha256'sı o koşununkidir. Hiçbir 176 ölçümü ondan önce
  koşmadı.
* **F0'ın `R_bant` alt-kapısı ISKALADI** (§2c) ve kurtarılmadı. Ölçüme
  devam edildi, çünkü `R_bant` bir *inşa* kapısı değil bir *ölçüm
  çıktısıdır*: onu bilmek için ölçmek gerekir. Beş inşa kapısının
  (G1–G5) hepsi tam geçmiştir. Hiçbir mühür bu kapıya dayandırılmadı.
* **F9'un literal metni dejeneredir** (§4) ve F1/F2'nin sayısal
  eşiklerinde zaten yazılı olan okuma uygulandı. Ölüm eşikleri
  değişmedi.
* **F3, F7 kuralıyla HÜKÜMSÜZ ilan edildi** (§3c). İlk yazımda `176f`
  F7'yi yalnız F1/F2/F4'e uyguluyordu; **koşudan önce** F3'e de
  uygulanacak biçimde düzeltildi (dondurulmuş F7'nin gereği).
  Düzeltme hükmü **zayıflattı**, güçlendirmedi: aksi hâlde dal (a)
  "YAŞADI" yazılacaktı.
* **Vekiller `kule = 0, duz = 0` ile ölçüldü** — λ/kesim ailesinin
  gazlarıyla aynı çağrı. `Hkeskin` `duz = 1` ile ölçülmüştür; bu bayrak
  yalnız `Ps2d`/`C4`'ü etkiler, `g_E`/`g_X`/`θ`/`M`'yi **etkilemez**
  (`166_T1.bant_agg` bayraktan bağımsızdır).
* **`c_gaz` vekillerde tek bantlıdır** (0.9256 / 0.9681), çünkü 167'nin
  hüküm filtresi `R_bant ≥ 0.98` istiyor ve vekilde yalnız en alt bant
  geçiyor. Bu sayı hiçbir hükümde kullanılmadı; bilgi olarak yazılır.
* **`174d`'nin kendi ön-mühürleri Ö1 ve Ö3 vekilde de öldü**
  (`ort(E) > 0`, en büyük 10 çizginin payı %6.7) — gerçek ve sentetik
  gazlarda olduğu gibi. Ö2 (bağımsız yol) ✓.
* **Boşa giden koşu yok**; smoke sınavı (`VFTEST`, 3000 tekne) koşuldu
  ve çıktıları silindi.

---

# HÜKÜM

## (i) ÖN-KAYITLI MÜHÜR **ÖLDÜ** — kurtarma yok

> KALEM'in cümlesi (*"Gerçeğin üçüncü-moment fazlası, kesim şekli
> (%39, E) düşüldükten sonra asal fazlarının kilidinin faturasıdır"*)
> **mühürlenmedi**. Ön-kayıtlı `F2` kesin biçimde tutmadı: faz kilidi
> sökülünce `g_E` **0.5837 → 0.7368** (+%26) ve `M` **0.2586 → 0.4241**
> (+%64) — ikisi de gerçeğin gittiği yönün **tersine**.

## (ii) AMA KİLİT ÖLÇÜLDÜ — ve devasa

> Zarfı bit-bit korunmuş, yalnız çizgi fazları karıştırılmış ve
> sadakatli zincirle sıfırdan kurulmuş iki gazda: `R_η` **1.26504 →
> 0.8085**, `Q_E` **−%68**, `μ̂²_E` **−%86**, `m3_çizgi` **−%77**,
> `skew(η_çizgi)` **−%55**. Bunlar λ ailesinin de kesim ailesinin de
> tamamen dışındadır. **Asal fazlarının kilidi, bu programda ilk kez
> NEDENSEL olarak ölçülmüştür.**

## (iii) ÖLÜMÜN SEBEBİ: `g_E` BİR KİLİT ÖLÇERİ DEĞİL, BİR ORANDIR

> `g_E ≡ 1 − Q_E/ρ_E` özdeşliğinde kilit **hem payı hem paydayı**
> taşır (`Q_E` −%68, `ρ_E` −%49) ve oran **yükselir**. Ön-kayıt
> `g_E`'yi tekdüze bir kilit ölçeri saymıştı; **değildir**. Bu, 175'in
> "%61 kilit" atfını çürütmez — onu **sınanamaz** kılar: 176'nın
> nedensel sınavı `g_E` para biriminde **ters işaretli** yanıt veriyor.

## (iv) TEK ÖDENEN FATURA: θ

> `θ` faz karıştırmasıyla **düşüyor** (0.893334 → 0.849404) ve düşüş,
> gerçek gazın fazlasının **%160'ı**; ödeme oranı `ω_θ = +0.716 ∈ (0,1]`.
> Ön-kayıtlı **beklenti (b) YANLIŞ ÇIKTI**: T-4'ün "Ç1 kilide kördür,
> θ çökmez" savı **öldü**. Ama (a) dalının literal eşiğine tohum
> gürültüsü içinde ulaşılamadı (θ: 0.8667 / 0.8321; saçılım 0.0345 >
> eşiğe uzaklık 0.0162) ⇒ **H-F1b HÜKÜMSÜZ**, üçüncü bir tohuma borç.
> **ΔM'nin %51.5'lik θ payı için "kilidin faturası" okuması ayakta
> kalan tek okumadır.**

## (v) H-F2 ÖLDÜ — ve bir SAHTE İSABET ortaya çıkardı

> `kilit^{HF2} = 0.2499` (hedef 0.61 ± 0.10). τ>0.70 açığının
> **%61.5'i GENLİK** payıdır, %13.5'i κ, yalnız **%25'i faz uyumu**.
> "Genlik %77 × faz %78" tipi kaba kalem 0.61'i **rastlantıyla**
> üretir; parametresiz ayrışım tersini söylüyor. **175'in %61'i ölçülü
> kilit metriklerinden türetilemez.**

## (vi) YENİ OLGU: NOMİNALİ AŞAN ÇİZGİ GÜCÜ KİLİDİN ÜRÜNÜDÜR

> `R_bant` (sadakat oranı, `√(Σ|c_ds|²_net/Σb_nom²)`): `Hkeskin`
> **1.22–1.40**, `son` **1.16–1.26**, `HA4` **1.33–6.48**; faz
> karıştırılınca **0.84–1.01**. Kurulmuş gazların nominali aşan çizgi
> gücü bir çözücü kusuru değil, **fazların kilidinin ürünüdür** — ve bu
> yüzden ön-kayıtlı `R_bant ≥ 0.98` kapısı, sınamak için üretilen
> nesneyi eliyor (§2c; ıska kurtarılmadı).

---

## Sıradaki adım (bu ölçümün işaret ettiği)

1. **ÜÇÜNCÜ TOHUM.** F3 tek başına tohum gürültüsünde kaldı
   (θ: 0.8667 / 0.8321). İki inşa daha (~10 dk/gaz + ~9 dk ölçüm)
   H-F1b'yi yaşat-ya-öldür kararına bağlar. **En ucuz, en yüksek
   getirili adım.**
2. **DOĞRU KİLİT PARA BİRİMİ.** `g_E` ters işaretli. Ama `Q_E`
   (−%68) ve `μ̂²_E` (−%86) tekdüze ve devasa çöküyor. 177 için
   ön-kayıtlı yeni bir fatura: gerçeğin `Q_E` fazlasının (`Hkeskin`
   0.91558 → `son` 0.87575, **−%4.3**) kilit payı, `ω_Q` ile ödenebilir
   mi? Bu, ΔM'yi `g_E` yerine `(Q_E, ρ_E)` çiftinde ayrıştırmayı
   gerektirir — ve ikisi de ölçülmüştür.
3. **`R_bant` FAZLASININ MEKANİZMASI** (§2c, vi). Aday: çizgi ↔ artık
   faz uyumu (kovaryans −0.0557 → −0.0237). Doğrudan sınav: vekilde
   `Kov(η_çiz, η_art)`'ın bant bant profili ile `R_bant` profilinin
   karşılaştırılması — veri **zaten diskte** (`176/K1_VF*.json`,
   `167/C_VF*.json`), yeni koşu gerekmez.
4. **KESİM × KİLİT ÇAPRAZ GAZI.** 176 yalnız `Hkeskin` zarfında
   karıştırdı. `HA4` zarfında bir vekil (VE1) kurulursa kesim ve kilit
   eksenlerinin **çarpıştığı** nokta ölçülür: `f_j`'nin kilit
   sökülmüş bir ailedeki değeri, ayrışımın parametresiz mi yoksa
   ailenin bir eseri mi olduğunu söyler.
5. **ARA KARIŞTIRMA.** `φ_q ~ U(0, 2πα)` ile `α ∈ [0,1]` bir
   **kilit merdiveni** kurar; gerçek gazın `α_eş`'i (λ_eş ve τ_c^eş'in
   yanına üçüncü koordinat) doğrudan ölçülür ve ΔM'nin %91.4'lük kilit
   payına bir SAYI verir. Beş inşa (~50 dk) yeterli.

---

## FİGÜR

`176_kilit_faturasi.png` (`176g_figur.py`) — dört panel:
**(a)** `R_η`: λ ailesi, kesim ailesi, ikiz/gerçek ve **vekiller**;
ön-kayıtlı bant, merkez ve ölüm eşiği; `Λ_R = +0.4565` oku.
**(b)** çarpan defteri: `Δlog(gerçek←ikiz)`, `Δlog(erfc←ikiz)` (kesim
ekseni) ve `Δlog(vekil←ikiz)` (kilidin toplam içeriği) — `g_E` ve `g_X`
kilit sökülünce **yukarı**, `θ` **aşağı**, `M` **yukarı**.
**(c)** DC kaçağının `⟨cosΔφ⟩` bant defteri: T-3'ün faz yarısının
ayakta kalışı ve en üst bantta gerçeğin ikiz ile vekil arasındaki yeri.
**(d)** H-F2 kaba kalemi: τ>0.70 açığının genlik/κ/faz payları
(%61.5 / %13.5 / %25.0) ve hedef bandı — sahte isabetin resmi.
