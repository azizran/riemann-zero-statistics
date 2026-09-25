# NOT 7 — Hakem Teftişi (25 Eylül 2026)

**Denetlenen taslak:** `qm_riemann/arxiv_comb_satellites.tex` (tamamı okundu; figürler ve kaynakça yer tutucu olduğu için denetlenmedi).
**Yöntem:** (a) Derivation 1 (durağan faz) ve Derivation 2 (BK kapalı-form) bağımsız olarak elle ve `sympy`/FFT ile yeniden türetildi/sayısal olarak sınandı (kod: bkz. Ek). (b) Taslaktaki her sayı, 7 paralel araştırmayla ilgili kaynak raporlarına (188, 190, 192, 193, 194, 195, 197, 198, KALEM_*, LITERATUR_TARAMASI, 196_bk_uydu_katsayilari.py) karşı satır satır karşılaştırıldı. (c) Aşırı iddia/eksiklik ve okunurluk taslak metninin kendi iç tutarlılığı açısından ayrıca denetlendi.

**Toplam bulgu: 26** (KRİTİK 6, ÖNEMLİ 13, KÜÇÜK 7), artı 3 bağımsız doğrulanmış türetim (hatasız).

---

## 0. Bağımsız türetim sonuçları (pozitif bulgular)

Bunlar hata İÇERMİYOR — teftişin dürüstlüğü için önce belirtiliyor:

- **Denklem (2)** (`δ_n = g_nL(m_n)/2π − 1 = Σ_q 2a_q sin(ω_qg_n/2)cos(ω_qm_n)`, satır 133-136): Denklem (1)'den elle yeniden türetildi (`sin A − sin B` özdeşliğiyle), işaret dahil **birebir doğrulandı**.
- **Derivation 1** (durağan faz, satır 227-243): `2πN̄(t)=tL(t)−t+7π/4` standart Riemann–von Mangoldt asimptotiğiyle (`N̄(T)=(T/2π)log(T/2πe)+7/8`) `sympy` ile **cebirsel olarak özdeş**; Jacobi–Anger işareti sayısal olarak `1e-16` hassasiyetle doğrulandı; `φ(t*) = t* − 7π/4` ve toplam faz sabitinin (`−7π/4 − π/4`) `mod 2π` tam olarak `0`'a sadeleştiği **sembolik olarak kanıtlandı**. Türetim matematiksel olarak sağlam.
- **Derivation 2** (BK kapalı-form, satır 359-372): `A_k`, `B_0,B_1,B_2` ve `ĝ_p(k)/ĝ_p(0)` formülleri hem sembolik (`sympy`, tam rasyonel aritmetik) hem sayısal (FFT, p=2,3,5,7,11) olarak **makine hassasiyetinde doğrulandı**. Sağlam bir türetim.

Bu üç doğrulama, taslağın matematiksel omurgasının sağlam olduğunu gösteriyor; aşağıdaki bulgular çoğunlukla **veri/kaynak eşleşmesi, çerçeveleme ve okunurlukla** ilgili.

---

## 1. KRİTİK bulgular

**K1 — satır 192.** "the kernel is 77% rank-one" nötr bir gözlem gibi sunuluyor. Kaynak (188_zeta_kimligi_RAPOR.md, K3(i)) bunun **ön-kayıtlı bir mekanizma testi (H-188a) olduğunu ve resmen ÖLDÜĞÜNÜ** gösteriyor: rank-1 payı 0.773±0.024, eşik ≥0.80 tutmadı; ayrıca corr(|v|,S)=0.127, gerekli eşik 0.50'nin çok altında. **Düzeltme:** "77%" rakamını korurken, bu sayının besleediği ön-kayıtlı eşiğin tutmadığını bir cümleyle açıkça belirtin (ör. "a pre-registered rank-one threshold of 0.80 was not met").

**K2 — satır 200-210 (Observation, Height independence).** Δω-ekseni/τ'-ekseni korelasyonu (0.980 vs 0.002) aynı "pre-registered" kutusuna konmuş, ama kaynak (190_tarak_evrensellik_RAPOR.md) bu korelasyon karşılaştırmasının **ön-kayıt kapsamı dışında, veri görüldükten SONRA yapılmış bir keşif** olduğunu açıkça yazıyor (ön-kayıt sha256=38d122a2…, 17:42:39 — korelasyon hesaplaması bundan sonra). **Düzeltme:** Bu iki sayıyı kutudan ayırın veya "(post-hoc)" ibaresiyle işaretleyin.

**K3 — satır 279-280.** "+log 7 ... magnitudes 1.5–1.9 times the prediction." Kaynağın (193_sinif_yasasi_RAPOR.md) kendi ham K1 tablosundan hesaplanan gerçek oranlar **1.84–2.44×**'tür, 1.5–1.9× değil. "1.5-1.9" rakamı kaynağın kendi HÜKÜM/MANŞET özetinde geçiyor ama bu özet kaynağın kendi ham verisiyle çelişiyor; taslak bu iç tutarsızlığı süzmeden devralmış. **Düzeltme:** Ham K1 tablosundan aralığı yeniden hesaplayıp düzeltin veya tutarsızlığı bir dipnotla açıklayın.

**K4 — satır 284-287.** "...they carry small negative residues, 9.4σ below the empty-window floor, which persist at all heights (Section~\ref{sec:height})." 9.4σ rakamı doğru (194_bos_pencere_RAPOR.md, H-194d) ama **"tüm yüksekliklerde kalıcı" iddiası desteksiz**: kaynak (194) yalnızca tek pencereyi (L=12.03) inceliyor, çok-yükseklik testi yok; taslağın kendi §7'si de bu μ(a)=0 rezidülerini hiç tartışmıyor (yalnız κ_p'den bahsediyor). Çapraz referans hem kaynakta hem hedef bölümde karşılıksız. **Düzeltme:** "which persist at all heights" ifadesini kaldırın veya gerçek çok-yükseklik kanıtını ekleyin/kaynak gösterin.

**K5 — satır 381-383.** "the factor of 2 explains why $+\log 2m$ is twice as bright as $+\log m$" — **niteleyici eksik ve iddia yanlış**. m çift olduğunda (ör. m=2 ⇒ +log4 = +log(2²)) bu satırdan üç önceki paragrafın kendi kuralına göre ($k\ge2$ için $f_p(k)=0$) uydu **tam karanlıktır**, "iki katı parlak" değil. Aynı paragrafın negatif-taraf ifadesi ("c(2/m)=2/m **for odd m**") bu niteleyiciyi doğru şekilde içeriyor; pozitif-taraf cümlesi aynı niteleyiciyi unutmuş — bu iç tutarsızlık, Derivation 2'nin "Three consequences" paragrafında somut bir hata. **Düzeltme:** "...twice as bright as $+\log m$ (for $m$ coprime to $2$)" şeklinde düzeltin.

**K6 — satır 397.** "The mirror law was tested in the band $\Delta\omega\in[-2.12,-1.58]$ ... which had never been computed." Kaynak (KALEM_AYNA_YASASI_25EYL2026.md) bandın yalnızca **$[-2.12,-1.775)$ kısmının hiç hesaplanmadığını**, $[-1.775,-1.58]$ kısmının ise daha önceki bir haritada (190b, $[-1.775,+3.55]$) **zaten hesaplanmış ama aranmamış/gösterilmemiş** olduğunu özenle ayırıyor. Taslak bu ayrımı sildi — bu, makalenin amiral gemisi kör testinin körlük iddiasını abartıyor. **Düzeltme:** "the sub-band $[-2.12,-1.775)$ had never been computed; $[-1.775,-1.58]$ had been computed earlier but not searched or analysed for this comb" şeklinde hassaslaştırın.

## 2. ÖNEMLİ bulgular

**Ö1 — satır 261-262.** "$-8.6°\pm3.9°$ apart ... exactly like the control $-\log 3$" — kontrolün gerçek değeri $-0.9°\pm5.5°$, $0.493/0.507$'dir; sayısal olarak özdeş değil (istatistiksel olarak ayırt edilemez düzeyde, ama "exactly like" abartılı). **Düzeltme:** "statistically indistinguishable from" ile değiştirin.

**Ö2 — satır 263-265.** "none of nine positions with $\mu(a)=0$ lights up" — kaynak (192_uydu_teorisi_RAPOR.md) 9 pozisyondan yalnız 3'ünün "temiz sıfır" olduğunu, 4'ünün ise yalnızca geniş bir taban tanımı sayesinde söndürülmüş sayıldığını ve bunların anlamlı negatif κ (−2.2σ ile −7.9σ) taşıdığını, bunu bizzat kendisinin "tasarım sınırı" diye işaretlediğini belirtiyor. **Düzeltme:** Bu nüansı bir yan cümleyle açıklayın.

**Ö3 — satır 280-281.** "$+\log 5$ ... shows the same stretch ($\approx1.3$)" — kaynak aralığı 1.3–1.8'dir (ham veriden yeniden hesap: 1.25–1.76), tek başına "≈1.3" alt ucu seçiyor; ayrıca bu aralık +log7'nin (K3'te düzeltilen) gerçek aralığıyla (1.84–2.44) "aynı" değil. **Düzeltme:** Tam aralığı verin, "the same stretch" ifadesini kaldırın.

**Ö4 — satır 318.** "allowed ones $+9.8$ to $+12.1\sigma$" — kaynak (195_uydu_karakteri_RAPOR.md) 24 izinli uydudan **22'sinin** bu aralıkta olduğunu, 2'sinin (χ₃/χ₅ₑ'de +log4, μ(4)=0) istisna olduğunu gösteriyor. **Düzeltme:** "22 of 24" olarak belirtin veya dipnotla istisnaları açıklayın.

**Ö5 — satır 320.** "$z=-55$ to $-82$" yalnızca birincil/seçilmiş uydu çiftleri altkümesini temsil ediyor; kaynağın tam B-ii tablosu daha geniş bir aralığa (ör. −25.0, −30.9) yayılıyor. **Düzeltme:** "for the primary satellite pairs" gibi bir niteleyici ekleyin.

**Ö6 — satır 348-350 (eq:BK).** $\Phi(\varepsilon)=\prod_p[1-(1-p^{i\varepsilon})^2/(p-1)^2]$ açık formülü KALEM_AYNA_YASASI kaynağında yazılı değil (kaynak yalnız $f=|\zeta|^2\Phi_{\text{off}}=\prod g_p$ veriyor); matematiksel olarak tutarlı ve türetilebilir ama kaynakta doğrudan yazılı olmayan bir formül sanki doğrudan alıntılanmış gibi sunuluyor. **Düzeltme:** Bir cümlelik ara-adım veya doğrudan kaynak/literatür ataması ekleyin.

**Ö7 — satır 385-393 ("A lesson in line shape").** 23-uydu kataloğu karşılaştırması veri-sonrasıdır (kör değildir) ama bu alt-bölüm, §7'deki analog durumun aksine ("A post-hoc look suggested..."), hiçbir yerde "post-hoc"/"not blind" ifadesi kullanmıyor; dürüstlük yalnız hemen ardından gelen "The blind test" alt-başlığıyla zımni kontrastla iletiliyor. **Düzeltme:** Alt-bölüm başına açık bir "(post-hoc, not blind)" notu ekleyin.

**Ö8 — satır 461-464 (Observation, Height dependence).** $\Lambda_{\text{last}}/\Lambda_{\text{alt}}$ oranının sabiti (γ=0) 7.9σ'da dışladığı raporlanıyor, ama kaynak aynı oranın γ=1 hipotezinden de 1.9σ saptığını gösteriyor — bu gerilim atlanmış. **Düzeltme:** Bir yarım cümleyle bu kalıntı gerilimi belirtin.

**Ö9 — satır 464-466.** "The half-integer ... families relax towards one in the same way ($0.79\to0.83\to0.92$...)" yalnız TEK yarı-tamsayı ailesini (R(5/2)) raporluyor; kaynakta (198_kappa_yukseklik_RAPOR.md) İKİNCİ bir aile (R(7/2): $0.75\to0.81\to0.88$) de var ve taslakta hiç anılmıyor. **Düzeltme:** Her iki aileyi de verin veya seçici raporlamanın gerekçesini açıklayın.

**Ö10 — Özet satır 36-40 vs. gövde satır 278-281.** Özet "the refined law ... passes a blind test" derken, gövdedeki +log7/+log5'in açıklanamayan büyüklük fazlalığını (K3, Ö3) hiç yansıtmıyor; yalnız özeti okuyan biri yasanın genel olarak doğrulandığını sanabilir. **Düzeltme:** Özete "...though it overshoots two other satellites by a still-unexplained 30–140%" türü bir niteleyici ekleyin.

**Ö11 — satır 166-169 (eq:K).** $C_j(q)$ hiçbir yerde açıkça tanımlanmıyor; yalnız düzyazıda "the slice's per-gap series is projected onto the pool lines within the block" deniyor. Dış okur bu niceliği önceki $c_q$ tanımına benzer bir formülle yeniden kuramaz. **Düzeltme:** $C_j(q) := 2\langle \delta_{n,j}^{\text{block}} e^{-i\omega_q m_n}\rangle$ türü açık bir formül ekleyin; "havuz"un $\mathcal P$ ile aynı küme olduğunu da belirtin.

**Ö12 — satır 144-146.** "mix, produced by the other window lines through the zero-dependent sampling grid" — mix'in neden sıfır olmadığının açıklaması tek cümleye sıkıştırılmış; dış okur için yetersiz. **Düzeltme:** Düzensiz örnekleme ızgarasında farklı frekansların neden tam ortogonal olmadığını (ve bunun T→∞'da kalıcı mı yoksa sonlu-pencere etkisi mi olduğunu) 1-2 cümleyle açın.

**Ö13 — satır 240 (Derivation 1) / provenance.** Sınıf-yasası kör testinin başarısı için kritik olan $-\pi/4$ durağan-faz terimi, iki en ilgili kaynakta (192_uydu_teorisi_RAPOR.md, KALEM_UYDU_TEORISI_23EYL2026.md) **yer almıyor** — bu iki kaynağın kendi türetimi yalnız $\phi^*\equiv 2\pi q'b/a - 7\pi/4 \pmod{2\pi}$'dir, $-\pi/4$ terimi hiç geçmiyor. Ben bunun standart durağan-faz teoreminin doğru bir parçası olduğunu bağımsız olarak doğruladım (bkz. §0), yani matematiksel olarak SORUN YOK; ama denetim listesindeki iki kaynakta bu düzeltmenin ilk nereden geldiği izlenemiyor. **Düzeltme:** $-\pi/4$ teriminin ilk türetildiği kalemi (192'den sonraki, muhtemelen 193 civarı) açıkça kaynak gösterin.

## 3. KÜÇÜK bulgular

**KU1 — satır 171-172.** "verified to $10^{-15}$" — kaynak (190_tarak_evrensellik_RAPOR.md) aslında $2.2\times10^{-16}$ gösteriyor; yanıltıcı değil ama gereksiz gevşek. Düzeltme: rakamı sıkılaştırın veya "to at least $10^{-15}$" yazın.

**KU2 — satır 122-125 ("alt" penceresi).** Tanımlayıcı parametreleri (t-aralığı, gap sayısı, L=9.34) verilen kaynak listesinde (188/190) izlenemiyor; muhtemelen listede olmayan erken bir kalemden geliyor. Düzeltme: doğrudan kaynağı dipnotla belirtin.

**KU3 — κ notasyon çakışması.** $\kappa(\omega')/\kappa_j$ (çekirdek, §2-4, işaretli gerçel sayı) ile $\kappa_p$ (§7, tahmine oranla ~1 civarı aktarım çarpanı) aynı simgeyi paylaşıyor; kavramsal olarak ilişkili ama yapısal olarak farklı nicelikler. Düzeltme: §7'deki niceliği $\eta_p$ veya $\lambda_p$ olarak yeniden adlandırın.

**KU4 — Özet uzunluğu.** ~373 kelime (doğrudan sayıldı) — math.NT için alışılmadık uzun/yoğun, beş bulguyu çok sayıda rakamla paketliyor. Düzeltme: özeti kısaltıp rakam yoğunluğunu gövdeye taşıyın.

**KU5 — "pre-registered" vs "blind".** Terimler açık bir operasyonel tanım olmadan, sanki eşanlamlıymış gibi kullanılıyor. Düzeltme: §2'de bir cümleyle ikisi arasındaki farkı (ör. "blind" = analiz kodunun veriden tamamen habersiz donması; "pre-registered" = eşiklerin veriden önce yazılı olarak dondurulması) tanımlayın.

**KU6 — Tablo 1 (satır 445-458).** Yalnız "alt" satırı "(blind)" diye işaretli; "last" penceresinin kısmi-görülmüşlük durumu yalnız alttaki gövde metninde açıklanıyor, tabloda değil. Düzeltme: tablo dipnotuna kısa bir not ekleyin.

**KU7 — $\tau(\chi)$ (satır 304, Gauss toplamı) ile $\tau_q/\tau'$ (satır 141, 194, frekans oranı) aynı harfi paylaşıyor.** Alanda $\tau(\chi)$ çok standart olduğu için düşük önem, ama bir dipnotla ayrım belirtilebilir.

---

## 4. Literatür konumlandırması (madde 5) — genel sonuç: ADİL

Landau–Gonek/BK/Conrey–Snaith/Rodgers/Zhang–Martelli–Torquato zinciri (satır 96-108) ölçülü ve doğru karakterize edilmiş; "we do not claim the mechanism is new" ifadesi literatür taramasının önerdiği çerçeveyle birebir örtüşüyor. Kanivets (Kan26) dipnotu (satır 513-518) adil ve doğru hedge'li ("may account for"); literatür taramasının kendi çekince listesinden en güçlü tek argümanı seçiyor — haksızlık değil, ama Ö7'nin de işaret ettiği gibi §6'daki "A lesson in line shape" ile aynı dürüstlük standardını tutarlı şekilde uygulaması iyi olur. companion1-4 atıflarında çelişki bulunamadı.

---

## Ek — Bağımsız doğrulama betikleri

`/private/tmp/.../scratchpad/verify_derivations.py` ve `verify_bk_closed_form.py` (bu oturumda çalıştırıldı, `/Users/ugursezen/Desktop/arin/deney/.venv/bin/python` ile): (1) $2\pi\bar N(t)=tL(t)-t+7\pi/4$ sympy ile standart asimptotikle özdeşliği, (2) Jacobi–Anger işareti sayısal spot-check, (3) durağan-faz sabitinin $-7\pi/4-\pi/4\equiv0\pmod{2\pi}$ sembolik ispatı, (4) $\hat g_p(k)/\hat g_p(0)$ kapalı-formunun $p\in\{2,3,5,7,11\}$ için hem FFT hem sembolik doğrulaması. Tüm dört test **geçti**; kodlar tekrar çalıştırılabilir haldedir (repoya yazılmadı, scratchpad'te).
