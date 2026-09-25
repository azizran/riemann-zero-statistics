# 118b ONARIM RAPORU (25 Eylül 2026) — χ₃ ve χ₈ₒ dosyalarında eksik sıfır çiftleri

**Kaynak:** `qm_riemann/195_uydu_karakteri_RAPOR.md` (K0a tablosu, satır ~75-115)
mpmath χ₃ adasında t≈246.30/246.42 ve χ₈ₒ adasında t≈242.55/242.69 civarında
dosyada OLMAYAN birer yakın sıfır çifti bulmuştu (L≈4.8; 195'in analiz
bölgesi L≥8.5'in çok altında). Bu rapor o iki çiftin yüksek hassasiyetli
konumlarını verir, dosyaların TÜM t aralığında başka eksik sıfır olup
olmadığını dört bağımsız yöntemle denetler ve onarılmış kopyaları üretir.

Git'e DOKUNULMADI. Orijinal `.npz` dosyaları DEĞİŞTİRİLMEDİ; onarılmış
kopyalar `_onarim` ekiyle yeni dosya adlarıyla yazıldı.

---

## 0) Hangi dosya hangi karakter?

`118b_yedi_ada_taperli_dokum.py`'deki `ADALAR` listesi ve `118a_taperli_L_motoru.py`
karakter tabloları (195 raporunun K0b tablosuyla çapraz doğrulandı):

| ada etiketi | dosya | q (=k) | a (parite: 0 çift/1 tek) | χ tablosu (0..q−1) |
|---|---|---|---|---|
| χ₃ | `118b_chi3_zeros.npz` | 3 | 1 (tek) | {0:0, 1:1, 2:−1} |
| χ₈ₒ | `118b_chi8o_zeros.npz` | 8 | 1 (tek) | {0:0,1:1,2:0,3:1,4:0,5:−1,6:0,7:−1} |

İkisi de gerçel, ilkel, tek (χ(−1)=−1) karakterler — 195'in K0b tablosuyla
(satır χ₃, χ₈ₒ) bit-bit aynı. `q`/`a` alanları npz içinde de doğrulandı
(`d["q"]`, `d["a"]`).

## 1) Dosya içerikleri (onarımdan ÖNCE)

| alan | `118b_chi3_zeros.npz` | `118b_chi8o_zeros.npz` |
|---|---|---|
| anahtarlar | zeros, bolgeler, q, a, T1, taper_c, nham, nB, nC, nD, neval, surlo, surhi, sure | (aynı) |
| q, a | 3, 1 | 8, 1 |
| T1 (hedef üst sınır) | 55000.0 | 42000.0 |
| sıfır sayısı (`zeros`) | 80205 | 65943 |
| t aralığı | [201.588914, 54999.622653] | [201.010137, 41999.583323] |
| sıralı mı / tekrar var mı | evet / hayır (min boşluk 4.45e-3) | evet / hayır (min boşluk 1.23e-2) |
| `bolgeler` (üretim sırasında bayraklanan şüpheli bölge) | BOŞ (0×2) | BOŞ (0×2) |

`bolgeler`in boş olması önemli: orijinal 118b boru hattının kendi K1/K2 tipi
sayım-sertifikası bu iki boşluğu HİÇ bayraklamamıştı (aşağıda §3'te bunun
nedeni tartışılıyor).

## 2) Bulunan eksik sıfırlar

**Yöntem:** L(1/2+it,χ) için gerçel Hardy tipi Z_χ(t) = e^{iθ(t)}·L(1/2+it,χ),
θ(t) = (t/2)log(q/π) + Im logΓ((s+a)/2) − arg(ε)/2 (98/118a motorunun AYNI
tanımı). L(1/2+it,χ) doğrudan Hurwitz-ζ toplamıyla (q^{-s}Σχ(r)ζ(s,r/q),
motorun yaklaşık-fonksiyonel-denklem toplamından TAMAMEN BAĞIMSIZ) mpmath
`mp.zeta`/`mp.loggamma` ile hesaplandı (dps=40). Her iki karakter de gerçel
ve ilkel olduğundan ε=1 (arg ε = 0) doğrudan hesaplanarak doğrulandı
(`|eps|=1.0` makine hassasiyetinde, tüm ondalıklarda 0'a yakın arg — bkz.
`hp.arg_eps_mp`), yani θ formülündeki −arg(ε)/2 terimi bu iki adada 0'dır.

Her iki bilinen boşlukta (chi3: idx 32, z∈[244.448, 248.395]; chi8o: idx 36,
z∈[241.500, 244.609]) adım 0.02 ile tarama önce İKİ işaret değişimi
bulundu, sonra her kök 100 iterasyonlu bisection ile (dps=40) rafine edildi:

| ada | t (yeni sıfır, ≥25 anlamlı hane) | Z(t) kalıntısı | Z'(t) (basit kök doğrulaması) | L=log(kt/2π) |
|---|---|---|---|---|
| χ₃ | 246.3028244345012325694137 | 3.0e-33 | +0.4966 | 4.7660 |
| χ₃ | 246.4148964765635420564502 | 1.3e-34 | −0.5087 | 4.7699 |
| χ₈ₒ | 242.5460897026670992725834 | −5.3e-33 | −0.7876 | 5.7302 |
| χ₈ₒ | 242.6908031431432587435340 | −5.4e-34 | +0.8392 | 5.7350 |

Z(t) kalıntıları dps=40'ın makine sıfırı mertebesinde (1e-33…1e-34) —
yani kökler istenen ≥1e-10 hassasiyetinin ÇOK ötesinde, pratikte dps=40
hassasiyetinde kesin. Z'(t)≠0 her ikisinde de basit (transversal) kök
olduğunu doğruluyor — teğet/çift kök değil.

Bu değerler 195 raporundaki yaklaşık değerlerle (246.303/246.415 ve
242.545/242.691) tam uyumlu.

## 3) Bütünlük denetimi — dosyaların TÜM t aralığında başka eksik sıfır var mı?

Dört bağımsız yöntem kullanıldı; hepsi TÜM dosya aralığında (yalnız 195'in
L≥8.5 penceresinde DEĞİL) çalıştırıldı:

### (A) N̄(t)−index tamsayı-sıçrama denetimi (Riemann–von Mangoldt tipi)

N̄_χ(t) = (t·log(kt/2π) − t)/(2π) + C_χ, C_χ = +1/8 (her ikisi de TEK
karakter). j0 (z₁'den önceki sıfır sayısı) üç yoldan hesaplandı:
- (i) ortalama N̄(mid)−i yuvarlaması,
- (ii) tek-nokta N̄(z₁)−½,
- (iii) BAĞIMSIZ: mpmath Hardy Z işaret değişimi sayımı, (0, z₁) aralığı,
  adım 0.02, dps=15 — 195'in yönteminden bağımsız yeniden üretildi.

| ada | j0 (i) | j0 (ii) | j0 (iii, bağımsız) | 195 raporundaki j0(iii) |
|---|---|---|---|---|
| χ₃ | 116 | 114 | **114** | 114 ✓ |
| χ₈ₒ | 147 | 145 | **145** | 145 ✓ |

d(i) = N̄(mid_i) − i − j0(iii) dizisinin 41-noktalık kayan medyanı, yuvarlanmış
tamsayı seviyesindeki HER sıçrama (dosyanın TÜM sıfır aralığında, 80204/65942
ara-nokta) arandı:

| ada | toplam sıçrama | benzersiz şüpheli boşluk | konum | seviye atlaması |
|---|---|---|---|---|
| χ₃ | 2 | **1** | z=[244.448267, 248.395359], boşluk=3.9471, L=4.7597 | +1 |
| χ₈ₒ | 2 | **1** | z=[241.500336, 244.609375], boşluk=3.1090, L=5.7284 | +1 |

Yani dosyanın TAMAMINDA (t≈201'den 55000/42000'e kadar) bu iki tanesi
DIŞINDA hiçbir integer-seviye sıçraması yok.

### (B) Tam-aralık "kaldırılmış dip" adayı taraması (118b'nin C adımının aynısı, tüm dosyaya uygulandı)

118b'nin dip-kurtarma algoritması (işaretsiz |Z| yerel minimumu, eşik 0.8,
komşulukta [-60,+60] ızgara noktası işaret DEĞİŞMİYOR şartı, mevcut sıfıra
0.05⟨g⟩'den yakınsa elenir) orijinal boru hattında yalnız bayraklanan
`bolgeler` içinde çalıştırılmıştı — bu iki dosyada `bolgeler` BOŞ olduğu için
bu adım hiç çalışmamıştı. Burada adım 0.005⟨g⟩ ızgarayla TÜM t aralığına
uygulandı:

| ada | ızgara noktası | Z hesap süresi | yerel minimum | eşik-altı + izole aday |
|---|---|---|---|---|
| χ₃ | 16 041 614 | 24.8 s | 80204 | **1** (t≈246.3575, \|Z\|_min≈1.60e-2) |
| χ₈ₒ | 13 189 231 | 28.1 s | 65944 | **1** (t≈242.6197, \|Z\|_min≈1.57e-2) |

Her iki adada da TEK aday bulundu ve bu aday tam olarak §2'de bulunan
çiftin ORTA NOKTASI — yani bu yöntem de dosyanın geri kalanında başka
bir "kaldırılmış dip" olmadığını gösteriyor.

### (C) Bağımsız ince tam-yeniden-tarama + eşleştirme

Taper motoru (`LmotorT`, 118a) `grid_frac=0.001` ile (üretimin ilk 0.03
ızgarasından 30× ince) t∈[200, T1] TAM aralığında SIFIRDAN yeniden
çalıştırıldı (chi3: 80203 sıfır; chi8o: 65941 sıfır) ve dosyayla injektif
en-yakın-komşu eşlemesi yapıldı (tolerans 0.25⟨g⟩_yerel):

- İnce taramada bulunup dosyada KARŞILIĞI OLMAYAN (= potansiyel yeni eksik
  sıfır): **0** (her iki adada da).
- Dosyada olup ince taramanın YAKALAYAMADIĞI: χ₃'te 2 (t≈1634.9113,
  1634.9662 — ayrım ≈0.055), χ₈ₒ'da 2 (t≈33452.4637, 33452.4789 — ayrım
  ≈0.015). Bunlar EKSİK DEĞİL: doğrudan mpmath nokta kontrolüyle
  (dps=20) o t değerlerinde gerçekten işaret değiştiren GERÇEK sıfırlar
  olduğu doğrulandı (Z değerleri sırasıyla +9.7e-3→+2.3e-5→−9.7e-3→
  +1.3e-4→+3.9e-2 ve −7.6e-3→+4.1e-4→+1.2e-3→−6.9e-4→−1.9e-2 — iki
  ayrı işaret değişimi net). Bunlar zaten 118b'nin B/C/D adımlarında
  (gerçek mpmath dip-kurtarmasıyla) doğru şekilde dosyaya girmiş dar
  çiftler; basit-taramanın bunları KAÇIRMASI aynı "kaldırılmış dip"
  mekanizmasının (taper Z'sinin bu dar çiftlerde işaret değiştirmeme
  eğilimi) ters yönde bir tezahürüdür, dosya hatası DEĞİLDİR.

### (D) t ≤ 1000 için motor-BAĞIMSIZ yoğun mpmath taraması

t∈[201, 1000] aralığında, adım 0.03⟨g⟩(t) (yerel ortalama boşluğun kesri,
dar çift ayrımının ~0.08-0.11⟨g⟩'sinden belirgin ince), dps=12, DOĞRUDAN
Hurwitz-ζ tabanlı Z_χ (motordan bağımsız) ile işaret değişimi sayıldı ve
dosyanın aynı [201,1000] aralığındaki sıfır sayısıyla karşılaştırıldı.

| ada | mpmath işaret değişimi | dosya sıfır sayısı [201,1000] | FARK | eşleşmeyen adaylar |
|---|---|---|---|---|
| χ₃ | 709 | 707 | **+2** | t≈246.297, t≈246.376 (§2'deki çiftin kaba parantezleri) |
| χ₈ₒ | 834 | 832 | **+2** | t≈242.530, t≈242.662 (§2'deki çiftin kaba parantezleri) |

Her iki adada da FARK tam olarak +2 ve eşleşmeyen adaylar §2'de bulunan
çiftlerin konumunda (kaba ızgara/dps=12 nedeniyle rafine değerlerden
~0.006-0.03 kayık, beklenen) — yani bu bağımsız, motor KULLANMAYAN dördüncü
yöntem de tam olarak AYNI iki sıfırı buluyor ve [201,1000] aralığında
BAŞKA fazlalık/eksik YOK (χ₈ₒ koşusu 1171 s sürdü, `dense_low_t.log`).

**SONUÇ:** Dört yöntem de birbirinden bağımsız biçimde AYNI sonuca varıyor:
her iki dosyada da §2'deki İKİ ÇİFT (toplam 4 sıfır) dışında, dosyanın TÜM
t aralığında (yalnız 195'in L≥8.5 üst bölgesinde değil) başka hiçbir eksik
sıfır YOK.

## 4) 195 analiz bölgesinin (L ≥ 8.5) etkilenmediğinin gerekçesi

Her iki eksik çift de L≈4.8-5.7 civarında (χ₃: L=4.766-4.770, χ₈ₒ:
L=5.730-5.735) — 195'in üst-bölge sınırı L_üst=8.5'in ÇOK altında. 195'in
K0a kapısı zaten bunu "kusur üst bölgede mi: HAYIR" olarak kaydetmişti;
burada bağımsız olarak da L hesabı doğrulandı. A ve B ölçümleri (195'in
çekirdek K̃ ve güç Ĝ istatistikleri) indeks n'yi HİÇ kullanmıyor — yalnız
mutlak sıfır konumları m_n, g_n kullanılıyor — bu yüzden n-kayması onları
zaten etkilemiyordu; C_χ kalibrasyon kapısı (K0a) ise kusur düzeltilince
|r|≤1e-5 veriyordu. Bu onarım, 195'in GEÇTİ hükmünü DEĞİŞTİRMİYOR; yalnızca
kaynak `.npz` dosyalarının kendisini (gelecekteki tüm kullanımlar için)
düzeltiyor.

## 5) Onarılmış dosyalar

Aynı anahtar kümesi korunarak (`zeros, bolgeler, q, a, T1, taper_c, nham,
nB, nC, nD, neval, surlo, surhi, sure`), yalnız `zeros` (2 yeni sıfır
eklendi, yeniden sıralandı, tekrar kaydı YOK — min boşluk kontrolü ile
doğrulandı) ve `nD` (tanımı gereği = len(zeros), yeni toplamı yansıtacak
şekilde güncellendi) değiştirildi. `bolgeler, q, a, T1, taper_c, nham, nB,
nC, neval, surlo, surhi, sure` alanları orijinal üretim koşusunun
istatistikleri olarak AYNEN korundu.

| dosya | eski n | yeni n | eklenen sıfırlar |
|---|---|---|---|
| `118b_chi3_zeros_onarim.npz` | 80205 | **80207** | 246.3028244345012, 246.4148964765635 |
| `118b_chi8o_zeros_onarim.npz` | 65943 | **65945** | 242.5460897026671, 242.6908031431433 |

### SHA-256

| dosya | sha256 |
|---|---|
| `118b_chi3_zeros.npz` (orijinal, DEĞİŞMEDİ) | `b11f797d4cdaf452bbf4557a99f19963808826d463695938827e8aa1e7f1e548` |
| `118b_chi8o_zeros.npz` (orijinal, DEĞİŞMEDİ) | `7ac99c3dde76d2bc310c659d36a044d9bd64c9fed88461245511500f4c72ebe6` |
| `118b_chi3_zeros_onarim.npz` (yeni) | `e8ec6da890fcbfc95f98262a09576ca0508ccbe881b72ef30271631c3ba8f79b` |
| `118b_chi8o_zeros_onarim.npz` (yeni) | `1259e0896ee8aa8be58a5db4b0ff6cb5ed7be55ba2bece28726700cec1dab4ba` |

## 6) Araçlar

Python: `/Users/ugursezen/Desktop/arin/deney/.venv/bin/python` (mpmath).
Kullanılan betikler (proje diziniyle karışmasın diye scratchpad'te tutuldu,
kalıcı değildir): `lib_setup.py` (118a'dan Lmotor/LmotorT/χ tablolarının
yan-etkisiz içe aktarımı), `hp.py` (bağımsız yüksek-hassasiyetli Hardy Z),
`nbar_check.py` (§3A), `dip_scan_full.py` (§3B), `match_fine.py` +
ince-tarama önbellekleri (§3C), `dense_low_t.py` (§3D).

Git'e DOKUNULMADI. Bu rapor ve onarılmış dosyalar dışında hiçbir mevcut
dosya değiştirilmedi.
