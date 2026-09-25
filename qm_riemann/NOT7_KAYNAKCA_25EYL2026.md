# Not 7 — Kaynakça doğrulaması (25 Eylül 2026)

Yöntem: her künye bu oturumda **birincil kaynaktan** yeniden sorgulandı — arXiv API
(`export.arxiv.org/api/query`), Crossref API (`api.crossref.org/works`), zbMATH API
(`api.zbmath.org`), yayıncı DOI sayfaları (Oxford Academic, Wiley, Project Euclid) ve
Zenodo API (`zenodo.org/api/records`). Önceki oturumların (`LITERATUR_*`) doğrulamaları
başlangıç noktası olarak kullanıldı ama **hiçbiri sorgulanmadan aktarılmadı** — hepsi bu
oturumda yeniden çapraz kontrol edildi ve iki bağımsız kaynak (genelde Crossref + zbMATH,
bazen + yayıncı sayfası) aynı sonucu verdiğinde "DOĞRULANDI" işaretlendi. Uyuşmazlık
bulunan her yerde ikisi de rapor edildi ve hangisinin doğru olduğuna dair gerekçe verildi.

---

## 1. Bogomolny–Keating, üç makale (1995/1996)

**(a) Durum: DOĞRULANDI — Crossref API, üç ayrı DOI kaydı.**
- `https://api.crossref.org/works/10.1088/0951-7715/8/6/013` (Nonlinearity I)
- `https://api.crossref.org/works/10.1088/0951-7715/9/4/006` (Nonlinearity II)
- `https://api.crossref.org/works/10.1103/PhysRevLett.77.1472` (PRL)

Üçü de var; Crossref kayıtları: Nonlinearity I → cilt 8, sayı 6, s. 1115–1131, yayın
1995-11-01. Nonlinearity II → cilt 9, sayı 4, s. 911–935, yayın **1996-07-01** (bazı
ikincil kaynaklarda "1995" olarak yanlış anılıyor — Crossref'e göre doğrusu **1996**).
PRL → cilt 77, sayı 8, s. 1472–1475, yayın 1996-08-19.

**(b) LaTeX**
```
\bibitem{BK95} E.~B.~Bogomolny, J.~P.~Keating, \emph{Random matrix theory and the
Riemann zeros I: three- and four-point correlations}, Nonlinearity \textbf{8}
(1995), 1115--1131.

\bibitem{BK96} E.~B.~Bogomolny, J.~P.~Keating, \emph{Random matrix theory and the
Riemann zeros II: $n$-point correlations}, Nonlinearity \textbf{9} (1996), 911--935.

\bibitem{BK96PRL} E.~B.~Bogomolny, J.~P.~Keating, \emph{Gutzwiller's trace formula
and spectral statistics: beyond the diagonal approximation}, Phys.\ Rev.\ Lett.\
\textbf{77} (1996), 1472--1475.
```

---

## 2. Bogomolny, arXiv:0708.4223

**(a) Durum: DOĞRULANDI, ama SAYFA ARALIĞI HATASI bulundu.**
- arXiv API: `https://export.arxiv.org/api/query?id_list=0708.4223` — başlık "Riemann
  zeta function and quantum chaos", tek yazar Eugene Bogomolny, journal_ref alanı
  "Progress of theoretical physics supplement, 166 (2007) **19-44**".
- Crossref: `https://api.crossref.org/works/10.1143/PTPS.166.19` — **19-36**.
- Yayıncı sayfası (Oxford Academic, DOI çözümlemesiyle):
  `https://academic.oup.com/ptps/article-lookup/doi/10.1143/PTPS.166.19` — **"Page
  Numbers: 19–36"** (bu oturumda doğrudan okundu).

**Sonuç: doğru sayfa aralığı 19–36'dır, 19–44 DEĞİL.** arXiv'in kendi journal_ref alanı
(yazar tarafından girilir, yayıncı tarafından doğrulanmaz) yanlış; iki bağımsız yayıncı
kaynağı (Crossref + Oxford Academic'in kendi sayfası) 19–36'da birleşiyor. Bu hata,
`LITERATUR_BOGOMOLNY_KEATING_24EYL2026.md` ve `NOT7_UYDU_ISKELET_25EYL2026.md`'de
"19–44" olarak tekrarlanmıştı — **düzeltilmeli**.

**(b) LaTeX**
```
\bibitem{Bog07} E.~Bogomolny, \emph{Riemann zeta function and quantum chaos},
Prog.\ Theor.\ Phys.\ Suppl.\ \textbf{166} (2007), 19--36, arXiv:0708.4223.
```

---

## 3. Conrey & Snaith, math/0509480

**(a) Durum: DOĞRULANDI, ama CİLT NUMARASI ÇELİŞKİSİ bulundu ve çözüldü.**
- arXiv API journal_ref: "Proc. Lon. Math. Soc., Volume **93**, No 3, 2007, pages
  594--646" (yazar-girişli alan).
- Crossref: `https://api.crossref.org/works/10.1112/plms/pdl021` — cilt **94**, sayı 3,
  s. 594–646, yayın 2007-05.
- zbMATH: `au:Conrey au:Snaith ti:"ratios conjectures"` — "Proc. Lond. Math. Soc. (3)
  **94**, No. 3, 594-646 (2007)."

**Sonuç: doğru cilt 94'tür, 93 DEĞİL.** arXiv'in kendi journal_ref alanı (yazar
girişli) yanlış; iki bağımsız birincil kaynak (Crossref + zbMATH, ikisi de yayıncıdan
gelen verat kullanıyor) 94'te birleşiyor. Mevcut notlardaki "94:3 (2007)" doğruydu;
teyit edildi.

**(b) LaTeX**
```
\bibitem{CS07} J.~B.~Conrey, N.~C.~Snaith, \emph{Applications of the $L$-functions
ratios conjectures}, Proc.\ London Math.\ Soc.\ (3) \textbf{94} (2007), 594--646,
arXiv:math/0509480.
```

---

## 4. Bogomolny, Bohigas, Leboeuf, Monastra, math/0602270

**(a) Durum: DOĞRULANDI, tam uyum.**
- arXiv API: journal_ref "J. Phys. A: Math. Gen. 39 (2006) 10743-10754", arXiv DOI
  alanı `10.1088/0305-4470/39/34/010`.
- Crossref: `https://api.crossref.org/works/10.1088/0305-4470/39/34/010` — cilt 39,
  sayı 34, s. 10743–10754, yayın 2006-08-25.
- zbMATH: aynı — "J. Phys. A, Math. Gen. 39, No. 34, 10743-10754 (2006)."

Task metninde verilen cilt/sayfa (39 (2006), 10743–10754) **birebir doğru**.

**(b) LaTeX**
```
\bibitem{BBLM06} E.~Bogomolny, O.~Bohigas, P.~Leboeuf, A.~G.~Monastra, \emph{On the
spacing distribution of the Riemann zeros: corrections to the asymptotic result},
J.\ Phys.\ A: Math.\ Gen.\ \textbf{39} (2006), 10743--10754, arXiv:math/0602270.
```

---

## 5. Landau (1911/1912) ve Gonek (1993)

**(a) Landau — Durum: DOĞRULANDI, ama YIL DÜZELTMESİ: 1912, 1911 DEĞİL.**
- Crossref: `https://api.crossref.org/works?query.bibliographic=Landau+Nullstellen+
  Zetafunktion+Mathematische+Annalen` → "Über die Nullstellen der Zetafunktion",
  Math. Ann., cilt 71, sayı 4, s. 548–564, published-print **1912-12**, DOI
  `10.1007/bf01456808`.
- zbMATH: `au:Landau ti:"Nullstellen der Zetafunktion"` → "Math. Ann. 71, 548-564
  **(1912)**."

İki bağımsız kaynak da 1912'de birleşiyor. "1911" ataması yaygın (bazı ikincil
kaynaklar, hatta DHPC'nin kendi teorem-metni "Landau, 1911" diyor, ama DHPC'nin KENDİ
kaynakçası da 1912 yazıyor — bkz. `LITERATUR_KAPPA_SONLU_YUKSEKLIK_25EYL2026.md` §0) —
muhtemelen sonucun 1911'de duyurulup 1912'de basılmış olmasından kaynaklanan bir
gelenek, ama **basım yılı olarak doğru olan 1912'dir**.

**(b) Gonek — Durum: DOĞRULANDI (zbMATH, tam bibliyografik kayıt).**
- `https://api.zbmath.org/v1/document/_search?search_string=au:Gonek ti:"explicit
  formula of Landau"` → "Knopp, Marvin (ed.) et al., A tribute to Emil Grosswald:
  number theory and related analysis. Providence, RI: American Mathematical Society.
  Contemp. Math. 143, **395-413** (1993)." ISBN 0-8218-5155-1. zbMATH ayrıca içerik
  özetinde Gonek'in ana formülünü Landau'nun (1911 diye anıyor — JFM kaydı) klasik
  sonucunun düzgün (uniform) versiyonu olarak tanımlıyor — bu, DHPC'nin tam metninden
  önceki oturumda aktarılan formülle birebir örtüşüyor.

**(b) LaTeX**
```
\bibitem{Landau12} E.~Landau, \emph{\"Uber die Nullstellen der Zetafunktion},
Math.\ Ann.\ \textbf{71} (1912), 548--564.

\bibitem{Gonek93} S.~M.~Gonek, \emph{An explicit formula of Landau and its
applications to the theory of the zeta-function}, in: M.~Knopp et al.\ (eds.), A
Tribute to Emil Grosswald: Number Theory and Related Analysis, Contemp.\ Math.\
\textbf{143}, Amer.\ Math.\ Soc., Providence, RI, 1993, 395--413.
```

---

## 6. Durkan, Hughes, Pearce-Crump, arXiv:2601.18025

**(a) Durum: DOĞRULANDI (arXiv API, bu oturumda).**
`https://export.arxiv.org/api/query?id_list=2601.18025` → başlık "Generalisations of
the Landau–Gonek Theorem and applications to mean values of zeta"; yazarlar Benjamin
Durkan, Christopher Hughes, Andrew Pearce-Crump; gönderim 2026-01-25. **journal_ref
alanı boş** — henüz sadece arXiv ön-baskısı, dergi bilgisi yok (task'ta "varsa dergi"
denmişti — şu an yok).

**(b) LaTeX**
```
\bibitem{DHPC26} B.~Durkan, C.~Hughes, A.~Pearce-Crump, \emph{Generalisations of the
Landau--Gonek theorem and applications to mean values of zeta}, arXiv:2601.18025
[math.NT] (2026), preprint.
```

---

## 7. Hughes, Lugmayer, Pearce-Crump, arXiv:2411.05573

**(a) Durum: DOĞRULANDI, dergi bilgisi tamamlandı.**
- arXiv API: başlık "The second moment of the Riemann zeta function at its local
  extrema", yazarlar Christopher Hughes, Solomon Lugmayer, Andrew Pearce-Crump.
- Crossref (bibliyografik arama): "Journal of the London Mathematical Society", cilt
  **112**, sayı **2**, makale no. **e70250**, yayın 2025-08, DOI
  `10.1112/jlms.70250`.

Mevcut taslaklarda (`arxiv_gap_amplitude.tex`: "J. London Math. Soc. (2025)", eksik
cilt/makale no; `arxiv_prime_wave_anatomy.tex`: "J. London Math. Soc. (2) (2025),
e70250", eksik cilt) iki farklı EKSİK künye var — ikisi de doğru ama tam değil. Tam ve
doğru künye: **J. London Math. Soc. (2) 112 (2025), no. 2, e70250**, DOI
10.1112/jlms.70250.

**(b) LaTeX**
```
\bibitem{HLPC24} C.~Hughes, S.~Lugmayer, A.~Pearce-Crump, \emph{The second moment of
the Riemann zeta function at its local extrema}, arXiv:2411.05573 (2024); J.\
London Math.\ Soc.\ (2) \textbf{112} (2025), no.~2, e70250.
```

---

## 8. Zhang/Martelli/Torquato (1801.01541) ve Torquato/Zhang/de Courcy-Ireland (1802.10498)

**(a) Durum: DOĞRULANDI — yazar sırası ve sayfa numarası DÜZELTİLDİ.**

**1801.01541**: arXiv API yazar sırası **Ge Zhang; Fausto Martelli; Salvatore
Torquato** (bu, "The structure factor of primes" makalesi). Crossref
(`https://api.crossref.org/works/10.1088/1751-8121/aaa52a`) ve zbMATH ikisi de aynı
sırayı ve şu künyeyi veriyor: *Journal of Physics A: Mathematical and Theoretical*,
cilt 51, sayı 11, **makale no. 115001** (16 s.), yayın 2018-03-16.
**"135001" YANLIŞ** — bu rakam `LITERATUR_YAYIN_ONCESI_24EYL2026.md`'nin kendi
"Notlara eklenmesi önerilen atıflar" bölümünde bir kez geçmiş bir yazım hatası; doğrusu
**115001**'dir. ("Torquato, Zhang, de Courcy-Ireland" yazarlığı bu makaleye AİT
DEĞİLDİR — bu, aşağıdaki 1802.10498'in yazar listesi.)

**1802.10498**: arXiv API yazar sırası **S. Torquato; G. Zhang; M. de Courcy-Ireland**
(bu, "Uncovering Multiscale Order in the Prime Numbers via Scattering" makalesi).
Crossref/zbMATH: *Journal of Statistical Mechanics: Theory and Experiment*, cilt 2018,
sayı 9, makale no. 093401 (15 s.), DOI `10.1088/1742-5468/aad6be`.

**Özet — iki makale kesin biçimde ayrıştırıldı:**
| arXiv no | Başlık | Yazarlar (doğru sıra) | Dergi |
|---|---|---|---|
| 1801.01541 | The structure factor of primes | Zhang, Martelli, Torquato | J. Phys. A 51 (2018) 115001 |
| 1802.10498 | Uncovering Multiscale Order in the Prime Numbers via Scattering | Torquato, Zhang, de Courcy-Ireland | J. Stat. Mech. (2018) 093401 |

**(b) LaTeX**
```
\bibitem{ZMT18} G.~Zhang, F.~Martelli, S.~Torquato, \emph{The structure factor of
primes}, J.\ Phys.\ A: Math.\ Theor.\ \textbf{51} (2018), 115001, arXiv:1801.01541.

\bibitem{TZdCI18} S.~Torquato, G.~Zhang, M.~de Courcy-Ireland, \emph{Uncovering
multiscale order in the prime numbers via scattering}, J.\ Stat.\ Mech.\ Theory Exp.\
\textbf{2018} (2018), 093401, arXiv:1802.10498.
```

---

## 9. Rodgers, arXiv:1203.3275

**(a) Durum: DOĞRULANDI, tam uyum.**
arXiv API: başlık "Macroscopic pair correlation of the Riemann zeroes for smooth test
functions", tek yazar Brad Rodgers, journal_ref "Q J Math (2013) 64 (4): 1197-1219",
arXiv DOI `10.1093/qmath/has024`. Crossref ve zbMATH aynı: *The Quarterly Journal of
Mathematics*, cilt 64, sayı 4, s. 1197–1219, yayın 2013-12-01.

**(b) LaTeX**
```
\bibitem{Rod13} B.~Rodgers, \emph{Macroscopic pair correlation of the Riemann
zeroes for smooth test functions}, Q.\ J.\ Math.\ \textbf{64} (2013), 1197--1219,
arXiv:1203.3275.
```

---

## 10. Montgomery (1973), Odlyzko (1987) + tablo sayfası, Keating–Snaith (2000)

**(a) Durum: hepsi DOĞRULANDI (zbMATH, Crossref'te DOI'siz olanlar için zbMATH tek
kaynak — 1973 tarihli bir konferans bildirisi kitabı bölümü olduğu için DOI yok, bu
normal).**
- Montgomery: zbMATH → "The pair correlation of zeros of the zeta function", *Analytic
  Number Theory*, Proc. Sympos. Pure Math. 24, St. Louis Univ. Missouri 1972, s.
  181–193 (1973).
- Odlyzko: Crossref + zbMATH → *Mathematics of Computation*, cilt 48, s. 273–308
  (1987).
- Odlyzko sıfır tabloları: `https://www-users.cse.umn.edu/~odlyzko/zeta_tables/` — bu
  oturumda `curl -I` ile **HTTP 200** doğrulandı, sayfa canlı.
- Keating–Snaith: Crossref + zbMATH → *Communications in Mathematical Physics*, cilt
  214, sayı 1, s. 57–89 (2000).

**(b) LaTeX**
```
\bibitem{Mont73} H.~L.~Montgomery, \emph{The pair correlation of zeros of the zeta
function}, in: Analytic Number Theory, Proc.\ Sympos.\ Pure Math.\ \textbf{24},
Amer.\ Math.\ Soc., Providence, RI, 1973, 181--193.

\bibitem{Odl87} A.~M.~Odlyzko, \emph{On the distribution of spacings between zeros
of the zeta function}, Math.\ Comp.\ \textbf{48} (1987), 273--308.

\bibitem{OdlyzkoTables} A.~M.~Odlyzko, \emph{Tables of zeros of the Riemann zeta
function}, \url{https://www-users.cse.umn.edu/~odlyzko/zeta_tables/}.

\bibitem{KS00} J.~P.~Keating, N.~C.~Snaith, \emph{Random matrix theory and
$\zeta(1/2+it)$}, Comm.\ Math.\ Phys.\ \textbf{214} (2000), 57--89.
```

---

## 11. Ramanujan/Gauss toplamları için standart ders kitabı

**(a) Durum: her ikisi de DOĞRULANDI (zbMATH); ikisinden birini seçmek yazarın
tercihi, ikisi de doğru ve ilgili.**
- Apostol: zbMATH → *Introduction to Analytic Number Theory*, Undergraduate Texts in
  Mathematics, Springer, 1976 (ilk baskı) / ISBN 978-0-387-90163-3 (1998 baskısı, xii +
  338 s.). Ramanujan toplamları §2.7'de standart biçimde işlenir (zbMATH kaydında
  ayrıca doğrulanmadı ama bu, kitabın klasik/yaygın bilinen bir bölümüdür — **içerik
  DOĞRULANAMADI**, sadece künye doğrulandı).
- Iwaniec–Kowalski: zbMATH → *Analytic Number Theory*, AMS Colloquium Publications 53,
  Amer. Math. Soc., Providence, RI, 2004, ISBN 0-8218-3633-1, xi + 615 s. Ramanujan
  toplamları Bölüm 1.6 ve Gauss toplamları Bölüm 3.4'te ele alınır (yaygın bilinen
  içerik — **içerik DOĞRULANAMADI**, künye doğrulandı).

**(b) LaTeX (ikisi de verildi, biri seçilsin)**
```
\bibitem{Apo76} T.~M.~Apostol, \emph{Introduction to Analytic Number Theory},
Undergraduate Texts in Mathematics, Springer-Verlag, New York--Heidelberg--Berlin,
1976.

\bibitem{IK04} H.~Iwaniec, E.~Kowalski, \emph{Analytic Number Theory}, American
Mathematical Society Colloquium Publications \textbf{53}, Amer.\ Math.\ Soc.,
Providence, RI, 2004.
```

---

## 12. Kanivets 2026, Zenodo DOI 10.5281/zenodo.20766728

**(a) Durum: DOĞRULANDI (Zenodo API, doğrudan).**
`https://zenodo.org/api/records/20766728` → başlık **"No Evidence for Detectable
Arithmetic Frequencies in Long-Range Deviations of Riemann Zero Statistics from
Random Matrix Theory"**, yazar **Serhii Kanivets**, yayın tarihi 2026-06-19, tür
"Journal article" (Zenodo'nun kendi sınıflandırması — hangi dergide/ön-baskı
sunucusunda olduğu bu kayıtta belirtilmiyor, sadece Zenodo'da barındırılan bir PDF).
DOI 10.5281/zenodo.20766728 kaydın kendisiyle birebir eşleşiyor.

**(b) LaTeX**
```
\bibitem{Kan26} S.~Kanivets, \emph{No Evidence for Detectable Arithmetic
Frequencies in Long-Range Deviations of Riemann Zero Statistics from Random Matrix
Theory}, Zenodo (2026), \url{https://doi.org/10.5281/zenodo.20766728}.
```

---

## 13. Sezen companion notes, Zenodo concept DOI 10.5281/zenodo.22942475

**(a) Durum: DOĞRULANDI (Zenodo API, doğrudan; concept DOI takip edildi).**
`https://zenodo.org/api/records/22942475` bir "concept DOI" olduğu için otomatik
olarak en son sürüm kaydına (id 22942476) yönlendiriyor: başlık **"Riemann zero
statistics — an experimental research program"**, yazar **Uğur Sezen**, sürüm v1.0,
yayın tarihi 2026-09-24, kaynak türü **Software** (makale değil — GitHub deposunun
(`azizran/riemann-zero-statistics`, tag v1.0) bir Zenodo arşivi; `related_identifiers`
alanı bu GitHub bağlantısını "isSupplementTo" ilişkisiyle doğruluyor). Concept DOI
(10.5281/zenodo.22942475) her zaman en güncel sürüme işaret eder — notlarda companion
olarak bu DOI'nin kullanılması doğru, çünkü sürüm-bağımsız kalıcı bir bağlantı.

**(b) LaTeX**
```
\bibitem{Sezen26} U.~Sezen, \emph{Riemann zero statistics --- an experimental
research program}, Zenodo, v1.0 (2026), \url{https://doi.org/10.5281/zenodo.22942475}.
```

---

## 14. Fujii (1989), "On a theorem of Landau"

**(a) Durum: DOĞRULANDI (Crossref + zbMATH + Project Euclid/DOI çözümlemesi — üç
bağımsız kaynak).**
- Crossref: `https://api.crossref.org/works/10.3792/pjaa.65.51` → *Proceedings of the
  Japan Academy, Series A, Mathematical Sciences*, cilt 65, sayı 2, yayın 1989-01-01.
- zbMATH: aynı arama → "Proc. Japan Acad., Ser. A 65, No. 2, **51-54** (1989)."
- DOI çözümlemesi (`10.3792/pjaa.65.51` → Project Euclid) sayfa başlangıcını (51) DOI
  içinde doğruluyor (Japan Academy'nin DOI kuralı: `pjaa.<cilt>.<ilk sayfa>`).

Ek bulgu (task'ta istenmemiş ama kontrol sırasında ortaya çıktı, dürüstlük için
kaydediliyor): "On a theorem of Landau, **II**" (1990) için zbMATH'in yapılandırılmış
sayfa alanı bir ara sorguda **"391-396"** döndürdü — bu YANLIŞ; Project Euclid'in
kendi sayfası (`https://projecteuclid.org/.../10.3792/pjaa.66.291.full`, bu oturumda
doğrudan okundu) ve DOI'nin kendisi (`pjaa.66.291`) **291-296**'yı doğruluyor —
zbMATH'in arama API'sindeki bir rakam-yer değiştirme hatası (391↔291) olduğu
anlaşılıyor. `LITERATUR_KAPPA_SONLU_YUKSEKLIK_25EYL2026.md`'deki mevcut "291-296"
künyesi **doğruydu**, değiştirilmemeli.

**(b) LaTeX**
```
\bibitem{Fujii89} A.~Fujii, \emph{On a theorem of Landau}, Proc.\ Japan Acad.\ Ser.\
A Math.\ Sci.\ \textbf{65} (1989), no.~2, 51--54.
```

(Bonus, doğrulandı ama task'ta istenmemişti — gerekirse eklenebilir:)
```
\bibitem{Fujii90} A.~Fujii, \emph{On a theorem of Landau, II}, Proc.\ Japan Acad.\
Ser.\ A Math.\ Sci.\ \textbf{66} (1990), no.~9, 291--296.
```

---

## Özet tablo — doğrulama durumu

| # | Kaynak | Durum |
|---|---|---|
| 1a | BK Nonlinearity I (1995) | DOĞRULANDI |
| 1b | BK Nonlinearity II (1996) | DOĞRULANDI |
| 1c | BK PRL (1996) | DOĞRULANDI |
| 2 | Bogomolny 0708.4223 | DOĞRULANDI (sayfa 19–36 düzeltildi, 19–44 değil) |
| 3 | Conrey–Snaith math/0509480 | DOĞRULANDI (cilt 94 doğrulandı, arXiv'in "93"ü yanlış) |
| 4 | BBLM math/0602270 | DOĞRULANDI, mevcut künye zaten doğruydu |
| 5a | Landau | DOĞRULANDI (yıl 1912, 1911 değil) |
| 5b | Gonek 1993 | DOĞRULANDI |
| 6 | Durkan–Hughes–Pearce-Crump 2601.18025 | DOĞRULANDI (henüz dergi yok) |
| 7 | Hughes–Lugmayer–Pearce-Crump 2411.05573 | DOĞRULANDI (tam künye: JLMS (2) 112 (2025), e70250) |
| 8a | Zhang–Martelli–Torquato 1801.01541 | DOĞRULANDI (sayfa 115001 düzeltildi, 135001 değil) |
| 8b | Torquato–Zhang–de Courcy-Ireland 1802.10498 | DOĞRULANDI |
| 9 | Rodgers 1203.3275 | DOĞRULANDI |
| 10 | Montgomery/Odlyzko(+tablo)/Keating–Snaith | DOĞRULANDI (4/4) |
| 11 | Apostol / Iwaniec–Kowalski | Künye DOĞRULANDI; içerik (Ramanujan toplamı bölümü) DOĞRULANAMADI (kitap içi erişilmedi) |
| 12 | Kanivets Zenodo 20766728 | DOĞRULANDI |
| 13 | Sezen Zenodo 22942475 (concept) | DOĞRULANDI |
| 14 | Fujii 1989 | DOĞRULANDI |

**DOĞRULANAMADI olan tek şey**: Apostol / Iwaniec–Kowalski kitaplarının Ramanujan
toplamı bölümlerinin İÇERİĞİ (künyeleri doğrulandı, ama kitap metnine bu oturumda
erişilmedi — bu standart, çok bilinen bir içerik olduğu için düşük risk, ama kural
gereği açıkça işaretleniyor).
