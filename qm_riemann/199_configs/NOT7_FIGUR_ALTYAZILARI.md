# Not 7 — figür altyazı taslakları (LaTeX, İngilizce)

Figürler `qm_riemann/fig_n7_*_en.{png,pdf}` (6.5 in genişlik, 300 dpi PNG + vektör PDF).
Betikler `199_configs/199fig_1.py … 199fig_5.py`; hepsi yalnız mühürlü json/npz okur, yeni
ölçüm yapmaz. Altyazıdaki her sayı betiklerin ekrana bastığı ya da assert ettiği değerlerle
ve ilgili raporla (190/193/195/197/198) birebir kontrol edildi. Her altyazının altında kısa
Türkçe kaynak/uyarı notu var (yayına girmez).

---

## Figür 1 — `fig_n7_positions_en`  (kaynak: 190)

```latex
\begin{figure}
  \centering
  \includegraphics[width=\columnwidth]{fig_n7_positions_en.pdf}
  \caption{\textbf{Satellite positions are universal in $\Delta\omega$.}
  (a)~Pooled cancellation-kernel density $\kappa$, normalized to its maximum on
  $-1.3\le\Delta\omega\le2.3$, versus $\Delta\omega=\omega'-L_{\rm loc}$ for two windows of
  zeros, $L=10.48$ (solid, bins of $0.025$) and $L=12.03$ (dashed, $\tau'$ bins of width
  $0.060$ in $\omega$); vertical lines mark satellites at $\log(a/b)$, solid for the
  pre-registered targets $0,\pm\log2,\pm\log3,\log6$ and dotted for the other satellites
  that are bright in both windows.
  In the $L=10.48$ window the pre-registered block-median peaks lie at $0.025$, $0.725$,
  $1.113$ and $1.800$ (targets $0,\log2,\log3,\log6$), each within $0.05$ of its target and
  away from the $\tau'$-fixed alternatives $0.604$, $0.957$, $1.562$, and the two profiles
  correlate at $0.980$.
  (b)~On the $\tau'$-scaled axis $\Delta\omega/L$ the same profiles no longer coincide
  (correlation $0.002$; ticks mark $\log n/L$ for $n=2,3,6$ in each window).
  The shaded strip $\Delta\omega<-1.07$ is only partly covered by the eight blocks of the
  $L=10.48$ window.}
  \label{fig:n7-positions}
\end{figure}
```

*Not:* profil kurgusu 190d/190e AYNEN; korelasyonlar yeniden hesaplanıp `K_kesif_190.json`
ile 1e-12'de assert edildi (0.97976 / 0.00170). Tepe medyanları HUKUM_190 H190a_tablo.
Noktalı uydular (±log 3/2, +log 5/4, +log 5, +log 10) 192 B kataloğunda iki pencerede de
"yanar" olanlardır; 190'da ±0.4 çiftinin log(3/2) etiketi "SINANMADI" idi, konum etiketi
192'den geliyor. Gölge: düşük pencerede Δω < 0.86·L − min L_b = −1.071 bölgesinde 8 bloğun
hepsi kapsamıyor (190d yoğunluğu kapsamayan blokları 0 sayar) — −log 3 tepesi bu yüzden
düşük pencerede görece küçük.

---

## Figür 2 — `fig_n7_classes_en`  (kaynak: 193)

```latex
\begin{figure}
  \centering
  \includegraphics[width=\columnwidth]{fig_n7_classes_en.pdf}
  \caption{\textbf{Residue-class law of the satellites.}
  Class shares $s_r=\kappa_r/\kappa_{\rm tot}$ of the lines $q'\equiv r\pmod a$ (points,
  8-block jackknife errors) against the blind prediction
  $s_r=\cos(2\pi rb/a)/\mu(a)$ (bars), window $L=12.03$.
  (a)~At $\Delta\omega=+\log10$ all four classes have the predicted sign and lie within
  $0.11$ of the prediction: $s_1=0.899\pm0.023$ and $s_9=0.902\pm0.017$ (predicted $0.809$),
  $s_3=-0.416\pm0.035$ and $s_7=-0.386\pm0.035$ (predicted $-0.309$).
  (b)~At $+\log7$ the signs are again correct ($|z|=4.9$--$6.2$ for $r=1,3,4,6$) but the
  magnitudes are $1.8$--$2.4$ times the prediction, while the controls $+\log3$
  ($0.462/0.538$) and $-\log3$ ($0.499/0.501$, no class structure expected) sit at
  $\tfrac12/\tfrac12$.}
  \label{fig:n7-classes}
\end{figure}
```

*Not:* değerler `K1_sinif.json`; öngörüler cos(2πrb/a)/μ(a)'dan yeniden hesaplanıp `s_pred`
ile assert edildi. +log10 farkları 0.090/0.107/0.077/0.093 (maks 0.107 → "within 0.11").
+log7'de r2/r5 (hükme girmeyen yan-KAYIT) de çizili: 0.332±0.096 ve 0.264±0.095 (öngörü
0.223). Oran ölçülen/öngörü: r1 2.31, r3 1.89, r4 1.84, r6 2.44 → altyazıda "1.8–2.4".

---

## Figür 3 — `fig_n7_characters_en`  (kaynak: 195)

```latex
\begin{figure}
  \centering
  \includegraphics[width=\columnwidth]{fig_n7_characters_en.pdf}
  \caption{\textbf{The character selects the satellites of Dirichlet $L$-function zeros.}
  (a)~$\mathrm{Re}\,\tilde K$ (symlog scale, jackknife errors) at nine satellites for five
  islands of zeros ($\chi_4,\chi_3,\chi_{5e},\chi_{8e},\chi_{8o}$; $L_\chi\ge8.5$) and for
  $\zeta$ (stars; low window on the same $\Delta L$ block grid); filled symbols are allowed
  satellites ($n=a/b$ has no prime factor $p\mid k$), open symbols forbidden ones.
  Allowed satellites with $\mu(a)\neq0$ carry the sign of $\zeta$ in both parities
  ($z=-55$ to $-82$ for the six primary pairs; $\chi_{8o}/\chi_{8e}=1.009\pm0.024$ at
  $+\log3$), whereas forbidden satellites leave only a small residual of opposite sign,
  $+0.0014$ to $+0.0063$, at the level of the $\mu(4)=0$ position $+\log4$ (positive also
  for $\zeta$, $\chi_3$ and $\chi_{5e}$).
  (b)~Phases: $\angle\tilde K=180^\circ\pm1.3^\circ$ for every allowed $\mu\neq0$ satellite,
  odd or even, and $0^\circ\pm7.2^\circ$ for every forbidden one.}
  \label{fig:n7-characters}
\end{figure}
```

*Not:* ada değerleri `B_cekirdek.json` (üst bölge), ζ `K0c.json` → `dusuk_izgara_tum`
(raporun K0c tablosundaki "K̃_ζ düşük (ΔL ızgarası)" satırı). Gösterilen büyüklük B
çekirdeğidir (Re K̃); A'nın güç z'leri çizilmedi (yasak z'leri negatif ve null kalibre
değil — rapor uyarısı). Betik denetimi: izinli maks |∠−180°| = 1.348°, maks |Im z| = 2.07;
yasak (n = 24) maks |∠| = 7.17°, Re K̃ ∈ [+0.00144, +0.00629]; ρ₈ HUKUM_195 ile assert.
−log5 ön-kayıtla B'den dışlandığı için yok.

---

## Figür 4 — `fig_n7_mirror_en`  (kaynak: 197)

```latex
\begin{figure}
  \centering
  \includegraphics[width=\columnwidth]{fig_n7_mirror_en.pdf}
  \caption{\textbf{Mirror law in the blind band.}
  (a)~Kernel $\kappa_\Sigma$ (32 blocks, window $L=10.48$, jackknife errors) in the
  pre-registered, previously unviewed band $-2.12\le\Delta\omega\le-1.58$, with the
  free-amplitude fit $M_{\rm BK}$ and the parameter-free Bogomolny--Keating prediction
  $A_x=s/x$ ($s$ from the calibration band $-1.30\le\Delta\omega\le-0.55$; same linear
  baseline and fixed quarter-family terms); solid and dotted grid lines mark integer and
  half-integer $x$ at $\Delta\omega=-\log x$.
  (b)~Amplitude ratios $R_x=A_x/(s/x)$: the integer family gives
  $\bar R_{\mathbb Z}=1.023\pm0.023$ ($\sigma_{\rm eff}$), whereas the line-density rival
  predicts $\bar R^g=0.706$ and its own fit (grey) returns $0.972$, more than
  $12\sigma_{\rm eff}$ away.
  The half-integer family follows the same $1/x$ envelope at a reduced level,
  $\rho=0.789\pm0.037$, and the prime power $x=8$ is not suppressed
  ($R_8=1.048\pm0.060$, $\psi=0.993\pm0.039$).}
  \label{fig:n7-mirror}
\end{figure}
```

*Not:* 197c_hukum.py importlib ile yüklendi (düzenlenmedi); `hesap(tam=False)` mühürlü
`harita_omega32.npz` üzerinde ön-kayıtlı uyumu yeniden koşar (M6 simülasyonu yok). s, tüm
R_x (değer ve jk se, iki model) ve R̄_Z HUKUM_197 ile 1e-10'da assert edildi. Hata
çubukları jk se; R̄_Z ve ρ σ_eff = f·σ_jk (f_Z = 1.565, f_ρ = 1.374). (R̄_Z^Rg − R̄^g)/σ_eff
= 12.3. Kör bantta 21 dilim; rms(veri − M_BK) 1.4e-4, rms(veri − BK öngörüsü) 5.0e-4
(medyan jk se 1.8e-4) — fark büyük ölçüde buçuklu tepelerden.

---

## Figür 5 — `fig_n7_height_en`  (kaynak: 198)

```latex
\begin{figure}
  \centering
  \includegraphics[width=\columnwidth]{fig_n7_height_en.pdf}
  \caption{\textbf{Height dependence of the positive-exponent factors.}
  $\kappa_p=A(+\log p)/(s\,c_{\rm BK}(p))$ for $p=2,3,5,7$ in three windows,
  $L_W=9.34$ (blind), $10.48$ (reference) and $12.03$, with jackknife errors.
  Solid: joint fit $1-\kappa_p=B_p(L/L_d)^{-\gamma}$ with $L_d=10.48$, giving
  $\hat\gamma=1.36\pm0.21$ ($\sigma_{\rm eff}$); dashed: finite-height hypothesis $H_S$
  ($\gamma=1$); dotted: constant structure $H_C$ ($\gamma=0$), both anchored at the
  pre-registered values $\kappa_p=0.795,\,0.720,\,0.608,\,0.530$.
  $\hat\gamma$ lies $1.7\sigma$ from $H_S$ and $6.4\sigma$ from $H_C$: the factors relax
  towards the Bogomolny--Keating value $\kappa_p=1$ as the height grows.}
  \label{fig:n7-height}
\end{figure}
```

*Not:* κ_p ± jk se, γ̂, σ_eff ve B_p `HUKUM_198.json`; H_S/H_C çapası κ_p^{197} =
A(+log p)/(s·c_BK) `HUKUM_197.json` İKİNCİL'den (0.7953/0.7205/0.6084/0.5296) — ön-kayıtlı
öngörü tablosunun çapası budur. 198c'nin kendi figürü H_S/H_C'yi 198'in düşük-pencere
ölçümüne (0.789/0.713/0.604/0.538, eşit-ΔL geometri) çapalıyordu; bu figür ön-kayıtlı çapayı
kullanır (fark ≤ 0.009, ≤ 0.6 se). |γ̂−1|/σ = 1.69, |γ̂|/σ = 6.36. Çapraz-yükseklik
ilişkisi için 198 ŞERHİ geçerli: bu L-bağımlılığı BK'nın düzeltmesi değil, gözlenebilirin
aktarım özelliği — altyazı "relax towards the BK value" diyerek yalnız ölçüleni söyler.
