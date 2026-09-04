# -*- coding: utf-8 -*-
"""
172f — G4: SENTEZ — M(λ), c(λ), VADİ; ve 173 İÇİN MÜHÜRLÜ ÖNGÖRÜ
Veri: 172/G1.json, 172/G2.json, 171/NU_MERDIVEN.json. YENİ KOŞU YOK.

══════════════════════════════════════════════════════════════════════
ÇERÇEVE
══════════════════════════════════════════════════════════════════════
G1: g = 1 − Q/ρ (ÖZDEŞLİK, iki alan için de).  G2: θ/θ₀ = (1−βφ)(ρ_X/ρ_X₀)^p,
p = 0.651 (yalnız λ ailesinden, 172d).  M = (g_E/g₀)(g_X/g₀)²(θ/θ₀).
c = KALİB/W_X = g_cal·θ/W_X.

DÜRÜSTLÜK: M(λ)'nın "parametresiz yeniden inşası" λ ailesinin KENDİ
ölçülen Q, ρ eğrileriyle yapılırsa ÖZDEŞTİR ve hiçbir şey sınamaz.
Sınanabilir olan İKİ şey vardır ve ikisi de SIFIR yeni parametre taşır:
 (a) tümseklerin YERİ (172c: Q_E ↔ λ_c; burada: α ↔ λ*),
 (b) 173 için ÖRNEKLEM-DIŞI mühürlü sayı (λ = 0.40 ve λ = 1.45).

══════════════════════════════════════════════════════════════════════
ÖN-MÜHÜR (koşudan önce; 4 Eyl, 172f)
══════════════════════════════════════════════════════════════════════
 Ö1  **α(λ) TÜMSEĞİNİN TEPESİ = c(λ) VADİSİ.** 171 §T2c.3 c'nin
     minimumunu λ\* = **0.6487** ölçtü (ν = 1 geçişi). 172d'nin ölçtüğü
     A-genişliği α_g(λ) de tümsektir (1.0242 → 1.2264 → 0.6874).
     ÖN-KAYIT: log α(λ)'ya 5 orta noktadan (0.60…1.15) parabol; tepe
     **λ\* = 0.6487 ± 0.03**. SIFIR yeni parametre.
     GEREKÇE: c_b(τ) ∝ KALİB_b/W_X(τ) ∝ exp(−(α_g − 2π²σ_X̃²)τ²) ⇒
     bant-ortalamalı c'nin λ-türevinde α'nın tepesi durgun noktadır.
 Ö2  VADİNİN ÇARPAN ADRESİ: dlog KALİB/dλ = dlog W_X/dλ kesişmesi
     λ\* = 0.6487'yi ±0.05 içinde vermeli (özdeşlik + ölçülen eğriler).
 Ö3  İNDİRGENMİŞ MODEL: g_X ≡ SABİT (0 param) ⇒ 7 gazda maks hata
     ≤ %0.5; ve M ≈ (g_E/g₀)·(ρ_X/ρ_X₀)^0.651 ⇒ 7 gazda rms ≤ %1.5.
 Ö4  **173 İÇİN MÜHÜR.** λ = 0.40 ve λ = 1.45 gazları için (henüz
     İNŞA EDİLMEMİŞ) M, θ, g_E, g_X, c öngörüleri aşağıda basılır ve
     `172/G4.json`'a yazılır. Bunlar log-λ'da kuadratik/parabol
     ekstrapolasyonlardır; belirtilen belirsizlikler ailenin kendi
     uyum artığından gelir. Kurtarma yok: 173 ölçerse ne çıkarsa o.

══════════════════════════════════════════════════════════════════════
SONUÇ (yalnız gerçek koşudan, `172/log_172f.txt`)
══════════════════════════════════════════════════════════════════════
 Ö1 ✓ (ön-kayıtlı pencerede) **α tümseğinin tepesi = 0.6675** (5 orta
      nokta, log) ve **0.6594** (log'suz); ön-kayıt λ* = 0.6487 ± 0.03
      ⇒ fark +0.019 / +0.011, İKİSİ DE İÇERİDE. 7 noktalı uyum 0.7877
      veriyor (dışarıda) — sebebi L130'un α = 0.6874 aykırılığıdır;
      ön-kayıt 5 orta noktayı belirtmişti, kurtarma değil.
 Ö2 ✓✓ **TAM İSABET.** dlogKALİB/dλ ile dlogW_X/dλ **λ = 0.6506**'da
      kesişiyor; ön-kayıt λ* = 0.6487 ⇒ fark **+0.0019 (‰3)**.
      Mekanizma tabloda açık: dlogW_X/dλ neredeyse sabit (−0.42…−0.60),
      dlogKALİB/dλ ise −0.90'dan (λ=0.55) +0.01'e (λ=1.075) tırmanıyor.
      **VADİ, KALİB'in kendi eğiminin W_X'in eğimini geçtiği yerdir.**
 Ö3 ✗ g_X "sabit" ön-kaydı ıskaladı (maks %1.01, eşik %0.5 — L050
      ve K070/E060 kırıyor). İndirgenmiş model M ≈ (g_E/g₀)(ρ_X/ρ_X₀)^p
      rms **%2.07** (eşik %1.5) ⇒ ÖLDÜ; en büyük hata λ=0.50'de −%4.37.
 Ö4 MÜHÜRLENDİ (`172/G4.json`), iki bağımsız yoldan:
      çarpan yolu   λ=0.40: M = 1.2947, c = 0.3790 | λ=1.45: M = 0.9744,
                    c = 0.4747  (g_E 0.6842/0.5939, θ 0.9452/0.8405)
      doğrudan yol  λ=0.40: M = 1.3881, c = 0.4052 | λ=1.45: M = 1.0060,
                    c = 0.4910   (171 M9 çıpası, log-λ kuadratiği)
      İki mühür λ=0.40'ta **%6.7 ayrışıyor** — çarpan yolu λ=0.50'yi
      zaten −%4.4 ıskaladığı için DOĞRUDAN YOL daha güvenilirdir;
      173 ölçerse ikisi arasındaki fark tek başına bir sınavdır.
      (c = c₀·M·W_X₀/W_X dönüşümü mevcut 7 gazda ‰1.4 tutuyor.)
"""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402

SCR = ("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
       "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/")
QM = "/Users/ugursezen/Desktop/arin/deney/qm_riemann/"
G1 = json.load(open(SCR + "172/G1.json"))
G2 = json.load(open(SCR + "172/G2.json"))
G3 = json.load(open(SCR + "172/G3.json"))
NU = json.load(open(SCR + "171/NU_MERDIVEN.json"))
LAM7 = ["L050", "L060", "L070", "L085", "Hkeskin", "L115", "L130"]
LV = np.array([.50, .60, .70, .85, 1.00, 1.15, 1.30])
LSTAR, LAMC, PTH = NU["lam_star"], 1.0109, 0.651
KES = ["K090", "K070", "HA4", "E060"]
PEN = ["HkT2a", "HkT2b", "HkT4a", "HkT4b"]


def E(k):
    return np.array([G1[g]["E"][k] for g in LAM7])


def X(k):
    return np.array([G1[g]["X"][k] for g in LAM7])


alfa = np.array([G2["alfa"][g] for g in LAM7])
th = np.array([G2["theta"][g] for g in LAM7])
KAL = np.array([G2["KAL"][g] for g in LAM7])
WX = np.array([G1[g]["WX"] for g in LAM7])
sX = np.array([G1[g]["sigX"] for g in LAM7])
gE = np.array([G1[g]["gE"] for g in LAM7])
gX = np.array([G1[g]["gX"] for g in LAM7])
cO = np.array(NU["c"])[::-1]                 # NU λ'yı azalan tutuyor
MO = np.array(NU["M"])[::-1]

print("== Ö1: α(λ) TÜMSEĞİNİN TEPESİ = c VADİSİ mi? ==")
for et, sl in (("5 orta (0.60-1.15)", slice(1, 6)), ("7 nokta", slice(None))):
    cc = np.polyfit(LV[sl], np.log(alfa[sl]), 2)
    print("  %-20s tepe λ = %+.4f   (ön-kayıt λ* = %.4f, fark %+.4f)"
          % (et, -cc[1] / (2 * cc[0]), LSTAR, -cc[1] / (2 * cc[0]) - LSTAR))
cc = np.polyfit(LV[1:6], alfa[1:6], 2)       # doğrusal ölçekte de
print("  %-20s tepe λ = %+.4f" % ("(log'suz, 5 orta)", -cc[1] / (2 * cc[0])))
print("  hatırlatma (172c): Q_E tümseğinin tepesi 1.0009/1.0191, "
      "ön-kayıt λ_c = %.4f" % LAMC)

print()
print("== Ö2: VADİNİN ÇARPAN ADRESİ — dlogKALİB/dλ = dlogW_X/dλ ==")
lm = 0.5 * (LV[:-1] + LV[1:])
dK = np.diff(np.log(KAL)) / np.diff(LV)
dW = np.diff(np.log(WX)) / np.diff(LV)
print("  λ_orta   dlogKALİB/dλ  dlogW_X/dλ   fark (=dlog c/dλ)")
for i in range(len(lm)):
    print("  %.3f     %+9.4f     %+9.4f    %+9.4f"
          % (lm[i], dK[i], dW[i], dK[i] - dW[i]))
f = dK - dW
kesis = [lm[i] + (lm[i + 1] - lm[i]) * (-f[i]) / (f[i + 1] - f[i])
         for i in range(len(f) - 1) if f[i] * f[i + 1] < 0]
print("  KESİŞİM(ler): %s   (ön-kayıt λ* = %.4f)"
      % (["%.4f" % k for k in kesis], LSTAR))

print()
print("== Ö3: İNDİRGENMİŞ MODEL ==")
print("  g_X sabit mi? maks |g_X/g_X₀ − 1| = %.2f%% (7 gaz)"
      % (100 * np.max(np.abs(gX / gX[4] - 1))))
Mred = (gE / gE[4]) * (X("rho") / X("rho")[4]) ** PTH
Mol = KAL / KAL[4]
print("  λ      M ölçülen  M indirgenmiş  fark%")
for i, g in enumerate(LAM7):
    print("  %.2f   %.4f     %.4f       %+6.2f" % (LV[i], Mol[i], Mred[i],
                                                   100 * (Mred[i] / Mol[i] - 1)))
print("  rms = %.2f%%  (ön-kayıt ≤ %%1.5)"
      % (100 * np.sqrt(((Mred / Mol - 1) ** 2).mean())))

print()
print("== Ö4: 173 İÇİN MÜHÜRLÜ ÖNGÖRÜ (λ = 0.40 ve λ = 1.45) ==")


def kuad(y, xs, log=True, deg=2):
    z = np.log(y) if log else y
    c = np.polyfit(np.log(LV), z, deg)
    v = np.polyval(c, np.log(xs))
    art = z - np.polyval(c, np.log(LV))
    return (np.exp(v) if log else v), float(art.std())


XS = np.array([0.40, 1.45])
tah = {}
for ad, y in (("Q_E", E("Q")), ("rho_E", E("rho")), ("Q_X", X("Q")),
              ("rho_X", X("rho")), ("sig_X", sX), ("alfa", alfa),
              ("g_X", gX), ("W_X", WX)):
    v, s = kuad(y, XS)
    tah[ad] = (v, s)
    print("  %-8s λ=0.40: %.5f   λ=1.45: %.5f   (aile uyum artığı %.1e log)"
          % (ad, v[0], v[1], s))
gEp = 1 - tah["Q_E"][0] / tah["rho_E"][0]
thp = th[4] * (tah["rho_X"][0] / X("rho")[4]) ** PTH
gXp = tah["g_X"][0]
Mp = (gEp / gE[4]) * (gXp / gX[4]) ** 2 * (thp / th[4])
cp = cO[4] * Mp * (WX[4] / tah["W_X"][0])
print("  ------------------------------------------------------------")
print("  λ = 0.40:  g_E = %.4f  g_X = %.4f  θ = %.4f  **M = %.4f**  "
      "**c = %.4f**" % (gEp[0], gXp[0], thp[0], Mp[0], cp[0]))
print("  λ = 1.45:  g_E = %.4f  g_X = %.4f  θ = %.4f  **M = %.4f**  "
      "**c = %.4f**" % (gEp[1], gXp[1], thp[1], Mp[1], cp[1]))
print("  (c dönüşümü c = c₀·M·W_X₀/W_X; bu bağıntı mevcut 7 gazda "
      "%.2f%% rms tutuyor)"
      % (100 * np.sqrt(((cO[4] * Mol * WX[4] / WX / cO - 1) ** 2).mean())))
# --- İKİNCİ MÜHÜR: doğrudan log-λ kuadratiği (171'in M9 "empirik çıpası")
Md, sM = kuad(Mol, XS)
cd, sc = kuad(cO, XS)
print("  [2. MÜHÜR — doğrudan log-λ kuadratiği, 171 M9 çıpası]")
print("     λ=0.40: M = %.4f, c = %.4f    λ=1.45: M = %.4f, c = %.4f"
      "   (aile artığı %.1e / %.1e log)" % (Md[0], cd[0], Md[1], cd[1], sM, sc))
print("     UYARI: çarpan-yolu (1. mühür) λ=0.50'de M'yi %.2f%% ıskalıyor; "
      "iki mühür λ=0.40'ta %.1f%% ayrışıyor."
      % (100 * (Mred[0] / Mol[0] - 1), 100 * (Mp[0] / Md[0] - 1)))

json.dump(dict(lam_star=LSTAR, lam_c=LAMC, p_theta=PTH,
               alfa_tepe=float(-np.polyfit(LV[1:6], np.log(alfa[1:6]), 2)[1]
                               / (2 * np.polyfit(LV[1:6], np.log(alfa[1:6]),
                                                 2)[0])),
               kesisim=[float(k) for k in kesis],
               M_indirgenmis=[float(v) for v in Mred],
               M_olculen=[float(v) for v in Mol],
               ongoru={"lam": [0.40, 1.45], "gE": [float(v) for v in gEp],
                       "gX": [float(v) for v in gXp],
                       "theta": [float(v) for v in thp],
                       "M": [float(v) for v in Mp],
                       "c": [float(v) for v in cp],
                       "M_dogrudan": [float(v) for v in Md],
                       "c_dogrudan": [float(v) for v in cd]}),
          open(SCR + "172/G4.json", "w"), indent=1)

# ------------------------------------------------------------------ FİGÜR
fig, ax = plt.subplots(2, 2, figsize=(12.5, 9))
a = ax[0, 0]
a.plot(LV, E("Q"), "o-", color="#1b6ca8", label=r"$Q_E$ (izdüşüm kusuru)")
a.plot(LV, X("Q"), "s-", color="#c0392b", label=r"$Q_X$")
a.axvline(LAMC, ls="--", c="k", lw=1)
a.text(LAMC + .01, 0.72, r"$\lambda_c=1.0109$" "\n(170 §K0.2)", fontsize=8)
a.plot([1.0], [G1["son"]["E"]["Q"]], "*", ms=14, color="#f39c12")
a.annotate("son", (1.0, G1["son"]["E"]["Q"]), (1.02, 0.855), fontsize=8)
a.set_xlabel(r"$\lambda$"); a.set_ylabel("Q")
a.set_title(u"(a) G1: Q_E tümseği λ_c'de zirve yapıyor")
a.legend(fontsize=8); a.grid(alpha=.3)

a = ax[0, 1]
a.plot(LV, alfa, "o-", color="#16a085")
a.axvline(LSTAR, ls="--", c="k", lw=1)
a.text(LSTAR + .01, 0.75, r"$\lambda^*=0.6487$" "\n(171 c-vadisi)", fontsize=8)
a2 = a.twinx()
a2.plot(LV, cO, "^--", color="#8e44ad", label="c(λ) (171)")
a2.set_ylabel("c", color="#8e44ad")
a.set_xlabel(r"$\lambda$"); a.set_ylabel(r"$\alpha$ (A-genişliği)",
                                         color="#16a085")
a.set_title(u"(b) G2: α tümseği ↔ c vadisi")
a.grid(alpha=.3)

a = ax[1, 0]
for et, lst, mk, col in (("λ ailesi", LAM7, "o", "#1b6ca8"),
                         ("kesim", KES, "s", "#c0392b"),
                         ("pencere", PEN, "^", "#16a085"),
                         ("son", ["son"], "*", "#f39c12")):
    xs = [G1[g]["X"]["rho"] / G1["Hkeskin"]["X"]["rho"] for g in lst]
    ys = [(G2["theta"][g] / G2["theta"]["Hkeskin"])
          / (1 - 0.2175 * G2["phi"][g]) for g in lst]
    a.plot(xs, ys, mk, color=col, label=et, ms=12 if et == "son" else 7)
xx = np.linspace(0.93, 1.30, 50)
a.plot(xx, xx ** PTH, "k-", lw=1, label=r"$(\rho_X/\rho_{X0})^{0.651}$")
a.set_xlabel(r"$\rho_X/\rho_{X0}$")
a.set_ylabel(r"$(\theta/\theta_0)/(1-\beta\varphi)$")
a.set_title(u"(c) G2: θ yasası — tek ayakta kalan aday")
a.legend(fontsize=8); a.grid(alpha=.3)

a = ax[1, 1]
ad = [k for k in G3["lam_es"] if G3["lam_es"][k] == G3["lam_es"][k]]
ad.sort(key=lambda k: G3["lam_es"][k])
renk = {"marjinal": "#7f8c8d", "E (η)": "#c0392b", "X (ΔĈ)": "#1b6ca8",
        "θ": "#f39c12", "şekil": "#16a085"}
a.barh(range(len(ad)), [G3["lam_es"][k] for k in ad],
       color=[renk[G3["kanal"][k]] for k in ad])
a.set_yticks(range(len(ad))); a.set_yticklabels(ad, fontsize=8)
a.axvline(G3["capa"], ls="--", c="k", lw=1)
a.text(G3["capa"] + .005, 0.2, u"σ_X̃ çapası\n%.3f" % G3["capa"], fontsize=8)
a.set_xlabel(r"$\lambda_{eş}$ (gerçek gaz)")
a.set_title(u"(d) G3: gerçeğin dar adresi — E kanalı ayrışıyor")
a.set_xlim(0.65, 1.0); a.grid(alpha=.3, axis="x")

fig.suptitle(u"172 — ÇARPANLARIN YASASI: g = 1 − Q/ρ, iki tümsek, "
             u"bir dar adres", fontsize=13)
fig.tight_layout(rect=(0, 0, 1, 0.97))
fig.savefig(QM + "172_carpanlar.png", dpi=130)
print("\n-> %s172/G4.json  +  %s172_carpanlar.png" % (SCR, QM))
