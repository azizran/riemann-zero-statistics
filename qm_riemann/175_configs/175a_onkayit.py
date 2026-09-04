# -*- coding: utf-8 -*-
"""
175a — ÖN-KAYIT (K2) + K1 İNŞA KAPILARININ DENETİMİ
====================================================
KALEM_KESIM_VE_KILIT_04EYL2026: fazlanın ayrışımı = kesim payı + kilit payı.

Bu betik ÖLÇÜMDEN ÖNCE koşar, kendi sha256'sını yazar ve dosyayı bir daha
yazmaz. İçinde ÜÇ şey vardır:
  (A) erfc zarfının TAM formülü — 152'nin raporundan ve 164_insa.merdiven
      kaynağından BİREBİR doğrulanır (iddia değil, kod denetimi);
  (B) K1 inşa kapılarının denetimi (erfc-ikiz = `HA4`, 164'ün SADAKATLİ
      çözücüsüyle zaten inşa edilmiştir — aşağıda dürüstlük beyanı);
  (C) DONDURULMUŞ KURAL: hangi ölçümden hangi öngörünün hangi formülle
      ve hangi ölüm eşiğiyle çıkacağı.

────────────────────────────────────────────────────────────────────────
DÜRÜSTLÜK BEYANI (ön-kayda aynen giriyor)
────────────────────────────────────────────────────────────────────────
1. **Erfc-ikiz YENİDEN İNŞA EDİLMEDİ.** KALEM "167_insa.main ile erfc
   zarflı gaz kur" diyor; ama o gaz ZATEN VAR: `HA4` (164_insa,
   `tip="sadakatli"`, `par={"tau_c":0.68,"delta":0.125,"c":-0.5}`),
   300 000 tekne, `maxF = 1.86e−09`, sıralılık TAM. 167_insa.main aynı
   çözücüyü çağırır ve DETERMİNİSTİKTİR: yeniden inşa bit-bit aynı
   diziyi verirdi. 10 dakikalık artıksız bir tekrar yerine mevcut gazın
   kapıları BAĞIMSIZ olarak yeniden denetlenir (aşağıda, z dosyasından).
2. **σ_ds ÖN-MÜHÜRÜ GEÇERSİZDİR ve öyle kaydedilir.** HA4'ün σ_ds²'si
   164'te YAYIMLANMIŞTIR (0.18861). Ölçülmüş bir sayı için "geniş bant
   ön-mühür" atmak sahte olurdu. KALEM'in bu maddesi bu koşuda
   UYGULANAMAZ — uzatmanın dört ıskası defterde kalır, beşincisi
   eklenmez.
3. **172'nin defteri OKUNMUŞTUR.** Ön-kayıt anında `172/G1.json` ve
   `G2.json`'un `HA4/K090/K070/E060` satırları (gE, gX, θ, KAL, μ̂²_E,
   π_E, ρ_E, σ'lar) BİLİNİYORDU; körlük iddiası YOKTUR. Ön-kayıt edilen
   şey KURALDIR ve öngörüler bu defterden TÜRETİLİR — sınanan şey,
   defterin (172, s_n sitesi, τ≤0.95, 8981 çizgi) 162 makinesinde
   (m_n sitesi, τ≤0.86, 3425 çizgi) ve 174d'nin DC-kaçağı defterinde
   TUTUP TUTMADIĞIDIR.
4. **ÖLÇÜLMEMİŞ olanlar (bu koşuda ilk kez ölçülecek):** kesim
   ailesinin (HA4, K090, K070, E060) hiçbirinde R_η / R_Ĉ / R_X,
   bant-bant defteri, üçüncü momentler (162 makinesi) YOK; hiçbirinde
   Σξ_q DC-kaçağı bant defteri (174d makinesi) YOK; κ'nın çizgi-tipine
   (asal ↔ kule) göre ayrışması hiçbir gazda YOK.
"""
import hashlib
import json
import math
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S175, S172, S167, S164, S155 = (SCR / "175", SCR / "172", SCR / "167",
                                SCR / "164", SCR / "155")
TWO_PI = 2 * np.pi
YOL = S175 / "ONKAYIT_K2.json"

TAU_C, DELTA = 0.68, 0.125
# kesim ailesi: ad -> ("keskin", τ_ust) | ("erfc", (τ_c, Δ))
KESIM = {"Hkeskin": ("keskin", 1.00), "K090": ("keskin", 0.90),
         "K070": ("keskin", 0.70), "HA4": ("erfc", (0.68, 0.125)),
         "E060": ("erfc", (0.60, 0.125))}
LAM7 = ["L050", "L060", "L070", "L085", "Hkeskin", "L115", "L130"]

print("=" * 74)
print("175a — ÖN-KAYIT (K2) + K1 İNŞA KAPILARI")
print("=" * 74)
if YOL.exists():
    print(f"  ** {YOL} ZATEN VAR — üzerine yazılmaz. **")
    d = json.load(open(YOL))
    print(f"  zaman={d['zaman']}  sha256={d['sha256'][:16]}")
    sys.exit(0)

zaman = time.strftime("%Y-%m-%d %H:%M:%S %z")
sha = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
print(f"  zaman  : {zaman}")
print(f"  sha256 : {sha}")

# ═════════════════════════════════════════════════════════════════════
# (A) ERFC ZARFININ TAM FORMÜLÜ — kaynak denetimi
# ═════════════════════════════════════════════════════════════════════
print("\n" + "=" * 74)
print("(A) ERFC ZARFI — 164_insa.merdiven kaynağıyla birebir denetim")
print("=" * 74)
sys.path.insert(0, str(QM / "164_configs"))
import importlib
I164 = importlib.import_module("164_insa")

d = np.load(QM / "128_odl_zeros6_2e6_zeros.npz")
Z = np.sort(np.asarray(d["zeros"], dtype=float))
zr = Z[len(Z) - 300000:]
t0, t1 = float(zr[0]), float(zr[-1])
L_hedef = float(np.log(0.5 * (t0 + t1) / TWO_PI))
om, a_pen, a_ham, w164 = I164.merdiven(L_hedef, TAU_C, DELTA, tau_ust=1.00)
tau = om / L_hedef
w_lit = np.array([0.5 * math.erfc((t - TAU_C) / DELTA) for t in tau])
fark = float(np.max(np.abs(w164 - w_lit)))
print(f"  L_hedef = {L_hedef:.12f}   nline = {len(om)}")
print("  FORMÜL:  w_q = ½·erfc( (τ_q − τ_c)/Δ ) ,  τ_c = %.2f , Δ = %.3f"
      % (TAU_C, DELTA))
print("           τ_q = ω_q/L_hedef ,  ω_q = log q ,  q = p^k (asal kuvvet)")
print("           a_q = 1/(π·k·√q)  ,  A_q = λ·a_q·w_q  (λ = 1.00)")
print("           S(t) = −Σ_q A_q sin(ω_q t) ,  N̄(z)+S(z) = n − ½")
print(f"  164_insa.merdiven ile maks |w_kod − w_formül| = {fark:.2e}   "
      f"{'✓' if fark == 0.0 else '✗'}")
print("  152 raporu (satır 37): 'A4 erfc 0.68/0.125 | a_q → a_q·½erfc("
      "(τ_q−0.68)/0.125) — S3'ü oturtan kesim'  ⇒ BİREBİR AYNI.")
print("  w(τ) örnekleri: " + "  ".join(
    "τ=%.2f→%.4f" % (t, 0.5 * math.erfc((t - TAU_C) / DELTA))
    for t in (0.50, 0.60, 0.68, 0.75, 0.80, 0.90)))

# kesim ekseni koordinatı (yalnız İNŞA'dan, ölçüm değil)
def wgt(ad):
    tip, p = KESIM[ad]
    if tip == "keskin":
        return (tau <= p + 1e-12).astype(float)
    return np.array([0.5 * math.erfc((t - p[0]) / p[1]) for t in tau])


KES = {}
print("\n  KESİM EKSENİ (inşa tarafı; τ̄_A := Σ A²τ / Σ A²):")
for ad in KESIM:
    A = a_ham * wgt(ad)
    KES[ad] = dict(tau_bar=float((A ** 2 * tau).sum() / (A ** 2).sum()),
                   rmsS=float(np.sqrt(0.5 * np.sum(A ** 2))),
                   P_hi=float((A[tau > 0.5] ** 2).sum() / (A ** 2).sum()),
                   tip=KESIM[ad][0], par=KESIM[ad][1])
    print("    %-8s %-7s τ̄_A = %.5f   rms S = %.5f   P(τ>0.5)/P = %.5f"
          % (ad, KESIM[ad][0], KES[ad]["tau_bar"], KES[ad]["rmsS"],
             KES[ad]["P_hi"]))

# ═════════════════════════════════════════════════════════════════════
# (B) K1 — İNŞA KAPILARI (mevcut erfc-ikiz HA4'ün bağımsız denetimi)
# ═════════════════════════════════════════════════════════════════════
print("\n" + "=" * 74)
print("(B) K1 — ERFC-İKİZ (HA4) İNŞA KAPILARI")
print("=" * 74)
ins = json.load(open(S164 / "insa_HA4.json"))
z = np.load(S164 / "z_HA4.npy")
z155 = np.load(S155 / "z_HA4.npy")
ayni = bool(np.array_equal(z, z155))
dz = np.diff(z)
mid = 0.5 * (z[:-1] + z[1:])
Lw = np.log(mid / TWO_PI)
ds = dz * Lw / TWO_PI - 1
L = float(Lw.mean())
K1G = dict(
    n=int(len(z)), n_hedef=300000, hucre_benzersiz=int(ins["hucre_benzersiz"]),
    maxF=float(ins["maxF"]), esik_maxF=1e-8, nF_asan=int(ins["nF_asan"]),
    sirali=bool(np.all(dz > 0)), min_dz=float(dz.min()),
    sigma_ds=float(np.std(ds)), sigma_ds2=float(np.var(ds)),
    sigma_ds2_164=float(ins["sigma_ds2"]), L=L, L_164=float(ins["L"]),
    z155_ayni=ayni, tip=ins["tip"], par=ins["par"], h=float(ins["h"]))
print("  ilk-kök hücre : %d / %d   %s"
      % (K1G["hucre_benzersiz"], 300000,
         "✓" if K1G["hucre_benzersiz"] == 300000 else "✗"))
print("  n             : %d   %s" % (K1G["n"], "✓" if K1G["n"] == 300000 else "✗"))
print("  maks|F|       : %.3e  (eşik 1e−8)  aşan tekne %d   %s"
      % (K1G["maxF"], K1G["nF_asan"], "✓" if K1G["maxF"] <= 1e-8 else "✗"))
print("  sıralılık     : %s   min Δz = %.6f   %s"
      % ("TAM" if K1G["sirali"] else "BOZUK", K1G["min_dz"],
         "✓" if K1G["sirali"] else "✗"))
print("  σ_ds (z'den)  : %.5f  (σ_ds² = %.5f ; 164 yayımı %.5f, fark %.1e)"
      % (K1G["sigma_ds"], K1G["sigma_ds2"], K1G["sigma_ds2_164"],
         abs(K1G["sigma_ds2"] - K1G["sigma_ds2_164"])))
print("  L (z'den)     : %.9f   (164: %.9f)" % (K1G["L"], K1G["L_164"]))
print("  155/z_HA4 ≡ 164/z_HA4 : %s" % ("✓ bit-bit aynı" if ayni else "✗ FARKLI"))
print("  çözücü        : %s   par=%s   h=%.4f" % (K1G["tip"], K1G["par"],
                                                  K1G["h"]))
# sağlık: R_bant (167_olcum'un C_HA4.json'undan; 169'un değişmez filtresi)
C = json.load(open(S167 / "C_HA4.json"))
bnt = [b for b in C["bant"] if b.get("olculdu")]
fil = [b for b in bnt if 0.52 <= b["lo"] <= 0.68 and b["SNR"] >= 3
       and b["tau_eff"] < 0.85]
K1G["R_bant"] = [float(b["R_bant"]) for b in fil]
K1G["R_bant_lo"] = [float(b["lo"]) for b in fil]
K1G["R_bant_min"] = float(min(K1G["R_bant"])) if fil else float("nan")
print("  SAĞLIK R_bant (169 filtresi, lo∈[0.52,0.68]): %s  ⇒ min %.4f "
      "(eşik 0.98)  %s"
      % (" ".join("%.4f" % r for r in K1G["R_bant"]), K1G["R_bant_min"],
         "✓" if K1G["R_bant_min"] >= 0.98 else "✗"))
K1G["gecti"] = bool(K1G["hucre_benzersiz"] == 300000 and K1G["n"] == 300000
                    and K1G["maxF"] <= 1e-8 and K1G["sirali"]
                    and K1G["R_bant_min"] >= 0.98)
print("  ⇒ K1 KAPILARI: %s" % ("**GEÇTİ**" if K1G["gecti"] else "**KALDI**"))

# ═════════════════════════════════════════════════════════════════════
# (C) DONDURULMUŞ KURAL — öngörüler ve ölüm eşikleri
# ═════════════════════════════════════════════════════════════════════
print("\n" + "=" * 74)
print("(C) DONDURULMUŞ KURAL")
print("=" * 74)
G1 = json.load(open(S172 / "G1.json"))
G2 = json.load(open(S172 / "G2.json"))


def V_O(g):                      # π := K/V_O  ⇒  V_O = K/π (özdeşlik)
    return G1[g]["E"]["K"] / G1[g]["E"]["pi"]


def ortE(g):                     # μ̂² = ort(E)²/V_O  ⇒  |ort E| = √(μ̂² V_O)
    return math.sqrt(G1[g]["E"]["mu2"] * V_O(g))


# --- P1: 174'ün K2-F köprüsü kesim gazlarında da tutar mı? ------------
# 174 K2-F: R_η(162 makinesi)/π_E(172 defteri) − 1 , 8 gazda −0.86…−1.83%
K2F = np.array([-1.19, -1.46, -0.86, -0.99, -1.12, -1.30, -1.58, -1.83])
rbar = float(np.mean(1 + K2F / 100))
rsd = float(np.std(1 + K2F / 100))
P1 = {g: dict(piE=G1[g]["E"]["pi"], R_eta_ongoru=G1[g]["E"]["pi"] * rbar)
      for g in KESIM}
print("\n  P1 — R_η ÖNGÖRÜSÜ (K2-F köprüsü: R_η = π_E × r̄, "
      "r̄ = %.5f ± %.5f, 174'ün 8 gazı)" % (rbar, rsd))
for g in KESIM:
    print("      %-8s π_E = %.5f  ⇒  R_η^ön = %.5f  (±%%1.5 bandı "
          "[%.4f, %.4f])" % (g, P1[g]["piE"], P1[g]["R_eta_ongoru"],
                             P1[g]["R_eta_ongoru"] * 0.985,
                             P1[g]["R_eta_ongoru"] * 1.015))
print("      ÖLÇÜT: dört kesim gazının hepsinde |R_η/π_E − 1| ≤ %2  ⇒ ✓;"
      " biri bile aşarsa köprü kesim ekseninde KIRIK.")

# --- P2: gerçek gaz KESİM AİLESİNİN İÇİNDE mi? -----------------------
Rson = 1.28829                        # 174b'nin ölçtüğü (yayımlı)
alt = min(P1[g]["R_eta_ongoru"] for g in KESIM)
ust = max(P1[g]["R_eta_ongoru"] for g in KESIM)
print("\n  P2 — AİLE-DIŞI FAZLANIN ADRESİ  (174: R_η(son)=1.2883, λ "
      "ailesinin tavanı 1.2745'in %1.08 ÜSTÜNDE)")
print("      öngörülen kesim ailesi aralığı: [%.4f, %.4f]  ∋ 1.28829 ? %s"
      % (alt, ust, "EVET (öngörü)" if alt <= Rson <= ust else "HAYIR"))
print("      ÖLÇÜT: ölçülen kesim ailesi R_η(son)'u KAPSIYORSA ⇒ "
      "**fazlanın ekseni KESİMDİR**; kapsamıyorsa H-K1'in η-kanalı ölür.")

# --- P3: DC kaçağı ort(E) — özdeşlikten öngörü -----------------------
print("\n  P3 — DC KAÇAĞI ort(E) ÖNGÖRÜSÜ (özdeşlik |ort E| = √(μ̂²_E·V_O),"
      " V_O = K_E/π_E; 174d'nin Ö2'si bunu 1.3e−13'te doğrulamıştı)")
P3 = {}
for g in ["son", "Hkeskin"] + [x for x in KESIM if x != "Hkeskin"]:
    P3[g] = dict(V_O=V_O(g), mu2=G1[g]["E"]["mu2"], ortE=ortE(g))
    print("      %-8s V_O = %.6f   μ̂²_E = %.6f   ⇒ ort(E)^ön = %+.6f"
          % (g, P3[g]["V_O"], P3[g]["mu2"], P3[g]["ortE"]))
print("      ÖLÇÜT (Ö-P3): 174d makinesinin Σ_q ξ_q'su bu değeri ‰5 içinde"
      " vermeli (işaret dahil, hepsi POZİTİF öngörülüyor).")

# --- P4: MÜHÜR KURALI (KALEM'in literal kuralı) ----------------------
# 174 K3b (yayımlı): Σξ bant bant, gerçek ↔ Hkeskin
X174 = {"son":     {(0.70, 0.80): +1.5384e-2, (0.80, 0.95): +9.7396e-3},
        "Hkeskin": {(0.70, 0.80): +2.0628e-2, (0.80, 0.95): +1.7196e-2}}
acik_keskin = sum(X174["son"].values()) - sum(X174["Hkeskin"].values())
print("\n  P4 — MÜHÜR KURALI (KALEM, literal)")
print("      174: τ>0.70 açığı (gerçek − keskin-ikiz) = %+.6e"
      % acik_keskin)
print("      KURAL: |Σξ(son) − Σξ(erfc-ikiz)|(τ>0.70) ≤ %.6e "
      "(= %%50) ⇒ **H-K1 YAŞAR**" % (0.5 * abs(acik_keskin)))
print("      BEKLENTİ (ön-kayıtlı, P3'ten): erfc-ikiz ALTTAN ıskalayacak —"
      " işaret DÖNECEK (açık +), büyüklük 0.8–1.6× keskin açığı.")
print("      Yani literal kuralın ÖLECEĞİNİ önceden yazıyoruz; ölüm "
      "bilgilendiricidir: gerçek gaz iki ikizin ARASINDADIR (164'ün "
      "sandviçi), tek bir kesim NOKTASI değil.")

# --- P5: PARAMETRESİZ AYRIŞIM — kesim kesri f ------------------------
print("\n  P5 — PARAMETRESİZ AYRIŞIM: KESİM KESRİ f")
print("      Her gözlenebilir y için  f(y) := [y(son) − y(Hkeskin)] /"
      " [y(HA4) − y(Hkeskin)]")
print("      (koordinat YOK, ara değer YOK, serbest parametre YOK.)")
print("      * Kesim ekseni gerçeği TEK BAŞINA taşıyorsa bütün f'ler AYNI")
print("        çıkmalıdır (tek-parametreli aile cümlesi).")
print("      * ÖLÇÜT (H-K1 GÜÇLÜ): DC-kaçağı bantlarının f_b'leri ±0.15")
print("        içinde sabit ⇒ açık SAF KESİMDİR.")
print("      * f'lerin yayılımı = **KİLİT PAYI** (kesimin taşımadığı).")
print("      * f < 0 ⇒ kesim ekseni o nicelikte TERS yönde ⇒ o pay kesime")
print("        yazılamaz (negatif kesim payı defterde ayrı satır).")
P5 = {}
for ad, f in (("π_E", lambda g: G1[g]["E"]["pi"]),
              ("ort(E)", ortE), ("g_E", lambda g: G1[g]["gE"]),
              ("g_X", lambda g: G1[g]["gX"]),
              ("θ(G2)", lambda g: G2["theta"][g]),
              ("θ(G1)", lambda g: G1[g]["th"]),
              ("M = KAL(G1)", lambda g: G1[g]["KAL"]),
              ("ρ_E", lambda g: G1[g]["E"]["rho"]),
              ("σ_ds", lambda g: G1[g]["sigds"])):
    ys, yh, ye = f("son"), f("Hkeskin"), f("HA4")
    P5[ad] = dict(son=ys, Hkeskin=yh, HA4=ye,
                  f=(ys - yh) / (ye - yh) if ye != yh else float("nan"))
    print("      f(%-11s) = %+7.4f    [son %.5f | keskin %.5f | erfc %.5f]"
          % (ad, P5[ad]["f"], ys, yh, ye))
print("      (bu satırlar 172 defterinden TÜRETİLDİ — okunmuş sayılar;")
print("       sınav, 162 makinesinin R_η'sı ve 174d'nin Σξ'si aynı f'yi")
print("       verecek mi? P1/P3 bunun ön-kaydıdır.)")

# --- P6: ΔM DEFTERİ --------------------------------------------------
dlogM = math.log(G1["son"]["KAL"] / G1["Hkeskin"]["KAL"])
dgE = math.log(G1["son"]["gE"] / G1["Hkeskin"]["gE"])
dgX = math.log(G1["son"]["gX"] / G1["Hkeskin"]["gX"])
dth = math.log(G1["son"]["th"] / G1["Hkeskin"]["th"])
print("\n  P6 — ΔM DEFTERİ (özdeşlik: M ≡ g_E·g_X²·θ ; çapa = KESKİN İKİZ)")
print("      Δlog M = %+.5f = Δlog g_E %+.5f + 2Δlog g_X %+.5f + Δlog θ "
      "%+.5f  (kalıntı %.1e)"
      % (dlogM, dgE, 2 * dgX, dth, abs(dlogM - dgE - 2 * dgX - dth)))
print("      AYRIŞIM (dondurulmuş formül, f* := f(R_η) ÖLÇÜLEN kesim kesri):")
print("        kesim payı := f*·[Δlog M](keskin→erfc)")
print("        kilit payı := [Δlog g_E + 2Δlog g_X](son) − f*·(aynısı)(erfc)")
print("        θ payı     := Δlog θ(son) − f*·Δlog θ(erfc)")
print("        kesim + kilit + θ ≡ Δlog M   (ÖZDEŞ, kalıntı ≤ 1e−12)")

# --- P7: K3 ÇARPIMSAL κ ---------------------------------------------
print("\n  P7 — K3: ÇARPIMSAL κ (H-K3)")
print("      TÜRETİM (ölçümden önce). s_n = m̄_n + δ(m̄_n),")
print("        δ(m̄) = ḡ Σ_c A_c cos(πτ_c) sin(ω_c m̄)   [143 G-yasası köşesi]")
print("      κ(ν) = ⟨e^{iνm̄}e^{iνδ}⟩ ; ⟨e^{iνm̄}cos(Ωm̄)⟩ = ½·1{Ω=±ν}:")
print("        1. mertebe: κ⁽¹⁾ = −π A_ν τ_ν cos(πτ_ν)        [174e'nin yasası]")
print("        2. mertebe: κ⁽²⁾ = (π²τ_ν²/2)·(S_ν − D_ν) ,")
print("           S_ν = Σ_{q₁q₂=Q, ikisi de merdivende} A₁A₂cos(πτ₁)cos(πτ₂)")
print("           D_ν = Σ_{q₁/q₂=Q, ikisi de merdivende} 2A₁A₂cos(πτ₁)cos(πτ₂)")
print("      MERDİVEN = ASAL KUVVETLER (q = p^k) ⇒ q₁q₂ = Q'nun merdiven")
print("      içi çözümü YALNIZ Q = p^k (k≥2) için vardır — **KULE**.")
print("      Q = p² için kapalı oran (sıfır parametre):")
print("        |κ⁽²⁾/κ⁽¹⁾| = 2τ_p cos²(πτ_p) / (|cos 2πτ_p| · W_pos(2τ_p))")
for tp in (0.25, 0.30, 0.35, 0.40, 0.45):
    W = math.exp(-2 * math.pi ** 2 * (2 * tp) ** 2 * 0.2730 ** 2)
    den = abs(math.cos(2 * math.pi * tp)) * W
    print("        τ_p=%.2f (τ_Q=%.2f): %s"
          % (tp, 2 * tp, ("%.3f" % (2 * tp * math.cos(math.pi * tp) ** 2 / den))
             if den > 1e-6 else "∞ (cos 2πτ_p = 0)"))
print("      ÖN-KAYITLI ÖLÇÜTLER:")
print("        K3-a: |κ|/|κ̂| oranı ÇİZGİ TİPİNE göre ayrılır. H-K3 doğruysa")
print("              medyan(kule) ≫ medyan(asal) ve medyan(asal) ≈ 1.")
print("              **medyan(asal) ≥ 1.5 (τ>0.60) ⇒ H-K3 (tam-rezonans")
print("              çarpımsal) ÖLÜR** — fazla asallarda da var demektir.")
print("        K3-b: kule çizgilerinde κ⁽²⁾ eklenince oran 1'e YAKLAŞMALI;")
print("              yaklaşmazsa (|1−oran| küçülmüyorsa) H-K3 ölür.")
print("        K3-c (RAKİP, ön-kayıtlı): fazla W_pos'un AŞIRI SÖNÜMÜ olabilir.")
print("              Sınav: κ̂₀ := −πA_ντ_ν cos(πτ_ν)  (W_pos YOK).")
print("              **τ>0.55'te medyan |κ|/|κ̂₀| ∈ [0.8,1.25] ⇒ rakip KAZANIR**")
print("              ve H-K3 gereksizdir.")

# --- P8: K4 θ DEFTERİ ------------------------------------------------
print("\n  P8 — K4: θ DEFTERİ")
print("      ÖN-KAYITLI BEKLENTİ: f(θ) < 0 (kesim ekseni θ'yı TERS yöne")
print("      taşır) ⇒ ΔM'nin %51.5'lik θ payı kesim şekliyle KAPANMAZ,")
print("      **büyür**. Ölçüt: f(θ) < 0 ⇒ H-K4 ÖLÜR (ve ölüm keskindir:")
print("      'θ'nın taşıyıcısı kesim değildir').")

ONK = dict(
    zaman=zaman, sha256=sha, betik="175a_onkayit.py",
    durustluk=[
        "erfc-ikiz YENİDEN İNŞA EDİLMEDİ: HA4 (164, sadakatli) kullanıldı; "
        "çözücü deterministik, yeniden inşa bit-bit aynı olurdu.",
        "σ_ds ön-mühürü GEÇERSİZ: HA4'ün σ_ds²'si 164'te yayımlanmıştır.",
        "172/G1+G2 defterinin HA4/K090/K070/E060 satırları ön-kayıt anında "
        "OKUNMUŞTU; körlük iddiası yoktur.",
        "ÖLÇÜLMEMİŞ: kesim ailesinde R_η/R_Ĉ/R_X, bant defteri, üçüncü "
        "momentler, Σξ DC-kaçağı defteri, κ'nın asal↔kule ayrışması."],
    zarf=dict(formul="w_q = 0.5*erfc((tau_q - tau_c)/delta)",
              tau_c=TAU_C, delta=DELTA, lam=1.00,
              a_q="1/(pi*k*sqrt(q)), q=p^k", A_q="lam*a_q*w_q",
              S="S(t) = -sum_q A_q sin(omega_q t)", seviye="N̄+S = n-1/2",
              L_hedef=L_hedef, nline=int(len(om)),
              kod_farki=fark),
    kesim_ekseni=KES, K1=K1G,
    kural=dict(
        P1=dict(ad="R_η = π_E × r̄ (K2-F köprüsü)", rbar=rbar, rsd=rsd,
                esik_pct=2.0, ongoru={g: P1[g]["R_eta_ongoru"] for g in P1}),
        P2=dict(ad="gerçek gaz kesim ailesinin İÇİNDE mi?", R_son=Rson,
                aralik_ongoru=[alt, ust]),
        P3=dict(ad="ort(E) özdeşlik öngörüsü", esik_permil=5.0,
                ongoru={g: P3[g]["ortE"] for g in P3}),
        P4=dict(ad="MÜHÜR: τ>0.70 açığı ≥%50 küçülür mü?",
                acik_keskin=acik_keskin, esik=0.5 * abs(acik_keskin),
                beklenti="ölecek; işaret dönecek (sandviç)"),
        P5=dict(ad="kesim kesri f(y)", tanim="[y(son)-y(Hk)]/[y(HA4)-y(Hk)]",
                sabitlik_esigi=0.15, defter={k: v for k, v in P5.items()}),
        P6=dict(ad="ΔM ayrışımı", dlogM=dlogM, dgE=dgE, dgX2=2 * dgX,
                dth=dth,
                formul="kesim=f*·ΔlogM(erfc); kilit=(E,X artığı); θ=θ artığı"),
        P7=dict(ad="çarpımsal κ", merdiven="asal kuvvetler ⇒ yalnız KULE",
                K3a_esik_medyan_asal=1.5, K3c_bant=[0.8, 1.25]),
        P8=dict(ad="θ defteri", beklenti="f(θ) < 0 ⇒ H-K4 ölür")),
)
S175.mkdir(parents=True, exist_ok=True)
YOL.write_text(json.dumps(ONK, indent=1, ensure_ascii=False))
print("\n-> %s   (sha256 %s)" % (YOL, sha[:16]))
