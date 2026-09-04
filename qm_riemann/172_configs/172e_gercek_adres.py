# -*- coding: utf-8 -*-
"""
172e — G3: GERÇEK GAZIN DAR ADRESİ (4 Eylül 2026)
Veri: 172/G1.json + 172/G2.json. YENİ KOŞU YOK.

══════════════════════════════════════════════════════════════════════
YÖNTEM (uyumdan ÖNCE)
══════════════════════════════════════════════════════════════════════
G1 ve G2, λ ekseninin BÜTÜN çarpanlarını ölçülür eğrilere bağladı:
Q_E(λ), ρ_E(λ), Q_X(λ), ρ_X(λ), g_E(λ), g_X(λ), θ(λ), α(λ).
"Gerçek gaz bir λ-gazı mıdır?" sorusu artık TEK BİR sınavdır:

  her büyüklük X için  λ_eş(X) := X(son)'u veren λ  (7-nokta merdiveninde
  tekdüze kolda ters çevirme). Gerçek gaz bir λ-gazıysa BÜTÜN λ_eş'ler
  aynı olmalıdır. **λ_eş'lerin YAYILIMI = gerçeğin fazlasının adresi.**

Artıklar ayrıca bir ÇAPA λ'sında okunur: σ_X̃'nin λ_eş'i (171 §T2b.2'nin
konvansiyonu, orada 0.9335).

══════════════════════════════════════════════════════════════════════
ÖN-MÜHÜR (koşudan önce; 4 Eyl, 172e)
══════════════════════════════════════════════════════════════════════
 Ö1  172a'nın beklentisi AYNEN sınanır: σ-çapasında ölçülen artıklar
     g_E ≈ 0 (|·| ≤ %0.5), (g_X)² ≈ 0 (≤ %0.5), θ ≈ **+%3** (2-4 aralığı).
 Ö2  λ_eş TUTARLILIĞI: gerçek gaz bir λ-gazıysa bütün λ_eş'ler ±0.05
     içinde uyuşur. (172c Ö7 zaten Q_E ↔ Q_X için 0.162 fark bulup
     bunu öldürdü; burada BÜTÜN büyüklüklere genişletiliyor.)
 Ö3  **162 KÖPRÜSÜ — YÖN ÖN-KAYITLI.** 162 §6: yıkıcı girişim η
     kanalındadır (Σ|c_η|²/2 ÷ Var(η) = 1.288, %29 sessizleşme),
     **Ĉ kanalında YOKTUR** (1.026). 165'in E alanı η kanalıdır,
     X alanı X̃ ≈ ΔĈ kanalıdır. Bellek fazlardaysa **E-kanalı
     büyüklüklerinin λ_eş'i X-kanalınınkinden KÜÇÜK** çıkmalıdır
     (η daha sessiz ⇒ daha küçük λ gibi görünür):
         **λ_eş(E-kanalı) < λ_eş(X-kanalı) ve < λ_eş(marjinaller)**
     Bu tek yönlü bir öngörüdür; ters çıkarsa 162 köprüsü ölür.
 Ö4  ŞEKİL BORCU YOK: 172d'nin τ-eğimi yasası son'da +%3.0 tuttu
     (ön-kayıtlı %20 eşiğinin çok içinde) ⇒ gerçek gazın fazlası
     ŞEKİLDE değil SEVİYEDE olmalı; burada yalnız alıntılanır.
 Ö5  ΔM'nin ÇARPAN DEFTERİ kapanmalı: log-artıkların toplamı
     ölçülen ΔlogM(son) ile ‰1 içinde uyuşmalı (özdeşlik denetimi).

══════════════════════════════════════════════════════════════════════
SONUÇ (yalnız gerçek koşudan, `172/log_172e.txt`)
══════════════════════════════════════════════════════════════════════
 Ö1 ✗/✓ σ-çapasında (λ_eş(σ_X̃) = 0.9363): g_E artığı **+%2.05**
      (ön-kayıt ≈0, |·|≤%0.5 ⇒ **ÖLDÜ**); (g_X)² **+%0.32** ✓;
      θ **+%2.53** (ön-kayıt +%3, 2-4 aralığı ⇒ ✓ ucundan).
      172a'nın "g_E artığı ≈ 0" beklentisi g_E–rE DOĞRUSU için doğruydu
      (+%0.25) ama λ-ÇAPASI için yanlış: gerçek gaz doğrunun ÜSTÜNDE
      ama doğru boyunca λ ≈ 0.76'ya KAYMIŞ durumda.
 Ö2 ✗ **KESİN ÖLÜM.** λ_eş yayılımı 0.69 … 0.97 (ön-kayıt ±0.05).
      **Gerçek gaz bir λ-gazı DEĞİLDİR.**
 Ö3 ✓ **162 KÖPRÜSÜ AYAKTA (yön ön-kayıtlıydı).**
      λ_eş(E kanalı, η) = **0.7768**  <  λ_eş(X kanalı, ΔĈ) = **0.9089**
      <  λ_eş(marjinaller) = **0.9354**.
      E kanalının beş büyüklüğü de (rE 0.742, μ̂²_E 0.749, g_E 0.758,
      Q_E 0.785, ρ_E 0.823) X kanalının hepsinden küçük. Bu, 162 §6'nın
      "yıkıcı girişim η kanalında (1.288), Ĉ kanalında YOK (1.026)"
      ölçümünün çarpan dilindeki KARŞILIĞIDIR.
 Ö4 ✓ 172d: son'un τ-eğimi yasayı +%3.0 ile sağlıyor ⇒ şekil borcu yok.
 Ö5 ✓ Çarpan defteri **3.3e−06** log farkla kapandı: g_E +%2.05,
      (g_X)² +%0.32, θ +%2.53  ⇒ toplam **+%4.97** = ölçülen ΔlogM
      (M(son) = 1.0561, λ-eğrisi 1.0061).
 DAR ADRES (σ-çapasında artık sıralaması): **μ̂²_E −%23.8** (model
      alanının DC kaçağı = sıfır tarağının asal rezonansları),
      rE −%6.5, α +%4.6, **Q_E −%3.6**, θ +%2.5, g_E +%2.1;
      X kanalı ise NEREDEYSE TAM SIFIR (ρ_X +%0.15, g_X +%0.16,
      Q_X −%0.23). Cebirsel sürücü Q_E'dir (ΔrE = −Δρ_E + 2ΔQ_E − Δμ̂²_E
      = +0.0155 − 0.0652 + 0.0127 = −0.0371 ✓).
"""
import json
import numpy as np

SCR = ("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
       "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/172/")
G1 = json.load(open(SCR + "G1.json"))
G2 = json.load(open(SCR + "G2.json"))
LAM7 = ["L050", "L060", "L070", "L085", "Hkeskin", "L115", "L130"]
LV = np.array([.50, .60, .70, .85, 1.00, 1.15, 1.30])

BUY = {   # ad -> (kanal, gazdan değer)
    "σ_ds": ("marjinal", lambda g: G1[g]["sigds"]),
    "σ_X̃": ("marjinal", lambda g: G1[g]["sigX"]),
    "σ_Ĉ": ("marjinal", lambda g: G1[g]["sigC"]),
    "rE": ("E (η)", lambda g: G1[g]["E"]["r"]),
    "Q_E": ("E (η)", lambda g: G1[g]["E"]["Q"]),
    "ρ_E": ("E (η)", lambda g: G1[g]["E"]["rho"]),
    "μ̂²_E": ("E (η)", lambda g: G1[g]["E"]["mu2"]),
    "g_E": ("E (η)", lambda g: G1[g]["gE"]),
    "rX": ("X (ΔĈ)", lambda g: G1[g]["X"]["r"]),
    "Q_X": ("X (ΔĈ)", lambda g: G1[g]["X"]["Q"]),
    "ρ_X": ("X (ΔĈ)", lambda g: G1[g]["X"]["rho"]),
    "g_X": ("X (ΔĈ)", lambda g: G1[g]["gX"]),
    "θ": ("θ", lambda g: G2["theta"][g]),
    "α (A-genişliği)": ("şekil", lambda g: G2["alfa"][g]),
}


def egri(f):
    return np.array([f(g) for g in LAM7])


def lam_es(y, v):
    """y(λ) eğrisinin v'yi verdiği λ (tekdüze kolda; kol otomatik)."""
    # tekdüze en uzun kolu seç: λ ≤ 1.15 (indeks 0..5) tercih
    for kol in (slice(0, 6), slice(0, 5), slice(None)):
        yy, ll = y[kol], LV[kol]
        d = np.diff(yy)
        if np.all(d > 0):
            if yy[0] <= v <= yy[-1]:
                return float(np.interp(v, yy, ll)), True
        elif np.all(d < 0):
            if yy[-1] <= v <= yy[0]:
                return float(np.interp(v, yy[::-1], ll[::-1])), True
    return float("nan"), False


print("== Ö2/Ö3: HER BÜYÜKLÜĞÜN KENDİ λ_eş'i ==")
print("%-16s %-10s %10s %10s %10s   %s"
      % ("büyüklük", "kanal", "son", "λ=1.00", "**λ_eş**", "tekdüze?"))
LE = {}
for ad, (kan, f) in BUY.items():
    y = egri(f)
    v = f("son")
    le, ok = lam_es(y, v)
    LE[ad] = (kan, le)
    print("%-16s %-10s %10.5f %10.5f %10s   %s"
          % (ad, kan, v, f("Hkeskin"),
             ("%.4f" % le) if le == le else "  —  ",
             "evet" if ok else "HAYIR (aralık dışı/tekdüze değil)"))

for kan in ("marjinal", "E (η)", "X (ΔĈ)", "θ", "şekil"):
    v = [LE[a][1] for a in LE if LE[a][0] == kan and LE[a][1] == LE[a][1]]
    if v:
        print("  %-10s λ_eş ortalaması = %.4f   [%.4f, %.4f]  (%d büyüklük)"
              % (kan, np.mean(v), min(v), max(v), len(v)))

mar = np.mean([LE[a][1] for a in ("σ_ds", "σ_X̃", "σ_Ĉ")])
eK = [LE[a][1] for a in ("rE", "Q_E", "ρ_E", "g_E") if LE[a][1] == LE[a][1]]
xK = [LE[a][1] for a in ("rX", "Q_X", "ρ_X") if LE[a][1] == LE[a][1]]
print("  **Ö3 SINAVI:** λ_eş(E) = %.4f  <  λ_eş(X) = %.4f ?  %s"
      % (np.mean(eK), np.mean(xK),
         "EVET (162 köprüsü ayakta)" if np.mean(eK) < np.mean(xK) else "HAYIR"))
print("  marjinaller λ_eş = %.4f (171 §T2b.2: 0.9335 σ_X̃ ile)" % mar)

print()
print("== Ö1/Ö5: σ-ÇAPASINDA (λ_eş = %.4f) ÇARPAN ARTIKLARI ==" % mar)
capa = LE["σ_X̃"][1]
print("  çapa = λ_eş(σ_X̃) = %.4f" % capa)
top = 0.0
print("%-16s %10s %10s %10s   %s"
      % ("çarpan", "ölçülen", "λ-eğrisi", "artık %", "ön-kayıt"))
ONK = {"g_E": "≈0 (|·|≤%0.5)", "(g_X)²": "≈0 (≤%0.5)", "θ": "+%3 (2-4)"}
for ad, f, us in (("g_E", lambda g: G1[g]["gE"], 1),
                  ("(g_X)²", lambda g: G1[g]["gX"], 2),
                  ("θ", lambda g: G2["theta"][g], 1)):
    y = egri(f)
    beklenen = float(np.interp(capa, LV, y))
    art = us * np.log(f("son") / beklenen)
    top += art
    print("%-16s %10.5f %10.5f %+10.2f   %s"
          % (ad, f("son") ** us, beklenen ** us, 100 * (np.exp(art) - 1),
             ONK.get(ad, "")))
Mson = G2["KAL"]["son"] / G2["KAL"]["Hkeskin"]
Mcapa = float(np.interp(capa, LV, np.array(
    [G2["KAL"][g] / G2["KAL"]["Hkeskin"] for g in LAM7])))
print("  ÇARPAN TOPLAMI  = %+.5f (log)   ⇒ %+.2f%%" % (top, 100 * (np.exp(top) - 1)))
print("  ÖLÇÜLEN ΔlogM   = %+.5f (log)   ⇒ %+.2f%%   (M(son)=%.4f, "
      "λ-eğrisi %.4f)" % (np.log(Mson / Mcapa), 100 * (Mson / Mcapa - 1),
                          Mson, Mcapa))
print("  Ö5 defter farkı = %.2e" % abs(top - np.log(Mson / Mcapa)))

print()
print("== DAR ADRES: ARTIKLARIN SIRALAMASI (σ-çapasında) ==")
sir = []
for ad, (kan, _le) in LE.items():
    f = BUY[ad][1]
    y = egri(f)
    bek = float(np.interp(capa, LV, y))
    sir.append((100 * (f("son") / bek - 1), ad, kan, f("son"), bek))
for a, ad, kan, v, b in sorted(sir, key=lambda z: -abs(z[0])):
    print("  %-16s %-10s ölç %10.5f  eğri %10.5f  artık %+7.2f%%"
          % (ad, kan, v, b, a))

json.dump(dict(capa=capa, lam_es={a: LE[a][1] for a in LE},
               kanal={a: LE[a][0] for a in LE},
               artik={ad: a for a, ad, _k, _v, _b in sir}),
          open(SCR + "G3.json", "w"), indent=1)
print("\n-> %sG3.json" % SCR)
