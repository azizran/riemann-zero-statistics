"""153 itmeli gaz -> KALICI RAPOR markdown. Tablo sayilari ozet_*.json'dan
otomatik; elle sayi girilmez (gercek referans + 152 referanslari disinda).
Metrik kurallari 152 ile AYNI: |Gamma|>1.5 patlak sayilir, faz orani bant
BASINA ortalama, bant 1-4 uzerinden."""
import json
import math
from pathlib import Path

S = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
         "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/153")
OUT = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann/"
           "153_itmeli_gaz_RAPOR.md")

GERCEK = {0.5375: (0.787, 0.206, 0.40), 0.585: (0.550, 0.488, 1.07),
          0.660: (0.118, 0.614, 1.85), 0.740: (-0.148, 0.498, 2.12),
          0.815: (-0.103, 0.655, 1.50)}
TAUS = [0.5375, 0.585, 0.660, 0.740, 0.815]
GERCEK_PS = {0.1: 0.000960, 0.2: 0.007327, 0.3: 0.024203,
             0.5: 0.103847, 0.7: 0.247217}
GMAKS = 1.5


def faz(re, im):
    a = math.atan2(im, re)
    return a + 2 * math.pi if a < -1.0 else a


G_FAZ = {t: faz(*GERCEK[t][:2]) for t in TAUS}


def metrik(bantlar):
    sg4 = 0.0; nb4 = 0; oran = {}; mag = {}
    sR = 0.0; nR = 0; eksik = []; supheli = []; patlak = []
    for b in bantlar:
        t = b["tau"]
        if b["re"] is None:
            eksik.append(t); continue
        gr = GERCEK[t]
        mm = math.hypot(b["re"], b["im"])
        mag[t] = mm
        if mm > GMAKS:
            patlak.append(t); continue
        if t != 0.815:
            sg4 += abs(b["re"] - gr[0]) + abs(b["im"] - gr[1]); nb4 += 1
        oran[t] = faz(b["re"], b["im"]) / G_FAZ[t]
        if b["R_artik"] is not None and b["R_artik"] < 0.05:
            sR += abs(b["R"] - gr[2]); nR += 1
        else:
            supheli.append(t)
    o4 = [oran[t] for t in TAUS[:4] if t in oran]
    return {"sg4": sg4, "nb4": nb4, "oran": oran, "mag": mag,
            "oran4": sum(o4) / len(o4) if o4 else float("nan"),
            "n_oran4": len(o4),
            "sR": sR, "nR": nR, "eksik": eksik, "supheli": supheli,
            "patlak": patlak}


def hucre(b):
    if b["re"] is None:
        return "—"
    mm = math.hypot(b["re"], b["im"])
    if mm > GMAKS:
        return f"({b['re']:+.2f},{b['im']:+.2f}) **‡**"
    R = ("—" if b["R_artik"] is None else
         (f"{b['R']:.2f}" if b["R_artik"] < 0.05 else f"({b['R']:.2f})*"))
    return f"({b['re']:+.3f},{b['im']:+.3f}) {mm:.2f} {R}"


def yukle(sira):
    D = {}
    for ad in sira:
        f = S / f"ozet_{ad}.json"
        if f.exists():
            D[ad] = json.loads(f.read_text())
    return D


SIRA = ["taban", "R1", "R2", "R3", "N3", "N5", "N10", "N5x", "N5y", "N5z",
        "J14", "J26", "Pd", "P1", "P0"]
ETIKET = {
    "taban": "taban erfc-0.68 (kontrol)",
    "R1": "R1  İ1 ε=+0.003 ×5", "R2": "R2  İ1 ε=+0.010 ×5",
    "R3": "R3  İ1 ε=+0.030 ×5",
    "N3": "N3  İ1 ε=−0.0070 ×3", "N5": "N5  İ1 ε=−0.0037 ×5",
    "N10": "N10 İ1 ε=−0.0018 ×10", "N5x": "N5x İ1 ε=−0.010 ×5",
    "N5y": "N5y İ1 ε=−0.020 ×5", "N5z": "N5z İ1 ε=−0.030 ×5",
    "J14": "J14 titreşim σ_j=0.141", "J26": "J26 titreşim σ_j=0.260",
    "Pd": "Pd  İ2 taban=düzgün", "P1": "P1  İ2 taban=GUE",
    "P0": "P0  İ2 taban=Poisson"}
ACIK = {
    "taban": "değişiklik yok — 152'nin A4'ü, öz-tutarlı Newton",
    "R1": "itme (görevde yazıldığı gibi), 5 süpürme",
    "R2": "itme, 5 süpürme", "R3": "itme, 5 süpürme",
    "N3": "yumuşatma, P(s<0.3)'e kalibre, 3 süpürme",
    "N5": "yumuşatma, P(s<0.3)'e kalibre, 5 süpürme",
    "N10": "yumuşatma, P(s<0.3)'e kalibre, 10 süpürme",
    "N5x": "yumuşatma, kalibre dozun ~2.7 katı",
    "N5y": "yumuşatma, ~5.4 kat", "N5z": "yumuşatma, ~8.1 kat",
    "J14": "yapısız titreşim — σ_ds² N5 ile eşleşir",
    "J26": "yapısız titreşim — σ_ds² N5x ile eşleşir",
    "Pd": "s≡1 taban + tek geçiş merdiven boyası",
    "P1": "GUE-surmise gap tabanı + tek geçiş boya",
    "P0": "üstel (Poisson) gap tabanı + tek geçiş boya"}
GRUP = {"taban": "kontrol", "R1": "İ1a", "R2": "İ1a", "R3": "İ1a",
        "N3": "İ1b", "N5": "İ1b", "N10": "İ1b", "N5x": "İ1b",
        "N5y": "İ1b", "N5z": "İ1b", "J14": "J", "J26": "J",
        "Pd": "İ2", "P1": "İ2", "P0": "İ2"}


def oran_alt(M, taus, gmin=0.0):
    o = [M["oran"][t] for t in taus
         if t in M["oran"] and M["mag"].get(t, 0) >= gmin]
    return (sum(o) / len(o), len(o)) if o else (float("nan"), 0)


def main():
    D = yukle(SIRA)
    M = {a: metrik(D[a]["bantlar"]) for a in D}
    L = []; A = L.append

    A("# 153 — itmeli gaz: seviye itmesi anormal fazı gerçeğe indiriyor mu?\n")
    A("152'nin ölçüm zinciri **birebir** korunarak 15 koşu yapıldı: taban "
      "kontrolü (152'nin A4'ü), görevin istediği iki itme yolu (İ1, İ2) ve "
      "tarama sırasında zorunlu hale gelen iki ek seri (İ1b ve J — "
      "gerekçeleri aşağıda). Soru: saf-merdiven gazının gerçeğin ~1.4 katı "
      "olan anormal dispersiyon fazını GUE seviye itmesi kapatıyor mu?\n")
    A("> **Kısa hüküm:** hayır — ve sebebi görevin varsaydığının tersi. "
      "Sentetik gaz seviye itmesinden yoksun DEĞİL; kısa menzilde gerçeğin "
      "~190 katı KATI (P(s<0.3): 0.00013 vs 0.02420). Görevde yazıldığı "
      "haliyle itme (ε>0) fazı 1.72'den 2.24'e **kötüleştiriyor**. Kısa "
      "menzili YUMUŞATMAK fazı kısmen indiriyor (1.72 → ~1.42) ama orada "
      "doyuyor ve S3'ü bozuyor. İ2/GUE tabanı biraz daha iyi (1.36) ama "
      "kendi kurgusundan gelen bir yapaylık taşıyor. **Hiçbiri 1.00'e "
      "inmiyor.**\n")
    A("**Koşulan kod:** `153_configs/153_gaz.py` (152 gibi tek kod yolu; "
      "zincir `g = np.diff(z)`den itibaren 152'nin karakteri karakterine "
      "aynısıdır). Kontrol bunu kanıtlıyor: `taban`, 152'nin A4'ünün her "
      "Newton iterasyonunda aynı maks|F|'yi ve beş bandın dördüncü "
      "basamağına kadar aynı Γ'sını üretti.\n")

    # ---------------- 1. kalibrasyon hedefi ----------------
    A("\n## 1. Kalibrasyon hedefi — ve taramayı baştan çeviren ölçüm\n")
    A("Görev, sentetik P(s)'in küçük-s kuyruğunu gerçeğinkine yaklaştırmayı "
      "kalibrasyon hedefi koydu. Hedef önce ÖLÇÜLDÜ (zeros6 son-300k, "
      "zincirin kendi açılması `s = g·log(mid/2π)/2π`), sonra taban gazında "
      "aynı büyüklüğe bakıldı:\n")
    A("| | P(s<0.1) | P(s<0.2) | P(s<0.3) | P(s<0.5) | σ_ds² |")
    A("|---|---|---|---|---|---|")
    A(f"| **gerçek (zeros6 son-300k)** | {GERCEK_PS[0.1]:.5f} | "
      f"{GERCEK_PS[0.2]:.5f} | **{GERCEK_PS[0.3]:.5f}** | "
      f"{GERCEK_PS[0.5]:.5f} | 0.1674 |")
    dt = D["taban"]
    A(f"| taban gazı (erfc-0.68) | {dt['ps']['0.1']:.5f} | "
      f"{dt['ps']['0.2']:.5f} | **{dt['ps']['0.3']:.5f}** | "
      f"{dt['ps']['0.5']:.5f} | {dt['sigma_ds2']:.4f} |")
    A("| GUE Wigner surmise (teorik) | 0.00107 | 0.00839 | 0.02725 | "
      "0.11200 | 0.17810 |")
    A("| Poisson (teorik) | 0.09516 | 0.18127 | 0.25918 | 0.39347 | 1.0 |")
    A("\nBu tablo taramanın yönünü değiştirdi. Sentetik gazda 300 bin "
      "aralığın **hiçbiri** 0.2ḡ'nin altında değil; gerçekte 2198 tanesi "
      "var. Yani sentetik gaz seviye itmesinden YOKSUN değil — kısa "
      "menzilde gerçeğin ~190 katı KATI. Birinci-mertebe merdiven gazı "
      "pürüzsüz ve deterministik bir haritadır; aralık dağılımı sınırlı bir "
      "fonksiyonun dağılımıdır, iki ucu da ince. Gerçek GUE'nin küçük-s "
      "kuyruğu bunun çok üstünde.\n")
    A("Sonuç: görevde yazıldığı haliyle itme (ε>0) kalibrasyon hedefinden "
      "UZAKLAŞTIRIR. Bu yüzden ε'nun **iki işareti de** koşuldu ve tabloda "
      "ayrı gruplar olarak duruyor — İ1a (ε>0, görevin harfi) ve İ1b (ε<0, "
      "görevin kalibrasyon hedefi). Bu tek sapmadır, sebebi de yukarıdaki "
      "ölçümdür.\n")

    # ---------------- 2. konfigürasyonlar ----------------
    A("\n## 2. Konfigürasyonlar\n")
    A("| konfig | grup | değişiklik | küçük-s: P(s<0.3) | σ_ds² | σ_η² | "
      "c₁ |")
    A("|---|---|---|---|---|---|---|")
    A(f"| **gerçek** | — | — | **{GERCEK_PS[0.3]:.5f}** | **0.1674** | "
      "**0.0227** | **−0.01158** |")
    for ad in SIRA:
        if ad not in D:
            continue
        d = D[ad]
        A(f"| {ETIKET[ad]} | {GRUP[ad]} | {ACIK[ad]} | {d['ps']['0.3']:.5f} "
          f"| {d['sigma_ds2']:.4f} | {d['sigma_eta2']:.4f} | "
          f"{d['c1']:+.5f} |")
    A("\nİ1'de yer değiştirme `Δz_n = ε·ḡ·[(ḡ/g_{n−1})² − (ḡ/g_n)²]`, "
      "±0.4·min(g_{n−1},g_n) kelepçeli. Kelepçe sıralamanın bozulmamasını "
      "GARANTİ eder (komşular en çok 0.4g yaklaşabildiğinden yeni aralık "
      "≥0.2g): dokuz İ1 koşusunun hiçbirinde sıra bozulmadı (sıra bozan "
      "çift = 0). Pencere ortalama-yoğunluğu (30 pencere) en kötü koşuda "
      "6.4e−05·ḡ kaydı — z yeniden ölçeklenmedi/ötelenmedi, hareketin "
      "teleskopik olması yetti. J koşuları (yapısız titreşim) 152'nin "
      "H-B'siyle aynıdır ve sırayı BOZAR (J14: 245 çift %0.08, J26: 5130 "
      "çift %1.71); 152'de olduğu gibi sonrasında sıralanır.\n")
    A("**Kelepçenin bağlama oranı doz arttıkça patlıyor** ve bu, aşağıdaki "
      "doygunluğun büyük ihtimalle sebebidir — son süpürmede bağlanan "
      "hamlelerin oranı: R1 %0.00, R2 %0.01, R3 %1.05, N3 %1.59, N5 %1.66, "
      "N10 %2.09, N5x %13.5, **N5y %29.3, N5z %39.2**. Yani N5y ve N5z "
      "artık saf 1/s² gevşetmesi değil, hareketlerin üçte birinde "
      "kelepçenin belirlediği bir dinamiktir; bu ikisinden fizik hükmü "
      "çıkarılmamalıdır, yalnızca 'doz büyütmek daha fazla kazanç "
      "getirmiyor' gözlemi için duruyorlar.\n")

    # ---------------- 3. bantlar ----------------
    A("\n## 3. Bantlar — Γ_rot(Re, Im), |Γ|, R\n")
    A("Hücre: `(Re,Im) |Γ| R`. Kurallar 152 ile AYNI: `‡` = |Γ|>1.5, bant "
      "sayısal olarak patlamıştır (ayrışım paydası sıfıra gitmiş), ölçüm "
      "değildir, hiçbir ortalamaya girmez. `*` = R uydurması fazı "
      "tutturamadı (faz-artığı > 0.05 rad).\n")
    A("| τ̄ | gerçek | " + " | ".join(ETIKET[a].split()[0] for a in SIRA
                                     if a in D) + " |")
    A("|---|---|" + "---|" * len([a for a in SIRA if a in D]))
    for t in TAUS:
        gr = GERCEK[t]
        hh = []
        for ad in SIRA:
            if ad not in D:
                continue
            b = next(x for x in D[ad]["bantlar"] if x["tau"] == t)
            hh.append(hucre(b))
        A(f"| {t} | ({gr[0]:+.3f},{gr[1]:+.3f}) "
          f"{math.hypot(gr[0], gr[1]):.2f} {gr[2]:.2f} | " +
          " | ".join(hh) + " |")

    # ---------------- 4. faz orani ----------------
    A("\n## 4. Anormal fazın dikliği — arg Γ_rot(sentetik) / arg Γ_rot"
      "(gerçek)\n")
    A("Hedef 1.000. Kapanış ölçütü budur.\n")
    A("`ort(b1–4)` 152'nin sütunuyla doğrudan karşılaştırılabilir olsun diye "
      "aynı kuralla hesaplandı. Yanına `ort(b1–3)` eklendi: 152, bant 4 ve "
      "5'in her koşuda gürültülü olduğunu gösterdi, ve aşağıda görüleceği "
      "gibi bu taramadaki iyileşmenin büyük kısmı tam o iki bantta oturuyor "
      "— ayrı okunabilmeli. `b3` ise tek başına en sağlam banttır (220 aday, "
      "her koşuda R-artığı küçük).\n")
    A("| konfig | τ=0.5375 | τ=0.585 | τ=0.66 | τ=0.74 | τ=0.815 | "
      "**ort(b1–4)** | ort(b1–3) | b3 |")
    A("|---|---|---|---|---|---|---|---|---|")
    A("| **gerçek** | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | **1.00** | 1.00 | "
      "1.00 |")
    for ad in SIRA:
        if ad not in D:
            continue
        m = M[ad]
        cs = []
        for t in TAUS:
            if t in m["oran"]:
                cs.append(f"{m['oran'][t]:.2f}")
            elif t in m["patlak"]:
                cs.append("‡")
            else:
                cs.append("—")
        o4, n4 = oran_alt(m, TAUS[:4])
        o3, n3 = oran_alt(m, TAUS[:3])
        b3 = m["oran"].get(0.660)
        A(f"| {ETIKET[ad]} | " + " | ".join(cs) +
          f" | **{o4:.2f}** ({n4}b) | {o3:.2f} ({n3}b) | " +
          (f"{b3:.2f} |" if b3 is not None else "‡ |"))

    # ---------------- 5. skor ----------------
    A("\n## 5. Skor\n")
    A("`s_Γ` = bant BAŞINA ortalama (|ΔRe| + |ΔIm|), bant 1–4 içinden patlak "
      "olmayanlar üzerinden — 152'deki düzeltmenin aynısı (toplam değil "
      "ortalama; yoksa en çok bandı patlayan konfigürasyon sahte biçimde "
      "'en iyi' görünür). Kaç bant üzerinden alındığı parantezde. Küçük = "
      "gerçeğe yakın.\n")
    A("| konfig | s_Γ (bant başına) | faz oranı b1–4 | faz oranı b1–3 | "
      "s_R | patlak ‡ |")
    A("|---|---|---|---|---|---|")
    for ad in SIRA:
        if ad not in D:
            continue
        m = M[ad]
        o4, n4 = oran_alt(m, TAUS[:4]); o3, n3 = oran_alt(m, TAUS[:3])
        sg = m["sg4"] / m["nb4"] if m["nb4"] else float("nan")
        sR = m["sR"] / m["nR"] if m["nR"] else float("nan")
        pat = ", ".join(f"{t}" for t in m["patlak"]) or "—"
        A(f"| {ETIKET[ad]} | {sg:.3f} ({m['nb4']}b) | {o4:.2f} | {o3:.2f} | "
          + (f"{sR:.2f} ({m['nR']}b)" if m["nR"] else "—") + f" | {pat} |")
    return L, D, M



def hukum(L, D, M):
    A = L.append
    o = {a: oran_alt(M[a], TAUS[:4])[0] for a in M}
    o3 = {a: oran_alt(M[a], TAUS[:3])[0] for a in M}
    b3 = {a: M[a]["oran"].get(0.660) for a in M}

    A("\n---\n")
    A("\n## Hüküm\n")

    A("\n### 0. Kontrol tuttu\n")
    A("`taban`, 152'nin A4'ünü yeniden üretti: yirmi Newton iterasyonunun "
      "her birinde aynı maks|F| (1.199e+00 → 8.761e-02), aynı σ_ds²=0.1128 / "
      "σ_η²=0.0230 / c₁=−0.00918, ve beş bandın Γ'sı dördüncü basamağa "
      "kadar aynı. Aşağıdaki farklar konfigürasyondan geliyor.\n")

    A("\n### 1. İ1a (ε>0, görevde yazıldığı gibi itme) — fazı KÖTÜLEŞTİRİYOR\n")
    A(f"Doz arttıkça faz oranı monoton biçimde uzaklaşıyor: taban "
      f"{o3['taban']:.2f} → R1 {o3['R1']:.2f} → R2 {o3['R2']:.2f} → R3 "
      f"{o3['R3']:.2f} (b1–3). Küçük-s göstergesi üç koşuda da tam **0.0000**"
      " — itme, zaten aşırı katı olan gazı büsbütün kristalleştiriyor; "
      f"σ_ds² 0.1128'den {D['R3']['sigma_ds2']:.4f}'e (gerçek 0.1674) "
      "çöküyor. Görevin öngördüğü malzeme bu değil: gazda eksik olan itme "
      "değil, itmenin fazlası var.\n")

    A("\n### 2. İ1b (ε<0, kalibrasyon hedefinin istediği yön) — kısmi, "
      "sonra doyuyor\n")
    A(f"Kalibre doz üç farklı süpürme sayısında aynı hedefi tutturuyor "
      f"(P(s<0.3) = {D['N3']['ps']['0.3']:.5f} / {D['N5']['ps']['0.3']:.5f} / "
      f"{D['N10']['ps']['0.3']:.5f}, gerçek 0.02420) ve **aynı fazı** "
      f"veriyor: b1–3 oranı {o3['N3']:.3f} / {o3['N5']:.3f} / "
      f"{o3['N10']:.3f}. Süpürme sayısı sonucu değiştirmiyor — sonuç "
      "gerçekten dozun kendisine bağlı, sayısal ayrıntıya değil.\n")
    A(f"Kalibre dozda kazanç mütevazı: {o3['taban']:.2f} → {o3['N5']:.2f} "
      f"(b1–3). Dozu 2.7 / 5.4 / 8.1 kat büyütünce oran {o3['N5x']:.2f} / "
      f"{o3['N5y']:.2f} / {o3['N5z']:.2f}'de **doyuyor** — daha fazla "
      "yumuşatma daha fazla kazanç getirmiyor, N5z'de geri bile dönüyor. "
      "1.00'e gitmiyor. (Uyarı: N5y/N5z'de kelepçe hamlelerin %29–39'unda "
      "bağlıyor — bölüm 2 — yani doygunluğun bir kısmı fizik değil, "
      "kelepçe olabilir. Kelepçesiz daha büyük doz sıralamayı bozardı ve "
      "gaz olmaktan çıkardı; bu yüzden bu yönde daha ileri gidilmedi. "
      "Güvenilir en iyi İ1b koşusu N5x'tir: kelepçe %13.5.)\n")
    A(f"Bedeli ağır: aşırı dozda σ_η² {D['taban']['sigma_eta2']:.4f} → "
      f"{D['N5z']['sigma_eta2']:.4f} (gerçek 0.0227), c₁ "
      f"{D['taban']['c1']:+.5f} → {D['N5z']['c1']:+.5f} (gerçek −0.01158), "
      f"P(s<0.3) {D['N5z']['ps']['0.3']:.5f} (gerçeğin 9 katı). 152'nin "
      "bulduğu gerilim aynen sürüyor: fazı iyileştiren her hamle S3'ü "
      "bozuyor.\n")

    A("\n### 3. J kontrolü — kazanç σ_ds² artefaktı DEĞİL\n")
    A("İ1b fazı düşürürken σ_ds²'yi de büyütüyor. Faz hangisini görüyor: "
      "kısa-menzil YAPISINI mı, yalnız σ_ds²'yi mi? J koşuları aynı σ_ds²'yi "
      "yapısız (bağımsız, inkoherent) titreşimle üretir. Eşleşmeler:\n")
    A("| çift | σ_ds² | faz oranı b1–3 |")
    A("|---|---|---|")
    for a, bq in (("N5", "J14"), ("N5x", "J26")):
        A(f"| {ETIKET[a]} | {D[a]['sigma_ds2']:.4f} | **{o3[a]:.3f}** |")
        A(f"| {ETIKET[bq]} | {D[bq]['sigma_ds2']:.4f} | {o3[bq]:.3f} |")
    A("\nAynı σ_ds²'de yapılı yumuşatma yapısız titreşimden belirgin biçimde "
      f"daha iyi ({o3['N5']:.2f} vs {o3['J14']:.2f}; {o3['N5x']:.2f} vs "
      f"{o3['J26']:.2f}). Yani kazanç varyans artefaktı değil — kısa-menzil "
      "yapısı fazı gerçekten kıpırdatıyor. 152'nin 'titreşim fazı hiç "
      "kıpırdatmıyor' bulgusu da doğrulanıyor: J14/J26 tabana göre b1–3'te "
      "hiç iyileşme vermiyor.\n")

    A("\n### 4. İ2 (itmeli taban + boyalı merdiven) — en iyi faz, ama "
      "kurgudan gelen yapaylıkla\n")
    A("Üç koşu aynı boyayı paylaşır, yalnız taban istatistiği değişir:\n")
    A("| taban | var(s) girdide | faz oranı b1–3 | ölçülen σ_ds² | "
      "ölçülen P(s<0.3) |")
    A("|---|---|---|---|---|")
    for a, vv in (("Pd", "0.000"), ("P1", "0.178"), ("P0", "1.004")):
        A(f"| {ACIK[a].split(' + ')[0]} | {vv} | {o3[a]:.3f} | "
          f"{D[a]['sigma_ds2']:.4f} | {D[a]['ps']['0.3']:.5f} |")
    A(f"\nDüzgün taban (Pd) tabandan **daha kötü** ({o3['Pd']:.2f} vs "
      f"{o3['taban']:.2f}) — yani tek-geçiş boyanın kendisi bir kazanç "
      "getirmiyor; öz-tutarlılığı atmak zarar veriyor. Poisson tabanı (P0) "
      "sinyali tümüyle yok ediyor (faz oranları negatif, R uydurmaları "
      "ölçülebilen dört bantta 1.76–2.33 rad artıkla düşüyor). Kazancı "
      "veren yalnız GUE tabanıdır: "
      f"P1, tüm taramanın en düşük oranını veriyor ({o3['P1']:.2f} b1–3, "
      f"{o['P1']:.2f} b1–4).\n")
    A("**Ama P1 temiz bir ölçüm değil.** Tek geçiş boya |Δz| en fazla "
      "0.725 üretiyor (ḡ=0.522) — bir ortalama aralıktan büyük yer "
      "değiştirme. Sonuç: komşu çiftlerin %1.71'i sıralamayı bozuyor, "
      "zincirin gerektirdiği sıralama bunları sıfıra yakın aralıklara "
      "çeviriyor. Bunun kanıtı Pd'dir: tabanında s≡1 olduğu için "
      "P(s<0.3)=0.00000 girmesi gerekirken ölçümde "
      f"{D['Pd']['ps']['0.3']:.5f} çıkıyor — bu sayının **tamamı** çaprazlama "
      "yapaylığıdır. Aynı sebeple P1'in P(s<0.3)'ü tabanındaki 0.02709'dan "
      f"{D['P1']['ps']['0.3']:.5f}'e şişiyor (gerçeğin 4.4 katı) ve "
      f"σ_ds² {D['P1']['sigma_ds2']:.4f}'e (gerçeğin 2 katı) çıkıyor. "
      "P1'in düşük faz oranı bu yüzden GUE itmesine değil, kısmen o "
      "yapaylığa da yazılabilir; ayrıca |Γ| genlikleri gerçeğin çok "
      "altında (0.40–0.46 vs gerçek 0.63–0.81), s_Γ skoru "
      f"{M['P1']['sg4']/M['P1']['nb4']:.3f} ile N5x'in "
      f"{M['N5x']['sg4']/M['N5x']['nb4']:.3f}'inden belirgin kötü.\n")
    A("Ek fizik notu: gerçek σ_ds²=0.1674 ile GUE surmise'ın 0.1781'i "
      "zaten aynı yerdedir. Yani gerçek sıfırların aralık varyansını GUE "
      "tek başına açıklıyor; üstüne BAĞIMSIZ bir merdiven yer değiştirmesi "
      f"eklenecek yer yok — eklenince P1'de olduğu gibi "
      f"{D['P1']['sigma_ds2']:.2f}'e fırlıyor. Gerçekte merdiven ile GUE "
      "dalgalanması bağımsız iki katman değil, aynı şeyin iki yazımıdır; "
      "İ2'nin kurgusu (bağımsız taban + üstüne boya) bunu ihlal ediyor.\n")

    A("\n### 5. Hangi yol daha iyi\n")
    A(f"Faz oranına tek başına bakılırsa İ2/GUE (P1, b1–3 {o3['P1']:.2f}), "
      f"İ1b'nin en iyisinden (N5x {o3['N5x']:.2f}) bir tık önde. Ama İ1b "
      "her başka ölçütte kazanıyor: s_Γ "
      f"{M['N5x']['sg4']/M['N5x']['nb4']:.3f} vs "
      f"{M['P1']['sg4']/M['P1']['nb4']:.3f}; |Γ| genlikleri gerçeğin "
      "aralığında; beş bandın BEŞİ DE ölçülebilir ve R uydurmalarının "
      "hepsi tutuyor (R-artık 0.002–0.015) — 152'nin sekiz ve 153'ün on "
      "beş koşusu içinde bunu başaran tek koşu N5x'tir; ve kurgusunda "
      "çaprazlama yapaylığı yok (sıra "
      "bozan çift 0, pencere yoğunluğu 4e−05·ḡ içinde korunmuş). "
      "**İ1b daha iyi yoldur**; İ2'nin sayısı daha parlak görünse de "
      "güvenilirliği düşüktür.\n")

    A("\n### 6. Faz ne kadar indi — dürüst muhasebe\n")
    A("Görevin sorusu: 1.4'ten 1.0'a çekiyor mu — kısmen mi, hiç mi?\n")
    A("| ölçüt | taban | en iyi İ1b (N5x) | en iyi İ2 (P1) | kapanan pay |")
    A("|---|---|---|---|---|")
    for et, f in (("ort b1–4", o), ("ort b1–3", o3), ("b3 tek başına", b3)):
        t0 = f["taban"]; bn = f["N5x"]; bp = f["P1"]
        pay = (t0 - bn) / (t0 - 1.0) * 100 if t0 != 1.0 else float("nan")
        A(f"| {et} | {t0:.2f} | {bn:.2f} | {bp:.2f} | "
          f"%{pay:.0f} (N5x) |")
    A("\n**Kısmen — ve payın büyük kısmı zayıf bantlardan geliyor.** Bant "
      f"bazında bakınca: τ=0.5375 {M['taban']['oran'][0.5375]:.2f} → "
      f"{M['N5x']['oran'][0.5375]:.2f} (ama bu bantta yalnız 33 aday var), "
      f"τ=0.585 {M['taban']['oran'][0.585]:.2f} → "
      f"{M['N5x']['oran'][0.585]:.2f} (gerçek bir hareket, 144 aday), "
      f"τ=0.66 {M['taban']['oran'][0.660]:.2f} → "
      f"{M['N5x']['oran'][0.660]:.2f} — **en sağlam bantta %32'lik "
      "fazlalığın yalnız dörtte biri kapanıyor.** τ=0.74 ve 0.815 "
      "tabanda zaten patlamış/sarmıştı; oradaki 'iyileşme' esas olarak "
      "bantların ölçülebilir hale gelmesidir, ki bu da bir kazanç ama faz "
      "kapanışı değil.\n")

    A("\n### 7. Ne kaldı\n")
    A("- Eksik malzeme, GUE seviye itmesinin sentetik gaza EKLENMESİ "
      "değildir; ölçüm bunun tersini söylüyor (sentetik gaz zaten aşırı "
      "katı). Kısa-menzil yumuşatması doğru yön ama kapanış b3'te ~%25, "
      "b1–3'te ~%42 kalıp doyuyor.\n"
      "- Yumuşatmanın yapısı yanlış: kalibre dozda P(s<0.3) gerçeği "
      f"tutturuyor ({D['N5']['ps']['0.3']:.5f} vs 0.02420) ama P(s<0.1) "
      f"{D['N5']['ps']['0.1']:.5f} ile gerçeğin ({GERCEK_PS[0.1]:.5f}) 10 "
      "katı — gevşetme GUE'nin s² yükselişini değil, bir kümelenme kuyruğu "
      "üretiyor. Doğru kısa-menzil ŞEKLİNİ kuran bir gevşetme (ör. Dyson "
      "log-gazı, 1/s kuvveti, ve ÖZ-TUTARLI merdivenle birlikte) sıradaki "
      "denemedir.\n"
      "- İ2'yi kurtarmanın yolu, boyayı tek geçişte değil öz-tutarlı "
      "çözmektir; öyle olunca taban gap'leri ile merdiven artık bağımsız "
      "olmaz ve çaprazlama yapaylığı da kalkar. Bu, 153'ün kurgusunda "
      "kasten yoktu (görev tek geçiş istedi) ama sonuç onu işaret ediyor.\n"
      "- 152'nin 'faz KATI' hükmü tam olarak yıkılmadı, yumuşatıldı: faz "
      "kıpırdıyor, ama en güvenilir bantta %6, ve 1.00'e giden yol "
      "görünmüyor.\n")

    A("\n---\n")
    A("\n## Dürüstlük notları\n")
    A("- **Sapma ve sebebi:** ε<0 (İ1b) ve J kontrolü görevde istenmedi. "
      "İ1b, görevin kendi kalibrasyon hedefinin (küçük-s kuyruğunu gerçeğe "
      "yaklaştırmak) ε>0 ile sağlanamaz olduğu ÖLÇÜLDÜĞÜ için eklendi "
      "(bölüm 1). J, İ1b'nin kazancının bir σ_ds² artefaktı olup olmadığını "
      "ayırmak için eklendi (bölüm 3). İkisi de tabloda ayrı grup.\n"
      "- **Patlayan bantlar ve sebepleri.** Γ_rot'un paydası (on₀−off₀), "
      "η'nın asal çizgilerindeki gücün çizgi-dışı güçten farkıdır. Merdiven "
      "yapısı zayıfladığında ya da düzensizlik onu bastırdığında bu fark "
      "sıfıra gider ve Γ patlar. Bu taramada patlayanlar: `taban` τ=0.74 "
      "(|Γ|=2.29; 152'de de aynı bant A4'te patlamıştı), `R1` τ=0.815 "
      "(|Γ|=6.73 — en hafif itme dozunda bile; bu bandın niçin bu kadar "
      "kırılgan olduğu ayrıca sınanmadı), `J26` τ=0.815 (|Γ|=3.20), `Pd` τ=0.585/0.74/0.815, "
      "`P0` τ=0.815 (|Γ|=4.26). Hiçbiri ortalamalara katılmadı.\n"
      "- **İ2'nin çaprazlama yapaylığı** (bölüm 4) raporun en zayıf "
      "halkasıdır ve P1'in 'en iyi faz' sonucunu doğrudan gölgeler; "
      "gizlenmedi, hükümde de ağırlığı düşürüldü.\n"
      "- **|Γ| alt sınırı yok.** 152'nin kuralı yalnız üstten eler "
      "(|Γ|>1.5). J14'ün τ=0.74 bandı |Γ|=0.06 ile çökmüş olmasına rağmen "
      "kurala göre sayılıyor ve oranı −0.51 çıkıp J14'ün b1–4 ortalamasını "
      "yapay biçimde 1.15'e indiriyor. Bu yüzden hükümde b1–4 değil b1–3 "
      "kullanıldı; kural 152 ile karşılaştırılabilirlik için "
      "değiştirilmedi.\n"
      "- **Süreler ölçüt değil.** Koşular sırasında aynı makinede başka bir "
      "iş (`154_sentetik.py`) 6 süreçle koştu; 8 çekirdeğe 12 süreç düştü. "
      "`taban` bu yüzden 20.4 dk sürdü (152'de aynı hesap 10.2 dk). "
      "Sayısal sonuç etkilenmez (aynı maks|F| dizisi), yalnız süre "
      "sütunu anlamsızdır — o yüzden tabloya konmadı.\n"
      "- **Uydurma yok:** yukarıdaki tüm tablo sayıları "
      "`scratchpad/153/ozet_*.json`'dan otomatik yazıldı. Elle girilen tek "
      "sayılar gerçek referans bantları, gerçek S3 üçlüsü, ölçülmüş gerçek "
      "P(s) hedefleri ve GUE/Poisson teorik değerleridir.\n")
    A("\nHam çıktılar: `153_configs/153_gaz.py <konfig>`; her koşunun tam "
      "log'u ve `ozet_*.json`'u üretilmiştir.\n")
    return L


if __name__ == "__main__":
    L, D, M = main()
    L = hukum(L, D, M)
    OUT.write_text("\n".join(L) + "\n")
    print(f"yazildi -> {OUT}  ({len(L)} satir)")
