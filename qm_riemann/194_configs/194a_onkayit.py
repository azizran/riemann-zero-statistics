# -*- coding: utf-8 -*-
"""
194a — K0: KURAL-ÖNCE ÖN-KAYIT (BOŞ PENCERE TABANI — kalıntı-sınıfından ve
uydudan bağımsız eksi taban hipotezi; KALEM_BOS_PENCERE_TABANI_24EYL2026 AYNEN)
================================================================================
ÖLÇÜMDEN ÖNCE donan (yalnız aritmetik — HİÇBİR veri/scratchpad-sonuç okunmaz):
 (i)  BOŞ PENCERE listesi: Δω ∈ [0.30, 2.30] içinde, a,b ≤ 12 sade kesirlerin
      log(a/b) konumlarının HEPSİNDEN (μ(a)≠0 olanlar VE μ(a)=0 olanlar)
      ≥ 0.06 uzak merkezler (boşlukların ortası; birbirinden ≥0.06 ayrık).
      ≥6 pencere çıkmazsa sınır a,b≤10'a iner (kaydedilir). Her boş pencere
      için en yakın kataloğ komşusu (sol+sağ, eşit uzaklıkta — orta nokta
      olduğu için) ve uzaklığı.
 (ii) μ=0 pencereleri KALEM'deki listeden AYNEN (Δω aralığı dışındakiler
      KAYIT — ölçülmez, yalnız not düşülür).
 (iii) H-194a..e eşikleri KALEM'den AYNEN.

Çıktı: scratchpad/194/ONKAYIT_194.json (sha256 = bu betiğin sha'sı + damga).
"""
import hashlib
import json
import subprocess
import time
from math import gcd, log
from pathlib import Path

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S194 = SCR / "194"
S194.mkdir(exist_ok=True)

DELTA = 0.03                 # uydu/pencere seçim yarı-genişliği (blok-yerel Δω)
UZAKLIK_ESIK = 0.06            # KALEM AYNEN: merkezden HER kataloğ noktasına >= bu kadar uzak
RANGE_LO, RANGE_HI = 0.30, 2.30


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def mu(n):
    n = int(n)
    if n == 1:
        return 1
    k, p, m = 0, 2, n
    while p * p <= m:
        if m % p == 0:
            m //= p
            if m % p == 0:
                return 0
            k += 1
        p += 1
    if m > 1:
        k += 1
    return (-1) ** k


def phi(n):
    n = int(n)
    res = n
    p = 2
    m = n
    while p * p <= m:
        if m % p == 0:
            while m % p == 0:
                m //= p
            res -= res // p
        p += 1
    if m > 1:
        res -= res // m
    return res


# ---------------- (i) BOŞ PENCERE: katalog + boşluk-ortası arama ----------------
def katalog(maxd):
    """a,b <= maxd, gcd(a,b)=1, a>b (pozitif log) tüm sade kesirlerin log(a/b)
    konumları — μ(a) DEĞERİNE BAKILMAKSIZIN (hem μ≠0 hem μ=0)."""
    pts = []
    for a in range(1, maxd + 1):
        for b in range(1, a):
            if gcd(a, b) == 1:
                pts.append({"pos": log(a / b), "a": a, "b": b, "mu_a": mu(a)})
    return sorted(pts, key=lambda x: x["pos"])


def bos_pencereler_bul(maxd):
    pts = katalog(maxd)
    merkezler = []
    for i in range(len(pts) - 1):
        p0, p1 = pts[i], pts[i + 1]
        bosluk = p1["pos"] - p0["pos"]
        if bosluk >= 2 * UZAKLIK_ESIK - 1e-12:
            mid = (p0["pos"] + p1["pos"]) / 2
            if RANGE_LO <= mid <= RANGE_HI:
                uzaklik = bosluk / 2
                merkezler.append({
                    "merkez": mid, "bosluk_genisligi": bosluk,
                    "en_yakin_komsu_sol": {**p0, "uzaklik": uzaklik},
                    "en_yakin_komsu_sag": {**p1, "uzaklik": uzaklik}})
    # birbirinden >=0.06 ayrık mı (komşu merkezler arası) — kontrol
    ms = sorted(m["merkez"] for m in merkezler)
    mutual = [ms[i + 1] - ms[i] for i in range(len(ms) - 1)]
    return merkezler, mutual, len(pts)


BOS12, MUTUAL12, NKAT12 = bos_pencereler_bul(12)
sinir_kullanilan = 12
dusme_notu = None
if len(BOS12) >= 6:
    BOS = BOS12
else:
    BOS10, MUTUAL10, NKAT10 = bos_pencereler_bul(10)
    sinir_kullanilan = 10
    dusme_notu = (f"a,b<=12 katalogunda ({NKAT12} nokta) yalnız {len(BOS12)} "
                  f"boş pencere çıktı (<6) -> sınır a,b<=10'a indirildi "
                  f"({NKAT10} nokta): {len(BOS10)} boş pencere.")
    BOS = BOS10
    MUTUAL = MUTUAL10
    if len(BOS10) < 6:
        dusme_notu += (f" UYARI: a,b<=10'da da 6'nın altında (N={len(BOS10)}); "
                       "KALEM'de bundan öte bir gevşetme adımı YOK — eşik "
                       "GEVŞETİLMEDİ, mevcut pencerelerle devam edilir.")
if dusme_notu is None:
    MUTUAL = MUTUAL12

MUTUAL_OK = all(d >= UZAKLIK_ESIK - 1e-9 for d in MUTUAL) if MUTUAL else True

BOS_OUT = []
for i, m in enumerate(BOS):
    BOS_OUT.append({
        "ad": f"bos{i+1}", "merkez": m["merkez"],
        "sinir_ab": sinir_kullanilan,
        "bosluk_genisligi": m["bosluk_genisligi"],
        "komsu_sol": {"a": m["en_yakin_komsu_sol"]["a"],
                      "b": m["en_yakin_komsu_sol"]["b"],
                      "pos": m["en_yakin_komsu_sol"]["pos"],
                      "mu_a": m["en_yakin_komsu_sol"]["mu_a"],
                      "uzaklik": m["en_yakin_komsu_sol"]["uzaklik"]},
        "komsu_sag": {"a": m["en_yakin_komsu_sag"]["a"],
                      "b": m["en_yakin_komsu_sag"]["b"],
                      "pos": m["en_yakin_komsu_sag"]["pos"],
                      "mu_a": m["en_yakin_komsu_sag"]["mu_a"],
                      "uzaklik": m["en_yakin_komsu_sag"]["uzaklik"]}})

# ---------------- (ii) μ=0 pencereleri (KALEM AYNEN) ----------------
MU0_TANIM = [
    ("+log(8/5)", 8, 5), ("+log(9/5)", 9, 5), ("+log(9/4)", 9, 4),
    ("+log(8/3)", 8, 3), ("+log4", 4, 1), ("+log(9/2)", 9, 2),
    ("+log8", 8, 1), ("+log9", 9, 1),
    ("+log(4/3)", 4, 3),   # KALEM'de listeli ama Δω aralığı DIŞINDA -> KAYIT
]
MU0_OUT = []
for ad, a, b in MU0_TANIM:
    assert mu(a) == 0, (ad, a, "mu(a) != 0 olmamalı")
    pos = log(a / b)
    ici = RANGE_LO <= pos <= RANGE_HI
    MU0_OUT.append({"ad": ad, "a": a, "b": b, "merkez": pos,
                    "durum": "aralik_ici" if ici else "KAYIT_aralik_disi"})

N_MU0_ICI = sum(1 for m in MU0_OUT if m["durum"] == "aralik_ici")

# ---------------- (iii) H-194a..e eşikleri (KALEM AYNEN) ----------------
H_194A = {
    "hipotez": "EKSİ TABAN VAR: boş pencerelerin ortalaması kappa_bar_bos",
    "bant": [-0.0035, -0.0015],
    "anlamlilik": ">=3 sigma negatif",
    "olum": "|kappa_bar_bos| < 2*se YA DA kappa_bar_bos > 0 (-> sızıntı rakibi kazanır)",
    "kurtarma": "YOK"}
H_194B = {
    "hipotez": "SABİTLİK: boş pencerelerde Delta-omega ile eğilim yok",
    "kosul": "doğrusal eğim sıfırla 2 sigma içinde",
    "statu": "KAYIT (eğim varsa biçimi raporlanır, ölüm yok)"}
H_194C = {
    "hipotez": "SINIFTAN BAĞIMSIZ: boş pencerelerde mod 3,5,10 sınıf payları eşit",
    "kosul": "s_r = 1/phi(a) +- [0.10*(1/phi(a)) + 2*se]",
    "statu": "cos deseni YOK beklentisi; eşiksiz hüküm ama bant kontrolü raporlanır"}
H_194D = {
    "hipotez": "ÇUKUR = TABAN: mu=0 pencerelerinin ortalaması kappa_bar_bos ile uyumlu",
    "kosul": "2 sigma içinde (sigma = sqrt(se_mu0^2+se_bos^2))",
    "olum": "kappa_bar_mu0, kappa_bar_bos'tan >=3 sigma daha derin (negatif) (-> çukurlar ayrı olgu)",
    "kurtarma": "YOK"}
H_194E = {
    "hipotez": "GERME KAPANIŞI: ÖLÇÜLEN kappa_bar_bos ile düzeltilmiş 193 payları "
               "(+log10,+log7,+log5) cos yasası öngörüsünün +-0.10 içinde",
    "duzeltme_formulu": "s_r' = (kappa_r - kappa_bar_bos/phi(a)) / (kappa_top - kappa_bar_bos)",
    "bant": 0.10,
    "olum": "+log10'da herhangi bir düzeltilmiş pay +-0.15 dışında",
    "kurtarma": "YOK"}

# ---------------- girdi dosyaları (denetim; okuma DEĞİL, yalnız sha) ----------
GIRDI_SHA = {
    "KALEM_BOS_PENCERE_TABANI_24EYL2026.md": sha(QM / "KALEM_BOS_PENCERE_TABANI_24EYL2026.md"),
    "193_sinif_yasasi_RAPOR.md": sha(QM / "193_sinif_yasasi_RAPOR.md"),
    "193_configs/193a_onkayit.py": sha(QM / "193_configs" / "193a_onkayit.py"),
    "193_configs/193b_sinif.py": sha(QM / "193_configs" / "193b_sinif.py")}

ONKAYIT = {
    "gorev": "194 — BOŞ PENCERE TABANI: kalıntı-sınıfından/uydudan bağımsız eksi taban mı?",
    "kalem": "KALEM_BOS_PENCERE_TABANI_24EYL2026.md AYNEN (commit 9e74e4a)",
    "girdi_sha256": GIRDI_SHA,
    "delta": DELTA,
    "gap_esik": UZAKLIK_ESIK,
    "aralik": [RANGE_LO, RANGE_HI],
    "sinir_ab_kullanilan": sinir_kullanilan,
    "sinir12_pencere_sayisi": len(BOS12),
    "dusme_notu": dusme_notu,
    "mutual_spacing_ok": bool(MUTUAL_OK),
    "mutual_diffs": MUTUAL,
    "bos_pencereler": BOS_OUT,
    "n_bos": len(BOS_OUT),
    "mu0_pencereler": MU0_OUT,
    "n_mu0_ici": N_MU0_ICI,
    "H_194a": H_194A, "H_194b": H_194B, "H_194c": H_194C,
    "H_194d": H_194D, "H_194e": H_194E,
    "phi_tablosu": {"3": phi(3), "5": phi(5), "7": phi(7), "10": phi(10)},
    "kurallar": ("TEK DALGA; ölümler kurtarmasız; eşik gevşetme yok; ölçülemeyen "
                 "'erişilemedi'; ön-kayıtlı kapılar (makine mührü dahil) "
                 "ATLANMAZ; git'e dokunulmaz; sonuç ORTAK TEFTİŞE."),
    "zaman": None,
    "sha256": None,
}

if __name__ == "__main__":
    try:
        ONKAYIT["zaman"] = subprocess.check_output(["date"], text=True).strip()
    except Exception:
        ONKAYIT["zaman"] = time.strftime("%a %b %d %H:%M:%S %Z %Y")
    s = hashlib.sha256(Path(__file__).resolve().read_bytes()).hexdigest()
    ONKAYIT["sha256"] = s
    yol = S194 / "ONKAYIT_194.json"
    if yol.exists():
        raise SystemExit(f"ONKAYIT_194 ZATEN VAR — üzerine yazılmaz: {yol}")
    json.dump(ONKAYIT, open(yol, "w"), indent=1, ensure_ascii=False)
    print("=" * 78)
    print("194a / K0 ÖN-KAYIT yazıldı")
    print("=" * 78)
    print(f"  sha256 = {s}")
    print(f"  damga  = {ONKAYIT['zaman']}")
    print(f"  sınır a,b<={sinir_kullanilan}  (a,b<=12 sonucu: {len(BOS12)} pencere)")
    if dusme_notu:
        print(f"  DÜŞME NOTU: {dusme_notu}")
    print(f"  mutual_spacing_ok = {MUTUAL_OK}  diffs={[round(d,4) for d in MUTUAL]}")
    print(f"\n  BOŞ PENCERELER (N={len(BOS_OUT)}):")
    for b in BOS_OUT:
        print(f"    {b['ad']}: merkez={b['merkez']:.5f}  bosluk={b['bosluk_genisligi']:.5f}  "
              f"sol={b['komsu_sol']['a']}/{b['komsu_sol']['b']}(pos={b['komsu_sol']['pos']:.5f},"
              f"mu={b['komsu_sol']['mu_a']},uzk={b['komsu_sol']['uzaklik']:.5f})  "
              f"sag={b['komsu_sag']['a']}/{b['komsu_sag']['b']}(pos={b['komsu_sag']['pos']:.5f},"
              f"mu={b['komsu_sag']['mu_a']},uzk={b['komsu_sag']['uzaklik']:.5f})")
    print(f"\n  MU=0 PENCERELERİ (N_ici={N_MU0_ICI}, toplam={len(MU0_OUT)}):")
    for m in MU0_OUT:
        print(f"    {m['ad']} (a={m['a']},b={m['b']}): merkez={m['merkez']:.5f}  {m['durum']}")
    print(f"\n  -> {yol}")
