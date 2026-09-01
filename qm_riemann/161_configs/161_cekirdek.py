"""
161 — ÇEKİRDEK: δ(τ_eff) = φ − (4π·τ_eff − 2π) tayfı ve kuadratik fit
====================================================================
Bu dosya HİÇBİR YENİ ÖLÇÜM YAPMAZ. 158'in kayıtlı çıktılarını
(scratchpad/158/F_<veri>_t<taban>.json — 46 koşu, bant düzeyinde
tau_eff / phi / sPhi_jk / phi_jk / teff_jk) ve gerektiğinde 155/157'nin
kayıtlı koşularını okur; hepsi TABLODAN-OKUNAN değil, ham JSON'dan
yeniden hesaplanandır.

KİMLİK (159 §2):
    φ_Γ = (4π·τ_eff − 2π) + δ          ← omurga τ'da TAM doğrusal
    δ   = arg[zp·conj(zc)·e^{−iA}]

Omurga model uzayında (derece ≥ 1) olduğu için, AYNI ağırlıkla yapılan
ağırlıklı polinom fitleri ÖZDEŞ olarak şu şekilde bağlıdır:
    c_φ = c_δ + (0, 4π, −2π)
yani
    b(φ) ≡ b(δ)                       (kuadratikte c₀ aynı)
    a    = c₁(φ) + 2c₀τ₀ = 4π + dδ/dτ|_{τ₀}      (ÖZDEŞ)
    τ₀   : φ'nin kökü ⇒ 4πτ₀ − 2π + δ(τ₀) = 0 ⇒ τ₀ = ½ − δ(τ₀)/4π (ÖZDEŞ)
         ve doğrusallaştırılmış hâli τ₀ = ½ − δ(½)/a  (159'un yazdığı)

FİT KONVANSİYONU — 158_analiz.py'nin `_fit`/`secim`/`olc`'siyle BİREBİR:
  apsis x = τ_eff · ağırlık w = 1/σ (varsayılan σ_φ jackknife)
  np.polyfit(x, y, derece, w=1/σ) · b = ½·y''(çapa)
  jackknife: 8 grup, her replikada (teff_jk, phi_jk) ile fit tekrarı

161'İN TEK EKLEMESİ: okuma ÇAPASI. 158 τ₀'da okuyordu (gaza bağlı nokta);
161 ÇAPAYI τ = ½'ye alıyor — omurganın KENDİ sıfırı, gazdan bağımsız,
konvansiyonsuz. Üçlü:
    δ(½) ,  dδ/dτ|½ ,  b
Türetilmişler (158'e köprü):
    a_½  = 4π + dδ/dτ|½        a_τ₀ = 4π + dδ/dτ|τ₀   (= 158'in a'sı, ÖZDEŞ)
    τ₀*  = ½ − δ(½)/a_½
"""
import json
from pathlib import Path

import numpy as np

TWO_PI = 2.0 * np.pi
FOUR_PI = 4.0 * np.pi

SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
SCR158, SCR157, SCR155 = SCR / "158", SCR / "157", SCR / "155"
SCR161 = SCR / "161"

TABAN = (0.28, 0.34, 0.40, 0.46, 0.52)
# 155/158'in korelasyon merdiveni sırası (τ₀ artan)
MERDIVEN = ("keskin", "A4", "J14", "J26", "N5", "N5z", "P1", "son", "orta")
GERCEK = ("son", "orta")
SENTETIK = ("keskin", "A4", "J14", "J26", "N5", "N5z", "P1")

# 158 §4 W-A tablosu — YALNIZ tutarlılık denetimi için alıntı (tablodan okunan)
REF158_WA = {
    ("keskin", 0.28): (0.4882, 11.667, -4.120), ("keskin", 0.34): (0.4885, 11.708, -3.843),
    ("keskin", 0.40): (0.4886, 11.660, -3.811), ("keskin", 0.46): (0.4886, 11.505, -3.909),
    ("keskin", 0.52): (0.4871, 10.654, +1.317),
    ("A4", 0.28): (0.4906, 13.094, +10.060), ("A4", 0.34): (0.4904, 13.125, +10.428),
    ("A4", 0.40): (0.4897, 13.054, +5.541), ("A4", 0.46): (0.4892, 13.019, -1.756),
    ("A4", 0.52): (0.4842, 10.773, +8.748),
    ("J14", 0.28): (0.4912, 13.009, +14.332), ("J14", 0.34): (0.4905, 13.063, +11.532),
    ("J14", 0.40): (0.4898, 12.987, +6.212), ("J14", 0.46): (0.4893, 12.876, +0.815),
    ("J14", 0.52): (0.4890, 12.136, +2.313),
    ("J26", 0.28): (0.4881, 12.817, +10.618), ("J26", 0.34): (0.4876, 12.912, +9.409),
    ("J26", 0.40): (0.4876, 12.834, +7.330), ("J26", 0.46): (0.4867, 12.424, +4.056),
    ("J26", 0.52): (0.4915, 13.381, -3.623),
    ("N5", 0.28): (0.4959, 12.465, +2.981), ("N5", 0.34): (0.4952, 12.621, +5.106),
    ("N5", 0.40): (0.4947, 12.599, +4.092), ("N5", 0.46): (0.4941, 12.426, -0.783),
    ("N5", 0.52): (0.4907, 10.818, +6.226),
    ("N5z", 0.28): (0.4941, 10.910, -4.558), ("N5z", 0.34): (0.4942, 11.161, +1.066),
    ("N5z", 0.40): (0.4949, 11.329, +5.761), ("N5z", 0.46): (0.4949, 10.785, +10.086),
    ("N5z", 0.52): (0.4861, 8.131, +21.047),
    ("P1", 0.28): (0.4971, 12.181, +7.920), ("P1", 0.34): (0.4986, 12.166, +8.638),
    ("P1", 0.40): (0.4996, 11.988, +6.835), ("P1", 0.46): (0.5026, 11.227, +18.000),
    ("P1", 0.52): (0.4890, 6.191, +40.940),
    ("son", 0.28): (0.4999, 10.790, -5.890), ("son", 0.34): (0.5043, 10.689, -6.877),
    ("son", 0.40): (0.5088, 10.615, -7.759), ("son", 0.46): (0.5120, 10.724, -7.502),
    ("son", 0.52): (0.5127, 11.013, -7.272),
    ("orta", 0.28): (0.4995, 10.814, -6.495), ("orta", 0.34): (0.5047, 10.695, -7.123),
    ("orta", 0.40): (0.5087, 10.656, -7.801), ("orta", 0.46): (0.5121, 10.753, -8.338),
    ("orta", 0.52): (0.5124, 10.845, -6.152),
}
# 159 §2c/§2d — δ tablosu (tablodan okunan; V-DENETİM'de karşılaştırılıyor)
REF159_DELTA_HALF = {"son": -0.0955, "orta": -0.0966, "keskin": +0.1352,
                     "A4": +0.1364, "P1": -0.0087}
REF159_DDELTA = {"son": -1.990, "orta": -1.981, "keskin": -0.882,
                 "A4": +0.584, "P1": -0.559}


# --------------------------------------------------------------- yükleme
def yukle(kok=SCR158, kalip="F_{v}_t{t}.json", veriler=None, tabanlar=TABAN):
    D = {}
    for v in (veriler or (MERDIVEN + ("P0",))):
        for t in tabanlar:
            for isim in (kalip.format(v=v, t=t), kalip.format(v=v, t=f"{t:g}")):
                p = kok / isim
                if p.exists():
                    D[(v, t)] = json.load(open(p))
                    break
    return D


def bant(d):
    return [b for b in d["bantlar"] if b.get("olculdu")]


def omurga(x):
    """Kinematik omurga: 4π·τ_eff − 2π (159 §2a)."""
    return FOUR_PI * np.asarray(x, float) - TWO_PI


def delta_bant(bs):
    """Bant listesinden (x, δ, σ_φ, σ_δ, δ_jk, teff_jk, φ, τ̄) döndür.

    σ_δ, δ'nın KENDİ jackknife'ıdır: δ_jk^(k) = φ_jk^(k) − omurga(τeff_jk^(k)).
    φ ve τ_eff aynı çizgilerden geldiği için güçlü korelelidirler; σ_δ
    dolayısıyla σ_φ'den farklıdır ve ayrı ölçülür.
    """
    x = np.array([b["tau_eff"] for b in bs], float)
    ph = np.array([b["phi"] for b in bs], float)
    sph = np.array([b["sPhi_jk"] for b in bs], float)
    tb = np.array([b["tau"] for b in bs], float)
    dl = ph - omurga(x)
    njk = min(len(b["phi_jk"]) for b in bs)
    PJ = np.array([b["phi_jk"][:njk] for b in bs], float)     # (nb, njk)
    TJ = np.array([b["teff_jk"][:njk] for b in bs], float)
    DJ = PJ - omurga(TJ)
    sdl = np.sqrt((njk - 1) / njk * np.sum((DJ - DJ.mean(1, keepdims=True))**2, 1))
    return dict(x=x, phi=ph, d=dl, s_phi=sph, s_d=sdl, DJ=DJ, TJ=TJ, PJ=PJ,
                tau=tb, njk=njk, absG=np.array([b["absG"] for b in bs], float))


# ------------------------------------------------------------------- fit
def _fit(xx, yy, ee, derece, capa=0.5):
    """158_analiz._fit ile aynı çekirdek; okuma çapası eklendi.

    Döndürür: c (katsayılar), çapada değer/türev/eğrilik, kök (τ₀), χ².
    """
    m = np.isfinite(xx) & np.isfinite(yy)
    xx, yy, ee = np.asarray(xx)[m], np.asarray(yy)[m], np.asarray(ee)[m]
    if len(xx) < derece + 2:
        return None
    ee = np.where(np.isfinite(ee) & (ee > 0), ee, np.nanmedian(ee))
    c = np.polyfit(xx, yy, derece, w=1.0 / ee)
    d1, d2 = np.polyder(c), np.polyder(c, 2)
    r = np.roots(c)
    r = np.array([z.real for z in r if abs(z.imag) < 1e-9])
    t0 = float(r[np.argmin(np.abs(r - xx.mean()))]) if len(r) else float("nan")
    chi2 = float(np.sum(((yy - np.polyval(c, xx)) / ee)**2))
    return dict(c=c, n=len(xx), chi2=chi2, t0=t0,
                # ÇAPADA (τ = capa) okunan üçlü
                v0=float(np.polyval(c, capa)),
                v1=float(np.polyval(d1, capa)),
                b=float(0.5 * np.polyval(d2, capa)),
                # τ₀'da okunan (158'in konvansiyonu)
                v1_t0=float(np.polyval(d1, t0)) if np.isfinite(t0) else float("nan"),
                b_t0=float(0.5 * np.polyval(d2, t0)) if np.isfinite(t0) else float("nan"),
                xmin=float(xx.min()), xmax=float(xx.max()))


def secim(P, pencere, t0=None):
    """158_analiz.secim ile BİREBİR aynı bant seçimi (τ̄ üzerinden)."""
    tb, x = P["tau"], P["x"]
    if pencere in ("A", "C"):
        m = (tb >= 0.43 - 1e-9) & (tb <= 0.61 + 1e-9)
    elif pencere == "D":
        m = (tb >= 0.43 - 1e-9) & (tb <= 0.59 + 1e-9)
    elif pencere == "B":
        m = (x - t0 >= -0.075) & (x - t0 <= 0.105)
    else:
        raise ValueError(pencere)
    return m


PENCERE_DERECE = {"A": 2, "D": 2, "B": 2, "C": 3}


def olc(P, pencere="A", agirlik="phi", capa=0.5, njack=8):
    """Bir (gaz, taban) koşusundan δ üçlüsü + φ fiti + jackknife.

    agirlik: 'phi' → w = 1/σ_φ (158/159 konvansiyonu, BİRİNCİL)
             'delta' → w = 1/σ_δ (δ'nın kendi jackknife'ı, sistematik)
    """
    derece = PENCERE_DERECE[pencere]
    ee_all = P["s_phi"] if agirlik == "phi" else P["s_d"]
    if pencere == "B":
        ilk = olc(P, "A", agirlik, capa, njack=0)
        if ilk is None:
            return None
        t0k = ilk["fphi"]["t0"]
        m = None
        for _ in range(4):
            m = secim(P, "B", t0k)
            f = _fit(P["x"][m], P["phi"][m], ee_all[m], derece, capa)
            if f is None:
                return None
            if abs(f["t0"] - t0k) < 1e-9:
                break
            t0k = f["t0"]
    else:
        m = secim(P, pencere)
    x, ee = P["x"][m], ee_all[m]
    fd = _fit(x, P["d"][m], ee, derece, capa)          # δ fiti  (161)
    fp = _fit(x, P["phi"][m], ee, derece, capa)        # φ fiti  (158)
    if fd is None or fp is None:
        return None
    a_half = FOUR_PI + fd["v1"]
    a_t0 = FOUR_PI + float(np.polyval(np.polyder(fd["c"]), fp["t0"]))
    out = dict(pencere=pencere, derece=derece, agirlik=agirlik, capa=capa,
               n=fd["n"], mask=m, fd=fd, fphi=fp,
               d_half=fd["v0"], dd=fd["v1"], b=fd["b"],
               a_half=a_half, a_t0=a_t0,
               t0_yildiz=0.5 - fd["v0"] / a_half if a_half != 0 else float("nan"),
               t0_phi=fp["t0"], a_phi=fp["v1_t0"], b_phi=fp["b_t0"],
               chi2=fd["chi2"], chi2dof=fd["chi2"] / max(fd["n"] - derece - 1, 1),
               xmin=fd["xmin"], xmax=fd["xmax"],
               kenar=fp["t0"] - fd["xmin"],       # <0 ⇒ τ₀ pencere DIŞINDA
               kenar_half=0.5 - fd["xmin"])       # <0 ⇒ ½ pencere DIŞINDA
    # ---- jackknife (158'in mimarisi: replika başına fit tekrarı)
    rep = {k: [] for k in ("d_half", "dd", "b", "t0_yildiz", "a_t0")}
    idx = np.where(m)[0]
    for k in range(min(njack, P["njk"])):
        xk, yk = P["TJ"][idx, k], P["DJ"][idx, k]
        pk = P["PJ"][idx, k]
        g = _fit(xk, yk, ee, derece, capa)
        gp = _fit(xk, pk, ee, derece, capa)
        if g is None or gp is None:
            continue
        ah = FOUR_PI + g["v1"]
        rep["d_half"].append(g["v0"]); rep["dd"].append(g["v1"])
        rep["b"].append(g["b"]); rep["t0_yildiz"].append(0.5 - g["v0"] / ah)
        rep["a_t0"].append(FOUR_PI + float(np.polyval(np.polyder(g["c"]), gp["t0"])))
    for q, vv in rep.items():
        v = np.array([z for z in vv if np.isfinite(z)], float)
        out["e_" + q] = (float(np.sqrt((len(v) - 1) / len(v)
                                       * np.sum((v - v.mean())**2)))
                         if len(v) >= 4 else float("nan"))
    return out


def ozet(vals):
    v = np.array([x for x in vals if np.isfinite(x)], float)
    if len(v) == 0:
        return float("nan"), float("nan"), float("nan")
    return (float(v.mean()),
            float(v.std(ddof=1)) if len(v) > 1 else 0.0,
            float(v.max() - v.min()))
