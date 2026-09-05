# -*- coding: utf-8 -*-
"""
178b — HESAP B: VEKİLLERİN θ_bw'si (SAF PYTHON, numpy YOK)
==========================================================
ÖN-MÜHÜR — bu betik 178a'nın DONMUŞ ön-kaydındaki formülün BİREBİR
UYGULAMASIDIR. Hiçbir formül, hiçbir eşik, hiçbir ağırlık burada
türetilmez; hepsi ONKAYIT.json'dan SAYI olarak okunur.

    ONKAYIT: .../scratchpad/178/ONKAYIT.json
             sha256 f81dac6d8c3acebfc61d49d698ab2bed92ead407c479531489b
                    40db7c6b7a1f6

Uygulanan (ön-kayıt "kestirici.formul" alanı, harf harf):
    theta_bw := M / (g_E * g_X_bw**2)
    g_X_bw   := SUM_b W_b * g_X_b
    g_X_b    := <X_b . x1> / <X_b . X_b>
    X_b      := sentez(s, w[lo<tau<=hi], y[lo<tau<=hi]) - mean
    x1       := Y.Xtil0
    M        := KALIB_u2(lo=0.60)   (DEĞİŞMEDİ)
    g_E      := artik.gE            (DEĞİŞMEDİ)
W_b, eşik_a, eşik_b, kilit_θ, f_θ, bant listesi: ONKAYIT'tan okunur.

══════════════════════════════════════════════════════════════════
NEDEN "HESAP B" — BAĞIMSIZ İKİNCİ YOL
══════════════════════════════════════════════════════════════════
178a çapaları numpy zinciriyle (165_cekirdek.sentez, np.dot) ölçtü.
178b AYNI cebiri saf Python'da, elle yazılmış döngülerle ve
math.fsum (tam toplam) ile kurar; numpy hiç import edilmez. Girdiler
DİSKTEKİ ÖNBELLEKLERDEN okunur — SIFIR YENİ GAZ KOŞUSU:

  s   = mid[1:]                     ← 155/eta_<ad>_t0.4_c4000.npz
  L                                 ← aynı dosya
  x1  = Xtil0 = X̃ − ⟨X̃⟩,
        X̃_n = (mid[n+1] − mid[n])·L/2π − 1     (159_cekirdek satır 127-128)
  w, y(=y_olc)                      ← 165/tayf_<ad>_t0.4_tm0.95.npz
        (Model165.__init__'in yazdığı tayf önbelleği; kaynak="olculen")
  tau = w / L                       (163_cekirdek.merdiven satır 118)
  M, g_E                            ← 176/G_<ad>.json  ("KAL", "gE")
                                      (çapalarda 167/C_<ad>.json)
  W_b, eşikler                      ← 178/ONKAYIT.json (DONMUŞ)

Vekil z verisi 155/ ile 176/ arasında BAYT BAYT ÖZDEŞTİR (cmp ile
denetlendi); eta önbelleği o z'den deterministik olarak üretilmiştir.

══════════════════════════════════════════════════════════════════
K-KİMLİK^B KAPISI (bu betiğin kendi kapısı)
══════════════════════════════════════════════════════════════════
Vekillere BAKILMADAN ÖNCE üç çapa (Hkeskin, son, HA4) bu saf-Python
makineyle yeniden ölçülür ve ONKAYIT'ın bant defteriyle (gX_b, pay,
payda, jk, s_b) karşılaştırılır. Tutmazsa hesap DURUR — vekil sayısı
üretilmez. Eşik: |Δ g_X,b| ≤ 1e-9 (saf-Python/numpy toplam sırası
farkının üstünde, hükmün duyarlılığının çok altında).

══════════════════════════════════════════════════════════════════
KURAL (ONKAYIT "kural" bloğu — DEĞİŞTİRİLMEDEN)
══════════════════════════════════════════════════════════════════
y_i = θ_bw(VF_i), i=1..4 (DÖRDÜ DE girer; eleme/ağırlıklama YASAK)
ȳ = ort(y), SAÇ_4 = 2·std(y,ddof=1)/√4 = std(y,ddof=1)
 (a) ȳ ≤ eşik_a VE |ȳ−eşik_a| ≥ SAÇ_4 ⇒ H-F1b^ba YAŞADI
 (b) ȳ ≥ eşik_b VE |ȳ−eşik_b| ≥ SAÇ_4 ⇒ H-F1b^ba ÖLDÜ
 aksi ⇒ HÜKÜMSÜZ ⇒ (n=4 tavanı + Ş2.6) KALICI HÜKÜMSÜZ
H2: Δ_i = log θ_bw(Hk) − log θ_bw(VF_i); ω_θ = kilit_θ/Δ̄;
    üst kenar Δ̄ ≥ kilit_θ, alt kenar Δ̄ > 0; KESİNLİK ⟺ marj ≥ SAÇ_4(Δ)

GİT'E DOKUNULMAZ. Bu betik yalnız HESAP_B.json yazar.
Kullanım: 178b_hesap_B.py   (argümansız)
"""
import array
import ast
import hashlib
import json
import math
import struct
import sys
import time
import zipfile
from pathlib import Path

SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S155, S165, S167, S176, S178 = (SCR / "155", SCR / "165", SCR / "167",
                                SCR / "176", SCR / "178")
ONKAYIT = S178 / "ONKAYIT.json"
OUT = S178 / "HESAP_B.json"
BU = Path(__file__).resolve()

TWO_PI = 2 * math.pi
TABAN, TAU_MAX = 0.4, 0.95
NJACK = 8
KIMLIK_TOL = 1e-9
CAPALAR = ("Hkeskin", "son", "HA4")
VEKILLER = ("VF1", "VF2", "VF3", "VF4")


# ---------------------------------------------------------------------
# 0. SAF-PYTHON .npy / .npz OKUYUCU (numpy YOK)
# ---------------------------------------------------------------------
def _npy_coz(ham):
    """Bir .npy bayt dizisini (baslik_sozlugu, veri_baytlari) olarak açar."""
    if ham[:6] != b"\x93NUMPY":
        raise SystemExit("npy sihirli sayısı yok")
    major = ham[6]
    if major == 1:
        uz = struct.unpack("<H", ham[8:10])[0]
        off = 10
    else:
        uz = struct.unpack("<I", ham[8:12])[0]
        off = 12
    bas = ast.literal_eval(ham[off:off + uz].decode("latin1"))
    if bas["fortran_order"]:
        raise SystemExit("fortran düzeni desteklenmiyor")
    return bas, ham[off + uz:]


def npz_oku(yol, ad):
    z = zipfile.ZipFile(str(yol))
    bas, veri = _npy_coz(z.read(ad + ".npy"))
    z.close()
    d = bas["descr"]
    if d == "<f8":
        a = array.array("d")
        a.frombytes(veri)
        return list(a) if bas["shape"] else a[0]
    if d == "<c16":
        a = array.array("d")
        a.frombytes(veri)
        return [a[i] for i in range(0, len(a), 2)], \
               [a[i] for i in range(1, len(a), 2)]
    if d == "<i8":
        a = array.array("q")
        a.frombytes(veri)
        return list(a) if bas["shape"] else a[0]
    raise SystemExit("bilinmeyen dtype: %s" % d)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


# ---------------------------------------------------------------------
# 1. GAZ ALANLARI — ÖNBELLEKTEN, SIFIR YENİ KOŞU
# ---------------------------------------------------------------------
def gaz_alanlari(ad):
    """s = mid[1:], x1 = Xtil0, w, y(reel/sanal), tau = w/L."""
    eta_yol = S155 / ("eta_%s_t%s_c4000.npz" % (ad, TABAN))
    tayf_yol = S165 / ("tayf_%s_t%s_tm%s.npz" % (ad, TABAN, TAU_MAX))
    L = npz_oku(eta_yol, "L")
    mid = npz_oku(eta_yol, "mid")
    N = len(mid) - 1
    s = array.array("d", mid[1:])
    kat = L / TWO_PI
    xt = array.array("d", bytes(8 * N))
    for n in range(N):
        xt[n] = (mid[n + 1] - mid[n]) * kat - 1.0
    ort = math.fsum(xt) / N
    x1 = array.array("d", bytes(8 * N))
    for n in range(N):
        x1[n] = xt[n] - ort
    w = npz_oku(tayf_yol, "w")
    yr, yi = npz_oku(tayf_yol, "y")
    tau = [v / L for v in w]
    return dict(ad=ad, L=L, N=N, s=s, x1=x1, w=w, yr=yr, yi=yi, tau=tau,
                Sx1=math.fsum(x1), eta_yol=str(eta_yol),
                tayf_yol=str(tayf_yol))


def sentez_ekle(out, s, wq, ar, ai):
    """out[n] += Re[amp·e^{i·wq·s_n}] — 165_cekirdek.sentez'in tek-çizgi hâli."""
    cos, sin = math.cos, math.sin
    for n in range(len(s)):
        a = wq * s[n]
        out[n] += cos(a) * ar - sin(a) * ai


def _jk_hata(v):
    """166_T1._jk biçimi: sqrt((n−1)/n · Σ(v−v̄)²)."""
    g = [x for x in v if x == x and abs(x) != float("inf")]
    n = len(g)
    if n < 4:
        return float("nan")
    m = math.fsum(g) / n
    return math.sqrt((n - 1) / n * math.fsum((x - m) ** 2 for x in g))


def bant_defteri(G, bant_lo, genislik):
    """Her bant için g_X,b, pay, payda, 8 jackknife eğimi ve s_b."""
    s, x1, N, Sx1 = G["s"], G["x1"], G["N"], G["Sx1"]
    w, yr, yi, tau = G["w"], G["yr"], G["yi"], G["tau"]
    defter = []
    for lo in bant_lo:
        hi = lo + genislik
        idx = [i for i in range(len(tau)) if tau[i] > lo and tau[i] <= hi + 1e-12]
        # q artan sıra = merdiven sırası (w = log q monoton) ⇒ idx zaten artan
        grup = [i % NJACK for i in range(len(idx))]
        t0 = time.time()
        Xg = [array.array("d", bytes(8 * N)) for _ in range(NJACK)]
        for k, i in enumerate(idx):
            sentez_ekle(Xg[grup[k]], s, w[i], yr[i], yi[i])
        Xb = array.array("d", bytes(8 * N))
        for j in range(NJACK):
            gj = Xg[j]
            for n in range(N):
                Xb[n] += gj[n]
        SXb = math.fsum(Xb)
        c = SXb / N
        A = math.fsum(a * b for a, b in zip(Xb, x1))
        Q = math.fsum(a * a for a in Xb)
        pay = A - c * Sx1
        payda = Q - N * c * c
        g_b = pay / payda if payda > 0 else float("nan")
        jk = []
        for j in range(NJACK):
            gj = Xg[j]
            sd = math.fsum(a - b for a, b in zip(Xb, gj))
            sdx = math.fsum((a - b) * u for a, b, u in zip(Xb, gj, x1))
            sdd = math.fsum((a - b) * (a - b) for a, b in zip(Xb, gj))
            cj = sd / N
            pj = sdx - cj * Sx1
            qj = sdd - N * cj * cj
            jk.append(pj / qj if qj > 0 else float("nan"))
        defter.append(dict(lo=float(lo), hi=float(hi), n=len(idx),
                           gX_b=g_b, pay=pay, payda=payda, jk=jk,
                           s_b=_jk_hata(jk), sure_s=time.time() - t0))
        print("    bant tau∈(%.2f,%.17g]  n=%3d  pay=%.10f  payda=%.10f  "
              "g_X,b=%.12f  s_b=%.3e  (%.0f s)"
              % (lo, hi, len(idx), pay, payda, g_b, defter[-1]["s_b"],
                 defter[-1]["sure_s"]))
        sys.stdout.flush()
    return defter


# ---------------------------------------------------------------------
# 2. ARİTMETİK YARDIMCILARI (saf Python)
# ---------------------------------------------------------------------
def ortalama(v):
    return math.fsum(v) / len(v)


def std1(v):
    n = len(v)
    m = ortalama(v)
    return math.sqrt(math.fsum((x - m) ** 2 for x in v) / (n - 1))


def sac_n(v):
    return 2.0 * std1(v) / math.sqrt(len(v))


def birlestir(defter, W):
    t = 0.0
    for k in range(len(W)):
        t += W[k] * defter[k]["gX_b"]
    return t


# ---------------------------------------------------------------------
def main():
    t_bas = time.time()
    print("=" * 74)
    print("178b — HESAP B: vekillerin θ_bw'si (SAF PYTHON, numpy YOK)")
    print("=" * 74)
    if "numpy" in sys.modules:
        raise SystemExit("numpy yüklenmiş — saf-Python koşulu bozuldu")
    print("  numpy sys.modules'te YOK  ✓")

    OK = json.load(open(ONKAYIT))
    print("  ÖN-KAYIT okundu : %s" % ONKAYIT)
    print("    zaman         : %s" % OK["zaman"])
    print("    sha256(betik) : %s" % OK["sha256"])
    print("    sha256(dosya) : %s" % sha(ONKAYIT))
    print("    hipotez       : %s" % OK["hipotez"])
    print("    bakış sayacı  : %d" % OK["bakis_sayaci"])

    BANT_LO = OK["kestirici"]["bant_kumesi"]["lo"]
    GEN = OK["kestirici"]["bant_kumesi"]["genislik"]
    W = OK["kestirici"]["agirlik"]["W"]
    esik_a = OK["esikler"]["esik_a"]
    esik_b = OK["esikler"]["esik_b"]
    kilit = OK["esikler"]["kilit_theta"]
    f_th = OK["esikler"]["f_theta"]
    th_Hk = OK["capalar"]["theta_bw"]["Hkeskin"]
    print("\n  DONMUŞ SAYILAR (ön-kayıttan, türetilmedi):")
    print("    bant lo listesi = %s   genişlik = %.2f" % (BANT_LO, GEN))
    for k in range(len(W)):
        print("    W_%d (lo=%.2f) = %.18f" % (k + 1, BANT_LO[k], W[k]))
    print("    Σ W_b        = %.18f" % math.fsum(W))
    print("    eşik_a       = %.18f" % esik_a)
    print("    eşik_b       = %.18f" % esik_b)
    print("    kilit_θ      = %.18f" % kilit)
    print("    f_θ          = %.18f" % f_th)
    print("    θ_bw(Hkeskin)= %.18f" % th_Hk)
    print("    P1 (D>0) = %s   P2 (kilit>0) = %s"
          % (OK["onkosullar"]["P1_D_pozitif"],
             OK["onkosullar"]["P2_kilit_pozitif"]))

    # ---- K-KİMLİK^B: çapaları saf-Python makineyle yeniden ölç --------
    print("\n" + "-" * 74)
    print("  K-KİMLİK^B KAPISI — çapalar saf-Python ile yeniden ölçülüyor")
    print("  (vekillere BAKILMADAN önce; tol |Δg_X,b| ≤ %.0e)" % KIMLIK_TOL)
    print("-" * 74)
    kimlik, kimlik_ok, en_kotu = {}, True, 0.0
    for ad in CAPALAR:
        print("  çapa %s:" % ad)
        G = gaz_alanlari(ad)
        print("    L=%.15f  N=%d  çizgi=%d  Σx1=%.3e"
              % (G["L"], G["N"], len(G["w"]), G["Sx1"]))
        D = bant_defteri(G, BANT_LO, GEN)
        ref = OK["capalar"]["bant_defteri"][ad]
        dg = [abs(D[k]["gX_b"] - ref[k]["gX_b"]) for k in range(len(D))]
        dp = [abs(D[k]["pay"] - ref[k]["pay"]) for k in range(len(D))]
        dq = [abs(D[k]["payda"] - ref[k]["payda"]) for k in range(len(D))]
        ds = [abs(D[k]["s_b"] - ref[k]["s_b"]) for k in range(len(D))]
        dn = [D[k]["n"] - ref[k]["n"] for k in range(len(D))]
        ok = max(dg) <= KIMLIK_TOL and max(dn) == 0 and min(dn) == 0
        kimlik_ok &= ok
        en_kotu = max(en_kotu, max(dg))
        kimlik[ad] = dict(max_dg=max(dg), max_dpay=max(dp), max_dpayda=max(dq),
                          max_ds_b=max(ds), n_farki=dn, gecti=bool(ok))
        print("    ONKAYIT ile fark: max|Δg_X,b|=%.3e  max|Δpay|=%.3e  "
              "max|Δpayda|=%.3e  max|Δs_b|=%.3e  n farkı=%s  %s"
              % (max(dg), max(dp), max(dq), max(ds), dn,
                 "✓" if ok else "✗ KAPI TUTMADI"))
        # yeni-para çapa yeniden üretimi (denetim)
        C = json.load(open(S167 / ("C_%s.json" % ad)))
        orta = [b for b in C["bant"] if abs(b["lo"] - 0.60) < 1e-9][0]
        gbw = birlestir(D, W)
        th = orta["KALIB_u2"] / (C["artik"]["gE"] * gbw * gbw)
        kimlik[ad].update(gX_bw_B=gbw, gX_bw_onkayit=OK["capalar"]["gX_bw"][ad],
                          theta_bw_B=th,
                          theta_bw_onkayit=OK["capalar"]["theta_bw"][ad],
                          d_theta=abs(th - OK["capalar"]["theta_bw"][ad]))
        print("    g_X,bw(B)=%.15f  (ön-kayıt %.15f)  Δ=%.3e"
              % (gbw, OK["capalar"]["gX_bw"][ad],
                 abs(gbw - OK["capalar"]["gX_bw"][ad])))
        print("    θ_bw(B)  =%.15f  (ön-kayıt %.15f)  Δ=%.3e"
              % (th, OK["capalar"]["theta_bw"][ad], kimlik[ad]["d_theta"]))
        sys.stdout.flush()
    if not kimlik_ok:
        raise SystemExit("ENGELLENDİ — K-KİMLİK^B tutmadı; saf-Python "
                         "makine 178a'nın makinesi değil. Vekile BAKILMADI.")
    print("  K-KİMLİK^B ✓  (en kötü |Δg_X,b| = %.3e)" % en_kotu)

    # ---- VEKİLLER ----------------------------------------------------
    print("\n" + "-" * 74)
    print("  VEKİLLER — θ_bw ölçümü (dördü de girer; eleme YASAK)")
    print("-" * 74)
    V = {}
    for ad in VEKILLER:
        print("  vekil %s:" % ad)
        G = gaz_alanlari(ad)
        print("    L=%.15f  N=%d  çizgi=%d  Σx1=%.3e"
              % (G["L"], G["N"], len(G["w"]), G["Sx1"]))
        D = bant_defteri(G, BANT_LO, GEN)
        kapi = all(b["payda"] > 0 and b["gX_b"] == b["gX_b"] for b in D)
        Gj = json.load(open(S176 / ("G_%s.json" % ad)))
        M, gE = Gj["KAL"], Gj["gE"]
        gbw = 0.0
        print("    bant-ağırlıklı birleşim (Σ W_b·g_X,b), terim terim:")
        for k in range(len(W)):
            terim = W[k] * D[k]["gX_b"]
            gbw += terim
            print("      lo=%.2f  W_b=%.15f × g_X,b=%.15f = %.15f  "
                  "(kısmi toplam %.15f)"
                  % (BANT_LO[k], W[k], D[k]["gX_b"], terim, gbw))
        th = M / (gE * gbw * gbw)
        kal = abs(th - M / (gE * gbw ** 2))
        print("    M   (KALIB_u2 lo=0.60, 176/G_%s.json) = %.18f" % (ad, M))
        print("    g_E (artik.gE, 176/G_%s.json)         = %.18f" % (ad, gE))
        print("    g_X,bw                                = %.18f" % gbw)
        print("    g_X (eski, global, donmuş)            = %.18f" % Gj["gX"])
        print("    g_E·g_X,bw²                           = %.18f"
              % (gE * gbw * gbw))
        print("    θ_bw = M/(g_E·g_X,bw²)                = %.18f" % th)
        print("    θ    (eski para, 176/G_%s.json)       = %.18f"
              % (ad, Gj["th"]))
        print("    özdeşlik kalıntısı                    = %.3e" % kal)
        print("    kapı (her bantta payda>0, θ_bw tanımlı): %s"
              % ("✓" if kapi else "✗ KAPI ÖLÜMÜ"))
        print("    R_bant_min (raporlanır, hükme dayanak DEĞİL) = %.6f"
              % Gj["R_bant_min"])
        V[ad] = dict(bant=D, M=M, gE=gE, gX_bw=gbw, theta_bw=th,
                     theta_eski=Gj["th"], gX_eski=Gj["gX"],
                     ozdeslik_kalinti=kal, kapi=bool(kapi),
                     R_bant_min=Gj["R_bant_min"])
        sys.stdout.flush()

    kapi_hepsi = all(V[a]["kapi"] for a in VEKILLER)
    print("\n  KAPI ÖLÜMÜ denetimi (dört vekil): %s"
          % ("✓ hiçbir vekilde payda ≤ 0 yok" if kapi_hepsi else "✗"))

    # ---- H1 ----------------------------------------------------------
    print("\n" + "=" * 74)
    print("  H1 — DONMUŞ KARAR KURALI (ONKAYIT 'kural' bloğu)")
    print("=" * 74)
    y = [V[a]["theta_bw"] for a in VEKILLER]
    for k, a in enumerate(VEKILLER):
        print("    θ_bw(%s) = %.18f" % (a, y[k]))
    ybar = ortalama(y)
    sdd = std1(y)
    S4 = sac_n(y)
    kat = max(y) - min(y)
    print("    ortalama ȳ            = %.18f" % ybar)
    print("    std(y, ddof=1)        = %.18f" % sdd)
    print("    SAÇ_4 = 2s/√4 = s     = %.18f" % S4)
    print("    SAÇ^kat = max−min     = %.18f  (BAĞLAYICI DEĞİL)" % kat)
    print("    eşik_a                = %.18f" % esik_a)
    print("    eşik_b                = %.18f" % esik_b)
    marj_a = abs(ybar - esik_a)
    marj_b = abs(ybar - esik_b)
    taraf_a = ybar <= esik_a
    taraf_b = ybar >= esik_b
    kesin_a = marj_a >= S4
    kesin_b = marj_b >= S4
    print("    dal (a): ȳ ≤ eşik_a ? %s   |ȳ−eşik_a| = %.18f   ≥ SAÇ_4 ? %s"
          % (taraf_a, marj_a, kesin_a))
    print("    dal (b): ȳ ≥ eşik_b ? %s   |ȳ−eşik_b| = %.18f   ≥ SAÇ_4 ? %s"
          % (taraf_b, marj_b, kesin_b))
    dal_a = bool(taraf_a and kesin_a)
    dal_b = bool(taraf_b and kesin_b)
    if dal_a:
        H1 = "YASADI"
    elif dal_b:
        H1 = "OLDU"
    else:
        H1 = "HUKUMSUZ"
    print("    dal (a) sonucu: %s" % ("SAĞLANDI ⇒ YAŞADI" if dal_a
                                      else "SAĞLANMADI"))
    print("    dal (b) sonucu: %s" % ("SAĞLANDI ⇒ ÖLDÜ" if dal_b
                                      else "SAĞLANMADI"))
    print("    H1 = %s" % H1)
    if H1 == "HUKUMSUZ":
        print("    ⇒ n=4 TAVAN + Ş2.6: KALICI HÜKÜMSÜZ (2. para biriminde "
              "de); 3. kestirici DENENMEZ")

    # ---- H2 ----------------------------------------------------------
    print("\n" + "=" * 74)
    print("  H2 — F9-θ satırı (ω_θ bandı)")
    print("=" * 74)
    P2 = OK["onkosullar"]["P2_kilit_pozitif"]
    print("    P2 (kilit_θ > 0) = %s" % P2)
    lg_Hk = math.log(th_Hk)
    dlt = [lg_Hk - math.log(V[a]["theta_bw"]) for a in VEKILLER]
    for k, a in enumerate(VEKILLER):
        w_i = kilit / dlt[k] if dlt[k] != 0 else float("nan")
        print("    Δ_%s = log θ_bw(Hk) − log θ_bw(%s) = %+.18f   "
              "(tohum ω = %+.12f)" % (a, a, dlt[k], w_i))
    dbar = ortalama(dlt)
    S4d = sac_n(dlt)
    katd = max(dlt) - min(dlt)
    om = kilit / dbar if dbar != 0 else float("nan")
    print("    Δ̄                     = %+.18f" % dbar)
    print("    std(Δ, ddof=1)        = %.18f" % std1(dlt))
    print("    SAÇ_4(Δ) = 2s/√4 = s  = %.18f" % S4d)
    print("    SAÇ^kat(Δ) = max−min  = %.18f  (BAĞLAYICI DEĞİL)" % katd)
    print("    kilit_θ               = %+.18f" % kilit)
    print("    ω_θ = kilit_θ / Δ̄     = %+.18f" % om)
    m_ust = abs(dbar - kilit)
    m_alt = abs(dbar - 0.0)
    ust_saglandi = dbar >= kilit
    alt_saglandi = dbar > 0.0
    ust_kesin = m_ust >= S4d
    alt_kesin = m_alt >= S4d
    print("    üst kenar: Δ̄ ≥ kilit_θ ? %s   marj = %.18f   KESİN ? %s"
          % (ust_saglandi, m_ust, ust_kesin))
    print("    alt kenar: Δ̄ > 0        ? %s   marj = %.18f   KESİN ? %s"
          % (alt_saglandi, m_alt, alt_kesin))
    if not P2:
        H2 = "HUKUMSUZ"
    elif ust_saglandi and alt_saglandi and ust_kesin and alt_kesin:
        H2 = "ODENEBILIR"
    elif ((not ust_saglandi) and ust_kesin) or \
         ((not alt_saglandi) and alt_kesin):
        H2 = "ODENEMEZ"
    else:
        H2 = "HUKUMSUZ"
    print("    ω_θ ∈ (0,1] bandında mı ? %s"
          % ("EVET" if (0.0 < om <= 1.0) else "HAYIR — BANDIN DIŞINDA"))
    print("    H2 = %s" % H2)

    # ---- bileşim -----------------------------------------------------
    anahtar = "%s|%s" % (H1, H2)
    cumle = OK["bilesim_tablosu"][anahtar]["cumle"]
    muhur = OK["bilesim_tablosu"][anahtar]["muhur_adayi"]
    print("\n" + "=" * 74)
    print("  BİLEŞİM (ön-kayıtta ÖNCEDEN yazılmış cümle — seçilmedi)")
    print("=" * 74)
    print("    anahtar     : %s" % anahtar)
    print("    cümle       : %s" % cumle)
    print("    mühür adayı : %s" % muhur)

    # ---- kayıt -------------------------------------------------------
    rec = dict(
        zaman=time.strftime("%Y-%m-%d %H:%M:%S %z"),
        betik=BU.name, betik_yolu=str(BU), sha256=sha(BU),
        onkayit_yolu=str(ONKAYIT), onkayit_sha256=sha(ONKAYIT),
        onkayit_betik_sha256=OK["sha256"],
        beyan="178a ön-kaydındaki formülün BİREBİR uygulaması; saf Python "
              "(numpy YOK), math.fsum ile tam toplam. Sıfır yeni gaz "
              "koşusu: girdiler 155/eta, 165/tayf, 176/G_*.json "
              "önbelleklerinden. Hiçbir eşik/ağırlık burada türetilmedi.",
        yol="B (saf Python, döngü)",
        kimlik_kapisi=dict(tol=KIMLIK_TOL, gecti=bool(kimlik_ok),
                           en_kotu_dg=en_kotu, ayrinti=kimlik),
        vekil={a: dict(theta_bw=V[a]["theta_bw"], gX_bw=V[a]["gX_bw"],
                       M=V[a]["M"], gE=V[a]["gE"],
                       theta_eski=V[a]["theta_eski"],
                       gX_eski=V[a]["gX_eski"],
                       ozdeslik_kalinti=V[a]["ozdeslik_kalinti"],
                       kapi=V[a]["kapi"], R_bant_min=V[a]["R_bant_min"],
                       bant=[{k: b[k] for k in
                              ("lo", "hi", "n", "gX_b", "pay", "payda",
                               "jk", "s_b")} for b in V[a]["bant"]])
               for a in VEKILLER},
        H1=dict(theta=y, ort=ybar, std_ddof1=sdd, sac_4=S4, sac_kat=kat,
                esik_a=esik_a, esik_b=esik_b, marj_a=marj_a, marj_b=marj_b,
                taraf_a=bool(taraf_a), taraf_b=bool(taraf_b),
                kesin_a=bool(kesin_a), kesin_b=bool(kesin_b),
                dal_a=dal_a, dal_b=dal_b, hukum=H1,
                kalici=bool(H1 == "HUKUMSUZ")),
        H2=dict(delta=dlt, delta_ort=dbar, sac_4_delta=S4d,
                sac_kat_delta=katd, kilit_theta=kilit, omega_theta=om,
                omega_tohum=[kilit / d for d in dlt],
                band_icinde=bool(0.0 < om <= 1.0),
                ust_saglandi=bool(ust_saglandi), ust_kesin=bool(ust_kesin),
                alt_saglandi=bool(alt_saglandi), alt_kesin=bool(alt_kesin),
                hukum=H2),
        bilesim=dict(anahtar=anahtar, cumle=cumle, muhur_adayi=muhur),
        kapi_olumu=dict(hepsi_saglam=bool(kapi_hepsi)),
        kural_kaynagi=OK["kural"],
        git="GİT'E DOKUNULMADI",
        sure_s=time.time() - t_bas)
    OUT.write_text(json.dumps(rec, indent=1, ensure_ascii=False))
    print("\n  -> %s" % OUT)
    print("  toplam süre: %.0f s" % rec["sure_s"])
    return rec


if __name__ == "__main__":
    main()
