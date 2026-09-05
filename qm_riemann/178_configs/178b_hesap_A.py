# -*- coding: utf-8 -*-
"""
178b — HESAP A: DÖRT VEKİLİN θ_bw'si ve DONMUŞ KURALIN UYGULANMASI
==================================================================
ÖN-MÜHÜR — ÖN-KAYIT FORMÜLÜNÜN BİREBİR UYGULAMASI
--------------------------------------------------
Bu betik HİÇBİR YENİ KARAR VERMEZ. Yaptığı tek şey, 178a'nın
mühürlediği ön-kaydı (SCRATCHPAD/178/ONKAYIT.json) diskten okumak ve
oradaki formülleri VF1..VF4 vekillerine harfiyen uygulamaktır:

  · KESTİRİCİ (ONKAYIT["kestirici"]["formul"], birebir):
        X_b     := sentez(s, w[lo<τ≤hi], y[lo<τ≤hi]) ,  X_b ← X_b − ⟨X_b⟩
        g_X,b   := ⟨X_b·x1⟩ / ⟨X_b·X_b⟩ ,   x1 := Y.Xtil0
        g_X,bw  := Σ_b W_b · g_X,b
        θ_bw    := M / ( g_E · g_X,bw² )
    M := KALİB_u2(lo=0.60)  ve  g_E := artik.gE  DEĞİŞMEZ.
  · AĞIRLIKLAR W_b: ÖN-KAYITTAN SAYI OLARAK OKUNUR. Bu betikte
    YENİDEN HESAPLANMAZ; vekilin hiçbir sayısı ağırlığa giremez
    (Ş3.2). Vekillerin jackknife'ı yalnız TANI olarak basılır.
  · BANT KÜMESİ: ONKAYIT["kestirici"]["bant_kumesi"]["lo"] =
    (0.52,0.56,0.60,0.64,0.68), genişlik 0.04, üyelik lo < τ_q ≤ hi.
  · EŞİKLER: ONKAYIT["esikler"]["esik_a"], ["esik_b"],
    ["kilit_theta"] — DİSKTEN, yeniden türetilmez (yalnız sağlama
    amacıyla θ_bw(Hk)²/θ_bw(son) ile karşılaştırılır).
  · ÖLÇÜM ZİNCİRİ: 165_cekirdek.gaz(ad,0.40,4000) → Model165(ad,
    0.40,4000,0.95,'olculen',Y=Y,tau_c=0.95).sec('olculen',0.95)
    .alanlar(kmax=3) — 178a'nın çapalarda koştuğu zincirin AYNISI.
  · KİMLİK KAPISI: her vekilde |g_E^yeni − g_E^disk| ≤ 1e−12,
    |g_X^yeni − g_X^disk| ≤ 1e−12 ve θ^eski = KALİB_u2/gcal'ın
    177'nin defterindeki sayıyla ≤ 1e−12 uyuşması denetlenir.
    Tutmazsa ölçüm makinesi aynı makine değildir ⇒ ENGELLENDİ.
  · HÜKÜM KURALI (ONKAYIT["kural"], birebir):
        SAÇ_n := 2·std(y, ddof=1)/√n ,  n = 4  ⇒  SAÇ_4 = std(ddof=1)
        (a) ȳ ≤ eşik_a VE |ȳ−eşik_a| ≥ SAÇ_4 ⇒ H-F1b^ba YAŞADI
        (b) ȳ ≥ eşik_b VE |ȳ−eşik_b| ≥ SAÇ_4 ⇒ H-F1b^ba ÖLDÜ
        aksi ⇒ HÜKÜMSÜZ (n=4 tavanı + Ş2.6 ⇒ KALICI HÜKÜMSÜZ)
        H2: Δ_i = log θ_bw(Hk) − log θ_bw(VF_i) ; ω_θ = kilit_θ/Δ̄ ;
        üst kenar Δ̄ ≥ kilit_θ, alt kenar Δ̄ > 0; KESİNLİK ⟺
        marj ≥ SAÇ_4(Δ). Sonuç cümlesi bilesim_tablosu'ndan SEÇİLİR
        (yeniden yazılmaz).
  · DÖRT TOHUM DA GİRER: eleme / ağırlıklama / dışlama YASAK.
    Bir vekilde payda ≤ 0 ya da θ_bw tanımsız kalırsa vekil
    DÜŞÜRÜLEMEZ; H1 ve H2 HÜKÜMSÜZ yazılır (Ş4.4).
  · YAZILAN TEK DOSYA: SCRATCHPAD/178/HESAP_A.json. Git'e
    DOKUNULMAZ. Eşik gevşetme YOKTUR; üçüncü kestirici DENENMEZ.

Kullanım: 178b_hesap_A.py     (argümansız)
"""
import hashlib
import importlib
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
for _p in ("167_configs", "166_configs", "165_configs", "163_configs",
           "160_configs", "159_configs", "155_configs", "154_configs"):
    sys.path.insert(0, str(QM / _p))
K = importlib.import_module("165_cekirdek")
ORT = importlib.import_module("167_ortak")
K.PENCERE.update(ORT.pencere_dict())

SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S167 = SCR / "167"
S176 = SCR / "176"
S178 = SCR / "178"
ONK = S178 / "ONKAYIT.json"
OUT = S178 / "HESAP_A.json"
BU = Path(__file__).resolve()

TABAN, CAP, TAU_MAX, TAU_C = 0.40, 4000, 0.95, 0.95
NJACK = 8
LO_M = 0.60
VEKILLER = ("VF1", "VF2", "VF3", "VF4")
KIMLIK_TOL = 1e-12


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def _jk(v):
    """166_T1._jk — birebir (yalnız TANI için)."""
    v = np.asarray([x for x in v if np.isfinite(x)], float)
    if len(v) < 4:
        return float("nan")
    return float(np.sqrt((len(v) - 1) / len(v) * np.sum((v - v.mean()) ** 2)))


def _egim(Xb, x1):
    """⟨X_b·x1⟩/⟨X_b·X_b⟩ — 178a._egim ile birebir."""
    Xb = Xb - Xb.mean()
    pay, payda = float(np.dot(Xb, x1)), float(np.dot(Xb, Xb))
    return (pay / payda if payda > 0 else float("nan")), pay, payda


def alanlar_ve_bantlar(ad, bant_lo, bant_gen):
    """178a.alanlar_ve_bantlar — körlük bekçisi dışında BİREBİR."""
    t0 = time.time()
    Y = K.gaz(ad, TABAN, CAP)
    Mo = K.Model165(ad, TABAN, CAP, TAU_MAX, "olculen", Y=Y, tau_c=TAU_C)
    Mo.sec("olculen", TAU_C).alanlar(kmax=3)
    A = Mo.artik
    x1 = Y.Xtil0
    tau, w, q = Mo.M["tau"], Mo.M["w"], Mo.M["q"]

    B = []
    for lo in bant_lo:
        hi = lo + bant_gen
        mb = Mo.msk & (tau > lo) & (tau <= hi + 1e-12)
        idx = np.flatnonzero(mb)
        idx = idx[np.argsort(q[idx], kind="stable")]        # q artan
        grup = np.arange(len(idx)) % NJACK                  # round-robin
        Xg = [K.sentez(Mo.s, w[idx[grup == j]], Mo.y[idx[grup == j]])
              for j in range(NJACK)]
        Xb = np.sum(Xg, axis=0)
        g_b, pay, payda = _egim(Xb, x1)
        jk = [_egim(Xb - Xg[j], x1)[0] for j in range(NJACK)]
        B.append(dict(lo=float(lo), hi=float(hi), n=int(len(idx)),
                      gX_b=g_b, pay=pay, payda=payda,
                      jk=[float(v) for v in jk], s_b=_jk(jk)))
    d = json.load(open(S167 / f"C_{ad}.json"))
    orta = [b for b in d["bant"] if abs(b["lo"] - LO_M) < 1e-9][0]
    return dict(ad=ad, gE=A["gE"], gX_global=A["gX"], gcal=A["gcal"],
                M=orta["KALIB_u2"], nline=A["nline"],
                gE_disk=d["artik"]["gE"], gX_disk=d["artik"]["gX"],
                gcal_disk=d["artik"]["gcal"],
                th_eski=orta["KALIB_u2"] / d["artik"]["gcal"],
                bant=B, sure_s=time.time() - t0)


def birlestir(B, W):
    return float(sum(w * b["gX_b"] for w, b in zip(W, B)))


def theta_bw(R, W):
    g = birlestir(R["bant"], W)
    return R["M"] / (R["gE"] * g * g), g


def main():
    if not ONK.exists():
        raise SystemExit("ÖN-KAYIT YOK — hesap yapılamaz.")
    P = json.load(open(ONK))

    W = list(P["kestirici"]["agirlik"]["W"])                  # DONMUŞ
    BANT_LO = tuple(P["kestirici"]["bant_kumesi"]["lo"])      # DONMUŞ
    BANT_GEN = float(P["kestirici"]["bant_kumesi"]["genislik"])
    ESIK_A = float(P["esikler"]["esik_a"])                    # DONMUŞ
    ESIK_B = float(P["esikler"]["esik_b"])                    # DONMUŞ
    KILIT = float(P["esikler"]["kilit_theta"])                # DONMUŞ
    TH_HK = float(P["capalar"]["theta_bw"]["Hkeskin"])
    TH_SON = float(P["capalar"]["theta_bw"]["son"])
    N_TAVAN = int(P["kural"]["n"])
    ESKI_TH = P["durustluk"]["eski_para"]["theta"]

    print("=" * 74)
    print("178b — HESAP A: VF1..VF4 θ_bw ; DONMUŞ KURALIN UYGULANMASI")
    print("=" * 74)
    print("  ön-kayıt : %s" % ONK)
    print("  sha256   : %s   (ön-kayıttaki: %s)"
          % (sha(P["betik_yolu"])[:16], P["sha256"][:16]))
    print("  ÖN-MÜHÜR : ön-kayıt formülünün BİREBİR uygulaması; hiçbir")
    print("             eşik/ağırlık burada türetilmez, hepsi DİSKTEN.")
    print("  W_b      = [%s]  (Σ=%.15f)"
          % (", ".join("%.6f" % v for v in W), sum(W)))
    print("  bantlar  : lo=%s  genişlik=%.2f  üyelik: lo<τ_q≤hi"
          % (list(BANT_LO), BANT_GEN))
    print("  eşik (a) = %.15f   eşik (b) = %.15f" % (ESIK_A, ESIK_B))
    print("  kilit_θ  = %+.15f   n = %d (TAVAN)\n" % (KILIT, N_TAVAN))

    # ---- ölçüm: dört vekil ------------------------------------------
    R = {}
    for v in VEKILLER:
        R[v] = alanlar_ve_bantlar(v, BANT_LO, BANT_GEN)
        r = R[v]
        print("  %-5s g_E=%.10f  g_X(global)=%.10f  M=%.10f  çizgi=%d "
              " (%.0f s)" % (v, r["gE"], r["gX_global"], r["M"],
                             r["nline"], r["sure_s"]))

    # ---- K-KİMLİK kapısı --------------------------------------------
    print("\n  K-KİMLİK KAPISI (tol %.0e):" % KIMLIK_TOL)
    kim, kimlik_ok = {}, True
    for v in VEKILLER:
        r = R[v]
        dE = abs(r["gE"] - r["gE_disk"])
        dX = abs(r["gX_global"] - r["gX_disk"])
        dT = abs(r["th_eski"] - ESKI_TH[v])
        ok = (dE <= KIMLIK_TOL) and (dX <= KIMLIK_TOL) and (dT <= 1e-12)
        kimlik_ok &= ok
        kim[v] = dict(d_gE=dE, d_gX=dX, d_theta_eski=dT, gecti=bool(ok),
                      th_eski=r["th_eski"], th_eski_177=ESKI_TH[v])
        print("    %-5s |Δg_E|=%.3e  |Δg_X|=%.3e  |Δθ_eski|=%.3e   %s"
              % (v, dE, dX, dT, "✓" if ok else "✗ KAPI TUTMADI"))
    if not kimlik_ok:
        raise SystemExit("ENGELLENDİ — K-KİMLİK kapısı tutmadı; ölçüm "
                         "makinesi aynı makine değil.")

    # ---- bant-bant tablo --------------------------------------------
    print("\n  BANT TABLOSU — g_X,b (satır: vekil, sütun: bant)")
    print("    %-6s %-5s" % ("", "n") +
          "".join("  τ∈(%.2f,%.2f]" % (lo, lo + BANT_GEN)
                  for lo in BANT_LO))
    for v in VEKILLER:
        B = R[v]["bant"]
        print("    %-6s %-5s" % (v, "") +
              "".join("  %+13.8f" % b["gX_b"] for b in B))
    print("    %-6s %-5s" % ("ikiz", "") +
          "".join("  %+13.8f" % g
                  for g in P["kestirici"]["agirlik"]["gX_b_ikiz"]))
    print("    %-6s %-5s" % ("n_çizgi", "") +
          "".join("  %13d" % b["n"] for b in R["VF1"]["bant"]))
    print("    %-6s %-5s" % ("W_b", "") +
          "".join("  %13.6f" % w for w in W))

    print("\n  BANT TABLOSU — pay ⟨X_b·x1⟩ / payda ⟨X_b·X_b⟩ / s_b(TANI)")
    for v in VEKILLER:
        for b in R[v]["bant"]:
            print("    %-5s τ∈(%.2f,%.2f]  n=%3d  pay=%12.5f  "
                  "payda=%12.5f  g_X,b=%+.8f  s_b(tanı)=%.3e"
                  % (v, b["lo"], b["hi"], b["n"], b["pay"], b["payda"],
                     b["gX_b"], b["s_b"]))

    # ---- KAPI ÖLÜMÜ denetimi (Ş4.4) ---------------------------------
    kapi_olum = []
    for v in VEKILLER:
        for b in R[v]["bant"]:
            if (not np.isfinite(b["gX_b"])) or b["payda"] <= 0:
                kapi_olum.append("%s τ∈(%.2f,%.2f]" % (v, b["lo"], b["hi"]))

    # ---- yeni para biriminde vekiller --------------------------------
    th, gbw, ind = {}, {}, {}
    for v in VEKILLER:
        th[v], gbw[v] = theta_bw(R[v], W)
        Bg = R[v]["bant"]
        ind[v] = dict(gX_global_donmus=R[v]["gX_global"],
                      gX_S_guc=sum(b["pay"] for b in Bg)
                      / sum(b["payda"] for b in Bg),
                      gX_bw=gbw[v],
                      ozdeslik_kalinti=abs(
                          th[v] - R[v]["M"] / (R[v]["gE"] * gbw[v] ** 2)))
        if not np.isfinite(th[v]):
            kapi_olum.append("%s θ_bw tanımsız" % v)

    print("\n  VEKİLLER — YENİ PARA BİRİMİ:")
    for v in VEKILLER:
        print("    %-5s g_X,bw=%.10f   θ_bw=%.10f   (eski θ=%.10f, "
              "eski g_X=%.10f)"
              % (v, gbw[v], th[v], ESKI_TH[v], R[v]["gX_global"]))
    print("\n  İNDİRGEME (Ş3.1) — g_X üç okuma:")
    for v in VEKILLER:
        i = ind[v]
        print("    %-5s global=%.8f  güç-ağırlıklı(S)=%.8f  bw=%.8f  "
              "|özdeşlik kalıntı|=%.2e"
              % (v, i["gX_global_donmus"], i["gX_S_guc"], i["gX_bw"],
                 i["ozdeslik_kalinti"]))

    # ---- DONMUŞ KURAL: H1 --------------------------------------------
    y = np.array([th[v] for v in VEKILLER], float)
    n = len(y)
    ort = float(y.mean())
    s = float(np.std(y, ddof=1))
    sac = 2.0 * s / np.sqrt(n)
    sac_kat = float(y.max() - y.min())
    marj_a = abs(ort - ESIK_A)
    marj_b = abs(ort - ESIK_B)
    taraf_a = bool(ort <= ESIK_A)
    taraf_b = bool(ort >= ESIK_B)
    kesin_a = bool(marj_a >= sac)
    kesin_b = bool(marj_b >= sac)
    yasadi = bool(taraf_a and kesin_a)
    oldu = bool(taraf_b and kesin_b)
    if kapi_olum:
        H1 = "HUKUMSUZ"
    elif yasadi:
        H1 = "YASADI"
    elif oldu:
        H1 = "OLDU"
    else:
        H1 = "HUKUMSUZ"

    print("\n" + "=" * 74)
    print("  H1 — DONMUŞ KARAR KURALI (n=%d TAVAN; dört tohum da girer)" % n)
    print("=" * 74)
    print("    θ_bw = [%s]" % ", ".join("%.10f" % v for v in y))
    print("    ȳ = %.15f   s(ddof=1) = %.15f" % (ort, s))
    print("    SAÇ_%d = 2s/√%d = %.15f   (SAÇ^kat = max−min = %.15f, "
          "BAĞLAYICI DEĞİL)" % (n, n, sac, sac_kat))
    print("    dal (a): ȳ ≤ %.15f ? %s   |ȳ−eşik_a| = %.15f  ≥ SAÇ ? %s"
          % (ESIK_A, taraf_a, marj_a, kesin_a))
    print("    dal (b): ȳ ≥ %.15f ? %s   |ȳ−eşik_b| = %.15f  ≥ SAÇ ? %s"
          % (ESIK_B, taraf_b, marj_b, kesin_b))
    print("    ⇒ H1 = %s" % H1)

    # ---- DONMUŞ KURAL: H2 --------------------------------------------
    P2 = bool(P["onkosullar"]["P2_kilit_pozitif"])
    D = np.array([np.log(TH_HK) - np.log(th[v]) for v in VEKILLER], float)
    Dort = float(D.mean())
    sD = float(np.std(D, ddof=1))
    sacD = 2.0 * sD / np.sqrt(n)
    sacD_kat = float(D.max() - D.min())
    omega = KILIT / Dort if Dort != 0 else float("inf")
    omega_i = [KILIT / d if d != 0 else float("inf") for d in D]
    ust_ok = bool(Dort >= KILIT)
    alt_ok = bool(Dort > 0.0)
    marj_ust = abs(Dort - KILIT)
    marj_alt = abs(Dort - 0.0)
    ust_kesin = bool(marj_ust >= sacD)
    alt_kesin = bool(marj_alt >= sacD)
    if (not P2) or kapi_olum:
        H2 = "HUKUMSUZ"
    elif ust_ok and alt_ok and ust_kesin and alt_kesin:
        H2 = "ODENEBILIR"
    elif ((not ust_ok) and ust_kesin) or ((not alt_ok) and alt_kesin):
        H2 = "ODENEMEZ"
    else:
        H2 = "HUKUMSUZ"

    print("\n" + "=" * 74)
    print("  H2 — F9'un θ satırı (kilit_θ = %+.15f)" % KILIT)
    print("=" * 74)
    print("    Δ_i = log θ_bw(Hk) − log θ_bw(VF_i):")
    for v, d_, w_ in zip(VEKILLER, D, omega_i):
        print("      %-5s Δ = %+.15f   ω_i = %+.10f" % (v, d_, w_))
    print("    Δ̄ = %+.15f   s(Δ) = %.15f   SAÇ_%d(Δ) = %.15f  "
          "(SAÇ^kat = %.15f)" % (Dort, sD, n, sacD, sacD_kat))
    print("    ω_θ = kilit_θ/Δ̄ = %.15f   (bant: 0 < ω ≤ 1)" % omega)
    print("    üst kenar Δ̄ ≥ kilit_θ ? %s   marj=%.15f ≥ SAÇ ? %s"
          % (ust_ok, marj_ust, ust_kesin))
    print("    alt kenar Δ̄ > 0      ? %s   marj=%.15f ≥ SAÇ ? %s"
          % (alt_ok, marj_alt, alt_kesin))
    print("    ⇒ H2 = %s" % H2)

    cumle = P["bilesim_tablosu"]["%s|%s" % (H1, H2)]["cumle"]
    muhur = P["bilesim_tablosu"]["%s|%s" % (H1, H2)]["muhur_adayi"]
    print("\n" + "=" * 74)
    print("  SONUÇ CÜMLESİ (ön-kayıttan SEÇİLDİ, yazılmadı):")
    print("    %s" % cumle)
    print("    mühür adayı: %s" % muhur)
    print("=" * 74)

    # ---- TANI: R_bant (hükme dayanak DEĞİL) --------------------------
    rb = {}
    for v in VEKILLER:
        g = json.load(open(S176 / f"G_{v}.json"))
        rb[v] = dict(R_bant_min=g.get("R_bant_min"),
                     F0_R_bant=g.get("F0_R_bant"),
                     R_bant_hukum=g.get("R_bant_hukum"))
    print("\n  TANI — R_bant (raporlanır; HİÇBİR HÜKME DAYANAK YAPILMAZ):")
    for v in VEKILLER:
        print("    %-5s R_bant_min=%.6f  F0=%s"
              % (v, rb[v]["R_bant_min"], rb[v]["F0_R_bant"]))

    # ---- kayıt -------------------------------------------------------
    tarih = subprocess.run(["date"], capture_output=True,
                           text=True).stdout.strip()
    rec = dict(
        zaman=time.strftime("%Y-%m-%d %H:%M:%S %z"), date_cikti=tarih,
        betik=BU.name, betik_yolu=str(BU), sha256=sha(BU),
        on_muhur="ön-kayıt formülünün birebir uygulaması",
        onkayit=dict(yol=str(ONK), sha256_onkayitta=P["sha256"],
                     sha256_178a_simdi=sha(P["betik_yolu"]),
                     hipotez=P["hipotez"], bakis_sayaci=P["bakis_sayaci"]),
        donmus_girdi=dict(W=W, bant_lo=list(BANT_LO), bant_gen=BANT_GEN,
                          esik_a=ESIK_A, esik_b=ESIK_B, kilit_theta=KILIT,
                          f_theta=P["esikler"]["f_theta"],
                          D_dlog=P["esikler"]["D_dlog_theta_son_Hk"],
                          E_dlog=P["esikler"]["E_dlog_theta_erfc_Hk"],
                          theta_bw_capalar=P["capalar"]["theta_bw"],
                          n=N_TAVAN),
        vekil=dict(
            theta_bw={v: th[v] for v in VEKILLER},
            gX_bw={v: gbw[v] for v in VEKILLER},
            gE={v: R[v]["gE"] for v in VEKILLER},
            M={v: R[v]["M"] for v in VEKILLER},
            gX_global={v: R[v]["gX_global"] for v in VEKILLER},
            nline={v: R[v]["nline"] for v in VEKILLER},
            bant_defteri={v: R[v]["bant"] for v in VEKILLER},
            sure_s={v: R[v]["sure_s"] for v in VEKILLER}),
        indirgeme=ind,
        kapilar=dict(K_KIMLIK=dict(tol=KIMLIK_TOL, gecti=bool(kimlik_ok),
                                   ayrinti=kim),
                     KAPI_OLUMU=dict(var=bool(kapi_olum), liste=kapi_olum)),
        H1=dict(theta=[float(v) for v in y], ort=ort, s_ddof1=s,
                sac_n=sac, sac_kat=sac_kat, n=n,
                esik_a=ESIK_A, esik_b=ESIK_B,
                taraf_a=taraf_a, kesin_a=kesin_a, marj_a=marj_a,
                taraf_b=taraf_b, kesin_b=kesin_b, marj_b=marj_b,
                hukum=H1,
                saglama_esik_a=TH_HK ** 2 / TH_SON),
        H2=dict(delta=[float(v) for v in D], delta_ort=Dort, s_ddof1=sD,
                sac_n=sacD, sac_kat=sacD_kat, kilit_theta=KILIT,
                omega_theta=float(omega),
                omega_tohum=[float(v) for v in omega_i],
                ust_ok=ust_ok, ust_kesin=ust_kesin, marj_ust=marj_ust,
                alt_ok=alt_ok, alt_kesin=alt_kesin, marj_alt=marj_alt,
                P2=P2, hukum=H2),
        sonuc=dict(anahtar="%s|%s" % (H1, H2), cumle=cumle,
                   muhur_adayi=bool(muhur)),
        tani=dict(R_bant=rb,
                  not_="R_bant ve vekil jackknife'ları TANI'dır; hükme "
                       "girmez. Ağırlıklar YALNIZ ikizden, ön-kayıttan "
                       "SAYI olarak alındı."),
        kural_metni=P["kural"],
        git="GİT'E DOKUNULMADI",
        yazilan_dosyalar=[str(OUT)])
    OUT.write_text(json.dumps(rec, indent=1, ensure_ascii=False))
    print("\n  -> %s" % OUT)
    print("  sha256(betik) = %s" % rec["sha256"])
    return rec


if __name__ == "__main__":
    main()
