"""
157 — ANALİZ: R_lad tayfı, τ₀_lad kimliği, taban çöküşü, kapalı-form yarışı
==========================================================================
Girdi : scratchpad/157/R_<veri>_t<taban>_<bantset>.json   (157_kos.py)
Çıktı : scratchpad/157/analiz_cikti.txt (ekrana da basılır)

BÖLÜMLER
  A   R_lad(τ) TAM TAYFI (iki pencere × dört taban × iki ızgara)
  A2  Ayrıntı tablosu: φ_Γ, arg M_emp, R_tam, R_lad, |ΔRe|, n_eff, ρ
  B   τ₀: φ, R_tam, R_lad sıfırları — KİMLİK sınavı + taban yasası
  C   TABAN ÇÖKÜŞÜ: R_lad tabanla ne kadar oynuyor; apsis kaydırılınca
      eğriler çöküyor mu (yani kalan taban etkisi SALT τ₀ kayması mı)?
  D   YARIŞ: ≤2 parametreli adaylar + φ-ilkel (h) ailesi
  E   Artık yapıları ve parametre hataları (Δχ²=1 profili)
"""
import json
from itertools import product
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/157")
VERI = ("son", "orta")
TABAN = (0.28, 0.34, 0.40, 0.46, 0.52, 0.58)
NEFF_ESIK = 3000.0      # 156'nın kuralı: ağırlık %1'den az etkin noktaya
                        # çökmüşse ortalama uydurma birkaç bağa dayanır
# 156'nın "R_lad taban-bağımsız" hükmü yalnız {0.46,0.52,0.58}'de sınanmıştı
TABAN_156 = (0.46, 0.52, 0.58)
TAU_UST = 0.76          # 154: τ>0.76 çözünmeyen-çizgi bölgesi → fite girmez
OUT = []


def yaz(*a):
    s = " ".join(str(x) for x in a)
    OUT.append(s)
    print(s, flush=True)


def yukle():
    D = {}
    for v, t, g in product(VERI, TABAN, ("kaba", "ince", "tau0")):
        for isim in (f"R_{v}_t{t}_{g}.json", f"R_{v}_t{t:g}_{g}.json"):
            p = SCR / isim
            if p.exists():
                D[(v, t, g)] = json.load(open(p))
                break
    return D


def bant(d):
    return [b for b in d["bantlar"] if b.get("olculdu")]


def sarmali(iz, R):
    """Kök R'ye giden yolda arg M SARIYOR mu?

    arg M_k(R) R=0'da arg M_emp'tir ve R büyüdükçe sürekli değişir; ama
    ±π'yi geçince numpy'ın arg'ı atlar. Kök bulucu bu atlamanın öbür
    yakasında sahte bir kök yakalayabilir (154'ün Denetim 3'ü). R=0'dan
    köke giden ızgara parçasında ardışık fark > 1 rad ise SARMA vardır:
    o bantta kanal ölçülen fazı gerçekte ÜRETEMİYOR.
    """
    if iz is None or not np.isfinite(R):
        return False
    Rg = np.asarray(iz["R"], float); ang = np.asarray(iz["arg"], float)
    i0 = int(np.argmin(np.abs(Rg)))
    j = int(np.argmin(np.abs(Rg - R)))
    lo, hi = (i0, j) if j >= i0 else (j, i0)
    seg = ang[lo:hi + 1]
    return bool(len(seg) > 1 and np.max(np.abs(np.diff(seg))) > 1.0)


def saglam(b, kan, neff=True):
    """Bant bu kanalda ÖLÇÜM sayılır mı? (faz tuttu, sarma yok, |Γ|
    patlak değil, ağırlık n_eff ≥ 3000 ile sağlıklı)"""
    k = b["kanal"][kan]
    return (k["artik"] < 0.02 and not b["patlak"]
            and not sarmali(k.get("izgara"), k["Rn"])
            and (not neff or k["neff"] >= NEFF_ESIK))


# ------------------------------------------------------------ kök yardımı
def kok_izgara(Rg, ang, phi):
    """arg M(R̃) = phi kökü: R=0'DAN DIŞA ilk kök (154 kuralı), doğrusal ara
    değerle. Kök yoksa None."""
    Rg = np.asarray(Rg, float); ang = np.asarray(ang, float)
    i0 = int(np.argmin(np.abs(Rg)))
    d = ang - phi
    kok = None
    for ileri in (True, False):
        rng = range(i0, len(Rg) - 1) if ileri else range(i0, 0, -1)
        aday = None
        for i in rng:
            j = i + 1 if ileri else i - 1
            if d[i] == 0.0:
                aday = float(Rg[i]); break
            if d[i] * d[j] < 0:
                aday = float(Rg[i] - d[i] * (Rg[j] - Rg[i]) / (d[j] - d[i]))
                break
        if aday is None:
            continue
        if kok is None or abs(aday) < abs(kok):
            kok = aday
    return kok


def sifir_fit(x, y, derece=2):
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if len(x) < derece + 1:
        return float("nan"), None, 0
    c = np.polyfit(x, y, derece)
    r = np.roots(c)
    r = np.array([z.real for z in r if abs(z.imag) < 1e-9])
    if len(r) == 0:
        return float("nan"), c, len(x)
    return float(r[np.argmin(np.abs(r - x.mean()))]), c, len(x)


def tau0_bir(bs, alan, ust, replik=None):
    """alan: 'phi' | 'dphi' (= φ_Γ − arg M_emp) | 'tam' | 'lad' | 'eta'

    'dphi' KANAL-BAĞIMSIZ ve TAHMİNCİ-BAĞIMSIZ tanımdır: R_k(τ)=0 tam
    olarak φ_Γ(τ) = arg M_emp(τ) noktasında olur (M_k(0)=M_emp, kanaldan
    bağımsız). Kanal başına 'tam'/'lad' fitleri ise R eğrisinin KENDİ
    eğriliğini taşır; kuadratik fitin kökü bu yüzden kanaldan kanala
    ~0.002 oynar — kimliğin ihlali değil, fit-formu farkıdır.
    """
    xs, ys = [], []
    for b in bs:
        if b["tau"] > ust + 1e-9:
            continue
        if replik is None:
            x = b["tau_eff"]
            y = (b["phi"] if alan == "phi" else
                 b["phi"] - b["argMe"] if alan == "dphi" else
                 b["kanal"][alan]["Rham"])
        else:
            x = b["teff_jk"][replik]
            if alan == "phi":
                y = b["phi_jk"][replik]
            elif alan == "dphi":
                y = b["phi_jk"][replik] - b["argMe"]
            else:
                y = (b["kanal"][alan]["R_jk"][replik]
                     * b["kanal"][alan]["olcek"])
        if np.isfinite(x) and np.isfinite(y):
            xs.append(x); ys.append(y)
    return sifir_fit(xs, ys)


def tau0_hatali(bs, alan, ust, njack=8):
    t0, c, n = tau0_bir(bs, alan, ust)
    rep = []
    for k in range(njack):
        tk, _, nk = tau0_bir(bs, alan, ust, replik=k)
        if np.isfinite(tk) and nk >= 3:
            rep.append(tk)
    if len(rep) >= 4:
        rep = np.array(rep)
        e = float(np.sqrt((len(rep) - 1) / len(rep)
                          * np.sum((rep - rep.mean())**2)))
    else:
        e = float("nan")
    return t0, e, n


# =====================================================================  A
def bolum_A(D):
    yaz("\n" + "=" * 96)
    yaz("A — R_lad(τ) TAM TAYFI   (ham R = X_lad'ın ağırlıktaki katsayısı; "
        "apsis τ_eff; ±jk)")
    yaz("=" * 96)
    for g in ("kaba", "ince"):
        for v in VERI:
            tb = [t for t in TABAN if (v, t, g) in D]
            if not tb:
                continue
            yaz(f"\n--- {g} ızgara ({'0.03' if g=='kaba' else '0.02'}) / "
                f"pencere {v} ---")
            yaz("  τ̄     τ_eff  " + "".join(
                f"|  R_lad@{t:.2f}  ±jk  " for t in tb))
            taus = sorted({b["tau"] for t in tb for b in bant(D[(v, t, g)])})
            for tau in taus:
                te, cells = float("nan"), ""
                for t in tb:
                    bb = [b for b in bant(D[(v, t, g)]) if b["tau"] == tau]
                    if not bb:
                        cells += "|      —           "
                        continue
                    b = bb[0]; k = b["kanal"]["lad"]
                    if np.isnan(te):
                        te = b["tau_eff"]
                    f = "*" if k["artik"] >= 0.02 else " "
                    cells += f"| {k['Rham']:+8.3f}{f} {k['sRham_jk']:5.3f} "
                yaz(f"  {tau:.4f} {te:.4f} {cells}")
            yaz("  (* = faz artığı ≥ 0.02 rad → R eşleşmesi TUTMADI)")


def bolum_A2(D, g="kaba", t=0.46):
    yaz("\n" + "=" * 96)
    yaz(f"A2 — AYRINTI (taban {t}, {g} ızgara): kanallar ve aşırı-belirleme")
    yaz("=" * 96)
    yaz("  pen  τ̄     τ_eff   |Γ|   φ_Γ     argMe   R_tam  ±jk   |ΔRe|t "
        " R_lad  ±jk   |ΔRe|l   n_eff_l    ρ")
    for v in VERI:
        if (v, t, g) not in D:
            continue
        for b in bant(D[(v, t, g)]):
            kt, kl = b["kanal"]["tam"], b["kanal"]["lad"]
            bay = ""
            if b["tau"] > TAU_UST:
                bay += " [τ>0.76 fit-dışı]"
            if kl["artik"] >= 0.02:
                bay += " [LAD FAZ TUTMADI]"
            if b["patlak"]:
                bay += " [|Γ|>1.5]"
            if kl["dRe"] >= 0.06:
                bay += " [ΔRe>0.06]"
            yaz(f"  {v:4s} {b['tau']:.4f} {b['tau_eff']:.4f} {b['absG']:.3f} "
                f"{b['phi']:+.4f} {b['argMe']:+.4f} {kt['Rham']:+6.3f} "
                f"{kt['sRham_jk']:5.3f} {kt['dRe']:6.4f} {kl['Rham']:+6.3f} "
                f"{kl['sRham_jk']:5.3f} {kl['dRe']:6.4f} {kl['neff']:9.0f} "
                f"{b['rho']:.4f}{bay}")


# =====================================================================  B
def bolum_B(D):
    yaz("\n" + "=" * 96)
    yaz("B — τ₀: φ_Γ, R_tam ve R_lad'ın SIFIRLARI   (KİMLİK SINAVI)")
    yaz("=" * 96)
    yaz("KİMLİK (analitik): R_k = 0 ⇒ w_k ≡ 1 ⇒ M_k(0) = ⟨e^{−iA·X_tam}⟩ =")
    yaz("M_emp, ve M_emp KANALDAN BAĞIMSIZDIR. Dolayısıyla her kanalın R'si")
    yaz("AYNI τ'da sıfırlanır: φ_Γ(τ) = arg M_emp(τ). Aşağıdaki Δ sütunu bu")
    yaz("kimliğin sayısal doğrulamasıdır (fit gürültüsü dışında 0 olmalı).")
    yaz("")
    yaz("  pen  taban  σ_η²    τ₀(φ)  ±jk     τ₀*(φ−argMe) ±jk    "
        "τ₀(R_tam) ±jk    τ₀(R_lad) ±jk    Δ(lad−tam)  n  aralık")
    kayit = {}
    for v in VERI:
        for t in TABAN:
            k = (v, t, "ince") if (v, t, "ince") in D else (v, t, "tau0")
            if k not in D:
                continue
            d = D[k]
            bs = bant(d)
            ust = 0.57 if len([b for b in bs if b["tau"] <= 0.57]) >= 5 \
                else 0.63
            a, ea, na = tau0_hatali(bs, "phi", ust)
            s_, es, _ = tau0_hatali(bs, "dphi", ust)
            b_, eb, _ = tau0_hatali(bs, "tam", ust)
            c_, ec, _ = tau0_hatali(bs, "lad", ust)
            kayit[(v, t)] = dict(seta=d["s_eta"], phi=a, ephi=ea, tam=b_,
                                 etam=eb, lad=c_, elad=ec, dphi=s_,
                                 edphi=es, n=na, ust=ust)
            yaz(f"  {v:4s} {t:.2f}  {d['s_eta']:.4f}  {a:.4f} {ea:.4f}   "
                f"  {s_:.4f} {es:.4f}      {b_:.4f} {eb:.4f}   "
                f"{c_:.4f} {ec:.4f}   {c_-b_:+.6f}  {na:2d}  τ̄≤{ust}")
    yaz("  NOT: τ₀* = (φ_Γ − arg M_emp)'in sıfırı = HER kanalın R'sinin")
    yaz("  TAM sıfır-geçişi (kimlik). Δ(lad−tam) ≠ 0 olması kimliğin")
    yaz("  ihlali değil, iki R eğrisinin KUADRATİK FİT formunun farkıdır.")

    yaz("\n  TABAN YASASI (155'in gerçek gazda ölçtüğü: τ₀ = 0.519 − "
        "0.19·σ_η²)")
    yaz("  pen  nicelik  dτ₀/dσ_η²   kesişim(σ_η²→0)  yayılım(taban ekseni)")
    for v in VERI:
        for ad in ("phi", "dphi", "tam", "lad"):
            xs = [kayit[(v, t)]["seta"] for t in TABAN if (v, t) in kayit]
            ys = [kayit[(v, t)][ad] for t in TABAN if (v, t) in kayit]
            xs = [x for x, y in zip(xs, ys) if np.isfinite(y)]
            ys = [y for y in ys if np.isfinite(y)]
            if len(xs) < 3:
                continue
            p = np.polyfit(xs, ys, 1)
            yaz(f"  {v:4s} {ad:7s} {p[0]:+8.3f}      {p[1]:.4f}"
                f"            {max(ys)-min(ys):.4f}")
    return kayit


def bolum_B2(D, taban=0.46):
    """τ₀'ın FİT KONVANSİYONUNA duyarlılığı: hangi tanım en dayanıklı?

    155, τ₀'ın taban konvansiyonuna asılı olduğunu ölçmüştü. Burada ikinci
    bir konvansiyon ekseni taranıyor: fit derecesi, fit aralığı, bant
    ızgarası. Soru: dört tanım (φ, φ−argM_emp, R_tam, R_lad) arasında
    hangisi bu eksende en az oynuyor?
    """
    yaz("\n" + "=" * 96)
    yaz(f"B2 — τ₀'ın FİT KONVANSİYONUNA duyarlılığı (taban {taban})")
    yaz("=" * 96)
    yaz("  pen  ızgara  τ̄≤   derece | τ₀(φ)   τ₀*     τ₀(R_tam) τ₀(R_lad) | n")
    kol = {a: [] for a in ("phi", "dphi", "tam", "lad")}
    for v in VERI:
        for g in ("ince", "kaba"):
            if (v, taban, g) not in D:
                continue
            bs = bant(D[(v, taban, g)])
            for ust in (0.55, 0.57, 0.60, 0.63):
                for deg in (1, 2):
                    r, n = [], 0
                    for alan in ("phi", "dphi", "tam", "lad"):
                        xs, ys = [], []
                        for b in bs:
                            if b["tau"] > ust + 1e-9:
                                continue
                            y = (b["phi"] if alan == "phi" else
                                 b["phi"] - b["argMe"] if alan == "dphi" else
                                 b["kanal"][alan]["Rham"])
                            xs.append(b["tau_eff"]); ys.append(y)
                        n = len(xs)
                        v0 = sifir_fit(xs, ys, deg)[0]
                        r.append(v0)
                        if np.isfinite(v0):
                            kol[alan].append(v0)
                    if n < deg + 2:
                        continue
                    yaz(f"  {v:4s} {g:5s} {ust:.2f}    {deg}   | "
                        + "  ".join(f"{x:.4f}" for x in r) + f"  | {n}")
    yaz("\n  KONVANSİYON YAYILIMI (taban SABİT 0.46; yalnız ızgara/aralık/"
        "derece değişiyor):")
    for alan, ad in (("phi", "τ₀(φ_Γ)"), ("dphi", "τ₀* = (φ−argM_emp) sıfırı"),
                     ("tam", "τ₀(R_tam)"), ("lad", "τ₀(R_lad)")):
        a = np.array(kol[alan], float); a = a[np.isfinite(a)]
        yaz(f"    {ad:28s} : {a.mean():.4f} ± {a.std(ddof=1):.4f}  "
            f"yayılım {a.max()-a.min():.4f}  (n={len(a)} konvansiyon)")
    yaz("  → R_lad'ın sıfırı bu eksende EN OYNAK olan; τ₀* en dayanıklı.")


# =====================================================================  C
def bolum_C(D, kayit, g="kaba"):
    yaz("\n" + "=" * 96)
    yaz("C — TABAN ÇÖKÜŞÜ: R_tam ve R_lad tabanla ne yapıyor?")
    yaz("=" * 96)
    yaz("Bant-bant tabanlar arası YAYILIM (maks−min) ve ortalamaya oranı.")
    yaz("156'nın '±%3' iddiası burada DÖRT tabanda sınanıyor.")
    for v in VERI:
        yaz(f"\n--- pencere {v} ---")
        yaz("  τ̄     τ_eff  |  R_tam ort   yayılım    %   |  R_lad ort  "
            " yayılım    %   | tabanlar")
        taus = sorted({b["tau"] for t in TABAN if (v, t, g) in D
                       for b in bant(D[(v, t, g)])})
        for tau in taus:
            vv = {"tam": [], "lad": []}
            tbs = {"tam": [], "lad": []}
            te = float("nan")
            for t in TABAN:
                if (v, t, g) not in D:
                    continue
                bb = [b for b in bant(D[(v, t, g)]) if b["tau"] == tau]
                if not bb:
                    continue
                b = bb[0]
                for kan in ("tam", "lad"):
                    if saglam(b, kan):
                        te = b["tau_eff"]
                        vv[kan].append(b["kanal"][kan]["Rham"])
                        tbs[kan].append(t)
            if len(vv["tam"]) < 2 and len(vv["lad"]) < 2:
                continue
            hc = []
            for kan in ("tam", "lad"):
                if len(vv[kan]) < 2:
                    hc.append("     —         —      — "); continue
                y_ = max(vv[kan]) - min(vv[kan]); m_ = float(np.mean(vv[kan]))
                hc.append(f"{m_:+7.3f} {y_:8.4f} {100*y_/abs(m_):6.1f}")
            yaz(f"  {tau:.4f} {te:.4f} | {hc[0]} ({len(vv['tam'])}) | "
                f"{hc[1]} ({len(vv['lad'])}) | "
                f"lad tabanları: {','.join(f'{x:.2f}' for x in tbs['lad'])}")

    # --- taban EKSENİ boyunca R: 156'nın dar penceresi vs geniş eksen ---
    yaz("\n  R(taban) SABİT τ'DA — 156'nın {0.46,0.52,0.58} penceresi vs")
    yaz("  geniş eksen {0.28…0.58}. (156'nın ±%3'ü DURAĞAN bir noktanın")
    yaz("  yerel düzlüğü mü, yoksa gerçek bir bağımsızlık mı?)")
    for v in VERI:
        taus = sorted({b["tau"] for t in TABAN if (v, t, g) in D
                       for b in bant(D[(v, t, g)])})
        taus = [x for x in taus if 0.58 < x <= TAU_UST]
        for kan in ("tam", "lad"):
            yaz(f"\n  --- {v} / R_{kan}(taban) ---")
            yaz("  taban σ_η²    σ_lad²  pay_lad " + "".join(
                f"  τ={x:.3f}" for x in taus))
            for t in TABAN:
                if (v, t, g) not in D:
                    continue
                d = D[(v, t, g)]
                row = ""
                for x in taus:
                    bb = [b for b in bant(d) if b["tau"] == x]
                    if not bb or bb[0]["kanal"][kan]["artik"] >= 0.02:
                        row += "        —"
                    else:
                        f = "" if saglam(bb[0], kan) else "⚠"
                        row += f" {bb[0]['kanal'][kan]['Rham']:+7.3f}{f:1s}"
                yaz(f"  {t:.2f}  {d['s_eta']:.4f}  {d['s_lad']:.4f}  "
                    f"{d['kov']['lad']/d['sA2']:+.4f}{row}")
            yaz("  (⚠ = n_eff < 3000: ağırlık %1'den az etkin bağa çökmüş)")
            for etiket, tbs, nf in (
                    ("156 penceresi {0.46,0.52,0.58}", TABAN_156, True),
                    ("geniş eksen {0.28…0.58}", TABAN, True),
                    ("geniş eksen, n_eff süzgeci YOK", TABAN, False)):
                sat = ""
                for x in taus:
                    vv = []
                    for t in tbs:
                        if (v, t, g) not in D:
                            continue
                        bb = [b for b in bant(D[(v, t, g)]) if b["tau"] == x]
                        if bb and saglam(bb[0], kan, nf):
                            vv.append(bb[0]["kanal"][kan]["Rham"])
                    sat += (f" {100*(max(vv)-min(vv))/abs(np.mean(vv)):6.1f}%"
                            f"/{len(vv)}"
                            if len(vv) >= 2 else "       — ")
                yaz(f"  yayılım% [{etiket:32s}]    {sat}")
            yaz("  (yayılım% / kaç taban; n_eff ≥ 3000 süzgeci ilk iki "
                "satırda AÇIK)")

    yaz("\n  ÇÖKÜŞ SINAVI — apsis (τ_eff − τ₀(taban)) yapılınca tabanlar arası")
    yaz("  yayılım küçülüyor mu? (τ₀ = aynı koşunun KENDİ R sıfırı, B'den)")
    yaz("  Sabit bir apsis ızgarasında doğrusal ara değerle; ⟨⟩ = ızgara ort.")
    yaz("  pen kanal | ham apsis (τ_eff)     | kaydırılmış (τ_eff−τ₀)")
    for v in VERI:
        for kan in ("tam", "lad"):
            eg = {}
            for t in TABAN:
                if (v, t, g) not in D or (v, t) not in kayit:
                    continue
                t0 = kayit[(v, t)][kan]
                pts = [(b["tau_eff"], b["kanal"][kan]["Rham"])
                       for b in bant(D[(v, t, g)])
                       if b["kanal"][kan]["artik"] < 0.02
                       and b["tau"] <= TAU_UST]
                if len(pts) < 4 or not np.isfinite(t0):
                    continue
                eg[t] = (np.array([p[0] for p in pts]),
                         np.array([p[1] for p in pts]), t0)
            if len(eg) < 3:
                continue

            def yayilim(kaydir):
                gl = (np.arange(0.09, 0.2401, 0.01) if kaydir
                      else np.arange(0.60, 0.7401, 0.01))
                sp, rl = [], []
                for xq in gl:
                    vv = []
                    for t, (xs, ys, t0) in eg.items():
                        xx = xs - t0 if kaydir else xs
                        if xq < xx.min() or xq > xx.max():
                            vv = []; break
                        vv.append(float(np.interp(xq, xx, ys)))
                    if len(vv) >= 3:
                        sp.append(max(vv) - min(vv))
                        rl.append((max(vv) - min(vv)) / abs(np.mean(vv)))
                return ((float(np.mean(sp)), float(np.mean(rl)))
                        if sp else (np.nan, np.nan))
            a0, r0 = yayilim(False)
            a1, r1 = yayilim(True)
            yaz(f"  {v:4s} {kan:3s}  | ⟨yayılım⟩={a0:.4f} ({100*r0:5.1f}%) "
                f"| ⟨yayılım⟩={a1:.4f} ({100*r1:5.1f}%)")

    # --- İLKEL nesnenin çöküşü: φ_Γ(τ − τ₀(taban)) tabanlar arası üst üste
    #     biniyor mu? (Taban konvansiyonu φ'yi SALT ÖTELİYOR mu?)
    yaz("\n  φ ÇÖKÜŞÜ — taban konvansiyonu ilkel φ_Γ'yı SALT ÖTELİYOR mu?")
    yaz("  φ(τ) = a(τ−τ₀) + b(τ−τ₀)² fitinin ŞEKİL katsayıları tabanla ne")
    yaz("  yapıyor? (τ₀ tabanla 0.013 kayıyor — a ve b de kayıyor mu?)")
    yaz("  pen  taban  τ₀(φ)    a       b     |  a/⟨a⟩   b/⟨b⟩")
    kuresel = []
    for v in VERI:
        kayd = []
        for t in TABAN:
            k = (v, t, "ince") if (v, t, "ince") in D else (v, t, "tau0")
            if k not in D:
                continue
            bs = [b for b in bant(D[k]) if b["tau"] <= 0.61]
            if len(bs) < 5:
                continue
            xx = np.array([b["tau_eff"] for b in bs])
            yy = np.array([b["phi"] for b in bs])
            ee = np.array([b["sPhi_jk"] for b in bs])
            ee = np.where(np.isfinite(ee) & (ee > 0), ee, np.nanmedian(ee))
            c = np.polyfit(xx, yy, 2, w=1.0 / ee)
            r = np.roots(c)
            r = np.array([z.real for z in r if abs(z.imag) < 1e-9])
            t0 = float(r[np.argmin(np.abs(r - xx.mean()))])
            kayd.append((t, t0, float(c[1] + 2 * c[0] * t0), float(c[0])))
            kuresel.append((v, t, t0, float(c[1] + 2 * c[0] * t0),
                            float(c[0])))
        if not kayd:
            continue
        am = np.mean([k[2] for k in kayd]); bm = np.mean([k[3] for k in kayd])
        for t, t0, a_, b_ in kayd:
            yaz(f"  {v:4s} {t:.2f}  {t0:.4f} {a_:7.3f} {b_:+7.3f}  |  "
                f"{a_/am:.4f}  {b_/bm:.4f}")
        yaz(f"  {v:4s} → τ₀ yayılımı {max(k[1] for k in kayd)-min(k[1] for k in kayd):.4f}"
            f" ({100*(max(k[1] for k in kayd)-min(k[1] for k in kayd))/np.mean([k[1] for k in kayd]):.2f}%)"
            f";  a yayılımı %{100*(max(k[2] for k in kayd)-min(k[2] for k in kayd))/am:.1f}"
            f";  b yayılımı %{100*(max(k[3] for k in kayd)-min(k[3] for k in kayd))/abs(bm):.1f}")
    if kuresel:
        aa = np.array([k[3] for k in kuresel])
        bb2 = np.array([k[4] for k in kuresel])
        t0s = np.array([k[2] for k in kuresel])
        yaz(f"\n  → BÜTÜN (pencere × taban) koşuları, n = {len(kuresel)}:")
        yaz(f"     a  = dφ/dτ|_τ₀ = {aa.mean():.3f} ± {aa.std(ddof=1):.3f} "
            f"(yayılım %{100*(aa.max()-aa.min())/aa.mean():.1f})")
        yaz(f"     b  = {bb2.mean():+.3f} ± {bb2.std(ddof=1):.3f} "
            f"(yayılım %{100*(bb2.max()-bb2.min())/abs(bb2.mean()):.1f})")
        yaz(f"     τ₀ = {t0s.mean():.4f} ± {t0s.std(ddof=1):.4f} "
            f"(yayılım {t0s.max()-t0s.min():.4f})")
        yaz("     Yani: taban konvansiyonu τ₀'ı 0.013 kaydırıyor ama EĞİMİ")
        yaz("     %2–4 içinde bırakıyor — φ'nin sıfır çevresindeki ŞEKLİ")
        yaz("     konvansiyondan bağımsız, YERİ değil. a·Δτ₀ = "
            f"{aa.mean()*(t0s.max()-t0s.min()):.3f} rad; φ'nin sabit τ'daki")
        yaz("     taban değişimi de bu mertebede (aşağıdaki denetim).")
        for v in VERI:
            vals = []
            for t in TABAN:
                k = (v, t, "ince") if (v, t, "ince") in D else (v, t, "tau0")
                if k not in D:
                    continue
                bb3 = [b for b in bant(D[k]) if abs(b["tau"] - 0.53) < 1e-9]
                if bb3:
                    vals.append(bb3[0]["phi"])
            if len(vals) >= 3:
                yaz(f"     denetim {v}: φ_Γ(τ̄=0.53) tabanla "
                    f"{min(vals):+.4f} … {max(vals):+.4f}, aralık "
                    f"{max(vals)-min(vals):.4f} rad")


# =====================================================================  D
def veri_kur(D, kanal="lad", g="kaba", taban_ana=0.46, taban_sis=TABAN,
             dRe_esik=None):
    P = []
    for v in VERI:
        if (v, taban_ana, g) not in D:
            continue
        for b in bant(D[(v, taban_ana, g)]):
            k = b["kanal"][kanal]
            if b["tau"] > TAU_UST or not saglam(b, kanal):
                continue
            if dRe_esik is not None and k["dRe"] >= dRe_esik:
                continue

            def yar_menzil(tbs):
                vs = []
                for t in tbs:
                    if (v, t, g) not in D:
                        continue
                    bb = [x for x in bant(D[(v, t, g)])
                          if x["tau"] == b["tau"]]
                    if bb and saglam(bb[0], kanal):
                        vs.append(bb[0]["kanal"][kanal]["Rham"])
                return (0.5 * (max(vs) - min(vs)) if len(vs) >= 2 else 0.0,
                        len(vs))
            s_tab, nb = yar_menzil(taban_sis)
            s_tab156, _ = yar_menzil(TABAN_156)
            oth = "orta" if v == "son" else "son"
            s_pen = 0.0
            if (oth, taban_ana, g) in D:
                bb = [x for x in bant(D[(oth, taban_ana, g)])
                      if x["tau"] == b["tau"]]
                if bb and bb[0]["kanal"][kanal]["artik"] < 0.02:
                    s_pen = 0.5 * abs(bb[0]["kanal"][kanal]["Rham"]
                                      - k["Rham"])
            sjk = k["sRham_jk"]
            if not np.isfinite(sjk) or sjk <= 0:
                sjk = 0.02 * abs(k["Rham"]) + 0.01
            sphi = b.get("sPhi_jk", float("nan"))
            P.append(dict(v=v, tau=b["tau"], x=b["tau_eff"], R=k["Rham"],
                          sjk=float(sjk), stab=float(s_tab),
                          stab156=float(s_tab156),
                          spen=float(s_pen), phi=b["phi"],
                          sphi=float(sphi), dRe=k["dRe"], neff=k["neff"],
                          olcek=k["olcek"], izgara=k.get("izgara"),
                          sA2=b["sA2"], nbase=nb))
    return P


# hata modelleri: (etiket, kullanılan bileşenler)
HATA_MODU = {
    "I":   ("yalnız istatistik (σ_jk)", ("sjk",)),
    "II":  ("istatistik + pencere (σ_jk ⊕ σ_pen)", ("sjk", "spen")),
    "III": ("+ 156'nın taban penceresi {0.46,0.52,0.58}",
            ("sjk", "spen", "stab156")),
    "IV":  ("+ TAM taban ekseni {0.28…0.58}  [dürüst konvansiyon bütçesi]",
            ("sjk", "spen", "stab")),
}


def sigma(P, mod="II"):
    bil = HATA_MODU[mod][1]
    return np.array([np.sqrt(sum(p[b]**2 for b in bil)) for p in P])


def _opt(chi2, th0, n=6):
    best = (chi2(np.array(th0, float)), np.array(th0, float))
    cur = best[1]
    for i in range(n):
        r = minimize(chi2, cur, method="Nelder-Mead",
                     options=dict(maxiter=40000, maxfev=40000,
                                  xatol=1e-11, fatol=1e-11))
        cur = r.x
        if r.fun < best[0]:
            best = (float(r.fun), np.array(r.x))
    return best


def adaylar(P, sg, tau0_olc):
    """Aday sözlüğü: ad -> (tanım, k, chi2_fn, th0, model_fn(th)->R dizisi)"""
    x = np.array([p["x"] for p in P])
    y = np.array([p["R"] for p in P])
    A2s = (2 * np.pi * x)**2 * np.array([p["sA2"] for p in P])
    olc = np.array([p["olcek"] for p in P])

    def mk(f):
        def chi2(th):
            try:
                m = f(th)
            except Exception:
                return 1e12
            if not np.all(np.isfinite(m)):
                return 1e12
            return float(np.sum(((y - m) / sg)**2))
        return chi2

    def h_model(th, kuad, b_sabit):
        a, t0 = th[0], th[1]
        b = (0.0 if not kuad else
             (b_sabit if b_sabit is not None else th[2]))
        out = np.empty(len(P))
        for i, p in enumerate(P):
            ph = a * (p["x"] - t0) + b * (p["x"] - t0)**2
            iz = p["izgara"]
            r = kok_izgara(iz["R"], iz["arg"], ph)
            if r is None:
                j = int(np.argmin(np.abs(np.array(iz["arg"]) - ph)))
                r = iz["R"][j] + 1e3          # kök yok → ağır ceza
            out[i] = r
        return out * olc

    def h_model_sabit_t0(th, kuad, b_sabit):
        return h_model([th[0], tau0_olc] + list(th[1:]), kuad, b_sabit)

    M = {}
    M["a"] = ("R = c  (sabit)", 1,
              lambda th: np.full_like(x, th[0]), [1.0])
    M["b"] = ("R = c·(τ−τ₀)  doğrusal", 2,
              lambda th: th[0] * (x - th[1]), [6.0, 0.51])
    M["c"] = (f"R = R∞(1−e^{{−(τ−τ₀)/w}}), τ₀={tau0_olc:.4f} SABİT (ölçülen)",
              2, lambda th: th[0] * (1 - np.exp(-(x - tau0_olc) / th[1])),
              [2.0, 0.12])
    M["c3"] = ("R = R∞(1−e^{−(τ−τ₀)/w})  [τ₀ serbest, 3 par]", 3,
               lambda th: th[0] * (1 - np.exp(-(x - th[1]) / th[2])),
               [2.0, 0.51, 0.12])
    M["g"] = ("R = c·(τ−τ₀)/τ²", 2,
              lambda th: th[0] * (x - th[1]) / x**2, [4.0, 0.515])
    M["g½"] = ("R = c·(τ−½)/τ²   [τ₀ ≡ ½ SABİT]", 1,
               lambda th: th[0] * (x - 0.5) / x**2, [3.0])
    M["e"] = ("R = k·A²σΔ²", 1, lambda th: th[0] * A2s, [0.5])
    M["h1"] = ("(h) φ = a(τ−τ₀) → R = argM⁻¹(φ)   [φ DOĞRUSAL]", 2,
               lambda th: h_model(th, False, None), [10.7, 0.510])
    M["h1½"] = ("(h) φ = a(τ−½) → R   [τ₀ ≡ ½ SABİT]", 1,
                lambda th: h_model([th[0], 0.5], False, None), [10.7])
    M["h1f"] = (f"(h) φ = a(τ−τ₀) → R, τ₀={tau0_olc:.4f} SABİT", 1,
                lambda th: h_model([th[0], tau0_olc], False, None), [10.7])
    M["h2b"] = ("(h) φ = a(τ−τ₀) − 6.5(τ−τ₀)²  [b=155'in değeri SABİT]", 2,
                lambda th: h_model(th, True, -6.5), [17.0, 0.509])
    M["h2"] = ("(h) φ = a(τ−τ₀) + b(τ−τ₀)²  [3 par]", 3,
               lambda th: h_model(th, True, None), [17.0, 0.509, -6.5])
    return {k: (v[0], v[1], mk(v[2]), v[3], v[2]) for k, v in M.items()}


def profil_hata(chi2, th, i, c2min, adim):
    """Δχ²=1 profili: i. parametre kaydırılır, diğerleri yeniden optimize."""
    out = []
    for yon in (+1, -1):
        d = adim * yon
        for _ in range(400):
            hedef = th.copy(); hedef[i] = th[i] + d
            if len(th) == 1:
                c2 = chi2(hedef)
            else:
                idx = [j for j in range(len(th)) if j != i]

                def f(u):
                    q = hedef.copy()
                    for j, uu in zip(idx, u):
                        q[j] = uu
                    return chi2(q)
                r = minimize(f, th[idx], method="Nelder-Mead",
                             options=dict(maxiter=8000, xatol=1e-10,
                                          fatol=1e-10))
                c2 = float(r.fun)
            if c2 - c2min >= 1.0:
                out.append(abs(d)); break
            d += adim * yon
        else:
            out.append(float("nan"))
    return float(np.nanmean(out))


def bolum_D(D, kayit, kanal="lad", g="kaba", taban_ana=0.46, baslik="",
            modlar=("I", "II", "III", "IV"), hata_mod=None, dRe_esik=None,
            ayrinti=True):
    yaz("\n" + "=" * 96)
    yaz(f"D — KAPALI-FORM YARIŞI {baslik}")
    yaz(f"    kanal={kanal}  ızgara={g}  ana taban={taban_ana}  "
        f"apsis=τ_eff  iki pencere ORTAK fit"
        + (f"  |ΔRe|<{dRe_esik} süzgeci" if dRe_esik else ""))
    yaz("=" * 96)
    P = veri_kur(D, kanal, g, taban_ana, dRe_esik=dRe_esik)
    if len(P) < 6:
        yaz("  yeterli nokta yok"); return None
    x = np.array([p["x"] for p in P]); y = np.array([p["R"] for p in P])
    t0m = float(np.nanmean([kayit[(v, taban_ana)][kanal] for v in VERI
                            if (v, taban_ana) in kayit
                            and np.isfinite(kayit[(v, taban_ana)][kanal])]))
    if not np.isfinite(t0m):
        t0m = 0.515
    yaz(f"  n = {len(P)} nokta; τ_eff ∈ [{x.min():.4f}, {x.max():.4f}]; "
        f"ölçülen τ₀({kanal}) = {t0m:.4f}")
    if ayrinti:
        yaz("\n  pen   τ̄     τ_eff     R_ölç    σ_jk   σ_pen  σ_tab156 "
            "σ_tab_tam  |ΔRe|   n_eff  nb")
        for p in P:
            yaz(f"  {p['v']:4s} {p['tau']:.4f} {p['x']:.4f} {p['R']:+8.3f} "
                f"{p['sjk']:7.4f} {p['spen']:6.4f} {p['stab156']:7.4f} "
                f"{p['stab']:8.4f}  {p['dRe']:.4f} {p['neff']:8.0f} "
                f"{p['nbase']:2d}")

    tum = {}
    for mod in modlar:
        sg = sigma(P, mod)
        AD = adaylar(P, sg, t0m)
        res = []
        for ad, (tanim, k, chi2, th0, mfn) in AD.items():
            c2, th = _opt(chi2, th0)
            res.append(dict(ad=ad, tanim=tanim, k=k, chi2=c2,
                            aic=c2 + 2 * k, th=th, chi2f=chi2, mfn=mfn))
        res.sort(key=lambda r: r["aic"])
        tum[mod] = (res, sg)
        yaz(f"\n  --- HATA MODELİ {mod}: {HATA_MODU[mod][0]} "
            f"(⟨σ⟩={np.mean(sg):.4f}) ---")
        yaz("  aday   k      χ²       AIC    ΔAIC   χ²/dof   parametreler")
        for r in res:
            pr = ", ".join(f"{u:+.4f}" for u in r["th"])
            yaz(f"  {r['ad']:5s} {r['k']:2d} {r['chi2']:10.2f} "
                f"{r['aic']:9.2f} {r['aic']-res[0]['aic']:7.2f} "
                f"{r['chi2']/max(len(P)-r['k'],1):7.2f}   [{pr}]   "
                f"{r['tanim']}")

    ana = hata_mod or modlar[-1]
    res, sg = tum[ana]
    yaz(f"\n  ARTIK YAPILARI (R_ölç − model), pencere '{P[0]['v']}':")
    yaz("  aday  " + " ".join(f"{p['x']:7.3f}" for p in P
                              if p['v'] == P[0]['v']))
    for r in tum["I"][0][:6]:
        m = r["mfn"](r["th"])
        vals = [f"{yy-mm:+7.3f}" for yy, mm, p in zip(y, m, P)
                if p['v'] == P[0]['v']]
        yaz(f"  {r['ad']:5s} " + " ".join(vals))

    if ayrinti:
        yaz(f"\n  PARAMETRE HATALARI (hata modeli {ana}; Δχ²=1 profili, "
            "diğer parametreler her adımda yeniden optimize)")
        for r in res[:5]:
            es = []
            for i in range(r["k"]):
                adim = max(abs(r["th"][i]) * 0.002, 1e-4)
                es.append(profil_hata(r["chi2f"], np.array(r["th"], float),
                                      i, r["chi2"], adim))
            yaz(f"  {r['ad']:5s} " + "  ".join(
                f"θ{i} = {t:+.4f} ± {e:.4f}"
                for i, (t, e) in enumerate(zip(r["th"], es))))
    return P, tum, t0m


# =====================================================================  E
def bolum_E(D, P, g="kaba", taban_ana=0.46):
    """φ'nin kendi şekli: doğrusal mı kuadratik mi (155'in b≈−6.5 sınavı)"""
    yaz("\n" + "=" * 96)
    yaz("E — İLKEL GÖZLENEBİLİR φ_Γ(τ_eff)'in ŞEKLİ (155'in b ≈ −6.5 sınavı)")
    yaz("=" * 96)
    yaz("  pen taban aralık          doğrusal: τ₀      a      χ²/dof | "
        "kuadratik: τ₀      a       b      χ²/dof | ΔAIC   hüküm")
    for v in VERI:
        for t in TABAN:
            for g2 in (("ince",) if (v, t, "ince") in D else ("tau0",)):
                if (v, t, g2) not in D:
                    continue
                bs = [b for b in bant(D[(v, t, g2)]) if b["tau"] <= 0.61]
                if len(bs) < 5:
                    continue
                xx = np.array([b["tau_eff"] for b in bs])
                yy = np.array([b["phi"] for b in bs])
                ee = np.array([b["sPhi_jk"] for b in bs])
                ee = np.where(np.isfinite(ee) & (ee > 0), ee, np.nanmedian(ee))
                out = []
                for deg in (1, 2):
                    c = np.polyfit(xx, yy, deg, w=1.0 / ee)
                    r = np.roots(c)
                    r = np.array([z.real for z in r if abs(z.imag) < 1e-9])
                    t0 = float(r[np.argmin(np.abs(r - xx.mean()))]) \
                        if len(r) else np.nan
                    c2 = float(np.sum(((yy - np.polyval(c, xx)) / ee)**2))
                    out.append((t0, c, c2, deg + 1))
                d1, d2 = out
                dof1 = max(len(xx) - 2, 1); dof2 = max(len(xx) - 3, 1)
                daic = (d2[2] + 2 * 3) - (d1[2] + 2 * 2)
                # kuadratiği (τ−τ₀) tabanında yaz: φ = a(τ−τ₀)+b(τ−τ₀)²
                t0q = d2[0]; cq = d2[1]
                a_q = float(cq[1] + 2 * cq[0] * t0q)
                b_q = float(cq[0])
                yaz(f"  {v:4s} {t:.2f} τ̄≤0.61 ({len(xx):2d} bant)  "
                    f"{d1[0]:.4f} {float(d1[1][0]):7.3f} {d1[2]/dof1:7.2f} | "
                    f"    {t0q:.4f} {a_q:7.3f} {b_q:+7.3f} {d2[2]/dof2:7.2f} "
                    f"| {daic:+6.2f}  "
                    f"{'EĞRİ' if daic < -2 else 'ayrılamaz'}")


if __name__ == "__main__":
    D = yukle()
    yaz(f"157 ANALİZ — {len(D)} koşu: "
        f"{sorted(set((k[0], k[1], k[2]) for k in D))}")
    bolum_A(D)
    bolum_A2(D, "kaba", 0.46)
    bolum_A2(D, "kaba", 0.52)
    kayit = bolum_B(D)
    bolum_B2(D, 0.46)
    bolum_B2(D, 0.40)
    bolum_C(D, kayit, "kaba")
    bolum_E(D, None)
    bolum_D(D, kayit, "lad", "kaba", 0.46, "[BİRİNCİL: R_lad]")
    bolum_D(D, kayit, "tam", "kaba", 0.46, "[KONTROL: R_tam]",
            ayrinti=False)
    bolum_D(D, kayit, "lad", "ince", 0.46,
            "[ROBUSTLUK 1: 0.02 ince ızgara]", modlar=("I", "IV"),
            ayrinti=False)
    bolum_D(D, kayit, "lad", "kaba", 0.52,
            "[ROBUSTLUK 2: ana taban 0.52]", modlar=("I", "IV"),
            ayrinti=False)
    bolum_D(D, kayit, "lad", "kaba", 0.40,
            "[ROBUSTLUK 3: ana taban 0.40]", modlar=("I", "IV"),
            ayrinti=False)
    bolum_D(D, kayit, "lad", "kaba", 0.46,
            "[ROBUSTLUK 4: yalnız aşırı-belirlemeyi GEÇEN bantlar "
            "(|ΔRe|<0.06)]", modlar=("I", "IV"), dRe_esik=0.06,
            ayrinti=False)
    (SCR / "analiz_cikti.txt").write_text("\n".join(OUT))
    print(f"\n-> {SCR/'analiz_cikti.txt'}")
