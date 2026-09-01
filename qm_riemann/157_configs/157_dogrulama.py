"""
157 — KOPYA-KAYMASI DENETİMİ
===========================
Üç ayrı düzeyde, 155/154/156'nın kayıtlı çıktılarına karşı:

 (V1) η ÖZDEŞLİĞİ. 156.zincir3'ün ayrıştırdığı η, 154.eta_zinciri'nin
      (155'in önbelleğe aldığı) η'sı ile aynı mı? ds = drift+lad+eta
      kapanışı ayrıca ölçülür.
 (V2) BANT DÖNGÜSÜ. 157'nin Γ, φ_Γ, τ_eff, ρ değerleri 155'in AYNI
      pencere/taban/bant için kaydettikleriyle aynı mı?
 (V3) R KÖKÜ. 157'nin R_tam'ı 155'in R'si ile aynı mı? (155 Rmin=−6,
      adım 0.02; 157 Rmin=−10, adım 0.05 + bisection — kök R=0'dan DIŞA
      arandığı için aralık/adım farkı kökü değiştirmemeli.)
 (V4) 156 REFERANSI. 156'nın std bantlarındaki R_lad(ham) değerleri
      157'nin aynı taban/pencere ölçümüyle uyumlu mu? (Bantlar farklı
      genişlikte olduğu için birebir değil, YÖN denetimi.)
"""
import importlib
import json
import sys
from pathlib import Path

import numpy as np

_C = Path(__file__).resolve().parent
sys.path.insert(0, str(_C))
sys.path.insert(0, str(_C.parent / "154_configs"))
sys.path.insert(0, str(_C.parent / "156_configs"))
K = importlib.import_module("157_cekirdek")
C154 = importlib.import_module("154_cekirdek")
C156 = importlib.import_module("156_cekirdek")

SCR = K.SCR
SCR155 = SCR.parent / "155"
SCR156 = SCR.parent / "156"
HERE = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")


def z_yukle(ad):
    d = np.load(HERE / "128_odl_zeros6_2e6_zeros.npz")
    Z = np.sort(np.asarray(d["zeros"], float))
    return Z[len(Z) - 300000:] if ad == "son" else Z[850000:1150000]


def V1():
    print("=" * 78)
    print("V1 — η ÖZDEŞLİĞİ (156.zincir3 vs 154.eta_zinciri / 155 önbelleği)")
    print("=" * 78)
    print("  pencere taban  maks|Δη|    Δσ_η²      Δc₁      ds-kapanış")
    for v in ("son", "orta"):
        for t in (0.34, 0.40, 0.46, 0.52):
            p155 = SCR155 / f"eta_{v}_t{t:g}_c4000.npz"
            p157 = SCR / f"kanal_{v}_t{t}_c4000.npz"
            if not (p155.exists() and p157.exists()):
                print(f"  {v:5s} {t:.2f}   (önbellek yok: "
                      f"{'155' if not p155.exists() else '157'})")
                continue
            a = np.load(p155); b = np.load(p157)
            de = float(np.max(np.abs(a["eta"] - b["eta"])))
            kap = float(np.max(np.abs(b["ds"] - (b["drift"] + b["lad"]
                                                 + b["eta"]))))
            print(f"  {v:5s} {t:.2f}  {de:.3e}  "
                  f"{abs(float(a['s_eta'])-float(b['s_eta'])):.3e}  "
                  f"{abs(float(a['c1'])-float(b['c1'])):.3e}  {kap:.3e}")


def V23():
    print("\n" + "=" * 78)
    print("V2/V3 — BANT DÖNGÜSÜ ve R KÖKÜ (155'in kayıtlı JSON'larına karşı)")
    print("=" * 78)
    ciftler = [("son", 0.40, "ince", "tau0_son_t0.4_ince.json"),
               ("orta", 0.40, "ince", "tau0_orta_t0.4_ince.json"),
               ("son", 0.46, "ince", "tau0_son_t0.46_ince.json"),
               ("son", 0.52, "ince", "tau0_son_t0.52_yuksek.json"),
               ("orta", 0.52, "ince", "tau0_orta_t0.52_yuksek.json")]
    for v, t, g, ref in ciftler:
        p = SCR / f"R_{v}_t{t}_{g}.json"
        q = SCR155 / ref
        if not (p.exists() and q.exists()):
            print(f"  [{v} t{t}] atlandı (dosya yok)")
            continue
        a = {b["tau"]: b for b in json.load(open(p))["bantlar"]
             if b.get("olculdu")}
        b_ = {b["tau"]: b for b in json.load(open(q))["bantlar"]
              if b.get("olculdu")}
        ort = sorted(set(a) & set(b_))
        if not ort:
            print(f"  [{v} t{t}] ortak bant yok")
            continue
        d = {k: [] for k in ("absG", "phi", "tau_eff", "rho", "R")}
        for tau in ort:
            d["absG"].append(abs(a[tau]["absG"] - b_[tau]["absG"]))
            d["phi"].append(abs(a[tau]["phi"] - b_[tau]["phi"]))
            d["tau_eff"].append(abs(a[tau]["tau_eff"] - b_[tau]["tau_eff"]))
            d["rho"].append(abs(a[tau]["rho"] - b_[tau]["rho"]))
            if a[tau]["kanal"]["tam"]["artik"] < 0.02 \
                    and b_[tau]["artik"] < 0.02:
                d["R"].append(abs(a[tau]["kanal"]["tam"]["Rham"]
                                  - b_[tau]["R"]))
        print(f"  [{v} taban {t}] {len(ort)} ortak bant ({ref}):  "
              f"maks Δ|Γ|={max(d['absG']):.2e}  Δφ={max(d['phi']):.2e}  "
              f"Δτ_eff={max(d['tau_eff']):.2e}  Δρ={max(d['rho']):.2e}  "
              f"ΔR_tam={max(d['R']) if d['R'] else float('nan'):.2e} "
              f"({len(d['R'])} bant)")


def V4():
    print("\n" + "=" * 78)
    print("V4 — 156 REFERANSI (R_lad ham). 156 bantları: (0.525,0.55),")
    print("     (0.55,0.62), (0.62,0.70), (0.70,0.78), (0.78,0.85);")
    print("     157'ninkiler 0.03 genişlikte → birebir DEĞİL, yön denetimi.")
    print("=" * 78)
    p = SCR156 / "k3_gercek_std.json"
    if not p.exists():
        print("  156 çıktısı yok"); return
    r = json.load(open(p))
    print("  156 (taban 0.52, son): " + "  ".join(
        f"τ={b['tau']:.4f}→R_lad={b['kanal']['lad']['Rham']:+.3f}"
        for b in r["bantlar"] if b.get("olculdu")))
    q = SCR / "R_son_t0.52_kaba.json"
    if q.exists():
        s = json.load(open(q))
        print("  157 (taban 0.52, son): " + "  ".join(
            f"τ={b['tau']:.4f}→R_lad={b['kanal']['lad']['Rham']:+.3f}"
            for b in s["bantlar"] if b.get("olculdu")))
    for ad, tb in (("t0.46", 0.46), ("t0.58", 0.58)):
        p2 = SCR156 / f"k3_gercek_gorev_{ad}.json"
        if p2.exists():
            r2 = json.load(open(p2))
            print(f"  156 (taban {tb}, görev bantları): " + "  ".join(
                f"τ={b['tau']:.4f}→R_lad={b['kanal']['lad']['Rham']:+.3f}"
                for b in r2["bantlar"] if b.get("olculdu")))


def V5():
    """R_k(τ)=0 kimliğinin DOĞRUDAN sınavı: M_k(0) kanaldan bağımsız mı?"""
    print("\n" + "=" * 78)
    print("V5 — KİMLİK: M_k(R=0) kanaldan bağımsız mı? (τ₀ özdeşliğinin")
    print("     analitik temeli; sayıyla doğrulanıyor)")
    print("=" * 78)
    z = z_yukle("son")[:120000]
    C = C156.zincir3(z, 0.46, 4000)
    X = {k: C156._bond(C[k] if k != "tam" else C["ds"])
         for k in ("tam", "lad", "eta", "drift")}
    Xt = X["tam"]
    for tau in (0.47, 0.51, 0.55, 0.66):
        A = 2 * np.pi * tau
        ex = np.exp(-1j * A * Xt)
        vals = []
        for k in ("tam", "lad", "eta", "drift"):
            M, _ = C154._M(0.0, A, X[k], ex)
            vals.append(M)
        print(f"  τ={tau:.2f}  " + "  ".join(
            f"{k}: {np.angle(m):+.10f}" for k, m in
            zip(("tam", "lad", "eta", "dri"), vals))
              + f"   maks fark = {max(abs(m-vals[0]) for m in vals):.2e}")


def V6():
    """Ağırlık OPERATÖRÜ tabana bağlı mı? (raporun analitik omurgası)

    X_tam = bond(ds) TABANDAN BAĞIMSIZDIR (ds regresyondan önce tanımlı).
    Dolayısıyla arg M_tam(R) eğrisi de, M_emp de tabandan bağımsız olmalı —
    yani R_tam, ilkel φ_Γ'nın SABİT ve tersinir bir yeniden-parametrizasyonu.
    X_lad ise TANIMI GEREĞİ tabana bağlıdır (hangi çizgiler merdivene
    sayılıyor?), yani R_lad φ'nin konvansiyonuna KENDİ konvansiyonunu ekler.
    """
    print("\n" + "=" * 78)
    print("V6 — AĞIRLIK OPERATÖRÜNÜN TABAN BAĞIMLILIĞI")
    print("=" * 78)
    for v in ("son", "orta"):
        ref = SCR / f"R_{v}_t0.52_kaba.json"
        if not ref.exists():
            continue
        R = json.load(open(ref))
        b_ = {b["tau"]: b for b in R["bantlar"] if b.get("olculdu")}
        for t in (0.28, 0.34, 0.40, 0.46, 0.58):
            p = SCR / f"R_{v}_t{t}_kaba.json"
            if not p.exists():
                continue
            a = {b["tau"]: b for b in json.load(open(p))["bantlar"]
                 if b.get("olculdu")}
            ort = sorted(set(a) & set(b_))
            if not ort:
                continue
            dt = max(float(np.max(np.abs(
                np.array(a[x]["kanal"]["tam"]["izgara"]["arg"])
                - np.array(b_[x]["kanal"]["tam"]["izgara"]["arg"]))))
                for x in ort)
            dme = max(abs(a[x]["argMe"] - b_[x]["argMe"]) for x in ort)
            dl = max(float(np.max(np.abs(
                np.array(a[x]["kanal"]["lad"]["izgara"]["arg"])
                - np.array(b_[x]["kanal"]["lad"]["izgara"]["arg"]))))
                for x in ort)
            print(f"  {v:4s} taban {t:.2f} vs 0.52 ({len(ort)} ortak bant): "
                  f"maks|Δ arg M_tam(R)| = {dt:.3e}   "
                  f"maks|Δ arg M_emp| = {dme:.3e}   "
                  f"maks|Δ arg M_lad(R̃)| = {dl:.3e}")


def V7():
    """(h) adayının sayısal temeli: JSON'a yazılan arg M ızgarasının tersi,
    koşuda bisection'la bulunan kökü geri veriyor mu?"""
    print("\n" + "=" * 78)
    print("V7 — IZGARA TERSİ ?= BİSECTION KÖKÜ  ((h) adayının hata payı)")
    print("=" * 78)
    Z = importlib.import_module("157_analiz")
    for v, t, g in (("son", 0.46, "kaba"), ("orta", 0.46, "kaba"),
                    ("son", 0.52, "kaba")):
        p = SCR / f"R_{v}_t{t}_{g}.json"
        if not p.exists():
            continue
        d = json.load(open(p))
        dm = {"tam": 0.0, "lad": 0.0}
        for b in [x for x in d["bantlar"] if x.get("olculdu")]:
            for kan in ("tam", "lad"):
                k = b["kanal"][kan]; iz = k["izgara"]
                r = Z.kok_izgara(iz["R"], iz["arg"], b["phi"])
                if r is not None:
                    dm[kan] = max(dm[kan], abs(r - k["Rn"]))
        print(f"  {v} taban {t}: maks|Δ R̃| tam={dm['tam']:.2e}  "
              f"lad={dm['lad']:.2e}   (σ_jk ≈ 5e−3 ile kıyasla ihmal "
              f"edilebilir)")


if __name__ == "__main__":
    V1(); V23(); V4(); V5(); V6(); V7()
