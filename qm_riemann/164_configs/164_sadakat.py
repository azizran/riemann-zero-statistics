"""
164 — SADAKAT DENETİMİ: gerçekleşen / nominal merdiven genliği (bant bant)
=========================================================================
163 §1c'nin yöntemi (`163_tani_merdiven.py`) — ÇİZGİ bazında:

    ds_n = g_n·log(m_n/2π)/2π − 1                    (ölçüm zincirinin ds'i)
    c_q  = 2·⟨(ds − ⟨ds⟩)·e^{−iω_q m}⟩                (çizgi izdüşümü)
    b_q  = 2·a_q^(pencereli)·sin(π τ_q) ,  τ_q = ω_q/L

    gerçekleşen/nominal = |c_q| / b_q        (hedef 1.00; faz arg c_q ≈ 0)

164'ÜN İKİ EKLEMESİ
  (1) BANT AGREGASYONU + GÜRÜLTÜ ÇIKARIMI. Yüksek τ'da b_q ölçüm gürültü
      tabanının altına iner (tek çizgi SNR < 1). Bu yüzden her çizgi için
      ölçüm zincirinin kendi ARA-NOKTA referansı da ölçülür
      (W' = w + gap/2, aynı 2.5·dres eleği) ve bant oranı
          R_bant = sqrt( (Σ|c_on|² − Σ|c_off|²) / Σ b² )
      ile gürültü-çıkarımlı verilir. Ham (çıkarımsız) sürüm de basılır.
  (2) NOMİNAL = PENCERELİ. Sentetik gazın inşasına GİREN genlik hangisiyse
      sadakat ona göre ölçülür (A4 için erfc ile çarpılmış a_q). 163
      penceresiz b_q kullanıyordu; her iki sütun da veriliyor, çünkü
      163'ün tablosuyla süreklilik gerekiyor.

GERÇEK VERİ KONTROLÜ: `son` penceresiz merdivene karşı ölçülür ve oranın
~1.00 kaldığı görülür (yöntemin kendisinin sınavı).

Kullanım:  164_sadakat.py <veri> [<veri> ...]
"""
import importlib
import json
import sys
import time
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
sys.path.insert(0, str(QM / "155_configs"))
sys.path.insert(0, str(QM / "154_configs"))
KOS155 = importlib.import_module("155_kos")
C154 = importlib.import_module("154_cekirdek")

SCR164 = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
              "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/164")
TWO_PI = 2 * np.pi

# gazın inşasında kullanılan erfc penceresi (yoksa None)
PENCERE = {"A4": (0.68, 0.125), "SA4": (0.68, 0.125), "eski_A4": (0.68, 0.125),
           "HA4": (0.68, 0.125), "SA1": (0.75, 0.10)}
# 163 §1c tablosunun çizgileri (birebir karşılaştırma için)
Q163 = (2, 3, 5, 7, 11, 101, 1009)
BANT = [(round(0.05 * i, 3), round(0.05 * (i + 1), 3)) for i in range(20)]
NMAX = 200          # bant başına örneklenen çizgi (tohum 21)


def erfc_w(tau, tc, dl):
    import math
    return np.array([0.5 * math.erfc((float(x) - tc) / dl) for x in tau])


def olc(veri, tau_ust=0.95, tohum=21):
    tb = time.time()
    z = np.sort(KOS155.veri_yukle(veri))
    g = np.diff(z)
    mid = 0.5 * (z[:-1] + z[1:])
    Lw = np.log(mid / TWO_PI)
    L = float(Lw.mean())
    ds = g * Lw / TWO_PI - 1
    dsm = ds - ds.mean()
    dres = TWO_PI / (mid[-1] - mid[0])
    qm = C154.pk_m(int(np.exp(1.0 * L)))
    allq = sorted(qm)
    allw = np.array([np.log(q) for q in allq])
    pen = PENCERE.get(veri)
    print(f"=== SADAKAT {veri}: n={len(z)}  L={L:.4f}  σ_ds²={np.var(ds):.4f}  "
          f"çizgi={len(allq)}  dres={dres:.3e}  "
          f"pencere={pen}", flush=True)

    rng = np.random.default_rng(tohum)
    sec = {}
    for lo, hi in BANT:
        if hi > tau_ust + 1e-9:
            continue
        idx = [i for i, w in enumerate(allw) if lo < w / L <= hi]
        if len(idx) > NMAX:
            idx = [idx[i] for i in rng.choice(len(idx), NMAX, replace=False)]
        sec[(lo, hi)] = idx
    # 163'ün çizgileri her hâlükârda ölçülsün
    ekstra = [i for i, q in enumerate(allq) if q in Q163]

    kayit = {}
    hedefler = sorted(set([i for v in sec.values() for i in v] + ekstra))
    for i in hedefler:
        q, w = allq[i], allw[i]
        koms = [allw[k] for k in (i - 1, i + 1)
                if 0 <= k < len(allw) and abs(allw[k] - w) > 1e-12]
        gap = min(abs(w - k) for k in koms)
        if gap < 2.5 * dres:
            continue
        tau = w / L
        a_ham = 1.0 / (np.pi * qm[q] * np.sqrt(q))
        wg = 1.0 if pen is None else float(erfc_w(np.array([tau]), *pen)[0])
        b_ham = 2 * a_ham * np.sin(np.pi * tau)
        b_pen = b_ham * wg
        cs = {}
        for W, et in ((w, "on"), (w + gap / 2, "off")):
            cw, sw = np.cos(W * mid), np.sin(W * mid)
            cs[et] = complex(2 * np.mean(dsm * cw), -2 * np.mean(dsm * sw))
        kayit[q] = dict(q=int(q), tau=float(tau), w=float(w), wg=wg,
                        b_ham=float(b_ham), b_pen=float(b_pen),
                        con=[cs["on"].real, cs["on"].imag],
                        cof=[cs["off"].real, cs["off"].imag])

    # ---- 163 tablosu (çizgi bazında)
    print("   q   τ       b_nom(pen)  |c_on|     oran_pen  oran_ham   arg",
          flush=True)
    for q in Q163:
        r = kayit.get(q)
        if r is None:
            continue
        c = complex(*r["con"])
        print(f"{q:5d} {r['tau']:.4f}  {r['b_pen']:9.5f}  {abs(c):9.5f}  "
              f"{abs(c)/r['b_pen']:8.3f}  {abs(c)/r['b_ham']:8.3f}  "
              f"{np.angle(c):+7.4f}", flush=True)

    # ---- bant tablosu
    print("\n  τ-bant      nq   Σb_pen^½    R_bant(çıkarımlı)  R_ham   "
          "R_bant(çıkarımsız)  SNR", flush=True)
    bant_out = []
    for (lo, hi), idx in sec.items():
        rs = [kayit[allq[i]] for i in idx if allq[i] in kayit]
        if not rs:
            continue
        Pon = sum(abs(complex(*r["con"]))**2 for r in rs)
        Pof = sum(abs(complex(*r["cof"]))**2 for r in rs)
        Bp = sum(r["b_pen"]**2 for r in rs)
        Bh = sum(r["b_ham"]**2 for r in rs)
        net = Pon - Pof
        R = float(np.sqrt(max(net, 0.0) / Bp)) if Bp > 0 else float("nan")
        Rh = float(np.sqrt(max(net, 0.0) / Bh)) if Bh > 0 else float("nan")
        Rraw = float(np.sqrt(Pon / Bp)) if Bp > 0 else float("nan")
        snr = float(Pon / Pof) if Pof > 0 else float("inf")
        fazlar = [float(np.angle(complex(*r["con"]))) for r in rs
                  if abs(complex(*r["con"])) > 3 * np.sqrt(Pof / len(rs))]
        bant_out.append(dict(lo=lo, hi=hi, nq=len(rs), R=R, R_ham=Rh,
                             R_ham_raw=Rraw, snr=snr,
                             Bp=float(np.sqrt(Bp)), Pon=float(Pon),
                             Pof=float(Pof),
                             faz_maks=float(np.max(np.abs(fazlar)))
                             if fazlar else float("nan"),
                             nfaz=len(fazlar)))
        print(f"  {lo:.2f}-{hi:.2f}  {len(rs):4d}  {np.sqrt(Bp):9.5f}   "
              f"{R:8.3f}          {Rh:7.3f}   {Rraw:8.3f}          "
              f"{snr:8.2f}", flush=True)

    out = dict(veri=veri, L=L, sigma_ds2=float(np.var(ds)), n=int(len(z)),
               pencere=pen, dres=float(dres), tau_ust=tau_ust,
               cizgi={str(k): v for k, v in kayit.items()},
               bant=bant_out, sure_s=time.time() - tb)
    SCR164.mkdir(parents=True, exist_ok=True)
    (SCR164 / f"sadakat_{veri}.json").write_text(json.dumps(out, indent=1))
    print(f"-> sadakat_{veri}.json  ({time.time()-tb:.0f}s)\n", flush=True)
    return out


if __name__ == "__main__":
    for v in sys.argv[1:]:
        olc(v)
