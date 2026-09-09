# -*- coding: utf-8 -*-
"""
183e — K2: ALAN-SURROGATE TABANI (SUR-B) — YENİDEN ÇÖZÜM YOK
=============================================================
ÖN-MÜHÜR: 183/ONKAYIT_183.json → SUR["SUR_B"] (veriden önce donduruldu).

MAKİNE KOPYALANMAZ: `162_cekirdek.Taban162` (taban 0.40, cap 4000,
τ_çizgi 0.86, NITER 0, BLOK 2000) ve `cizgi_kur` / `cizgi_cikar`
AYNEN kullanılır — 174b/176d'nin ta kendisi.

SUR-B (ön-kayıttan, birebir):
    η°  := cizgi_kur(m, |c_η|·e^{iψ}, w) + η_artık(VF)
    c°  := cizgi_cikar(m, η°, w)
    R°_η := Σ|c°|²/2 / Var(η°) ;  Kov-oranı° := 1 − R°_η
Aynı işlem X̃ (=ΔĈ) kanalında; m3°_çizgi := m3(η°_çizgi, X̃°_çizgi).

DİKKAT (raporda aynen yazılacak): SUR-B yalnız ÇİZGİ FAZLARINI
rastgeleleştirir. η'nın artık bileşeninin çizgi gösterimi yoktur ve
yeniden çözmeden surrogate'lanamaz; bu yüzden artık AYNEN devralınır.
Yani SUR-B "çizgi↔artık faz eşleşmesi silinmiş" tabandır.

ÜREME KAPISI: surrogate'tan ÖNCE, aynı kodun ölçtüğü R_η, P, Var(η),
Kov(artık,η), Kov(çizgi,artık), m3_çizgi değerleri 176/K1_<gaz>.json
ile karşılaştırılır (|Δ| basılır).

Kullanım: 183e_surrogate_taban.py VF1 VF2 VF3 VF4
"""
import importlib
import json
import sys
import time
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
for _p in ("162_configs", "160_configs", "159_configs", "156_configs",
           "155_configs", "154_configs"):
    sys.path.insert(0, str(QM / _p))
C162 = importlib.import_module("162_cekirdek")

SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S183, S176 = SCR / "183", SCR / "176"
TABAN, CAP = 0.40, 4000


def _m3(a, b):
    a = a - a.mean()
    b = b - b.mean()
    sa, sb = a.std(), b.std()
    return float((a * b * b).mean() / (sa * sb * sb)) if sa * sb > 0 else np.nan


def defter(x, c, m, w):
    """Bir kanalın çizgi/artık defteri (174b'nin kanal_defteri'nin çekirdeği,
    bant bölmesi olmadan — toplam çizgi alanı tek seferde kurulur)."""
    xciz = C162.cizgi_kur(m, [c], w)[0]
    art = x - xciz
    x0 = x - x.mean()
    V = float(np.var(x))
    P = float(np.sum(np.abs(c) ** 2) / 2.0)
    return dict(Var=V, P=P, R=P / V, Var_cizgi=float(np.var(xciz)),
                Var_artik=float(np.var(art)),
                kov_art_x=float(np.mean((art - art.mean()) * x0)),
                kov_ciz_art=P - float(np.var(xciz)),
                kovoran=float(np.mean((art - art.mean()) * x0)) / V), xciz, art


def kos(gaz, npsi=1):
    t0 = time.time()
    print("=" * 78, flush=True)
    print(f"183e — SUR-B alan-surrogate tabanı: {gaz}  (162 makinesi aynen)",
          flush=True)
    print("=" * 78, flush=True)
    T = C162.Taban162(gaz, taban=TABAN, cap=CAP)
    Y = T.Y
    m, w = Y.mid, T.w
    print(f"  [Taban162 {time.time()-t0:.0f}s]  nline={len(w)}  N={len(Y.eta)}",
          flush=True)

    # ---------------- ÖLÇÜLEN (üreme kapısı) --------------------------
    d_eta, eta_ciz, eta_art = defter(Y.eta, T.c_eta, m, w)
    s = np.ascontiguousarray(m[1:])
    x1 = Y.Xtil0
    c_X = C162.cizgi_cikar(s, [x1], w)[0]
    d_X, X_ciz, X_art = defter(x1, c_X, s, w)
    m3_ciz = _m3(eta_ciz[1:], X_ciz)
    ref = json.load(open(S176 / f"K1_{gaz}.json"))
    print(f"  ÖLÇÜLEN : R_η={d_eta['R']:.5f}  Kov-oranı={d_eta['kovoran']:+.5f}"
          f"  Kov(çiz,art)={d_eta['kov_ciz_art']:+.6f}  m3_çiz={m3_ciz:+.5f}",
          flush=True)
    ur = dict(R=abs(d_eta["R"] - ref["eta"]["R"]),
              P=abs(d_eta["P"] - ref["eta"]["P"]),
              Var=abs(d_eta["Var"] - ref["eta"]["Var"]),
              kov_art=abs(d_eta["kov_art_x"] - ref["eta"]["kov_art_x"]),
              kov_ciz_art=abs(d_eta["kov_ciz_art"]
                              - (ref["eta"]["P"] - ref["eta"]["Var_cizgi"])),
              m3=abs(m3_ciz - ref["m3"]["cizgi"]))
    print(f"  ÜREME KAPISI |Δ| vs 176/K1_{gaz}.json: "
          + "  ".join(f"{k}={v:.2e}" for k, v in ur.items()), flush=True)

    # ---------------- SUR-B -------------------------------------------
    sonuc = []
    for j in range(npsi):
        tj = time.time()
        rng = np.random.default_rng(183000 + 100 * int(gaz.replace("VF", "")) + j)
        ce = np.abs(T.c_eta) * np.exp(1j * rng.uniform(0, 2 * np.pi, len(w)))
        cx = np.abs(c_X) * np.exp(1j * rng.uniform(0, 2 * np.pi, len(w)))
        eta_s = C162.cizgi_kur(m, [ce], w)[0] + eta_art
        x1_s = C162.cizgi_kur(s, [cx], w)[0] + X_art
        cse = C162.cizgi_cikar(m, [eta_s], w)[0]
        csx = C162.cizgi_cikar(s, [x1_s], w)[0]
        de, ec_s, ea_s = defter(eta_s, cse, m, w)
        dx, xc_s, xa_s = defter(x1_s, csx, s, w)
        m3s = _m3(ec_s[1:], xc_s)
        print(f"  SUR-B[{j}] R°_η={de['R']:.5f}  Kov-oranı°={de['kovoran']:+.5f}"
              f"  Kov(çiz,art)°={de['kov_ciz_art']:+.6f}  m3°_çiz={m3s:+.5f}"
              f"  R°_X̃={dx['R']:.5f}   [{time.time()-tj:.0f}s]", flush=True)
        sonuc.append(dict(psi=j, eta=de, Xtil=dx, m3_cizgi=m3s))

    rec = dict(gaz=gaz, olculen=dict(eta=d_eta, Xtil=d_X, m3_cizgi=m3_ciz),
               ureme_kapisi=ur, surrogate=sonuc,
               nline=int(len(w)), N=int(len(Y.eta)), sure_s=time.time() - t0)
    S183.mkdir(parents=True, exist_ok=True)
    p = S183 / f"SURB_{gaz}.json"
    p.write_text(json.dumps(rec, indent=1, ensure_ascii=False, default=float))
    print(f"  -> {p}  ({(time.time()-t0)/60:.1f} dk)", flush=True)
    return rec


if __name__ == "__main__":
    gz = sys.argv[1:] or ["VF1", "VF2", "VF3", "VF4"]
    for g in gz:
        kos(g)
