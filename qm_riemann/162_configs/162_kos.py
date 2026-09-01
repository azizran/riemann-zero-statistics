"""
162 — KOŞU SÜRÜCÜSÜ.   Kullanım: 162_kos.py <veri> <taban> <plan>

plan: virgülle ayrılmış çeşit listesi. Her çeşit:
      gercek            → 160'ın ta kendisi (bit düzeyi denetim için)
      V0                → boru hattı kontrolü (faz karıştırma YOK)
      VG<tohum>         → ANA Gauss vekil (çizgi + süreklilik)
      VL<tohum>         → yalnız çizgi fazları karıştırılmış
      VC<tohum>         → yalnız süreklilik karıştırılmış
Örn:  162_kos.py son 0.40 gercek,V0,VG101,VG102,VG103,VG104,VG105,VL101,VC101

Her çeşit ÖNCE t1 ızgarasında (τ∈0.44–0.80), G158 kümesindekiler ayrıca
158'in ızgarasında (a/b için) koşulur. Taban162 bir KEZ kurulur.
"""
import importlib
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
C = importlib.import_module("162_cekirdek")

IZG = {"t1": C.IZGARA_T1, "g158": C.IZGARA158}
G158 = ("gercek", "VG101", "VG102", "VG103")     # b/a için ek ızgara


if __name__ == "__main__":
    veri, taban, plan = sys.argv[1], float(sys.argv[2]), sys.argv[3]
    cesitler = plan.split(",")
    t0 = time.time()
    C.SCR.mkdir(parents=True, exist_ok=True)
    print(f"=== 162 / veri={veri} / taban={taban} / plan={plan} ===", flush=True)
    T = C.Taban162(veri, taban)
    print(f"    taban kuruldu ({(time.time()-t0)/60:.1f} dk)", flush=True)

    for ad in cesitler:
        tA = time.time()
        if ad == "gercek":
            Y = T.Y
        elif ad == "V0":
            Y = T.vekil_yerel(None, "V0")
        else:
            Y = T.vekil_yerel(int(ad[2:]), ad[:2])
        izgaralar = ["t1"] + (["g158"] if ad in G158 else [])
        for gad in izgaralar:
            bant = [b for b in IZG[gad] if b[0] >= taban - 1e-9]
            r = C.olc162(Y, f"{veri}-{ad}-{gad}", bant)
            r.update(veri=veri, taban=taban, izgara=gad, cesit=ad,
                     sure_s=time.time() - tA)
            if ad != "gercek":
                r["v162"] = {k: v for k, v in Y._162.items()}
            r["taban162"] = dict(pay_eta=T.pay_eta, pay_C=T.pay_C,
                                 kapanis=T.kapanis, m_kimlik=T.m_kimlik,
                                 nq=len(T.q), Var_Chat=float(T.Chat.var()))
            p = C.SCR / f"S_{veri}_t{taban}_{gad}_{ad}.json"
            p.write_text(json.dumps(r, indent=1))
            print(f"-> {p.name}  ({(time.time()-tA)/60:.1f} dk)", flush=True)
        del Y
    print(f"\nTOPLAM {(time.time()-t0)/60:.1f} dk", flush=True)
