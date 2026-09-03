"""
170 — K3 (TAÇ): GERÇEK GAZIN +3.06σ FAZLASI T-YASASIYLA KAPANIYOR MU?
=====================================================================
Yeni ölçüm YOK: `scratchpad/170/{K0.json, K1_T.json}` ve
`scratchpad/169/K2b_*.json` (170e'nin koşturduğu `169_k2b.main`) okunur.

H-D2 (KALEM): 162'ye göre kilitli asal fazlar η'yı inkoheran toplamdan
%29 SESSİZLEŞTİRİR ⇒ gerçek gazda birikmiş-faz σ_b'leri sadakatli
sentetikten KÜÇÜK ⇒ daha az sarılma ⇒ T daha BÜYÜK ⇒ c daha YÜKSEK.
İşaret bu yönde mi, büyüklük tutuyor mu?

Sınavlar:
  (K3-a) σ_Ĉ(son) ↔ σ_Ĉ(Hkeskin): 162'nin %29'u görünüyor mu?
  (K3-b) T-yasası fazlanın ne kadarını kapatıyor?
         beklenen c(son)/c(Hk) = ΠT(son)/ΠT(Hk)
  (K3-c) TUTARLILIK: aynı σ_Ĉ→c kuralı λ ekseninde de geçerli mi?
         (λ dizisinden ÖLÇÜLEN dc/dσ_Ĉ ile `son` öngörülür)
  (K3-d) korelatör düzeyinde: `son`un bacak kırpma aktarımları
         Hkeskin'inkinden BÜYÜK mü?

ÖN-MÜHÜR (koşudan ÖNCE):
  * σ_Ĉ(son)/σ_Ĉ(Hk) = 0.27303/0.27768 = 0.9833 — yani **−%1.7**,
    162'nin %29'u DEĞİL ⇒ (K3-a) büyüklükte tutmayacak.
  * ΠT(son)/ΠT(Hk) = 0.78140/0.77072 = 1.0139 (+%1.4); ölçülen
    c oranı 0.4122/0.4035 = 1.0216 (+%2.2) ⇒ T-yasası fazlanın ~2/3'ünü
    kapatır, artık ≈ +1σ'ya iner. **İŞARET doğru, büyüklük eksik.**
  * (K3-c) λ ekseninden ölçülen dc/dσ_Ĉ POZİTİF ve dik; `son`un daha
    küçük σ_Ĉ'si oradan c'nin DÜŞMESİNİ ister ⇒ fazla BÜYÜR (≥ +4σ).
    İki kural çelişecek.

Çıktı: scratchpad/170/K3.json
"""
import json
import sys
from pathlib import Path

import numpy as np

SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
C_HIP = 4 / np.pi ** 2
NB = {"100_E": 1, "010_Xa": 1, "001_Xb": 1,
      "011_XaXb": 2, "110_EXa": 2, "111_hepsi": 3}


def main():
    K0 = json.load(open(SCR / "170/K0.json"))["gaz"]
    KT = json.load(open(SCR / "170/K1_T.json"))["gaz"]
    h, s = K0["Hkeskin"], K0["son"]
    th, ts = KT["Hkeskin"], KT["son"]
    print("=" * 96)
    print("K3 (TAÇ) — GERÇEK GAZIN FAZLASI ve H-D2")
    print("=" * 96)

    d = s["c_WX"] - h["c_WX"]
    sd = float(np.hypot(s["s_tot"], h["s_tot"]))
    print(f"\n  ÖLÇÜLEN FAZLA:  c(son) − c(Hkeskin) = {s['c_WX']:.4f} − "
          f"{h['c_WX']:.4f} = {d:+.4f} ± {sd:.4f}  = **{d/sd:+.2f}σ**")
    print(f"     (169'un '+3.06σ'sı 4/π²'ye göreydi: c(son) − 4/π² = "
          f"{s['c_WX']-C_HIP:+.4f} = {(s['c_WX']-C_HIP)/s['s_tot']:+.2f}σ)")

    print("\n  (K3-a) 162'nin %29 SESSİZLEŞMESİ MARJİNALLERDE GÖRÜNÜYOR MU?")
    for k, et in (("sigds", "σ_ds"), ("sigX", "σ_X̃"), ("sigC", "σ_Ĉ")):
        print(f"     {et:5s}: son/Hkeskin = {s[k]/h[k]:.4f}  "
              f"({100*(s[k]/h[k]-1):+.2f}%)")
    print("     162'nin beklentisi: −%29.  ⇒ marjinal düzeyde GÖRÜNMÜYOR.")

    print("\n  (K3-b) T-YASASI FAZLANIN NE KADARINI KAPATIYOR?")
    rt = ts["PiT"] / th["PiT"]
    pred = h["c_WX"] * rt
    print(f"     Π T(son)/Π T(Hk) = {ts['PiT']:.5f}/{th['PiT']:.5f} = "
          f"{rt:.4f}  ({100*(rt-1):+.2f}%)")
    print(f"     ⇒ öngörülen c(son) = {pred:.4f} ;  ölçülen "
          f"{s['c_WX']:.4f} ;  artık = {s['c_WX']-pred:+.4f} = "
          f"**{(s['c_WX']-pred)/sd:+.2f}σ**")
    print(f"     kapanan pay: {100*(rt-1)/(s['c_WX']/h['c_WX']-1):.0f}% "
          f"(+3.0σ → {(s['c_WX']-pred)/sd:+.2f}σ)")

    print("\n  (K3-c) TUTARLILIK: aynı kural λ ekseninde de geçerli mi?")
    lam = ["Hkeskin", "L085", "L070", "L060"]
    x = np.array([K0[a]["sigC"] for a in lam])
    y = np.array([K0[a]["c_WX"] for a in lam])
    p = np.polyfit(x, y, 1)
    pl = float(np.polyval(p, s["sigC"]))
    print(f"     λ dizisinden ÖLÇÜLEN doğru: c = {p[0]:+.4f}·σ_Ĉ "
          f"{p[1]:+.4f}   (σ_Ĉ: " +
          " ".join(f"{v:.4f}" for v in x) + ")")
    print(f"     `son`un σ_Ĉ = {s['sigC']:.5f} ⇒ o doğrudan c(son) = "
          f"{pl:.4f} ; ölçülen {s['c_WX']:.4f} ⇒ artık "
          f"{s['c_WX']-pl:+.4f} = **{(s['c_WX']-pl)/s['s_tot']:+.2f}σ**")
    print("     ⇒ λ ekseninin kuralı fazlayı KAPATMIYOR, BÜYÜTÜYOR.")

    print("\n  (K3-d) KORELATÖR DÜZEYİ: `son`un bacak kırpma aktarımları")
    tab = {}
    for g in ("Hkeskin", "L085", "L070", "L060", "son"):
        p2 = SCR / f"169/K2b_{g}.json"
        if not p2.exists():
            continue
        o = json.loads(p2.read_text())["ortalama"]
        tab[g] = {n: float(np.mean([o[a] for a in o if NB[a] == n]))
                  for n in (1, 2, 3)}
    print(f"     {'gaz':9s} {'1 bacak':>9s} {'2 bacak':>9s} {'3 bacak':>9s} "
          f"{'(1b)^4':>8s} {'ölçülen c':>10s}")
    for g, v in tab.items():
        print(f"     {g:9s} {v[1]:9.4f} {v[2]:9.4f} {v[3]:9.4f} "
              f"{v[1]**4:8.4f} {K0[g]['c_WX']:10.4f}")
    print(f"     {'Gauss':9s} {np.sqrt(2/np.pi):9.4f} {2/np.pi:9.4f} "
          f"{(2/np.pi)**1.5:9.4f} {C_HIP:8.4f}")

    print("\n  (K3-e) ÖLÇÜLEN 3-BACAK AKTARIMININ 4-BACAĞA UZATIMI "
          "(ρ₃^{4/3}) ile c ÖNGÖRÜSÜ")
    print(f"     {'gaz':9s} {'ρ₃':>8s} {'ρ₃^{4/3}':>9s} {'ölçülen c':>10s} "
          f"{'ρ₃^{4/3}/c':>10s} | {'Hk-göreli öngörü':>17s} {'ölçüm':>8s} "
          f"{'artık σ':>9s}")
    e = {}
    if "Hkeskin" in tab:
        base3, basec = tab["Hkeskin"][3], K0["Hkeskin"]["c_WX"]
        for g, v in tab.items():
            ext = v[3] ** (4.0 / 3.0)
            pr = basec * (v[3] / base3) ** (4.0 / 3.0)
            sg = float(np.hypot(K0[g]["s_tot"], K0["Hkeskin"]["s_tot"]))
            z = (K0[g]["c_WX"] - pr) / sg if g != "Hkeskin" else 0.0
            e[g] = dict(rho3=v[3], ext=ext, c=K0[g]["c_WX"], pred=pr, z=z)
            print(f"     {g:9s} {v[3]:8.4f} {ext:9.4f} "
                  f"{K0[g]['c_WX']:10.4f} {ext/K0[g]['c_WX']:10.4f} | "
                  f"{pr:17.4f} {K0[g]['c_WX']:8.4f} {z:+8.2f}σ")
        print("     ⇒ `son` satırı H-D2'nin NİCEL sınavıdır (aynı λ, "
              "yalnız faz kilidi farklı); λ satırları kontroldür.")

    print("\n  (K3-f) BÜTÜN λ DİZİSİNİN TOPLU TABLOSU (ΠT kaynakları "
          "K1_T.json / ONKAYIT_*.json)")
    ok60 = json.load(open(SCR / "170/ONKAYIT_L060.json"))
    p115 = SCR / "170/ONKAYIT_L115.json"
    ok115 = json.load(open(p115)) if p115.exists() else None
    print(f"     {'gaz':9s} {'λ':>5s} {'σ_Ĉ':>8s} {'Π_b T(σ_b)':>11s} "
          f"{'ΠT/ΠT₀':>8s} | {'KALİB':>7s} {'W_X':>7s} {'c_WX':>7s} "
          f"{'±σ_tot':>7s}")
    sira = [("L115", 1.15), ("Hkeskin", 1.00), ("L085", 0.85),
            ("L070", 0.70), ("L060", 0.60), ("son", None)]
    toplu = []
    for g, l in sira:
        if g not in K0 or "c_WX" not in K0[g]:
            continue
        if g in KT:
            pt = KT[g]["PiT"]
        elif g == "L060":
            pt = ok60["ongoru"]["P1a"]
        elif g == "L115" and ok115:
            pt = (ok115["ongoru"]["Q1"] / h["c_WX"]) * th["PiT"]
        else:
            pt = float("nan")
        k = K0[g]
        toplu.append(dict(gaz=g, lam=l, sigC=k["sigC"], PiT=pt,
                          KALIB=k["KALIB"], WX=k["W_X"], c=k["c_WX"],
                          s=k["s_tot"]))
        print(f"     {g:9s} {('—' if l is None else '%.2f'%l):>5s} "
              f"{k['sigC']:8.5f} {pt:11.5f} {pt/th['PiT']:8.4f} | "
              f"{k['KALIB']:7.4f} {k['W_X']:7.4f} {k['c_WX']:7.4f} "
              f"{k['s_tot']:7.4f}")

    out = dict(toplu=toplu, fazla=d, s_fazla=sd, z_fazla=d / sd, PiT_oran=rt,
               c_pred_T=pred, z_T=(s["c_WX"] - pred) / sd,
               c_pred_lam=pl, z_lam=(s["c_WX"] - pl) / s["s_tot"],
               bacak=tab, uzatma=e, dogru=[float(p[0]), float(p[1])])
    json.dump(out, open(SCR / "170/K3.json", "w"), indent=1, default=float)
    print(f"\n-> {SCR}/170/K3.json")


if __name__ == "__main__":
    main()
