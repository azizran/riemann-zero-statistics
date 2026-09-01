"""
156 — NAKİL SINAVI (transplant): ×1.4'ün kaynağı pay-yapısı mı?

Mantık: eğer sağkalım-ağırlık bağlaşımı kanal k'de yaşayan EVRENSEL bir
bağlaşım sabitiyse (R_k soğurmanın kendi özelliği, gazın değil), o zaman
GERÇEKTE ölçülen ham R_k, SENTETİK gazın kendi ampirik dağılımına
uygulandığında sentetiğin ÖLÇÜLEN fazını vermelidir. Vermiyorsa o kanal
evrensel değildir.

Bu, Gauss yaklaşımı KULLANMAZ: her iki gazda da ağırlık ortalaması
kendi 300k'lık ampirik örneği üzerinden alınır. Tek varsayım, R_k'nin
gazdan bağımsız olduğudur — hipotezin kendi iddiası.

Ayrıca: iki gazın faz-eşlenmiş (R_l, R_e) SIRTLARININ ham koordinatlarda
KESİŞİMİ aranır. Kesişim varsa, tek bir evrensel çift her iki gazın
fazını da açıklıyor demektir; o çiftin Re'yi de tutup tutmadığı ayrıca
sınanır.

Kullanım: 156_nakil.py
"""
import importlib
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
C3 = importlib.import_module("156_cekirdek")

HERE = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
OUT = SCR / "156"
TWO_PI = 2 * np.pi


def kanallar(ad):
    p = OUT / f"kanal_{ad}.npy"
    if p.exists():
        M = np.load(p)
        return dict(tam=M[0], lad=M[1], eta=M[2], dri=M[3])
    if ad == "gercek":
        d = np.load(HERE / "128_odl_zeros6_2e6_zeros.npz")
        Z = np.sort(np.asarray(d["zeros"], dtype=float))
        z = Z[len(Z) - 300000:]
    else:
        z = np.sort(np.load(SCR / f"154/z_{ad}.npy"))
    C = C3.zincir3(z)
    X = {k: C3._bond(C[v]) for k, v in
         (("tam", "ds"), ("lad", "lad"), ("eta", "eta"), ("dri", "drift"))}
    np.save(p, np.vstack([X["tam"], X["lad"], X["eta"], X["dri"]]))
    return X


def M_of(A, x_lin, ex):
    x = -A * x_lin
    w = np.exp(x - x.max())
    sw = w.sum()
    return (w @ ex) / sw, float(sw * sw / (w @ w))


def main():
    G = {a: json.load(open(OUT / f"k3_{a}_std.json"))
         for a in ("gercek", "A4")}
    X = {a: kanallar(a) for a in ("gercek", "A4")}
    ex = {}
    cikti = []

    bant_g = {b["tau"]: b for b in G["gercek"]["bantlar"] if b.get("olculdu")}
    bant_s = {b["tau"]: b for b in G["A4"]["bantlar"] if b.get("olculdu")}

    print("=== NAKİL SINAVI: gerçekte ölçülen ham R_k → sentetik A4 gazı ===")
    print("(φ_öng = nakil öngörüsü; φ_ölç = sentetiğin ölçülen fazı; "
          "hedef φ_öng ≈ φ_ölç)\n")
    for tau in sorted(set(bant_g) & set(bant_s)):
        bg, bs = bant_g[tau], bant_s[tau]
        A = TWO_PI * tau
        for a in ("gercek", "A4"):
            ex[a] = np.exp(-1j * A * X[a]["tam"])
        pat = "  [SENTETİK BANT PATLAK |Γ|>1.5]" if bs["patlak"] else ""
        print(f"τ={tau:.4f}   φ_gerçek={bg['phi']:+.4f}  "
              f"φ_A4(ölç)={bs['phi']:+.4f}  oran="
              f"{bs['phi']/bg['phi']:+.3f}{pat}")
        sat = dict(tau=tau, phi_g=bg["phi"], phi_s=bs["phi"],
                   oran_olcu=bs["phi"] / bg["phi"], patlak_s=bs["patlak"],
                   nakil={})
        for k in ("tam", "lad", "eta"):
            kg = bg["kanal"][k]
            if kg["artik"] >= 0.02:
                print(f"    {k:3s}: gerçekte faz üretilemedi → nakil YOK")
                sat["nakil"][k] = None
                continue
            Rham = kg["Rham"]
            Mv, neff = M_of(A, Rham * X["A4"][k], ex["A4"])
            php = float(np.angle(Mv))
            # ters yön: sentetikte ölçülen R → gerçek gaz
            ks = bs["kanal"][k]
            ters = None
            if ks["artik"] < 0.02:
                Mv2, ne2 = M_of(A, ks["Rham"] * X["gercek"][k], ex["gercek"])
                ters = dict(phi=float(np.angle(Mv2)), Re=float(Mv2.real),
                            neff=ne2)
            print(f"    {k:3s}: R_ham(gerçek)={Rham:+9.3f} → φ_öng="
                  f"{php:+.4f}  (ölç {bs['phi']:+.4f}, öng/ölç="
                  f"{php/bs['phi']:+.3f})  ReM={Mv.real:+.4f} "
                  f"(ReΓ {bs['Gre']:+.4f})  n_eff={neff:.0f}"
                  + (f"   | ters: R_ham(A4)={ks['Rham']:+.3f} → φ_öng"
                     f"={ters['phi']:+.4f} (gerçek ölç {bg['phi']:+.4f})"
                     if ters else ""))
            sat["nakil"][k] = dict(R_ham_g=Rham, phi_ong=php,
                                   Re_ong=float(Mv.real), neff=neff,
                                   ters=ters)
        # iki-kanallı: gerçeğin en iyi çifti → sentetik
        e2 = bg.get("en_iyi2")
        if e2:
            xl = e2["Rl_ham"] * X["A4"]["lad"] + e2["Re_ham"] * X["A4"]["eta"]
            Mv, neff = M_of(A, xl, ex["A4"])
            print(f"    2K : (R_l,R_e)_gerçek=({e2['Rl_ham']:+.3f},"
                  f"{e2['Re_ham']:+.3f}) → φ_öng={np.angle(Mv):+.4f}  "
                  f"ReM={Mv.real:+.4f}  n_eff={neff:.0f}")
            sat["nakil"]["2K"] = dict(Rl=e2["Rl_ham"], Re=e2["Re_ham"],
                                      phi_ong=float(np.angle(Mv)),
                                      Re_ong=float(Mv.real), neff=neff)
        # ---- sırt kesişimi (ham koordinat) ----
        kes = sirt_kesisim(bg, bs)
        if kes:
            print(f"    SIRT KESİŞİMİ (ham): R_l={kes['Rl']:+.4f}  "
                  f"R_e={kes['Re']:+.4f}")
            for a, b in (("gercek", bg), ("A4", bs)):
                xl = kes["Rl"] * X[a]["lad"] + kes["Re"] * X[a]["eta"]
                Mv, neff = M_of(A, xl, ex[a])
                print(f"        {a:7s}: φ={np.angle(Mv):+.4f} (ölç "
                      f"{b['phi']:+.4f}, artık {abs(np.angle(Mv)-b['phi']):.4f})"
                      f"  ReM={Mv.real:+.4f} (ReΓ {b['Gre']:+.4f}, |ΔRe|="
                      f"{abs(Mv.real-b['Gre']):.4f})  n_eff={neff:.0f}")
                kes[a] = dict(phi=float(np.angle(Mv)), Re=float(Mv.real),
                              neff=neff, dphi=float(abs(np.angle(Mv)-b['phi'])),
                              dRe=float(abs(Mv.real - b["Gre"])))
            sat["kesisim"] = kes
        else:
            print("    SIRT KESİŞİMİ: yok (iki sırt ham düzlemde kesişmiyor)")
            sat["kesisim"] = None
        print()
        cikti.append(sat)

    (OUT / "nakil.json").write_text(json.dumps(cikti, indent=1))
    print(f"-> {OUT / 'nakil.json'}")


def sirt_kesisim(bg, bs):
    """İki faz-eşlenmiş sırtın ham (R_l, R_e) düzlemindeki kesişimi.

    Sırtlar R_l'ye göre sıralı eğrilerdir; ortak R_l aralığında
    d(R_l) = R_e^gerçek(R_l) − R_e^A4(R_l) fonksiyonunun işaret değiştirdiği
    yerde doğrusal ara değerle kesişim bulunur.
    """
    def egri(b):
        p = [(q["Rl_ham"], q["Re_ham"]) for q in b.get("sirt", [])
             if q["artik"] < 0.02]
        return sorted(p)
    A_, B_ = egri(bg), egri(bs)
    if len(A_) < 2 or len(B_) < 2:
        return None
    lo = max(A_[0][0], B_[0][0]); hi = min(A_[-1][0], B_[-1][0])
    if hi <= lo:
        return None
    xs = np.linspace(lo, hi, 2001)
    ya = np.interp(xs, [p[0] for p in A_], [p[1] for p in A_])
    yb = np.interp(xs, [p[0] for p in B_], [p[1] for p in B_])
    d = ya - yb
    sg = np.where(np.sign(d[:-1]) * np.sign(d[1:]) < 0)[0]
    if len(sg) == 0:
        return None
    i = int(sg[0])
    t = d[i] / (d[i] - d[i + 1])
    Rl = float(xs[i] + t * (xs[i + 1] - xs[i]))
    Re = float(ya[i] + t * (ya[i + 1] - ya[i]))
    return dict(Rl=Rl, Re=Re, ortak_Rl=[float(lo), float(hi)])


if __name__ == "__main__":
    main()
