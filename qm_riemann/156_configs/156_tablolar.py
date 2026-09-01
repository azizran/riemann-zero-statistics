"""
156 — rapor tablolarını JSON'lardan üretir (elle sayı girilmez).

Kullanım: 156_tablolar.py > tablolar.md
"""
import json
from pathlib import Path

import numpy as np

SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/156")
NEFF_TABAN = 3000          # 300k'nın %1'i: altındaki uydurma "birkaç nokta"
GAZ = [("gercek", "gerçek"), ("A4", "A4"), ("keskin", "keskin")]
GOREV = (0.585, 0.66, 0.74)


def yuk(dosya):
    return json.load(open(SCR / dosya))


def en_iyi2(b, taban=None):
    ok = [p for p in b.get("sirt", []) if p["artik"] < 0.02
          and (taban is None or p["neff"] >= taban)]
    return min(ok, key=lambda p: p["dRe"]) if ok else None


def main():
    D = {a: yuk(f"k3_{a}_std.json") for a, _ in GAZ}

    print("### T0 — kanal istatistikleri ve PAY YAPISI\n")
    print("| büyüklük | gerçek (son-300k) | sentetik A4 | sentetik keskin |")
    print("|---|---|---|---|")
    for ad, k in (("σ_ds² (nokta)", "s_ds"), ("σ_lad² (nokta)", "s_lad"),
                  ("σ_η² (nokta)", "s_eta"), ("c₁(η)", "c1")):
        print(f"| {ad} | " + " | ".join(f"{D[a][k]:.5f}" for a, _ in GAZ) + " |")
    for ad, k in (("σΔ² = Var(X_tam)", "tam"), ("Var(X_lad)", "lad"),
                  ("Var(X_eta)", "eta"), ("Var(X_drift)", "dri")):
        print(f"| {ad} | " + " | ".join(f"{D[a]['var'][k]:.5f}"
                                        for a, _ in GAZ) + " |")
    for k, ad in (("lad", "Cov(X_lad, X_tam)"), ("eta", "Cov(X_eta, X_tam)")):
        print(f"| {ad} | " + " | ".join(f"{D[a]['kov'][k]:.5f}"
                                        for a, _ in GAZ) + " |")
    for k, ad in (("lad", "**pay** Cov(X_lad,X_tam)/σΔ²"),
                  ("eta", "**pay** Cov(X_eta,X_tam)/σΔ²"),
                  ("dri", "**pay** Cov(X_drift,X_tam)/σΔ²")):
        print(f"| {ad} | " + " | ".join(
            f"{D[a]['kov'][k]/D[a]['var']['tam']:+.4f}" for a, _ in GAZ) + " |")
    for k, ad in (("lad", "korel(X_lad, X_tam)"),
                  ("eta", "korel(X_eta, X_tam)")):
        print(f"| {ad} | " + " | ".join(
            f"{D[a]['kov'][k]/np.sqrt(D[a]['var'][k]*D[a]['var']['tam']):+.4f}"
            for a, _ in GAZ) + " |")
    print("| ölçek s_lad = σΔ²/Cov(X_lad,X_tam) | " + " | ".join(
        f"{D[a]['olcek']['lad']:.3f}" for a, _ in GAZ) + " |")
    print("| ölçek s_eta = σΔ²/Cov(X_eta,X_tam) | " + " | ".join(
        f"{D[a]['olcek']['eta']:.3f}" for a, _ in GAZ) + " |")

    for a, ad in GAZ:
        print(f"\n\n### T1-{ad} — tek-kanallı eşleşmeler\n")
        print("| τ̄ | \\|Γ\\| | φ_Γ | ReΓ | kanal | R_ham | R̃ | faz artığı | "
              "ReM | \\|ΔRe\\| | n_eff | hüküm |")
        print("|---|---|---|---|---|---|---|---|---|---|---|---|")
        for b in D[a]["bantlar"]:
            if not b.get("olculdu"):
                print(f"| {b['tau']:.4f} | — | — | — | — | — | — | — | — | "
                      "— | — | ölçülemedi |")
                continue
            ilk = True
            for k in ("tam", "lad", "eta"):
                c = b["kanal"][k]
                tut = c["artik"] < 0.02
                hk = ("faz üretilemiyor" if not tut
                      else ("✓" if c["dRe"] < 0.06 else "✗ Re"))
                if tut and c["neff"] < NEFF_TABAN:
                    hk += " ⚠n_eff"
                sol = (f"| {b['tau']:.4f} | {b['absG']:.3f}"
                       f"{' ‡' if b['patlak'] else ''} | {b['phi']:+.4f} | "
                       f"{b['Gre']:+.4f} " if ilk else "|  |  |  |  ")
                ilk = False
                print(sol + f"| {k} | {c['Rham']:+.3f} | {c['Rn']:+.3f} | "
                      f"{c['artik']:.4f} | {c['ReM']:+.4f} | {c['dRe']:.4f} | "
                      f"{c['neff']:.0f} | {hk} |")
        print("\n‡ = |Γ|>1.5, bant sayısal olarak patlak (152/153 kuralı); "
              "hüküm çıkarılmaz. ⚠n_eff = ağırlık "
              f"{NEFF_TABAN}'den az etkin noktaya çökmüş (uydurma birkaç "
              "noktaya dayanıyor).")

    print("\n\n### T2 — iki-kanallı sırt (w = e^{−A(R_l·X_lad + R_e·X_eta)})\n")
    print("| gaz | τ̄ | sırt nokta | serbest en iyi (R_l,R_e)_ham | \\|ΔRe\\| |"
          f" n_eff | n_eff≥{NEFF_TABAN} kısıtlı (R_l,R_e)_ham | \\|ΔRe\\| | "
          "n_eff |")
    print("|---|---|---|---|---|---|---|---|---|")
    for a, ad in GAZ:
        for b in D[a]["bantlar"]:
            if not b.get("olculdu"):
                continue
            s = [p for p in b.get("sirt", []) if p["artik"] < 0.02]
            f = (lambda p: (f"({p['Rl_ham']:+.3f}, {p['Re_ham']:+.3f}) | "
                            f"{p['dRe']:.4f} | {p['neff']:.0f}")
                 if p else "— | — | —")
            print(f"| {ad}{' ‡' if b['patlak'] else ''} | {b['tau']:.4f} | "
                  f"{len(s)} | {f(en_iyi2(b))} | {f(en_iyi2(b, NEFF_TABAN))} |")

    print("\nSırttaki en iyi nokta SIĞ bir minimumdur; kabul edilebilir BÖLGE "
          f"(faz artığı<0.02, |ΔRe|<0.06, n_eff≥{NEFF_TABAN}) daha bilgilidir:"
          "\n")
    print("| gaz | τ̄ | kabul edilebilir R_l (ham) | R_e (ham) | nokta |")
    print("|---|---|---|---|---|")
    for a, ad in GAZ:
        for b in D[a]["bantlar"]:
            if not b.get("olculdu"):
                continue
            ok = [p for p in b.get("sirt", []) if p["artik"] < 0.02
                  and p["dRe"] < 0.06 and p["neff"] >= NEFF_TABAN]
            if not ok:
                print(f"| {ad}{' ‡' if b['patlak'] else ''} | {b['tau']:.4f} "
                      "| — | — | 0 |")
                continue
            rl = [p["Rl_ham"] for p in ok]; re = [p["Re_ham"] for p in ok]
            print(f"| {ad}{' ‡' if b['patlak'] else ''} | {b['tau']:.4f} | "
                  f"{min(rl):+.3f} … {max(rl):+.3f} | {min(re):+.3f} … "
                  f"{max(re):+.3f} | {len(ok)} |")

    # ---- taban dayanıklılığı ----
    print("\n\n### T3 — TABAN DAYANIKLILIĞI (gerçek gaz, regresyon tabanı "
          "0.46 / 0.52 / 0.58)\n")
    print("lad/eta ayrımı tabanın nereye çizildiğine bağlıdır; hüküm tabana "
          "bağlıysa fizik değil konvansiyondur.\n")
    T = {"0.46": yuk("k3_gercek_gorev_t0.46.json"),
         "0.52": D["gercek"],
         "0.58": yuk("k3_gercek_gorev_t0.58.json")}
    print("| τ̄ | taban | φ_Γ | R_tam | \\|ΔRe\\|_tam | R_lad | \\|ΔRe\\|_lad "
          "| R_eta | \\|ΔRe\\|_eta | 2K en iyi (R_l,R_e)_ham |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    for tau in GOREV:
        for tb in ("0.46", "0.52", "0.58"):
            b = next((x for x in T[tb]["bantlar"] if x.get("olculdu")
                      and abs(x["tau"] - tau) < 1e-6), None)
            if b is None:
                continue
            c = b["kanal"]
            def g(k, alan, fmt="{:+.3f}"):
                return ("—" if c[k]["artik"] >= 0.02
                        else fmt.format(c[k][alan]))
            e = en_iyi2(b, NEFF_TABAN)
            es = (f"({e['Rl_ham']:+.3f}, {e['Re_ham']:+.3f})" if e else "—")
            print(f"| {tau:.4f} | {tb} | {b['phi']:+.4f} | "
                  f"{g('tam','Rham')} | {g('tam','dRe','{:.4f}')} | "
                  f"{g('lad','Rham')} | {g('lad','dRe','{:.4f}')} | "
                  f"{g('eta','Rham')} | {g('eta','dRe','{:.4f}')} | {es} |")

    # ---- görev bantları özeti ----
    print("\n\n### T4 — görev bantları özeti (hangi kanal kazanıyor)\n")
    print("| τ̄ | gaz | R_tam | R_lad | R_eta | \\|ΔRe\\|_tam | \\|ΔRe\\|_lad |"
          " \\|ΔRe\\|_eta | en iyi kanal |")
    print("|---|---|---|---|---|---|---|---|---|")
    for tau in GOREV:
        for a, ad in GAZ:
            b = next((x for x in D[a]["bantlar"] if x.get("olculdu")
                      and abs(x["tau"] - tau) < 1e-6), None)
            if b is None:
                continue
            c = b["kanal"]
            gec = {k: (c[k]["dRe"] if c[k]["artik"] < 0.02 else None)
                   for k in ("tam", "lad", "eta")}
            uy = [k for k in gec if gec[k] is not None]
            iyi = min(uy, key=lambda k: gec[k]) if uy else "—"
            def g(k, alan, fmt="{:+.3f}"):
                return ("—" if c[k]["artik"] >= 0.02
                        else fmt.format(c[k][alan]))
            print(f"| {tau:.4f} | {ad}{' ‡' if b['patlak'] else ''} | "
                  f"{g('tam','Rham')} | {g('lad','Rham')} | {g('eta','Rham')} |"
                  f" {g('tam','dRe','{:.4f}')} | {g('lad','dRe','{:.4f}')} | "
                  f"{g('eta','dRe','{:.4f}')} | {iyi} |")

    # ---- faz kapasitesi ----
    print("\n\n### T4b — FAZ KAPASİTESİ: kanal TEK BAŞINA hangi fazı "
          "üretebilir?\n")
    print("Sarma-açılmış arg M(R̃) menzili (`156_kapasite.py`). İkinci menzil "
          f"n_eff ≥ {NEFF_TABAN} kısıtı altındadır: ağırlığı çökertmeden "
          "erişilebilen faz. Ölçülen φ_Γ menzilin dışındaysa o kanal fazı "
          "HİÇBİR R ile üretemez.\n")
    print("| gaz | τ̄ | φ_Γ | kanal | serbest menzil | kök? | n_eff≥ kısıtlı "
          "menzil | kök? | η tavanı / φ_Γ |")
    print("|---|---|---|---|---|---|---|---|---|")
    for k in yuk("kapasite.json"):
        pay = (f"{k['arg_hi_s']/k['phi']:.0%}"
               if k["kanal"] == "eta" and k["phi"] > 0 else "")
        print(f"| {k['gaz']} | {k['tau']:.4f} | {k['phi']:+.4f} | "
              f"{k['kanal']} | [{k['arg_lo']:+.3f}, {k['arg_hi']:+.3f}] | "
              f"{'var' if k['kok_var'] else '**YOK**'} | "
              f"[{k['arg_lo_s']:+.3f}, {k['arg_hi_s']:+.3f}] | "
              f"{'var' if k['kok_var_s'] else '**YOK**'} | {pay} |")

    # ---- nakil ----
    print("\n\n### T5 — NAKİL SINAVI: gerçekte ölçülen ham R_k → A4 gazına\n")
    N = yuk("nakil.json")
    print("`öngörü/φ_A4` = nakil sentetiğin ölçülen fazını ne kadar tutuyor "
          "(hedef 1.00). `öngörü/φ_gerçek` = bağlaşım EVRENSEL olsaydı "
          "sentetik/gerçek faz oranı ne çıkardı — bunu ölçülen `oran` "
          "sütunuyla karşılaştırın.\n")
    print("| τ̄ | φ_gerçek | φ_A4 (ölç) | ölçülen oran | kanal | "
          "R_ham(gerçek) | φ_öngörü(A4) | öngörü/φ_A4 | öngörü/φ_gerçek | "
          "n_eff |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    for s in N:
        ilk = True
        for k in ("tam", "lad", "eta"):
            n = s["nakil"].get(k)
            sol = (f"| {s['tau']:.4f}{' ‡' if s['patlak_s'] else ''} | "
                   f"{s['phi_g']:+.4f} | {s['phi_s']:+.4f} | "
                   f"{s['oran_olcu']:+.3f} " if ilk else "|  |  |  |  ")
            ilk = False
            if n is None:
                print(sol + f"| {k} | — | — | — | — | — |")
                continue
            print(sol + f"| {k} | {n['R_ham_g']:+.3f} | {n['phi_ong']:+.4f} | "
                  f"{n['phi_ong']/s['phi_s']:+.3f} | "
                  f"{n['phi_ong']/s['phi_g']:+.3f} | {n['neff']:.0f} |")

    # ---- ortak çift ----
    print("\n\n### T6 — ORTAK ÇİFT: tek bir (R_l,R_e) iki gazı birden "
          "açıklıyor mu?\n")
    try:
        O = yuk("ortak_cift.json")
    except FileNotFoundError:
        print("(ortak_cift.json yok)")
        return
    print("| τ̄ | ortak faz çözümü sayısı | en iyi (R_l,R_e)_ham | "
          "\\|ΔRe\\|_gerçek | \\|ΔRe\\|_A4 | n_eff(min) | hüküm |")
    print("|---|---|---|---|---|---|---|")
    for s in O:
        if not s.get("cozum"):
            print(f"| {s['tau']:.4f}{' ‡' if s.get('patlak') else ''} | 0 | — "
                  "| — | — | — | ortak faz çözümü yok |")
            continue
        e = s["en_iyi"]
        print(f"| {s['tau']:.4f}{' ‡' if s.get('patlak') else ''} | "
              f"{len(s['cozum'])} | ({e['Rl']:+.3f}, {e['Re']:+.3f}) | "
              f"{e['gercek']['dRe']:.4f} | {e['A4']['dRe']:.4f} | "
              f"{e['neff_min']:.0f} | {s['hukum']} |")


if __name__ == "__main__":
    main()
