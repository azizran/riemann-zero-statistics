# -*- coding: utf-8 -*-
"""
180a — K0: DONMUŞ ÖN-KAYIT (NEDENSEL DEFTER)
=============================================
ÖN-MÜHÜR. Bu betik HİÇBİR yüzleşme niceliği hesaplanmadan ÖNCE koşar;
kendi sha256'sını ve `date` zaman damgasını `180/ONKAYIT_K0.json`'a yazar
ve dosya bir daha DEĞİŞTİRİLMEZ. Sonraki her betik (180b…) bu dosyayı
okur, sha'sını doğrular ve yalnız burada dondurulmuş formülleri uygular.

NE DONDURULUR
-------------
(1) Vekil ailelerinin TANIMI ve inşa kapıları (VS = gerçek-zarflı karışık).
(2) Ayrışım formülleri — İKİ KONVANSİYON (doğrusal VE log), ikisi de
    raporlanır; hiçbiri "seçilmez".
(3) Esas para birimi HAM defter: {KALİB(M), Q_E, ρ_E, Q_X, ρ_X}.
    Oran dili {g_E, g_X, θ, g_cal} YALNIZ EK olarak raporlanır.
(4) H-180a (çifte-sayım hakemi) karar kuralı ve İKİ DAL YORUMU — sayısal
    eşiklerle, ölçümden önce.
(5) Tablo şablonları (rapor bunları birebir doldurur).
(6) ÇAPALAR: hepsi ölçüm önbelleğinden okunur (yeni hesap değil).

BU BETİKTE HESAPLANMAYAN (ve hesaplanmaması ŞART olan):
    ZARF payı, KİLİT payı, bunların oranları, VS ailesinin herhangi bir
    defter satırı, VF ailesinin 169-K2b kırpma aktarımı. Bunların hepsi
    180b/180c/180d'de, bu dosya diske yazıldıktan SONRA ölçülür.

Kullanım: 180a_onkayit.py
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
for _p in ("172_configs", "167_configs", "164_configs"):
    sys.path.insert(0, str(QM / _p))
B172 = importlib.import_module("172b_gE_yasasi")
I164 = importlib.import_module("164_insa")

SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S180 = SCR / "180"
S167 = SCR / "167"
S169 = SCR / "169"
BU = Path(__file__).resolve()
TWO_PI = 2 * np.pi

HAM = ("M", "Q_E", "rho_E", "Q_X", "rho_X")        # ESAS para birimi
ORAN = ("g_E", "g_X", "theta", "g_cal")            # YALNIZ EK


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def defter(g):
    """Bir gazın HAM + ORAN defter satırı (172b makinesi, aynen)."""
    v = B172.oku(g)
    E = B172.kusur(v["varE_mod"], v["varE_olc"], v["varE_res"], v["korE"],
                   v["gE"])
    X = B172.kusur(v["varX_mod"], v["varX_olc"], v["varX_res"], v["korX"],
                   v["gX"])
    return dict(M=v["KAL"], Q_E=E["Q"], rho_E=E["rho"], pi_E=E["pi"],
                mu2_E=E["mu2"], r_E=E["r"], Q_X=X["Q"], rho_X=X["rho"],
                pi_X=X["pi"], g_E=v["gE"], g_X=v["gX"], theta=v["th"],
                g_cal=v["gcal"], W_X=v["WX"], N=v["N"], nline=v["nline"])


# =====================================================================
# 1. VEKİL AİLELERİNİN TANIMI
# =====================================================================
# VF (zaten kurulu, 176/177): Hkeskin ZARFI (a_q = 1/(πk√q), τ_q ≤ 1.00),
#     φ_q ~ U(0,2π), tohum 1,2,3,4  →  VF1..VF4.
# VS (180'de kurulur): GERÇEK gazın (son) ölçülen çizgi-gücü profiliyle
#     ölçeklenmiş zarf; aynı ω_q, aynı çözücü, φ_q ~ U(0,2π), tohum 1,2.
#
# ZARF TANIMI (dondurulmuş, parametresiz):
#     r(τ) := R_bant(son, τ_eff) / R_bant(Hkeskin, τ_eff)
#     — 167 bant defterinden, ÖLÇÜLEN çizgi gücünün nominale oranı;
#       iki gazın aynı bantlarındaki oranı saf sayıdır.
#     A_q(VS) := a_q · r(τ_q),  τ_q = ω_q / L_hedef
#     r(τ) τ_eff ızgarasında DOĞRUSAL interpolasyon; ölçülen aralığın
#     dışında SABİT TUTMA (clamp). Bu seçim burada donuyor.
#
# BİLİNEN SINIR (ölçümden önce yazıldı, süslenmeyecek):
#     R_bant'ın kendisi 176-HÜKÜM(vi)'ye göre kilit tarafından da
#     beslenir (karışık vekillerde 0.84–1.01'e düşüyor). Dolayısıyla
#     r(τ) SAF bir zarf niceliği değildir; "gerçeğin ölçülen çizgi-gücü
#     profilinin nominal merdivene taşınması"dır. ZARF PAYI bu tanıma
#     GÖREdir ve raporda böyle adlandırılacaktır. Bu, çapa-BAĞIMSIZLIĞI
#     bozmaz (hiçbir üçüncü referans gaza başvurulmuyor), ama "saf zarf"
#     iddiasını sınırlar.
#
# İNŞA KAPILARI (VS, 176-F0 ile aynı + biri yeni):
#   V1  sha256(A_VS1) == sha256(A_VS2)   (genlik dizisi bit-bit AYNI)
#   V2  maks|A_VS − a_Hk·r(τ)| == 0.0    (tarif birebir uygulandı)
#   V3  φ≡0 sağlaması: 0.0 / 0.0
#   V4  ilk-kök hücre 300000 / 300000
#   V5  maks|F| ≤ 1e−8, aşan tekne 0
#   V6  sıralılık TAM
#   Biri tutmazsa gaz ÖLÇÜLMEZ.
#
# ÖLÇÜM ZİNCİRİ: 167_olcum.kos(ad, 0.40, 0.95, 0, 0) + 176c defteri —
#   VF ile BİREBİR aynı çağrı. (KALİB_u2 b_nom'a bağlı değildir; b_nom
#   yalnız R_bant'a girer, o da yalnız SAĞLAMA olarak kullanılacaktır.)
#   SAĞLAMA (kapı değil, tanı): R_bant(VS)/R_bant(VF) ≈ r(τ) beklenir.

# =====================================================================
# 2. AYRIŞIM FORMÜLLERİ — İKİ KONVANSİYON, İKİSİ DE RAPORLANIR
# =====================================================================
FORMUL = {
 "aile": {
   "son": "gerçek ζ gazı (son 300000 sıfır)",
   "Hk": "Hkeskin — keskin ikiz (nominal zarf, φ≡0)",
   "VF": "Hk-zarflı karışık vekiller VF1..VF4 (⟨·⟩_VF = aritmetik ort.)",
   "VS": "son-zarflı karışık vekiller VS1,VS2 (⟨·⟩_VS = aritmetik ort.)",
 },
 "DOGRUSAL": {
   "delta":  "ΔN := N(son) − N(Hk)",
   "zarf":   "ZARF_N := ⟨N⟩_VS − ⟨N⟩_VF        [ÇAPASIZ, doğrudan ölçülür]",
   "kilit":  "KİLİT_N := ΔN − ZARF_N            [ARTIK: kilit + zarf×kilit]",
   "pay":    "p_zarf := ZARF_N/ΔN ,  p_kilit := KİLİT_N/ΔN ,  toplam ≡ 1",
 },
 "LOG": {
   "delta":  "ΔlogN := log N(son) − log N(Hk)",
   "zarf":   "ZARFlog_N := ⟨log N⟩_VS − ⟨log N⟩_VF",
   "kilit":  "KİLİTlog_N := ΔlogN − ZARFlog_N",
   "pay":    "p_zarf^log := ZARFlog/ΔlogN ,  p_kilit^log := 1 − p_zarf^log",
 },
 "esas_para_birimi": list(HAM),
 "ek_para_birimi": list(ORAN),
 "hata_cubugu": {
   "tanim": "se(⟨N⟩_A) = sd(N_A, ddof=1)/√n_A ; n_VS=2 ⇒ se=|N1−N2|/2",
   "zarf":  "se(ZARF_N) = √(se_VS² + se_VF²)",
   "kilit": "se(KİLİT_N) = se(ZARF_N)  (ΔN tohumsuz, tek gaz çifti)",
   "pay":   "se(p_zarf) = se(ZARF_N)/|ΔN|",
 },
 "isaret_kurali":
   "ZARF ve ΔN aynı işaretli değilse pay > 1 veya < 0 çıkar; bu bir "
   "hata değil ÖLÇÜMDÜR ve aynen yazılır (176'nın g_E dersi).",
}

# =====================================================================
# 3. H-180a — ÇİFTE-SAYIM HAKEMİ: KARAR KURALI ve İKİ DAL
# =====================================================================
GAUSS3 = float((2.0 / np.pi) ** 1.5)      # 3 bacak Gauss = (2/π)^{3/2}
GAUSS1 = float(np.sqrt(2.0 / np.pi))
GAUSS2 = float(2.0 / np.pi)

H180A = {
 "olculen": "169_k2b.main(g, 0.40, 0.95) → ortalama['111_hepsi'] =: ρ₃(g) "
            "(3 bacak kırpma aktarımı, lo∈[0.52,0.68] bantlarının ort.)",
 "gazlar": ["VF1", "VF2", "VF3", "VF4"],
 "capalar": {"rho3_Hkeskin": 0.5284, "rho3_son": 0.5356,
             "GAUSS3": GAUSS3, "GAUSS2": GAUSS2, "GAUSS1": GAUSS1},
 "istatistik": {
   "rho3bar_VF": "⟨ρ₃⟩ over VF1..VF4",
   "se_VF": "sd(ρ₃, ddof=1)/√4",
   "C": "ÇÖKÜŞ KESRİ  C := (ρ₃(Hk) − ⟨ρ₃⟩_VF) / (ρ₃(Hk) − GAUSS3)",
   "z_G": "z_G := (⟨ρ₃⟩_VF − GAUSS3)/se_VF",
 },
 "DAL_A": {
   "kosul": "C ≥ 0.70 VE |z_G| ≤ 3",
   "hukum": "Karışık vekilde kırpma aktarımı Gauss'a ÇÖKTÜ ⇒ ölçülen "
            "kırpma fazlası KİLİDİN BİR GÖRÜNÜMÜdür; %93 kırpma-kapanışı "
            "ile %91.4 kilit defteri AYNI mekanizmanın iki defteridir; "
            "ÇİFTE SAYIM YOKTUR (tek nedensel zincir: kilit → kırpma → ΔM).",
 },
 "DAL_B": {
   "kosul": "C ≤ 0.30",
   "hukum": "Kırpma fazlası kilit sökülünce AYAKTA KALDI ⇒ İKİ AYRI "
            "BİLEŞEN vardır; %184 gerçek çifte sayımdır ve paylar aynı "
            "defterde yeniden bölünmelidir.",
 },
 "HUKUMSUZ": {
   "kosul": "0.30 < C < 0.70, ya da (C ≥ 0.70 ama |z_G| > 3), ya da "
            "VF tohumları arasında ρ₃ saçılımı |ρ₃(Hk) − GAUSS3|'ün "
            "yarısını aşarsa (çözünürlük yok).",
   "hukum": "HÜKÜMSÜZ yazılır; kurtarma yapılmaz, süslenmez.",
 },
 "ek_kayit": "1 ve 2 bacak oranları (100_E, 010_Xa, 001_Xb, 011_XaXb, "
             "110_EXa) da tabloya yazılır; karar YALNIZ ρ₃'e bağlıdır.",
}

# =====================================================================
# 4. TABLO ŞABLONLARI (rapor bunları birebir doldurur)
# =====================================================================
SABLON = {
 "T1_insa": ["gaz", "tohum", "sha256(A)", "V1", "V2", "V3", "V4(hücre)",
             "V5(maks|F|)", "V6(sıralılık)", "σ_ds", "L", "süre"],
 "T2_defter": ["gaz", "M", "Q_E", "ρ_E", "Q_X", "ρ_X",
               "| g_E", "g_X", "θ", "g_cal"],
 "T3_ayrisim_dogrusal": ["nicelik N", "N(son)", "N(Hk)", "ΔN",
                         "⟨N⟩_VS", "⟨N⟩_VF", "ZARF_N ± se",
                         "KİLİT_N ± se", "p_zarf ± se", "p_kilit ± se"],
 "T4_ayrisim_log": ["nicelik N", "ΔlogN", "ZARFlog ± se", "KİLİTlog ± se",
                    "p_zarf^log ± se", "p_kilit^log ± se"],
 "T5_H180a": ["gaz", "ρ₁(E)", "ρ₁(Xa)", "ρ₁(Xb)", "ρ₂(XaXb)", "ρ₂(EXa)",
              "ρ₃(hepsi)"],
 "T6_eski_capa": ["kestirici", "eski (çapa-bağıl) değer", "konvansiyon",
                  "180'in çapasız karşılığı", "fark"],
 "T7_K4_capraz": ["nicelik", "176 karıştırma çöküşü (Hk zarfı)",
                  "180 ΔN(son↔Hk)", "ZARF payı", "tutarlı mı"],
}

# =====================================================================
# 5. ESKİ ÇAPA-BAĞIL SAYILAR (179/S4 doğrulamasından, aynen)
# =====================================================================
ESKI = {
 "f_b_tau070_080": 0.3975, "f_b_tau080_095": 0.3871,
 "f_gE_dogrusal": 0.3828, "f_gE_log": 0.3903, "f_mu2E": 0.3968,
 "f_theta_G1_dogrusal": -0.1685, "f_theta_G1_log": -0.1495,
 "f_theta_G2_dogrusal": -0.1430, "f_toplamDC_dogrusal": 0.4957,
 "mansetkesim_pay": 0.0857, "manset_kilit_pay": 0.9143,
 "manset_kesim_mutlak": 0.0049292, "manset_kilit_mutlak": 0.0526124,
 "DlogM_son_Hk": 0.057542,
 "not": "Hepsi (Hkeskin, HA4) ÇAPA üçlüsüne ve konvansiyon seçimine "
        "bağıldır — 179/S4 §4.4. 180 bunları YENİDEN ÜRETMEZ, yalnız "
        "çapasız sayının yanına koyar.",
}


def main():
    S180.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    dmg = subprocess.run(["date"], capture_output=True, text=True).stdout.strip()
    print("=" * 74, flush=True)
    print("180a — K0: DONMUŞ ÖN-KAYIT (NEDENSEL DEFTER)", flush=True)
    print("=" * 74, flush=True)
    print(f"  zaman (`date`): {dmg}", flush=True)

    # ---------------- ÇAPALAR: önbellekten okunur ---------------------
    capa = {}
    for g in ("son", "Hkeskin", "HA4", "VF1", "VF2", "VF3", "VF4"):
        capa[g] = defter(g)
        d = capa[g]
        print(f"  ÇAPA {g:8s} M={d['M']:.6f}  Q_E={d['Q_E']:.6f}  "
              f"ρ_E={d['rho_E']:.6f}  Q_X={d['Q_X']:.6f}  "
              f"ρ_X={d['rho_X']:.6f}  | g_E={d['g_E']:.6f} θ={d['theta']:.6f}",
              flush=True)

    # K2b çapaları (169 önbelleği)
    k2b = {}
    for g in ("Hkeskin", "son", "L060", "L070", "L085"):
        p = S169 / f"K2b_{g}.json"
        if p.exists():
            k2b[g] = json.load(open(p))["ortalama"]
    print(f"  ÇAPA ρ₃: Hkeskin={k2b['Hkeskin']['111_hepsi']:.4f}  "
          f"son={k2b['son']['111_hepsi']:.4f}  GAUSS3={GAUSS3:.6f}", flush=True)

    # ---------------- ZARF PROFİLİ r(τ) — inşa girdisi ----------------
    bs = json.load(open(S167 / "C_son.json"))["bant"]
    bh = json.load(open(S167 / "C_Hkeskin.json"))["bant"]
    ts = np.array([b["tau_eff"] for b in bs if b.get("olculdu")])
    rs = np.array([b["R_bant"] for b in bs if b.get("olculdu")])
    th_ = np.array([b["tau_eff"] for b in bh if b.get("olculdu")])
    rh = np.array([b["R_bant"] for b in bh if b.get("olculdu")])
    assert len(ts) == len(th_) and np.allclose(ts, th_, atol=2e-3), \
        "bant ızgaraları uyuşmuyor"
    tg = 0.5 * (ts + th_)
    rr = rs / rh
    print("\n  ZARF PROFİLİ r(τ) = R_bant(son)/R_bant(Hkeskin):", flush=True)
    for a, b, c, e in zip(tg, rs, rh, rr):
        print(f"    τ_eff={a:.4f}  R_son={b:.5f}  R_Hk={c:.5f}  "
              f"r={e:.5f}", flush=True)

    # A_q(VS) dizisi — bit-bit dondurulur
    d = np.load(QM / "128_odl_zeros6_2e6_zeros.npz")
    Z = np.sort(np.asarray(d["zeros"], dtype=float))
    zr = Z[len(Z) - 300000:]
    L_hedef = float(np.log(0.5 * (float(zr[0]) + float(zr[-1])) / TWO_PI))
    om, a_h, _, _ = I164.merdiven(L_hedef, None, None)
    tau = om / L_hedef
    r_q = np.interp(tau, tg, rr)          # dışarıda SABİT TUTMA (clamp)
    A_vs = a_h * r_q
    np.save(S180 / "A_sonzarf.npy", A_vs)
    np.save(S180 / "r_tau_profil.npy", np.vstack([tg, rr]))
    sha_A = hashlib.sha256(np.ascontiguousarray(A_vs).tobytes()).hexdigest()
    sha_Ah = hashlib.sha256(np.ascontiguousarray(a_h).tobytes()).hexdigest()
    print(f"\n  merdiven: {len(om)} çizgi, L_hedef={L_hedef:.9f}", flush=True)
    print(f"  sha256(A_VS)      = {sha_A}", flush=True)
    print(f"  sha256(A_Hkeskin) = {sha_Ah}", flush=True)
    print(f"  Σ|A_VS| = {np.abs(A_vs).sum():.9f}   (Hk: "
          f"{np.abs(a_h).sum():.9f})", flush=True)
    print(f"  ΣA_VS²  = {np.sum(A_vs**2):.9f}   (Hk: "
          f"{np.sum(a_h**2):.9f})   oran = "
          f"{np.sum(A_vs**2)/np.sum(a_h**2):.6f}", flush=True)
    print(f"  r_q: min={r_q.min():.5f} maks={r_q.max():.5f}  "
          f"τ>0.7766 clamp'li çizgi sayısı = {int((tau > tg[-1]).sum())} "
          f"/ {len(tau)}  (ΣA² payı = "
          f"{float(np.sum(A_vs[tau>tg[-1]]**2)/np.sum(A_vs**2)):.4f})",
          flush=True)

    # clamp duyarlılığı — TANI (karar niceliği değil)
    r_lin = np.interp(tau, tg, rr)
    sl = float(np.polyfit(tg[-4:], rr[-4:], 1)[0])
    r_ext = np.where(tau > tg[-1], rr[-1] + sl * (tau - tg[-1]), r_lin)
    r_ext = np.where(tau < tg[0], rr[0], r_ext)
    print(f"  [tanı] doğrusal-uzatmalı alternatif: ΣA² oranı "
          f"{float(np.sum((a_h*r_ext)**2)/np.sum(A_vs**2)):.6f} "
          f"(clamp'e göre)", flush=True)

    rec = dict(
        gorev=180, ad="NEDENSEL DEFTER — K0 ÖN-KAYIT",
        zaman=dmg, zaman_unix=time.time(), betik=BU.name, sha256=sha(BU),
        kalem="KALEM_NEDENSEL_DEFTER_07EYL2026.md",
        sha_kalem=sha(QM / "KALEM_NEDENSEL_DEFTER_07EYL2026.md"),
        formul=FORMUL, H180a=H180A, sablon=SABLON, eski_capa=ESKI,
        ham=list(HAM), oran=list(ORAN),
        capa_defter=capa, capa_k2b=k2b,
        zarf=dict(tanim="r(τ)=R_bant(son)/R_bant(Hkeskin), τ_eff ızgarasında"
                        " doğrusal interpolasyon, dışarıda SABİT TUTMA",
                  tau_eff=tg.tolist(), R_son=rs.tolist(), R_Hk=rh.tolist(),
                  r=rr.tolist(), L_hedef=L_hedef, nline=int(len(om)),
                  sha_A_VS=sha_A, sha_A_Hkeskin=sha_Ah,
                  sumA2_VS=float(np.sum(A_vs ** 2)),
                  sumA2_Hk=float(np.sum(a_h ** 2)),
                  sumabsA_VS=float(np.abs(A_vs).sum()),
                  clamp_cizgi=int((tau > tg[-1]).sum()),
                  clamp_A2_payi=float(np.sum(A_vs[tau > tg[-1]] ** 2)
                                      / np.sum(A_vs ** 2)),
                  clamp_duyarlilik_A2=float(np.sum((a_h * r_ext) ** 2)
                                            / np.sum(A_vs ** 2))),
        insa_kapilari=["V1 sha(A_VS1)==sha(A_VS2)",
                       "V2 maks|A_VS − a_Hk·r(τ)| == 0.0",
                       "V3 φ≡0 sağlaması 0.0/0.0",
                       "V4 ilk-kök hücre 300000/300000",
                       "V5 maks|F| ≤ 1e−8, aşan 0", "V6 sıralılık TAM"],
        tohumlar=dict(VS=[1, 2], VF=[1, 2, 3, 4]),
        sure_s=time.time() - t0)
    OUT = S180 / "ONKAYIT_K0.json"
    if OUT.exists():
        raise SystemExit(f"ÖN-KAYIT ZATEN VAR: {OUT} — üzerine yazılmaz.")
    OUT.write_text(json.dumps(rec, indent=1, ensure_ascii=False,
                              default=float))
    print(f"\n-> {OUT}", flush=True)
    print(f"   betik sha256 = {rec['sha256']}", flush=True)
    print(f"   KALEM  sha256 = {rec['sha_kalem']}", flush=True)
    print("   ÖN-KAYIT KAPANDI. Bundan sonra yalnız ÖLÇÜM.", flush=True)


if __name__ == "__main__":
    main()
