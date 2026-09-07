# -*- coding: utf-8 -*-
"""
181a — K0: n=4 TOHUM UZANTISININ DONMUŞ ÖN-KAYDI
=================================================
ÖN-MÜHÜR. Bu betik **VS3/VS4 doğmadan önce** koşar ve hiçbir yeni
yüzleşme niceliği hesaplamaz: yalnız (a) 180a'nın donmuş ön-kaydını ve
180f'nin K3 defterini OKUR, (b) n=4'te hata çubuğunun ne olacağını
ARİTMETİKLE önceden yazar, (c) hangi HÜKÜMSÜZ satırın hangi tohum
saçılımı eşiğinde OKUNUR hâle geleceğini sayısallaştırır.

GEVŞETME YOK — eşikler 180a'nındır:
  * ayrışım formülleri (doğrusal + log), hata çubuğu tanımı
    se(⟨N⟩_A)=sd(ddof=1)/√n_A, se(ZARF)=√(se_VS²+se_VF²),
    se(KİLİT)=se(ZARF), se(p)=se(ZARF)/|ΔN|  — 180a'dan AYNEN;
  * işaret kuralı (pay>1 veya <0 bir ölçümdür) — 180a'dan AYNEN;
  * inşa kapıları V1..V6 — 180a'dan AYNEN (ilk-kök 300000/300000,
    maks|F| ≤ 1e−8, sıralılık TAM, sha(A) = ön-kayıttaki).

TEK YENİ ŞEY (ve bu yüzden burada, sayılardan ÖNCE donuyor):
180a "okunur / HÜKÜMSÜZ" ayrımını 180 raporunda SÖZEL uyguladı
(§4.1). 181 n=4'e çıkacağı için bu ayrımın SAYISAL hâli gerekiyor.
Uydurmuyoruz: 180'in kendi partisyonunu birebir üreten en dar ölçüt
seçiliyor ve KALİBRASYON burada, veri gelmeden gösteriliyor:

    O1 (OKUNUR ÖLÇÜTÜ):  |ZARF_N| ≥ 2 · se(ZARF_N)
    O2 (yalnız RAPORLANIR, karar vermez): se(p_zarf)

180'in partisyonu (rapor §4.1'den, HÜKÜMSÜZ = {Q_X, ρ_X, θ, g_X,
g_cal}; okunur = {M, Q_E, ρ_E, g_E}) O1 ile İKİ KONVANSİYONDA DA
birebir yeniden üretilmelidir; üretmezse bu betik ÖLÜR (SystemExit) ve
181 ön-kayıtsız koşmaz.

ÖN-KAYITLI ÖNGÖRÜLER (n=4, koşullu — koşul açıkça yazılıyor):
  H1: tohum popülasyonunun sd'si n=2 kestirimindeki değerinde kalırsa
      se_VS(n=4) = sd_VS(n=2)/2 = se_VS(n=2)/√2,
      se_ZARF(n=4) = √(se_VS(n=2)²/2 + se_VF²)      [se_VF DEĞİŞMEZ]
  H2: bir satır, ZARF merkezi yerinde kalırsa, ancak ve ancak
      sd_VS(n=4) ≤ sd*_VS := 2·√((|ZARF|/2)² − se_VF²) olduğunda
      OKUNUR olur. Sağdaki karekökün içi negatifse (|ZARF| < 2·se_VF)
      satır **VS tohumu eklemekle ULAŞILAMAZ**dır — n ne olursa olsun.
  H3: manşet (M, log) için öngörülen yeni hata çubuğu ve iki
      konvansiyonda öngörülen p_zarf ± se.

Bu öngörüler KEHANET DEĞİL, ÖLÇÜT'tür: gerçek n=4 sayıları geldiğinde
hüküm yalnız O1'e göre verilir; öngörü tutmazsa tutmadığı yazılır.

KAPSAM (donduruluyor): kırpma hakemi (K2 / VF ailesi / 169_k2b) 181'in
DIŞINDADIR. 180'in H-180a HÜKÜMSÜZ'ü aynen durur; VF ailesi ve ρ₃
sayılarına dokunulmaz.

Kullanım: 181a_onkayit_n4.py
"""
import hashlib
import json
import math
import time
from pathlib import Path

import numpy as np

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S180, S181 = SCR / "180", SCR / "181"

HAM = ("M", "Q_E", "rho_E", "Q_X", "rho_X")
ORAN = ("g_E", "g_X", "theta", "g_cal")
AD = {"M": "M (=KALİB_u2)", "Q_E": "Q_E", "rho_E": "ρ_E", "Q_X": "Q_X",
      "rho_X": "ρ_X", "g_E": "g_E", "g_X": "g_X", "theta": "θ",
      "g_cal": "g_cal"}

# --- 180 raporu §4.1'in SÖZEL partisyonu (kalibrasyon hedefi) ---------
P180_OKUNUR = {"M", "Q_E", "rho_E", "g_E"}
P180_HUKUMSUZ = {"Q_X", "rho_X", "theta", "g_X", "g_cal"}
K_SIGMA = 2.0            # O1 eşiği — aşağıda kalibre ediliyor


def main():
    t0 = time.time()
    S181.mkdir(parents=True, exist_ok=True)
    me = Path(__file__).resolve()
    sha_me = hashlib.sha256(me.read_bytes()).hexdigest()
    OK180 = json.load(open(S180 / "ONKAYIT_K0.json"))
    K3 = json.load(open(S180 / "K3_DEFTER.json"))
    A = np.load(S180 / "A_sonzarf.npy")
    sha_A = hashlib.sha256(np.ascontiguousarray(A).tobytes()).hexdigest()

    print("=" * 78, flush=True)
    print("181a — K0 ÖN-KAYIT: VS ailesinin n=2 → n=4 UZANTISI", flush=True)
    print(f"    zaman = {time.ctime()}", flush=True)
    print(f"    betik sha256 = {sha_me}", flush=True)
    print(f"    180a ön-kayıt sha = {OK180['sha256']}", flush=True)
    print(f"    180 K3 defteri     = {K3['zaman']}", flush=True)
    print(f"    A_sonzarf.npy sha  = {sha_A}", flush=True)
    print(f"      ön-kayıttaki sha = {OK180['zarf']['sha_A_VS']}", flush=True)
    zarf_ok = (sha_A == OK180["zarf"]["sha_A_VS"])
    print(f"      ZARF BİT-BİT AYNI: {'✓' if zarf_ok else '✗'}", flush=True)
    if not zarf_ok:
        raise SystemExit("A_sonzarf.npy 180a'nınki değil — 181 koşmaz.")

    # ---------------------------------------------------- O1 kalibrasyonu
    print("\n" + "=" * 78)
    print("  K0.1 — OKUNUR ÖLÇÜTÜNÜN KALİBRASYONU (veri gelmeden)")
    print(f"  O1: |ZARF| ≥ {K_SIGMA:g}·se(ZARF).  180 §4.1 partisyonunu iki "
          "konvansiyonda da birebir üretmeli.")
    print("=" * 78)
    print("  nicelik    konv.     |ZARF|/se   O1      180 §4.1   uyum")
    kalib_ok = True
    for k in HAM + ORAN:
        for konv in ("dogrusal", "log"):
            r = K3[konv][k]
            t = abs(r["ZARF"]) / r["se_ZARF"]
            o1 = "OKUNUR" if t >= K_SIGMA else "HÜKÜMSÜZ"
            p180 = "OKUNUR" if k in P180_OKUNUR else "HÜKÜMSÜZ"
            u = (o1 == p180)
            kalib_ok &= u
            print(f"  {AD[k]:10s} {konv:9s} {t:9.3f}   {o1:9s} {p180:9s} "
                  f"{'✓' if u else '✗'}", flush=True)
    print(f"\n  KALİBRASYON: {'✓ birebir' if kalib_ok else '✗ TUTMADI'}",
          flush=True)
    if not kalib_ok:
        raise SystemExit("O1, 180'in partisyonunu üretmiyor — ön-kayıt YOK.")

    # ------------------------------------------- n=4 öngörüleri (H1..H3)
    print("\n" + "=" * 78)
    print("  K0.2 — n=4 ÖNGÖRÜLERİ (H1: sd_VS n=2 kestiriminde kalırsa;")
    print("         H2: sd*_VS = 2√((|ZARF|/2)² − se_VF²), ZARF merkezi sabit)")
    print("=" * 78)
    print("  nicelik    konv.    se_ZARF(2) se_ZARF(4)ö  se_p(2) se_p(4)ö  "
          "|Z|/se(4)ö  öngörü      sd_VS(2)   sd*_VS     durum")
    ONG = {}
    for k in HAM + ORAN:
        for konv in ("dogrusal", "log"):
            r = K3[konv][k]
            seVS2, seVF, Z, D = r["VS_se"], r["VF_se"], r["ZARF"], r["delta"]
            sdVS2 = r["VS_sd"]
            seZ2 = r["se_ZARF"]
            seVS4 = sdVS2 / 2.0                       # = seVS2/√2
            seZ4 = math.sqrt(seVS4 ** 2 + seVF ** 2)
            t4 = abs(Z) / seZ4
            ong = "OKUNUR" if t4 >= K_SIGMA else "HÜKÜMSÜZ"
            ic = (abs(Z) / K_SIGMA) ** 2 - seVF ** 2
            if ic <= 0.0:
                sdstar, durum = float("nan"), "ULAŞILAMAZ (|ZARF|<2·se_VF)"
            else:
                sdstar = 2.0 * math.sqrt(ic)
                durum = ("şimdiden geçer" if sdVS2 <= sdstar
                         else "sd %.2f× küçülmeli" % (sdVS2 / sdstar))
            simdi = "OKUNUR" if abs(Z) / seZ2 >= K_SIGMA else "HÜKÜMSÜZ"
            print(f"  {AD[k]:10s} {konv:9s} {seZ2:9.6f} {seZ4:11.6f} "
                  f"{r['se_p']:8.4f} {seZ4/abs(D):8.4f} {t4:10.3f}  "
                  f"{ong:9s} {sdVS2:10.6f} {sdstar:10.6f}  {durum}",
                  flush=True)
            ONG[f"{konv}/{k}"] = dict(
                nicelik=k, konvansiyon=konv, delta=D, ZARF=Z,
                se_ZARF_n2=seZ2, se_p_n2=r["se_p"], t_n2=abs(Z) / seZ2,
                simdiki_hukum=simdi, sd_VS_n2=sdVS2, se_VF=seVF,
                se_VS_n2=seVS2, se_VS_n4_ongoru=seVS4,
                se_ZARF_n4_ongoru=seZ4, se_p_n4_ongoru=seZ4 / abs(D),
                t_n4_ongoru=t4, ongorulen_hukum=ong,
                sd_yildiz_VS=(None if ic <= 0 else sdstar),
                ulasilabilir=bool(ic > 0), durum=durum)

    # ------------------------------------------------------------- manşet
    mlog = ONG["log/M"]
    mlin = ONG["dogrusal/M"]
    qlog = ONG["log/Q_E"]
    print("\n" + "=" * 78)
    print("  K0.3 — MANŞETİN ÖNGÖRÜLEN HÂLİ (H3)")
    print("=" * 78)
    print(f"  180 manşeti  : M log zarf %{100*K3['log']['M']['p_zarf']:.1f} "
          f"± %{100*K3['log']['M']['se_p']:.1f}")
    print(f"  n=4 öngörüsü : ZARF merkezi değişmezse ± "
          f"%{100*mlog['se_p_n4_ongoru']:.1f}  (H1 altında)")
    print(f"  Q_E log      : ± %{100*K3['log']['Q_E']['se_p']:.1f} → ± "
          f"%{100*qlog['se_p_n4_ongoru']:.1f} (öngörü)")
    print(f"  M doğrusal   : ± %{100*K3['dogrusal']['M']['se_p']:.1f} → ± "
          f"%{100*mlin['se_p_n4_ongoru']:.1f} (öngörü)")

    yeni = [v for v in ONG.values()
            if v["simdiki_hukum"] == "HÜKÜMSÜZ"
            and v["ongorulen_hukum"] == "OKUNUR"]
    ulasilamaz = [v for v in ONG.values()
                  if v["simdiki_hukum"] == "HÜKÜMSÜZ" and not v["ulasilabilir"]]
    print("\n  ÖN-KAYITLI BEKLENTİ — n=4'te OKUNUR olması beklenen yeni "
          "satırlar (H1):")
    print("   ", ", ".join(f"{v['nicelik']}/{v['konvansiyon']}"
                           for v in yeni) or "(yok)")
    print("  ÖN-KAYITLI BEKLENTİ — VS tohumu eklemekle ULAŞILAMAZ satırlar "
          "(H2, |ZARF| < 2·se_VF):")
    print("   ", ", ".join(f"{v['nicelik']}/{v['konvansiyon']}"
                           for v in ulasilamaz) or "(yok)")

    # -------------------------------------------------------------- kapılar
    OUT = dict(
        gorev=181, ad="VS3/VS4 TOHUMLARI — K0 ÖN-KAYIT (n=4 uzantısı)",
        zaman=time.ctime(), zaman_unix=time.time(),
        betik="181a_onkayit_n4.py", sha256=sha_me,
        miras=dict(onkayit_180a_sha=OK180["sha256"],
                   K3_defter_zaman=K3["zaman"],
                   sha_A_sonzarf=sha_A,
                   sha_A_onkayit=OK180["zarf"]["sha_A_VS"],
                   formul=OK180["formul"]),
        tohumlar=dict(VS_eski=[1, 2], VS_yeni=[3, 4], VS_n4=[1, 2, 3, 4],
                      VF=[1, 2, 3, 4]),
        insa_kapilari=dict(
            V1="sha256(A_VS3)=sha256(A_VS4)=sha_A_sonzarf (BİT-BİT)",
            V2="maks|A_VS − a_Hk·r(τ)| == 0.0",
            V3="φ≡0 sağlaması maks|ΔS| = maks|ΔS'| = 0.0",
            V4="ilk-kök hücre 300000/300000 (benzersiz)",
            V5="maks|F| ≤ 1e−8, aşan tekne 0",
            V6="sıralılık TAM (min Δz > 0)",
            not_="Biri tutmazsa o gaz ÖLÇÜLMEZ; kurtarma yok."),
        okunurluk=dict(
            O1=f"|ZARF_N| ≥ {K_SIGMA:g}·se(ZARF_N)  [KARAR VEREN ÖLÇÜT]",
            O2="se(p_zarf) — yalnız raporlanır, karar vermez",
            kalibrasyon="180 §4.1 partisyonu iki konvansiyonda da birebir "
                        "üretildi (bu betik, veri gelmeden)",
            p180_okunur=sorted(P180_OKUNUR),
            p180_hukumsuz=sorted(P180_HUKUMSUZ), k_sigma=K_SIGMA),
        H1="se_VS(n=4)=sd_VS(n=2)/2; se_ZARF(n=4)=√(se_VS(n=2)²/2+se_VF²); "
           "se_VF DEĞİŞMEZ (VF zaten n=4)",
        H2="sd*_VS = 2√((|ZARF|/2)²−se_VF²); iç negatifse satır VS tohumu "
           "eklemekle ULAŞILAMAZ",
        H3=dict(manset_180="M log zarf payı %81.5 ± %11.3",
                manset_n4_ongoru_se_p=mlog["se_p_n4_ongoru"],
                QE_log_n4_ongoru_se_p=qlog["se_p_n4_ongoru"],
                M_dogrusal_n4_ongoru_se_p=mlin["se_p_n4_ongoru"]),
        ongoruler=ONG,
        beklenen_yeni_okunur=[f"{v['nicelik']}/{v['konvansiyon']}"
                              for v in yeni],
        ulasilamaz=[f"{v['nicelik']}/{v['konvansiyon']}" for v in ulasilamaz],
        kapsam_disi=dict(
            kirpma_hakemi="169_k2b / ρ₃ / C / z_G / VF ailesi — 181'in "
                          "DIŞINDA. 180'in H-180a HÜKÜMSÜZ hükmü aynen "
                          "durur; hiçbir sayısı yeniden hesaplanmaz.",
            W_pos="179 ŞART ③ — 181'de de ele alınmıyor, borç duruyor."),
        isaret_kurali=OK180["formul"]["isaret_kurali"],
        sure_s=time.time() - t0)
    p = S181 / "ONKAYIT_181_K0.json"
    p.write_text(json.dumps(OUT, indent=1, ensure_ascii=False, default=float))
    print(f"\n-> {p}  ({time.time()-t0:.2f} s)", flush=True)
    print("ÖN-KAYIT MÜHÜRLENDİ — VS3/VS4 henüz inşa EDİLMEDİ.", flush=True)


if __name__ == "__main__":
    main()
