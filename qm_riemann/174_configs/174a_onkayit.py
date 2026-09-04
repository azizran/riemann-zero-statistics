# -*- coding: utf-8 -*-
"""
174a — K2'NİN ÖN-KAYDI: KURAL, K1 KOŞULMADAN ÖNCE DONDURULUYOR
================================================================
Görev 174 (GERÇEĞİN İMZASI). Bu betik HİÇBİR ÖLÇÜM YAPMAZ; yalnız
K2'nin eşleme KURALINI, hedeflerini ve ölüm ölçütlerini zaman damgalı
olarak dondurur ve kendi sha256'sını kaydeder. 174b (K1) bundan SONRA
koşar.

══════════════════════════════════════════════════════════════════════
DÜRÜSTLÜK BEYANI — neyi biliyorum, neyi bilmiyorum (ön-kayıt anı)
══════════════════════════════════════════════════════════════════════
BİLİYORUM (defterden, göreve verilen ve okunan):
  • 172e/G3'ün ölçtüğü hedefler: λ_eş(E kanalı) = 0.7768, λ_eş(X) =
    0.9089, λ_eş(marjinal) = 0.9354, λ_eş(g_E) = 0.7581, λ_eş(θ) =
    0.6916; σ-çapası λ_çapa = 0.9363; ΔM defteri g_E +%2.05,
    (g_X)² +%0.32, θ +%2.53, toplam +%4.97.
  • 167/C_*.json'un `artik` blokları (Var(η), Var(E_mod), g_E, ...) ve
    172/G1.json (K, S, μ̂², ρ, Q) — yani 172'nin BÜTÜN ÇARPAN DEFTERİ.
  • 162'nin gerçek gazda ölçtüğü iki sayı: R_η = 1.288, R_Ĉ = 1.026.
BİLMİYORUM (bu ön-kayıt anında ölçülmemiş):
  • R_η ve R_Ĉ'nin **Hkeskin (sadakatli ikiz)** ve **L085 (kontrol)**
    değerleri — 162 bunları hiç ölçmedi (162 §"sıradaki adım" 2 bunu
    açıkça borç yazıyor). 162'nin "keskin"i 152'nin gazıdır, Hkeskin
    DEĞİLDİR.
  • R'nin bant-bant dağılımı hiçbir gazda.
  • X̃ (=ΔĈ) kanalının girişim oranı hiçbir gazda (162 yalnız Ĉ'yi
    ölçtü).
  • Çizgi alanlarının ÜÇÜNCÜ momentleri hiçbir gazda (162 yalnız
    gerçek gazın η çarpıklığını, −0.151, kaydetti).

══════════════════════════════════════════════════════════════════════
KÖPRÜNÜN CEBRİ (ön-kayıtta TÜRETİLİYOR — ölçümden önce)
══════════════════════════════════════════════════════════════════════
165_cekirdek.Model165.alanlar:  E = Σ_q Re[hp_q e^{iω_q s}],
hp_q = 2⟨e1 e^{−iω_q s}⟩ ,  e1 = η_{n+1} − ort.  O halde ÖZDEŞ olarak

    K := ⟨E·e1⟩ = Σ_q Re[ conj(hp_q)·⟨e1 e^{iω_q s}⟩ ] = Σ_q |hp_q|²/2

yani **172'nin π := K/V_O'su, 162'nin girişim oranının ta kendisidir**
(V_O = Var(η)). Ve 172b'nin özdeşliği g_E = 1 − Q/ρ = π/ρ olduğundan

    ***  g_E = R_η / ρ_E  ***          (R_η := Σ|c_η|²/2 ÷ Var(η))

Bu bir YASA DEĞİL, CEBİRDİR. Ampirik içerik iki yerdedir:
 (i) 162'nin R'si BAŞKA bir çizgi kümesinde (τ ≤ 0.86, 3425 çizgi),
     BAŞKA bir sitede (m_n, η_n) ölçülür; 172'nin π'si τ ≤ 0.95'te
     8981 çizgide, s_n = m_{n+1} sitesinde. İkisinin AYNI sayıyı
     vermesi zorunlu değildir ve sınanmamıştır.
 (ii) R'nin λ merdiveni boyunca nasıl aktığı hiç ölçülmedi.

══════════════════════════════════════════════════════════════════════
K2 KURALLARI — DONDURULDU
══════════════════════════════════════════════════════════════════════
Girdi: 174b'nin (K1) ölçtüğü R_η(gaz), R_Ĉ(gaz), R_X(gaz); 162
makinesiyle, AYNI pencere (son 300k), AYNI taban (0.40), AYNI cap
(4000), AYNI çizgi evreni (τ ≤ 0.86).

 K2-A  λ_eş ÖNGÖRÜSÜ.
       λ_eş^ön(E) := R_η(λ) eğrisinin R_η(son)'u verdiği λ.
       Birincil eğri: 172e'nin tanımıyla AYNI — 7-nokta λ merdiveni
       (0.50,0.60,0.70,0.85,1.00,1.15,1.30) üzerinde tekdüze kolda
       doğrusal ters çevirme; KOŞULURSA. Yedek/asgari: KALEM'in
       istediği iki nokta (L085, Hkeskin) üzerinde LOG-DOĞRUSAL
       ters çevirme.
       HEDEF: 0.7768 (E kanalı ortalaması); ikinci hedef 0.7581 (g_E).
       ÖLÇÜT: |λ_eş^ön − 0.7768| ≤ 0.05 ⇒ TAM İSABET;
              ≤ 0.10 ⇒ KISMİ; > 0.10 ⇒ K2-A ÖLDÜ.

 K2-B  g_E PAYI ÖNGÖRÜSÜ (parametresiz).
       Δlog g_E^ön := log R_η(son) − log R_η(λ_çapa=0.9363),
       ρ_E DÜZELTMESİ YOK (ön-kayıtlı varsayım: ρ_E λ boyunca düz).
       HEDEF: +2.05% (log 0.02029).
       ÖLÇÜT: |Δ^ön − 2.05| ≤ 0.5 puan ⇒ TAM; ≤ 1.0 ⇒ KISMİ;
              > 1.0 puan ⇒ K2-B ÖLDÜ.
       (Ayrıca ρ-düzeltmeli sürüm de yazılır; ama HÜKÜM düzeltmesiz
        sürüme göre verilir — parametresizlik ölçütü budur.)

 K2-C  θ PAYI ÖNGÖRÜSÜ.
       162: Ĉ kanalında girişim YOK (1.026). Ön-kayıt: R_Ĉ üç gazda
       ±%2 içinde AYNI çıkacak ⇒ bu köprü θ'ya HİÇBİR ŞEY söylemez ⇒
       **Δθ^ön = 0.00%**.  HEDEF: +2.53%. Yani K2-C, ÖNCEDEN, θ'nın
       bu köprüyle kapanMAyacağını yazıyor.
       (Eğer R_Ĉ ya da R_X üç gazda ±%2'yi aşan bir λ-eğilimi
        gösterirse ön-kayıt bu maddede yanılmış olur ve yazılır.)

 K2-D  ΔM KAPANIŞI ve H-G3 MÜHRÜ.
       kapanış := (Δlog g_E^ön) / Δlog M(ölçülen = 0.04851 log/+4.97%).
       ÖN-KAYITLI BEKLENTİ: kapanış ≈ %41 (2.05/4.97) — yani
       **H-G3 (taç) BU KÖPRÜDEN MÜHÜRLENMEZ** (eşik %70).
       Ayakta kalması beklenen H-G1'dir: fazlanın g_E payı girişim
       eksikliğinin parametresiz faturasıdır.
       H-G3'ün mühürlenmesi ancak K2-E aşağıdaki gibi çıkarsa mümkün.

 K2-E  ÜÇÜNCÜ MOMENT KÖPRÜSÜ (θ için tek aday, tek yönlü ön-kayıt).
       θ tanımı gereği bir ÜÇÜNCÜ moment nesnesidir
       (θ = KALİB_u2/g_cal, KALİB ∝ ⟨e1·x1²·e^{−iWs}⟩); girişim oranı
       İKİNCİ moment nesnesidir. Rastgele fazlı çizgi alanının bütün
       üçüncü momentleri ÖZDEŞ SIFIRDIR (Rice/Gauss) ⇒ ölçülen üçüncü
       moment saf faz-kilidi ölçüsüdür.
       Ölçülecek: m3 := ⟨e1·x1²⟩/(σ_e1·σ_x1²) ve çizgi-alanı eşleniği.
       ÖN-KAYITLI TEK YÖNLÜ ÖNGÖRÜ:
           λ_eş(m3) < λ_eş(R_η)          [çünkü λ_eş(θ)=0.6916 <
                                          λ_eş(E)=0.7768]
       Ters çıkarsa K2-E ÖLÜR. Nicel hedef: λ_eş(m3) ≈ 0.69 (±0.10).

 K2-F  KÖPRÜ KİMLİĞİ DENETİMİ (162 ↔ 172).
       R_η(162 makinesi, τ≤0.86, site m_n) ile π_E(172 defteri,
       τ≤0.95, site s_n) her gazda karşılaştırılır.
       ÖN-KAYIT: |R_η/π_E − 1| ≤ %5 (üç gazda da).  Aşarsa köprünün
       "aynı nesne" okuması zayıflar ve yazılır.

 K2-G  X KANALI SAĞLAMASI (K4).
       Aynı eşleme X kanalına uygulanır: λ_eş^ön(X) := R_X(λ)
       eğrisinin R_X(son)'u verdiği λ.  HEDEF: 0.9089.
       ÖLÇÜT: |λ_eş^ön(X) − 0.9089| ≤ 0.05 ⇒ TAM; ≤ 0.10 ⇒ KISMİ.

══════════════════════════════════════════════════════════════════════
ÖLÜM KURALI: hiçbir kural koşudan sonra gevşetilmez; ıskalayan madde
"ÖLDÜ" yazılır ve raporda aynen kalır.
══════════════════════════════════════════════════════════════════════
"""
import datetime
import hashlib
import json
import pathlib

SCR = pathlib.Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
                   "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/174")
SELF = pathlib.Path(__file__).resolve()

HEDEF = dict(
    lam_es_E=0.7768, lam_es_gE=0.7581, lam_es_X=0.9089,
    lam_es_marj=0.9354, lam_es_theta=0.6916, lam_capa=0.9363,
    d_gE_pct=2.05, d_gX2_pct=0.32, d_theta_pct=2.53, dM_pct=4.97,
    R_eta_162_gercek=1.288, R_Chat_162_gercek=1.026,
)

KURAL = {
 "K2-A": dict(ad="λ_eş öngörüsü (E kanalı)",
              tanim="λ_eş^ön(E) = R_η(λ) eğrisinin R_η(son)'u verdiği λ; "
                    "birincil 7-nokta merdiveni (172e ile aynı ters "
                    "çevirme), yedek 2-nokta (L085,Hkeskin) log-doğrusal",
              hedef=0.7768, tam=0.05, kismi=0.10),
 "K2-B": dict(ad="g_E payı (parametresiz)",
              tanim="Δlog g_E^ön = log R_η(son) − log R_η(0.9363); "
                    "ρ_E düzeltmesi YOK",
              hedef=2.05, tam=0.5, kismi=1.0, birim="%"),
 "K2-C": dict(ad="θ payı",
              tanim="R_Ĉ üç gazda ±%2 içinde aynı ⇒ Δθ^ön = 0",
              ongoru=0.0, hedef=2.53, birim="%"),
 "K2-D": dict(ad="ΔM kapanışı / H-G3",
              tanim="kapanış = Δlog g_E^ön / 0.04851",
              beklenen_pct=41.0, muhur_esigi_pct=70.0),
 "K2-E": dict(ad="üçüncü moment köprüsü (θ adayı)",
              tanim="m3 = ⟨e1 x1²⟩/(σ_e1 σ_x1²); tek yönlü: "
                    "λ_eş(m3) < λ_eş(R_η)",
              hedef=0.6916, tolerans=0.10),
 "K2-F": dict(ad="köprü kimliği 162↔172",
              tanim="|R_η/π_E − 1| ≤ %5, üç gazda", esik_pct=5.0),
 "K2-G": dict(ad="X kanalı sağlaması (K4)",
              tanim="λ_eş^ön(X) = R_X(λ) eğrisinin R_X(son)'u verdiği λ",
              hedef=0.9089, tam=0.05, kismi=0.10),
}

if __name__ == "__main__":
    SCR.mkdir(parents=True, exist_ok=True)
    src = SELF.read_bytes()
    h = hashlib.sha256(src).hexdigest()
    ts = datetime.datetime.now().astimezone().isoformat(timespec="seconds")
    rec = dict(gorev=174, kapi="K2", zaman=ts, betik=str(SELF),
               sha256=h, hedef=HEDEF, kural=KURAL,
               bilinmeyen=["R_η(Hkeskin)", "R_η(L085)", "R_Ĉ(Hkeskin)",
                           "R_Ĉ(L085)", "R_X(hepsi)", "bant-bant defteri",
                           "üçüncü momentler (hepsi)"],
               bilinen=["172/G1.json", "172/G3.json", "167/C_*.json",
                        "162 raporu: R_η=1.288, R_Ĉ=1.026 (yalnız gerçek)"])
    p = SCR / "ONKAYIT_K2.json"
    if p.exists():
        eski = json.loads(p.read_text())
        print("!! ONKAYIT_K2.json ZATEN VAR — üzerine YAZILMAZ.")
        print("   ilk kayıt:", eski["zaman"], eski["sha256"][:16])
        raise SystemExit(0)
    p.write_text(json.dumps(rec, indent=1, ensure_ascii=False))
    print("=" * 70)
    print("174a — K2 ÖN-KAYDI DONDURULDU")
    print("=" * 70)
    print("zaman   :", ts)
    print("sha256  :", h)
    print("dosya   :", p)
    for k, v in KURAL.items():
        print(f"  {k}  {v['ad']}")
        print(f"        {v['tanim']}")
    print("\nBİLİNMEYENLER (ön-kayıt anında ölçülmemiş):")
    for b in rec["bilinmeyen"]:
        print("   -", b)
