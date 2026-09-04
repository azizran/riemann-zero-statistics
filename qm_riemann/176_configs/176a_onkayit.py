# -*- coding: utf-8 -*-
"""
176a — K2 ÖN-KAYIT: KARIŞTIRMA SINAVININ ÇÖKÜŞ ÖNGÖRÜLERİ
==========================================================
Bu betik HİÇBİR 176 ölçümü koşmadan ÖNCE koşar; kendi sha256'sını ve
zaman damgasını `176/ONKAYIT_K2.json`'a yazar ve dosya bir daha
YAZILMAZ (varlık kontrolü). Okuduğu her şey ÖNBELLEKTEN gelir
(172/G1.json, 174/K1_*.json, 174/K3_*.json, 175/K2.json) — yani
174/175'in yayımlanmış defteridir.

══════════════════════════════════════════════════════════════════
DÜRÜSTLÜK BEYANI (ön-kayıtta aynen durur)
══════════════════════════════════════════════════════════════════
1. KÖRLÜK İDDİASI YOKTUR. 174/175'in bütün sayıları (ΔM = +%5.92,
   kesim %39 / kilit %61, θ'nın %51.5'i, τ>0.70 faz uyumu %78,
   genlik %77, Kov(çizgi,artık) = −0.0435) bu tayfa tarafından
   OKUNMUŞTUR. Ön-kayıt edilen şey KURALDIR: hangi ölçümden hangi
   öngörünün hangi formülle ve hangi ÖLÜM EŞİĞİYLE çıkacağı.
2. ÖN-KAYIT ANINDA ÖLÇÜLMEMİŞ OLANLAR (hepsi 176'da ilk kez ölçülür):
   faz-karıştırılmış vekil gazın HİÇBİR niceliği. Böyle bir gaz
   (merdiven fazları rastgele + sadakatli zincirle YENİDEN İNŞA)
   bu programda hiç kurulmamıştır. 162'nin vekilleri KONUM-UZAYINDA,
   kurulmuş gazın ALANLARININ fazını karıştırır; 176'nınki
   MERDİVENİN fazını karıştırıp gazı sıfırdan çözer. İkisi farklı
   nesnelerdir.
3. VEKİLİN ZARFI GERÇEĞİN ZARFIDIR. Gerçek ζ gazının merdiveni
   nominal olarak S(t) = −Σ_q a_q sin(ω_q t), a_q = 1/(πk√q),
   τ_q ≤ 1.00 — yani `Hkeskin`'in ta kendisi. Vekil bu a_q dizisini
   BİREBİR (bit-bit) korur, yalnız her çizgiye rastgele bir φ_q
   takar. Kesim koordinatı böylece Hkeskin'e ÇİVİLENİR; ölçülen her
   fark yalnız KİLİDİN farkıdır.
4. ÖN-KAYIT SONRASI HİÇBİR EŞİK DEĞİŞTİRİLMEZ. Ölümler kurtarmasız
   yazılır.

══════════════════════════════════════════════════════════════════
TÜRETİMLER (ölçümden ÖNCE, kalemle)
══════════════════════════════════════════════════════════════════
T-1  VEKİL MERDİVENİ.  S_v(t) = −Σ_q A_q sin(ω_q t + φ_q),
     φ_q ~ U(0,2π) bağımsız.  |A_q| birebir Hkeskin'inki.
     Güç tayfı ÖZDEŞ (rms S, rms S', Σa_qω_q hepsi aynı);
     değişen tek şey çizgiler-arası BAĞIL FAZ.

T-2  R_η'NİN KİLİT İÇERİĞİ.  174-K1d'nin özdeşliği:
        P ≡ Kov(η_çizgi, η) ,  R_η = P/Var(η) ,
        Var(η) = Var(η_çiz) + 2Kov(η_çiz,η_art) + Var(η_art).
     R_η > 1 ⟺ Kov(η_çizgi, η_artık) < 0 (Hkeskin: −0.05574).
     Artık, merdivenin τ>0.86 kuyruğu + doğrusal-olmayan terimlerdir;
     kilitli gazda bunlar çizgilerle AYNI faz kökünden doğar.
     Fazlar rastgele olunca bu kovaryans YAPISAL olarak kaybolur
     (162-V5: rastgele fazlı bir alanda çıkarım yansızdır, 1.0045).
     Kov → 0 limitinde, Hkeskin'in KENDİ sayılarıyla:
        R_η(vekil) ≈ P/(P + Var_artık) = 0.085933/(0.085933+0.037731)
                   = 0.6949
     Bu bir NOKTA öngörüsü değil bir MERKEZDİR (vekil yeniden
     kurulduğu için Var_artık da değişir); bant aşağıda.

T-3  DC KAÇAĞI KİLİDE KÖRDÜR (ters yönlü öngörü).
     174d/175g'nin ölçtüğü yasa: κ(ω_Q) = −π A_Q τ_Q cos(πτ_Q).
     Faz taşıyıcısıyla birlikte yazılırsa, d(t) = ḡΣA_q sin(ω_q t+φ_q)
     yer değiştirmesinden
        κ(ω_Q) ≈ −π τ_Q A_Q cos(πτ_Q) · e^{−iφ_Q}   (arg κ = −φ_Q + c₂)
     ve ölçülen çizgi katsayısı hp_Q ∝ b_Q e^{+iφ_Q} (arg hp = +φ_Q + c₁).
     Dolayısıyla
        Δφ_Q := arg hp_Q + arg κ(ω_Q) = c₁ + c₂   —  φ_Q'DAN BAĞIMSIZ.
     ⇒ DC kaçağı ξ_q = |hp||κ|cosΔφ bir ÖZ-REZONANSTIR ve faz
     karıştırmasına KÖRDÜR. ÖNGÖRÜ: ort(E)(vekil) POZİTİF kalır ve
     Hkeskin'in +0.062139'undan ±%30'dan fazla sapmaz; ⟨cosΔφ⟩
     0.50–0.80 bantlarında ≥ 0.90 kalır. Bu türetim yanlışsa ölçüm
     onu öldürür ve öyle yazılır.

T-4  θ VE Ç1 KANALI (H-F1b'nin ters yönlü türetimi).
     θ = KALİB_u2 / (g_E g_X²) bir ORANDIR: ölçülen üçüncü moment /
     model alanının üçüncü momenti. 165'in Ç1 kanalı (`kanal1`,
     "çift-iptal, ν ≡ 0 aritmetiksiz") q₂ = q₃ eşleşmesiyle çalışır;
     fazı arg hp_Q'dan başka hiçbir çizginin fazını taşımaz. Yani Ç1
     de KİLİDE KÖRDÜR. Üçüncü momentin baskın kanalı Ç1 ise θ
     ÇÖKMEZ. ÖN-KAYITLI BEKLENTİ: **(b) θ çökmez** ⇒ H-F1b ölür ve
     θ "üçüncü eksen" olarak açık kalır. (a) çıkarsa T-4 ölür.

T-5  ÜÇÜNCÜ MOMENTLERİN RICE SIFIRI.  Rastgele fazlı bir çizgi
     alanının bütün üçüncü momentleri özdeş sıfırdır (174 §1e).
     Vekilin ÇİZGİ alanları (η_çizgi, X_çizgi) tam olarak böyledir
     ⇒ m3_çizgi ve skew(η_çizgi) çökmelidir. Buna karşılık ds'nin
     çarpıklığı ÇÖKMEK ZORUNDA DEĞİLDİR (162-V4: konum-uzayı
     vekilinde skew(ds) 0.282 → 0.534 BÜYÜDÜ; öz-tutarlılığın
     doğrusal-olmayanlığı fazdan bağımsız çarpıklık üretir).

══════════════════════════════════════════════════════════════════
DONDURULMUŞ KURAL (ölçüm bunlara göre yargılanır)
══════════════════════════════════════════════════════════════════
F0  İNŞA KAPILARI: ilk-kök hücre 300000/300000; maks|F| ≤ 1e−8;
    sıralılık TAM; ZARF ÖZDEŞLİĞİ maks|A_q(vekil) − A_q(Hkeskin)| = 0
    (tam sıfır, makine hassasiyeti değil); ölçüm tarafında
    R_bant ≥ 0.98 (169 filtresi, lo ∈ [0.52,0.68]).
    Biri tutmazsa vekil ÖLÇÜLMEZ.

F1  R_η(vekil):  merkez 0.6949, BANT [0.60, 1.05].
    GEREK ŞART (kilit faturasının alt sınırı):
      Λ_R := R_η(Hkeskin) − R_η(vekil) ≥ (1−f(R_η))·[R_η(son)−R_η(Hk)]
           = 0.2555 × 0.0232450 = 0.0059391
    ÖLÜM: R_η(vekil) ≥ 1.15  ⇒ "girişim oranı = faz kilidi" okuması ÖLÜR.

F2  g_E(vekil): BANT [0.25, 0.57].
    GEREK ŞART: Λ_E := log g_E(Hkeskin) − log g_E(vekil) ≥ 0.0149737
      ( = 0.61 × Δlog g_E(son←Hkeskin) = 0.61 × 0.0245470 )
    ÖLÜM: g_E(vekil) ≥ g_E(Hkeskin) = 0.5837010  ⇒ H-F1'in E-kanalı ÖLÜR.

F3  θ(vekil): iki dal, ikisi de ön-kayıtlı.
    (a) H-F1b YAŞAR:  θ(vekil) ≤ 0.8656146  ( = θ(Hk)·e^{−0.0315180} )
    (b) H-F1b ÖLÜR :  θ(vekil) ≥ θ(Hkeskin) = 0.8933338
    ara bölge: KISMİ (ne yaşar ne ölür; sayı yazılır).
    ÖN-KAYITLI BEKLENTİ: **(b)** — T-4'ün Ç1 argümanı yüzünden.

F4  ΔM: Δlog M(vekil ← Hkeskin) BANT [−1.00, −0.06].
    GEREK ŞART: Λ_M := −Δlog M(vekil←Hk) ≥ 0.0575417 (bütün ΔM).
    ÖLÜM: Δlog M(vekil←Hk) ≥ 0.

F5  ÜÇÜNCÜ MOMENTLER (T-5): |m3_çizgi(vekil)| ≤ 0.040 (Hkeskin 0.0802)
    ve |skew(η_çizgi)(vekil)| ≤ 0.130 (Hkeskin 0.2501). skew(ds) için
    ÖN-KAYIT YOKTUR (T-5'in son cümlesi) — yalnız yazılır.

F6  DC KAÇAĞI (T-3, TERS YÖNLÜ): ort(E)(vekil) > 0 ve
    |ort(E)(vekil)/0.0621390 − 1| ≤ 0.30; ⟨cosΔφ⟩(vekil) 0.50–0.80
    bantlarında ≥ 0.90. Tutmazsa T-3 ÖLÜR (ve o zaman kilit DC
    kaçağını da taşıyor demektir — H-F1 için LEHTE bir sürpriz).

F7  TOHUM SAÇILIMI: en az iki tohum. Her F-maddesi için iki tohumun
    ORTALAMASI yargılanır; tohumlar-arası fark |y₁−y₂| HATA ÇUBUĞUDUR.
    Eğer |y₁−y₂| ölçülen değerin eşiğe uzaklığından BÜYÜKSE, o madde
    "HÜKÜMSÜZ (tohum gürültüsü)" yazılır — ne yaşar ne ölür.

F8  H-F2 (kilit payının ilk-ilke öngörüsü) — FORMÜL ŞİMDİ DONDURULUR,
    DEĞERİ HENÜZ HESAPLANMADI. 174d'nin özdeşliği
       ξ_b = Σ|hp| · ⟨|κ|⟩_w · ⟨cosΔφ⟩ ,
       amp_b := Σ_{q∈b}|hp_q| ,  res_b := Σ|hp||κ| / Σ|hp| ,
       faz_b := ξ_b / Σ|hp||κ|
    τ > 0.70 bölgesinde (açığın %83'ü), gerçek ↔ keskin ikiz:
       Δ_amp = log[amp(son)/amp(Hk)] , Δ_res = log[res/res] ,
       Δ_faz = log[faz/faz]
       **kilit_payı^{HF2} := Δ_faz / (Δ_amp + Δ_res + Δ_faz)**
    SINAV: |kilit^{HF2} − 0.61| ≤ 0.10 ⇒ H-F2 YAŞAR; değilse ÖLÜR
    (kurtarma yok, aday ölürse ölür).

F9  NİHAİ AYRIŞIM TABLOSU (K3) — formül donduruldu:
       ΔlogM(son←Hk) ≡ Δlog g_E + 2Δlog g_X + Δlog θ   (özdeş)
       her çarpan için  kesim payı := f_j·ΔlogX_j(erfc←Hk)
                        kilit payı := ΔlogX_j(son←Hk) − kesim payı
       f_j := ΔlogX_j(son←Hk)/ΔlogX_j(erfc←Hk)  (parametresiz)
    ve VEKİL SÜTUNU: ΔlogX_j(Hk←vekil) = gazın TOPLAM kilit içeriği.
    "Kilit faturası ödenebilir mi?" ölçütü:
       ödeme oranı  ω_j := kilit payı_j / ΔlogX_j(Hk←vekil)
    ω_j ∈ (0, 1] ise fatura ödenebilir; ω_j > 1 ya da ω_j < 0 ise
    kilit o çarpanı taşıyamaz.

MÜHÜR KURALI (dondurulmuş):
  * F1 ✓ ve F2 ✓ ve F3(a) ve F4 ✓  ⇒ **MÜHÜR-TAM**: "Gerçeğin
    üçüncü-moment fazlası, kesim şekli (%39, E) düşüldükten sonra
    asal fazlarının kilidinin faturasıdır."
  * F1 ✓ ve F2 ✓ ama F3(b)         ⇒ **MÜHÜR-E**: cümle YALNIZ E
    kanalı için mühürlenir; θ "üçüncü eksen" olarak keskin bir
    negatifle açık kalır.
  * F2 ✗                            ⇒ **ÖLÜM**: kilit faturası ödemiyor.
"""
import hashlib
import json
import math
import sys
import time
from pathlib import Path

QM = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")
S176 = SCR / "176"
OUT = S176 / "ONKAYIT_K2.json"

BU = Path(__file__).resolve()


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    S176.mkdir(parents=True, exist_ok=True)
    if OUT.exists():
        d = json.load(open(OUT))
        print("ÖN-KAYIT ZATEN VAR — ÜZERİNE YAZILMAZ.")
        print("  zaman : %s" % d["zaman"])
        print("  sha256: %s" % d["sha256"])
        return d

    G1 = json.load(open(SCR / "172" / "G1.json"))
    K1 = {g: json.load(open(SCR / "174" / f"K1_{g}.json"))
          for g in ("son", "Hkeskin", "HA4")}

    Hk, sn, ha = G1["Hkeskin"], G1["son"], G1["HA4"]
    dgE = math.log(sn["gE"] / Hk["gE"])
    dgX = math.log(sn["gX"] / Hk["gX"])
    dth = math.log(sn["th"] / Hk["th"])
    dM = math.log(sn["KAL"] / Hk["KAL"])
    egE = math.log(ha["gE"] / Hk["gE"])
    egX = math.log(ha["gX"] / Hk["gX"])
    eth = math.log(ha["th"] / Hk["th"])
    eM = math.log(ha["KAL"] / Hk["KAL"])

    Rs, Rh, Ra = (K1[g]["eta"]["R"] for g in ("son", "Hkeskin", "HA4"))
    fR = (Rs - Rh) / (Ra - Rh)
    P_h = K1["Hkeskin"]["eta"]["P"]
    Va_h = K1["Hkeskin"]["eta"]["Var_artik"]
    R_merkez = P_h / (P_h + Va_h)

    KILIT_PAY = 0.61          # 175'in ölçtüğü E-kanalı kilit payı
    lamE_gerek = KILIT_PAY * dgE
    lamR_gerek = (1.0 - fR) * (Rs - Rh)
    th_a = Hk["th"] * math.exp(-dth)

    kural = {
        "F0": {
            "ad": "inşa kapıları",
            "ilk_kok": "300000/300000", "maxF": 1e-8, "sirali": "TAM",
            "zarf_farki": 0.0, "R_bant": 0.98,
        },
        "F1": {
            "ad": "R_η(vekil)",
            "merkez": R_merkez, "bant": [0.60, 1.05],
            "gerek_Lambda_R": lamR_gerek,
            "olum": "R_η(vekil) >= 1.15",
            "capa": {"Hkeskin": Rh, "son": Rs, "HA4": Ra, "f(R_η)": fR},
        },
        "F2": {
            "ad": "g_E(vekil)",
            "bant": [0.25, 0.57],
            "gerek_Lambda_E": lamE_gerek,
            "olum": "g_E(vekil) >= %.7f" % Hk["gE"],
            "capa": {"Hkeskin": Hk["gE"], "son": sn["gE"], "HA4": ha["gE"],
                     "dlog_gE(son<-Hk)": dgE, "dlog_gE(erfc<-Hk)": egE},
        },
        "F3": {
            "ad": "θ(vekil)",
            "dal_a_yasar": th_a, "dal_b_olur": Hk["th"],
            "beklenti": "b (T-4: Ç1 kilide kördür)",
            "capa": {"Hkeskin": Hk["th"], "son": sn["th"], "HA4": ha["th"],
                     "dlog_th(son<-Hk)": dth, "dlog_th(erfc<-Hk)": eth},
        },
        "F4": {
            "ad": "Δlog M(vekil<-Hkeskin)",
            "bant": [-1.00, -0.06],
            "gerek_Lambda_M": dM,
            "olum": "Δlog M(vekil<-Hk) >= 0",
            "capa": {"Hkeskin": Hk["KAL"], "son": sn["KAL"], "HA4": ha["KAL"],
                     "dlogM(son<-Hk)": dM, "dlogM(erfc<-Hk)": eM,
                     "dlog_gX(son<-Hk)": dgX, "dlog_gX(erfc<-Hk)": egX},
        },
        "F5": {
            "ad": "üçüncü momentlerin Rice sıfırı",
            "m3_cizgi_esik": 0.040, "skew_eta_ciz_esik": 0.130,
            "capa": {"m3_cizgi(Hk)": K1["Hkeskin"]["m3"]["cizgi"],
                     "skew_eta_ciz(Hk)": K1["Hkeskin"]["m3"]["skew_eta_ciz"],
                     "skew_ds(Hk)": K1["Hkeskin"]["m3"]["skew_ds"]},
            "not": "skew(ds) için ön-kayıt YOKTUR (162-V4)",
        },
        "F6": {
            "ad": "DC kaçağı kilide kördür (T-3)",
            "ortE_capa": 0.0621390, "bagil_tolerans": 0.30,
            "cos_esik": 0.90, "cos_bantlari": [[0.50, 0.60], [0.60, 0.70],
                                               [0.70, 0.80]],
            "olum": "tutmazsa T-3 ölür",
        },
        "F7": {"ad": "tohum saçılımı", "min_tohum": 2,
               "kural": "|y1-y2| > |y_ort - esik| ise HÜKÜMSÜZ"},
        "F8": {
            "ad": "H-F2 kaba kalem",
            "formul": "kilit^HF2 = D_faz/(D_amp+D_res+D_faz), tau>0.70",
            "hedef": KILIT_PAY, "tolerans": 0.10,
            "not": "DEĞER ÖN-KAYIT ANINDA HESAPLANMADI",
        },
        "F9": {
            "ad": "nihai ayrışım tablosu",
            "formul": ("ΔlogM = Δlog g_E + 2Δlog g_X + Δlog θ; "
                       "kesim_j = f_j·Δlog_j(erfc<-Hk); "
                       "kilit_j = Δlog_j(son<-Hk) − kesim_j; "
                       "ω_j = kilit_j / Δlog_j(Hk<-vekil)"),
        },
        "MUHUR": {
            "TAM": "F1 ✓ ve F2 ✓ ve F3(a) ve F4 ✓",
            "E": "F1 ✓ ve F2 ✓ ama F3(b)",
            "OLUM": "F2 ✗",
            "cumle": ("Gerçeğin üçüncü-moment fazlası, kesim şekli "
                      "(%39, E) düşüldükten sonra asal fazlarının "
                      "kilidinin faturasıdır."),
        },
    }

    rec = dict(
        zaman=time.strftime("%Y-%m-%d %H:%M:%S %z"),
        betik=BU.name, sha256=sha(BU),
        vekil_recete=dict(
            merdiven="S_v(t) = -Sum_q A_q sin(w_q t + phi_q)",
            A_q="birebir Hkeskin: a_q = 1/(pi k sqrt(q)), tau_q <= 1.00, lam=1",
            phi="rng.uniform(0,2pi), np.random.default_rng(tohum)",
            cozucu="164_insa.coz_sadakatli (izgara braketi + korumali Newton, "
                   "SIRALI ILK-KOK), h=0.015, nz=300000, c=-0.5",
            tohumlar=[1, 2],
            adlar=["VF1", "VF2"],
        ),
        kural=kural,
        durustluk=[
            "korluk iddiasi yoktur; 174/175 defteri okunmustur",
            "on-kayit aninda olculmemis olan: vekil gazin HICBIR niceligi",
            "vekilin zarfi Hkeskin'in zarfidir (bit-bit)",
            "esikler bir daha degistirilmez; olumler kurtarmasiz yazilir",
        ],
    )
    OUT.write_text(json.dumps(rec, indent=1, ensure_ascii=False))

    print("=" * 74)
    print("176a — K2 ÖN-KAYIT YAZILDI")
    print("=" * 74)
    print("  zaman : %s" % rec["zaman"])
    print("  sha256: %s" % rec["sha256"])
    print("  -> %s" % OUT)
    print()
    print("  ÇAPALAR (önbellekten):")
    print("    R_η   : Hkeskin %.5f | son %.5f | HA4 %.5f   f(R_η) = %+.5f"
          % (Rh, Rs, Ra, fR))
    print("    g_E   : Hkeskin %.6f | son %.6f | HA4 %.6f" % (Hk["gE"],
          sn["gE"], ha["gE"]))
    print("    θ     : Hkeskin %.6f | son %.6f | HA4 %.6f" % (Hk["th"],
          sn["th"], ha["th"]))
    print("    M     : Hkeskin %.6f | son %.6f | HA4 %.6f" % (Hk["KAL"],
          sn["KAL"], ha["KAL"]))
    print("    Δlog(son←Hk): g_E %+.7f  2g_X %+.7f  θ %+.7f  M %+.7f"
          % (dgE, 2 * dgX, dth, dM))
    print("    Δlog(erfc←Hk): g_E %+.7f  2g_X %+.7f  θ %+.7f  M %+.7f"
          % (egE, 2 * egX, eth, eM))
    print()
    print("  DONDURULMUŞ ÖNGÖRÜ BANTLARI:")
    print("    F1  R_η(vekil)  merkez %.4f   bant [0.60, 1.05]   "
          "gerek Λ_R ≥ %.7f   ölüm ≥ 1.15" % (R_merkez, lamR_gerek))
    print("    F2  g_E(vekil)  bant [0.25, 0.57]   gerek Λ_E ≥ %.7f   "
          "ölüm ≥ %.6f" % (lamE_gerek, Hk["gE"]))
    print("    F3  θ(vekil)    (a) ≤ %.7f YAŞAR | (b) ≥ %.7f ÖLÜR   "
          "[beklenti: b]" % (th_a, Hk["th"]))
    print("    F4  ΔlogM(vekil←Hk) bant [−1.00, −0.06]  gerek Λ_M ≥ %.7f"
          % dM)
    print("    F6  ort(E)(vekil) ≈ +0.0621390 ±%%30, ⟨cosΔφ⟩ ≥ 0.90 "
          "(T-3: kilide kör)")
    print("    F8  kilit^HF2 = Δ_faz/(Δ_amp+Δ_res+Δ_faz) → |·−0.61| ≤ 0.10")
    print("=" * 74)
    return rec


if __name__ == "__main__":
    main()
