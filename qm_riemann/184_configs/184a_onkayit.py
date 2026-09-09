# -*- coding: utf-8 -*-
"""
184a — K0: DONMUŞ ÖN-KAYIT (ZARFIN ANATOMİSİ)
==============================================
ÖN-MÜHÜR. Bu betik HİÇBİR zarf niceliği ölçülmeden ÖNCE koşar; kendi
sha256'sını ve `date` damgasını 184/ONKAYIT_184.json'a yazar; dosya bir
daha DEĞİŞTİRİLMEZ. 184b… bu dosyayı okur, sha'sını doğrular, yalnız
burada donmuş tanım/form/eşikleri uygular.

────────────────────────────────────────────────────────────────────────
NE DONDURULUR
────────────────────────────────────────────────────────────────────────
(1) â_q ÖLÇÜM TANIMI (132c özdeşliği (★) + 162/174 makinesi):
    ds_n = Σ_Q 2 a_Q sin(ω_Q g_n/2) cos(ω_Q m_n)          (★, 132c)
      a_Q = 1/(π m √Q)  (Q=p^m),  ω_Q = log Q,  m_n = ½(z_n+z_{n+1}) HAM
      orta nokta,  g_n = z_{n+1}−z_n,  τ_Q = ω_Q/L.
    ÖLÇÜM:  c_q(ds) = 2⟨ds_n e^{−i ω_q m_n}⟩   (155/159/162 ile BİT-BİT
      aynı ham izdüşüm: `ds`, `mid` = eta_onbellek'ten; kalıntı yok).
    ETKİN GENLİK:  â_q := |c_q(ds)|.
    NOMİNAL (analitik, öz-tutarlı bastırmasız):
      a_q^eff := 2 a_q sin(π τ_q).
    ORAN:  w_q := â_q / a_q^eff.
    İkiz-kontrol: w_twin(τ) ≈ 1 (mütevazı τ-düşüşüyle; 132c Y2: %84–108).

(2) SAF ZARF (iki dilde, ikisi de raporlanır):
    r(τ)  := â_q(gerçek) / â_q(ikiz) = w_gerçek/w_ikiz   [ORAN — nominal
             TAM sadeleşir; 175'in "%77" dili; K2'nin BİRİNCİL nesnesi].
    D(τ)  := w_gerçek(τ) − w_ikiz(τ)                     [FARK — KALEM-K1
             sözü "w(gerçek)−w(ikiz) = saf zarf"].

(3) BANT IZGARASI (τ):  KENAR = [.45 .50 .55 .60 .65 .70 .75 .80 .86].
    Bant apsisi = o banttaki çizgilerin nominal-genlik-ağırlıklı ⟨τ⟩.
    Bant değeri (ORAN dili):
      r_b = (Σ_{q∈b} â_gerçek) / (Σ_{q∈b} â_ikiz)   [genlik-ağırlıklı].
    Bant değeri (w dili): w_b = (Σ â) / (Σ a_q^eff), gerçek ve ikiz ayrı.
    KUYRUK bandı ayrıca: τ>0.70 tek toplu bant (açığın %83'ü orada).

(4) JACKKNIFE (hata çubuğu): N=299999 nokta 8 BİTİŞİK bloğa bölünür.
    İzdüşümün blok-kısmi toplamları TEK GEÇİŞTE toplanır; leave-one-out
      c_q^(−k) = 2(Σ_hepsi − Σ_blok_k)/(N − n_k);
    bant nicelikleri her (−k) için yeniden kurulur; jackknife varyansı
      se² = (K−1)/K · Σ_k (x_(−k) − x̄)² ,  K=8.

(5) HİPOTEZLER — formlar + ÖLÜM EŞİKLERİ (r(τ)'ye fit; τ∈[.45,.86]):
    H-Z1 (erfc):     r(τ) = ½ erfc((τ − τ_c)/Δ)            [par: τ_c, Δ]
    H-Z2 (güç):      r(τ) = 1 − c·τ^α       (c,α>0; r∈[0,1] kırpılır)
    H-Z3 (Gauss-DW): r(τ) = exp(−(2π τ)² σ²/2)             [par: σ]
    En küçük kareler = jackknife-se ile ters-varyans ağırlıklı; dof =
      n_bant − n_par.
    SURVIVAL:  χ²/dof ≤ 2.0.       ÖLÜM:  χ²/dof > 2.0 (kurtarmasız).
    KUYRUK SINAVI (τ>0.70 bantları): standardize artık ortalaması
      |⟨R⟩|/se_R ≤ 2.0 → geçer;  > 2.0 → "kuyrukta ölür" (ayrı yazılır,
      açığın %83'ü orada).
    KIYAS: iki+ form survival ise Δ(χ²/dof) < 0.3 → ikisi de yaşar
      (zorla seçim YOK); aksi halde düşük χ²/dof kazanır.
    PARAMETRE KIYASLARI (ölüm sonrası, hayatta kalanlar için):
      H-Z1 τ_c ↔ 152'nin 0.68'i ve ikiz (τ_c→∞) limiti;
      H-Z3 σ ↔ A(τ) kalibrasyonunun λ_eff = 0.99'u (K4).

(6) K3 KÖPRÜ KARARI (ileri-hesap, parametresiz): ölçülen r(τ) profili,
    zarfın güç-payı Σ(â_gerçek² − â_ikiz²)/… üzerinden 180'in ΔM zarf-payını
    (~%80, log defteri) YENİDEN ÜRETİYOR mu? TUTMA EŞİĞİ: 180'in
    %70–88 aralığını (181: %87.6±8.5; 180: %81.5±11.3) kesişiyorsa
    "zarf kimliği MÜHÜRLENİR"; kesmezse köprü açık yazılır. (Formül 184d'de
    dondurulur — bu ön-kayıt yalnız ARALIĞI ve YÖNÜ sabitler.)

BU BETİKTE HESAPLANMAYAN (ŞART): hiçbir â_q, w, r, D, χ², σ, τ_c —
hepsi 184b…'de, bu dosya yazıldıktan SONRA.
Kullanım: 184a_onkayit.py
"""
import hashlib
import json
import subprocess
import time
from pathlib import Path

SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
           "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/184")
BU = Path(__file__).resolve()

KENAR = [0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.86]
KUYRUK_TAU = 0.70
NJACK = 8
GAZLAR = {"gercek": "son", "ikiz": "keskin"}   # veri_yukle adları
TAU_CIZGI = 0.86
TABAN = 0.40
CAP = 4000

HIPOTEZ = {
    "H-Z1_erfc":  {"form": "0.5*erfc((tau-tau_c)/Delta)", "par": ["tau_c", "Delta"]},
    "H-Z2_guc":   {"form": "1 - c*tau**alpha",             "par": ["c", "alpha"]},
    "H-Z3_gaussDW": {"form": "exp(-(2*pi*tau)**2 * sigma**2 / 2)", "par": ["sigma"]},
}
ESIK = {
    "survival_chi2_dof": 2.0,
    "kuyruk_z": 2.0,
    "kiyas_dchi2": 0.3,
    "K3_zarf_pay_araligi": [0.70, 0.88],   # 180/181 log defteri
}

DONMUS = {
    "aq_tanim": "a_q = 1/(pi*m*sqrt(q)), q=p^m, m=us",
    "aq_eff_nominal": "2*a_q*sin(pi*tau_q)",
    "chat_q": "c_q(ds) = 2*<ds_n * exp(-i*w_q*mid_n)>, mid=ham orta nokta",
    "ahat_q": "|c_q(ds)|",
    "w_q": "ahat_q / aq_eff_nominal",
    "r_tau": "ahat_gercek/ahat_ikiz  (nominal sadelesir; BIRINCIL K2 nesnesi)",
    "D_tau": "w_gercek - w_ikiz",
    "kanal": "ds (132c ozdesligi *; taban-usti tau>=0.45'te eta ile ozdes)",
    "kaynak_onbellek": "155/eta_<gaz>_t0.4_c4000.npz (mid,ds,L)",
    "capraz_kontrol": "176/K1c_<gaz>.npz c_eta ile |c_eta(q)|~|c_ds(q)| (tau>=0.45)",
}


def main():
    SCR.mkdir(parents=True, exist_ok=True)
    try:
        damga = subprocess.check_output(["date"], text=True).strip()
    except Exception:
        damga = time.strftime("%a %b %d %H:%M:%S %Z %Y")
    sha = hashlib.sha256(BU.read_bytes()).hexdigest()
    onk = dict(
        gorev="184 — ZARFIN ANATOMISI (K0 on-kayit)",
        zaman=damga, sha256=sha, betik=str(BU),
        donmus_tanim=DONMUS, kenar=KENAR, kuyruk_tau=KUYRUK_TAU,
        njack=NJACK, gazlar=GAZLAR, taban=TABAN, cap=CAP, tau_cizgi=TAU_CIZGI,
        hipotez=HIPOTEZ, esik=ESIK,
    )
    p = SCR / "ONKAYIT_184.json"
    p.write_text(json.dumps(onk, indent=1, ensure_ascii=False))
    print("=" * 70)
    print("184a / K0 ON-KAYIT  yazildi")
    print("=" * 70)
    print(f"  zaman : {damga}")
    print(f"  sha256: {sha}")
    print(f"  -> {p}")
    print(f"  bant kenar : {KENAR}")
    print(f"  hipotezler : {list(HIPOTEZ)}")
    print(f"  esikler    : {ESIK}")


if __name__ == "__main__":
    main()
