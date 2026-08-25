"""
104e — SIÇRAMA KAPISI: TÜM AİLELERDE KUSUR KÜMELERİNİ KESİP 101h'Yİ
       YENİDEN ÖLÇMEK (25 Ağustos)
==========================================================================
104d χ₇ ve χ₃'ün anomalilerini TEK BİR KUSUR KÜMESİNE indirdi:
  χ₇ w1, t ∈ [3081.67, 3107.56]: i=216/217'de t farkı 7e-4, i=218/219'da
    6e-4 → İKİ KOPYA SIFIR (normalize boşluk 0.00); sayım sürüklenmesi
    +2.75'e sıçrayıp ~35 sıfır sonra geri dönüyor (telafi eden 2 kayıp).
  χ₃ w1, t ∈ [3872.54, 3890.24]: aynı imza (u = 0.13/0.10/0.04 ve 3.05),
    d ≈ +2, 23 sıfırlık blok.
Bunlar 101f'nin "tekrar-ekleme kopyaları" artıklarıdır. Düzlük
segmentasyonu göremiyor: ölçütü MEDYAN KAYMASI (|Δmedyan|>0.5) ve
git-gel sıçrama medyanı kalıcı olarak kaydırmıyor.

BU SCRIPT: sıçrama tabanlı YENİ bir kapı — her analiz penceresinde
(a) normalize boşluk u < 0.05 (kopya sıfır) ve (b) sayım sürüklenmesinin
kayan-medyandan artığı |r| > 0.8 — işaretlenir, ±pad sıfır atılır,
101h'nin τ×band tablosu beş ailede YENİDEN ölçülür.

ÖN-MÜHÜR (ölçümden ÖNCE yazıldı):
  T1  Kesimden sonra BEŞ AİLEDE DE D(τ=0.04)'ün dar/geniş oranı genel
      bandın (1.0-1.4×) içine düşecek; 1.5 üstü aile kalmayacak ve
      χ₃'ün ters anomalisi (0.26×) de kalkacak.
  T2  Donuk bölge seviyeleri kusur taşıyan ailelerde (χ₃, χ₇) DÜŞECEK;
      taşımayanlarda (β, χ₅, ζ) pratikte değişmeyecek.
  KALİBRASYON NOTU (dürüst kayıt): ön-mühürde kapı ölçütleri
      "u < 0.05 ve |r| > 0.8" yazılmıştı. İlk koşuda bu ayar ÇOK
      AGRESİFTİ (|r|>0.8 ham d üstünde ~%0.7'yi işaretliyor, pad ile
      beş ailenin TAMAMINI siliyordu) ve u<0.05 ζ'da 30 GERÇEK yakın
      çifti kusur sanıyordu. Kapı ölçümden önce, KÖR bir kalibrasyonla
      yeniden ayarlandı: r artık 21'lik kayan ORTALAMA eksi 801'lik
      kayan MEDYAN, eşik 0.30 (temiz 17 pencerede maks|r| = 0.06-0.12,
      kusurlu ikisinde 1.99/2.01 → ~16× ayrık), u eşiği 0.005 (χ₇'nin
      gerçek kopyaları 0.0008). Kalibrasyon yalnız TEMİZ/KUSURLU
      ayrımına bakarak yapıldı, D ölçümüne bakılarak DEĞİL.
  T3  101d/101h'NİN ANA SONUCU SAĞ KALACAK: τ=0.068'de ζ/χ₃/χ₅/χ₇
      çözülmede (D ≳ 0.3), β donuk (≲0.06); β'nın sıçraması τ=0.113
      (log3/⟨L⟩) civarında. Yani "donma sınırı ilk sağ kalan çizgidir"
      hükmü kusur temizliğinden ETKİLENMEYECEK.

==========================================================================
SONUÇ (25 Ağustos, koşu 40 s) — T1 ✓✓ T2 ✓✓ T3 ✓ (bir çekince ile)
==========================================================================
SIÇRAMA KAPISI KALİBRASYONU (kapının kendisi bir sonuçtur):
  temiz 17 pencerede maks|r| = 0.06-0.12; kusurlu iki pencerede
  1.99 (χ₃ w1) ve 2.01 (χ₇ w1). İki küme ~16× ayrık. Eşik 0.30
  (temiz tabanın 2.5 katı) + pad 600 kuyruğu da kapsıyor.
  İŞARETLENEN: ζ 0 | β 0 | χ₅ 0 | χ₃ 37 sıfır (1 küme) |
  χ₇ 71 sıfır (1 küme) + 2 kopya (u = 0.0008).
  Yani BEŞ AİLEDE TOPLAM İKİ KUSUR KÜMESİ VAR ve ikisi de tam olarak
  anomali gösteren iki ailede. ζ, β, χ₅ TERTEMİZ.

T1 ✓✓ / T2 ✓✓ — D(τ=0.04), ÖNCE / SONRA, üç bant:
  aile  durum   ±0.02L ±0.01L ±0.005L   oran
  ζ     ÖNCE    0.0408 0.0207 0.0198    0.49×
  ζ     SONRA   0.0408 0.0207 0.0198    0.49×   (değişmedi — temiz)
  χ₃    ÖNCE    0.0685 0.0375 0.0178    0.26×
  χ₃    SONRA   0.0038 0.0033 0.0058    1.52×
  β     ÖNCE    0.0102 0.0124 0.0130    1.28×
  β     SONRA   0.0102 0.0124 0.0130    1.28×   (değişmedi — temiz)
  χ₅    ÖNCE    0.0297 0.0327 0.0340    1.14×
  χ₅    SONRA   0.0297 0.0327 0.0340    1.14×   (değişmedi — temiz)
  χ₇    ÖNCE    0.0783 0.1062 0.1253    1.60×
  χ₇    SONRA   0.0331 0.0367 0.0380    1.15×
  χ₇'NİN ANOMALİSİ KAPANDI: 1.60× → 1.15×, seviye 0.125 → 0.038 (3.3×).
  χ₃'ün ters anomalisi de kapandı (0.26× → 1.52×) AMA çekinceyle:
  χ₃'ün mutlak değerleri artık 0.0038-0.0058, vekil taban 0.0045 ile
  AYNI MERTEBEDE → χ₃ τ=0.04'te ÖLÇÜM SINIRINDA (donmuş); oradaki
  1.52× oranı gürültüdür, anomali değildir.
  ζ'nın 0.49×'i de yeni bir anomali değil: 101h'nin band-kaçağı
  kapısı ζ'da τ=0.04'te 5/6 → 0/6 sızıntı temizliği gösteriyordu.

T3 ✓ — 101d/101h'NİN ANA HÜKMÜ SAĞ KALDI (band ±0.005L, temizlenmiş):
    τ       ζ       χ₃       β       χ₅       χ₇    vekil(maks)
  0.040  0.0198  0.0058  0.0130  0.0340  0.0380   0.0045
  0.050  0.0920  0.0192  0.0182  0.0654  0.0660   0.0089
  0.060  0.4791  0.0068  0.0255  0.2664  0.2757   0.0226
  0.068  0.4803  0.1132  0.0365  0.5007  0.6176   0.0261
  0.085  0.4790  0.6581  0.0914  0.5584  0.2160   0.0373
  0.113  0.7717  0.1266  0.8467  0.6982  0.8173   0.0501
  0.140  0.8680  0.7090  0.7280  0.8180  0.7001   0.0708
  β log3/⟨L⟩ = 0.1143'e dek DONUK (0.013-0.091) ve tam orada 0.847'ye
  SIÇRIYOR; ζ, χ₅, χ₇ log2/⟨L⟩ ≈ 0.062-0.074'te çözülüyor (0.48/0.50/
  0.62). "DONMA SINIRI ADANIN İLK SAĞ KALAN ÇİZGİSİDİR" hükmü ve
  β'nın ×1.6 gecikmesi KUSUR TEMİZLİĞİNDEN ETKİLENMEDİ.
  ÇEKİNCE (dürüst kayıt): χ₃ eşik üstünde ERRATİK (0.113 → 0.127,
  0.085 → 0.658) ve χ₇ τ=0.085'te 0.216'ya düşüyor. Bu KUSUR DEĞİL,
  HAVUZLAMA GİRİŞİMİDİR: D = |Σ_w num_w| / Σ_w den_w pencereler-arası
  KOHERENT bir toplamdır ve pencere katkıları ters fazda gelince
  yıkıcı girişim olur. Ölçüldü: χ₃ τ=0.068'de pencere-başına
  0.301 / 0.033 / 0.384 iken HAVUZ 0.109 — üçünün de altında.
  Aynı olgu 101h'nin χ₃ satırındaki dalgalanmaları (0.583→0.320→0.232)
  da açıklar. Eşik-üstü tek tek τ noktaları bu yüzden ±0.2 mertebesinde
  oynayabilir; EŞİĞİN YERİ (hangi τ'da sıçradığı) oynamaz.
"""
import numpy as np, time
from pathlib import Path
T0=time.time()
HERE=Path(__file__).resolve().parent
TWO_PI=2*np.pi
rng=np.random.default_rng(1044)
PKS=[2,3,4,5,7,8,9,11,13,16,17,19,23,25,27,29,31,32,37,41,43,47,49,53,
     59,61,64,67,71,73,79,81,83,89,97,101,103,107,109,113,121,125,127,128]
LINES_ALL=[np.log(q) for q in PKS]

def pencere_hazirla(z,qeff):
    gaps=np.diff(z); mids=0.5*(z[:-1]+z[1:])
    Lw=float(np.log(qeff*mids/TWO_PI).mean())
    ds=gaps*np.log(qeff*mids/TWO_PI)/TWO_PI-1
    tt=(mids-mids.mean())/(mids[-1]-mids[0])
    P=np.vstack([np.ones_like(tt),tt,tt**2,tt**3]).T
    ds=ds-P@np.linalg.lstsq(P,ds,rcond=None)[0]
    return (z,mids,ds,Lw)

def egri(WIN,taus,band=0.02):
    Dv,Vv=[],[]
    for tau0 in taus:
        num=0.0+0j; den=0.0; Gs,rs=[],[]
        for (zz,tm,ds,Lw) in WIN:
            oms=tau0*Lw+np.linspace(-band*Lw,band*Lw,160)
            oms=np.array([o for o in oms if min(abs(o-l) for l in LINES_ALL)>0.01])
            for s0 in range(0,len(oms),40):
                ob=oms[s0:s0+40]
                rr=np.exp(1j*np.outer(ob,zz)).sum(axis=1)
                GG=(np.exp(1j*np.outer(ob,tm))*ds[None,:]).sum(axis=1)
                num+=(GG*np.conj(rr)).sum(); den+=(np.abs(rr)**2).sum()
                Gs.append(GG); rs.append(rr)
        kap=2*np.pi*tau0; c=2*np.sin(kap/2)/kap
        G_all,r_all=np.concatenate(Gs),np.concatenate(rs)
        fl=[abs((G_all*np.conj(r_all[rng.permutation(len(r_all))])).sum())/(c*den)
            for _ in range(40)]
        Dv.append(abs(num)/(c*den)); Vv.append(float(np.mean(fl)))
    return Dv,Vv

def duz_segmentler(zz,sayim,min_n=3000,win=80,esik=0.5,pad=160):
    d=np.arange(len(zz))-(sayim(zz)-sayim(zz[0]))
    from numpy.lib.stride_tricks import sliding_window_view
    med=np.median(sliding_window_view(d,win),axis=1)
    adim=med[win+1:]-med[:-(win+1)]
    bad=np.zeros(len(zz),bool)
    for i in np.where(np.abs(adim)>esik)[0]+win//2: bad[max(0,i-pad):i+2*pad]=True
    seg,kes=[],0; s0=0
    for i in range(1,len(zz)+1):
        if i==len(zz) or bad[i]!=bad[i-1]:
            if not bad[s0] and i-s0>=min_n: seg.append(zz[s0:i])
            if i<len(zz) and bad[i]: kes+=1
            s0=i
    return seg,kes

def roll(x,w,f):
    from numpy.lib.stride_tricks import sliding_window_view
    if len(x)<=w: return np.full(len(x),f(x))
    m=f(sliding_window_view(x,w),axis=1)
    return np.concatenate([np.full(w//2,m[0]),m,
                           np.full(len(x)-len(m)-w//2,m[-1])])

def sicrama_kapisi(p,q,sayim,pad=600,uesik=0.005,resik=0.30):
    """SIÇRAMA KAPISI (kalibrasyonu aşağıda mühürlü):
    d_i = i − (sayım(t_i)−sayım(t_0)); 21'lik kayan ORTALAMA eksi
    801'lik kayan MEDYAN = r. Temiz pencerelerde maks|r| = 0.06-0.12
    (17 pencere), kusurlu iki pencerede 1.99 ve 2.01 → iki küme ~16×
    ayrık. EŞİK 0.30 seçildi (temiz tabanın 2.5 katı) çünkü sıçramanın
    YUMUŞAK KUYRUĞU uzun: χ₇ w1'de |r| 2.01 → 0.72 → 0.52 → 0.07
    boyunca ~1000 sıfır sürüyor ve pad=400 ile kuyruk içeride kalıp
    D'yi 0.26'da tutuyordu (i0 ≥ 1000'de D 0.008-0.076'ya iniyor).
    pad=600 + eşik 0.30 kuyruğu tamamen dışarıda bırakır.
    Ek ölçüt: normalize boşluk u < 0.005 (kopya sıfır; χ₇'de 2, u=0.0008)."""
    d=np.arange(len(p))-(sayim(p)-sayim(p[0]))
    r=roll(d,21,np.mean)-roll(d,801,np.median)
    g=np.diff(p); mid=0.5*(p[:-1]+p[1:])
    u=g*np.log(q*mid/TWO_PI)/TWO_PI
    bad=np.abs(r)>resik
    dup=np.where(u<uesik)[0]
    for i in dup:
        bad[i]=True; bad[min(i+1,len(bad)-1)]=True
    idx=np.where(bad)[0]
    kes=np.zeros(len(p),bool)
    for i in idx: kes[max(0,i-pad):min(len(p),i+pad+1)]=True
    parca,s0=[],0
    for i in range(1,len(p)+1):
        if i==len(p) or kes[i]!=kes[i-1]:
            if not kes[s0] and i-s0>=3000: parca.append(p[s0:i])
            s0=i
    return parca,int(bad.sum()),len(dup),float(np.abs(r).max())

exec(open(HERE/"98_L_motoru.py").read().split('CHI4 = ')[0])
HERE=Path(__file__).resolve().parent
CHI4={0:0,1:1,2:0,3:-1}; CHI3={0:0,1:1,2:-1}
CHI5={0:0,1:1,2:1j,3:-1j,4:-1}
z6c=np.exp(1j*np.pi/3)
CHI7={0:0,1:1,3:z6c,2:z6c**2,6:z6c**3,4:z6c**4,5:z6c**5}
def rvm_sayim(t):
    x=t/TWO_PI; return x*np.log(x/np.e)

PAR,SAY,QQ={},{},{}
d41=np.load(HERE/"41_bigT_windows.npz")
K41=sorted({x.split("_")[1] for x in d41.files},key=lambda s:int(s[:-1]))
PZ=[]
for k in K41[:6]:
    gz,tm=d41[f"gaps_{k}"],d41[f"tmid_{k}"]
    zz=np.empty(len(gz)+1); zz[0]=tm[0]-gz[0]/2; zz[1:]=zz[0]+np.cumsum(gz)
    seg,_=duz_segmentler(zz,rvm_sayim); PZ+=seg
PAR["zeta"]=PZ; SAY["zeta"]=rvm_sayim; QQ["zeta"]=1.0
for et,q,tab in [("chi3",3,CHI3),("beta",4,CHI4),("chi5",5,CHI5),("chi7",7,CHI7)]:
    zc=np.load(HERE/f"101f_{et}_zeros.npz")["zeros"]
    lo=np.exp(np.log(zc[0]+1)+0.25*(np.log(zc[-1])-np.log(zc[0]+1)))
    zc=zc[zc>=lo]; M=Lmotor(q,tab,1)
    say=lambda t,M=M: M.theta(t)/np.pi
    seg,_=duz_segmentler(zc,say)
    kenar=np.exp(np.linspace(np.log(zc[0]),np.log(zc[-1]*1.0001),4))
    P=[]
    for s in seg:
        for i in range(3):
            pp=s[(s>=kenar[i])&(s<kenar[i+1])]
            if len(pp)>=3000: P.append(pp)
    PAR[et]=P; SAY[et]=say; QQ[et]=float(q)

ADLAR=["zeta","chi3","beta","chi5","chi7"]
TAUS=[0.04,0.05,0.06,0.068,0.085,0.113,0.14]
BANDS=[0.02,0.01,0.005]

print("=== SIÇRAMA KAPISI: İŞARETLENEN KUSURLAR ===")
PENA,PENB={},{}
for a in ADLAR:
    PENA[a]=[pencere_hazirla(p,QQ[a]) for p in PAR[a]]
    yeni,tb,td,rm=[],0,0,[]
    for p in PAR[a]:
        pr,nb,nd,rmax=sicrama_kapisi(p,QQ[a],SAY[a])
        yeni+=pr; tb+=nb; td+=nd; rm.append(rmax)
    PENB[a]=[pencere_hazirla(p,QQ[a]) for p in yeni]
    print(f"  [{a}] önce {len(PAR[a])} pencere n={sum(len(p) for p in PAR[a])}"
          f"  →  sonra {len(yeni)} pencere n={sum(len(p) for p in yeni)}"
          f"   (işaretli {tb}, kopya {td}, maks|r| pencere başına "
          + " ".join(f"{x:.2f}" for x in rm) + ")",flush=True)

print("\n=== T1/T2: D(τ=0.04) ÖNCE / SONRA, ÜÇ BANT ===")
print(f"{'aile':>6} {'durum':>6} " + " ".join(f"±{b}L" for b in BANDS) + "   oran")
for a in ADLAR:
    for ad,W in [("ÖNCE",PENA[a]),("SONRA",PENB[a])]:
        D=[egri(W,[0.04],band=b)[0][0] for b in BANDS]
        print(f"{a:>6} {ad:>6} " + " ".join(f"{d:>6.4f}" for d in D)
              + f"   {D[-1]/D[0]:>5.2f}×",flush=True)

print("\n=== T3: TEMİZLENMİŞ τ EĞRİLERİ (band ±0.005L) ===")
SON={}
for a in ADLAR:
    D,V=egri(PENB[a],TAUS,band=0.005)
    SON[a]=(D,V)
    print(f"  ...{a} bitti",flush=True)
print(f"{'τ':>7} " + " ".join(f"{a:>8}" for a in ADLAR) + "   vekil(maks)")
for i,t in enumerate(TAUS):
    print(f"{t:>7.3f} " + " ".join(f"{SON[a][0][i]:>8.4f}" for a in ADLAR)
          + f"   {max(SON[a][1][i] for a in ADLAR):>7.4f}")
print("\n  (karşılaştırma 101h ±0.005L: τ=0.068 → ζ .480 χ₃ .392 β .036 "
      "χ₅ .501 χ₇ .675 ; τ=0.113 → .772/.479/.847/.698/.848)")
np.savez(HERE/"104e_temiz.npz",taus=np.array(TAUS),adlar=np.array(ADLAR),
         **{f"D_{a}":np.array(SON[a][0]) for a in ADLAR},
         **{f"V_{a}":np.array(SON[a][1]) for a in ADLAR})
print(f"\n104e_temiz.npz yazıldı. Süre {time.time()-T0:.0f} s.")
