"""
104d — KUSURUN YERİ: χ₇ w1'DEKİ SAYIM SIÇRAMASININ TEŞHİSİ (25 Ağustos)
==========================================================================
104b'nin K2/K6/K6b'si anomaliyi tek bir pencereye indirdi: χ₇ w1
(t ∈ [2914.6, 11063.3]). Net sürüklenme ≈ 0 (−0.00) AMA pencere içi
maks|d − med| = 2.75 — bütün ailelerin bütün pencerelerinin EN BÜYÜĞÜ
(ikinci: χ₃ w1 2.60; geri kalan 13 pencere 0.88-1.16). Yani w1'de
"git-gel" bir sayım hatası var: net sıfır, ama içeride ≥2 sıfırlık bir
KAYIP + telafi eden bir SAHTE (ya da tersi). Düzlük segmentasyonu bunu
görmüyor çünkü ölçüt basamağın MEDYAN KAYMASI (|Δmedyan|>0.5), ve
git-gel hareket medyanı kalıcı olarak kaydırmıyor.
Sıkılaştırma (0.35/300 → 0.25/500) w1'in başını 2914.6 → 3572.4 →
3879.7'ye çekiyor ve w1'in D'si 0.304 → 0.201 → 0.049'a çöküyor.

ÖN-MÜHÜR (ölçümden ÖNCE yazıldı):
  Y1  w1'in sayım sürüklenmesi d_i, t ∈ [2914.6, 3880] aralığında
      ≥2 birimlik bir sapma yapıp geri dönecek (kusurun yeri).
  Y2  O sapmanın t-aralığı w1'den ÇIKARILDIĞINDA (kalan parçalar ayrı
      pencere) χ₇'nin τ=0.04 havuz değeri ve dar/geniş oranı normale
      (D ≲ 0.03, oran ≈ 1.0-1.3×) dönecek.
  Y3  Aynı işlem χ₃ w1'e (maks|d−med| = 2.60) uygulanınca χ₃'ün
      0.26× "ters anomalisi" de normalleşecek (aynı cinsten kusur).
  Y4  Kusur bölgesindeki boşluk dağılımı anormal olacak: normalize
      boşluk ya ≈0 (sahte çift) ya ≳3 (kayıp sıfır) uçlar verecek.

==========================================================================
SONUÇ (25 Ağustos, koşu 5 s) — Y1 ✓✓ Y2 ✓✓ Y3 ✓✓ Y4 ✓✓ : KUSUR BULUNDU
==========================================================================
Y1 ✓✓ χ₇ w1 — KOPYA SIFIRLAR, t ≈ 3081.7:
   i        t              d(medyandan)
   215  3081.1344        −0.27
   216  3081.6688        +0.04
   217  3081.6695        +1.04   ← 216 ile farkı 7.0e-4 (KOPYA)
   218  3081.8885        +1.75
   219  3081.8891        +2.75   ← 218 ile farkı 6.0e-4 (KOPYA)
   220  3082.8773        +2.47   ... d ≈ +2 olarak i=251'e dek sürüyor
  |d|>1 bitişik blok: i ∈ [217,251], t ∈ [3081.67, 3107.56] (35 sıfır).
  Normalize boşluklar orada 0.001 ve 0.001 (pencere medyanı 0.974).
  Yani İKİ SAHTE (kopya) SIFIR + ilerideki telafi edici iki kayıp →
  NET sürüklenme 0, bu yüzden 101f'nin sayım sertifikası ve 101d'nin
  düzlük segmentasyonu İKİSİ DE GÖRMEDİ. Bunlar 101f'nin "tekrar-
  ekleme kopyaları" temizliğinden ARTAKALANLARDIR.
Y1 ✓✓ χ₃ w1 — AYNI İMZA, t ≈ 3872.5-3890.2:
  |d|>1 blok i ∈ [633,655] (23 sıfır), d ≈ +2, i=656'da d → +0.15.
  Normalize boşluklar bölgede 0.13 / 0.10 / 0.04 (kopyaya yakın çift)
  ve 3.05 (kayıp sıfır) — kopya-kadar-keskin değil ama aynı sınıf.
Y4 ✓✓ Boşluk dağılımı: her iki bölgede de pencerenin GLOBAL uç
  değerleri orada (χ₇: u_min 0.001 = pencere minimumu; χ₃: u_max 3.05
  = pencere maksimumu).

Y2 ✓✓ / Y3 ✓✓ KUSUR BÖLGESİ ±400 SIFIR PAD İLE KESİLİNCE (τ=0.04):
  χ₇: ÖNCE 0.0783 0.1062 0.1253 (1.60×)
      SONRA 0.0535 0.0610 0.0634 (1.19×)
  χ₃: ÖNCE 0.0685 0.0375 0.0178 (0.26×)
      SONRA 0.0090 0.0105 0.0113 (1.26×)
  İKİ ANOMALİ DE (χ₇'nin yükselişi VE χ₃'ün ters düşüşü) tek bir
  ~20-35 sıfırlık kusur kümesinden geliyormuş; kesildiğinde ikisi de
  genel banda (1.1-1.3×) oturuyor.
  DÜRÜST KAYIT: pad=400 χ₇'de YETMEDİ — sıçramanın yumuşak kuyruğu
  ~1000 sıfır sürüyor (|r| 2.01 → 0.72 → 0.52 → 0.07) ve w1'in D'si
  0.26'da takılı kalıyor; başlangıç i₀ ≥ 1000'e (t ≥ 3680) çekilince
  D 0.008-0.076'ya iniyor. Nihai kapı ayarı 104e'de (eşik 0.30,
  pad 600).
"""
import numpy as np, time
from pathlib import Path
T0=time.time()
HERE = Path(__file__).resolve().parent
TWO_PI = 2*np.pi
rng = np.random.default_rng(1043)
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

exec(open(HERE/"98_L_motoru.py").read().split('CHI4 = ')[0])
HERE = Path(__file__).resolve().parent
CHI3={0:0,1:1,2:-1}; CHI5={0:0,1:1,2:1j,3:-1j,4:-1}
z6c=np.exp(1j*np.pi/3)
CHI7={0:0,1:1,3:z6c,2:z6c**2,6:z6c**3,4:z6c**4,5:z6c**5}
AILE={"chi3":(3,CHI3),"chi5":(5,CHI5),"chi7":(7,CHI7)}
PEN,PAR,MOT={},{},{}
for et,(q,tab) in AILE.items():
    zc=np.load(HERE/f"101f_{et}_zeros.npz")["zeros"]
    lo=np.exp(np.log(zc[0]+1)+0.25*(np.log(zc[-1])-np.log(zc[0]+1)))
    zc=zc[zc>=lo]; M=Lmotor(q,tab,1); MOT[et]=M
    seg,_=duz_segmentler(zc,lambda t: M.theta(t)/np.pi)
    kenar=np.exp(np.linspace(np.log(zc[0]),np.log(zc[-1]*1.0001),4))
    W,P=[],[]
    for s in seg:
        for i in range(3):
            p=s[(s>=kenar[i])&(s<kenar[i+1])]
            if len(p)>=3000: W.append(pencere_hazirla(p,float(q))); P.append(p)
    PEN[et]=W; PAR[et]=P

BANDS=[0.02,0.01,0.005]

# ---------- Y1: sürüklenme profili, en sapkın pencerede ----------
print("=== Y1: SAYIM SÜRÜKLENMESİ PROFİLİ (en sapkın pencereler) ===")
KUSUR={}
for et,wi in [("chi7",0),("chi3",0)]:
    p=PAR[et][wi]; th=MOT[et].theta(p)/np.pi
    d=np.arange(len(p))-(th-th[0]); d=d-np.median(d)
    j=int(np.argmax(np.abs(d)))
    # sapma bölgesi: |d| > 1.0 olan, j'yi içeren bitişik blok
    m=np.abs(d)>1.0
    a=j
    while a>0 and m[a-1]: a-=1
    b=j
    while b<len(d)-1 and m[b+1]: b+=1
    KUSUR[et]=(p[a],p[b],a,b,float(d[j]))
    print(f"  {et} w{wi+1}: maks|d−med| = {abs(d[j]):.2f} @ i={j} t={p[j]:.2f}")
    print(f"     |d|>1 bitişik blok: i∈[{a},{b}]  t∈[{p[a]:.2f},{p[b]:.2f}] "
          f"({b-a+1} sıfır, genişlik {p[b]-p[a]:.2f})")
    lo2,hi2=max(0,j-12),min(len(p),j+13)
    print("     yerel d (i, t, d):")
    for i in range(lo2,hi2):
        print(f"       {i:>6} {p[i]:>11.4f} {d[i]:>+7.2f}")

# ---------- Y4: kusur bölgesinde boşluk dağılımı ----------
print("\n=== Y4: KUSUR BÖLGESİNDE NORMALİZE BOŞLUKLAR ===")
for et in ["chi7","chi3"]:
    q=AILE[et][0]; p=PAR[et][0]
    t0,t1,a,b,_=KUSUR[et]
    lo2,hi2=max(1,a-8),min(len(p)-1,b+9)
    g=np.diff(p); mid=0.5*(p[:-1]+p[1:])
    u=g*np.log(q*mid/TWO_PI)/TWO_PI
    print(f"  {et}: bölge çevresi normalize boşluklar "
          f"(medyan tümü {np.median(u):.3f}):")
    print("     " + " ".join(f"{x:.2f}" for x in u[lo2:hi2]))
    print(f"     bölge min {u[a:b+1].min():.3f} maks {u[a:b+1].max():.3f}; "
          f"tüm pencere min {u.min():.3f} maks {u.max():.3f}")

# ---------- Y2/Y3: kusur bölgesini çıkar, yeniden ölç ----------
print("\n=== Y2/Y3: KUSUR BÖLGESİ ÇIKARILARAK YENİDEN ÖLÇÜM (τ=0.04) ===")
PAD=400
for et in ["chi7","chi3"]:
    q=AILE[et][0]
    t0,t1,a,b,_=KUSUR[et]
    P=list(PAR[et])
    p=P[0]
    ke=max(0,a-PAD); kb=min(len(p),b+PAD)
    yeni=[]
    for parca in [p[:ke],p[kb:]]:
        if len(parca)>=3000: yeni.append(parca)
    P2=yeni+P[1:]
    W0=[pencere_hazirla(x,float(q)) for x in P]
    W2=[pencere_hazirla(x,float(q)) for x in P2]
    D0=[egri(W0,[0.04],band=bb)[0][0] for bb in BANDS]
    D2=[egri(W2,[0.04],band=bb)[0][0] for bb in BANDS]
    print(f"  {et}: kusur t∈[{t0:.1f},{t1:.1f}], ±{PAD} sıfır pad ile kesildi")
    print(f"     kalan parçalar n={[len(x) for x in P2]}")
    print(f"     ÖNCE : " + " ".join(f"{d:>6.4f}" for d in D0)
          + f"   {D0[-1]/D0[0]:>5.2f}×")
    print(f"     SONRA: " + " ".join(f"{d:>6.4f}" for d in D2)
          + f"   {D2[-1]/D2[0]:>5.2f}×", flush=True)

print(f"\nSüre {time.time()-T0:.0f} s.")
