"""200-C: LMFDB (Platt) sıfır dosyası okuyucu — YALNIZ KONUMLAR (Z/M hesaplamaz).
Biçim: 8 bayt blok sayısı; her blok: 32 bayt başlık (t0, t1: double; Nt0, Nt1: uint64 birikimli sayım),
ardından (Nt1−Nt0) adet 13 baytlık kayıt: küçük-uçlu 104-bit tamsayı × 2^-101 = bir önceki sıfıra olan aralık
(bloğun ilk kaydı t0'dan). Sıfırlar (taban = ilk bloğun t0'ı, TAMSAYI olması beklenir) + ofset olarak döner."""
import numpy as np, struct
def oku(yol, n_max=None):
    with open(yol,'rb') as f:
        nblok=struct.unpack('<Q',f.read(8))[0]
        taban=None; ofs=[]; say=0; bloklar=[]
        for _ in range(nblok):
            t0,t1,N0,N1=struct.unpack('<ddQQ',f.read(32)); m=N1-N0
            raw=np.frombuffer(f.read(13*m),dtype=np.uint8).reshape(m,13)
            lo=raw[:,:8].copy().view('<u8').ravel().astype(np.float64)
            hi=np.zeros(m,dtype=np.uint64)
            for k in range(5): hi|=raw[:,8+k].astype(np.uint64)<<np.uint64(8*k)
            gap=hi.astype(np.float64)*2.0**-37+lo*2.0**-101
            if taban is None: taban=t0
            ofs.append((t0-taban)+np.cumsum(gap)); bloklar.append((t0,t1,N0,N1)); say+=m
            if n_max and say>=n_max: break
    z=np.concatenate(ofs)
    return taban, (z[:n_max] if n_max else z), bloklar
