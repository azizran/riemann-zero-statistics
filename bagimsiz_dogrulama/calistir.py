# -*- coding: utf-8 -*-
"""
SUDO'SUZ KOŞTURUCU — 152-187 scriptlerini diskteki baytlarına DOKUNMADAN çalıştırır.
Yerelleştirilmiş kaynak, script'in YANINDA geçici bir .kosu.py dosyasına yazılır ve
`__main__` olarak koşulur (multiprocessing spawn'ın `_init`'i bulabilmesi için).
Orijinal dosya diskte değişmez → sha256 ön-kayıt mühürleri korunur.

Kullanım:  python3 calistir.py <script.py> [argümanlar...]
"""
import sys, pathlib, runpy, os

QM  = pathlib.Path("/Users/ugur/Desktop/Deney/qm_riemann")
SCR = QM / "scratchpad"
DISI_QM  = "/Users/ugursezen/Desktop/arin/deney/qm_riemann"
DISI_SCR = ("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
            "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad")

orijinal = pathlib.Path(sys.argv[1]).resolve()
kaynak = orijinal.read_text(encoding="utf-8")
n1 = kaynak.count(DISI_QM); n2 = kaynak.count(DISI_SCR)
kaynak = kaynak.replace(DISI_QM, str(QM)).replace(DISI_SCR, str(SCR))
gecici = orijinal.with_name(orijinal.stem + ".kosu.py")
gecici.write_text(kaynak, encoding="utf-8")
print(f"[calistir] {orijinal.name}: {n1} QM + {n2} scratchpad yolu yerelleştirildi → {gecici.name}", flush=True)
SCR.mkdir(parents=True, exist_ok=True)
sys.argv = [str(gecici)] + sys.argv[2:]
try:
    runpy.run_path(str(gecici), run_name="__main__")
finally:
    try: gecici.unlink()
    except OSError: pass
