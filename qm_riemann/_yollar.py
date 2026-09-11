# -*- coding: utf-8 -*-
"""
Taşınabilir yol yardımcısı (yeni scriptler için).

Eski scriptler mutlak yollar kullanıyor:
    QM  = Path("/Users/ugursezen/Desktop/arin/deney/qm_riemann")
    SCR = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/<UUID>/scratchpad")

Bu dosya aynı nesneleri taşınabilir biçimde verir:
    QM  = depo içindeki qm_riemann/ dizini
    SCR = ara çıktıların dizini (ortam değişkeniyle taşınabilir)

Kullanım:
    from _yollar import QM, SCR
    ya da
    QM_KOK=/başka/yer  QM_SCRATCH=/başka/scratchpad python3 script.py

NOT: 152-187 arku scriptleri sha256 ile dondurulmuş ön-kayıtlara sahip;
onları DÜZENLEMEK mührü kırar. Bu yüzden eski dosyalara dokunulmadı,
onun yerine 09_yol_koprular.sh ile sembolik köprü kuruluyor.
"""
import os
from pathlib import Path

QM = Path(os.environ.get("QM_KOK", Path(__file__).resolve().parent)).resolve()
SCR = Path(os.environ.get("QM_SCRATCH", QM / "scratchpad")).resolve()

def scr(gorev, *parcalar):
    """SCR/'NNN'/... yolunu üretir (klasör yoksa oluşturmaz)."""
    return SCR.joinpath(str(gorev), *parcalar)

__all__ = ["QM", "SCR", "scr"]
