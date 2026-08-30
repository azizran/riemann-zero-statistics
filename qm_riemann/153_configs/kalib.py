"""I1 kalibrasyonu: z_taban uzerinde (eps, nsup) tarayip SADECE kucuk-s
gostergesini ve var(ds)'yi olcer. Tam olcum zinciri KOSULMAZ (ucuz tarama)."""
import numpy as np
import sys
from pathlib import Path
sys.path.insert(0, "/Users/ugursezen/Desktop/arin/deney/qm_riemann/153_configs")
import importlib.util
spec = importlib.util.spec_from_file_location(
    "g153", "/Users/ugursezen/Desktop/arin/deney/qm_riemann/153_configs/153_gaz.py")
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)

TWO_PI = 2 * np.pi
S = Path("/private/tmp/claude-501/-Users-ugursezen-Desktop-arin/"
         "71b922d7-9b7c-485d-9927-3b0db425e2a2/scratchpad/153")
z0 = np.load(S / "z_taban.npy")


def stat(z):
    g = np.diff(z); mid = 0.5 * (z[:-1] + z[1:]); Lw = np.log(mid / TWO_PI)
    ds = g * Lw / TWO_PI - 1; s = ds + 1
    return np.var(ds), np.mean(s < 0.3), np.mean(s < 0.1), np.mean(s < 0.5)


def sessiz(*a, **k):
    pass


print("GERCEK           : var(ds)=0.16744  P(s<0.3)=0.02420  "
      "P(s<0.1)=0.00096  P(s<0.5)=0.10385")
v, p3, p1, p5 = stat(z0)
print(f"taban (eps=0)    : var(ds)={v:.5f}  P(s<0.3)={p3:.5f}  "
      f"P(s<0.1)={p1:.5f}  P(s<0.5)={p5:.5f}")
print()
for nsup in (3, 5, 10):
    for eps in (-0.1, -0.03, -0.01, -0.003, 0.003, 0.01, 0.03, 0.1):
        try:
            z = m.itme_supur(z0.copy(), eps, nsup, log=sessiz)
        except SystemExit as e:
            print(f"eps={eps:+.3f} nsup={nsup:2d} : DUSTU ({e})")
            continue
        v, p3, p1, p5 = stat(z)
        print(f"eps={eps:+.3f} nsup={nsup:2d} : var(ds)={v:.5f}  "
              f"P(s<0.3)={p3:.5f}  P(s<0.1)={p1:.5f}  P(s<0.5)={p5:.5f}",
              flush=True)
    print()
