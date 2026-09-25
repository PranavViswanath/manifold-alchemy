"""Redox-ladder window per element at pH 0.3: scan U, record the MP-Pourbaix stable entry; the solid-oxide window
[U_red, U_ox] and the channels (reductive to lower-valent ion / oxidative to higher-valent oxyanion) with electron counts."""
import json, time, pathlib, numpy as np
from mp_api.client import MPRester
from pymatgen.analysis.pourbaix_diagram import PourbaixDiagram
ROOT = pathlib.Path(r"C:\Users\pvisw\Documents\CS\research\interp\sci disc\mlips"); key = open(ROOT / ".mp_api_key").read().strip()
ELS = ["Ag","Al","Bi","Co","Cr","Cu","Fe","Ge","Hf","In","Mn","Mo","Nb","Ni","Pb","Pd","Re","Sb","Si","Sn","Ta","Ti","V","W","Zn","Zr","Ru","Ir","Pt","Rh","Os","Au"]
OUT = DATA / "pbx_window.json"; res = json.load(open(OUT)) if OUT.exists() else {}
def client():
    for t in range(10):
        try: return MPRester(key)
        except Exception as ex: time.sleep(10)
m = client(); pH = 0.3; Us = np.round(np.arange(0.4, 2.61, 0.02), 2)
for el in ELS:
    if el in res: continue
    for t in range(6):
        try: ents = m.get_pourbaix_entries([el]); break
        except Exception as ex: print("retry", el, str(ex)[:60], flush=True); time.sleep(10); m = client()
    else: continue
    pd = PourbaixDiagram(ents, filter_solids=True)
    seq = []
    for U in Us:
        e = pd.get_stable_entry(pH, float(U)); seq.append((float(U), e.name, e.phase_type))
    # compress into segments
    segs = []
    for U, n, p in seq:
        if segs and segs[-1][1] == n: segs[-1][2] = U
        else: segs.append([U, n, U, p])
    res[el] = segs; json.dump(res, open(OUT, "w"), indent=1)
    print(f"{el:3s} " + " | ".join(f"{a:.2f}-{b:.2f} {n} ({p[0]})" for a, n, b, p in segs), flush=True)
print("\nSOLID WINDOW containing or nearest to U=1.68 V, pH 0.3:")
for el, segs in res.items():
    sol = [s for s in segs if s[3] == "Solid" and el != s[1].replace("(s)", "")]
    print(el, [(s[1], s[0], s[2]) for s in sol])
