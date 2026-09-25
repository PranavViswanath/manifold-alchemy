import json, pathlib, time
from mp_api.client import MPRester
ROOT = pathlib.Path(r"C:\Users\pvisw\Documents\CS\research\interp\sci disc\mlips"); key = open(ROOT / ".mp_api_key").read().strip()
win = json.load(open(DATA / "pbx_window.json")); out = {}
with MPRester(key) as m:
    for el, segs in win.items():
        # operating oxide = solid segment containing 1.68 V (excluding MO4/M2O7-type volatile oxides), else None
        f = None
        for lo, name, hi, ph in segs:
            if ph != "Solid" or name == el + "(s)": continue
            hi2 = 2.6 if hi >= 2.58 else hi
            if lo <= 1.68 <= hi2: f = name.replace("(s)", "")
        if f is None: out[el] = dict(oxide=None, gap=None); continue
        try:
            docs = m.materials.summary.search(formula=f, fields=["material_id", "band_gap", "energy_above_hull", "is_metal"])
            docs = sorted(docs, key=lambda d: d.energy_above_hull or 0)
            out[el] = dict(oxide=f, gap=float(docs[0].band_gap), ehull=float(docs[0].energy_above_hull or 0), n=len(docs))
        except Exception as ex: out[el] = dict(oxide=f, gap=None, err=str(ex)[:60])
        print(el, out[el], flush=True)
json.dump(out, open(DATA / "gaps.json", "w"), indent=1)
