import json, numpy as np, re
from scipy.interpolate import RBFInterpolator
from pymatgen.core.periodic_table import Element
from pymatgen.core import Composition
import pathlib; DATA = pathlib.Path(__file__).parent / "data"
out = json.load(open(DATA / "twodoor_chart.json")); win = json.load(open(DATA / "pbx_window.json")); gaps = json.load(open(DATA / "gaps.json"))
def volatile(name, el):
    c = Composition(name.replace("(s)", "")); 
    try: return c["O"] / c[el] > 2.5
    except Exception: return False
def margin(el):
    best = -1.0
    for lo, name, hi, ph in win[el]:
        if ph != "Solid" or name == el + "(s)" or volatile(name, el): continue
        hi2 = 2.6 if hi >= 2.58 else hi
        if lo <= 1.68 <= hi2: best = max(best, min(1.68 - lo, hi2 - 1.68))
        else: best = max(best, -min(abs(1.68 - lo), abs(1.68 - hi2)))
    return min(best, 1.0)
def cond(el):
    g = gaps.get(el, {}); 
    if g.get("oxide") is None or volatile(g["oxide"], el): return -1.0
    return 1.0 if (g["gap"] <= 0.3 or el == "Pt") else -1.0
els = [e for e in win if e not in ("Al", "Si")]  # keep p/d elements with (group,row) defined
gp = np.array([[Element(e).group, Element(e).row] for e in els], float)
MARG = RBFInterpolator(gp, np.array([margin(e) for e in els]), kernel="thin_plate_spline", smoothing=0.05)
COND = RBFInterpolator(gp, np.array([cond(e) for e in els]), kernel="thin_plate_spline", smoothing=0.05)
print("real-element terms: el margin cond D surf")
for e in ["Ru", "Ir", "Pt", "Pd", "Rh", "Os", "Mn", "Co", "Ni", "Pb", "Sn", "Ti", "Ge", "Mo", "W", "Ag"]:
    r = out.get(f"real_{e}") or out.get(f"{Element(e).group:.1f}_{Element(e).row:.1f}")
    print(f"{e:3s} {margin(e):+.2f} {cond(e):+.0f}  D {r['D']:+.2f} surf {r['surf']:+.2f}" if r else f"{e:3s} {margin(e):+.2f} {cond(e):+.0f}")
COST = {"Ir": (9, 6), "Ru": (8, 5), "Pt": (10, 6)}
apex = out["real_Ru"]["D"] + (1.6 - 1.28)
def J(g, p, wD=1.0, lam=1.5, wc=3.0, ws=0.3, mu=2.0):
    r = out[f"{g:.1f}_{p:.1f}"]; pen = sum(np.exp(-((g-cg)**2 + (p-cp)**2) / 0.3) for cg, cp in COST.values())
    return -wD * abs(r["D"] - apex) + lam * float(MARG([[g, p]])[0]) + wc * float(COND([[g, p]])[0]) + ws * r["surf"] - mu * pen
G = np.arange(3, 14.01, 0.5); P = np.arange(4, 6.01, 0.5)
# nearest element by chart coords
tbl = {(Element(e).group, Element(e).row): e for e in [str(x) for x in ["Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","La","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb"]]}
def near(g, p): return min(tbl, key=lambda k: (k[0]-g)**2 + (k[1]-p)**2)
for wD in (1.0, 0.3):
    grid = {(float(g), float(p)): J(g, p, wD=wD) for g in G for p in P}
    top = sorted(grid.items(), key=lambda kv: -kv[1])[:8]
    print(f"\nwD={wD} TOP:", [(k, round(v, 2), tbl[near(*k)]) for k, v in top])
    for start in [(8.0, 5.0), (9.0, 6.0), (7.0, 4.0)]:
        g, p = start; path = [tbl[near(g, p)]]
        for _ in range(15):
            cands = [(g+dg, p+dp) for dg in (-0.5, 0, 0.5) for dp in (-0.5, 0, 0.5) if 3 <= g+dg <= 14 and 4 <= p+dp <= 6]
            g2, p2 = max(cands, key=lambda c: J(*c, wD=wD))
            if (g2, p2) == (g, p): break
            g, p = g2, p2; path.append((g, p, tbl[near(g, p)]))
        print(f"  walk from {start}:", path)
# dump objective grid + component terms + walk for the animation
dump = dict(J={f"{g:.1f},{p:.1f}": float(J(g, p, wD=1.0)) for g in G for p in P},
            margin={f"{g:.1f},{p:.1f}": float(MARG([[g, p]])[0]) for g in G for p in P},
            cond={f"{g:.1f},{p:.1f}": float(COND([[g, p]])[0]) for g in G for p in P},
            D={f"{g:.1f},{p:.1f}": out[f"{g:.1f}_{p:.1f}"]["D"] for g in G for p in P}, apex=apex,
            real={e: dict(margin=margin(e), cond=cond(e)) for e in els})
g, p = 8.0, 5.0; path = [(g, p)]
for _ in range(15):
    cands = [(g+dg, p+dp) for dg in (-0.5, 0, 0.5) for dp in (-0.5, 0, 0.5) if 3 <= g+dg <= 14 and 4 <= p+dp <= 6]
    g2, p2 = max(cands, key=lambda c: J(*c)); 
    if (g2, p2) == (g, p): break
    g, p = g2, p2; path.append((g, p))
dump["walk"] = path; json.dump(dump, open(DATA / "twodoor_anim.json", "w"))
print("dumped", path)
