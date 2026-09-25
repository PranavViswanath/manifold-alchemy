"""Two-door host walk over the (group, period) chart. All cations virtual. Rutile MO2(110), O-terminated cus row, unrelaxed
(fixed RuO2 lattice: navigation only). Per chart point: E_cov, E_OH (H on the site's O*), E_vac (site + O* + 1 lattice O
removed), E_bulk (virtual MO2 rutile 6-atom cell). Real-element sanity at Ru/Ir/Pd/Pt. Then objective + gradient walk."""
import sys, json, time, pathlib, numpy as np, torch
sys.path.insert(0, r"C:\Users\pvisw\Documents\CS\research\interp\sci disc\mlips\code")
from uma_manifold import VirtualUMA, Chart, set_sites, TABLES
from ase import Atoms
from scipy.interpolate import RBFInterpolator
from pymatgen.core.periodic_table import Element
ROOT = pathlib.Path(r"C:\Users\pvisw\Documents\CS\research\interp\sci disc\mlips"); HERE = pathlib.Path(__file__).parent
t0 = time.time(); log = lambda *a: print(f"[{time.time()-t0:6.0f}s]", *a, flush=True)
d = json.load(open(DATA / "ruo2_110_slab.json"))
slab0 = Atoms(numbers=d["numbers"], positions=d["positions"], cell=d["cell"], pbc=True); SITE = d["i"]; UNIT_O = d["unit_O"]
CATS = [int(k) for k in np.where(slab0.numbers == 44)[0]]; CUS = [24,25,26,27,28,29,30,31]
cov = slab0.copy()
for i in CUS: cov += Atoms("O", positions=[slab0.positions[i] + [0, 0, 1.75]])
n0 = len(slab0); ads = {i: n0 + k for k, i in enumerate(CUS)}
oh = cov.copy(); oh += Atoms("H", positions=[cov.positions[ads[SITE]] + [0.62, 0, 0.75]])
rm = sorted([SITE, ads[SITE], UNIT_O[0]], reverse=True); vac = cov.copy(); del vac[rm]
keep = [i for i in range(len(cov)) if i not in rm]; remap = {o: n for n, o in enumerate(keep)}
a, c, u = 4.49, 3.11, 0.306
bulk = Atoms("Ru2O4", scaled_positions=[[0,0,0],[.5,.5,.5],[u,u,0],[1-u,1-u,0],[.5+u,.5-u,.5],[.5-u,.5+u,.5]], cell=[a,a,c], pbc=True)
U = VirtualUMA("cuda"); U.prime("oc22"); W = {k: U.W[k] for k in TABLES}; chart = Chart(W)
def Ev(atoms, cat_idx, vec):
    V = U.default_V(atoms); set_sites(V, cat_idx, vec); return U.energy(atoms, V)
def desc(vec):
    Ecov = Ev(cov, CATS, vec); Eoh = Ev(oh, CATS, vec); Evac = Ev(vac, [remap[i] for i in CATS if i in remap], vec); Eb = Ev(bulk, [0, 1], vec) / 2
    return dict(D=Ecov - Eoh, surf=Evac + Eb - Ecov)   # D lacks +1/2 E_H2 (constant); surf reference-free
# door-2 margin from MP windows: margin(U=1.68) = min(U - U_red, U_ox - U) for the solid oxide window containing/nearest 1.68; negative if outside
win = json.load(open(DATA / "pbx_window.json"))
def margin(el):
    segs = [s for s in win.get(el, []) if s[3] == "Solid" and s[1] != el + "(s)"]
    if not segs: return -1.0
    best = -1.0
    for lo, name, hi, _ in segs:
        hi2 = 2.6 if hi >= 2.58 else hi
        if lo <= 1.68 <= hi2: best = max(best, min(1.68 - lo, hi2 - 1.68))
        else: best = max(best, -min(abs(1.68 - lo), abs(1.68 - hi2)))
    return float(min(best, 1.0))
gp, mg = [], []
for el in win:
    try: e = Element(el); gp.append([e.group, e.row]); mg.append(margin(el))
    except Exception: pass
MARG = RBFInterpolator(np.array(gp, float), np.array(mg), kernel="thin_plate_spline", smoothing=0.05)
COST = {"Ir": (9, 6), "Ru": (8, 5), "Pt": (10, 6)}
out_p = DATA / "twodoor_chart.json"; out = json.load(open(out_p)) if out_p.exists() else {}
for M in ["Ru", "Ir", "Pd", "Pt", "Rh", "Mn"]:
    g, p = chart.coords(M); r = desc(chart.vectors(g, p)); log(f"real {M} ({g},{p}) D {r['D']:+.2f} surf {r['surf']:+.2f} margin {margin(M):+.2f}")
    out[f"real_{M}"] = r
for g in np.arange(3, 14.01, 0.5):
    for p in np.arange(4, 6.01, 0.5):
        k = f"{g:.1f}_{p:.1f}"
        if k in out: continue
        out[k] = desc(chart.vectors(float(g), float(p))); json.dump(out, open(out_p, "w"), indent=1)
        log(f"chart ({g:4.1f},{p:3.1f}) D {out[k]['D']:+.2f} surf {out[k]['surf']:+.2f} margin {float(MARG([[g,p]])[0]):+.2f} nearest {chart.nearest(g,p)}")
json.dump(out, open(out_p, "w"), indent=1); log("landscape done")
# objective: D calibrated so that real Ir/Ru sit where DFT puts them is out of scope; use apex = D(Ru) + (1.6 - 1.28) from the relaxed covered run
D_ru = out["real_Ru"]["D"]; apex = D_ru + (1.6 - 1.28)
def J(g, p, lam=1.5, mu=2.0):
    k = f"{round(g*2)/2:.1f}_{round(p*2)/2:.1f}"; r = out[k]
    pen = sum(np.exp(-((g-cg)**2 + (p-cp)**2) / 0.3) for cg, cp in COST.values())
    return -abs(r["D"] - apex) + lam * float(MARG([[g, p]])[0]) + 0.3 * r["surf"] - mu * pen
grid = {(round(g,1), round(p,1)): J(g, p) for g in np.arange(3, 14.01, 0.5) for p in np.arange(4, 6.01, 0.5)}
top = sorted(grid.items(), key=lambda kv: -kv[1])[:8]
log("TOP objective points:", [(k, round(v, 2), chart.nearest(*k)[0]) for k, v in top])
g, p = 8.0, 5.0; path = [(g, p, chart.nearest(g, p)[0])]
for step in range(12):
    cands = [(g+dg, p+dp) for dg in (-0.5, 0, 0.5) for dp in (-0.5, 0, 0.5) if 3 <= g+dg <= 14 and 4 <= p+dp <= 6]
    g2, p2 = max(cands, key=lambda c: J(*c))
    if (g2, p2) == (g, p): break
    g, p = g2, p2; path.append((g, p, chart.nearest(g, p)[0]))
log("WALK from Ru:", path)
json.dump(dict(top=[(k, v, chart.nearest(*k)[0]) for k, v in top], walk=path, apex=apex), open(DATA / "twodoor_walk.json", "w"), indent=1)
