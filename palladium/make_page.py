"""Self-contained plotly page of the two-door host walk (RuO2 -> Pd/Rh ridge). Output results/p25_manifold/lila/host_walk.html"""
import json, numpy as np, pathlib, matplotlib
from matplotlib import cm
from scipy.interpolate import RBFInterpolator, CubicSpline
from pymatgen.core.periodic_table import Element
from ase import Atoms
HERE = pathlib.Path(__file__).parent; OUT = HERE / "data"
A = json.load(open(OUT / "twodoor_anim.json"))
pts = np.array([[float(a) for a in k.split(",")] for k in A["J"]]); vals = np.array(list(A["J"].values()))
rbf = RBFInterpolator(pts, vals, kernel="thin_plate_spline", smoothing=1e-3)
gs = np.linspace(3, 14, 89); ps = np.linspace(4, 6, 41); G, P = np.meshgrid(gs, ps); Z = rbf(np.c_[G.ravel(), P.ravel()]).reshape(G.shape)
fz = lambda g, p: float(rbf(np.array([[g, p]]))[0])
SHOW = "Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Hf Ta W Re Os Ir Pt Au Hg Tl Pb".split()
beads = [dict(s=s, g=float(Element(s).group), p=float(Element(s).row)) for s in SHOW]
for b in beads: b["z"] = fz(b["g"], b["p"])
W = np.array(A["walk"], float); t = np.linspace(0, 1, len(W)); cs = CubicSpline(t, W, axis=0); tt = np.linspace(0, 1, 90)
path = cs(tt); pz = [fz(g, p) for g, p in path]
turbo = [[i / 15, matplotlib.colors.to_hex(cm.turbo(i / 15))] for i in range(16)]
base = json.load(open(OUT / "ruo2_110_slab.json"))
slab = Atoms(numbers=base["numbers"], positions=base["positions"], cell=base["cell"], pbc=True)
zz = slab.positions[:, 2]; idx = np.where(zz > zz.max() - 3.6)[0]
atoms = [dict(x=float(slab.positions[i, 0]), y=float(slab.positions[i, 1]), z=float(zz[i]), O=bool(slab.numbers[i] == 8)) for i in idx]
ELCOL = {"Ru": "#8A96A3", "Rh": "#5B8FF9", "Pd": "#1E8449", "Ir": "#7B6D8D", "Pt": "#B8B8B8", "Tc": "#3B9AB2", "Os": "#4A6572", "Fe": "#B5651D", "Co": "#3A6EA5", "Ni": "#5E9F5E"}
data = dict(gs=gs.tolist(), ps=ps.tolist(), Z=Z.tolist(), beads=beads, path=path.tolist(), pz=pz, steps=len(W) - 1, turbo=turbo,
            atoms=atoms, cell=slab.cell.lengths().tolist(), elcol=ELCOL)
HTML = r"""<!doctype html><html><head><meta charset="utf-8"><title>Two-door walk</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>body{margin:0;background:#FAF6EE;font-family:Georgia,serif;color:#2B2B2B;overflow:hidden}
#title{text-align:center;font-size:30px;font-weight:bold;padding:16px 0 0}
#wrap{display:flex;height:calc(100vh - 150px)} #surf{flex:2.1} #right{flex:1;display:flex;flex-direction:column;align-items:center}
#mol{width:100%;flex:1} #eltitle{font-size:28px;font-weight:bold;margin-top:6px}
#cap{text-align:center;font-size:23px;padding:6px 0 10px;min-height:34px}
#ctl{text-align:center;font-size:15px;color:#6B6B6B} button{font:15px Georgia;padding:4px 12px;margin:0 4px}</style></head><body>
<div id="title">The right objective walks UMA's periodic table from RuO&#8322; to palladium</div>
<div id="wrap"><div id="surf"></div><div id="right"><div id="eltitle">RuO&#8322; (110) surface</div><div id="mol"></div></div></div>
<div id="cap"></div><div id="ctl"><button onclick="play()">&#9654; play</button><button onclick="reset()">&#8635; reset</button> &nbsp; drag to rotate &middot; scroll to zoom</div>
<script>
const D = DATA_JSON; const INK="#2B2B2B", GREEN="#1E8449", GOLD="#FFD23F";
const surf = {type:"surface", x:D.gs, y:D.ps, z:D.Z, colorscale:D.turbo, showscale:false, opacity:1.0, lighting:{ambient:0.9,diffuse:0.45,specular:0.05,roughness:0.9,fresnel:0.25}, lightposition:{x:-4,y:-6,z:10}, contours:{z:{show:true,usecolormap:true,project:{z:true},width:1}}};
const beads = {type:"scatter3d", mode:"markers+text", x:D.beads.map(b=>b.g), y:D.beads.map(b=>b.p), z:D.beads.map(b=>b.z+0.08), text:D.beads.map(b=>b.s), textposition:"top center", textfont:{size:13,color:INK}, marker:{size:4,color:INK}, hoverinfo:"text"};
const pathTr = {type:"scatter3d", mode:"lines", x:[], y:[], z:[], line:{color:"white",width:12}, hoverinfo:"skip"};
const pathIn = {type:"scatter3d", mode:"lines", x:[], y:[], z:[], line:{color:INK,width:5}, hoverinfo:"skip"};
const head = {type:"scatter3d", mode:"markers", x:[D.path[0][0]], y:[D.path[0][1]], z:[D.pz[0]+0.12], marker:{size:13,color:"white",line:{color:INK,width:3}}, hoverinfo:"skip"};
const stars = {type:"scatter3d", mode:"text", x:[], y:[], z:[], text:[], textfont:{size:44,color:GOLD}, hoverinfo:"skip"};
const lay3 = {paper_bgcolor:"#FAF6EE", margin:{l:0,r:0,t:0,b:0}, scene:{uirevision:"keep", bgcolor:"#FAF6EE", xaxis:{title:"group",color:"#6B6B6B",showbackground:false,autorange:"reversed"}, yaxis:{title:"period",color:"#6B6B6B",showbackground:false,tickvals:[4,5,6]}, zaxis:{title:"anode objective (higher = better)",color:"#6B6B6B",showbackground:false}, aspectratio:{x:2.2,y:1.1,z:0.9}, camera:{eye:{x:-0.9,y:-1.75,z:1.25},up:{x:0,y:0,z:1}}}, showlegend:false};
Plotly.newPlot("surf",[surf,beads,pathTr,pathIn,head,stars],lay3,{displayModeBar:false,responsive:true});
const O = D.atoms.filter(a=>a.O), M = D.atoms.filter(a=>!a.O); const zmin=Math.min(...D.atoms.map(a=>a.z)), zmax=Math.max(...D.atoms.map(a=>a.z));
function shade(hex,z){const d=0.55+0.45*(z-zmin)/(zmax-zmin+1e-6); const c=hex.match(/\w\w/g).map(h=>Math.round(parseInt(h,16)*d)); return "rgb("+c.join(",")+")";}
function molTraces(el,found){ const col = found?GREEN:(D.elcol[el]||"#4C78A8");
  return [{type:"scatter",mode:"markers",x:O.map(a=>a.x),y:O.map(a=>a.y),marker:{size:22,color:O.map(a=>shade("#E2574C",a.z)),line:{color:"#FAF6EE",width:1.5}},hoverinfo:"skip"},
          {type:"scatter",mode:"markers+text",x:M.map(a=>a.x),y:M.map(a=>a.y),text:M.map(_=>el),textfont:{size:13,color:"white",family:"Georgia"},textposition:"middle center",marker:{size:38,color:M.map(a=>shade(col,a.z)),line:{color:INK,width:1.5}},hoverinfo:"skip"}]; }
const layM = {paper_bgcolor:"#FAF6EE",plot_bgcolor:"#FAF6EE",margin:{l:10,r:10,t:10,b:10},xaxis:{visible:false,range:[-1,D.cell[0]+1]},yaxis:{visible:false,range:[-1,D.cell[1]+1.5],scaleanchor:"x"},showlegend:false};
let curEl=null;
function mol(el,found){ if(el===curEl&&!found) return; curEl=el; Plotly.react("mol",molTraces(el,found),layM,{displayModeBar:false,responsive:true});
  const t=document.getElementById("eltitle"); t.textContent = el+"O₂ (110) surface"+(found?"  ✓":""); t.style.color = found?GREEN:INK; }
function nearest(g,p){let b=D.beads[0],d=1e9;for(const q of D.beads){const dd=Math.hypot(q.g-g,q.p-p);if(dd<d){d=dd;b=q}}return b.s}
let k=0, timer=null; const N=D.path.length;
function frame(){ const done=k>=N-1; const g=D.path[k][0],p=D.path[k][1]; const cur=nearest(g,p);
  const xs=D.path.slice(0,k+1).map(q=>q[0]), ys=D.path.slice(0,k+1).map(q=>q[1]), zs=D.pz.slice(0,k+1).map(z=>z+0.1);
  Plotly.restyle("surf",{x:[xs,xs,[g]],y:[ys,ys,[p]],z:[zs,zs,[D.pz[k]+0.12]]},[2,3,4]);
  const step=Math.min(D.steps,Math.floor(k/N*D.steps)+1);
  document.getElementById("cap").textContent = done? "Three steps from RuO₂ to the Pd / Rh ridge.  Lila Sciences, Sept 2026: 2,942 experiments → PdOₓ, 1,000 h" : "Step "+step+": near "+cur;
  mol(done?"Pd":cur,done);
  if(done){const pd=D.beads.find(b=>b.s=="Pd"),rh=D.beads.find(b=>b.s=="Rh"); Plotly.restyle("surf",{x:[[pd.g,rh.g]],y:[[pd.p,rh.p]],z:[[pd.z+0.3,rh.z+0.3]],text:[["★","★"]]},[5]); clearInterval(timer);timer=null;}
  else k++; }
function play(){ if(timer) return; if(k>=N-1) reset(); timer=setInterval(frame,35); }
function reset(){ k=0; if(timer){clearInterval(timer);timer=null} Plotly.restyle("surf",{x:[[],[],[]],y:[[],[],[]],z:[[],[],[]]},[2,3,5]); Plotly.restyle("surf",{x:[[D.path[0][0]]],y:[[D.path[0][1]]],z:[[D.pz[0]+0.12]]},[4]); curEl=null; mol("Ru",false);
  document.getElementById("cap").textContent="Objective: near the activity apex · cation solid at 1.68 V · conducts · not Ir / Ru / Pt"; }
reset();
</script></body></html>"""
open(HERE / "host_walk.html", "w", encoding="utf-8").write(HTML.replace("DATA_JSON", json.dumps(data)))
print("wrote", HERE / "host_walk.html")
