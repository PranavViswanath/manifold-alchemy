"""Minimal tools for walking UMA's element manifold.

UMA (fairchem, uma-s-1p2) turns each element into learned 128-d vectors in four lookup tables:
  sphere_embedding       initial node state           indexed by atom
  source_embedding       one end of every bond message indexed by edge source atom
  target_embedding       other end                     indexed by edge target atom
  composition_embedding  mean-pooled -> expert routing indexed by atom
`VirtualUMA` lets any atom carry ANY 128-d vector in all four tables, so element identity becomes continuous.
Energies are the model's raw per-atom-summed energy (no per-element reference offsets) so virtual atoms are well defined;
every quantity we use is an energy difference, where references cancel anyway.
"""
from __future__ import annotations
import numpy as np, torch, warnings
from scipy.interpolate import RBFInterpolator
from ase import Atoms
from ase.data import chemical_symbols
from pymatgen.core.periodic_table import Element
warnings.filterwarnings("ignore")

NORM = 1.342                      # UMA-S-1.2 energy normaliser: raw head output x NORM = eV
TABLES = ["sph", "src", "tgt", "comp"]


class VirtualUMA:
    """Load uma-s-1p2 and expose element identity as per-atom vectors that can be overridden."""

    def __init__(self, device="cuda"):
        from fairchem.core import pretrained_mlip, FAIRChemCalculator
        from fairchem.core.units.mlip_unit.api.inference import InferenceSettings
        self.Calc = FAIRChemCalculator
        settings = InferenceSettings(tf32=False, activation_checkpointing=False, merge_mole=False, compile=False, execution_mode="general")
        self.pred = pretrained_mlip.get_predict_unit("uma-s-1p2", device=device, inference_settings=settings)
        m = self.pred.model
        while hasattr(m, "module"):
            m = m.module
        self.bb = m.backbone
        self.head = m.output_heads["energyandforcehead"].head.energy_block
        self.mods = {"sph": self.bb.sphere_embedding, "src": self.bb.source_embedding, "tgt": self.bb.target_embedding, "comp": self.bb.composition_embedding}
        self.W = {k: self.mods[k].weight.detach().clone().to(device) for k in TABLES}   # the four element tables [100,128]
        self.dev = device; self.task = None; self.calc = None
        self.S = {"V": None, "ei": None, "E_raw": None, "net": None}
        self._ei_cache = {}
        self._install_hooks()

    # -- hooks: replace the four lookups by per-atom vectors when S["V"] is set -------------------------------------
    def _install_hooks(self):
        S = self.S
        def emb_hook(name):
            def h(mod, inp, out):
                V = S["V"]
                if V is None or V.get(name) is None:
                    return out
                v = V[name].to(out.device, out.dtype)
                if name in ("sph", "comp"):
                    return v                                        # [N,128]
                ei = S["ei"].to(out.device)
                return v[ei[0] if name == "src" else ei[1]]        # [E,128]
            return h
        for k in TABLES:
            self.mods[k].register_forward_hook(emb_hook(k))
        def ei_hook(mod, args, kw, out):                             # capture the graph on real-element passes
            if S["V"] is None:
                S["ei"] = (kw["edge_index"] if "edge_index" in kw else args[2]).detach()
        self.bb.blocks[0].edge_wise.register_forward_hook(ei_hook, with_kwargs=True)
        def energy_hook(mod, inp, out):                              # raw energy from the final scalar channels
            if S["net"] is not None:
                S["E_raw"] = float(S["net"](out["node_embedding"][:, 0, :]).sum())
        self.bb.register_forward_hook(energy_hook)

    def prime(self, task: str):
        """Select the level of theory (omat: bulk PBE; oc22: oxide surfaces; oc20: metal surfaces; omol: molecules)."""
        if self.task == task:
            return
        self.S["net"] = None
        a = Atoms("MgO", positions=[[0, 0, 0], [2.1, 0, 0]], cell=[8, 8, 8], pbc=True); a.info.update(charge=0, spin=1)
        a.calc = self.Calc(self.pred, task_name=task); a.get_potential_energy()
        eb = self.head
        self.S["net"] = torch.nn.Sequential(eb[0].merged_linear_layer(), torch.nn.SiLU(), eb[2].merged_linear_layer(), torch.nn.SiLU(), eb[4].merged_linear_layer())
        self.task = task; self.calc = self.Calc(self.pred, task_name=task)

    def _run(self, atoms):
        a = atoms.copy(); a.info.update(charge=0, spin=1)
        self.calc.reset()                                            # ASE caches identical atoms; embeddings changed
        a.calc = self.calc; a.get_potential_energy()
        return self.S["E_raw"] * NORM

    def default_V(self, atoms):
        """The real element vectors of each atom, as overridable per-atom tensors."""
        Z = torch.tensor(atoms.numbers, device=self.dev)
        return {k: self.W[k][Z].clone() for k in TABLES}

    def energy(self, atoms, V=None):
        """Energy in eV (reference-free). V = per-atom vectors from default_V, possibly edited with set_sites."""
        key = (atoms.numbers.tobytes(), np.round(atoms.positions, 5).tobytes(), np.round(atoms.cell.array, 5).tobytes())
        if V is None:
            self.S["V"] = None; E = self._run(atoms); self._ei_cache[key] = self.S["ei"]; return E
        if key not in self._ei_cache:
            while len(self._ei_cache) >= 64:
                self._ei_cache.pop(next(iter(self._ei_cache)))
            self.S["V"] = None; self._run(atoms); self._ei_cache[key] = self.S["ei"]
        self.S["ei"] = self._ei_cache[key]; self.S["V"] = V
        try:
            return self._run(atoms)
        finally:
            self.S["V"] = None


def set_sites(V, sites, vecs):
    """Give atoms `sites` the element vectors `vecs` (dict table -> 128-vector). Edits V in place."""
    for k in V:
        v = vecs[k]
        v = torch.as_tensor(np.asarray(v) if not torch.is_tensor(v) else v, dtype=V[k].dtype, device=V[k].device)
        V[k][list(sites)] = v
    return V


def delete_rows(V, idx):
    keep = [i for i in range(next(iter(V.values())).shape[0]) if i not in set(idx)]
    return {k: V[k][keep] for k in V}


def element_vectors(W, symbol):
    from ase.data import atomic_numbers
    return {k: W[k][atomic_numbers[symbol]] for k in TABLES}


# ---------------------------------------------------------------------------------------------------------------------
class Chart:
    """The continuous periodic table: a smooth map (group, period) -> element vectors, through the real elements.
    Thin-plate spline per table (Goodfire manifold-steering recipe: known concept coordinates, spline through centroids)."""

    def __init__(self, W, zmax=83):
        self.els = []
        for z in range(3, zmax + 1):
            e = Element(chemical_symbols[z])
            if 58 <= z <= 71 or e.group == 18:          # skip lanthanides and noble gases (little bonding data)
                continue
            self.els.append((z, float(e.group), float(e.row)))
        self.Z = np.array([z for z, _, _ in self.els]); self.gp = np.array([[g, p] for _, g, p in self.els])
        self.W = {k: (W[k].cpu().numpy() if hasattr(W[k], "cpu") else np.asarray(W[k])) for k in W}
        self.rbf = {k: RBFInterpolator(self.gp, self.W[k][self.Z], kernel="thin_plate_spline") for k in self.W}

    def vectors(self, g, p):
        return {k: self.rbf[k](np.array([[g, p]], float))[0] for k in self.W}

    def coords(self, symbol):
        e = Element(symbol); return float(e.group), float(e.row)

    def nearest(self, g, p):
        d = np.hypot(self.gp[:, 0] - g, self.gp[:, 1] - p); i = int(np.argmin(d))
        return chemical_symbols[self.Z[i]], float(d[i])


def walk(objective, start, chart, steps=8, h=0.25, step=0.5, decay=1.0, bounds=((2, 14), (4, 6)), log=print):
    """Gradient ascent of `objective(g, p)` on the chart by central finite differences. Returns the path."""
    g, p = chart.coords(start); path = [(g, p)]
    for it in range(steps):
        f0 = objective(g, p)
        dg = (objective(g + h, p) - objective(g - h, p)) / (2 * h)
        dp = (objective(g, p + h) - objective(g, p - h)) / (2 * h)
        n = np.hypot(dg, dp) + 1e-9; st = step * decay ** it
        g = float(np.clip(g + st * dg / n, *bounds[0])); p = float(np.clip(p + st * dp / n, *bounds[1]))
        path.append((g, p))
        log(f"step {it + 1}: f {f0:+.3f}  ->  (group {g:.2f}, period {p:.2f})  nearest {chart.nearest(g, p)[0]}")
    return path
