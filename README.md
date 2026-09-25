# Traveling UMA's manifold

**The project:** Discover new materials with desirable chemical properties by traveling the smooth periodic table manifold inside MLIPs.

![Walking UMA's manifold to a stabler green-hydrogen catalyst](media/manifold_walk.gif)

**Background:** To screen new materials for many applications (drug design, cheaper batteries, carbon capture, etc), we usually simulate chemical reactions with expensive computational methods like DFT, before synthesizing in a lab. However, the search space for possible materials across all elements at many possible sites and configurations is explosive, and most don't have the properties we want.

For example, take the problem of water splitting to produce hydrogen. This process relies on a catalyst (in this example, ruthenium oxide) that holds onto oxygen long enough to separate it from hydrogen in H₂O. However, during the reaction, RuO₂ degrades too quickly because its own oxygen atoms get consumed, and the material disintegrates. The question becomes, can you put (dope) a different atom next to the ruthenium to hold its oxygen atoms in place?

To answer this, traditionally you'd just simulate each of the 118 elements at that position, and rank oxygen removal energy. Try tin? Underbinds. Cobalt? Overpotential. Nickel? Destabilizes. For single sites, this is slow but acceptable computationally. But for most materials problems, with hundreds of sites in different geometries and long horizon reactions, this gets too expensive.

Here's where machine learning interatomic potentials come in. Models like UMA (Universal Model for Atoms) are neural networks trained on hundreds of millions of quantum chemistry calculations to predict the energy of any arrangement of atoms. They do in seconds what traditional methods do in hours.

To work across all of chemistry, UMA has to represent each element as a vector - a point in a high-dimensional space. These vectors weren't designed. They were learned. The model adjusted them until the energies came out right across molecules, crystals, surfaces, catalysts, everything. When you plot where the model placed the elements, you find it independently rediscovered the periodic table! The axes correspond to atomic row, size, and electronegativity. It built its own map of chemistry from nothing but energy.

And unlike Mendeleev's table, this map is continuous. There's a point for every position between real elements, not just at the 90 elements nature gives us. The energy varies smoothly across this surface. Which means it has a slope. Which means you can do calculus on it.

Instead of testing every element one by one, you can place a virtual atom at a site, ask "which direction on this surface improves my objective," and follow the gradient. The walk moves continuously through the model's learned element space, stepping in whichever direction makes the oxygen hold on tighter, until the slope flattens and you've reached the peak. Snap to the nearest real element. That's our answer.

I tested this on the RuO₂ water-splitting problem. The walk started at ruthenium, followed the slope, and arrived at tungsten and tantalum - the same dopants that experimentalists identified over a decade of synthesis and testing. The walk found them in 40 energy evaluations. It was never told any chemistry.

The curvature of the surface (the second derivative/Hessian) turns out to encode even more. From one structure and a few hundred evaluations, it predicts the full ranking of every single-element substitution across 445 validated sites and every co-doping pair. The gradient pointed toward fluorine as the best anion substitution - stepping off the oxide island entirely - which matches the experimental finding that fluorine-doped RuO₂ lasts over 1400 hours.

Three correct predictions, from traveling a surface that a neural network learned from energies and nobody told it was there. Our job: find where the manifold takes us for unsolved materials challenges!

## Reproduce it

**Colab (no setup):** [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/PranavViswanath/manifold-alchemy/blob/main/walkthrough_colab.ipynb) — `walkthrough_colab.ipynb` is fully standalone: it installs fairchem, asks for your Hugging Face token, fetches the two data files from this repo, and runs on a T4.

**Locally:**
```
pip install -r requirements.txt        # needs a Hugging Face token with access to facebook/UMA
jupyter notebook walkthrough.ipynb     # 43 short cells, ~10 min on a laptop GPU with the cached landscape
```

| file | what it is |
|---|---|
| `walkthrough_colab.ipynb` | the same walkthrough, standalone for Colab (library inlined, data fetched from GitHub) |
| `walkthrough.ipynb` | load UMA → open its element tables → plot its periodic table → virtual atoms → the continuous chart → the RuO₂ problem → the objective → the gradient → the landscape → the walk → confirm against exact substitution and the literature |
| `uma_manifold.py` | the whole method in ~150 lines: `VirtualUMA` (any atom carries any element vector), `Chart` (the continuous periodic table), `walk` (gradient ascent on it) |
| `data/ruo2_110_slab.json` | the RuO₂(110) surface, with the active Ru site and dopant positions |
| `data/landscape_RuO2.json` | the cached objective over the continuous table (set `RECOMPUTE = True` to redo it, ~20 min) |
| `media/` | the animation and its final frame |
