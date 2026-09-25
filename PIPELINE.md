# The walk as a repeatable instrument: what went wrong before the palladium result, and the protocol that follows

The RuO₂→Pd result took one day. Most of that day was spent on things that did not work. Each failure below became a rule.
The protocol at the end is the rule set in order. Point the walk at the next problem by following it top to bottom.

## Failures, in the order they happened

**1. The objective was inherited from the previous problem.** The RuO₂ dopant walk optimised lattice-oxygen binding, so the
host walk started with the same quantity. It ranked Pd below Ru (Pd's lattice O is 1.3 eV weaker) while experiment has PdOₓ
outliving RuOₓ four to one. A day earlier the same descriptor had failed to order eight measured oxide stabilities and we had
not asked why. Rule: derive the objective from how the *incumbent* fails, not from what worked last time. Enumerate every
loss channel (here: anion door and cation door) before computing anything.

**2. The activity descriptor was computed on the wrong surface.** Unrelaxed O* and OH* on a bare cus row put TiO₂ at the
volcano apex and Ir/Ru on the strong-binding leg. At OER potential the cus row is O-covered; on the covered, relaxed surface
the ordering is sane. Rule: compute descriptors in the operating surface state and relax; check that the known champions land
where the field puts them before trusting the rest of the table.

**3. Tabulated inputs had bugs that flattered the wrong elements.** Materials Project lists RuO₄ and OsO₄ as "solids", so the
first window scan gave ruthenium a positive stability margin. Its GGA band gaps call SnO₂ a near-conductor. Its aqueous-ion
data make Ta and Nb dissolve when real Ta₂O₅ is passive. Rule: before a tabulated term enters the objective, check it against
four textbook facts (Ru dissolves, Ir survives, MnO₂'s window closes near 1.7 V, SnO₂ insulates). Exclude volatile MO₄-type
"solids" by O:M ratio.

**4. A physical requirement was weighted instead of imposed.** With conductivity as a soft term, germanium won the objective:
GeO₂ is perfectly Pourbaix-stable and perfectly useless as an anode. Rule: requirements (conducts, stays solid, not the
incumbent) are hard constraints; only genuine trade-offs (activity vs margin) get weights.

**5. The endpoint model was the incumbent's frame.** Rutile PdO₂ is a metastable high-pressure phase; evaluating Pd in the
RuO₂ frame gave a surface dissolution term 2 eV weaker than Ru's, which is a statement about the frame, not about palladium.
Rule: the chart walk may use one frame for navigation, but every endpoint is re-evaluated in its own operando phase
(here Pd(IV)-terminated PdO(101)).

**6. A sweep ran to completion before a three-element sanity check.** The first 29-element host sweep cost ten minutes and
produced a descriptor that was unphysical (failure 2). Rule: run Ru, Ir and one known failure first; only then the table.

**7. A descriptor was degenerate by construction.** In the dopant sweep the "leached" state removes the dopant entirely, so
it is identical for every dopant. That turned out to be informative (any sacrificial dopant leaves the same defect) but it was
not designed; it was noticed after 29 runs. A second bug deleted the dopant site itself in the leaching step. Rule: before a
sweep, write the expected value of each descriptor for two limiting cases by hand and confirm the code reproduces them on one
element.

**8. A retrodiction target that could not falsify anything.** Two hours went into planning a dopant retrodiction against
Lila's In/Mn result before checking the ground truth: the lead is 99.91 at% Pd, lifetimes are n = 1, and their own screen does
not separate the five Pd systems. Rule: audit the ground truth (n, protocol, what actually differs between samples) before
designing the test. If it cannot distinguish hypotheses, do not build a test on it.

**9. A rule was carried across hosts.** High-valence dopants (Ta, W, Nb) stabilise RuO₂. In PdO they destabilise the Pd site,
because Pd is 2+ in that lattice and a 5+ neighbour drains it. Rule: dopant rules are host-specific; recompute, never transfer.

**10. Circularity risk in the retrodiction.** The Ir/Ru/Pt penalty was chosen knowing the answer. It is defensible (the design
goal is a non-incumbent anode) but it had to be argued after the fact. Rule: write the objective, its constraints, the pass
criterion and what is tabulated versus learned *before* running, in a file, and keep it.

**11. Infrastructure that hid failure.** Duplicate OpenMP runtimes killed two GPU jobs silently (fix: `KMP_DUPLICATE_LIB_OK=TRUE`);
piping a long job through `tail` hid all progress until the end; a GIF render competed with the GPU job for CPU and took
twenty minutes; shell heredocs with embedded quotes failed to parse twice. Rule: log to a file with per-line flush, check the
first line of output within a minute, never share the machine between a render and a sweep, and write scripts with a file
tool rather than a shell heredoc.

## The protocol

0. **State the design goal as a sentence** including what is excluded (here: "an acid OER anode that is not Ir, Ru or Pt").
1. **Read how the incumbent fails.** List every loss channel with its governing quantity and whether that quantity is
   computable from energies (UMA) or must be tabulated (redox couples, band gaps, cost). Find the most recent paper that
   decomposes the failure into terms; today that was Gauthier & Maraschin 2026.
2. **Write the objective file before computing.** Terms, which are constraints and which are trade-offs, weights, the operating
   point (pH, U, T), the surface state, the pass criterion, and the known cases it must reproduce.
3. **Validate every tabulated term** on four textbook facts. Exclude artefacts explicitly and record them.
4. **Validate every computed term on three elements** (incumbent, the stable reference, one known failure) in the operating
   surface state, relaxed, before running the table.
5. **Run the landscape** with all sites virtual for navigation; **walk** from the incumbent with greedy or gradient steps;
   report the top of the landscape as well as the walk, and walks from two other starts.
6. **Re-evaluate every endpoint in its own operando phase**, then exact substitution, then DFT anchors on two points.
7. **Write the honest attribution**: which terms decided the landing, which were learned versus tabulated, what the model
   gets wrong in that region (UMA over-binds O on Ru/Ir by ~0.3 eV).
8. **Write the falsification tests** as single experiments before announcing anything. A result without a kill test is a plot.

Applied to the next challenge, steps 0 to 4 are a day of reading and a few hundred forward passes. Step 5 is minutes. Step 6 is
the compute. Everything after is writing.
