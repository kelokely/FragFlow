# FragFlow: depth-resolved cosolvent mapping in lipid bilayers

Data, simulation inputs and figure code for

> Nael, M. A.; Alakonda, L. M.; Elokely, K. M. *FragFlow: depth-resolved cosolvent mapping recovers
> cryptic and allosteric sites inside lipid bilayers.*

Cosolvent molecular dynamics converts probe density to free energy against a reference density. In
water that reference is a single scalar, the bulk probe density. In a bilayer it cannot be: benzene
partitions into the acyl core, so the whole membrane clears a bulk-referenced threshold before a
protein is present. FragFlow replaces the scalar with a depth-resolved reference measured in the
same box and the same frames,

    ddG(r) = -RT ln [ rho_cav(r) / rho_ref(zeta(r)) ],

where zeta = 2*depth/<D_PP> is the reduced membrane depth, zero at the bilayer midplane and +/-1 at
the phosphate planes. This archive holds everything needed to check that claim against the numbers
in the paper, and to run the three systems again from the coordinates they started from.

## Layout

    spec/       ALGORITHMS.md      the five procedures, in full, sufficient to reimplement
                DATA_MANIFEST.md   what every file in data/ is, and which figure or table it feeds
    data/       source data for the figures and tables, plus the cavity inventories,
                per-replicate channel values, residence statistics and lipid-competition events
    scripts/    fig4_reference.py, fig5_specificity.py
    inputs/     aac_1okc_popc, c5ar1_5o9h, trek1_apo: run parameters, topology, index groups,
                starting coordinates and the include-topology files each system needs

## Systems

| directory | protein | organism | PDB | bilayer |
|---|---|---|---|---|
| `aac_1okc_popc` | mitochondrial ADP/ATP carrier | *Bos taurus* | 1OKC | POPC |
| `c5ar1_5o9h` | complement C5a receptor 1 | *Homo sapiens* | 5O9H | POPC |
| `trek1_apo` | TREK-1 | *Mus musculus* | 6CQ8 | POPC |

The panel is not all one organism, and the three architectures are a mitochondrial carrier, a
class-A G-protein-coupled receptor and a two-pore-domain potassium channel. All three are deposited
here in the apo state, each with 91 benzene probes placed in the aqueous slab only.

Production used GROMACS with the CHARMM36m protein and lipid parameters and the TIP3P water model,
at 303.15 K under semi-isotropic pressure coupling. `production.mdp` in each directory carries every
run setting. The starting coordinates in `start.gro` are the equilibrated systems at the point
production began, so a run reproduced from them starts where the reported trajectories start.

Every `#include` in the carrier's and TREK-1's topologies resolves inside its own directory. The
receptor's topology is the one exception: its first line is
`#include "charmm36-jul2021.ff/forcefield.itp"`, which names the published CHARMM36 July 2021
GROMACS force-field distribution rather than a file shipped here. Put that distribution beside
`topol.top`, or point the include at `toppar/forcefield.itp`, which carries the same parameters.

One setting needs saying out loud. `nsteps` is not the production length in every directory, because
the carrier was produced by extension rather than in one shot: its `production.mdp` carries
`continuation = yes` and an `nsteps` of 500 000, which is a single 1 ns chunk. Reported production
per replicate is 200 ns for the carrier, 150 ns for the receptor and 200 ns for TREK-1, and the
first 50 ns of each is discarded before analysis. To reproduce a replicate, set `nsteps` to the
target length or extend the run to it; do not take `nsteps` as given.

## Regenerating the figures

Requires Python with NumPy and Matplotlib.

    cd scripts
    python fig4_reference.py      # Figure 4 and Data File S1
    python fig5_specificity.py    # Figure 5 and Data File S2

Both read only from `data/`. `fig5_specificity.py` verifies its table against `fp2_frozen.json` and
its site depths against `site_depth_frame.json` before it draws, and stops rather than draw a figure
from numbers that disagree with the measurement. Run it without `python -O`: the checks are
assertions, and it refuses to run when they have been stripped.

## Reimplementing the method

`spec/ALGORITHMS.md` gives the depth-resolved reference, the per-frame membrane frame, the
depth-matched control, the density estimator and the pooling rule, each with its frozen constants.
Every quantity reported in the paper follows from those five procedures applied to the trajectories
the inputs here produce. `spec/DATA_MANIFEST.md` closes the loop: for any number in the manuscript
it names the file the number comes from and how to recompute it.

## Reference densities and replicates

Production is eight replicates for AAC and C5aR1 and four for TREK-1 apo. The replicate is the unit
of replication throughout: intervals quoted on a value are computed across replicates, and where a
frame-level interval is also given it is labeled as a statement of precision on the pooled estimate
rather than of uncertainty on the value. The two are not interchangeable, and the manuscript reports
both wherever they differ.

## License

Data, inputs and documentation are released under CC BY 4.0. The scripts in `scripts/` are released
under the MIT License. See `LICENSE`.

The CHARMM36 stream and parameter files under `inputs/*/toppar/` are redistributed under the terms
set by their authors and carry their own headers; they are included so each system is self-contained.
