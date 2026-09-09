# Algorithms S1 to S5: numerical specification of FragFlow

This file is a plain-text copy of **Section S11 of the Supporting Information**, provided for convenience;
the typeset version in the Supporting Information is authoritative.
FragFlow is defined by this specification, which gives the five procedures that are not fully determined
by the Methods text, in a form sufficient for independent reimplementation. Every constant is the
frozen value used for all results in the paper.

Notation: `rho` is a probe number density, `zeta` the reduced membrane depth, `RT = 0.602` kcal/mol
at 303.15 K. Coordinates are in Angstroem. Vectors are 3-component.

---

## Algorithm S1. The depth-resolved reference, rho_ref(zeta)

The quantity that replaces the scalar bulk denominator. It is measured in the same box, in the same
frames, as the numerator it divides.

```
INPUT   trajectory frames F, probe set P, lipid phosphorus set Q, protein heavy atoms A
CONST   REF_CUT      = 25.0     # A, minimum lateral distance from protein for a "far-field" probe
        N_BINS       = 22       # zeta bins spanning [-1.5, +1.5]
        DISCARD_PS   = 50000    # partitioning is still transient before this
        RT           = 0.602    # kcal/mol at 303.15 K

FOR each frame f in F with time(f) >= DISCARD_PS:

    # --- membrane frame for THIS frame, on raw coordinates -------------------------------
    z_q      <- z coordinates of Q in frame f
    z_mid    <- circular mean of z_q about the box z-period          # NOT the arithmetic mean:
                                                                     # the bilayer straddles the
                                                                     # periodic boundary
    upper    <- { q in Q : z_q > z_mid },  lower <- { q in Q : z_q < z_mid }
    D_PP(f)  <- mean(z of upper) - mean(z of lower)

    # --- assign every probe a depth and a far-field flag ---------------------------------
    FOR each probe p in P:
        c_p        <- minimum-image centroid of the heavy atoms of p
        zeta(p,f)  <- 2 * (z of c_p - z_mid) / D_PP(f)
        lateral    <- min over a in A of the xy-plane distance |c_p - a|
        far(p,f)   <- ( lateral >= REF_CUT )

    # --- accumulate the far-field histogram ----------------------------------------------
    FOR each probe p with far(p,f):
        b <- bin index of zeta(p,f) in N_BINS bins over [-1.5, +1.5]
        counts[b] += 1
    FOR each bin b:
        volume[b] += (volume of the slab of bin b lying beyond REF_CUT of the protein, in frame f)

rho_ref[b] <- counts[b] / volume[b]          for every bin b
```

**Two points that are not free choices.** The midplane is a *circular* mean because the bilayer
straddles the periodic boundary; an arithmetic mean places it in the water. And `volume[b]` is the
far-field slab volume in that frame, not a fixed constant, because the box fluctuates under
semi-isotropic coupling.

**Scoring a voxel** then uses the reference at the voxel's own depth:

```
ddG(r) = -RT * ln( rho_cav(r) / rho_ref[ bin of zeta(r) ] )
```

---

## Algorithm S2. Per-voxel density with a Gamma-Jeffreys posterior

A voxel visited once must not be scored as if it were measured. The estimator returns a posterior
mean of `ln rho` and its variance, so a sparsely sampled voxel carries a wide interval rather than a
confident wrong value.

```
INPUT   n_eff   effective probe-observation count in the voxel
        E_eff   effective exposure (voxel volume x sampled time)
CONST   detection threshold tau = -2.0 kcal/mol
        count gate n_obs >= 1

ln_rho   <- digamma(n_eff + 0.5) - ln(E_eff)
Var      <- polygamma(1, n_eff + 0.5)

IF n_obs = 0:
    RETURN "insufficient data"        # NOT a value. The posterior is finite here and describes
                                      # the prior, not the site.
```

**Pooling across replicates** uses DerSimonian-Laird random effects with empirical-Bayes shrinkage,
with the **replicate** as the unit. A frame-level interval on the same quantity is narrower by 7x to
15x and does not describe uncertainty on the value.

---

## Algorithm S3. Cavity detection and the zeta-extent percolation test

```
INPUT   protein heavy atoms A for each frame
CONST   LATTICE     = 1.0    # A
        FREE_R      = 3.0    # A, a lattice point is free if no protein heavy atom is nearer
        ENCLOSE_N   = 12     # lattice steps searched along each axis
        AXES        = 7      # 3 Cartesian + 4 body diagonals
        MIN_AXES    = 5      # a point is a pocket point if enclosed on this many axes
        OCC_MIN     = 0.20   # retain points holding pocket status in >= this fraction of frames
        VOL_MIN     = 60     # A^3

FOR each frame:
    mark lattice points free / occupied on a LATTICE grid
    FOR each free point, FOR each of the AXES directions:
        enclosed <- protein encountered within ENCLOSE_N steps in BOTH directions
    pocket point  <- enclosed on >= MIN_AXES axes

retain points with occupancy >= OCC_MIN ; label connected regions ; discard volume < VOL_MIN

# --- percolation, diagnosed by depth extent and NOT by volume -------------------------------
FOR each labeled region R:
    IF min(zeta over R) < -1 AND max(zeta over R) > +1:
        R has crossed both phosphate planes: it is a channel or a lipid-filled groove.
        Split R by iterative erosion until the pieces separate; relabel.
```

**A volume-only criterion was tested first and rejected on measurement**: it admitted a 2389 A^3
merged region that separated into four pieces at the first erosion step, the merge being a
single-voxel neck.

**Known limitation, stated rather than fixed:** the test is one-sided. A groove leaving the bilayer
through a single leaflet passes it. Tightening the rule to reject exit through either leaflet
changes none of the numbers reported in this work.

---

## Algorithm S4. Residence events, and what counts as one

```
INPUT   probe centroids per frame, aligned to a common reference
CONST   CLUSTER_EPS   = 1.5    # A, DBSCAN neighborhood on the density grid
        CLUSTER_MIN   = 5      # minimum grid points per density cluster
        SHELL_CUTOFF  = 8.0    # A, analysis shell around the protein
        CONTACT       = 4.0    # A, probe centroid to site center
        NO GAP BRIDGING

sites <- DBSCAN(probe centroid density, eps=CLUSTER_EPS, min_samples=CLUSTER_MIN)
         restricted to the SHELL_CUTOFF shell

FOR each (probe p, site s):
    occupied(p,s,t) <- distance( centroid(p,t), center(s) ) <= CONTACT
    events          <- maximal runs of consecutive sampled frames with occupied = true
                       A SINGLE ABSENT FRAME TERMINATES A RUN.
```

**No gap bridging is a measured choice, not a convention.** The smallest inter-event gap observed on
the carrier was 400 ps and never 200 ps, across 934 gaps; bridging at 600 ps moves the fast-component
weight from 86.5 % to 76.4 %, so part of what bridging recovers is boundary rattle rather than
continued residence.

**Which of these constants matters.** Sweeping one at a time: `CONTACT` is **inert** across 3.5 to
4.5 A, returning 152 sites and 1670 events at every value. `CLUSTER_EPS` moves the site count 152 to
20 over a 1 A change, and `SHELL_CUTOFF` moves it 82 to 203. A residence timescale may therefore be
quoted without its parameters; **a count of sites may not**.

---

## Algorithm S5. Lipid treated as a competing ligand

Used for the extra-helical site analysis. The pocket is defined from the protein so that the lipid
cannot define its own site.

```
INPUT   protein wall residues W of the site (recorded, not re-derived per frame)
CONST   PROBE_R = 5.0    # A, the same radius the occupancy channel uses

FOR each frame f:
    center(f)  <- centroid of the heavy atoms of W in frame f     # tracks the pocket as it moves
    present(f) <- { lipid residue L : any heavy atom of L within PROBE_R of center(f) }
    probe(f)   <- count of probe heavy atoms within PROBE_R of center(f)

events <- maximal contiguous runs per lipid residue id, no gap bridging (as Algorithm S4)
```

Reported per site: occupancy fraction, event count, dwell distribution, the number of distinct
lipids visiting, the split of contacting atoms between acyl chain and headgroup, and the
within-replicate correlation between lipid and probe occupancy.
