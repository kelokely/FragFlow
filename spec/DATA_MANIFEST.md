# File manifest

For **FragFlow**, depth-resolved cosolvent mapping in lipid bilayers.

The method is specified in full in **Section S11 of the Supporting Information (Algorithms S1 to S5)**,
in a form sufficient for independent reimplementation. `ALGORITHMS.md` in this archive is a plain-text
copy of that section, provided for convenience.

Labels below match the Supporting Information exactly. Typeset tables in the Supporting Information are
**Table S1** to **Table S8**; the machine-readable sources for the figures and those tables are
**Data File S1** to **Data File S5**; the remaining data files are **File S1** to **File S12**.

The archive holds four directories. `data/` is the manifest below. `spec/` is this file and the
algorithm specification. `scripts/` regenerates Data File S1 and Data File S2, and with them
Figures 4 and 5, from `data/`. `inputs/` is the simulation input set for the three systems: run
parameters, topology, index groups, starting coordinates and the stream files each system needs.

## Source data for the figures and tables

| file | referred to as | contents |
|---|---|---|
| `fig4_reference.csv` | Data File S1 | far-field reference profile per replicate; probe transfer ladder |
| `fig5_specificity.csv` | Data File S2 | firing rate against threshold, both references, all three systems (Table S3) |
| `fig6_c5ar1_channels.csv` | Data File S3 | C5aR1 cavities, four channels, pooled values and intervals (Table S4) |
| `fig11_cav35_channels.csv` | Data File S4 | AAC cavities, same four channels (Table S5) |
| `fig12_residence.csv` | Data File S5 | residence mixture decomposition, per replicate and pooled (Tables S6, S7) |

## Cavity inventories and site characterization

| file | referred to as | contents |
|---|---|---|
| `cavities_c5ar1.json` | File S1 | C5aR1 cavity inventory: volume, occupancy, zeta, extent, centroid, distance to ligand |
| `cavities_aac.json` | File S2 | AAC cavity inventory, same fields |
| `cavities_trek.json` | File S3 | TREK-1 cavity inventory, same fields |
| `cav35_characterization.json` | File S4 | cav35 wall residues, per-replicate occupancy, lining composition |

## Per-replicate channel values

| file | referred to as | contents |
|---|---|---|
| `e30_channels_c5ar1.json` | File S5 | occupancy and hydration per cavity per replicate, n = 8; underlies every replicate-level comparison |

## Residence

| file | referred to as | contents |
|---|---|---|
| `trek_residence.json` | File S6 | TREK-1 re-derived site set and per-site dwell statistics |
| `pi_c_pairs.json` | File S7 | all density-matched site pairs with their dwell ratios, per replicate |

## Detector sensitivity

| file | referred to as | contents |
|---|---|---|
| `engine_sweep.json` | File S8 | site and event counts under the one-at-a-time parameter sweep |

## Lipid treated as a competing ligand

| file | referred to as | contents |
|---|---|---|
| `lipid_occupancy_summary.json` | File S9 | per-site lipid and probe occupancy, dwell, contact character, per-replicate correlation |
| `c5ar1_9p2_lipid_events.csv` | File S10 | every lipid entry and escape event at 9P2: replicate, lipid id, start, end, duration |
| `aac_cav35_lipid_events.csv` | File S11 | the same at cav39/cav35 |
| `site_depth_frame.json` | File S12 | reduced site depth in the membrane frame: per replicate, per site, both poolings, for all three systems |

## A note on two occupancy conventions

Both are over all **eight** production replicates; they differ in how the eight are combined, not in
how many are used. `fig6_c5ar1_channels.csv` and `fig11_cav35_channels.csv` pool frames across the
eight and divide once, so a replicate contributes in proportion to its frame count.
`e30_channels_c5ar1.json` (File S5) reports each replicate separately, and quantities quoted in the
main text as a site occupancy are the mean of those eight per-replicate means, which weights every
replicate equally. At the 9P2 site the two give 3.49 and 3.48 probe per frame respectively. The
difference is small here and need not be, because per-replicate occupancy at this site ranges from
0.12 to 7.57. The file a number comes from is stated wherever it is quoted.

## A note on the depth column

The `zeta` column of `fig6_c5ar1_channels.csv` and `fig11_cav35_channels.csv` is the reduced site
depth in the frozen membrane-frame convention: the circular-mean phosphate midplane taken on raw
coordinates, minimum-image depths, zeta = 2d/thickness, averaged over all eight production
replicates. It is read from a single artifact rather than typed, and that artifact is validated: its
rep1 values reproduce two independently recorded depths to 0.001.

This replaces a column measured along the **aligned** z. Superposing each production
frame onto the analysis reference carries the protein's tilt in the bilayer over to the bilayer, so
that axis is not the membrane normal; the mean tilt is 6.0 deg on C5aR1 and 10.9 deg on AAC. The
same site therefore reads -0.14 on the old axis and -0.07 on this one. The C5aR1 bilayer also
straddles the periodic boundary in every replicate, which an arithmetic midplane would place inside
the acyl core -- the circular mean is what makes the number well defined there at all.

No delta-delta-G in this archive changes with the convention, and none did: the depth-matched
control is placed at the same z as the cavity it is matched to and within the same frame, so every
ratio is taken inside one frame and the depth enters only as a label.

## Reproducing a quoted number

Every quantity in the manuscript is recomputable from these files. For example, the replicate-level
occupancy of the 9P2 site is the mean of `e30_channels_c5ar1.json -> occupancy -> cav3` over its
eight entries; its reduced depth is `2 * (z of the cav3 centroid in cavities_c5ar1.json - z_mid) /
D_PP`, with `z_mid` and `D_PP` as measured in Methods.

## Intermediates used by the figure scripts

| file | contents |
|---|---|
| `fp2_frozen.json` | firing rates of the frozen detection rule, both references, four populations, three systems; the artifact `fig5_specificity.py` asserts against before it draws |
| `fig1c_boxes.json` | mean equilibrated box height and the resulting geometric no-preference baseline for each of the three systems in panel (c) of Figure 4 |
| `perrep/{apo,holo}_rep{1..4}.npz` | far-field probe depths per TREK-1 replicate, the source of panel (a) of Figure 4 |

## Regenerating the figure source data

From `scripts/`, with Python, NumPy and Matplotlib available:

```
python fig4_reference.py      # writes fig4_reference.pdf/.png/.csv   (Data File S1)
python fig5_specificity.py    # writes fig5_specificity.pdf/.png/.csv (Data File S2)
```

Both read only from `data/` and write into the working directory. `fig5_specificity.py` checks its
transcribed table against `fp2_frozen.json` and its site depths against `site_depth_frame.json`
before drawing, and stops rather than draw numbers that disagree with the measurement.
