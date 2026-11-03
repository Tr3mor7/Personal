# gromacs-tiling: tile a peptide/bilayer system 3×3 in GROMACS

## The problem

Peptide–membrane simulations usually start from one small patch, for example
**one peptide in a bilayer of 500 lipids (a 1:500 ratio)**. One peptide can't
tell you how peptides behave together: whether they cluster, insert
cooperatively or form pores. You need several peptides at the same dilute
ratio. Tiling the patch 3×3 in the membrane plane gives
**9 peptides among 4,500 lipids**: the ratio is still 1:500, and the
equilibrated local structure of every copy is kept.

GROMACS can copy the coordinates with `gmx_mpi genconf -nbox 3 3 1`. It writes
the **whole** system nine times, back to back:

```
[peptide, lipids…, ions, water]  [peptide, lipids…, ions, water]  … ×9
```

The rest of the setup no longer matches that file:

* **Topology.** Multiplying each `[ molecules ]` count by 9 (`PROA 9`, `POPC 4500`, …)
  describes *all peptides, then all lipids, …*. That's not the order in the file,
  so `grompp` fails with "atom names don't match". Worse, if the names happen to
  line up, atoms quietly get the wrong parameters.
* **Index.** The original `index.ndx` only covers the atoms of the first copy, so
  groups like `SOLU`/`MEMB`/`SOLV` used in `tc-grps` and `comm-grps` are wrong.

`tile_system.py` fixes both, and checks the result, in one step.

## What it does

1. Runs `gmx_mpi genconf -nbox 3 3 1`, or takes a `genconf` output you already have.
2. Checks that the tiled file is exactly nine copies of the original, atom by atom.
3. Writes the coordinates in a consistent order:
   * `--order grouped` (default): each molecule type is gathered into one
     contiguous block, so the topology simply has every count ×9.
   * `--order copy`: keeps genconf's order and repeats the `[ molecules ]` list nine times.
4. Writes a topology whose `[ molecules ]` section matches the coordinates. Relative
   `#include` lines are rewritten so they still find your force-field files.
5. Writes an index file in which every original group covers all nine copies, so
   your `.mdp` files work unchanged.

## Requirements

* **Python 3.6 or newer.** Only the standard library is used; there's nothing to install.
* **GROMACS**, with **`gmx_mpi`** on your `PATH`. It's used for `genconf`, and for
  `editconf` when the input is a PDB file. You don't need GROMACS if you pass an
  existing `genconf` output with `--tiled-gro`.

## Inputs and outputs

| | File | Notes |
|---|---|---|
| **In** | coordinates (`.gro` or `.pdb`) | the untiled, equilibrated system; it must contain the box (for a PDB, a `CRYST1` line) |
| **In** | topology (`.top`) | the matching untiled topology |
| **In** | index (`.ndx`) | optional |
| **In** | `genconf` output (`.gro`) | optional; left blank, the script runs `gmx_mpi genconf` itself |
| **Out** | `<name>_tiled.gro` | tiled coordinates in a consistent molecule order |
| **Out** | `topol_tiled.top` | topology with the matching `[ molecules ]` section |
| **Out** | `index_tiled.ndx` | every original group, covering all nine copies |

All outputs go into an output directory you choose (by default, the current one).

## How to run it

Interactively. It asks for each path in turn, with tab completion and defaults
in brackets, and asks before overwriting anything:

```bash
python3 tile_system.py
```

Or all at once with flags (`''` skips the index or the existing `genconf` output):

```bash
python3 tile_system.py --gro system.gro --top topol.top --ndx index.ndx --tiled-gro '' \
    --out-dir Tiled_3x3 --out-gro tiled.gro --out-top tiled.top --out-ndx tiled.ndx
```

Other options:

* `--nbox NX NY NZ` sets the tiling (default `3 3 1`; don't tile a bilayer in z).
* `--order copy` keeps genconf's atom order instead of grouping molecule types.
* `-I DIR` adds a search path for `#include` files. If a molecule type's definition
  can't be found (for example, water inside a force-field directory), the script
  assumes one residue per molecule and prints a note.

Then use the three tiled files in your usual `grompp` commands: the tiled `.gro`
for `-c` (and for `-r` wherever you used the original coordinates as the
restraint reference), the tiled `.top` for `-p`, and the tiled `.ndx` for `-n`.
Your `.mdp` files stay as they are, because the index group names don't change.

## Tests

```bash
python3 -m unittest discover -s tests
```

The tests need neither GROMACS nor real simulation data. They build a small
synthetic peptide/bilayer/water system and put a fake `gmx_mpi` on `PATH` whose
`genconf` copies the system the same way GROMACS does. They check that:

* the `[ molecules ]` list matches the output `.gro` atom by atom, in both orderings;
* the index groups cover all nine copies;
* a wrong `-nbox` is rejected;
* the interactive prompts work, and `genconf` is run as `gmx_mpi`.

## AI disclosure

This code was written with Claude (Anthropic's AI assistant). The project and its
design were directed and reviewed by Gurnoor Hayer.
