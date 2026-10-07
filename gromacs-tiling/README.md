# Tiling a peptide/bilayer system 3×3 in GROMACS

All GROMACS commands here run through **`gmx_mpi`**.

`gmx_mpi genconf -nbox 3 3 1` copies the **whole** system nine times, back to back:

```
[peptide, POPC…, ions, water]  [peptide, POPC…, ions, water]  … ×9
```

Two things break after that:

* **Topology.** Multiplying each `[ molecules ]` count by 9 (`PROA 9`, `POPC 1152`, …)
  gives the order *all peptides, then all lipids, …*. That doesn't match the file,
  so grompp fails with "atom names don't match", or it quietly gives atoms the wrong
  parameters.
* **Index.** The original `index.ndx` only covers the first copy's atom numbers.

`tile_system.py` fixes both. It runs genconf and checks that the result is nine
exact copies. Then it either reorders the atoms so each molecule type sits in one
contiguous block (`--order grouped`, the default, so counts ×9 are correct) or keeps
genconf's order and repeats the `[ molecules ]` list nine times (`--order copy`).
Last, it writes an index file where every original group (e.g. `SOLU`, `MEMB`, `SOLV`)
covers all nine copies, so your `.mdp` `tc-grps`/`comm-grps` work unchanged.
It only needs the Python standard library.

## Usage

```bash
python3 tile_system.py
```

It asks for each file path in turn. You can press Tab to complete paths:

* **Inputs:** the original coordinates, `.top` and `.ndx` (optional, leave blank to
  skip), and an existing `genconf` output. Coordinates can be `.gro` or PDB (e.g.
  CHARMM-GUI's `step5_input.pdb`). PDB content is spotted whatever the file extension,
  converted with `gmx_mpi editconf`, and the tiled output is always written as `.gro`. Leave that blank and it runs
  `gmx_mpi genconf -nbox 3 3 1` itself. Paths that don't exist are asked for again.
* **Outputs:** the tiled `.gro`, `.top` and `.ndx`. Defaults are shown in brackets.
  If a file already exists, it asks before overwriting.

If you'd rather give paths as flags, it uses them and asks only for the rest
(pass `''` to skip the index or genconf output without being asked):

```bash
python3 tile_system.py --gro step5_input.gro --top topol.top --ndx index.ndx --tiled-gro '' \
    --out-gro tiled.gro --out-top tiled.top --out-ndx tiled.ndx
```

* The output topology can go in any directory. Relative `#include "toppar/..."`
  lines are rewritten so they still point at the original files.
* `--nbox NX NY NZ` changes the tiling (default `3 3 1`; don't tile a bilayer in z).
  `--order copy` keeps genconf's atom order instead of grouping molecule types.
* `-I DIR` adds search paths for `#include` files. If a moleculetype can't be found
  (such as water inside a force-field directory), the script assumes one residue per
  molecule and prints a note.

## Using the tiled files

Your `.mdp` files stay as they are: index group names don't change, so `tc-grps`,
`comm-grps` and the rest still match. In your existing `grompp` commands, swap in
the three tiled files: the tiled `.gro` for `-c` (and for `-r` wherever you used the
original coordinates as the restraint reference), the tiled `.top` for `-p`, and the
tiled `.ndx` for `-n`. Make any new index groups from the tiled `.gro`.

## Tests

```bash
python3 -m unittest discover -s tests
```

The tests put a fake `gmx_mpi` on `PATH` whose `genconf` copies the whole system the
same way GROMACS does. They check that the expanded `[ molecules ]` list matches the output `.gro`
atom by atom in both orderings, that index groups cover all nine copies, and that a
wrong `-nbox` is rejected. They also feed answers to the prompts and
check that `genconf` is run as `gmx_mpi`.
