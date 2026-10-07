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

### Full pipeline (tile → EM → equilibration → production)

```bash
bash run_tiled.sh
```

First the script asks for every path it needs. You can press Tab to complete paths:

* **Inputs:** original `.gro`, `.top`, `.ndx` (optional, leave blank to skip),
  the minimization `.mdp`, the number of equilibration steps and each one's `.mdp`,
  and the production `.mdp`. Paths that don't exist are asked for again.
* **Outputs:** tiled `.gro`, `.top`, `.ndx`, and an output name (no extension) for
  the minimization, each equilibration step and the production run. Defaults are
  shown in brackets. If a file already exists, it asks before overwriting.
* An optional MPI launcher for `mdrun` (e.g. `mpirun -np 4`) and extra `mdrun`
  flags (e.g. `-ntomp 8 -nb gpu`).

Then it runs `gmx_mpi genconf`, `tile_system.py`, and `gmx_mpi grompp`/`gmx_mpi mdrun`
for each step. It stops early if `gmx_mpi` isn't on `PATH`, so load your GROMACS module first.

### Tiling only

```bash
python3 tile_system.py
```

It asks for the original `.gro`, `.top` and `.ndx`, and for an existing genconf output.
Leave that blank and it runs `gmx_mpi genconf -nbox 3 3 1` itself. Then it asks for
the output `.gro`, `.top` and `.ndx`.

If you'd rather give paths as flags, it uses them and asks only for the rest
(pass `''` to skip the index or genconf output without being asked):

```bash
python3 tile_system.py --gro step5_input.gro --top topol.top --ndx index.ndx --tiled-gro '' \
    --out-gro tiled.gro --out-top tiled.top --out-ndx tiled.ndx
```

* The output topology can go in any directory. Relative `#include "toppar/..."`
  lines are rewritten so they still point at the original files.
* `--nbox NX NY NZ` changes the tiling (default `3 3 1`). `--order copy` keeps
  genconf's atom order instead of grouping molecule types.
* `-I DIR` adds search paths for `#include` files. If a moleculetype can't be found
  (such as water inside a force-field directory), the script assumes one residue per
  molecule and prints a note.

## Things to check in the simulation setup

* **Tile only in x and y** (`-nbox 3 3 1`). Tiling in z stacks bilayers.
* **Position restraints:** the `-r` reference must be the *tiled* structure
  (`tiled.gro` for EM, then the minimized structure). The original `-r` file has
  too few atoms.
* **Velocities:** all nine copies start with identical velocities. Set
  `gen_vel = yes` (with a random `gen_seed = -1`) in the first equilibration step so
  the copies don't stay correlated. Don't use `-t` with a checkpoint from the
  untiled system.
* The `.mdp` groups keep the same names, so `tc-grps`, `comm-grps` and
  semi-isotropic pressure coupling need no changes. The box is 3× larger in x/y, so
  expect about 9× the cost per step. Consider rechecking `-ntomp`/`-ntmpi` or the GPU
  settings.
* Custom groups you add later with `gmx_mpi make_ndx` should be made from `tiled.gro`.

## Tests

```bash
python3 -m unittest discover -s tests
```

The tests put a fake `gmx_mpi` on `PATH`. Its `genconf` copies the whole system the
same way GROMACS does, and its `grompp`/`mdrun` only create the files the next step
needs. They check that the expanded `[ molecules ]` list matches the output `.gro`
atom by atom in both orderings, that index groups cover all nine copies, and that a
wrong `-nbox` is rejected. They also feed answers to the prompts of both scripts
and check that every GROMACS call is a `gmx_mpi` call with the paths that were entered.
