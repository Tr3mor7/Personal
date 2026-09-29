# Tiling a peptide/bilayer system 3×3 in GROMACS

`gmx genconf -nbox 3 3 1` copies the **whole** system nine times, back to back:

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

Tiling only:

```bash
python3 tile_system.py --gro step5_input.gro --top topol.top --ndx index.ndx \
    --nbox 3 3 1 --out-gro tiled.gro --out-top tiled.top --out-ndx tiled.ndx
```

* Write `--out-top` next to the original `topol.top`: `#include "toppar/..."` paths
  are relative to the topology's location.
* Already ran genconf? Pass its output with `--tiled-gro genconf_out.gro`.
* Use `--gmx gmx_mpi` (or similar) if your GROMACS binary has another name.
* `-I DIR` adds search paths for `#include` files. If a moleculetype can't be found
  (such as water inside a force-field directory), the script assumes one residue per
  molecule and prints a note.

The full pipeline (tile → EM → equilibration → production) is in `run_tiled.sh`. Edit
the file names at the top (CHARMM-GUI names are the defaults) and run
`bash run_tiled.sh`.

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
* Custom groups you add later with `gmx make_ndx` should be made from `tiled.gro`.

## Tests

```bash
python3 -m unittest discover -s tests
```

The tests use a fake `gmx genconf` that copies the whole system the same way
GROMACS does. They check that the expanded `[ molecules ]` list matches the output
`.gro` atom by atom in both orderings, that index groups cover all nine copies, and
that a wrong `-nbox` is rejected.
