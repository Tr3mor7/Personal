# m3backmap: Martini 3 → all-atom backmapping

`m3backmap` takes a Martini 3 coarse-grained structure (`.gro` or `.pdb`) and
rebuilds an all-atom structure that you can minimise and simulate with an
atomistic force field (set up for CHARMM36 naming by default).

It needs only Python ≥ 3.8 and NumPy.

| CG species | How it is rebuilt | Output |
|---|---|---|
| **Proteins** (martinize2 `BB`, `SC1`…`SC5` beads) | Backbone from a BB-trace fragment library, side chains from ideal internal coordinates with chi angles fitted to the SC beads | Heavy atoms + C-terminal `OXT` (hydrogens are added by `pdb2gmx -ignh`) |
| **Lipids**: POPC, DOPC, DPPC, POPE, DOPE built in | Fragment-wise fit of an all-atom template, one fragment per bead | All atoms, CHARMM36 names and atom order |
| **Any other molecule** | The same template fit, from a `.map` file you make with `make-map` | All atoms |
| **Water** `W` / `SW` / `TW` | 4 / 3 / 2 water molecules per bead | TIP3P (or SPC/E, TIP4P, OPC3 geometry) |
| **Ions** `NA`, `CL`, `K`, `CA`, `MG`, … | 1 ion per bead, plus optional hydration waters | CHARMM (`SOD`, `CLA`, …) or AMBER names |

After placement, a geometric relaxation restores ideal bond lengths, angles,
planarity and chirality while keeping every bead's atoms centred on the CG
bead.

## Quick start

```bash
# from this directory (or: pip install .  → gives the `m3backmap` command)
python -m m3backmap backmap -c cg.gro -o aa.gro --prune-water 0.2
```

This writes:

* `aa.gro`: the full system (protein heavy atoms, everything else all-atom),
* `aa_protein.pdb`: protein only, with chain IDs, ready for `pdb2gmx`,
* `aa_nonprotein.gro`: everything else,
* a `[ molecules ]` block on stderr, in the same order as the output file.

Then finish with GROMACS:

```bash
gmx pdb2gmx -f aa_protein.pdb -o aa_protein_h.gro -p topol.top -ff charmm36-jul2022 -water tip3p -ignh
python -m m3backmap merge -f aa_protein_h.gro aa_nonprotein.gro -o system.gro
# edit topol.top: #include the lipid/ion itps and append the [ molecules ] lines
gmx grompp -f gromacs/em.mdp -c system.gro -p topol.top -o em.tpr
gmx mdrun -deffnm em
# then a short position-restrained run: gromacs/eq_posre.mdp (1 fs, h-bonds constraints)
```

`gromacs/finalize.sh <base> [topol.top]` runs the same steps for you.

### Useful options

| Option | Meaning |
|---|---|
| `-m X.map …`, `-d dir …` | extra templates (files or directories); they override the built-in ones |
| `--itp lipid.itp …` | reorder each residue's atoms to match an itp `[ atoms ]` section (matched by name) |
| `--prune-water R` | delete waters whose O is closer than `R` nm to a solute heavy atom (0.2 is a good start) |
| `--water tip3p\|spce\|tip4p\|opc3`, `--water-resname SOL` | water geometry and residue name |
| `--ions charmm\|amber\|gromos`, `--ion-waters N` | ion naming, and waters added around each ion (Martini 3 ion beads include their hydration shell) |
| `--relax-iter N`, `--no-relax` | relaxation iterations (default 500) or none |
| `--bb-break D` | BB–BB distance (nm) that starts a new chain (default 0.6) |
| `--keep-order` | keep the CG residue order (by default, molecules of one type are grouped together) |
| `--skip-unknown` | drop residues with no rule instead of stopping |
| `--seed` | random seed (orientation of water, ions and 1–2 bead molecules) |

## Adding a molecule (`make-map`)

You need one all-atom structure of the molecule, using the atom names of your
force field (from CHARMM-GUI, ATB, LigParGen or an equilibrated simulation), and
the Martini 3 mapping written as an index-like file:

```text
; beads.ndx  –  bead names must match your Martini itp
[ NH3 ]
N C12
[ PO4 ]
C11 P O11 O12 O13 O14
...
```

```bash
python -m m3backmap make-map -f pope_aa.pdb -b beads.ndx -n POPE -o POPE.map
python -m m3backmap backmap -c cg.gro -m POPE.map -o aa.gro
```

Hydrogens may be left out of `beads.ndx`; each one follows its bonded heavy
atom. Any other atom you leave out goes to the nearest bead. An atom listed
under several beads counts towards each of them, as in Martini mappings that
share atoms.

The `.map` format is plain text: `[ molecule ]` (residue names),
`[ martini ]` (bead order), `[ atoms ]` (`name x y z bead [bead…]`, in nm; this
order is the output order) and an optional `[ bonds ]` section. When
`[ bonds ]` is missing, bonds are inferred from distances.

## How it works

**Protein backbone.** Martini 3 BB beads are the centre of N, CA, C and O. A
library of about 40,000 five-residue windows is built at start-up (about
0.3 s) from ideal-geometry chains that sample the helix, strand, PPII, turn
and left-handed regions of the Ramachandran plot. For each residue, the BB
window i−2…i+2 is compared with the library on distances and signed
dihedrals; the 24 best candidates are then superposed exactly, and the N, CA,
C and O of the best fit are kept. The idea is the same as BBQ, but trained
on BB beads instead of CA atoms. Chain ends use truncated windows.

**Side chains.** CB is placed with ideal L-geometry. The other heavy atoms
come from Engh–Huber internal coordinates. The chi angles come from a rotamer
grid, followed by local refinement, chosen to minimise the distance between
the SC bead positions and the centres of the corresponding atoms, plus a
clash penalty against the nearby backbone.

**Other molecules.** For every bead, the template's bead and its bonded
neighbour beads (≥ 3 points) are superposed onto the CG positions, and the
rotation is applied to the bead's own atoms, centred on the CG bead. Flexible
molecules, such as bent lipid tails, therefore follow the CG conformation.
One- and two-bead molecules are oriented at random about the bead or the
bead axis. Molecules split across periodic boundaries are made whole along
the bead graph.

**Relaxation.** All 1-2 and 1-3 distances are constrained, plus 1-4
distances across planar bonds (peptide bonds, C=C, esters) and rigid
aromatic rings. Signed-volume constraints keep chiral centres (CA, lipid C2,
and others) and keep sp² centres planar. During most of the iterations,
bead-centre restraints keep each fragment on its bead; they are then switched
off, so the final geometry satisfies the constraints.

## Validation

`tools/validate.py` forward-maps a real all-atom structure to Martini 3,
backmaps it, and compares the result with the original heavy atoms. The test
structure was the YiiP transporter dimer (564 residues) with 221 POPE lipids,
CHARMM36, from the MDAnalysis test data. The POPE template was the one built
into this package.

| group | heavy atoms | RMSD | within 1 Å | within 2 Å |
|---|---:|---:|---:|---:|
| protein backbone | 2256 | 0.22 Å | 99.2 % | 99.8 % |
| protein side chains | 2116 | 0.69 Å | 86.5 % | 96.5 % |
| POPE | 10829 | 1.12 Å | 61.0 % | 94.2 % |

No CA or lipid C2 centre had the wrong chirality, and all cis double bonds
stayed cis. The whole run (about 32,000 atoms) takes about 8 s. A system with
about 480,000 atoms, mostly water, takes about 16 s.

Unit tests (round trips with PBC wrapping, geometry, CLI, `make-map`, itp
reordering):

```bash
python -m unittest discover -s tests -v
```

## Limitations and notes

* **This gives a starting structure, not an equilibrated one.** Always run
  energy minimisation and a short position-restrained MD.
* The protein output has heavy atoms only, with PDB names (`ILE CD1`, `OXT`).
  `pdb2gmx -ignh` adds hydrogens and translates names for the force field;
  CHARMM36 and AMBER ports handle these names. Disulfides are detected by
  `pdb2gmx` (`-ss`). Histidine protonation can be chosen with `-his`.
* Unusual residues (ligands, PTMs, non-standard amino acids) need a `.map`
  template. Protein residues are recognised by residue name together with a
  `BB` bead.
* The side-chain bead definitions are in `m3backmap/protein.py`
  (`M3_SIDECHAIN_MAP`) and follow the martinize2 Martini 3 mappings. Edit them
  if your CG model differs.
* The lipid templates use the Martini 3 lipid bead names (`NC3`/`NH3`, `PO4`,
  `GL1`, `GL2`, `C1A`…`C4B`, with `D` for unsaturated beads). If your Martini
  lipid topology uses other names or bead counts, use `make-map` with your own
  bead definitions. `tools/build_pc_templates.py` regenerates the built-in
  templates.
* Side chains with a single SC bead and several chi angles (Gln, Glu, Arg, Lys
  ends) are underdetermined by the CG model, so expect about 1 Å deviations
  there.
* Water is placed on a grid of 4 molecules per bead. Use `--prune-water` to
  remove waters that overlap the solute, and update the `SOL` count printed in
  `[ molecules ]`.

## Layout

```
m3backmap/
  cli.py          command line (backmap, make-map, merge)
  backmap.py      driver: classification, chains, water/ion handling, pruning, itp reordering
  protein.py      backbone library, side-chain builder, Martini 3 protein mapping
  templates.py    .map files, fragment fitting, make-map
  relax.py        constraint-projection relaxation
  solvent.py      water models and ions
  hydrogens.py    geometric H placement (template building)
  structure.py    .gro / .pdb I/O
  geometry.py     superposition, NeRF, PBC helpers
  data/*.map      built-in lipid templates
tools/            template builder, validation script
gromacs/          em.mdp, eq_posre.mdp, finalize.sh
tests/            unittest suite
```

## AI disclosure

This code was written with Claude (Anthropic's AI assistant). The project and its design were directed and reviewed by Gurnoor Hayer.
