#!/usr/bin/env python3
"""Validate backmapping on a real atomistic structure.

The atomistic structure is forward-mapped to Martini 3 beads (protein with the
built-in Martini 3 mapping, other molecules with the .map templates), then
backmapped, and the result is compared with the original heavy atoms.

    python tools/validate.py -f system_aa.gro [-m POPE.map ...] [-o backmapped.gro]
"""

from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from m3backmap import protein as P  # noqa: E402
from m3backmap.backmap import Options, backmap  # noqa: E402
from m3backmap.geometry import dihedral  # noqa: E402
from m3backmap.structure import Atom, Residue, System, read_structure, write_structure  # noqa: E402
from m3backmap.templates import load_maps  # noqa: E402

# atomistic force-field names -> PDB names used by m3backmap
RENAME = {("ILE", "CD"): "CD1", ("*", "OT1"): "O", ("*", "O1"): "O", ("*", "OC1"): "O",
          ("*", "OT2"): "OXT", ("*", "O2"): "OXT", ("*", "OC2"): "OXT"}


def aa_names(res: Residue):
    rn = P.canonical(res.resname)
    out = {}
    for a in res.atoms:
        n = RENAME.get((rn, a.name), RENAME.get(("*", a.name), a.name))
        out[n] = a.xyz
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-f", required=True)
    ap.add_argument("-m", "--maps", nargs="*", default=[])
    ap.add_argument("-o", "--out")
    ap.add_argument("--cg-out", help="write the forward-mapped CG structure")
    ap.add_argument("--no-relax", action="store_true")
    a = ap.parse_args()

    aa = read_structure(a.f)
    maps = load_maps(files=a.maps)
    cg_res, refs = [], []
    for r in aa.residues:
        rn = P.canonical(r.resname)
        if rn and r.get("CA") is not None:
            pos = aa_names(r)
            if not all(n in pos for n in P.atom_order(rn)):
                print(f"skip incomplete {r.resname}{r.resnr}", file=sys.stderr)
                continue
            beads = P.forward_map_residue(rn, pos)
            cg_res.append(Residue(r.resnr, r.resname, [Atom(r.resnr, r.resname, b, x) for b, x in beads.items()]))
            refs.append(("protein", rn, pos))
        elif r.resname in maps:
            t = maps[r.resname]
            pos = {at.name: at.xyz for at in r.atoms}
            if not all(n in pos for n in t.atom_names):
                print(f"skip {r.resname}{r.resnr}: atom names differ from template", file=sys.stderr)
                continue
            xyz = np.array([pos[n] for n in t.atom_names])
            # forward mapping must not be broken by PBC: atoms are assumed whole
            cg_res.append(Residue(r.resnr, r.resname, [Atom(r.resnr, r.resname, b, xyz[m].mean(0))
                                                       for b, m in zip(t.beads, t.members)]))
            refs.append(("map", r.resname, pos))
    cg = System("forward mapped", cg_res, aa.box)
    if a.cg_out:
        write_structure(cg, a.cg_out)
    t0 = time.time()
    out, rep = backmap(cg, maps, Options(verbose=False, relax=not a.no_relax))
    dt = time.time() - t0
    if a.out:
        write_structure(out, a.out)

    box = np.diag(out.box_vectors())
    stats = {}
    chir_bad = 0
    for (kind, rn, pos), r in zip(refs, out.residues):
        for at in r.atoms:
            if at.name.startswith("H") or at.name not in pos:
                continue
            d = at.xyz - pos[at.name]
            if np.all(box > 0):
                d -= np.round(d / box) * box
            key = ("protein backbone" if at.name in P.BB_ATOMS else "protein side chain") if kind == "protein" else rn
            stats.setdefault(key, []).append(np.sum(d * d))
        if kind == "protein" and r.get("CB") is not None:
            imp = np.degrees(dihedral(r.get("C").xyz, r.get("N").xyz, r.get("CA").xyz, r.get("CB").xyz))
            chir_bad += abs(imp + 122.5) > 25
    print(f"backmapped {len(cg_res)} residues / {cg.n_atoms} beads -> {out.n_atoms} atoms in {dt:.1f} s")
    print(f"{'group':<22s}{'atoms':>8s}{'RMSD (A)':>10s}{'<1A':>8s}{'<2A':>8s}")
    for k, v in stats.items():
        v = np.sqrt(np.array(v)) * 10
        print(f"{k:<22s}{len(v):8d}{np.sqrt(np.mean(v ** 2)):10.2f}{np.mean(v < 1) * 100:7.1f}%{np.mean(v < 2) * 100:7.1f}%")
    print(f"CA centres with wrong chirality: {chir_bad}")


if __name__ == "__main__":
    main()
