#!/usr/bin/env python3
"""Build all-atom phosphatidylcholine/-ethanolamine templates (CHARMM36 atom names) and the
corresponding Martini 3 mapping, written as m3backmap .map files.

The geometry is generated from ideal internal coordinates (all-trans tails,
cis double bonds) and hydrogens are added geometrically. The global
conformation does not matter much because backmapping fits each bead's
fragment independently and the relaxation/EM fixes the junctions.

    python tools/build_pc_templates.py            # writes m3backmap/data/*.map
"""

from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from m3backmap.geometry import place_atom  # noqa: E402
from m3backmap.hydrogens import add_hydrogens  # noqa: E402
from m3backmap.templates import MapTemplate  # noqa: E402

DEG = np.pi / 180


def zbuild(zmat):
    """zmat rows: (name, a, b, c, bond[A], angle, torsion); first 3 rows special."""
    pos = {}
    n0 = zmat[0][0]
    pos[n0] = np.zeros(3)
    n1, _, _, c1, b1, _, _ = zmat[1]
    pos[n1] = pos[c1] + np.array([b1 / 10, 0, 0])
    n2, _, b2ref, c2, b2, ang2, _ = zmat[2]
    v = pos[b2ref] - pos[c2]
    v /= np.linalg.norm(v)
    perp = np.array([0.0, 1.0, 0.0])
    a = ang2 * DEG
    pos[n2] = pos[c2] + b2 / 10 * (np.cos(a) * v + np.sin(a) * perp)
    for name, ra, rb, rc, bond, ang, tor in zmat[3:]:
        pos[name] = place_atom(pos[ra], pos[rb], pos[rc], bond / 10, ang * DEG, tor * DEG)
    return pos


def tail(prefix, first, last, ref_a, ref_b, ref_c, double_bonds=()):
    """Acyl chain atoms C{prefix}{first+1..last}; the first one is placed from
    the reference atoms (ref_a, ref_b, ref_c = carbon ``first``)."""
    rows = []
    prev3, prev2, prev1 = ref_a, ref_b, ref_c
    for k in range(first + 1, last + 1):
        name = f"C{prefix}{k}"
        bond, ang, tor = 1.53, 112.0, 180.0
        if (k - 1) in double_bonds:      # this atom is the second atom of a double bond
            bond, ang = 1.34, 123.0
            tor = 180.0
        if k in double_bonds:            # first atom of the double bond (sp2)
            ang = 112.0
        if (k - 2) in double_bonds:      # atom after the double bond: cis
            ang = 123.0
            tor = 0.0
        rows.append((name, prev3, prev2, prev1, bond, ang, tor))
        prev3, prev2, prev1 = prev2, prev1, name
    return rows


def pc_lipid(sn1_len, sn2_len, sn1_db=(), sn2_db=(), head="PC"):
    """Return the z-matrix of a PC (or PE) lipid; tails given as number of acyl
    carbons after the ester carbon."""
    z = [
        ("P", None, None, None, 0, 0, 0),
        ("O11", None, None, "P", 1.60, 0, 0),
        ("O12", None, "O11", "P", 1.60, 104.0, 0),
        ("O13", "O12", "O11", "P", 1.48, 108.0, 120.0),
        ("O14", "O12", "O11", "P", 1.48, 108.0, -120.0),
        # choline
        ("C11", "O11", "P", "O12", 1.43, 120.0, 180.0),
        ("C12", "P", "O12", "C11", 1.52, 110.0, 180.0),
        ("N", "O12", "C11", "C12", 1.50, 115.0, 180.0),
        ("C13", "C11", "C12", "N", 1.49, 109.5, 180.0),
        ("C14", "C11", "C12", "N", 1.49, 109.5, 60.0),
        ("C15", "C11", "C12", "N", 1.49, 109.5, -60.0),
        # glycerol
        ("C1", "O12", "P", "O11", 1.43, 120.0, 180.0),
        ("C2", "P", "O11", "C1", 1.52, 110.0, 180.0),
        ("C3", "O11", "C1", "C2", 1.52, 110.0, 180.0),
        ("O21", "O11", "C1", "C2", 1.45, 109.0, -60.0),
        ("C21", "C1", "C2", "O21", 1.34, 117.0, 150.0),
        ("O22", "C2", "O21", "C21", 1.21, 123.0, 0.0),
        ("C22", "C2", "O21", "C21", 1.51, 111.0, 180.0),
        ("O31", "C1", "C2", "C3", 1.45, 109.0, 180.0),
        ("C31", "C2", "C3", "O31", 1.34, 117.0, 180.0),
        ("O32", "C3", "O31", "C31", 1.21, 123.0, 0.0),
        ("C32", "C3", "O31", "C31", 1.51, 111.0, 180.0),
    ]
    z += tail("2", 2, sn2_len + 1, "O21", "C21", "C22", sn2_db)
    z += tail("3", 2, sn1_len + 1, "O31", "C31", "C32", sn1_db)
    if head == "PE":
        z = [r for r in z if r[0] not in ("C13", "C14", "C15")]
    return z


def bonds_of(zmat):
    b = []
    for name, _, _, c, *_ in zmat:
        if c is not None:
            b.append((c, name))
    return b


def build(name, sn1_len, sn2_len, sn1_db, sn2_db, beads, head="PC"):
    z = pc_lipid(sn1_len, sn2_len, sn1_db, sn2_db, head)
    pos = zbuild(z)
    # enforce R configuration at C2 (sn-glycero-3-phosphocholine)
    c2 = pos["C2"]
    trip = np.dot(pos["O21"] - c2, np.cross(pos["C1"] - c2, pos["C3"] - c2))
    if trip > 0:
        z = [(r[0], r[1], r[2], r[3], r[4], r[5], -r[6]) if r[0] == "O21" else r for r in z]
        pos = zbuild(z)
        c2 = pos["C2"]
        trip = np.dot(pos["O21"] - c2, np.cross(pos["C1"] - c2, pos["C3"] - c2))
    assert trip < 0, "C2 chirality"
    heavy_bonds = bonds_of(z)
    nb = {n: [] for n in pos}
    for a, b in heavy_bonds:
        nb[a].append(b)
        nb[b].append(a)

    # hydrogens (CHARMM36 names)
    H = {}
    owner = {}

    def addh(heavy, names, planar=False):
        h = add_hydrogens(pos, nb, heavy, names, planar=planar)
        H.update(h)
        for n in names:
            owner[n] = heavy

    if head == "PE":
        addh("N", ["HN1", "HN2", "HN3"])
    else:
        for m in ("C13", "C14", "C15"):
            addh(m, [f"H{m[1:]}A", f"H{m[1:]}B", f"H{m[1:]}C"])
    addh("C12", ["H12A", "H12B"])
    addh("C11", ["H11A", "H11B"])
    addh("C1", ["HA", "HB"])
    addh("C2", ["HS"])
    addh("C3", ["HX", "HY"])
    for prefix, length, dbs, (s1, s2, s3) in (("2", sn2_len, sn2_db, ("R", "S", "T")),
                                              ("3", sn1_len, sn1_db, ("X", "Y", "Z"))):
        last = length + 1
        dbatoms = set()
        for d in dbs:
            dbatoms |= {d, d + 1}
        for k in range(2, last + 1):
            c = f"C{prefix}{k}"
            if k in dbatoms:
                hn = f"H{k}1" if prefix == "2" else f"H{k}X"
                addh(c, [hn], planar=True)
            elif k == last:
                addh(c, [f"H{k}{s1}", f"H{k}{s2}", f"H{k}{s3}"])
            else:
                addh(c, [f"H{k}{s1}", f"H{k}{s2}"])
    pos.update(H)

    # CHARMM-like atom order
    if head == "PE":
        order = ["N", "HN1", "HN2", "HN3"]
    else:
        order = ["N", "C13", "H13A", "H13B", "H13C", "C14", "H14A", "H14B", "H14C",
                 "C15", "H15A", "H15B", "H15C"]
    order += ["C12", "H12A", "H12B", "C11", "H11A", "H11B",
              "P", "O13", "O14", "O11", "O12", "C1", "HA", "HB", "C2", "HS", "O21", "C21", "O22"]

    def chain_atoms(prefix, first, last):
        out = []
        for k in range(first, last + 1):
            c = f"C{prefix}{k}"
            out.append(c)
            out += [h for h in H if owner[h] == c]
        return out

    order += chain_atoms("2", 2, 2)
    order += ["C3", "HX", "HY", "O31", "C31", "O32"]
    order += chain_atoms("3", 2, 2)
    order += chain_atoms("2", 3, sn2_len + 1)
    order += chain_atoms("3", 3, sn1_len + 1)
    assert sorted(order) == sorted(pos), set(pos) ^ set(order)

    # bead assignment (heavy atoms; hydrogens follow their carbon)
    bead_of = {}
    for b, atoms in beads.items():
        for a in atoms:
            bead_of.setdefault(a, []).append(b)
    for h, c in owner.items():
        bead_of[h] = [bead_of[c][0]]
    missing = [a for a in order if a not in bead_of]
    assert not missing, missing

    xyz = np.array([pos[a] for a in order])
    xyz -= xyz.mean(0)
    idx = {a: i for i, a in enumerate(order)}
    bonds = [(idx[a], idx[b]) for a, b in heavy_bonds] + [(idx[owner[h]], idx[h]) for h in H]
    bonds = sorted((min(i, j), max(i, j)) for i, j in bonds)
    t = MapTemplate([name], list(beads), order, xyz, [bead_of[a] for a in order], bonds)
    return t


def chain_beads(prefix, letter, length, dbs, unsat_bead=2):
    """Split acyl carbons 2..length+1 into 4 Martini beads."""
    carbons = list(range(2, length + 2))
    n = len(carbons)
    # 4 beads; distribute as evenly as possible, extra carbons in the middle beads
    sizes = [n // 4] * 4
    for k in range(n % 4):
        sizes[[1, 2, 0, 3][k]] += 1
    out, s = {}, 0
    for b in range(4):
        cs = carbons[s:s + sizes[b]]
        s += sizes[b]
        kind = "D" if any(c in dbs or c - 1 in dbs for c in cs) else "C"
        out[f"{kind}{b + 1}{letter}"] = [f"C{prefix}{c}" for c in cs]
    return out


def pc_beads(sn1_len, sn2_len, sn1_db, sn2_db, head="PC"):
    beads = {
        "NC3": ["N", "C12", "C13", "C14", "C15"],
        "PO4": ["C11", "P", "O11", "O12", "O13", "O14"],
        "GL1": ["C2", "O21", "C21", "O22"],
        "GL2": ["C1", "C3", "O31", "C31", "O32"],
    }
    if head == "PE":
        beads = {"NH3": ["N", "C12"], **{k: v for k, v in beads.items() if k != "NC3"}}
    beads.update(chain_beads("2", "A", sn2_len, sn2_db))
    beads.update(chain_beads("3", "B", sn1_len, sn1_db))
    return beads


LIPIDS = {
    # name: (head, sn-1 acyl carbons after C=O, sn-2 ..., sn-1 double bonds, sn-2 double bonds)
    "POPC": ("PC", 15, 17, (), (9,)),
    "DOPC": ("PC", 17, 17, (9,), (9,)),
    "DPPC": ("PC", 15, 15, (), ()),
    "POPE": ("PE", 15, 17, (), (9,)),
    "DOPE": ("PE", 17, 17, (9,), (9,)),
}


def main():
    outdir = os.path.join(os.path.dirname(HERE), "m3backmap", "data")
    os.makedirs(outdir, exist_ok=True)
    for name, (head, l1, l2, d1, d2) in LIPIDS.items():
        beads = pc_beads(l1, l2, d1, d2, head)
        t = build(name, l1, l2, d1, d2, beads, head)
        path = os.path.join(outdir, f"{name}.map")
        t.write(path, comment=(f"{name}: CHARMM36 atom names, Martini 3 beads.\n"
                               "Generated by tools/build_pc_templates.py (ideal geometry)."))
        print(f"{path}: {len(t.atom_names)} atoms, beads {' '.join(t.beads)}")


if __name__ == "__main__":
    main()
