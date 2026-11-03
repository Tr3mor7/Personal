"""Generic template-based backmapping driven by ``.map`` files.

Map file format (plain text, ``;`` starts a comment)::

    [ molecule ]
    POPC            ; residue name(s) this template applies to
    [ martini ]
    NC3 PO4 GL1 GL2 C1A D2A C3A C4A C1B C2B C3B C4B
    [ atoms ]
    ; name     x       y       z    (nm)  bead(s)
    N        0.102   0.405  1.993   NC3
    ...
    [ bonds ]       ; optional - inferred from distances when absent
    N C12
    ...

The atoms are written in the order given in ``[ atoms ]``. An atom may list
several beads: it then contributes to the centre of every listed bead (as in
Martini mappings that share atoms) and is placed with the first one.

Placement: every bead owns a fragment of atoms. For each bead the template
positions of the bead and its bonded neighbour beads are superposed onto the
CG positions; the rotation is applied to the bead's fragment, which is centred
on the CG bead. Molecules of one or two beads are oriented randomly about the
bead / bead axis.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .geometry import kabsch, random_rotation, rotation_aligning, axis_rotation, normalize

COVALENT_RADII = {"H": 0.031, "C": 0.076, "N": 0.071, "O": 0.066, "P": 0.107,
                  "S": 0.105, "F": 0.057, "B": 0.084, "I": 0.139}
TWO_LETTER = {"CL": 0.102, "BR": 0.120, "SI": 0.111, "SE": 0.120}


def element_radius(name: str) -> float:
    n = name.lstrip("0123456789").upper()
    if n[:2] in TWO_LETTER:
        return TWO_LETTER[n[:2]]
    return COVALENT_RADII.get(n[:1], 0.076)


def is_hydrogen(name: str) -> bool:
    return name.lstrip("0123456789").upper().startswith("H")


def infer_bonds(names: Sequence[str], xyz: np.ndarray, tol: float = 1.25) -> List[Tuple[int, int]]:
    r = np.array([element_radius(n) for n in names])
    h = np.array([is_hydrogen(n) for n in names])
    d = np.linalg.norm(xyz[:, None] - xyz[None], axis=-1)
    cut = tol * (r[:, None] + r[None])
    bonds = []
    n = len(names)
    for i in range(n):
        for j in range(i + 1, n):
            if h[i] and h[j]:
                continue
            if 0.04 < d[i, j] < cut[i, j]:
                bonds.append((i, j))
    # hydrogens keep only their closest partner
    keep = []
    for i in range(n):
        if h[i]:
            cand = [(d[a, b], (a, b)) for (a, b) in bonds if i in (a, b)]
            if cand:
                keep.append(min(cand)[1])
    heavy_bonds = [(a, b) for (a, b) in bonds if not h[a] and not h[b]]
    return sorted(set(heavy_bonds) | set(keep))


@dataclass
class MapTemplate:
    resnames: List[str]
    beads: List[str]
    atom_names: List[str]
    xyz: np.ndarray                       # (n_atoms, 3) nm
    atom_beads: List[List[str]]           # beads per atom (first = owner)
    bonds: List[Tuple[int, int]] = field(default_factory=list)
    source: str = ""

    # derived -------------------------------------------------------------
    def __post_init__(self):
        if not self.bonds:
            self.bonds = infer_bonds(self.atom_names, self.xyz)
        self.bead_index = {b: k for k, b in enumerate(self.beads)}
        for a, bl in zip(self.atom_names, self.atom_beads):
            for b in bl:
                if b not in self.bead_index:
                    raise ValueError(f"{self.source}: atom {a} refers to unknown bead {b}")
        self.owner = np.array([self.bead_index[bl[0]] for bl in self.atom_beads])
        self.members = [[i for i, bl in enumerate(self.atom_beads) if b in bl] for b in self.beads]
        for b, m in zip(self.beads, self.members):
            if not m:
                raise ValueError(f"{self.source}: bead {b} has no atoms")
        self.bead_xyz = np.array([self.xyz[m].mean(0) for m in self.members])
        nb = len(self.beads)
        adj = np.zeros((nb, nb), bool)
        for i, j in self.bonds:
            for bi in self.atom_beads[i]:
                for bj in self.atom_beads[j]:
                    if bi != bj:
                        adj[self.bead_index[bi], self.bead_index[bj]] = True
                        adj[self.bead_index[bj], self.bead_index[bi]] = True
        for bl in self.atom_beads:
            for x in bl:
                for y in bl:
                    if x != y:
                        adj[self.bead_index[x], self.bead_index[y]] = True
        self.bead_adj = adj
        # spanning tree over the bead graph used to make molecules whole
        self.tree: List[Tuple[int, int]] = []
        seen = {0}
        queue = [0]
        while len(seen) < nb:
            if not queue:  # disconnected: attach the nearest unseen bead
                rest = [j for j in range(nb) if j not in seen]
                j = min(rest, key=lambda r: min(np.linalg.norm(self.bead_xyz[r] - self.bead_xyz[s]) for s in seen))
                p = min(seen, key=lambda s: np.linalg.norm(self.bead_xyz[j] - self.bead_xyz[s]))
                self.tree.append((j, p))
                seen.add(j)
                queue.append(j)
                continue
            p = queue.pop(0)
            for j in np.nonzero(adj[p])[0]:
                if int(j) not in seen:
                    seen.add(int(j))
                    queue.append(int(j))
                    self.tree.append((int(j), p))
        self.anchors = [self._anchors(k) for k in range(nb)]

    def _anchors(self, k: int) -> List[int]:
        nb = len(self.beads)
        if nb <= 3:
            return list(range(nb))
        shell = [k]
        seen = {k}
        frontier = [k]
        while len(shell) < 3 and frontier:
            nxt = []
            for f in frontier:
                for j in np.nonzero(self.bead_adj[f])[0]:
                    if j not in seen:
                        seen.add(int(j))
                        nxt.append(int(j))
            # take the whole shell to keep the fit symmetric
            shell.extend(sorted(nxt))
            frontier = nxt
        if len(shell) < 3:  # disconnected bead graph: fall back to nearest beads
            d = np.linalg.norm(self.bead_xyz - self.bead_xyz[k], axis=1)
            for j in np.argsort(d):
                if int(j) not in shell:
                    shell.append(int(j))
                if len(shell) >= 3:
                    break
        return shell

    def make_whole(self, cg: np.ndarray, box: np.ndarray) -> np.ndarray:
        from .geometry import minimum_image
        cg = np.array(cg, float, copy=True)
        if not np.any(box):
            return cg
        for child, parent in self.tree:
            cg[child] = cg[parent] + minimum_image(cg[child] - cg[parent], box)
        return cg

    # placement -----------------------------------------------------------
    def place(self, cg: np.ndarray, rng: np.random.Generator) -> np.ndarray:
        """Return atom coordinates (n_atoms, 3) for CG bead coordinates (n_beads, 3)."""
        nb = len(self.beads)
        out = np.zeros_like(self.xyz)
        if nb == 1:
            R = random_rotation(rng)
            return (self.xyz - self.bead_xyz[0]) @ R.T + cg[0]
        if nb == 2:
            R1 = rotation_aligning(self.bead_xyz[1] - self.bead_xyz[0], cg[1] - cg[0])
            R2 = axis_rotation(cg[1] - cg[0], rng.uniform(0, 2 * np.pi))
            R = R2 @ R1
            for k in range(2):
                m = self.owner == k
                out[m] = (self.xyz[m] - self.bead_xyz[k]) @ R.T + cg[k]
            return out
        for k in range(nb):
            anc = self.anchors[k]
            w = np.where(np.array(anc) == k, 2.0, 1.0)
            R, _ = kabsch(self.bead_xyz[anc], cg[anc], w)
            m = self.owner == k
            out[m] = (self.xyz[m] - self.bead_xyz[k]) @ R.T + cg[k]
        return out

    # io ------------------------------------------------------------------
    def write(self, path: str, comment: str = "") -> None:
        lines = []
        if comment:
            lines += [f"; {c}" for c in comment.splitlines()]
        lines += ["[ molecule ]", " ".join(self.resnames), "", "[ martini ]", " ".join(self.beads), "",
                  "[ atoms ]", ";  name        x         y         z   bead(s)"]
        for n, x, bl in zip(self.atom_names, self.xyz, self.atom_beads):
            lines.append(f"{n:<6s} {x[0]:9.4f} {x[1]:9.4f} {x[2]:9.4f}   {' '.join(bl)}")
        lines += ["", "[ bonds ]"]
        for i, j in self.bonds:
            lines.append(f"{self.atom_names[i]} {self.atom_names[j]}")
        with open(path, "w") as fh:
            fh.write("\n".join(lines) + "\n")


def _sections(path: str) -> Dict[str, List[List[str]]]:
    sec: Dict[str, List[List[str]]] = {}
    cur = None
    with open(path) as fh:
        for raw in fh:
            line = raw.split(";")[0].strip()
            if not line:
                continue
            if line.startswith("[") and line.endswith("]"):
                cur = line.strip("[] ").lower()
                sec.setdefault(cur, [])
                continue
            if cur is None:
                raise ValueError(f"{path}: data outside a section: {line}")
            sec[cur].append(line.split())
    return sec


def read_map(path: str) -> MapTemplate:
    sec = _sections(path)
    for s in ("molecule", "martini", "atoms"):
        if s not in sec:
            raise ValueError(f"{path}: missing [ {s} ] section")
    resnames = [t for row in sec["molecule"] for t in row]
    beads = [t for row in sec["martini"] for t in row]
    names, xyz, ab = [], [], []
    for row in sec["atoms"]:
        if len(row) < 5:
            raise ValueError(f"{path}: bad atom line {' '.join(row)}")
        names.append(row[0])
        xyz.append([float(v) for v in row[1:4]])
        ab.append(row[4:])
    idx = {n: i for i, n in enumerate(names)}
    bonds = []
    for row in sec.get("bonds", []):
        bonds.append((idx[row[0]], idx[row[1]]))
    return MapTemplate(resnames, beads, names, np.array(xyz), ab, bonds, source=path)


def default_map_dirs() -> List[str]:
    return [os.path.join(os.path.dirname(__file__), "data")]


def load_maps(dirs: Sequence[str] = (), files: Sequence[str] = ()) -> Dict[str, MapTemplate]:
    """Load all ``*.map`` files. Later directories/files override earlier ones."""
    maps: Dict[str, MapTemplate] = {}
    paths = []
    for d in list(default_map_dirs()) + list(dirs):
        if os.path.isdir(d):
            paths += [os.path.join(d, f) for f in sorted(os.listdir(d)) if f.endswith(".map")]
    paths += list(files)
    for p in paths:
        t = read_map(p)
        for r in t.resnames:
            maps[r] = t
    return maps


# ---------------------------------------------------------------------------
# make-map: build a template from an atomistic structure + bead definitions
# ---------------------------------------------------------------------------

def read_bead_definitions(path: str) -> Tuple[List[str], Dict[str, List[str]]]:
    """Read an index-like file::

        [ NC3 ]
        N C12 C13 C14 C15
        [ PO4 ]
        ...
    """
    order, defs, cur = [], {}, None
    with open(path) as fh:
        for raw in fh:
            line = raw.split(";")[0].strip()
            if not line:
                continue
            if line.startswith("["):
                cur = line.strip("[] ")
                order.append(cur)
                defs[cur] = []
            else:
                defs[cur] += line.split()
    return order, defs


def make_map(atom_names: Sequence[str], xyz: np.ndarray, bead_order: List[str],
             bead_atoms: Dict[str, List[str]], resnames: List[str],
             hydrogens_follow_parent: bool = True) -> MapTemplate:
    names = list(atom_names)
    beads_of: Dict[str, List[str]] = {n: [] for n in names}
    for b in bead_order:
        for a in bead_atoms[b]:
            if a not in beads_of:
                raise ValueError(f"bead {b}: atom {a} not found in structure")
            beads_of[a].append(b)
    bonds = infer_bonds(names, xyz)
    if hydrogens_follow_parent:
        for i, n in enumerate(names):
            if beads_of[n] or not is_hydrogen(n):
                continue
            partners = [b if a == i else a for a, b in bonds if i in (a, b)]
            if partners and beads_of[names[partners[0]]]:
                beads_of[n] = [beads_of[names[partners[0]]][0]]
    missing = [n for n in names if not beads_of[n]]
    if missing:
        # assign unmapped atoms to the nearest bead centre
        centres = {b: xyz[[names.index(a) for a in bead_atoms[b]]].mean(0) for b in bead_order}
        for n in missing:
            p = xyz[names.index(n)]
            beads_of[n] = [min(bead_order, key=lambda b: np.linalg.norm(centres[b] - p))]
    return MapTemplate(resnames, bead_order, names, np.asarray(xyz, float),
                       [beads_of[n] for n in names], bonds)
