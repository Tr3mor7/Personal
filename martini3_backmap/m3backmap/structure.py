"""Minimal structure I/O for GROMACS .gro and PDB files (coordinates in nm)."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional

import numpy as np


@dataclass
class Atom:
    resnr: int
    resname: str
    name: str
    xyz: np.ndarray  # nm
    chain: str = ""


@dataclass
class Residue:
    resnr: int
    resname: str
    atoms: List[Atom] = field(default_factory=list)
    chain: str = ""

    @property
    def names(self) -> List[str]:
        return [a.name for a in self.atoms]

    def coords(self) -> np.ndarray:
        return np.array([a.xyz for a in self.atoms], dtype=float).reshape(-1, 3)

    def get(self, name: str) -> Optional[Atom]:
        for a in self.atoms:
            if a.name == name:
                return a
        return None


@dataclass
class System:
    title: str
    residues: List[Residue]
    box: np.ndarray  # 3 (rectangular) or 9 (triclinic, gro order) values in nm

    @property
    def n_atoms(self) -> int:
        return sum(len(r.atoms) for r in self.residues)

    def box_vectors(self) -> np.ndarray:
        """Return the 3x3 matrix of box vectors (rows) in nm."""
        b = np.asarray(self.box, dtype=float)
        if b.size == 3:
            return np.diag(b)
        if b.size == 9:
            # gro order: v1(x) v2(y) v3(z) v1(y) v1(z) v2(x) v2(z) v3(x) v3(y)
            v = np.zeros((3, 3))
            v[0, 0], v[1, 1], v[2, 2] = b[0], b[1], b[2]
            v[0, 1], v[0, 2] = b[3], b[4]
            v[1, 0], v[1, 2] = b[5], b[6]
            v[2, 0], v[2, 1] = b[7], b[8]
            return v
        raise ValueError(f"Unsupported box specification: {b}")


# --------------------------------------------------------------------------
# reading
# --------------------------------------------------------------------------

def _group_residues(atoms: List[Atom]) -> List[Residue]:
    residues: List[Residue] = []
    prev = None
    for a in atoms:
        key = (a.resnr, a.resname, a.chain)
        if prev is None or key != prev:
            residues.append(Residue(a.resnr, a.resname, [], a.chain))
            prev = key
        residues[-1].atoms.append(a)
    return residues


def read_gro(path: str) -> System:
    with open(path) as fh:
        lines = fh.read().splitlines()
    title = lines[0].strip()
    n = int(lines[1].strip())
    atoms: List[Atom] = []
    # Detect precision from the distance between decimal points.
    first = lines[2]
    try:
        p1 = first.index(".", 20)
        p2 = first.index(".", p1 + 1)
        width = p2 - p1
    except ValueError:
        width = 8
    for line in lines[2:2 + n]:
        resnr = int(line[0:5])
        resname = line[5:10].strip()
        name = line[10:15].strip()
        x = float(line[20:20 + width])
        y = float(line[20 + width:20 + 2 * width])
        z = float(line[20 + 2 * width:20 + 3 * width])
        atoms.append(Atom(resnr, resname, name, np.array([x, y, z])))
    box = np.array([float(v) for v in lines[2 + n].split()])
    return System(title, _group_residues(atoms), box)


def read_pdb(path: str) -> System:
    atoms: List[Atom] = []
    box = np.zeros(3)
    title = path
    with open(path) as fh:
        for line in fh:
            rec = line[:6]
            if rec.startswith("CRYST1"):
                a, b, c = (float(line[6:15]), float(line[15:24]), float(line[24:33]))
                box = np.array([a, b, c]) / 10.0
            elif rec.startswith("TITLE"):
                title = line[10:].strip()
            elif rec in ("ATOM  ", "HETATM"):
                name = line[12:16].strip()
                resname = line[17:21].strip()
                chain = line[21].strip()
                resnr = int(line[22:26])
                xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]) / 10.0
                atoms.append(Atom(resnr, resname, name, xyz, chain))
            elif rec.startswith("ENDMDL"):
                break
    return System(title, _group_residues(atoms), box)


def read_structure(path: str) -> System:
    if path.lower().endswith(".pdb"):
        return read_pdb(path)
    return read_gro(path)


# --------------------------------------------------------------------------
# writing
# --------------------------------------------------------------------------

def write_gro(system: System, path: str) -> None:
    out = [system.title or "backmapped", f"{system.n_atoms:5d}"]
    idx = 0
    for res in system.residues:
        for a in res.atoms:
            idx += 1
            x, y, z = a.xyz
            out.append(f"{res.resnr % 100000:5d}{res.resname[:5]:<5s}{a.name[:5]:>5s}"
                       f"{idx % 100000:5d}{x:8.3f}{y:8.3f}{z:8.3f}")
    b = np.asarray(system.box, dtype=float)
    out.append(" ".join(f"{v:10.5f}" for v in b))
    with open(path, "w") as fh:
        fh.write("\n".join(out) + "\n")


def _pdb_name(name: str) -> str:
    # PDB convention: names shorter than 4 chars start in column 14 unless the
    # element symbol has two letters; this is good enough for organic molecules.
    if len(name) >= 4:
        return name[:4]
    return " " + name.ljust(3)


def write_pdb(system: System, path: str) -> None:
    out = [f"TITLE     {system.title}"]
    bv = system.box_vectors() * 10.0
    a, b, c = np.linalg.norm(bv, axis=1)
    if a > 0 and b > 0 and c > 0:
        alpha = np.degrees(np.arccos(np.dot(bv[1], bv[2]) / (b * c)))
        beta = np.degrees(np.arccos(np.dot(bv[0], bv[2]) / (a * c)))
        gamma = np.degrees(np.arccos(np.dot(bv[0], bv[1]) / (a * b)))
        out.append(f"CRYST1{a:9.3f}{b:9.3f}{c:9.3f}{alpha:7.2f}{beta:7.2f}{gamma:7.2f} P 1           1")
    idx = 0
    prev_chain = None
    for res in system.residues:
        if prev_chain is not None and res.chain != prev_chain:
            out.append("TER")
        prev_chain = res.chain
        for at in res.atoms:
            idx += 1
            x, y, z = at.xyz * 10.0
            element = at.name.lstrip("0123456789")[:1]
            out.append(f"ATOM  {idx % 100000:5d} {_pdb_name(at.name)} {res.resname[:4]:<4s}"
                       f"{(res.chain or ' ')[:1]}{res.resnr % 10000:4d}    "
                       f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}")
    out.append("TER")
    out.append("END")
    with open(path, "w") as fh:
        fh.write("\n".join(out) + "\n")


def write_structure(system: System, path: str) -> None:
    if path.lower().endswith(".pdb"):
        write_pdb(system, path)
    else:
        write_gro(system, path)
