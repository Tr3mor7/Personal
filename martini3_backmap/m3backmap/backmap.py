"""Backmapping driver: CG system -> all-atom system."""

from __future__ import annotations

import string
import sys
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from . import protein as P
from . import solvent as S
from .geometry import minimum_image, place_atom
from .relax import Relaxer, chiral_quads, graph_constraints
from .structure import Atom, Residue, System
from .templates import MapTemplate

DEG = np.pi / 180.0


@dataclass
class Options:
    water_model: str = "tip3p"
    water_resname: str = "SOL"
    ion_names: str = "charmm"
    ion_waters: int = 0
    relax: bool = True
    relax_iter: int = 500
    max_bb_distance: float = 0.6     # nm; larger BB-BB gap = chain break
    prune_water: float = 0.0         # nm; 0 = keep all waters
    skip_unknown: bool = False
    group_molecules: bool = True     # group residues of the same type (simpler [ molecules ])
    seed: int = 1
    verbose: bool = True


@dataclass
class Report:
    n_cg: int = 0
    n_aa: int = 0
    counts: Dict[str, int] = field(default_factory=dict)
    molecules: List[Tuple[str, int]] = field(default_factory=list)
    chains: List[Tuple[str, int]] = field(default_factory=list)
    skipped: Dict[str, int] = field(default_factory=dict)
    pruned_waters: int = 0
    warnings: List[str] = field(default_factory=list)


def _log(opts: Options, msg: str):
    if opts.verbose:
        print(msg, file=sys.stderr)


# ---------------------------------------------------------------------------
# protein helpers
# ---------------------------------------------------------------------------

def _ideal_peptide_distances() -> Dict[Tuple[str, str], float]:
    X = P.build_backbone(np.array([-60.0, -60.0]), np.array([-45.0, -45.0]))
    N0, CA0, C0, O0 = X[0]
    N1, CA1, C1, O1 = X[1]
    d = lambda a, b: float(np.linalg.norm(a - b))
    # proline CD placed on N1 (cis to CA0 across the peptide bond)
    CD = place_atom(CA0, C0, N1, 0.1473, 128.4 * DEG, 0.0)
    out = {
        ("C", "+N"): d(C0, N1), ("CA", "+N"): d(CA0, N1), ("O", "+N"): d(O0, N1),
        ("C", "+CA"): d(C0, CA1), ("CA", "+CA"): d(CA0, CA1), ("O", "+CA"): d(O0, CA1),
        ("C", "+CD"): d(C0, CD), ("CA", "+CD"): d(CA0, CD), ("O", "+CD"): d(O0, CD),
    }
    return out


_PEPTIDE = None

RING_CLOSURES = {
    "PHE": [("CE2", "CZ")], "TYR": [("CE2", "CZ")],
    "TRP": [("NE1", "CE2"), ("CZ3", "CH2")], "HIS": [("CE1", "NE2")], "PRO": [("CD", "N")],
}
RIGID_GROUPS = {
    "PHE": ["CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "TYR": ["CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"],
    "TRP": ["CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
    "HIS": ["CB", "CG", "ND1", "CD2", "CE1", "NE2"],
}
# ideal proline ring values replacing the (approximately closed) built ring
PRO_IDEAL = {("CD", "N"): 0.1473, ("CA", "CD"): 0.2430, ("CG", "N"): 0.2343}


def _residue_bonds(res: str, names: List[str]) -> List[Tuple[str, str]]:
    b = [("N", "CA"), ("CA", "C"), ("C", "O")]
    if res != "GLY":
        b.append(("CA", "CB"))
    for e in P.SIDECHAINS[res]:
        b.append((e[3], e[0]))
    b += RING_CLOSURES.get(res, [])
    if "OXT" in names:
        b.append(("C", "OXT"))
    return b


def _backmap_chain(residues: List[Residue], box: np.ndarray, rng, opts: Options, relaxer: Relaxer,
                   offset: int, chain_id: str, report: Report) -> List[Residue]:
    global _PEPTIDE
    if _PEPTIDE is None:
        _PEPTIDE = _ideal_peptide_distances()
    L = len(residues)
    bb = np.array([r.get("BB").xyz for r in residues])
    # unwrap along the chain
    for i in range(1, L):
        bb[i] = bb[i - 1] + minimum_image(bb[i] - bb[i - 1], box)
    sc_targets = []
    for r, b in zip(residues, bb):
        t = {}
        for a in r.atoms:
            if a.name != "BB":
                t[a.name] = b + minimum_image(a.xyz - r.get("BB").xyz, box)
        sc_targets.append(t)
    backbone = P.rebuild_backbone(bb, rng)
    out: List[Residue] = []
    all_names, all_xyz, res_of = [], [], []
    for i, r in enumerate(residues):
        rn = P.canonical(r.resname)
        N, CA, C, O = backbone[i]
        lo, hi = max(0, i - 2), min(L, i + 3)
        env = [backbone[j, k] for j in range(lo, hi) for k in range(4) if not (j == i and k < 3)]
        expected = set(P.M3_SIDECHAIN_MAP[rn])
        missing = expected - set(sc_targets[i])
        if missing:
            report.warnings.append(f"{r.resname}{r.resnr}: missing beads {sorted(missing)}")
        sc = P.fit_sidechain(rn, N, CA, C, sc_targets[i], np.array(env))
        pos = {"N": N, "CA": CA, "C": C, "O": O, **sc}
        if i == L - 1:
            pos["OXT"] = P.build_oxt(N, CA, C, O)
        names = [n for n in P.atom_order(rn) if n in pos] + (["OXT"] if "OXT" in pos else [])
        res = Residue(r.resnr, r.resname, [Atom(r.resnr, r.resname, n, pos[n], chain_id) for n in names], chain_id)
        out.append(res)
        for n in names:
            all_names.append(n)
            all_xyz.append(pos[n])
            res_of.append(i)
    if opts.relax:
        _protein_constraints(out, all_names, np.array(all_xyz), res_of, bb, sc_targets, relaxer, offset)
    return out


def _protein_constraints(residues, names, xyz, res_of, bb, sc_targets, relaxer: Relaxer, offset: int):
    index = {}
    for k, (n, i) in enumerate(zip(names, res_of)):
        index[(i, n)] = k
    bonds = []
    planar = [False] * len(names)
    for i, r in enumerate(residues):
        rn = P.canonical(r.resname)
        rnames = [a.name for a in r.atoms]
        for a, b in _residue_bonds(rn, rnames):
            if (i, a) in index and (i, b) in index:
                bonds.append((index[(i, a)], index[(i, b)]))
        if i + 1 < len(residues):
            bonds.append((index[(i, "C")], index[(i + 1, "N")]))
        planar[index[(i, "C")]] = True
        planar[index[(i, "N")]] = True
    I, J, D = graph_constraints(names, xyz, bonds, planar=planar, rigid_rings=False)
    D = D.copy()
    for k, (a, b) in enumerate(zip(I, J)):
        ra, rb = res_of[a], res_of[b]
        na, nb = names[a], names[b]
        if ra != rb:
            if ra > rb:
                na, nb = nb, na
            key = (na, "+" + nb)
            if key in _PEPTIDE:
                D[k] = _PEPTIDE[key]
        elif P.canonical(residues[ra].resname) == "PRO":
            key = tuple(sorted((na, nb)))
            if key in PRO_IDEAL:
                D[k] = PRO_IDEAL[key]
    relaxer.add_constraints(I, J, D, offset)
    relaxer.add_chirality(*chiral_quads(xyz, bonds), offset=offset)
    # rigid aromatic groups
    for i, r in enumerate(residues):
        rn = P.canonical(r.resname)
        grp = [index[(i, n)] for n in RIGID_GROUPS.get(rn, []) if (i, n) in index]
        gi, gj = [], []
        for x in range(len(grp)):
            for y in range(x + 1, len(grp)):
                gi.append(grp[x])
                gj.append(grp[y])
        if gi:
            gi, gj = np.array(gi), np.array(gj)
            relaxer.add_constraints(gi, gj, np.linalg.norm(xyz[gi] - xyz[gj], axis=1), offset)
    # bead restraints
    for i, r in enumerate(residues):
        rn = P.canonical(r.resname)
        relaxer.add_group([index[(i, n)] for n in P.BB_ATOMS], bb[i], offset)
        for bead, atoms in P.M3_SIDECHAIN_MAP[rn].items():
            if bead in sc_targets[i]:
                relaxer.add_group([index[(i, n)] for n in atoms], sc_targets[i][bead], offset)


# ---------------------------------------------------------------------------
# main entry point
# ---------------------------------------------------------------------------

def classify(res: Residue, maps: Dict[str, MapTemplate]) -> str:
    names = res.names
    if P.canonical(res.resname) and "BB" in names:
        return "protein"
    if res.resname in maps:
        return "map"
    if S.is_water(res.resname, names):
        return "water"
    if len(names) == 1 and S.ion_key(res.resname, names[0]):
        return "ion"
    return "unknown"


def backmap(cg: System, maps: Dict[str, MapTemplate], opts: Options = Options()) -> Tuple[System, Report]:
    rng = np.random.default_rng(opts.seed)
    box = cg.box_vectors()
    report = Report(n_cg=cg.n_atoms)
    relaxer = Relaxer()
    out: List[Residue] = []
    n_atoms = 0
    kinds = [classify(r, maps) for r in cg.residues]
    unknown = sorted({r.resname for r, k in zip(cg.residues, kinds) if k == "unknown"})
    if unknown and not opts.skip_unknown:
        raise SystemExit(
            "No backmapping rule for residue(s): " + ", ".join(unknown) +
            "\nProvide a .map template (see `m3backmap make-map -h`) via --maps, "
            "or pass --skip-unknown to drop them.")
    template_constraints: Dict[int, tuple] = {}
    chain_letters = iter(string.ascii_uppercase + string.ascii_lowercase + string.digits)
    resid = 0
    i = 0
    n_res = len(cg.residues)
    while i < n_res:
        res, kind = cg.residues[i], kinds[i]
        if kind == "protein":
            # collect a chain: consecutive protein residues with close BB beads
            chain = [res]
            j = i + 1
            while j < n_res and kinds[j] == "protein":
                d = np.linalg.norm(minimum_image(cg.residues[j].get("BB").xyz - chain[-1].get("BB").xyz, box))
                if d > opts.max_bb_distance:
                    break
                chain.append(cg.residues[j])
                j += 1
            cid = next(chain_letters, "X")
            built = _backmap_chain(chain, box, rng, opts, relaxer, n_atoms, cid, report)
            for r in built:
                n_atoms += len(r.atoms)
                out.append(r)
            report.chains.append((cid, len(chain)))
            report.molecules.append((f"Protein_chain_{cid}", 1))
            report.counts["protein residues"] = report.counts.get("protein residues", 0) + len(chain)
            i = j
            continue
        if kind == "map":
            t = maps[res.resname]
            names = res.names
            missing = [b for b in t.beads if b not in names]
            if missing:
                raise SystemExit(f"{res.resname}{res.resnr}: beads {missing} expected by {t.source} not found")
            raw = np.array([res.get(b).xyz for b in t.beads])
            cgx = t.make_whole(raw, box)
            xyz = t.place(cgx, rng)
            resid += 1
            out.append(Residue(res.resnr, res.resname,
                               [Atom(res.resnr, res.resname, n, x) for n, x in zip(t.atom_names, xyz)]))
            if opts.relax:
                key = id(t)
                if key not in template_constraints:
                    template_constraints[key] = (graph_constraints(t.atom_names, t.xyz, t.bonds),
                                                 chiral_quads(t.xyz, t.bonds))
                cons, quads = template_constraints[key]
                relaxer.add_constraints(*cons, offset=n_atoms)
                relaxer.add_chirality(*quads, offset=n_atoms)
                for k, m in enumerate(t.members):
                    relaxer.add_group(m, cgx[k], n_atoms)
            n_atoms += len(t.atom_names)
            _add_molecule(report, res.resname)
        elif kind == "water":
            for a in res.atoms:
                for mol in S.place_waters(a.xyz, S.WATER_BEADS[a.name.upper()], opts.water_model, rng):
                    out.append(_water_residue(mol, opts))
                    n_atoms += len(mol)
                    _add_molecule(report, opts.water_resname)
        elif kind == "ion":
            a = res.atoms[0]
            key = S.ion_key(res.resname, a.name)
            rname, aname = S.ION_NAMES[opts.ion_names][key]
            out.append(Residue(res.resnr, rname, [Atom(res.resnr, rname, aname, a.xyz.copy())]))
            n_atoms += 1
            _add_molecule(report, rname)
            for mol in S.place_waters(a.xyz, opts.ion_waters, opts.water_model, rng, exclude_center=True) if opts.ion_waters else []:
                out.append(_water_residue(mol, opts))
                n_atoms += len(mol)
                _add_molecule(report, opts.water_resname)
        else:
            report.skipped[res.resname] = report.skipped.get(res.resname, 0) + 1
        i += 1

    if opts.relax and relaxer.ci:
        _log(opts, f"Relaxing geometry ({opts.relax_iter} iterations)...")
        flat = np.array([a.xyz for r in out for a in r.atoms])
        flat = relaxer.run(flat, n_iter=opts.relax_iter, verbose=opts.verbose)
        k = 0
        for r in out:
            for a in r.atoms:
                a.xyz = flat[k]
                k += 1

    if opts.prune_water > 0:
        out, n = prune_waters(out, box, opts.prune_water, opts.water_resname)
        report.pruned_waters = n
        _recount(report, out, opts)

    if opts.group_molecules:
        first: Dict[str, int] = {}
        keys = []
        for r in out:
            k = "\0protein" if r.chain else r.resname
            first.setdefault(k, len(first))
            keys.append(first[k])
        order = sorted(range(len(out)), key=lambda k: keys[k])  # stable
        out = [out[k] for k in order]
        _recount(report, out, opts)

    # sequential residue numbering for non-protein residues
    seq = 0
    for r in out:
        seq += 1
        if not r.chain:
            r.resnr = seq
            for a in r.atoms:
                a.resnr = seq
    aa = System(f"Backmapped from: {cg.title}", out, np.array(cg.box, float))
    report.n_aa = aa.n_atoms
    return aa, report


def _water_residue(mol: np.ndarray, opts: Options) -> Residue:
    names = S.WATER_MODELS[opts.water_model][2]
    return Residue(0, opts.water_resname, [Atom(0, opts.water_resname, n, x) for n, x in zip(names, mol)])


def _add_molecule(report: Report, name: str, n: int = 1):
    if report.molecules and report.molecules[-1][0] == name:
        report.molecules[-1] = (name, report.molecules[-1][1] + n)
    else:
        report.molecules.append((name, n))
    report.counts[name] = report.counts.get(name, 0) + n


def _recount(report: Report, residues: List[Residue], opts: Options):
    mols: List[Tuple[str, int]] = []
    counts: Dict[str, int] = {}
    seen_chains = []
    for r in residues:
        if r.chain:
            if r.chain not in seen_chains:
                seen_chains.append(r.chain)
                mols.append((f"Protein_chain_{r.chain}", 1))
            continue
        if mols and mols[-1][0] == r.resname:
            mols[-1] = (r.resname, mols[-1][1] + 1)
        else:
            mols.append((r.resname, 1))
        counts[r.resname] = counts.get(r.resname, 0) + 1
    if "protein residues" in report.counts:
        counts["protein residues"] = report.counts["protein residues"]
    report.molecules = mols
    report.counts = counts


def prune_waters(residues: List[Residue], box: np.ndarray, cutoff: float, water_resname: str):
    """Remove water molecules whose oxygen is within ``cutoff`` of a non-water
    heavy atom. Periodic for rectangular boxes (triclinic: non-periodic)."""
    solute = np.array([a.xyz for r in residues if r.resname != water_resname
                       for a in r.atoms if not a.name.lstrip("0123456789").startswith("H")]).reshape(-1, 3)
    widx = [k for k, r in enumerate(residues) if r.resname == water_resname]
    if len(solute) == 0 or not widx:
        return residues, 0
    oxy = np.array([residues[k].atoms[0].xyz for k in widx])
    clash = near_any(oxy, solute, cutoff, box)
    drop = {widx[k] for k in np.nonzero(clash)[0]}
    return [r for k, r in enumerate(residues) if k not in drop], len(drop)


def near_any(query: np.ndarray, ref: np.ndarray, cutoff: float, box: np.ndarray) -> np.ndarray:
    """Boolean mask: query points with at least one ref point closer than cutoff."""
    rect = np.allclose(box, np.diag(np.diag(box))) and np.all(np.diag(box) > 0)
    if rect:
        L = np.diag(box)
        nc = np.maximum((L // cutoff).astype(int), 1)
        cs = L / nc
        q, r = query % L, ref % L
    else:
        lo = np.minimum(query.min(0), ref.min(0))
        q, r = query - lo, ref - lo
        cs = np.full(3, cutoff)
        nc = (np.maximum(q.max(0), r.max(0)) // cutoff).astype(int) + 3
    rc = np.minimum((r // cs).astype(int), nc - 1)
    qc = np.minimum((q // cs).astype(int), nc - 1)
    key = (rc[:, 0] * nc[1] + rc[:, 1]) * nc[2] + rc[:, 2]
    order = np.argsort(key)
    skey, sref = key[order], r[order]
    out = np.zeros(len(q), bool)
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for dz in (-1, 0, 1):
                c = qc + (dx, dy, dz)
                if rect:
                    c %= nc
                    ok = np.ones(len(q), bool)
                else:
                    ok = np.all((c >= 0) & (c < nc), axis=1)
                k = (c[:, 0] * nc[1] + c[:, 1]) * nc[2] + c[:, 2]
                start = np.searchsorted(skey, k, "left")
                stop = np.searchsorted(skey, k, "right")
                n = np.where(ok, stop - start, 0)
                for m in range(int(n.max()) if len(n) else 0):
                    sel = np.nonzero((n > m) & ~out)[0]
                    if not len(sel):
                        break
                    d = sref[start[sel] + m] - q[sel]
                    if rect:
                        d -= np.round(d / L) * L
                    out[sel] |= np.einsum("ij,ij->i", d, d) < cutoff * cutoff
    return out


# ---------------------------------------------------------------------------
# reorder residues to match .itp atom order
# ---------------------------------------------------------------------------

def read_itp_orders(paths: Sequence[str]) -> Dict[str, List[str]]:
    """Return {residue name: atom name order} for single-residue moleculetypes."""
    orders: Dict[str, List[str]] = {}
    for path in paths:
        sec = None
        cur: List[Tuple[str, str]] = []

        def flush():
            if cur:
                resn = {r for r, _ in cur}
                if len(resn) == 1:
                    orders[resn.pop()] = [n for _, n in cur]

        with open(path) as fh:
            for raw in fh:
                line = raw.split(";")[0].strip()
                if not line or line.startswith("#"):
                    continue
                if line.startswith("["):
                    sec = line.strip("[] ").lower()
                    if sec == "moleculetype":
                        flush()
                        cur = []
                    continue
                if sec == "atoms":
                    f = line.split()
                    if len(f) >= 5:
                        cur.append((f[3], f[4]))
        flush()
    return orders


def reorder_residues(system: System, orders: Dict[str, List[str]]) -> List[str]:
    msgs = []
    for r in system.residues:
        order = orders.get(r.resname)
        if not order:
            continue
        have = {a.name: a for a in r.atoms}
        if set(have) != set(order):
            missing = set(order) - set(have)
            extra = set(have) - set(order)
            msg = f"{r.resname}: atom names differ from itp (missing {sorted(missing)[:5]}, extra {sorted(extra)[:5]})"
            if msg not in msgs:
                msgs.append(msg)
            continue
        r.atoms = [have[n] for n in order]
    return msgs
