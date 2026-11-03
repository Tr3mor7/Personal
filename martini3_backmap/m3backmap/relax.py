"""Geometric relaxation of a backmapped structure.

Fragment-based placement leaves distorted bonds at fragment junctions. This
module iteratively projects distance constraints (1-2, 1-3 and selected 1-4
pairs, i.e. bonds, angles and planar groups) while restraining every bead's
atoms to stay centred on the CG bead. It is not a force field minimisation -
it just provides a much better starting point for one.
"""

from __future__ import annotations

from collections import defaultdict
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np

from .templates import is_hydrogen


def scatter_add(index: np.ndarray, values: np.ndarray, n: int) -> np.ndarray:
    """Fast equivalent of np.add.at(out, index, values) for (m, 3) values."""
    return np.stack([np.bincount(index, weights=values[:, k], minlength=n) for k in range(3)], axis=1)


def _adjacency(n: int, bonds: Iterable[Tuple[int, int]]) -> List[List[int]]:
    adj: List[List[int]] = [[] for _ in range(n)]
    for i, j in bonds:
        adj[i].append(j)
        adj[j].append(i)
    return adj


def small_rings(adj: List[List[int]], max_size: int = 6) -> List[List[int]]:
    rings = set()
    n = len(adj)

    def dfs(start, cur, path):
        if len(path) > max_size:
            return
        for nb in adj[cur]:
            if nb == start and len(path) >= 3:
                rings.add(tuple(sorted(path)))
            elif nb not in path and nb > start:
                dfs(start, nb, path + [nb])

    for s in range(n):
        if len(adj[s]) >= 2:
            dfs(s, s, [s])
    # keep minimal rings only (drop rings that are unions of smaller ones)
    out = sorted(rings, key=len)
    minimal: List[List[int]] = []
    for r in out:
        rs = set(r)
        if not any(set(m) < rs for m in minimal):
            minimal.append(list(r))
    return minimal


def is_planar_atom(name: str, n_neighbors: int, has_h: bool) -> bool:
    e = name.lstrip("0123456789").upper()[:1]
    if e == "C":
        return n_neighbors == 3 if has_h else n_neighbors <= 3
    return False


def graph_constraints(names: Sequence[str], xyz: np.ndarray, bonds: Sequence[Tuple[int, int]],
                      planar: Sequence[bool] | None = None, rigid_rings: bool = True):
    """Return (i, j, d0) arrays of constraints derived from a reference
    geometry ``xyz``: bonds, angles (1-3), torsions about bonds between two
    planar atoms (1-4) and all pairs within small rings."""
    n = len(names)
    adj = _adjacency(n, bonds)
    has_h = any(is_hydrogen(x) for x in names)
    if planar is None:
        planar = [is_planar_atom(names[i], len(adj[i]), has_h) for i in range(n)]
    pairs = set()
    for i, j in bonds:
        pairs.add((min(i, j), max(i, j)))
    for c in range(n):
        for a in adj[c]:
            for b in adj[c]:
                if a < b:
                    pairs.add((a, b))
    for j, k in bonds:
        if planar[j] and planar[k]:
            for i in adj[j]:
                if i == k:
                    continue
                for l in adj[k]:
                    if l != j and l != i:
                        pairs.add((min(i, l), max(i, l)))
    if rigid_rings:
        for ring in small_rings(adj):
            members = set(ring)
            for r in ring:  # include ring substituents
                members.update(adj[r])
            ml = sorted(members)
            for x in range(len(ml)):
                for y in range(x + 1, len(ml)):
                    pairs.add((ml[x], ml[y]))
    if not pairs:
        return np.zeros(0, int), np.zeros(0, int), np.zeros(0)
    p = np.array(sorted(pairs))
    d0 = np.linalg.norm(xyz[p[:, 0]] - xyz[p[:, 1]], axis=1)
    return p[:, 0], p[:, 1], d0


def chiral_quads(xyz: np.ndarray, bonds: Sequence[Tuple[int, int]], planar_cut: float = 0.3):
    """Signed-volume constraints (centre, a, b, c, V0) for every atom with 3 or
    4 neighbours. Nearly planar centres get V0 = 0 (enforces planarity),
    tetrahedral centres keep their reference handedness."""
    adj = _adjacency(len(xyz), bonds)
    out = []
    for c, nb in enumerate(adj):
        if len(nb) not in (3, 4):
            continue
        a, b, d = nb[:3]
        u, v, w = xyz[a] - xyz[c], xyz[b] - xyz[c], xyz[d] - xyz[c]
        V = float(np.dot(u, np.cross(v, w)))
        norm = np.linalg.norm(u) * np.linalg.norm(v) * np.linalg.norm(w)
        if norm < 1e-12:
            continue
        if abs(V) / norm < planar_cut:
            if len(nb) == 3:
                out.append((c, a, b, d, 0.0))
        else:
            out.append((c, a, b, d, V))
    if not out:
        return np.zeros((0, 4), int), np.zeros(0)
    q = np.array([o[:4] for o in out], int)
    return q, np.array([o[4] for o in out])


class Relaxer:
    """Collects constraints/bead groups for a whole system and relaxes it."""

    def __init__(self):
        self.ci: List[np.ndarray] = []
        self.cj: List[np.ndarray] = []
        self.cd: List[np.ndarray] = []
        self.g_atoms: List[np.ndarray] = []
        self.g_target: List[np.ndarray] = []
        self.quads: List[np.ndarray] = []
        self.qv: List[np.ndarray] = []

    def add_chirality(self, quads, v0, offset: int = 0):
        if len(quads):
            self.quads.append(np.asarray(quads, int) + offset)
            self.qv.append(np.asarray(v0, float))

    def add_constraints(self, i, j, d0, offset: int = 0):
        if len(i):
            self.ci.append(np.asarray(i) + offset)
            self.cj.append(np.asarray(j) + offset)
            self.cd.append(np.asarray(d0, float))

    def add_group(self, atoms: Sequence[int], target: np.ndarray, offset: int = 0):
        self.g_atoms.append(np.asarray(atoms, int) + offset)
        self.g_target.append(np.asarray(target, float))

    def run(self, x: np.ndarray, n_iter: int = 200, omega: float = 1.2,
            bead_strength: float = 0.2, verbose: bool = False) -> np.ndarray:
        if not self.ci:
            return x
        full = np.array(x, float, copy=True)
        # work only on atoms that take part in a constraint (e.g. not waters)
        used = np.unique(np.concatenate(self.ci + self.cj + [g for g in self.g_atoms] +
                                        [q.ravel() for q in self.quads]))
        remap = np.full(len(full), -1, int)
        remap[used] = np.arange(len(used))
        x = full[used]
        I = remap[np.concatenate(self.ci)]
        J = remap[np.concatenate(self.cj)]
        D = np.concatenate(self.cd)
        IJ = np.concatenate([I, J])
        cnt = np.bincount(IJ, minlength=len(x)).astype(float)
        cnt[cnt == 0] = 1.0
        if self.g_atoms:
            G = remap[np.concatenate(self.g_atoms)]
            gid = np.concatenate([np.full(len(a), k) for k, a in enumerate(self.g_atoms)])
            gsize = np.array([len(a) for a in self.g_atoms], float)
            T = np.array(self.g_target)
            gcnt = np.bincount(G, minlength=len(x)).astype(float)
            gcnt[gcnt == 0] = 1.0
        if self.quads:
            Q = remap[np.concatenate(self.quads)]
            QV = np.concatenate(self.qv)
            qcnt = np.bincount(Q[:, 0], minlength=len(x)).astype(float)
            qcnt[qcnt == 0] = 1.0
        for it in range(n_iter):
            if self.quads:
                c = x[Q[:, 0]]
                u, v, w = x[Q[:, 1]] - c, x[Q[:, 2]] - c, x[Q[:, 3]] - c
                V = np.einsum("ij,ij->i", u, np.cross(v, w))
                # dV/dc = -(v x w + w x u + u x v)
                g = -(np.cross(v, w) + np.cross(w, u) + np.cross(u, v))
                g2 = np.einsum("ij,ij->i", g, g)
                step = np.where(g2 > 1e-16, (QV - V) / np.maximum(g2, 1e-16), 0.0)[:, None] * g
                x += 0.5 * scatter_add(Q[:, 0], step, len(x)) / qcnt[:, None]
            d = x[J] - x[I]
            L = np.linalg.norm(d, axis=1)
            L = np.where(L < 1e-6, 1e-6, L)
            corr = (0.5 * (L - D) / L)[:, None] * d
            delta = scatter_add(IJ, np.concatenate([corr, -corr]), len(x))
            x += omega * delta / cnt[:, None]
            # bead restraints are annealed away so that the final structure
            # satisfies the geometric constraints
            frac = it / max(1, n_iter - 1)
            strength = bead_strength * float(np.clip((0.85 - frac) / 0.25, 0.0, 1.0))
            if self.g_atoms and strength > 0:
                cog = scatter_add(gid, x[G], len(self.g_atoms)) / gsize[:, None]
                err = T - cog
                x += strength * scatter_add(G, err[gid], len(x)) / gcnt[:, None]
            if verbose and (it % 50 == 0 or it == n_iter - 1):
                dev = np.abs(L - D)
                print(f"  relax iter {it:4d}: max |d-d0| = {dev.max():.4f} nm, rms = {np.sqrt(np.mean(dev**2)):.4f} nm")
        full[used] = x
        return full
