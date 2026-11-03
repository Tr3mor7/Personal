"""Geometric hydrogen placement for building atomistic templates."""

from __future__ import annotations

from typing import Dict, List, Sequence

import numpy as np

from .geometry import normalize, place_atom

DEG = np.pi / 180.0
TET = 109.47


def add_hydrogens(pos: Dict[str, np.ndarray], neighbors: Dict[str, List[str]],
                  heavy: str, h_names: Sequence[str], planar: bool = False,
                  bond: float | None = None) -> Dict[str, np.ndarray]:
    """Return positions for the hydrogens ``h_names`` bonded to ``heavy``.

    ``neighbors`` gives the heavy-atom neighbours of each heavy atom. ``planar``
    marks sp2 centres. Bond lengths default to 0.109 nm (C-H) / 0.101 nm (N-H)
    / 0.097 nm (O-H)."""
    if bond is None:
        bond = {"N": 0.101, "O": 0.097, "S": 0.134}.get(heavy[0], 0.109)
    X = pos[heavy]
    nb = neighbors[heavy]
    u = [normalize(pos[n] - X) for n in nb]
    nh = len(h_names)
    out = {}
    if nh == 0:
        return out
    if len(nb) >= 3 and nh == 1:
        d = normalize(-np.sum(u, axis=0))
        out[h_names[0]] = X + bond * d
    elif len(nb) == 2 and nh == 1:
        # sp2 (in-plane bisector) or sp3 methine with 2 heavy neighbours
        d = normalize(-(u[0] + u[1]))
        out[h_names[0]] = X + bond * d
    elif len(nb) == 2 and nh == 2:
        b = normalize(-(u[0] + u[1]))
        n = normalize(np.cross(u[0], u[1]))
        half = 0.5 * TET * DEG
        out[h_names[0]] = X + bond * (np.cos(half) * b + np.sin(half) * n)
        out[h_names[1]] = X + bond * (np.cos(half) * b - np.sin(half) * n)
    elif len(nb) == 1:
        a = nb[0]
        # reference atom for the dihedral: a neighbour of the neighbour
        refs = [r for r in neighbors.get(a, []) if r != heavy]
        if refs:
            ref = pos[refs[0]]
        else:
            ref = pos[a] + np.array([0.0, 0.0, 0.1]) if abs(u[0][2]) < 0.9 else pos[a] + np.array([0.1, 0, 0])
        angle = 120.0 if planar else TET
        if nh == 1:
            tors = [180.0]
        elif nh == 2:
            tors = [0.0, 180.0] if planar else [60.0, -60.0]
        else:
            tors = [180.0, 60.0, -60.0]
        for name, t in zip(h_names, tors):
            out[name] = place_atom(ref, pos[a], X, bond, angle * DEG, t * DEG)
    else:
        raise ValueError(f"Cannot place {nh} H on {heavy} with {len(nb)} heavy neighbours")
    return out
