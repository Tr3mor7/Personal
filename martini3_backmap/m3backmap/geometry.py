"""Geometry helpers: superposition, internal coordinates, PBC handling."""

from __future__ import annotations

import numpy as np


def normalize(v: np.ndarray) -> np.ndarray:
    n = np.linalg.norm(v, axis=-1, keepdims=True)
    return v / np.where(n < 1e-12, 1.0, n)


def kabsch(P: np.ndarray, Q: np.ndarray, w: np.ndarray | None = None):
    """Optimal rotation R and translation t so that ``P @ R.T + t`` ~= Q.

    P, Q: (N, 3). w: optional (N,) weights. Proper rotations only.
    """
    P = np.asarray(P, float)
    Q = np.asarray(Q, float)
    if w is None:
        w = np.ones(len(P))
    w = np.asarray(w, float) / np.sum(w)
    pc = (w[:, None] * P).sum(0)
    qc = (w[:, None] * Q).sum(0)
    H = ((P - pc) * w[:, None]).T @ (Q - qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d if d != 0 else 1.0])
    R = Vt.T @ D @ U.T
    t = qc - pc @ R.T
    return R, t


def rmsd(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.sum((a - b) ** 2, axis=-1))))


def random_rotation(rng: np.random.Generator) -> np.ndarray:
    q = rng.normal(size=4)
    q /= np.linalg.norm(q)
    a, b, c, d = q
    return np.array([
        [a * a + b * b - c * c - d * d, 2 * (b * c - a * d), 2 * (b * d + a * c)],
        [2 * (b * c + a * d), a * a - b * b + c * c - d * d, 2 * (c * d - a * b)],
        [2 * (b * d - a * c), 2 * (c * d + a * b), a * a - b * b - c * c + d * d],
    ])


def rotation_aligning(u: np.ndarray, v: np.ndarray) -> np.ndarray:
    """Rotation matrix that rotates unit vector u onto unit vector v."""
    u = normalize(u)
    v = normalize(v)
    c = float(np.dot(u, v))
    if c > 1 - 1e-10:
        return np.eye(3)
    if c < -1 + 1e-10:
        # 180 degrees about any axis perpendicular to u
        axis = normalize(np.cross(u, [1.0, 0, 0]) if abs(u[0]) < 0.9 else np.cross(u, [0, 1.0, 0]))
        return 2 * np.outer(axis, axis) - np.eye(3)
    k = np.cross(u, v)
    K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    return np.eye(3) + K + K @ K * (1 / (1 + c))


def axis_rotation(axis: np.ndarray, angle: float) -> np.ndarray:
    axis = normalize(np.asarray(axis, float))
    x, y, z = axis
    c, s = np.cos(angle), np.sin(angle)
    C = 1 - c
    return np.array([
        [c + x * x * C, x * y * C - z * s, x * z * C + y * s],
        [y * x * C + z * s, c + y * y * C, y * z * C - x * s],
        [z * x * C - y * s, z * y * C + x * s, c + z * z * C],
    ])


def place_atom(a, b, c, bond, angle, torsion):
    """NeRF: position of D such that |CD| = bond, angle BCD = angle and
    dihedral ABCD = torsion (angles in radians). Fully broadcasting:
    a, b, c are (..., 3); bond/angle/torsion scalars or (...,) arrays."""
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    c = np.asarray(c, float)
    bc = normalize(c - b)
    n = normalize(np.cross(b - a, bc))
    m = np.cross(n, bc)
    bond = np.asarray(bond, float)[..., None]
    angle = np.asarray(angle, float)[..., None]
    torsion = np.asarray(torsion, float)[..., None]
    return (c - bond * np.cos(angle) * bc
            + bond * np.sin(angle) * np.cos(torsion) * m
            + bond * np.sin(angle) * np.sin(torsion) * n)


def dihedral(a, b, c, d) -> np.ndarray:
    b0 = np.asarray(a) - np.asarray(b)
    b1 = normalize(np.asarray(c) - np.asarray(b))
    b2 = np.asarray(d) - np.asarray(c)
    v = b0 - np.sum(b0 * b1, -1, keepdims=True) * b1
    w = b2 - np.sum(b2 * b1, -1, keepdims=True) * b1
    x = np.sum(v * w, -1)
    y = np.sum(np.cross(b1, v) * w, -1)
    return np.arctan2(y, x)


def minimum_image(d: np.ndarray, box: np.ndarray) -> np.ndarray:
    """Apply the minimum image convention to displacement(s) d (..., 3).

    ``box`` is a 3x3 matrix of box vectors (rows), triclinic allowed
    (GROMACS lower-triangular convention)."""
    d = np.array(d, float, copy=True)
    if not np.any(box):
        return d
    for i in (2, 1, 0):
        L = box[i, i]
        if L <= 0:
            continue
        shift = np.round(d[..., i] / L)
        d -= shift[..., None] * box[i]
    return d


def make_whole(coords: np.ndarray, box: np.ndarray, order=None) -> np.ndarray:
    """Unwrap a set of particles so that each is the minimum image of the
    previous one (sequential unwrapping along ``order``)."""
    coords = np.array(coords, float, copy=True)
    if len(coords) < 2 or not np.any(box):
        return coords
    idx = list(range(len(coords))) if order is None else list(order)
    for k in range(1, len(idx)):
        i, j = idx[k - 1], idx[k]
        coords[j] = coords[i] + minimum_image(coords[j] - coords[i], box)
    return coords


def make_whole_tree(coords: np.ndarray, box: np.ndarray) -> np.ndarray:
    """Unwrap a small cluster of particles relative to the first one using a
    nearest-neighbour spanning tree (robust for branched molecules)."""
    coords = np.array(coords, float, copy=True)
    n = len(coords)
    if n < 2 or not np.any(box):
        return coords
    placed = [0]
    remaining = set(range(1, n))
    while remaining:
        best = None
        for i in placed:
            for j in remaining:
                d = np.linalg.norm(minimum_image(coords[j] - coords[i], box))
                if best is None or d < best[0]:
                    best = (d, i, j)
        _, i, j = best
        coords[j] = coords[i] + minimum_image(coords[j] - coords[i], box)
        placed.append(j)
        remaining.remove(j)
    return coords
