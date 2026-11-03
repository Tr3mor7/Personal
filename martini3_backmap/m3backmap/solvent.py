"""Martini 3 water and ion beads -> atomistic water / ions."""

from __future__ import annotations

from typing import Dict, List, Tuple

import numpy as np

from .geometry import random_rotation

DEG = np.pi / 180.0

# Martini 3 water beads: regular (4:1), small (3:1), tiny (2:1)
WATER_BEADS = {"W": 4, "SW": 3, "TW": 2, "WN": 4, "PW": 4}
WATER_RESNAMES = {"W", "SW", "TW", "PW", "WN", "SOL"}

# Martini 3 ion bead name -> element key
ION_KEYS = {
    "NA": "NA", "NA+": "NA", "SOD": "NA", "TNA": "NA",
    "K": "K", "K+": "K", "POT": "K", "TK": "K",
    "CL": "CL", "CL-": "CL", "CLA": "CL", "TCL": "CL",
    "CA": "CA", "CA+": "CA", "CAL": "CA",
    "MG": "MG", "MG+": "MG",
    "BR": "BR", "BR-": "BR", "I": "I", "I-": "I",
}
ION_RESNAMES = {"ION", "NA", "CL", "K", "CA", "MG", "NA+", "CL-", "SOD", "CLA", "POT", "CAL", "BR", "I"}

ION_NAMES = {
    # force field family -> element key -> (resname, atomname)
    "charmm": {"NA": ("SOD", "SOD"), "K": ("POT", "POT"), "CL": ("CLA", "CLA"),
               "CA": ("CAL", "CAL"), "MG": ("MG", "MG"), "BR": ("BR", "BR"), "I": ("I", "I")},
    "amber": {"NA": ("NA", "NA"), "K": ("K", "K"), "CL": ("CL", "CL"),
              "CA": ("CA", "CA"), "MG": ("MG", "MG"), "BR": ("BR", "BR"), "I": ("I", "I")},
    "gromos": {"NA": ("NA", "NA"), "K": ("K", "K"), "CL": ("CL", "CL"),
               "CA": ("CA", "CA"), "MG": ("MG", "MG"), "BR": ("BR", "BR"), "I": ("I", "I")},
}

WATER_MODELS = {
    # name: (OH bond nm, HOH angle deg, atom names, virtual site distance or None)
    "tip3p": (0.09572, 104.52, ("OW", "HW1", "HW2"), None),
    "spce": (0.1, 109.47, ("OW", "HW1", "HW2"), None),
    "spc": (0.1, 109.47, ("OW", "HW1", "HW2"), None),
    "tip4p": (0.09572, 104.52, ("OW", "HW1", "HW2", "MW"), 0.015),
    "opc3": (0.09789, 109.47, ("OW", "HW1", "HW2"), None),
}


def is_water(resname: str, bead_names: List[str]) -> bool:
    return resname.upper() in WATER_RESNAMES and all(b.upper() in WATER_BEADS for b in bead_names)


def ion_key(resname: str, bead_name: str) -> str | None:
    for k in (bead_name.upper(), resname.upper()):
        if k in ION_KEYS:
            return ION_KEYS[k]
    return None


def water_molecule(model: str, rng: np.random.Generator) -> np.ndarray:
    """Randomly oriented water with O at the origin: (n_sites, 3)."""
    bond, ang, names, vs = WATER_MODELS[model]
    a = ang * DEG / 2
    h1 = bond * np.array([np.sin(a), np.cos(a), 0.0])
    h2 = bond * np.array([-np.sin(a), np.cos(a), 0.0])
    sites = [np.zeros(3), h1, h2]
    if vs is not None:
        sites.append(np.array([0.0, vs, 0.0]))
    R = random_rotation(rng)
    return np.array(sites) @ R.T


def cluster_offsets(n: int, rng: np.random.Generator, spacing: float = 0.28) -> np.ndarray:
    """Positions of ``n`` water oxygens around a bead centre (centred at 0)."""
    if n == 1:
        pts = np.zeros((1, 3))
    elif n == 2:
        pts = np.array([[0, 0, 0.5], [0, 0, -0.5]]) * spacing
    elif n == 3:
        r = spacing / np.sqrt(3)
        pts = np.array([[r * np.cos(t), r * np.sin(t), 0] for t in (0, 2 * np.pi / 3, 4 * np.pi / 3)])
    else:
        tet = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], float)
        pts = tet[:n] * spacing / (2 * np.sqrt(2))
        if n > 4:
            extra = rng.normal(size=(n - 4, 3))
            pts = np.vstack([pts, extra / np.linalg.norm(extra, axis=1)[:, None] * spacing])
    pts = pts - pts.mean(0)
    return pts @ random_rotation(rng).T


def place_waters(center: np.ndarray, n: int, model: str, rng: np.random.Generator,
                 exclude_center: bool = False) -> List[np.ndarray]:
    """List of water molecules (each (n_sites,3)) around ``center``.

    With ``exclude_center`` the waters are put on a shell of radius ~0.24 nm
    (used for hydration shells of ions)."""
    if exclude_center:
        v = rng.normal(size=(n, 3))
        offs = v / np.linalg.norm(v, axis=1)[:, None] * 0.24
    else:
        offs = cluster_offsets(n, rng)
    return [center + o + water_molecule(model, rng) for o in offs]
