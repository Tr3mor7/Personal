"""Protein backmapping for Martini 3 (martinize2-style bead names).

Backbone
    A fragment library is generated from ideal-geometry peptide chains that
    sample the populated regions of the Ramachandran plot. For every residue the
    window of BB beads i-2..i+2 is compared with the library (distance +
    signed-dihedral features, then exact superposition RMSD on the best
    candidates) and the N/CA/C/O atoms of the best fragment are superposed onto
    the CG window. This is the BBQ idea, but trained directly on Martini-3
    BB beads (centre of geometry of N, CA, C, O) instead of C-alpha atoms.

Side chains
    CB is built from N/CA/C with ideal L-geometry. The remaining heavy atoms are
    built from ideal internal coordinates; the chi angles are chosen from a
    rotamer grid (then refined) to reproduce the SC bead positions, with a mild
    clash penalty against the nearby backbone.

Hydrogens are not built: run ``gmx pdb2gmx -ignh`` on the protein output, which
adds them with the naming of your chosen force field.
"""

from __future__ import annotations

import itertools
from functools import lru_cache
from typing import Dict, List, Sequence, Tuple

import numpy as np

from .geometry import kabsch, place_atom, dihedral, normalize, random_rotation

DEG = np.pi / 180.0

# Ideal backbone geometry (nm, degrees) - Engh & Huber
B_N_CA, B_CA_C, B_C_N, B_C_O = 0.1458, 0.1525, 0.1329, 0.1231
A_N_CA_C, A_CA_C_N, A_C_N_CA, A_CA_C_O = 111.2, 116.2, 121.7, 120.5

AMINO_ACIDS = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}
# Alternative residue names (protonation states, force-field specific names)
ALIASES = {
    "HSD": "HIS", "HSE": "HIS", "HSP": "HIS", "HID": "HIS", "HIE": "HIS",
    "HIP": "HIS", "HISD": "HIS", "HISE": "HIS", "HISH": "HIS", "HIS1": "HIS",
    "HIS2": "HIS", "CYX": "CYS", "CYM": "CYS", "CYS2": "CYS", "ASH": "ASP",
    "ASPP": "ASP", "GLH": "GLU", "GLUP": "GLU", "LYN": "LYS", "LSN": "LYS",
    "ARGN": "ARG", "LYSH": "LYS",
}


def canonical(resname: str) -> str | None:
    r = resname.upper()
    if r in AMINO_ACIDS:
        return r
    return ALIASES.get(r)


# ---------------------------------------------------------------------------
# Side-chain internal coordinates (heavy atoms, PDB names)
# Each entry: (atom, ref_a, ref_b, ref_c, bond[nm], angle[deg], torsion)
# torsion is ("chi", k, offset) -> chi_k + offset, or ("fix", value)
# ---------------------------------------------------------------------------

def _c(k, off=0.0):
    return ("chi", k, off)


def _f(v):
    return ("fix", v)


SIDECHAINS: Dict[str, List[tuple]] = {
    "GLY": [],
    "ALA": [],
    "SER": [("OG", "N", "CA", "CB", 0.1417, 110.8, _c(0))],
    "CYS": [("SG", "N", "CA", "CB", 0.1808, 113.8, _c(0))],
    "VAL": [("CG1", "N", "CA", "CB", 0.1527, 110.7, _c(0)),
            ("CG2", "N", "CA", "CB", 0.1527, 110.4, _c(0, 120.0))],
    "THR": [("OG1", "N", "CA", "CB", 0.1430, 109.2, _c(0)),
            ("CG2", "N", "CA", "CB", 0.1530, 111.1, _c(0, -120.0))],
    "ILE": [("CG1", "N", "CA", "CB", 0.1527, 110.7, _c(0)),
            ("CG2", "N", "CA", "CB", 0.1527, 110.4, _c(0, -120.0)),
            ("CD1", "CA", "CB", "CG1", 0.1520, 114.0, _c(1))],
    "LEU": [("CG", "N", "CA", "CB", 0.1530, 116.1, _c(0)),
            ("CD1", "CA", "CB", "CG", 0.1524, 110.3, _c(1)),
            ("CD2", "CA", "CB", "CG", 0.1525, 110.6, _c(1, 120.0))],
    "MET": [("CG", "N", "CA", "CB", 0.1520, 113.7, _c(0)),
            ("SD", "CA", "CB", "CG", 0.1810, 112.7, _c(1)),
            ("CE", "CB", "CG", "SD", 0.1790, 100.6, _c(2))],
    "PRO": [("CG", "N", "CA", "CB", 0.1500, 104.5, _c(0)),
            ("CD", "CA", "CB", "CG", 0.1510, 105.5, _c(1))],
    "PHE": [("CG", "N", "CA", "CB", 0.1500, 113.9, _c(0)),
            ("CD1", "CA", "CB", "CG", 0.1390, 120.0, _c(1)),
            ("CD2", "CA", "CB", "CG", 0.1390, 120.0, _c(1, 180.0)),
            ("CE1", "CB", "CG", "CD1", 0.1390, 120.0, _f(180.0)),
            ("CE2", "CB", "CG", "CD2", 0.1390, 120.0, _f(180.0)),
            ("CZ", "CG", "CD1", "CE1", 0.1390, 120.0, _f(0.0))],
    "TYR": [("CG", "N", "CA", "CB", 0.1510, 113.8, _c(0)),
            ("CD1", "CA", "CB", "CG", 0.1390, 120.8, _c(1)),
            ("CD2", "CA", "CB", "CG", 0.1390, 120.8, _c(1, 180.0)),
            ("CE1", "CB", "CG", "CD1", 0.1390, 121.2, _f(180.0)),
            ("CE2", "CB", "CG", "CD2", 0.1390, 121.2, _f(180.0)),
            ("CZ", "CG", "CD1", "CE1", 0.1390, 119.6, _f(0.0)),
            ("OH", "CD1", "CE1", "CZ", 0.1360, 119.9, _f(180.0))],
    "TRP": [("CG", "N", "CA", "CB", 0.1500, 114.1, _c(0)),
            ("CD1", "CA", "CB", "CG", 0.1370, 127.1, _c(1)),
            ("CD2", "CA", "CB", "CG", 0.1430, 126.6, _c(1, 180.0)),
            ("NE1", "CB", "CG", "CD1", 0.1380, 110.2, _f(180.0)),
            ("CE2", "CB", "CG", "CD2", 0.1400, 107.2, _f(180.0)),
            ("CE3", "CB", "CG", "CD2", 0.1400, 133.9, _f(0.0)),
            ("CZ2", "CG", "CD2", "CE2", 0.1400, 122.3, _f(180.0)),
            ("CZ3", "CG", "CD2", "CE3", 0.1390, 118.8, _f(180.0)),
            ("CH2", "CD2", "CE2", "CZ2", 0.1370, 117.5, _f(0.0))],
    "HIS": [("CG", "N", "CA", "CB", 0.1500, 113.7, _c(0)),
            ("ND1", "CA", "CB", "CG", 0.1380, 122.7, _c(1)),
            ("CD2", "CA", "CB", "CG", 0.1360, 131.0, _c(1, 180.0)),
            ("CE1", "CB", "CG", "ND1", 0.1320, 109.0, _f(180.0)),
            ("NE2", "CB", "CG", "CD2", 0.1370, 107.0, _f(180.0))],
    "ASP": [("CG", "N", "CA", "CB", 0.1520, 113.0, _c(0)),
            ("OD1", "CA", "CB", "CG", 0.1250, 119.2, _c(1)),
            ("OD2", "CA", "CB", "CG", 0.1250, 118.2, _c(1, 180.0))],
    "ASN": [("CG", "N", "CA", "CB", 0.1520, 112.6, _c(0)),
            ("OD1", "CA", "CB", "CG", 0.1230, 120.9, _c(1)),
            ("ND2", "CA", "CB", "CG", 0.1330, 116.5, _c(1, 180.0))],
    "GLU": [("CG", "N", "CA", "CB", 0.1520, 113.8, _c(0)),
            ("CD", "CA", "CB", "CG", 0.1520, 113.3, _c(1)),
            ("OE1", "CB", "CG", "CD", 0.1250, 119.0, _c(2)),
            ("OE2", "CB", "CG", "CD", 0.1250, 118.1, _c(2, 180.0))],
    "GLN": [("CG", "N", "CA", "CB", 0.1520, 113.8, _c(0)),
            ("CD", "CA", "CB", "CG", 0.1520, 112.8, _c(1)),
            ("OE1", "CB", "CG", "CD", 0.1230, 120.9, _c(2)),
            ("NE2", "CB", "CG", "CD", 0.1330, 116.5, _c(2, 180.0))],
    "LYS": [("CG", "N", "CA", "CB", 0.1520, 113.8, _c(0)),
            ("CD", "CA", "CB", "CG", 0.1520, 111.8, _c(1)),
            ("CE", "CB", "CG", "CD", 0.1520, 111.7, _c(2)),
            ("NZ", "CG", "CD", "CE", 0.1490, 111.9, _c(3))],
    "ARG": [("CG", "N", "CA", "CB", 0.1520, 113.8, _c(0)),
            ("CD", "CA", "CB", "CG", 0.1520, 111.8, _c(1)),
            ("NE", "CB", "CG", "CD", 0.1460, 111.7, _c(2)),
            ("CZ", "CG", "CD", "NE", 0.1330, 124.8, _c(3)),
            ("NH1", "CD", "NE", "CZ", 0.1330, 120.3, _f(0.0)),
            ("NH2", "CD", "NE", "CZ", 0.1330, 119.6, _f(180.0))],
}

_SP3 = [-65.0, 180.0, 65.0]
_HALF = [-60.0, -30.0, 0.0, 30.0, 60.0, 90.0]           # 180-periodic sp2 groups
_FULL = list(np.arange(-150.0, 181.0, 30.0))
_RING = [0.0, 30.0, 60.0, 90.0, 120.0, 150.0]

ROTAMER_GRID: Dict[str, List[Sequence[float]]] = {
    "SER": [_SP3], "CYS": [_SP3], "VAL": [_SP3], "THR": [_SP3],
    "ILE": [_SP3, _SP3], "LEU": [_SP3, _SP3], "MET": [_SP3, _SP3, _SP3],
    "PHE": [_SP3, _RING], "TYR": [_SP3, _RING], "TRP": [_SP3, _FULL],
    "HIS": [_SP3, _FULL], "ASP": [_SP3, _HALF], "ASN": [_SP3, _FULL],
    "GLU": [_SP3, _SP3, _HALF], "GLN": [_SP3, _SP3, _FULL],
    "LYS": [_SP3, _SP3, _SP3, _SP3], "ARG": [_SP3, _SP3, _SP3, [180.0, -90.0, 90.0]],
}
# Proline ring puckers (chi1, chi2): endo / exo
PRO_ROTAMERS = [(30.0, -35.0), (-30.0, 40.0)]

# ---------------------------------------------------------------------------
# Martini 3 protein mapping (heavy atoms; centre of geometry). These follow the
# martinize2 Martini 3 mapping files. Edit here if your CG model differs.
# ---------------------------------------------------------------------------
BB_ATOMS = ("N", "CA", "C", "O")

M3_SIDECHAIN_MAP: Dict[str, Dict[str, Tuple[str, ...]]] = {
    "GLY": {},
    "ALA": {"SC1": ("CB",)},
    "CYS": {"SC1": ("CB", "SG")},
    "SER": {"SC1": ("CB", "OG")},
    "THR": {"SC1": ("CB", "OG1", "CG2")},
    "VAL": {"SC1": ("CB", "CG1", "CG2")},
    "LEU": {"SC1": ("CB", "CG", "CD1", "CD2")},
    "ILE": {"SC1": ("CB", "CG1", "CG2", "CD1")},
    "MET": {"SC1": ("CB", "CG", "SD", "CE")},
    "PRO": {"SC1": ("CB", "CG", "CD")},
    "ASN": {"SC1": ("CB", "CG", "OD1", "ND2")},
    "ASP": {"SC1": ("CB", "CG", "OD1", "OD2")},
    "GLN": {"SC1": ("CB", "CG", "CD", "OE1", "NE2")},
    "GLU": {"SC1": ("CB", "CG", "CD", "OE1", "OE2")},
    "LYS": {"SC1": ("CB", "CG", "CD"), "SC2": ("CE", "NZ")},
    "ARG": {"SC1": ("CB", "CG", "CD"), "SC2": ("NE", "CZ", "NH1", "NH2")},
    "HIS": {"SC1": ("CB", "CG"), "SC2": ("ND1", "CE1"), "SC3": ("CD2", "NE2")},
    "PHE": {"SC1": ("CB", "CG"), "SC2": ("CD1", "CE1", "CZ"), "SC3": ("CD2", "CE2", "CZ")},
    "TYR": {"SC1": ("CB", "CG"), "SC2": ("CD1", "CE1"), "SC3": ("CD2", "CE2"),
            "SC4": ("CZ", "OH")},
    "TRP": {"SC1": ("CB", "CG"), "SC2": ("CD1", "NE1"), "SC3": ("CD2", "CE2"),
            "SC4": ("CE3", "CZ3"), "SC5": ("CZ2", "CH2")},
}


def atom_order(resname: str) -> List[str]:
    return list(BB_ATOMS) + ([] if resname == "GLY" else ["CB"]) + [e[0] for e in SIDECHAINS[resname]]


# ---------------------------------------------------------------------------
# Ideal chain building (used for the fragment library and for tests)
# ---------------------------------------------------------------------------

def build_backbone(phi, psi, omega=None) -> np.ndarray:
    """Build ideal backbone(s) from dihedrals (degrees).

    phi, psi: (..., L). Returns (..., L, 4, 3) coordinates of N, CA, C, O (nm).
    """
    phi = np.asarray(phi, float) * DEG
    psi = np.asarray(psi, float) * DEG
    omega = np.full_like(phi, np.pi) if omega is None else np.asarray(omega, float) * DEG
    shape = phi.shape[:-1]
    L = phi.shape[-1]
    X = np.zeros(shape + (L, 4, 3))
    N0 = np.zeros(shape + (3,))
    CA0 = N0 + np.array([B_N_CA, 0, 0])
    ang = (180.0 - A_N_CA_C) * DEG
    C0 = CA0 + B_CA_C * np.array([np.cos(ang), np.sin(ang), 0.0])
    X[..., 0, 0, :], X[..., 0, 1, :], X[..., 0, 2, :] = N0, CA0, C0
    for i in range(L):
        N, CA, C = X[..., i, 0, :], X[..., i, 1, :], X[..., i, 2, :]
        if i + 1 < L:
            Nn = place_atom(N, CA, C, B_C_N, A_CA_C_N * DEG, psi[..., i])
            CAn = place_atom(CA, C, Nn, B_N_CA, A_C_N_CA * DEG, omega[..., i])
            Cn = place_atom(C, Nn, CAn, B_CA_C, A_N_CA_C * DEG, phi[..., i + 1])
            X[..., i + 1, 0, :], X[..., i + 1, 1, :], X[..., i + 1, 2, :] = Nn, CAn, Cn
            X[..., i, 3, :] = place_atom(Nn, CA, C, B_C_O, A_CA_C_O * DEG, np.pi)
        else:
            # last residue: O placed as if followed by a trans residue
            Nn = place_atom(N, CA, C, B_C_N, A_CA_C_N * DEG, psi[..., i])
            X[..., i, 3, :] = place_atom(Nn, CA, C, B_C_O, A_CA_C_O * DEG, np.pi)
    return X


def build_cb(N, CA, C):
    """Ideal CB for L-amino acids (broadcasting over leading dims)."""
    # Ideal-geometry construction (coefficients from trRosetta, in Angstrom).
    b = (CA - N) * 10.0
    c = (C - CA) * 10.0
    a = np.cross(b, c)
    return CA + (-0.58273431 * a + 0.56802827 * b - 0.54067466 * c) / 10.0


def build_sidechain(resname: str, N, CA, C, chis) -> Dict[str, np.ndarray]:
    """Build side-chain heavy atoms for (M,) sets of chi angles.

    N, CA, C: (3,) or (M,3). chis: (M, nchi) in degrees.
    Returns a dict atom -> (M, 3)."""
    chis = np.atleast_2d(np.asarray(chis, float))
    M = chis.shape[0]
    pos = {"N": np.broadcast_to(N, (M, 3)), "CA": np.broadcast_to(CA, (M, 3)),
           "C": np.broadcast_to(C, (M, 3))}
    if resname == "GLY":
        return {}
    pos["CB"] = build_cb(pos["N"], pos["CA"], pos["C"])
    out = {"CB": pos["CB"]}
    for name, ra, rb, rc, bond, ang, tors in SIDECHAINS[resname]:
        if tors[0] == "chi":
            t = (chis[:, tors[1]] + tors[2]) * DEG
        else:
            t = np.full(M, tors[1] * DEG)
        pos[name] = place_atom(pos[ra], pos[rb], pos[rc], bond, ang * DEG, t)
        out[name] = pos[name]
    return out


def n_chi(resname: str) -> int:
    ks = [e[6][1] for e in SIDECHAINS[resname] if e[6][0] == "chi"]
    return (max(ks) + 1) if ks else 0


# ---------------------------------------------------------------------------
# Backbone fragment library
# ---------------------------------------------------------------------------

# (phi, psi, sigma, weight, typical segment length)
_SS_STATES = [
    ((-63.0, -42.0), 8.0, 0.40, (4, 18)),    # alpha helix
    ((-120.0, 130.0), 15.0, 0.25, (3, 10)),  # beta strand
    ((-70.0, 145.0), 12.0, 0.12, (2, 6)),    # PPII
    ((-90.0, 0.0), 20.0, 0.10, (1, 3)),      # bridge / turn
    ((60.0, 40.0), 12.0, 0.05, (1, 2)),      # left-handed helix
    ((-75.0, -20.0), 15.0, 0.08, (2, 4)),    # 3-10 / turns
]

WINDOW = 2  # beads on each side


def _sample_dihedrals(rng, n_chains, L):
    phi = np.zeros((n_chains, L))
    psi = np.zeros((n_chains, L))
    w = np.array([s[2] for s in _SS_STATES])
    w /= w.sum()
    for c in range(n_chains):
        i = 0
        while i < L:
            s = rng.choice(len(_SS_STATES), p=w)
            (p0, q0), sig, _, (lmin, lmax) = _SS_STATES[s]
            seg = rng.integers(lmin, lmax + 1)
            j = min(L, i + seg)
            phi[c, i:j] = p0 + rng.normal(0, sig, j - i)
            psi[c, i:j] = q0 + rng.normal(0, sig, j - i)
            i = j
    return phi, psi


@lru_cache(maxsize=2)
def backbone_library(n_chains: int = 400, length: int = 104, seed: int = 2024):
    """Return (bb_windows (M,5,3), atoms (M,4,3)) centred on the central bead."""
    rng = np.random.default_rng(seed)
    phi, psi = _sample_dihedrals(rng, n_chains, length)
    X = build_backbone(phi, psi)                     # (C, L, 4, 3)
    BB = X.mean(axis=2)                              # (C, L, 3) COG of N CA C O
    win, atoms = [], []
    for i in range(WINDOW, length - WINDOW):
        w = BB[:, i - WINDOW:i + WINDOW + 1]         # (C, 5, 3)
        win.append(w - BB[:, i:i + 1])
        atoms.append(X[:, i] - BB[:, i:i + 1])
    win = np.concatenate(win, 0)
    atoms = np.concatenate(atoms, 0)
    return win, atoms


def _features(win: np.ndarray, mask: Tuple[bool, ...]) -> np.ndarray:
    """Distance + signed-dihedral features of a (M, 5, 3) window restricted
    to the available positions in ``mask``."""
    idx = [k for k, m in enumerate(mask) if m]
    P = win[:, idx]
    feats = []
    for a, b in itertools.combinations(range(len(idx)), 2):
        if idx[b] - idx[a] >= 2:
            feats.append(np.linalg.norm(P[:, a] - P[:, b], axis=-1))
    for s in range(len(idx) - 3):
        if idx[s + 3] - idx[s] == 3:
            d = dihedral(P[:, s], P[:, s + 1], P[:, s + 2], P[:, s + 3])
            feats.append(0.15 * np.sin(d))
            feats.append(0.15 * np.cos(d))
    if not feats:
        return np.zeros((len(win), 1))
    return np.stack(feats, -1)


def _batched_kabsch_rmsd(P: np.ndarray, Q: np.ndarray):
    """P: (K, n, 3) candidates, Q: (n, 3) target -> rotations (K,3,3), rmsd (K,)."""
    pc = P.mean(1, keepdims=True)
    qc = Q.mean(0, keepdims=True)
    P0 = P - pc
    Q0 = Q - qc
    H = np.einsum("kni,nj->kij", P0, Q0)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(np.einsum("kji,klj->kil", Vt, U)))
    D = np.ones((len(P), 3))
    D[:, 2] = np.where(d == 0, 1.0, d)
    R = np.einsum("kji,kj,klj->kil", Vt, D, U)       # V D U^T
    fit = np.einsum("kni,kji->knj", P0, R)
    r = np.sqrt(np.mean(np.sum((fit - Q0) ** 2, -1), -1))
    return R, r, pc[:, 0], qc[0]


def rebuild_backbone(bb: np.ndarray, rng: np.random.Generator, n_candidates: int = 24) -> np.ndarray:
    """Reconstruct N, CA, C, O for a chain from its (L, 3) BB bead trace.

    Returns (L, 4, 3)."""
    L = len(bb)
    out = np.zeros((L, 4, 3))
    lib_win, lib_atoms = backbone_library()
    if L == 1:
        R = random_rotation(rng)
        out[0] = lib_atoms[0] @ R.T + bb[0]
        return out
    # group residues by window availability pattern
    groups: Dict[Tuple[bool, ...], List[int]] = {}
    for i in range(L):
        mask = tuple(0 <= i + k < L for k in range(-WINDOW, WINDOW + 1))
        groups.setdefault(mask, []).append(i)
    for mask, residues in groups.items():
        idx = [k for k, m in enumerate(mask) if m]
        lib_f = _features(lib_win, mask)
        q_win = np.zeros((len(residues), 2 * WINDOW + 1, 3))
        for r, i in enumerate(residues):
            for k in idx:
                q_win[r, k] = bb[i + k - WINDOW] - bb[i]
        q_f = _features(q_win, mask)
        lib_n2 = np.sum(lib_f ** 2, 1)
        for start in range(0, len(residues), 256):
            qf = q_f[start:start + 256]
            d2 = np.sum(qf ** 2, 1)[:, None] + lib_n2[None] - 2 * qf @ lib_f.T
            k = min(n_candidates, len(lib_f))
            cand = np.argpartition(d2, k - 1, axis=1)[:, :k]
            for r, i in enumerate(residues[start:start + 256]):
                c = cand[r]
                P = lib_win[c][:, idx]
                Q = q_win[start + r][idx]
                R, rm, pc, qc = _batched_kabsch_rmsd(P, Q)
                best = int(np.argmin(rm))
                atoms = (lib_atoms[c[best]] - pc[best]) @ R[best].T + qc
                out[i] = atoms + bb[i]
    return out


# ---------------------------------------------------------------------------
# Side-chain fitting
# ---------------------------------------------------------------------------

def _sidechain_score(resname, built, targets, env, w_clash=20.0):
    M = next(iter(built.values())).shape[0]
    score = np.zeros(M)
    for bead, atoms in M3_SIDECHAIN_MAP[resname].items():
        if bead not in targets:
            continue
        cog = np.mean([built[a] for a in atoms], axis=0)
        score += np.sum((cog - targets[bead]) ** 2, -1)
    if env is not None and len(env):
        far = [a for a in built if a != "CB"]
        if far:
            P = np.stack([built[a] for a in far], 1)                 # (M, n, 3)
            d = np.linalg.norm(P[:, :, None, :] - env[None, None], axis=-1)
            score += w_clash * np.sum(np.clip(0.28 - d, 0, None) ** 2, axis=(1, 2))
    return score


def fit_sidechain(resname: str, N, CA, C, targets: Dict[str, np.ndarray], env=None):
    """Return dict of side-chain atom positions best matching the SC beads."""
    if resname == "GLY":
        return {}
    nc = n_chi(resname)
    if nc == 0:
        return {k: v[0] for k, v in build_sidechain(resname, N, CA, C, np.zeros((1, 0))).items()}
    if resname == "PRO":
        grid = np.array(PRO_ROTAMERS)
    else:
        grid = np.array(list(itertools.product(*ROTAMER_GRID[resname])), float)
    if not any(b in targets for b in M3_SIDECHAIN_MAP[resname]):
        # no SC information at all: take the most common rotamer
        best = grid[0]
    else:
        sc = _sidechain_score(resname, build_sidechain(resname, N, CA, C, grid), targets, env)
        best = grid[int(np.argmin(sc))]
        best_score = float(np.min(sc))
        start = best.copy()
        if resname != "PRO":
            for step in (20.0, 10.0, 5.0):
                for _ in range(4):
                    offs = np.array(list(itertools.product((-step, 0.0, step), repeat=nc)))
                    cand = best[None] + offs
                    cand = np.clip(cand, start - 45.0, start + 45.0)
                    s = _sidechain_score(resname, build_sidechain(resname, N, CA, C, cand), targets, env)
                    j = int(np.argmin(s))
                    if s[j] < best_score - 1e-9:
                        best, best_score = cand[j], float(s[j])
                    else:
                        break
    built = build_sidechain(resname, N, CA, C, best[None])
    return {k: v[0] for k, v in built.items()}


def build_oxt(N, CA, C, O):
    """Second carboxylate oxygen of the C-terminus."""
    return place_atom(O, CA, C, 0.1250, 117.0 * DEG, np.pi)


# ---------------------------------------------------------------------------
# Forward mapping (useful for testing / validation)
# ---------------------------------------------------------------------------

def forward_map_residue(resname: str, atoms: Dict[str, np.ndarray]) -> Dict[str, np.ndarray]:
    beads = {"BB": np.mean([atoms[a] for a in BB_ATOMS], axis=0)}
    for bead, names in M3_SIDECHAIN_MAP[resname].items():
        beads[bead] = np.mean([atoms[a] for a in names], axis=0)
    return beads
