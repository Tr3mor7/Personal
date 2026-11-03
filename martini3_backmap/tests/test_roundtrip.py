"""Round-trip tests: atomistic -> Martini 3 (forward mapping) -> backmapped.

Run with:  python -m unittest discover -s tests -v
"""

import os
import sys
import tempfile
import unittest

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "tools"))

from m3backmap import protein as P  # noqa: E402
from m3backmap.backmap import Options, backmap  # noqa: E402
from m3backmap.geometry import dihedral, kabsch  # noqa: E402
from m3backmap.structure import Atom, Residue, System, read_gro, write_gro  # noqa: E402
from m3backmap.templates import load_maps  # noqa: E402
from m3backmap.cli import main as cli_main  # noqa: E402

SEQ = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWUAVTTFSMNGQ"
THREE = {"A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN", "E": "GLU", "G": "GLY",
         "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO", "S": "SER",
         "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL"}


def make_protein(rng, n=60):
    seq = [THREE[c] for c in SEQ if c in THREE][:n]
    phi, psi = P._sample_dihedrals(rng, 1, n)
    X = P.build_backbone(phi[0], psi[0])
    residues = []
    for i, rn in enumerate(seq):
        N, CA, C, O = X[i]
        nc = P.n_chi(rn)
        if rn == "PRO":
            chis = np.array([P.PRO_ROTAMERS[rng.integers(2)]])
        else:
            chis = np.array([[rng.choice([-65.0, 180.0, 65.0]) + rng.normal(0, 8) for _ in range(nc)]])
        sc = {k: v[0] for k, v in P.build_sidechain(rn, N, CA, C, chis).items()}
        pos = {"N": N, "CA": CA, "C": C, "O": O, **sc}
        residues.append((rn, pos))
    return residues


def cg_protein(residues, offset):
    out = []
    for i, (rn, pos) in enumerate(residues):
        beads = P.forward_map_residue(rn, pos)
        out.append(Residue(i + 1, rn, [Atom(i + 1, rn, b, x + offset) for b, x in beads.items()]))
    return out


def cg_lipid(template, xyz, resnr):
    atoms = []
    for k, b in enumerate(template.beads):
        atoms.append(Atom(resnr, template.resnames[0], b, xyz[template.members[k]].mean(0)))
    return Residue(resnr, template.resnames[0], atoms)


def flexed_popc(rng):
    """POPC with random gauche defects in the tails (different from the template)."""
    import build_pc_templates as B
    z = B.pc_lipid(15, 17, (), (9,))
    z2 = []
    for row in z:
        name = row[0]
        if name.startswith(("C2", "C3")) and len(name) >= 3 and name[2:].isdigit() and int(name[2:]) >= 4 \
                and row[6] == 180.0 and rng.random() < 0.25:
            row = row[:6] + (rng.choice([-65.0, 65.0]),)
        z2.append(row)
    pos = B.zbuild(z2)
    return pos


class RoundTrip(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        rng = np.random.default_rng(11)
        cls.maps = load_maps()
        cls.box = np.array([8.0, 8.0, 8.0])
        cls.protein = make_protein(rng)
        residues = cg_protein(cls.protein, np.array([4.0, 4.0, 4.0]) - 0.0)
        # lipids
        t = cls.maps["POPC"]
        cls.lipids = []
        for k in range(3):
            pos = flexed_popc(rng)
            heavy = [n for n in t.atom_names if not n.startswith("H")]
            xyz = np.array([pos[n] for n in heavy])
            # put into template order with hydrogens at their carbons (only heavy atoms are compared)
            full = np.array([pos.get(n, pos[_parent(t, n)]) for n in t.atom_names])
            shift = np.array([1.0 + 2.0 * k, 7.5, 1.0])  # near the box edge -> crosses PBC
            full = full + shift
            cls.lipids.append({n: full[i] for i, n in enumerate(t.atom_names)})
            residues.append(cg_lipid(t, full, 100 + k))
        # wrap lipid beads into the box to exercise PBC unwrapping
        for r in residues[-3:]:
            for a in r.atoms:
                a.xyz = a.xyz % cls.box
        # water and ions
        for k in range(20):
            residues.append(Residue(200 + k, "W", [Atom(200 + k, "W", "W", rng.uniform(0, 8, 3))]))
        residues.append(Residue(300, "ION", [Atom(300, "ION", "NA", np.array([1.0, 1.0, 1.0]))]))
        residues.append(Residue(301, "CL", [Atom(301, "CL", "CL", np.array([2.0, 1.0, 1.0]))]))
        cls.cg = System("test", residues, cls.box)
        cls.aa, cls.report = backmap(cls.cg, cls.maps, Options(verbose=False, seed=3))

    def test_counts(self):
        rep = self.report
        self.assertEqual(rep.counts["SOL"], 80)
        self.assertEqual(rep.counts["POPC"], 3)
        self.assertEqual(rep.counts["SOD"], 1)
        self.assertEqual(rep.counts["CLA"], 1)
        self.assertEqual(rep.chains, [("A", len(self.protein))])

    def test_protein_rmsd(self):
        ref, got, bb_ref, bb_got = [], [], [], []
        prot = [r for r in self.aa.residues if r.chain]
        for (rn, pos), r in zip(self.protein, prot):
            for a in r.atoms:
                if a.name in pos:
                    ref.append(pos[a.name] + 4.0)
                    got.append(a.xyz)
                    if a.name in P.BB_ATOMS:
                        bb_ref.append(pos[a.name] + 4.0)
                        bb_got.append(a.xyz)
        ref, got = np.array(ref), np.array(got)
        rmsd_all = np.sqrt(np.mean(np.sum((ref - got) ** 2, 1)))
        rmsd_bb = np.sqrt(np.mean(np.sum((np.array(bb_ref) - np.array(bb_got)) ** 2, 1)))
        print(f"\n  protein backbone RMSD {rmsd_bb * 10:.2f} A, all heavy atoms {rmsd_all * 10:.2f} A")
        self.assertLess(rmsd_bb, 0.05)
        self.assertLess(rmsd_all, 0.15)

    def test_protein_geometry(self):
        prot = [r for r in self.aa.residues if r.chain]
        cn = [np.linalg.norm(prot[i + 1].get("N").xyz - prot[i].get("C").xyz) for i in range(len(prot) - 1)]
        self.assertLess(max(abs(np.array(cn) - 0.1329)), 0.02)
        # L-chirality: improper C-N-CA-CB around -122 deg
        for r in prot:
            if r.get("CB") is not None:
                d = np.degrees(dihedral(r.get("C").xyz, r.get("N").xyz, r.get("CA").xyz, r.get("CB").xyz))
                self.assertLess(abs(d + 122.5), 20.0, f"{r.resname}{r.resnr} chirality {d}")
        self.assertIsNotNone(prot[-1].get("OXT"))

    def test_lipid_rmsd(self):
        lip = [r for r in self.aa.residues if r.resname == "POPC"]
        for ref, r in zip(self.lipids, lip):
            names = [a.name for a in r.atoms if not a.name.startswith("H")]
            got = np.array([r.get(n).xyz for n in names])
            exp = np.array([ref[n] for n in names])
            # compare modulo a periodic shift of the whole molecule
            shift = np.round((got.mean(0) - exp.mean(0)) / self.box) * self.box
            exp = exp + shift
            rmsd = np.sqrt(np.mean(np.sum((got - exp) ** 2, 1)))
            print(f"  POPC heavy-atom RMSD {rmsd * 10:.2f} A")
            self.assertLess(rmsd, 0.12)
            # bonds sane after relaxation
            t = self.maps["POPC"]
            xyz = np.array([a.xyz for a in r.atoms])
            d = np.array([np.linalg.norm(xyz[i] - xyz[j]) for i, j in t.bonds])
            self.assertLess(d.max(), 0.19)
            self.assertGreater(d.min(), 0.085)

    def test_cli_and_gro_io(self):
        with tempfile.TemporaryDirectory() as tmp:
            cgf = os.path.join(tmp, "cg.gro")
            write_gro(self.cg, cgf)
            out = os.path.join(tmp, "aa.gro")
            rc = cli_main(["backmap", "-c", cgf, "-o", out, "-q", "--prune-water", "0.2"])
            self.assertEqual(rc, 0)
            s = read_gro(out)
            self.assertGreater(s.n_atoms, 600)
            self.assertTrue(os.path.exists(os.path.join(tmp, "aa_protein.pdb")))
            self.assertTrue(os.path.exists(os.path.join(tmp, "aa_nonprotein.gro")))


def _parent(t, hname):
    i = t.atom_names.index(hname)
    for a, b in t.bonds:
        if i in (a, b):
            return t.atom_names[b if a == i else a]


if __name__ == "__main__":
    unittest.main(verbosity=2)
