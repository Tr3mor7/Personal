"""Tests for map files, make-map, itp reordering and small helpers."""

import os
import sys
import tempfile
import unittest

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from m3backmap.backmap import Options, backmap, near_any, read_itp_orders, reorder_residues  # noqa: E402
from m3backmap.cli import main as cli_main  # noqa: E402
from m3backmap.structure import Atom, Residue, System, read_gro, write_gro, read_pdb, write_pdb  # noqa: E402
from m3backmap.templates import load_maps, read_map  # noqa: E402


class Tools(unittest.TestCase):

    def test_builtin_templates(self):
        maps = load_maps()
        for name, n in (("POPC", 134), ("DOPC", 138), ("DPPC", 130)):
            self.assertIn(name, maps)
            t = maps[name]
            self.assertEqual(len(t.atom_names), n)
            self.assertEqual(len(t.bonds), n - 1)          # acyclic molecule
            self.assertEqual(len(t.beads), 12)

    def test_make_map_roundtrip(self):
        t = load_maps()["POPC"]
        with tempfile.TemporaryDirectory() as tmp:
            aa = System("popc", [Residue(1, "POPC", [Atom(1, "POPC", n, x) for n, x in zip(t.atom_names, t.xyz)])],
                        np.array([5.0, 5.0, 5.0]))
            write_gro(aa, os.path.join(tmp, "popc.gro"))
            with open(os.path.join(tmp, "beads.ndx"), "w") as fh:
                for k, b in enumerate(t.beads):
                    heavy = [n for n, bl in zip(t.atom_names, t.atom_beads) if b in bl and not n.startswith("H")]
                    fh.write(f"[ {b} ]\n{' '.join(heavy)}\n")
            out = os.path.join(tmp, "X.map")
            self.assertEqual(cli_main(["make-map", "-f", os.path.join(tmp, "popc.gro"), "-b",
                                       os.path.join(tmp, "beads.ndx"), "-n", "XPC", "-o", out]), 0)
            t2 = read_map(out)
            self.assertEqual(t2.resnames, ["XPC"])
            self.assertEqual(t2.atom_beads, t.atom_beads)   # hydrogens followed their carbons
            self.assertEqual(len(t2.bonds), len(t.bonds))

    def test_itp_reorder(self):
        with tempfile.TemporaryDirectory() as tmp:
            itp = os.path.join(tmp, "m.itp")
            with open(itp, "w") as fh:
                fh.write("[ moleculetype ]\nMOL 3\n[ atoms ]\n"
                         "1 C 1 MOL C2 1 0 12\n2 C 1 MOL C1 1 0 12\n3 O 1 MOL O1 1 0 16\n")
            s = System("t", [Residue(1, "MOL", [Atom(1, "MOL", n, np.zeros(3)) for n in ("C1", "C2", "O1")])],
                       np.zeros(3))
            msgs = reorder_residues(s, read_itp_orders([itp]))
            self.assertEqual(msgs, [])
            self.assertEqual(s.residues[0].names, ["C2", "C1", "O1"])

    def test_near_any_periodic(self):
        box = np.diag([3.0, 3.0, 3.0])
        q = np.array([[0.05, 1.0, 1.0], [1.5, 1.5, 1.5]])
        r = np.array([[2.95, 1.0, 1.0]])
        self.assertEqual(list(near_any(q, r, 0.2, box)), [True, False])

    def test_unknown_residue(self):
        cg = System("t", [Residue(1, "XYZ", [Atom(1, "XYZ", "B1", np.zeros(3))])], np.array([3.0, 3, 3]))
        with self.assertRaises(SystemExit):
            backmap(cg, load_maps(), Options(verbose=False))
        aa, rep = backmap(cg, load_maps(), Options(verbose=False, skip_unknown=True))
        self.assertEqual(aa.n_atoms, 0)
        self.assertEqual(rep.skipped, {"XYZ": 1})

    def test_pdb_io(self):
        s = System("t", [Residue(5, "ALA", [Atom(5, "ALA", "CA", np.array([0.1, 0.2, 0.3]), "A")], "A")],
                   np.array([4.0, 4.0, 4.0]))
        with tempfile.TemporaryDirectory() as tmp:
            p = os.path.join(tmp, "x.pdb")
            write_pdb(s, p)
            s2 = read_pdb(p)
            np.testing.assert_allclose(s2.residues[0].atoms[0].xyz, [0.1, 0.2, 0.3], atol=1e-4)
            np.testing.assert_allclose(s2.box, [4.0, 4.0, 4.0], atol=1e-4)


if __name__ == "__main__":
    unittest.main(verbosity=2)
