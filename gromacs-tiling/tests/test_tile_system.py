"""Tests for tile_system.py on a small synthetic peptide/bilayer/water system.

A fake `gmx_mpi` on PATH stands in for GROMACS: its `genconf -nbox` works the way
GROMACS does (whole-system copies back to back, shifted by the box vectors). Every call
is logged.
Run with:  python3 -m unittest discover -s tests
"""

import os
import subprocess
import sys
import tempfile
import textwrap
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
TOOL = os.path.join(HERE, "..", "tile_system.py")

FAKE_GMX = r'''#!/usr/bin/env python3
import os, sys
args = sys.argv[1:]
with open(os.environ["FAKE_GMX_LOG"], "a") as log:
    log.write(" ".join(["gmx_mpi"] + args) + "\n")
if args[0] == "editconf":  # PDB (Angstrom, CRYST1 box) -> .gro (nm)
    f = args[args.index("-f") + 1]; o = args[args.index("-o") + 1]
    assert f.endswith(".pdb"), "editconf needs a .pdb extension to read PDB"
    rows, box = [], None
    for l in open(f):
        if l.startswith("CRYST1"):
            box = [float(l[6:15]) / 10, float(l[15:24]) / 10, float(l[24:33]) / 10]
        elif l.startswith(("ATOM", "HETATM")):
            rows.append(f"{int(l[22:26]):5d}{l[17:21].strip():<5}{l[12:16].strip():>5}{len(rows) + 1:5d}"
                        f"{float(l[30:38]) / 10:8.3f}{float(l[38:46]) / 10:8.3f}{float(l[46:54]) / 10:8.3f}")
    with open(o, "w") as fh:
        fh.write("converted\n" + f"{len(rows):5d}\n" + "\n".join(rows) + "\n")
        fh.write(f"{box[0]:10.5f}{box[1]:10.5f}{box[2]:10.5f}\n")
    sys.exit(0)
assert args[0] == "genconf"
f = args[args.index("-f") + 1]; o = args[args.index("-o") + 1]
i = args.index("-nbox"); nx, ny, nz = map(int, args[i + 1:i + 4])
lines = open(f).read().splitlines()
n = int(lines[1]); atoms = lines[2:2 + n]; box = [float(x) for x in lines[2 + n].split()[:3]]
out = []
for a in range(nx):
    for b in range(ny):
        for c in range(nz):
            for l in atoms:
                x, y, z = float(l[20:28]) + a * box[0], float(l[28:36]) + b * box[1], float(l[36:44]) + c * box[2]
                out.append(f"{l[:20]}{x:8.3f}{y:8.3f}{z:8.3f}{l[44:]}")
with open(o, "w") as fh:
    fh.write(lines[0] + "\n" + f"{len(out):5d}\n" + "\n".join(out) + "\n")
    fh.write(f"{box[0]*nx:10.5f}{box[1]*ny:10.5f}{box[2]*nz:10.5f}\n")
'''

# PROA: 2 residues, 5 atoms | POPC: 3 atoms x2 | SOL: 3 atoms x3 (itp not reachable) | CLA x1
ATOMS = [
    (1, "ALA", "N"), (1, "ALA", "CA"), (1, "ALA", "C"),
    (2, "GLY", "N"), (2, "GLY", "CA"),
    (3, "POPC", "P"), (3, "POPC", "C1"), (3, "POPC", "C2"),
    (4, "POPC", "P"), (4, "POPC", "C1"), (4, "POPC", "C2"),
    (5, "SOL", "OW"), (5, "SOL", "HW1"), (5, "SOL", "HW2"),
    (6, "SOL", "OW"), (6, "SOL", "HW1"), (6, "SOL", "HW2"),
    (7, "SOL", "OW"), (7, "SOL", "HW1"), (7, "SOL", "HW2"),
    (8, "CLA", "CLA"),
]
PROA = ["N", "CA", "C", "N", "CA"]
POPC = ["P", "C1", "C2"]


def read(path):
    with open(path) as f:
        return f.read()


def itp(name, names):
    rows = "\n".join(f"{k + 1:5d}  XX  1  RES  {nm:<4}  {k + 1}  0.0  1.0" for k, nm in enumerate(names))
    return f"[ moleculetype ]\n; name nrexcl\n{name} 3\n\n[ atoms ]\n{rows}\n"


class TileTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        d = self.d = self.tmp.name
        os.makedirs(os.path.join(d, "toppar"))
        os.makedirs(os.path.join(d, "bin"))
        with open(os.path.join(d, "bin", "gmx_mpi"), "w") as f:
            f.write(FAKE_GMX)
        os.chmod(os.path.join(d, "bin", "gmx_mpi"), 0o755)
        self.env = dict(os.environ, PATH=os.path.join(d, "bin") + os.pathsep + os.environ["PATH"],
                        FAKE_GMX_LOG=os.path.join(d, "gmx.log"))
        with open(os.path.join(d, "toppar", "PROA.itp"), "w") as f:
            f.write(itp("PROA", PROA))
        with open(os.path.join(d, "toppar", "POPC.itp"), "w") as f:
            f.write(itp("POPC", POPC))
        with open(os.path.join(d, "toppar", "CLA.itp"), "w") as f:
            f.write(itp("CLA", ["CLA"]))
        with open(os.path.join(d, "topol.top"), "w") as f:
            f.write(textwrap.dedent("""\
                #include "missing_ff.ff/forcefield.itp"
                #include "missing_ff.ff/tip3p.itp"
                #include "toppar/PROA.itp"
                #include "toppar/POPC.itp"
                #include "toppar/CLA.itp"

                [ system ]
                ; Name
                Peptide in POPC

                [ molecules ]
                ; Compound  #mols
                PROA   1
                POPC   2
                SOL    3
                CLA    1
                """))
        with open(os.path.join(d, "in.gro"), "w") as f:
            f.write("test system\n%5d\n" % len(ATOMS))
            for k, (r, rn, an) in enumerate(ATOMS):
                f.write(f"{r:5d}{rn:<5}{an:>5}{k + 1:5d}{0.1 * k:8.3f}{0.2:8.3f}{0.3:8.3f}\n")
            f.write("   4.00000   5.00000   6.00000\n")
        with open(os.path.join(d, "index.ndx"), "w") as f:
            f.write("[ SOLU ]\n1 2 3 4 5\n[ MEMB ]\n6 7 8 9 10 11\n[ SOLV ]\n"
                    + " ".join(str(i) for i in range(12, 22)) + "\n")

    def tearDown(self):
        self.tmp.cleanup()

    def run_tool(self, order):
        d = self.d
        subprocess.run([sys.executable, TOOL,
                        "--gro", os.path.join(d, "in.gro"), "--top", os.path.join(d, "topol.top"),
                        "--ndx", os.path.join(d, "index.ndx"), "--tiled-gro", "", "--nbox", "3", "3", "1",
                        "--order", order, "--out-dir", d, "--out-gro", "out.gro",
                        "--out-top", os.path.join(d, "out.top"), "--out-ndx", os.path.join(d, "out.ndx")],
                       check=True, capture_output=True, text=True, env=self.env)
        return self.read_outputs()

    def read_outputs(self):
        d = self.d
        lines = read(os.path.join(d, "out.gro")).splitlines()
        n = int(lines[1])
        atoms = [(l[5:10].strip(), l[10:15].strip()) for l in lines[2:2 + n]]
        top = read(os.path.join(d, "out.top"))
        mols = [tuple(l.split()) for l in top.split("[ molecules ]")[1].splitlines()
                if l.strip() and not l.startswith(";")]
        ndx = {}
        for block in read(os.path.join(d, "out.ndx")).split("[ ")[1:]:
            name, rest = block.split(" ]", 1)
            ndx[name] = [int(x) for x in rest.split()]
        return lines, atoms, top, [(m, int(c)) for m, c in mols], ndx

    def expand(self, mols, names):
        out = []
        for m, c in mols:
            out += [(m, a) for a in names[m]] * c
        return out

    def check_consistent(self, atoms, mols):
        # Every [ molecules ] entry, expanded atom by atom, must match the .gro atom names.
        names = {"PROA": PROA, "POPC": POPC, "SOL": ["OW", "HW1", "HW2"], "CLA": ["CLA"]}
        self.assertEqual([a for _, a in self.expand(mols, names)], [a for _, a in atoms])

    def test_grouped(self):
        lines, atoms, top, mols, ndx = self.run_tool("grouped")
        self.assertEqual(mols, [("PROA", 9), ("POPC", 18), ("SOL", 27), ("CLA", 9)])
        self.assertEqual(len(atoms), 21 * 9)
        self.check_consistent(atoms, mols)
        self.assertEqual(ndx["SOLU"], list(range(1, 46)))
        self.assertEqual(ndx["MEMB"], list(range(46, 100)))
        self.assertEqual(ndx["SOLV"], list(range(100, 190)))
        self.assertIn("Peptide in POPC (tiled 3x3x1)", top)
        self.assertIn('#include "toppar/PROA.itp"', top)
        self.assertEqual(lines[-1].split(), ["12.00000", "15.00000", "6.00000"])
        # residues renumbered sequentially: 2 per peptide, then one per molecule
        resids = [int(l[0:5]) for l in lines[2:-1]]
        self.assertEqual(resids[:10], [1, 1, 1, 2, 2, 3, 3, 3, 4, 4])
        self.assertEqual(resids[-1], 9 * 2 + 18 + 27 + 9)
        # peptide copy 2 sits one box length along y (genconf loops z fastest, then y, then x)
        self.assertAlmostEqual(float(lines[2 + 5][28:36]) - float(lines[2][28:36]), 5.0)

    def test_copy_order(self):
        _, atoms, _, mols, ndx = self.run_tool("copy")
        self.assertEqual(mols, [("PROA", 1), ("POPC", 2), ("SOL", 3), ("CLA", 1)] * 9)
        self.check_consistent(atoms, mols)
        self.assertEqual(ndx["SOLU"], [c * 21 + i for c in range(9) for i in range(1, 6)])

    def write_pdb(self, path):
        with open(path, "w") as f:
            f.write("REMARK    GENERATED BY CHARMM-GUI\n")
            f.write("CRYST1   40.000   50.000   60.000  90.00  90.00  90.00 P 1           1\n")
            for k, (r, rn, an) in enumerate(ATOMS):
                f.write(f"ATOM  {k + 1:5d} {an:<4} {rn:<4}X{r:4d}    "
                        f"{k:8.3f}{2.0:8.3f}{3.0:8.3f}  1.00  0.00\n")
            f.write("END\n")

    def test_pdb_input(self):
        d = self.d
        # PDB content in files named .pdb and .gro (as from a renamed CHARMM-GUI step5_input.pdb)
        for name in ("in.pdb", "pdb_named.gro"):
            self.write_pdb(os.path.join(d, name))
            open(os.path.join(d, "gmx.log"), "w").close()
            r = subprocess.run([sys.executable, TOOL, "--gro", os.path.join(d, name),
                                "--top", os.path.join(d, "topol.top"), "--ndx", os.path.join(d, "index.ndx"),
                                "--tiled-gro", "", "--out-dir", d, "--out-gro", "out.gro",
                                "--out-top", os.path.join(d, "out.top"), "--out-ndx", os.path.join(d, "out.ndx")],
                               input="", capture_output=True, text=True, env=self.env)
            self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
            self.assertIn("is in PDB format", r.stdout)
            log = [l.split()[1] for l in read(os.path.join(d, "gmx.log")).splitlines()]
            self.assertEqual(log, ["editconf", "genconf"])
            lines, atoms, _, mols, ndx = self.read_outputs()
            self.assertEqual(mols, [("PROA", 9), ("POPC", 18), ("SOL", 27), ("CLA", 9)])
            self.check_consistent(atoms, mols)
            self.assertEqual(ndx["SOLV"], list(range(100, 190)))
            self.assertEqual(lines[-1].split(), ["12.00000", "15.00000", "6.00000"])

    def test_not_a_gro_file(self):
        d = self.d
        with open(os.path.join(d, "junk.gro"), "w") as f:
            f.write("title\nnot a number\n1 2 3\n")
        r = subprocess.run([sys.executable, TOOL, "--gro", os.path.join(d, "junk.gro"),
                            "--top", os.path.join(d, "topol.top"), "--ndx", "", "--tiled-gro", "",
                            "--out-dir", d, "--out-gro", "o.gro", "--out-top", "o.top"],
                           input="", capture_output=True, text=True, env=self.env)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("is not a .gro file", r.stderr)
        self.assertNotIn("Traceback", r.stderr)

    def test_rejects_wrong_nbox(self):
        d = self.d
        subprocess.run(["gmx_mpi", "genconf", "-f", os.path.join(d, "in.gro"),
                        "-o", os.path.join(d, "t.gro"), "-nbox", "2", "2", "1"], check=True, env=self.env)
        r = subprocess.run([sys.executable, TOOL, "--gro", os.path.join(d, "in.gro"),
                            "--top", os.path.join(d, "topol.top"), "--ndx", "", "--tiled-gro", os.path.join(d, "t.gro"),
                            "--out-dir", d, "--out-gro", "o.gro", "--out-top", "o.top"],
                           input="", capture_output=True, text=True, env=self.env)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("expected 9 x 21", r.stderr)


    def test_genconf_uses_gmx_mpi(self):
        self.run_tool("grouped")
        log = read(os.path.join(self.d, "gmx.log")).splitlines()
        self.assertEqual(len(log), 1)
        self.assertTrue(log[0].startswith("gmx_mpi genconf -f "))
        self.assertTrue(log[0].endswith("-nbox 3 3 1"))

    def test_prompts_for_every_path(self):
        d = self.d
        run = os.path.join(d, "new", "run")             # does not exist yet
        answers = "\n".join([
            "nope.gro",                         # missing file: asked again
            os.path.join(d, "in.gro"),
            os.path.join(d, "topol.top"),
            os.path.join(d, "index.ndx"),
            "",                                 # no genconf output yet: run gmx_mpi genconf
            run,                                # output directory
            "",                                 # yes, create it
            "",                                 # accept default names: in_tiled.gro,
            "big.top",                          #   big.top
            "",                                 #   index_tiled.ndx
        ]) + "\n"
        r = subprocess.run([sys.executable, TOOL], input=answers, check=True,
                           capture_output=True, text=True, env=self.env)
        for label in ("Original (untiled) coordinate file", "Original topology file",
                      "Original index file", "Existing gmx_mpi genconf output",
                      "Output directory for the generated files", "does not exist; create it?",
                      "Output tiled coordinate file name", "Output tiled topology file name",
                      "Output tiled index file name"):
            self.assertIn(label, r.stdout)
        self.assertIn("not found: nope.gro", r.stdout)
        self.assertEqual(sorted(os.listdir(run)), ["big.top", "in_tiled.gro", "index_tiled.ndx"])
        os.replace(os.path.join(run, "in_tiled.gro"), os.path.join(d, "out.gro"))
        os.replace(os.path.join(run, "index_tiled.ndx"), os.path.join(d, "out.ndx"))
        top = read(os.path.join(run, "big.top"))
        # topology written to another directory: local includes now point back to toppar/
        self.assertIn('#include "../../toppar/PROA.itp"', top)
        self.assertIn('#include "missing_ff.ff/forcefield.itp"', top)
        with open(os.path.join(d, "out.top"), "w") as f:
            f.write(top)
        _, atoms, _, mols, ndx = self.read_outputs()
        self.assertEqual(mols, [("PROA", 9), ("POPC", 18), ("SOL", 27), ("CLA", 9)])
        self.check_consistent(atoms, mols)
        self.assertEqual(ndx["SOLU"], list(range(1, 46)))

    def test_out_dir_flag(self):
        d = self.d
        out = os.path.join(d, "flagdir")
        subprocess.run([sys.executable, TOOL, "--gro", os.path.join(d, "in.gro"),
                        "--top", os.path.join(d, "topol.top"), "--ndx", os.path.join(d, "index.ndx"),
                        "--tiled-gro", "", "--out-dir", out, "--out-gro", "a.gro", "--out-top", "a.top",
                        "--out-ndx", "a.ndx"], input="", check=True, capture_output=True, text=True,
                       env=self.env)
        self.assertEqual(sorted(os.listdir(out)), ["a.gro", "a.ndx", "a.top"])


if __name__ == "__main__":
    unittest.main()
