#!/usr/bin/env python3
"""Tile a GROMACS system with `gmx_mpi genconf -nbox` and rebuild the topology and index.

`gmx_mpi genconf -nbox NX NY NZ` writes NX*NY*NZ back-to-back copies of the whole input
system: [copy 1: peptide, lipids, water, ions][copy 2: peptide, lipids, ...] ...
A topology whose [ molecules ] counts are simply multiplied by 9 does NOT match that
atom order, and grompp will fail (or, worse, silently mis-assign parameters when atom
names happen to line up). An index file from the original system only covers copy 1.

This script:
  1. runs `gmx_mpi genconf` (or takes an already tiled .gro),
  2. checks the tiled file is N exact copies of the original,
  3. writes the coordinates in the requested molecule order
       grouped (default): every [ molecules ] block gathered across copies, so each
                          count is multiplied by N (e.g. PROA 9, POPC 1152, ...)
       copy:              genconf order kept, the [ molecules ] list repeated N times
  4. writes the topology with a matching [ molecules ] section,
  5. writes an index file where every original group covers all N copies.

Every input and output path is asked for interactively; a path already given on the
command line (e.g. --gro step5_input.gro) is used without asking.
Only the Python standard library is needed.
"""

import argparse
import os
import re
import subprocess
import sys
import tempfile

GMX = "gmx_mpi"  # every GROMACS command is run through gmx_mpi

# ----------------------------------------------------------------------------- gro


class Gro:
    def __init__(self, title, atom_lines, box_line):
        self.title = title
        self.atom_lines = atom_lines
        self.box_line = box_line

    @property
    def natoms(self):
        return len(self.atom_lines)

    def resid(self, i):
        return self.atom_lines[i][0:5].strip()

    def resname(self, i):
        return self.atom_lines[i][5:10].strip()

    def atomname(self, i):
        return self.atom_lines[i][10:15].strip()


def looks_like_pdb(path):
    """True for PDB content: no atom count on line 2, but ATOM/HETATM/CRYST1 records."""
    with open(path) as f:
        lines = f.read().splitlines()
    try:
        int(lines[1].split()[0])
        return False  # a .gro file
    except (IndexError, ValueError):
        return any(l.startswith(("ATOM  ", "HETATM", "CRYST1")) for l in lines)


def read_gro(path):
    with open(path) as f:
        lines = f.read().splitlines()
    try:
        natoms = int(lines[1].split()[0])
    except (IndexError, ValueError):
        sys.exit(f"error: {path} is not a .gro file (line 2 must be the number of atoms)")
    if len(lines) < natoms + 3:
        sys.exit(f"error: {path} declares {natoms} atoms but is truncated")
    return Gro(lines[0], lines[2:2 + natoms], lines[2 + natoms])


def write_gro(path, title, atom_lines, box_line, residue_breaks=()):
    """Write atoms, renumbering residues and atoms sequentially (wrapping at 100000).

    A new residue starts where the residue number/name changes or at any position in
    residue_breaks (the start of each copied block).
    """
    out = [title, f"{len(atom_lines):5d}"]
    resnr = 0
    prev_key = None
    for n, line in enumerate(atom_lines):
        key = (line[0:5], line[5:10])
        if key != prev_key or n in residue_breaks:
            resnr += 1
        prev_key = key
        out.append(f"{resnr % 100000:5d}{line[5:15]}{(n + 1) % 100000:5d}{line[20:]}")
    out.append(box_line)
    with open(path, "w") as f:
        f.write("\n".join(out) + "\n")


# ----------------------------------------------------------------------------- topology

DIRECTIVE = re.compile(r"^\s*\[\s*([A-Za-z_]+)\s*\]")
INCLUDE = re.compile(r'^\s*#include\s+["<](.+?)[">]')


def strip_comment(line):
    return line.split(";", 1)[0].strip()


def include_dirs(top_path, extra):
    dirs = [os.path.dirname(os.path.abspath(top_path))] + list(extra)
    for env in ("GMXLIB",):
        if os.environ.get(env):
            dirs += os.environ[env].split(os.pathsep)
    if os.environ.get("GMXDATA"):
        dirs.append(os.path.join(os.environ["GMXDATA"], "top"))
    return dirs


def scan_moleculetypes(path, dirs, found, visited, unresolved):
    """Collect {moleculetype: [atom names]} from path and everything it #includes."""
    path = os.path.abspath(path)
    if path in visited:
        return
    visited.add(path)
    here = os.path.dirname(path)
    section = None
    current = None
    with open(path) as f:
        for raw in f:
            inc = INCLUDE.match(raw)
            if inc:
                name = inc.group(1)
                for d in [here] + dirs:
                    cand = os.path.join(d, name)
                    if os.path.isfile(cand):
                        scan_moleculetypes(cand, dirs, found, visited, unresolved)
                        break
                else:
                    unresolved.add(name)
                continue
            line = strip_comment(raw)
            if not line or line.startswith("#"):
                continue
            m = DIRECTIVE.match(line)
            if m:
                section = m.group(1).lower()
                continue
            if section == "moleculetype":
                current = line.split()[0]
                found[current] = []
                section = "moleculetype_done"
            elif section == "atoms" and current is not None:
                fields = line.split()
                if len(fields) >= 5:
                    found[current].append(fields[4])


def split_top(top_path):
    """Return (lines before [ molecules ], molecule entries, lines after the entries)."""
    with open(top_path) as f:
        lines = f.read().splitlines()
    start = None
    for i, raw in enumerate(lines):
        m = DIRECTIVE.match(strip_comment(raw))
        if m and m.group(1).lower() == "molecules":
            start = i
    if start is None:
        sys.exit(f"error: no [ molecules ] section in {top_path} (it must be in the main .top)")
    entries = []
    end = start + 1
    while end < len(lines):
        line = strip_comment(lines[end])
        if DIRECTIVE.match(line):
            break
        if line and not line.startswith("#"):
            name, count = line.split()[:2]
            entries.append((name, int(count)))
        end += 1
    return lines[:start + 1], entries, lines[end:]


def molecule_blocks(gro, entries, moltypes):
    """Split the original system into [ molecules ] blocks: (name, count, start, size)."""
    blocks = []
    pos = 0
    for name, count in entries:
        names = moltypes.get(name)
        if names:
            per_mol = len(names)
        else:
            # Moleculetype not readable (e.g. water from a force-field dir that is not on
            # the include path): assume one residue per molecule, as for water/ions/lipids.
            if pos >= gro.natoms:
                sys.exit(f"error: [ molecules ] lists more atoms than the .gro ({name})")
            per_mol = 1
            while (pos + per_mol < gro.natoms
                   and gro.resid(pos + per_mol) == gro.resid(pos)
                   and gro.resname(pos + per_mol) == gro.resname(pos)):
                per_mol += 1
            print(f"note: moleculetype {name} not found in the included .itp files; "
                  f"using {per_mol} atoms/molecule from residue {gro.resname(pos)} in the .gro")
        size = per_mol * count
        if pos + size > gro.natoms:
            sys.exit(f"error: [ molecules ] needs more atoms than the {gro.natoms} in the .gro "
                     f"(ran out at {name})")
        if names:
            got = [gro.atomname(pos + k) for k in range(per_mol)]
            if got != names:
                print(f"warning: first {name} molecule in the .gro has atom names that differ "
                      f"from its .itp; check the [ molecules ] order")
        blocks.append((name, count, pos, size))
        pos += size
    if pos != gro.natoms:
        sys.exit(f"error: [ molecules ] accounts for {pos} atoms but the .gro has {gro.natoms}")
    return blocks


# ----------------------------------------------------------------------------- index


def read_ndx(path):
    groups = []
    with open(path) as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            m = re.match(r"^\[\s*(.+?)\s*\]$", line)
            if m:
                groups.append((m.group(1), []))
            elif groups:
                groups[-1][1].extend(int(x) for x in line.split())
    return groups


def write_ndx(path, groups):
    with open(path, "w") as f:
        for name, atoms in groups:
            f.write(f"[ {name} ]\n")
            for k in range(0, len(atoms), 15):
                f.write(" ".join(f"{a:>6d}" for a in atoms[k:k + 15]) + "\n")
            f.write("\n")


# ----------------------------------------------------------------------------- prompts


def enable_path_completion():
    """Tab-complete file paths at the prompts (where readline is available)."""
    try:
        import glob
        import readline
    except ImportError:
        return
    readline.set_completer_delims(" \t\n")
    readline.set_completer(lambda text, state: (glob.glob(os.path.expanduser(text) + "*") + [None])[state])
    readline.parse_and_bind("tab: complete")


def ask(prompt):
    try:
        return input(prompt).strip()
    except EOFError:
        sys.exit(f"\nerror: no answer for: {prompt.rstrip(': ')}")


def ask_dir(label, default):
    """Prompt for a directory; offer to create it if it does not exist yet."""
    while True:
        path = os.path.expanduser(ask(f"{label} [{default}]: ") or default)
        if os.path.isdir(path):
            return os.path.abspath(path)
        if os.path.exists(path):
            print(f"  not a directory: {path}")
            continue
        if ask(f"  {path} does not exist; create it? [Y/n]: ").lower() in ("", "y", "yes"):
            os.makedirs(path)
            return os.path.abspath(path)


def ask_path(label, default=None, must_exist=True, optional=False, base=None):
    """Prompt until a usable path is entered. Blank returns the default (or None if optional).

    With base, a relative answer (or the default) is taken inside that directory.
    """
    hint = f" [{default}]" if default else (" (blank to skip)" if optional else "")
    while True:
        answer = ask(f"{label}{hint}: ")
        path = os.path.expanduser(answer) if answer else default
        if path and base:
            path = os.path.join(base, path)
        if not path:
            if optional:
                return None
            print("  a path is required")
            continue
        if must_exist and not os.path.isfile(path):
            print(f"  not found: {path}")
            continue
        if not must_exist:
            parent = os.path.dirname(os.path.abspath(path))
            if not os.path.isdir(parent):
                print(f"  directory does not exist: {parent}")
                continue
            if os.path.exists(path):
                if ask(f"  {path} exists; overwrite? [y/N]: ").lower() not in ("y", "yes"):
                    continue
        return path


def fill_paths(a):
    """Ask for every input/output path not given on the command line."""
    enable_path_completion()
    print("Enter the file paths (Tab completes paths).")
    if not a.gro:
        a.gro = ask_path("Original (untiled) coordinate file, .gro or .pdb")
    if not a.top:
        a.top = ask_path("Original topology file, .top")
    if a.ndx is None:
        a.ndx = ask_path("Original index file, .ndx", optional=True) or ""
    if a.tiled_gro is None:
        a.tiled_gro = ask_path(f"Existing {GMX} genconf output, .gro (blank to run {GMX} genconf now)",
                               optional=True) or ""
    if not a.out_dir:
        a.out_dir = ask_dir("Output directory for the generated files", os.getcwd())
    os.makedirs(a.out_dir, exist_ok=True)
    # output names given as flags are also taken relative to the output directory
    for attr in ("out_gro", "out_top", "out_ndx"):
        if getattr(a, attr):
            setattr(a, attr, os.path.join(a.out_dir, os.path.expanduser(getattr(a, attr))))
    stem = os.path.splitext(os.path.basename(a.gro))[0]
    if not a.out_gro:
        a.out_gro = ask_path("Output tiled coordinate file name, .gro", f"{stem}_tiled.gro",
                             must_exist=False, base=a.out_dir)
    if not a.out_top:
        a.out_top = ask_path("Output tiled topology file name, .top", "topol_tiled.top",
                             must_exist=False, base=a.out_dir)
    if a.ndx and not a.out_ndx:
        a.out_ndx = ask_path("Output tiled index file name, .ndx", "index_tiled.ndx",
                             must_exist=False, base=a.out_dir)


def relocate_include(raw, src_dir, dst_dir):
    """Rewrite a relative #include so it still resolves from the output topology's directory."""
    m = INCLUDE.match(raw)
    if not m or os.path.isabs(m.group(1)) or src_dir == dst_dir:
        return raw
    target = os.path.join(src_dir, m.group(1))
    if not os.path.isfile(target):
        return raw  # e.g. a force field found through GMXLIB; leave it to grompp
    return raw.replace(m.group(1), os.path.relpath(target, dst_dir), 1)


# ----------------------------------------------------------------------------- main


def run_gmx(args):
    cmd = [GMX] + args
    print("running:", " ".join(cmd))
    try:
        subprocess.run(cmd, check=True)
    except FileNotFoundError:
        sys.exit(f"error: {GMX} not found on PATH (load your GROMACS module first)")
    except subprocess.CalledProcessError as e:
        sys.exit(f"error: {GMX} {args[0]} failed (exit code {e.returncode}); see its output above")


def load_coords(path, tmp, tag):
    """Read a .gro file; PDB content (whatever the extension) is converted with editconf.

    Returns the coordinates and the path of a .gro file holding them.
    """
    if not looks_like_pdb(path):
        return read_gro(path), path
    print(f"note: {path} is in PDB format; converting it to .gro with {GMX} editconf")
    pdb = os.path.join(tmp, f"{tag}.pdb")  # editconf picks the format from the extension
    gro = os.path.join(tmp, f"{tag}.gro")
    with open(path) as src, open(pdb, "w") as dst:
        dst.write(src.read())
    run_gmx(["editconf", "-f", pdb, "-o", gro])
    return read_gro(gro), gro


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gro", help="original (untiled) coordinates, .gro (PDB is converted)")
    p.add_argument("--top", help="original topology, .top")
    p.add_argument("--ndx", help="original index file; pass '' to go without one")
    p.add_argument("--nbox", nargs=3, type=int, default=[3, 3, 1], metavar=("NX", "NY", "NZ"),
                   help="number of copies in x y z (default: 3 3 1; do not tile a bilayer in z)")
    p.add_argument("--tiled-gro", help=f"output of `{GMX} genconf -nbox` if already run; "
                   "pass '' to have genconf run for you")
    p.add_argument("--order", choices=["grouped", "copy"], default="grouped",
                   help="grouped: each molecule type contiguous (default); copy: genconf order")
    p.add_argument("-I", "--include", action="append", default=[],
                   help="extra directory to search for #include files (repeatable)")
    p.add_argument("--out-dir", help="directory to write the generated files into (created if needed)")
    p.add_argument("--out-gro", help="tiled coordinates to write, .gro (relative to --out-dir)")
    p.add_argument("--out-top", help="tiled topology to write, .top (relative to --out-dir)")
    p.add_argument("--out-ndx", help="tiled index file to write, .ndx (relative to --out-dir)")
    a = p.parse_args(argv)
    fill_paths(a)
    for path in filter(None, (a.gro, a.top, a.ndx, a.tiled_gro)):
        if not os.path.isfile(path):
            sys.exit(f"error: not found: {path}")

    ncopies = a.nbox[0] * a.nbox[1] * a.nbox[2]
    with tempfile.TemporaryDirectory() as tmp:
        orig, orig_gro = load_coords(a.gro, tmp, "original")
        box = [float(x) for x in orig.box_line.split()[:3]]
        if min(box) <= 0:
            sys.exit(f"error: {a.gro} has no periodic box ({orig.box_line.strip()!r}); "
                     "genconf needs the box to place the copies (for a PDB, add a CRYST1 line)")
        if a.tiled_gro:
            tiled, _ = load_coords(a.tiled_gro, tmp, "tiled")
        else:
            raw = os.path.join(tmp, "genconf.gro")
            run_gmx(["genconf", "-f", orig_gro, "-o", raw, "-nbox"] + [str(n) for n in a.nbox])
            tiled = read_gro(raw)

    # genconf must have produced ncopies back-to-back copies of the original atom list.
    n = orig.natoms
    if tiled.natoms != ncopies * n:
        sys.exit(f"error: tiled .gro has {tiled.natoms} atoms, expected {ncopies} x {n} = {ncopies * n}; "
                 f"was it made with -nbox {' '.join(map(str, a.nbox))}?")
    for c in range(ncopies):
        for i in range(n):
            if tiled.atomname(c * n + i) != orig.atomname(i) or tiled.resname(c * n + i) != orig.resname(i):
                sys.exit(f"error: tiled atom {c * n + i + 1} is {tiled.resname(c * n + i)}/"
                         f"{tiled.atomname(c * n + i)}, expected {orig.resname(i)}/{orig.atomname(i)} "
                         f"(copy {c + 1}); the tiled file is not in genconf copy order")

    head, entries, tail = split_top(a.top)
    moltypes, unresolved = {}, set()
    dirs = include_dirs(a.top, a.include)
    scan_moleculetypes(a.top, dirs, moltypes, set(), unresolved)
    if unresolved:
        print("note: could not open #include files (fine if they only hold force-field "
              "parameters): " + ", ".join(sorted(unresolved)))
    blocks = molecule_blocks(orig, entries, moltypes)

    # new_index[c][i]: 0-based position in the output of original atom i in copy c
    new_index = [[0] * n for _ in range(ncopies)]
    order = []
    breaks = set()
    if a.order == "grouped":
        for name, count, start, size in blocks:
            for c in range(ncopies):
                breaks.add(len(order))
                for i in range(start, start + size):
                    new_index[c][i] = len(order)
                    order.append(c * n + i)
        mol_lines = [(name, count * ncopies) for name, count, _, _ in blocks]
    else:
        for c in range(ncopies):
            breaks.add(len(order))
            for i in range(n):
                new_index[c][i] = len(order)
                order.append(c * n + i)
        mol_lines = []
        for _ in range(ncopies):
            for name, count in entries:
                if mol_lines and mol_lines[-1][0] == name:
                    mol_lines[-1] = (name, mol_lines[-1][1] + count)
                else:
                    mol_lines.append((name, count))

    nbox_txt = "x".join(map(str, a.nbox))
    write_gro(a.out_gro, f"{orig.title.strip()} (tiled {nbox_txt})",
              [tiled.atom_lines[k] for k in order], tiled.box_line, breaks)

    width = max(len(nm) for nm, _ in mol_lines)
    body = ["; Compound        #mols"] + [f"{nm:<{max(width, 16)}} {cnt:>8d}" for nm, cnt in mol_lines]
    out_head = []
    in_system = False
    src_dir = os.path.dirname(os.path.abspath(a.top))
    dst_dir = os.path.dirname(os.path.abspath(a.out_top))
    for raw in head:
        raw = relocate_include(raw, src_dir, dst_dir)
        m = DIRECTIVE.match(strip_comment(raw))
        if m:
            in_system = m.group(1).lower() == "system"
            out_head.append(raw)
        elif in_system and strip_comment(raw):
            out_head.append(f"{raw.rstrip()} (tiled {nbox_txt})")
            in_system = False
        else:
            out_head.append(raw)
    with open(a.out_top, "w") as f:
        f.write("\n".join(out_head + body
                           + ([""] + [relocate_include(t, src_dir, dst_dir) for t in tail]
                              if tail else [])) + "\n")

    if a.ndx:
        groups = []
        for name, atoms in read_ndx(a.ndx):
            bad = [x for x in atoms if not 1 <= x <= n]
            if bad:
                sys.exit(f"error: index group {name} has atom {bad[0]} outside 1..{n}")
            groups.append((name, sorted(new_index[c][x - 1] + 1 for c in range(ncopies) for x in atoms)))
        write_ndx(a.out_ndx, groups)

    print(f"\ntiled {n} atoms x {ncopies} copies = {len(order)} atoms ({a.order} order)")
    print(f"box: {orig.box_line.split()[:3]} -> {tiled.box_line.split()[:3]}")
    print("[ molecules ]")
    for nm, cnt in mol_lines:
        print(f"  {nm:<{width}} {cnt}")
    print(f"wrote into {a.out_dir}:")
    for path in filter(None, (a.out_gro, a.out_top, a.out_ndx if a.ndx else None)):
        print(f"  {os.path.relpath(path, a.out_dir)}")


if __name__ == "__main__":
    main()
