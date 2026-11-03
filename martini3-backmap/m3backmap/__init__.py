"""m3backmap - backmap Martini 3 coarse-grained structures to all-atom resolution.

Pipeline
--------
1. Read a Martini 3 ``.gro`` (or ``.pdb``) file and split it into residues.
2. Proteins: rebuild the backbone from the BB-bead trace with a fragment
   library (BBQ-like quadrilateral lookup), then build every side chain from
   ideal internal coordinates, choosing the chi angles that best reproduce the
   SC bead positions.
3. Water / ions: place atomistic molecules on the bead positions.
4. Everything else (lipids, small molecules...): fragment-wise fit of an
   atomistic template described in a ``.map`` file.
5. Optional geometric relaxation that restores ideal bond lengths/angles
   while keeping each bead's atoms centred on the CG bead.

The output is meant to be energy-minimised (and briefly equilibrated with
position restraints) with an atomistic force field afterwards.
"""

__version__ = "0.1.0"
