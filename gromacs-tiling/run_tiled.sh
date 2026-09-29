#!/usr/bin/env bash
# Tile a peptide/bilayer system 3x3x1 with gmx genconf, rebuild topology + index,
# then minimize, equilibrate and run production on the tiled system.
#
# Edit the variables below (defaults follow CHARMM-GUI Membrane Builder file names),
# then run from the directory holding your system files:
#     bash run_tiled.sh
set -euo pipefail

GMX=${GMX:-gmx}                       # e.g. "gmx_mpi" or "gmx"
MDRUN_FLAGS=${MDRUN_FLAGS:-}          # e.g. "-ntomp 8 -nb gpu"
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

GRO=step5_input.gro                   # original, untiled coordinates
TOP=topol.top                         # original topology
NDX=index.ndx                         # original index file
NBOX="3 3 1"                          # copies in x y z (never tile a bilayer in z)

EM_MDP=step6.0_minimization.mdp
EQ_MDPS=(step6.1_equilibration.mdp step6.2_equilibration.mdp step6.3_equilibration.mdp
         step6.4_equilibration.mdp step6.5_equilibration.mdp step6.6_equilibration.mdp)
PROD_MDP=step7_production.mdp

OUT=tiled                             # prefix for all tiled files

# 1. Tile coordinates and write the matching topology and index -------------------
python3 "$SCRIPT_DIR/tile_system.py" \
    --gmx "$GMX" --gro "$GRO" --top "$TOP" --ndx "$NDX" --nbox $NBOX \
    --out-gro "$OUT.gro" --out-top "$OUT.top" --out-ndx "$OUT.ndx"

# 2. Energy minimization -----------------------------------------------------------
# -r: position-restraint reference must be the tiled coordinates, not the original.
$GMX grompp -f "$EM_MDP" -c "$OUT.gro" -r "$OUT.gro" -p "$OUT.top" -n "$OUT.ndx" \
    -o "${OUT}_em.tpr"
$GMX mdrun -v -deffnm "${OUT}_em" $MDRUN_FLAGS

# 3. Equilibration (restraints referenced to the minimized structure) --------------
prev="${OUT}_em"
for mdp in "${EQ_MDPS[@]}"; do
    name="${OUT}_$(basename "$mdp" .mdp)"
    cpt=()
    [[ -f "$prev.cpt" ]] && cpt=(-t "$prev.cpt")
    $GMX grompp -f "$mdp" -c "$prev.gro" -r "${OUT}_em.gro" -p "$OUT.top" -n "$OUT.ndx" \
        "${cpt[@]}" -o "$name.tpr"
    $GMX mdrun -v -deffnm "$name" $MDRUN_FLAGS
    prev="$name"
done

# 4. Production -----------------------------------------------------------------------
name="${OUT}_$(basename "$PROD_MDP" .mdp)"
$GMX grompp -f "$PROD_MDP" -c "$prev.gro" -t "$prev.cpt" -p "$OUT.top" -n "$OUT.ndx" \
    -o "$name.tpr"
$GMX mdrun -v -deffnm "$name" $MDRUN_FLAGS
