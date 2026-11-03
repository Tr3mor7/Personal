#!/usr/bin/env bash
# Turn m3backmap output into a minimised all-atom system with GROMACS.
#
#   bash finalize.sh <base> [topol.top]
#
# <base> is the -o name given to m3backmap without extension, i.e. the script
# expects <base>_protein.pdb (if there is protein) and <base>_nonprotein.gro.
#
# Environment variables:
#   FF     force field directory name for pdb2gmx   (default charmm36-jul2022)
#   WATER  water model for pdb2gmx                  (default tip3p)
#   GMX    gromacs binary                           (default gmx)
set -euo pipefail
base=${1:?usage: finalize.sh <base> [topol.top]}
top=${2:-topol.top}
FF=${FF:-charmm36-jul2022}
WATER=${WATER:-tip3p}
GMX=${GMX:-gmx}
here=$(cd "$(dirname "$0")" && pwd)
PY=${PYTHON:-python3}
export PYTHONPATH="$here/..:${PYTHONPATH:-}"

parts=()
if [[ -f ${base}_protein.pdb ]]; then
  echo ">> pdb2gmx: adding hydrogens and building the protein topology"
  $GMX pdb2gmx -f ${base}_protein.pdb -o ${base}_protein_h.gro -p protein.top \
       -i posre.itp -ff "$FF" -water "$WATER" -ignh
  parts+=(${base}_protein_h.gro)
fi
[[ -f ${base}_nonprotein.gro ]] && parts+=(${base}_nonprotein.gro)
if [[ ${#parts[@]} -eq 0 ]]; then parts=(${base}.gro); fi

$PY -m m3backmap merge -f "${parts[@]}" -o ${base}_system.gro

if [[ ! -f $top ]]; then
  cat <<MSG

${base}_system.gro is ready. Create $top now:
  * start from protein.top (written by pdb2gmx) if you have protein,
  * #include the lipid / ion itp files of your force field,
  * append the non-protein [ molecules ] lines printed by m3backmap.
Then run:  bash $0 $base $top
MSG
  exit 0
fi

echo ">> energy minimisation"
$GMX grompp -f "$here/em.mdp" -c ${base}_system.gro -p "$top" -o ${base}_em.tpr -maxwarn 1
$GMX mdrun -deffnm ${base}_em
echo ">> restrained equilibration"
$GMX grompp -f "$here/eq_posre.mdp" -c ${base}_em.gro -r ${base}_em.gro -p "$top" -o ${base}_eq.tpr -maxwarn 1
$GMX mdrun -deffnm ${base}_eq
echo "Done: ${base}_eq.gro"
