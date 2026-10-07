#!/usr/bin/env bash
# Tile a peptide/bilayer system 3x3x1 with gmx_mpi genconf, rebuild topology + index,
# then minimize, equilibrate and run production on the tiled system.
#
# Every GROMACS command runs through gmx_mpi. The script asks for the path of every
# input file and every output file before anything runs (Tab completes paths):
#     bash run_tiled.sh
set -euo pipefail

GMX=gmx_mpi
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

command -v "$GMX" >/dev/null || { echo "error: $GMX not found on PATH (load your GROMACS module)" >&2; exit 1; }

# ask_in VAR "label" [optional] -- existing input file; "optional" allows blank
ask_in() {
    local __var=$1 label=$2 optional=${3:-} path
    while true; do
        read -r -e -p "$label${optional:+ (blank to skip)}: " path
        path=${path/#\~/$HOME}
        if [[ -z $path && -n $optional ]]; then printf -v "$__var" ''; return; fi
        if [[ -f $path ]]; then printf -v "$__var" '%s' "$path"; return; fi
        echo "  not found: ${path:-<blank>}"
    done
}

# ask_out VAR "label" default -- output path; checks the directory, confirms overwrite
ask_out() {
    local __var=$1 label=$2 default=$3 path ok
    while true; do
        read -r -e -p "$label [$default]: " path
        path=${path/#\~/$HOME}
        path=${path:-$default}
        if [[ ! -d $(dirname -- "$path") ]]; then echo "  directory does not exist: $(dirname -- "$path")"; continue; fi
        if [[ -e $path ]]; then
            read -r -p "  $path exists; overwrite? [y/N]: " ok
            [[ $ok == [yY]* ]] || continue
        fi
        printf -v "$__var" '%s' "$path"
        return
    done
}

echo "== Input files =="
ask_in GRO "Original (untiled) coordinate file, .gro"
ask_in TOP "Original topology file, .top"
ask_in NDX "Original index file, .ndx" optional
ask_in EM_MDP "Energy minimization parameter file, .mdp"
read -r -p "Number of equilibration steps [6]: " NEQ
NEQ=${NEQ:-6}
[[ $NEQ =~ ^[0-9]+$ ]] || { echo "error: not a number: $NEQ" >&2; exit 1; }
EQ_MDPS=()
for ((i = 1; i <= NEQ; i++)); do
    ask_in mdp "Equilibration step $i parameter file, .mdp"
    EQ_MDPS+=("$mdp")
done
ask_in PROD_MDP "Production parameter file, .mdp"

echo "== Output files =="
ask_out OUT_GRO "Tiled coordinate file, .gro" tiled.gro
ask_out OUT_TOP "Tiled topology file, .top" tiled.top
OUT_NDX=
[[ -n $NDX ]] && ask_out OUT_NDX "Tiled index file, .ndx" tiled.ndx
# mdrun -deffnm takes a name without extension; each step writes NAME.tpr/.gro/.edr/.log/...
ask_out EM_OUT "Energy minimization output name (no extension)" tiled_em
EQ_OUTS=()
for ((i = 1; i <= NEQ; i++)); do
    ask_out name "Equilibration step $i output name (no extension)" "tiled_eq$i"
    EQ_OUTS+=("$name")
done
ask_out PROD_OUT "Production output name (no extension)" tiled_md

read -r -e -p "MPI launcher for mdrun, e.g. 'mpirun -np 4' (blank for none): " LAUNCH
read -r -e -p "Extra mdrun flags, e.g. '-ntomp 8 -nb gpu' (blank for none): " MDRUN_FLAGS

ndx_flag=()
[[ -n $NDX ]] && ndx_flag=(-n "$OUT_NDX")

mdrun() {  # mdrun DEFFNM
    $LAUNCH "$GMX" mdrun -v -deffnm "$1" $MDRUN_FLAGS
}

# 1. Tile coordinates and write the matching topology and index -----------------------
python3 "$SCRIPT_DIR/tile_system.py" --gro "$GRO" --top "$TOP" --ndx "$NDX" --tiled-gro "" \
    --nbox 3 3 1 --out-gro "$OUT_GRO" --out-top "$OUT_TOP" --out-ndx "$OUT_NDX"

# 2. Energy minimization -----------------------------------------------------------------
# -r: position-restraint reference must be the tiled coordinates, not the original.
"$GMX" grompp -f "$EM_MDP" -c "$OUT_GRO" -r "$OUT_GRO" -p "$OUT_TOP" ${ndx_flag[@]+"${ndx_flag[@]}"} \
    -o "$EM_OUT.tpr"
mdrun "$EM_OUT"

# 3. Equilibration (restraints referenced to the minimized structure) --------------------
prev=$EM_OUT
for ((i = 0; i < NEQ; i++)); do
    name=${EQ_OUTS[i]}
    cpt=()
    [[ -f $prev.cpt ]] && cpt=(-t "$prev.cpt")
    "$GMX" grompp -f "${EQ_MDPS[i]}" -c "$prev.gro" -r "$EM_OUT.gro" -p "$OUT_TOP" ${ndx_flag[@]+"${ndx_flag[@]}"} \
        ${cpt[@]+"${cpt[@]}"} -o "$name.tpr"
    mdrun "$name"
    prev=$name
done

# 4. Production ----------------------------------------------------------------------------
cpt=()
[[ -f $prev.cpt ]] && cpt=(-t "$prev.cpt")
"$GMX" grompp -f "$PROD_MDP" -c "$prev.gro" -p "$OUT_TOP" ${ndx_flag[@]+"${ndx_flag[@]}"} ${cpt[@]+"${cpt[@]}"} \
    -o "$PROD_OUT.tpr"
mdrun "$PROD_OUT"
