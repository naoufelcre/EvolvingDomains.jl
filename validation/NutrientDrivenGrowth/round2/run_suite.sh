#!/bin/bash
# Campaign II / Experiment 5: full run matrix.  Each run writes output/<tag>/.
set -u
cd "$(dirname "$0")"
export EVOLVINGDOMAINS_CUTCELLS_LIBRARY="/Users/naoufel/Documents/Scientific stack/EvolvingDomains-Redesign/validation/lib/cutcells/lib/libevolvingdomains_cut_cells.dylib"
PROJ="/Users/naoufel/Documents/Scientific stack/EvolvingDomains-Redesign/validation/env"
J=(julia --project="$PROJ" --startup-file=no)

run() { # ncells tfinal method dtmax tag
  echo "### $(date +%H:%M:%S) start $5 : n=$1 t=$2 method=$3 dtmax=$4"
  "${J[@]}" run.NutrientGrowthRefined.jl run "$1" "$2" "$3" "$4" "$5" > "output/$5.log" 2>&1
  echo "### $(date +%H:%M:%S) done  $5 : exit=$?"
}

mkdir -p output
run 32  0.5 A 5.0e-3 n32_A
run 32  0.5 B 5.0e-3 n32_B
run 64  0.5 A 5.0e-3 n64_A
run 64  0.5 B 5.0e-3 n64_B
run 96  0.5 A 5.0e-3 n96_A
run 96  0.5 B 5.0e-3 n96_B
run 128 0.5 A 5.0e-3 n128_A
run 128 0.5 B 5.0e-3 n128_B
run 64  0.5 A 1.0e-3 n64_A_dt1e-3
run 64  0.5 B 1.0e-3 n64_B_dt1e-3
run 64  0.5 A 2.5e-3 n64_A_dt2p5e-3
echo "### suite complete $(date +%H:%M:%S)"
