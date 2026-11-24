#!/usr/bin/env bash
# Full presentation case matrix -- REAL UNITS (concrete slab, MN & m).
# Runner defaults already set: L0=20, H0=10, D=23, nu=0.2, f=0.01.
# Stiffeners: 30x50 -> EI=100, GIt=10 ; 30x100 -> EI=825, GIt=82.5 (GIt=0.1*EI).
# Normalized+decay step (nstep) -> scale-invariant, one value fits all cases.
set -e
cd "$(dirname "$0")"
FF="/c/Program Files (x86)/FreeFem++/FreeFem++.exe"
NS="-nstep 0.4 -iters 60"
run(){ echo "==== $1 ===="; "$FF" -nw runner.edp "${@:2}" 2>&1 | grep -E "DONE" || true; }

# -- reference: method demo (hinged slab, off-centred start -> optimum) --------
run ref       -name ref       -nV 2 -nH 0 -bc0 1 -bc1 1 -bc2 1 -bc3 1 -EI 100 -GIt 10  -asym 1 -grid 1 -nstep 0.4 -iters 80

# -- boundary conditions (30x50, nV=2) ----------------------------------------
run bcClamped -name bcClamped -nV 2 -nH 0 -bc0 0 -bc1 0 -bc2 0 -bc3 0 -EI 100 -GIt 10  -grid 1 $NS
run bcOneway  -name bcOneway  -nV 2 -nH 0 -bc0 0 -bc1 2 -bc2 0 -bc3 2 -EI 100 -GIt 10  -grid 1 $NS
run bcHinged  -name bcHinged  -nV 2 -nH 0 -bc0 1 -bc1 1 -bc2 1 -bc3 1 -EI 100 -GIt 10  -grid 1 $NS

# -- number of stiffeners (hinged base, 30x50) --------------------------------
run n0        -name n0        -nV 0 -nH 0 -bc0 1 -bc1 1 -bc2 1 -bc3 1 -EI 100 -GIt 10  -grid 1 -iters 2
run n1        -name n1        -nV 1 -nH 0 -bc0 1 -bc1 1 -bc2 1 -bc3 1 -EI 100 -GIt 10  $NS
run n3        -name n3        -nV 3 -nH 0 -bc0 1 -bc1 1 -bc2 1 -bc3 1 -EI 100 -GIt 10  -grid 1 $NS
run n2h1      -name n2h1      -nV 2 -nH 1 -bc0 1 -bc1 1 -bc2 1 -bc3 1 -EI 100 -GIt 10  -grid 1 $NS

# -- cross-section: 30x100 (= torsion ON) ; and torsion OFF --------------------
run cs825     -name cs825     -nV 2 -nH 0 -bc0 1 -bc1 1 -bc2 1 -bc3 1 -EI 825 -GIt 82.5 -grid 1 $NS
run torOff    -name torOff    -nV 2 -nH 0 -bc0 1 -bc1 1 -bc2 1 -bc3 1 -EI 825 -GIt 0    -grid 1 $NS

echo "==== ALL DONE ===="
for p in ref bcClamped bcOneway bcHinged n0 n1 n3 n2h1 cs825 torOff; do
  echo "---- $p ----"; cat "runs/${p}_pos.txt" 2>/dev/null || echo "MISSING";
done
