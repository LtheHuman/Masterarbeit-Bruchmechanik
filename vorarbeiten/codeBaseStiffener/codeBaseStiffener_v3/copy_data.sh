#!/usr/bin/env bash
# Copy the result files the presentation needs into Latex/StiffenedPlate/data/
set -e
cd "$(dirname "$0")"
DST="Latex/StiffenedPlate/data"
mkdir -p "$DST"
for n in ref n0 bcClamped bcHinged bcOneway n1 n3 n2h1 cs825 torOff; do
  cp -f "runs/${n}_obj.txt" "$DST/" 2>/dev/null || true
  cp -f "runs/${n}_pos.txt" "$DST/" 2>/dev/null || true
done
# clean CSV + convert compliance to kJ (MN*m -> kJ is *1000) for readable axes
for f in "$DST"/*_obj.txt; do
  awk -F'[,[:space:]]+' 'NR==1{print "iteration,C_kJ";next} NF>=2{printf "%d,%.6f\n",$1,$2*1000}' "$f" > "$f.t" && mv "$f.t" "$f"
done

# deflection grids -- copied, downsampled 4x, and scaled to mm (col 3 *1000)
ds(){ awk '
  function flush(){ if(kept) print ""; kept=0; col=0 }
  /^[[:space:]]*$/ { rownum++; flush(); next }
  { if(col==0) keeprow=(rownum%2==0); col++;
    if(keeprow && (col%2==1)){ printf "%s %s %g\n",$1,$2,$3*1000; kept=1 } }' "$1" > "$2"; }
for n in ref n0 bcClamped bcHinged bcOneway n3 n2h1 cs825 torOff; do
  [ -f "runs/${n}_w_init.dat"  ] && ds "runs/${n}_w_init.dat"  "$DST/${n}_w_init.dat"
  [ -f "runs/${n}_w_final.dat" ] && ds "runs/${n}_w_final.dat" "$DST/${n}_w_final.dat"
done
echo "copied + downsampled to $DST:"; ls "$DST"
