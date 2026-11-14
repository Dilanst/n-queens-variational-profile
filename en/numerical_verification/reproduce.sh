#!/usr/bin/env bash
# Reproduces the numerical values of [A] and [P]; the references are in the README.
#   bash reproduce.sh          everything (about four and a half minutes)
#   bash reproduce.sh quick    short check: n = 10..13 and grids of 64 and 128 points
#                              (about 15 seconds)
set -euo pipefail
cd "$(dirname "$0")"
mode="${1:-full}"

echo "== 1. Exact counts (Table 1 of [A])"
cd counts
cc -O2 -o layer_count layer_count.c
python3 brute_force.py > brute_force.txt
for n in 4 5 6; do for h in 0 1 2; do ./layer_count $n $h; done; done > small_counter.txt
if diff brute_force.txt small_counter.txt > /dev/null; then
  echo "The counter agrees with the brute-force search for n = 4, 5, 6."
else
  echo "ERROR: the counter does not agree with the brute-force search."; exit 1
fi
if [ "$mode" = "quick" ]; then sizes="10 11 12 13"; else sizes="10 11 12 13 14 15"; fi
: > new_counts.txt
for n in $sizes; do for h in 0 1 2; do ./layer_count $n $h >> new_counts.txt; done; done
if diff new_counts.txt <(head -n "$(wc -l < new_counts.txt)" counts.txt) > /dev/null; then
  echo "The recomputed counts agree with counts.txt."
else
  echo "ERROR: the recomputed counts do not agree with counts.txt."; exit 1
fi
python3 table1.py
cd ..

echo "== 2. Variational profile (Section 9 of [P])"
cd profile
if [ "$mode" = "quick" ]; then python3 solve_profile.py 64 128; else python3 solve_profile.py 64 128 256 | tee profile_table.txt; fi
