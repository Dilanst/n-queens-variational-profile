"""Table 1 of [A]: f_h(n) = (n-h+1) N_{h-1}(n) / (h^2 N_h(n)) for h = 1, 2 and n = 10..15.

Reads counts.txt (lines "n h N_h(n)", produced by layer_count), computes the exact fractions
and their rounding to six decimal places, compares them with the values printed in Table 1 of [A]
and writes table1_output.csv.
"""
import os
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
PRINTED = {10: ("0.127393", "0.200891"), 11: ("0.111437", "0.181943"),
           12: ("0.123894", "0.170782"), 13: ("0.122402", "0.165452"),
           14: ("0.109393", "0.158389"), 15: ("0.114223", "0.152947")}

N = {}
with open(os.path.join(HERE, "counts.txt")) as f:
    for line in f:
        if line.strip():
            n, h, c = map(int, line.split())
            N[(n, h)] = c

rows, ok = [], True
for n in sorted(PRINTED):
    for h in (1, 2):
        fr = Fraction((n - h + 1) * N[(n, h - 1)], h * h * N[(n, h)])
        dec = f"{float(fr):.6f}"
        agrees = dec == PRINTED[n][h - 1]
        ok &= agrees
        rows.append(f"{n},{h},{N[(n, h)]},{fr.numerator}/{fr.denominator},{dec},{PRINTED[n][h - 1]},{agrees}")

with open(os.path.join(HERE, "table1_output.csv"), "w") as f:
    f.write("n,h,N_h,f_h_fraction,f_h_decimal,printed,agrees\n" + "\n".join(rows) + "\n")
print("\n".join(rows))
print("All 12 entries agree with Table 1." if ok else "THERE ARE DIFFERENCES with Table 1.")
