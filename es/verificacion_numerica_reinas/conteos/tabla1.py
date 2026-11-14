"""Tabla 1 de [A]: f_h(n) = (n-h+1) N_{h-1}(n) / (h^2 N_h(n)) para h = 1, 2 y n = 10..15.

Lee conteos.txt (líneas «n h N_h(n)», producidas por conteo_capas), calcula las fracciones exactas
y su redondeo a seis decimales, las compara con los valores impresos en la Tabla 1 de [A] y escribe
tabla1_salida.csv.
"""
import os
from fractions import Fraction

AQUI = os.path.dirname(os.path.abspath(__file__))
IMPRESOS = {10: ("0.127393", "0.200891"), 11: ("0.111437", "0.181943"),
            12: ("0.123894", "0.170782"), 13: ("0.122402", "0.165452"),
            14: ("0.109393", "0.158389"), 15: ("0.114223", "0.152947")}

N = {}
with open(os.path.join(AQUI, "conteos.txt")) as f:
    for linea in f:
        if linea.strip():
            n, h, c = map(int, linea.split())
            N[(n, h)] = c

filas, ok = [], True
for n in sorted(IMPRESOS):
    for h in (1, 2):
        fr = Fraction((n - h + 1) * N[(n, h - 1)], h * h * N[(n, h)])
        dec = f"{float(fr):.6f}"
        coincide = dec == IMPRESOS[n][h - 1]
        ok &= coincide
        filas.append(f"{n},{h},{N[(n, h)]},{fr.numerator}/{fr.denominator},{dec},{IMPRESOS[n][h - 1]},{coincide}")

with open(os.path.join(AQUI, "tabla1_salida.csv"), "w") as f:
    f.write("n,h,N_h,f_h_fraccion,f_h_decimal,impreso,coincide\n" + "\n".join(filas) + "\n")
print("\n".join(filas))
print("Las 12 entradas coinciden con la Tabla 1." if ok else "HAY DIFERENCIAS con la Tabla 1.")
