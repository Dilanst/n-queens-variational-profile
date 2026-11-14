#!/usr/bin/env bash
# Reproduce los valores numéricos de [A] y [P].
#   bash reproducir.sh          todo (unos 5 minutos)
#   bash reproducir.sh rapido   comprobación corta (n = 10..13 y mallas 64 y 128; unos 15 segundos)
set -euo pipefail
cd "$(dirname "$0")"
modo="${1:-completo}"

echo "== 1. Conteos exactos (Tabla 1 de [A])"
cd conteos
cc -O2 -o conteo_capas conteo_capas.c
python3 fuerza_bruta.py > fuerza_bruta.txt
for n in 4 5 6; do for h in 0 1 2; do ./conteo_capas $n $h; done; done > contador_pequeno.txt
diff fuerza_bruta.txt contador_pequeno.txt > /dev/null && echo "Contador = fuerza bruta para n = 4, 5, 6."
if [ "$modo" = "rapido" ]; then tamanos="10 11 12 13"; else tamanos="10 11 12 13 14 15"; fi
: > conteos_nuevos.txt
for n in $tamanos; do for h in 0 1 2; do ./conteo_capas $n $h >> conteos_nuevos.txt; done; done
if diff conteos_nuevos.txt <(head -n "$(wc -l < conteos_nuevos.txt)" conteos.txt) > /dev/null; then
  echo "Los conteos recalculados coinciden con conteos.txt."
else
  echo "ATENCIÓN: los conteos recalculados no coinciden con conteos.txt."; exit 1
fi
python3 tabla1.py
cd ..

echo "== 2. Perfil variacional (Sección 9 de [P])"
cd perfil
if [ "$modo" = "rapido" ]; then python3 perfil.py 64 128; else python3 perfil.py 64 128 256 | tee perfil_tabla.txt; fi
