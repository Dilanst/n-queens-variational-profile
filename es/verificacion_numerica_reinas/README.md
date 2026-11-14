# Verificación numérica

Estos programas reproducen los valores numéricos de dos artículos:

- **[A]** D. S. Mejia, *Completabilidad de las últimas capas en el problema de las n reinas*, 2026. Repositorio: [n-queens-last-layers](https://github.com/Dilanst/n-queens-last-layers).
- **[P]** D. S. Mejia, *El perfil variacional de las colocaciones parciales de reinas: existencia, ecuaciones de Euler y regularidad hasta el tablero lleno*, 2026. DOI: [10.5281/zenodo.24220042](https://doi.org/10.5281/zenodo.24220042). Repositorio: [n-queens-variational-profile](https://github.com/Dilanst/n-queens-variational-profile).

## Cómo reproducirlos

Hace falta un compilador de C (gcc o clang) y Python 3 con NumPy (`pip install -r requirements.txt`). Desde esta carpeta:

```
bash reproducir.sh          # todo, unos 4 minutos y medio
bash reproducir.sh rapido   # comprobación corta, unos 15 segundos
```

El guion compila el contador y lo contrasta con una búsqueda por fuerza bruta. Después recalcula los conteos y los compara con los guardados; si algo no coincide, lo avisa y se detiene. Por último resuelve el perfil e imprime sus valores.

## Equipo y tiempos

Los tiempos de esta sección se midieron en este equipo:

- Ubuntu 24.04 (Linux)
- Procesador Intel Xeon a 2,1 GHz, con 2 núcleos
- 8 GB de RAM
- GCC 13.3, Python 3.13 y NumPy 2.5

Los programas no trabajan en paralelo: cada uno usa un solo núcleo.

| Paso | Tiempo |
|---|---|
| Compilar el contador y compararlo con la fuerza bruta | 1 s |
| Conteos para n = 10, …, 15 | 3 min 34 s, de los que n = 15, h = 2 se lleva 2 min 37 s |
| Perfil en las mallas de 64, 128 y 256 puntos | 54 s, de los que la malla de 256 se lleva 43 s |
| Total (`bash reproducir.sh`) | Unos 4 min 30 s |

Los resultados que se obtienen son idénticos, byte a byte, a los archivos guardados.

## 1. Conteos exactos (Tabla 1 de [A])

- **N_h(n)** es el número de colocaciones de n − h reinas que no se atacan en el tablero n × n.
- **f_h(n) = (n − h + 1) N_{h−1}(n) / (h² N_h(n))** es la fracción de la ecuación (1.1) de [A].

| Archivo | Contenido |
|---|---|
| `conteos/conteo_capas.c` | Contador: recorre todas las colocaciones fila a fila. Uso: `./conteo_capas n h` |
| `conteos/fuerza_bruta.py` | Control para n = 4, 5, 6: prueba todos los subconjuntos de casillas |
| `conteos/tabla1.py` | Calcula f_h(n) como fracción exacta, la redondea a seis decimales y la compara con la Tabla 1 de [A] |
| `conteos/conteos.txt` | N₀, N₁ y N₂ para n = 10, …, 15 |
| `conteos/tabla1_salida.csv` | Lo que escribe `tabla1.py` |

| n | N₀ | N₁ | N₂ | f₁(n) | f₂(n) |
|---|---|---|---|---|---|
| 10 | 724 | 56 832 | 636 524 | 905/7104 ≈ 0,127393 | 31968/159131 ≈ 0,200891 |
| 11 | 2 680 | 264 544 | 3 634 976 | 3685/33068 ≈ 0,111437 | 41335/227186 ≈ 0,181943 |
| 12 | 14 200 | 1 375 368 | 22 146 756 | 7100/57307 ≈ 0,123894 | 630377/3691126 ≈ 0,170782 |
| 13 | 73 712 | 7 828 784 | 141 952 580 | 59891/489299 ≈ 0,122402 | 5871588/35488145 ≈ 0,165452 |
| 14 | 365 596 | 46 788 624 | 960 059 292 | 639793/5848578 ≈ 0,109393 | 745407/4706173 ≈ 0,158389 |
| 15 | 2 279 184 | 299 307 416 | 6 849 285 616 | 4273470/37413427 ≈ 0,114223 | 261893989/1712321404 ≈ 0,152947 |

Las 12 fracciones redondeadas coinciden con la Tabla 1 de [A]. N₀ es el número de soluciones del problema de las n reinas, y para n = 4, 5 y 6 el contador da lo mismo que la fuerza bruta.

## 2. Perfil variacional (Sección 9 de [P] y último valor de la Tabla 1 de [A])

`perfil/perfil.py` toma una malla de N × N casillas y, para cada masa 1 − x, resuelve el sistema de Euler reescalado en su versión discreta:

- ν_ij = b r(i) r(j) u(i − j) u(i + j − N + 1);
- r = 1/(x + bK) y u = 1/(1 + bJ);
- la media de r vale 1.

Aquí K y J son las integrales de fibra discretas, y **F(x) = (1 − x)/b**.

Los perfiles de filas y de columnas coinciden por simetría, y también los de las dos familias de diagonales. El programa usa esa simetría y al final comprueba las cuatro ecuaciones sin ella; todos los residuos quedan por debajo de 10⁻¹³.

El sistema se resuelve por iteración de punto fijo en doble precisión. F′(0) se obtiene con una diferencia centrada de paso 10⁻⁴, ∫₀¹ log F con la regla de Simpson de paso 0,005, y s_N es el valor del funcional discreto en masa uno.

```
python3 perfil/perfil.py    # mallas de 64, 128 y 256 puntos
```

| Malla | F(0) | F′(0) | ∫₀¹ log F | s_N en masa uno |
|---|---|---|---|---|
| 64 | 0,08086911055 | 0,5157082143 | −0,9441617323 | −1,9441617315 |
| 128 | 0,08090292111 | 0,5157547443 | −0,9440411935 | −1,9440411927 |
| 256 | 0,08091137290 | 0,5157663798 | −0,9440110655 | −1,9440110647 |
| Richardson (128, 256) | 0,0809141902 | 0,51577026 | −0,944001023 | −1,944001022 |

Comparación con [P]:

- **F(0)** coincide en las once cifras decimales en las tres mallas.
- **F′(0)** difiere en 5·10⁻⁹, porque [P] la calcula por diferenciación implícita.
- **∫₀¹ log F** difiere en 7·10⁻⁹, por la cuadratura en x.
- **Las extrapolaciones** coinciden con las que imprime [P]: 0,080914; 0,515770; −0,944001.

Dos comprobaciones más:

- **La identidad ∫₀¹ log F = 1 + s*(1) de [P].** El valor del funcional en masa uno, calculado directamente, coincide con ∫₀¹ log F − 1 hasta 10⁻⁹.
- **El intervalo de Nobel, Agrawal y Boyd.** El valor extrapolado, −s*(1) ≈ 1,944001022, cae dentro del intervalo [1,944000752; 1,944001082] que dan P. Nobel, A. Agrawal y S. Boyd en *Computing tighter bounds on the n-queens constant via Newton's method*, Optim. Lett. 17 (2023) 1229–1240.

Son valores numéricos: no llevan una cota de error certificada.

## Licencia

[MIT](LICENSE).
