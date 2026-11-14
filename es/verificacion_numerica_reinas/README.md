# Verificación numérica: últimas capas y perfil variacional del problema de las n reinas

Programas que reproducen los valores numéricos de dos artículos:

- **[A]** D. S. Mejia, *Completabilidad de las últimas capas en el problema de las n reinas*, 2026.
- **[P]** D. S. Mejia, *El perfil variacional de las colocaciones parciales de reinas: existencia, ecuaciones de Euler y regularidad hasta el tablero lleno*, 2026.

Los programas son deterministas: no hacen muestreo ni necesitan GPU. Se escribieron a partir de las definiciones de los artículos.

> **English summary.** This repository reproduces the numerical values reported in [A] and [P]:
> - exhaustive counts of partial n-queens placements, which give Table 1 of [A];
> - a deterministic solver for the discretized Euler equations of the variational profile, which gives the values in Section 9 of [P].
>
> Run `bash reproducir.sh` (C compiler, Python 3 and NumPy required).

## Reproducir todo

```
bash reproducir.sh          # todo, unos 5 minutos
bash reproducir.sh rapido   # comprobación corta, unos 15 segundos
```

Requisitos: un compilador de C (gcc o clang) y Python 3 con NumPy (`pip install -r requirements.txt`).

## 1. Conteos exactos (Tabla 1 de [A])

- **N_h(n)** es el número de colocaciones de n − h reinas que no se atacan en el tablero n × n.
- **f_h(n) = (n − h + 1) N_{h−1}(n) / (h² N_h(n))** es la ecuación (1.1) de [A].

| Archivo | Contenido |
|---|---|
| `conteos/conteo_capas.c` | Enumeración completa por filas con máscaras de bits. Uso: `./conteo_capas n h` |
| `conteos/fuerza_bruta.py` | Control independiente para n = 4, 5, 6: prueba todos los subconjuntos de casillas |
| `conteos/tabla1.py` | Fracciones exactas f_h(n), redondeo a seis decimales y comparación con la Tabla 1 de [A] |
| `conteos/conteos.txt` | N₀, N₁ y N₂ para n = 10, …, 15 |
| `conteos/tabla1_salida.csv` | Salida de `tabla1.py` |

| n | N₀ | N₁ | N₂ | f₁(n) | f₂(n) |
|---|---|---|---|---|---|
| 10 | 724 | 56 832 | 636 524 | 905/7104 ≈ 0,127393 | 31968/159131 ≈ 0,200891 |
| 11 | 2 680 | 264 544 | 3 634 976 | 3685/33068 ≈ 0,111437 | 41335/227186 ≈ 0,181943 |
| 12 | 14 200 | 1 375 368 | 22 146 756 | 7100/57307 ≈ 0,123894 | 630377/3691126 ≈ 0,170782 |
| 13 | 73 712 | 7 828 784 | 141 952 580 | 59891/489299 ≈ 0,122402 | 5871588/35488145 ≈ 0,165452 |
| 14 | 365 596 | 46 788 624 | 960 059 292 | 639793/5848578 ≈ 0,109393 | 745407/4706173 ≈ 0,158389 |
| 15 | 2 279 184 | 299 307 416 | 6 849 285 616 | 4273470/37413427 ≈ 0,114223 | 261893989/1712321404 ≈ 0,152947 |

- **Comprobaciones.**
  - Las 12 entradas de la Tabla 1 de [A] coinciden con estos redondeos.
  - N₀ es el número de soluciones de n reinas.
  - Para n = 4, 5, 6 el contador coincide con la fuerza bruta.
- **Tiempo.** El cálculo completo tarda unos 3,5 minutos en un núcleo; casi todo se va en n = 15, h = 2.

## 2. Perfil variacional (Sección 9 de [P]; último valor de la Tabla 1 de [A])

**Discretización.** Se toma una malla de N × N casillas. Para cada masa 1 − x se resuelve el sistema de Euler reescalado discreto:

- ν_ij = b r(i) r(j) u(i − j) u(i + j − N + 1);
- r = 1/(x + bK) y u = 1/(1 + bJ);
- la media de r vale 1.

Aquí K y J son las integrales de fibra discretas, y **F(x) = (1 − x)/b**.

Por simetría coinciden los perfiles de filas y columnas, y también las dos familias diagonales. El programa impone esa simetría y comprueba al final las cuatro ecuaciones sin ella. Todos los residuos son menores que 10⁻¹³.

**Método.**
- **Sistema:** iteración de punto fijo en doble precisión.
- **F′(0):** diferencia centrada con paso 10⁻⁴.
- **∫₀¹ log F:** regla de Simpson con paso 0,005.
- **s_N:** valor discreto del funcional en masa uno.

```
python3 perfil/perfil.py               # mallas 64, 128 y 256, alrededor de un minuto
```

| Malla | F(0) | F′(0) | ∫₀¹ log F | s_N en masa uno |
|---|---|---|---|---|
| 64 | 0,08086911055 | 0,5157082143 | −0,9441617323 | −1,9441617315 |
| 128 | 0,08090292111 | 0,5157547443 | −0,9440411935 | −1,9440411927 |
| 256 | 0,08091137290 | 0,5157663798 | −0,9440110655 | −1,9440110647 |
| Richardson (128, 256) | 0,0809141902 | 0,51577026 | −0,944001023 | −1,944001022 |

**Comparación con [P].**
- **F(0).** Coincide en las once cifras decimales en las tres mallas.
- **F′(0).** Difiere en 5·10⁻⁹, porque [P] usa diferenciación implícita.
- **∫₀¹ log F.** Difiere en 7·10⁻⁹, por la cuadratura en x.
- **Extrapolaciones.** Coinciden con las que imprime [P]: 0,080914; 0,515770; −0,944001.

**Dos comprobaciones más.**
- **Identidad ∫₀¹ log F = 1 + s*(1).** El valor del funcional en masa uno, calculado directamente, coincide con ∫₀¹ log F − 1 hasta 10⁻⁹, como afirma esa identidad de [P].
- **Comparación con Nobel–Agrawal–Boyd.** El valor extrapolado, −s*(1) ≈ 1,944001022, está dentro del intervalo [1,944000752; 1,944001082] de P. Nobel, A. Agrawal y S. Boyd, *Computing tighter bounds on the n-queens constant via Newton's method*, Optim. Lett. 17 (2023) 1229–1240.

Son valores numéricos, sin certificación por intervalos.

## Entorno de referencia

Python 3.13, NumPy 2.5 y GCC 13 en Linux.

## Licencia

MIT (ver `LICENSE`).
