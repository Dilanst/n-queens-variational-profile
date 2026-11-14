# Numerical verification

These programs reproduce the numerical values of two papers:

- **[A]** D. S. Mejia, *Completability of the final layers in the n-queens problem*, 2026. Repository: [n-queens-last-layers](https://github.com/Dilanst/n-queens-last-layers).
- **[P]** D. S. Mejia, *The variational profile of partial queen placements: existence, Euler equations, and regularity up to the full board*, 2026. Repository: [n-queens-variational-profile](https://github.com/Dilanst/n-queens-variational-profile).

## How to reproduce them

You need a C compiler (gcc or clang) and Python 3 with NumPy (`pip install -r requirements.txt`). From this folder:

```
bash reproduce.sh          # everything, about four and a half minutes
bash reproduce.sh quick    # short check, about 15 seconds
```

The script compiles the counter and checks it against a brute-force search. It then recomputes the counts and compares them with the stored ones; if anything differs, it says so and stops. Finally, it solves for the profile and prints its values.

## Machine and running times

The times in this section were measured on this machine:

- Ubuntu 24.04 (Linux)
- Intel Xeon processor at 2.1 GHz, with 2 cores
- 8 GB of RAM
- GCC 13.3, Python 3.13 and NumPy 2.5

The programs do not run in parallel: each one uses a single core.

| Step | Time |
|---|---|
| Compiling the counter and comparing it with the brute-force search | 1 s |
| Counts for n = 10, …, 15 | 3 min 34 s, of which n = 15, h = 2 takes 2 min 37 s |
| Profile on the grids of 64, 128 and 256 points | 54 s, of which the 256-point grid takes 43 s |
| Total (`bash reproduce.sh`) | About 4 min 30 s |

The results are byte-for-byte identical to the stored files.

## 1. Exact counts (Table 1 of [A])

- **N_h(n)** is the number of placements of n − h non-attacking queens on the n × n board.
- **f_h(n) = (n − h + 1) N_{h−1}(n) / (h² N_h(n))** is the fraction in equation (1.1) of [A].

| File | Contents |
|---|---|
| `counts/layer_count.c` | Counter: goes through all the placements row by row. Usage: `./layer_count n h` |
| `counts/brute_force.py` | Check for n = 4, 5, 6: tries every subset of squares |
| `counts/table1.py` | Computes f_h(n) as an exact fraction, rounds it to six decimal places and compares it with Table 1 of [A] |
| `counts/counts.txt` | N₀, N₁ and N₂ for n = 10, …, 15 |
| `counts/table1_output.csv` | What `table1.py` writes |

| n | N₀ | N₁ | N₂ | f₁(n) | f₂(n) |
|---|---|---|---|---|---|
| 10 | 724 | 56,832 | 636,524 | 905/7104 ≈ 0.127393 | 31968/159131 ≈ 0.200891 |
| 11 | 2,680 | 264,544 | 3,634,976 | 3685/33068 ≈ 0.111437 | 41335/227186 ≈ 0.181943 |
| 12 | 14,200 | 1,375,368 | 22,146,756 | 7100/57307 ≈ 0.123894 | 630377/3691126 ≈ 0.170782 |
| 13 | 73,712 | 7,828,784 | 141,952,580 | 59891/489299 ≈ 0.122402 | 5871588/35488145 ≈ 0.165452 |
| 14 | 365,596 | 46,788,624 | 960,059,292 | 639793/5848578 ≈ 0.109393 | 745407/4706173 ≈ 0.158389 |
| 15 | 2,279,184 | 299,307,416 | 6,849,285,616 | 4273470/37413427 ≈ 0.114223 | 261893989/1712321404 ≈ 0.152947 |

The 12 rounded fractions agree with Table 1 of [A]. N₀ is the number of solutions of the n-queens problem, and for n = 4, 5 and 6 the counter gives the same values as the brute-force search.

## 2. Variational profile (Section 9 of [P] and the last value in Table 1 of [A])

`profile/solve_profile.py` takes a grid of N × N squares and, for each mass 1 − x, solves the discrete version of the rescaled Euler system:

- ν_ij = b r(i) r(j) u(i − j) u(i + j − N + 1);
- r = 1/(x + bK) and u = 1/(1 + bJ);
- the mean of r is 1.

Here K and J are the discrete fiber integrals, and **F(x) = (1 − x)/b**.

The row and column profiles coincide by symmetry, and so do those of the two diagonal families. The program uses this symmetry and at the end checks the four equations without it; all the residuals are below 10⁻¹³.

The system is solved by fixed-point iteration in double precision. F′(0) is obtained by a centered difference with step 10⁻⁴, ∫₀¹ log F by Simpson's rule with step 0.005, and s_N is the value of the discrete functional at mass one.

```
python3 profile/solve_profile.py    # grids of 64, 128 and 256 points
```

| Grid | F(0) | F′(0) | ∫₀¹ log F | s_N at mass one |
|---|---|---|---|---|
| 64 | 0.08086911055 | 0.5157082143 | −0.9441617323 | −1.9441617315 |
| 128 | 0.08090292111 | 0.5157547443 | −0.9440411935 | −1.9440411927 |
| 256 | 0.08091137290 | 0.5157663798 | −0.9440110655 | −1.9440110647 |
| Richardson (128, 256) | 0.0809141902 | 0.51577026 | −0.944001023 | −1.944001022 |

Comparison with [P]:

- **F(0)** agrees with [P] to all eleven decimal places on the three grids.
- **F′(0)** differs by 5·10⁻⁹, because [P] computes it by implicit differentiation.
- **∫₀¹ log F** differs by 7·10⁻⁹, because of the quadrature in x.
- **The extrapolations** agree with those printed in [P]: 0.080914, 0.515770, −0.944001.

Two further checks:

- **The identity ∫₀¹ log F = 1 + s*(1) of [P].** The value of the functional at mass one, computed directly, agrees with ∫₀¹ log F − 1 up to 10⁻⁹.
- **The interval of Nobel, Agrawal and Boyd.** The extrapolated value, −s*(1) ≈ 1.944001022, lies inside the interval [1.944000752, 1.944001082] given by P. Nobel, A. Agrawal and S. Boyd in *Computing tighter bounds on the n-queens constant via Newton's method*, Optim. Lett. 17 (2023) 1229–1240.

These are numerical values: they carry no certified error bound.

## License

[MIT](LICENSE).
