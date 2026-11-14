/* conteo_capas.c: cuenta N_h(n), el número de colocaciones de n-h reinas que no se atacan
   en el tablero n x n ([A], ecuación (1.1) y Tabla 1).

   Recorre el tablero fila a fila con máscaras de bits: cada fila recibe una reina o queda
   vacía, hasta dejar exactamente h filas vacías. Las diagonales ocupadas siguen avanzando
   también al pasar por una fila vacía.

   Compilar:  cc -O2 -o conteo_capas conteo_capas.c
   Uso:       ./conteo_capas n h        (n <= 32)
   Salida:    n h N_h(n)
*/
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

static int n, H;
static uint64_t lleno;
static unsigned long long total;

static void fila(int r, uint64_t cols, uint64_t d1, uint64_t d2, int vacias) {
    int quedan = n - r;                          /* filas por decidir, incluida la r */
    if (H - vacias > quedan) return;             /* ya no caben las filas vacías que faltan */
    uint64_t libres = lleno & ~(cols | d1 | d2);
    if (r == n - 1) {                            /* última fila */
        if (vacias == H) total += (unsigned long long)__builtin_popcountll(libres);
        else if (vacias == H - 1) total += 1ULL;
        return;
    }
    if (vacias < H)                              /* dejar la fila r vacía */
        fila(r + 1, cols, (d1 << 1) & lleno, d2 >> 1, vacias + 1);
    if (H - vacias == quedan) return;            /* todas las restantes deben quedar vacías */
    while (libres) {                             /* poner una reina en cada casilla libre */
        uint64_t bit = libres & (~libres + 1);
        libres ^= bit;
        fila(r + 1, cols | bit, ((d1 | bit) << 1) & lleno, (d2 | bit) >> 1, vacias);
    }
}

int main(int argc, char **argv) {
    if (argc != 3) { fprintf(stderr, "uso: %s n h\n", argv[0]); return 1; }
    n = atoi(argv[1]); H = atoi(argv[2]);
    if (n < 1 || n > 32 || H < 0 || H > n) { fprintf(stderr, "fuera de rango\n"); return 1; }
    lleno = (n == 64) ? ~0ULL : ((1ULL << n) - 1);
    total = 0;
    fila(0, 0, 0, 0, 0);
    printf("%d %d %llu\n", n, H, total);
    return 0;
}
