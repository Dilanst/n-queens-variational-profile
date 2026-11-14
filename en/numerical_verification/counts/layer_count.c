/* layer_count.c: computes N_h(n), the number of placements of n-h non-attacking queens
   on the n x n board ([A], equation (1.1) and Table 1).

   It scans the board row by row with bit masks: each row receives a queen or is left
   empty, until exactly h rows are empty. The occupied diagonals keep shifting when they
   pass through an empty row as well.

   Build:   cc -O2 -o layer_count layer_count.c
   Usage:   ./layer_count n h        (n <= 32)
   Output:  n h N_h(n)
*/
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

static int n, H;
static uint64_t full;
static unsigned long long total;

static void row(int r, uint64_t cols, uint64_t d1, uint64_t d2, int empty) {
    int left = n - r;                            /* rows still to decide, row r included */
    if (H - empty > left) return;                /* the empty rows still needed no longer fit */
    uint64_t open = full & ~(cols | d1 | d2);
    if (r == n - 1) {                            /* last row */
        if (empty == H) total += (unsigned long long)__builtin_popcountll(open);
        else if (empty == H - 1) total += 1ULL;
        return;
    }
    if (empty < H)                               /* leave row r empty */
        row(r + 1, cols, (d1 << 1) & full, d2 >> 1, empty + 1);
    if (H - empty == left) return;               /* all the remaining rows must stay empty */
    while (open) {                               /* put a queen on each open square */
        uint64_t bit = open & (~open + 1);
        open ^= bit;
        row(r + 1, cols | bit, ((d1 | bit) << 1) & full, (d2 | bit) >> 1, empty);
    }
}

int main(int argc, char **argv) {
    if (argc != 3) { fprintf(stderr, "usage: %s n h\n", argv[0]); return 1; }
    n = atoi(argv[1]); H = atoi(argv[2]);
    if (n < 1 || n > 32 || H < 0 || H > n) { fprintf(stderr, "out of range\n"); return 1; }
    full = (n == 64) ? ~0ULL : ((1ULL << n) - 1);
    total = 0;
    row(0, 0, 0, 0, 0);
    printf("%d %d %llu\n", n, H, total);
    return 0;
}
