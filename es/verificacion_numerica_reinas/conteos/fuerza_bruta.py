"""Control del contador: cuenta N_h(n) probando todos los subconjuntos de casillas (n <= 6)."""
import itertools
def contar(n, q):
    casillas = [(r, c) for r in range(n) for c in range(n)]
    total = 0
    for S in itertools.combinations(casillas, q):
        if (len({r for r, c in S}) == q and len({c for r, c in S}) == q and
                len({r - c for r, c in S}) == q and len({r + c for r, c in S}) == q):
            total += 1
    return total
if __name__ == "__main__":
    for n in range(4, 7):
        for h in range(3):
            print(n, h, contar(n, n - h))
