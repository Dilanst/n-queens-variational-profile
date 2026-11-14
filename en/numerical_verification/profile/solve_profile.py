"""Discrete variational profile.

Reproduces the values in the table of Section 9 of [P] (F(0), F'(0) and the integral of log F)
and the value of F(0) in Table 1 of [A]. The references [A] and [P] are in the README.

Discretization: a grid of N x N squares. At mass 1 - x it solves the rescaled Euler system
    nu_ij = b r(i) r(j) u(i-j) u(i+j-N+1),
    r = 1/(x + b K),   u = 1/(1 + b J),   mean(r) = 1,
with K(i) = (1/N) sum_j r(j) u(i-j) u(i+j-N+1) and J(k) = m(k)/(b u(k)), where m(k) is 1/N times
the sum of nu over the diagonal k (the two diagonal families coincide by symmetry). F(x) = (1-x)/b.
By the symmetry of the problem, the row profile and the column profile coincide, and so do the two
diagonal families; the program imposes this and checks at the end that the four equations hold
without that symmetry. At x = 0 it also computes the discrete value s_N of the functional.

The system is solved by fixed-point iteration in double precision.
Usage:  python solve_profile.py            (grids of 64, 128 and 256 points)
        python solve_profile.py 64 128     (only those grids)
"""
import sys
import numpy as np

def solve(N, x, tol=1e-13, maxit=200000, start=None):
    i = np.arange(N)
    D = (i[:, None] - i[None, :]) + (N - 1)            # diagonal i - j      -> 0..2N-2
    A = (i[:, None] + i[None, :])                      # antidiagonal i + j  -> 0..2N-2
    if start is None:
        r, u, b = np.ones(N), np.full(2 * N - 1, 0.5), 4.0
    else:
        r, u, b = start[0].copy(), start[1].copy(), start[2]
    for it in range(maxit):
        P = np.outer(r, r) * u[D] * u[A]                # nu / b
        K = P.sum(axis=1) / N / r
        m = b * np.bincount(D.ravel(), weights=P.ravel(), minlength=2 * N - 1) / N
        res = max(np.abs(r * (x + b * K) - 1).max(), np.abs(u + m - 1).max(), abs(r.mean() - 1))
        if res < tol:
            break
        u = u / (u + m)                                 # diagonals: u = 1/(1 + bJ)
        P = np.outer(r, r) * u[D] * u[A]
        K = P.sum(axis=1) / N / r
        r = np.sqrt(r / (x + b * K))                    # rows: geometric mean, avoids oscillation
        c = r.mean()
        r, b = r / c, b * c * c                         # mean(r) = 1
    nu = b * np.outer(r, r) * u[D] * u[A]
    # check without symmetry: rows, columns and the two diagonal families
    row_mass, col_mass = nu.sum(axis=1) / N, nu.sum(axis=0) / N
    mD = np.bincount(D.ravel(), weights=nu.ravel(), minlength=2 * N - 1) / N
    mA = np.bincount(A.ravel(), weights=nu.ravel(), minlength=2 * N - 1) / N
    chk = max(np.abs(row_mass - (1 - x * r)).max(), np.abs(col_mass - (1 - x * r)).max(),
              np.abs(mD - (1 - u)).max(), np.abs(mA - (1 - u)).max())
    s = -3.0 - (nu * np.log(nu)).sum() / N**2 - 2.0 * (u * np.log(u)).sum() / N
    return dict(F=(1 - x) / b, b=b, s=s, it=it, res=res, chk=chk, mass=nu.sum() / N**2,
                state=(r, u, b))

def table(N, step=0.005, h=1e-4):
    R0 = solve(N, 0.0)
    # F'(0): centered difference (the discrete system is analytic in x around 0)
    Fp = (solve(N, h, start=R0['state'])['F'] - solve(N, -h, start=R0['state'])['F']) / (2 * h)
    # integral of log F over [0,1] by Simpson's rule; F(1) = 1 (zero mass), log F(1) = 0
    xs = np.linspace(0.0, 1.0, int(round(1 / step)) + 1)
    logF, st, worst = [], R0['state'], 0.0
    for x in xs:
        if x == 1.0:
            logF.append(0.0); continue
        R = solve(N, x, start=st)
        st = R['state']; worst = max(worst, R['res'], R['chk'])
        logF.append(np.log(R['F']))
    logF = np.array(logF)
    hx = xs[1] - xs[0]
    integral = hx / 3 * (logF[0] + logF[-1] + 4 * logF[1:-1:2].sum() + 2 * logF[2:-1:2].sum())
    return dict(N=N, F0=R0['F'], Fp0=Fp, intlogF=integral, s=R0['s'], mass=R0['mass'],
                res=max(R0['res'], R0['chk'], worst))

if __name__ == "__main__":
    grids = [int(a) for a in sys.argv[1:]] or [64, 128, 256]
    results = []
    for N in grids:
        T = table(N)
        results.append(T)
        print(f"N={N:4d}  F(0)={T['F0']:.11f}  F'(0)={T['Fp0']:.10f}  int log F={T['intlogF']:.10f}"
              f"  s_N(1)={T['s']:.10f}  max residual={T['res']:.1e}")
    if len(results) >= 2:
        a, c = results[-2], results[-1]
        rich = lambda k: c[k] + (c[k] - a[k]) / 3
        print(f"Richardson ({a['N']},{c['N']}):  F(0)={rich('F0'):.10f}  F'(0)={rich('Fp0'):.8f}"
              f"  int log F={rich('intlogF'):.9f}  s_N(1)={rich('s'):.9f}")
