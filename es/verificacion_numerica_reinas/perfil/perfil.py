"""Verificación numérica independiente del perfil variacional.

Reproduce los valores numéricos de [P] (tabla de la Sección 9: F(0), F'(0) e integral de log F)
y el valor F(0) de la Tabla 1 de [A]. Las referencias [A] y [P] están en el README.

Discretización: malla de N x N casillas. En masa 1 - x se resuelve el sistema de Euler reescalado
    nu_ij = b r(i) r(j) u(i-j) u(i+j-N+1),
    r = 1/(x + b K),   u = 1/(1 + b J),   media(r) = 1,
con K(i) = (1/N) sum_j r(j) u(i-j) u(i+j-N+1) y J(k) = m(k)/(b u(k)), donde m(k) es (1/N) por la
suma de nu sobre la diagonal k (las dos familias diagonales coinciden por simetría). F(x) = (1-x)/b.
Por la simetría del problema, el perfil de filas y el de columnas coinciden, y también las dos
familias diagonales; el programa lo impone y comprueba al final que las cuatro ecuaciones se cumplen
sin esa simetría. En x = 0 calcula también el valor discreto del funcional s_N.

Iteración de punto fijo determinista en doble precisión: sin muestreo ni GPU.
Uso:  python perfil.py            (tabla completa, mallas 64, 128 y 256)
      python perfil.py 64 128     (sólo esas mallas)
"""
import sys
import numpy as np

def resolver(N, x, tol=1e-13, maxit=200000, ini=None):
    i = np.arange(N)
    D = (i[:, None] - i[None, :]) + (N - 1)            # diagonal i - j      -> 0..2N-2
    A = (i[:, None] + i[None, :])                      # antidiagonal i + j  -> 0..2N-2
    if ini is None:
        r, u, b = np.ones(N), np.full(2 * N - 1, 0.5), 4.0
    else:
        r, u, b = ini[0].copy(), ini[1].copy(), ini[2]
    for it in range(maxit):
        P = np.outer(r, r) * u[D] * u[A]                # nu / b
        K = P.sum(axis=1) / N / r
        m = b * np.bincount(D.ravel(), weights=P.ravel(), minlength=2 * N - 1) / N
        res = max(np.abs(r * (x + b * K) - 1).max(), np.abs(u + m - 1).max(), abs(r.mean() - 1))
        if res < tol:
            break
        u = u / (u + m)                                 # diagonales: u = 1/(1 + bJ)
        P = np.outer(r, r) * u[D] * u[A]
        K = P.sum(axis=1) / N / r
        r = np.sqrt(r / (x + b * K))                    # filas: media geométrica, evita oscilar
        c = r.mean()
        r, b = r / c, b * c * c                         # media(r) = 1
    nu = b * np.outer(r, r) * u[D] * u[A]
    # comprobación sin simetría: filas, columnas y las dos familias diagonales
    fil, col = nu.sum(axis=1) / N, nu.sum(axis=0) / N
    mD = np.bincount(D.ravel(), weights=nu.ravel(), minlength=2 * N - 1) / N
    mA = np.bincount(A.ravel(), weights=nu.ravel(), minlength=2 * N - 1) / N
    chk = max(np.abs(fil - (1 - x * r)).max(), np.abs(col - (1 - x * r)).max(),
              np.abs(mD - (1 - u)).max(), np.abs(mA - (1 - u)).max())
    s = -3.0 - (nu * np.log(nu)).sum() / N**2 - 2.0 * (u * np.log(u)).sum() / N
    return dict(F=(1 - x) / b, b=b, s=s, it=it, res=res, chk=chk, masa=nu.sum() / N**2,
                estado=(r, u, b))

def tabla(N, paso=0.005, h=1e-4):
    R0 = resolver(N, 0.0)
    # F'(0): diferencia centrada (el sistema discreto es analítico en x alrededor de 0)
    Fp = (resolver(N, h, ini=R0['estado'])['F'] - resolver(N, -h, ini=R0['estado'])['F']) / (2 * h)
    # integral de log F en [0,1] por Simpson; F(1) = 1 (masa nula), log F(1) = 0
    xs = np.linspace(0.0, 1.0, int(round(1 / paso)) + 1)
    logF, est, peor = [], R0['estado'], 0.0
    for x in xs:
        if x == 1.0:
            logF.append(0.0); continue
        R = resolver(N, x, ini=est)
        est = R['estado']; peor = max(peor, R['res'], R['chk'])
        logF.append(np.log(R['F']))
    logF = np.array(logF)
    hx = xs[1] - xs[0]
    integral = hx / 3 * (logF[0] + logF[-1] + 4 * logF[1:-1:2].sum() + 2 * logF[2:-1:2].sum())
    return dict(N=N, F0=R0['F'], Fp0=Fp, intlogF=integral, s=R0['s'], masa=R0['masa'],
                res=max(R0['res'], R0['chk'], peor))

if __name__ == "__main__":
    mallas = [int(a) for a in sys.argv[1:]] or [64, 128, 256]
    filas = []
    for N in mallas:
        T = tabla(N)
        filas.append(T)
        print(f"N={N:4d}  F(0)={T['F0']:.11f}  F'(0)={T['Fp0']:.10f}  int log F={T['intlogF']:.10f}"
              f"  s_N(1)={T['s']:.10f}  residuo max={T['res']:.1e}")
    if len(filas) >= 2:
        a, c = filas[-2], filas[-1]
        rich = lambda k: c[k] + (c[k] - a[k]) / 3
        print(f"Richardson ({a['N']},{c['N']}):  F(0)={rich('F0'):.10f}  F'(0)={rich('Fp0'):.8f}"
              f"  int log F={rich('intlogF'):.9f}  s_N(1)={rich('s'):.9f}")
