import numpy as np
from numba import njit

"""Explicit Runge-Kutta methods"""

# Bogacki-Shampine
RK23 = np.array([
    [0, 0, 0, 0],
    [1 / 2, 0, 0, 0],
    [0, 3 / 4, 0, 0],
    [2 / 9, 1 / 3, 4 / 9, 0],  # order 3
    [7 / 24, 1 / 4, 1 / 3, 1 / 8]  # order 2
])
# Dormand-Prince
RK45 = np.array([
    [0, 0, 0, 0, 0, 0, 0],
    [1 / 5, 0, 0, 0, 0, 0, 0],
    [3 / 40, 9 / 40, 0, 0, 0, 0, 0],
    [44 / 45, -56 / 15, 32 / 9, 0, 0, 0, 0],
    [19372 / 6561, -25360 / 2187, 64448 / 6561, -212 / 729, 0, 0, 0],
    [9017 / 3168, -355 / 33, 46732 / 5247, 49 / 176, -5103 / 18656, 0, 0],
    [35 / 384, 0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84, 0],  # order 5
    [5179 / 57600, 0, 7571 / 16695, 393 / 640, -92097 / 339200, 187 / 2100, 1 / 40]  # order 4
])


@njit
def adaptive_RK(fun, y0, args=(), dt=1., tol=1e-10, tableau=RK45, order=5):
    k = np.zeros((y0.size, tableau.shape[1]))
    stages = np.zeros((tableau.shape[0], y0.size))

    stages[0] = y0
    for i, coeff in enumerate(tableau[1:]):
        k[:, i] = fun(stages[i], *args)
        stages[i + 1] = y0 + dt * np.sum(coeff * k, axis=1)

    E = np.linalg.norm(stages[-2] - stages[-1], np.inf) + 1e-20
    maxerr = tol * (1 + np.linalg.norm(y0, np.inf))
    q = min(0.8 * (maxerr / E) ** (1 / (order + 1)), 4.)
    dt_new = dt * q
    if E > maxerr:
        return adaptive_RK(fun, y0, args, dt_new, tol, tableau, order)
    return stages[-2], dt, dt_new
