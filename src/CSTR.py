import numpy as np
import scipy.optimize as opt
from numba import njit
from typing import Iterable
from numpy.typing import NDArray

import ode_solve


class CSTR:
    def __init__(self, reactions: Iterable[dict], C0: dict,
                 T0, v, VR, rho, Cp, Tc=300, UA=500, initial=None):
        self.C0 = np.array(list(C0.values()), float)
        self.components = list(C0.keys())
        self.reactions = {'coeff': [], 'exp': [], 'k0': [], 'Ea': [], 'dH': []}
        shape = (1, len(self.components))
        for r in reactions:
            r_mat = np.zeros((2, len(self.components)))
            for i, name in enumerate(self.components):
                r_mat[:, i] = r[name] if name in r else [0, 0]
            idx_base_comp = self.components.index(list(r.keys())[0])
            self.reactions['coeff'].append((r_mat[0] / abs(r_mat[0, idx_base_comp])).reshape(shape))
            self.reactions['exp'].append(r_mat[1].reshape(shape))
            self.reactions['k0'].append(np.full(shape, r['k0']))
            self.reactions['Ea'].append(np.full(shape, r['Ea']))
            dH = np.zeros(shape)
            dH[0, idx_base_comp] = r['dH']
            self.reactions['dH'].append(dH)
        for key, value in self.reactions.items():
            self.reactions[key] = np.array(value)  # shape (n_reac, 1, n_comp)
        self.reac_args = tuple(self.reactions.values())

        self.VR = float(VR)
        self.v = float(v)
        self.T0 = float(T0)
        self.Cpv = float(rho * Cp)
        self.Tc = float(Tc)
        self.UA = float(UA)
        self.initial = initial if initial is not None else {}
        self.mgrid = None
        self.C_2d = None
        self.Q_gen_2d = None

    def C_balance_eql(self, C, T, sv):
        return self.VR * sv * (self.C0 - C) + self.VR * rate(C, T, *self.reac_args[:-1]).sum(axis=0)

    def solve_steady_state(self):
        T = np.linspace(270, 1000)
        k0, Ea = self.reactions['k0'][0, 0, 0], self.reactions['Ea'][0, 0, 0]
        k = k0 * np.exp(-Ea / (8.31446261815324e-3 * T))
        T_max = T[np.argmax(k > 1000)]

        T = np.linspace(270, T_max, 100)
        sv = np.linspace(0, 2, 20)
        C = np.zeros((T.size * sv.size, self.C0.size))
        x0 = self.C0.reshape((1, -1)).repeat(T.size, axis=0)
        C[:T.size] = np.zeros_like(x0)

        for i, svi in enumerate(sv[1:], start=1):
            sol = opt.root(self.C_balance_eql, x0, args=(T, svi), method='krylov')
            C[i * T.size: (i + 1) * T.size] = sol.x

        T_repeated = T.reshape((1, -1)).repeat(sv.size, axis=0).reshape(-1)
        Q_gen = self.Q_gen(C, T_repeated)
        self.C_2d = C.reshape((self.C0.size, sv.size, T.size))
        self.Q_gen_2d = Q_gen.reshape((sv.size, T.size))
        self.mgrid = np.meshgrid(sv, T, indexing='ij')

        return T_max

    def get_2d_data(self, name):
        sv = self.mgrid[0][:, :1]
        T = self.mgrid[1][0]
        match name:
            case 'Tc':
                return (sv * self.VR * self.Cpv * (T - self.T0) - self.Q_gen_2d) / self.UA + T
            case 'Q_gen':
                return self.Q_gen_2d
            case _:
                ...  # TODO add contour plot for concentration

    def Q_gen_steady(self, T, sv):
        x0 = self.C0.reshape((1, -1)).repeat(T.size, axis=0)
        if sv > 0:
            sol = opt.root(self.C_balance_eql, x0, args=(T, sv), method='krylov')
            C = sol.x
        else:
            C = np.zeros_like(x0)
        return self.Q_gen(C, T)

    def Q_rem(self, T):
        return self.v * self.Cpv * (T - self.T0) + self.UA * (T - self.Tc)

    def Q_gen(self, C, T):
        return Q_gen(C, T, self.VR, *self.reac_args)

    def step(self, Ci, Ti, t_max, tol, method='RK45'):
        args = self.v, self.C0, self.T0, self.VR, self.Cpv, self.Tc, self.UA, *self.reac_args
        C = np.array([Ci])
        T = np.array([Ti])
        t = np.zeros(1)
        n_step = 0
        dt = 1.0
        match method:
            case 'RK23':
                tableau = ode_solve.RK23
                order = 3
            case 'RK45':
                tableau = ode_solve.RK45
                order = 5
            case _:
                raise Exception(f'Unknown method {method}.')

        while (t_remain := t_max - t[-1]) > 0 and n_step < 200:
            dt = min(dt, t_remain)
            y0 = np.concatenate((C[-1], T[-1:]))
            y, dt, dt_new = ode_solve.adaptive_RK(dydt, y0, args, dt, tol, tableau, order)
            C_new, T_new = y[:-1], y[-1:]
            C = np.append(C, C_new.reshape((1, -1)), axis=0)
            T = np.append(T, T_new)
            t = np.append(t, t[-1] + dt)
            dt = dt_new
            n_step += 1

        return C[1:], T[1:], t[1:]


@njit(cache=True)
def rate(C, T, coeff, exp, k0, Ea) -> NDArray:  # shape (n_reac, -1, n_comp)
    C = C.reshape((1, -1, C.shape[-1]))
    T = T.reshape((1, -1, 1))
    f_T = k0 * np.exp(-Ea / (8.31446261815324e-3 * T))
    C_exp = C ** exp
    f_C = np.zeros_like(C_exp)
    for i in range(C_exp.shape[0]):
        for j in range(C_exp.shape[1]):
            f_C[i, j] = np.cumprod(C_exp[i, j])
    f_C = f_C[:, :, -1:]
    return coeff * f_T * f_C


@njit(cache=True)
def Q_gen(C, T, VR, coeff, exp, k0, Ea, dH) -> NDArray:  # shape (-1,)
    return (dH * VR * rate(C, T, coeff, exp, k0, Ea)).sum(axis=0).sum(axis=1)


@njit(cache=True)
def dydt(y, v, C0, T0, VR, Cpv, Tc, UA, coeff, exp, k0, Ea, dH) -> NDArray:  # shape (n_comp + 1,)
    C, T = y[:-1], y[-1:]
    dCdt = (C0 - C) * v / VR + rate(C, T, coeff, exp, k0, Ea).sum(axis=0).reshape(-1)
    dTdt = Q_gen(C, T, VR, coeff, exp, k0, Ea, dH) / (Cpv * VR) + (T0 - T) * v / VR + UA * (Tc - T) / (Cpv * VR)
    return np.concatenate((dCdt, dTdt))
