import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
from ode_solve import adaptive_RK

"""first-order reaction A→P"""
k1 = lambda T: 8e15 * np.exp(-133.024 / 8.314e-3 / T)  # 1/s
dH1 = -2e8  # J/kmol

VR = 2.37  # m3
CA0 = 3.0  # kmol/m3
v = 0.01  # m3/s  1.905e-3
T0 = 335.15  # K
rho = 1000  # kg/m3
Cp = 4000  # J/(kg*K)
Tc = 369.15  # K
UA = 80000  # J/(s*K)


def Q_gen(CA, T):
    return -dH1 * k1(T) * CA * VR


def Q_rem(T):
    return (T - T0) * v * rho * Cp + UA * (T - Tc)


def dydt(y):
    CA, T = y
    kCA = (CA0 - CA) * v / VR - k1(T) * CA
    kT = (Q_gen(CA, T) - Q_rem(T)) / (VR * rho * Cp)
    return np.array([kCA, kT])


t_max = 2000
dt = 1.0

n = 10
T_lim = 300, 500
initials = np.array([[T0, CA0]]).repeat(n, axis=0)
initials = np.append(initials, np.zeros((n, 2)), axis=0)
initials[:, 0] = np.array([
    300, 365, 368, 371, 380, 400, 420, 440, 460, 480,
    305, 325, 345, 365, 385, 400, 408.61, 408.7, 409.5, 412,
], dtype=float)

paths = []
Q_gen_list = []

for val in initials:
    t = np.zeros(1)
    T = np.array([val[0]])
    CA = np.array([val[1]])
    Q_gen_transient = Q_gen(CA[-1], T[-1])

    while t[-1] < t_max:
        y0 = np.array([CA[-1], T[-1]])
        y, dt, dt_new = adaptive_RK(dydt, y0, dt=dt, tol=1e-10)
        CA_new, T_new = y

        t = np.append(t, t[-1] + dt)
        CA = np.append(CA, CA_new)
        T = np.append(T, T_new)
        Q_gen_transient = np.append(Q_gen_transient, Q_gen(CA_new, T_new))

        if T_new > T_lim[1] + 50: break
        print(f'\rt = {t[-1]:.2f} s, CA = {CA[-1]:.3f} kmol/m3, '
              f'T = {T[-1]:.2f} K, dt = {dt:.2e} s', end='')
        dt = dt_new

    xy = np.concatenate((T, CA)).reshape((2, -1))
    paths.append(xy)
    Q_gen_list.append(Q_gen_transient)

"""plots"""
fig_QT = plt.figure(1, layout='tight')
fig_CT = plt.figure(2, layout='tight')
ax_QT = fig_QT.add_subplot()
ax_CT = fig_CT.add_subplot()

Ts = np.linspace(*T_lim, 100)
CAs = lambda T: CA0 / (1 + k1(T) * VR / v)

ax_QT.plot(Ts, Q_gen(CAs(Ts), Ts), color='#FF731D', label='$Q_{gen}$')
ax_QT.plot(Ts, Q_rem(Ts), color='#5F9DF7', label='$Q_{rem}$')
ax_QT.plot(paths[-1][0], Q_gen_list[-1], ':', color='darkcyan', label='$Q_{gen}$ (transient)')
ax_QT.legend(loc='lower right')
ax_QT.set(xlabel='T (K)', ylabel='Q (J/s)',
          xlim=T_lim, ylim=(0, 1e7))
ax_QT.ticklabel_format(axis='y', scilimits=[-3, 3])

for xy in paths:
    T, CA = xy
    ax_CT.plot(T, CA, linewidth=0.8, color='#79B4B7')
    idx = round(len(T) / 4.3)
    ax_CT.annotate("", xy=(T[idx + 1], CA[idx + 1]), xytext=(T[idx], CA[idx]),
                   arrowprops=dict(width=0, headwidth=4, headlength=8, color='#79B4B7'))

residual = lambda T: Q_gen(CAs(T), T) - Q_rem(T)
saddle_point = opt.root(residual, x0=380).x[0]
stable_point_1 = opt.root(residual, x0=350).x[0]
stable_point_2 = opt.root(residual, x0=400).x[0]

ax_CT.scatter(saddle_point, CAs(saddle_point), marker='+', color='tab:blue', zorder=2, s=80)
ax_CT.scatter(stable_point_1, CAs(stable_point_1), marker='*', color='tab:blue', zorder=2, s=80)
ax_CT.scatter(stable_point_2, CAs(stable_point_2), marker='*', color='tab:blue', zorder=2, s=80)
ax_CT.annotate("stable points",
               xy=(stable_point_1 - 3, CAs(stable_point_1) - 0.01),
               xytext=(310, 2.0),
               fontsize='large',
               arrowprops=dict(facecolor='black',
                               width=0, headwidth=4, headlength=8,
                               connectionstyle='arc3,rad=-0.2'))
ax_CT.annotate("", xy=(stable_point_2 - 3, CAs(stable_point_2) + 0.01),
               xytext=(350, 1.95),
               arrowprops=dict(facecolor='black',
                               width=0, headwidth=4, headlength=8,
                               connectionstyle='arc3,rad=0.2'))
ax_CT.annotate("saddle point",
               xy=(saddle_point + 5, CAs(saddle_point) - 0.04),
               fontsize='large')
ax_CT.set(xlabel='T (K)', ylabel='$C_A$ (mol/L)', xlim=T_lim, ylim=(-0.1, CA0 + 0.1))

# plt.savefig(f'figures/phase_plane_2.jpg', dpi=600, bbox_inches='tight')
plt.show()
