import os

import numpy as np
import matplotlib.pyplot as plt
from ode_solve import adaptive_RK

"""first-order reaction A→P"""
k1 = lambda T: 1.2e9 * np.exp(-72.7475 / 8.314e-3 / T)  # 1/s
dH1 = -5e7  # J/kmol

VR = 0.1  # m3
CA0 = 1.0  # kmol/m3
CP0 = 0.0
v = 1.667e-3  # m3/s
T0 = 340  # K
rho = 1000  # kg/m3
Cp = 239  # J/(kg*K)
Tc = 305  # K
UA = 900  # J/(s*K)


def Q_gen(CA, T):
    return -dH1 * k1(T) * CA * VR


def Q_rem(T, T0):
    return (T - T0) * v * rho * Cp + UA * (T - Tc)


def dydt(y, T0):
    CA, CP, T = y
    kCA = (CA0 - CA) * v / VR - k1(T) * CA
    kCP = (CP0 - CP) * v / VR + k1(T) * CA
    kT = (Q_gen(CA, T) - Q_rem(T, T0)) / (VR * rho * Cp)
    return np.array([kCA, kCP, kT])


t_max = 1000
dt = 1.0

"""initial values"""
t = np.zeros(1)
CA = np.array([0.9])
CP = np.array([0.1])
T = np.array([322.5])
T0_arr = np.array([340.])
Q_gen_transient = Q_gen(CA[-1], T[-1])
Q_rem_transient = Q_rem(T[-1], T0_arr[-1])

while t[-1] < t_max:
    y0 = np.array([CA[-1], CP[-1], T[-1]])
    y, dt, dt_new = adaptive_RK(dydt, y0, args=(T0_arr[-1],), dt=dt, tol=1e-8)
    CA_new, CP_new, T_new = y

    t = np.append(t, t[-1] + dt)
    CA = np.append(CA, CA_new)
    CP = np.append(CP, CP_new)
    T = np.append(T, T_new)
    Q_gen_transient = np.append(Q_gen_transient, Q_gen(CA_new, T_new))
    Q_rem_transient = np.append(Q_rem_transient, Q_rem(T_new, T0_arr[-1]))
    T0_arr = np.append(T0_arr, T0_arr[-1] + 0.5 / (10 / dt))  # T0 increases 0.5 K every 10 seconds

    print(f'\rt = {t[-1]:.2f} s, CA = {CA[-1]:.3f} kmol/m3, CP = {CP[-1]:.3f} kmol/m3, '
          f'T = {T[-1]:.2f} K, dt = {dt:.2e} s', end='')
    dt = dt_new

"""plots"""

text_color = '#333333'
plt.rcParams.update({
    'text.color': text_color,
    'axes.labelcolor': text_color,
    'axes.edgecolor': text_color,
    'xtick.color': text_color,
    'ytick.color': text_color,
})

fig_steady = plt.figure(1, figsize=(6.4, 4.8), layout='tight')
fig_transient = plt.figure(2, figsize=(6.4, 7.2), layout='tight')
ax_QT = fig_steady.add_subplot()
ax_Ct = fig_transient.add_subplot(311)
ax_Tt = fig_transient.add_subplot(312)
ax_Qt = fig_transient.add_subplot(313)

Tx = np.linspace(250, 1000, 100)
Q_gen_steady = -dH1 * VR * k1(Tx) * CA0 / (1 + k1(Tx) * VR / v)
Q_rem_steady = Q_rem(Tx, T0_arr[-1])

ax_QT.plot(Tx, Q_gen_steady, color='#FF731D', label='$Q_{gen}$')
ax_QT.plot(Tx, Q_rem_steady, color='#5F9DF7', label='$Q_{rem}$')
ax_QT.plot(T, Q_gen_transient, ':', color='darkcyan', label='$Q_{gen}$ (transient)')
ax_QT.legend(loc='lower right')
ax_QT.set(xlabel='T (K)', ylabel='Q (J/s)',
          xlim=(250, 500), ylim=(0, 1.75e5))
ax_QT.ticklabel_format(axis='y', scilimits=[-3, 3])

ax_Tt.plot(t, T, color='#F45050', label='$T_r$')
ax_Tt.plot(t, T0_arr, color='#0079FF', label='$T_0$')
ax_Tt.legend(loc='upper left')
ax_Tt.set(xlabel='t (s)', ylabel='T (K)',
          xlim=(0, t_max))

ax_Ct.plot(t, CA, color='#FD8A8A', label='$c_A$')
ax_Ct.plot(t, CP, color='#9EA1D4', label='$c_P$')
ax_Ct.legend(loc='upper left')
ax_Ct.set(xlabel='t (s)', ylabel='Concentration (kmol/m$^3$)',
          xlim=(0, t_max))

ax_Qt.semilogy(t, Q_gen_transient, color='#FF731D', label='$Q_{gen}$')
ax_Qt.semilogy(t, Q_rem_transient, color='#5F9DF7', label='$Q_{rem}$')
ax_Qt.legend(loc='upper left')
ax_Qt.set(xlabel='t (s)', ylabel='Q (J/s)')
ax_Qt.ticklabel_format(axis='x', useOffset=False)

DIR = 'first_order'
os.makedirs(DIR, exist_ok=True)
for i, name in enumerate(['QT', 'CTQ-t'], start=1):
    plt.figure(i)
    plt.savefig(f'{DIR}/{name}.jpg', dpi=600, bbox_inches='tight')

plt.show()
